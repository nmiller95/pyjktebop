import os
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import corner

plt.style.use('jktebop.mplstyle')


def wrap_phase(phase_array):
    """ Convert phase array from (0,1) to (-0.5,0.5) """
    for i, ph in enumerate(phase_array):
        if ph >= 0.5:
            phase_array[i] = ph - 1
    return phase_array


def extend_array(array, phase=False):
    """ Extend/copy arrays so secondary eclipses at phase 0.5 plot correctly """
    if phase:
        return np.hstack([array-1, array, array+1])
    else:
        return np.hstack([array, array, array])


def make_bins(ph_array, mag_array, n_bin=1000, min_in_bin=8):
    """Convert phase and magnitude array into binned phase and magnitude arrays"""
    n_in_bin, bin_edges = np.histogram(ph_array, bins=n_bin, range=(-1.5, 1.5))
    bin_indices = np.digitize(ph_array, bin_edges)

    ph_bin = np.zeros(n_bin)
    mag_bin = np.zeros(n_bin)
    for i in range(len(n_in_bin)):
        if n_in_bin[i] > 0:
            ph_bin[i] = np.mean(ph_array[bin_indices == i + 1])
            mag_bin[i] = np.median(mag_array[bin_indices == i + 1])

    ph_bin = ph_bin[n_in_bin > min_in_bin]
    mag_bin = mag_bin[n_in_bin > min_in_bin]

    return ph_bin, mag_bin


class JKTEBOP:
    def __init__(self, target_name, recompile=True, compiler='gfortran', custom_compile_command=None):
        """
        Initialise and run JKTEBOP.

        All input/output files must use the following naming convention:
         - LC only: {target_name}.in/dat/par/out/fit
         - LC + RV: {target_name}.in/par/out/fit,
                    {target_name}-phot.dat,
                    {target_name}-rv1.dat/out,
                    {target_name}-rv2.dat/out

        Parameters
        ----------
        target_name : str
            Name of target, used in filenames (see filename convention)
        recompile : bool, optional
            Whether to recompile JKTEBOP each time
        compiler : str, optional
            Name of compiler to use (e.g. 'gfortran')
        custom_compile_command : str, optional
            Own command to use to compile JKTEBOP if e.g. gfortran is not working.
        """
        self.compiler = compiler
        self.target_name = target_name
        self.custom_compile_command = custom_compile_command

        self.clean_old_results()
        if recompile:
            self.compile_jktebop()
        self.run_jktebop()

    def compile_jktebop(self):
        """ Compile jktebop.f """
        if not self.custom_compile_command:
            # -std=legacy suppresses warnings about arithmetic IF statements
            os.system(f'{self.compiler} -o jktebop jktebop.f -std=legacy')
        else:
            os.system(self.custom_compile_command)

    def clean_old_results(self):
        """ Remove pre-existing output files. JKTEBOP doesn't run if these exist. """
        for ext in [".out", ".fit", ".par", "-rv1.out", "-rv2.out"]:
            try:
                os.system(f'rm -rf {self.target_name}{ext}')
            except FileNotFoundError:
                pass

    def run_jktebop(self):
        """ Run jktebop executable """
        os.system(f'./jktebop {self.target_name}.in')


class TASK3:
    def __init__(self, target_name, rv=False, v_shift=0):
        """ Read the JKTEBOP results from TASK3 output files """
        self.target_name = target_name
        self.v_shift = v_shift
        # Load light curve observations and fit from files
        self.o_phase, self.o_mag, self.o_c = self.data_from_output_file()
        self.c_phase, self.c_mag = self.model_from_fit_file()
        self.ph2 = self.secondary_phase_from_param_file()
        # If applicable, load RV observations and fit (as tuple)
        self.rv_results = self.get_rv_results() if rv else None

    def data_from_output_file(self):
        """ Extract data from the .out file as np.arrays and extend beyond +/-0.5 phase """
        output_file = Table.read(f'{self.target_name}.out', format='ascii')  # Phase-folded light curve
        data_phase = extend_array(wrap_phase(np.array(output_file['PHASE'])), phase=True)
        data_mag = extend_array(np.array(output_file['MAGNITUDE'] - self.v_shift))
        data_residual = extend_array(np.array(output_file['(O-C)']))
        return data_phase, data_mag, data_residual

    def model_from_fit_file(self):
        """Extract fit from the .fit file as np.arrays and extend beyond +/-0.5 phase """
        fit_file = Table.read(f'{self.target_name}.fit', format='ascii')
        fit_phase = wrap_phase(np.array(fit_file['PHASE']))
        fit_mag = np.array(fit_file['MAGNITUDE'] - self.v_shift)

        # Sort, since we'll be plotting as a line instead of points, and extend beyond phase +/-0.5
        fit_phase_sort = extend_array(np.array([fit_phase for fit_phase, _ in sorted(zip(fit_phase, fit_mag))]), True)
        fit_mag_sort = extend_array(np.array([fit_mag for _, fit_mag in sorted(zip(fit_phase, fit_mag))]))
        return fit_phase_sort, fit_mag_sort

    def secondary_phase_from_param_file(self):
        """ Extract phase of secondary eclipse from .par file """
        with open(f'{self.target_name}.par') as f:
            for line in f.readlines():
                if 'Phase of secondary eclipse' in line:
                    ph2 = float(line[43:])
        if ph2 >= 0.5:
            ph2 -= 1
        return ph2

    def get_rv_results(self):
        """ Retrieve model and data from RV output files """
        # Load observed data and errors
        output_file_1 = Table.read(f'{self.target_name + '-rv1'}.out', format='ascii')
        output_file_2 = Table.read(f'{self.target_name + '-rv2'}.out', format='ascii')
        data_phase_1 = extend_array(wrap_phase(np.array(output_file_1['PHASE'])), phase=True)
        data_phase_2 = extend_array(wrap_phase(np.array(output_file_2['PHASE'])), phase=True)
        data_rv1 = extend_array(np.array(output_file_1['RV_STAR_A']))
        data_rv2 = extend_array(np.array(output_file_2['RV_STAR_B']))
        data_residual_1 = extend_array(np.array(output_file_1['(O-C)']))
        data_residual_2 = extend_array(np.array(output_file_2['(O-C)']))

        # Load model RV curves
        fit_file = Table.read(f'{self.target_name}.fit', format='ascii')
        fit_phase = wrap_phase(np.array(fit_file['PHASE']))
        fit_rv1 = np.array(fit_file['RV_A'])
        fit_rv2 = np.array(fit_file['RV_B'])
        # Sort and extend beyond +/-0.5 phase
        fit_phase_sort = extend_array(np.array([fit_phase for fit_phase, _ in sorted(zip(fit_phase, fit_rv1))]), True)
        fit_rv1_sort = extend_array(np.array([fit_mag for _, fit_mag in sorted(zip(fit_phase, fit_rv1))]))
        fit_rv2_sort = extend_array(np.array([fit_mag for _, fit_mag in sorted(zip(fit_phase, fit_rv2))]))

        return (data_phase_1, data_phase_2, data_rv1, data_rv2, data_residual_1, data_residual_2,
                fit_phase_sort, fit_rv1_sort, fit_rv2_sort)

    def plot_lightcurve(self, save=True, y_buffer=0.003):
        """
        Plot data and model of the entire phase-folded light curve

        Parameters
        ----------
        save: bool, optional
            Whether to save the plot to file as <target>-lightcurve.png
        y_buffer: float, optional
            Number to add to the min/max of the observed data to determine the y axis limits
        """
        figure, axes = plt.subplots(2, sharex='col', figsize=(12, 6), gridspec_kw={'height_ratios': [5, 2]})
        plt.subplots_adjust(hspace=0)
        axes[0].scatter(self.o_phase, self.o_mag, color='C1', s=10)
        axes[0].plot(self.c_phase, self.c_mag, color='C0')
        axes[0].set(ylabel='Differential Magnitude', ylim=(max(self.o_mag) + y_buffer, min(self.o_mag) - y_buffer))
        axes[1].scatter(self.o_phase, self.o_c, color='C1', s=10)
        # Centre the plot: primary eclipse first, one full cycle
        x_mid = 0.5 * self.ph2 if self.ph2 > 0 else 0.5 * (self.ph2 + 1)
        axes[1].set(xlabel='Phase', ylabel='Residual', xlim=(x_mid - 0.5, x_mid + 0.5))
        figure.align_labels()
        figure.tight_layout()
        if save:
            plt.savefig(f'{target}-lightcurve.png')
        figure.show()

    def plot_eclipses(self, save=True, ecl_width=0.016, n_bin=1000, y_buffer=0.005, marker_size=20):
        """
        Plot data and model centered on the primary and secondary eclipses

        Parameters
        ----------
        save: bool, optional
            Whether to save the plot to file as <target>-eclipses.png
        ecl_width: float, optional
            Width of the eclipses to restrict the x axes, in phase units
        n_bin: int, optional
            Number of bins to use when binning the residual
        y_buffer: float, optional
            Number to add to the min/max of the observed data to determine the y axis limits
        marker_size: int, optional
            Size of the ax.scatter markers
        """
        # Make binned residuals to show structure
        binned_phase, binned_residual = make_bins(self.o_phase, self.o_c, n_bin=n_bin)

        # Create figure
        figure, axes = plt.subplots(2, 2, figsize=(6, 6), sharex='col', sharey='row',
                                    gridspec_kw={'height_ratios': [5, 2]})
        figure.subplots_adjust(hspace=0.05, wspace=0.05)
        # Primary eclipse
        axes[0, 0].scatter(self.o_phase, self.o_mag, color='C1', s=marker_size)
        axes[0, 0].plot(self.c_phase, self.c_mag, color='C0')
        # Primary residual
        axes[1, 0].scatter(self.o_phase, self.o_c, color='C1', s=marker_size)
        axes[1, 0].scatter(binned_phase, binned_residual, color='C2', s=marker_size)
        # Secondary eclipse
        axes[0, 1].scatter(self.o_phase, self.o_mag, color='C1', s=marker_size)
        axes[0, 1].plot(self.c_phase, self.c_mag, color='C0')
        # Secondary residual
        axes[1, 1].scatter(self.o_phase, self.o_c , color='C1', s=marker_size)
        axes[1, 1].scatter(binned_phase, binned_residual, color='C2', s=marker_size)
        # Axis options - eclipse panels (upper)
        axes[0, 0].set(ylabel='Differential Magnitude', xlim=(-ecl_width, ecl_width),
                       ylim=(max(self.o_mag) + y_buffer, min(self.o_mag) - y_buffer))
        axes[1, 0].set(xlabel='Phase', ylabel='Residual')
        # Axis options - residual panels (lower)
        axes[0, 1].set(xlim=(self.ph2 - ecl_width, self.ph2 + ecl_width))
        axes[1, 1].set(xlabel='Phase', ylim=(min(self.o_c) * 1.02, max(self.o_c) * 1.02))
        axes[1, 0].xaxis.set_minor_locator(MultipleLocator(0.005))
        axes[1, 1].xaxis.set_minor_locator(MultipleLocator(0.005))

        # Align and save
        figure.align_labels()
        figure.tight_layout()
        figure.show()
        if save:
            figure.savefig(f'{target}-eclipses.png')

    def plot_rv_curve(self, save=True, marker_size=60):
        """
        Plot data and model of the entire phase-folded radial velocity curve

        Parameters
        ----------
        save: bool, optional
            Whether to save the plot to file as <target>-rvs.png
        marker_size: int, optional
            Size of the ax.scatter markers
        """
        if not self.rv_results:
            raise RuntimeError('RV results not available! Make sure you set rv=True when initialising TASK3')
        o_rv_phase_1, o_rv_phase_2, o_rv1, o_rv2, o_c_rv_1, o_c_rv_2, c_rv_phase, c_rv1, c_rv2 = self.rv_results

        figure, axes = plt.subplots(2, sharex='col', figsize=(10, 6), gridspec_kw={'height_ratios': [4, 2]})
        axes[0].scatter(o_rv_phase_1, o_rv1, color='C3', s=marker_size)
        axes[0].scatter(o_rv_phase_2, o_rv2, color='C1', s=marker_size)
        axes[0].plot(c_rv_phase, c_rv1, color='C3')
        axes[0].plot(c_rv_phase, c_rv2, color='C1')
        axes[0].set(ylabel='RV [km/s]', ylim=(min(min(c_rv1), min(c_rv2)) - 5, max(max(c_rv1), max(c_rv2)) + 5))
        axes[1].scatter(o_rv_phase_1, o_c_rv_1, color='C3', s=marker_size)
        axes[1].scatter(o_rv_phase_2, o_c_rv_2, color='C1', s=marker_size)
        axes[1].set(xlabel='Phase', ylabel='Residual', xlim=(-0.5, 0.5))
        # Align and save
        figure.align_labels()
        figure.tight_layout()
        figure.show()
        if save:
            figure.savefig(f'{target}-rvs.png')

    def plot_rv_lc(self, save=True, y_buffer=0.003, rv_marker_size=60):
        """
        Plot data and model of the entire phase-folded light curve and radial velocity curve

        Parameters
        ----------
        save: bool, optional
            Whether to save the plot to file as <target>-lc-rvs-combo.png
        y_buffer: float, optional
            Number to add to the min/max of the observed light curve data to determine the y axis limits
        rv_marker_size: int, optional
            Size of the ax.scatter markers for the RV panel
        """
        if not self.rv_results:
            raise RuntimeError('RV results not available! Make sure you set rv=True when initialising TASK3')
        o_rv_phase_1, o_rv_phase_2, o_rv1, o_rv2, o_c_rv_1, o_c_rv_2, c_rv_phase, c_rv1, c_rv2 = self.rv_results

        figure, axes = plt.subplots(2, sharex='col', figsize=(10, 6), gridspec_kw={'height_ratios': [3, 3]})
        # Light curve panel
        axes[0].scatter(self.o_phase, self.o_mag, color='C4', s=2)
        axes[0].plot(self.c_phase, self.c_mag, color='C5')
        axes[0].set(ylabel='Magnitude', ylim=(max(self.o_mag) + y_buffer, min(self.o_mag) - y_buffer))
        # RV panel
        axes[1].scatter(o_rv_phase_1, o_rv1, color='C3', s=rv_marker_size, zorder=0)
        axes[1].scatter(o_rv_phase_2, o_rv2, color='C1', s=rv_marker_size, zorder=0)
        axes[1].plot(c_rv_phase, c_rv1, color='C3')
        axes[1].plot(c_rv_phase, c_rv2, color='C1')
        axes[1].set(ylabel='RV [km/s]', ylim=(min(min(c_rv1), min(c_rv2)) - 5, max(max(c_rv1), max(c_rv2)) + 5),
                    xlabel='Phase', xlim=(-0.6, 0.6))
        # Align and save
        figure.align_labels()
        figure.tight_layout()
        figure.show()
        if save:
            figure.savefig(f'{target}-lc-rvs-combo.png')


class TASK8:
    def __init__(self, target_name):
        self.target_name = target_name
        self.samples = self.load_results()

    def load_results(self):
        """ Retrieve TASK8 Monte Carlo results from .fit file """
        return Table.read(f'{self.target_name}.fit', format='ascii')

    def print_results(self):
        """ Print mean and standard deviation of the MC samples from TASK8 """
        print(f"TASK8 Monte Carlo results for {self.target_name}. Number of MC iterations: {len(self.samples)}")
        print("Parameter       Mean         Stdev")
        print("---------       -------      -------")
        for col in list(self.samples.itercols())[2:]:
            mean = np.mean(np.array(col))
            stdev = np.std(np.array(col))
            print(f"{col.name.ljust(10)} {mean:>12.5f} {stdev:>12.5f}")

    def plot_corner(self, n_params=13, save=False):
        """
        Make corner plot to display MC samples from TASK8

        Parameters
        ----------
        n_params: int
            Number of free parameters used in MC fit. The .fit file presents derived as well as fitted parameters!
        save: bool, optional
            Whether to save the plot to file as <target>-corner.png
        """
        if n_params >= 15:
            print("Caution: Using a large number of parameters may overload your RAM!")
        samples = self.samples[self.samples.colnames[2:]]
        samples = np.array([samples[col] for col in samples.colnames]).T
        figure = corner.corner(samples[0:n_params+1], labels=self.samples.colnames[2:])
        figure.show()
        if save:
            figure.savefig(f'{target}-corner.png')

    def make_task2_with_best_fit(self):
        """ Generate a TASK2 or TASK3 file for best parameters, so the user (or a later function) can rerun and plot """
        pass


if __name__ == "__main__":
    # PyJKTEBOP currently supports TASK3 and TASK8 of JKTEBOP, with/without RVs
    # Future improvements:
    # - support TASK9 (Residual Permutation)
    # - automatically take TASK8 best results, put into TASK2/TASK3 file to make best model and plot
    # - class-based structure with command-line interface

    target = 'hd161701'  # Name of target / JKTEBOP files
    task = 3  # 3, 8, or 9
    rerun = False  # True = run JKTEBOP. False = plot existing results
    rv_fit = True  # True = using RVs. False = standard LC only. Expects extensions -rv1, -rv2, -phot on .dat files

    if rerun:
        j = JKTEBOP(target)

    if task == 3:
        t = TASK3(target, rv=True)
        t.plot_lightcurve()
        t.plot_eclipses()
        t.plot_rv_curve()
        t.plot_rv_lc()

    if task == 8:
        t = TASK8(target)
        t.print_results()
        t.plot_corner(n_params=5)
        t.make_task2_with_best_fit()  # TODO
