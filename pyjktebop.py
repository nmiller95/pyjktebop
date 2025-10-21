import os
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import corner

plt.style.use('jktebop.mplstyle')


def compile_jktebop():
    """Use bash command to compile jktebop.f"""
    os.system('gfortran -o jktebop jktebop.f')


def clean_old_results(target_name: str):
    """ Remove pre-existing output files. JKTEBOP doesn't run if these exist. """
    for ext in ['out', 'fit', 'par']:
        try:
            os.system(f'rm -rf {target_name}.{ext}')
        except FileNotFoundError:
            pass


def run_jktebop(target_name: str):
    """Run jktebop executable"""
    os.system(f'./jktebop {target_name}.in')


def wrap_phase(phase_array):
    """ Convert phase array from (0,1) to (-0.5,0.5) """
    for i, ph in enumerate(phase_array):
        if ph >= 0.5:
            phase_array[i] = ph - 1
    return phase_array


def data_from_output_file(target_name: str, v_shift=0):
    """ Extract data from the .out file as np.arrays """
    output_file = Table.read(f'{target_name}.out', format='ascii')  # Phase-folded light curve
    data_phase = wrap_phase(np.array(output_file['PHASE']))
    data_mag = np.array(output_file['MAGNITUDE'] - v_shift)
    data_residual = np.array(output_file['(O-C)'])

    # Extend the arrays to wrap properly if ph2=0.5
    data_phase = extend_array(data_phase, True)
    data_mag = extend_array(data_mag)
    data_residual = extend_array(data_residual)
    return data_phase, data_mag, data_residual


def model_from_fit_file(target_name: str):
    """Extract fit from the .fit file as np.arrays """
    fit_file = Table.read(f'{target_name}.fit', format='ascii')
    fit_phase = wrap_phase(np.array(fit_file['PHASE']))
    fit_mag = np.array(fit_file['MAGNITUDE'])

    # Sort, since we'll be plotting as a line instead of points
    fit_phase_sort = [fit_phase for fit_phase, _ in sorted(zip(fit_phase, fit_mag))]
    fit_mag_sort = [fit_mag for _, fit_mag in sorted(zip(fit_phase, fit_mag))]

    # Extend the arrays to wrap properly if ph2=0.5
    fit_phase_sort = extend_array(np.array(fit_phase_sort), True)
    fit_mag_sort = extend_array(np.array(fit_mag_sort))
    return fit_phase_sort, fit_mag_sort


def secondary_phase_from_param_file(target_name: str) -> float:
    """Extract phase of secondary eclipse from .par file """
    with open(f'{target_name}.par') as f:
        for line in f.readlines():
            if 'Phase of secondary eclipse' in line:
                ph2 = float(line[43:])
    if ph2 >= 0.5:
        ph2 -= 1
    return ph2


def plot_lightcurve(data_phase, data_mag, data_residual, fit_phase, fit_mag, ph2, save=True, y_axis_buffer=0.003):
    """
    Plot data and model of the entire phase-folded light curve

    Parameters
    ----------
    data_phase: np.array
        Observed phase
    data_mag: np.array
        Observed magnitudes
    data_residual: np.array
        Observed minus calculated magnitudes
    fit_phase: np.array
        Model phase
    fit_mag: np.array
        Model magnitudes
    ph2: float
        Phase of secondary eclipse
    save: bool, optional
        Whether to save the plot to file as <target>-lightcurve.png
    y_axis_buffer: float, optional
        Number to add to the min/max of the observed data to determine the y axis limits
    """
    figure, axes = plt.subplots(2, sharex='col', figsize=(12, 6), gridspec_kw={'height_ratios': [5, 2]})
    plt.subplots_adjust(hspace=0)
    axes[0].scatter(data_phase, data_mag, color='C1', s=10)
    axes[0].plot(fit_phase, fit_mag, color='C0')
    axes[0].set(ylabel='Differential Magnitude', ylim=(max(data_mag) + y_axis_buffer, min(data_mag) - y_axis_buffer))
    axes[1].scatter(data_phase, data_residual, color='C1', s=10)
    # Centre the plot: primary eclipse first, one full cycle
    x_mid = 0.5 * ph2 if ph2 > 0 else 0.5 * (ph2 + 1)
    axes[1].set(xlabel='Phase', ylabel='Residual', xlim=(x_mid-0.5, x_mid+0.5))
    figure.align_labels()
    figure.tight_layout()
    if save:
        plt.savefig(f'{target}-lightcurve.png')
    figure.show()


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


def plot_eclipses(data_phase, data_mag, data_residual, fit_phase, fit_mag, ph2, save=True,
                  ecl_width=0.016, n_bin=1000, y_axis_buffer=0.005, marker_size=20):
    """
    Plot data and model centered on the primary and secondary eclipses

    Parameters
    ----------
    data_phase: np.array
        Observed phase
    data_mag: np.array
        Observed magnitudes
    data_residual: np.array
        Observed minus calculated magnitudes
    fit_phase: np.array
        Model phase
    fit_mag: np.array
        Model magnitudes
    ph2: float
        Phase of secondary eclipse
    save: bool, optional
        Whether to save the plot to file as <target>-eclipses.png
    ecl_width: float, optional
        Width of the eclipses to restrict the x axes, in phase units
    n_bin: int, optional
        Number of bins to use when binning the residual
    y_axis_buffer: float, optional
        Number to add to the min/max of the observed data to determine the y axis limits
    marker_size: int, optional
        Size of the ax.scatter markers
    """
    # Make binned residuals to show structure
    binned_phase, binned_residual = make_bins(data_phase, data_residual, n_bin=n_bin)

    # Create figure
    figure, axes = plt.subplots(2, 2, figsize=(6, 6), sharex='col', sharey='row',
                                gridspec_kw={'height_ratios': [5, 2]})
    figure.subplots_adjust(hspace=0.05, wspace=0.05)
    # Primary eclipse
    axes[0, 0].scatter(data_phase, data_mag, color='C1', s=marker_size)
    axes[0, 0].plot(fit_phase, fit_mag, color='C0')
    # Primary residual
    axes[1, 0].scatter(data_phase, data_residual, color='C1', s=marker_size)
    axes[1, 0].scatter(binned_phase, binned_residual, color='C2', s=marker_size)
    # Secondary eclipse
    axes[0, 1].scatter(data_phase, data_mag, color='C1', s=marker_size)
    axes[0, 1].plot(fit_phase, fit_mag, color='C0')
    # Secondary residual
    axes[1, 1].scatter(data_phase, data_residual , color='C1', s=marker_size)
    axes[1, 1].scatter(binned_phase, binned_residual, color='C2', s=marker_size)
    # Axis options - eclipse panels (upper)
    axes[0, 0].set(ylabel='Differential Magnitude', xlim=(-ecl_width, ecl_width),
                   ylim=(max(data_mag) + y_axis_buffer, min(data_mag) - y_axis_buffer))
    axes[1, 0].set(xlabel='Phase', ylabel='Residual')
    # Axis options - residual panels (lower)
    axes[0, 1].set(xlim=(ph2 - ecl_width, ph2 + ecl_width))
    axes[1, 1].set(xlabel='Phase', ylim=(min(data_residual) * 1.02, max(data_residual) * 1.02))
    axes[1, 0].xaxis.set_minor_locator(MultipleLocator(0.005))
    axes[1, 1].xaxis.set_minor_locator(MultipleLocator(0.005))

    # Align and save
    figure.align_labels()
    figure.tight_layout()
    if save:
        figure.savefig(f'{target}-eclipses.png')
    figure.show()


def get_rv_results(target_name):
    """Retrieve model and data from RV output files"""
    # Grab data and errors
    output_file_1 = Table.read(f'{target_name+'-rv1'}.out', format='ascii')
    output_file_2 = Table.read(f'{target_name+'-rv2'}.out', format='ascii')
    data_phase_1 = wrap_phase(np.array(output_file_1['PHASE']))
    data_phase_2 = wrap_phase(np.array(output_file_2['PHASE']))
    data_rv1 = np.array(output_file_1['RV_STAR_A'])
    data_rv2 = np.array(output_file_2['RV_STAR_B'])
    data_residual_1 = np.array(output_file_1['(O-C)'])
    data_residual_2 = np.array(output_file_2['(O-C)'])

    # Grab model
    fit_file = Table.read(f'{target_name}.fit', format='ascii')
    fit_phase = wrap_phase(np.array(fit_file['PHASE']))
    fit_rv1 = np.array(fit_file['RV_A'])
    fit_rv2 = np.array(fit_file['RV_B'])
    fit_phase_sort = [fit_phase for fit_phase, _ in sorted(zip(fit_phase, fit_rv1))]
    fit_rv1_sort = [fit_mag for _, fit_mag in sorted(zip(fit_phase, fit_rv1))]
    fit_rv2_sort = [fit_mag for _, fit_mag in sorted(zip(fit_phase, fit_rv2))]

    return (data_phase_1, data_phase_2, data_rv1, data_rv2, data_residual_1, data_residual_2,
            fit_phase_sort, fit_rv1_sort, fit_rv2_sort)


def plot_rv_results(*params, save=True):
    """Plots RV curve and combined RV-LC plot"""
    (o_rv_phase_1, o_rv_phase_2, o_rv1, o_rv2, o_c_rv_1, o_c_rv_2, c_rv_phase, c_rv1, c_rv2, data_phase,
     data_mag, data_residual, fit_phase, fit_mag) = params

    # PLOT 1 : RVs only
    figure, axes = plt.subplots(2, sharex='col', figsize=(10, 6), gridspec_kw={'height_ratios': [4, 2]})
    plt.subplots_adjust(hspace=0)
    axes[0].scatter(o_rv_phase_1, o_rv1, color='#8B0000', s=60)
    axes[0].scatter(o_rv_phase_2, o_rv2, color='C1', s=60)
    axes[0].plot(c_rv_phase, c_rv1, color='#8B0000')
    axes[0].plot(c_rv_phase, c_rv2, color='C1')
    axes[0].set(ylabel='RV [km/s]', ylim=(min(o_rv1) - 5, max(o_rv2) + 5))
    axes[1].scatter(o_rv_phase_1, o_c_rv_1, color='#8B0000', s=60)
    axes[1].scatter(o_rv_phase_2, o_c_rv_2, color='C1', s=60)
    axes[1].set(xlabel='Phase', ylabel='Residual', xlim=(-0.5, 0.5))
    figure.align_labels()
    figure.tight_layout()
    figure.show()
    if save:
        figure.savefig(f'{target}-rvs.png')

    # PLOT 2 : EVERYTHING
    figure, axes = plt.subplots(2, sharex='col', figsize=(10, 6), gridspec_kw={'height_ratios': [3, 3]})
    plt.subplots_adjust(hspace=0)
    # Light curve panel
    axes[0].scatter(data_phase, data_mag, color='#b3b3b3', s=2)
    axes[0].plot(fit_phase, fit_mag, color='#0c0f19')
    axes[0].set(ylabel='Magnitude', ylim=(max(data_mag) + 0.05, min(data_mag) - 0.05),)
    # RV panel
    axes[1].scatter(o_rv_phase_1, o_rv1, color='#8B0000', s=60, zorder=0)
    axes[1].scatter(o_rv_phase_2, o_rv2, color='C1', s=60, zorder=0)
    axes[1].plot(c_rv_phase, c_rv1, color='#8B0000')
    axes[1].plot(c_rv_phase, c_rv2, color='C1')
    axes[1].set(ylabel='RV [km/s]', ylim=(min(o_rv1) - 5, max(o_rv2) + 5), xlabel='Phase', xlim=(-0.5, 0.5))
    figure.align_labels()
    figure.tight_layout()
    figure.show()
    if save:
        figure.savefig(f'{target}-lc-rvs-combo.png')


def get_plot_mc_results(target_name, make_corner=True):
    """Retrieve TASK8 Monte Carlo results from .fit file, print and plot results"""
    results = Table.read(f'{target_name}.fit', format='ascii')
    print(f"TASK8 Monte Carlo results for {target_name}. Number of MC iterations: {len(results)}")
    print("Parameter       Mean         Stdev")
    print("---------       -------      -------")
    for col in list(results.itercols())[2:]:
        mean = np.mean(np.array(col))
        stdev = np.std(np.array(col))
        print(f"{col.name.ljust(10)} {mean:>12.5f} {stdev:>12.5f}")

    # Prepare table object for corner plot
    if make_corner:
        results_data = results[results.colnames[2:]]
        results_data = np.array([results_data[col] for col in results_data.colnames]).T
        figure = corner.corner(results_data, labels=results.colnames[2:])
        figure.show()


if __name__ == "__main__":
    # PyJKTEBOP currently supports TASK3 and TASK8 of JKTEBOP, with/without RVs
    # Future improvements:
    # - support TASK9 (Residual Permutation)
    # - automatically take TASK8 best results, put into TASK2/TASK3 file to make best model and plot
    # - class-based structure with command-line interface
    # - rebase master -> main
    # - update RV plotting to publication standard

    target = ''  # Name of target / JKTEBOP files
    # task2_file = 'llaqr_best'
    task = 3  # 3, 8, or 9
    rerun = True  # True = run JKTEBOP. False = plot existing results
    rv_fit = False  # True = using RVs. False = standard LC only. Expects extensions -rv1, -rv2, -phot on .dat files
    make_corner_plot = False  # For TASK8. Will probably break your RAM if you have too many free parameters
    mag_shift = 0  # -5.7265  # Constant to shift magnitude by in plots

    # Procedure to set up and run JKTEBOP
    if rerun:
        # Clean out old result files that prevent JKTEBOP from running
        clean_old_results(target)
        if rv_fit:
            clean_old_results(target + '-phot')
            clean_old_results(target + '-rv1')
            clean_old_results(target + '-rv2')

        # Compile and run JKTEBOP
        compile_jktebop()
        run_jktebop(target)

    if task == 3:
        # Extract light curve data and fit information from output files
        o_phase, o_mag, o_c = data_from_output_file(target, v_shift=mag_shift)
        c_phase, c_mag = model_from_fit_file(target)
        c_mag = np.array(c_mag) - mag_shift
        phase_2 = secondary_phase_from_param_file(target)

        # Make light curve plots
        plot_lightcurve(o_phase, o_mag, o_c, c_phase, c_mag, phase_2)
        plot_eclipses(o_phase, o_mag, o_c, c_phase, c_mag, phase_2)

        # Extract and plot RVs
        if rv_fit:
            rv_params = get_rv_results(target)
            p = *rv_params, o_phase, o_mag, o_c, c_phase, c_mag
            plot_rv_results(*p)

    elif task == 8:
        # Prints mean and stdev of MC results. Can make a corner plot, but may break RAM if too many params
        get_plot_mc_results(target, make_corner=make_corner_plot)

    # elif task == 2:
    #     # Extract light curve data and fit information from output files
    #     o_phase, o_mag, o_c = data_from_output_file(target, v_shift=mag_shift)
    #     c_phase, c_mag = model_from_fit_file(task2_file)
    #     c_mag = np.array(c_mag) - mag_shift
    #     phase_2 = 0.3175753908  # secondary_phase_from_param_file(target)
    #
    #     # Make light curve plots
    #     plot_lightcurve(o_phase, o_mag, o_c, c_phase, c_mag)
    #     plot_eclipses(o_phase, o_mag, o_c, c_phase, c_mag, phase_2)