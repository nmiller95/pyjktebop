import os
from matplotlib import pyplot as plt
import numpy as np
from astropy.table import Table

plt.style.use('jktebop.mplstyle')


def compile_jktebop():
    """Use bash command to compile jktebop.f"""
    os.system('gfortran -o jktebop jktebop.f')


def run_jktebop(target_name):
    """Use bash commands to run jktebop executable
    :type target_name: str
    """
    # JKTEBOP doesn't run if output files already exist
    os.system(f'rm -rf {target_name}.out')
    os.system(f'rm -rf {target_name}.fit')
    os.system(f'rm -rf {target_name}.par')
    # Run jktebop executable
    os.system(f'./jktebop {target_name}.in')


def wrap_phase(phase_array):
    """Function to convert phase array from (0,1) to (-0.5,0.5)"""
    for i, ph in enumerate(phase_array):
        if ph >= 0.5:
            phase_array[i] = ph - 1
    return phase_array


def data_from_output_file(target_name, v_shift=0):
    """Extract correctly formatted arrays from the .out file
    :type target_name: str
    :type v_shift: float
    :returns: np.arrays of phase, magnitude, residual
    """
    output_file = Table.read(f'{target_name}.out', format='ascii')  # Phase-folded light curve
    data_phase = wrap_phase(np.array(output_file['PHASE']))
    data_mag = np.array(output_file['MAGNITUDE'] - v_shift)
    data_residual = np.array(output_file['(O-C)'] - v_shift)
    return data_phase, data_mag, data_residual


def model_from_fit_file(target_name):
    """Extract correctly formatted arrays from the .fit file
        :type target_name: str
        :returns: np.arrays of phase, magnitude
        """
    fit_file = Table.read(f'{target_name}.fit', format='ascii')
    fit_phase = wrap_phase(np.array(fit_file['PHASE']))
    fit_mag = np.array(fit_file['MAGNITUDE'])
    fit_phase_sort = [fit_phase for fit_phase, _ in sorted(zip(fit_phase, fit_mag))]
    fit_mag_sort = [fit_mag for _, fit_mag in sorted(zip(fit_phase, fit_mag))]
    return fit_phase_sort, fit_mag_sort


def secondary_phase_from_param_file(target_name):
    """Extract correctly formatted phase of secondary eclipse
    :type target_name: str
    :returns: phase of secondary eclipse as float
    """
    with open(f'{target_name}.par') as f:
        for line in f.readlines():
            if 'Phase of secondary eclipse' in line:
                ph2 = float(line[43:])
    if ph2 >= 0.5:
        ph2 -= 1
    return ph2


def plot_lightcurve(data_phase, data_mag, data_residual, fit_phase, fit_mag):
    """Plots data and model of the entire phase-folded light curve"""
    figure, axes = plt.subplots(2, sharex='col', figsize=(8, 6),
                                gridspec_kw={'height_ratios': [4, 2]})
    axes[0].scatter(data_phase, data_mag, color='#a1a1a1')
    axes[0].plot(fit_phase, fit_mag)
    axes[0].set(ylabel='Magnitude',
                ylim=(max(data_mag) + 0.05, min(data_mag) - 0.05))
    axes[1].scatter(data_phase, data_residual, color='C1')
    axes[1].set(xlabel='Phase', ylabel='Residual')
    figure.align_labels()
    figure.show()


def plot_eclipses(data_phase, data_mag, data_residual, fit_phase, fit_mag, ph2, save=True):
    """Plots data and model centered on the primary and secondary eclipses"""
    figure, axes = plt.subplots(2, 2, figsize=(8, 6), sharex='col', sharey='row',
                                gridspec_kw={'height_ratios': [4, 2]})
    figure.subplots_adjust(hspace=0.05)
    # Primary eclipse
    axes[0, 0].scatter(data_phase, data_mag, color='C1')
    axes[0, 0].plot(fit_phase, fit_mag, color='#b1b1b1')
    axes[0, 0].set(ylabel='Magnitude', xlim=(-0.015, 0.015),
                   ylim=(max(data_mag) + 0.05, min(data_mag) - 0.05))
    # Primary residual
    axes[1, 0].scatter(data_phase, data_residual, color='C1')
    axes[1, 0].set(xlabel='Phase', ylabel='Residual')
    # Secondary eclipse
    axes[0, 1].scatter(data_phase, data_mag, color='C1')
    axes[0, 1].plot(fit_phase, fit_mag, color='#b1b1b1')
    axes[0, 1].set(xlim=(ph2 - 0.015, ph2 + 0.015))
    # Secondary residual
    axes[1, 1].scatter(data_phase, data_residual, color='C1')
    axes[1, 1].set(xlabel='Phase')
    figure.align_labels()
    figure.show()
    if save:
        figure.savefig(f'{target}-eclipses.png')


if __name__ == "__main__":

    target = 'TESS-j052'

    # Compile and run JKTEBOP
    compile_jktebop()
    run_jktebop(target)

    # Extract information from output files
    o_phase, o_mag, o_c = data_from_output_file(target)
    c_phase, c_mag = model_from_fit_file(target)
    phase_2 = secondary_phase_from_param_file(target)

    # Make plots
    plot_lightcurve(o_phase, o_mag, o_c, c_phase, c_mag)
    plot_eclipses(o_phase, o_mag, o_c, c_phase, c_mag, phase_2)
