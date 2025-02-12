import os
from matplotlib import pyplot as plt
import numpy as np
from astropy.table import Table

plt.style.use('jktebop.mplstyle')


def compile_jktebop():
    """Use bash command to compile jktebop.f"""
    os.system('gfortran -o jktebop jktebop.f')


def clean_old_results(target_name):
    """Use bash commands to clean output files so jktebop can run
    :type target_name: str"""
    # JKTEBOP doesn't run if output files already exist
    os.system(f'rm -rf {target_name}.out')
    os.system(f'rm -rf {target_name}.fit')
    os.system(f'rm -rf {target_name}.par')


def run_jktebop(target_name):
    """Use bash commands to run jktebop executable
    :type target_name: str
    """
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
    plt.subplots_adjust(hspace=0)
    axes[0].scatter(data_phase, data_mag, color='#92a6e5')
    axes[0].plot(fit_phase, fit_mag, color='#0c0f19')
    axes[0].set(ylabel='Magnitude',
                ylim=(max(data_mag) + 0.05, min(data_mag) - 0.05))
    axes[1].scatter(data_phase, data_residual, color='#92a6e5')
    axes[1].set(xlabel='Phase', ylabel='Residual', xlim=(-0.5, 0.5))
    figure.align_labels()
    figure.tight_layout()
    figure.show()


def plot_eclipses(data_phase, data_mag, data_residual, fit_phase, fit_mag, ph2, save=True):
    """Plots data and model centered on the primary and secondary eclipses"""
    figure, axes = plt.subplots(2, 2, figsize=(8, 6), sharex='col', sharey='row',
                                gridspec_kw={'height_ratios': [4, 2]})
    figure.subplots_adjust(hspace=0.05, wspace=0.05)
    # Primary eclipse
    axes[0, 0].scatter(data_phase, data_mag, color='#92a6e5')
    axes[0, 0].plot(fit_phase, fit_mag, color='#0c0f19')
    axes[0, 0].set(ylabel='Magnitude', xlim=(-0.012, 0.012),
                   ylim=(max(data_mag) + 0.05, min(data_mag) - 0.05))
    # Primary residual
    axes[1, 0].scatter(data_phase, data_residual, color='#92a6e5')
    axes[1, 0].set(xlabel='Phase', ylabel='Residual')
    # Secondary eclipse
    axes[0, 1].scatter(data_phase, data_mag, color='#92a6e5')
    axes[0, 1].plot(fit_phase, fit_mag, color='#0c0f19')
    axes[0, 1].set(xlim=(ph2 - 0.012, ph2 + 0.012))
    # Secondary residual
    axes[1, 1].scatter(data_phase, data_residual, color='#92a6e5')
    axes[1, 1].set(xlabel='Phase')
    figure.align_labels()
    figure.tight_layout()
    figure.show()
    if save:
        figure.savefig(f'{target}-eclipses.png')


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
    figure, axes = plt.subplots(2, sharex='col', figsize=(8, 6), gridspec_kw={'height_ratios': [4, 2]})
    plt.subplots_adjust(hspace=0)
    axes[0].scatter(o_rv_phase_1, o_rv1, color='#8B0000', s=60)
    axes[0].scatter(o_rv_phase_2, o_rv2, color='C1', s=60)
    axes[0].plot(c_rv_phase, c_rv1, color='#8B0000')
    axes[0].plot(c_rv_phase, c_rv2, color='C1')
    axes[0].set(ylabel='RV [km/s]', ylim=(max(o_rv2) + 5, min(o_rv1) - 5))
    axes[1].scatter(o_rv_phase_1, o_c_rv_1, color='#8B0000', s=60)
    axes[1].scatter(o_rv_phase_2, o_c_rv_2, color='C1', s=60)
    axes[1].set(xlabel='Phase', ylabel='Residual', xlim=(-0.5, 0.5))
    figure.align_labels()
    figure.tight_layout()
    figure.show()
    if save:
        figure.savefig(f'{target}-rvs.png')

    # PLOT 2 : EVERYTHING
    figure, axes = plt.subplots(2, sharex='col', figsize=(8, 6), gridspec_kw={'height_ratios': [3, 3]})
    plt.subplots_adjust(hspace=0)
    # Light curve panel
    axes[0].scatter(data_phase, data_mag, color='#b3b3b3', s=2)
    axes[0].plot(fit_phase, fit_mag, color='#0c0f19')
    axes[0].set(ylabel='Magnitude', ylim=(max(data_mag) + 0.05, min(data_mag) - 0.05),)
    # RV panel
    axes[1].scatter(o_rv_phase_1, o_rv1, color='#8B0000', s=60)
    axes[1].scatter(o_rv_phase_2, o_rv2, color='C1', s=60)
    axes[1].plot(c_rv_phase, c_rv1, color='#8B0000')
    axes[1].plot(c_rv_phase, c_rv2, color='C1')
    axes[1].set(ylabel='RV [km/s]', ylim=(max(o_rv2) + 5, min(o_rv1) - 5), xlabel='Phase', xlim=(-0.5, 0.5))
    figure.align_labels()
    figure.tight_layout()
    figure.show()
    if save:
        figure.savefig(f'{target}-lc-rvs-combo.png')


if __name__ == "__main__":

    target = 'llaqr'
    rv_fit = True

    # Clean out old result files that prevent JKTEBOP from running
    clean_old_results(target)
    if rv_fit:
        clean_old_results(target + '-phot')
        clean_old_results(target + '-rv1')
        clean_old_results(target + '-rv2')

    # Compile and run JKTEBOP
    compile_jktebop()
    run_jktebop(target)

    # Extract information from output files
    o_phase, o_mag, o_c = data_from_output_file(target)
    c_phase, c_mag = model_from_fit_file(target)
    phase_2 = secondary_phase_from_param_file(target)

    if rv_fit:
        rv_params = get_rv_results(target)
        p = *rv_params, o_phase, o_mag, o_c, c_phase, c_mag
        plot_rv_results(*p)

    # Make plots
    plot_lightcurve(o_phase, o_mag, o_c, c_phase, c_mag)
    plot_eclipses(o_phase, o_mag, o_c, c_phase, c_mag, phase_2)


