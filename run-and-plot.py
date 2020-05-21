import os
from matplotlib import pyplot as plt
from astropy.table import Table
plt.style.use('custom')

os.system('gfortran -o jktebop jktebop.f')  # Compile JKTEBOP
os.system('rm -rf j052.out')  # JKTEBOP doesn't like overwriting existing files
os.system('rm -rf j052.fit')
os.system('rm -rf j052.par')
os.system('./jktebop j052.in')  # Run JKTEBOP

# Extract the phase of secondary eclipse from the parameter file
with open('j052.par') as f:
    for line in f.readlines():
        if 'Phase of secondary eclipse' in line:
            phase2 = float(line[43:])

# Read in the JKTEBOP output files
out_file = Table.read('j052.out', format='ascii')  # Phase-folded light curve
fit_file = Table.read('j052.fit', format='ascii')  # Best fit model

# Plot the phase-folded data and model
fig, ax = plt.subplots(2, sharex=True, figsize=(8, 6),
                       gridspec_kw={'height_ratios': [4, 2]})
fig.subplots_adjust(hspace=0.05)
ax[0].scatter(out_file['PHASE'], out_file['MAGNITUDE']-9.5, color='#b3cde0')
ax[0].plot(fit_file['PHASE'], fit_file['MAGNITUDE']+1, color='#011f4b')
ax[0].set(ylabel='Magnitude', ylim=(1.45, 0.95))
ax[1].scatter(out_file['PHASE'], out_file['(O-C)']-10.5, color='#005b96')
ax[1].set(xlabel='Phase', ylabel='Residual')
fig.show()
