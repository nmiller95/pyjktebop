import os
from matplotlib import pyplot as plt
from astropy.table import Table

os.system('gfortran -o jktebop jktebop.f')  # Compile JKTEBOP
os.system('rm -rf j052.out')  # JKTEBOP doesn't like overwriting existing files
os.system('rm -rf j052.fit')
os.system('rm -rf j052.par')
os.system('./jktebop j052.in')  # Run JKTEBOP

# Read in the JKTEBOP output files
out_file = Table.read('j052.out', format='ascii')
fit_file = Table.read('j052.fit', format='ascii')

# Plot the phase-folded data and model
fig, ax = plt.subplots(figsize=(10, 5))
ax.scatter(out_file['PHASE'], out_file['MAGNITUDE'], color='orange')
ax.plot(fit_file['PHASE'], fit_file['MAGNITUDE']+10.5, color='black')
ax.set(xlabel='Phase', ylabel='Magnitude', ylim=(10.95, 10.45))
fig.show()
