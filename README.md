# PyJKTEBOP: A simple JKTEBOP python interface

Compiles and runs the eclipsing binary fitting code JKTEBOP via a python script. Returns plots of the light curve (and radial velocity curve) and fit. Currently TASK3 (Levenberg-Marquardt minimisation) is supported, TASK 8 (Monte Carlo) and TASK9 (residual permutation) will follow in the near future.

This interface uses JKTEBOP (version 43); for information visit [John's website](http://www.astro.keele.ac.uk/~jkt/codes/jktebop.html).


Compilation
-----------

JKTEBOP is written in FORTRAN 77, so you will need to check you have a suitable compiler on your system.

You may need to modify the os.system() commands to suit your system. Quoting directly from JKTEBOP's documentation:

> On a UNIX system and depending on your compiler use one of these commands:
> * g77 -o jktebop jktebop.f
> * g95 -o jktebop jktebop.f
> * gfortran -o jktebop jktebop.f
> * ifort -o jktebop jktebop.f
>
> For Windows systems please refer to your FORTRAN compiler user manual (or alternatively get rid of Windows and use Linux instead)."


Running
-------

In `pyjktebop.py` you should adjust the setup parameters `target`, `rerun`, `rv_fit` and `mag_shift` to suit your needs. PyJKTEBOP assumes your JKTEBOP files have the same base name, i.e. setting `target='llaqr'` implies your files have names 'llaqr.in', 'llaqr.dat', etc. 

To modify your fit parameters, change the values directly in the JKTEBOP configuration (.in) file. 

#### Simultaneous LC-RV fitting
PyJKTEBOP currently supports the JKTEBOP TASK3 routine with a light curve and optionally radial velocity data. Make sure your input files are labelled the same way as the example LL Aqr files in the JKTEBOP tarfile (see the [JKTEBOP website](https://www.astro.keele.ac.uk/jkt/codes/jktebop.html) to download the latest version): as before, with extensions -phot.dat, -rv1.dat and -rv2.dat. Read the JKTEBOP documentation for how to set up a fit with RVs.
