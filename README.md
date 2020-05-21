# JKTEBOP python interface

Compiles and runs the eclipsing binary fitting code JKTEBOP via a python script. Returns plots of the light curve and fit.

This interface uses JKTEBOP (version 34); for information visit [John's website](http://www.astro.keele.ac.uk/~jkt/codes/jktebop.html).


Compilation
-----------

JKTEBOP is written in FORTRAN 77, so you will need to check you have a suitable compiler on your system.

You may need to modify the os.system() commands to suit your system. From JKTEBOP's documentation:

> On a UNIX system and depending on your compiler use one of these commands:
> * g77 -o jktebop jktebop.f
> * g95 -o jktebop jktebop.f
> * gfortran -o jktebop jktebop.f
> * ifort -o jktebop jktebop.f
>
> For Windows systems please refer to your FORTRAN compiler user manual (or alternatively get rid of Windows and use Linux instead)."


Running
-------

To modify your fit parameters, change the values directly in the .in file. Then simply run the code in python and a plot should appear with the results.
