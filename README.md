# Simple python interface for light curve fitting with JKTEBOP

Compiles and runs the eclipsing binary fitting code JKTEBOP via a python script. Returns nicely formatted plots of the light curve and fit.

<img src="https://nikkehmiller.files.wordpress.com/2021/02/example-eclipses.png" alt="Example output" width="400"/>

This interface uses JKTEBOP (version 40); for information visit [John Taylor's website](http://www.astro.keele.ac.uk/~jkt/codes/jktebop.html).

Installation
------------

Clone the repository to your local machine with
```
git clone https://github.com/nmiller95/pyjktebop.git
```

Setup
-----

### Compilation

JKTEBOP is written in FORTRAN 77, so you will need to check you have a suitable compiler installed on your system.

You may need to modify the `os.system()` commands in `run-and-plot.py` to suit your system (currently it is set to `gfortran`). From JKTEBOP's documentation:

> On a UNIX system and depending on your compiler use one of these commands:
> * `g77 -o jktebop jktebop.f`
> * `g95 -o jktebop jktebop.f`
> * `gfortran -o jktebop jktebop.f`
> * `ifort -o jktebop jktebop.f`
>
> For Windows systems please refer to your FORTRAN compiler user manual (or alternatively get rid of Windows and use Linux instead)."

### Light curve data format

Your data should be in a `.dat` file with no header in the format: time (days), magnitude, error (optional).

### Input parameter file

You will need to generate an `.in` file that contains estimated fit parameters (you can adjust these later) in terminal: 
```
gfortran -o jktebop jktebop.f  # your compile command
./jktebop newfile
```
The python interface is currently optimised for TASK3. I would recommend always leaving 'light scale factor' as a free parameter.

Running
-------

Run the python script in terminal with
```
python3 run-and-plot.py
```

Future fixes & improvements
---------------------------
 
* Support for other JKTEBOP tasks
