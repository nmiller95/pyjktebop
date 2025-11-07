###############
Getting started
###############

Prepare JKTEBOP input files
===========================
You can generate a blank input file with the following commands:

.. prompt:: bash
    gfortran -o jktebop.f jktebop
    ./jktebop newfile
PyJKTEBOP will assume that your file naming structure is the following:

* ``<target>.in`` - configuration file
* ``<target>-phot.dat`` - photometry data, in HJD and normalised magnitudes (out of eclipse = 0)
* ``<target>.out``, ``<target>.par``, ``<target>.fit`` - output files (N.B. there's no ``.out`` for TASK8)

If you're also fitting radial velocities, these files should be named as follows:

* ``<target>-rv1.dat``, ``<target>-rv2.dat`` - radial velocity data, in HJD and km/s
* ``<target>-rv1.out``, ``<target>-rv2.out`` - output files

See the full `JKTEBOP <https://www.astro.keele.ac.uk/jkt/codes/jktebop.html>`__ documentation and example files
for complete instructions on preparing input files.

Running the command-line version
================================
This option is best for pre-written routines. For more flexibility, consider scripting.
Here are a few example use cases of the command-line version of PyJKTEBOP.
If running the script as written in your IDE, make sure to set the run configuration to suit your needs.

.. note::
    When running JKTEBOP for the first time, make sure to set the flag ``-e`` to ensure it compiles.

Run TASK3 light curve fit; plot and save
----------------------------------------

.. prompt:: bash

    python pyjktebop.py target 3 -rps
Plot and save an existing TASK3 light curve and RV fit
------------------------------------------------------

.. prompt:: bash

    python pyjktebop.py target 3 -vps
Plot and save existing TASK8 light curve results
------------------------------------------------
Here, we also specify the number of parameters to show in the corner plot (default = 8).
These are counted from the leftmost column of the ``<target>.fit`` output file.

.. prompt:: bash

    python pyjktebop.py target 8 -ps --n_corner_params 10
Scripting
=========
For more flexibility, you can import PyJKTEBOP as a package (but make sure to place ``jktebop.f`` in the same directory).
Future versions may allow more flexibility in this.

.. code-block:: python

    from pyjktebop import JKTEBOP, TASK3, TASK8
    # Compile and run JKTEBOP
    JKTEBOP('target', recompile=True, compiler='ifort')

    # Read in TASK3 results
    t = TASK3('target', rv=True)
    t.plot_eclipses(ecl_width=0.02, n_bin=5000)
    t.plot_lightcurve(save=False, y_buffer=0.005)
    t.plot_rv_curve(marker_size=30)
    t.plot_rv_lc(y_buffer=0.005, rv_marker_size=30)

    # Read in and plot a different target's existing TASK8 results
    t8 = TASK8('target2', ld_a='cub', ld_b='sqrt')
    t8.plot_corner(n_params=5, save=True)
    t3 = t8.model_from_best_fit(recompile=False)
    t3.plot_eclipses()
See the full documentation for information on each class/function available.