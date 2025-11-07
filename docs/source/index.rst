############################################
PyJKTEBOP: A simple Python JKTEBOP interface
############################################

**PyJKTEBOP** is a simple, flexible Python interface for using the light curve fitting code
`JKTEBOP <https://www.astro.keele.ac.uk/jkt/codes/jktebop.html>`__.

The JKTEBOP fortran code fits models to light curves and radial velocities of **detached eclipsing binary stars**.
The original code does not have plotting capability, which is where PyJKTEBOP comes in!

PyJKTEBOP currently uses JKTEBOP version 43.

Features
========
- Compile and run JKTEBOP with Python
- Run ``TASK3`` (Levenberg-Marquardt minimisation) with optional RVs
- Run ``TASK8`` (Monte Carlo) and make a corner plot from the samples
- Automatically take best parameters from ``TASK8``, generate and plot the best model using ``TASK3``
- Customisable, auto-generated, publication standard plots via the included matplotlib style file
- Scriptable or command-line interface

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage
   pyjktebop