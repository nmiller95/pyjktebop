JKTEBOP                                          John Southworth (2014/05/02)
=======                                          ============================

This tarball contains the JKTEBOP (version 34) code, which is used to model
the light and radial velocity curves of detached eclipsing binary star systems.
It also contains all input and output files for three example systems:

1) The transiting planetary system WASP-4. One transit light curve is given.
   This is the simplest dataset to understand

2) The eclipsing binary system WW Aurigae. One V-band light curve is given.
   This is a simple example of modelling an eclipsing binary system.

3) The eclipsing binary system LL Aquarii. For this system there is the
   SuperWASP light curve plus the radial velocity measurements from Griffin
   (2013, The Observatory, vol 133, page 156). This is the most complicated
   of the three examples as it includes light curves, radial velocities, use
   of the "chif" command and also times of minimum light.

For information visit:  http://www.astro.keele.ac.uk/~jkt/codes/jktebop.html
Limited support can be obtained by emailing me at:  astro.js[at].keele.ac.uk

JKTEBOP is monolithic and written in (almost) standard FORTRAN 77.  Input and
output are done using only text files.  JKTEBOP should therefore be very easy
to compile and run.


Compilation
-----------

On a UNIX system and depending on your compiler use one of these commands:
  g77 -o jktebop jktebop.f
  g95 -o jktebop jktebop.f
  gfortran -o jktebop jktebop.f
  ifort -o jktebop jktebop.f

I now use ifort on my work desktop and gfortran on my personal laptop.
JKTEBOP should therefore work well with these two codes.

For Windows systems please refer to your FORTRAN compiler user manual (or
alternatively get rid of Windows and use Linux instead).


Running JKTEBOP
---------------

You will need an input file giving the initial parameter values and other
information. This can be obtained by running JKTEBOP with this command:
  jktebop newfile

Once you have an input file (e.g. wasp4.in) run JKTEBOP with this command:
  jktebop wasp4.in

All output will be to text files - JKTEBOP does not create any graphics.
I recommend using GNUPLOT or IDL to make plots of the light curve fits etc,
but there are many (probably hundreds) of alternative packages out there.


Files included
--------------

readme.txt        This file
jktebop.f         JKTEBOP version 34 source code

wasp4.in          Example input file for the transiting planet WASP-4
wasp4.dat         Light curve datafile to be fitted
wasp4.par         Output parameter file for WASP-4
wasp4.out         Output light curve file (contains residuals)
wasp4.fit         Output fit file (contains the best fit for plotting)

wwaur.in          Example input file for the eclipsing binary WW Aurigae
wwaur-V.dat       Light curve datafile to be fitted
wwaur-V.par       Output parameter file for WW Aurigae
wwaur-V.out       Output light curve file (contains residuals)
wwaur-V.fit       Output fit file (contains the best fit for plotting)

llaqr.in          Example input file for the eclipsing binary LL Aquarii
llaqr-phot.dat    Light curve datafile to be fitted
llaqr-rv1.dat     Radial velocity measurements for the primary star
llaqr-rv2.dat     Radial velocity measurements for the secondary star
llaqr.par         Output parameter file for LL Aquarii
llaqr.out         Output light curve file (contains residuals)
llaqr.fit         Output fit file (contains the best fit for plotting)



