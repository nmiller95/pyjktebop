!=======================================================================
!     PROGRAM JKTEBOP       John Southworth (taylorsouthworth~gmail.com)
!                           Astrophysics Group, Keele University, UK
!-----------------------------------------------------------------------
! PARAMETERS OF THE MODEL:
! V(1) = central surface brightness ratio (starB/starA)
! V(2) = sum of the fractional radii (radii divided by semimajor axis)
! V(3) = ratio of stellar radii (starB/starA)
! V(4) = linear limb darkening coefficient for star A
! V(5) = linear limb darkening coefficient for star B
! V(6) = orbital inclination (degrees)
! V(7) = e cos(omega) OR ecentricity
! V(8) = e sin(omega) OR omega(degrees)
! V(9) = gravity darkening coefficient of star A
! V(10) = gravity darkening coefficient of star B
! V(11) = reflected light for star A
! V(12) = reflected light for star A
! V(13) = mass ratio (starB/starA) for the light curve calculation
! V(14) = tidal lead/lag angle (degrees)
! V(15) = third light in units where (LA + LB + LC = 1)
! V(16) = phase correction factor (i.e. phase of primary eclipse)
! V(17) = light scaling factor (magnitudes)
! V(18) = integration ring size (degrees)
! V(19) = orbital period (days)
! V(20) = ephemeris timebase (days)
! V(21) = limb darkening coefficient 2 for star A
! V(22) = limb darkening coefficient 3 for star A
! V(23) = limb darkening coefficient 4 for star A
! V(24) = limb darkening coefficient 2 for star B
! V(25) = limb darkening coefficient 3 for star B
! V(26) = limb darkening coefficient 4 for star B
! V(27) = velocity amplitude of star A (km/s)
! V(28) = velocity amplitude of star B (km/s)
! V(29) = systemic velocity of star A (km/s)
! V(30) = systemic velocity of star B (km/s)
! V(31-90) = twenty lots of sine curve parameters  [T0,period,amplitude]
! V(91-540) = 50 polynomial pars [pivot,Tstart,Tend,const,x,x2,x3,x4,x5]
! VEXTRA(1) = star A fractional radius
! VEXTRA(2) = star B fractional radius
! VEXTRA(3) = light ratio starB/starA
! VEXTRA(4) = eccentricity
! VEXTRA(5) = periastron longitude
! VEXTRA(6) = impact par (pri ecl)
! VEXTRA(7) = impact par (sec ecl)
! VEXTRA(8) = reduced chi-squared
! VEXTRA(9) = orbital semimajor axis
! VEXTRA(10) = mass ratio from RVs
! VEXTRA(11) = mass of star A (Msun)
! VEXTRA(12) = mass of star B (Msun)
! VEXTRA(13) = radius of star A (Rsun)
! VEXTRA(14) = radius of star B (Rsun)
! VEXTRA(15) = logg of star A (c.g.s.)
! VEXTRA(16) = logg of star B (c.g.s.)
! VEXTRA(17) = density of star A (solar units)
! VEXTRA(18) = density of star B (solar units)
!-----------------------------------------------------------------------
! If Claret's four-parameter limb darkening law is being used then:
! I(mu)/I(0) = 1 - a*[1-mu^0.5] - b*[1-mu] - c*[1-mu^1.5] - d*[1-mu^2]
! The coefficients a,b,c,d for star A will be in V(4),V(21),V(22),V(23),
! and these coefficients for star B will be in V(5),V(24),V(25),V(26).
!-----------------------------------------------------------------------
! The identity of star A versus that of star B is not important for the
! code as it will run happily either way. But primary eclipse is defined
! to be at phase zero (if the phase correction parameter is zero).
! By convention the primary eclipse is deeper than the secondary eclipse
! and star A is the star which is eclipsed during primary eclipse (i.e.
! is at inferior conjunction). This normally means that star A is the
! HOTTER star (i.e. has the higher surface brightness), but exceptions
! are possible in a small fraction of cases with eccentric orbits.
!-----------------------------------------------------------------------
! JKTEBOP VERSIONS:
! Version 1:  Simplex minimisation algorithm and new input / output used
! Version 2:  Monte Carlo simulation and parametr perturbation algorithm
! Version 3:  Adjustment to Monte Carlo LD coeffs and input/output files
! Version 4:  Solves for sum of radii and convergence criterion modified
! Version 5:  Added TASK0 to find LD and GD coeffs.  Minor modifications
! Version 6:  Reflection and  scale factor  can all be fixed or adjusted
! Version 7:  Can use (e,w) or (ecosw,esinw).    SFACT modified to be in
!             magnitudes; observ'nal errors found; spherical star option
! Version 8:  Bootstrapping error analysis added and the output modified
! Version 9:  Command-line arguments allowed, Monte Carlo without param-
!             eter kicking option and fitting for period and Tzero added
! Version 10: Can now use 99999 datapoints. Whole code now in magnitudes
! Version 11: Bug fixes, tasks renumbered,  added sigma clipping and the
!             global fit procedures, but not thoroughly tested these yet
! Version 12: Nonlinear limb darkening law, fitting for times of minimum
!             light, and FMAX corrections included (from Alvaro Gimenez)
! Version 13: Removed  BILINEAR  and modified  TASK1  to just call JKTLD
!             Modified input file and arrays  for the diff types of data
! Version 14: Fixed the requirement for inputtd INTRING to be an integer
!             Fixed formal errors (when observational ones not supplied)
!             Sorted out numerical derivatives problem when  i ~ 90  deg
! Version 15: Added TASK 9 and cubic LD law, modified simulation output,
!             made MRQMIN a factor of 3 faster, and working on red noise
! Version 16: Added ability to include sine perturbations on parameters.
! Version 17: Added ability to specify a third light and its uncertainty
! Version 18: Made possible to fit for r1 and r2 instead of r1+r2 and k.
! Version 19: Added VARY=3 flag + finer fit phasing if r1 or r2 is small
! Version 20: Polynomial, optimisation check, TASK4 debug, DV centralise
! Version 21: Added input of ecosw and esinw observational  constraints.
! Version 22: Added input of  e and omega  as observational constraints.
! Version 23: (did not) fix a minor bug with calculation of third light.
! Version 24: Moved to gfortran compiler.  Adjusted date_and_time usage.
!             Fixed bug with third light. Added ECQUADPHASES subroutine.
! Version 25: Added numerical integration over long exposure times.
! Version 26: Fixed bug with TASK 5 and modified TASK 5 output slightly.
! Version 27: Fixed TASK 8 bug, converted all to real*8, 999999 datapnt.
! Version 28: Corrected ECQUADPHASES. All write(*) or print* => write(6)
! Version 29: Modified to have a specific MAXDATA and split DATA arrays.
! Version 30: Changed max iterations to 100; check for this in TASK 789.
! Version 31: Added RVs, four-parameter LD, calculation of impact param-
!             eter and physical properties of stars.  Modified input and
!             output. Added "chif" function to adjust the errorbar sizes
! Version 32: Added OCCULTQUAD routine as alternative to the EBOP model.
! Version 33: Extended to 9 polynomials, applied to subsets of the data.
! Version 34: Can now put polynomials and sines on e, w, ecosw and esinw
! Version 35: Can now put polynomials and sines on T0 to allow for TTVs.
! Version 36: Can now put polynomials and sines on orbital parameters.
! Version 37: Fixed inability to have negative velocity amplitudes, and
!             modified to the IAU 2015 resolution Rsun and GMsun values.
! Version 38: Added eclipse durations, reviewed gravity darkening and
!             added printout of the most discrepant light curve points.
! Version 39: Allowed for input and output files of 50 characters length
! Version 40: Fixed problem with ellipsoidal effect and nonlinear LD.
! Version 41: Extended to having 20 sines & 50 polynomials. Minor fixes.
! Version 42: Implemented the power-2 LD law (using c,alpha *or* h1,h2),
!             and changed task 2 to allow input reflection coefficients.
! Version 43: Implemented combination LDs for quad, log, sqrt and cubic,
!             checks for data outside limits or NaN, added nofit option.
! Last modified: 2022/12/09
!-----------------------------------------------------------------------
! POSSIBLE MODIFICATIONS IN FUTURE:
! 1) Extend to WD2003 and WINK
! 2) Port to F90 or F95 to have long lines, modules, and improved output
! 3) Incorporate change of omega (apsidal motion)
! 4) Include a light-time effect
! 5) Allow for multiple light and radial velocity curves
!-----------------------------------------------------------------------
! MISCELLANEOUS NOTES:
! 1) Phase shift has been redefined compared to original EBOP so that it
!    corresponds directly to the phase of primary minimum.
! 2) MRQMIN adjusts coeffs only if VARY (called 'ia' in MRQMIN) is 1
! 3) If VARY=2 then the parameter is fixed during the initial fit but is
!    perturbed by a set amount (flat distribution) for later analyses.
! 4) If VARY(11) and / or  VARY(12) are "-1" then V(11) and/or V(12) are
!    calculated from the system geometry;  if set to 0 they are fixed at
!    the input value and if 1 are freely adjusted to best fit.
! 5) If the mass ratio is <= 0 then both stars are assumed to be spheres
!    The mass ratio is used only in the calculation of the light curve -
!    it does not have any effect on the radial velocities.
! 6) If ecosw > 5.0 then (ecosw,esinw) will be taken to be (10+e,omega)
!    and fitting will occur using e and omega as parameters. e and omega
!    can be strongly correlated, but this option is useful if e is known
!    but omega isn't; this can happen for EBs exhibiting apsidal motion.
! 7) Observational errors are looked for in the  input light curve file.
!    If they are not found then  equal weight  is given  to each  point.
! 8) Nonlinear LD is now supported for the two-coefficient  logarithmic,
!    quadratic, square-root and power-2 laws, and for the "Claret" four-
!    parameter law. The type of law must be specified on input.  Star B
!    can also be forced to have the same coefficients as star A. If the
!    Claret four-parameter law is used, it must be used for both stars.
! 9) Fitting for times of minimum light is directly possible.  The cycle
!    numbers and times are inputted  on lines  immediately below all the
!    parameter lines in the input file.
! 10) EBOP results are symmetric about 90 degrees for inclination (which
!     means i=89 gives same answer as i=91),  which causes problems with
!     numerical derivativs when i>89.9. In this case some extra is added
!     to the numerical derivative to keep the solution slightly below 90
! 11) If input (rA+rB) is negative, (rA+rB,k) is interpreted as (-rA,rB)
! 12) If  VARY=3  then the parameter is optimised during all fits but is
!     not perturbed by a set amount in later analyses (eg. Monte Carlo).
! 13) If the maximum number of datapoints is too small for you,  you can
!     now change it to your choice of value.  Look for these quantities:
!     MAXDATA and DATAFORM near the top of the main program.
!-----------------------------------------------------------------------
! TASK NUMBERS AND PURPOSES:
! (1) This outputs LD coefficients for a given Teff, logg, [M/H], Vmicro
! (2) This outputs a model light curve for fixed input parameters.
! (3) This fits a model to an observed light curve and outputs results.
! (4) This fits a model, rejects discrepant observations, and refits.
! (5) This does a pseudo-global minimisation by perturbing input params.
! (6) This investigates how different parameters vary around best fit.
! (7) This conducts bootstrapping simulations to find robust errors.
! (8) This conducts Monte Carlo simulations to find robust errors.
! (9) This conducts residual permutations to deal with correlated noise.
!-----------------------------------------------------------------------
! LANGUAGE:  JKTEBOP is written in FORTRAN 77, using several extensions
!   to the ANSI standard:   ==   <=   <   >   >=   /=   !   enddo  endif
!   Until version 24 it was only ever compiled in g77.
! gfortran compiler:  this compiles successfully, executes correctly and
!   is the compiler I use for JKTEBOP. I currently use gfortran from gcc
!   gcc version 5.4.0 running on Linux Mint 18.3
! g77 compiler: this is obsolete and should not be used. Use gfortran.
! g95 compiler: this is available and should work but I no longer use it
! Intel-Fortran compiler: up to JKTEBOP version 35 this was periodically
!   checked and found to work well.  I no longer use this compiler as my
!   license has expired. JKTEBOP version 26 and earlier must be compiled
!   with the -r8 command-line flag.    For JKTEBOP version 29 onwards, I
!   found that implicit-length arrays caused a  segmentation fault  that
!   could be fixed by adding the flag "-heap-arrays" at compilation.
!-----------------------------------------------------------------------
! INPUT/OUTPUT FILE UNIT NUMBERS:
! 50  output new file (used only in subroutine NEWFILE)
! 60  input parameter file (used only in subroutine INPUT)
! 61  all input data files (used only in subroutine READDATA)
! 62  output parameter file
! 63  output light curve data file
! 64  output best-fit file
! 65  output file for RVs of star A
! 66  output file for RVs of star B
! 70-79  used for debugging output
!-----------------------------------------------------------------------
! TYPES OF DATA:
! DTYPE=1 for light curve data points
! DTYPE=2 for light ratio
! DTYPE=3 for times of minimum light
! DTYPE=4 for third light
! DTYPE=5 for e*cos(omega) or e
! DTYPE=6 for e*sin(omega) or omega
! DTYPE=7 for radial velocities (RVs) of star A
! DTYPE=8 for radial velocities of star B
!-----------------------------------------------------------------------
! TYPES OF LIMB DARKENING:
! LDTYPE=1: linear LD law
! LDTYPE=2: quadratic LD law
! LDTYPE=3: logarithmic LD law
! LDTYPE=4: square-root LD law
! LDTYPE=5: cubic LD law
! LDTYPE=6: power-2 law
! LDTYPE=9: Claret four-parameter LD law
! LDTYPE=-2: quadratic LD law using u+0.6v and u-0.6v
! LDTYPE=-3: logarithmic LD law using u+0.6v and u-0.6v
! LDTYPE=-4: square-root LD law using u+0.6v and u-0.6v
! LDTYPE=-5: cubic LD law using u+0.3v and u-0.3v
! LDTYPE=-6: power-2 law re-parameterised into h1 and h2
! LDTYPE=-1: star 2 has same LD as star 1 (possible for LDTYPE(2) only)
!=======================================================================
!=======================================================================
      PROGRAM JKTEBOP
      implicit none

            ! It is now possible to specify the  maximum number of data-
            ! points you want the code to deal with. The advantage of a
            ! small number of datapoints is that the code requires less
            ! memory and may run slightly faster.   The disadvantage is
            ! that you cannot fit large datasets.  I recommend that you
            ! opt for either 99,999 or 999,999 datapoints.  Specify the
            ! number using the "parameter" statement below. The largest
            ! possible MAXDATA probably depends on the Fortran compiler
            ! you are using.
            ! If you play around with MAXDATA then you should also modi-
            ! fy the DATAFORM parameter as well. DATAFORM specifies the
            ! number of characters required to output NDATA (the number
            ! of datapoints) to screen or to file. The largest value of
            ! NDATA is limited by MAXDATA.   For MAXDATA=99999 you want
            ! DATAFORM='5',   whereas for MAXDATA=999999 you would need
            ! DATAFORM='6' or greater.
            ! The reason to make this easy to specify  is that the user
            ! has more control over the exact spacing of the diagnostic
            ! output,  which may be useful if that output is being read
            ! by another code  which needs numbers to be in a  specific
            ! place in the output.
      integer MAXDATA
      parameter (MAXDATA = 999999)
      character*1 DATAFORM
      parameter (DATAFORM = '6')

      integer VERSION               ! Version of this code
      integer TASK                  ! Task to perform (between 0 and 5)
      character INFILE*50           ! Name of the input parameter file
      real*8  V(540)                ! Fundamental EBOP model parameters
      integer VARY(540)             ! Params fixed (0) or variable (1)
      integer LDTYPE(2)             ! Type of LD law for each star
      real*8  DATX(MAXDATA)         ! Data abscissa (time or phase)
      real*8  DATY(MAXDATA)         ! Data ordinate (measured value)
      real*8  DATERR(MAXDATA)       ! Data errorbars
      integer DTYPE(MAXDATA)        ! Type of data (between 1 and 6)
      integer NDATA                 ! Number of data points
      integer NMIN                  ! Number of times of minimum light
      integer NLR                   ! Number of observed light ratios
      integer NL3                   ! Number of obsd third light values
      integer NECW                  ! Number of observed e*cos(omega)'s
      integer NESW                  ! Number of observed e*sin(omega)'s
      integer NSINE                 ! Number of sine curves to include
      integer PSINE(20)             ! Which parameter each sine acts on
      integer NPOLY                 ! Number of polynomials to include
      integer PPOLY(50)             ! Which parameter each poly acts on
      integer NSIM                  ! Number of simulations or refits
      real*8  SIGMA                 ! Sigma rejection value for task 4
      integer NUMINT                ! Number of numerical integrations
      real*8  NINTERVAL             ! Time interval for numerical ints.
      integer i,ERROR               ! Loop counter and error flag

      VERSION = 43

      ERROR = 0
      do i = 1,MAXDATA
        DATX(i) = 0.0d0
        DATY(i) = 0.0d0
        DATERR(i) = 0.0d0
        DTYPE(i) = 0
      end do
      do i = 1,540
        V(i) = -100.0d0
        VARY(i) = 0
      end do
      V(27) = -1.0d10               ! Need very low values to avoid ...
      V(28) = -1.0d10               ! ... ignoring negative values.
      NDATA = 0
      NMIN = 0
      NLR = 0
      NL3 = 0
      NECW = 0
      NESW = 0
      NSIM = 0
      NSINE = 0
      NPOLY = 0
      do i = 1,20
        PSINE(i) = 0
      end do
      do i = 1,50
        PPOLY(i) = 0
      end do
      LDTYPE(1) = -100
      LDTYPE(2) = -100
      NUMINT = 1
      NINTERVAL = 0.0d0

      V(14) = 0.0d0                 ! Tidal angle (leave at 0 normally)
      VARY(14) = 0                  ! Don't vary tidal lead/lag angle
      VARY(18) = 0                  ! Don't vary integration ring size


            ! First output the code name and version, then check for
            ! command-line arguments. If one is found, then take it to
            ! be the input file name.

      write(6,*) " "
      write(6,'(A10,I2,A13,A55)') "JKTEBOP  v",VERSION,"  John Southw",
     &        "orth (Keele University, UK, taylorsouthworth~gmail.com)"

      if ( iargc() == 1 ) then
        CALL GETARG (1,INFILE)
      else
        write(6,'(A39,A41)') "A package for modelling the light curve",
     &                     "s of well-detached eclipsing binary stars"
        write(6,'(A39,A41)') "Task 1  outputs limb and gravity darken",
     &                     "ing coefficients for given Teff and log g"
        write(6,'(A39,A41)') "Task 2  outputs one model light curve c",
     &                     "alculated using a set of input parameters"
        write(6,'(A39,A41)') "Task 3  finds the best fit of the model",
     &                     " to observations  (formal errorbars only)"
        write(6,'(A39,A41)') "Task 4  finds the best fit to the obser",
     &                     "vations, sigma clips, and refits the data"
        write(6,'(A39,A41)') "Task 5  finds global best fit, by pertu",
     &                     "rbing parameters and refitting many times"
        write(6,'(A39,A41)') "Task 6  fits observations and finds goo",
     &                     "dness of fit for several parameter values"
        write(6,'(A39,A41)') "Task 7  finds robust reliable errors by",
     &                     " analysing with a bootstrapping algorithm"
        write(6,'(A39,A41)') "Task 8  finds robust errors by analysin",
     &                     "g with a Monte Carlo simulation algorithm"
        write(6,'(A39,A41)') "Task 9  finds robust errors using Monte",
     &                     " Carlo given significant correlated noise"
        write(6,*) " "
        write(6,'(A39,A41)') "Usage:  'jktebop  [inputfile]' to under",
     &                     "take one of these tasks                  "
        write(6,'(A39,A41)') "Usage:  'jktebop   newfile'    to outpu",
     &                     "t an empty input file for a given task   "
        write(6,'(A39,A41)') "Usage:  'jktebop     1'        to under",
     &                     "take Task 1 (limb darkening coefficients)"
        write(6,*) " "
        stop
      end if

            ! Check the MAXDATA and DATAFORM values. These are the two
            ! parameters which are intended to be changeable by users.

      if ( MAXDATA < 999 ) then
        write(6,'(A39,A41)') "## Warning: maximum number of datapoint",
     &                     "s is less than 999, which is very small. "
      end if
      if ( MAXDATA > 9999999 ) then
        write(6,'(A39,A41)') "## Warning: maximum number of datapoint",
     &                     "s is more than 9,999,999, which is huge. "
      end if

      if ( DATAFORM /= '1' .and. DATAFORM /= '2' .and. DATAFORM /= '3'
     &     .and. DATAFORM /= '4' .and. DATAFORM /= '5' .and.
     &     DATAFORM /= '6' .and. DATAFORM /= '7' .and. DATAFORM /= '8'
     &     .and. DATAFORM /= '9' ) then
        write (6,'(A39,A41)') "### ERROR: parameter DATAFORM must be '",
     &                      "1' '2' '3' '4' '5' '6' '7' '8' or '9'.   "
        write(6,*) " "
        stop
      end if

      if ( MAXDATA > 99 .and. (DATAFORM =='1' .or. DATAFORM=='2') ) then
        write(6,'(A39,A41)') "## Warning: DATAFORM is too small for M",
     &                     "AXDATA. Some format statements may break."
      end if
      if ( MAXDATA > 999 .and. (DATAFORM == '1' .or. DATAFORM == '2'
     &                                      .or. DATAFORM == '3') ) then
        write(6,'(A39,A41)') "## Warning: DATAFORM is too small for M",
     &                     "AXDATA. Some format statements may break."
      end if
      if ( MAXDATA > 9999 .and. (DATAFORM == '1' .or. DATAFORM == '2'
     &                 .or. DATAFORM == '3' .or. DATAFORM == '4') ) then
        write(6,'(A39,A41)') "## Warning: DATAFORM is too small for M",
     &                     "AXDATA. Some format statements may break."
      end if
      if ( MAXDATA > 99999 .and. (DATAFORM == '1' .or. DATAFORM == '2'
     &  .or. DATAFORM=='3' .or. DATAFORM=='4' .or. DATAFORM=='5') ) then
        write(6,'(A39,A41)') "## Warning: DATAFORM is too small for M",
     &                     "AXDATA. Some format statements may break."
      end if
      if ( MAXDATA > 999999 .and. (DATAFORM == '1' .or. DATAFORM == '2'
     &                 .or. DATAFORM == '3' .or. DATAFORM == '4'
     &                 .or. DATAFORM == '5' .or. DATAFORM == '6') ) then
        write(6,'(A39,A41)') "## Warning: DATAFORM is too small for M",
     &                     "AXDATA. Some format statements may break."
      end if
      if ( MAXDATA>9999999 .and. (DATAFORM == '1' .or. DATAFORM == '2'
     &    .or. DATAFORM == '3' .or. DATAFORM == '4' .or. DATAFORM == '5'
     &                 .or. DATAFORM == '6' .or. DATAFORM == '7') ) then
        write(6,'(A39,A41)') "## Warning: DATAFORM is too small for M",
     &                     "AXDATA. Some format statements may break."
      end if

            ! Now to the main part. Check the command-line input to see
            ! what to do, and then call the appropriate subroutines.

      if ( INFILE == "newfile" .or. INFILE == "NEWFILE" ) then
        CALL NEWFILE ()
        write(6,*) " "
        stop
      else if ( INFILE == "1" ) then
        CALL TASK1 ()
        STOP
      end if

      CALL INPUT (VERSION,TASK,INFILE,V,VARY,LDTYPE,DATX,DATY,DATERR,
     &           DTYPE,NDATA,MAXDATA,DATAFORM,NLR,NMIN,NSIM,SIGMA,NSINE,
     &           PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL,ERROR)
      if ( ERROR /= 0 ) then
        write(6,*) " "
        stop
      end if

      if ( TASK == 2 ) CALL TASK2 (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY)

      if ( TASK == 3 .or. TASK == 4 ) CALL TASK34 (TASK,V,VARY,LDTYPE,
     &          DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,NMIN,
     &     SIGMA,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      if ( TASK == 6 ) CALL TASK6 (V,VARY,LDTYPE,DATX,DATY,DATERR,DTYPE,
     &                      NDATA,MAXDATA,DATAFORM,NLR,NMIN,NSINE,PSINE,
     &                       NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      if ( TASK == 5 .or. TASK == 7 .or. TASK == 8 .or. TASK == 9)
     &  CALL TASK5789 (TASK,V,VARY,LDTYPE,DATX,DATY,DATERR,DTYPE,NDATA,
     &                       MAXDATA,DATAFORM,NLR,NMIN,NSIM,NSINE,PSINE,
     &                       NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      write(6,*) " "

      END PROGRAM JKTEBOP
!=======================================================================
!=======================================================================
      SUBROUTINE NEWFILE ()                ! Makes all empty input files
      implicit none
      integer TASK                         ! Task the file is wanted for
      character*50 OUTFILE                 ! Name of output file to make
      integer ERROR                        ! Loop counter and error flag

      ERROR = 0
      write(6,'(A39,$)')   "Enter name of input file to create >>  "
      read (*,*) OUTFILE
      write(6,'(A39,$)')   "Enter task number (between 2 and 9) >> "
      read (*,*,iostat=ERROR) TASK

      if ( ERROR /= 0 ) then
        write(6,'(A40)') "### ERROR: did not understand input.    "
        return
      else if ( TASK < 2 .or. TASK > 9 ) then
        write(6,'(A40)') "### ERROR: task integer is out of range."
        return
      else
        CALL OPENFILE (50,"new"," output   ",OUTFILE,ERROR)
        if ( ERROR /= 0 ) return
      end if

      write(50,100)"Task to do (from 2 to 9)   Integ. ring size (deg)  "
      write(50,100)"Sum of fractional radii    Ratio of the radii      "
      write(50,100)"Orbital inclination (deg)  Mass ratio of system    "
      write(50,100)"ecosw or eccentricity      esinw or periastron long"
      write(50,100)"Gravity darkening (star A) Grav darkening (star B) "
      write(50,100)"Surface brightness ratio   Amount of third light   "
      write(50,100)"LD law type for star A     LD law type for star B  "
      write(50,100)"LD star A (linear coeff)   LD star B (linear coeff)"
      write(50,100)"LD star A (nonlin coeff)   LD star B (nonlin coeff)"
      write(50,100)"Reflection effect star A   Reflection effect star B"
      write(50,100)"Phase of primary eclipse   Light scale factor (mag)"

      if ( TASK >= 3 .and. TASK <= 9 ) then
        write(50,101)"Orbital period of eclipsing binary system (days) "
        write(50,101)"Reference time of primary minimum (JD/HJD/BJD)   "
      end if

      if ( TASK == 2 ) then
        write(50,101)"Output file name (continuous character string)   "
      else if ( TASK == 4 ) then
        write(50,101)"Sigma value to reject discrepant observations    "
      else if ( TASK == 5 ) then
        write(50,101)"Number of refits for perturbed initial parameters"
      else if ( TASK == 7 ) then
        write(50,101)"Number of bootstrapping simulations to do        "
      else if ( TASK == 8 ) then
        write(50,101)"Number of Monte Carlo simulations to do          "
      end if

      if ( TASK >= 3 .and. TASK <= 9 ) then
      write(50,103)"Adjust RADII SUM    or  RADII RATIO    (0, 1, 2, 3)"
      write(50,103)"Adjust INCLINATION  or  MASSRATIO      (0, 1, 2, 3)"
      write(50,103)"Adjust ECOSW-or-E   or  ESINW-or-OMEGA (0, 1, 2, 3)"
      write(50,103)"Adjust GRAVDARK1    or  GRAVDARK2      (0, 1, 2, 3)"
      write(50,103)"Adjust SURFBRIGHT2  or  THIRDLIGHT     (0, 1, 2, 3)"
      write(50,103)"Adjust LD-LIN starA or  LD-LIN starB   (0, 1, 2, 3)"
      write(50,103)"Adjust LD-NONLIN A  or  LD-NONLIN B    (0, 1, 2, 3)"
      write(50,103)"Adjust REFLECTION A or  REFLECTION B   (-1,0,1,2,3)"
      write(50,103)"Adjust PHASESHIFT   or  SCALE FACTOR   (0, 1, 2, 3)"
      write(50,103)"Adjust PERIOD       or  T(pri.ecl.)    (0, 1, 2, 3)"
      write(50,100)"Name of file containing light curve                "
      write(50,100)"Name of output parameter file                      "
      end if

      if ( TASK == 3 .or. TASK == 4 ) then
        write(50,101)"Name of output light curve file                  "
        write(50,101)"Name of output model light curve fit file        "
      end if

      if ( TASK == 5 .or. TASK == 7 .or. TASK == 8 .or. TASK == 9 ) then
        write(50,101)"Name of output file of individual fits           "
      end if

      write (50,*) " "
      write (50,*) " "

      write (50,102) "# Enter the appropriate numbers on the l",
     &               "eft-hand side of each line of this file."
      write (50,102) "# Most of the lines require two numerica",
     &               "l parameters separated by spaces.       "
      write (50,*) " "

      write (50,102) "# Put a negative number for the mass rat",
     &               "io to force the stars to be spherical.  "
      write (50,102) "# The mass ratio will then be irrelevant",
     &               " (it is only used to get deformations). "
      write (50,*) " "

      if ( TASK == 2 ) then
        write (50,102) "# The reflection coefficients will be ta",
     &                 "ken from the input file. If you want the"
        write (50,102) "# code to calculate them itself (which w",
     &                 "ill be the usual situation), put a large"
        write (50,102) "# negative number (e.g. -100) for each r",
     &                 "eflection coefficient in the input file."
        write (50,*) " "
      end if

      write (50,102) "# To fit for rA and rB instead of (rA+rB",
     &               ") and k, give a negative value for r1+r2"
      write (50,102) "# Then (rA+rB) will be interpreted to me",
     &               "an rA,  and k will be interpreted as rB."
      write (50,102) "# The adjustment indicators will similar",
     &               "ly refer to rA,rB rather than (rA+rB),k."
      write (50,*) " "

      write (50,102) "# If e < 10 then e and omega will be ass",
     &               "umed to be e*cos(omega) and e*sin(omega)"
      write (50,102) "# If e >= 10 then e and omega will be as",
     &               "sumed to be (e+10) and omega (degrees). "
      write (50,102) "# The first option is in general better ",
     &               "unless eccentricity is larger or fixed. "
      write (50,*) " "

      write (50,102) "# The possible entries for the type of l",
     &               "imb darkening law are 'lin' (for linear)"
      write (50,102) "# 'quad' (quadratic), 'log' (logarithmic",
     &               "), 'sqrt' (square-root), 'cub' (cubic), "
      write (50,102) "# 'pow2' (power-2), qrep (quadratic re-p",
     &               "arameterisation) lrep (logarithmic re-  "
      write (50,102) "# parameterisation, srep (square-root re",
     &               "-parameterisation, crep (cubic re-param-"
      write (50,102) "# eterisation, and 'h1h2' (power-2 with ",
     &               "the h1 and h2 parameterisation).        "
      write (50,102) "# Put 'same' for star B to force its coe",
     &               "fficients to be equal those of star A.  "
      write (50,*) " "

      if ( TASK >= 3 .and. TASK <= 9 ) then
        write (50,102) "# For each adjustable parameter the adju",
     &                 "stment integer can be 0  (parameter will"
        write (50,102) "# be fixed at the input file value),  1 ",
     &                 "(parameter will be freely adjusted),  2 "
        write (50,102) "# (parameter will be fixed for initial f",
     &                 "it but perturbed during later analysis)."
        write (50,102) "# or  3 (adjusted in initial fit but not",
     &                 " perturbed during Monte Carlo analysis)."
        write (50,*) " "
        write (50,102) "# When fitting a light curve  the reflec",
     &                 "tion coefficients can be calculated from"
        write (50,102) "# the system geometry  (put -1 for the a",
     &                 "djustment integers),  held fixed (put 0)"
        write (50,102) "# or freely adjusted to fit the light cu",
     &                 "rve (put 1) - useful for close binaries."
        write (50,*) " "
      end if

      if ( TASK == 7 .or. TASK == 8 ) then
        write (50,102) "# When doing Monte Carlo or bootstrappin",
     &                 "g the starting parameters for each simu-"
        write (50,102) "# lation are perturbed by a pre-defined ",
     &                 "amount to avoid biassing the results. If"
        write (50,102) "# this is not wanted, put a minus sign i",
     &                 "n front of the number of simulations.   "
        write (50,*) " "
      end if

      write (50,102) "# FOUR-PARAMETER LIMB DARKENING:   you c",
     &               "an alternatively put '4par' for the limb"
      write (50,102) "# darkening in which case the input file",
     &               " format differs a bit. Change the lines:"
      write (50,102) "#                   LD star A (linear co",
     &               "eff)   LD star B (linear coeff)         "
      write (50,102) "#                   LD star A (nonlin co",
     &               "eff)   LD star B (nonlin coeff)         "
      write (50,102) "# to the following lines and put in the ",
     &               "information at the line starts as usual:"
      write (50,102) "#                   LD star A (coefficie",
     &               "nt 1)  LD star B (coefficient 1)        "
      write (50,102) "#                   LD star A (coefficie",
     &               "nt 2)  LD star B (coefficient 2)        "
      write (50,102) "#                   LD star A (coefficie",
     &               "nt 3)  LD star B (coefficient 3)        "
      write (50,102) "#                   LD star A (coefficie",
     &               "nt 4)  LD star B (coefficient 4)        "
      write (50,102) "# You also need to change the lines for ",
     &               "the adjustment parameters from:         "
      write (50,102) "#                   Adjust LD-lin1  or  ",
     &               "LD-lin2            (0, 1, 2, 3)         "
      write (50,102) "#                   Adjust LD-nonlin1  o",
     &               "r  LD-nonlin2      (0, 1, 2, 3)         "
      write (50,102) "# to the following lines and put in the ",
     &               "information at the line starts as usual:"
      write (50,102) "#                   Adjust LDcoeff-A1   ",
     &               "or  LDcoeff-B1     (0, 1, 2, 3)         "
      write (50,102) "#                   Adjust LDcoeff-A2   ",
     &               "or  LDcoeff-B2     (0, 1, 2, 3)         "
      write (50,102) "#                   Adjust LDcoeff-A3   ",
     &               "or  LDcoeff-B3     (0, 1, 2, 3)         "
      write (50,102) "#                   Adjust LDcoeff-A4   ",
     &               "or  LDcoeff-B4     (0, 1, 2, 3)         "
      write (50,102) "# Remember not to include the '#' symbol",
     &               ": it is used only to comment lines out. "
      write (50,102) "# Reference for the '4par' law: Claret (",
     &               "2000, A&A, 363, 1081)                   "
      write (50,*) " "

      write (50,102) "# TIMES OF MINIMUM LIGHT: add a line bel",
     &               "ow the parameter line to input each one:"
      write (50,102) "#   TMIN  [cycle]  [time]  [error]      ",
     &               "                                        "
      write (50,102) "# where [cycle] is cycle number (integer",
     &               " for primary minimum  or integer+0.5 for"
      write (50,102) "# secondary minimum), [time] and [error]",
     &               " are the observed time and uncertainty. "
      write (50,*) " "

      write (50,102) "# LIGHT RATIO: add a line below the para",
     &               "meter line to input each observed one:  "
      write (50,102) "#   LRAT  [time]  [light_ratio]  [error]",
     &               "                                        "
      write (50,102) "# where [time] is the time when the spec",
     &               "troscopic light ratio was measured,     "
      write (50,102) "# [light_ratio] is its value and [error]",
     &               " is its measurement uncertainty.        "
      write (50,*) " "

      write (50,102) "# MEASURED THIRD LIGHT VALUE:   include ",
     &               "as observed constraint by adding a line:"
      write (50,102) "#   THDL  [value]  [uncertainty]        ",
     &               "                                        "
      write (50,102) "# which gives the third light measuremen",
     &               "t and its observational uncertainty.    "
      write (50,*) " "

      write (50,102) "# MEASURED orbital shape parameters (dep",
     &               "ending on the value of eccentricity):   "
      write (50,102) "#   ECSW  [value]  [uncertainty]   (inte",
     &               "rpreted as either e*cos(omega) or e)    "
      write (50,102) "#   ENSW  [value]  [uncertainty]   (inte",
     &               "rpreted as either e*sin(omega) or omega)"
      write (50,*) " "

      write (50,102) "# SINE AND POLYNOMIAL FITTING:   the par",
     &               "ameters of sine curves or polynomials of"
      write (50,102) "# order 5) can be included.  You can hav",
     &               "e up to nine sines and five polynomials,"
      write (50,102) "# each acting on a specific parameter. T",
     &               "he information for each one is specified"
      write (50,102) "# by an additional line below the main i",
     &               "nput file parameters.     Line format:  "
      write (50,102) "#   SINE  [par]  [T0]  [P]  [amp]  [vary",
     &               "(T0)]  [vary(P)]  [vary(amp)]           "
      write (50,'(A40,A40,A40,A40)')
     &               "#   POLY  [par]  [pivot]  [const]  [x]  ",
     &               "[x^2]  [x^3] [x^4]  [x^5]  [vary(const)]",
     &               "  [vary(x)]  [vary(x^2)]  [vary(x^3)]   ",
     &               "[vary(x^4)]   [vary(x^5)]               "
      write (50,102) "# where the required parameters are give",
     &               "n inside square brackets. [T0] is a time"
      write (50,102) "# of zero phase (JD/HJD/BJD), [P] is per",
     &               "iod (days), [amp] is amplitude and [x^n]"
      write (50,102) "# are the coefficients of the polynomial",
     &               ".  Each parameter has a [vary()] which  "
      write (50,102) "# is 0, 1, 2 or 3 to indicate how the pa",
     &               "rameter is treated.                     "
      write (50,102) "# [par] indicates the parameter to apply",
     &               " it to: J, r1, r2, i, L3, sf, L1, L2, e,"
      write (50,102) "# w, ec, es, T0, KA, KB, VA, VB. Here sf",
     &               " is the scale factor, L1 is the light   "
      write (50,102) "# from star A, L2 is the light from star",
     &               " B, ec is ecosw, es is esinw, T0 is the "
      write (50,102) "# reference time of minimum, KA and KB a",
     &               "re the velocity amplitudes of the two   "
      write (50,102) "# stars, and VA and VB are the systemtic",
     &               " velocities of the two stars.           "
      write (50,102) "# Note that the independent parameter is",
     &               " always time (JD/HJD/BJD or phase).     "
      write (50,102) "# If you want to apply a polynomial to o",
     &               "nly part of the data, in a specific time"
      write (50,102) "# interval, then use the following line ",
     &               "to give the extra information:          "
      write (50,'(A40,A40,A40,A40)')
     &               "#   POLY  [par]  [pivot]  [const] [x] [x",
     &               "^2] [x^3] [x^4] [x^5]  [vary(const)] [va",
     &               "ry(x)] [vary(x^2)] [vary(x^3)] [vary(x^4",
     &               ")] [vary(x^5)]  [start-time]  [end-time]"
      write (50,102) "# JKTEBOP will check for the two extra n",
     &               "umbers and use them automatically.      "
      write (50,*) " "

      write (50,102) "# NUMERICAL INTEGRATION:  long exposure ",
     &               "times can be split up into NUMINT points"
      write (50,102) "# occupying a total time interval of NIN",
     &               "TERVAL (seconds) by including this line:"
      write (50,102) "#   NUMI  [numint]  [ninterval]         ",
     &               "                                        "
      write (50,*) " "

      write (50,102) "# FITTING FOR RADIAL VELOCITIES:    the ",
     &               "observed RVs should be in separate files"
      write (50,102) "# for the two stars and the data should ",
     &               "be in the same format as the light curve"
      write (50,102) "# data. Then add a line below the main i",
     &               "nput parameters for each RV file:       "
      write (50,102) "#   RV1  [infile]  [outfile]  [K]  [Vsys",
     &               "]  [vary(K)]  [vary(Vsys)]              "
      write (50,102) "#   RV2  [infile]  [outfile]  [K]  [Vsys",
     &               "]  [vary(K)]  [vary(Vsys)]              "
      write (50,102) "# where RV1 is for primary star velociti",
     &               "es, RV2 is for secondary star velocities"
      write (50,102) "# [infile] is the input data file, [outf",
     &               "ile] is the output data file, [K] is the"
      write (50,102) "# velocity amplitude of the star (km/s),",
     &               " [Vsys] is its systemic velocity (km/s),"
      write (50,102) "# and [vary(K)] and [vary(Vsys)] are 0 t",
     &               "o fix and 1 to fit for these quantities."
      write (50,102) "# The mass ratio parameter is not used f",
     &               "or the RVs, only for the light curve.   "
      write (50,102) "# If you want to fix the systemic veloci",
     &               "ty for star B to that for star A, simply"
      write (50,102) "# set vary(Vsys) for star B to be equal ",
     &               "to -1                                   "
      write (50,*) " "

      write (50,102) "# ERROR BARS IN THE DATAFILES: whenever ",
     &               "JKTEBOP reads in a set of data from file"
      write (50,102) "# it checks to see if there are three nu",
     &               "mbers on the first line.  If so, JKTEBOP"
      write (50,102) "# assumes that the datafile contains thr",
     &               "ee columns (TIME, OBSERVATION, ERRORBAR)"
      write (50,102) "# and reads the data in accordingly. If ",
     &               "it can only find two numbers, it assumes"
      write (50,102) "# that these represent TIME, OBSERVATION",
     &               "  and that error bars are not available."
      write (50,102) "# If errorbars are not available, or if ",
     &               "they are too large or too small, JKTEBOP"
      write (50,102) "# can iteratively scale them until a red",
     &               "uced chi-squared of 1.0 is obtained.  To"
      write (50,102) "# use this option put the word 'chif' on",
     &               " a line on its own below the main set of"
      write (50,102) "# parameters.  Warning: this is availabl",
     &               "e only for light curves with ten or more"
      write (50,102) "# datapoints and for RV curves with five",
     &               " or more datapoints. Be careful!        "
      write (50,*) " "

      write (50,102) "# Specify 'nofit' at the start of a line",
     &               " to fix all the parameters at the input "
      write (50,102) "# values. This can be useful for debuggi",
     &               "ng when a fit is poor or fails.         "
      write (50,*) " "

      write (50,*) " "

      close (50,status="keep")
      write(6,'(A29,I1,A21,A50)') "An empty input file for task ",TASK,
     &                                   " has been written to ",OUTFILE

100   FORMAT (18X,A51)
101   FORMAT (18X,A49)
102   FORMAT (A40,A40)
! 1021  FORMAT (A40,A40,A40,A40)
103   FORMAT (" 0  0",13X,A51)

      END SUBROUTINE NEWFILE
!=======================================================================
!=======================================================================
      SUBROUTINE INPUT (VERSION,TASK,INFILE,V,VARY,LDTYPE,DATX,DATY,
     &                  DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,NMIN,
     &                  NSIM,SIGMA,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,
     &                  NESW,NUMINT,NINTERVAL,STATUS)

            ! This subroutine reads in all the input file and input data
            ! It opens several files for reading later, but the names of
            ! the files are not retained. See top for file unit numbers.

      implicit none
      integer MAXDATA                     ! IN: max number of datapoints
      character DATAFORM*1                ! IN: MAXDATA character length
      integer VERSION,TASK                ! IN:  Code version Task numbr
      character INFILE*50                 ! IN:  Name of the  input file
      real*8 V(540)                       ! OUT: Photometric  parameters
      integer VARY(540)                   ! OUT: Par adjustment ingeters
      integer LDTYPE(2)                   ! OUT: Type of LD law to adopt
      real*8 DATX(MAXDATA)                ! OUT: Data independent varble
      real*8 DATY(MAXDATA)                ! OUT: Data dependent variable
      real*8 DATERR(MAXDATA)              ! OUT: Data errorbars
      integer DTYPE(MAXDATA)              ! OUT: Type of data (1 2 or 3)
      integer NDATA,NLR,NMIN              ! OUT: Amount of types of data
      integer NSIM                        ! OUT: Number  of  simulations
      real*8 SIGMA                        ! OUT: Sigma  rejection  value
      integer NSINE,NL3                   ! OUT: Number of sines and L3s
      integer NECW,NESW                   ! OUT: Number of ecosw & esinw
      integer PSINE(20)                   ! OUT: Which par for each sine
      integer NPOLY,PPOLY(50)             ! OUT: Similar for polynomials
      integer NUMINT                      ! OUT: Number of numerical int
      real*8  NINTERVAL                   ! OUT: Time interval of numint
      character OBSFILE*50                ! LOCAL: Input lightcurve file
      character LCFILE*50                 ! LOCAL: Output lghtcurve file
      character PARFILE*50                ! LOCAL: Output parameter file
      character FITFILE*50                ! LOCAL: Output model fit file
      character RV1OBSFILE*50             ! LOCAL: Input RV file, star A
      character RV2OBSFILE*50             ! LOCAL: Input RV file, star B
      character RV1OUTFILE*50             ! LOCAL: Output RV file star A
      character RV2OUTFILE*50             ! LOCAL: Output RV file star B
      integer i,j,k,ERROR,STATUS          ! LOCAL: Counter & error flags
      real*8 LP,LS                        ! LOCAL: Light from  each star
      real*8 ARRAY(MAXDATA),SELLECT       ! LOCAL: Help array & FUNCTION
      character DTDATE*8,DTTIME*10        ! LOCAL: runtime time and date
!!      character DTDATE*8,DTTIME*9    ! form needed for g77 compilation
      character LD1*4,LD2*4               ! LOCAL: type of LD law to use
      character CHARHELP200*200           ! LOCAL: variable to read in
      character CHARHELP4*4               ! LOCAL: variable to read in
      real*8 GETMODEL,MAG                 ! LOCAL: function, value ret'd
      character WHICHPAR*2                ! LOCAL: par name for sin/poly
      integer IARRAY(6),IHELP             ! LOCAL: for reading stuff in
      real*8 R1,R2                        ! LOCAL: frac radii of stars
      integer CHIFUDGE                    ! LOCAL: indicator for option
      real*8 POLYPIV,POLYSTART,POLYSTOP   ! LOCAL: polynomial parameters
      real*8 POLYCONST,POLYX1,POLYX2      ! LOCAL: polynomial parameters
      real*8 POLYX3,POLYX4,POLYX5         ! LOCAL: polynomial parameters
      integer NDATA_LC,NDATA_NOTLC        ! LOCAL: number of datapoints
      integer NOFIT                       ! LOCAL: flag to fit nothing

      CHIFUDGE = 0                        ! LOCAL: to force chi^2_nu = 1
      NOFIT = 0                           ! LOCAL: set to 0 to ignore

      ERROR = 0
      STATUS = 0
      CALL OPENFILE (60,"old","input file",INFILE,ERROR)
      if ( ERROR /= 0 ) then
        write(6,*) " "
        STOP
      end if

!-----------------------------------------------------------------------

      read (60,*,iostat=ERROR) TASK,V(18)
      if ( ERROR /= 0 ) then
        write(6,'(A39,A41)') "### ERROR: cannot read first line of in",
     &                     "put parameter file (TASK and INTRING).   "
        STATUS = 1
        write(6,'(A28,I5)') "### Error message returned: ",ERROR
      else if ( TASK < 0 .or. TASK > 9 ) then
        write(6,'(A39,A41)') "### ERROR: task integer does not corres",
     &                     "pond to a value between 1 and 9.         "
        STATUS = 1
      else if ( V(18) < 0.1d0 .or. V(18) > 10.0d0 ) then
        if ( V(18) < -1.01d0 .or. V(18) > -0.99d0 ) then
          write(6,'(A37,A43)') "### ERROR: integration ring size must",
     &                   " be between 0.1 and 10.0 degrees.          "
          STATUS = 1
        end if
      end if
      if ( STATUS /= 0 ) return

      if ( TASK == 1 ) write(6,'(A25,A55)')"Task 1  outputs limb and ",
     &        "gravity darkening coefficients for given Teff and log g"
      if ( TASK == 2 ) write(6,'(A25,A55)')"Task 2  outputs one model",
     &        " light curve calculated using a set of input parameters"
      if ( TASK == 3 ) write(6,'(A25,A55)')"Task 3  finds the best fi",
     &        "t of the model to observations  (formal errorbars only)"
      if ( TASK == 4 ) write(6,'(A25,A55)')"Task 4  finds the best fi",
     &        "t to the observations, sigma clips, and refits the data"
      if ( TASK == 5 ) write(6,'(A25,A55)')"Task 5  finds global best",
     &        " fit, by perturbing parameters and refitting many times"
      if ( TASK == 6 ) write(6,'(A25,A55)')"Task 6  fits observations",
     &        " and finds goodness of fit for several parameter values"
      if ( TASK == 7 ) write(6,'(A25,A55)')"Task 7  finds robust reli",
     &        "able errors by analysing with a bootstrapping algorithm"
      if ( TASK == 8 ) write(6,'(A25,A55)')"Task 8  finds robust erro",
     &        "rs by analysing with a Monte Carlo simulation algorithm"
      if ( TASK == 9 ) write(6,'(A25,A55)')"Task 9  finds robust erro",
     &        "rs using Monte Carlo given significant correlated noise"

      CALL READFF (60,"RADII SUM ",-0.8d0, 0.8d0, V( 2),
     &                "RADIIRATIO", 0.0d0, 1.0d2, V( 3),STATUS)

      if ( V(3) < 0.0 ) then
        write (6,'(A39,A41)') "### Warning: the ratio of the radii is ",
     &                      "less than zero. If you are trying to use "
        write (6,'(A39,A41)') "### the option of inputting rA and rB r",
     &                      "ather than (rA+rB) and k=rB/rA, then you "
        write (6,'(A39,A41)') "### need to put a minus sign in front o",
     &                      "f the (rA+rB) value, not the value for k."
      end if

      CALL READFF (60,"INCLNATION",50.0d0, 1.4d2, V( 6),
     &                "MASS RATIO",-1.0d9, 1.0d3, V(13),STATUS)
      CALL READFF (60,"ECCENTRCTY",-1.0d0,11.0d0, V( 7),
     &                "OMEGA     ",-3.6d2, 3.6d2, V( 8),STATUS)
      CALL READFF (60,"GRAVDARK-A",-1.0d1, 1.0d1, V( 9),
     &                "GRAVDARK-B",-1.0d1, 1.0d1, V(10),STATUS)
      CALL READFF (60,"SURF-BRT-B", 0.0d0, 1.0d4, V( 1),
     &                "THIRDLIGHT", 0.0d0, 1.0d1, V(15),STATUS)

!------------------------------------------------------------------------

      read (60,*) LD1,LD2

      if ( LD1 == "lin"  ) LDTYPE(1) = 1
      if ( LD1 == "quad" ) LDTYPE(1) = 2
      if ( LD1 == "log"  ) LDTYPE(1) = 3
      if ( LD1 == "sqrt" ) LDTYPE(1) = 4
      if ( LD1 == "cub"  ) LDTYPE(1) = 5
      if ( LD1 == "pow2" ) LDTYPE(1) = 6
      if ( LD1 == "4par" ) LDTYPE(1) = 9
      if ( LD1 == "qrep" ) LDTYPE(1) = -2
      if ( LD1 == "lrep" ) LDTYPE(1) = -3
      if ( LD1 == "srep" ) LDTYPE(1) = -4
      if ( LD1 == "crep" ) LDTYPE(1) = -5
      if ( LD1 == "h1h2" ) LDTYPE(1) = -6

      if ( LD2 == "lin"  ) LDTYPE(2) = 1
      if ( LD2 == "quad" ) LDTYPE(2) = 2
      if ( LD2 == "log"  ) LDTYPE(2) = 3
      if ( LD2 == "sqrt" ) LDTYPE(2) = 4
      if ( LD2 == "cub"  ) LDTYPE(2) = 5
      if ( LD2 == "pow2" ) LDTYPE(2) = 6
      if ( LD2 == "4par" ) LDTYPE(2) = 9
      if ( LD2 == "qrep" ) LDTYPE(2) = -2
      if ( LD2 == "lrep" ) LDTYPE(2) = -3
      if ( LD2 == "srep" ) LDTYPE(2) = -4
      if ( LD2 == "crep" ) LDTYPE(2) = -5
      if ( LD2 == "h1h2" ) LDTYPE(2) = -6
      if ( LD2 == "same" ) LDTYPE(2) = -1

      if ( LDTYPE(1) < -98 ) then
        write(6,'(A39,A41)') "### ERROR: limb darkening law for star ",
     &                     "A should be one of: lin, quad, log, sqrt,"
        write(6,'(A39,A41)') "### cub, pow2, 4par, qrep, lrep, srep, ",
     &                     "crep, or h1h2.                           "
        STATUS = 1
      end if
      if ( LDTYPE(2) < -98 ) then
        write(6,'(A39,A41)') "### ERROR: limb darkening law for star ",
     &                     "B should be one of: lin, quad, log, sqrt,"
        write(6,'(A39,A41)') "### cub, pow2, 4par, qrep, lrep, srep, ",
     &                     "crep, h1h2, or sam.                      "
        STATUS = 1
      end if
      if ( STATUS /= 0 ) return

      if ( LDTYPE(1)==9 .and. (LDTYPE(2)/=9 .and. LDTYPE(2)/=-1) ) then
        write(6,'(A39,A41)') "### ERROR: star A has four-parameter l ",
     &                     "imb darkening so star B must have too.   "
        STATUS = 1
      end if
      if ( LDTYPE(1) /= 9 .and. LDTYPE(2) == 9 ) then
        write(6,'(A39,A41)') "### ERROR: star B has four-parameter l ",
     &                     "imb darkening so star A must have too.   "
        STATUS = 1
      end if
      if ( STATUS /= 0 ) return

      if ( LDTYPE(1) /= 9 ) then
        CALL READFF (60,"LD-lin-A  ",-1.0d0, 2.0d0, V( 4),
     &                  "LD-lin-B  ",-1.0d0, 2.0d0, V( 5),STATUS)
        if ( LDTYPE(2) == -1 ) V(5) = V(4)
        CALL READFF (60,"LDnonlin-A",-1.0d0, 2.0d0, V(21),
     &                  "LDnonlin-B",-1.0d0, 2.0d0, V(24),STATUS)
        if ( LDTYPE(2) == -1 ) V(24) = V(21)
      end if

      if ( LDTYPE(1) == 9 ) then
        CALL READFF (60,"LDC-A-1   ",-9.0d0, 9.0d0, V( 4),
     &                  "LDC-B-1   ",-9.0d0, 9.0d0, V( 5),STATUS)
        CALL READFF (60,"LDC-A-2   ",-9.0d0, 9.0d0, V(21),
     &                  "LDC-B-2   ",-9.0d0, 9.0d0, V(24),STATUS)
        CALL READFF (60,"LDC-A-3   ",-9.0d0, 9.0d0, V(22),
     &                  "LDC-B-3   ",-9.0d0, 9.0d0, V(25),STATUS)
        CALL READFF (60,"LDC-A-4   ",-9.0d0, 9.0d0, V(23),
     &                  "LDC-B-4   ",-9.0d0, 9.0d0, V(26),STATUS)
        if ( LDTYPE(2) == -1 ) then
          V(5) = V(4)
          V(24) = V(21)
          V(25) = V(22)
          V(26) = V(23)
        end if
      end if

!------------------------------------------------------------------------

      CALL READFF (60,"REFLECTN-A",-1.0d9, 1.0d9, V(11),
     &                "REFLECTN-B",-1.0d9, 1.0d9, V(12),STATUS)
      CALL READFF (60,"PHASESHIFT",-1.0d0, 1.0d0, V(16),
     &                "SCALEFACTR",-1.0d3, 1.0d3, V(17),STATUS)

      if ( STATUS /= 0 ) return

      if ( TASK >= 3 .and. TASK <= 9 ) then
        CALL READF (60, "PERIOD    ", 0.0d0, 1.0d6, V(19),STATUS)
        CALL READF (60, "TIME-ZERO ",-3.0d7, 3.0d7, V(20),STATUS)
      end if

      if ( TASK == 2 ) then
        CALL READCHAR30 (60,"OUTFILE   ",PARFILE,STATUS)
        close (60)
        CALL OPENFILE (62,"new","lightcurve",PARFILE,STATUS)
        return
      end if

      if ( STATUS /= 0 ) return

      if ( TASK == 5 .or. TASK == 7 .or. TASK == 8 ) then
        read (60,*,iostat=ERROR) NSIM
        if ( ERROR /= 0 ) then
          if (TASK==5) write(6,'(A24,A56)') "### ERROR reading number",
     &       " of perturbed initial parameter sets to refit data with."
          if (TASK==7) write(6,'(A24,A56)') "### ERROR reading the nu",
     &       "mber of bootstrapping simulations to do.                "
          if (TASK==8) write(6,'(A24,A56)') "### ERROR reading the nu",
     &       "mber of Monte Carlo simulations to do.                  "
          STATUS = 1
        end if
        if ( abs(NSIM) < 8 .or. abs(NSIM) > 99999 ) then
          write(6,'(A37,A43)') "### ERROR: number of simulations/sets",
     &                    " to do must be between +/-8 and +/-99999.  "
          STATUS = 1
          return
        end if
      else if ( TASK == 4 ) then
        read (60,*,iostat=ERROR) SIGMA
        if ( ERROR /= 0 ) then
          write(6,'(A37,A43)') "### ERROR reading the sigma number fo",
     &                    "r rejection of discrepant data.            "
          STATUS = 1
        end if
      end if

      if ( STATUS /= 0 ) return

      if ( TASK >= 3 .and. TASK <= 9 ) then
        CALLREAD2(60,"adj(rA+rB)","adj(rB/rA)",VARY( 2),VARY( 3),STATUS)
        CALLREAD2(60,"adj(_INC_)","adj(MB/MA)",VARY( 6),VARY(13),STATUS)
        CALLREAD2(60,"adj(_ECC_)","adj(OMEGA)",VARY( 7),VARY( 8),STATUS)
        CALLREAD2(60,"adj(_GD-A)","adj(_GD-B)",VARY( 9),VARY(10),STATUS)
        CALLREAD2(60,"adj(_SB-B)","adj(_L_3_)",VARY( 1),VARY(15),STATUS)

        if ( LDTYPE(1) /= 9 ) then
        CALLREAD2(60,"adj(LD-l1)","adj(LD-l1)",VARY( 4),VARY( 5),STATUS)
        CALLREAD2(60,"adj(LD-n2)","adj(LD-n2)",VARY(21),VARY(24),STATUS)
        else
        CALLREAD2(60,"adj(LD-A1)","adj(LD-B1)",VARY( 4),VARY( 5),STATUS)
        CALLREAD2(60,"adj(LD-A2)","adj(LD-B2)",VARY(21),VARY(24),STATUS)
        CALLREAD2(60,"adj(LD-A3)","adj(LD-B3)",VARY(22),VARY(25),STATUS)
        CALLREAD2(60,"adj(LD-A4)","adj(LD-B4)",VARY(23),VARY(26),STATUS)
        end if

        CALLREAD2(60,"adj(REFLA)","adj(REFLB)",VARY(11),VARY(12),STATUS)
        CALLREAD2(60,"adj(PSHFT)","adj(SFACT)",VARY(16),VARY(17),STATUS)
        CALLREAD2(60,"adj(PERIOD","adj(TZERO)",VARY(19),VARY(20),STATUS)
        CALL READCHAR30 (60,"OBSFILE   ",OBSFILE,STATUS)
        CALL READCHAR30 (60,"PARAMFILE ",PARFILE,STATUS)
        if ( STATUS /= 0 ) return
      end if

      if ( LDTYPE(1) == 1 ) VARY(21) = 0              ! linear LD star A
      if ( LDTYPE(2) == 1 ) VARY(24) = 0              ! linear LD star B
      if ( LDTYPE(2) == -1 ) then               ! same LD for both stars
        VARY(5) = 0
        VARY(24) = 0
        VARY(25) = 0
        VARY(26) = 0
      end if

      if ( VARY(20) == 1 .and. VARY(16) == 1 ) then
        VARY(16) = 0
        write(6,'(A39,A41)') ">> TZERO and PSHIFT cannot both be adju",
     &                      "sted so adj(PSHIFT) has been set to zero."
      end if

      if ( V(13) < 0 .and. VARY(13) /= 0 ) then
        write(6,'(A39,A41)') ">> The photometric mass ratio is equal t",
     &                      "o or less than 0 so cannot be fitted.    "
        VARY(13) = 0
        write(6,'(A39,A41)') ">> Therefore its adjustment integer has ",
     &                      "been set to zero.                        "
      end if

      if(TASK==3.or.TASK==4.or.TASK==5.or.TASK==7.or.TASK==8.or.TASK==9)
     &                  CALL READCHAR30 (60,"LC FILE   ",LCFILE, STATUS)
      if ( TASK == 3 .or. TASK==4)
     &                  CALL READCHAR30 (60,"FITFILE   ",FITFILE,STATUS)


!-----------------------------------------------------------------------
            ! Read in any additional information from the input file.
            ! TMIN: time of minimum light
            ! LRAT: spectroscopic light ratio
            ! THDL: third light measurement
            ! SINE: sine curve parameters
            ! POLY: polynomial parameters
            ! ECSW: either ecosw or eccentricity (depending on V(7))
            ! ESNW: either esinw or w (omega)    (depending on V(7))
            ! RV1: radial velocities of star A
            ! RV2: radial velocities of star B

      NLR = 0
      NMIN = 0
      do i = 1,9999
        read (60,'(A200)',iostat=ERROR) CHARHELP200
        if ( ERROR /= 0 ) exit
        CHARHELP4 = '#   '

        read (CHARHELP200,*,iostat=ERROR) CHARHELP4
!-----------------------------------------------------------------------
        if ( CHARHELP4 == "CHIF" .or. CHARHELP4 == "chif" ) then
          CHIFUDGE = 1
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "TMIN" .or. CHARHELP4 == "tmin" ) then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,(ARRAY(j),j=1,3)
          if ( ERROR /= 0 ) then
            write(6,'(A32,A28,I3)') "### ERROR reading in data for ti",
     &                             "me of minimum light, number ",NMIN+1
            STATUS = 1
          else
            if ( NDATA >= MAXDATA ) then
              write (6,'(A33,A47)') "### ERROR: maximum number of data",
     &                  "points exceeded when reading a TMIN line.     "
              STATUS = 1
            else
              NMIN = NMIN + 1
              NDATA = NDATA + 1
              DATX(NDATA) = ARRAY(1)
              DATY(NDATA) = ARRAY(2)
              DATERR(NDATA) = ARRAY(3)
              DTYPE(NDATA) = 3
            end if
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "LRAT" .or. CHARHELP4 == "lrat" ) then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,(ARRAY(j),j=1,3)
          if ( ERROR /= 0 ) then
            write(6,'(A32,A18,I3)') "### ERROR reading in data for li",
     &                                        "ght ratio, number ",NLR+1
            STATUS = 1
          else
            if ( NDATA >= MAXDATA ) then
              write (6,'(A33,A47)') "### ERROR: maximum number of data",
     &                  "points exceeded when reading an LRAT line.    "
              STATUS = 1
            else
              NLR = NLR + 1
              NDATA = NDATA + 1
              DATX(NDATA) = ARRAY(1)
              DATY(NDATA) = ARRAY(2)
              DATERR(NDATA) = ARRAY(3)
              DTYPE(NDATA) = 2
            end if
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "THDL" .or. CHARHELP4 == "thdl") then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,(ARRAY(j),j=1,2)
          if ( ERROR /= 0 ) then
            write(6,'(A32,A24,I3)') "### ERROR reading in data for th",
     &                               "ird light value, number ",NL3+1
            STATUS = 1
          else
            if ( NDATA >= MAXDATA ) then
              write (6,'(A33,A47)') "### ERROR: maximum number of data",
     &                  "points exceeded when reading a THDL line.     "
              STATUS = 1
            else
              NL3 = NL3 + 1
              NDATA = NDATA + 1
              DATX(NDATA) = 0.0d0            ! dummy value as it is unused
              DATY(NDATA) = ARRAY(1)
              DATERR(NDATA) = ARRAY(2)
              DTYPE(NDATA) = 4
            end if
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "SINE" .or. CHARHELP4 == "sine" ) then
          read (CHARHELP200,*,iostat=ERROR)  CHARHELP4, WHICHPAR,
     &                                (ARRAY(j),j=1,3),(IARRAY(j),j=1,3)
          if ( ERROR /= 0 ) then
            write(6,'(A47,I3)')
     &         "### ERROR reading in data for sine wave number ",NSINE+1
            STATUS = 1
          else
            if ( NSINE >= 20 ) then
              write (6,'(A33,A47)') "## Warning: maximum number of sin",
     &                 "e waves is 20. Any extra ones will be ignored. "
            else
              NSINE = NSINE + 1
              V(28+3*NSINE) = ARRAY(1)
              V(29+3*NSINE) = ARRAY(2)
              V(30+3*NSINE) = ARRAY(3)
              VARY(28+3*NSINE) = IARRAY(1)
              VARY(29+3*NSINE) = IARRAY(2)
              VARY(30+3*NSINE) = IARRAY(3)
              if ( WHICHPAR == "J"  ) PSINE(NSINE) = 1
              if ( WHICHPAR == "r1" ) PSINE(NSINE) = 2
              if ( WHICHPAR == "r2" ) PSINE(NSINE) = 3
              if ( WHICHPAR == "i"  ) PSINE(NSINE) = 6
              if ( WHICHPAR == "e"  ) PSINE(NSINE) = 7
              if ( WHICHPAR == "ec" ) PSINE(NSINE) = 7
              if ( WHICHPAR == "w"  ) PSINE(NSINE) = 8
              if ( WHICHPAR == "es" ) PSINE(NSINE) = 8
              if ( WHICHPAR == "L3" ) PSINE(NSINE) = 15
              if ( WHICHPAR == "sf" ) PSINE(NSINE) = 17
              if ( WHICHPAR == "T0" ) PSINE(NSINE) = 20
              if ( WHICHPAR == "KA" ) PSINE(NSINE) = 27
              if ( WHICHPAR == "KB" ) PSINE(NSINE) = 28
              if ( WHICHPAR == "VA" ) PSINE(NSINE) = 29
              if ( WHICHPAR == "VB" ) PSINE(NSINE) = 30
              if ( WHICHPAR == "L1" ) PSINE(NSINE) = -1
              if ( WHICHPAR == "L2" ) PSINE(NSINE) = -2

              if ( (WHICHPAR == "w" .or. WHICHPAR == "e") .and.
     &             V(7) < 10.0d0 ) then
                write(6,'(A33,A47)')"### ERROR: the orbital shape para",
     &                 "meters have been given as ecosw and esinw, but"
                write(6,'(A33,A47)')"### a sine wave in eccentricity o",
     &                 "r omega has been specified. To specify 'e' or "
                write(6,'(A33,A47)')"### 'w' you must change the param",
     &                 "eters of the fit from ecosw,esinw to e,omega. "
                STATUS = 1
              end if

              if ( (WHICHPAR == "ec" .or. WHICHPAR == "es") .and.
     &             V(7) >= 10.0d0 ) then
                write(6,'(A33,A47)')"### ERROR: the orbital shape para",
     &                 "meters have been given as e and omega, but a  "
                write(6,'(A33,A47)')"### sine wave in ecosw or esinw h",
     &                 "as been specified. To specify 'ec' or 'es' you"
                write(6,'(A33,A47)')"### must change the parameters of",
     &                 " the fit from ecosw and esinw to e and omega. "
                STATUS = 1
              end if

              if ( PSINE(NSINE) == 0 ) then
                write(6,'(A29,A40,A4)') "### ERROR: a valid sine wave ",
     &           "parameter has not been specified. It is ",PSINE(NSINE)
                write (6,'(A27,A42,A42)') "### and needs to be one of ",
     &                     "these parameters: J, r1, r2, i, e, w, ec, ",
     &                     "es, L3, sf, T0, L1, L2, KA, KB, VA, or VB."
                STATUS = 1
              else
                write (6,'(A28,I2,A30,A2,A1)')
     &                             ">> Read parameters for sine ",NSINE,
     &                     ', to be applied to parameter "',WHICHPAR,'"'
!                 write (62,'(A25,I1,A30,A2,A1)')
!      &                                "Read parameters for sine ",NSINE,
!      &                     ', to be applied to parameter "',WHICHPAR,'"'
              end if
            end if
          end if
!-----------------------------------------------------------------------
            ! Polynomial input can be in two forms: with and without an
            ! time interval over which to apply the polynomial. If the
            ! time interval is not given then the polynomial is applied
            ! to all datapoints. Examples:
! poly  sf  55364.66542  0.0 0.0 0.0 0.0 0.0 0.0  1 1 0 0 0 0
! poly  sf  55364.66542  0.0 0.0 0.0 0.0 0.0 0.0  1 1 0 0 0 0  55364.0 55365.0

        else if ( CHARHELP4 == "POLY" .or. CHARHELP4 == "poly" ) then

          read (CHARHELP200,*,iostat=ERROR)  CHARHELP4,WHICHPAR,POLYPIV,
     &                     POLYCONST,POLYX1,POLYX2,POLYX3,POLYX4,POLYX5,
     &                             (IARRAY(j),j=1,6),POLYSTART,POLYSTOP

          if ( ERROR /= 0 ) then
            POLYSTART = -1.0d10     ! set very low value for first data
            POLYSTOP = 1.0d10       ! set very low value for last data
            read (CHARHELP200,*,iostat=ERROR)  CHARHELP4,WHICHPAR,
     &                    POLYPIV,POLYCONST,POLYX1,POLYX2,POLYX3,POLYX4,
     &                                          POLYX5,(IARRAY(j),j=1,6)
          end if

          if ( ERROR /= 0 ) then
            write(6,'(A48,I3)')
     &        "### ERROR reading in data for polynomial number ",NPOLY+1
            write(6,*) " "
            write(6,'(A36,A44)') "Here are example polynomial input li",
     &                   "nes for fits to the total system light:     "
            write(6,'(A36,A44)') "poly  sf  55364.66542  0.0 0.0 0.0 0",
     &                   ".0 0.0 0.0  1 1 0 0 0 0                     "
            write(6,'(A36,A44)') "poly  sf  55364.66542  0.0 0.0 0.0 0",
     &                   ".0 0.0 0.0  1 1 0 0 0 0  55364.0 55365.0    "
            write(6,'(A36,A44)') "which I have used before. Remember t",
     &                   "o specify a good pivot point for your data: "
            write(6,'(A36,A44)') "a good pivot point is normally somew",
     &                   "here near the midpoint of your data in time."
            STATUS = 1
          else

            if ( NPOLY >= 50 ) then
              write (6,'(A33,A47)') "## Warning: maximum number of pol",
     &                "ynomials is 50. Any extra ones will be ignored."
            else

                  ! Add 1 to the number of polynomials.
                  ! j is the array index of a give parameter
                  ! Put the parameters and vary's into V and VARY

              NPOLY = NPOLY + 1
              j = 91 + 9*(NPOLY-1)
              V(j) = POLYPIV
              V(j+1) = POLYSTART
              V(j+2) = POLYSTOP
              V(j+3) = POLYCONST
              V(j+4) = POLYX1
              V(j+5) = POLYX2
              V(j+6) = POLYX3
              V(j+7) = POLYX4
              V(j+8) = POLYX5
              VARY(j) = 0
              VARY(j+1) = 0
              VARY(j+2) = 0
              VARY(j+3) = IARRAY(1)
              VARY(j+4) = IARRAY(2)
              VARY(j+5) = IARRAY(3)
              VARY(j+6) = IARRAY(4)
              VARY(j+7) = IARRAY(5)
              VARY(j+8) = IARRAY(6)
              if ( WHICHPAR == "J"  ) PPOLY(NPOLY) = 1
              if ( WHICHPAR == "r1" ) PPOLY(NPOLY) = 2
              if ( WHICHPAR == "r2" ) PPOLY(NPOLY) = 3
              if ( WHICHPAR == "i"  ) PPOLY(NPOLY) = 6
              if ( WHICHPAR == "e"  ) PPOLY(NPOLY) = 7
              if ( WHICHPAR == "ec" ) PPOLY(NPOLY) = 7
              if ( WHICHPAR == "w"  ) PPOLY(NPOLY) = 8
              if ( WHICHPAR == "es" ) PPOLY(NPOLY) = 8
              if ( WHICHPAR == "L3" ) PPOLY(NPOLY) = 15
              if ( WHICHPAR == "sf" ) PPOLY(NPOLY) = 17
              if ( WHICHPAR == "T0" ) PPOLY(NPOLY) = 20
              if ( WHICHPAR == "KA" ) PPOLY(NPOLY) = 27
              if ( WHICHPAR == "KB" ) PPOLY(NPOLY) = 28
              if ( WHICHPAR == "VA" ) PPOLY(NPOLY) = 29
              if ( WHICHPAR == "VB" ) PPOLY(NPOLY) = 30
              if ( WHICHPAR == "L1" ) PPOLY(NPOLY) = -1
              if ( WHICHPAR == "L2" ) PPOLY(NPOLY) = -2

              if ( (WHICHPAR == "w" .or. WHICHPAR == "e") .and.
     &             V(7) < 10.0d0 ) then
                write(6,'(A33,A47)')"### ERROR: the orbital shape para",
     &                 "meters have been given as ecosw and esinw, but"
                write(6,'(A33,A47)')"### a polynomial in eccentricity ",
     &                 "or omega has been specified. To specify 'e' or"
                write(6,'(A33,A47)')"### 'w' you must change the param",
     &                 "eters of the fit from ecosw,esinw to e,omega. "
                STATUS = 1
              end if

              if ( (WHICHPAR == "ec" .or. WHICHPAR == "es") .and.
     &             V(7) >= 10.0d0 ) then
                write(6,'(A33,A47)')"### ERROR: the orbital shape para",
     &                 "meters have been given as e and omega, but a  "
                write(6,'(A33,A47)')"### polynomial in ecosw or esinw ",
     &                 "has been specified. To go for 'ec' or 'es' you"
                write(6,'(A33,A47)')"### must change the parameters of",
     &                 " the fit from ecosw and esinw to e and omega. "
                STATUS = 1
              end if

              if ( PPOLY(NPOLY) == 0 ) then
                write(6,'(A30,A40,A4)')"### ERROR: a valid polynomial ",
     &           "parameter has not been specified. It is ",PPOLY(NPOLY)
                write (6,'(A27,A42,A42)') "### and needs to be one of ",
     &                     "these parameters: J, r1, r2, i, e, w, ec, ",
     &                     "es, L3, sf, T0, L1, L2, KA, KB, VA, or VB."
                STATUS = 1
              else
                if ( POLYSTART < -9.9d9 ) then
                  write (6,'(A34,I2,A30,A2,A1)')
     &                       ">> Read parameters for polynomial ",NPOLY,
     &                     ', to be applied to parameter "',WHICHPAR,'"'
                else
                  write (6,'(A34,I2,A30,A2,A28)')
     &                       ">> Read parameters for polynomial ",NPOLY,
     &                        ', to be applied to parameter "',WHICHPAR,
     &                                    '" over a given time interval'
                end if
              end if
            end if
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "ECSW" .or. CHARHELP4 == "ecsw") then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,(ARRAY(j),j=1,2)
          if ( ERROR /= 0 ) then
            if ( V(7) < 5.0d0 ) then
              write(6,'(A30,A33,I3)') "### ERROR reading in data for ",
     &                        "e*cos(omega) observation, number ",NECW+1
            else
              write(6,'(A30,A33,I3)') "### ERROR reading in data for ",
     &                        "eccentricity observation, number ",NECW+1
            end if
            STATUS = 1
          else
            NECW = NECW + 1
            NDATA = NDATA + 1
            DATX(NDATA) = 0.0d0            ! dummy value as it is unused
            DATY(NDATA) = ARRAY(1)
            DATERR(NDATA) = ARRAY(2)
            DTYPE(NDATA) = 5
          end if
          if ( ARRAY(1) >= 10.0d0 .and. V(7) < 10.0d0 ) then
            write(6,'(A35,A45)') "## Warning: the ecosw constraint val",
     &                 "ue is greater than or equal to 10, indicating "
            write(6,'(A35,A45)') "## that it represents e, but the mai",
     &                 "n input value for this quantity is less than  "
            write(6,'(A35,A45)') "## ten, indicating that it represent",
     &                 "s e*cos(omega). This inconsistency could cause"
            write(6,'(A35,A45)') "## a poor fit so the main input quan",
     &                 "tity has been changed from e*cos(omega) to e. "
            V(7) = V(7) + 10.0d0
          endif
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "ESNW" .or. CHARHELP4 == "esnw") then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,(ARRAY(j),j=1,2)
          if ( ERROR /= 0 ) then
            if ( V(7) < 5.0d0 ) then
              write(6,'(A30,A33,I3)') "### ERROR reading in data for ",
     &                        "e*sin(omega) observation, number ",NECW+1
            else
              write(6,'(A30,A41,I3)') "### ERROR reading in data for ",
     &                "periastron longitude observation, number ",NECW+1
            end if
            STATUS = 1
          else
            NESW = NESW + 1
            NDATA = NDATA + 1
            DATX(NDATA) = 0.0d0            ! dummy value as it is unused
            DATY(NDATA) = ARRAY(1)
            DATERR(NDATA) = ARRAY(2)
            DTYPE(NDATA) = 6
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "NUMI" .or. CHARHELP4 == "numi") then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,NUMINT,NINTERVAL
          if ( ERROR /= 0 ) then
            write(6,'(A35,A45)') "### ERROR reading in instructions f",
     &                  "or numerical integration.                    "
            STATUS = 1
          else
            if ( NUMINT < 2 ) then
              write(6,'(A35,A45)') "### ERROR reading instructions for",
     &                  " numerical integration: NUMINT is less than 2"
              STATUS = 1
            end if
            if ( NUMINT > 10 .and. NUMINT < 100 ) then
              write(6,'(A35,A45)') "### Warning: NUMINT is more than 1",
     &                  "0. This will slow JKTEBOP down significantly."
            endif
            if ( NUMINT >= 100 ) then
              write(6,'(A35,A45)') "### ERROR with instructions for nu",
     &                  "merical integration: NUMINT is more than 100."
              STATUS = 1
            end if
         end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "RV1" .or. CHARHELP4 == "rv1" ) then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,RV1OBSFILE,
     &                          RV1OUTFILE,V(27),V(29),VARY(27),VARY(29)
          if ( ERROR /= 0 ) then
            write(6,'(A35,A45)') "### ERROR reading in instructions f",
     &                  "or radial velocities of star A.              "
            STATUS = 1
          end if
          if ( VARY(27) < 0 .or. VARY(27) > 3 ) then
            write(6,'(A35,A45)') "### ERROR reading in vary(K) for sta",
     &                  "r A. The value should be 0, 1, 2 or 3.       "
            STATUS = 1
          end if
          if ( VARY(29) < 0 .or. VARY(29) > 3 ) then
            write(6,'(A35,A45)') "### ERROR reading in vary(Vsys) for ",
     &                  "star A. The value should be 0, 1, 2 or 3.    "
            STATUS = 1
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "RV2" .or. CHARHELP4 == "rv2" ) then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,RV2OBSFILE,
     &                          RV2OUTFILE,V(28),V(30),VARY(28),VARY(30)
          if ( ERROR /= 0 ) then
            write(6,'(A35,A45)') "### ERROR reading in instructions f",
     &                  "or radial velocities of star B.              "
            STATUS = 1
          end if
          if ( VARY(28) < 0 .or. VARY(28) > 3 ) then
            write(6,'(A35,A45)') "### ERROR reading in vary(K) for sta",
     &                  "r B. The value should be 0, 1, 2 or 3.       "
            STATUS = 1
          end if
          if ( VARY(30) < -1 .or. VARY(30) > 3 ) then
            write(6,'(A35,A45)') "### ERROR reading in vary(Vsys) for ",
     &                  "star B. The value should be -1, 0, 1, 2 or 3."
            STATUS = 1
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "NOFI" .or. CHARHELP4 == "nofi" ) then
          NOFIT = 1
!-----------------------------------------------------------------------
        end if
      end do

      close (unit=60,status="keep")
      if ( STATUS /= 0 ) return

!-----------------------------------------------------------------------
            ! Open the output files and start writing results to them.

      if ( TASK >= 3 ) then

        CALL OPENFILE(62,"new","parameter ",PARFILE,STATUS)
        if (TASK==5) CALL OPENFILE(63,"new","refit     ",LCFILE, STATUS)
        if (TASK==7) CALL OPENFILE(63,"new","bootstrap ",LCFILE, STATUS)
        if (TASK==8) CALL OPENFILE(63,"new","simulation",LCFILE, STATUS)
        if (TASK==9) CALL OPENFILE(63,"new","simulation",LCFILE, STATUS)

        if (TASK == 3 .or. TASK == 4) then
          CALL OPENFILE (63,"new","LC output ",LCFILE, STATUS)
          if ( V(27) >= -1.0d9 ) CALL OPENFILE (65,"new","RV1 output",
     &                                                RV1OUTFILE,STATUS)
          if ( V(28) >= -1.0d9 ) CALL OPENFILE (66,"new","RV2 output",
     &                                                RV2OUTFILE,STATUS)
          CALL OPENFILE (64,"new","model fit ",FITFILE,STATUS)
        end if
        if ( STATUS /= 0 ) return

        write (62,'(A38,A42)') "======================================",
     &                     "=========================================="
        write (62,'(A38,A42)') "JKTEBOP  output results               ",
     &                     "    John Southworth (astro.js~keele.ac.uk)"
        write (62,'(A38,A42)') "======================================",
     &                     "=========================================="

        CALL date_and_time(DTDATE,DTTIME)
        write (62,'(A32,I2,6X,A24,A2,":",A2,1X,A2,"/",A2,"/",A4,A8)')
     &                       "Version number of JKTEBOP code: ",VERSION,
     &               "Time and date at start: ",DTTIME(1:2),DTTIME(3:4),
     &                               DTDATE(7:8),DTDATE(5:6),DTDATE(1:4)
        write (62,*)   " "

        write(62,'(A32,I1)')  "Task to perform:                ", TASK
        write(62,'(A32,A50)') "Input parameter file:           ", INFILE
        write(62,'(A32,A50)') "Output parameter file:          ",PARFILE
        if ( TASK == 3 .or. TASK == 4 ) write (62,'(A32,A50)')
     &                        "Output model fit:               ",FITFILE
        write(62,'(A32,A50)') "Input light curve data file:    ",OBSFILE
        if ( V(27) >= -1.0d9 ) write (62,'(A32,A50)')
     &                     "Input RV data file for star A:  ",RV1OBSFILE
        if ( V(28) >= -1.0d9 ) write (62,'(A32,A50)')
     &                     "Input RV data file for star B:  ",RV2OBSFILE

        if ( TASK == 3 .or. TASK == 4 ) then
          write(62,'(A32,A50)')"Output light curve data file:   ",LCFILE
          if ( V(27) >= -1.0d9 ) write(62,'(A32,A50)')
     &                     "Output RV datafile for star A:  ",RV1OUTFILE
          if ( V(28) >= -1.0d9 ) write(62,'(A32,A50)')
     &                     "Output RV datafile for star B:  ",RV2OUTFILE
        else if ( TASK == 5 ) then
          write(62,'(A32,A50)')"Perturbed refit output file:    ",LCFILE
        else if ( TASK == 7 ) then
          write(62,'(A32,A50)')"Bootstrapping output file:      ",LCFILE
        else if ( TASK == 8 ) then
          write(62,'(A32,A50)')"Monte Carlo simulation file:    ",LCFILE
        else if ( TASK == 9 ) then
          write(62,'(A32,A50)')"Residual-permutation sim file:  ",LCFILE
        end if

        write (62,*)   " "
      end if

!-----------------------------------------------------------------------
            ! Read in the photometric and radial velocity data
            ! Enter flag into the DATERR array for the case of "chif"

      if ( TASK >= 3 .and. TASK <= 9 ) then
        NDATA_NOTLC = NDATA
        CALL READDATA (1,OBSFILE,DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,
     &                 DATAFORM,STATUS)
        if ( STATUS /= 0 ) return
        NDATA_LC = NDATA - NDATA_NOTLC

        if ( CHIFUDGE == 1 ) V(14) = -99.5d0

            ! Have to set NSIM for TASK 9 here, as it is used to specify
            ! the sizes of some arrays in the subroutine TASK5789.  This
            ! is based on the assumption that the number of LC datapoint
            ! is greater than the number of RV1 or RV2 datapoints.

        if ( TASK == 9 ) NSIM = NDATA_LC - 1

        if ( V(27) >= -1.0d9 ) CALL READDATA (7,RV1OBSFILE,DATX,DATY,
     &                       DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,STATUS)

        if ( V(28) >= -1.0d9 ) CALL READDATA (8,RV2OBSFILE,DATX,DATY,
     &                       DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,STATUS)
      end if

!-----------------------------------------------------------------------
            ! Output various bits of useful information to the output
            ! parameter file and also to standard output.

      if ( V(18) > -1.01d0 .and. V(18) < -0.99d0 ) then
        write (6,'(A38,A42)')  ">> INTRING = -1  so the Mandel & Agol ",
     &                     "OCCULTSMALL subroutine will be used.      "
        write (62,'(A38,A42)') ">> INTRING = -1  so the Mandel & Agol ",
     &                     "OCCULTSMALL subroutine will be used.      "
      end if

      if ( NMIN > 0 ) write(6,'(A7,I3,A43)') ">> Read",NMIN,
     &                     " times of minimum light from the input file"
      if ( NMIN > 0 ) write (62,'(A4,I3,A43)') "Read",NMIN,
     &                     " times of minimum light from the input file"

      if ( NLR > 0 )  write(6,'(A7,I3,A47)') ">> Read",NLR,
     &                 " spectroscopic light ratios from the input file"
      if ( NLR > 0 )  write (62,'(A4,I3,A47)') "Read",NLR,
     &                 " spectroscopic light ratios from the input file"

      if ( NL3 > 0 )  write(6,'(A7,I3,A45)') ">> Read",NL3,
     &                   " third light measurements from the input file"
      if ( NL3 > 0 )  write (62,'(A4,I3,A45)') "Read",NL3,
     &                   " third light measurements from the input file"

      if ( NECW > 0 ) then
        if ( V(7) < 5.0d0 ) write(6,'(A7,I3,A46)') ">> Read",NECW,
     &                  " e*cos(omega) measurements from the input file"
        if ( V(7) < 5.0d0 ) write (62,'(A4,I3,A46)') "Read",NECW,
     &                  " e*cos(omega) measurements from the input file"
        if ( V(7) > 5.0d0 ) write(6,'(A7,I3,A46)') ">> Read",NECW,
     &                  " eccentricity measurements from the input file"
        if ( V(7) > 5.0d0 ) write (62,'(A4,I3,A46)') "Read",NECW,
     &                  " eccentricity measurements from the input file"
      end if

      if ( NESW > 0 ) then
        if ( V(7) < 5.0d0 ) write(6,'(A7,I3,A46)') ">> Read",NESW,
     &                  " e*sin(omega) measurements from the input file"
        if ( V(7) < 5.0d0 ) write (62,'(A4,I3,A46)') "Read",NESW,
     &                  " e*sin(omega) measurements from the input file"
        if ( V(7) > 5.0d0 ) write(6,'(A7,I3,A46)') ">> Read",NESW,
     &                  " periast-long measurements from the input file"
        if ( V(7) > 5.0d0 ) write (62,'(A4,I3,A46)') "Read",NESW,
     &                  " periast-long measurements from the input file"
      end if

      if ( NSINE > 0 ) write(6,'(A7,I2,A39)') ">> Read ",NSINE,
     &                         " sine wave datasets from the input file"
      if ( NSINE > 0 ) write (62,'(A5,I2,A39)') "Read ",NSINE,
     &                         " sine wave datasets from the input file"

      if ( NPOLY > 0 ) then
        write(6,'(A8,I2,A40)') ">> Read ",NPOLY,
     &                        " polynomial datasets from the input file"
        write (62,'(A5,I2,A40)') "Read ",NPOLY,
     &                        " polynomial datasets from the input file"
        do i = 1,NPOLY
          write (62,'(A11,I2,A36,$)') "Polynomial ",i,
     &                            " will be applied to this parameter: "
          if ( PPOLY(i) ==  1 ) write (62,'(A2)') "J "
          if ( PPOLY(i) ==  2 ) write (62,'(A2)') "r1"
          if ( PPOLY(i) ==  3 ) write (62,'(A2)') "r2"
          if ( PPOLY(i) ==  6 ) write (62,'(A2)') "i "
          if ( PPOLY(i) ==  7 ) write (62,'(A2)') "e "
          if ( PPOLY(i) ==  7 ) write (62,'(A2)') "ec"
          if ( PPOLY(i) ==  8 ) write (62,'(A2)') "w "
          if ( PPOLY(i) ==  8 ) write (62,'(A2)') "es"
          if ( PPOLY(i) == 15 ) write (62,'(A2)') "L3"
          if ( PPOLY(i) == 17 ) write (62,'(A2)') "sf"
          if ( PPOLY(i) == 20 ) write (62,'(A2)') "T0"
          if ( PPOLY(i) == 27 ) write (62,'(A2)') "KA"
          if ( PPOLY(i) == 28 ) write (62,'(A2)') "KB"
          if ( PPOLY(i) == 29 ) write (62,'(A2)') "VA"
          if ( PPOLY(i) == 30 ) write (62,'(A2)') "VB"
          if ( PPOLY(i) == -1 ) write (62,'(A2)') "L1"
          if ( PPOLY(i) == -2 ) write (62,'(A2)') "L2"
        end do
      end if

      if ( NUMINT > 1 ) write(6,'(A23,A43)') ">> Read instructions fo",
     &                    "r numerical integration from the input file"
      if ( NUMINT > 1 ) write (62,'(A30,I4,A18,F7.2,A8)')
     &                        "Numerical integration invoked:",NUMINT,
     &                        " samples covering ",NINTERVAL," seconds"

      write (62,*) " "


!-----------------------------------------------------------------------
            ! Check the polynomial time intervals to see if any contain
            ! no datapoints

      if ( NPOLY >= 1 ) then
        do i = 1,NPOLY
          k = 91 + 9*(i-1)
          IHELP = 0
          do j = 1,NDATA
            if ( DATX(j)>V(k+1) .and. DATX(j)<V(k+2) ) IHELP = IHELP + 1
          end do
          if ( IHELP < 1 ) then
            write(6,'(A1)') " "
            write(6,'(A30,A47,I1,A1)') "### ERROR: there are no datapo",
     &           "ints in the time interval given for polynomial ",i,"."
            write(6,'(A36,A44)') "### This will cause the fitting rout",
     &                    "ine to fail. Check the input data and file. "
            write(6,'(A1)') " "
            STATUS = 1
          end if
        end do
      end if

!-----------------------------------------------------------------------

      if ( DATERR(1) < 0.0d0 ) then
        if ( NMIN>0 .or. NLR>0 .or. NL3>0 .or. NECW>0 .or. NESW>0 ) then
          write(6,*) " "
          write(6,'(A37,A43)') "## Warning: no observational errors h",
     &                    "ave been found but additional information  "
          write(6,'(A37,A43)') "## (times of minimum light, spectrosc",
     &                    "opic light ratios, e or omega values) have "
          write(6,'(A37,A43)') "## been entered. Proper uncertainties",
     &                    " are therefore needed for all input data.  "
          write(6,*) " "
        end if
      end if

      if ( NL3 > 0 .and. VARY(15) == 0 ) then
        write(6,*) " "
        write(6,'(A39,A41)') "## Warning: a constraint on third light",
     &                     " value has been given but third light is "
        write(6,'(A39,A41)') "## not a fitted parameter, so the const",
     &                     "raint will not be taken into account.    "
        write(6,*) " "
      endif

      if ( (NECW>0.or.NESW>0 ) .and. VARY(7)==0 .and. VARY(8)==0 ) then
        write(6,*) " "
        write(6,'(A39,A41)') "## Warning: a constraint on orbital ecc",
     &                     "entricity or periastron longitude has    "
        write(6,'(A39,A41)') "## been given, but both parameters are ",
     &                     "fixed so the constraint will be ignored. "
        write(6,*) " "
      endif

      if ( V(13) < 0.0d0 ) write (62,'(A19,A61)') "Mass ratio is below",
     &   " zero: the stellar shapes will be forced to be spherical.    "
      if ( V(7) < 5.0d0 )  write (62,'(A19,A61)') "Eccentricity is bel",
     &   "ow 5:  [e,omega] are taken to be [e*cos(omega),e*sin(omega)]."
      if ( V(7) > 5.0d0 )  write (62,'(A19,A61)') "Eccentricity is mor",
     &   "e than 5:  e and omega will be taken to be (e+10) and omega. "
      if ( V(2) < 0.0d0 )  write (62,'(A19,A61)') "r1+r2 is less than ",
     &   "zero so (r1+r2) and k will be interpreted to mean -r1 and r2"

      if ( V(7) < 5.0d0 ) then
        if (V(7)<-1.0d0.or.V(7)>1.0d0.or.V(8)<-1.0d0.or.V(8)>1.0d0) then
          write(6,'(A36,A44)')  "### ERROR: e*sin(omega) or e*cos(ome",
     &                    "ga) are unphysical: must be between -1 and 1"
          STATUS = 1
          return
        end if
      end if

                    ! Sort out radii as that there are two entry options

      if ( V(2) >= 0.0d0 ) then
        R1 = V(2) / (1.0d0 + V(3))
        R2 = V(2) / (1.0d0 + (1.0d0/V(3)))
      else
        R1 = abs(V(2))
        R2 = V(3)
      end if

                    ! Find reflection coefficients if they are not fixed

      if ( TASK >= 3 .and. TASK <= 9 ) then

        MAG = GETMODEL(V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,V(20),1,LP,
     &                                              LS,NUMINT,NINTERVAL)
        if ( VARY(11) == -1 ) V(11) = 0.4d0 * (LS/(1.0d0-V(15))) * R1**2
        if ( VARY(12) == -1 ) V(12) = 0.4d0 * (LP/(1.0d0-V(15))) * R2**2

        do i = 11,12
100       FORMAT (A4,I2,A37,I1,A36)
          if ( VARY(i) == -1 ) write (62,100)
     &                 "Adj(",i,") = -1 so reflection effect for star ",
     &                      i-10," is calculated from system geometry. "
          if ( VARY(i) == 0 )  write (62,100)
     &                 "Adj(",i,") = -1 so reflection effect for star ",
     &                      i-10," is fixed at the input file value.   "
          if ( VARY(i) == 1 )  write (62,100)
     &                 "Adj(",i,") = +1 so reflection effect for star ",
     &                      i-10," is freely adjusted to the best fit. "
          if ( VARY(i) == 2 )  write (62,100)
     &                 "Adj(",i,") = +1 so reflection effect for star ",
     &                      i-10," is set to input value but perturbed."
        end do

        if ( VARY(i) == -1 .and. NLR > 0 ) then
          write(6,'(A37,A43)') "## Warning: if there is a spectroscop",
     &                    "ic light ratio it's best to directly adjust"
          write(6,'(A37,A43)') "## reflection coefficients, rather th",
     &                    "an calculate them from the system geometry."
        end if
        write (62,*) " "

                  ! Now to avoid problems with partial derivatives later

        if ( V(7)>5.0d0 .and. VARY(7)+VARY(8)/=0)  V(7)=max(V(7),10.001)
        if ( V(7) < 6.0d0 .and. V(7) >= 1.0d0 )    V(7) = 0.0d0
        if ( V(11) < 0.001d0 .and. VARY(11) == 1 ) V(11) = 0.001d0
        if ( V(12) < 0.001d0 .and. VARY(12) == 1 ) V(12) = 0.001d0

            ! Find a good starting value for the light scale factor

        if ( VARY(17) /= 0 ) then
          do i = 1,NDATA
            ARRAY(i) = DATY(i)
          end do
          V(17) = sellect(ARRAY,NDATA,int(0.2d0*NDATA))
        end if
      end if

            ! If LD of star B is forced to be same as star A, then set
            ! the (unused) LD of star B to be fixed to avoid problems.

      if ( LDTYPE(2) == -1 ) then
        V(5) = V(4)
        VARY(5) = 0
      end if

            ! If no parameters should be fitted (because 'nofit' has
            ! been specified) then put all the VARY integers to zero.

      if ( NOFIT == 1) then
        do i = 1,540
          VARY(i) = 0
        end do
        write (6,*) " "
        write (6,'(A38,A42)') "## Warning: 'nofit' has been specified",
     &                    " so no parameters will be fitted.         "
        write (6,*) " "
      end if

      END SUBROUTINE INPUT
!=======================================================================
!=======================================================================
      SUBROUTINE READFF (UNIT,NAME1,LOWER1,UPPER1,VALUE1,
     &                        NAME2,LOWER2,UPPER2,VALUE2,STATUS)
            ! This reads in two double precision numbers and checks if
            ! they are within given bounds. STATUS is set to 1 if error
            ! occurs and keeps its previous value if no error occurs.
      implicit none
      integer UNIT                        ! IN: unit number to read from
      character NAME1*10,NAME2*10         ! IN: names of the  parameters
      real*8 LOWER1,UPPER1,LOWER2,UPPER2  ! IN: allowed ranges of values
      real*8 VALUE1,VALUE2                ! OUT:values of the parameters
      integer STATUS                      ! OUT:set to 1 if error occurs
      integer ERROR                       ! LOCAL:  error flag for input

      ERROR = 0
      read (UNIT,*,iostat=ERROR) VALUE1,VALUE2
      if ( ERROR /= 0 ) then
        write(6,'(A33,A10,A5,A10)')
     &           "### ERROR reading the parameters ",NAME1," and ",NAME2
        STATUS = 1
      end if

      if ( VALUE1 < LOWER1 .or. VALUE1 > UPPER1 ) then
        write(6,'(A24,A10,A4,F12.5)')
     &                    "### ERROR: the value of ",NAME1," is ",VALUE1
        write(6,'(A21,F12.5,A4,F12.5)')
     &                      "The allowed range is ",LOWER1," to ",UPPER1
        STATUS = 1
      end if

      if ( VALUE2 < LOWER2 .or. VALUE2 > UPPER2 ) then
        write(6,'(A24,A10,A4,F12.5)')
     &                    "### ERROR: the value of ",NAME2," is ",VALUE2
        write(6,'(A21,F12.5,A4,F12.5)')
     &                      "The allowed range is ",LOWER2," to ",UPPER2
        STATUS = 1
      end if

      END SUBROUTINE READFF
!-----------------------------------------------------------------------
      SUBROUTINE READF (UNIT,NAME,LOWER,UPPER,VALUE,STATUS)
            ! This reads one double precision number and checks if it is
            ! within given bounds. STATUS is set to 1 if an error occurs
      implicit none
      integer UNIT                  ! IN: number of unit to read from
      character NAME*10             ! IN: name of parameter
      real*8 LOWER,UPPER            ! IN: allowed range of values
      real*8 VALUE                  ! OUT: value of the parameter
      integer STATUS                ! IN/OUT: set to 1 if error occurs
      integer ERROR                 ! LOCAL: error flag for input

      ERROR = 0
      read (UNIT,*,iostat=ERROR) VALUE
      if ( ERROR /= 0 ) then
        write(6,'(A32,A10)') "### ERROR reading the parameter ",NAME
        STATUS = 1
      end if
      if ( VALUE < LOWER .or. VALUE > UPPER ) then
        write(6,'(A24,A10,A4,F12.5)')
     &                      "### ERROR: the value of ",NAME," is ",VALUE
        write(6,'(A21,F12.5,A4,F12.5)')
     &                        "The allowed range is ",LOWER," to ",UPPER
        STATUS = 1
      end if

      END SUBROUTINE READF
!-----------------------------------------------------------------------
      SUBROUTINE READCHAR30 (UNIT,NAME,VALUE,STATUS)
            ! This reads in a 30-char value. STATUS is set to 1 if an
            ! error occurs and keeps its previous value if no error.
      implicit none
      integer UNIT                  ! IN: number of unit to read from
      character NAME*10             ! IN: name of character to be read
      character VALUE*50            ! OUT: value of the character
      integer STATUS                ! IN/OUT: set to 1 if error occurs
      integer ERROR                 ! LOCAL: error flag for input

      ERROR = 0
      read (UNIT,*,iostat=ERROR) VALUE
      if ( ERROR /= 0 ) then
        write(6,'(A29,A50)') "### ERROR reading the string ",NAME
        STATUS = 1
      end if
      if ( VALUE == "#" ) then
        write(6,'(A39,A41)') "### ERROR: you cannot use '#' as an out",
     &                      "put file name.                           "
        STATUS = 1
      end if

      END SUBROUTINE READCHAR30
!-----------------------------------------------------------------------
      SUBROUTINE READ2 (UNIT,NAME1,NAME2,VALUE1,VALUE2,STATUS)
            ! This reads two integers on one line and tests if both are
            ! between -1 and 3. If not then STATUS is set to 1.
            ! If there is no error then STATUS keeps its previous value.
      implicit none
      integer UNIT                  ! IN: number of unit to read from
      character NAME1*10,NAME2*10   ! IN: names of characters to read
      integer VALUE1,VALUE2         ! OUT: values of the integers
      integer STATUS                ! IN/OUT: set to 1 if error occurs
      integer ERROR                 ! LOCAL: error flag for input

      ERROR = 0
      read (UNIT,*,iostat=ERROR) VALUE1,VALUE2
      if ( ERROR /= 0 ) then
        write(6,'(A31,A10,A5,A10)')
     &            "### ERROR reading the integers ",NAME1," and ",NAME2
        STATUS = 1
      end if
      if ( VALUE1 < -1 .or. VALUE1 > 3 ) then
        write(6,'(A30,A10,A4,I6,A28)') "### ERROR: adjustment integer ",
     &                NAME1," is ",VALUE1," and should be -1, 0, 1 or 2"
        STATUS = 1
      end if
      if ( VALUE2 < -1 .or. VALUE2 > 3 ) then
        write(6,'(A30,A10,A4,I6,A28)') "### ERROR: adjustment integer ",
     &                NAME2," is ",VALUE2," and should be -1, 0, 1 or 2"
        STATUS = 1
      end if

      END SUBROUTINE READ2
!=======================================================================
      SUBROUTINE READDATA (WHICHDTYPE,OBSFILE,DATX,DATY,DATERR,DTYPE,
     &                                    NDATA,MAXDATA,DATAFORM,STATUS)

            ! This subroutine opens an input data file, reads all the
            ! data, and closes the file again. Different types of data
            ! are specified using WHICHDTYPE (in) and DTYPE (out).

      implicit none
      integer MAXDATA               ! IN: maximum number of datapoints
      character DATAFORM*1          ! IN: character length of MAXDATA
      character*50 OBSFILE          ! IN: Name of input light curve file
      integer WHICHDTYPE            ! IN: which type of data it is
      real*8 DATX(MAXDATA)          ! OUT: Data independent variable
      real*8 DATY(MAXDATA)          ! OUT: Data dependent variable
      real*8 DATERR(MAXDATA)        ! OUT: Data errorbars
      integer DTYPE(MAXDATA)        ! OUT: Type of data (here all "1")
      integer NDATA                 ! OUT: Number of datapoints in total
      integer STATUS                ! IN/OUT: set to 1 if there is error
      integer NDATAIN               ! IN:    Number of datapoints before
      integer i,j,ERROR,ERRFLAG     ! LOCAL: Loop counter + error flags
      integer N_NAN,N_LOW,N_HIGH    ! LOCAL: number of NaN, low or high
      integer N_NANERR              ! LOCAL: number of NaNs in errorbars
      character*200 CHARHELP        ! LOCAL: Helper character string
      real*8 HELP1,HELP2,HELP3      ! LOCAL: Helper variables

      NDATAIN = NDATA

      ERROR = 0

      if ( WHICHDTYPE == 1 ) then
        CALL OPENFILE (61,"old","LC data   ",OBSFILE,ERROR)
      else if ( WHICHDTYPE == 7 ) then
        CALL OPENFILE (61,"old","RV A data ",OBSFILE,ERROR)
      else if ( WHICHDTYPE == 8 ) then
        CALL OPENFILE (61,"old","RV B data ",OBSFILE,ERROR)
      else
        write (6,'(A36,A40,I4)') "### ERROR: wrong type of data specif",
     &             "ied in subroutine READDATA:  WHICHDTYPE=",WHICHDTYPE
        STOP
      end if

      if ( ERROR /= 0 ) then
        STATUS = 1
        return
      end if

      read (61,'(A200)',iostat=ERROR) CHARHELP
      rewind (61)

      read (CHARHELP,*,iostat=ERROR) HELP1,HELP2,HELP3
      if ( ERROR /= 0 ) then
        read (CHARHELP,*,iostat=ERROR) HELP1,HELP2
        if ( ERROR /= 0 ) then
          write(6,'(A48,A50)')
     &        "### ERROR: cannot understand first line of file ",OBSFILE
          STATUS = 1
          return
        end if
        ERRFLAG = 0                  ! There are no observational errors
      else
        ERRFLAG = 1                     ! There are observational errors
      end if

      do i = NDATAIN+1,MAXDATA
        if ( ERRFLAG == 0 ) then
          read (61,*,iostat=ERROR) DATX(i),DATY(i)
          DATERR(i) = -1.0d0
        else if ( ERRFLAG == 1 ) then
          read (61,*,iostat=ERROR) DATX(i),DATY(i),DATERR(i)
          if ( ERROR == 0 .and. DATERR(i) <= 0.0d0 ) then
            j = i-NDATAIN
            write(6,'(A44,I'//DATAFORM//',A9,A50)')
     &                   "### ERROR: found errorbar <= 0 for datapoint",
     &                                             j," in file ",OBSFILE
            write(6,'(A21,F16.7,1X,F15.7,1X,F15.7)')
     &                "### Data values are: ", DATX(j),DATY(j),DATERR(j)
            STOP
          end if
        end if
        if ( ERROR /= 0 ) exit
        DTYPE(i) = WHICHDTYPE
        NDATA = NDATA + 1
      end do

      if ( NDATA < 5 ) then
        write(6,'(A30)') "### ERROR: too few data to fit"
        STATUS = 1
      end if

      if ( ERRFLAG == 1 ) then
        write (62,'(A5,I'//DATAFORM//',A39,A50)') "Read ",NDATA-NDATAIN,
     &                 " datapoints (with errorbars) from file ",OBSFILE
        write(6,'(A8,I'//DATAFORM//',A36,A50)')">> Read ",NDATA-NDATAIN,
     &                    " datapoints (with errors) from file ",OBSFILE
      else if  ( ERRFLAG == 0 ) then
        write (62,'(A5,I'//DATAFORM//',A37,A50)') "Read ",NDATA-NDATAIN,
     &                   " datapoints (no errorbars) from file ",OBSFILE
        write(6,'(A8,I'//DATAFORM//',A37,A50)')">> Read ",NDATA-NDATAIN,
     &                   " datapoints (no errorbars) from file ",OBSFILE
      end if

      close (61)

            ! Time to check if the input data make sense. Search for
            ! three things: values less than -100 or more than +100
            ! (for the magnitudes, not the RVs), and values or errorbars
            ! that are NaN. The NaN values can be caught by requiring
            ! datapoints to equal themselves - NaN by definition equal
            ! no value (including themselves) so fail this test.
            ! Added on 2022/12/05 for version 43 of the code.

      N_NAN = 0
      N_NANERR = 0
      N_LOW = 0
      N_HIGH = 0
      do i = NDATAIN+1,MAXDATA
        if ( DATY(i) <= -100.0 ) N_LOW = N_LOW + 1
        if ( DATY(i) >= 100.0 ) N_HIGH = N_HIGH + 1
        if ( DATY(i) /= DATY(i) ) N_NAN = N_NAN + 1
        if ( DATERR(i) /= DATERR(i) ) N_NANERR = N_NANERR + 1
      end do

      if ( WHICHDTYPE == 1 ) then
        if ( N_LOW >= 1 ) then
          write (6,'(A1)') " "
          write (6,'(A12,I'//DATAFORM//',A30,A30)') "## Warning: ",
     &                          N_LOW, " datapoints are -100 or less. ",
     &                                 " This very likely means that  "
          write (6,'(A37,A43)') ">> the data are not in magnitude unit",
     &                    "s. Please check input and output CAREFULLY."
          write (6,'(A1)') " "
        end if
        if ( N_HIGH >= 1 ) then
          write (6,'(A1)') " "
          write (6,'(A12,I'//DATAFORM//',A30,A30)') "## Warning: ",
     &                         N_HIGH, " datapoints are 100 or more.  ",
     &                                 " This very likely means that  "
          write (6,'(A37,A43)') ">> the data are not in magnitude unit",
     &                    "s. Please check input and output CAREFULLY."
          write (6,'(A1)') " "
        end if
      endif

      if ( N_NAN >= 1 ) then
        write (6,'(A17,I'//DATAFORM//',A30,A50)') "### ERROR: found ",
     &                  N_NAN, " datapoints to be NaN in file ", OBSFILE
        STATUS = 1
      end if
      if ( N_NANERR >= 1 ) then
        write (6,'(A17,I'//DATAFORM//',A30,A50)') "### ERROR: found ",
     &                N_NANERR, " errorbars to be NaN in file ", OBSFILE
        STATUS = 1
      end if

      END SUBROUTINE READDATA
!=======================================================================
      SUBROUTINE OPENFILE (UNIT,STATE,FILE,FILENAME,STATUS)
            ! Opens a file. STATUS = 1 if the action was not
            ! successful and left unchanged if the action was successful
      implicit none
      integer UNIT                  ! IN: unit number to open file to
      character*50 FILENAME         ! IN: name of datafile to open
      character*10 FILE             ! IN: identifier of this datafile
      character*3 STATE             ! IN: "old" or "new"
      integer STATUS                ! OUT: indicates success of opening
      integer ERROR                 ! LOCAL: error flag for opening file

      ERROR = 0

      OPEN (unit=UNIT,file=FILENAME,status=STATE,iostat=ERROR)

      if ( ERROR /= 0 ) then
        write(6,'(A18,A3,1X,A10,A8,A50)')
     &              "### ERROR opening ", STATE,FILE," file:  ",FILENAME
        STATUS = 1
      end if

      if (STATE=="new".and.ERROR==0)  write(6,'(A10,A3,1X,A10,A8,A50)')
     &                       ">> Opened ",STATE,FILE," file:  ",FILENAME

      END SUBROUTINE OPENFILE
!=======================================================================
!=======================================================================
      SUBROUTINE TASK1 ()
                                    ! This task outputs  limb  darkening
            ! coefficients for given Teff and logg and [M/H] and Vmicro.
            ! It simply interfaces with the  JKTLD  code, which performs
            ! bilinear interpolation in Teff,logg for given [M/H],Vmicro
            ! Usage:   jktld  <Teff>  <logg>  <M/H>  <Vmicro>  <outfile>

      implicit none
      real*8 TEFF,LOGG              ! IN: Parameters  to  interpolate to
      real*8 MOH,VMICRO             ! IN: Other parameters forthe tables
      character*20 CTEFF,CLOGG      ! LOCAL: character version of values
      character*20 CMOH,CMICRO      ! LOCAL: character version of values
      character*50 OUTFILE          ! LOCAL: name of output file to make

      write(6,'(A40,$)') "Enter the effective temperature (K)  >> "
      read (*,*) TEFF
      write(6,'(A40,$)') "Enter the surface gravity (log cm/s) >> "
      read (*,*) LOGG
      write(6,'(A40,$)') "Enter the metal abundance  ([M/H])   >> "
      read (*,*) MOH
      write(6,'(A40,$)') "Enter the microturbulence velocity   >> "
      read (*,*) VMICRO
      write(6,'(A40,$)') "Enter the output file name to create >> "
      read (*,*) OUTFILE

      if ( TEFF < 3500.0d0 ) then
        write(6,'(A39,A41)') "## Warning: a Teff below 3500 K is out ",
     &                      "of range of most of the LD coeff tables. "
      else if ( TEFF < 2000.0d0 ) then
        write(6,'(A39,A41)') "## Warning: a Teff below 2000 K is out ",
     &                      "of range all the LD coefficient tables.  "
      else if ( TEFF > 50000.0d0 ) then
        write(6,'(A39,A41)') "## Warning: a Teff above 50000 K is out",
     &                      " of range of all LD coefficient tables.  "
      end if
      if ( LOGG > 5.0d0 .or. LOGG < 0.0d0 )
     &  write(6,'(A39,A41)') "## Warning: log(g)s outside the range 0",
     &                      ".0 to 5.0 are not covered in the tables. "
      if ( MOH /= 0.0d0 .and. VMICRO /= 2.0d0 )
     &  write(6,'(A39,A41)') "## Warning: for [M/H] /= 0 and Vmicro /",
     &                      "= 2 you will probably get no results.    "

      write (CTEFF,'(F20.8)') TEFF
      write (CLOGG,'(F20.8)') LOGG
      write (CMOH,'(F20.8)') MOH
      write (CMICRO,'(F20.8)') VMICRO

      CALL SYSTEM ( "jktld " // CTEFF // " " //  CLOGG // " " // CMOH
     &                              // " " // CMICRO // " " // OUTFILE )

      END SUBROUTINE TASK1
!=======================================================================
      SUBROUTINE TASK2 (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY)
                                          ! Produces a model light curve
      implicit none
      real*8 V(540)                       ! IN: light  curve  parameters
      integer NSINE,PSINE(20)             ! IN: number sines, which pars
      integer NPOLY,PPOLY(50)             ! IN: number polys, which pars
      integer VARY(540)                   ! IN: parameters vary or fixed
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      integer i                           ! LOCAL: loop counter
      real*8 MAG,LP,LS                    ! LOCAL: EBOP/GETMODEL  output
      real*8 PHASE                        ! LOCAL:  phase for evaluation
      real*8 GETMODEL                     ! FUNCTION: evaluate the model
      real*8 HJD                          ! LOCAL: time to calculate for
      real*8 R1,R2                        ! LOCAL: radii of the 2 stars
      integer NPHASE                      ! LOCAL: number of phases todo

            ! NSINE=0 and NPOLY=0 and NUMINT=1 but don't need to do this

      V(19) = 1.0d0           ! Set period to 1.0
      V(20) = 0.0d0           ! Set Tzero to 0.0
      LP = 0.0d0
      LS = 0.0d0

      if ( V(2) >= 0.0d0 ) then
        R1 = V(2) / (1.0d0 + V(3))
        R2 = V(2) / (1.0d0 + (1.0d0/V(3)))
      else
        R1 = abs(V(2))
        R2 = V(3)
      end if

            ! Sort out the reflection coefficients.

      MAG=GETMODEL(V,VARY,LDTYPE,0,PSINE,0,PPOLY,V(20),1,LP,LS,1,0.0d0)

      if ( V(11) < -99.9 ) then
        V(11) = 0.4d0 * (LS/(1.0d0-V(15))) * R1**2
        write(6,'(A39,A41)') ">> The reflection coefficient for star ",
     &                     "1 was calculated from the system geometry"
      else
        write(6,'(A39,A41)') ">> The reflection coefficient used for ",
     &                     "star 1 was that given in the input file. "
      end if

      if ( V(12) < -99.9 ) then
        V(12) = 0.4d0 * (LP/(1.0d0-V(15))) * R2**2
        write(6,'(A39,A41)') ">> The reflection coefficient for star ",
     &                     "1 was calculated from the system geometry"
      else
        write(6,'(A39,A41)') ">> The reflection coefficient used for ",
     &                     "star 2 was that given in the input file. "
      end if

            ! We want 10000 phase points in most cases, but if one
            ! component is very small we want a factor of 10 more.

      NPHASE = 10001
      if ( R1 < 0.01d0 .or. R2 < 0.01d0 ) NPHASE = 100001

            ! Now calculate the light curve as a function of orbital
            ! phase and write each point into the output file.

      write(62,'(A47)')"#  PHASE  MAGNITUDE    L1         L2         L3"
      do i = 1,NPHASE
        PHASE = (i-1) / dble(NPHASE-1)
        HJD = V(20) + PHASE * V(19)
        MAG =GETMODEL(V,VARY,LDTYPE,0,PSINE,0,PPOLY,HJD,1,LP,LS,1,0.0d0)
        write (62,'(F8.6,4(1X,F10.6))') PHASE,MAG,LP,LS,V(15)
      end do
      close (62)

      END SUBROUTINE TASK2
!=======================================================================
      SUBROUTINE ERRFUDGE (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,
     &                   NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,IFAIL,NSINE,
     &                 PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

            ! This subroutine  modifies the error bars  of the LC and RV
            ! datasets to force the reduced chi-squared for each dataset
            ! versus the best fit to equal 1.0. It is optionally invoked
            ! if the command  "chif" (i.e. chi^2 fudge)  is specified in
            ! the input file. It only occurs for light curve datasets of
            ! 10 or more points and RV datasets of 5 or more points. The
            ! procedure occurs iteratively, to cope with the possibility
            ! that adjusting the errorbars for multiple datasets changes
            ! the best fit significantly.

      implicit none
      real*8 DATX(MAXDATA)                ! IN: Data independent varible
      real*8 DATY(MAXDATA)                ! IN: Data dependent variable
      real*8 DATERR(MAXDATA)              ! IN/OUT: Data errorbars
      integer DTYPE(MAXDATA)              ! IN: type of each datapoint
      integer NDATA,MAXDATA               ! IN: number, max number of dp
      integer NLR,NMIN,NL3,NECW,NESW      ! IN: numbers of types of data
      real*8 V(540)                       ! IN: light  curve  parameters
      integer VARY(540)                   ! IN: parameters vary or fixed
      real*8 VERR(540)                    ! SUB: parameter formal errors
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      integer ITER,IFAIL                  ! IN: number iterations + flag
      integer NSINE,PSINE(20)             ! IN: Number of sines and pars
      integer NPOLY,PPOLY(50)             ! IN: Number of polys and pars
      integer NUMINT                      ! IN: Number of numerical ints
      real*8  NINTERVAL                   ! IN: Time interval for numint
      real*8 CHISQ                        ! OUT: chi-square of model fit
      real*8 CHISQIN                      ! LOCAL: local chi-squared val
      real*8 CHISQLC,CHISQRV1,CHISQRV2    ! LOCAL: chisqred for datasets
      real*8 RES2SUMLC,RES2SUMRV1,RES2SUMRV2 ! LOCAL: sum squared resids
      integer NLC,NRV1,NRV2               ! LOCAL: number of dataset dps
      integer i,j                         ! LOCAL: loop counters
      real*8 MAG,LP,LS                    ! LOCAL: EBOP, GETMODEL output
      real*8 GETMODEL                     ! FUNCTION: evaluate the model

            ! First set up each iteration by zeroing the sum variables

      do i = 1,9         ! Maximum of 9 iterations to avoid endless loop
        NLC = 0
        NRV1 = 0
        NRV2 = 0
        CHISQLC = 0.0d0
        CHISQRV1 = 0.0d0
        CHISQRV2 = 0.0d0
        RES2SUMLC  = 0.0d0
        RES2SUMRV1 = 0.0d0
        RES2SUMRV2 = 0.0d0

            ! Evaluate the model to obtain the sum of the squared resid-
            ! uals and the sum of the chi-squared, for each of the three
            ! possible datasets.

        do j = 1,NDATA
          MAG = GETMODEL(V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,DATX(j),
     &                                DTYPE(j),LP,LS,NUMINT,NINTERVAL)
          CHISQIN = ( (DATY(j)-MAG) / DATERR(j) )**2
          if ( DTYPE(j) == 1 ) then
            RES2SUMLC  = RES2SUMLC  + (DATY(j)-MAG)**2
            CHISQLC = CHISQLC + CHISQIN
            NLC = NLC + 1
          end if
          if ( DTYPE(j) == 7 ) then
            RES2SUMRV1  = RES2SUMRV1  + (DATY(j)-MAG)**2
            CHISQRV1 = CHISQRV1 + CHISQIN
            NRV1 = NRV1 + 1
          end if
          if ( DTYPE(j) == 8 ) then
            RES2SUMRV2  = RES2SUMRV2  + (DATY(j)-MAG)**2
            CHISQRV2 = CHISQRV2 + CHISQIN
            NRV2 = NRV2 + 1
          end if
        end do

            ! Calculate the overall rms and the sqrt(CHISQRED) for each.

        RES2SUMLC = sqrt(RES2SUMLC / dble(NLC))
        RES2SUMRV1 = sqrt(RES2SUMRV1 / dble(NRV1))
        RES2SUMRV2 = sqrt(RES2SUMRV2 / dble(NRV2))
        CHISQLC = sqrt(CHISQLC / dble(NLC))
        CHISQRV1 = sqrt(CHISQRV1 / dble(NRV1))
        CHISQRV2 = sqrt(CHISQRV2 / dble(NRV2))

            ! Apply them to the errorbars if there are sufficient points
            ! for each datasets individually. The rms is used as the
            ! errorbar if there were no errorbars. If errorbars already
            ! exist then they are scaled by sqrt(CHISQRED) so chisqred=1

        do j = 1,NDATA
          if ( DTYPE(j) == 1 .and. NLC >= 10 ) then
            if ( DATERR(j) <= 0.0d0 ) DATERR(j) = RES2SUMLC
            if ( DATERR(j) > 0.0d0 ) DATERR(j) = DATERR(j) * CHISQLC
          else if ( DTYPE(j) == 7 .and. NRV1 >= 5 ) then
            if ( DATERR(j) <= 0.0d0 ) DATERR(j) = RES2SUMRV1
            if ( DATERR(j) > 0.0d0 ) DATERR(j) = DATERR(j) * CHISQRV1
          else if ( DTYPE(j) == 8 .and. NRV2 >= 5 ) then
            if ( DATERR(j) <= 0.0d0 ) DATERR(j) = RES2SUMRV2
            if ( DATERR(j) > 0.0d0 ) DATERR(j) = DATERR(j) * CHISQRV2
          end if
        end do

            ! Now we have new errorbars we recalculate the fit.
            ! Then print out the CHISQ values and continue.

        CALL FITEBOP (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,V,
     &                VARY,LDTYPE,ITER,CHISQ,VERR,IFAIL,NSINE,PSINE,
     &                NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

        if ( NLC >= 10 .or. NRV1 >= 5 .or. NRV2 >= 5 ) then
          write (6,'(A18,I1,A38,$)') ">> Done iteration ",i,
     &                          " to adjust errorbars. Chisqred values:"
          if ( NLC >= 10 ) write (6,'(1X,F7.3,$)') CHISQLC**2
          if ( NRV1 >= 5 ) write (6,'(1X,F7.3,$)') CHISQRV1**2
          if ( NRV2 >= 5 ) write (6,'(1X,F7.3,$)') CHISQRV2**2
          write (6,'(A1)') " "
        end if

            ! Check if the relevant chisq values are sufficiently close
            ! to 1.0 for good errorbars; exit loop if this is so.

        if ( NLC >= 10 .and. NRV1 >= 5 .and. NRV2 >= 5 .and.
     &       CHISQLC > 0.999 .and. CHISQLC < 1.001 .and.
     &       CHISQRV1 > 0.999 .and. CHISQRV1 < 1.001 .and.
     &       CHISQRV2 > 0.999 .and. CHISQRV2 < 1.001 ) exit

        if ( NLC >= 10 .and. NRV1 >= 5 .and. NRV2 < 5 .and.
     &       CHISQLC > 0.999 .and. CHISQLC < 1.001 .and.
     &       CHISQRV1 > 0.999 .and. CHISQRV1 < 1.001 ) exit

        if ( NLC >= 10 .and. NRV1 < 5 .and. NRV2 >= 5 .and.
     &       CHISQLC > 0.999 .and. CHISQLC < 1.001 .and.
     &       CHISQRV2 > 0.999 .and. CHISQRV2 < 1.001 ) exit

        if ( NLC >= 10 .and. NRV1 < 5 .and. NRV2 < 5 .and.
     &       CHISQLC > 0.999 .and. CHISQLC < 1.001 ) exit

      end do

      END SUBROUTINE ERRFUDGE
!=======================================================================
      SUBROUTINE TASK34 (TASK,V,VARY,LDTYPE,DATX,DATY,DATERR,DTYPE,
     &                   NDATA,MAXDATA,DATAFORM,NLR,NMIN,SIGMA,NSINE,
     &                 PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

            ! Find the best-fitting light curve parameters for the data.
            ! If TASK=4 it iteratively rejects all datapoints which are
            ! > SIGMA sigma from the best fit and refits the remainder.

      implicit none
      integer MAXDATA                     ! IN: max number of datapoints
      character DATAFORM*1                ! IN: MAXDATA character length
      integer TASK                        ! IN: which task to do
      real*8 V(540)                       ! IN: light  curve  parameters
      integer VARY(540)                   ! IN: parameters vary or fixed
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      real*8 DATX(MAXDATA)                ! IN: Data independent varible
      real*8 DATY(MAXDATA)                ! IN: Data dependent variable
      real*8 DATERR(MAXDATA)              ! IN: Data errorbars
      integer DTYPE(MAXDATA)              ! IN: type of each datapoint
      integer NDATA,NLR,NMIN              ! IN: number of  types of data
      integer NSINE,NL3                   ! IN: Numbers of sines and L3s
      integer NECW,NESW                   ! IN: Numbers of ecosw & esinw
      integer PSINE(20)                   ! IN: Which par for each sine
      integer NPOLY,PPOLY(50)             ! IN: Similar for polynomials
      integer NUMINT                      ! IN: Number of numerical ints
      real*8  NINTERVAL                   ! IN: Time interval for numint
      real*8 SIGMA                        ! IN:  std. dev. for  clipping
      real*8 CHISQ                        ! SUB: chi-square of model fit
      real*8 VERR(540)                    ! SUB: parameter formal errors
      real*8 OMC(MAXDATA)                 ! LOCAL: (O-C) residual values
      real*8 SIG                          ! LOCAL: rms of the O-C values
      real*8 RESIDSQSUM                   ! LOCAL: sum of resid. squares
      integer ITER,IFAIL                  ! LOCAL: iter number & success
      real*8 MAG,LP,LS                    ! LOCAL: EBOP/GETMODEL  output
      integer KEEP(MAXDATA),ACOUNT        ! LOCAL: Datapoint bookkeeping
      integer i,j                         ! LOCAL: loop counter
      integer NREJ                        ! LOCAL: no of rejected values
      real*8 GETMODEL                     ! FUNCTION: evaluate the model

      do i = 1,MAXDATA
        OMC(i) = 0.0d0
        KEEP(i) = 0
      end do
      SIG = 0.0d0
      NREJ = 0

            ! Output initial values and then find the best fit.
            ! Optionally modify the errorbars to force chi^2_red = 1
            ! If TASK=3 then output various results and finish.

      CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,
     &             NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,0,NSINE,PSINE,
     &             NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      CALL FITEBOP (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,V,
     &              VARY,LDTYPE,ITER,CHISQ,VERR,IFAIL,NSINE,PSINE,
     &              NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      if ( V(14) > -100.0d0 .and. V(14) < -99.0d0 ) then
        if ( IFAIL == 0 .and. ITER < 100 ) then
          write(6,'(A33,$)') ">> Best fit has been found after "
          if ( ITER < 100 ) write(6,'(I2,$)') ITER
          if ( ITER >= 100 ) write(6,'(I3,$)') ITER
          write(6,'(A12)') " iterations."
        else
         write(6,'(A43,$)')"## Warning: a good fit was not found after "
          if ( ITER < 100 ) write(6,'(I2,$)') ITER
          if ( ITER >= 100 ) write(6,'(I3,$)') ITER
          write(6,'(A12)') " iterations."
        end if
        CALL ERRFUDGE (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,V,
     &                 VARY,LDTYPE,ITER,CHISQ,VERR,IFAIL,NSINE,PSINE,
     &                 NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      end if

      if ( TASK == 3 ) then
        CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,
     &               NMIN,V,VARY,LDTYPE,ITER,CHISQ, VERR,1,NSINE,PSINE,
     &               NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

        if ( IFAIL == 0 .and. ITER < 100 ) then
          write(6,'(A33,$)') ">> Best fit has been found after "
          if ( ITER < 100 ) write(6,'(I2,$)') ITER
          if ( ITER >= 100 ) write(6,'(I3,$)') ITER
          write(6,'(A12)') " iterations."
        else
         write(6,'(A43,$)')"## Warning: a good fit was not found after "
          if ( ITER < 100 ) write(6,'(I2,$)') ITER
          if ( ITER >= 100 ) write(6,'(I3,$)') ITER
          write(6,'(A12)') " iterations."
        end if

        return
      end if

            ! If TASK=4 then output the best fit and continue.

      CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,
     &             NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,2,NSINE,PSINE,
     &             NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      write(6,'(A40,A40)')  ">> Best fit to all data points has been ",
     &                       "found and written to the parameter file."

            ! Now, for nine iterations, calculate the scatter of the
            ! observations (if there are no observational errors). Then
            ! go through the data checking for ones with O-C values
            ! greater than SIGMA sigmas. Finish if an iteration does not
            ! result in a rejected datapoint. Remove the rejected datap-
            ! oints then refit the data. Note that the light ratios and
            ! minimum times always have uncertainties.

      do j = 1,9

            ! First calculate the rms of the light curve residuals if
            ! there are no observational uncertainties.

        if ( DATERR(1) < 0.0d0 ) then
          RESIDSQSUM = 0.0d0
          do i = 1,NDATA
            if ( DTYPE(i) == 1 ) then
              MAG = GETMODEL (V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                        DATX(i),DTYPE(i),LP,LS,NUMINT,NINTERVAL)
              OMC(i) = MAG - DATY(i)
              RESIDSQSUM = RESIDSQSUM + OMC(i)**2
            end if
          end do
          SIG = sqrt( RESIDSQSUM / dble(NDATA-NLR-NMIN) )
        else
          do i = 1,NDATA
            if ( DTYPE(i) == 1 ) then
              MAG = GETMODEL (V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                        DATX(i),DTYPE(i),LP,LS,NUMINT,NINTERVAL)
              OMC(i) = MAG - DATY(i)
            end if
          end do
        end if

          ! Now put the array indices of the good datapoints into KEEP

        ACOUNT = 0
        do i = 1,NDATA                      ! if no observational errors
          if ( DTYPE(i) == 1 .and. DATERR(i) < 0.0d0 ) then
            if ( abs(OMC(i)) <= SIG*SIGMA ) then
              ACOUNT = ACOUNT + 1
              KEEP(ACOUNT) = i
            end if
          else                            ! if ob'l errors were supplied
            if ( abs(OMC(i)/DATERR(i)) <= SIGMA ) then
              ACOUNT = ACOUNT + 1
              KEEP(ACOUNT) = i
            end if
          end if
        end do

            ! Now keep only those datapoints which are specified by KEEP
            ! Have to recount number of light ratios and minimum times

        do i = 1,ACOUNT
          DATX(i) = DATX(KEEP(i))
          DATY(i) = DATY(KEEP(i))
          DATERR(i) = DATERR(KEEP(i))
          DTYPE(i) = DTYPE(KEEP(i))
        end do

        NREJ = NDATA - ACOUNT
        write(6,'(A13,I1,A1,I5,A8,I'//DATAFORM//',A47)')
     &                      ">> Iteration ",j,":",NREJ," out of ",NDATA,
     &                 " datapoints have been rejected from the dataset"

            ! Now recount the numbers and types of data and refit them.

        NDATA = ACOUNT
        NLR = 0
        NMIN = 0
        do i = 1,NDATA
          if ( DTYPE(i) == 2 ) NLR = NLR + 1
          if ( DTYPE(i) == 3 ) NMIN = NMIN + 1
        end do

        CALL FITEBOP (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,V,
     &                VARY,LDTYPE,ITER,CHISQ,VERR,IFAIL,NSINE,PSINE,
     &                NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
        if ( NREJ < 1 ) exit
      end do

      CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,
     &             NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,1,NSINE,PSINE,
     &             NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      write(6,'(A32,I3,A45)') ">> Best fit has been found from ",ITER,
     &                 " iterations and written to the parameter file"

      END SUBROUTINE TASK34
!=======================================================================
      SUBROUTINE TASK6 (V,VARY,LDTYPE,DATX,DATY,DATERR,DTYPE,NDATA,
     &                  MAXDATA,DATAFORM,NLR,NMIN,NSINE,PSINE,NPOLY,
     &                  PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

            ! Fits model parameters to an observed light curve. Then
            ! investigates each adjustable parameter by fixing it at a
            ! range of values round the best fit and seeing what the
            ! effect is on the other parameters.

      implicit none
      integer MAXDATA                     ! IN: max number of datapoints
      character DATAFORM*1                ! IN: MAXDATA character length
      real*8 V(540)                       ! IN: The light curve params
      integer VARY(540)                   ! IN: Par adjustment integers
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      real*8 DATX(MAXDATA)                ! IN: Data independent varible
      real*8 DATY(MAXDATA)                ! IN: Data dependent variables
      real*8 DATERR(MAXDATA)              ! IN: Data errorbars
      integer DTYPE(MAXDATA)              ! IN: Type of each datapoint
      integer NDATA,NLR,NMIN              ! IN: Numbers of datapoints
      integer NSINE,NL3                   ! IN: Numbers of sines and L3s
      integer NECW,NESW                   ! IN: Numbers of ecosw & esinw
      integer PSINE(20)                   ! IN: Which par for each sine
      integer NPOLY,PPOLY(50)             ! IN: Similar for polynomials
      integer NUMINT                      ! IN: Number of numerical ints
      real*8  NINTERVAL                   ! IN: Time interval for numint
      integer NUMVARY                     ! LOCAL: Number of vary params
      integer VWHERE(540)                 ! LOCAL: Which parameters vary
      real*8 VSTORE(540)                  ! LOCAL: Store best-fit params
      integer VARYFLAG                    ! LOCAL: Store a VARY integer
      real*8 VERR(540),CHISQ              ! LOCAL: Output from FITEBOP
      integer i,j,k,IFAIL,ITER            ! LOCAL: Loop counters etc

      CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,
     &             NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,0,NSINE,PSINE,
     &             NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      CALL FITEBOP (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,V,
     &              VARY,LDTYPE,ITER,CHISQ,VERR,IFAIL,NSINE,PSINE,
     &              NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      if ( V(14) > -100.0d0 .and. V(14) < -99.0d0 ) then
        CALL ERRFUDGE (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,V,
     &                 VARY,LDTYPE,ITER,CHISQ,VERR,IFAIL,NSINE,PSINE,
     &                 NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      end if

      CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,
     &             NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,2,NSINE,PSINE,
     &             NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      write(6,'(A40,A40)') ">> Fitting process completed. Best fit f",
     &                     "ound and output to the parameter file.  "

      NUMVARY = 0
      j = 1
      do i = 1,540
        if ( i /= 11 .and. i /= 12 ) then
          if ( VARY(i) /= 0 ) then
            NUMVARY = NUMVARY + 1
            VWHERE(j) = i
            j = j + 1
          end if
        end if
      end do

      do i = 1,540
        VSTORE(i) = V(i)
      end do

            ! Parametrs with VARY=1 are adjustable and those with VARY=2
            ! are fixed when finding the best fit but are perturbed like
            ! adjustable parameters here.
            ! Foreach perturbed parameter store its value and adjustment
            ! integer, then step through various values whilst adjusting
            ! all adjustable parameters to find the best fit.   Start at
            ! the best-fit value  and gradually step away from it in one
            ! direction. Then do the same for the other direction.

      CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,
     &             NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,0,NSINE,PSINE,
     &             NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      CALL FITEBOP (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,V,
     &              VARY,LDTYPE,ITER,CHISQ,VERR,IFAIL,NSINE,PSINE,
     &              NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      do i = 1,NUMVARY
        j = VWHERE(i)
        do k = 1,540
          V(k) = VSTORE(k)
        end do
        j = VWHERE(i)
        if ( j /= 16 .and. j /= 17 ) then
          VARYFLAG = VARY(j)
          VARY(j) = 0
          CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,
     &                 NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,100+j,
     &           NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

          do k = 0,80
            if ( j == 6 ) V(j) = VSTORE(j) - 0.1d0*k
            if ( j /= 6 ) V(j) = VSTORE(j) * (1.0d0 - k/40.0d0)
            CALL FITEBOP (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,
     &                  V,VARY,LDTYPE,ITER,CHISQ,VERR,IFAIL,NSINE,PSINE,
     &                       NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
            CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,
     &               NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,200+j,NSINE,
     &                 PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
          end do

          do k = 1,540
            V(k) = VSTORE(k)
          end do

          do k = 1,540
            if ( j == 6 ) V(j) = VSTORE(j) + 0.1d0*k
            if ( j /= 6 ) V(j) = VSTORE(j) * (1.0d0 + k/40.0d0)
            if ( j /= 6 .or. V(6) < 90.1d0 ) then
              CALL FITEBOP (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,
     &                   NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,IFAIL,NSINE,
     &                 PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
              CALL OUTPUT(DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,
     &               NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,200+j,NSINE,
     &                 PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
            end if
          end do

          VARY(j) = VARYFLAG
        end if

      end do

      END SUBROUTINE TASK6
!=======================================================================
      SUBROUTINE TASK5789 (TASK,V,VARY,LDTYPE,DATX,DATY,DATERR,DTYPE,
     &                     NDATA,MAXDATA,DATAFORM,NLR,NMIN,NSIM,NSINE,
     &                 PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
!
!   This big subroutine executes four of the JKTEBOP tasks:
!
! TASK 5:  this has some characteristics of a global search. Once a best
!   fit has been found the params are perturbed and refitted many times.
!   Each result is outputted  and the best one is given as the solution.
!
! TASK 7:  this performs a bootstrapping error analysis. Once a best fit
!   is found many new datasets are made by copying (with replacement) at
!   random from the actual light curve. Each one is fitted to find 1-sig
!   errors. All light ratios and times of minimum are kept each time.
!
! TASK 8:  this does a Monte Carlo simulation error analysis.  Once best
!   fit has been found, the model is evaluated  at the phases of the obs
!   to make a synthetic light curve. Then, for each simulation, Gaussian
!   noise is added and the data are refitted to find the 1-sigma spread.
!
! TASK 9:  this does a Monte Carlo simulation error analysis  similar to
!   TASK 8,  but instead of adding simulated Gaussian noise to the data,
!   the residuals of the fit are used.    For each simulation the set of
!   residuals is moved  and applied to the model values for the adjacent
!   datapoints. Residuals which drop off the end of the data are wrapped
!   around to the start of the dataset.   Thus the number of simulations
!   done is (NDATA - 1). This has the advantage that correlated noise in
!   the residuals is retained and affects the quality of the fit.
!
! For TASKS 7, 8, 9,  the starting pars are perturbed each time to avoid
!   finding unrealistically small errors in a local minimum.   For TASKS
!   7 and 8 this is optional  (won't be done if NSIM is less than zero).
!   Can also be turned off for individual parameters by putting VARY=3
!
! Monte Carlo simulations are an excellent error indicator,   but do not
!   properly account for systematic errors.     The residual permutation
!   (TASK 9) is better. Bootstrapping is the best, but is likely to give
!   pesimistic results due to the loss of time sampling in the simulated
!   datasets. Correct solution to problem: model and remove systematics.
!
! For information on bootstrapping, Monte Carlo, and minimisation algor-
! ithms, read Numerical Recipes in Fortran 77 (Press et al 1993) chap.15

      implicit none
      integer MAXDATA               ! IN: maximum number of datapoints
      character DATAFORM*1          ! IN: character length of MAXDATA
      integer TASK                  ! IN: The task to undertake
      real*8 V(540)                 ! IN: The photometric parameters
      integer VARY(540)             ! IN: Parameter adjustment integers
      integer LDTYPE(2)             ! IN: Type of LD law for each star
      real*8 DATX(MAXDATA)          ! IN: Data independent variables
      real*8 DATY(MAXDATA)          ! IN: Data dependent variables
      real*8 DATERR(MAXDATA)        ! IN: Data errorbars
      integer DTYPE(MAXDATA)        ! IN: Types of observational data
      integer NDATA,NLR,NMIN        ! IN: Numbers of different datatypes
      integer NSINE,NL3             ! IN: Numbers of sines and L3s
      integer NECW,NESW             ! IN: Numbers of ecosw and esinw
      integer PSINE(20)             ! IN: Which parameter for each sine
      integer NPOLY,PPOLY(50)       ! IN: Similar for the polynomials
      integer NSIM                  ! IN: Number of simulations to do
      integer NUMINT                ! IN: Number of numerical integrat's
      real*8  NINTERVAL             ! IN: Time interval numerical integ.
      character PERTURB*1           ! LOCAL: Perturb params ('y' or 'n')
      real*8 VERR(540)              ! SUB:   Param  formal uncertainties
      real*8 DV(540)                ! LOCAL: Perturbations for paramters
      real*8 VEXTRA(18)             ! SUB:   Extra (dependent) parametrs
      real*8 VALL(540,abs(NSIM))    ! SUB:   All the best-fitting params
      real*8 VALLEXTRA(18,abs(NSIM))! SUB:   All extra (dependent) pars
      integer NVARY                 ! LOCAL: Number of varying parametrs
      integer ITER                  ! SUB:   Numberof fitting iterations
      real*8 INX(MAXDATA)           ! LOCAL: Synthetic obs to be fitted
      real*8 INY(MAXDATA)           ! LOCAL: Synthetic obs to be fitted
      real*8 INERR(MAXDATA)         ! LOCAL: Synthetic obs to be fitted
      real*8 RESID(MAXDATA)         ! LOCAL: Residuals of the best fit
      real*8 INV(540)               ! SUB:   Parameters of  synthetic LC
      real*8 INVEXTRA(18)           ! SUB:   Dependent parameter array
      real*8 CHISQ,ERRSIZELC        ! SUB:   Chi-squared of the best fit
      real*8 ERRSIZERV1,ERRSIZERV2  ! SUB:   Errorbar sizes for datasets
      integer SEEDSTART,SEED        ! FUNCTION: Start randomnumber maker
      real*8 RANDOMG                ! FUNCTION: Gaussian  random  number
      real*8 RANDOM                 ! FUNCTION: Flat-distrib  random num
      real*8 GETMODEL               ! FUNCTION: Gets  model  predictions
      real*8 LP,LS                  ! LOCAL: EBOP/GETMODEL output values
      real*8 HELP1                  ! LOCAL: Useful variable storage
      integer i,j,k,m,ERROR         ! LOCAL: Loop counters + error flags
      character*5 NAME5(540)        ! LOCAL: short names of parameters
      character*5 NAMEXTRA5(18)     ! LOCAL: short names of extra params
      character*19 NAME19(540)      ! LOCAL: long names of parameters
      character*19 NAMEXTRA19(18)   ! LOCAL: long names of extra params
      integer DOPRINT(540)          ! OUT:   1 to print fitted parameter
      integer DOPRINTEXTRA(18)      ! OUT:   1 to print dependent param
      real*8 FITCHISQ               ! LOCAL: chisq of original best fit
      integer ILC(MAXDATA)          ! LOCAL: array indices of phot data
      integer IRV1(MAXDATA)         ! LOCAL: array indices of RV1 data
      integer IRV2(MAXDATA)         ! LOCAL: array indices of RV2 data
      integer ITMIN(MAXDATA)        ! LOCAL: array i of times of minimum
      real*8 HELP2,HELP3            ! LOCAL: helper variables
      integer IMOD                  ! LOCAL: helper variable
      real*8 MODELVAL(MAXDATA)      ! LOCAL: model value at data times
      integer ILCFIRST,ILCLAST      ! LOCAL: where LC data are in array
      integer IRV1FIRST,IRV1LAST    ! LOCAL: where RV1 data are in array
      integer IRV2FIRST,IRV2LAST    ! LOCAL: where RV2 data are in array
      integer ITMINFIRST,ITMINLAST  ! LOCAL: where tmin data arein array
      integer NLC,NRV1,NRV2,NTMIN   ! LOCAL: numbers of types of data

      CALL GETNAME(V,NAME5,NAMEXTRA5,NAME19,NAMEXTRA19)

            ! First get DV values (how much to perturb parameters by).
            ! Set 10x larger than the numerical derivative intervals.
            ! Check whether perturbations should be applied.
            ! Set perturbations a factor of 50 larger for TASK 5.

      if ( NSIM < 0 ) then
        PERTURB = "n"
        NSIM = abs(NSIM)
      else
        PERTURB = "y"
      end if
      if (task == 9) perturb="n"

      CALL GET_DV(V,DV,NPOLY,PPOLY)
      do i = 1,540
        DV(i) = DV(i) * 10.0d0
      end do

      if ( TASK == 5 ) then
        do i = 1,540
          DV(i) = DV(i) * 5.0d0
        end do
      end if

            ! Now find the best fit to the light curve and output it.

      CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,
     &             NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,0,NSINE,PSINE,
     &             NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      CALL FITEBOP (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,V,
     &              VARY,LDTYPE,ITER,CHISQ,VERR,ERROR,NSINE,PSINE,
     &              NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      if ( V(14) > -100.0d0 .and. V(14) < -99.0d0 ) then
        write(6,'(A39,A41)')  ">> Fitting process completed. Now inves",
     &                     "tigating the errorbars for the datasets.  "
        CALL ERRFUDGE (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,V,
     &                 VARY,LDTYPE,ITER,CHISQ,VERR,ERROR,NSINE,PSINE,
     &                 NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      end if

      FITCHISQ = CHISQ

      CALL GETEXTRAS (V,VARY,NDATA,LDTYPE,CHISQ,NSINE,PSINE,NPOLY,PPOLY,
     &                                          NUMINT,NINTERVAL,VEXTRA)
      if ( DATERR(1) <= 0.0d0 ) VEXTRA(8) = -999.0d0

      CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,
     &             NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,2,NSINE,PSINE,
     &             NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      write(6,'(A40,A40)')  ">> Fitting process completed. Best fit i",
     &                      "s found and outputted to parameter file."

            ! Need to set up arrays for each type of data for the resid-
            ! ual permutation simulations. This allows the residuals to
            ! be permuted within each dataset, rather than over all data
            ! sets, which makes much more mathematical sense. The photo-
            ! metric and RV datapoints are permuted. Light ratios, tmin,
            ! ecosw, esinw and L3 values are not.

      do i = 1,NDATA
        MODELVAL(i) = GETMODEL (V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                          DATX(i),DTYPE(i),LP,LS,NUMINT,NINTERVAL)
        RESID(i) = DATY(i) - MODELVAL(i)
      end do

      do i = 1,MAXDATA
        ILC(i) = 0
        ITMIN(i) = 0
        IRV1(i) = 0
        IRV2(i) = 0
      end do

      NLC = 0
      NTMIN = 0
      NRV1 = 0
      NRV2 = 0
      ILCFIRST = MAXDATA
      ILCLAST = 0
      ITMINFIRST = MAXDATA
      ITMINLAST = 0
      IRV1FIRST = MAXDATA
      IRV1LAST = 0
      IRV2FIRST = MAXDATA
      IRV2LAST = 0

      do i = 1,NDATA
        if ( DTYPE(i) == 1 ) then
          NLC = NLC + 1
          ILC(NLC) = i
          if ( ILCFIRST > i ) ILCFIRST = i
          if ( ILCLAST < i ) ILCLAST = i
        end if
        if ( DTYPE(i) == 3 ) then
          NTMIN = NTMIN + 1
          ITMIN(NTMIN) = i
          if ( ITMINFIRST > i ) ITMINFIRST = i
          if ( ITMINLAST < i ) ITMINLAST = i
        end if
        if ( DTYPE(i) == 7 ) then
          NRV1 = NRV1 + 1
          IRV1(NRV1) = i
          if ( IRV1FIRST > i ) IRV1FIRST = i
          if ( IRV1LAST < i ) IRV1LAST = i
        end if
        if ( DTYPE(i) == 8 ) then
          NRV2 = NRV2 + 1
          IRV2(NRV2) = i
          if ( IRV2FIRST > i ) IRV2FIRST = i
          if ( IRV2LAST < i ) IRV2LAST = i
        end if
      end do

            ! Write the column headings  for those individual simulation
            ! results to be output later.   The complication is that the
            ! proliferation of parameters (now 540 fit pars and 18 depen-
            ! dent pars) yields extensive output, of which a significant
            ! amount is unimportant because not all parameters are being
            ! varied at that point. Thus restrict the output to the pars
            ! which are being fitted or are otherwise changing.

      CALL GETPRINTSIM (V,VARY,DOPRINT,DOPRINTEXTRA)

      write(63,'(A12,$)') "#   N  ITER "
      do i = 1,540
        if ( DOPRINT(i) /= 0 ) write (63,'(7X,A5,A7,$)') NAME5(i),
     &                                                   "       "
      end do
      do i = 1,18
        if (DOPRINTEXTRA(i) /= 0) write(63,'(7X,A5,A7,$)') NAMEXTRA5(i),
     &                                                     "       "
      end do

      write (63,*) " "

            ! Start the random number generator. Calculate the number of
            ! variable params.  Store original data and best-fit params.

      SEED = SEEDSTART ()

      NVARY = 0
      do i = 1,540
        if ( VARY(i) == 1 ) NVARY = NVARY + 1
      end do

      do i = 1,540
        INV(i) = V(i)
      end do

            ! If the inputted NSIM is below zero no perturbations should
            ! be applied to the initial  light curve parameter estimates
            ! prior to each MC simulation. This is useful if convergence
            ! only happens for a  very narrow range of parameter values.

      if ( TASK == 5 ) write (62,'(A23,A57)') "Task 5: initial paramet",
     &      "ers will be widely perturbed before each fit is done.    "
      if ( TASK == 7 .or. TASK == 8) then
        if (PERTURB == "n") write (62,'(A18,A60)') "The number of simu",
     &   "lations is below zero: the parameters will not be perturbed."
        if (PERTURB == "y") write (62,'(A18,A60)') "The number of simu",
     &   "lations is above zero: the parameters will be perturbed.    "
        write (62,*) " "
      end if
      if ( TASK == 9 ) write (62,'(A23,A57)') "Task 9: initial paramet",
     &      "ers will be perturbed before each fit is performed.      "

      if ( TASK == 5 ) write(6,'(A46,I6)')
     &             ">> Number of perturbed parameter refits to do:",NSIM
      if ( TASK == 7 ) write(6,'(A45,I6)')
     &              ">> Number of bootstrapping simulations to do:",NSIM
      if ( TASK == 8 ) write(6,'(A43,I6)')
     &                ">> Number of Monte Carlo simulations to do:",NSIM
      if ( TASK == 9 ) write(6,'(A41,I6)')
     &                  ">> Number of residual permutations to do:",NSIM

            ! For Monte Carlo simulations we must evaluate the best-fit-
            ! ting model at the observed phases. If there were no errors
            ! in the input light or radial velocity curve, calculate the
            ! rms of the residuals  and  use to scale the Gaussian noise
            ! added to the simulated data.

      if ( TASK == 8 ) then

        HELP1 = 0.0d0
        HELP2 = 0.0d0
        HELP3 = 0.0d0
        do i = 1,NDATA
          if ( DTYPE(i) == 1 ) HELP1 = HELP1 + RESID(i)**2
          if ( DTYPE(i) == 7 ) HELP2 = HELP2 + RESID(i)**2
          if ( DTYPE(i) == 8 ) HELP3 = HELP3 + RESID(i)**2
        end do
        ERRSIZELC = sqrt(HELP1 / dble(NLC-NVARY))
        ERRSIZERV1 = sqrt(HELP2 / dble(NRV1))
        ERRSIZERV2 = sqrt(HELP3 / dble(NRV2))

        do i = 1,NDATA
          if ( DTYPE(i) == 1 ) INERR(i) = ERRSIZELC
          if ( DTYPE(i) == 7 ) INERR(i) = ERRSIZERV1
          if ( DTYPE(i) == 8 ) INERR(i) = ERRSIZERV2
        end do

        if ( DATERR(ILCFIRST) <= 0.0d0 ) then
          write (62,'(A35,A45)') "No observational errors were suppli",
     &                 "ed with the input light curve.  Noise of size"
          write (62,'(F6.3,A12,A62)')abs(ERRSIZELC)*1.d3," mmag (stand",
     &  "ard error of best fit) will be added to the synthetic datasets"
          write (6,'(A17,F6.3,A50)') ">> Noise of size ",abs(ERRSIZELC)*
     &         1.d3," mmag will be added to the synthetic light curves."
          VEXTRA(8) = 1.0d0
          do i = 1,NDATA
            if ( DTYPE(i) == 1 ) DATERR(i) = ERRSIZELC
          end do
        else
          write (62,'(A36,A44)') "Observational errors were supplied w",
     &                   "ith the input light curve.  These have been "
          write (62,'(A36,A44)') "assumed to be correct and used to se",
     &                   "t the size of the simulated Gaussian noise. "
        end if

        if ( NRV1 >= 1 ) then
          if ( DATERR(IRV1FIRST) <= 0.0d0 ) then
            write (62,'(A34,A46)') "No observational errors were suppl",
     &                  "ied with the input star A RVs.  Noise of size "
            write (62,'(F6.3,A13,A60)') abs(ERRSIZERV1)," km/s (standa",
     &   "rd error of best fit) will be added to the synthetic datasets"
            do i = 1,NDATA
              if ( DTYPE(i) == 7 ) DATERR(i) = ERRSIZERV1
            end do
          else
            write (62,'(A34,A46)') "Observational errors were supplied",
     &                 " with the input star A RVs.   These have been "
            write (62,'(A34,A46)') "assumed to be correct and used to ",
     &                 "set the size of the simulated Gaussian noise. "
          end if
        end if

        if ( NRV2 >= 1 ) then
          if ( DATERR(IRV2FIRST) <= 0.0d0 ) then
            write (62,'(A34,A46)') "No observational errors were suppl",
     &                  "ied with the input star B RVs.  Noise of size "
            write (62,'(F6.3,A15,A58)') abs(ERRSIZERV2)," km/s (standa",
     &   "rd error of best fit) will be added to the synthetic datasets"
            do i = 1,NDATA
              if ( DTYPE(i) == 8 ) DATERR(i) = ERRSIZERV2
            end do
          else
            write (62,'(A34,A46)') "Observational errors were supplied",
     &                 " with the input star B RVs.   These have been "
            write (62,'(A34,A46)') "assumed to be correct and used to ",
     &                 "set the size of the simulated Gaussian noise. "
          end if
        end if

      end if

!-----------------------------------------------------------------------

      write(6,'(A13,$)') ">> Completed:"

      do i = 1,abs(NSIM)             ! Loop over number of simulations

500     continue                     ! Enter here if previous sim failed

            ! TASK 5: use actual dataset for every simulation

        if ( TASK == 5 ) then
          do j = 1,NDATA
            INX(j) = DATX(j)
            INY(j) = DATY(j)
            INERR(j) = DATERR(j)
          end do
        end if

            ! TASK 7: randomly sample the observations with replacement
            ! to create a new light curve to fit. Don't do this with the                        ! Need to deal with the RVs here
            ! light ratios or minimum times etc (at end of DATA array) -
            ! these should be preserved as they are.

        if ( TASK == 7 ) then
          do j = 1,NDATA
            if ( DTYPE(j) == 1 ) then
              do m = 1,100000
                k = 1 + int( random(SEED) * (dble(NDATA)-0.000001d0) )
                if ( DTYPE(k) == 1 ) exit
              end do
              INX(j) = DATX(k)
              INY(j) = DATY(k)
              INERR(j) = DATERR(k)
            else
              INX(j) = DATX(j)
              INY(j) = DATY(j)
              INERR(j) = DATERR(j)
            end if
          end do
        end if

            ! TASK 8: create a simulated light curve by adding Gaussian
            ! noise to the best-fit model light curve.

        if ( TASK == 8 ) then
          do j = 1,NDATA
            INX(j) = DATX(j)
            INY(j) = MODELVAL(j) + randomg(SEED,0.0d0,DATERR(j))
            INERR(j) = DATERR(j)
          end do
        end if

            ! TASK 9: create simulated dataset (not for the times of
            ! minimum and spectroscopic light ratios) by taking the
            ! model fit for each datapoint and adding the residual of
            ! the best fit from the datapoint distant by index i.
            ! For the RV datasets, which I am assuming have fewer data-
            ! points than the LC dataset, the residuals are cycled
            ! multiple times by using the modulus function.

          if ( TASK == 9 ) then

          do j = 1,NDATA
            INX(j) = DATX(j)
            INY(j) = DATY(j)        ! needed for data not being permuted
            INERR(j) = DATERR(j)
          end do

          do j = 1,NLC-i
            INY(ILC(j)) = MODELVAL(ILC(j)) + RESID(ILC(j+i))
          end do
          do j = NLC-i+1,NLC
            INY(ILC(j)) = MODELVAL(ILC(j)) + RESID(ILC(j+i-NLC))
          end do

          if ( NRV1 >= 1 ) then
            IMOD = mod(i,NRV1)
            do j = 1,NRV1-IMOD
              INY(IRV1(j)) = MODELVAL(IRV1(j)) + RESID(IRV1(j+IMOD))
            end do
            do j = NRV1-IMOD+1,NRV1
              INY(IRV1(j)) = MODELVAL(IRV1(j)) +RESID(IRV1(j+IMOD-NRV1))
            end do
          end if

          if ( NRV2 >= 1 ) then
            IMOD = mod(i,NRV2)
            do j = 1,NRV2-IMOD
              INY(IRV2(j)) = MODELVAL(IRV2(j)) + RESID(IRV2(j+IMOD))
            end do
            do j = NRV2-IMOD+1,NRV2
              INY(IRV2(j)) = MODELVAL(IRV2(j)) +RESID(IRV2(j+IMOD-NRV2))
            end do
          end if

        end if

            ! Now perturb the initial values of the fitted parameters
            ! by 10 times the numerical derivative sizes
            ! Return to the original best-fit parameters without
            ! perturbation if PERTURB=="n", to avoid getting piles
            ! of rubbish results if the data are a bit dodgy.
            ! If DV(i) > 0 then is it a multiplicative perturbation.
            ! If DV(i) < 0 then is it an additive perturbation.

        do j = 1,540
          if ( (VARY(j)==1 .and. PERTURB=="y") .or. VARY(j)==2 ) then
            INV(j) = V(j) + (2.0d0*random(SEED)-1.0d0) * abs(DV(j))
          else
            INV(j) = V(j)
          end if
        end do

        if ( INV(6) > 90.0d0 ) GOTO 500         ! Avoid inclination > 90

            ! Now fit the simulated light curve, then catch some failed
            ! iterations by checking for some of the  common  problems.
            ! Don't catch TASK9 problems - this actually causes further
            ! attempts to fit *exactly* the same data, which inevitably
            ! yield the same result, and so the code gets stuck redoing
            ! the same failed fit over and over again.

        CALL FITEBOP (INX,INY,INERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,INV,
     &                VARY,LDTYPE,ITER,CHISQ,VERR,ERROR,NSINE,PSINE,
     &                NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

        if ( ERROR /= 0 )  GOTO 500

            ! Now catch *very* bad fits where k < 0, or the linear LD
            ! coefficients are < -5 or > +5, or the inclination is < 0
            ! or > 180, or the CHISQ is really awful. Don't record them
            ! except in the case of TASK9 (see above for the reason).

        if ( TASK /= 9 ) then
          if ( INV(3) <= 0.0d0 ) then
         write(6,'(A32,E12.5)')" Retry: V(3) is less than zero: ",INV(3)
            GOTO 500
          end if
          if ( INV(4) < -5.0d0 .or. INV(4) > 5.0d0 ) then
            write(6,'(A25,E12.5)') " Retry: V(4) is haywire: ",INV(4)
            GOTO 500
          end if
          if ( INV(5) < -5.0d0 .or. INV(5) > 5.0d0 ) then
            write(6,'(A25,E12.5)') " Retry: V(5) is haywire: ",INV(5)
            GOTO 500
          end if
          if ( INV(6) < 0.0d0 .or. INV(6) > 180.0d0 ) then
            write(6,'(A25,E12.5)') " Retry: V(6) is haywire: ",INV(6)
            GOTO 500
          end if
          if ( FITCHISQ > 0.0d0 .and. CHISQ/FITCHISQ > 10.0d0 ) then
            write(6,'(A17,E12.6,A1,E12.6,A20)') " CHISQ mismatch (",
     &                         CHISQ,",",FITCHISQ,"): reject iteration."
            GOTO 500
          end if
          if ( ITER >= 100 ) then
            write(6,'(A32,I4,A11)') " Retry: failed to converge after",
     &                                                ITER," iterations"
            GOTO 500
          end if
        end if

            ! Now store the final fitted parameters and calculate the
            ! the dependent parameters for storage as well.

        do j = 1,540
          VALL(j,i) = INV(j)
        end do

        CALL GETEXTRAS (INV,VARY,NDATA,LDTYPE,CHISQ,NSINE,PSINE,NPOLY,
     &                                  PPOLY,NUMINT,NINTERVAL,INVEXTRA)
        do j = 1,18
          VALLEXTRA(j,i) = INVEXTRA(j)
        end do

            ! Write the relevant results to the big output file.

        write(63,'(I5,1X,I3,$)') i,ITER
        do j = 1,540
          if ( DOPRINT(j) /= 0 ) write(63,'(1X,f18.10,$)') INV(j)
        end do
        do j = 1,18
          if (DOPRINTEXTRA(j) /= 0) write(63,'(1X,f18.10,$)')INVEXTRA(j)
        end do
        write (63,*) " "

            ! Output the simulation number to screen if required.

        if ( i < 10 ) write(6,'(1X,I1,$)') i
        if ( abs((dble(i)/10.0d0)-int(dble(i)/10.0d0)) < 0.0001d0 ) then
          if ( i >= 10 .and. i < 100 )   write(6,'(1X,I2,$)') i
          if ( i >= 100 .and. i < 1000 )  write(6,'(1X,I3,$)') i
          if ( i >= 1000 .and. i < 10000 ) write(6,'(1X,I4,$)') i
          if ( i >= 10000 )                 write(6,'(1X,I5,$)') i
        end if

            ! Check that a halt has not been requested.

        CALL STOPCHECK (ERROR)
        if ( ERROR /= 0 ) exit

      end do
      NSIM = i - 1

      write(6,*) " "

!-----------------------------------------------------------------------

            ! TASK 5: pick out the result with the lowest chi-square
            ! value, put its parameters into V, refit the result, and
            ! print the final quantities to the output file.

      if ( TASK == 5 ) then
        j = 0
        HELP1 = 1.0d8
        do i = 1,NSIM
          if ( VALLEXTRA(8,i) < HELP1 ) then
            j = i
            HELP1 = VALLEXTRA(8,i)
          end if
        end do

        do i = 1,540
          V(i) = VALL(i,j)
        end do
        CHISQ = HELP1
        CALL FITEBOP (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,V,
     &                VARY,LDTYPE,ITER,CHISQ,VERR,ERROR,NSINE,PSINE,
     &                NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

        write(62,'(A44)') "                                            "
        write(62,'(A44)') "-----------------------------------------   "
        write(62,'(A44)') "Overall best fit from all parameter sets:   "
        write(62,'(A44)') "-----------------------------------------   "
        write(62,'(A44)') "                                            "
        CALL OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,NLR,
     &               NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,3,NSINE,PSINE,
     &               NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      end if

            ! TASKS 7,8,9: call OUTPUTSIM to write out 1-sigma uncertai-
            ! nties in each of the fitted and dependent parameters.

      if ( TASK == 7 .or. TASK == 8 .or. TASK == 9 ) then

        write (62,*) " "
        write (62,'(A38,A56)') "======================================",
     &    "==========================================================="
        if ( TASK == 7 ) write (62,'(A13,I6,A27)')
     &          "Results from ",abs(NSIM)," bootstrapping simulations:"
        if ( TASK == 8 ) write (62,'(A13,I6,A25)')
     &            "Results from ",abs(NSIM)," Monte Carlo simulations:"
        if ( TASK == 9 ) write (62,'(A13,I6,A25)')
     &                 "Results from ",NSIM," residual-shifted fits:  "
        write (62,'(A38,A56)') "======================================",
     &    "==========================================================="
        write (62,*) " "


            ! Print out the 1-sigma results for each variable parameter

        write (62,'(A38,A53)') "              Best fit        one_sigm",
     &          "a         median         +68.27%         -68.27%     "
        CALL OUTPUTSIM (V,VARY,VEXTRA,VALL,VALLEXTRA,NSIM,0.6827d0)
        write (62,*) " "

            ! Also print out the 2-sigma results for reference

        if ( int(dble(NSIM)*(1.d0-0.9545d0)*0.5d0) > 1 ) then
          write (62,'(A36,A55)') "              Best fit        two_si",
     &        "gma         median         +95.45%         -95.45%     "
          CALL OUTPUTSIM (V,VARY,VEXTRA,VALL,VALLEXTRA,NSIM,0.9545d0)
          write (62,*) " "
        end if

            ! Also print out the 3-sigma results for reference

        if ( int(dble(NSIM)*(1.d0-0.9973d0)*0.5d0) > 1 ) then
          write (62,'(A36,A55)') "              Best fit        three_",
     &        "sigma       median         +99.73%         -99.73%     "
          CALL OUTPUTSIM (V,VARY,VEXTRA,VALL,VALLEXTRA,NSIM,0.9973d0)
          write (62,*) " "
        end if

      end if

      END SUBROUTINE TASK5789
!=======================================================================
!=======================================================================
      SUBROUTINE STOPCHECK (STATUS)  ! Checks if the file "jktebop.stop"
                                     ! exists and if so puts STATUS=1000
      implicit none
      integer ERROR,STATUS

      STATUS = 0
      ERROR = 100
      open (99,file="jktebop.stop",status="old",iostat=ERROR)
      close (99,status="delete")
      if ( ERROR == 0 ) STATUS = 1000

      END SUBROUTINE STOPCHECK
!=======================================================================
      SUBROUTINE OUTPUTSIM (V,VARY,VEXTRA,VALL,VALLEXTRA,NSIM,CONFINT)

            ! Given a set of Monte Carlo or bootstrapping simulation
            ! results, this subroutine calculates the confindence inter-
            ! val for each of the adjusted and dependent parameters.

      implicit none
      integer NSIM                  ! IN: The number of simulations done
      real*8 V(540)                 ! IN: The photometric model paramtrs
      integer VARY(540)             ! IN: Which parametrs fixed/adjusted
      real*8 VEXTRA(18)             ! IN: Extra  (dependent)  parameters
      real*8 VALL(540,abs(NSIM))    ! IN: All the best-fitting parametrs
      real*8 VALLEXTRA(18,abs(NSIM))! IN: All the extra (dependent) pars
      real*8 CONFINT                ! IN: Fractional confidence interval
      real*8 ARRAY(NSIM)            ! LOCAL: All results for one paramtr
      real*8 HELP1,HELP2,HELP3      ! LOCAL: Useful storage for variabls
      real*8 CONFERR,PERCENT        ! LOCAL: Error estimates for a param
      character*5 MESSAGE           ! LOCAL: Warning message
      integer i,j,k                 ! LOCAL: Loop counters
      real*8 SELLECT                ! FUNCTION: Gives nth array value
      character*5 NAME5(540),NAMEXTRA5(18)
      character*19 NAME19(540),NAMEXTRA19(18)
      integer DOPRINT(540),DOPRINTEXTRA(18)

      CALL GETNAME(V,NAME5,NAMEXTRA5,NAME19,NAMEXTRA19)

            ! Firstly construct an array  specifying which parameters to
            ! print out.  The ones of interest are the fitted parameters
            ! and some of the  extra parameters  depending on situation.
            ! Also specify DOPRINT=2 for those pars where the percentage
            ! value should also be printed out.

      CALL GETPRINTSIM (V,VARY,DOPRINT,DOPRINTEXTRA)

            ! Sort out the radii for the two possibilities of their pars

      if ( V(2) < 0.0d0 ) then
        do i = 1,NSIM
          VALL(2,i) = abs(VALL(2,i))
        end do
        V(2) = abs(V(2))
      end if

            ! Now for each adjusted parameter (and the other calculated
            ! quantities),   put the NSIM parameter evaluations into an
            ! array and select the median and 1_sigma bounds (using the
            ! median +/- 0.5_sigma) and output to the parameter file.

      do j = 1,558      ! Npar + Nextra = 540 + 18 = 566

        if ( j <= 540 ) then
          do k = 1,NSIM
            ARRAY(k) = mod(VALL(j,k),1.0d3)
          end do
        else
          do k = 1,NSIM
            ARRAY(k) = mod(VALLEXTRA(j-540,k),1.0d3)
          end do
        end if

        HELP1 = SELLECT(ARRAY,NSIM,int(dble(NSIM)*(1.d0-CONFINT)*0.5d0))
        HELP2 = SELLECT(ARRAY,NSIM,int(dble(NSIM)*0.5d0))
        HELP3 = SELLECT(ARRAY,NSIM,int(dble(NSIM)*(1.d0+CONFINT)*0.5d0))

        HELP1 = abs(HELP2-HELP1)        ! Minus CONFINT sigma error
        HELP3 = abs(HELP3-HELP2)        ! Plus  CONFINT sigma error

        MESSAGE = "     "               ! Warn if errors are asymmetric
        if ( HELP1<0.67d0*HELP3.or. HELP1>1.5d0*HELP3) MESSAGE = "skew?"
        if ( HELP1<0.5d0*HELP3 .or. HELP1>2.d0*HELP3 ) MESSAGE = "SKEW!"

        PERCENT = 50.d0 * abs( (abs(HELP1)+abs(HELP3)) / HELP2 )
        CONFERR = 0.5d0 * abs(HELP3+HELP1)

        if ( CONFERR > 0.0d0 ) then
          if ( j <= 540 ) then
            if (DOPRINT(j)==1) write(62,101) j,NAME5(j),mod(V(j),1.0d3),
     &                      CONFERR,mod(HELP2,1.0d3),HELP3,HELP1,MESSAGE
            if (DOPRINT(j)==2) write(62,100) j,NAME5(j),mod(V(j),1.0d3),
     &              CONFERR,mod(HELP2,1.0d3),HELP3,HELP1,PERCENT,MESSAGE
          else
            k = j - 540

            if ( DOPRINTEXTRA(k) == 1 )     write (62,103) NAMEXTRA5(k),
     &                    mod(VEXTRA(k),1.0d3),CONFERR,mod(HELP2,1.0d3),
     &                                               HELP3,HELP1,MESSAGE
            if ( DOPRINTEXTRA(k) == 2 )     write (62,102) NAMEXTRA5(k),
     &                    mod(VEXTRA(k),1.0d3),CONFERR,mod(HELP2,1.0d3),
     &                                       HELP3,HELP1,PERCENT,MESSAGE
          end if
        end if

      end do

!       if ( NAME(38) == "r1+r2" ) V(2) = abs(V(2))

100   FORMAT (I3,1X,A5,1X,F14.10," +/-",F14.10,2X,
     &               F14.10," +",F14.10," -",F14.10," (",F5.2,"%)  ",A5)
101   FORMAT (I3,1X,A5,1X,F14.10," +/-",F14.10,2X,
     &                    F14.10," +",F14.10," -",F14.10,11X,A5)
102   FORMAT    (4X,A5,1X,F14.10," +/-",F14.10,2X,
     &               F14.10," +",F14.10," -",F14.10," (",F5.2,"%)  ",A5)
103   FORMAT    (4X,A5,1X,F14.10," +/-",F14.10,2X,
     &                    F14.10," +",F14.10," -",F14.10,11X,A5)

      END SUBROUTINE OUTPUTSIM
!=======================================================================
      SUBROUTINE GETPRINTSIM (V,VARY,DOPRINT,DOPRINTEXTRA)
            ! This subroutine determines which fitted variable and depe-
            ! ndent variables to print out  for each simulation, for the
            ! tasks 5, 7, 8 and 9. Itis called from subroutine TASK5789.
      implicit none
      real*8 V(540)                       ! IN:  Fitted parameters
      integer VARY(540)                   ! IN:  Par adjustment integers
      integer DOPRINT(540)                ! OUT: 0 or 1 or 2 (see below)
      integer DOPRINTEXTRA(18)            ! OUT: 0 or 1 or 2 (see below)
      integer i                           ! LOCAL: loop counter

            ! 0: do not print out this parameter
            ! 1: print out this parameter
            ! 2: print out and give percentage error in OUTPUTSIM

      do i = 1,540
        DOPRINT(i) = 0
      end do
      do i = 1,18
        DOPRINTEXTRA(i) = 0
      end do

      do i = 1,540
        if ( VARY(i) /= 0 ) then
          DOPRINT(i) = 1
          if ( i==2 .or. i==3 .or. i==27 .or. i==28 ) DOPRINT(i) = 2
        end if
      end do

      if ( VARY(2) /= 0 .or. VARY(3) /= 0 ) then
        DOPRINTEXTRA(1) = 2                                ! rA+rB or rA
        DOPRINTEXTRA(2) = 2                                    ! k or rB
      end if

      if ( VARY(1)/=0 .or. VARY(2)/=0 .or. VARY(3)/=0 .or. VARY(5)/=0
     &     .or. VARY(6)/=0 .or. VARY(21)/=0 .or. VARY(22)/=0 ) then
        if ( V(1) > 0.0d0 ) DOPRINTEXTRA(3) = 1                  ! lB/lA
      end if

      if ( VARY(7) /= 0 .or. VARY(8) /= 0 ) then
        DOPRINTEXTRA(4) = 1                                 ! ecosw or e
        DOPRINTEXTRA(5) = 1                                 ! esinw or w
      end if

      if ( VARY(2) /= 0 .or. VARY(3) /= 0 .or. VARY(6) /= 0
     &                  .or. VARY(7) /= 0 .or. VARY(8) /= 0 ) then
        DOPRINTEXTRA(6) = 1                                     ! b(pri)
        DOPRINTEXTRA(7) = 1                                     ! b(sec)
      end if

      DOPRINTEXTRA(8) = 1

      if ( V(27) > -1.0d9 .and. V(28) > -1.0d9 ) then
        DOPRINTEXTRA(9) = 1                             ! semimajor axis
        DOPRINTEXTRA(10) = 1                                ! mass ratio
        do i = 11,14
          DOPRINTEXTRA(i) = 2                           ! M1, M2, R1, R2
        end do
        do i = 15,18
          DOPRINTEXTRA(i) = 1                         ! logg and density
        end do
      end if

      END SUBROUTINE GETPRINTSIM
!=======================================================================
      SUBROUTINE OUTPUT (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,DATAFORM,
     &                NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,VERR,WHAT,NSINE,
     &                 PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

            ! This subroutine writes information to the output files.
            ! Precisely what to write is given by parameter WHAT.
            ! If WHAT=0 (all TASKs) output initial pars and useful stuff
            ! If WHAT=1 (TASKs 3 and 4), output the final fitted params,
            !   useful quantities, residual curve, best-fit model curve.
            ! If WHAT=2 (TASKs 4, 5, 6, 7, 8, 9), output the best fitted
            !   parameters, useful quantities and various messages.
            ! If WHAT=3 (TASK 5), output best fit pars and useful stuff.
            ! If WHAT=101-122,  then this has been invoked from TASK6 to
            !   output the name of the parameter=(WHAT-100).
            ! If WHAT=201-222,  then this has come from TASK6  to output
            !   fit results for a certain value of parameter=(WHAT-200).

      implicit none
      integer MAXDATA                     ! IN: max number of datapoints
      character DATAFORM*1                ! IN: MAXDATA character length
      integer WHAT                        ! IN: Indicates what to output
      real*8 DATX(MAXDATA)                ! IN: Data independent varible
      real*8 DATY(MAXDATA)                ! IN: Data dependent variables
      real*8 DATERR(MAXDATA)              ! IN: Data errorbars
      integer DTYPE(MAXDATA)              ! IN: Type of  each data point
      real*8 V(540),VERR(540)             ! IN: Parameters  and  errors
      integer VARY(540)                   ! IN: Par adjustment  integers
      integer NDATA,NLR,NMIN              ! IN: Numbers  of  data points
      integer NSINE,NL3                   ! IN: Numbers of sines and L3s
      integer NECW,NESW                   ! IN: Numbers of ecosw & esinw
      integer PSINE(20)                   ! IN: Which par  for each sine
      integer NPOLY,PPOLY(50)             ! IN: Similar  for polynomials
      integer NUMINT                      ! IN: Number of numerical ints
      real*8  NINTERVAL                   ! IN: Time interval for numint
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      integer ITER                        ! IN: number of fit iterations
      real*8 CHISQ                        ! IN:   chi-squared of the fit
      integer NVARY,VWHERE(540)           ! LOCAL: which parameters vary
      character*5 NAME5(540)              ! GETIN: short parameter names
      character*5 NAMEXTRA5(18)           ! GETIN: short extra par names
      character*19 NAME19(540)            ! GETIN: long  parameter names
      character*19 NAMEXTRA19(18)         ! GETIN: long  extra par names
      real*8 ECC,OMEGA,ECOSW,ESINW        ! LOCAL: orbital    parameters
      real*8 R1,R2,R12,RRAT               ! LOCAL: fractional star radii
      real*8 RESIDSQSUM,SE                ! LOCAL: sum resids^2, std.err
      real*8 A,B,EPSLN1,EPSLN2,EPSLN      ! LOCAL: stellar shape  params
      real*8 MAG,LP,LS,HJD                ! LOCAL: LIGHT subroutine pars
      real*8 PHASE,HELP                   ! LOCAL: some useful variables
      real*8 OMC(MAXDATA)                 ! LOCAL: observed minus calc'd
      real*8 LIMBRIGHT                    ! LOCAL: total limb brightness
      real*8 UNINTMAG                     ! LOCAL: unintegrated magntude
      integer NPHASE                      ! LOCAL: number phases to plot
      character MESSAGE2*19               ! LOCAL: message to go to file
      real*8 ECQPHASES(6)                 ! LOCAL: important orbitphases
      integer i,j,k,ERROR,STATUS          ! LOCAL: counters & errorflags
      real*8 GETPHASE                     ! FUNCTION: calc orbital phase
      real*8 GETMODEL                     ! FUNCTION: calcs model output
      real*8 CHISQLC,CHISQLR              ! LOCAL: chisq of LCs, and LRs
      real*8 CHISQMIN,CHISQL3             ! LOCAL: chisq of tmin and L3s
      real*8 CHISQECW,CHISQESW            ! LOCAL: chisq of ecosw, esinw
      real*8 CHISQRV1,CHISQRV2            ! LOCAL: chisq of the  RV sets
      real*8 RES2SUMLC,RES2SUMLR          ! LOCAL: sum-squared of resids
      real*8 RES2SUMMIN,RES2SUML3         ! LOCAL: sum-squared of resids
      real*8 RES2SUMECW,RES2SUMESW        ! LOCAL: sum-squared of resids
      real*8 RES2SUMRV1,RES2SUMRV2        ! LOCAL: sum-squared of resids
      integer NLC,NRV1,NRV2               ! LOCAL: number of LCs and RVs
      real*8 RV1,RV2                      ! LOCAL: RVs of the two stars
      real*8 VEXTRA(18)                   ! GETIN: additional quantities
      integer NOPRINT                     ! LOCAL: flags what to output
      real*8 DEG2RAD,PI                   ! LOCAL: useful constants
      real*8 SINI,COS2I                   ! LOCAL: sin and cos-sq of inc
      real*8 ETERM,ETERM2,BTERM           ! LOCAL: equation calculations
      real*8 PRIDUR,SECDUR                ! LOCAL: durations of eclipses
      real*8 OMC_HI1,OMC_HI2,OMC_LO1      ! LOCAL: OMC extremum helper..
      real*8 OMC_LO2,SIG_HI1,SIG_HI2      ! LOCAL:   .. finder variables
      real*8 SIG_LO1,SIG_LO2,SIGMA        ! LOCAL: and sigma value.
      integer iOMC_HI1,iOMC_HI2,iOMC_LO1  ! LOCAL: OMC extremum helper..
      integer iOMC_LO2,iSIG_HI1,iSIG_HI2  ! LOCAL:  .. finder array ..
      integer iSIG_LO1,iSIG_LO2           ! LOCAL:  .. indices.

      DEG2RAD = 45.0d0 / atan(1.0d0)
      PI = atan(1.0d0) * 4.0d0

      ERROR = 0
      STATUS = 0
      NLC = 0
      NRV1 = 0
      NRV2 = 0
      do i = 1,NDATA
        if ( DTYPE(i) == 1 ) NLC = NLC + 1
        if ( DTYPE(i) == 7 ) NRV1 = NRV1 + 1
        if ( DTYPE(i) == 8 ) NRV2 = NRV2 + 1
      end do

      CALL GETNAME(V,NAME5,NAMEXTRA5,NAME19,NAMEXTRA19)

      CALL GETEXTRAS (V,VARY,NDATA,LDTYPE,CHISQ,NSINE,PSINE,NPOLY,PPOLY,
     &                                          NUMINT,NINTERVAL,VEXTRA)

            ! Find the stellar radii for both options (rA,rB), (rA+rB,k)

      if ( V(2) < 0.0d0 ) then
        R1 = abs(V(2))
        R2 = V(3)
        R12 = abs(V(2)) + V(3)
        RRAT = V(3) / abs(V(2))
      else
        R1 = abs(V(2)) / (1.0d0 + V(3))
        R2 = abs(V(2)) / (1.0d0 + (1.0d0/V(3)))
        R12 = abs(V(2))
        RRAT = V(3)
      end if

            ! Find the light contributions and reflection coefficients
            ! for the two stars. Do this at phase of first quadrature.

      call ECQUADPHASES(V(7),V(8),V(16),ECQPHASES)
      MAG = GETMODEL (V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                V(20)+V(19)*ECQPHASES(2),1,LP,LS,NUMINT,NINTERVAL)
!       LTOTAL = LP + LS + V(15)
      if ( VARY(11) == -1 )  V(11) = 0.4d0 * (LS/(1.0d0-V(15))) * R1**2
      if ( VARY(12) == -1 )  V(12) = 0.4d0 * (LP/(1.0d0-V(15))) * R2**2

            ! Find the number of variable parameters

      NVARY = 0
      do i = 1,540
        if ( VARY(i) == 1 .or. VARY(i) == 3 ) NVARY = NVARY + 1
      end do

            ! Output the values of the parameters and various other info

      if ( WHAT == 0 ) then
        write(62,107) "Initial values of the parameters:               "
      else if ( WHAT == 1 .or. WHAT == 2 ) then
        write(62,*)   " "
        write(62,107) "---------------------------------------         "
        write(62,110) ITER," iterations of EBOP completed     "
        write(62,107) "---------------------------------------         "
        write(62,*)   " "
        write(62,105)   "Warning: the uncertainties given below a",
     &                  "re only formal errors and are probably a"
        write(62,105)   "bit optimistic.  Use TASKs 6, 7 and 8 to",
     &                  " get reliable errors for the parameters."
        write(62,105)   "They also assume the noise is only Poiss",
     &                  "onian: use TASK9 if red noise is strong."
        if ( DATERR(1) < 0.0d0 ) then
          write(62,105) "No observational errors; so parameter er",
     &                  "rors have been scaled by sqrt(chi^2_red)"
        end if
        write(62,*)   " "
        if ( V(6) >= 89.9d0 .and. VARY(6) == 1 ) then
          write(62,105) "## Warning: the light curve solution has",
     &                  " an inclination of 89.9 degrees or more."
          write(62,105) "The minimisation algorithm cannot go bey",
     &                  "ond 90 degrees, so the solution has been"
          write(62,105) "kept artificially slightly below this va",
     &                  "lue to preserve its numerical stability."
          write (62,*) " "
        end if
        write(62,107) "Final values of the parameters:                 "
      end if

!-----------------------------------------------------------------------

      if ( WHAT == 0 .or. WHAT == 1 .or. WHAT == 2 .or. WHAT == 3 ) then
        if ( V(6) > 90.0d0 )  V(6) = 180.0d0 - V(6)
        if ( V(6) < 0.0d0 ) V(6) = V(6) + 180.0d0
        if ( LDTYPE(2) == -1 ) then
          V(5) = V(4)
          VERR(5) = VERR(4)
        end if

        do i = 1,540
          MESSAGE2 = "                   "
          if ( WHAT == 0 .and. (VARY(i) == 1 .or. VARY(i) == 3) )
     &                                  MESSAGE2 = "   (adjusted)      "
          if ( WHAT == 1 .and. (VARY(i) == 0 .or. VARY(i) == 2) )
     &                                  MESSAGE2 = "   (fixed)         "
          if ( i==11 .and. VARY(i)==-1) MESSAGE2 = "   (from geometry) "
          if ( i==12 .and. VARY(i)==-1) MESSAGE2 = "   (from geometry) "
          if ( i==30 .and. VARY(i)==-1) MESSAGE2 = "   (same as star A)"

            ! Now check for specific parameters which do not need to be
            ! outputted. Check every parameter included in the fit.

            ! V(14) is tidal lead-lag angle - JKTEBOP ignores this.
            ! V(27) and V(29) are velocity amplitudes
            ! V(21) to V(26) are nonlinear LD coefficients
            ! V(27) to V(30) are velocity parameters (K1,K2,Vsys1,Vsys2)
            ! V(31-90) are sine curve parameters (T0, period, amplitude)
            ! V(91-540) are poly (pivot,Tstart,Tend,const,x,x2,x3,x4,x5)

          if ( i <= 30+3*NSINE .or. (i>=91 .and. i<=90+9*NPOLY)) then
            NOPRINT = 0

            if ( i == 14 ) NOPRINT = 1
            if ( (i == 27 .or. i == 29) .and. V(27) < -1.d9) NOPRINT = 1
            if ( (i == 28 .or. i == 30) .and. V(28) < -1.d9) NOPRINT = 1
            if ( LDTYPE(1) == 1 .and. i >= 21 .and. i <= 23) NOPRINT = 1
            if ( LDTYPE(1) /= 9 .and. i >= 22 .and. i <= 23) NOPRINT = 1
            if ( LDTYPE(2) == 1 .and. i >= 24 .and. i <= 26) NOPRINT = 1
            if ( LDTYPE(2) >= 2 .and. LDTYPE(2) <= 6 .and. i >= 25 .and.
     &                                             i <= 26 ) NOPRINT = 1
            if ( LDTYPE(2) >= -6 .and. LDTYPE(2) <= -2 .and. i >= 25
     &                                       .and. i <= 26 ) NOPRINT = 1
            if ( LDTYPE(2) == -1 ) then
              if ( LDTYPE(1)==1 .and. i >= 24 .and. i <= 26) NOPRINT = 1
              if ( LDTYPE(1)/=9 .and. i >= 25 .and. i <= 26) NOPRINT = 1
            end if
            do j = 1,NPOLY
              k = 91 + 9*(j-1)
              if ( i == k+1 .and. V(k+1) <= -9.9d9 ) NOPRINT = 1
              if ( i == k+2 .and. V(k+2) >= 9.9d9 ) NOPRINT = 1
            end do

                    ! Now we have winnowed the unwanted output, print out
                    ! the results for each of the V parameters.

            if ( NOPRINT == 0 ) then
              if ( WHAT == 0 ) then
                write(62,113) i,NAME19(i),V(i),MESSAGE2
              else
                if ( VARY(i) /= 1 .and. VARY(i) /= 3 ) then
                  write(62,113) i,NAME19(i),V(i),MESSAGE2
                else
                  write(62,114) i,NAME19(i),V(i),VERR(i)
                end if
              end if
            end if
          end if

        end do

116   FORMAT ("Limb darkening law for primary star:       ",A19)
117   FORMAT ("Limb darkening law for secondary star:     ",A19)

        write (62,*) " "
        if ( LDTYPE(1) ==  1 ) write (62,116) "linear             "
        if ( LDTYPE(1) ==  2 ) write (62,116) "quadratic          "
        if ( LDTYPE(1) ==  3 ) write (62,116) "logarithmic        "
        if ( LDTYPE(1) ==  4 ) write (62,116) "square-root        "
        if ( LDTYPE(1) ==  5 ) write (62,116) "cubic              "
        if ( LDTYPE(1) ==  6 ) write (62,116) "power-2            "
        if ( LDTYPE(1) ==  9 ) write (62,116) "four-parameter     "
        if ( LDTYPE(1) == -2 ) write (62,116) "quadratic re-par   "
        if ( LDTYPE(1) == -3 ) write (62,116) "logarithmic re-par "
        if ( LDTYPE(1) == -4 ) write (62,116) "square-root re-par "
        if ( LDTYPE(1) == -5 ) write (62,116) "cubic re-parameter "
        if ( LDTYPE(1) == -6 ) write (62,116) "power-2 (h1 and h2)"

        if ( LDTYPE(2) ==  1 ) write (62,116) "linear             "
        if ( LDTYPE(2) ==  2 ) write (62,116) "quadratic          "
        if ( LDTYPE(2) ==  3 ) write (62,116) "logarithmic        "
        if ( LDTYPE(2) ==  4 ) write (62,116) "square-root        "
        if ( LDTYPE(2) ==  5 ) write (62,116) "cubic              "
        if ( LDTYPE(2) ==  6 ) write (62,116) "power-2            "
        if ( LDTYPE(2) ==  9 ) write (62,116) "four-parameter     "
        if ( LDTYPE(2) == -2 ) write (62,116) "quadratic re-par   "
        if ( LDTYPE(2) == -3 ) write (62,116) "logarithmic re-par "
        if ( LDTYPE(2) == -4 ) write (62,116) "square-root re-par "
        if ( LDTYPE(2) == -5 ) write (62,116) "cubic re-parameter "
        if ( LDTYPE(2) == -6 ) write (62,116) "power-2 (h1 and h2)"
        if ( LDTYPE(2) == -1 ) write (62,117) "same as primary    "

            ! Now check to see if the limb darkening coefficients are
            ! physically possible or reasonable, and print a warning
            ! if they are not. Remember that EBOP can happily function
            ! with physically impossible parameters (e.g. LDlin > 1.0)

        if ( LDTYPE(1) == 1 ) then
          LIMBRIGHT = 1.0d0 - V(4)
        else if ( LDTYPE(1) >= 2 .and. LDTYPE(1) <= 5 ) then
          LIMBRIGHT = 1.0d0 - V(4) - V(21)
        else if ( LDTYPE(1) >= -5 .and. LDTYPE(1) <= -2 ) then
          LIMBRIGHT = 0.5d0   !!! SORT THIS OUT LATER !!!
        else if ( LDTYPE(1) == 6 ) then
          LIMBRIGHT = V(4)
         else if ( LDTYPE(1) == -6 ) then
          LIMBRIGHT = 1.0d0 - V(4) + V(21)
        else if ( LDTYPE(1) == 9 ) then
          LIMBRIGHT = 1.0d0 - V(4) - V(21) - V(22) - V(23)
        else
          write (6,*) '### ERROR: LDTYPE(1) is wrong: ',LDTYPE(1)
          STOP
        end if
        if ( LIMBRIGHT > 1.0d0 .and. LDTYPE(1) /= 9 ) then
          write (62,*) " "
          write (62,105) "## Warning: the total limb darkening at ",
     &                   "the limb of star A is less than 0.0, so "
          write (62,105) "## the limb darkening coefficient(s) for",
     &                   " this star are physically unrealistic.  "
        else if ( LIMBRIGHT < 0.0d0 .and. LDTYPE(1) /= 9 ) then
          write (62,*) " "
          write (62,105) "## Warning: the total limb darkening at ",
     &                   "the limb of star A is greater than 1.0, "
          write (62,105) "## so the limb darkening coefficient(s) ",
     &                   "for this star are physically unlikely.  "
        end if

        if ( LDTYPE(2) == 1 ) then
          LIMBRIGHT = 1.0d0 - V(5)
        else if ( LDTYPE(2) >= 2 .and. LDTYPE(2) <= 5 ) then
          LIMBRIGHT = 1.0d0 - V(5) - V(24)
        else if ( LDTYPE(2) >= -5 .and. LDTYPE(2) <= -2 ) then
          LIMBRIGHT = 0.5d0   !!! SORT THIS OUT LATER !!!
        else if ( LDTYPE(2) == 6 ) then
          LIMBRIGHT = V(5)
         else if ( LDTYPE(2) == -6 ) then
          LIMBRIGHT = 1.0d0 - V(5) + V(24)
        else if ( LDTYPE(2) == 9 ) then
          LIMBRIGHT = 1.0d0 - V(5) - V(24) - V(25) - V(26)
        else
          if ( LDTYPE(2) /= -1 ) then
            write (6,*) '### ERROR: LDTYPE(2) is wrong: ',LDTYPE(2)
            STOP
          end if
        end if
        if ( LIMBRIGHT < 0.0d0 .and. LDTYPE(2) /= 9  ) then
          write (62,105) "## Warning: the total limb darkening at ",
     &                   "the limb of star B is less than 0.0, so "
          write (62,105) "## the limb darkening coefficient(s) for",
     &                   " this star are physically unrealistic.  "
        else if ( LIMBRIGHT > 1.0d0 .and. LDTYPE(2) /= 9  ) then
          write (62,105) "## Warning: the total limb darkening at ",
     &                   "the limb of star B is more than 1.0, so "
          write (62,105) "## the limb darkening coefficient(s) for",
     &                   " this star are physically unlikely.     "
        end if

        write (62,*) " "

        if ( WHAT /= 0 ) then
        write(62,102)"Phase of primary eclipse:           ",ECQPHASES(1)
        write(62,102)"Phase of first quadrature:          ",ECQPHASES(2)
        write(62,102)"Phase of secondary eclipse:         ",ECQPHASES(3)
        write(62,102)"Phase of second quadrature:         ",ECQPHASES(4)
        write(62,102)"Phase of periastron:                ",ECQPHASES(5)
        write(62,102)"Phase of apastron:                  ",ECQPHASES(6)
        write (62,*) " "
        endif

            ! Print radii of the stars and check if they are too large
            ! for EBOP to provide a good approximation to their shapes.

        if ( V(2) >= 0.0d0 ) then
          write (62,102) "Fractional primary radius:          ",R1
          write (62,102) "Fractional secondary radius:        ",R2
        else
          write (62,102) "Sum of the fractional radii:        ",R12
          write (62,102) "Ratio of the radii:                 ",RRAT
        end if

        if ( R1 >= 0.3d0 .or. R2 >= 0.3d0 ) then
          write (62,*) " "
          write (62,105) "## Warning: average radius is greater th",
     &                   "an 0.3, so the radii may be wrong by 5%!"
          write (62,105) "## See North & Zahn (2004, New Astronomy",
     &                   " Review, 48, 741) for further details.  "
          write (62,105) "## Do not use EBOP: go for a more sophis",
     &                   "ticated model which uses Roche geometry."
          write (62,*) " "
        else if ( R2 >= 0.25d0 .or. R2 >= 0.25d0 )  then
          write (62,*) " "
          write (62,105) "## Warning: average radius is greater th",
     &                   "an 0.25 so the radii may be wrong by 1%."
          write (62,105) "## See North & Zahn (2004, New Astronomy",
     &                   " Review, 48, 741) for further details.  "
          write (62,*) " "
        end if

        write (62,'(A26,F7.4,A3,2X,F17.10)')
     &         "Stellar light ratio (phase",ECQPHASES(2),"): ",LS/LP
        write (62,102) "Primary contribut'n to system light:", LP
        write (62,102) "Secondary contrib'n to system light:", LS
        write (62,*) " "

        write(62,102) "Impact parameter (primary eclipse): ",VEXTRA(6)
        write(62,102) "Impact paramtr (secondary eclipse): ",VEXTRA(7)

            ! If the ratio of the radii is small then calculate the
            ! durations of primary and secondary eclipse using equation
            ! 15 from Kipping (2010MNRAS.407..301K), as expanded by
            ! Matson et al (2016AJ....151..139M), equations 25 and 26.

        CALL GETEOMEGA (V(7),V(8),ECC,OMEGA,ECOSW,ESINW)

        if ( RRAT <= 0.2 ) then
          SINI = sin( V(6) / DEG2RAD )
          COS2I = cos( V(6) / DEG2RAD )**2
          ETERM = (1.0d0 - ECC**2) / (1.0d0 + ESINW)
          ETERM2 = (1.0d0 - ECC**2)**1.5d0 / (1.0d0 + ESINW)**2.0d0
          BTERM = sqrt(R12**2 - ETERM**2 * COS2I) / ETERM / SINI
          PRIDUR = V(19) / PI * ETERM2 * asin(BTERM)
          ETERM = (1.0d0 - ECC**2) / (1.0d0 - ESINW)
          ETERM2 = (1.0d0 - ECC**2)**1.5d0 / (1.0d0 - ESINW)**2.0d0
          BTERM = sqrt(R12**2 - ETERM**2 * COS2I) / ETERM / SINI
          SECDUR = V(19) / PI * ETERM2 * asin(BTERM)
          write(62,102) "Primary eclipse duration:           ",PRIDUR
          write(62,102) "Secondary eclipse duration:         ",SECDUR
        endif

        if ( V(7) > 5.0d0 ) then
          write (62,102) "Eccentricity * cos(omega):          ",ECOSW
          write (62,102) "Eccentricity * sin(omega):          ",ESINW
        else if ( V(7) /= 0.0d0 .or. V(8) /= 0.0d0 ) then
          write (62,102) "Orbital eccentricity e:             ",ECC
          write (62,102) "Periastron longitude omega (degree):",OMEGA
        end if
        write (62,*) " "

        if ( VARY(11) /= -1 .or. VARY(12) /= -1 ) then
          HELP = 0.4d0 * LS * R1**2
          write (62,102) "Geometric reflection coeff (star A):",HELP
          HELP = 0.4d0 * LP * R2**2
          write (62,102) "Geometric reflection coeff (star B):",HELP
        end if

        CALL BIAX (R1,abs(V(13)),A,B,EPSLN1)
        CALL BIAX (R2,1.0d0/abs(V(13)),A,B,EPSLN2)
        if ( V(13) <= 0.0d0 ) EPSLN1 = 0.0d0
        if ( V(13) <= 0.0d0 ) EPSLN2 = 0.0d0

        write (62,102) "Oblateness of the primary star:     ",EPSLN1
        write (62,102) "Oblateness of the secondary star:   ",EPSLN2
        if ( (WHAT == 1 .or. WHAT == 2) .and. V(13) < -0.000001d0 ) then
          CALL BIAX(R1,abs(V(13)),A,B,EPSLN)
          write(62,102)"Expected oblateness of primary:     ",EPSLN
          CALL BIAX(R2,1.0d0/abs(V(13)),A,B,EPSLN)
          write(62,102)"Expected oblateness of secondary:   ",EPSLN
        end if

        if ( EPSLN1 > 0.04d0 .or. EPSLN2 > 0.04d0 )  then
          write (62,*) " "
          write (62,105) "## Warning: oblateness is above the reco",
     &                   "mmended maximum value for EBOP of 0.04. "
          write (62,105) "## See Popper & Etzel (1981, AJ, 86, 102",
     &                   ") for justification and further details."
        end if
        write (62,*) " "

        if ( WHAT >= 1 .and. NRV1 >= 1 .and. NRV2 >= 1 ) then
          write(62,102)"Orbital semimajor axis (Rsun):      ",VEXTRA( 9)
          write(62,102)"Orbital semimajor axis (AU):        ",
     &                                               VEXTRA( 9)/214.94d0
          write(62,102)"Mass ratio (star B / star A):       ",VEXTRA(10)
          write(62,102)"Mass of star A (Msun)               ",VEXTRA(11)
          write(62,102)"Mass of star B (Msun)               ",VEXTRA(12)
          write(62,102)"Radius of star A (Rsun)             ",VEXTRA(13)
          write(62,102)"Radius of star B (Rsun)             ",VEXTRA(14)
          write(62,102)"Log surface gravity of star A (cgs):",VEXTRA(15)
          write(62,102)"Log surface gravity of star B (cgs):",VEXTRA(16)
          write(62,102)"Density of star A (solar units):    ",VEXTRA(17)
          write(62,102)"Density of star B (solar units):    ",VEXTRA(18)
          write (62,*) " "
        end if

      end if

!-----------------------------------------------------------------------
            ! Output the observations, fitted model values and residuals
            ! for the case of WHAT=1.    For other values of WHAT, these
            ! quantities are not outputted but the CHISQ is still needed
            ! so this loop is still executed.

      if ( WHAT == 1 .or. WHAT == 2 .or. WHAT == 3 ) then
        if ( WHAT == 1 ) then
          write (63,'(A35,A35)') "#     TIME       MAGNITUDE    ERROR",
     &                           "      PHASE        MODEL      (O-C)"
          if ( NRV1 >= 1 )
     &    write (65,'(A35,A35)') "#     TIME       RV_STAR_A    ERROR",
     &                           "      PHASE        MODEL      (O-C)"
          if ( NRV2 >= 1 )
     &    write (66,'(A35,A35)') "#     TIME       RV_STAR_B    ERROR",
     &                           "      PHASE        MODEL      (O-C)"
        end if

        OMC_HI1 = 0.0d0     ! Highest O-C value for the lightcrve data
        OMC_HI2 = 0.0d0     ! Second-highest O-C value for the LC data
        OMC_LO1 = 0.0d0     ! Lowest O-C value for the lightcurve data
        OMC_LO2 = 0.0d0     ! Second-lowest O-C value  for the LC data
        SIG_HI1 = 0.0d0     ! Highest sigma value for the lightcv data
        SIG_HI2 = 0.0d0     ! Second-highest sigmavalue forthe LC data
        SIG_LO1 = 0.0d0     ! Lowest sigma value for the lightcrv data
        SIG_LO2 = 0.0d0     ! Second-lowest sigma value forthe LC data
        iOMC_HI1 = 0        ! Data array index value for OMC_HI1
        iOMC_HI2 = 0        ! Data array index value for OMC_HI2
        iOMC_LO1 = 0        ! Data array index value for OMC_LO1
        iOMC_LO2 = 0        ! Data array index value for OMC_LO2
        iSIG_HI1 = 0        ! Data array index value for SIG_HI1
        iSIG_HI2 = 0        ! Data array index value for SIG_HI2
        iSIG_LO1 = 0        ! Data array index value for SIG_LO1
        iSIG_LO2 = 0        ! Data array index value for SIG_LO2

        RES2SUMLC  = 0.0d0
        RES2SUMLR  = 0.0d0
        RES2SUMMIN = 0.0d0
        RES2SUML3  = 0.0d0
        RES2SUMECW = 0.0d0
        RES2SUMESW = 0.0d0
        RES2SUMRV1 = 0.0d0
        RES2SUMRV2 = 0.0d0
        CHISQLC  = 0.0d0
        CHISQLR  = 0.0d0
        CHISQMIN = 0.0d0
        CHISQL3  = 0.0d0
        CHISQECW = 0.0d0
        CHISQESW = 0.0d0
        CHISQRV1 = 0.0d0
        CHISQRV2 = 0.0d0

        do i = 1,NDATA
          MAG = GETMODEL(V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,DATX(i),
     &                                  DTYPE(i),LP,LS,NUMINT,NINTERVAL)
          OMC(i) = DATY(i)-MAG
          if ( WHAT == 1 ) then
            if ( DTYPE(i) == 1 ) write (63,104) DATX(i),DATY(i),
     &                DATERR(i),GETPHASE(DATX(i),V(19),V(20)),MAG,OMC(i)
            if ( DTYPE(i) == 7 ) write (65,1047) DATX(i),DATY(i),
     &                DATERR(i),GETPHASE(DATX(i),V(19),V(20)),MAG,OMC(i)
            if ( DTYPE(i) == 8 ) write (66,1047) DATX(i),DATY(i),
     &                DATERR(i),GETPHASE(DATX(i),V(19),V(20)),MAG,OMC(i)
          end if
          if ( DTYPE(i)==1 ) RES2SUMLC  = RES2SUMLC  + OMC(i)**2
          if ( DTYPE(i)==2 ) RES2SUMLR  = RES2SUMLR  + OMC(i)**2
          if ( DTYPE(i)==3 ) RES2SUMMIN = RES2SUMMIN + OMC(i)**2
          if ( DTYPE(i)==4 ) RES2SUML3  = RES2SUML3  + OMC(i)**2
          if ( DTYPE(i)==5 ) RES2SUMECW = RES2SUMECW + OMC(i)**2
          if ( DTYPE(i)==6 ) RES2SUMESW = RES2SUMESW + OMC(i)**2
          if ( DTYPE(i)==7 ) RES2SUMRV1 = RES2SUMRV1 + OMC(i)**2
          if ( DTYPE(i)==8 ) RES2SUMRV2 = RES2SUMRV2 + OMC(i)**2
          if ( DTYPE(i)==1 ) CHISQLC = CHISQLC + (OMC(i)/DATERR(i))**2
          if ( DTYPE(i)==2 ) CHISQLR = CHISQLR + (OMC(i)/DATERR(i))**2
          if ( DTYPE(i)==3 ) CHISQMIN= CHISQMIN+ (OMC(i)/DATERR(i))**2
          if ( DTYPE(i)==4 ) CHISQL3 = CHISQL3 + (OMC(i)/DATERR(i))**2
          if ( DTYPE(i)==5 ) CHISQECW= CHISQECW+ (OMC(i)/DATERR(i))**2
          if ( DTYPE(i)==6 ) CHISQESW= CHISQESW+ (OMC(i)/DATERR(i))**2
          if ( DTYPE(i)==7 ) CHISQRV1= CHISQRV1+ (OMC(i)/DATERR(i))**2
          if ( DTYPE(i)==8 ) CHISQRV2= CHISQRV2+ (OMC(i)/DATERR(i))**2

                  ! Now keep track of the two highest and two lowest
                  ! values of (O-C) in the light curve data. Retain
                  ! their values and also the datapoint numbers.

          if ( WHAT == 1 .or. WHAT == 2 .or. WHAT == 3 ) then
           if ( DTYPE(i) == 1 ) then
              if ( OMC(i) > OMC_HI1 ) then
                OMC_HI2 = OMC_HI1
                iOMC_HI2 = iOMC_HI1
                OMC_HI1 = OMC(i)
                iOMC_HI1 = i
              else if ( OMC(i) > OMC_HI2 ) then
                OMC_HI2 = OMC(i)
                iOMC_HI2 = i
              endif
              if ( OMC(i) < OMC_LO1 ) then
                OMC_LO2 = OMC_LO1
                iOMC_LO2 = iOMC_LO1
                OMC_LO1 = OMC(i)
                iOMC_LO1 = i
              else if ( OMC(i) < OMC_LO2 ) then
                OMC_LO2 = OMC(i)
                iOMC_LO2 = i
              endif
            endif

                  ! Now keep track of the two highest and two lowest
                  ! sigma values in the light curve data. Retain
                  ! their values and also the datapoint numbers.

            if ( DATERR(1) > 0.0d0 ) then
              if ( DTYPE(i) == 1 ) then
                SIGMA = OMC(i) / DATERR(i)
                if ( SIGMA > SIG_HI1 ) then
                  SIG_HI2 = SIG_HI1
                  iSIG_HI2 = iSIG_HI1
                  SIG_HI1 = SIGMA
                  iSIG_HI1 = i
                else if ( SIGMA > SIG_HI2 ) then
                  SIG_HI2 = SIGMA
                  iSIG_HI2 = i
                endif
                if ( SIGMA < SIG_LO1 ) then
                  SIG_LO2 = SIG_LO1
                  iSIG_LO2 = iSIG_LO1
                  SIG_LO1 = SIGMA
                  iSIG_LO1 = i
                else if ( SIGMA < SIG_LO2 ) then
                  SIG_LO2 = SIGMA
                  iSIG_LO2 = i
                endif
              endif
            endif
          endif

        end do

        CHISQ = CHISQLC + CHISQLR + CHISQMIN + CHISQL3 + CHISQECW
     &          + CHISQESW + CHISQRV1 + CHISQRV2

        write(62,1022)"Total number of datapoints:         ",NDATA
        write(62,1022)"Total number of fitted parameters:  ",NVARY
        write(62,1022)"Total number of degrees of freedom: ",NDATA-NVARY

        if ( DATERR(1) > 0.0d0 ) then
          write (62,102) "Total chisq of the fit:             ",CHISQ
          write (62,102) "Reduced chisq of the fit:           ",
     &                                           CHISQ/dble(NDATA-NVARY)
        end if

        write (62,1022) "Total number of LC datapoints:      ",NLC
        write (62,102) "rms of the LC residuals (mmag):     ",
     &                            sqrt(RES2SUMLC / dble(NLC)) * 1000.0d0
        if ( DATERR(1) > 0.0d0 ) write (62,102)
     &          "Reduced chisq for the LC data:      ",CHISQLC/dble(NLC)

        if ( NLR >= 1 ) then
          write (62,102) "rms of the light ratio residuals:   ",
     &                                       sqrt(RES2SUMLR / dble(NLR))
          write (62,102) "Reduced chisq for the light ratios: ",
     &                                                 CHISQLR/dble(NLR)
        end if

        if ( NL3 >= 1 ) then
          write (62,102) "rms of the third light residuals:   ",
     &                                       sqrt(RES2SUML3 / dble(NL3))
          write (62,102) "Reduced chisq for the third light:  ",
     &                                                 CHISQL3/dble(NL3)
        end if

        if ( NMIN >= 1 ) then
          write (62,102) "rms of the tmin residuals (d):      ",
     &                                     sqrt(RES2SUMMIN / dble(NMIN))
          write (62,102) "Reduced chisq for the tmin data:    ",
     &                                               CHISQMIN/dble(NMIN)
        end if

        if ( NECW >= 1 ) then
          if ( V(7) > 9.99 ) then
            write (62,102) "rms of the eccentricity residuals:  ",
     &                                     sqrt(RES2SUMECW / dble(NECW))
            write (62,102) "Reduced chisq for the eccentricity: ",
     &                                               CHISQECW/dble(NECW)
          else
            write (62,102) "rms of the ecosw residuals:         ",
     &                                     sqrt(RES2SUMECW / dble(NECW))
            write (62,102) "Reduced chisq for the ecosw values: ",
     &                                               CHISQECW/dble(NECW)
          end if
        end if

        if ( NESW >= 1 ) then
          if ( V(7) > 9.99 ) then
            write (62,102) "rms of the omega residuals:         ",
     &                                     sqrt(RES2SUMESW / dble(NESW))
            write (62,102) "Reduced chisq for the omega values: ",
     &                                               CHISQESW/dble(NESW)
          else
            write (62,102) "rms of the esinw residuals:         ",
     &                                     sqrt(RES2SUMESW / dble(NESW))
            write (62,102) "Reduced chisq for the esinw values: ",
     &                                               CHISQESW/dble(NESW)
          end if
        end if

        if ( NRV1 >= 1 ) then
          write (62,1022) "Total number of RVs for star A:     ",NRV1
          write (62,102) "rms of star A RV residuals (km/s):  ",
     &                                     sqrt(RES2SUMRV1 / dble(NRV1))
          write (62,102) "Reduced chisq for the RVs of star A:",
     &                                               CHISQRV1/dble(NRV1)
        end if

        if ( NRV2 >= 1 ) then
          write (62,1022) "Total number of RVs for star B:     ",NRV2
          write (62,102) "rms of star B RV residuals (km/s):  ",
     &                                     sqrt(RES2SUMRV2 / dble(NRV2))
          write (62,102) "Reduced chisq for the RVs of star B:",
     &                                               CHISQRV2/dble(NRV2)
        end if

        write (62,*) " "

1567  FORMAT        (I7,2X,F14.6,3X,F10.8,2X,F10.5,1X,F10.5,1X,F10.5)

        if ( WHAT == 1 .or. WHAT == 2 .or. WHAT == 3 ) then
          write (62,'(A36,A44)') "The two highest and two lowest O-C v",
     &                   "alues in the light curve:                   "
          write (62,'(A36,A44)') " Number        Time       Phase     ",
     &                   "    Magnitude    Model      (O-C)           "
          write(62,1567)iOMC_LO1,DATX(iOMC_LO1),GETPHASE(DATX(iOMC_LO1),
     &        V(19),V(20)),DATY(iOMC_LO1),DATY(iOMC_LO1)-OMC_LO1,OMC_LO1
          write(62,1567)iOMC_LO2,DATX(iOMC_LO2),GETPHASE(DATX(iOMC_LO2),
     &        V(19),V(20)),DATY(iOMC_LO2),DATY(iOMC_LO2)-OMC_LO2,OMC_LO2
          write(62,1567)iOMC_HI2,DATX(iOMC_HI2),GETPHASE(DATX(iOMC_HI2),
     &        V(19),V(20)),DATY(iOMC_HI2),DATY(iOMC_HI2)-OMC_HI2,OMC_HI2
          write(62,1567)iOMC_HI1,DATX(iOMC_HI1),GETPHASE(DATX(iOMC_HI1),
     &        V(19),V(20)),DATY(iOMC_HI1),DATY(iOMC_HI1)-OMC_HI1,OMC_HI1
          write (62,*) " "

        if ( DATERR(1) > 0.0d0 ) then
            write (62,'(A34,A46)') "The two highest and two lowest res",
     &                 "idual values in the light curve:              "
            write (62,'(A34,A46)') " Number        Time       Phase   ",
     &                 "      Magnitude    Model      sigma           "
          write(62,1567)iOMC_LO1,DATX(iOMC_LO1),GETPHASE(DATX(iOMC_LO1),
     &        V(19),V(20)),DATY(iOMC_LO1),DATY(iOMC_LO1)-OMC_LO1,OMC_LO1
          write(62,1567)iOMC_LO2,DATX(iOMC_LO2),GETPHASE(DATX(iOMC_LO2),
     &        V(19),V(20)),DATY(iOMC_LO2),DATY(iOMC_LO2)-OMC_LO2,OMC_LO2
          write(62,1567)iOMC_HI2,DATX(iOMC_HI2),GETPHASE(DATX(iOMC_HI2),
     &        V(19),V(20)),DATY(iOMC_HI2),DATY(iOMC_HI2)-OMC_HI2,OMC_HI2
          write(62,1567)iOMC_HI1,DATX(iOMC_HI1),GETPHASE(DATX(iOMC_HI1),
     &        V(19),V(20)),DATY(iOMC_HI1),DATY(iOMC_HI1)-OMC_HI1,OMC_HI1
            write (62,*) " "
          endif
        endif

      end if

!-----------------------------------------------------------------------
            ! Now output the best fit on a grid of orbital phases

      if ( WHAT == 1 ) then
        NPHASE = 10001
        if ( R1 < 0.01d0 .or. R2 < 0.01d0 ) NPHASE = 100001
        write (64,'(A38,A37)') "# PHASE  MAGNITUDE    L1       L2     ",
     &                          "  L3     UNINTMAG    RV_A      RV_B  "

        do i = 1,NPHASE
          PHASE = (i-1) / dble(NPHASE-1)
          HJD = V(20) + PHASE * V(19)
          MAG = GETMODEL(V,VARY,LDTYPE,0,PSINE,0,PPOLY,HJD,1,LP,LS,
     &                                                 NUMINT,NINTERVAL)
          UNINTMAG = GETMODEL(V,VARY,LDTYPE,0,PSINE,0,PPOLY,HJD,1,LP,LS,
     &                                                      1,NINTERVAL)
          if ( V(27) >= -1.0d9 ) then
            RV1 = GETMODEL(V,VARY,LDTYPE,0,PSINE,0,PPOLY,HJD,7,LP,LS,
     &                                                 NUMINT,NINTERVAL)
          else
            RV1 = 0.0d0
          end if
          if ( V(28) >= -1.0d9 ) then
          RV2 = GETMODEL(V,VARY,LDTYPE,0,PSINE,0,PPOLY,HJD,8,LP,LS,
     &                                                 NUMINT,NINTERVAL)
          else
            RV2 = 0.0d0
          end if

          write (64,12) PHASE,MAG,LP,LS,V(15),UNINTMAG,RV1,RV2
12        FORMAT (F7.5,1X,F9.6,1X,3(1X,F8.5),1X,F9.6,2(1X,F9.4))
        end do

      end if

!-----------------------------------------------------------------------
      ! Now output the results and residuals for the direct observationl
      ! constraints on the times of minimum, spectroscopic light ratios,
      ! third light, ecosw and esinw.

      if ( WHAT == 1 .or. WHAT == 2 .or. WHAT == 3  ) then
        if ( NLR>0 .or. NMIN>0 .or. NL3>0 .or. NECW>0 .or. NESW>0) then
          write (62,'(A36,A36)') "------------------------------------",
     &                           "-------------------------------     "
          write (62,'(A36,A36)') "Results for the additional observed ",
     &                           "quantities                          "
          write (62,'(A36,A36)') "------------------------------------",
     &                           "-------------------------------     "
          write (62,'(A36,A36)') "Type      Time/Cycle     Measurement",
     &                           "     Error      Model     Sigma     "
          do i = 1,NDATA
            if ( DTYPE(i) >= 2 .and. DTYPE(i) <= 6 ) then
              HELP = GETMODEL(V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                        DATX(i),DTYPE(i),LP,LS,NUMINT,NINTERVAL)
!               HELP2 = DATX(i)
            end if

            if ( DTYPE(i) == 2 ) write (62,120) "LRAT", DATX(i),
     &               DATY(i),DATERR(i),HELP,(HELP-DATY(i))/DATERR(i)
            if ( DTYPE(i) == 4 ) write (62,120) "THDL", DATX(i),
     &               DATY(i),DATERR(i),HELP,(HELP-DATY(i))/DATERR(i)
            if ( DTYPE(i) == 3 ) write (62,120) "TMIN", DATX(i),
     &               DATY(i),DATERR(i),HELP,(HELP-DATY(i))/DATERR(i)

            if ( DTYPE(i) == 5 ) then
              if ( V(7) > 5.0d0 ) write (62,120) "ECC ", DATX(i),
     &                            DATY(i),DATERR(i),HELP,
     &                                 (HELP-DATY(i))/DATERR(i)
              if ( V(7) < 5.0d0 ) write (62,120) "ECSW", DATX(i),
     &               DATY(i),DATERR(i),HELP,(HELP-DATY(i))/DATERR(i)
            end if

            if ( DTYPE(i) == 6 ) then
              if ( V(7) > 5.0d0 ) write (62,120) "OMGA", DATX(i),
     &               DATY(i),DATERR(i),HELP,(HELP-DATY(i))/DATERR(i)
              if ( V(7) < 5.0d0 ) write (62,120) "ESNW", DATX(i),
     &               DATY(i),DATERR(i),HELP,(HELP-DATY(i))/DATERR(i)
            end if

          end do
          write (62,'(A36,A36)') "------------------------------------",
     &                           "-------------------------------     "
          write (62,*) " "
        end if

      end if

!--------------------------TASK6----------------------------------------
            ! If WHAT is between 100 and 122 then output the name of the
            ! parameter which is to be investigated next.
            ! If WHAT is between 200 and 222 then output various results
            ! from the last EBOP fit on one line for comparison.

      if ( WHAT > 100 .and. WHAT < 187 ) then
        j = 0
        do i = 1,540
          if ( VARY(i) == 1 ) then
            j = j + 1
            VWHERE(j) = i
          end if
        end do
        write (62,*) " "
        write (62,'(A41,A19)')
     &      "Now investigating the fitting parameter: ",NAME19(WHAT-100)
        write (62,'(A4,3X,37(A5,7X))')   "ITER",NAME5(WHAT-100)," rms ",
     &                  (NAME5(VWHERE(k)),k=1,j)," r_1 "," r_2 ","L2/L1"
        write(6,'(A44,A19)')
     &  ">> Now investigating the fitting parameter:  ",NAME19(WHAT-100)
      end if

!-----------------------------------------------------------------------

      if ( WHAT > 200 .and. WHAT < 287 ) then
        RESIDSQSUM = 0.0d0
        do i = 1,NDATA
          MAG = GETMODEL(V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,DATX(i),
     &                                         1,LP,LS,NUMINT,NINTERVAL)
          RESIDSQSUM = RESIDSQSUM + (DATY(i)-MAG)**2
        end do
        SE = sqrt(RESIDSQSUM / NDATA) * 1000.0d0
        j = 0
        do i = 1,540
          if ( VARY(i) == 1 ) then
            j = j + 1
            VWHERE(j) = i
          end if
        end do
       call ECQUADPHASES(V(7),V(8),V(16),ECQPHASES)
       write(62,'(I3,23(1X,F10.7))')   ITER, mod(V(WHAT-200),1.0d3), SE,
     &                          (mod(V(VWHERE(k)),1.0d3),k=1,j), R1, R2,
     &                   GETMODEL(V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                V(20)+V(19)*ECQPHASES(2),2,LP,LS,NUMINT,NINTERVAL)
      end if
!-----------------------------------------------------------------------

! 100   FORMAT (I2,2X,A19,2X,F18.10)
! 101   FORMAT (I2,2X,A19,2X,F18.10,A12)
102   FORMAT (A36,2X,F17.10)
! 1021  FORMAT (A36,2X,E17.10)
1022  FORMAT (A36,2X,I17)
104   FORMAT (F14.6,1X,F11.6,1X,F10.6,3X,F10.8,1X,F10.6,1X,F10.6)
1047  FORMAT (F14.6,1X,F11.5,1X,F10.6,3X,F10.8,1X,F10.5,1X,F10.5)
105   FORMAT (A40,A40)
107   FORMAT (A48)
! 108   FORMAT (F14.6,F10.6)
110   FORMAT (I6,A34)
! 111   FORMAT (I2,2X,A18,2X,F18.10,A17)
! 112   FORMAT (A39,F7.4)
113   FORMAT (I3,2X,A19,2X,F18.10,A19)
114   FORMAT (I3,2X,A19,2X,F18.10," +/- ",F14.10)
! 115   FORMAT (A36,3X,I10,2X,I4)
! 118   FORMAT (F8.1,1X,F13.5,1X,F8.5,2X,F13.5,1X,F6.2)
! 119   FORMAT (F13.5,1X,F7.4,1X,F6.4,2X,F7.4,F7.2)
120   FORMAT (A4,3X,F14.6,1X,F14.6,1X,F10.6,1X,F12.6,F7.2)

      END SUBROUTINE OUTPUT
!=======================================================================
      SUBROUTINE GET_DV (V,DV,NPOLY,PPOLY)
            ! This subroutine produces the DV values: the intervals for
            ! evaluating the numerical derivatives for each variable.
            ! The solutions are calculated for V+DV and V-DV, so the DV
            ! values are half the size of the actual numerical interval.
      implicit none
      real*8 V(540)                       ! IN: photometric parameters
      integer NPOLY                       ! IN: number of polynomials
      integer PPOLY(50)                   ! IN: parameters of the polys
      real*8 DV(540)                      ! OUT: DV values
      real*8 DVD(540)                     ! LOCAL: basic DV values
      real*8 VAL                          ! LOCAL: value of parameter
      integer i,j                         ! LOCAL: loop counters

      data DVD/ 0.01d0,   0.01d0,   0.01d0,   -0.01d0, ! J rs k LDA1
     &         -0.01d0,  -0.01d0, -0.001d0,  -0.001d0, ! LDB1 i e w
     &         -0.02d0,  -0.02d0,  0.001d0,   0.001d0, ! GD1,2 ref1,2
     &          0.01d0,  -1.0d0,   -0.01d0, -0.0001d0, ! q tangl L3 dPhi
     &        -0.001d0,   0.5d0,    0.1d-6, -0.0001d0, ! sfact ring P T0
     &         -0.01d0,  -0.01d0,  -0.01d0,            ! LDA2 LDA3 LDA4
     &         -0.01d0,  -0.01d0,  -0.01d0,            ! LDB2 LDB3 LDB4
     &          0.01d0,   0.01d0,   -0.1d0,    -0.1d0, ! KA KB Vsys Vsys
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 1 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 2 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 3 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 4 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 5 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 6 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 7 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 8 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 9 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine10 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine11 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine12 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine13 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine14 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine15 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine16 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine17 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine18 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine19 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine20 T0,p,amp
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 1
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 2
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 3
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 4
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 5
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 6
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 7
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 8
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 9
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 10
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 11
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 12
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 13
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 14
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 15
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 16
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 17
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 18
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 19
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 20
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 21
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 22
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 23
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 24
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 25
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 26
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 27
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 28
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 29
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 30
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 31
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 32
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 33
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 34
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 35
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 36
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 37
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 38
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 39
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 40
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 41
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 42
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 43
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 44
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 45
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 46
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 47
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 48
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  ! poly 49
     &        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/  ! poly 50

            ! The DVD values are the basic information from which the DV
            ! values are calculated, and essentially represent the order
            ! of magnitude expected for each parameter (based on experi-
            ! ence) and the way in which the DV values should be done.
            ! If DVD(i) > 0 then is it a multiplicative perturbation.
            ! If DVD(i) < 0 then is it an additive perturbation.

      do i = 1,540
        if ( DVD(i) >= 0.0d0 ) then
          DV(i) = abs(DVD(i)*V(i))
        else
          DV(i) = abs(DVD(i))
        end if

            ! The sine timesbases need special treatment: set them to
            ! one hundredth of the period of the sine wave in question.

        if ( i==31 .or. i==34 .or. i==37 .or. i==40 .or. i==43 .or.
     &       i==46 .or. i==49 .or. i==52 .or. i==55 .or. i==58 .or.
     &       i==61 .or. i==64 .or. i==67 .or. i==70 .or. i==73 .or.
     &       i==76 .or. i==79 .or. i==82 .or. i==85 .or. i==88 ) then
          DV(i) = V(i+1) / 100.0d0
        end if
      end do

            ! Ths polynomials need special treatment because the type
            ! and size of derivative wanted depends on the parameter
            ! the poynomial is applied to.

!               if ( WHICHPAR == "J"  ) PPOLY(NPOLY) = 1
!               if ( WHICHPAR == "r1" ) PPOLY(NPOLY) = 2
!               if ( WHICHPAR == "r2" ) PPOLY(NPOLY) = 3
!               if ( WHICHPAR == "i"  ) PPOLY(NPOLY) = 6
!               if ( WHICHPAR == "e"  ) PPOLY(NPOLY) = 7
!               if ( WHICHPAR == "ec" ) PPOLY(NPOLY) = 7
!               if ( WHICHPAR == "w"  ) PPOLY(NPOLY) = 8
!               if ( WHICHPAR == "es" ) PPOLY(NPOLY) = 8
!               if ( WHICHPAR == "L3" ) PPOLY(NPOLY) = 15
!               if ( WHICHPAR == "sf" ) PPOLY(NPOLY) = 17
!               if ( WHICHPAR == "L1" ) PPOLY(NPOLY) = -1
!               if ( WHICHPAR == "L2" ) PPOLY(NPOLY) = -2

      do i = 1,NPOLY
        j = 91 + 9*(i-1)
        DV(j)   = 0.0001d0         ! Timebase - will not be adjusted
        DV(j+1) = 0.0001d0         ! Timebase - will not be adjusted
        DV(j+2) = 0.0001d0         ! Timebase - will not be adjusted

            ! First fix the additive numerical intervals

        if ( PPOLY(i) ==  1 .or. PPOLY(i) ==  6 .or. PPOLY(i) ==  7 .or.
     &       PPOLY(i) ==  7 .or. PPOLY(i) == 15 .or. PPOLY(i) == 17 .or.
     &       PPOLY(i) == 20 .or. PPOLY(i) == 27 .or. PPOLY(i) == 28 .or.
     &       PPOLY(i) == 29 .or. PPOLY(i) == 30 .or. PPOLY(i) == -1 .or.
     &       PPOLY(i) == -2 ) then
          DV(j+3) = 0.001d0
          DV(j+4) = 0.001d0
          DV(j+5) = 0.0001d0
          DV(j+6) = 0.0001d0
          DV(j+7) = 0.00001d0
          DV(j+8) = 0.00001d0
        end if

            ! omega needs special treatment because apsidal motion (the
            ! rate of change of omega) is very slow.

        if ( PPOLY(i) == 8 ) then
          DV(j+3) = 0.00001d0
          DV(j+4) = 0.00001d0
          DV(j+5) = 0.00001d0
          DV(j+6) = 0.00001d0
          DV(j+7) = 0.00001d0
          DV(j+8) = 0.00001d0
        end if

            ! Finally fix the multiplicative numerical intervals

        if ( PPOLY(i) == 2 .or. PPOLY(i) == 3 ) then
          if ( PPOLY(i) == 2 ) VAL = V(2)
          if ( PPOLY(i) == 3 ) VAL = V(3)
          DV(j+3) = 0.001d0 * VAL
          DV(j+4) = 0.001d0 * VAL
          DV(j+5) = 0.001d0 * VAL
          DV(j+6) = 0.0001d0 * VAL
          DV(j+7) = 0.0001d0 * VAL
          DV(j+8) = 0.0001d0 * VAL
        end if
      end do

      END SUBROUTINE GET_DV
!=======================================================================
      SUBROUTINE GETNAME (V,NAME5,NAMEXTRA5,NAME19,NAMEXTRA19)
      implicit none
      real*8 V(540)                              ! IN: fitted parameters
      character*5 NAME5(540),NAMEXTRA5(18)       ! OUT: short names
      character*19 NAME19(540),NAMEXTRA19(18)    ! OUT: long names
      character*5 N5(540),NE5(18)                ! LOCAL: short names
      character*19 N19(540),NE19(18)             ! LOCAL: long names
      integer i                                  ! LOCAL: loop counter

            ! I appear to need to define local variables (N5, NE5, N19,
            ! NE19) and then put them in the output arrays (NAME5, NAME-
            ! XTRA5, NAME19, NAMEXTRA19). ifort doesn't allow the output
            ! arrays to be filled using a data statement. What a pain.

      data N5 /         "    J","rA+rB","    k","LD_A1","LD_B1","  inc",
     &                  "ecosw","esinw"," GD_A"," GD_B","reflA","reflB",
     &                  "qphot","T_lag","  L_3","phcor","sfact","iring",
     &                  "P_orb","  T_0","LD_A2","LD_A3","LD_A4","LD_B2",
     &                  "LD_B3","LD_B4","  K_A","  K_B","VsysA","VsysB",
     &          "sin1T","sin1P","sin1A", "sin2T","sin2P","sin2A",
     &          "sin3T","sin3P","sin3A", "sin4T","sin4P","sin4A",
     &          "sin5T","sin5P","sin5A", "sin6T","sin6P","sin6A",
     &          "sin7T","sin7P","sin7A", "sin8T","sin8P","sin8A",
     &          "sin9T","sin9P","sin9A", "sn10T","sn10P","sn10A",
     &          "sn11T","sn11P","sn11A", "sn12T","sn12P","sn12A",
     &          "sn13T","sn13P","sn13A", "sn14T","sn14P","sn14A",
     &          "sn15T","sn15P","sn15A", "sn16T","sn16P","sn16A",
     &          "sn17T","sn17P","sn17A", "sn18T","sn18P","sn18A",
     &          "sn19T","sn19P","sn19A", "sn20T","sn20P","sn20A",
     & "p01pv","p01ts","p01tf","p01co","p01x ","p01x2","p01x3","p01x4",
     & "p01x5","p02pv","p02ts","p02tf","p02co","p02x ","p02x2","p02x3",
     & "p02x4","p02x5","p03pv","p03ts","p03tf","p03co","p03x ","p03x2",
     & "p03x3","p03x4","p03x5","p04pv","p04ts","p04tf","p04co","p04x ",
     & "p04x2","p04x3","p04x4","p04x5","p05pv","p05ts","p05tf","p05co",
     & "p05x ","p05x2","p05x3","p05x4","p05x5","p05pv","p05ts","p05tf",
     & "p06co","p06x ","p06x2","p06x3","p06x4","p06x5","p06pv","p06ts",
     & "p07tf","p07co","p07x ","p07x2","p07x3","p07x4","p07x5","p07pv",
     & "p08ts","p08tf","p08co","p08x ","p08x2","p08x3","p08x4","p08x5",
     & "p09pv","p09ts","p09tf","p09co","p09x ","p09x2","p09x3","p09x4",
     & "p09x5","p10pv","p10ts","p10tf","p10co","p10x ","p10x2","p10x3",
     & "p10x4","p10x5","p11pv","p11ts","p11tf","p11co","p11x ","p11x2",
     & "p11x3","p11x4","p11x5","p12pv","p12ts","p12tf","p12co","p12x ",
     & "p12x2","p12x3","p12x4","p12x5","p13pv","p13ts","p13tf","p13co",
     & "p13x ","p13x2","p13x3","p13x4","p13x5","p14pv","p14ts","p14tf",
     & "p14co","p14x ","p14x2","p14x3","p14x4","p14x5","p15pv","p15ts",
     & "p15tf","p15co","p15x ","p15x2","p15x3","p15x4","p15x5","p16pv",
     & "p16ts","p16tf","p16co","p16x ","p16x2","p16x3","p16x4","p16x5",
     & "p17pv","p17ts","p17tf","p17co","p17x ","p17x2","p17x3","p17x4",
     & "p17x5","p18pv","p18ts","p18tf","p18co","p18x ","p18x2","p18x3",
     & "p18x4","p18x5","p19pv","p19ts","p19tf","p19co","p19x ","p19x2",
     & "p19x3","p19x4","p19x5","p20pv","p20ts","p20tf","p20co","p20x ",
     & "p20x2","p20x3","p20x4","p20x5","p21pv","p21ts","p21tf",
     & "p21co","p21x ","p21x2","p21x3","p21x4","p21x5","p22pv","p22ts",
     & "p22tf","p22co","p22x ","p22x2","p22x3","p22x4","p22x5","p23pv",
     & "p23ts","p23tf","p23co","p23x ","p23x2","p23x3","p23x4","p23x5",
     & "p24pv","p24ts","p24tf","p24co","p24x ","p24x2","p24x3","p24x4",
     & "p24x5","p25pv","p25ts","p25tf","p25co","p25x ","p25x2","p25x3",
     & "p25x4","p25x5","p26pv","p26ts","p26tf","p26co","p26x ","p26x2",
     & "p26x3","p26x4","p26x5","p27pv","p27ts","p27tf","p27co","p27x ",
     & "p27x2","p27x3","p27x4","p27x5","p28pv","p28ts","p28tf","p28co",
     & "p28x ","p28x2","p28x3","p28x4","p28x5","p29pv","p29ts","p29tf",
     & "p29co","p29x ","p29x2","p29x3","p29x4","p29x5","p30pv","p30ts",
     & "p30tf","p30co","p30x ","p30x2","p30x3","p30x4","p30x5","p31pv",
     & "p31ts","p31tf","p31co","p31x ","p31x2","p31x3","p31x4","p31x5",
     & "p32pv","p32ts","p32tf","p32co","p32x ","p32x2","p32x3","p32x4",
     & "p32x5","p33pv","p33ts","p33tf","p33co","p33x ","p33x2","p33x3",
     & "p33x4","p33x5","p34pv","p34ts","p34tf","p34co","p34x ","p34x2",
     & "p34x3","p34x4","p34x5","p35pv","p35ts","p35tf","p35co","p35x ",
     & "p35x2","p35x3","p35x4","p35x5","p36pv","p36ts","p36tf","p36co",
     & "p36x ","p36x2","p36x3","p36x4","p36x5","p37pv","p37ts","p37tf",
     & "p37co","p37x ","p37x2","p37x3","p37x4","p37x5","p38pv","p38ts",
     & "p38tf","p38co","p38x ","p38x2","p38x3","p38x4","p38x5","p39pv",
     & "p39ts","p39tf","p39co","p39x ","p39x2","p39x3","p39x4","p39x5",
     & "p40pv","p40ts","p40tf","p40co","p40x ","p40x2","p40x3","p40x4",
     & "p40x5","p41pv","p41ts","p41tf","p41co","p41x ","p41x2","p41x3",
     & "p41x4","p41x5","p42pv","p42ts","p42tf","p42co","p42x ","p42x2",
     & "p42x3","p42x4","p42x5","p43pv","p43ts","p43tf","p43co","p43x ",
     & "p43x2","p43x3","p43x4","p43x5","p44pv","p44ts","p44tf","p44co",
     & "p44x ","p44x2","p44x3","p44x4","p44x5","p45pv","p45ts","p45tf",
     & "p45co","p45x ","p45x2","p45x3","p45x4","p45x5","p46pv","p46ts",
     & "p46tf","p46co","p46x ","p46x2","p46x3","p46x4","p46x5","p47pv",
     & "p47ts","p47tf","p47co","p47x ","p47x2","p47x3","p47x4","p47x5",
     & "p48pv","p48ts","p48tf","p48co","p48x ","p48x2","p48x3","p48x4",
     & "p48x5","p49pv","p49ts","p49tf","p49co","p49x ","p49x2","p49x3",
     & "p49x4","p49x5","p50pv","p50ts","p50tf","p50co","p50x ","p50x2",
     & "p50x3","p50x4","p50x5" /

      data NE5 /  "  r_A","  r_B","LB/LA","    e","omega","b_pri",
     &                  "b_sec","Rchi2","    a","    q","  M_A","  M_B",
     &                  "  R_A","  R_B","loggA","loggB","rho_A","rho_B"/

      data N19 /
     &"Surf. bright. ratio","Sum of frac radii  ","Ratio of the radii ",
     &"Limb darkening A1  ","Limb darkening B1  ","Orbit inclination  ",
     &"ecc * cos(omega)   ","ecc * sin(omega)   ","Grav darkening A   ",
     &"Grav darkening B   ","Reflected light A  ","Reflected light B  ",
     &"Phot mass ratio    ","Tide lead/lag angle","Third light (L_3)  ",
     &"Phase correction   ","Light scale factor ","Integration ring   ",
     &"Orbital period (P) ","Ephemeris timebase ","Limb darkening A2  ",
     &"Limb darkening A3  ","Limb darkening A4  ","Limb darkening B2  ",
     &"Limb darkening B3  ","Limb darkening B4  ","RV amplitude star A",
     &"RV amplitude star B","Systemic RV star A ","Systemic RV star B ",
     &"Sine01 timebase    ","Sine01 period      ","Sine01 amplitude   ",
     &"Sine02 timebase    ","Sine02 period      ","Sine02 amplitude   ",
     &"Sine03 timebase    ","Sine03 period      ","Sine03 amplitude   ",
     &"Sine04 timebase    ","Sine04 period      ","Sine04 amplitude   ",
     &"Sine05 timebase    ","Sine05 period      ","Sine05 amplitude   ",
     &"Sine06 timebase    ","Sine06 period      ","Sine06 amplitude   ",
     &"Sine07 timebase    ","Sine07 period      ","Sine07 amplitude   ",
     &"Sine08 timebase    ","Sine08 period      ","Sine08 amplitude   ",
     &"Sine09 timebase    ","Sine09 period      ","Sine09 amplitude   ",
     &"Sine10 timebase    ","Sine10 period      ","Sine10 amplitude   ",
     &"Sine11 timebase    ","Sine11 period      ","Sine11 amplitude   ",
     &"Sine12 timebase    ","Sine12 period      ","Sine12 amplitude   ",
     &"Sine13 timebase    ","Sine13 period      ","Sine13 amplitude   ",
     &"Sine14 timebase    ","Sine14 period      ","Sine14 amplitude   ",
     &"Sine15 timebase    ","Sine15 period      ","Sine15 amplitude   ",
     &"Sine16 timebase    ","Sine16 period      ","Sine16 amplitude   ",
     &"Sine17 timebase    ","Sine17 period      ","Sine17 amplitude   ",
     &"Sine18 timebase    ","Sine18 period      ","Sine18 amplitude   ",
     &"Sine19 timebase    ","Sine19 period      ","Sine19 amplitude   ",
     &"Sine20 timebase    ","Sine20 period      ","Sine20 amplitude   ",
     &"Poly01 pivot       ","Poly01 time start  ","Poly01 time finish ",
     &"Poly01 constant    ","Poly01 coeff of x  ","Poly01 coeff of x^2",
     &"Poly01 coeff of x^3","Poly01 coeff of x^4","Poly01 coeff of x^5",
     &"Poly02 pivot       ","Poly02 time start  ","Poly02 time finish ",
     &"Poly02 constant    ","Poly02 coeff of x  ","Poly02 coeff of x^2",
     &"Poly02 coeff of x^3","Poly02 coeff of x^4","Poly02 coeff of x^5",
     &"Poly03 pivot       ","Poly03 time start  ","Poly03 time finish ",
     &"Poly03 constant    ","Poly03 coeff of x  ","Poly03 coeff of x^2",
     &"Poly03 coeff of x^3","Poly03 coeff of x^4","Poly03 coeff of x^5",
     &"Poly04 pivot       ","Poly04 time start  ","Poly04 time finish ",
     &"Poly04 constant    ","Poly04 coeff of x  ","Poly04 coeff of x^2",
     &"Poly04 coeff of x^3","Poly04 coeff of x^4","Poly04 coeff of x^5",
     &"Poly05 pivot       ","Poly05 time start  ","Poly05 time finish ",
     &"Poly05 constant    ","Poly05 coeff of x  ","Poly05 coeff of x^2",
     &"Poly05 coeff of x^3","Poly05 coeff of x^4","Poly05 coeff of x^5",
     &"Poly06 pivot       ","Poly06 time start  ","Poly06 time finish ",
     &"Poly06 constant    ","Poly06 coeff of x  ","Poly06 coeff of x^2",
     &"Poly06 coeff of x^3","Poly06 coeff of x^4","Poly06 coeff of x^5",
     &"Poly07 pivot       ","Poly07 time start  ","Poly07 time finish ",
     &"Poly07 constant    ","Poly07 coeff of x  ","Poly07 coeff of x^2",
     &"Poly07 coeff of x^3","Poly07 coeff of x^4","Poly07 coeff of x^5",
     &"Poly08 pivot       ","Poly08 time start  ","Poly08 time finish ",
     &"Poly08 constant    ","Poly08 coeff of x  ","Poly08 coeff of x^2",
     &"Poly08 coeff of x^3","Poly08 coeff of x^4","Poly08 coeff of x^5",
     &"Poly09 pivot       ","Poly09 time start  ","Poly09 time finish ",
     &"Poly09 constant    ","Poly09 coeff of x  ","Poly09 coeff of x^2",
     &"Poly09 coeff of x^3","Poly09 coeff of x^4","Poly09 coeff of x^5",
     &"Poly10 pivot       ","Poly10 time start  ","Poly10 time finish ",
     &"Poly10 constant    ","Poly10 coeff of x  ","Poly10 coeff of x^2",
     &"Poly10 coeff of x^3","Poly10 coeff of x^4","Poly10 coeff of x^5",
     &"Poly11 pivot       ","Poly11 time start  ","Poly11 time finish ",
     &"Poly11 constant    ","Poly11 coeff of x  ","Poly11 coeff of x^2",
     &"Poly11 coeff of x^3","Poly11 coeff of x^4","Poly11 coeff of x^5",
     &"Poly12 pivot       ","Poly12 time start  ","Poly12 time finish ",
     &"Poly12 constant    ","Poly12 coeff of x  ","Poly12 coeff of x^2",
     &"Poly12 coeff of x^3","Poly12 coeff of x^4","Poly12 coeff of x^5",
     &"Poly13 pivot       ","Poly13 time start  ","Poly13 time finish ",
     &"Poly13 constant    ","Poly13 coeff of x  ","Poly13 coeff of x^2",
     &"Poly13 coeff of x^3","Poly13 coeff of x^4","Poly13 coeff of x^5",
     &"Poly14 pivot       ","Poly14 time start  ","Poly14 time finish ",
     &"Poly14 constant    ","Poly14 coeff of x  ","Poly14 coeff of x^2",
     &"Poly14 coeff of x^3","Poly14 coeff of x^4","Poly14 coeff of x^5",
     &"Poly15 pivot       ","Poly15 time start  ","Poly15 time finish ",
     &"Poly15 constant    ","Poly15 coeff of x  ","Poly15 coeff of x^2",
     &"Poly15 coeff of x^3","Poly15 coeff of x^4","Poly15 coeff of x^5",
     &"Poly16 pivot       ","Poly16 time start  ","Poly16 time finish ",
     &"Poly16 constant    ","Poly16 coeff of x  ","Poly16 coeff of x^2",
     &"Poly16 coeff of x^3","Poly16 coeff of x^4","Poly16 coeff of x^5",
     &"Poly17 pivot       ","Poly17 time start  ","Poly17 time finish ",
     &"Poly17 constant    ","Poly17 coeff of x  ","Poly17 coeff of x^2",
     &"Poly17 coeff of x^3","Poly17 coeff of x^4","Poly17 coeff of x^5",
     &"Poly18 pivot       ","Poly18 time start  ","Poly18 time finish ",
     &"Poly18 constant    ","Poly18 coeff of x  ","Poly18 coeff of x^2",
     &"Poly18 coeff of x^3","Poly18 coeff of x^4","Poly18 coeff of x^5",
     &"Poly19 pivot       ","Poly19 time start  ","Poly19 time finish ",
     &"Poly19 constant    ","Poly19 coeff of x  ","Poly19 coeff of x^2",
     &"Poly19 coeff of x^3","Poly19 coeff of x^4","Poly19 coeff of x^5",
     &"Poly20 pivot       ","Poly20 time start  ","Poly20 time finish ",
     &"Poly20 constant    ","Poly20 coeff of x  ","Poly20 coeff of x^2",
     &"Poly20 coeff of x^3","Poly20 coeff of x^4","Poly20 coeff of x^5",
     &"Poly21 pivot       ","Poly21 time start  ","Poly21 time finish ",
     &"Poly21 constant    ","Poly21 coeff of x  ","Poly21 coeff of x^2",
     &"Poly21 coeff of x^3","Poly21 coeff of x^4","Poly21 coeff of x^5",
     &"Poly22 pivot       ","Poly22 time start  ","Poly22 time finish ",
     &"Poly22 constant    ","Poly22 coeff of x  ","Poly22 coeff of x^2",
     &"Poly22 coeff of x^3","Poly22 coeff of x^4","Poly22 coeff of x^5",
     &"Poly23 pivot       ","Poly23 time start  ","Poly23 time finish ",
     &"Poly23 constant    ","Poly23 coeff of x  ","Poly23 coeff of x^2",
     &"Poly23 coeff of x^3","Poly23 coeff of x^4","Poly23 coeff of x^5",
     &"Poly24 pivot       ","Poly24 time start  ","Poly24 time finish ",
     &"Poly24 constant    ","Poly24 coeff of x  ","Poly24 coeff of x^2",
     &"Poly24 coeff of x^3","Poly24 coeff of x^4","Poly24 coeff of x^5",
     &"Poly25 pivot       ","Poly25 time start  ","Poly25 time finish ",
     &"Poly25 constant    ","Poly25 coeff of x  ","Poly25 coeff of x^2",
     &"Poly25 coeff of x^3","Poly25 coeff of x^4","Poly25 coeff of x^5",
     &"Poly26 pivot       ","Poly26 time start  ","Poly26 time finish ",
     &"Poly26 constant    ","Poly26 coeff of x  ","Poly26 coeff of x^2",
     &"Poly26 coeff of x^3","Poly26 coeff of x^4","Poly26 coeff of x^5",
     &"Poly27 pivot       ","Poly27 time start  ","Poly27 time finish ",
     &"Poly27 constant    ","Poly27 coeff of x  ","Poly27 coeff of x^2",
     &"Poly27 coeff of x^3","Poly27 coeff of x^4","Poly27 coeff of x^5",
     &"Poly28 pivot       ","Poly28 time start  ","Poly28 time finish ",
     &"Poly28 constant    ","Poly28 coeff of x  ","Poly28 coeff of x^2",
     &"Poly28 coeff of x^3","Poly28 coeff of x^4","Poly28 coeff of x^5",
     &"Poly29 pivot       ","Poly29 time start  ","Poly29 time finish ",
     &"Poly29 constant    ","Poly29 coeff of x  ","Poly29 coeff of x^2",
     &"Poly29 coeff of x^3","Poly29 coeff of x^4","Poly29 coeff of x^5",
     &"Poly30 pivot       ","Poly30 time start  ","Poly30 time finish ",
     &"Poly30 constant    ","Poly30 coeff of x  ","Poly30 coeff of x^2",
     &"Poly30 coeff of x^3","Poly30 coeff of x^4","Poly30 coeff of x^5",
     &"Poly31 pivot       ","Poly31 time start  ","Poly31 time finish ",
     &"Poly31 constant    ","Poly31 coeff of x  ","Poly31 coeff of x^2",
     &"Poly31 coeff of x^3","Poly31 coeff of x^4","Poly31 coeff of x^5",
     &"Poly32 pivot       ","Poly32 time start  ","Poly32 time finish ",
     &"Poly32 constant    ","Poly32 coeff of x  ","Poly32 coeff of x^2",
     &"Poly32 coeff of x^3","Poly32 coeff of x^4","Poly32 coeff of x^5",
     &"Poly33 pivot       ","Poly33 time start  ","Poly33 time finish ",
     &"Poly33 constant    ","Poly33 coeff of x  ","Poly33 coeff of x^2",
     &"Poly33 coeff of x^3","Poly33 coeff of x^4","Poly33 coeff of x^5",
     &"Poly34 pivot       ","Poly34 time start  ","Poly34 time finish ",
     &"Poly34 constant    ","Poly34 coeff of x  ","Poly34 coeff of x^2",
     &"Poly34 coeff of x^3","Poly34 coeff of x^4","Poly34 coeff of x^5",
     &"Poly35 pivot       ","Poly35 time start  ","Poly35 time finish ",
     &"Poly35 constant    ","Poly35 coeff of x  ","Poly35 coeff of x^2",
     &"Poly35 coeff of x^3","Poly35 coeff of x^4","Poly35 coeff of x^5",
     &"Poly36 pivot       ","Poly36 time start  ","Poly36 time finish ",
     &"Poly36 constant    ","Poly36 coeff of x  ","Poly36 coeff of x^2",
     &"Poly36 coeff of x^3","Poly36 coeff of x^4","Poly36 coeff of x^5",
     &"Poly37 pivot       ","Poly37 time start  ","Poly37 time finish ",
     &"Poly37 constant    ","Poly37 coeff of x  ","Poly37 coeff of x^2",
     &"Poly37 coeff of x^3","Poly37 coeff of x^4","Poly37 coeff of x^5",
     &"Poly38 pivot       ","Poly38 time start  ","Poly38 time finish ",
     &"Poly38 constant    ","Poly38 coeff of x  ","Poly38 coeff of x^2",
     &"Poly38 coeff of x^3","Poly38 coeff of x^4","Poly38 coeff of x^5",
     &"Poly39 pivot       ","Poly39 time start  ","Poly39 time finish ",
     &"Poly39 constant    ","Poly39 coeff of x  ","Poly39 coeff of x^2",
     &"Poly39 coeff of x^3","Poly39 coeff of x^4","Poly39 coeff of x^5",
     &"Poly40 pivot       ","Poly40 time start  ","Poly40 time finish ",
     &"Poly40 constant    ","Poly40 coeff of x  ","Poly40 coeff of x^2",
     &"Poly40 coeff of x^3","Poly40 coeff of x^4","Poly40 coeff of x^5",
     &"Poly41 pivot       ","Poly41 time start  ","Poly41 time finish ",
     &"Poly41 constant    ","Poly41 coeff of x  ","Poly41 coeff of x^2",
     &"Poly41 coeff of x^3","Poly41 coeff of x^4","Poly41 coeff of x^5",
     &"Poly42 pivot       ","Poly42 time start  ","Poly42 time finish ",
     &"Poly42 constant    ","Poly42 coeff of x  ","Poly42 coeff of x^2",
     &"Poly42 coeff of x^3","Poly42 coeff of x^4","Poly42 coeff of x^5",
     &"Poly43 pivot       ","Poly43 time start  ","Poly43 time finish ",
     &"Poly43 constant    ","Poly43 coeff of x  ","Poly43 coeff of x^2",
     &"Poly43 coeff of x^3","Poly43 coeff of x^4","Poly43 coeff of x^5",
     &"Poly44 pivot       ","Poly44 time start  ","Poly44 time finish ",
     &"Poly44 constant    ","Poly44 coeff of x  ","Poly44 coeff of x^2",
     &"Poly44 coeff of x^3","Poly44 coeff of x^4","Poly44 coeff of x^5",
     &"Poly45 pivot       ","Poly45 time start  ","Poly45 time finish ",
     &"Poly45 constant    ","Poly45 coeff of x  ","Poly45 coeff of x^2",
     &"Poly45 coeff of x^3","Poly45 coeff of x^4","Poly45 coeff of x^5",
     &"Poly46 pivot       ","Poly46 time start  ","Poly46 time finish ",
     &"Poly46 constant    ","Poly46 coeff of x  ","Poly46 coeff of x^2",
     &"Poly46 coeff of x^3","Poly46 coeff of x^4","Poly46 coeff of x^5",
     &"Poly47 pivot       ","Poly47 time start  ","Poly47 time finish ",
     &"Poly47 constant    ","Poly47 coeff of x  ","Poly47 coeff of x^2",
     &"Poly47 coeff of x^3","Poly47 coeff of x^4","Poly47 coeff of x^5",
     &"Poly48 pivot       ","Poly48 time start  ","Poly48 time finish ",
     &"Poly48 constant    ","Poly48 coeff of x  ","Poly48 coeff of x^2",
     &"Poly48 coeff of x^3","Poly48 coeff of x^4","Poly48 coeff of x^5",
     &"Poly49 pivot       ","Poly49 time start  ","Poly49 time finish ",
     &"Poly49 constant    ","Poly49 coeff of x  ","Poly49 coeff of x^2",
     &"Poly49 coeff of x^3","Poly49 coeff of x^4","Poly49 coeff of x^5",
     &"Poly50 pivot       ","Poly50 time start  ","Poly50 time finish ",
     &"Poly50 constant    ","Poly50 coeff of x  ","Poly50 coeff of x^2",
     &"Poly50 coeff of x^3","Poly50 coeff of x^4","Poly50 coeff of x^5"/

      data NE19 /
     &"Frac radius star A ","Frac radius star B ","Light ratio (B/A)  ",
     &"       Eccentricity","Periastronlongitude","Impact par (pri ec)",
     &"Impact par (sec ec)","Reduced chi-squared","Semimajor axis Rsun",
     &"Mass ratio from RVs","Star A mass    Msun","Star B mass    Msun",
     &"Star A radius  Rsun","Star B radius  Rsun","Star A logg   c.g.s",
     &"Star B logg   c.g.s","Star A density  sun","Star B density  sun"/

      do i = 1,540
        NAME5(i) = N5(i)
        NAME19(i) = N19(i)
      end do

      do i = 1,18
        NAMEXTRA5(i) = NE5(i)
        NAMEXTRA19(i) = NE19(i)
      end do

      if ( V(2) < 0.0d0 ) then                ! if fitting for radii and
        NAME5(2) = "  r_A"                    ! not their sum and ratio
        NAME5(3) = "  r_B"
        NAMEXTRA5(1) = "rA+rB"
        NAMEXTRA5(2) = "    k"
        NAME19(2) = "Frac radius star A "
        NAME19(3) = "Frac radius star B "
        NAMEXTRA19(1) = "Sum of frac radii  "
        NAMEXTRA19(2) = "Ratio of the radii "
      end if

      if ( V(7) > 5.0d0 ) then               ! if fitting for e and w
        NAME5(7) = "    e"                   ! not for ecosw and esinw
        NAME5(8) = "omega"
        NAMEXTRA5(4) = "ecosw"
        NAMEXTRA5(5) = "esinw"
        NAME19(7) = "Eccentricity       "
        NAME19(8) = "Periastronlongitude"
        NAMEXTRA19(4) = "ecc * cos(omega)   "
        NAMEXTRA19(5) = "ecc * sin(omega)   "
      end if

      END SUBROUTINE GETNAME
!=======================================================================
      SUBROUTINE GETEXTRAS (V,VARY,NDATA,LDTYPE,CHISQ,NSINE,PSINE,NPOLY,
     &                                    PPOLY,NUMINT,NINTERVAL,VEXTRA)

            ! This subroutine calculates various dependent quantities of
            ! interest from the fitted quantities.

      implicit none
      real*8 V(540)                       ! IN:  Parameters of the fit
      integer VARY(540)                   ! IN:  Which parameters fitted
      integer NDATA                       ! IN:  Number of datapoints
      integer LDTYPE(2)                   ! IN:  Type of LD law for star
      real*8 CHISQ                        ! IN:  Chi-squared of the fit
      integer NSINE,PSINE(20)             ! IN:  Which par for each sine
      integer NPOLY,PPOLY(50)             ! IN:  Which par for each poly
      integer NUMINT                      ! IN:  Number numerical integs
      real*8  NINTERVAL                   ! IN:  Time interval num.integ
      real*8 VEXTRA(18)                   ! OUT:  Dependent quantities
      integer NVARY                       ! LOCAL: number of fitted pars
      real*8 PI,DEG2RAD                   ! LOCAL: useful constants
      real*8 ACONST,MCONST                ! LOCAL: orbital constants
      real*8 EFAC                         ! LOCAL: (1-e^2) quantity
      real*8 SINI,COSI                    ! LOCAL: useful quantities
      real*8 R1,R2                        ! LOCAL: fractional radii
      real*8 RSUN,GMSUN,GSUN              ! LOCAL: physical constants
      real*8 LP,LS                        ! LOCAL: EBOP/GETMODEL outputs
      real*8 K1,K2,KSUM                   ! LOCAL: velocity amplitudes
      real*8 ECQPHASES(6)                 ! SUB:   useful orbital phases
      integer i                           ! LOCAL: loop counter
      real*8 GETMODEL                     ! FUNCTION: gets model value
      real*8 ECC,OMEGA,ECOSW,ESINW        ! LOCAL: eccentricity values

            ! VEXTRA(1) = star A fractional radius
            ! VEXTRA(2) = star B fractional radius
            ! VEXTRA(3) = light ratio starB/starA
            ! VEXTRA(4) = eccentricity
            ! VEXTRA(5) = periastron longitude
            ! VEXTRA(6) = impact par (pri ecl)
            ! VEXTRA(7) = impact par (sec ecl)
            ! VEXTRA(8) = reduced chi-squared
            ! VEXTRA(9) = orbital semimajor axis
            ! VEXTRA(10) = mass ratio from RVs
            ! VEXTRA(11) = mass of star A (Msun)
            ! VEXTRA(12) = mass of star B (Msun)
            ! VEXTRA(13) = radius of star A (Rsun)
            ! VEXTRA(14) = radius of star B (Rsun)
            ! VEXTRA(15) = logg of star A (c.g.s.)
            ! VEXTRA(16) = logg of star B (c.g.s.)
            ! VEXTRA(17) = density of star A (solar units)
            ! VEXTRA(18) = density of star B (solar units)

      do i = 1,18
        VEXTRA(i) = -999.0d0
      end do

      PI = atan(1.0d0) * 4.0d0
      DEG2RAD = 45.0d0 / atan(1.0d0)
      RSUN = 6.957d8                      ! The IAU 2015 resolution ....
      GMSUN = 1.3271244d20                !   .... B3 nominal quantities

      GSUN = log10(GMSUN*100.0d0/RSUN**2)
      SINI = sin( V(6) / DEG2RAD )
      COSI = cos( V(6) / DEG2RAD )
      MCONST = ( 1.0d9 * 86400.0d0 ) / ( 2.0d0 * PI * GMSUN )
      ACONST = ( 1.0d3 * 86400.0d0 ) / ( 2.0d0 * PI * RSUN )

            ! Calculate the reduced chi-squared (need number of fit par)

      NVARY = 0
      do i = 1,540
        if ( VARY(i) /= 0 ) NVARY = NVARY + 1
      end do
      VEXTRA(8) = CHISQ

            ! First calculate r_A and r_B or their combination depending
            ! on whether (r_A+r_B,k) or (r_A,r_B) were the fitted ones.

      if ( V(2) >= 0.0d0 ) then
        VEXTRA(1) = V(2) / (1.0d0 + V(3))                   ! r_A
        VEXTRA(2) = V(2) / (1.0d0 + (1.0d0/V(3)))           ! r_B
        R1 = VEXTRA(1)
        R2 = VEXTRA(2)
      else
        VEXTRA(1) = abs(V(2)) + V(3)                        ! rA+rB
        VEXTRA(2) = V(3) / abs(V(2))                        ! k=r_B/r_A
        R1 = abs(V(2))
        R2 = V(3)
      end if

            ! Get the light ratio at the time of primary quadrature

      call ECQUADPHASES(V(7),V(8),V(16),ECQPHASES)
      VEXTRA(3) = GETMODEL (V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                V(20)+V(19)*ECQPHASES(2),2,LP,LS,NUMINT,NINTERVAL)

            ! Get (e,w) or (ecosw,esinw) depending on which was fitted
            ! Also calculate the eccentricity factor for use below,

      CALL GETEOMEGA (V(7),V(8),ECC,OMEGA,ECOSW,ESINW)

      if ( V(7) > 5.0d0 ) then
        VEXTRA(4) = ECOSW
        VEXTRA(5) = ESINW
      else
        VEXTRA(4) = ECC
        VEXTRA(5) = OMEGA
      end if

      EFAC = 1.0d0 - ECC**2

            ! calculate impact parameters using these equations:
            ! b(pri) = (cos i / rA) * (1-e^2 / 1+esinw)
            ! b(sec) = (cos i / rA) * (1-e^2 / 1-esinw)
            ! modified from Winn (http://arxiv.org/abs/1001.2010).

      VEXTRA(6) = EFAC / (1.0d0 + ESINW) * COSI / (R1)
      VEXTRA(7) = EFAC / (1.0d0 - ESINW) * COSI / (R1)

      ! calculate these quantities if possible (if K1 and K2 are known):
      !   VEXTRA(9) = semimajor axis      VEXTRA(14) = radius of star B
      !   VEXTRA(10) = mass ratio         VEXTRA(15) = logg of star A
      !   VEXTRA(11) = mass of star A     VEXTRA(16) = logg of star B
      !   VEXTRA(12) = mass of star B     VEXTRA(17) = density of star A
      !   VEXTRA(13) = radius of star A   VEXTRA(18) = density of star B

      if ( V(27) >= -1.0d9 .and. V(28) >= -1.0d9 ) then
        K1 = V(27)
        K2 = V(28)
        if ( K1 < 0.0d0 .and. K2 < 0.0d0 ) then
          write (6,'(A37,A43)') "## Warning: velocity amplitudes were ",
     &                    "both negative so have been swapped for the "
          write (6,'(A37,A43)') "## calculation of the physical proper",
     &                    "ties of the stars. Check the input data.   "
          K1 = 0.0d0 - V(28)
          K2 = 0.0d0 - V(27)
        end if
        if ( K1 < 0.0d0 ) write (6,'(A58,E15.8)')
     &   "## Warning: the velocity amplitude of star A is negative: ",K1
        if ( K2 < 0.0d0 ) write (6,'(A58,E15.8)')
     &   "## Warning: the velocity amplitude of star B is negative: ",K2

        KSUM = K1 + K2
        VEXTRA(9) = ACONST * sqrt(EFAC) * KSUM * V(19) / SINI
        VEXTRA(10) = K1 / K2
        VEXTRA(11) = MCONST * EFAC**1.5d0 * KSUM**2 *K2* V(19) / SINI**3
        VEXTRA(12) = MCONST * EFAC**1.5d0 * KSUM**2 *K1* V(19) / SINI**3
        VEXTRA(13) = VEXTRA(9) * R1
        VEXTRA(14) = VEXTRA(9) * R2
        VEXTRA(15) = GSUN + log10(VEXTRA(11)) -(2.0d0*log10(VEXTRA(13)))
        VEXTRA(16) = GSUN + log10(VEXTRA(12)) -(2.0d0*log10(VEXTRA(14)))
        VEXTRA(17) = VEXTRA(11) / VEXTRA(13)**3
        VEXTRA(18) = VEXTRA(12) / VEXTRA(14)**3
      end if

      END SUBROUTINE GETEXTRAS
!=======================================================================
      DOUBLEPRECISION FUNCTION GETMODEL (V,VARY,LDTYPE,NSINE,PSINE,NPOLY
     &                        ,PPOLY,TIME,DTYPE1,LA,LB,NUMINT,NINTERVAL)
            ! Output a predicted model value according to the parameters
            ! in array V. Precise meaning of the value depends on DTYPE.
            ! DTYPE1=1  it outputs an EBOP magnitude for given time
            ! DTYPE1=2  it outputs a light ratio for the given time
            ! DTYPE1=3  outputs a time of eclipse for the given =CYCLE=
            ! DTYPE1=4  it simply outputs the third light value
            ! DTYPE1=5  it outputs e or e*cos(omega)
            ! DTYPE1=6  it outputs omega or e*sin(omega)
            ! DTYPE1=7  it outputs the RV of star A
            ! DTYPE1=8  it outputs the RV of star B
      implicit none
      real*8 V(540)                  ! IN: Photometric parameters
      integer VARY(540)              ! IN: Which parameters are fitted
      integer LDTYPE(2)             ! IN: LD law type for the two stars
      real*8 TIME                   ! IN: The given TIME, PHASE or CYCLE
      integer DTYPE1                ! IN: 1-8 depending on wanted result
      integer NSINE,PSINE(20)       ! IN: number and parameters of sines
      integer NPOLY,PPOLY(50)       ! IN: number and parameters of polys
      integer NUMINT                ! IN: Number of numerical integratns
      real*8 NINTERVAL              ! IN: Time interval for integrations
      real*8 LA,LB                  ! OUT: Light produced by each star
      real*8 FMAG,LP,LS             ! LOCAL: LIGHT subroutine output
      real*8 ECC,OMEGA,ECOSW,ESINW  ! LOCAL: orbital shape parameters
      real*8 GETMIN                 ! FUNCTION: returns time of minimum
!       real*8 GETPHASE               ! FUNCTION: returns orbital phase
      real*8 GETRV                  ! FUNCTION: returns the RV of a star
      real*8 FMAGSUM,LASUM,LBSUM
      integer i
      real*8 TIMEIN

      GETMODEL = 0.0d0

      if ( DTYPE1 == 1 ) then
        if ( NUMINT == 1 ) then
         CALL LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIME,FMAG,LP,LS)
          LA = LP
          LB = LS
          GETMODEL = FMAG

        else if ( NUMINT > 1 ) then
          FMAGSUM = 0.0d0
          LASUM = 0.0d0
          LBSUM = 0.0d0
          CALL LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIME,FMAG,LP,LS)
          do i = 1,NUMINT
            TIMEIN  =  TIME  -  NINTERVAL / 86400.0d0 / dble(NUMINT) *
     &        ( dble(NUMINT) - 2.0d0*dble(i) + 1.0d0 ) / 2.0d0
          CALL LIGHT(V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIMEIN,FMAG,LP,LS)
            FMAGSUM = FMAGSUM + FMAG
            LASUM = LASUM + LP
            LBSUM = LBSUM + LS
          end do
          FMAG = FMAGSUM / NUMINT
          LA = LASUM / NUMINT
          LB = LBSUM / NUMINT
          GETMODEL = FMAG
        else
          write(6,*)"NUMINT is less than 1 in function GETMODEL. Abort."
          write(6,*)"NUMINT =    ", NUMINT
          write(6,*)"NINTERVAL = ", NINTERVAL
          stop
        end if

      else if ( DTYPE1 == 2 ) then
        CALL LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIME,FMAG,LP,LS)
        LA = LP
        LB = LS
        GETMODEL = LS / LP

      else if ( DTYPE1 == 3 ) then
        GETMODEL = GETMIN (V(20),V(19),V(7),V(8),TIME)

      else if ( DTYPE1 == 4 ) then
        GETMODEL = V(15)

      else if ( DTYPE1 == 5 .or. DTYPE1 == 6 ) then
        CALL GETEOMEGA (V(7),V(8),ECC,OMEGA,ECOSW,ESINW)
        if ( V(7) > 5.0d0 ) then
          if ( DTYPE1 == 5 )  GETMODEL = ECC
          if ( DTYPE1 == 6 )  GETMODEL = OMEGA
        else
          if ( DTYPE1 == 5 )  GETMODEL = ECOSW
          if ( DTYPE1 == 6 )  GETMODEL = ESINW
        end if

      else if ( DTYPE1 == 7 ) then
        GETMODEL = GETRV (TIME,V,NSINE,PSINE,NPOLY,PPOLY,'A')

      else if ( DTYPE1 == 8 ) then
        if ( VARY(30) == -1 ) then
          GETMODEL = GETRV (TIME,V,NSINE,PSINE,NPOLY,PPOLY,'C')
        else
          GETMODEL = GETRV (TIME,V,NSINE,PSINE,NPOLY,PPOLY,'B')
        end if

      else
        GETMODEL = -100.0d0
        write(6,*) "### ERROR: wrong datatype asked for in GETMODEL: ",
     &              DTYPE1
        STOP
      end if

      END FUNCTION GETMODEL
!=======================================================================
      DOUBLEPRECISION FUNCTION GETPHASE (HJD,PERIOD,TZERO)
            ! Returns phase from given time and orbital ephemeris
      implicit none
      real*8 HJD,PERIOD,TZERO

      GETPHASE = (HJD - TZERO) / PERIOD
      GETPHASE = GETPHASE - int(GETPHASE)
      if ( GETPHASE < 0.0d0 ) GETPHASE = GETPHASE + 1.0d0

      END FUNCTION GETPHASE
!=======================================================================
      SUBROUTINE ECQUADPHASES (ECCIN,OMEGAIN,PSHIFT,PHASES)
            ! Calculates orbital phases of the eclipses and quadratures.
            ! PHASES(1) is the phase of primary eclipse
            ! PHASES(2) is the phase of secondary eclipse
            ! PHASES(3) is the phase of *photometric* primary quadrature
            ! PHASES(4) is the phase of *photometric* secndry quadrature
            ! PHASES(5) is the phase of periastron
            ! PHASES(6) is the phase of apastron
      implicit none
      real*8 ECCIN,OMEGAIN          ! IN: eccentricity and peri.long.
      real*8 PSHIFT                 ! IN: time of primary eclipse
      real*8 PHASES(6)              ! OUT: phases of eclipses and quads
      real*8 PI,DEG2RAD             ! LOCAL: constants
      real*8 ECC,OMEGA              ! LOCAL: values of ecc and peri.long
      real*8 ECOSW,ESINW            ! LOCAL: values of combination terms
      real*8 EFAC                   ! LOCAL: useful eccentricity factor
      real*8 TERM1,TERM2            ! LOCAL: calculation helper varibles
      real*8 PHASEDIFF              ! LOCAL: diff of prim and sec minima
      real*8 OMEGARAD               ! LOCAL: periastron longitude in rad
      real*8 TRUEANOM               ! LOCAL: true anomaly at primary ecl
      real*8 SQRTE,TANEDIV2         ! LOCAL: useful intermediate numbers
      real*8 ECCANOM                ! LOCAL: eccentric anomaly at pr ecl
      real*8 MEANANOM               ! LOCAL: mean anomaly at primary ecl

      PI = atan(1.0d0) * 4.0d0
      DEG2RAD = 45.0d0 / atan(1.0d0)

            ! Get actual eccentricity and periastron longitude values
            ! from the input, which could be (e+10,w) or (ecosw,esinw)

      if ( ECCIN > 10.0d0 ) then
        ECC = ECCIN - 10.0d0
        OMEGA = OMEGAIN / DEG2RAD
        ECOSW = ECC * cos(OMEGA)
        ESINW = ECC * sin(OMEGA)
      else
        ECC = sqrt(ECCIN**2 + OMEGAIN**2)
        OMEGA = atan2(OMEGAIN,ECCIN)
        ECOSW = ECCIN
        ESINW = OMEGAIN
      end if

      EFAC = sqrt(1.0d0 - ECC**2)

!        write(6,*)" "
!        write(6,'(a16,2(f13.8))')"pi,deg2rad      ",pi,deg2rad
!        write(6,'(a16,2(f13.8))')"eccin,omegain   ",eccin,omegain
!        write(6,'(a16,2(f13.8))')"ecc,omega       ",ecc,omega
!        write(6,'(a16,2(f13.8))')"ecosw,esinw     ",ecosw,esinw
!        write(6,'(a16,2(f13.8))')"efac,pshift     ",efac,pshift

! The equation for phase difference comes from Hilditch (2001) page 238
!  equation 5.66, originally credited to the monograph by Kopal (1959).

      TERM1 = 2.0d0 * atan( ECOSW / EFAC )
      TERM2 = 2.0d0 * ECOSW * EFAC / (1.0d0 - ESINW**2)
      PHASEDIFF = ( PI + TERM1 + TERM2 ) / ( 2.0d0 * PI )

      PHASES(1) = PSHIFT
      PHASES(2) = PSHIFT + PHASEDIFF/2.0d0
      PHASES(3) = PSHIFT + PHASEDIFF
      PHASES(4) = PSHIFT + PHASEDIFF/2.0d0 + 0.50d0

! This set of equations for calculating the phase of periastron was
! obtained using the MCMC code of Andrew Cameron as a starting point.

      OMEGARAD = OMEGA * PI / 180.0d0
      TRUEANOM = (PI / 2.0d0) - OMEGA
      SQRTE = sqrt( (1.0d0 + ECC) / (1.0d0 - ECC) )
      TANEDIV2 = tan(TRUEANOM / 2.0d0) / SQRTE
      ECCANOM = 2.0d0 * atan( TANEDIV2 )
      MEANANOM = ECCANOM - (ECC * sin(ECCANOM))
      PHASES(5) = PHASES(1) - (MEANANOM / ( PI * 2.0d0))
      PHASES(6) = PHASES(5) + 0.5d0

!       write(6,*)" "
!       write(6,'(A22,2(F14.8))')"Peri long (rad) =    ", OMEGA
!       write(6,'(A22,2(F14.8))')"True anomaly =       ", TRUEANOM
!       write(6,'(A22,2(F14.8))')"Eccentric anomaly =  ", ECCANOM
!       write(6,'(A22,2(F14.8))')"Mean anomaly =       ", MEANANOM
!       write(6,'(A22,2(F14.8))')"Phase of periastron =", PERIPHASE

      END SUBROUTINE ECQUADPHASES
!=======================================================================
      DOUBLEPRECISION FUNCTION GETRV (TIME,V,NSINE,PSINE,NPOLY,PPOLY,
     $                                WHICHSTAR)
            ! Returns a radial velocity of a star at the given time TIME
            ! WHICHSTAR='A' returns an RV for star A
            ! WHICHSTAR='B' returns an RV for star B
            ! WHICHSTAR='C' returns RV for star B with Vsys from star A
      implicit none
      real*8 TIME                   ! IN: time to calculate the RV for
      real*8 V(540)                 ! IN: fitted parameters
      integer NSINE,PSINE(20)       ! IN: number and parameters of sines
      integer NPOLY,PPOLY(50)       ! IN: number and parameters of polys
      character*1 WHICHSTAR         ! IN: 'A' or 'B' to specify the star
      real*8 PI,DEG2RAD             ! LOCAL: useful constants
      real*8 PHASES(6)              ! LOCAL: useful orbital phases
      real*8 TRUEANOM               ! LOCAL: true anomaly
      real*8 MEANANOM               ! LOCAL: mean anomaly
      real*8 ECCANOM                ! LOCAL: eccentric anomaly
      real*8 KA,KB,VA,VB            ! LOCAL: orbital parameters
      real*8 KVAL                   ! LOCAL: velocity amplitude (km/s)
      real*8 VSYS                   ! LOCAL: systemic velocity (km/s)
      real*8 ECC,OMEGA,ECOSW,ESINW  ! LOCAL: orbital shape variables
      real*8 TANTD2,PHIANOM,GUESS   ! LOCAL: helper variables
      real*8 SINTERM,POLYTERM       ! LOCAL: sine and poly term values
      real*8 SINT,SINP,SINA         ! LOCAL: sine parameter values
      integer i,j                   ! LOCAL: loop counters
      real*8 GETPHASE               ! FUNCTION: calculate orbital phase

      PI = atan(1.0d0) * 4.0d0
      DEG2RAD = 45.0d0 / atan(1.0d0)

      CALL GETEOMEGA (V(7),V(8),ECC,OMEGA,ECOSW,ESINW)

            ! Modification of 2015/04/24:
            ! Now need to apply sine and polynomial functions to the
            ! orbit parameters if required. Copy them from the V array
            ! to avoid modifying the V array. Then check each sine and
            ! poly to see what parameter it applies to, and apply the
            ! term if relevant. All are additive in this case.

      KA = V(27)
      KB = V(28)
      VA = V(29)
      VB = V(30)

      if ( NSINE > 0 ) then
        do i = 1,NSINE
          SINT = V(28+i*3)      ! sine reference time
          SINP = V(29+i*3)      ! sine period
          SINA = V(30+i*3)      ! sine amplitude
          SINTERM = SINA * sin( 2.0d0 * PI * (TIME-SINT) / SINP )
          if ( PSINE(i) == 27 )  KA = KA + SINTERM
          if ( PSINE(i) == 28 )  KB = KB + SINTERM
          if ( PSINE(i) == 29 )  VA = VA + SINTERM
          if ( PSINE(i) == 30 )  VB = VB + SINTERM

          if ( PSINE(i) == 7 .or. PSINE(i) == 8 ) then
            if ( V(7) >= 9.99d0 ) then
              if ( PSINE(i) == 7 ) ECC = ECC + SINTERM
              if ( PSINE(i) == 8 ) OMEGA = OMEGA + SINTERM
              ECOSW = ECC * cos(OMEGA/DEG2RAD)
              ESINW = ECC * sin(OMEGA/DEG2RAD)
            endif
            if ( V(7) < 9.99d0 ) then
              if ( PSINE(i) == 7 ) ECOSW = ECOSW + SINTERM
              if ( PSINE(i) == 8 ) ESINW = ESINW + SINTERM
              CALL GETEW (ECOSW,ESINW,ECC,OMEGA)
            endif
          endif
        end do
      end if

      if ( NPOLY > 0 ) then
        do i = 1,NPOLY
          j = 91 + 9*(i-1)
          POLYTERM = V(j+3) + V(j+4)*(TIME-V(j))
     &                      + V(j+5)*((TIME-V(j))**2)
     &                      + V(j+6)*((TIME-V(j))**3)
     &                      + V(j+7)*((TIME-V(j))**4)
     &                      + V(j+8)*((TIME-V(j))**5)
          if ( TIME < V(j+1) ) POLYTERM = 0.0d0
          if ( TIME > V(j+2) ) POLYTERM = 0.0d0
          if ( PPOLY(i) == 27 )  KA = KA + POLYTERM
          if ( PPOLY(i) == 28 )  KB = KB + POLYTERM
          if ( PPOLY(i) == 29 )  VA = VA + POLYTERM
          if ( PPOLY(i) == 30 )  VB = VB + POLYTERM

          if ( PPOLY(i) == 7 .or. PPOLY(i) == 8 ) then
            if ( V(7) >= 9.99d0 ) then
              if ( PPOLY(i) == 7 ) ECC = ECC + POLYTERM
              if ( PPOLY(i) == 8 ) OMEGA = OMEGA + POLYTERM
              ECOSW = ECC * cos(OMEGA/DEG2RAD)
              ESINW = ECC * sin(OMEGA/DEG2RAD)
            endif
            if ( V(7) < 9.99d0 ) then
              if ( PPOLY(i) == 7 ) ECOSW = ECOSW + POLYTERM
              if ( PPOLY(i) == 8 ) ESINW = ESINW + POLYTERM
              CALL GETEW (ECOSW,ESINW,ECC,OMEGA)
            endif
          endif
        end do
      end if

            ! First calculate the mean anomaly at this time using the
            ! time of periastron from the ECQUADPHASES subroutine.

      CALL ECQUADPHASES (V(7),V(8),V(16),PHASES)
      PHIANOM = getphase(TIME,V(19),V(20)) - PHASES(5)
      MEANANOM = 2.0d0 * PI * mod(PHIANOM,1.0d0)

            ! Now calculate the true anomaly iteratively. This approach
            ! uses the MCMC code of Andrew Cameron as a starting point.

      GUESS = MEANANOM
      do i = 1,100
         ECCANOM = MEANANOM + ECC*sin(GUESS)
         if ( abs(ECCANOM-GUESS) < 1.0d-5 ) exit
         GUESS = ECCANOM
         if (i >= 100) write (6,'(A11,A54,A31,4(f11.6))') "## Warning:",
     &         " calculation of true anomaly in function GETRV did not",
     &         " converge. e, w, ecosw, esinw =", ECC,OMEGA,ECOSW,ESINW
      end do

      TANTD2 = sqrt((1.0d0+ECC)/(1.0d0-ECC)) * tan(ECCANOM/2.0d0)
      TRUEANOM = atan(TANTD2) * 2.0d0

            ! Now decide which orbit is wanted and do the calculation

      if ( WHICHSTAR == 'A' ) then
        KVAL = KA
        VSYS = VA
      else if ( WHICHSTAR == 'B' ) then
        KVAL = 0.0d0 - KB
        VSYS = VB
      else if ( WHICHSTAR == 'C' ) then
        KVAL = 0.0d0 - KB
        VSYS = VA
      else
        write (6,'(A39,A41)') "### ERROR: the value of WHICHSTAR passe",
     &                      "d to functon GETRV is not 'A', 'B' or 'C'"
        STOP
      end if

      GETRV = KVAL * (cos(TRUEANOM + OMEGA/DEG2RAD) + ECOSW) + VSYS

      END FUNCTION GETRV
!=======================================================================
      DOUBLEPRECISION FUNCTION GETMIN (TZERO,PERIOD,ECCIN,OMEGAIN,CICLE)
            ! Returns time of minimum for given cycle and ephemeris.  If
            ! the orbit's circular then the cycle number is used without
            ! restriction so can refer to any phase. If the orbit is ec-
            ! centric then the cycle number should be integer (indicates
            ! primary minimum) or half-integer (secondary minimum).
      implicit none
      real*8 TZERO,PERIOD           ! IN: reference time, orbital period
      real*8 ECCIN,OMEGAIN          ! IN: orbital (e,w) or (ecosw,esinw)
      real*8 CICLE                  ! IN: cycle number of minimum to use
      real*8 CICLEFRAC              ! LOCAL: fraction part of cycle nmbr
      real*8 ECC,OMEGA              ! LOCAL: eccentricity and peri.long.
      real*8 ECOSW,ESINW            ! LOCAL: eccentr'y combination terms
      real*8 PSEP,PHASES(6)         ! LOCAL: phase sep and useful phases
      real*8 PI,DEG2RAD             ! LOCAL: useful variables

      PI = atan(1.0d0) * 4.0d0
      DEG2RAD = 45.0d0 / atan(1.0d0)

            ! First must deal with the possibility that e and omega are
            ! actually e*cos(omega) and e*sin(omega)

      CALL GETEOMEGA (ECCIN,OMEGAIN,ECC,OMEGA,ECOSW,ESINW)

            ! If orbit is circular then simply use the orbital ephemeris
            ! If orbit is eccentric then call ECQUADPHASES to calculate
            ! the phase difference between primary and secondary minima.

      if ( abs(ECC) < 1.0d-7 ) then
        GETMIN = TZERO  +  PERIOD * CICLE
      else
        CICLEFRAC = mod(CICLE,1.0d0)

        if ( ( CICLEFRAC >= 0.0d0 .and. CICLEFRAC < 0.001d0 ) .or.
     &       ( CICLEFRAC > 0.999d0 .and. CICLEFRAC <= 1.0d0 ) ) then
          GETMIN = TZERO + PERIOD * CICLE
        else if ((CICLEFRAC > -0.501d0 .and. CICLEFRAC < -0.499d0) .or.
     &           (CICLEFRAC > 0.499d0 .and. CICLEFRAC < 0.501d0) ) then
          CALL ECQUADPHASES (ECCIN,OMEGAIN,0.0d0,PHASES)
          PSEP = PHASES(3) - PHASES(1)
          if ( PSEP < 0.0d0 ) PSEP = PSEP + 1.0d0
          GETMIN = TZERO + PERIOD * (CICLE-0.50d0+PSEP)
        else
          write(6,'(A37,A43)') "### ERROR: found a cycle number which",
     &                    " is not integer or half-integer. Abort.    "
          write(6,'(A18,F20.10)') "### Cycle number =",CICLE
          write(6,*) " "
          stop
        end if

      end if

      END FUNCTION GETMIN
!=======================================================================
      SUBROUTINE GETEOMEGA (EIN,OMEGAIN,ECC,OMEGA,ECOSW,ESINW)
            ! Calculates eccentricity and omega, and combinations ecosw
            ! and esinw from input which could be (e,w) or (ecosw,esinw)
      implicit none
      real*8 EIN,OMEGAIN            ! IN: could be  e,w  or  ecosw,esinw
      real*8 ECC                    ! OUT: orbital eccentricity
      real*8 OMEGA                  ! OUT: periastron longitude (degree)
      real*8 ECOSW,ESINW            ! OUT: e*cos(omega) and e*sin(omega)
      real*8 DEG2RAD                ! LOCAL: convert degrees  to radians

      DEG2RAD = 45.0d0 / atan(1.0d0)

      if ( EIN > 5.0d0 ) then               ! if input is e+10 and omega
        ECC = EIN - 10.0d0
        OMEGA = OMEGAIN
        ECOSW = (EIN-10.0d0) * cos(OMEGAIN/DEG2RAD)
        ESINW = (EIN-10.0d0) * sin(OMEGAIN/DEG2RAD)

      else                                 ! if input is ecosw and esinw
        ECC = sqrt(EIN**2 + OMEGAIN**2)
        OMEGA = atan2(OMEGAIN,EIN) * DEG2RAD
        if ( OMEGA < 0.0d0 ) OMEGA = OMEGA + 360.0d0
        if ( OMEGA > 360.0d0 ) OMEGA = OMEGA - 360.0d0
        ECOSW = EIN
        ESINW = OMEGAIN
      end if

      END SUBROUTINE GETEOMEGA
!=======================================================================
!=======================================================================
      SUBROUTINE FITEBOP (DATX,DATY,DATERR,DTYPE,NDATA,MAXDATA,NLR,NMIN,
     &                    V,VARY,LDTYPE,ITER,CHISQ,VERR,ERROR,NSINE,
     &                 PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
            ! This subroutine calls the  MRQMIN  algorithm (Press et al,
            ! 1992, Numerical recipes in FORTRAN 77, p.678)  which finds
            ! the best-fitting EBOP model for the data using the
            ! Levenberg-Marquardt optimisation method.
            ! Unfortunately I have had to use a COMMON block here to
            ! avoid passing parameters through MRQMIN and related code.
      implicit none
      integer MAXDATA                     ! IN: max number of datapoints
      real*8 DATX(MAXDATA)                ! IN: Data independent varible
      real*8 DATY(MAXDATA)                ! IN: Data dependent variables
      real*8 DATERR(MAXDATA)              ! IN: Data errorbars
      integer DTYPE(MAXDATA)              ! IN: Types of datapoints
      integer NDATA,NLR,NMIN              ! IN: Numbers of datapoints
      real*8 V(540)                       ! IN: EBOP parameter values
      integer VARY(540)                   ! IN: Whether parameters vary
      integer LDTYPE(2)                   ! IN: Type of LD law for stars
      integer NL3,NECW,NESW               ! IN: Number of L3,ecosw,esinw
      integer NSINE,PSINE(20)             ! IN: number sines, which pars
      integer NPOLY,PPOLY(50)             ! IN: number polys, which pars
      integer NUMINT                      ! IN: Number of numerical ints
      real*8 NINTERVAL                    ! IN: Time interval for numint
      integer ITER                        ! OUT: Number of iterations
      real*8 CHISQ                        ! OUT: Reduced chi-squared
      real*8 VERR(540)                    ! OUT: Formal parameter errors
      integer ERROR                       ! OUT: Whether fit successful
      real*8 SIG(MAXDATA)                 ! SUB: abs(DATERR)
      real*8 LAMBDA                       ! LOCAL: Marquardt lambda
      real*8 COVAR(540,540),ALPHA(540,540)! LOCAL: Covariance etc matrix
      real*8 OCHISQ                       ! LOCAL: Previous chi-squared
      integer i,NVARY                     ! LOCAL: helpful integers
      integer CLDTYPE(2)                  ! COMMON: for EBOP subroutine
      integer CNSINE,CPSINE(20)           ! COMMON: for EBOP subroutine
      integer CNPOLY,CPPOLY(50)           ! COMMON: for EBOP subroutine
      real*8 CNINTERVAL                   ! COMMON: for EBOP subroutine
      integer CNUMINT                     ! COMMON: for EBOP subroutine

      common / FOREBOP / CNINTERVAL,CLDTYPE,CNSINE,CPSINE,CNPOLY,CPPOLY,
     &                                                           CNUMINT
      ERROR = 0

      CLDTYPE(1) = LDTYPE(1)
      CLDTYPE(2) = LDTYPE(2)
      CNSINE = NSINE
      CNPOLY = NPOLY
      do i = 1,20
        CPSINE(i) = PSINE(i)
      end do
      do i = 1,50
        CPPOLY(i) = PPOLY(i)
      end do
      CNUMINT = NUMINT
      CNINTERVAL = NINTERVAL

      do i = 1,NDATA
        SIG(i) = abs(DATERR(i))
      end do

      NVARY = 0
      do i = 1,540
        if ( VARY(i) == 1 .or. VARY(i) == 3 ) NVARY = NVARY + 1
      end do

            ! Now find the best fit using MRQMIN. This requires an init-
            ! ial call with LAMBDA less than zero to initialise things.
            ! Then iterate to the best fit and assume this has been got
            ! once LAMBDA > 10^7. If no observational errors have been
            ! supplied then calculate them and send in once more. And
            ! finally set LAMBDA = 0.0 to get useful parameters out.

      if ( NVARY == 0 ) then
        ITER = 0
      else

        LAMBDA = -1.0d0
        OCHISQ = 1.0d10
        CALL MRQMIN (DATX,DATY,SIG,DTYPE,NDATA,V,VARY,540,COVAR,ALPHA,
     &                                           540,CHISQ,LAMBDA,ERROR)
        do i = 1,199
          CALL MRQMIN (DATX,DATY,SIG,DTYPE,NDATA,V,VARY,540,COVAR,ALPHA,
     &                                           540,CHISQ,LAMBDA,ERROR)
          if ( LAMBDA >= 1.0d10 .and. abs(OCHISQ/CHISQ) < 1.001 ) exit
          if ( ERROR /= 0 ) exit
          OCHISQ = CHISQ
        end do
        ITER = i

      end if

      LAMBDA = 0.0d0
      CALL MRQMIN (DATX,DATY,SIG,DTYPE,NDATA,V,VARY,540,COVAR,ALPHA,540,
     &                                               CHISQ,LAMBDA,ERROR)

            ! Now record the formal errors outputted by MRQMIN

      do i = 1,540
        VERR(i) = sqrt(COVAR(i,i))
        if ( DATERR(1) < 0.0d0 )
     &            VERR(i) =  VERR(i) * sqrt( CHISQ / dble(NDATA-NVARY) )
      end do

      if ( V(6) > 90.0d0 ) V(6) = 180.0d0 - V(6)

      if ( LDTYPE(2) == -1 ) then
        V(5) = V(4)
        VERR(5) = VERR(4)
        V(24) = V(21)
        VERR(24) = VERR(21)
      end if

      if ( DATERR(1) > 0.0d0 )  CHISQ = CHISQ / dble(NDATA-NVARY)
      if ( DATERR(1) <= 0.0d0 ) CHISQ = -1.0d0

      if ( VARY(30) == -1 ) then
        V(30) = V(29)
        VERR(30) = VERR(29)
      end if

      END SUBROUTINE FITEBOP
!=======================================================================
      SUBROUTINE EBOP (DTYPE1,X,V,Y,DYDA,NCOEFFS,VARY,DODYDA)
            ! This evaluates the model value (Y) for one datapoint (X)
            ! INDEX is the datapoint number and DTYPE1 is its type.
            ! Optionally (if DODYDA = 'y') it also calculates the numer-
            ! ical derivatives of X with respect to the variable params.
      implicit none
      integer DTYPE1                ! IN: type of the datapoint
      real*8 X                      ! IN: Time to  calculate results for
      real*8 V(540)                 ! IN: The   photometric   parameters
      integer NCOEFFS               ! IN: Total number of adjustd params
      integer VARY(540)             ! IN:  Which params  being  adjusted
      character*1 DODYDA            ! IN: 'y' or 'n' todothe derivatives
      real*8 Y                      ! OUT:  Output result for input time
      real*8 DYDA(540)              ! OUT:   The numerical differentials
      real*8 LP,LS                  ! LOCAL: Light produced by each star
      real*8 OUT1,OUT2              ! LOCAL: Help in finding derivatives
      real*8 STORE                  ! LOCAL: Help in finding derivatives
      real*8 DV(540)                ! LOCAL: Amount to perturb params by
      real*8 R1,R2                  ! LOCAL: fractional radii
      integer i                     ! LOCAL: Loop counter
      real*8 GETMODEL               ! FUNCTION: Returns model evaluation
      integer LDTYPE(2)             ! IN/COMMON: LD law types  for stars
      integer NSINE,PSINE(20)       ! IN/COMMON: numbr sines, which pars
      integer NPOLY,PPOLY(50)       ! IN/COMMON: numbr polys, which pars
      integer NUMINT                ! IN/COMMON: numbr of numerical ints
      real*8 NINTERVAL              ! IN/COMMON: timeinterval for numint

      common / FOREBOP / NINTERVAL,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,NUMINT

      CALL GET_DV(V,DV,NPOLY,PPOLY)

            ! First get the model prediction for this datapoint. And use
            ! this call to GETMODEL to get the reflection coefficients

      if ( V(2) >= 0.0d0 ) then
        R1 = V(2) / (1.0d0 + V(3))
        R2 = V(2) / (1.0d0 + (1.0d0/V(3)))
      else
        R1 = abs(V(2))
        R2 = V(3)
      end if

      Y = GETMODEL (V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,X,DTYPE1,
     &                                           LP,LS,NUMINT,NINTERVAL)
      if ( VARY(11) == -1 )   V(11) = 0.4d0 * (LS/(1.0d0-V(15))) * R1**2
      if ( VARY(12) == -1 )   V(12) = 0.4d0 * (LP/(1.0d0-V(15))) * R2**2

            ! Now for each adjustable parameter work out the adjustment
            ! to make to its value for calculating partial derivatives.
            ! NOTE: for the sine reference time the size ofthe numerical
            ! interval depends on the period of the sine: this must be
            ! treated separately to avoid numerical intervals which are
            ! either too large or too small.

      if ( DODYDA == 'y' ) then

        do i = 1,NCOEFFS
          if ( VARY(i) == 1 .or. VARY(i) == 3 ) then

            STORE = V(i)
            V(i) = STORE + DV(i)
            OUT1 = GETMODEL (V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,X,
     &                                    DTYPE1,LP,LS,NUMINT,NINTERVAL)
            V(i) = STORE - DV(i)
            OUT2 = GETMODEL (V,VARY,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,X,
     &                                    DTYPE1,LP,LS,NUMINT,NINTERVAL)
            V(i) = STORE
            DYDA(i) = (OUT1 - OUT2) / (2.0d0 * DV(i))
          else
            DYDA(i) = 0.0d0
          end if
        end do

        if ( V(6) > 89.9d0) then
          DYDA(6) = DYDA(6) + (V(6)-89.9d0)**2
        end if

      else
        do i = 1,NCOEFFS
          DYDA(i) = 0.0d0
        end do
      end if

      END SUBROUTINE EBOP
!=======================================================================
!=======================================================================

!=======================================================================
!================     SUBROUTINES ORIGINALLY FROM EBOP    ==============
!=======================================================================
      SUBROUTINE BIAX (R,Q,A,B,EPS)
            ! EBOP subroutine to calculate biaxial ellipsoid dimensions
            ! and oblateness for each star after Chandrasekhar (1933).
      implicit none
      real*8 R,Q,A,B,EPS

      if ( Q <= 0.0d0 )  then
        A = R
        B = R
        EPS = 0.0d0
      else
        A = R * ( 1.0d0 + (1.0d0 + 7.0d0*Q)/6.0d0 * R**3.0d0)
        B = R * ( 1.0d0 + (1.0d0 - 2.0d0*Q)/6.0d0 * R**3.0d0)
        EPS = (A - B) / A
        B = ( (1.0d0 - EPS) * R**3.0d0) ** (1.0d0/3.0d0)
        A = B / (1.0d0 - EPS)
      end if

      END SUBROUTINE BIAX
!=======================================================================
      SUBROUTINE GETEW (ECOSW,ESINW,E,W)
            ! EBOP subroutine to calculate e and w from e(cos)w e(sin)w
      implicit none
      real*8 ECOSW,ESINW,E,W

      if ( ECOSW == 0.0d0  .and.  ESINW == 0.0d0 ) then
        E = 0.0d0
        W = 0.0d0
      else
        W = atan2( ESINW,ECOSW )
        E = sqrt( ESINW*ESINW + ECOSW*ECOSW )
        W = W * 180.0d0 / 3.1415926536d0
      end if

      END SUBROUTINE GETEW
!=======================================================================
      SUBROUTINE LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,HJD,FMAG,LP,LS)
            ! This is the subroutine where a single brightness is calcu-
            ! lated for a specific time given a set of system parametrs.
            ! The EBOP code (Etzel, 1981psbs.conf..111E) implements the
            ! NDE model (Nelson & Davis, 1972ApJ...174..617N) to do so.
            ! The LIGHT, BIAX and GETEW subroutines are the only pieces
            ! of code from EBOP that still remains in JKTEBOP. They are
            ! used with permission of Dr Paul B. Etzel (see Southworth,
            ! Smalley, Maxted, Claret & Etzel, 2005MNRAS.363..529S).
      implicit real*8 (a-h,o-z)
      real*8 V(540),HJD,GETPHASE
      real*8 LP,LS,LECL,LE
      real*8 LD1U,LD2U            ! linear LD coeff for each star
      real*8 LD1Q,LD2Q            ! quadratic LD coeff for each star
      real*8 LD1L,LD2L            ! logarithmic LD coeff for each star
      real*8 LD1S,LD2S            ! square-root LD coeff for each star
      real*8 LD1C,LD2C            ! cubic LD coeff for each star
      real*8 LD1POW1,LD2POW1      ! power-2 LD coeff for each star
      real*8 LD1POW2,LD2POW2      ! power-2 LD coeff for each star
      real*8 LD11,LD21            ! 4-par LD coeff 1 for each star
      real*8 LD12,LD22            ! 4-par LD coeff 2 for each star
      real*8 LD13,LD23            ! 4-par LD coeff 3 for each star
      real*8 LD14,LD24            ! 4-par LD coeff 4 for each star
      real*8 LDU,LDQ,LDS,LDL,LDC  ! LD coeffs for the star in question
      real*8 LDPOW1,LDPOW2        ! LD coeffs for the star in question
      real*8 LD1,LD2,LD3,LD4      ! LD coeffs for the star in question
      real*8 LD1UEQ,LD2UEQ        ! equivalent linear LD coefficients
      integer LDTYPE(2)           ! LD law type for both stars
      integer NSINE,PSINE(20)             ! IN: number sines, which pars
      integer NPOLY,PPOLY(50)             ! IN: number polys, which pars
      real*8 PHASE,SINT,SINP,SINA,SINTERM,POLYTERM,LPMULT,LSMULT
      real*8 ECC,OMEGA,ECOSW,ESINW        ! LOCAL: eccentricity values
      real*8 PERIOD,TZERO                 ! Orbital period and ref time

      DEG2RAD = 45.0d0 / atan(1.0d0)

            ! First check if we are using Mandel & Agol instead.

      if ( V(18) > -1.01d0 .and. V(18) < -0.99d0 ) then
        CALL OCCULTSMALL(V,HJD,FMAG)
        LP = 1.0d0
        LS = 0.0d0
!         write(71,*)fmag
        return
      end if

            ! Now get on with the interface to the EBOP model.

      LPMULT = 1.0d0
      LSMULT = 1.0d0
      PI = 3.1415926536d0
      TWOPI = 6.28318531d0
      RAD = 0.0174532925d0

C
C        DETERMINE PRIMARY AND SECONDARY BIAXIAL DIMENSIONS
C        USE SPHERICAL RADII FOR THE COMPUTATION OF ECLIPSE FUNCTIONS
C        USE OBLATENESSES FOR THE COMPUTATION OF THE OUTSIDE ECLIPSE
C        PHOTOMETRIC VARIATIONS WITH LIMB AND GRAVITY DARKENING
C

      if ( V(2) >= 0.0d0 ) then
        RP = V(2) / ( 1.0d0 + V(3) )
        RS = V(2) / ( 1.0d0 + (1.0d0/V(3)) )
      else
        RP = abs( V(2) )
        RS = V(3)
      end if

      CALL GETEOMEGA (V(7),V(8),ECC,OMEGA,ECOSW,ESINW)

      BS     = V(1)
      FI     = V(6)
      YP     = V(9)
      YS     = V(10)
      SP     = V(11)
      SS     = V(12)
      Q      = V(13)
      TANGL  = 0.0d0!V(14)
      EL     = V(15)
      DPH    = 1.0d0 - V(16)
      SFACT  = V(17)
      DGAM   = V(18)
      PERIOD = V(19)
      TZERO  = V(20)

            ! Now put the fractional radii into RP and RS.
            ! There are two options here: the default that the numbers
            ! supplied are actually (r1+r2) and k, and the alternative
            ! (indicated by a negative r1+r2) that they are -r1 and r2

      if ( V(2) >= 0.0d0 ) then
        RP = V(2) / ( 1.0d0 + V(3) )
        RS = V(2) / ( 1.0d0 + (1.0d0/V(3)) )
      else
        RP = abs( V(2) )
        RS = V(3)
      end if

            ! Limb darkening now.
            ! First set all the possible LD coefficients to zero

      LD1U = 0.0d0              ! linear terms
      LD2U = 0.0d0
      LD1L = 0.0d0              ! log terms
      LD2L = 0.0d0
      LD1S = 0.0d0              ! sqrt terms
      LD2S = 0.0d0
      LD1Q = 0.0d0              ! quadratic terms
      LD2Q = 0.0d0
      LD1C = 0.0d0              ! cubic terms
      LD2C = 0.0d0

      LD1POW1 = 0.0d0           ! power-2 law coefficients
      LD1POW2 = 0.0d0
      LD2POW1 = 0.0d0
      LD2POW2 = 0.0d0

      LD11 = 0.0d0              ! Claret 4-par law ^(1/2) terms
      LD21 = 0.0d0
      LD12 = 0.0d0              ! Claret 4-par law ^1 terms
      LD22 = 0.0d0
      LD13 = 0.0d0              ! Claret 4-par law ^(3/2) terms
      LD23 = 0.0d0
      LD14 = 0.0d0              ! Claret 4-par law ^2 terms
      LD24 = 0.0d0

            ! Now insert the LD coefficients into the correct variables
            ! for star A. If we are using the re-parameterisations then
            ! we need to do the **opposite*( transformation to get back
            ! to the original LD coefficients as used below.

      if ( LDTYPE(1) == 1 ) then
        LD1U = V(4)
      else if ( LDTYPE(1) == 2 ) then
        LD1U = V(4)
        LD1Q = V(21)
      else if ( LDTYPE(1) == 3 ) then
        LD1U = V(4)
        LD1L = V(21)
      else if ( LDTYPE(1) == 4 ) then
        LD1U = V(4)
        LD1S = V(21)
      else if ( LDTYPE(1) == 5 ) then
        LD1U = V(4)
        LD1C = V(21)
      else if ( LDTYPE(1) == 6 ) then
        LD1POW1 = V(4)
        LD1POW2 = V(21)
      else if ( LDTYPE(1) == 9 ) then
        LD11 = V(4)
        LD12 = V(21)
        LD13 = V(22)
        LD14 = V(23)
      else if ( LDTYPE(1) == -2 ) then
        LD1U = (V(4) + V(21)) / 2.0d0
        LD1Q = (V(4) - V(21)) / 2.0d0 / 0.6d0
      else if ( LDTYPE(1) == -3 ) then
        LD1U = (V(4) + V(21)) / 2.0d0
        LD1L = (V(4) - V(21)) / 2.0d0 / 0.6d0
      else if ( LDTYPE(1) == -4 ) then
        LD1U = (V(4) + V(21)) / 2.0d0
        LD1S = (V(4) - V(21)) / 2.0d0 / 0.6d0
      else if ( LDTYPE(1) == -5 ) then
        LD1U = (V(4) + V(21)) / 2.0d0
        LD1S = (V(4) - V(21)) / 2.0d0 / 0.3d0
      else if ( LDTYPE(1) == -6 ) then
        LD1POW1 = 1.0d0 - V(4) + V(21)
        LD1POW2 = log10(LD1POW1/V(21)) / log10(2.0d0)    ! to get log_2
      else
        write(6,*) " "
        write(6,*) "### ERROR: wrong LDTYPE(1) in subroutine LIGHT."
        write(6,*) " "
        stop
      endif

            ! Now insert the LD coefficients into the correct variables
            ! for star B. See the description of this above for star A.
            ! The differences are the variable names and array indices,
            ! and that star B can be fixed to the same LDCs as star A.

      if ( LDTYPE(2) == 1 ) then
        LD2U = V(5)
      else if ( LDTYPE(2) == 2 ) then
        LD2U = V(5)
        LD2Q = V(24)
      else if ( LDTYPE(2) == 3 ) then
        LD2U = V(5)
        LD2L = V(24)
      else if ( LDTYPE(2) == 4 ) then
        LD2U = V(5)
        LD2S = V(24)
      else if ( LDTYPE(2) == 5 ) then
        LD2U = V(5)
        LD2C = V(24)
      else if ( LDTYPE(2) == 6 ) then
        LD2POW1 = V(5)
        LD2POW2 = V(24)
      else if ( LDTYPE(2) == 9 ) then
        LD21 = V(5)
        LD22 = V(24)
        LD23 = V(25)
        LD24 = V(26)
      else if ( LDTYPE(2) == -2 ) then
        LD2U = (V(5) + V(24)) / 2.0d0
        LD2Q = (V(5) - V(24)) / 2.0d0 / 0.6d0
      else if ( LDTYPE(2) == -3 ) then
        LD2U = (V(5) + V(24)) / 2.0d0
        LD2L = (V(5) - V(24)) / 2.0d0 / 0.6d0
      else if ( LDTYPE(2) == -4 ) then
        LD2U = (V(5) + V(24)) / 2.0d0
        LD2S = (V(5) - V(24)) / 2.0d0 / 0.6d0
      else if ( LDTYPE(2) == -5 ) then
        LD2U = (V(5) + V(24)) / 2.0d0
        LD2S = (V(5) - V(24)) / 2.0d0 / 0.3d0
      else if ( LDTYPE(2) == -6 ) then
        LD2POW1 = 1.0d0 - V(5) + V(24)
        LD2POW2 = log10(LD2POW1/V(24)) / log10(2.0d0)    ! to get log_2
      else if ( LDTYPE(2) == -1 ) then
        LD2U    = LD1U
        LD2Q    = LD1Q
        LD2L    = LD1L
        LD2S    = LD1S
        LD2C    = LD1C
        LD2POW1 = LD1POW1
        LD2POW2 = LD1POW2
        LD21    = LD11
        LD22    = LD12
        LD23    = LD13
        LD24    = LD14
      else
        write(6,*) " "
        write(6,*) "### ERROR: wrong LDTYPE(2) in subroutine LIGHT."
        write(6,*) " "
        stop
      endif

!        print*," "
!        print*,"RP   ",RP
!        print*,"RS   ",RS
!        print*,"ecosw",ECOSW
!        print*,"esinw",ESINW
!        print*,"J    ",BS
!        print*,"i    ",FI
!        print*,"GD1  ",YP
!        print*,"GD2  ",YS
!        print*,"refl1",SP
!        print*,"refl2",SS
!        print*,"Q    ",Q
!        print*,"TANGL",TANGL
!        print*,"L3   ",EL
!        print*,"Pshft",DPH
!        print*,"SFACT",SFACT
!        print*,"intrg",DGAM
!        print*,"Porb ",PERIOD
!        print*,"Tpri ",TZERO
!        print*,"LD1U    ",LD1U
!        print*,"LD2U    ",LD2U
!        print*,"LD1L    ",LD1L
!        print*,"LD2L    ",LD2L
!        print*,"LD1S    ",LD1S
!        print*,"LD2S    ",LD2S
!        print*,"LD1Q    ",LD1Q
!        print*,"LD2Q    ",LD2Q
!        print*,"LD1C    ",LD1C
!        print*,"LD2C    ",LD2C
!        print*,"LD11    ",LD11
!        print*,"LD21    ",LD21
!        print*,"LD12    ",LD12
!        print*,"LD22    ",LD22
!        print*,"LD13    ",LD13
!        print*,"LD23    ",LD23
!        print*,"LD14    ",LD14
!        print*,"LD24    ",LD24
!        print*,"LD1POW1 ",LD1POW1
!        print*,"LD1POW2 ",LD1POW2
!        print*,"LD1H1   ",LD1H1
!        print*,"LD1H2   ",LD1H2
!        print*,"LD2POW1 ",LD2POW1
!        print*,"LD2POW2 ",LD2POW2
!        print*,"LD2H1   ",LD2H1
!        print*,"LD2H2   ",LD2H2
!        stop

      if ( NSINE > 0 ) then
        do i = 1,NSINE
          SINT = V(31+(i-1)*3)      ! sine reference time
          SINP = V(32+(i-1)*3)      ! sine period
          SINA = V(33+(i-1)*3)      ! sine amplitude
          SINTERM = SINA * sin( TWOPI * (HJD-SINT) / SINP )
          if ( PSINE(i) ==  1 )  BS = BS * (1.0d0+SINTERM)
          if ( PSINE(i) ==  2 )  RP = RP * (1.0d0+SINTERM)
          if ( PSINE(i) ==  3 )  RS = RS * (1.0d0+SINTERM)
          if ( PSINE(i) ==  6 )  FI = FI + SINTERM
          if ( PSINE(i) == 15 )  EL = EL * (1.0d0+SINTERM)
          if ( PSINE(i) == 17 )  SFACT = SFACT + SINTERM
          if ( PSINE(i) == 20 )  TZERO = TZERO + SINTERM
          if ( PSINE(i) == -1 )  LPMULT = LPMULT * (1.0d0+SINTERM)
          if ( PSINE(i) == -2 )  LSMULT = LSMULT * (1.0d0+SINTERM)

          if ( PSINE(i) == 7 .or. PSINE(i) == 8 ) then
            if ( V(7) >= 9.99d0 ) then
              if ( PSINE(i) == 7 ) ECC = ECC + SINTERM
              if ( PSINE(i) == 8 ) OMEGA = OMEGA + SINTERM
              ECOSW = ECC * cos(OMEGA/DEG2RAD)
              ESINW = ECC * sin(OMEGA/DEG2RAD)
            endif
            if ( V(7) < 9.99d0 ) then
              if ( PSINE(i) == 7 ) ECOSW = ECOSW + SINTERM
              if ( PSINE(i) == 8 ) ESINW = ESINW + SINTERM
              CALL GETEW (ECOSW,ESINW,ECC,OMEGA)
            endif
          endif
        end do
      end if

            ! Apply the polynomial. First find the V array elements
            ! containing the polynomial quantities. Then calculate
            ! the size of the polynomial evaluated at the time of the
            ! datapoint. Add it to the appropriate parameter ONLY if
            ! the HJD is within the interval given by V(j+6) to V(j+7).

      if ( NPOLY > 0 ) then
        do i = 1,NPOLY
          j = 91 + 9*(i-1)
          POLYTERM = V(j+3) + V(j+4)*(HJD-V(j))
     &                      + V(j+5)*((HJD-V(j))**2)
     &                      + V(j+6)*((HJD-V(j))**3)
     &                      + V(j+7)*((HJD-V(j))**4)
     &                      + V(j+8)*((HJD-V(j))**5)
          if ( HJD < V(j+1) ) POLYTERM = 0.0d0
          if ( HJD > V(j+2) ) POLYTERM = 0.0d0
          if ( PPOLY(i) ==  1 )  BS = BS + POLYTERM
          if ( PPOLY(i) ==  2 )  RP = RP + POLYTERM
          if ( PPOLY(i) ==  3 )  RS = RS + POLYTERM
          if ( PPOLY(i) ==  6 )  FI = FI + POLYTERM
          if ( PPOLY(i) == 15 )  EL = EL + POLYTERM
          if ( PPOLY(i) == 17 )  SFACT = SFACT + POLYTERM
          if ( PPOLY(i) == 20 )  TZERO = TZERO + POLYTERM
          if ( PPOLY(i) == -1 )  LPMULT = LPMULT + POLYTERM
          if ( PPOLY(i) == -2 )  LSMULT = LSMULT + POLYTERM

          if ( PPOLY(i) == 7 .or. PPOLY(i) == 8 ) then
            if ( V(7) >= 9.99d0 ) then
              if ( PPOLY(i) == 7 ) ECC = ECC + POLYTERM
              if ( PPOLY(i) == 8 ) OMEGA = OMEGA + POLYTERM
              ECOSW = ECC * cos(OMEGA/DEG2RAD)
              ESINW = ECC * sin(OMEGA/DEG2RAD)
            endif
            if ( V(7) < 9.99d0 ) then
              if ( PPOLY(i) == 7 ) ECOSW = ECOSW + POLYTERM
              if ( PPOLY(i) == 8 ) ESINW = ESINW + POLYTERM
              CALL GETEW (ECOSW,ESINW,ECC,OMEGA)
            endif
          endif
        end do
      end if

      PHASE = GETPHASE(HJD,PERIOD,TZERO)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !         RS=RP*RATIO
      if ( Q <= 0.0d0 ) then
        CALL BIAX (RP,0.0d0,RPA,RPB,EP)
        CALL BIAX (RS,0.0d0,RSA,RSB,ES)
      else
        CALL BIAX (RP,Q,RPA,RPB,EP)
        CALL BIAX (RS,1.0d0/Q,RSA,RSB,ES)
      end if

C
C        CORRECT THE OBSERVED PHASE FOR ANY EPOCH ERROR IN EPHEMERIS
C
      THETA=PHASE+DPH
C
      SINI  = SIN(FI*RAD)
      SINI2 = SINI*SINI
      COSI2 = 1.0d0  - SINI2
C
C        TRANSLATE TIDAL LEAD/LAG ANGLE TO RADIANS
      TANGR=TANGL*RAD
C
C     EQUATION 9
C        CONVERT PHASE TO RADIANS
      FMN=THETA*TWOPI
C
C        GET CURRENT VALUES OF E, AND W
      CALL GETEW (ECOSW,ESINW,E,W)
C
C        TEST FOR CIRCULAR ORBIT
      IF (E)   17,20,17
   20 COSVW=COS(FMN)
      SINVW=SIN(FMN)
      RV=1.0d0
      GO TO 25
C
C        SOLUTION OF KEPLER'S EQUATION BY DIFFERENTIAL CORRECTIONS
C        (NON-ZERO ECCENTRICITY ONLY . . . )
C
C     EQUATION 6
C
   17 OMEGA = 450.0d0  - W
   23 IF (OMEGA - 360.0d0)         22,21,21
   21 OMEGA = OMEGA - 360.0d0
      GO TO 23
   22 OMEGA = OMEGA*RAD
C        SINE AND COSINE OF OMEGA
      COSW=COS(OMEGA)
      SINW=SIN(OMEGA)
C
C        COMPUTE MEAN ANOMALY CORRECTION TO PHASE
C        CORRESPONDING TO V=OMEGA=90-W
C        AT WHICH PHASE COS(V-OMEGA)=1
      E0=ATAN2(SQRT(1.0d0-E*E)*SINW,COSW+E)
C
C        MEAN ANOMALY OF MID-PRIMARY ECLIPSE
      FMA0=E0-E*SIN(E0)
C
C        MEAN ANOMALY
      FMA=FMN+FMA0
C     FIRST APPROXIMATION OF ECCENTRIC ANOMALY
      EA=FMA+E*SIN(FMA)
C
      DO 10 J=1,15
C        EVALUATE SINE AND COSINE OF ECCENTRIC ANOMALY
      SINE=SIN(EA)
      COSE=COS(EA)
      DENOM=1.0d0-E*COSE
      DISC=FMA-EA+E*SINE
      EA=EA+DISC/DENOM
C        TEST FOR CONVERGENCE
      IF (ABS(DISC) - 2.0d-5)     15,15,10
   10 CONTINUE
C
C
C        EVALUATE SINE AND COSINE OF TRUE ANOMALY
   15 COSV=(COSE-E)/DENOM
      SINV=SINE*SQRT(1.0d0-E*E)/DENOM
C
C        RADIUS VECTOR
      RV = (1.0d0-E*E)/(1.0d0+E*COSV)
C
C        THE PHOTOMETRIC PHASE ARGUMENT IN TERMS OF ORBIT PARAMETERS
C        VW = V-OMEGA
      COSVW=COSV*COSW+SINV*SINW
      SINVW=SINV*COSW-COSV*SINW
C
   25 COS2=COSVW*COSVW
      SIN2=1.0d0-COS2
C
      CSVWT=COS(TANGR)*COSVW-SIN(TANGR)*SINVW
C
C
C        PHOTOMETRIC EFFECTS
C
C

!----------------------------------------------------------------------!
! Modifications for version 40, May 2020, by JKT                       !
! and modifications for version 42, March 2022, by JKT                 !
! and modifications for version 43, September 2022, by JKT             !
!----------------------------------------------------------------------!
! Now to sort out the limb darkening and gravity darkening normalisation
! plus the ellipsoidal effect. One problem is that I only have an equat-
! ion for flux normalisation with gravity darkening in the case of the
! linear LD law. Solution: convert all the more complex laws into an
! equivalent for the linear law that gives the same disc-integrated LD.
! Some of my derivations were checked against Heller (arXiv:1901.01730).

      FMAXP = 0.0d0           ! Disc-integrated flux for the primary
      FMAXS = 0.0d0           ! Disc-integrated flux for the secondary
      DELTP = 0.0d0           ! Change in brightness due to oblateness
      DELTS = 0.0d0           ! Change in brightness due to oblateness
      SHORT = 0.0d0           ! Foreshortening

            ! Remember that all unused LD coeffs will be set to zero
            ! so we can do all the LD laws in one equation.
            ! First do the calculations for spherical stars.

      if ( Q < 0.0d0 ) then
        FMAXP = 1.0d0 - LD1U/3.0d0 - LD1Q/6.0d0 - LD1C/10.0d0
     &                - LD1S/5.0d0 + LD1L*2.0d0/9.0d0 - LD11/5.0d0
     &                - LD12/3.0d0 - LD13*3.0d0/7.0d0 - LD14/2.0d0
     &                - LD1POW1 + 2.0d0*LD1POW1/(LD1POW2+2.0d0)

        FMAXS = 1.0d0 - LD2U/3.0d0 - LD2Q/6.0d0 - LD2C/10.0d0
     &                - LD2S/5.0d0 + LD2L*2.0d0/9.0d0 - LD21/5.0d0
     &                - LD22/3.0d0 - LD23*3.0d0/7.0d0 - LD24/2.0d0
     &                - LD2POW1 + 2.0d0*LD2POW1/(LD2POW2+2.0d0)
        DELTP = 0.0d0
        DELTS = 0.0d0
        SHORT = 0.0d0
      end if

            ! Now for distorted stars.
            ! EP and ES are the oblatenesses of the two stars from BIAX
            ! YP and YS are gravity darkening coefficients for the stars

      if ( Q >= 0.0d0 ) then

            ! LD1UEQ and LD2UEQ are LD coeffs for the linear law that
            ! give the same overall LD as for the nonlinear laws.

        LD1UEQ = LD1U + LD1Q*0.5d0 + LD1C*0.3d0 + LD1S*0.6d0 -LD1L/1.5d0
     &               + LD11*0.6d0 + LD12 + LD13*9.0d0/7.0d0 + LD14*1.5d0
     &               + 3.0d0*LD1POW1*(1.0d0-2.0d0/(LD1POW2+2.0d0))
        LD2UEQ = LD2U + LD2Q*0.5d0 + LD2C*0.3d0 + LD2S*0.6d0 -LD2L/1.5d0
     &               + LD21*0.6d0 + LD22 + LD23*9.0d0/7.0d0 + LD24*1.5d0
     &               + 3.0d0*LD2POW1*(1.0d0-2.0d0/(LD2POW2+2.0d0))

            ! Get normalised fluxes accounting for LD, GD and oblateness

        FMAXP = ((1.0d0-LD1UEQ)+0.66666666667d0*LD1UEQ*(1.0d0+0.2d0*EP))
     &                                   *(1.0d0+3.0d0*YP*EP)/(1.0d0-EP)
        FMAXS = ((1.0d0-LD2UEQ)+0.66666666667d0*LD2UEQ*(1.0d0+0.2d0*ES))
     &                                   *(1.0d0+3.0d0*YS*ES)/(1.0d0-ES)

            ! Get the total change in brightness due to oblateness (not
            ! dependent on time) and the amount of this (is dep on time)

        DELTP = (15.0d0+LD1UEQ) / (15.0d0-5.0d0*LD1UEQ) * (1.0d0+YP)*EP
        DELTS = (15.0d0+LD1UEQ) / (15.0d0-5.0d0*LD1UEQ) * (1.0d0+YS)*ES
        SHORT = SINI2 * CSVWT * CSVWT
!         print'(7(f10.6))',LD1UEQ,LD2UEQ,FMAXP,FMAXS,DELTP,DELTS,SHORT
      end if

!----------------------------------------------------------------------
! Complete original code replaced by the code above:
! C
! C        PHOTOMETRIC EFFECTS
! C
! C        TEST FOR SIMPLE CASE OF TWO SPHERICAL STARS
!       IF (EP .EQ. 0.  .AND.  ES .EQ. 0.)   GO TO 26
! C
! C        EITHER OR BOTH STARS ARE OBLATE
! C
!       FMAXP=((1.0E0-UP)+0.666666667E0*UP*(1.0E0+0.2E0*EP))
!      1      *(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)
!       FMAXS=((1.0E0-US)+0.666666667E0*US*(1.0E0+0.2E0*ES))
!      1      *(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)
! C        CHANGE IN INTENSITY RATIO DUE TO OBLATENESS RELATED VARIABLES
! C        FROM QUADRATURE TO MINIMUM
! C        FACE ON TO END ON
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES
! C        FORE-SHORTENING FUNCTION OF OBLATENESS
!       SHORT=SINI2*CSVWT*CSVWT
!       GO TO 27
! C
! C        BOTH STARS ARE SPHERICAL
! C
!    26 FMAXP=1.0E0-UP/3.0E0
!       FMAXS=1.0E0-US/3.0E0
!       DELTP=0.0E0
!       DELTS=0.0E0
!       SHORT=0.0
!----------------------------------------------------------------------

C
C        UN-NORMALIZED BRIGHTNESS OF STELLAR COMPONENTS AT QUADRATURE
   27 OP=PI*RPB*RPB*FMAXP
      OS=PI*RSB*RSB*FMAXS*BS
C        THE NORMALIZING FACTOR
      OTOT=OP+OS
C        BRIGHTNESS CONTRIBUTION FROM EACH COMPONENT
      LP=OP/OTOT*(1.0d0-DELTP*SHORT)
      LS=OS/OTOT*(1.0d0-DELTS*SHORT)
C
C        REFLECTION AND RERADIATION EQUATION
      IF (SP .EQ. 0.0d0  .AND.  SS .EQ. 0.0d0)   GO TO 28
      HEAT=SINI*COSVW
      HEAT2=0.5d0+0.5d0*HEAT*HEAT
      DLP=SP*(HEAT2+HEAT)
      DLS=SS*(HEAT2-HEAT)
      GO TO 29
   28 DLP=0.0d0
      DLS=0.0d0
C
C        WHICH ECLIPSE COULD THIS BE
   29 IF (COSVW)         40,40,30
C
C     PRIMARY ECLIPSE
C
   30 R1 = RP
      R2 = RS
!----------------------------------------------------------------------!
! JKT mod (2006/08/10): the line these replaced was      UU = UP       !
!----------------------------------------------------------------------!
      LDU = LD1U                                                       !
      LDQ = LD1Q                                                       !
      LDL = LD1L                                                       !
      LDS = LD1S                                                       !
      LDC = LD1C                                                       !
      LD1 = LD11                                                       !
      LD2 = LD12                                                       !
      LD3 = LD13                                                       !
      LD4 = LD14                                                       !
      LDPOW1 = LD1POW1                                                 ! modified 2022/03/22
      LDPOW2 = LD1POW2                                                 ! modified 2022/03/22
!----------------------------------------------------------------------!
      LE=LP
      DLE=DLP
      GO TO 60
C
C
C     SECONDARY ECLIPSE
C
   40 R1 = RS
      R2 = RP
!-----------------------------------------------------------------------
! JKT mod (2006/08/10): the line these replaced was      UU = US       !
!----------------------------------------------------------------------!
      LDU = LD2U                                                       !
      LDL = LD2L                                                       !
      LDS = LD2S                                                       !
      LDQ = LD2Q                                                       !
      LDC = LD2C                                                       !
      LD1 = LD21                                                       !
      LD2 = LD22                                                       !
      LD3 = LD23                                                       !
      LD4 = LD24                                                       !
      LDPOW1 = LD2POW1                                                 ! modified 2022/03/22
      LDPOW2 = LD2POW2                                                 ! modified 2022/03/22
!----------------------------------------------------------------------!
      LE=LS
      DLE=DLS
C
   60 SUMS = 0.0d0
      ALAST = 0.0d0
      AREA=0.0d0
C
C     EQUATION  5
C
      DD = SINVW*SINVW + COSVW*COSVW*COSI2
!----------------------------------------------------------------------!
      IF (DD .LE. 1.0d-8)  DD=0.0d0
!----------------------------------------------------------------------!
! The limiting value in the statement above used to be 1.0d-6 but      !
! Willie Torres and I independently identified this as a problem for   !
! very long-period EBs. We found that the eclipse attained a weird     !
! boxy shape near the centre, and if the parameters were pushed        !
! further the eclipse vanished entirely. Fix: a new limiting value     !
! of 1.0d-8. The change was made on 2014/02/10 by JKT (JKTEBOP v33).   !
!----------------------------------------------------------------------!
      DD = DD*RV*RV
      D = SQRT(ABS(DD))
      R22 = R2*R2
C
C     EQUATION 17
C
      GAMN = 90.01d0*RAD
      DGAMA = DGAM*RAD
      DGM = DGAMA/2.0d0
      RK = 0.0d0
      GAM = 0.0d0
   50 GAM = GAM + DGAMA
C        HAS LIMIT OF INTEGRATION BEEN REACHED
      IF (GAM - GAMN)              48,48,49
C
   48 RR = R1*SIN(GAM)
      R12 = RR*RR
C
      AA = 0.0d0
C        ARE THE PROJECTED DISKS CONCENTRIC
      IF (D)                       405,406,405
  406 IF (RR - R2)                 230,230,403
  403 IF (RK - R2)                 404, 49, 49
  404 AA = PI*R22
      GO TO 215
C        TEST FOR NO ECLIPSE
  405 IF (D-R1-R2)                 240,216,216
  216 SUMS = 0.0d0
      GO TO 49
C        DECIDE WHICH AREA EQUATIONS FOR NON-CONCENTRIC ECLIPSE
  240 IF (D-RR-R2)                 245,215,215
  245 IF (D-R2+RR)                 230,230,250
  250 IF (R1-R2)                   255,255,280
  255 IF (DD-R22+R12)              205,210,210
  280 IF (D-RR+R2)                 290,260,260
  260 IF (RR-R2)                   255,255,265
  265 IF (DD-R12+R22)              270,210,210
C
C     EQUATION 12
C
  270 S1 = ABS((R12 - R22 - DD)*0.5d0/D)
      A1 = ABS(R2-S1)
      B2 = ABS(RR-S1-D  )
      AA=PI*R22-(R22*ACOS((R2-A1)/R2)
     1   - (R2-A1)*SQRT(2.0d0*R2*A1-A1*A1))
     2   +R12*ACOS((RR-B2)/RR)-(RR-B2)*SQRT(2.0d0*RR*B2-B2*B2)
      GO TO 215
C
  290 IF (R1 - R2 - D)             260,260,295
  295 IF (RK - R2 - D)             300,215,215
  300 RR = R2 + D
      R12 = RR*RR
      GAMN = 0.0d0
      GO TO 260
C
  230 AA = PI*R12
      GO TO 215
C
C     EQUATION 10
C
  205 S = ABS((R12 - R22 + DD)*0.5d0/D)
      A = ABS(RR-S)
      B1 = ABS(R2-S-D)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0d0*RR*A - A*A)
      AB1 = R22*ACOS((R2-B1)/R2) - (R2-B1)*SQRT(2.0d0*R2*B1-B1*B1)
      AA = PI*R12 - A1 + AB1
      GO TO 215
C
C     EQUATION 1
C
  210 S = ABS((R12 - R22 + DD)*0.5d0/D)
      A = ABS(RR-S)
      B = ABS(S-D+R2)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0d0*RR*A - A*A)
      AA1 = R22*ACOS((R2-B)/R2) - (R2-B)*SQRT(2.0d0*R2*B - B*B)
      AA = A1 + AA1
C
  215 DAREA = AA - ALAST
!----------------------------------------------------------------------!
! JKT modification (2006/09/10). The removed line was:                 !
!     SUM = SUM + DAREA*(1.0d0 - UU + UU*COS(GAM-DGM))                 !
!----------------------------------------------------------------------!
      COSGAM = cos(GAM-DGM)                                            !
      SUMS = SUMS + DAREA*(1.0d0 - LDU*(1.0d0-COSGAM)                  !
     &          - LDL*COSGAM*log(COSGAM) - LDS*(1.0d0-sqrt(COSGAM))    !
     &          - LDQ*(1.0d0-COSGAM)**2 - LDC*(1.0d0-COSGAM)**3        !
     &          - LD1*(1.0d0-COSGAM**0.5d0) - LD2*(1.0d0-COSGAM)       !
     &          - LD3*(1.0d0-COSGAM**1.5d0) - LD4*(1.0d0-COSGAM**2)    !
     &          - LDPOW1*(1.0d0-COSGAM**LDPOW2))                       ! JKT modification 2022/03/22
!       write(99,'(30(f15.8))')sums,COSGAM,darea,
!      &   COSGAM**LDPOW2,LDPOW1*(1.0d0-COSGAM**LDPOW2),
!      &    darea*(1.0-LDPOW1*(1.0d0-COSGAM**LDPOW2))
!       write(98,'(20(f10.5))')cosgam,darea,ldu,ldl,lds,ldq,ldc,ld1,ld2,
!      &     ld3,ld4,ldpow1,ldpow2
!       print*,"a",darea,SUMS,1.0d0 - LDU*(1.0d0-COSGAM)
!      &          - LDL*COSGAM*log(COSGAM) - LDS*(1.0d0-sqrt(COSGAM)
!      &          - LDQ*(1.0d0-COSGAM)**2 - LDC*(1.0d0-COSGAM)**3)
!      &          - LD1*(1.0d0-COSGAM**0.5d0) - LD2*(1.0d0-COSGAM)
!      &          - LD3*(1.0d0-COSGAM**1.5d0) - LD4*(1.0d0-COSGAM**2)
!       print*,"b",LD1*(1.0d0-COSGAM**0.5d0),LD2*(1.0d0-COSGAM),
!      &           LD3*(1.0d0-COSGAM**1.5d0),LD4*(1.0d0-COSGAM**2)
!       print*,"l",LDU,LDL,LDS,LDQ,LDC,LD1,LD2,LD3,LD4
!----------------------------------------------------------------------!
      ALAST = AA
      AREA = AREA + DAREA
C
      RK = RR
      GO TO 50
C
C        LIGHT LOSS FROM ECLIPSE
C
   49 ADISK = PI*R1*R1
!----------------------------------------------------------------------!
! JKT modification (2006/9/10).  See 1992A+A...259..227D for more info.!
! The removed line was:          ALPHA = SUM/(ADISK*(1.0-UU/3.0))      !
!----------------------------------------------------------------------!
      ALPHA = 1.0d0 - LDU/3.0d0 + LDL*2.0d0/9.0d0                      !
     &        - LDS/5.0d0 - LDQ/6.0d0 - LDC/10.0d0                     !
     &        - LD11/5.0d0 - LD12/3.0d0 - LD13*3.0/7.0d0 - LD14/2.0d0  !
     &        - LDPOW1 + 2.0d0*LDPOW1/(LDPOW2+2.0d0)                   ! JKT modification 2022/03/22
      ALPHA = SUMS/(ADISK*ALPHA)                                       !
!----------------------------------------------------------------------!
      LECL = ALPHA*LE
      AREA = AREA/ADISK
      REFL=DLP+DLS-AREA*DLE
C
C        THEORETICAL INTENSITY WITH THIRD LIGHT AND QUADRATURE
C        SCALE FACTOR APPLIED
C
!----------------------------------------------------------------------!
! This is the original line from EBOP:
!----------------------------------------------------------------------!
!      FLITE = ((LP+LS-LECL+REFL)*(1.0d0-EL)+EL)*SFACT
!----------------------------------------------------------------------!

      LP = LP * LPMULT               ! sine/poly applied to star A light
      LS = LS * LSMULT               ! sine/poly applied to star B light
      FLITE = ((LP+LS-LECL+REFL)*(1.0d0-EL)+EL)
      FMAG = -2.5d0 * log10(FLITE) + SFACT

      LP = LP * (1.0d0-EL)           ! account for third light *AFTER*
      LS = LS * (1.0d0-EL)           ! FLITE and FMAG have been found

      END
!=======================================================================
!=======================================================================

!=======================================================================
!=======================================================================
      SUBROUTINE OCCULTSMALL(V,HJD,MAG)
      implicit none
            ! This subroutine calculates the flux expected for a transit
            ! of a planet across a star, using the small-planet approxi-
            ! mation and the Claret four-parameterlimb darkening law.
            ! Reference: Mandel & Agol (2002ApJ...580L.171M)
            ! Taken from the website of Eric Agol on 2013/06/18
            ! http://www.astro.washington.edu/users/agol/

      real*8 V(540)                  ! IN: EBOP model parameters
      real*8 HJD                    ! IN: time of measurement
      real*8 MAG                    ! OUT: predicted brightness in mags
      real*8 ECC,OMEGA,ECOSW,ESINW  ! LOCAL: orbital shape parameters
      real*8 EFAC,SINI              ! LOCAL: useful parameters
      real*8 R1,R2                  ! LOCAL: fractional radii
      real*8 PHASES(6)              ! LOCAL: useful orbital phases
      real*8 TRUEANOM               ! LOCAL: true anomaly
      real*8 MEANANOM               ! LOCAL: mean anomaly
      real*8 ECCANOM                ! LOCAL: eccentric anomaly
      real*8 TANTD2,PHIANOM,GUESS   ! LOCAL: helper variables
      integer i                     ! LOCAL: loop counter
      real*8 GETPHASE               ! FUNCTION: calculate orbital phase
      real*8 SEPR                   ! LOCAL: star-planet separation
      real*8 PHANGLE                ! LOCAL: phase angle

      real*8 p                ! in: ratio of planet to stellar radius
      real*8 c1,c2,c3,c4      ! in: nonlinear limb darkening coeffs
      real*8 z                ! in: orbital phase
      real*8 mu               ! out: flux relative to unobscured source
      real*8 pi,norm,i1,x,tmp ! local: various parameters
      real*8 iofr             ! function: from Mandel & Agol
!             real*8 z                ! in: impact par (ignore planet radius)

      PI = atan(1.0d0) * 4.0d0

            ! Firstly set the contribution from reflected light to zero

      V(11) = 0.0d0
      V(12) = 0.0d0

            ! Transform the EBOP pars to those needed by OCCULTSMALL

      if ( V(2) >= 0.0d0 ) then
        R1 = V(2) / (1.0d0 + V(3))
        R2 = V(2) / (1.0d0 + (1.0d0/V(3)))
      else
        R1 = abs(V(2))
        R2 = V(3)
      end if
      p = R2 / R1

      c1 = V(4)
      c2 = V(21)
      c3 = V(22)
      c4 = V(23)

!       print'(A7,f15.5)', "HJD:   ", HJD
!       print'(A7,2(f10.5))', "E,OMEGA", ecc,omega
!       print'(A7,f10.5)', "R1:    ", R1
!       print'(A7,f10.5)', "R2:    ", R2
!       print'(A7,f10.5)', "p:     ", p
!       print'(A7,2(f10.5))', "c1,c2: ", c1,c2
!       print'(A7,2(f10.5))', "c3,c4: ", c3,c4

            ! Need the impact parameter for this orbital phase (which is
            ! not the same as the "impact parameter" referred to in gen-
            ! eral which is actually the minimum impact parameter).
            ! Get this from the mean anomaly, then the eccentric anomaly
            ! and then the true anomaly

      CALL GETEOMEGA (V(7),V(8),ECC,OMEGA,ECOSW,ESINW)
      EFAC = 1.0d0 - ECC**2
      SINI = sin( V(6) / (45.0d0/atan(1.0d0)) )

      CALL ECQUADPHASES (V(7),V(8),V(16),PHASES)
      PHIANOM = getphase(HJD,V(19),V(20)) - PHASES(5)
      MEANANOM = 2.0d0 * PI * mod(PHIANOM,1.0d0)

      GUESS = MEANANOM
      do i = 1,100
         ECCANOM = MEANANOM + ECC*sin(GUESS)
         if ( abs(ECCANOM-GUESS) < 1.0d-5 ) exit
         GUESS = ECCANOM
         if (i >= 100) write (6,'(A26,A54)') "## Warning: calculation ",
     &         "of mean anomaly in OCCULTSMALL did not converge.      "
      end do

      SEPR = 1.0d0 - ECC*cos(ECCANOM)

      TANTD2 = sqrt((1.0d0+ECC)/(1.0d0-ECC)) * tan(ECCANOM/2.0d0)
      TRUEANOM = atan(TANTD2) * 2.0d0

      PHANGLE = acos( SINI * cos(TRUEANOM + OMEGA - PI/2.0d0) )

      if ( cos(PHANGLE) > 0.0d0 ) then
        z = SEPR * sin(PHANGLE) / R1
      else
        z = 1.0d0 / R1
      end if

            ! Here comes the original code from OCCULTSMALL:

      pi=acos(-1.d0)
      norm=pi*(1.d0-c1/5.d0-c2/3.d0-3.d0*c3/7.d0-c4/2.d0)
      i1=1.d0-c1-c2-c3-c4

      mu=1.d0
      if(z.gt.1.d0-p.and.z.lt.1.d0+p) then
        x=1.d0-(z-p)**2
        tmp=(1.d0-c1*(1.d0-0.8d0*x**0.25d0)
     &           -c2*(1.d0-2.d0/3.d0*x**0.5d0)
     &           -c3*(1.d0-4.d0/7.d0*x**0.75d0)
     &           -c4*(1.d0-0.5d0*x))
        mu=1.d0-tmp*(p**2*acos((z-1.d0)/p)
     &      -(z-1.d0)*sqrt(p**2-(z-1.d0)**2))/norm
      endif
      if(z.le.1.d0-p.and.z.ne.0.d0) then
        mu=1.d0-pi*p**2*iofr(c1,c2,c3,c4,z,p)/norm
      endif
      if(z.eq.0.d0) then
        mu=1.d0-pi*p**2/norm
      endif

            ! Finally, transform the output flux into a magnitude
            ! Account for scale factor and third light here.

      MAG = -2.5d0 * log10( (mu*(1.0d0-V(15))) + V(15) )  +  V(17)

!       write(72,'(20(f15.7,1x))') hjd,phianom,meananom,eccanom,sepr,
!      &trueanom,phangle,z,mu,mag

      END SUBROUTINE OCCULTSMALL
!-----------------------------------------------------------------------
      DOUBLEPRECISION FUNCTION IOFR(c1,c2,c3,c4,r,p)     ! Mandel & Agol
      implicit none
      real*8 r,p,c1,c2,c3,c4,sig1,sig2
      sig1=sqrt(sqrt(1.d0-(r-p)**2))
      sig2=sqrt(sqrt(1.d0-(r+p)**2))
      iofr=1.d0-c1*(1.d0+(sig2**5-sig1**5)/5.d0/p/r)
     &         -c2*(1.d0+(sig2**6-sig1**6)/6.d0/p/r)
     &         -c3*(1.d0+(sig2**7-sig1**7)/7.d0/p/r)
     &         -c4*(p**2+r**2)
      END FUNCTION IOFR
!=======================================================================
!=======================================================================

!=======================================================================
!=================     NUMERICAL RECIPES SUBROUTINES     ===============
!=======================================================================
      SUBROUTINE MRQMIN (x,y,sig,DTYPE,ndata,a,ia,ma,covar,alpha,nca,
     &                   chisq,alamda,ifail)
      implicit none
      integer NDATA,MA,NCA          ! IN: NDATA, numcoeffs, maxcoeffs
      real*8 X(ndata),Y(ndata)      ! IN: data to be fitted
      real*8 SIG(ndata)             ! IN: data errors in y
      integer DTYPE(ndata)          ! IN: types of the datapoints
      real*8 A(ma)                  ! IN: coefficients
      integer IA(ma)                ! IN: adjust (1) or fix (0) coeffs
      real*8 COVAR(nca,nca)         ! OUT: curvature matrix
      real*8 ALPHA(nca,nca)         ! OUT: covariance matrix
      real*8 ALAMDA                 ! IN/OUT: Marquardt lambda factor
      real*8 CHISQ                  ! OUT: chi-squared of the fit

      integer MMAX,j,k,l,m,mfit,ifail
      parameter (MMAX = 540)
      real*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit

!       write(6,*)"a"
!       write(6,*)a
!       write(6,*)"ia"
!       write(6,*)ia
!       write(6,*)ndata,x(1),x(ndata)
!       write(6,*)ndata,y(1),y(ndata)
!       write(6,*)ndata,sig(1),sig(ndata)
!       write(6,*)"nca,chisq,alamda,ifail"
!       write(6,*)nca,chisq,alamda,ifail

      if(alamda.lt.0.0d0)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).eq.1) mfit=mfit+1
11      continue
        alamda=0.0010d0
        OCHISQ = 1.0d10
        call mrqcof(x,y,sig,DTYPE,ndata,a,ia,ma,alpha,beta,nca,OCHISQ,
     &              chisq)
        OCHISQ = CHISQ
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif

      j=0
      do 14 l=1,ma
        if(ia(l).eq.1) then
          j=j+1
          k=0
          do 13 m=1,ma
            if(ia(m).eq.1) then
              k=k+1
              covar(j,k)=alpha(j,k)
            endif
13        continue
          covar(j,j)=alpha(j,j)*(1.0d0+alamda)
          da(j)=beta(j)
        endif
14    continue
      call gaussj(covar,mfit,nca,da,1,1,ifail)
      if ( ifail /= 0 ) return

      if(alamda.eq.0.0d0)then
        call covsrt(covar,nca,ma,ia,mfit)
        return
      endif

      j=0
      do 15 l=1,ma
        if(ia(l).eq.1) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,DTYPE,ndata,atry,ia,ma,covar,da,nca,OCHISQ,
     &            chisq)
      if(chisq.lt.ochisq)then
        alamda=0.1d0*alamda
        ochisq=chisq
        j=0
        do 17 l=1,ma
          if(ia(l).eq.1) then
            j=j+1
            k=0
            do 16 m=1,ma
              if(ia(m).eq.1) then
                k=k+1
                alpha(j,k)=covar(j,k)
              endif
16          continue
            beta(j)=da(j)
            a(l)=atry(l)
          endif
17      continue
      else
        alamda=10.0d0*alamda
        chisq=ochisq
      endif
      return
      END
!-----------------------------------------------------------------------
      SUBROUTINE mrqcof (x,y,sig,DTYPE,ndata,a,ia,ma,alpha,beta,nalp,
     &                                                     OCHISQ,chisq)
      implicit none
      integer ma,nalp,NDATA,ia(ma),MMAX,mfit,i,j,k,l,m
      parameter ( MMAX = 540 )
      real*8 OCHISQ,CHISQ
      real*8 a(ma),alpha(nalp,nalp),beta(ma)
      real*8 wt,dyda(MMAX)
      real*8 X(NDATA),Y(NDATA),SIG(NDATA)
      real*8 YMOD(NDATA),DY(NDATA),SIG2(NDATA)
      integer DTYPE(ndata)

            ! Modified by JKT, 16/05/2007
            ! Now if CHISQ > OCHISQ it doesn't waste computing time in
            ! calculating the derivatives (as these are not used)

      mfit=0
      do j = 1,ma
        if ( ia(j)==1 ) mfit = mfit+1
      end do

      do j = 1,mfit
        do k = 1,j
          alpha(j,k) = 0.0d0
        end do
        beta(j) = 0.0d0
      end do

      CHISQ = 0.0d0
      do i = 1,NDATA

        CALL EBOP (DTYPE(i),x(i),a,ymod(i),dyda,ma,ia,'n')
!       CALL EBOP (INDEX,X,V,Y,DYDA,NCOEFFS,VARY,DODYDA)

        SIG2(i) = 1.0d0 / (SIG(i)*SIG(i))
        DY(i) = Y(i) - YMOD(i)
        CHISQ = CHISQ + DY(i)*DY(i)*SIG2(i)
      end do
      if ( CHISQ > OCHISQ ) return

      do i = 1,NDATA
        CALL EBOP (DTYPE(i),x(i),a,ymod(i),dyda,ma,ia,'y')

        j=0
        do l=1,ma
          if(ia(l).eq.1) then
            j=j+1
            wt=dyda(l)*sig2(i)
            k=0
            do m=1,l
              if(ia(m).eq.1) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
            end do
            beta(j)=beta(j)+dy(i)*wt
          endif
        end do
      end do

      do j = 2,mfit
        do k = 1,j-1
          alpha(k,j) = alpha(j,k)
        end do
      end do

      END
!-----------------------------------------------------------------------
      SUBROUTINE GAUSSJ (a,n,np,b,m,mp,ifail)
      implicit none
      integer m,mp,n,np,NMAX,ifail
      real*8 a(np,np),b(np,mp),big,dum,pivinv
      parameter ( NMAX = 540 )
      integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      ifail=0
      irow=0
      icol=0
!       write(6,*)a
!       write(6,*)" "
      do j=1,n
        ipiv(j)=0
      end do
      do i=1,n
        big=0.0d0
        do j=1,n
          if(ipiv(j).ne.1)then
            do k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(6,'(A25,A27,I3,A1)') "### PROBLEM A: singular m",
     &                               "atrix in gaussj (parameter:",k,")"
                ifail = 1
                return
              endif
            end do
          endif
        end do
!         print'(i5,$)',icol
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
          end do
          do l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
          end do
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.0d0) then
          write(6,'(A31,A21,I3,A1)') "### PROBLEM B: singular matrix ",
     &                                  "in gaussj (parameter:",ICOL,")"
          ifail = 1
          return
        end if
        pivinv=1.0d0/a(icol,icol)
        a(icol,icol)=1.0d0
        do l=1,n
          a(icol,l)=a(icol,l)*pivinv
        end do
        do l=1,m
          b(icol,l)=b(icol,l)*pivinv
        end do
        do ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.0d0
            do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
            end do
            do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
            end do
          endif
        end do
      end do
      do l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
          end do
        endif
      end do
      return
      END
!-----------------------------------------------------------------------
      SUBROUTINE COVSRT (covar,npc,ma,ia,mfit)
      implicit none
      integer ma,mfit,npc,ia(ma),i,j,k
      real*8 covar(npc,npc),swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.0d0
          covar(j,i)=0.0d0
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).eq.1)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END
!=======================================================================
!=======================================================================
!======================================================================
      INTEGER FUNCTION SEEDSTART ()
            ! This uses the Fortran-intrinsic function SYSTEM_CLOCK to
            ! generate an integer SEED dependent on real time.
            ! SYSTEM_CLOCK returns the present time, the number of
            ! times a second this changes (100 on my computer) and the
            ! maximum value of present time (after which present time
            ! starts again from zero) (2147483647 on my computer).
            ! SEED is outputted as a four-figure integer.
      implicit none
      integer a,b                   ! unused constant values

      CALL SYSTEM_CLOCK (SEEDSTART,a,b)
                        ! SEEDSTART becomes the current clock count

      END FUNCTION SEEDSTART
!=====================================================================
      DOUBLEPRECISION FUNCTION RANDOM (SEED)
            ! Minimal random number generator from Park and Miller.
            ! Return a random deviate between 0.0 and 1.0 with a flat
            ! ditribution. Its period is (2**31 - 1) = 2.15e9
            ! The constants are the best found by Park and Miller.
            ! Set and reset SEED to any integer value to inititalize
            ! the sequence and do not change it thereafter - the code
            ! updates SEED itself.
            ! Does not work if SEED=0, unless it bit-swaps.

      implicit none
      integer SEED                  ! All-important seed value
      integer IA,IM,IR,IQ,MASK      ! Various constants
      real*8 AM                     ! Another constant
      integer K                     ! Yet another constant

      IA = 16807
      IM = 2147483647               ! 2**31 ie largest integer the ..
      AM = 1.0d0 / IM               ! .. computer can handle.
      IQ = 127773
      IR = 2836
      MASK = 123459876

      SEED = IEOR(SEED,MASK)              ! IEOR is integer exclusive-or
      K = SEED / IQ                       !  bit-swapping, the non-nume-
      SEED = IA * (SEED - K*IQ) - K*IR    !  ric operation needed by all
      if (SEED < 0.0d0) SEED = SEED + IM  !  random number generstors.
      RANDOM = AM * SEED
      SEED = IEOR(SEED,MASK)              ! Seed is updated here

      END FUNCTION RANDOM
!=======================================================================
      DOUBLEPRECISION FUNCTION RANDOMG (SEED,MEAN,SD)
            ! This produces a random number with a Gaussian distribution
            ! It uses RANDOM to generate a random number distribution
            ! with zero mean and unit variance then scales the result
            ! with SD and offsets it with MEAN to produce the effect.
      implicit none
      integer SEED                  ! Seeding integer
      real*8 MEAN,SD                ! Desired mean and S.D. (variance)
      integer ISET                  ! LOCAL: integer variable
      real*8 FAC,GSET,RSQ,V1,V2     ! LOCAL: real variables
      real*8 RANDOM            ! FUNCTION: flat-distrib random generator
      SAVE ISET,GSET
      integer i

      if ( SEED < 0 ) ISET = 0      ! Reinitialise
      if ( ISET == 0 ) then
        do i = 1,100000
          V1 = 2.0d0 * RANDOM(SEED) - 1.0d0
          V2 = 2.0d0 * RANDOM(SEED) - 1.0d0
          RSQ = V1**2 + V2**2
          if ( RSQ /= 0.0d0 .and. RSQ <= 1.0d0 ) exit
        end do
        FAC = SQRT( -2.0d0 * log(RSQ) / RSQ )
        GSET = V1 * FAC
        RANDOMG = V2 * FAC
        ISET = 1
      else
        RANDOMG = GSET
        ISET = 0
      end if
      RANDOMG = ( RANDOMG * SD ) + MEAN

      END FUNCTION RANDOMG
!=====================================================================
      DOUBLEPRECISION FUNCTION SELLECT (ARRAY,NUM,K)
            ! Returns the Kth smallest value in ARRAY(NUM).
            ! ARRAY is simply sorted during this procedure.
      implicit none
      integer NUM                   ! IN: size of the array
      real*8 ARRAY(NUM)             ! IN: Input array
      integer K                     ! OUT: value to find
      real*8 STORE                  ! LOCAL: storage variable
      integer i,j,TAG               ! LOCAL: loop counters and a flag

      TAG = 0
      do i = 1,NUM
        STORE = 10.0**10.0
        do j = i,NUM
            if ( ARRAY(j) < STORE )  then
            STORE = ARRAY(j)
            TAG = j
          end if
        end do
        ARRAY(TAG) = ARRAY(i)
        ARRAY(i) = STORE
      end do
      SELLECT = ARRAY(K)

      END FUNCTION SELLECT
!=======================================================================
      DOUBLEPRECISION FUNCTION SIGMA (ARRAY,NUM)
            ! Returns the standard deviation of array ARRAY
      implicit none
      integer NUM                   ! IN: size of the array
      real*8 ARRAY(NUM)             ! IN: array
      real*8 MEAN,VAR,SUMS,SUMSQ    ! LOCAL: properties of ARRAY values
      integer i                     ! LOCAL: loop counter

      SIGMA = 0.0d0
      SUMS = 0.0d0
      SUMSQ = 0.0d0
      do i = 1,NUM
        SUMS = SUMS + ARRAY(i)
        SUMSQ = SUMSQ + ARRAY(i)**2
      end do

      MEAN = SUMS / NUM
      VAR = (SUMSQ / NUM) - MEAN**2
      if ( VAR > 0.0d0 )  SIGMA = sqrt( VAR )

      END FUNCTION SIGMA
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
