! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Global Data Module File
! 
! Generated by KPP-2.1 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : r2smh_Global.f90
! Time                 : Fri Mar 11 14:25:20 2016
! Working directory    : /global/scratch/azare/Josh/WRFV3/chem/KPP/mechanisms/r2smh
! Equation file        : r2smh.kpp
! Output root filename : r2smh
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE r2smh_Global

  USE r2smh_Parameters, ONLY: dp, NSPEC, NVAR, NFIX, NREACT
  PUBLIC
  SAVE


! Declaration of global variables

! C - Concentration of all species
  REAL(kind=dp) :: C(NSPEC)
! VAR - Concentrations of variable species (global)
  REAL(kind=dp) :: VAR(NVAR)
! FIX - Concentrations of fixed species (global)
  REAL(kind=dp) :: FIX(NFIX)
! VAR, FIX are chunks of array C
      EQUIVALENCE( C(1),VAR(1) )
      EQUIVALENCE( C(218),FIX(1) )
! RCONST - Rate constants (global)
  REAL(kind=dp) :: RCONST(NREACT)
! TIME - Current integration time
  REAL(kind=dp) :: TIME
! SUN - Sunlight intensity between [0,1]
  REAL(kind=dp) :: SUN
! TEMP - Temperature
  REAL(kind=dp) :: TEMP
! RTOLS - (scalar) Relative tolerance
  REAL(kind=dp) :: RTOLS
! TSTART - Integration start time
  REAL(kind=dp) :: TSTART
! TEND - Integration end time
  REAL(kind=dp) :: TEND
! DT - Integration step
  REAL(kind=dp) :: DT
! ATOL - Absolute tolerance
  REAL(kind=dp) :: ATOL(NSPEC)
! RTOL - Relative tolerance
  REAL(kind=dp) :: RTOL(NSPEC)
! STEPMIN - Lower bound for integration step
  REAL(kind=dp) :: STEPMIN
! STEPMAX - Upper bound for integration step
  REAL(kind=dp) :: STEPMAX
! CFACTOR - Conversion factor for concentration units
  REAL(kind=dp) :: CFACTOR

! INLINED global variable declarations

! INLINED global variable declarations


END MODULE r2smh_Global

