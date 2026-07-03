! ***********************************************************************
! colors_def.f90
!
! SED_Model's counterpart to MESA colors/public/colors_def.f90.
!
! MESA's colors_def.f90 defines the Colors_General_Info handle type plus
! the cache/inlist bookkeeping built around it. SED_Model has no handle:
! every kernel below receives its arrays directly from the caller, which
! is what a clean f2py binding requires (see fortran/cc_api.f90). There
! is therefore nothing to port for the handle itself -- what *does* port
! directly is the shared real kind and the physical constants used
! throughout the private kernels, so that is what this file provides.
!
! Copyright (C) 2025 Niall Miller
! LGPL-3.0-or-later
! ***********************************************************************

module colors_def
   implicit none
   private

   integer, parameter :: dp = kind(1.0d0)

   ! Photometric / bolometric constants (cgs, wavelength in Angstroms)
   real(dp), parameter :: CLIGHT_CM_S  = 2.99792458d10
   real(dp), parameter :: AB_FNU_ZP    = 3.631d-20      ! 3631 Jy in erg/s/cm^2/Hz
   real(dp), parameter :: ST_FLAM_ZP   = 3.63d-9        ! flat f_lambda zp erg/s/cm^2/A
   real(dp), parameter :: BAD_MAG      = -99.9d0        ! failure sentinel for magnitudes
   real(dp), parameter :: TINY_VALUE   = 1.0d-10

   public :: dp
   public :: CLIGHT_CM_S, AB_FNU_ZP, ST_FLAM_ZP, BAD_MAG, TINY_VALUE

end module colors_def
