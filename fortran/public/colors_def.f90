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
   real(dp), parameter :: BAD_MAG      = -99.9d0        ! failure sentinel for band magnitudes
   ! MESA colors/private/bolometric.f90 returns 99.0 as the bolometric
   ! magnitude when the integrated flux is not positive.
   real(dp), parameter :: BOL_MAG_NONPOSITIVE_FLUX = 99.0d0
   real(dp), parameter :: TINY_VALUE   = 1.0d-300        ! positivity floor for linear interpolation

   ! Bolometric zero-point: irradiance corresponding to M_bol = 0 (IAU 2015
   ! Resolution B2), i.e. a source with L = Lsun at 10 pc has M_bol = 4.74.
   ! Mirrors f_bol_zero_point in MESA const_def.
   real(dp), parameter :: LSUN_ERG_S   = 3.828d33
   real(dp), parameter :: MBOL_SUN     = 4.74d0
   real(dp), parameter :: PC_CM        = 3.0856775814913673d18
   real(dp), parameter :: PI_DP        = 3.141592653589793d0
   real(dp), parameter :: F_BOL_ZERO_POINT = &
      LSUN_ERG_S * 10.0_dp**(0.4_dp*MBOL_SUN) / (4.0_dp*PI_DP*(10.0_dp*PC_CM)**2)

   public :: dp
   public :: CLIGHT_CM_S, AB_FNU_ZP, ST_FLAM_ZP, BAD_MAG, BOL_MAG_NONPOSITIVE_FLUX, TINY_VALUE
   public :: F_BOL_ZERO_POINT

end module colors_def
