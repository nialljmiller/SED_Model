! ***********************************************************************
! colors_lib.f90
!
! SED_Model's counterpart to MESA colors/public/colors_lib.f90.
!
! MESA's colors_lib.f90 is the single module other MESA code is meant
! to `use`: it aggregates the private colors/ modules and adds handle
! lifecycle calls (colors_init, alloc_colors_handle, ...) on top. Those
! lifecycle calls have no SED_Model counterpart -- there is no handle --
! so this file keeps only the aggregation role, re-exporting the
! private kernels below under one module name, matching MESA's file
! and its purpose.
!
! NOT part of the compiled cc_api extension. numpy's f2py (at least
! under the numpy<2.0 this project's pyproject.toml pins) cannot wrap
! a module that only re-exports procedures from other modules -- it
! tries to build a Python getter for every public name as if it were
! a module variable and crashes with KeyError('void') the moment it
! hits a re-exported subroutine. fortran/cc_api.f90 therefore `use`s
! the private kernels directly instead of going through this file (see
! its header for the confirmed root cause and how it was verified).
! This file is kept for source-level parity with MESA and for any
! plain-Fortran caller that wants one aggregated `use colors_lib` --
! it compiles and behaves correctly on its own, it's simply excluded
! from the Makefile/meson.build/setup.py source lists that build the
! Python extension.
!
! Copyright (C) 2025 Niall Miller
! LGPL-3.0-or-later
! ***********************************************************************

module colors_lib
   use colors_def,    only: dp
   use hermite_interp, only: hermite_interp_vector
   use linear_interp,  only: trilinear_interp_vector
   use colors_utils,   only: dilute_flux, trapezoidal_integration, simpson_integration, &
                              find_containing_cell, find_nearest_point, find_interval, &
                              interp_filter_onto_sed
   use synthetic,       only: calculate_synthetic_flux, magnitude, &
                               compute_vega_zero_point, compute_ab_zero_point, compute_st_zero_point
   use bolometric,      only: bolometric_flux, bolometric_magnitude, calculate_bolometric_phot

   implicit none

   private

   public :: dp
   public :: hermite_interp_vector, trilinear_interp_vector
   public :: dilute_flux, trapezoidal_integration, simpson_integration
   public :: find_containing_cell, find_nearest_point, find_interval
   public :: interp_filter_onto_sed
   public :: calculate_synthetic_flux, magnitude
   public :: compute_vega_zero_point, compute_ab_zero_point, compute_st_zero_point
   public :: bolometric_flux, bolometric_magnitude, calculate_bolometric_phot

end module colors_lib
