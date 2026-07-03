! ***********************************************************************
! synthetic.f90
!
! SED_Model's counterpart to MESA colors/private/synthetic.f90.
!
! calculate_synthetic_flux and the three compute_*_zero_point functions
! match MESA's names and (for the zero-points) MESA's function-return
! convention: they return the zero-point directly and use -1.0 as the
! failure sentinel, exactly as MESA's colors/private/synthetic.f90
! does, rather than an ierr output. `magnitude` has no MESA counterpart
! -- MESA inlines the -2.5*log10(flux/zp) conversion directly inside
! its single fused calculate_synthetic function, whereas SED_Model
! keeps flux computation and magnitude conversion as separate steps
! (an existing SED_Model design choice, not something introduced here).
!
! MESA's calculate_synthetic also owns filter-onto-SED interpolation
! (via knn_interp), optional SED CSV output, and results-directory
! bookkeeping. None of that is reproduced: filter interpolation is
! colors_utils' interp_filter_onto_sed (linear, not KNN -- see that
! file for why), and SED_Model's I/O lives in sed_model/io.py rather
! than in the Fortran layer.
!
! Copyright (C) 2025 Niall Miller
! LGPL-3.0-or-later
! ***********************************************************************

module synthetic
   use colors_def, only: dp, CLIGHT_CM_S, AB_FNU_ZP, ST_FLAM_ZP, BAD_MAG
   use colors_utils, only: simpson_integration, interp_filter_onto_sed
   implicit none

   private
   public :: calculate_synthetic_flux, magnitude
   public :: compute_vega_zero_point, compute_ab_zero_point, compute_st_zero_point

contains

   ! Photon-counting in-band flux from an already-convolved SED:
   !   F_band = integral(convolved_flux * lambda) / integral(filter * lambda)
   ! `convolved_flux` is obs_flux * filter_on_sed_grid, computed by the
   ! caller -- matches MESA's calling convention in synthetic.f90.
   subroutine calculate_synthetic_flux(wavelengths, convolved_flux, filter_on_sed_grid, synthetic_flux)
      real(dp), intent(in)  :: wavelengths(:), convolved_flux(:), filter_on_sed_grid(:)
      real(dp), intent(out) :: synthetic_flux

      real(dp) :: integrated_flux, integrated_filter

      call simpson_integration(wavelengths, convolved_flux * wavelengths, integrated_flux)
      call simpson_integration(wavelengths, filter_on_sed_grid * wavelengths, integrated_filter)

      if (integrated_filter > 0.0_dp) then
         synthetic_flux = integrated_flux / integrated_filter
      else
         synthetic_flux = -1.0_dp
      end if
   end subroutine calculate_synthetic_flux


   subroutine magnitude(band_flux, zero_point, mag, ierr)
      real(dp), intent(in)  :: band_flux, zero_point
      real(dp), intent(out) :: mag
      integer,  intent(out) :: ierr

      ierr = 0
      if (band_flux > 0.0_dp .and. zero_point > 0.0_dp) then
         mag = -2.5_dp * log10(band_flux / zero_point)
      else
         mag = BAD_MAG
         ierr = 1
      end if
   end subroutine magnitude


   ! Vega zero-point -- called once at filter-load time.
   real(dp) function compute_vega_zero_point(vega_wave, vega_flux, filt_wave, filt_trans)
      real(dp), intent(in) :: vega_wave(:), vega_flux(:)
      real(dp), intent(in) :: filt_wave(:), filt_trans(:)

      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: filt_on_vega(:)
      integer :: ierr

      allocate (filt_on_vega(size(vega_wave)))
      call interp_filter_onto_sed(filt_wave, filt_trans, vega_wave, filt_on_vega, ierr)

      call simpson_integration(vega_wave, vega_flux * filt_on_vega * vega_wave, int_flux)
      call simpson_integration(vega_wave, filt_on_vega * vega_wave, int_filter)

      if (int_filter > 0.0_dp) then
         compute_vega_zero_point = int_flux / int_filter
      else
         compute_vega_zero_point = -1.0_dp
      end if

      deallocate (filt_on_vega)
   end function compute_vega_zero_point


   ! AB zero-point -- f_nu = 3631 Jy converted to f_lambda on the filter grid.
   real(dp) function compute_ab_zero_point(filt_wave, filt_trans)
      real(dp), intent(in) :: filt_wave(:), filt_trans(:)

      integer  :: i
      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: f_ab(:)

      allocate (f_ab(size(filt_wave)))
      do i = 1, size(filt_wave)
         if (filt_wave(i) > 0.0_dp) then
            f_ab(i) = AB_FNU_ZP * (CLIGHT_CM_S * 1.0d8) / filt_wave(i)**2
         else
            f_ab(i) = 0.0_dp
         end if
      end do

      call simpson_integration(filt_wave, f_ab * filt_trans * filt_wave, int_flux)
      call simpson_integration(filt_wave, filt_trans * filt_wave, int_filter)

      if (int_filter > 0.0_dp) then
         compute_ab_zero_point = int_flux / int_filter
      else
         compute_ab_zero_point = -1.0_dp
      end if

      deallocate (f_ab)
   end function compute_ab_zero_point


   ! ST zero-point -- flat f_lambda = 3.63e-9 erg/s/cm^2/A.
   real(dp) function compute_st_zero_point(filt_wave, filt_trans)
      real(dp), intent(in) :: filt_wave(:), filt_trans(:)

      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: f_st(:)

      allocate (f_st(size(filt_wave)))
      f_st = ST_FLAM_ZP

      call simpson_integration(filt_wave, f_st * filt_trans * filt_wave, int_flux)
      call simpson_integration(filt_wave, filt_trans * filt_wave, int_filter)

      if (int_filter > 0.0_dp) then
         compute_st_zero_point = int_flux / int_filter
      else
         compute_st_zero_point = -1.0_dp
      end if

      deallocate (f_st)
   end function compute_st_zero_point

end module synthetic
