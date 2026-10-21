! ***********************************************************************
! bolometric.f90
!
! SED_Model's counterpart to MESA colors/private/bolometric.f90.
!
! calculate_bolometric_phot matches MESA's name and role: integrate a
! diluted SED and convert the result to a bolometric magnitude. It is
! built from bolometric_flux / bolometric_magnitude, which are existing
! SED_Model steps (already separate in the pre-restructure kernel, kept
! separate here since sed_model/forward.py surfaces the flux and the
! magnitude as two distinct ForwardResult fields).
!
! Not ported: MESA's calculate_bolometric owns the interpolation-method
! dispatch (Hermite / Hermite_bounded / Linear / KNN) and the
! interpolation-radius diagnostic, both driven by the Colors_General_Info
! handle. In SED_Model that dispatch already happens in Python
! (sed_model/forward.py's interp_method argument), so it has no
! Fortran-layer counterpart here.
!
! Copyright (C) 2025 Niall Miller
! LGPL-3.0-or-later
! ***********************************************************************

module bolometric
   use colors_def, only: dp, BAD_MAG, F_BOL_ZERO_POINT
   use colors_utils, only: simpson_integration, is_finite
   implicit none

   private
   public :: bolometric_flux, bolometric_magnitude, calculate_bolometric_phot

contains

   ! Integrate the diluted SED over all wavelengths.  As in MESA, samples with
   ! a non-finite or non-positive wavelength, or a non-finite or negative flux,
   ! are zeroed before integrating.
   subroutine bolometric_flux(sed_wave, obs_flux, bol_flux, ierr)
      real(dp), intent(in)  :: sed_wave(:), obs_flux(:)
      real(dp), intent(out) :: bol_flux
      integer,  intent(out) :: ierr

      real(dp), allocatable :: clean_flux(:)
      integer :: i

      ierr = 0
      if (size(sed_wave) /= size(obs_flux) .or. size(sed_wave) < 2) then
         bol_flux = 0.0_dp
         ierr = 1
         return
      end if

      allocate (clean_flux(size(obs_flux)))
      clean_flux = obs_flux
      do i = 1, size(sed_wave)
         if (.not. is_finite(sed_wave(i)) .or. sed_wave(i) <= 0.0_dp .or. &
             .not. is_finite(clean_flux(i)) .or. clean_flux(i) < 0.0_dp) then
            clean_flux(i) = 0.0_dp
         end if
      end do

      call simpson_integration(sed_wave, clean_flux, bol_flux)
      if (bol_flux <= 0.0_dp) ierr = 1
   end subroutine bolometric_flux


   ! Mag_bol = -2.5 * log10(F_bol / F_bol_0)  [F_bol in erg/s/cm^2]
   ! F_bol_0 is the irradiance of an M_bol = 0 source (IAU 2015 B2), as in
   ! MESA's const_def f_bol_zero_point.
   subroutine bolometric_magnitude(bol_flux, bol_mag, ierr)
      real(dp), intent(in)  :: bol_flux
      real(dp), intent(out) :: bol_mag
      integer,  intent(out) :: ierr

      ierr = 0
      if (bol_flux > 0.0_dp) then
         bol_mag = -2.5_dp * log10(bol_flux / F_BOL_ZERO_POINT)
      else
         bol_mag = BAD_MAG
         ierr = 1
      end if
   end subroutine bolometric_magnitude


   ! Combined entry point -- matches MESA's calculate_bolometric_phot
   ! name and role (integrate, then convert to a magnitude in one call).
   subroutine calculate_bolometric_phot(sed_wave, obs_flux, bol_mag, bol_flux, ierr)
      real(dp), intent(in)  :: sed_wave(:), obs_flux(:)
      real(dp), intent(out) :: bol_mag, bol_flux
      integer,  intent(out) :: ierr

      integer :: ierr2
      call bolometric_flux(sed_wave, obs_flux, bol_flux, ierr)
      call bolometric_magnitude(bol_flux, bol_mag, ierr2)
      if (ierr2 /= 0 .and. ierr == 0) ierr = ierr2
   end subroutine calculate_bolometric_phot

end module bolometric
