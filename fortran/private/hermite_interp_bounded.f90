! ***********************************************************************
! hermite_interp_bounded.f90
!
! SED_Model's counterpart to MESA colors/private/hermite_interp_bounded.f90.
!
! Hermite interpolation can undershoot below zero near sharp spectral
! features.  This wrapper accepts the Hermite SED when it is finite and
! physically usable, projects negligible negative undershoot (below 0.1%
! of the positive SED, both by integral and by peak) onto f_lambda >= 0,
! and otherwise falls back to the positivity-preserving trilinear
! interpolant.
!
! Where MESA calls mesa_error when both interpolants are unusable, this
! routine returns ierr = 2 and lets the Python layer raise.
!   ierr = 0  Hermite accepted (possibly after repair)
!   ierr = 1  query clamped to the grid, or linear fallback used
!   ierr = 2  neither interpolant returned a valid SED
!
! Copyright (C) 2026 Niall Miller
! LGPL-3.0-or-later
! ***********************************************************************

module hermite_interp_bounded
   use colors_def, only: dp
   use colors_utils, only: simpson_integration, is_finite
   use hermite_interp, only: hermite_interp_vector
   use linear_interp, only: trilinear_interp_vector
   implicit none

   private
   public :: hermite_bounded_interp_vector

   ! Accept negligible Hermite undershoot when both its integrated magnitude
   ! and its peak amplitude are below 0.1% of the positive SED.
   real(dp), parameter :: negative_flux_fraction_tol = 1.0d-3

contains

   subroutine hermite_bounded_interp_vector(x_val, y_val, z_val, &
                                            x_grid, y_grid, z_grid, &
                                            f_values_4d, n_lambda, wavelengths, &
                                            result_flux, ierr)
      real(dp), intent(in)  :: x_val, y_val, z_val
      real(dp), intent(in)  :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in)  :: f_values_4d(:,:,:,:)
      integer,  intent(in)  :: n_lambda
      real(dp), intent(in)  :: wavelengths(n_lambda)
      real(dp), intent(out) :: result_flux(n_lambda)
      integer,  intent(out) :: ierr

      real(dp) :: integrated_flux
      logical  :: interp_ok
      integer  :: ierr_h

      call hermite_interp_vector(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 f_values_4d, n_lambda, result_flux, ierr_h)
      ierr = ierr_h

      interp_ok = valid_sed(wavelengths, result_flux, .false.)
      if (interp_ok) call repair_negative_flux(wavelengths, result_flux, interp_ok)

      if (interp_ok) then
         call simpson_integration(wavelengths, result_flux, integrated_flux)
         interp_ok = is_finite(integrated_flux) .and. integrated_flux > 0.0_dp
      end if

      if (interp_ok) return

      ! Hermite was unusable: fall back to the positivity-preserving linear
      ! interpolant.
      call trilinear_interp_vector(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                   f_values_4d, n_lambda, result_flux, ierr_h)
      ierr = 1

      interp_ok = ierr_h < 2 .and. valid_sed(wavelengths, result_flux, .true.)
      if (interp_ok) then
         call simpson_integration(wavelengths, result_flux, integrated_flux)
         interp_ok = is_finite(integrated_flux) .and. integrated_flux > 0.0_dp
      end if

      if (.not. interp_ok) ierr = 2
   end subroutine hermite_bounded_interp_vector


   subroutine repair_negative_flux(wavelengths, fluxes, repair_ok)
      real(dp), intent(in)    :: wavelengths(:)
      real(dp), intent(inout) :: fluxes(:)
      logical,  intent(out)   :: repair_ok

      real(dp), allocatable :: flux_part(:)
      real(dp) :: negative_integral, positive_integral
      real(dp) :: negative_peak, positive_peak
      integer  :: n

      repair_ok = .false.
      n = size(fluxes)

      if (n < 2 .or. size(wavelengths) /= n) return

      ! Nothing needs repairing.  Zero flux is already physically safe.
      if (.not. any(fluxes < 0.0_dp)) then
         repair_ok = any(fluxes > 0.0_dp)
         return
      end if

      ! A useful SED must retain a positive interior.
      if (.not. any(fluxes > 0.0_dp)) return

      allocate (flux_part(n))

      flux_part = max(0.0_dp, -fluxes)
      call simpson_integration(wavelengths, flux_part, negative_integral)
      negative_peak = maxval(flux_part)

      flux_part = max(0.0_dp, fluxes)
      call simpson_integration(wavelengths, flux_part, positive_integral)
      positive_peak = maxval(flux_part)

      if (.not. all(is_finite([negative_integral, positive_integral, &
                                    negative_peak, positive_peak]))) return
      if (positive_integral <= 0.0_dp .or. positive_peak <= 0.0_dp) return

      if (negative_integral/positive_integral > negative_flux_fraction_tol) return
      if (negative_peak/positive_peak > negative_flux_fraction_tol) return

      ! The undershoot is negligible; project only the unphysical samples onto
      ! the physical constraint f_lambda >= 0.
      where (fluxes < 0.0_dp) fluxes = 0.0_dp

      repair_ok = all(fluxes >= 0.0_dp) .and. any(fluxes > 0.0_dp)
   end subroutine repair_negative_flux


   logical function valid_sed(wavelengths, fluxes, require_positive_flux)
      real(dp), intent(in) :: wavelengths(:), fluxes(:)
      logical,  intent(in) :: require_positive_flux

      valid_sed = size(wavelengths) >= 2 .and. size(wavelengths) == size(fluxes)
      if (.not. valid_sed) return

      valid_sed = all(is_finite(wavelengths)) .and. all(is_finite(fluxes))
      if (.not. valid_sed) return

      valid_sed = all(wavelengths > 0.0_dp) .and. &
                  all(wavelengths(2:) > wavelengths(:size(wavelengths) - 1))
      if (.not. valid_sed) return

      if (require_positive_flux) valid_sed = all(fluxes > 0.0_dp)
   end function valid_sed

end module hermite_interp_bounded
