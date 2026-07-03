! ***********************************************************************
! colors_utils.f90
!
! SED_Model's counterpart to MESA colors/private/colors_utils.f90.
!
! Only the subset of MESA's colors_utils that SED_Model actually needs
! is ported here: grid-location helpers, flux dilution, and the two
! integration rules. Everything else in MESA's colors_utils.f90 (SED
! file/lookup-table loading, the stencil cache, grid-to-lookup mapping,
! path resolution) exists to support MESA's disk-backed, handle-cached
! loading strategy, which SED_Model does not use -- SED_Tools already
! hands SED_Model a fully-built in-memory flux cube, so that machinery
! has no SED_Model counterpart and is intentionally left out.
!
! Names and argument order match MESA's colors_utils.f90 exactly
! (assumed-shape arrays, no explicit size arguments) so this file reads
! as a direct diff against the MESA original rather than a rewrite.
! The one addition is an `ierr` output on find_interval-adjacent logic
! is NOT made -- ierr handling lives one layer up, in fortran/cc_api.f90,
! which is the only file in this tree that has no MESA counterpart at
! all (it exists solely to give f2py an explicit-shape, no-derived-type
! surface to bind against).
!
! Copyright (C) 2025 Niall Miller
! LGPL-3.0-or-later
! ***********************************************************************

module colors_utils
   use colors_def, only: dp, TINY_VALUE
   implicit none

   private
   public :: find_interval, find_containing_cell, find_nearest_point
   public :: dilute_flux
   public :: trapezoidal_integration, simpson_integration
   public :: interp_filter_onto_sed

contains

   ! Locate the interval in sorted array x that contains val.
   ! Returns i such that x(i) <= val <= x(i+1), and fractional position
   ! t = (val - x(i)) / (x(i+1) - x(i)) in [0,1], clamped to [1, n-1].
   subroutine find_interval(x, val, i, t)
      real(dp), intent(in)  :: x(:), val
      integer,  intent(out) :: i
      real(dp), intent(out) :: t

      integer  :: n, lo, hi, mid
      real(dp) :: denom

      n = size(x)

      if (n <= 1) then
         i = 1; t = 0.0_dp; return
      end if

      lo = 1; hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (val >= x(mid)) then
            lo = mid
         else
            hi = mid
         end if
      end do

      i = max(1, min(lo, n - 1))

      denom = x(i + 1) - x(i)
      if (abs(denom) > 0.0_dp) then
         t = (val - x(i)) / denom
      else
         t = 0.0_dp
      end if
      t = max(0.0_dp, min(1.0_dp, t))
   end subroutine find_interval


   subroutine find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                   i_x, i_y, i_z, t_x, t_y, t_z)
      real(dp), intent(in)  :: x_val, y_val, z_val
      real(dp), intent(in)  :: x_grid(:), y_grid(:), z_grid(:)
      integer,  intent(out) :: i_x, i_y, i_z
      real(dp), intent(out) :: t_x, t_y, t_z

      call find_interval(x_grid, x_val, i_x, t_x)
      call find_interval(y_grid, y_val, i_y, t_y)
      call find_interval(z_grid, z_val, i_z, t_z)
   end subroutine find_containing_cell


   ! Return the nearest grid index for a single sorted axis.
   integer function find_nearest_1d(grid, val) result(idx)
      real(dp), intent(in) :: grid(:), val
      integer  :: i
      real(dp) :: t
      call find_interval(grid, val, i, t)
      idx = merge(i, i + 1, t < 0.5_dp)
   end function find_nearest_1d


   subroutine find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
      real(dp), intent(in)  :: x_val, y_val, z_val
      real(dp), intent(in)  :: x_grid(:), y_grid(:), z_grid(:)
      integer,  intent(out) :: i_x, i_y, i_z

      i_x = find_nearest_1d(x_grid, x_val)
      i_y = find_nearest_1d(y_grid, y_val)
      i_z = find_nearest_1d(z_grid, z_val)
   end subroutine find_nearest_point


   ! Apply (R/d)^2 dilution to convert surface flux to observed flux.
   subroutine dilute_flux(surface_flux, R, d, calibrated_flux)
      real(dp), intent(in)  :: surface_flux(:)
      real(dp), intent(in)  :: R, d
      real(dp), intent(out) :: calibrated_flux(:)
      calibrated_flux = surface_flux * (R / d)**2
   end subroutine dilute_flux


   subroutine trapezoidal_integration(x, y, result)
      real(dp), intent(in)  :: x(:), y(:)
      real(dp), intent(out) :: result
      integer :: i, n
      n = size(x)
      result = 0.0_dp
      do i = 1, n - 1
         result = result + 0.5_dp * (y(i) + y(i+1)) * (x(i+1) - x(i))
      end do
   end subroutine trapezoidal_integration


   ! Composite Simpson's rule.
   !
   ! Odd n: full composite Simpson over all n points.
   ! Even n: composite Simpson over the first n-1 points (odd), then a single
   !         trapezoid panel for the last pair.  This is O(h^4) accurate
   !         everywhere except the last panel, which is O(h^2) — far better
   !         than falling back to trapezoid over the entire array.
   subroutine simpson_integration(x, y, result)
      real(dp), intent(in)  :: x(:), y(:)
      real(dp), intent(out) :: result

      integer  :: i, n
      real(dp) :: h, s

      n = size(x)

      if (n < 2) then
         result = 0.0_dp; return
      end if

      if (mod(n, 2) == 1) then
         s = 0.0_dp
         do i = 1, n - 2, 2
            h = x(i+2) - x(i)
            s = s + h / 6.0_dp * (y(i) + 4.0_dp*y(i+1) + y(i+2))
         end do
         result = s
      else
         ! Simpson on points 1..n-1, trapezoid on the last pair n-1..n
         s = 0.0_dp
         do i = 1, n - 3, 2
            h = x(i+2) - x(i)
            s = s + h / 6.0_dp * (y(i) + 4.0_dp*y(i+1) + y(i+2))
         end do
         s = s + 0.5_dp * (y(n-1) + y(n)) * (x(n) - x(n-1))
         result = s
      end if
   end subroutine simpson_integration


   ! Linear interpolation of a filter transmission curve onto an
   ! arbitrary wavelength grid, clamped to 0 outside the filter's own
   ! wavelength range.
   !
   ! MESA gets this behaviour from knn_interp's interpolate_array,
   ! called from synthetic.f90. SED_Model does not carry knn_interp, so
   ! this small linear-interpolation utility lives here instead -- it is
   ! the one routine in this file with no direct MESA colors_utils name
   ! to mirror.
   ! Interpolate filter transmission onto the SED wavelength grid using a
   ! two-pointer walk instead of a binary search per SED point.  Both grids
   ! are assumed to be sorted ascending, so the filter pointer lo only ever
   ! advances — O(n_sed + n_filt) rather than O(n_sed * log n_filt).
   subroutine interp_filter_onto_sed(filt_wave, filt_trans, sed_wave, filt_on_sed, ierr)
      real(dp), intent(in)  :: filt_wave(:), filt_trans(:)
      real(dp), intent(in)  :: sed_wave(:)
      real(dp), intent(out) :: filt_on_sed(:)
      integer,  intent(out) :: ierr

      integer  :: i, nf, lo
      real(dp) :: t, denom

      ierr = 0
      nf   = size(filt_wave)
      lo   = 1

      do i = 1, size(sed_wave)
         if (sed_wave(i) <= filt_wave(1) .or. sed_wave(i) >= filt_wave(nf)) then
            filt_on_sed(i) = 0.0_dp
            cycle
         end if

         do while (lo < nf - 1 .and. filt_wave(lo + 1) < sed_wave(i))
            lo = lo + 1
         end do

         denom = filt_wave(lo + 1) - filt_wave(lo)
         if (abs(denom) > 0.0_dp) then
            t = (sed_wave(i) - filt_wave(lo)) / denom
         else
            t = 0.0_dp
         end if
         filt_on_sed(i) = max(0.0_dp, filt_trans(lo) + t*(filt_trans(lo + 1) - filt_trans(lo)))
      end do
   end subroutine interp_filter_onto_sed

end module colors_utils
