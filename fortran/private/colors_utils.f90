! ***********************************************************************
! colors_utils.f90
!
! SED_Model's counterpart to MESA colors/private/colors_utils.f90.
!
! Only the subset of MESA's colors_utils that SED_Model actually needs
! is ported here: grid-location helpers, flux dilution, the two
! integration rules, and the filter-onto-SED interpolation. Everything
! else in MESA's colors_utils.f90 (SED file/lookup-table loading, the
! stencil cache, grid-to-lookup mapping, path resolution) exists to
! support MESA's disk-backed, handle-cached loading strategy, which
! SED_Model does not use -- SED_Tools already hands SED_Model a
! fully-built in-memory flux cube, so that machinery has no SED_Model
! counterpart and is intentionally left out.
!
! The numerics track MESA: find_interval collapses dummy axes
! (all 0 / 999 / -999) and clamps out-of-range values to the boundary
! cell, simpson_integration is the non-uniform-spacing Simpson rule,
! find_nearest_point is a true nearest-node search, and
! interp_filter_onto_sed reproduces knn_interp's interpolate_array
! followed by synthetic.f90's zero_filter_outside_support.
!
! Copyright (C) 2025 Niall Miller
! LGPL-3.0-or-later
! ***********************************************************************

module colors_utils
   use colors_def, only: dp
   implicit none

   private
   public :: find_interval, find_containing_cell, find_nearest_point
   public :: find_bracket_index
   public :: dilute_flux
   public :: trapezoidal_integration, simpson_integration
   public :: interp_filter_onto_sed, zero_filter_outside_support
   public :: is_finite

contains

   ! True for ordinary finite values; false for NaN and +/-Infinity.  Written
   ! without the ieee_arithmetic intrinsic module, which some toolchains
   ! (e.g. gfortran on musl/Alpine) do not provide.  Relies on IEEE comparison
   ! semantics, so the sources must not be built with -ffast-math.
   elemental logical function is_finite(x)
      real(dp), intent(in) :: x
      is_finite = (x == x) .and. (abs(x) <= huge(x))
   end function is_finite


   ! Find the interval in a sorted array containing `val`.
   ! Returns i such that x(i) <= val <= x(i+1) and fractional position t in
   ! [0,1].  Dummy axes (all 0, 999 or -999) and single-point axes collapse
   ! to i=1, t=0.  Values outside the axis clamp to the boundary cell.
   subroutine find_interval(x, val, i, t)
      real(dp), intent(in)  :: x(:), val
      integer,  intent(out) :: i
      real(dp), intent(out) :: t

      integer :: n, lo, hi, mid
      logical :: dummy_axis

      n = size(x)

      if (n <= 1) then
         i = 1; t = 0.0_dp; return
      end if

      dummy_axis = all(x == 0.0_dp) .or. all(x == 999.0_dp) .or. all(x == -999.0_dp)
      if (dummy_axis) then
         i = 1; t = 0.0_dp; return
      end if

      if (val <= x(1)) then
         i = 1; t = 0.0_dp; return
      else if (val >= x(n)) then
         i = n - 1; t = 1.0_dp; return
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

      i = lo
      if (abs(x(i + 1) - x(i)) < 1.0e-30_dp) then
         t = 0.0_dp   ! degenerate interval -- no interpolation needed
      else
         t = (val - x(i)) / (x(i + 1) - x(i))
      end if
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


   ! Nearest grid node on each axis (absolute distance).
   subroutine find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
      real(dp), intent(in)  :: x_val, y_val, z_val
      real(dp), intent(in)  :: x_grid(:), y_grid(:), z_grid(:)
      integer,  intent(out) :: i_x, i_y, i_z

      i_x = minloc(abs(x_val - x_grid), 1)
      i_y = minloc(abs(y_val - y_grid), 1)
      i_z = minloc(abs(z_val - z_grid), 1)
   end subroutine find_nearest_point


   ! Returns idx such that grid(idx) <= val < grid(idx+1), clamped to bounds.
   subroutine find_bracket_index(grid, val, idx)
      real(dp), intent(in)  :: grid(:), val
      integer,  intent(out) :: idx

      integer :: n, lo, hi, mid

      n = size(grid)
      if (n < 2) then
         idx = 1; return
      end if

      if (val <= grid(1)) then
         idx = 1; return
      else if (val >= grid(n)) then
         idx = n - 1; return
      end if

      lo = 1; hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (val >= grid(mid)) then
            lo = mid
         else
            hi = mid
         end if
      end do
      idx = lo
   end subroutine find_bracket_index


   ! Apply the (R/d)^2 dilution factor to convert surface flux to observed flux.
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


   ! Simpson's rule for non-uniform spacing.  Panels are taken in pairs of
   ! intervals; when n is even the final interval is closed with a trapezoid.
   subroutine simpson_integration(x, y, result)
      real(dp), intent(in)  :: x(:), y(:)
      real(dp), intent(out) :: result

      integer  :: i, n
      real(dp) :: s, h1, h2, f0, f1, f2

      n = size(x)
      s = 0.0_dp

      if (n < 2 .or. size(y) /= n) then
         result = 0.0_dp; return
      end if

      do i = 1, n - 2, 2
         h1 = x(i+1) - x(i)
         h2 = x(i+2) - x(i+1)
         f0 = y(i); f1 = y(i+1); f2 = y(i+2)

         s = s + (h1 + h2) / 6.0_dp * ( &
             f0 * (2.0_dp*h1 - h2) / h1 + &
             f1 * (h1 + h2)**2 / (h1*h2) + &
             f2 * (2.0_dp*h2 - h1) / h2)
      end do

      if (mod(n, 2) == 0) then
         s = s + 0.5_dp * (x(n) - x(n-1)) * (y(n) + y(n-1))
      end if

      result = s
   end subroutine simpson_integration


   ! Linear interpolation of a filter transmission curve onto an arbitrary
   ! wavelength grid.  Matches knn_interp's interpolate_array: values outside
   ! the tabulated range are first held at the end values, then
   ! zero_filter_outside_support enforces compact support.
   ! Both grids must be sorted ascending; the filter pointer only advances,
   ! so the cost is O(n_sed + n_filt).
   subroutine interp_filter_onto_sed(filt_wave, filt_trans, sed_wave, filt_on_sed, ierr)
      real(dp), intent(in)  :: filt_wave(:), filt_trans(:)
      real(dp), intent(in)  :: sed_wave(:)
      real(dp), intent(out) :: filt_on_sed(:)
      integer,  intent(out) :: ierr

      integer  :: i, nf, lo
      real(dp) :: t, denom

      ierr = 0
      nf   = size(filt_wave)

      if (nf < 2 .or. size(filt_trans) /= nf .or. size(sed_wave) < 1) then
         filt_on_sed = 0.0_dp
         ierr = 1
         return
      end if

      lo = 1
      do i = 1, size(sed_wave)
         if (sed_wave(i) <= filt_wave(1)) then
            filt_on_sed(i) = filt_trans(1)
            cycle
         else if (sed_wave(i) >= filt_wave(nf)) then
            filt_on_sed(i) = filt_trans(nf)
            cycle
         end if

         do while (lo < nf - 1 .and. filt_wave(lo + 1) <= sed_wave(i))
            lo = lo + 1
         end do

         denom = filt_wave(lo + 1) - filt_wave(lo)
         if (abs(denom) > 0.0_dp) then
            t = (sed_wave(i) - filt_wave(lo)) / denom
         else
            t = 0.0_dp
         end if
         filt_on_sed(i) = filt_trans(lo) + t * (filt_trans(lo + 1) - filt_trans(lo))
      end do

      call zero_filter_outside_support(sed_wave, filt_wave, filt_on_sed)
   end subroutine interp_filter_onto_sed


   ! Transmission outside the supplied passband must be zero, and negative
   ! interpolated values are unphysical.  Mirrors MESA's synthetic.f90.
   subroutine zero_filter_outside_support(target_wave, filt_wave, filt_on_target)
      real(dp), intent(in)    :: target_wave(:), filt_wave(:)
      real(dp), intent(inout) :: filt_on_target(:)

      real(dp) :: fmin, fmax
      integer  :: i

      fmin = minval(filt_wave)
      fmax = maxval(filt_wave)

      do i = 1, size(target_wave)
         if (target_wave(i) < fmin .or. target_wave(i) > fmax) then
            filt_on_target(i) = 0.0_dp
         else if (filt_on_target(i) < 0.0_dp) then
            filt_on_target(i) = 0.0_dp
         end if
      end do
   end subroutine zero_filter_outside_support

end module colors_utils
