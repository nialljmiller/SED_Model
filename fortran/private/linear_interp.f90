! ***********************************************************************
! linear_interp.f90
!
! SED_Model's counterpart to MESA colors/private/linear_interp.f90.
!
! trilinear_interp_vector below matches MESA's name, argument order,
! and assumed-shape interface. As with hermite_interp.f90, MESA's
! construct_sed_linear / construct_sed_from_files handle-and-cache
! wrapper has no SED_Model counterpart (SED_Tools always hands over a
! preloaded cube), and `ierr` is appended as an optional output for the
! same reason described in hermite_interp.f90.
!
! Copyright (C) 2025 Niall Miller
! LGPL-3.0-or-later
! ***********************************************************************

module linear_interp
   use colors_def, only: dp, TINY_VALUE
   use colors_utils, only: find_containing_cell, find_nearest_point, is_finite
   implicit none

   private
   public :: trilinear_interp_vector

contains

   ! Trilinear interpolation over all wavelengths.  Matches MESA's
   ! trilinear_interp_vector: out-of-range queries clamp to the boundary
   ! cell, singleton axes are degenerate, and any non-finite or
   ! non-positive interpolant falls back to the nearest grid node.
   ! ierr = 1 flags a nearest-node fallback; ierr = 2 means even the
   ! nearest node held an unusable (non-finite or negative) flux.
   subroutine trilinear_interp_vector(x_val, y_val, z_val, &
                                      x_grid, y_grid, z_grid, &
                                      f_values_4d, n_lambda, result_flux, ierr)
      real(dp), intent(in)  :: x_val, y_val, z_val
      real(dp), intent(in)  :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in)  :: f_values_4d(:,:,:,:)   ! (nx, ny, nz, n_lambda)
      integer,  intent(in)  :: n_lambda
      real(dp), intent(out) :: result_flux(n_lambda)
      integer,  intent(out), optional :: ierr

      integer  :: nx, ny, nz, lam
      integer  :: i_x, i_y, i_z, n_x, n_y, n_z
      integer  :: ix_max, iy_max, iz_max, ox, oy, oz
      real(dp) :: t_x, t_y, t_z
      real(dp) :: w_x(0:1), w_y(0:1), w_z(0:1)
      real(dp) :: accum
      integer  :: ierr_local

      nx = size(x_grid); ny = size(y_grid); nz = size(z_grid)
      ierr_local = 0

      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      call setup_axis(nx, x_val, x_grid, i_x, t_x, ix_max, w_x)
      call setup_axis(ny, y_val, y_grid, i_y, t_y, iy_max, w_y)
      call setup_axis(nz, z_val, z_grid, i_z, t_z, iz_max, w_z)

      call find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                              n_x, n_y, n_z)

      do lam = 1, n_lambda
         accum = 0.0_dp
         do oz = 0, iz_max
            do oy = 0, iy_max
               do ox = 0, ix_max
                  accum = accum + w_x(ox)*w_y(oy)*w_z(oz) &
                          * f_values_4d(i_x+ox, i_y+oy, i_z+oz, lam)
               end do
            end do
         end do

         if (.not. is_finite(accum) .or. accum <= 0.0_dp) then
            accum = f_values_4d(n_x, n_y, n_z, lam)
            ierr_local = max(ierr_local, 1)
         end if

         if (.not. is_finite(accum) .or. accum < 0.0_dp) then
            ierr_local = 2
            accum = TINY_VALUE
         end if
         result_flux(lam) = max(TINY_VALUE, accum)
      end do

      if (present(ierr)) ierr = ierr_local
   end subroutine trilinear_interp_vector


   ! Per-axis interpolation weights.  Singleton axes get a single node;
   ! queries beyond the axis are pinned to the boundary node.
   subroutine setup_axis(n, x_val, grid, i, t, i_max, w)
      integer,  intent(in)    :: n
      real(dp), intent(in)    :: x_val
      real(dp), intent(in)    :: grid(:)
      integer,  intent(inout) :: i
      real(dp), intent(inout) :: t
      integer,  intent(out)   :: i_max
      real(dp), intent(out)   :: w(0:1)

      if (n <= 1) then
         i = 1; t = 0.0_dp; i_max = 0
         w(0) = 1.0_dp; w(1) = 0.0_dp
         return
      end if

      if (x_val <= grid(1)) then
         i = 1; t = 0.0_dp
      else if (x_val >= grid(n)) then
         i = n - 1; t = 1.0_dp
      else
         if (i < 1) i = 1
         if (i >= n) i = n - 1
         t = max(0.0_dp, min(1.0_dp, t))
      end if

      i_max = 1
      w(0) = 1.0_dp - t
      w(1) = t
   end subroutine setup_axis

end module linear_interp
