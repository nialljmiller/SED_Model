! ***********************************************************************
! hermite_interp.f90
!
! SED_Model's counterpart to MESA colors/private/hermite_interp.f90.
!
! This is the file the restructure was for: hermite_interp_vector below
! has the exact same name, argument order, and assumed-shape array
! interface as MESA's hermite_interp_vector (x_val, y_val, z_val,
! x_grid, y_grid, z_grid, f_values_4d, n_lambda, result_flux), and the
! basis-function / finite-difference-tangent algorithm is unchanged
! from the existing SED_Model kernel. Dropping this file into either
! tree and pointing its `use colors_utils` line at that tree's
! colors_utils module is a same-file swap for the numerical core.
!
! Two things could not be ported verbatim, both structural rather than
! numerical:
!   - MESA's construct_sed_hermite / construct_sed_from_files wrap this
!     routine with the Colors_General_Info handle: stencil caching and
!     a disk-loading fallback for when the flux cube isn't preloaded.
!     SED_Model always receives a preloaded cube from SED_Tools, so
!     those wrappers have no SED_Model counterpart and are left out.
!   - `ierr` is appended as an OPTIONAL output. MESA's call sites never
!     pass it (out-of-grid queries silently fall back to the nearest
!     grid point); SED_Model's forward.py surfaces that fallback as a
!     `clamped` diagnostic, which is the one existing SED_Model
!     behaviour this restructure had to preserve. Because it's
!     optional, MESA-style call sites that omit it still work.
!
! Copyright (C) 2025 Niall Miller
! LGPL-3.0-or-later
! ***********************************************************************

module hermite_interp
   use colors_def, only: dp, TINY_VALUE
   use colors_utils, only: find_containing_cell, find_nearest_point
   implicit none

   private
   public :: hermite_interp_vector

contains

   ! Vectorised Hermite tensor interpolation over all wavelengths.
   ! The cell location (i_x, i_y, i_z, t_x, t_y, t_z) depends only on
   ! (x_val, y_val, z_val) and the grids -- not on wavelength -- so it
   ! is computed once and reused across all n_lambda samples.
   subroutine hermite_interp_vector(x_val, y_val, z_val, &
                                    x_grid, y_grid, z_grid, &
                                    f_values_4d, n_lambda, result_flux, ierr)
      real(dp), intent(in)  :: x_val, y_val, z_val
      real(dp), intent(in)  :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in)  :: f_values_4d(:,:,:,:)   ! (nx, ny, nz, n_lambda)
      integer,  intent(in)  :: n_lambda
      real(dp), intent(out) :: result_flux(n_lambda)
      integer,  intent(out), optional :: ierr

      integer  :: nx, ny, nz
      integer  :: i_x, i_y, i_z, ix, iy, iz, lam
      real(dp) :: t_x, t_y, t_z, dx, dy, dz, f_sum
      real(dp) :: df_dx, df_dy, df_dz
      real(dp) :: h_x(2), hx_d(2), h_y(2), hy_d(2), h_z(2), hz_d(2)
      integer  :: ierr_local

      nx = size(x_grid); ny = size(y_grid); nz = size(z_grid)
      ierr_local = 0

      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      ! Outside grid: fall back to nearest point
      if (i_x < 1 .or. i_x >= nx .or. &
          i_y < 1 .or. i_y >= ny .or. &
          i_z < 1 .or. i_z >= nz) then
         call find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
         result_flux = f_values_4d(i_x, i_y, i_z, :)
         ierr_local = 1   ! flag: clamped to boundary
         if (present(ierr)) ierr = ierr_local
         return
      end if

      dx = x_grid(i_x + 1) - x_grid(i_x)
      dy = y_grid(i_y + 1) - y_grid(i_y)
      dz = z_grid(i_z + 1) - z_grid(i_z)

      ! Precompute Hermite basis (same for all wavelengths)
      h_x  = [h00(t_x), h01(t_x)]
      hx_d = [h10(t_x), h11(t_x)]
      h_y  = [h00(t_y), h01(t_y)]
      hy_d = [h10(t_y), h11(t_y)]
      h_z  = [h00(t_z), h01(t_z)]
      hz_d = [h10(t_z), h11(t_z)]

      ! Wavelength loop: only the derivative reads and sum vary
      do lam = 1, n_lambda
         f_sum = 0.0_dp
         do iz = 0, 1
            do iy = 0, 1
               do ix = 0, 1
                  call compute_derivatives_4d( &
                     f_values_4d, nx, ny, nz, n_lambda, &
                     i_x+ix, i_y+iy, i_z+iz, lam, &
                     x_grid, y_grid, z_grid, &
                     df_dx, df_dy, df_dz)

                  f_sum = f_sum &
                     + h_x(ix+1)  * h_y(iy+1)  * h_z(iz+1)  * f_values_4d(i_x+ix, i_y+iy, i_z+iz, lam) &
                     + hx_d(ix+1) * h_y(iy+1)  * h_z(iz+1)  * dx * df_dx &
                     + h_x(ix+1)  * hy_d(iy+1) * h_z(iz+1)  * dy * df_dy &
                     + h_x(ix+1)  * h_y(iy+1)  * hz_d(iz+1) * dz * df_dz
               end do
            end do
         end do
         result_flux(lam) = max(TINY_VALUE, f_sum)
      end do

      if (present(ierr)) ierr = ierr_local
   end subroutine hermite_interp_vector


   subroutine compute_derivatives_4d(f, nx, ny, nz, n_lambda, &
                                     i, j, k, lam, &
                                     x_grid, y_grid, z_grid, &
                                     df_dx, df_dy, df_dz)
      integer,  intent(in)  :: nx, ny, nz, n_lambda
      real(dp), intent(in)  :: f(nx, ny, nz, n_lambda)
      integer,  intent(in)  :: i, j, k, lam
      real(dp), intent(in)  :: x_grid(nx), y_grid(ny), z_grid(nz)
      real(dp), intent(out) :: df_dx, df_dy, df_dz

      if (nx == 1) then
         df_dx = 0.0_dp
      else if (i > 1 .and. i < nx) then
         df_dx = (f(i+1,j,k,lam) - f(i-1,j,k,lam)) / (x_grid(i+1) - x_grid(i-1))
      else if (i == 1) then
         df_dx = (f(i+1,j,k,lam) - f(i,j,k,lam))   / (x_grid(i+1) - x_grid(i))
      else
         df_dx = (f(i,j,k,lam)   - f(i-1,j,k,lam)) / (x_grid(i)   - x_grid(i-1))
      end if

      if (ny == 1) then
         df_dy = 0.0_dp
      else if (j > 1 .and. j < ny) then
         df_dy = (f(i,j+1,k,lam) - f(i,j-1,k,lam)) / (y_grid(j+1) - y_grid(j-1))
      else if (j == 1) then
         df_dy = (f(i,j+1,k,lam) - f(i,j,k,lam))   / (y_grid(j+1) - y_grid(j))
      else
         df_dy = (f(i,j,k,lam)   - f(i,j-1,k,lam)) / (y_grid(j)   - y_grid(j-1))
      end if

      if (nz == 1) then
         df_dz = 0.0_dp
      else if (k > 1 .and. k < nz) then
         df_dz = (f(i,j,k+1,lam) - f(i,j,k-1,lam)) / (z_grid(k+1) - z_grid(k-1))
      else if (k == 1) then
         df_dz = (f(i,j,k+1,lam) - f(i,j,k,lam))   / (z_grid(k+1) - z_grid(k))
      else
         df_dz = (f(i,j,k,lam)   - f(i,j,k-1,lam)) / (z_grid(k)   - z_grid(k-1))
      end if
   end subroutine compute_derivatives_4d


   ! Hermite basis functions
   pure real(dp) function h00(t)
      real(dp), intent(in) :: t
      h00 = (1.0_dp + 2.0_dp*t) * (1.0_dp - t)**2
   end function h00

   pure real(dp) function h10(t)
      real(dp), intent(in) :: t
      h10 = t * (1.0_dp - t)**2
   end function h10

   pure real(dp) function h01(t)
      real(dp), intent(in) :: t
      h01 = t**2 * (3.0_dp - 2.0_dp*t)
   end function h01

   pure real(dp) function h11(t)
      real(dp), intent(in) :: t
      h11 = t**2 * (t - 1.0_dp)
   end function h11

end module hermite_interp
