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
! Numerics follow the current MESA kernel: non-uniform three-point
! derivative estimates (weighted by the neighbouring spacings), singleton
! axes treated as degenerate dimensions, polynomial Hermite basis, and no
! positivity floor (negative undershoot is handled by
! hermite_interp_bounded.f90).
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
   use colors_def, only: dp
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
      integer  :: ix_max, iy_max, iz_max
      real(dp) :: t_x, t_y, t_z, dx, dy, dz
      real(dp) :: val, df_dx, df_dy, df_dz
      real(dp) :: h_x(2), hx_d(2), h_y(2), hy_d(2), h_z(2), hz_d(2)
      real(dp) :: wx, wy, wz, wxd, wyd, wzd

      nx = size(x_grid); ny = size(y_grid); nz = size(z_grid)
      if (present(ierr)) ierr = 0

      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      ! Singleton axes are valid degenerate dimensions, not an out-of-grid
      ! condition.
      ix_max = 1; iy_max = 1; iz_max = 1
      if (nx == 1) then
         i_x = 1; t_x = 0.0_dp; ix_max = 0
      end if
      if (ny == 1) then
         i_y = 1; t_y = 0.0_dp; iy_max = 0
      end if
      if (nz == 1) then
         i_z = 1; t_z = 0.0_dp; iz_max = 0
      end if

      ! Outside a non-degenerate axis: fall back to the nearest node
      if ((nx > 1 .and. (i_x < 1 .or. i_x >= nx)) .or. &
          (ny > 1 .and. (i_y < 1 .or. i_y >= ny)) .or. &
          (nz > 1 .and. (i_z < 1 .or. i_z >= nz))) then
         call find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
         result_flux = f_values_4d(i_x, i_y, i_z, :)
         if (present(ierr)) ierr = 1   ! flag: clamped to boundary
         return
      end if

      ! Cell spacing; degenerate axes have no derivative contribution.
      dx = 0.0_dp; dy = 0.0_dp; dz = 0.0_dp
      if (nx > 1) dx = x_grid(i_x + 1) - x_grid(i_x)
      if (ny > 1) dy = y_grid(i_y + 1) - y_grid(i_y)
      if (nz > 1) dz = z_grid(i_z + 1) - z_grid(i_z)

      ! Hermite basis (same for all wavelengths)
      h_x  = [h00(t_x), h01(t_x)]
      hx_d = [h10(t_x), h11(t_x)]
      h_y  = [h00(t_y), h01(t_y)]
      hy_d = [h10(t_y), h11(t_y)]
      h_z  = [h00(t_z), h01(t_z)]
      hz_d = [h10(t_z), h11(t_z)]

      ! Stencil loop: the weights are invariant over lambda, so lambda is
      ! innermost.
      result_flux = 0.0_dp
      do iz = 0, iz_max
         wz = h_z(iz+1); wzd = hz_d(iz+1)
         do iy = 0, iy_max
            wy = h_y(iy+1); wyd = hy_d(iy+1)
            do ix = 0, ix_max
               wx = h_x(ix+1); wxd = hx_d(ix+1)
               do lam = 1, n_lambda
                  val = f_values_4d(i_x+ix, i_y+iy, i_z+iz, lam)

                  call compute_derivatives_4d( &
                     f_values_4d, nx, ny, nz, n_lambda, &
                     i_x+ix, i_y+iy, i_z+iz, lam, &
                     x_grid, y_grid, z_grid, &
                     df_dx, df_dy, df_dz)

                  result_flux(lam) = result_flux(lam) &
                     + wx*wy*wz*val &
                     + wxd*wy*wz*dx*df_dx &
                     + wx*wyd*wz*dy*df_dy &
                     + wx*wy*wzd*dz*df_dz
               end do
            end do
         end do
      end do
   end subroutine hermite_interp_vector


   ! Derivative estimates on a non-uniform grid.  Interior nodes use the
   ! spacing-weighted average of the one-sided slopes; boundary nodes use a
   ! one-sided difference; singleton axes have zero derivative.
   subroutine compute_derivatives_4d(f, nx, ny, nz, n_lambda, &
                                     i, j, k, lam, &
                                     x_grid, y_grid, z_grid, &
                                     df_dx, df_dy, df_dz)
      integer,  intent(in)  :: nx, ny, nz, n_lambda
      real(dp), intent(in)  :: f(nx, ny, nz, n_lambda)
      integer,  intent(in)  :: i, j, k, lam
      real(dp), intent(in)  :: x_grid(nx), y_grid(ny), z_grid(nz)
      real(dp), intent(out) :: df_dx, df_dy, df_dz

      real(dp) :: h_minus, h_plus, d_minus, d_plus

      if (nx == 1) then
         df_dx = 0.0_dp
      else if (i > 1 .and. i < nx) then
         h_minus = x_grid(i)   - x_grid(i-1)
         h_plus  = x_grid(i+1) - x_grid(i)
         d_minus = (f(i,j,k,lam)   - f(i-1,j,k,lam)) / h_minus
         d_plus  = (f(i+1,j,k,lam) - f(i,j,k,lam))   / h_plus
         df_dx = (h_plus*d_minus + h_minus*d_plus) / (h_minus + h_plus)
      else if (i == 1) then
         df_dx = (f(i+1,j,k,lam) - f(i,j,k,lam)) / (x_grid(i+1) - x_grid(i))
      else
         df_dx = (f(i,j,k,lam) - f(i-1,j,k,lam)) / (x_grid(i) - x_grid(i-1))
      end if

      if (ny == 1) then
         df_dy = 0.0_dp
      else if (j > 1 .and. j < ny) then
         h_minus = y_grid(j)   - y_grid(j-1)
         h_plus  = y_grid(j+1) - y_grid(j)
         d_minus = (f(i,j,k,lam)   - f(i,j-1,k,lam)) / h_minus
         d_plus  = (f(i,j+1,k,lam) - f(i,j,k,lam))   / h_plus
         df_dy = (h_plus*d_minus + h_minus*d_plus) / (h_minus + h_plus)
      else if (j == 1) then
         df_dy = (f(i,j+1,k,lam) - f(i,j,k,lam)) / (y_grid(j+1) - y_grid(j))
      else
         df_dy = (f(i,j,k,lam) - f(i,j-1,k,lam)) / (y_grid(j) - y_grid(j-1))
      end if

      if (nz == 1) then
         df_dz = 0.0_dp
      else if (k > 1 .and. k < nz) then
         h_minus = z_grid(k)   - z_grid(k-1)
         h_plus  = z_grid(k+1) - z_grid(k)
         d_minus = (f(i,j,k,lam)   - f(i,j,k-1,lam)) / h_minus
         d_plus  = (f(i,j,k+1,lam) - f(i,j,k,lam))   / h_plus
         df_dz = (h_plus*d_minus + h_minus*d_plus) / (h_minus + h_plus)
      else if (k == 1) then
         df_dz = (f(i,j,k+1,lam) - f(i,j,k,lam)) / (z_grid(k+1) - z_grid(k))
      else
         df_dz = (f(i,j,k,lam) - f(i,j,k-1,lam)) / (z_grid(k) - z_grid(k-1))
      end if
   end subroutine compute_derivatives_4d


   ! Hermite basis functions
   pure real(dp) function h00(t)
      real(dp), intent(in) :: t
      h00 = 2.0_dp*t**3 - 3.0_dp*t**2 + 1.0_dp
   end function h00

   pure real(dp) function h10(t)
      real(dp), intent(in) :: t
      h10 = t**3 - 2.0_dp*t**2 + t
   end function h10

   pure real(dp) function h01(t)
      real(dp), intent(in) :: t
      h01 = -2.0_dp*t**3 + 3.0_dp*t**2
   end function h01

   pure real(dp) function h11(t)
      real(dp), intent(in) :: t
      h11 = t**3 - t**2
   end function h11

end module hermite_interp
