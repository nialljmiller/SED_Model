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
   use colors_utils, only: find_containing_cell, find_nearest_point
   implicit none

   private
   public :: trilinear_interp_vector

contains

   subroutine trilinear_interp_vector(x_val, y_val, z_val, &
                                      x_grid, y_grid, z_grid, &
                                      f_values_4d, n_lambda, result_flux, ierr)
      real(dp), intent(in)  :: x_val, y_val, z_val
      real(dp), intent(in)  :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in)  :: f_values_4d(:,:,:,:)   ! (nx, ny, nz, n_lambda)
      integer,  intent(in)  :: n_lambda
      real(dp), intent(out) :: result_flux(n_lambda)
      integer,  intent(out), optional :: ierr

      integer  :: nx, ny, nz
      integer  :: i_x, i_y, i_z, lam
      real(dp) :: t_x, t_y, t_z
      real(dp) :: c000, c001, c010, c011, c100, c101, c110, c111
      real(dp) :: c00, c01, c10, c11, c0, c1, interp_val
      integer  :: ierr_local

      nx = size(x_grid); ny = size(y_grid); nz = size(z_grid)
      ierr_local = 0

      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      if (i_x < 1 .or. i_x >= nx .or. &
          i_y < 1 .or. i_y >= ny .or. &
          i_z < 1 .or. i_z >= nz) then
         call find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
         result_flux = f_values_4d(i_x, i_y, i_z, :)
         ierr_local = 1
         if (present(ierr)) ierr = ierr_local
         return
      end if

      do lam = 1, n_lambda
         c000 = f_values_4d(i_x,   i_y,   i_z,   lam)
         c100 = f_values_4d(i_x+1, i_y,   i_z,   lam)
         c010 = f_values_4d(i_x,   i_y+1, i_z,   lam)
         c110 = f_values_4d(i_x+1, i_y+1, i_z,   lam)
         c001 = f_values_4d(i_x,   i_y,   i_z+1, lam)
         c101 = f_values_4d(i_x+1, i_y,   i_z+1, lam)
         c011 = f_values_4d(i_x,   i_y+1, i_z+1, lam)
         c111 = f_values_4d(i_x+1, i_y+1, i_z+1, lam)

         c00 = c000*(1.0_dp - t_x) + c100*t_x
         c01 = c001*(1.0_dp - t_x) + c101*t_x
         c10 = c010*(1.0_dp - t_x) + c110*t_x
         c11 = c011*(1.0_dp - t_x) + c111*t_x

         c0 = c00*(1.0_dp - t_y) + c10*t_y
         c1 = c01*(1.0_dp - t_y) + c11*t_y

         interp_val = c0*(1.0_dp - t_z) + c1*t_z
         result_flux(lam) = max(TINY_VALUE, interp_val)
      end do

      if (present(ierr)) ierr = ierr_local
   end subroutine trilinear_interp_vector

end module linear_interp
