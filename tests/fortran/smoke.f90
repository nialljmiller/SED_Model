! Smoke test for the standalone Fortran kernels.  Reproduces the cases from
! the MESA colors unit test (colors/test/src/test_colors.f90) that do not
! need a flux cube.  Exits non-zero on the first failure.
program fortran_smoke
   use colors_lib
   implicit none

   integer :: nfail
   nfail = 0

   call check_simpson_nonuniform()
   call check_filter_support()
   call check_nonuniform_hermite()
   call check_bounded_fallback()
   call check_bolometric_zero_point()

   if (nfail > 0) then
      write (*, '(a,i0,a)') 'FAILED: ', nfail, ' check(s)'
      stop 1
   end if
   write (*, '(a)') 'all Fortran kernel checks passed'

contains

   subroutine expect(label, got, want, tol)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: got, want, tol
      if (abs(got - want) > tol*max(1.0_dp, abs(want))) then
         write (*, '(a,": got ",es24.16," want ",es24.16)') label, got, want
         nfail = nfail + 1
      end if
   end subroutine expect

   subroutine check_simpson_nonuniform()
      real(dp) :: x(5), r
      x = [0.0_dp, 1.0_dp, 3.0_dp, 4.0_dp, 7.0_dp]
      call simpson_integration(x, x*x, r)
      call expect('simpson_nonuniform', r, 343.0_dp/3.0_dp, 1.0d-12)
   end subroutine check_simpson_nonuniform

   subroutine check_filter_support()
      real(dp) :: sed_wave(5), sed_flux(5), filt_wave(3), trans(3), conv(5), on_sed(5)
      real(dp) :: band, zp
      integer :: ierr
      sed_wave = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      sed_flux = [100.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 100.0_dp]
      filt_wave = [2.0_dp, 3.0_dp, 4.0_dp]
      trans = 1.0_dp
      call interp_filter_onto_sed(filt_wave, trans, sed_wave, on_sed, ierr)
      conv = sed_flux*on_sed
      call calculate_synthetic_flux(sed_wave, conv, on_sed, band)
      zp = compute_vega_zero_point(sed_wave, sed_flux, filt_wave, trans)
      call expect('filter_support_band_flux', band, 2.0_dp, 1.0d-12)
      call expect('filter_support_vega_zp', zp, 2.0_dp, 1.0d-12)
   end subroutine check_filter_support

   subroutine check_nonuniform_hermite()
      real(dp) :: x(4), one(1), cube(4, 1, 1, 3), flux(3)
      integer :: lam
      x = [0.0_dp, 1.0_dp, 3.0_dp, 6.0_dp]
      one = 0.0_dp
      do lam = 1, 3
         cube(:, 1, 1, lam) = 100.0_dp + x*x
      end do
      call hermite_interp_vector(1.5_dp, 0.0_dp, 0.0_dp, x, one, one, cube, 3, flux)
      call expect('nonuniform_hermite', flux(1), 102.25_dp, 1.0d-12)
   end subroutine check_nonuniform_hermite

   subroutine check_bounded_fallback()
      real(dp) :: x(3), one(1), wave(3), cube(3, 1, 1, 3), flux(3)
      integer :: lam, ierr
      x = [0.0_dp, 1.0_dp, 2.0_dp]
      one = 0.0_dp
      wave = [1.0_dp, 2.0_dp, 4.0_dp]
      do lam = 1, 3
         cube(:, 1, 1, lam) = [1.0_dp, 1.0_dp, 100.0_dp]
      end do
      call hermite_bounded_interp_vector(0.5_dp, 0.0_dp, 0.0_dp, x, one, one, &
                                         cube, 3, wave, flux, ierr)
      call expect('bounded_fallback_flux', flux(1), 1.0_dp, 1.0d-12)
      call expect('bounded_fallback_ierr', real(ierr, dp), 1.0_dp, 0.0_dp)
   end subroutine check_bounded_fallback

   subroutine check_bolometric_zero_point()
      real(dp), parameter :: lsun = 3.828d33, pc = 3.0856775814913673d18
      real(dp), parameter :: pi = 3.141592653589793_dp
      real(dp) :: wave(3), flux(3), solar, mag, bol_flux
      integer :: ierr
      solar = lsun/(4.0_dp*pi*(10.0_dp*pc)**2)
      wave = [1.0_dp, 2.0_dp, 3.0_dp]
      flux = solar/2.0_dp
      call calculate_bolometric_phot(wave, flux, mag, bol_flux, ierr)
      call expect('bolometric_zero_point', mag, 4.74_dp, 1.0d-10)
   end subroutine check_bolometric_zero_point

end program fortran_smoke
