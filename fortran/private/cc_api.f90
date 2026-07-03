! ***********************************************************************
! fortran/private/cc_api.f90  —  NOT COMPILED
!
! This file is a legacy reference copy that routes through colors_lib.f90
! (the public aggregation module).  It is NOT included in the Makefile or
! meson.build build sources and is therefore never compiled.
!
! The file that IS compiled is fortran/cc_api.f90 (one level up), which
! uses the private kernels directly.  That workaround exists because the
! pinned f2py cannot wrap a module that only re-exports procedures from
! other modules (see fortran/cc_api.f90 header and MIGRATION.md).
!
! This copy is kept for reference only.  All changes to the compiled API
! must go into fortran/cc_api.f90, NOT this file.
! ***********************************************************************

module cc_api
   use colors_lib, only: dp, &
      hermite_interp_vector, &
      trilinear_interp_vector, &
      lib_dilute_flux    => dilute_flux, &
      lib_trapz          => trapezoidal_integration, &
      lib_simpson        => simpson_integration, &
      lib_interp_filter  => interp_filter_onto_sed, &
      lib_synthetic_flux => calculate_synthetic_flux, &
      lib_magnitude      => magnitude, &
      lib_bol_phot       => calculate_bolometric_phot, &
      lib_vega_zp        => compute_vega_zero_point, &
      lib_ab_zp          => compute_ab_zero_point, &
      lib_st_zp          => compute_st_zero_point

   implicit none
   private

   public :: interp_sed_hermite
   public :: interp_sed_linear
   public :: dilute_flux
   public :: synthetic_magnitude
   public :: bolometric
   public :: vega_zero_point
   public :: ab_zero_point
   public :: st_zero_point
   public :: trapz
   public :: simpson

contains

   ! -------------------------------------------------------------------------
   ! interp_sed_hermite
   !
   ! Interpolate a full SED at (teff, logg, meta) using cubic Hermite
   ! tensor interpolation on the preloaded flux cube. Forwards straight
   ! to hermite_interp_vector in fortran/private/hermite_interp.f90.
   ! -------------------------------------------------------------------------
   subroutine interp_sed_hermite(teff, logg, meta, &
                                  teff_grid, nt, &
                                  logg_grid, nl, &
                                  meta_grid, nm, &
                                  flux_cube, nw, &
                                  result_flux, ierr)
      !f2py intent(in)  :: teff, logg, meta
      !f2py intent(in)  :: teff_grid, logg_grid, meta_grid
      !f2py intent(in)  :: nt, nl, nm, nw
      !f2py intent(in)  :: flux_cube
      !f2py intent(out) :: result_flux
      !f2py intent(out) :: ierr
      real(dp), intent(in)  :: teff, logg, meta
      integer,  intent(in)  :: nt, nl, nm, nw
      real(dp), intent(in)  :: teff_grid(nt), logg_grid(nl), meta_grid(nm)
      real(dp), intent(in)  :: flux_cube(nt, nl, nm, nw)
      real(dp), intent(out) :: result_flux(nw)
      integer,  intent(out) :: ierr

      call hermite_interp_vector(teff, logg, meta, &
                                 teff_grid, logg_grid, meta_grid, &
                                 flux_cube, nw, result_flux, ierr)
   end subroutine interp_sed_hermite


   ! -------------------------------------------------------------------------
   ! interp_sed_linear
   !
   ! Same interface as interp_sed_hermite but using trilinear interpolation.
   ! -------------------------------------------------------------------------
   subroutine interp_sed_linear(teff, logg, meta, &
                                 teff_grid, nt, &
                                 logg_grid, nl, &
                                 meta_grid, nm, &
                                 flux_cube, nw, &
                                 result_flux, ierr)
      !f2py intent(in)  :: teff, logg, meta
      !f2py intent(in)  :: teff_grid, logg_grid, meta_grid
      !f2py intent(in)  :: nt, nl, nm, nw
      !f2py intent(in)  :: flux_cube
      !f2py intent(out) :: result_flux
      !f2py intent(out) :: ierr
      real(dp), intent(in)  :: teff, logg, meta
      integer,  intent(in)  :: nt, nl, nm, nw
      real(dp), intent(in)  :: teff_grid(nt), logg_grid(nl), meta_grid(nm)
      real(dp), intent(in)  :: flux_cube(nt, nl, nm, nw)
      real(dp), intent(out) :: result_flux(nw)
      integer,  intent(out) :: ierr

      call trilinear_interp_vector(teff, logg, meta, &
                                   teff_grid, logg_grid, meta_grid, &
                                   flux_cube, nw, result_flux, ierr)
   end subroutine interp_sed_linear


   ! -------------------------------------------------------------------------
   ! dilute_flux
   !
   ! Apply (R/d)^2 dilution to convert surface flux to observed flux.
   ! -------------------------------------------------------------------------
   subroutine dilute_flux(surface_flux, nw, R, d, observed_flux)
      !f2py intent(in)  :: surface_flux, nw, R, d
      !f2py intent(out) :: observed_flux
      integer,  intent(in)  :: nw
      real(dp), intent(in)  :: surface_flux(nw), R, d
      real(dp), intent(out) :: observed_flux(nw)
      call lib_dilute_flux(surface_flux, R, d, observed_flux)
   end subroutine dilute_flux


   ! -------------------------------------------------------------------------
   ! synthetic_magnitude
   !
   ! Compute a synthetic magnitude in a single filter from a diluted SED.
   !
   ! Steps performed internally:
   !   1. Interpolate filter transmission onto the SED wavelength grid
   !   2. Photon-counting in-band flux integration
   !   3. m = -2.5 log10(F_band / F_zp)
   !
   ! Parameters
   ! ----------
   ! sed_wave(nw)        : SED wavelength grid (A)
   ! obs_flux(nw)        : diluted (observer-frame) SED flux (erg/s/cm^2/A)
   ! filt_wave(nf)       : filter wavelength grid (A)
   ! filt_trans(nf)      : filter transmission [0,1]
   ! zero_point          : precomputed photometric zero-point
   ! mag                 : output magnitude
   ! band_flux           : output in-band flux (before zero-point)
   ! ierr                : 0=ok, 1=integration failure, 2=non-positive flux
   ! -------------------------------------------------------------------------
   subroutine synthetic_magnitude(sed_wave, obs_flux, nw, &
                                   filt_wave, filt_trans, nf, &
                                   zero_point, &
                                   mag, band_flux, ierr)
      !f2py intent(in)  :: sed_wave, obs_flux, nw
      !f2py intent(in)  :: filt_wave, filt_trans, nf
      !f2py intent(in)  :: zero_point
      !f2py intent(out) :: mag, band_flux, ierr
      integer,  intent(in)  :: nw, nf
      real(dp), intent(in)  :: sed_wave(nw), obs_flux(nw)
      real(dp), intent(in)  :: filt_wave(nf), filt_trans(nf)
      real(dp), intent(in)  :: zero_point
      real(dp), intent(out) :: mag, band_flux
      integer,  intent(out) :: ierr

      real(dp) :: filt_on_sed(nw)
      integer  :: ierr2

      call lib_interp_filter(filt_wave, filt_trans, sed_wave, filt_on_sed, ierr)
      if (ierr /= 0) then
         mag = -99.9_dp; band_flux = -1.0_dp; return
      end if

      call lib_synthetic_flux(sed_wave, obs_flux * filt_on_sed, filt_on_sed, band_flux)
      if (band_flux <= 0.0_dp) then
         mag = -99.9_dp; ierr = 1; return
      end if

      call lib_magnitude(band_flux, zero_point, mag, ierr2)
      if (ierr2 /= 0) ierr = ierr2
   end subroutine synthetic_magnitude


   ! -------------------------------------------------------------------------
   ! bolometric
   !
   ! Compute bolometric flux and magnitude from a diluted SED.
   ! -------------------------------------------------------------------------
   subroutine bolometric(sed_wave, obs_flux, nw, bol_flux, bol_mag, ierr)
      !f2py intent(in)  :: sed_wave, obs_flux, nw
      !f2py intent(out) :: bol_flux, bol_mag, ierr
      integer,  intent(in)  :: nw
      real(dp), intent(in)  :: sed_wave(nw), obs_flux(nw)
      real(dp), intent(out) :: bol_flux, bol_mag
      integer,  intent(out) :: ierr

      call lib_bol_phot(sed_wave, obs_flux, bol_mag, bol_flux, ierr)
   end subroutine bolometric


   ! -------------------------------------------------------------------------
   ! Zero-point subroutines (called once at init from Python)
   !
   ! colors_lib's compute_*_zero_point functions return the MESA-style
   ! function value with -1.0 as the failure sentinel; ierr here is
   ! derived from that sentinel to keep the existing ierr-based Python
   ! contract unchanged.
   ! -------------------------------------------------------------------------

   subroutine vega_zero_point(vega_wave, vega_flux, nv, &
                               filt_wave, filt_trans, nf, &
                               zp, ierr)
      !f2py intent(in)  :: vega_wave, vega_flux, nv
      !f2py intent(in)  :: filt_wave, filt_trans, nf
      !f2py intent(out) :: zp, ierr
      integer,  intent(in)  :: nv, nf
      real(dp), intent(in)  :: vega_wave(nv), vega_flux(nv)
      real(dp), intent(in)  :: filt_wave(nf), filt_trans(nf)
      real(dp), intent(out) :: zp
      integer,  intent(out) :: ierr

      zp = lib_vega_zp(vega_wave, vega_flux, filt_wave, filt_trans)
      ierr = merge(0, 1, zp > 0.0_dp)
   end subroutine vega_zero_point


   subroutine ab_zero_point(filt_wave, filt_trans, nf, zp, ierr)
      !f2py intent(in)  :: filt_wave, filt_trans, nf
      !f2py intent(out) :: zp, ierr
      integer,  intent(in)  :: nf
      real(dp), intent(in)  :: filt_wave(nf), filt_trans(nf)
      real(dp), intent(out) :: zp
      integer,  intent(out) :: ierr

      zp = lib_ab_zp(filt_wave, filt_trans)
      ierr = merge(0, 1, zp > 0.0_dp)
   end subroutine ab_zero_point


   subroutine st_zero_point(filt_wave, filt_trans, nf, zp, ierr)
      !f2py intent(in)  :: filt_wave, filt_trans, nf
      !f2py intent(out) :: zp, ierr
      integer,  intent(in)  :: nf
      real(dp), intent(in)  :: filt_wave(nf), filt_trans(nf)
      real(dp), intent(out) :: zp
      integer,  intent(out) :: ierr

      zp = lib_st_zp(filt_wave, filt_trans)
      ierr = merge(0, 1, zp > 0.0_dp)
   end subroutine st_zero_point


   ! -------------------------------------------------------------------------
   ! Standalone integration (exposed for testing / Python use)
   ! -------------------------------------------------------------------------

   subroutine trapz(x, y, n, result)
      !f2py intent(in)  :: x, y, n
      !f2py intent(out) :: result
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: x(n), y(n)
      real(dp), intent(out) :: result
      call lib_trapz(x, y, result)
   end subroutine trapz


   subroutine simpson(x, y, n, result)
      !f2py intent(in)  :: x, y, n
      !f2py intent(out) :: result
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: x(n), y(n)
      real(dp), intent(out) :: result
      call lib_simpson(x, y, result)
   end subroutine simpson

end module cc_api
