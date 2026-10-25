"""
tests/test_inverse.py
=====================
Tests for sed_model.inverse (photometry -> stellar parameters).

Strategy
--------
Rather than testing against external reference values, we test the
inverse model self-consistently:

  1. Run the forward model at known parameters to produce synthetic
     magnitudes.
  2. Perturb those magnitudes by a small known noise level.
  3. Run the inverse model and verify the posterior median recovers
     the true parameters within a physically reasonable tolerance.

MCMC chains are kept deliberately short (n_steps=300, n_burn=100,
n_walkers=16) so the test suite completes in under ~30 seconds.
The tolerances are set wide enough to accommodate shot noise from the
short chain while still being meaningful.

All tests that require the Fortran extension or data files skip
gracefully when those are not present.
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
from _data_paths import KURUCZ_DIR, JOHNSON_DIR, VEGA_SED

R_SUN  = 6.957e10
D_10PC = 3.0857e19

# ---------------------------------------------------------------------------
# Skip guards
# ---------------------------------------------------------------------------

def _cc_api_available():
    try:
        from sed_model import cc_api  # noqa: F401
        return True
    except ImportError:
        return False

requires_fortran = pytest.mark.skipif(
    not _cc_api_available(),
    reason="cc_api Fortran extension not built",
)
requires_grid = pytest.mark.skipif(
    not os.path.isdir(KURUCZ_DIR),
    reason=f"Kurucz2003all grid not found at {KURUCZ_DIR}",
)
requires_filters = pytest.mark.skipif(
    not os.path.isdir(JOHNSON_DIR),
    reason=f"Johnson filter directory not found at {JOHNSON_DIR}",
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def grid():
    from sed_model import load_grid
    return load_grid(KURUCZ_DIR)


@pytest.fixture(scope="module")
def filters():
    from sed_model import load_filters_from_instrument_dir
    vega = VEGA_SED if os.path.isfile(VEGA_SED) else None
    return load_filters_from_instrument_dir(JOHNSON_DIR, vega_sed_path=vega)


@pytest.fixture(scope="module")
def solar_obs(grid, filters):
    """Synthetic solar observations with 0.01 mag Gaussian noise."""
    from sed_model import run_forward
    TRUE_TEFF = 5778.0
    TRUE_LOGG = 4.44
    TRUE_META = 0.0

    result = run_forward(
        teff=TRUE_TEFF, logg=TRUE_LOGG, meta=TRUE_META,
        R=R_SUN, d=D_10PC,
        grid=grid, filters=filters,
        mag_system="AB",
    )

    rng = np.random.default_rng(42)
    filter_names = [f.name for f in filters]
    obs_mag = np.array([result.magnitudes[n] for n in filter_names])
    obs_err = np.full(len(filter_names), 0.01)
    obs_mag += rng.normal(0, obs_err)

    return {
        "true_teff": TRUE_TEFF,
        "true_logg": TRUE_LOGG,
        "true_meta": TRUE_META,
        "obs_magnitudes":    obs_mag,
        "obs_uncertainties": obs_err,
        "filter_names":      filter_names,
    }


# ---------------------------------------------------------------------------
# Smoke test — does run_inverse return an InverseResult?
# ---------------------------------------------------------------------------

@requires_fortran
@requires_grid
@requires_filters
class TestInverseSmoke:

    def test_returns_inverse_result(self, grid, filters, solar_obs):
        from sed_model import run_inverse, InverseResult
        result = run_inverse(
            obs_magnitudes=solar_obs["obs_magnitudes"],
            obs_uncertainties=solar_obs["obs_uncertainties"],
            filter_names=solar_obs["filter_names"],
            R=R_SUN, d=D_10PC,
            grid=grid, filters=filters,
            mag_system="AB",
            n_walkers=16, n_steps=200, n_burn=50,
            progress=False, seed=0,
        )
        assert isinstance(result, InverseResult)

    def test_samples_shape(self, grid, filters, solar_obs):
        from sed_model import run_inverse
        result = run_inverse(
            obs_magnitudes=solar_obs["obs_magnitudes"],
            obs_uncertainties=solar_obs["obs_uncertainties"],
            filter_names=solar_obs["filter_names"],
            R=R_SUN, d=D_10PC,
            grid=grid, filters=filters,
            mag_system="AB",
            n_walkers=16, n_steps=200, n_burn=50,
            n_thin=1, progress=False, seed=0,
        )
        # (n_steps - n_burn) * n_walkers samples, 3 parameters
        expected_samples = (200 - 50) * 16
        assert result.samples.shape == (expected_samples, 3)
        assert result.log_prob.shape == (expected_samples,)

    def test_acceptance_fraction_reasonable(self, grid, filters, solar_obs):
        from sed_model import run_inverse
        result = run_inverse(
            obs_magnitudes=solar_obs["obs_magnitudes"],
            obs_uncertainties=solar_obs["obs_uncertainties"],
            filter_names=solar_obs["filter_names"],
            R=R_SUN, d=D_10PC,
            grid=grid, filters=filters,
            mag_system="AB",
            n_walkers=16, n_steps=300, n_burn=100,
            progress=False, seed=0,
        )
        mean_af = float(np.mean(result.acceptance_fraction))
        assert 0.05 < mean_af < 0.99, (
            f"Suspicious acceptance fraction: {mean_af:.3f}"
        )


# ---------------------------------------------------------------------------
# Posterior recovery — do we get back roughly the right parameters?
# ---------------------------------------------------------------------------

@requires_fortran
@requires_grid
@requires_filters
class TestPosteriorRecovery:
    """
    Uses solar parameters which are well inside the grid.
    With 0.01 mag per-filter noise and 7 Johnson filters,
    the posterior should be narrow and centred near the truth.
    Tolerances are deliberately generous to accommodate short chains.
    """

    # These are physically motivated tolerances, not tight numerical ones.
    TEFF_TOL = 500.0    # K — grid spacing is ~250 K so 2 grid steps
    LOGG_TOL = 0.5      # dex — grid spacing is 0.5 dex
    META_TOL = 0.5      # dex — grid spacing is ~0.375 dex

    @pytest.fixture(scope="class")
    def posterior(self, grid, filters, solar_obs):
        from sed_model import run_inverse
        return run_inverse(
            obs_magnitudes=solar_obs["obs_magnitudes"],
            obs_uncertainties=solar_obs["obs_uncertainties"],
            filter_names=solar_obs["filter_names"],
            R=R_SUN, d=D_10PC,
            grid=grid, filters=filters,
            mag_system="AB",
            n_walkers=16, n_steps=500, n_burn=150,
            p0_centre={
                'teff': solar_obs["true_teff"],
                'logg': solar_obs["true_logg"],
                'meta': solar_obs["true_meta"],
            },
            p0_scatter=0.03,
            progress=False, seed=42,
        )

    def test_teff_recovered(self, posterior, solar_obs):
        median_teff = float(np.median(posterior.samples[:, 0]))
        err = abs(median_teff - solar_obs["true_teff"])
        assert err < self.TEFF_TOL, (
            f"Teff: median={median_teff:.1f}, true={solar_obs['true_teff']:.1f}, "
            f"err={err:.1f} K > tol={self.TEFF_TOL} K"
        )

    def test_logg_recovered(self, posterior, solar_obs):
        median_logg = float(np.median(posterior.samples[:, 1]))
        err = abs(median_logg - solar_obs["true_logg"])
        assert err < self.LOGG_TOL, (
            f"logg: median={median_logg:.3f}, true={solar_obs['true_logg']:.3f}, "
            f"err={err:.3f} > tol={self.LOGG_TOL}"
        )

    def test_meta_recovered(self, posterior, solar_obs):
        median_meta = float(np.median(posterior.samples[:, 2]))
        err = abs(median_meta - solar_obs["true_meta"])
        assert err < self.META_TOL, (
            f"[M/H]: median={median_meta:.3f}, true={solar_obs['true_meta']:.3f}, "
            f"err={err:.3f} > tol={self.META_TOL}"
        )

    def test_samples_within_grid(self, posterior, grid):
        """All posterior samples must lie within grid bounds."""
        t = posterior.samples[:, 0]
        l = posterior.samples[:, 1]
        m = posterior.samples[:, 2]
        assert np.all(t >= grid.teff_bounds[0]) and np.all(t <= grid.teff_bounds[1])
        assert np.all(l >= grid.logg_bounds[0]) and np.all(l <= grid.logg_bounds[1])
        assert np.all(m >= grid.meta_bounds[0]) and np.all(m <= grid.meta_bounds[1])

    def test_log_prob_finite(self, posterior):
        assert np.all(np.isfinite(posterior.log_prob))

    def test_posterior_not_degenerate(self, posterior):
        """Posterior should have non-zero spread in all three parameters."""
        for i, name in enumerate(["Teff", "logg", "[M/H]"]):
            std = float(np.std(posterior.samples[:, i]))
            assert std > 0, f"Posterior for {name} is a delta function (std=0)"


# ---------------------------------------------------------------------------
# InverseResult API
# ---------------------------------------------------------------------------

@requires_fortran
@requires_grid
@requires_filters
class TestInverseResultAPI:

    @pytest.fixture(scope="class")
    def result(self, grid, filters, solar_obs):
        from sed_model import run_inverse
        return run_inverse(
            obs_magnitudes=solar_obs["obs_magnitudes"],
            obs_uncertainties=solar_obs["obs_uncertainties"],
            filter_names=solar_obs["filter_names"],
            R=R_SUN, d=D_10PC,
            grid=grid, filters=filters,
            mag_system="AB",
            n_walkers=16, n_steps=200, n_burn=50,
            progress=False, seed=0,
        )

    def test_summary_keys(self, result):
        s = result.summary()
        for key in ("teff", "logg", "meta"):
            assert key in s
            for sub in ("median", "lo", "hi", "lower_1sigma", "upper_1sigma"):
                assert sub in s[key]

    def test_summary_ordering(self, result):
        s = result.summary()
        for key in ("teff", "logg", "meta"):
            assert s[key]["lo"] <= s[key]["median"] <= s[key]["hi"]

    def test_map_estimate_returns_finite_params(self, result):
        map_params = result.map_estimate()
        assert set(map_params) == {"teff", "logg", "meta"}
        for name, value in map_params.items():
            assert np.isfinite(value), f"MAP {name} is not finite"

    def test_print_summary_runs(self, result):
        result.print_summary()  # should not raise

    def test_save_and_load_roundtrip(self, result):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "posterior.npz")
            result.save(path)
            assert os.path.isfile(path)

            from sed_model import InverseResult
            loaded = InverseResult.load(path)

        np.testing.assert_array_equal(result.samples, loaded.samples)
        np.testing.assert_array_equal(result.log_prob, loaded.log_prob)
        np.testing.assert_array_equal(result.obs_magnitudes, loaded.obs_magnitudes)
        np.testing.assert_array_equal(result.obs_uncertainties, loaded.obs_uncertainties)
        assert result.filter_names == loaded.filter_names
        assert result.R == loaded.R
        assert result.d == loaded.d
        assert result.mag_system == loaded.mag_system
        assert result.n_walkers == loaded.n_walkers
        assert result.n_steps == loaded.n_steps
        assert result.n_burn == loaded.n_burn
        assert result.n_thin == loaded.n_thin

    def test_save_adds_npz_extension(self, result):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "posterior")  # no extension
            result.save(path)
            assert os.path.isfile(path + ".npz")

    def test_thinning_reduces_samples(self, grid, filters, solar_obs):
        from sed_model import run_inverse
        r1 = run_inverse(
            obs_magnitudes=solar_obs["obs_magnitudes"],
            obs_uncertainties=solar_obs["obs_uncertainties"],
            filter_names=solar_obs["filter_names"],
            R=R_SUN, d=D_10PC, grid=grid, filters=filters,
            n_walkers=16, n_steps=200, n_burn=50, n_thin=1,
            progress=False, seed=0,
        )
        r5 = run_inverse(
            obs_magnitudes=solar_obs["obs_magnitudes"],
            obs_uncertainties=solar_obs["obs_uncertainties"],
            filter_names=solar_obs["filter_names"],
            R=R_SUN, d=D_10PC, grid=grid, filters=filters,
            n_walkers=16, n_steps=200, n_burn=50, n_thin=5,
            progress=False, seed=0,
        )
        assert r5.samples.shape[0] < r1.samples.shape[0]


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

@requires_fortran
@requires_grid
@requires_filters
class TestInverseValidation:

    def test_mismatched_lengths_raises(self, grid, filters):
        from sed_model import run_inverse
        filter_names = [filters[0].name]
        with pytest.raises(ValueError, match="filter_names"):
            run_inverse(
                obs_magnitudes=[5.0, 6.0],   # 2 mags
                obs_uncertainties=[0.01, 0.01],
                filter_names=filter_names,   # 1 name
                R=R_SUN, d=D_10PC,
                grid=grid, filters=filters,
                n_walkers=8, n_steps=10, n_burn=2,
                progress=False,
            )

    def test_zero_uncertainty_raises(self, grid, filters):
        from sed_model import run_inverse
        with pytest.raises(ValueError, match="positive"):
            run_inverse(
                obs_magnitudes=[5.0],
                obs_uncertainties=[0.0],
                filter_names=[filters[0].name],
                R=R_SUN, d=D_10PC,
                grid=grid, filters=filters,
                n_walkers=8, n_steps=10, n_burn=2,
                progress=False,
            )

    def test_unknown_filter_name_raises(self, grid, filters):
        from sed_model import run_inverse
        with pytest.raises(ValueError, match="not in the supplied filters"):
            run_inverse(
                obs_magnitudes=[5.0],
                obs_uncertainties=[0.01],
                filter_names=["NonExistentFilter_XYZ"],
                R=R_SUN, d=D_10PC,
                grid=grid, filters=filters,
                n_walkers=8, n_steps=10, n_burn=2,
                progress=False,
            )


# ---------------------------------------------------------------------------
# Post-0.3.1 inverse improvements: fast model, global init, systematics,
# priors, upper limits, diagnostics
# ---------------------------------------------------------------------------

def _mags(grid, filters, teff, logg, meta, d=D_10PC, av=0.0):
    from sed_model import run_forward, make_extinction_model
    ext = make_extinction_model(enabled=True, law="fitzpatrick99", a_v=av) if av else None
    r = run_forward(teff=teff, logg=logg, meta=meta, R=R_SUN, d=d, grid=grid,
                    filters=filters, mag_system="AB", extinction=ext)
    names = [f.name for f in filters]
    return names, np.array([r.magnitudes[n] for n in names])


@requires_fortran
@requires_grid
@requires_filters
class TestInverseCore:

    def test_fast_model_matches_run_forward(self, grid, filters):
        from sed_model import make_extinction_model
        from sed_model._inv_core import _FastModel
        ext = make_extinction_model(enabled=True, law="fitzpatrick99", a_v=0.0)
        names = [f.name for f in filters]
        fm = _FastModel(grid, filters, names, R_SUN, "AB", "hermite_bounded", ext)
        for (t, l, m, av) in [(5800.0, 4.4, 0.0, 0.0), (7123.0, 3.7, -0.8, 1.3), (3650.0, 4.9, 0.45, 0.4)]:
            _, ref = _mags(grid, filters, t, l, m, av=av)
            np.testing.assert_allclose(fm.mags(t, l, m, av, D_10PC), ref, atol=1e-9, rtol=0)

    def test_grid_holes_are_unsupported(self, grid, filters):
        from sed_model._inv_core import _FastModel
        fm = _FastModel(grid, filters, [filters[0].name], R_SUN, "AB", "hermite_bounded", None)
        assert fm.supported(5800.0, 4.4, 0.0)
        # hot + low gravity corner of the Kurucz grid is an all-zero SED
        assert not fm.supported(40000.0, 1.0, 0.0)

    def test_auto_init_recovers_hot_star_from_far_away(self, grid, filters):
        """The legacy mid-box ball starts at ~27000 K; the global search must not care."""
        from sed_model import run_inverse
        names, m = _mags(grid, filters, 9500.0, 4.2, -0.4)
        res = run_inverse(m, np.full(m.size, 0.02), names, R=R_SUN, d=D_10PC,
                          grid=grid, filters=filters, n_walkers=16, n_steps=400,
                          n_burn=100, progress=False, verbose=False, seed=3)
        s = res.summary()
        assert abs(s["teff"]["median"] - 9500.0) < 400.0
        assert res.init_info["method"] == "global-search+LM"
        assert res.converged is not None and res.n_eff is not None
        assert res.chi2_best is not None and res.n_data == m.size

    def test_legacy_ball_path_still_works(self, grid, filters):
        from sed_model import run_inverse
        names, m = _mags(grid, filters, 5800.0, 4.4, 0.0)
        res = run_inverse(m, np.full(m.size, 0.02), names, R=R_SUN, d=D_10PC,
                          grid=grid, filters=filters, n_walkers=16, n_steps=200,
                          n_burn=50, progress=False, verbose=False, seed=1,
                          init="ball", p0_centre={"teff": 5800.0, "logg": 4.4, "meta": 0.0})
        assert res.init_info["method"] == "ball"

    def test_nan_band_is_dropped_with_warning(self, grid, filters):
        from sed_model import run_inverse
        names, m = _mags(grid, filters, 5800.0, 4.4, 0.0)
        m = m.copy()
        m[2] = np.nan
        with pytest.warns(UserWarning, match="non-finite"):
            res = run_inverse(m, np.full(m.size, 0.02), names, R=R_SUN, d=D_10PC,
                              grid=grid, filters=filters, n_walkers=16, n_steps=100,
                              n_burn=20, progress=False, verbose=False, seed=1)
        assert len(res.filter_names) == len(names) - 1

    def test_fit_jitter_adds_parameter(self, grid, filters):
        from sed_model import run_inverse
        names, m = _mags(grid, filters, 5800.0, 4.4, 0.0)
        res = run_inverse(m, np.full(m.size, 0.01), names, R=R_SUN, d=D_10PC,
                          grid=grid, filters=filters, n_walkers=16, n_steps=150,
                          n_burn=50, progress=False, verbose=False, seed=1,
                          fit_jitter=True, jitter_max=0.1)
        assert res.param_names == ["teff", "logg", "meta", "sigma_sys"]
        assert res.samples.shape[1] == 4
        assert np.all((res.samples[:, 3] >= 0) & (res.samples[:, 3] <= 0.1))

    def test_sys_floor_and_upper_limit_likelihood(self, grid, filters):
        from sed_model import fit_params_from_grid
        from sed_model._inv_core import _FastModel, _Problem
        names, m = _mags(grid, filters, 5800.0, 4.4, 0.0)
        fp = fit_params_from_grid(grid, d_cm=D_10PC)
        fm = _FastModel(grid, filters, names, R_SUN, "AB", "hermite_bounded", None)
        th = np.array([5800.0, 4.4, 0.0])
        obs = m + 0.1
        p0 = _Problem(fm, fp, obs, np.full(m.size, 0.01))
        p1 = _Problem(fm, fp, obs, np.full(m.size, 0.01), sys_floor=0.1)
        assert p1.log_like(th) > p0.log_like(th)          # floor tolerates the offset
        up = np.zeros(m.size, bool)
        up[0] = True
        # band 0 is a limit; the model is fainter than the limit (ok) or brighter (violated)
        obs_ok, obs_bad = m.copy(), m.copy()
        obs_ok[0] -= 0.5
        obs_bad[0] += 0.5
        ok = _Problem(fm, fp, obs_ok, np.full(m.size, 0.01), upper=up)
        bad = _Problem(fm, fp, obs_bad, np.full(m.size, 0.01), upper=up)
        assert ok.log_like(th) > bad.log_like(th) + 100.0

    def test_prior_parsing(self):
        from sed_model._inv_core import parse_priors
        p = parse_priors({"d": ("parallax", 2.0, 0.1), "a_v": ("halfnormal", 0.5)}, ["teff", "a_v", "d"])
        pc = 3.085677581491367e18
        assert p["d"](500.0 * pc) == pytest.approx(0.0)       # parallax exactly 2 mas
        assert p["d"](250.0 * pc) < -100.0
        assert p["a_v"](0.0) > p["a_v"](1.0)
        with pytest.raises(ValueError):
            parse_priors({"logg": ("normal", 4.0, 0.1)}, ["teff"])

    def test_parallax_prior_constrains_distance(self, grid, filters):
        from sed_model import run_inverse, fit_params_from_grid
        pc = 3.085677581491367e18
        names, m = _mags(grid, filters, 5800.0, 4.4, 0.0, d=100 * pc)
        fp = fit_params_from_grid(grid, teff=(5000, 6500), logg=(3.5, 5.0), meta=(-1.0, 0.5),
                                  d_cm=(20 * pc, 500 * pc))
        res = run_inverse(m, np.full(m.size, 0.02), names, R=R_SUN, grid=grid, filters=filters,
                          fit_params=fp, n_walkers=16, n_steps=400, n_burn=150,
                          progress=False, verbose=False, seed=2,
                          priors={"d": ("parallax", 10.0, 0.2)})
        dmed = np.median(res.samples[:, 3]) / pc
        assert abs(dmed - 100.0) < 5.0

    def test_auto_stop_and_auto_burn(self, grid, filters):
        from sed_model import run_inverse
        names, m = _mags(grid, filters, 5800.0, 4.4, 0.0)
        res = run_inverse(m, np.full(m.size, 0.02), names, R=R_SUN, d=D_10PC,
                          grid=grid, filters=filters, n_walkers=16, n_steps=3000,
                          n_burn=None, auto_stop=True, progress=False, verbose=False, seed=5)
        assert res.n_steps <= 3000
        assert 0 <= res.n_burn < res.n_steps

    def test_result_roundtrip_new_fields(self, grid, filters):
        from sed_model import run_inverse, InverseResult
        names, m = _mags(grid, filters, 5800.0, 4.4, 0.0)
        res = run_inverse(m, np.full(m.size, 0.02), names, R=R_SUN, d=D_10PC,
                          grid=grid, filters=filters, n_walkers=16, n_steps=120,
                          n_burn=30, progress=False, verbose=False, seed=5)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "r.npz")
            res.save(path)
            back = InverseResult.load(path)
        assert back.converged == res.converged
        assert back.chi2_best == pytest.approx(res.chi2_best)
        assert back.n_data == res.n_data
        np.testing.assert_allclose(back.n_eff, res.n_eff)
        assert back.init_info["method"] == res.init_info["method"]

    def test_truth_on_prior_bounds_does_not_break_sampler(self, grid, filters):
        """Truth pinned to a grid corner: init walkers must be distinct and finite
        (identical walkers make differential-evolution proposals NaN)."""
        from sed_model import run_inverse
        names, m = _mags(grid, filters, 3600.0, 5.0, 0.5)
        res = run_inverse(m, np.full(m.size, 0.005), names, R=R_SUN, d=D_10PC,
                          grid=grid, filters=filters, n_walkers=16, n_steps=200,
                          n_burn=50, progress=False, verbose=False, seed=11)
        assert np.all(np.isfinite(res.samples)) and np.all(np.isfinite(res.log_prob))
        assert np.median(res.samples[:, 1]) > 4.0

    def test_multimodal_init_distinct_walkers(self, grid, filters):
        from sed_model import fit_params_from_grid
        from sed_model._inv_core import _FastModel, _Problem, initial_walkers
        names, m = _mags(grid, filters, 6000.0, 4.5, -0.3)
        fp = fit_params_from_grid(grid, d_cm=D_10PC)
        fm = _FastModel(grid, filters, names, R_SUN, "AB", "hermite_bounded", None)
        pr = _Problem(fm, fp, m, np.full(m.size, 0.02))
        p0, info = initial_walkers(pr, 16, np.random.default_rng(0))
        assert p0.shape == (16, 3)
        assert len({tuple(r) for r in p0}) == 16          # no duplicate walkers
        assert np.all(np.isfinite([pr(t) for t in p0]))
        assert info["n_modes"] >= 1 and sum(info["mode_walkers"]) == 16
