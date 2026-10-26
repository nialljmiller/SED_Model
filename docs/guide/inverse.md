# Inverse Model

The inverse model recovers a posterior over stellar parameters from observed broadband photometry, using [`emcee`](https://emcee.readthedocs.io/) ensemble MCMC.

## Statistical model

- **Prior** — flat within the bounds of each free parameter (from [`FitParams`](parameters.md)); \(-\infty\) outside.
- **Likelihood** — Gaussian per filter, summed:

\[
\ln \mathcal{L} = -\frac{1}{2} \sum_i \left[ \frac{(m_i^{\rm obs} - m_i^{\rm model})^2}{\sigma_i^2} + \ln\left(2\pi \sigma_i^2\right) \right]
\]

Model magnitudes come from the **same Fortran kernels, called in the same order, as `run_forward`** (interpolate, dilute, extinguish, convolve), minus the parts that do not affect magnitudes (bolometric integral, nearest-grid diagnostic, rebuilding the extinction model every call). A test checks agreement with `run_forward` to 1e-9 mag. Any exception inside the evaluation, a non-finite model magnitude, or a query that touches an *empty* grid node (Kurucz-type grids store unphysical corners, e.g. hot + low gravity, as all-zero SEDs that the forward model turns into finite magnitudes near 770) returns \(-\infty\) so the walker simply moves away.

Optional extensions, all off by default:

- `sys_floor=0.03` adds a systematic error floor (mag) in quadrature to every uncertainty.
- `fit_jitter=True` samples one extra common term `sigma_sys` (flat on `[0, jitter_max]`), added in quadrature to all uncertainties; it appears as the last column of `samples`.
- `priors={'d': ('parallax', 2.0, 0.05), 'a_v': ('halfnormal', 0.5), ...}` multiplies in extra priors (`normal`, `halfnormal`, `parallax` and `edsd` for distance, or any callable returning a log prior).
- `upper_limits=[False, True, ...]` flags bands whose entry is a magnitude lower limit; they contribute \(\ln\Phi((m^{\rm model}-m^{\rm obs})/\sigma)\).
- Non-finite magnitudes or uncertainties are dropped (with a warning) instead of poisoning the likelihood.

## Basic usage

```python
from sed_model import run_inverse

posterior = run_inverse(
    obs_magnitudes=[5.03, 4.17],
    obs_uncertainties=[0.01, 0.02],
    filter_names=["G", "J"],
    R=6.957e10,           # cm — always required, not sampled
    d=3.086e19,           # cm — convenience: builds default FitParams with d fixed
    grid=grid,
    filters=filters,
    mag_system="AB",
    n_walkers=32, n_steps=2000, n_burn=500, n_thin=1,
    seed=42,
    progress=True,
)
```

With no `fit_params` supplied, a default is built from the grid: `Teff`, `logg`, `[M/H]` free across the full grid, `Av` fixed at 0, distance fixed at `d` (default 1 pc). If `fit_params` **is** supplied, the `d` keyword is ignored — set distance through `fit_params`.

## Choosing what to sample

```python
from sed_model import fit_params_from_grid
from sed_model.params import PC_TO_CM

# Av free with flat prior over [0, 3] mag
params = fit_params_from_grid(grid, a_v=(0.0, 3.0), d_cm=500 * PC_TO_CM)

# Av and distance both free — expect a correlated Av–d posterior
params = fit_params_from_grid(
    grid,
    a_v=(0.0, 2.0),
    d_cm=(100 * PC_TO_CM, 2000 * PC_TO_CM),
)

posterior = run_inverse(..., fit_params=params)
```

When `Av` is active (free, or fixed > 0) and no `extinction` model is supplied, a default Fitzpatrick (1999) model is created automatically; pass `extinction_law=...` or a full `ExtinctionModel` to control the prescription. See [Extinction](extinction.md).

## Walker initialisation

By default (`init='auto'`) the sampler does not start from a ball at the middle of the prior box (which, for the full Kurucz grid, is a 27 000 K, logg 2.5 point that lies in an empty corner of the grid). Instead it

1. evaluates the posterior on a stratified (Latin hypercube) draw over the prior box, log-uniform in Teff and distance;
2. refines the best few distinct candidates with Levenberg-Marquardt on the magnitude residuals (pure NumPy);
3. estimates each optimum's covariance from the Gauss-Newton Jacobian and its relative posterior mass by importance sampling, and
4. starts the walkers in those ellipsoids, giving every plausible mode a viable sub-ensemble.

With `init='auto'` the ensemble uses a differential-evolution move mixture (`moves='de'`, including mode-hopping jumps), which rebalances walkers between modes and decorrelates about twice as fast as emcee's stretch move. Pass `moves='stretch'` for the old behaviour.

To use the legacy Gaussian ball, give a starting point or ask for it:

```python
posterior = run_inverse(
    ...,
    p0_centre={"teff": 5800, "logg": 4.4, "meta": 0.0},   # implies init='ball'
    p0_scatter=0.02,                                      # fraction of each range
)
```

`n_burn=None` picks the burn-in automatically (twice the autocorrelation time, at most half the chain), and `auto_stop=True` samples in blocks and stops once the chain exceeds 50 autocorrelation times.

## Input validation

`run_inverse` raises `ValueError` before any sampling when:

- `obs_magnitudes` and `obs_uncertainties` differ in length, or `filter_names` does not match them;
- any uncertainty is ≤ 0 (the likelihood would be undefined);
- a filter name is not among the supplied filters;
- `n_walkers < 6`, `n_walkers < 2 × n_free`, or an integer `n_burn >= n_steps`;
- every parameter is fixed (nothing to sample).

## Diagnostics

After sampling, a `UserWarning` is emitted if any walker's acceptance fraction falls below 0.1 or above 0.9. The result carries

- `autocorr_time` (integrated autocorrelation time per parameter, `None` if it could not be estimated),
- `converged` (`True` when the post-burn-in chain is longer than 50 autocorrelation times for every parameter) and `n_eff` (effective number of independent samples per parameter),
- `chi2_best` / `n_data` (data chi-square of the maximum-posterior sample and the number of detected bands), and
- `init_info` (how the walkers were started: method, number of modes, walkers per mode, likelihood evaluations).

`posterior.credible_interval('meta', 0.68, method='hpd')` returns the shortest interval containing the requested mass; unlike the equal-tailed intervals in `summary()`, it extends to the bound when the posterior mode sits on a prior bound.

!!! note "Calibration caveats"
    Intervals are only as good as the model: with real photometry the atmosphere library, filter curves and zero points contribute systematic errors that the Gaussian likelihood does not know about (use `sys_floor` or `fit_jitter`). Fitting independent MIST bolometric-correction photometry (different atmospheres, filter curves and Vega zero point) with the Kurucz grid, the nominal 95 % Teff intervals covered the truth 39 % of the time; `fit_jitter=True` raised that to 85 % (0.4.0 benchmark: 46 % plain, 62 % with `sys_floor=0.03`, 89 % with `fit_jitter`). `fit_jitter` makes the intervals honest but does not remove the offset itself: Teff stays about 0.5 % low in that test. Weakly constrained parameters (logg from a few broad bands, Teff-Av degeneracy for few bands) have posteriors dominated by the flat prior and can be genuinely multimodal; check `converged`, and consider informative priors via `priors=`.

!!! note "What to expect"
    Injection-and-recovery on the Kurucz grid (16 stars per cell, 0.005 and 0.02 mag noise, Johnson, optical+NIR and wide band sets, 0.4.1 defaults):

    | Free parameters | Teff bias | Teff scatter (MAD) | Teff 68 % / 95 % coverage |
    |---|---|---|---|
    | Teff, logg, [M/H], Av | +0.1 % | 1.2 % | 0.75 / 0.97 |
    | Teff, logg, [M/H], Av, distance | +0.6 % | 4.9 % | 0.57 / 0.90 |

    With Av and distance both free the default 2000 steps are usually too few, and only about 5 % of those fits report `converged`; the coverage is still reasonable, so the flag is conservative, but increase `n_steps` for science use. Fits to only three bands with Av and distance both free are effectively unconstrained and should not be trusted. logg and [M/H] are poorly constrained by broadband photometry; an informative [M/H] prior (a width of 0.3 dex in the test) cut the [M/H] scatter from 0.60 to 0.21 dex.

!!! warning "Fixed radius"
    `R` is not sampled, so it feeds directly into Teff through L ∝ R²T⁴. A radius 14 % too small biases Teff about 7 % high. Fitting the apparent UBVRIJHK magnitudes of Vega at 7.68 pc, with R = 2.36 R☉ (its polar radius) the result was 10 400 K; with R = 2.71 R☉ (from the interferometric angular diameter) it was 9 500 K, against about 9 600 K from the literature. Use a radius that matches the measurement.

## The result

`run_inverse` returns an [`InverseResult`](results-io.md#inverseresult) holding the flattened post-burn-in chain (`samples`, shape `(n_samples, n_free)` with columns in `param_names` order), per-sample `log_prob`, the observations, sampler settings, acceptance fractions, and the fixed-parameter values used. See [Results and I/O](results-io.md) for summaries, MAP estimates, and `.npz` persistence.
