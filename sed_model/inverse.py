"""
sed_model.inverse
======================
Inverse model: observed magnitudes → posterior on stellar parameters.

The inverse model is the mirror of the forward model.  Both share the same
``FitParams`` object (from ``sed_model.params``) which declares each
physical parameter as fixed or free.  The MCMC samples only the free
parameters; fixed ones are threaded to the forward model unchanged.

This bidirectionality is concrete, not cosmetic:
  - The forward model maps ``FitParams.unpack(theta)`` → magnitudes.
  - The inverse model maps magnitudes → ``FitParams`` + posterior on theta.
  - The same ``FitParams`` instance flows in both directions.

Parameter modes
---------------
Every physical quantity — Teff, logg, [M/H], Av, distance — is described
by a ``ParamSpec`` inside a ``FitParams``.  Two behaviours:

``fixed(value)``
    Not sampled.  Passed unchanged to the forward model at every step.

``free(lo, hi)``
    Sampled by MCMC with a flat prior over [lo, hi].

``bounded(lo, hi)``
    Alias for ``free`` — use it to be explicit about hard limits.

The simplest use::

    from sed_model import run_inverse
    from sed_model.params import fit_params_from_grid, PC_TO_CM

    params = fit_params_from_grid(grid)                # Teff/logg/meta free
    result = run_inverse(obs_mags, obs_errs, filter_names,
                         R=R_sun, fit_params=params,
                         grid=grid, filters=filters)

To fit with extinction free::

    params = fit_params_from_grid(grid, a_v=(0.0, 3.0))   # Av free 0–3 mag
    ext    = ExtinctionModel(enabled=True, law='fitzpatrick99', a_v=0.0)
    result = run_inverse(..., fit_params=params, extinction=ext)

To fit with distance free (expect Av–d degeneracy in the posterior)::

    params = fit_params_from_grid(
        grid,
        a_v=(0.0, 2.0),
        d_cm=(100*PC_TO_CM, 2000*PC_TO_CM),
    )
"""

from __future__ import annotations

import warnings
from typing import Optional, TYPE_CHECKING

import numpy as np

from .grid import AtmosphereGrid
from .filters import Filter
from .forward import run_forward
from .params import FitParams, fit_params_from_grid, PC_TO_CM
from .io import InverseResult
from ._inv_core import _FastModel, _Problem, initial_walkers

if TYPE_CHECKING:
    from .sed_extinction import ExtinctionModel


# chi2_best / n_data above this triggers the fit_jitter=True hint
_POOR_FIT_CHI2_PER_DATUM = 5.0


# ---------------------------------------------------------------------------
# Prior
# ---------------------------------------------------------------------------

def _log_prior(theta: np.ndarray, fit_params: FitParams) -> float:
    """Flat prior within bounds for all free parameters; -inf outside."""
    return 0.0 if fit_params.in_prior(theta) else -np.inf


# ---------------------------------------------------------------------------
# Likelihood
# ---------------------------------------------------------------------------

def _log_likelihood(
    theta: np.ndarray,
    obs_mag: np.ndarray,
    obs_err: np.ndarray,
    filter_order: list,
    R: float,
    fit_params: FitParams,
    grid: AtmosphereGrid,
    filters: list,
    mag_system: str,
    interp_method: str,
    extinction: Optional["ExtinctionModel"],
) -> float:
    """Gaussian log-likelihood summed over filters.

    Calls ``run_forward`` in FitParams mode so both directions share
    the same parameter unpacking logic.
    """
    try:
        result = run_forward(
            fit_params=fit_params,
            theta=theta,
            R=R,
            grid=grid,
            filters=filters,
            mag_system=mag_system,
            interp_method=interp_method,
            extinction=extinction,
        )
    except Exception:
        return -np.inf

    ll = 0.0
    for i, fname in enumerate(filter_order):
        m_model = result.magnitudes.get(fname)
        if m_model is None or not np.isfinite(m_model):
            return -np.inf
        sigma2 = obs_err[i] ** 2
        ll += -0.5 * ((obs_mag[i] - m_model) ** 2 / sigma2
                      + np.log(2.0 * np.pi * sigma2))
    return ll


# ---------------------------------------------------------------------------
# Posterior
# ---------------------------------------------------------------------------

def _log_posterior(
    theta: np.ndarray,
    obs_mag, obs_err, filter_order,
    R, fit_params, grid, filters,
    mag_system, interp_method, extinction,
) -> float:
    lp = _log_prior(theta, fit_params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + _log_likelihood(
        theta, obs_mag, obs_err, filter_order,
        R, fit_params, grid, filters,
        mag_system, interp_method, extinction,
    )


def _named_moves(emcee, name: str):
    name = name.lower()
    if name == "stretch":
        return None
    if name == "de":
        return [(emcee.moves.DEMove(), 0.7),
                (emcee.moves.DEMove(gamma0=1.0), 0.15),   # mode-hopping jumps
                (emcee.moves.StretchMove(), 0.15)]
    raise ValueError("moves string must be 'stretch' or 'de'")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_inverse(
    obs_magnitudes:    list,
    obs_uncertainties: list,
    filter_names:      list,
    R:                 float,
    grid:              AtmosphereGrid,
    filters:           list,
    # FitParams — the shared parameter spec
    fit_params:        Optional[FitParams] = None,
    # convenience: if fit_params not supplied, these are used to build one
    d:                 Optional[float] = None,
    extinction_law:    Optional[str]   = None,
    # extinction model (required if a_v is free or if you want reddening)
    extinction:        Optional["ExtinctionModel"] = None,
    # sampler settings
    mag_system:        str   = "AB",
    interp_method:     str   = "hermite_bounded",
    n_walkers:         int   = 32,
    n_steps:           int   = 2000,
    n_burn:            Optional[int] = 500,
    n_thin:            int   = 1,
    p0_centre:         Optional[dict] = None,
    p0_scatter:        float = 0.02,
    seed:              Optional[int]  = None,
    progress:          bool  = True,
    verbose:           bool  = True,
    # --- options added after 0.3.1 (all default to a sensible/legacy-safe value) ---
    init:              str   = "auto",
    sys_floor:         float = 0.0,
    fit_jitter:        bool  = False,
    jitter_max:        float = 0.2,
    priors:            Optional[dict] = None,
    upper_limits:      Optional[list] = None,
    auto_stop:         bool  = False,
    moves=None,
    pool=None,
) -> InverseResult:
    """Infer stellar parameters from observed broadband photometry.

    Parameters
    ----------
    obs_magnitudes : array-like, shape (n_filters,)
        Observed magnitudes in the same order as ``filter_names``.  Non-finite
        entries (NaN/inf) are treated as missing and dropped, with a warning.
    obs_uncertainties : array-like, shape (n_filters,)
        1-sigma uncertainties.  Must all be > 0.
    filter_names : list of str
        Filter names matching elements of ``filters``.
    R : float
        Stellar radius in cm (always required; not currently sampled).
    grid : AtmosphereGrid
    filters : list of Filter
    fit_params : FitParams or None
        Full parameter specification.  If None, one is built from ``grid``
        using the ``d`` convenience argument (default: 1 pc, fixed).
    d : float or None
        Distance in cm.  Used only when ``fit_params`` is None to build a
        default ``FitParams`` with fixed distance.  If ``fit_params`` is
        supplied this argument is ignored — set distance via fit_params.
    extinction_law : str or None
        Convenience shorthand used to build a default extinction model when
        both ``fit_params`` and ``extinction`` are None and Av is active.
        Ignored (with a warning) when ``fit_params`` is supplied.
    extinction : ExtinctionModel or None
        Dust extinction model.  Required if Av is free in fit_params.
        Ignored if Av is fixed at 0.
    mag_system : {'AB', 'Vega', 'ST'}
    interp_method : {'hermite_bounded', 'hermite', 'linear'}
    n_walkers, n_steps, n_burn, n_thin : int
        emcee settings.  ``n_steps`` is the (maximum) chain length;
        ``n_burn=None`` chooses the burn-in automatically (2 integrated
        autocorrelation times, at most half the chain).
    p0_centre : dict or None
        Optional starting point keyed by parameter name, e.g.
        ``{'teff': 5800, 'logg': 4.4, 'meta': 0.0}``.  When given, walkers
        start in a small Gaussian ball around it (the legacy behaviour,
        equivalent to ``init='ball'``).
    p0_scatter : float
        Width of initial ball as fraction of each parameter's range
        (only used by ``init='ball'`` / ``p0_centre``).
    seed : int or None
    progress : bool
        Show the emcee progress bar.
    verbose : bool
        Print the parameter summary and extinction configuration before
        sampling.  Set to False for batch runs.
    init : {'auto', 'ball'}
        Walker initialisation.  ``'auto'`` (default) does a stratified global
        search over the prior box (log-uniform in Teff and distance), refines
        the best distinct candidates with Levenberg-Marquardt, and starts the
        walkers in Laplace-approximation ellipsoids around the optima
        (several modes are populated in proportion to their approximate mass).
        With ``p0_centre`` supplied, or ``init='ball'``, the legacy Gaussian
        ball (centre = p0_centre or the middle of the prior box) is used.
    sys_floor : float
        Systematic error floor in mag, added in quadrature to every
        uncertainty (default 0 = off).
    fit_jitter : bool
        Fit one extra common noise term ``sigma_sys`` (mag, flat prior on
        [0, ``jitter_max``]) added in quadrature to all uncertainties.  It
        appears as the last column of ``samples`` / entry of ``param_names``.
    priors : dict or None
        Extra priors on free parameters, e.g.
        ``{'d': ('parallax', 2.0, 0.1)}`` (mas) or ``{'a_v': ('halfnormal', 0.5)}``.
        See ``sed_model._inv_core.parse_priors`` for the supported forms; a
        callable ``f(value) -> log prior`` is also accepted.
    upper_limits : array-like of bool or None
        Mark bands whose ``obs_magnitudes`` entry is a magnitude *lower limit*
        (flux upper limit).  Their likelihood is Phi((m_model - m_obs)/sigma).
    auto_stop : bool
        Sample in blocks and stop early once the chain is >50 autocorrelation
        times long and the estimate has stabilised; warns if ``n_steps`` is
        reached first.
    moves : emcee move, list, str or None
        Passed to ``emcee.EnsembleSampler``.  ``None`` (default) means the
        differential-evolution mixture ``'de'`` (mode-hopping jumps; handles
        multimodal/degenerate posteriors) when ``init='auto'``, and emcee's
        stretch move for the legacy ball.  ``'stretch'`` forces the latter.
    pool : multiprocessing-style pool or None
        Evaluate the posterior in parallel.

    Returns
    -------
    InverseResult
        Includes ``converged``, ``n_eff``, ``chi2_best``, ``n_data`` and
        ``init_info`` diagnostics in addition to the chain.
    """
    try:
        import emcee
    except ImportError as exc:
        raise ImportError(
            "emcee is required for the inverse model. "
            "Install with: pip install emcee"
        ) from exc

    # ------------------------------------------------------------------
    # Input validation
    # ------------------------------------------------------------------
    obs_mag = np.asarray(obs_magnitudes,    dtype=np.float64)
    obs_err = np.asarray(obs_uncertainties, dtype=np.float64)

    if obs_mag.shape != obs_err.shape:
        raise ValueError(
            "obs_magnitudes and obs_uncertainties must have the same length."
        )
    if len(filter_names) != len(obs_mag):
        raise ValueError(
            "filter_names must have the same length as obs_magnitudes."
        )
    if not np.all(obs_err > 0):
        raise ValueError(
            "All obs_uncertainties must be positive (> 0). "
            "Zero or negative uncertainties make the likelihood undefined."
        )
    if n_walkers < 6:
        raise ValueError("n_walkers must be >= 6 (2 × n_dim).")
    if n_burn is not None and n_burn >= n_steps:
        raise ValueError("n_burn must be < n_steps.")
    if init not in ("auto", "ball"):
        raise ValueError("init must be 'auto' or 'ball'.")
    if sys_floor < 0:
        raise ValueError("sys_floor must be >= 0.")

    available = {f.name for f in filters}
    missing   = set(filter_names) - available
    if missing:
        raise ValueError(
            f"Filter names not in the supplied filters: {missing}. "
            f"Available: {available}"
        )

    upper = (np.zeros(obs_mag.size, bool) if upper_limits is None
             else np.asarray(upper_limits, dtype=bool))
    if upper.shape != obs_mag.shape:
        raise ValueError("upper_limits must have the same length as obs_magnitudes.")

    finite = np.isfinite(obs_mag) & np.isfinite(obs_err)
    if not finite.all():
        warnings.warn(
            f"Dropping {int((~finite).sum())} band(s) with non-finite magnitude or "
            f"uncertainty: {[n for n, ok in zip(filter_names, finite) if not ok]}",
            UserWarning, stacklevel=2,
        )
        obs_mag, obs_err, upper = obs_mag[finite], obs_err[finite], upper[finite]
        filter_names = [n for n, ok in zip(filter_names, finite) if ok]
    if obs_mag.size == 0 or upper.all():
        raise ValueError("No detected (finite, non-upper-limit) bands to fit.")

    # ------------------------------------------------------------------
    # Build FitParams if not supplied
    # ------------------------------------------------------------------
    if fit_params is None:
        d_cm = d if d is not None else PC_TO_CM
        fit_params = fit_params_from_grid(grid, d_cm=d_cm)
    elif extinction_law is not None:
        warnings.warn(
            "extinction_law is ignored when fit_params is provided. "
            "Configure the extinction law via the extinction parameter instead.",
            UserWarning,
            stacklevel=2,
        )

    n_dim_phys = fit_params.n_free
    if n_dim_phys == 0:
        raise ValueError(
            "All parameters are fixed — nothing to sample. "
            "Make at least one parameter free in fit_params."
        )
    n_dim = n_dim_phys + (1 if fit_jitter else 0)
    if n_walkers < 2 * n_dim:
        raise ValueError(
            f"n_walkers ({n_walkers}) must be >= 2 × n_free ({2*n_dim}). "
            f"Free parameters: {fit_params.free_names}"
        )

    # ------------------------------------------------------------------
    # Build a default extinction model if Av is active but none provided.
    # ------------------------------------------------------------------
    av_active = fit_params.a_v.is_free or (fit_params.a_v.is_fixed and fit_params.a_v.value > 0.0)
    if extinction is None and av_active:
        from .sed_extinction import make_extinction_model
        extinction = make_extinction_model(
            enabled=True,
            law=extinction_law or 'fitzpatrick99',
            a_v=0.0 if fit_params.a_v.is_free else fit_params.a_v.value,
        )

    # ------------------------------------------------------------------
    # Optionally log what is being sampled vs fixed
    # ------------------------------------------------------------------
    if verbose:
        print(fit_params.summary())
        if extinction is not None and getattr(extinction.config, 'enabled', False):
            cfg = extinction.config
            av_mode = "free (sampled)" if fit_params.a_v.is_free else f"fixed={fit_params.a_v.value:.3f}"
            print(
                f"[inverse] Extinction: law={cfg.law}, Rv={cfg.r_v:.2f}, "
                f"Av {av_mode}"
            )
        else:
            print("[inverse] Extinction: disabled")

    # ------------------------------------------------------------------
    # Posterior object (lean photometry evaluator; same kernels as run_forward)
    # ------------------------------------------------------------------
    model = _FastModel(grid, filters, list(filter_names), R, mag_system,
                       interp_method, extinction)
    problem = _Problem(model, fit_params, obs_mag, obs_err, upper=upper,
                       priors=priors, sys_floor=sys_floor,
                       fit_jitter=fit_jitter, jitter_max=jitter_max)

    # ------------------------------------------------------------------
    # Initial walker positions
    # ------------------------------------------------------------------
    rng = np.random.default_rng(seed)
    use_ball = (init == "ball") or (p0_centre is not None)
    if use_ball:
        p0 = fit_params.initial_ball(
            n_walkers, centre=p0_centre or None, scatter=p0_scatter, rng=rng
        )
        if fit_jitter:
            extra = np.clip(np.abs(rng.normal(0.0, 0.1 * jitter_max, n_walkers)) + 1e-4,
                            1e-4, jitter_max)
            p0 = np.column_stack([p0, extra])
        init_info = dict(method="ball")
    else:
        p0, init_info = initial_walkers(problem, n_walkers, rng)

    # ------------------------------------------------------------------
    # Run MCMC
    # ------------------------------------------------------------------
    if moves is None and not use_ball:
        # Differential-evolution moves mix between modes and decorrelate ~2x
        # faster than the plain stretch move (see CHANGELOG); the legacy ball
        # path keeps emcee's default stretch move.
        moves = "de"
    if isinstance(moves, str):
        moves = _named_moves(emcee, moves)
    if pool is None:
        sampler = emcee.EnsembleSampler(n_walkers, n_dim, problem.batch,
                                        moves=moves, vectorize=True)
    else:
        sampler = emcee.EnsembleSampler(n_walkers, n_dim, problem,
                                        moves=moves, pool=pool)

    done = 0
    state = p0
    converged_auto = False
    if auto_stop:
        block = max(100, n_steps // 20)
        last_tau = None
        while done < n_steps:
            k = min(block, n_steps - done)
            state = sampler.run_mcmc(state, k, progress=progress, skip_initial_state_check=done > 0)
            done += k
            try:
                tau_now = sampler.get_autocorr_time(tol=0, quiet=True)
            except Exception:
                continue
            tmax = float(np.nanmax(tau_now))
            if np.isfinite(tmax) and done > 50.0 * tmax:
                if last_tau is not None and abs(tmax - last_tau) / tmax < 0.05:
                    converged_auto = True
                    break
            if np.isfinite(tmax):
                last_tau = tmax
        if not converged_auto:
            warnings.warn(
                f"auto_stop: chain did not reach 50 autocorrelation times within "
                f"n_steps={n_steps}; increase n_steps or check for multimodality.",
                UserWarning, stacklevel=2,
            )
    else:
        sampler.run_mcmc(state, n_steps, progress=progress)
        done = n_steps

    acc = sampler.acceptance_fraction
    if np.any(acc < 0.1) or np.any(acc > 0.9):
        warnings.warn(
            f"Some walkers have extreme acceptance fractions "
            f"(min={acc.min():.2f}, max={acc.max():.2f}). "
            f"Consider adjusting n_walkers or p0_scatter.",
            UserWarning, stacklevel=2,
        )

    try:
        tau = sampler.get_autocorr_time(tol=0, quiet=True)
        if not np.all(np.isfinite(tau)):
            tau = None
    except Exception:
        tau = None

    if n_burn is None:
        n_burn_eff = 0 if tau is None else int(min(done // 2, np.ceil(2.0 * np.nanmax(tau))))
    else:
        n_burn_eff = int(min(n_burn, done - 1))

    flat_samples  = sampler.get_chain(discard=n_burn_eff, thin=n_thin, flat=True)
    flat_log_prob = sampler.get_log_prob(discard=n_burn_eff, thin=n_thin, flat=True)

    converged = None
    n_eff = None
    if tau is not None:
        n_post = done - n_burn_eff
        converged = bool(n_post > 50.0 * float(np.nanmax(tau)))
        n_eff = n_walkers * n_post / np.maximum(tau, 1.0)

    ib = int(np.argmax(flat_log_prob))
    chi2_best = float(problem.chi2(flat_samples[ib]))

    if converged is False and not (auto_stop and not converged_auto):
        warnings.warn(
            "The chain is not converged (shorter than 50 autocorrelation times "
            "after burn-in), so the posterior may be unreliable. Increase "
            "n_steps, or set auto_stop=True to sample until it is.",
            UserWarning, stacklevel=2,
        )
    if (not fit_jitter and problem.n_data > 0 and np.isfinite(chi2_best)
            and chi2_best / problem.n_data > _POOR_FIT_CHI2_PER_DATUM):
        warnings.warn(
            f"Poor fit: chi2_best/n_data = {chi2_best / problem.n_data:.1f} "
            f"(> {_POOR_FIT_CHI2_PER_DATUM:g}) with fit_jitter=False. The quoted "
            "uncertainties probably under-estimate the model/photometry "
            "mismatch; consider fit_jitter=True (or sys_floor).",
            UserWarning, stacklevel=2,
        )
    init_info = dict(init_info)
    init_info["n_evals"] = int(problem.n_evals)

    fixed_params = {name: getattr(fit_params, name).value for name in fit_params.fixed_names}

    return InverseResult(
        samples=flat_samples,
        log_prob=flat_log_prob,
        filter_names=list(filter_names),
        obs_magnitudes=obs_mag,
        obs_uncertainties=obs_err,
        R=float(R),
        d=float(fit_params.d.value if fit_params.d.is_fixed else np.nan),
        mag_system=mag_system,
        n_walkers=n_walkers,
        n_steps=done,
        n_burn=n_burn_eff,
        n_thin=n_thin,
        acceptance_fraction=acc,
        autocorr_time=tau,
        param_names=list(problem.names),
        fixed_params=fixed_params,
        converged=converged,
        n_eff=n_eff,
        chi2_best=chi2_best,
        n_data=problem.n_data,
        init_info=init_info,
    )
