"""
sed_model._inv_core
===================
Internals of the inverse model: a lean photometry-only forward evaluator, the
log-posterior object, prior parsing, and the global-search + Levenberg-
Marquardt walker initialiser.  Nothing here changes the forward numerics:
``_FastModel`` calls exactly the same Fortran kernels, in the same order, as
``run_forward`` (and is checked against it in the test-suite); it just skips
the work that is not needed for magnitudes (full ``ForwardResult``
construction, bolometric integral, nearest-grid-distance diagnostic,
rebuilding the extinction model at every call).
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from ._cc_ext import get_cc_api
from .params import FitParams, PC_TO_CM

_LN2PI = math.log(2.0 * math.pi)
_SQRT2 = math.sqrt(2.0)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _log_ndtr(x: np.ndarray) -> np.ndarray:
    """log of the standard-normal CDF, stable in the lower tail."""
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    out = np.empty_like(x)
    for k, v in enumerate(x):
        if v > -20.0:
            out[k] = math.log(max(0.5 * math.erfc(-v / _SQRT2), 1e-320))
        else:  # asymptotic expansion
            out[k] = -0.5 * v * v - math.log(-v) - 0.5 * _LN2PI + math.log1p(-1.0 / (v * v))
    return out


def grid_valid_nodes(grid) -> np.ndarray:
    """Boolean (nT, nL, nM) mask of grid nodes that hold a real spectrum.

    Kurucz-type grids are not a full cube: unphysical corners (e.g. hot and
    low gravity) are stored as all-zero SEDs, which the forward model turns
    into finite but absurd magnitudes (~770) rather than an error.
    """
    return np.asarray(grid.flux.sum(axis=-1) > 0.0)


def _axis_nodes(axis: np.ndarray, x: float):
    """Nodes of the bracketing cell that actually contribute at x."""
    n = len(axis)
    i0 = int(np.searchsorted(axis, x, side="right")) - 1
    i0 = min(max(i0, 0), n - 2)
    span = axis[i0 + 1] - axis[i0]
    f = (x - axis[i0]) / span if span > 0 else 0.0
    if f < 1e-9:
        return (i0,)
    if f > 1.0 - 1e-9:
        return (i0 + 1,)
    return (i0, i0 + 1)


# ---------------------------------------------------------------------------
# photometry-only forward model
# ---------------------------------------------------------------------------

class _FastModel:
    """Synthetic magnitudes for a fixed list of filters; same kernels as run_forward."""

    def __init__(self, grid, filters, filter_order, R, mag_system, interp_method, extinction):
        self.grid = grid
        self.R = float(R)
        self.method = interp_method.lower()
        if self.method not in ("hermite_bounded", "hermite", "linear"):
            raise ValueError(
                f"Unknown interp_method '{interp_method}'. "
                "Choose 'hermite_bounded', 'hermite' or 'linear'."
            )
        byname = {f.name: f for f in filters}
        self.bands = [
            (np.ascontiguousarray(byname[n].wavelengths, dtype=np.float64),
             np.ascontiguousarray(byname[n].transmission, dtype=np.float64),
             float(byname[n].zero_point(mag_system)))
            for n in filter_order
        ]
        self.n_bands = len(self.bands)
        self.valid = grid_valid_nodes(grid)
        self._has_holes = not bool(self.valid.all())

        self.ext = extinction if (extinction is not None
                                  and getattr(extinction.config, "enabled", False)) else None
        self.ext_k = None
        if self.ext is not None:
            # Av-linear extinction curves (all built-in laws) can be evaluated once.
            try:
                k1 = np.asarray(self.ext.with_av(1.0).extinction_curve(grid.wavelengths))
                k2 = np.asarray(self.ext.with_av(0.7).extinction_curve(grid.wavelengths))
                if np.allclose(k2, 0.7 * k1, rtol=1e-9, atol=1e-12):
                    self.ext_k = k1
            except Exception:
                self.ext_k = None
        self._cc = get_cc_api(required=True)

    # the compiled module is not picklable -> re-acquire after unpickling
    def __getstate__(self):
        d = dict(self.__dict__)
        d["_cc"] = None
        return d

    def __setstate__(self, d):
        self.__dict__.update(d)
        self._cc = get_cc_api(required=True)

    # ------------------------------------------------------------------
    def supported(self, teff, logg, meta) -> bool:
        """False when the query touches an empty (all-zero) grid node."""
        if not self._has_holes:
            return True
        g = self.grid
        it = _axis_nodes(g.teff_grid, teff)
        il = _axis_nodes(g.logg_grid, logg)
        im = _axis_nodes(g.meta_grid, meta)
        return bool(self.valid[np.ix_(it, il, im)].all())

    def mags(self, teff, logg, meta, av, d) -> Optional[np.ndarray]:
        """Magnitudes in ``filter_order``; None if the SED cannot be built."""
        cc = self._cc
        g = self.grid
        tq, lq, mq = g.clamp(teff, logg, meta)
        if self.method == "hermite_bounded":
            flux, ierr = cc.interp_sed_hermite_bounded(
                tq, lq, mq, g.teff_grid, g.logg_grid, g.meta_grid, g.flux_fortran, g.wavelengths)
            if ierr == 2:
                return None
        elif self.method == "hermite":
            flux, ierr = cc.interp_sed_hermite(
                tq, lq, mq, g.teff_grid, g.logg_grid, g.meta_grid, g.flux_fortran)
        else:
            flux, ierr = cc.interp_sed_linear(
                tq, lq, mq, g.teff_grid, g.logg_grid, g.meta_grid, g.flux_fortran)
        flux = cc.dilute_flux(flux, self.R, float(d))
        if self.ext is not None:
            if self.ext_k is not None:
                a_lam = av * self.ext_k if av != 0.0 else None
            else:
                a_lam = np.asarray(self.ext.with_av(av).extinction_curve(g.wavelengths))
            if a_lam is not None:
                flux = flux * 10.0 ** (-0.4 * a_lam)
        out = np.empty(self.n_bands)
        for i, (w, t, zp) in enumerate(self.bands):
            mag, _bf, ie = cc.synthetic_magnitude(g.wavelengths, flux, w, t, zp)
            out[i] = mag if ie == 0 else np.nan
        return out


# ---------------------------------------------------------------------------
# priors
# ---------------------------------------------------------------------------

def parse_priors(priors, free_names):
    """Turn the user ``priors`` dict into {name: callable(value) -> logp}.

    Accepted specs per parameter name (units are the parameter's own, i.e. K,
    dex, mag, and **cm** for ``d``; the 'parallax' and 'edsd' forms for ``d``
    take mas / pc for convenience):

    ``('normal', mu, sigma)``
    ``('halfnormal', sigma)``        (parameter >= 0, e.g. a_v)
    ``('parallax', plx_mas, err_mas)``   d only; Gaussian likelihood on parallax
    ``('edsd', L_pc)``               d only; exponentially decreasing space
                                     density prior  p(d) ~ d^2 exp(-d/L)
    ``callable(value) -> log prior``
    """
    out = {}
    if not priors:
        return out
    for name, spec in priors.items():
        if name not in free_names:
            raise ValueError(
                f"prior given for '{name}' which is not a free parameter "
                f"(free: {free_names})"
            )
        if callable(spec):
            out[name] = spec
            continue
        kind, *a = spec
        kind = str(kind).lower()
        if kind == "normal":
            mu, sg = map(float, a)
            out[name] = (lambda v, mu=mu, sg=sg: -0.5 * ((v - mu) / sg) ** 2)
        elif kind == "halfnormal":
            sg = float(a[0])
            out[name] = (lambda v, sg=sg: -0.5 * (v / sg) ** 2)
        elif kind == "parallax":
            if name != "d":
                raise ValueError("'parallax' prior only applies to d")
            plx, err = map(float, a)
            out[name] = (lambda v, plx=plx, err=err:
                         -0.5 * ((1000.0 * PC_TO_CM / v - plx) / err) ** 2)
        elif kind == "edsd":
            if name != "d":
                raise ValueError("'edsd' prior only applies to d")
            L = float(a[0]) * PC_TO_CM
            out[name] = (lambda v, L=L: 2.0 * math.log(v) - v / L)
        else:
            raise ValueError(f"unknown prior kind '{kind}'")
    return out


# ---------------------------------------------------------------------------
# log-posterior
# ---------------------------------------------------------------------------

class _Problem:
    """Callable log-posterior over theta = [free FitParams..., (sigma_sys)]."""

    def __init__(self, model: _FastModel, fit_params: FitParams, obs_mag, obs_err,
                 upper=None, priors=None, sys_floor=0.0, fit_jitter=False, jitter_max=0.2):
        self.model = model
        self.fp = fit_params
        self.obs = np.asarray(obs_mag, dtype=np.float64)
        self.err = np.asarray(obs_err, dtype=np.float64)
        self.upper = np.zeros(self.obs.size, bool) if upper is None else np.asarray(upper, bool)
        self.det = ~self.upper
        self.n_data = int(self.det.sum())
        self.sys_floor = float(sys_floor)
        self.fit_jitter = bool(fit_jitter)
        self.jitter_max = float(jitter_max)

        self.free_names = list(fit_params.free_names)
        self.names = self.free_names + (["sigma_sys"] if self.fit_jitter else [])
        self.n_dim = len(self.names)
        lo, hi = [], []
        for p in self.free_names:
            sp = getattr(fit_params, p)
            lo.append(sp.lo); hi.append(sp.hi)
        if self.fit_jitter:
            lo.append(0.0); hi.append(self.jitter_max)
        self.lo = np.array(lo, float)
        self.hi = np.array(hi, float)
        # position of each physical parameter in theta, or its fixed value
        self._slot = {}
        k = 0
        for p in ("teff", "logg", "meta", "a_v", "d"):
            sp = getattr(fit_params, p)
            if sp.is_free:
                self._slot[p] = ("free", k); k += 1
            else:
                self._slot[p] = ("fixed", float(sp.value))
        self.priors = [(self.names.index(n), f) for n, f in parse_priors(priors, self.names).items()]
        # log-space coordinates for the optimiser / initial draw
        self.is_log = np.array([n in ("teff", "d") for n in self.names])
        self.n_evals = 0

    # --------------------------------------------------------------
    def _phys(self, theta):
        v = []
        for p in ("teff", "logg", "meta", "a_v", "d"):
            kind, x = self._slot[p]
            v.append(float(theta[x]) if kind == "free" else x)
        return v

    def sigma_eff2(self, theta):
        s2 = self.err ** 2 + self.sys_floor ** 2
        if self.fit_jitter:
            s2 = s2 + float(theta[-1]) ** 2
        return s2

    def model_mags(self, theta):
        teff, logg, meta, av, d = self._phys(theta)
        if not self.model.supported(teff, logg, meta):
            return None
        try:
            return self.model.mags(teff, logg, meta, av, d)
        except Exception:
            return None

    def in_bounds(self, theta):
        return bool(np.all(np.isfinite(theta)) and np.all(theta >= self.lo) and np.all(theta <= self.hi))

    def log_prior(self, theta):
        if not self.in_bounds(theta):
            return -np.inf
        lp = 0.0
        for k, f in self.priors:
            lp += f(float(theta[k]))
        return lp

    def log_like(self, theta, m=None):
        if m is None:
            m = self.model_mags(theta)
        if m is None or not np.all(np.isfinite(m)):
            return -np.inf
        s2 = self.sigma_eff2(theta)
        det = self.det
        ll = -0.5 * float(np.sum((self.obs[det] - m[det]) ** 2 / s2[det] + np.log(2.0 * np.pi * s2[det])))
        if not det.all():
            up = self.upper
            ll += float(np.sum(_log_ndtr((m[up] - self.obs[up]) / np.sqrt(s2[up]))))
        return ll

    def chi2(self, theta):
        m = self.model_mags(theta)
        if m is None:
            return np.inf
        s2 = self.sigma_eff2(theta)
        return float(np.sum((self.obs[self.det] - m[self.det]) ** 2 / s2[self.det]))

    def __call__(self, theta):
        self.n_evals += 1
        theta = np.asarray(theta, dtype=np.float64)
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_like(theta)

    def batch(self, thetas):
        return np.array([self(t) for t in thetas])

    # --------------------------------------------------------------
    # residual vector for Levenberg-Marquardt (data + Gaussian priors)
    def residuals(self, theta):
        m = self.model_mags(theta)
        n_extra = len(self.priors)
        if m is None or not np.all(np.isfinite(m)):
            return np.full(self.n_data + int(self.upper.sum()) + n_extra, 1e3)
        s = np.sqrt(self.sigma_eff2(theta))
        r = [(self.obs[self.det] - m[self.det]) / s[self.det]]
        if not self.det.all():
            up = self.upper
            r.append(np.maximum(0.0, self.obs[up] - m[up]) / s[up])
        for k, f in self.priors:
            r.append(np.array([math.sqrt(max(-2.0 * f(float(theta[k])), 0.0))]))
        return np.concatenate(r)


# ---------------------------------------------------------------------------
# initialisation: global search -> LM refinement -> Laplace ball(s)
# ---------------------------------------------------------------------------

def _to_w(problem, theta):
    return np.where(problem.is_log, np.log(np.maximum(theta, 1e-300)), theta)


def _from_w(problem, w):
    return np.where(problem.is_log, np.exp(w), w)


def _lm(problem, u0, wlo, wspan, n_iter=60, frozen=None):
    """Levenberg-Marquardt in unit-box coordinates u (theta = f(wlo + u*wspan))."""
    n = u0.size

    def res(u):
        return problem.residuals(_from_w(problem, wlo + np.clip(u, 0, 1) * wspan))

    u = np.clip(u0, 0, 1)
    r = res(u)
    c = float(r @ r)
    lam = 1e-2
    h = 2e-4
    J = None
    for _ in range(n_iter):
        J = np.empty((r.size, n))
        for j in range(n):
            up = u.copy()
            step = h if u[j] + h <= 1.0 else -h
            up[j] += step
            J[:, j] = (res(up) - r) / step
        if frozen is not None:
            J[:, frozen] = 0.0
        A = J.T @ J
        g = J.T @ r
        improved = False
        for _try in range(12):
            D = np.diag(np.maximum(np.diag(A), 1e-9))
            try:
                step = np.linalg.solve(A + lam * D, -g)
            except np.linalg.LinAlgError:
                lam *= 10
                continue
            un = np.clip(u + step, 0.0, 1.0)
            rn = res(un)
            cn = float(rn @ rn)
            if cn < c:
                dc = c - cn
                u, r, c = un, rn, cn
                lam = max(lam / 3.0, 1e-9)
                improved = True
                break
            lam *= 4.0
        if not improved or dc < 1e-6 * max(c, 1.0) and np.max(np.abs(step)) < 1e-6:
            break
        if dc < 1e-8:
            break
    # final Jacobian at the optimum for the covariance estimate
    J = np.empty((r.size, n))
    for j in range(n):
        up = u.copy()
        step = h if u[j] + h <= 1.0 else -h
        up[j] += step
        J[:, j] = (res(up) - r) / step
    if frozen is not None:
        J[:, frozen] = 0.0
    return u, c, J


def initial_walkers(problem: _Problem, n_walkers: int, rng: np.random.Generator,
                    n_prior: Optional[int] = None, n_starts: int = 16, n_kicks: int = 3,
                    mult: int = 2500):
    """Global search + local refinement; returns (p0, info)."""
    n = problem.n_dim
    wlo = _to_w(problem, problem.lo)
    whi = _to_w(problem, problem.hi)
    wspan = whi - wlo
    if n_prior is None:
        n_prior = int(max(2000, mult * n))

    # stratified (Latin-hypercube) draw over the prior box (log-uniform in teff, d)
    U = (np.argsort(rng.random((n_prior, n)), axis=0) + rng.random((n_prior, n))) / n_prior
    frozen = None
    if problem.fit_jitter:               # hold sigma_sys fixed during the search
        U[:, -1] = 0.05
        frozen = np.zeros(n, bool); frozen[-1] = True
    lp = np.full(n_prior, -np.inf)
    for i in range(n_prior):
        th = _from_w(problem, wlo + U[i] * wspan)
        lp[i] = problem(th)
    order = np.argsort(-lp)
    order = order[np.isfinite(lp[order])]
    if order.size == 0:
        raise ValueError(
            "The likelihood is not finite anywhere in the prior box "
            "(bad bounds, or the whole box lies in an empty region of the grid)."
        )

    # greedy distinct starting points
    starts = []
    for i in order:
        if all(np.max(np.abs(U[i] - U[j])) > 0.05 for j in starts):
            starts.append(i)
        if len(starts) >= n_starts:
            break

    modes = []
    results = [_lm(problem, U[i], wlo, wspan, frozen=frozen) for i in starts]
    # kick the best optimum a few times to escape shallow local minima of the
    # (piecewise-cubic) chi-square surface
    if results and n_kicks > 0:
        ib = int(np.argmin([r[1] for r in results]))
        for _ in range(n_kicks):
            u0 = np.clip(results[ib][0] + 0.03 * rng.standard_normal(n), 0.0, 1.0)
            if frozen is not None:
                u0[frozen] = results[ib][0][frozen]
            results.append(_lm(problem, u0, wlo, wspan, frozen=frozen))
    for u, c, J in results:
        th = _from_w(problem, wlo + u * wspan)
        lpo = problem(th)
        if not np.isfinite(lpo):
            continue
        if any(np.max(np.abs(u - m["u"])) < 0.02 for m in modes):
            continue
        # covariance in u space from Gauss-Newton:  -2 ln L = sum r^2
        A = J.T @ J
        if frozen is not None:
            A[frozen, frozen] = 1.0 / 0.05 ** 2  # nominal width for the frozen jitter dim
        try:
            cov = np.linalg.inv(A + 1e-10 * np.eye(n) * max(np.trace(A), 1.0))
        except np.linalg.LinAlgError:
            cov = np.diag(np.full(n, 0.05 ** 2))
        ev, V = np.linalg.eigh(0.5 * (cov + cov.T))
        ev = np.clip(ev, (2e-6) ** 2, 0.15 ** 2)
        cov = (V * ev) @ V.T
        modes.append(dict(u=u, lp=lpo, cov=cov, chol=(V * np.sqrt(ev)), logdet=float(np.sum(np.log(ev))), chi2=c))
    if not modes:
        # fall back to best random draws
        for i in order[:3]:
            modes.append(dict(u=U[i], lp=lp[i], cov=np.eye(n) * 0.05 ** 2, chol=np.eye(n) * 0.05,
                              logdet=n * math.log(0.05 ** 2), chi2=float('nan')))

    modes.sort(key=lambda m: -m["lp"])
    best = modes[0]["lp"]
    keep = [m for m in modes if m["lp"] > best - 10.0]

    # Relative posterior mass of each mode by importance sampling from its
    # (box-clipped) Gaussian.  Density is measured in theta (Jacobian of the
    # log-coordinates included) so that modes are comparable.
    lj_idx = np.where(problem.is_log)[0]
    n_is = 120
    for m in keep:
        z = rng.standard_normal((n_is, n))
        U_is = m["u"] + z @ m["chol"].T
        ok = np.all((U_is > 0.0) & (U_is < 1.0), axis=1)
        logq = -0.5 * np.sum(z ** 2, axis=1) - 0.5 * m["logdet"] - 0.5 * n * _LN2PI
        terms = []
        for k in range(n_is):
            if not ok[k]:
                continue
            th = _from_w(problem, wlo + U_is[k] * wspan)
            lpk = problem(th)
            if np.isfinite(lpk):
                lpk += float(np.sum(np.log(th[lj_idx])))
                terms.append(lpk - logq[k])
        if terms:
            t = np.array(terms)
            tm = t.max()
            m["logZ"] = float(tm + np.log(np.sum(np.exp(t - tm)) / n_is))
        else:
            m["logZ"] = -np.inf
    logZ = np.array([m["logZ"] for m in keep])
    if not np.isfinite(logZ).any():
        logZ = np.array([m["lp"] + 0.5 * m["logdet"] for m in keep])
    w = np.exp(logZ - np.max(logZ[np.isfinite(logZ)]))
    w[~np.isfinite(logZ)] = 0.0
    w /= w.sum()

    # Share walkers by mass, but give every plausible mode a viable
    # sub-ensemble; differential-evolution moves then rebalance them.
    min_w = max(n + 3, int(round(0.15 * n_walkers)))
    max_modes = max(1, n_walkers // min_w)
    sel = np.argsort(-logZ)[:max_modes]
    keep = [keep[i] for i in sorted(sel)]
    w = np.array([w[i] for i in sorted(sel)]); w = w / w.sum()
    counts = np.maximum(np.floor(w * n_walkers).astype(int), min_w if len(keep) > 1 else 0)
    while counts.sum() > n_walkers:
        counts[np.argmax(counts)] -= 1
    counts[np.argmax(w)] += n_walkers - counts.sum()

    pos = []
    eps = 1e-6
    for m, cnt in zip(keep, counts):
        for _ in range(int(cnt)):
            th = None
            for scale in (1.0, 0.5, 0.2, 0.05):
                for _try in range(8):
                    u = m["u"] + scale * (m["chol"] @ rng.standard_normal(n))
                    if np.all(u > eps) and np.all(u < 1.0 - eps):
                        t_ = _from_w(problem, wlo + u * wspan)
                        if np.isfinite(problem(t_)):
                            th = t_
                            break
                if th is not None:
                    break
            if th is None:  # tiny, distinct offset from the mode itself
                u = np.clip(m["u"] + 1e-5 * rng.standard_normal(n), eps, 1.0 - eps)
                th = _from_w(problem, wlo + u * wspan)
            pos.append(th)
    p0 = np.array(pos)
    if problem.fit_jitter:  # start the jitter at small, spread values
        p0[:, -1] = np.clip(np.abs(rng.normal(0.0, 0.1 * problem.jitter_max, len(p0))) + 1e-4,
                            1e-4, problem.jitter_max)
    info = dict(
        method="global-search+LM",
        n_prior_draws=int(n_prior),
        n_modes=int(sum(1 for c in counts if c > 0)),
        mode_walkers=[int(c) for c in counts if c > 0],
        mode_weights=[round(float(x), 4) for x in w],
        best_chi2=float(modes[0]["chi2"]),
        best_theta=[float(x) for x in _from_w(problem, wlo + modes[0]["u"] * wspan)],
    )
    return p0, info
