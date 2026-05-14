"""
sed_model.io
=================
Result containers and persistence for SED_Model outputs.

InverseResult
    Holds the full emcee sampler output, summary statistics, and the
    input observations.  Can be saved to / loaded from a NumPy .npz
    archive so posteriors are portable without pickling the sampler.

ForwardResult serialisation
    Thin helpers to write a ForwardResult SED to a two-column CSV.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import numpy as np

from .forward import ForwardResult


# ---------------------------------------------------------------------------
# InverseResult
# ---------------------------------------------------------------------------

@dataclass
class InverseResult:
    """Output of the inverse (photometry → stellar parameters) inference.

    Attributes
    ----------
    samples : ndarray, shape (n_samples, 3)
        Flattened posterior samples after burn-in.
        Columns: Teff (K), logg, [M/H].
    log_prob : ndarray, shape (n_samples,)
        Log-posterior value for each sample.
    filter_names : list of str
        Names of filters used in the fit, in order.
    obs_magnitudes : ndarray, shape (n_filters,)
        Observed magnitudes supplied by the user.
    obs_uncertainties : ndarray, shape (n_filters,)
        Per-filter magnitude uncertainties (1-sigma).
    R : float
        Stellar radius used (cm).
    d : float
        Distance used (cm).
    mag_system : str
        Photometric system ('Vega', 'AB', 'ST').
    n_walkers : int
        Number of emcee walkers.
    n_steps : int
        Total steps per walker (including burn-in).
    n_burn : int
        Steps discarded as burn-in.
    n_thin : int
        Thinning factor applied to the chain.
    acceptance_fraction : ndarray, shape (n_walkers,)
        Per-walker acceptance fraction.
    autocorr_time : ndarray or None, shape (3,)
        Estimated integrated autocorrelation time for each parameter.
        None if estimation failed (chain too short).
    """

    samples:             np.ndarray        # (n_samples, n_free)
    log_prob:            np.ndarray        # (n_samples,)
    filter_names:        list[str]
    obs_magnitudes:      np.ndarray        # (n_filters,)
    obs_uncertainties:   np.ndarray        # (n_filters,)
    R:                   float
    d:                   float
    mag_system:          str
    n_walkers:           int
    n_steps:             int
    n_burn:              int
    n_thin:              int
    acceptance_fraction: np.ndarray        # (n_walkers,)
    autocorr_time:       Optional[np.ndarray] = None
    param_names:         list[str] = field(default_factory=lambda: ["teff", "logg", "meta"])
    fixed_params:        dict[str, float] = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------

    def summary(self) -> dict[str, dict[str, float]]:
        """Return median and 1-sigma credible interval for each sampled parameter."""
        if len(self.param_names) != self.samples.shape[1]:
            raise ValueError(
                f"param_names has {len(self.param_names)} entries but samples has "
                f"{self.samples.shape[1]} columns"
            )

        result = {}
        for i, label in enumerate(self.param_names):
            col = self.samples[:, i]
            med = float(np.percentile(col, 50))
            lo  = float(np.percentile(col, 15.865))
            hi  = float(np.percentile(col, 84.135))
            result[label] = {
                "median":       med,
                "lo":           lo,
                "hi":           hi,
                "lower_1sigma": med - lo,
                "upper_1sigma": hi - med,
            }
        return result

    def map_estimate(self) -> dict[str, float]:
        """Return the maximum a-posteriori sample as a {parameter: value} dict."""
        idx = int(np.argmax(self.log_prob))
        s = self.samples[idx]
        return {name: float(s[i]) for i, name in enumerate(self.param_names)}

    def print_summary(self) -> None:
        """Print a formatted parameter summary for arbitrary free parameters."""
        s = self.summary()
        map_params = self.map_estimate()
        labels = {
            "teff": "Teff   [K]",
            "logg": "logg      ",
            "meta": "[M/H]     ",
            "a_v":  "Av    [mag]",
            "d":    "d      [cm]",
        }

        print(f"\n{'─'*60}")
        print("  SED_Model  —  Posterior Summary")
        print(f"{'─'*60}")
        print(f"  Filters : {', '.join(self.filter_names)}")
        print(f"  System  : {self.mag_system}")
        print(f"  Sampled : {', '.join(self.param_names)}")
        if self.fixed_params:
            fixed = ', '.join(f"{k}={v:.6g}" for k, v in self.fixed_params.items())
            print(f"  Fixed   : {fixed}")
        print(f"  Samples : {len(self.samples)} "
              f"(walkers={self.n_walkers}, steps={self.n_steps}, "
              f"burn={self.n_burn}, thin={self.n_thin})")
        af_mean = float(np.mean(self.acceptance_fraction))
        print(f"  Mean acceptance fraction : {af_mean:.3f}")
        if self.autocorr_time is not None:
            tau = ' / '.join(f"{x:.1f}" for x in np.ravel(self.autocorr_time))
            print(f"  Autocorr times : {tau}")

        print(f"{'─'*60}")
        for param in self.param_names:
            v = s[param]
            label = labels.get(param, param)
            print(f"  {label:<12s} :  {v['median']:>10.3f}"
                  f"  +{v['upper_1sigma']:.3f}  -{v['lower_1sigma']:.3f}")

        map_text = ', '.join(f"{k}={v:.6g}" for k, v in map_params.items())
        print(f"\n  MAP : {map_text}")
        print(f"{'─'*60}\n")

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: str | Path) -> None:
        """Save the posterior to a compressed NumPy archive (.npz).

        Parameters
        ----------
        path : str or Path
            Output file.  The ``.npz`` extension is added if absent.
        """
        path = Path(path)
        if path.suffix != ".npz":
            path = path.with_suffix(".npz")

        meta_dict = {
            "filter_names":    self.filter_names,
            "R":               self.R,
            "d":               self.d,
            "mag_system":      self.mag_system,
            "n_walkers":       self.n_walkers,
            "n_steps":         self.n_steps,
            "n_burn":          self.n_burn,
            "n_thin":          self.n_thin,
            "param_names":     self.param_names,
            "fixed_params":    self.fixed_params,
        }

        arrays = dict(
            samples=self.samples,
            log_prob=self.log_prob,
            obs_magnitudes=self.obs_magnitudes,
            obs_uncertainties=self.obs_uncertainties,
            acceptance_fraction=self.acceptance_fraction,
            meta_json=np.array([json.dumps(meta_dict)]),
        )
        if self.autocorr_time is not None:
            arrays["autocorr_time"] = self.autocorr_time

        np.savez_compressed(path, **arrays)

    @classmethod
    def load(cls, path: str | Path) -> "InverseResult":
        """Load an InverseResult saved with :meth:`save`.

        Parameters
        ----------
        path : str or Path
            Path to a ``.npz`` file produced by :meth:`save`.
        """
        path = Path(path)
        if path.suffix != ".npz":
            path = path.with_suffix(".npz")

        data = np.load(path, allow_pickle=False)
        meta = json.loads(str(data["meta_json"][0]))

        autocorr = (
            data["autocorr_time"]
            if "autocorr_time" in data
            else None
        )

        return cls(
            samples=data["samples"],
            log_prob=data["log_prob"],
            filter_names=meta["filter_names"],
            obs_magnitudes=data["obs_magnitudes"],
            obs_uncertainties=data["obs_uncertainties"],
            R=float(meta["R"]),
            d=float(meta["d"]),
            mag_system=meta["mag_system"],
            n_walkers=int(meta["n_walkers"]),
            n_steps=int(meta["n_steps"]),
            n_burn=int(meta["n_burn"]),
            n_thin=int(meta["n_thin"]),
            acceptance_fraction=data["acceptance_fraction"],
            autocorr_time=autocorr,
            param_names=meta.get("param_names", ["teff", "logg", "meta"]),
            fixed_params=meta.get("fixed_params", {}),
        )


# ---------------------------------------------------------------------------
# ForwardResult  serialisation helpers
# ---------------------------------------------------------------------------

def save_sed(result: ForwardResult, path: str | Path) -> None:
    """Write the observed-frame SED from a ForwardResult to a CSV.

    Columns: wavelength (Å), surface_flux, observed_flux (erg/s/cm²/Å).

    Parameters
    ----------
    result : ForwardResult
    path : str or Path
    """
    path = Path(path)
    header = (
        f"# SED_Model SED output\n"
        f"# Teff={result.teff:.1f} K  logg={result.logg:.3f}  "
        f"[M/H]={result.meta:.3f}\n"
        f"# R={result.R:.4e} cm  d={result.d:.4e} cm\n"
        f"# wavelength_AA, surface_flux_erg_s_cm2_AA, observed_flux_erg_s_cm2_AA"
    )
    data = np.column_stack([
        result.wavelengths,
        result.surface_flux,
        result.observed_flux,
    ])
    np.savetxt(path, data, delimiter=",", header=header, comments="")


def save_magnitudes(result: ForwardResult, path: str | Path) -> None:
    """Write the synthetic magnitudes from a ForwardResult to a CSV.

    Columns: filter_name, magnitude, band_flux (erg/s/cm²/Å).

    Parameters
    ----------
    result : ForwardResult
    path : str or Path
    """
    path = Path(path)
    lines = [
        "# SED_Model synthetic photometry",
        f"# Teff={result.teff:.1f} K  logg={result.logg:.3f}  "
        f"[M/H]={result.meta:.3f}  system={result.mag_system}",
        "# filter_name, magnitude, band_flux_erg_s_cm2_AA",
    ]
    for name, mag in result.magnitudes.items():
        bf = result.band_fluxes.get(name, float("nan"))
        lines.append(f"{name},{mag:.6f},{bf:.6e}")

    path.write_text("\n".join(lines) + "\n")
