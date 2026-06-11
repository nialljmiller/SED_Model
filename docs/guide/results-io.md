# Results and I/O

## `ForwardResult` outputs

Two CSV writers persist forward-model products:

```python
from sed_model import save_sed, save_magnitudes

save_sed(result, "forward_sed.csv")            # wavelength, surface_flux, observed_flux
save_magnitudes(result, "forward_magnitudes.csv")  # filter, magnitude, band_flux + metadata
```

Both files carry a `#`-commented header recording the input parameters (`teff`, `logg`, `meta`, `R`, `d`, `mag_system`, ...) so the file is self-describing.

## `InverseResult`

Returned by `run_inverse`; a frozen container with:

| Field | Description |
|---|---|
| `samples` | `(n_samples, n_free)` flattened post-burn-in chain; columns follow `param_names` |
| `log_prob` | `(n_samples,)` log-posterior per sample |
| `param_names` | the free parameters, canonical order (e.g. `['teff', 'logg', 'meta']`) |
| `fixed_params` | `{name: value}` for the parameters that were held fixed |
| `obs_magnitudes`, `obs_uncertainties`, `filter_names` | the data that were fitted |
| `R`, `d`, `mag_system` | physical configuration |
| `n_walkers`, `n_steps`, `n_burn`, `n_thin` | sampler settings |
| `acceptance_fraction` | per-walker acceptance fractions |
| `autocorr_time` | integrated autocorrelation time per free parameter, or `None` |

### Summaries

```python
s = posterior.summary()
# {'teff': {'median': ..., 'lo': ..., 'hi': ...,
#           'lower_1sigma': ..., 'upper_1sigma': ...}, ...}
```

`lo`/`hi` are the 15.865 / 84.135 percentiles (±1σ equivalent); `lower_1sigma`/`upper_1sigma` are the distances from the median to those percentiles.

```python
posterior.print_summary()
```

prints a formatted block including the sampled and fixed parameters, the median ± asymmetric 1σ for each free parameter, mean acceptance fraction, and autocorrelation times.

### MAP estimate

```python
map_params = posterior.map_estimate()
# {'teff': ..., 'logg': ..., 'meta': ...}  — the sample with highest log_prob
```

Returns a **dict keyed by parameter name**, so it works for any free-parameter configuration (3 parameters or 5).

### Persistence

```python
posterior.save("posterior.npz")          # .npz appended if missing
loaded = InverseResult.load("posterior.npz")
```

The `.npz` stores the chain arrays plus a JSON-encoded metadata blob (observations, parameter names, fixed values, sampler settings), so a loaded result round-trips exactly — `summary()`, `map_estimate()`, and `print_summary()` all work on the reloaded object.
