"""
SED_Model API usage example.

This script demonstrates the intended high-level workflow:

1. Load an atmosphere grid.
2. Load filters.
3. Run the forward model.
4. Save synthetic SED/magnitudes.
5. Run the inverse model on synthetic observations.
6. Print posterior summaries.

Edit the paths below to point at SED_Tools-generated data products.
"""

from pathlib import Path

from sed_model import (
    load_grid,
    load_filters_from_instrument_dir,
    run_forward,
    run_inverse,
    save_sed,
    save_magnitudes,
    fit_params_from_grid,
    make_extinction_model,
    RSUN_TO_CM,
    PC_TO_CM,
)


GRID_DIR = Path("data/stellar_models/Kurucz2003all")
FILTER_DIR = Path("data/filters/Generic/Johnson")

R_SUN = RSUN_TO_CM
D_10_PC = 10.0 * PC_TO_CM


def main() -> None:
    grid = load_grid(GRID_DIR)
    filters = load_filters_from_instrument_dir(FILTER_DIR)

    print(grid)
    print(f"Loaded {len(filters)} filters:")
    for filt in filters:
        print(f"  {filt.name}")

    forward = run_forward(
        teff=5778.0,
        logg=4.44,
        meta=0.0,
        R=R_SUN,
        d=D_10_PC,
        grid=grid,
        filters=filters,
        mag_system="Vega",
        interp_method="hermite",
    )

    print("\nForward synthetic magnitudes:")
    for name, mag in sorted(forward.magnitudes.items()):
        print(f"  {name:8s} {mag:10.5f}")

    save_sed(forward, "forward_sed.csv")
    save_magnitudes(forward, "forward_magnitudes.csv")

    # Build a tiny synthetic inverse problem from the forward result.
    filter_names = ["B", "V"]
    obs_magnitudes = [forward.magnitudes[name] for name in filter_names]
    obs_uncertainties = [0.02, 0.02]

    params = fit_params_from_grid(
        grid,
        a_v=0.0,
        d_cm=D_10_PC,
    )

    posterior = run_inverse(
        obs_magnitudes=obs_magnitudes,
        obs_uncertainties=obs_uncertainties,
        filter_names=filter_names,
        R=R_SUN,
        grid=grid,
        filters=filters,
        fit_params=params,
        mag_system="Vega",
        interp_method="hermite",
        n_walkers=24,
        n_steps=400,
        n_burn=100,
        n_thin=2,
        seed=42,
        progress=False,
    )

    print("\nInverse posterior summary:")
    posterior.print_summary()
    posterior.save("posterior_samples.npz")

    # Example with extinction in the forward model.
    extinction = make_extinction_model(
        enabled=True,
        law="fitzpatrick99",
        a_v=0.3,
        r_v=3.1,
    )

    reddened = run_forward(
        teff=5778.0,
        logg=4.44,
        meta=0.0,
        R=R_SUN,
        d=D_10_PC,
        grid=grid,
        filters=filters,
        mag_system="Vega",
        extinction=extinction,
    )

    print("\nForward magnitudes with Av=0.3:")
    for name, mag in sorted(reddened.magnitudes.items()):
        print(f"  {name:8s} {mag:10.5f}")


if __name__ == "__main__":
    main()
