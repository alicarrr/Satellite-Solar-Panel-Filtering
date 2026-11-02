"""Run with simple_panel_properties.py inside Simple Optimiser Test.

optimise_filters.py and constants.py are in the parent folder.

Uses the existing optimiser unchanged, prints the checks and results, and saves
four PDF figures and a results pickle in simple_validation/.
Set solve_new_filter=False to replot the saved results without solving again.
"""
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

from pathlib import Path
import sys
import pickle

# Find the shared modules one folder above this test folder.
sys.path.insert(1, str(Path(__file__).resolve().parent.parent))

import numpy as np
import casadi as ca
import optimise_filters as optimiser
import simple_panel_properties as panel
from constants import h, c, k, nm, minute


# Run controls. False loads the saved solutions instead of solving again.
solve_new_filter = True  # True solves and saves; False loads the saved pickle.
make_plots = True        # Save and display figures after solving or loading.
results_file = panel.figure_dir / "validation_results.pkl"

# Solver settings for this validation. Physical inputs are in the panel file.
starts = 5
random_seed = 0
temperature_weight = 1.0
smoothness_weight = 1e-3
slope_per_nm = 1e-2
slope_per_log_wavelength = 5.0


def blackbody_band_power(temperature):
    """Independent Planck calculation using the existing physical constants."""
    wavelength_m = panel.wavelength_nm[:, None] * nm
    bandwidth_m = panel.bandwidth_nm[:, None] * nm
    temperature = np.atleast_1d(temperature)[None, :]
    exponent = h * c / (wavelength_m * k * temperature)
    return 2 * np.pi * h * c**2 * bandwidth_m / wavelength_m**5 / np.expm1(exponent)


def check_solution(filters, temperature, weights, solver_cost=None, slope_limits=False):
    """Recompute the cost, heat balance and temperature history from the solution."""
    alpha, rho, tau = [np.asarray(filters[name]) for name in ("alpha", "rho", "tau")]
    illumination = tau / (1 - rho * panel.panel_reflectivity)
    system_alpha = alpha * (1 + panel.panel_reflectivity * illumination) + panel.panel_absorptivity * illumination
    system_rho = rho + tau * panel.panel_reflectivity * illumination
    power = np.dot(panel.incident_power_w, illumination * panel.electrical_conversion)
    absorbed = np.dot(panel.incident_power_w, system_alpha)
    emitted = panel.area_m2 * (system_alpha @ blackbody_band_power(temperature))
    capacity = panel.mass_kg * panel.heat_capacity_j_per_kg_k
    residual = capacity * np.diff(temperature) - (absorbed - power - emitted[:-1]) * panel.interval_s[:-1]

    # Replay the returned filter independently using the same Euler time step.
    replay = np.empty_like(temperature)
    replay[0] = panel.initial_temperature_k
    for index in range(temperature.size - 1):
        emission = panel.area_m2 * np.dot(system_alpha, blackbody_band_power(replay[index])[:, 0])
        replay[index + 1] = replay[index] + panel.interval_s[index] * (absorbed - power - emission) / capacity

    energy = power * np.sum(panel.interval_s)
    incident_energy = np.sum(panel.incident_power_w) * np.sum(panel.interval_s)
    tracking = np.mean(((temperature - panel.target_temperature_k) / panel.temperature_reference_k)**2)
    roughness = sum(np.sum(np.diff(value)**2) for value in (alpha, rho, tau)) / alpha.size
    cost = -weights[0] * energy / incident_energy + weights[1] * tracking + weights[2] * roughness
    assert np.allclose(alpha + rho + tau, 1, atol=1e-7, rtol=0), "Filter conservation failed"
    assert np.allclose(system_alpha + system_rho, 1, atol=1e-7, rtol=0), "System conservation failed"
    assert all(np.min(value) >= -1e-7 and np.max(value) <= 1 + 1e-7 for value in (alpha, rho, tau))
    assert np.min(temperature) >= panel.temperature_min_k - 1e-7
    assert np.max(temperature) <= panel.temperature_max_k + 1e-7
    assert abs(temperature[0] - panel.initial_temperature_k) < 1e-7
    assert np.max(np.abs(residual)) < 1e-5, "Thermal balance failed"
    assert np.max(np.abs(replay - temperature)) < 1e-5, "Temperature replay failed"
    if solver_cost is not None:
        assert abs(cost - solver_cost) < 1e-7, "Independent cost check failed"
    if slope_limits:
        allowed_step = np.minimum(slope_per_nm * np.diff(panel.wavelength_nm),
                                  slope_per_log_wavelength * np.diff(np.log(panel.wavelength_nm)))
        assert all(np.all(np.abs(np.diff(value)) <= allowed_step + 1e-7) for value in (alpha, rho, tau))
    return {"filters": filters, "temperature": temperature, "power": power,
            "absorbed": absorbed, "emitted": emitted, "energy": energy, "cost": cost,
            "rmse": np.sqrt(np.mean((temperature - panel.target_temperature_k)**2)),
            "max_error": np.max(np.abs(temperature - panel.target_temperature_k)),
            "balance_error": np.max(np.abs(residual)), "replay_error": np.max(np.abs(replay - temperature)),
            "optics_error": np.max(np.abs(alpha + rho + tau - 1)), "system_alpha": system_alpha}


def solve_validation():
    optimiser.t_set = panel.target_temperature_k
    optimiser.t_ref_scale = panel.temperature_reference_k
    optimiser.t_min, optimiser.t_max = panel.temperature_min_k, panel.temperature_max_k
    optimiser.w_elec = 1.0
    optimiser.internal_heat_load_w = 0.0
    count = panel.wavelength_nm.size
    incident = np.broadcast_to(panel.incident_power_w, (1, panel.time_s.size, count)).copy()
    data = optimiser.SimData(
        lam_nm=panel.wavelength_nm, dlam_nm=panel.bandwidth_nm,
        lam_edges_nm=panel.wavelength_edges_nm, t=panel.time_s, dt=panel.interval_s,
        q_in=incident, material=np.array([0]), labels=["panel"],
        mass=np.array([panel.mass_kg]), heat_cap=np.array([panel.heat_capacity_j_per_kg_k]),
        area=np.array([panel.area_m2]), T0=np.array([panel.initial_temperature_k]),
        native_lam_nm=panel.wavelength_nm, native_idx=np.arange(count),
        native_q_spectrum=panel.incident_power_w * np.sum(panel.interval_s))

    # Supply a synthetic lookup in memory so the existing property builder
    # performs the wafer-transmission-to-backing-absorption calculation itself.
    temperature_grid = np.linspace(panel.temperature_min_k, panel.temperature_max_k, 5)
    lookup = {"temperature_C": temperature_grid - 273.15, "wavelength_nm": panel.wavelength_nm,
              "total_absorptivity": np.tile(panel.wafer_absorptivity, (5, 1)),
              "reflectivity": np.tile(panel.wafer_reflectivity, (5, 1)),
              "transmissivity": np.tile(panel.wafer_transmissivity, (5, 1)),
              "electrical_absorptivity_mpp": np.tile(panel.electrical_conversion, (5, 1)),
              "valid_optics_mask": np.ones(count, dtype=bool)}
    original_loader = optimiser.load_spectral_lookup
    optimiser.load_spectral_lookup = lambda: lookup
    optimiser.panel_spectral_source, optimiser.lookup_material = {0: "lookup"}, 0
    try:
        backed = optimiser.build_panel_properties(data, verbose=False)
    finally:
        optimiser.load_spectral_lookup = original_loader

    # Prescribe the equivalent opaque surface independently for comparison.
    temperature = ca.MX.sym("temperature")
    conversion = ca.Function("fixed_conversion", [temperature], [ca.DM(panel.electrical_conversion)])
    opaque = optimiser.PanelProperties(
        alpha_abs={0: panel.panel_absorptivity}, eps_emit={0: panel.panel_absorptivity},
        rho_abs={0: panel.panel_reflectivity}, rho_emit={0: panel.panel_reflectivity},
        alpha_elec={0: conversion}, source={0: "opaque equivalent"},
        t_bounds=(panel.temperature_min_k, panel.temperature_max_k), lookup_info={})
    for field in ("alpha_abs", "eps_emit", "rho_abs", "rho_emit"):
        assert np.allclose(getattr(backed, field)[0], getattr(opaque, field)[0], atol=1e-12, rtol=0)
    for temperature in np.linspace(panel.temperature_min_k, panel.temperature_max_k, 9):
        assert np.allclose(np.asarray(backed.alpha_elec[0](temperature)).ravel(),
                           panel.electrical_conversion, atol=1e-12, rtol=0)
    baseline_t, baseline_power = optimiser.simulate_unfiltered(data, backed)
    opaque_t, opaque_power = optimiser.simulate_unfiltered(data, opaque)
    assert np.allclose(baseline_t, opaque_t, atol=1e-9, rtol=0)
    assert np.allclose(baseline_power, opaque_power, atol=1e-11, rtol=0)

    wafer_heat = np.dot(panel.incident_power_w, panel.wafer_thermal_absorptivity)
    backing_heat = np.dot(panel.incident_power_w, panel.wafer_transmissivity)
    initial_emission = panel.area_m2 * np.dot(panel.panel_absorptivity, blackbody_band_power(panel.initial_temperature_k)[:, 0])
    recovered_heat = (panel.mass_kg * panel.heat_capacity_j_per_kg_k
                      * (baseline_t[0, 1] - baseline_t[0, 0]) / panel.interval_s[0] + initial_emission)
    assert backing_heat > 0
    assert abs(recovered_heat - wafer_heat - backing_heat) < 1e-10

    cases, payloads = {}, {}
    for name, properties, weight_t, weight_s in (
        ("energy_only", backed, 0.0, 0.0),
        ("combined", backed, temperature_weight, smoothness_weight),
        ("opaque_equivalent", opaque, temperature_weight, smoothness_weight),
    ):
        optimiser.w_temp, optimiser.w_smooth = weight_t, weight_s
        use_slopes = name != "energy_only"
        optimiser.d_smooth = {0: {key: slope_per_nm for key in ("alpha", "rho", "tau")}} if use_slopes else None
        optimiser.d_smooth_log = {0: {key: slope_per_log_wavelength for key in ("alpha", "rho", "tau")}} if use_slopes else None
        print(f"\n{name}", flush=True)
        problem = optimiser.build_problem(data, properties, verbose=False,
                                          solver_options=optimiser.ipopt_options_quiet.copy())
        best, trials = optimiser.multi_start(problem, starts=starts, seed=random_seed)
        weights = (1.0, weight_t, weight_s)
        for trial in trials:
            if trial.success:
                check_solution(trial.filters[0], trial.temperature["panel"], weights, trial.J, use_slopes)
        cases[name] = check_solution(best.filters[0], best.temperature["panel"], weights, best.J, use_slopes)
        successful_costs = [trial.J for trial in trials if trial.success]
        cases[name]["successful_starts"] = len(successful_costs)
        cases[name]["attempted_starts"] = len(trials)
        cases[name]["cost_spread"] = np.ptp(successful_costs)
        payloads[name] = optimiser.solution_payload(best, trials, problem, properties)
        payloads[name]["model_assumptions"]["electrical"] = "Temperature-independent electrical step from simple_panel_properties.py"
        payloads[name]["settings"]["n_filter_bands"] = count
        payloads[name]["settings"]["band_spacing"] = "piecewise logarithmic"

    clear_filter = {"alpha": np.zeros(count), "rho": np.zeros(count), "tau": np.ones(count)}
    baseline = check_solution(clear_filter, baseline_t[0], (1.0, temperature_weight, smoothness_weight))
    combined, equivalent = cases["combined"], cases["opaque_equivalent"]
    energy_bound = np.dot(panel.incident_power_w, panel.electrical_conversion) * np.sum(panel.interval_s)
    assert abs(cases["energy_only"]["energy"] - energy_bound) < 1e-3
    assert abs(combined["energy"] - energy_bound) < 1e-3
    assert combined["max_error"] < 1e-3
    assert abs(combined["energy"] - equivalent["energy"]) < 1e-3
    assert abs(combined["cost"] - equivalent["cost"]) < 1e-7
    assert np.max(np.abs(combined["temperature"] - equivalent["temperature"])) < 1e-3

    # Check that the unchanged emission function remains nonlinear with T.
    temperatures = panel.target_temperature_k + np.array([-10, 0, 10])
    emission = []
    for temperature in temperatures:
        bands = optimiser.planck_band_power(float(temperature), ca.DM(panel.wavelength_nm), ca.DM(panel.bandwidth_nm))
        emission.append(panel.area_m2 * np.dot(combined["system_alpha"], np.asarray(bands).ravel()))
    independent = panel.area_m2 * (combined["system_alpha"] @ blackbody_band_power(temperatures))
    assert np.allclose(emission, independent, atol=1e-10, rtol=0)
    assert emission[2] > emission[1] > emission[0] and emission[2] - 2 * emission[1] + emission[0] > 0
    print("All checks passed.")

    # Keep the combined case in the normal optimiser format, with the other
    # validation cases and original plot inputs stored in additional fields.
    payload = payloads["combined"].copy()
    payload["validation_payloads"] = payloads
    payload["validation_cases"] = cases
    payload["validation_unfiltered"] = baseline
    payload["validation_backing"] = (wafer_heat, backing_heat, recovered_heat)
    payload["validation_emission"] = {"temperature": temperatures, "power": np.array(emission)}
    payload["validation_properties"] = {name: getattr(panel, name) for name in (
        "wavelength_edges_nm", "wafer_absorptivity", "wafer_reflectivity",
        "wafer_transmissivity", "electrical_conversion", "wafer_thermal_absorptivity",
        "panel_absorptivity", "panel_reflectivity", "panel_transmissivity")}
    return payload


def print_results(payload):
    """Print the saved results in either solve mode or plot-only mode."""
    cases, baseline = payload["validation_cases"], payload["validation_unfiltered"]
    combined = cases["combined"]
    wafer_heat, backing_heat, recovered_heat = payload["validation_backing"]
    temperatures = payload["validation_emission"]["temperature"]
    emission = payload["validation_emission"]["power"]
    print(f"\nBacking heating: {wafer_heat:.2f} W wafer + {backing_heat:.2f} W backing = {recovered_heat:.2f} W total")
    print("\nFixed-filter emission:")
    for temperature, power in zip(temperatures, emission):
        print(f"  {temperature:.2f} K: {power:.6f} W")
    print(f"\n{'Case':<20} {'Energy (J)':>12} {'Final T (K)':>13} {'RMSE (K)':>13} {'Cost':>13}")
    for name, result in [("unfiltered", baseline)] + list(cases.items()):
        print(f"{name:<20} {result['energy']:12.6f} {result['temperature'][-1]:13.6f} {result['rmse']:13.6g} {result['cost']:13.8f}")
    print(f"Maximum target error: {combined['max_error']:.3e} K")
    print(f"Energy-balance residual: {combined['balance_error']:.3e} J")
    print(f"Temperature replay error: {combined['replay_error']:.3e} K")
    print(f"Optical conservation error: {combined['optics_error']:.3e}")

    for name, case in cases.items():
        print(f"{name}: {case['successful_starts']}/{case['attempted_starts']} starts converged, cost spread {case['cost_spread']:.3e}")


def plot_validation(payload):
    """Regenerate the figures using only the data stored in the pickle."""
    panel.plot_properties(payload["validation_properties"])
    cases, baseline = payload["validation_cases"], payload["validation_unfiltered"]
    combined = cases["combined"]
    fig, ax = plt.subplots(figsize=(7, 4))
    edges_nm = payload["lam_edges_nm"]
    ax.axvspan(400, 900, color="grey", alpha=0.15)
    for key, label in (("alpha", r"$\alpha_f(\lambda)$"), ("rho", r"$\rho_f(\lambda)$"),
                       ("tau", r"$\tau_f(\lambda)$")):
        ax.stairs(combined["filters"][key], edges_nm, baseline=None, label=label)
    ax.set(xscale="log", xlim=(edges_nm[0], edges_nm[-1]), ylim=(-0.02, 1.02),
           xlabel=r"Wavelength, $\lambda$ (nm)", ylabel="Fraction of incident power")
    ax.legend(loc="center left")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(panel.figure_dir / "simple_panel_filter.pdf")

    fig, axes = plt.subplots(3, 1, figsize=(7, 7), sharex=True)
    minutes = payload["t"] / minute
    axes[0].plot(minutes, baseline["temperature"], "C0--", label="Unfiltered")
    axes[0].plot(minutes, cases["energy_only"]["temperature"], "C1:", label="Energy only")
    axes[0].plot(minutes, combined["temperature"], "C2-", label="Combined")
    axes[0].axhline(payload["settings"]["t_set"], color="k", linestyle=":", label=r"Target, $T_{\mathrm{set}}$")
    axes[0].set_ylabel(r"Temperature, $T(t)$ (K)")
    axes[1].plot(minutes, np.full(minutes.size, baseline["power"]), "C0--", label="Unfiltered")
    axes[1].plot(minutes, np.full(minutes.size, combined["power"]), "C2-", label="Combined")
    axes[1].set(ylabel="Electrical power,\n" + r"$\dot{Q}_{\mathrm{elec}}(t)$ (W)", ylim=(0, 3))
    axes[2].plot(minutes, baseline["emitted"], "C0--", label="Unfiltered")
    axes[2].plot(minutes, combined["emitted"], "C2-", label="Combined")
    axes[2].set(xlabel=r"Time, $t$ (min)",
                ylabel="Thermal emission,\n" + r"$\dot{Q}_{\mathrm{emit}}(t)$ (W)")
    for ax in axes:
        ax.set_xlim(0, minutes[-1])
        ax.legend()
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(panel.figure_dir / "simple_panel_response.pdf")
    plt.show()


def main():
    if solve_new_filter:
        payload = solve_validation()
        results_file.parent.mkdir(parents=True, exist_ok=True)
        with results_file.open("wb") as stream:
            pickle.dump(payload, stream, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Saved {results_file}")
    else:
        if not results_file.is_file():
            raise FileNotFoundError("No saved validation results. Set solve_new_filter=True and run once first.")
        with results_file.open("rb") as stream:
            payload = pickle.load(stream)
        print(f"Loaded {results_file}")
    print_results(payload)
    if make_plots:
        plot_validation(payload)


if __name__ == "__main__":
    main()
