"""Optimise an ideal spectral coating using orbital incident-power data.

The mathematical model is kept in the GOVERNING EQUATIONS section. All file
handling, unit conversion, validation, plotting and CasADi conversion are kept
outside that section.

For incident power P_i in wavelength bin i, the one-pass coating/panel model is

    P_reflected,i = P_i (R_c,i + T_c,i R_p,i)
    P_absorbed,i  = P_i (A_c,i + T_c,i A_p,i)
    P_transmitted,i = P_i T_c,i T_p,i
    A_c,i = 1 - R_c,i - T_c,i

so reflected + absorbed + transmitted equals incident power. Panel reflection
is assumed to leave after one interaction; repeated coating/panel reflections
and interference are not modelled. Electrical power is part of panel-absorbed
power, not a fourth term in the optical energy partition.

Optimisation assumptions
------------------------
1. A coating spectrum is static, uniform, angle-independent and shared by all
   panel types assigned to the same coating group.
2. Reflected and transmitted radiation leaves the modelled satellite system;
   it does not strike another face or an internal component.
3. Irradiance is uniform across each face. Active-area fraction changes only
   electrical output, so a solar panel's optics must describe its whole face.
4. Every solar panel operates continuously at the saved maximum-power point.
5. Extracted electrical power is useful output and is not dissipated back into
   the satellite. All other absorbed optical power becomes heat.
6. Temperature is fixed by the property file. Thermal emission, conduction and
   temperature feedback are outside this optimisation.
7. The electrical-retention constraint applies to total mission energy, not to
   every individual panel or time sample.
8. Time integration is trapezoidal and IPOPT returns a local optimum only.
"""

import pickle
import warnings
from pathlib import Path

import casadi as ca
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from constants import c, h, nm, q

# ============================= USER INPUTS =============================

simulation_pickle_file = Path(__file__).resolve().parent / "Simulation Data" / "simulation_data.pkl"
panel_property_file = Path(__file__).resolve().parent / "Solar Panel Properties" / "panel_properties.pkl"
results_directory = Path(__file__).resolve().parent / "Optimisation Results"
results_pickle_file = results_directory / "panel_coating_results.pkl"

# Meaning of each panelPowers sample in the simulation pickle:
# "power_per_bin"              = W/bin over the complete panel
# "irradiance_per_bin"         = W/m²/bin incident on the panel
# "spectral_irradiance_per_nm" = W/m²/nm at the wavelength-bin centre
simulation_data_kind = "power_per_bin"

# Set either value to None once it is included in the simulation pickle.
manual_panel_types = np.array([0, 0, 0, 0, 1, 1])
manual_panel_areas_m2 = np.full(6, 0.01)

# Panels mapped to one name use the same coating spectrum.
coating_group_by_type = {0: "solar", 1: "aluminium"}

# The current density uses whole-face irradiance. This fraction then scales the
# solar-cell area producing electricity. Optical properties must describe the
# complete face, including inactive area, if this value is below one.
active_area_fraction_by_type = {0: 1.0}

minimum_electrical_fraction = 0.99

# This tiny secondary cost selects the no-coating state R=0, T=1 in wavelength
# bins with negligible physical influence, preventing arbitrary spectral tails.
no_coating_preference_weight = 1.0e-4

ipopt_maximum_iterations = 2000
ipopt_tolerance = 1.0e-9
ipopt_print_level = 5
show_plots = True
save_plots = True
save_results_to_file = True


# ======================== GOVERNING EQUATIONS =========================

def absorptivity(reflectivity, transmissivity):
    """A(lambda) = 1 - R(lambda) - T(lambda)."""
    return 1.0 - reflectivity - transmissivity


def split_incident_power(incident_power, coating_reflectivity, coating_transmissivity, panel):
    """Apply the one-pass spectral optical-power equations stated above."""
    coating_absorption = absorptivity(coating_reflectivity, coating_transmissivity)
    power_reaching_panel = multiply(incident_power, coating_transmissivity)
    coating_reflected = multiply(incident_power, coating_reflectivity)
    coating_absorbed = multiply(incident_power, coating_absorption)
    panel_reflected = multiply(power_reaching_panel, panel["reflectivity"])
    panel_absorbed = multiply(power_reaching_panel, panel["absorptivity"])
    transmitted = multiply(power_reaching_panel, panel["transmissivity"])
    return {"reflected": coating_reflected + panel_reflected, "coating_absorbed": coating_absorbed, "panel_absorbed": panel_absorbed, "absorbed": coating_absorbed + panel_absorbed, "transmitted": transmitted}


def photogenerated_current_density(incident_power, coating_transmissivity, eqe, wavelengths_m, panel_area_m2):
    """J_ph = q/A sum_i[P_i T_c,i EQE_i lambda_i/(h c)]."""
    collected_photons_per_second = multiply(multiply(multiply(incident_power, coating_transmissivity), eqe), wavelengths_m) / (h * c)
    return q * sum_spectrum(collected_photons_per_second) / panel_area_m2


def electrical_power(mpp_power_density_w_m2, panel_area_m2, active_area_fraction):
    """P_el = A_panel f_active (P_mpp/A_active)."""
    return panel_area_m2 * active_area_fraction * mpp_power_density_w_m2


def heat_power(absorbed_power, electrical_power_output):
    """P_heat = P_absorbed - P_el."""
    return absorbed_power - electrical_power_output


def integrate_energy(power_samples, times_s):
    """E = sum_k 0.5 (P_k + P_k+1) (t_k+1 - t_k)."""
    energy = 0
    for time_index in range(len(times_s) - 1):
        dt = times_s[time_index + 1] - times_s[time_index]
        energy += 0.5 * (power_samples[time_index] + power_samples[time_index + 1]) * dt
    return energy


def objective_function(heat_energy_j, incident_energy_j, coating_deviation):
    """J = E_heat/E_incident + w mean[R_c² + (1 - T_c)²]."""
    return heat_energy_j / incident_energy_j + no_coating_preference_weight * coating_deviation


def minimum_electrical_energy(baseline_electrical_energy_j):
    """E_el >= eta_min E_el,no-coating."""
    return minimum_electrical_fraction * baseline_electrical_energy_j


def add_coating_bounds(opti, reflectivity, transmissivity):
    """Apply 0 <= R_c, 0 <= T_c and R_c + T_c <= 1 at every wavelength."""
    opti.subject_to(opti.bounded(0.0, reflectivity, 1.0))
    opti.subject_to(opti.bounded(0.0, transmissivity, 1.0))
    opti.subject_to(reflectivity + transmissivity <= 1.0)


# ===================== NUMPY/CASADI CONVERSIONS =======================

def is_casadi(value):
    return isinstance(value, (ca.DM, ca.MX, ca.SX))


def multiply(left, right):
    """Use elementwise multiplication for either numeric or symbolic spectra."""
    if is_casadi(left) or is_casadi(right):
        return ca.times(left, right)
    return np.asarray(left) * np.asarray(right)


def sum_spectrum(values):
    """Implement the sum over wavelength with the appropriate array backend."""
    if is_casadi(values):
        return ca.sum1(values)
    return float(np.sum(values))


def sum_squares(values):
    if is_casadi(values):
        return ca.sumsqr(values)
    return float(np.sum(np.asarray(values) ** 2))


def as_casadi_spectrum(values):
    return ca.DM(np.asarray(values, dtype=float).reshape(-1, 1))


# ===================== DATA LOADING AND CHECKS ========================

def calculate_bin_widths_nm(wavelengths_nm):
    """Convert centre samples to bin widths only when input is per nm."""
    edges_nm = np.empty(wavelengths_nm.size + 1)
    edges_nm[1:-1] = 0.5 * (wavelengths_nm[:-1] + wavelengths_nm[1:])
    edges_nm[0] = wavelengths_nm[0] - 0.5 * (wavelengths_nm[1] - wavelengths_nm[0])
    edges_nm[-1] = wavelengths_nm[-1] + 0.5 * (wavelengths_nm[-1] - wavelengths_nm[-2])
    return np.diff(edges_nm)


def one_value_per_panel(value, panel_count, name, dtype):
    array = np.asarray(value, dtype=dtype)
    if array.ndim == 0:
        array = np.full(panel_count, array.item(), dtype=dtype)
    if array.shape != (panel_count,):
        raise ValueError(f"{name} must contain one value per panel.")
    return array


def convert_simulation_values_to_power_per_bin(values, wavelengths_nm, panel_areas_m2):
    """Perform the declared unit conversion before any governing equation."""
    if simulation_data_kind == "power_per_bin":
        return values
    if simulation_data_kind == "irradiance_per_bin":
        return values * panel_areas_m2[None, :, None]
    if simulation_data_kind == "spectral_irradiance_per_nm":
        widths_nm = calculate_bin_widths_nm(wavelengths_nm)
        return values * panel_areas_m2[None, :, None] * widths_nm[None, None, :]
    raise ValueError("simulation_data_kind is not one of the three documented choices.")


def load_simulation_data():
    saved = pd.read_pickle(simulation_pickle_file)
    times_s = np.asarray(saved["simulation_times"], dtype=float)
    wavelengths_nm = np.asarray(saved["wavelengths"], dtype=float)
    panel_labels = list(saved["panelLabels"])
    frames = [saved["panelPowers"][label] for label in panel_labels]

    for label, frame in zip(panel_labels, frames):
        frame_wavelengths_nm = np.asarray(frame.columns, dtype=float)
        frame_times_s = np.asarray(frame.index, dtype=float)
        if frame_wavelengths_nm.shape != wavelengths_nm.shape or not np.allclose(frame_wavelengths_nm, wavelengths_nm, rtol=0.0, atol=1.0e-9):
            raise ValueError(f"Panel {label} columns do not match the saved wavelengths.")
        if frame_times_s.shape != times_s.shape or not np.allclose(frame_times_s, times_s, rtol=1.0e-12, atol=1.0e-12):
            raise ValueError(f"Panel {label} rows do not match the saved simulation times.")

    saved_values = np.stack([frame.to_numpy(dtype=float) for frame in frames], axis=1)
    time_count, panel_count, wavelength_count = saved_values.shape
    type_source = manual_panel_types if manual_panel_types is not None else saved.get("panelTypes")
    area_source = manual_panel_areas_m2 if manual_panel_areas_m2 is not None else saved.get("panelAreas_m2")
    if type_source is None or area_source is None:
        raise ValueError("Panel types and areas must be supplied manually or stored in the simulation pickle.")
    panel_types = one_value_per_panel(type_source, panel_count, "panel types", int)
    panel_areas_m2 = one_value_per_panel(area_source, panel_count, "panel areas", float)

    if times_s.shape != (time_count,) or wavelengths_nm.shape != (wavelength_count,):
        raise ValueError("Simulation dimensions are inconsistent.")
    if times_s.size < 2 or np.any(np.diff(times_s) <= 0.0):
        raise ValueError("Simulation times must contain at least two increasing values.")
    if wavelengths_nm.size < 2 or np.any(np.diff(wavelengths_nm) <= 0.0):
        raise ValueError("Wavelengths must contain at least two increasing values.")
    if np.any(~np.isfinite(saved_values)) or np.any(panel_areas_m2 <= 0.0):
        raise ValueError("Incident data must be finite and panel areas must be positive.")
    tolerance = 1.0e-12 * max(1.0, float(np.max(np.abs(saved_values))))
    if np.min(saved_values) < -tolerance:
        raise ValueError("The simulation contains materially negative incident power.")
    saved_values = np.maximum(saved_values, 0.0)
    incident_power = convert_simulation_values_to_power_per_bin(saved_values, wavelengths_nm, panel_areas_m2)
    return {"times_s": times_s, "wavelengths_nm": wavelengths_nm, "panel_labels": panel_labels, "panel_types": panel_types, "panel_areas_m2": panel_areas_m2, "incident_power_w_per_bin": incident_power, "simulation_data_kind": simulation_data_kind}


def align_panel_wavelengths(panel, source_wavelengths_nm, target_wavelengths_nm):
    """Resample saved spectra onto the simulation grid outside the equations."""
    panel = panel.copy()
    if source_wavelengths_nm.shape == target_wavelengths_nm.shape and np.allclose(source_wavelengths_nm, target_wavelengths_nm, rtol=0.0, atol=1.0e-9):
        panel["wavelengths_nm"] = target_wavelengths_nm.copy()
        return panel

    source_step_nm = np.max(np.diff(source_wavelengths_nm))
    lower_extrapolation_nm = max(source_wavelengths_nm[0] - target_wavelengths_nm[0], 0.0)
    upper_extrapolation_nm = max(target_wavelengths_nm[-1] - source_wavelengths_nm[-1], 0.0)
    if lower_extrapolation_nm > source_step_nm + 1.0e-9 or upper_extrapolation_nm > source_step_nm + 1.0e-9:
        raise ValueError(f"Panel wavelengths contain {source_wavelengths_nm.size} bins from {source_wavelengths_nm[0]:g} to {source_wavelengths_nm[-1]:g} nm, but the simulation contains {target_wavelengths_nm.size} bins from {target_wavelengths_nm[0]:g} to {target_wavelengths_nm[-1]:g} nm. Run panel_initialise.py with use_simulation_wavelengths = True and the same simulation_pickle_file used by panel_optimise.py.")
    if lower_extrapolation_nm > 0.0 or upper_extrapolation_nm > 0.0:
        warnings.warn("The simulation extends by at most one property-grid step. The nearest panel value is used at that endpoint; rerun panel_initialise.py to remove this compatibility fallback.", stacklevel=2)

    panel["reflectivity"] = np.interp(target_wavelengths_nm, source_wavelengths_nm, panel["reflectivity"])
    panel["transmissivity"] = np.interp(target_wavelengths_nm, source_wavelengths_nm, panel["transmissivity"])
    panel["absorptivity"] = absorptivity(panel["reflectivity"], panel["transmissivity"])
    panel["emissivity"] = np.interp(target_wavelengths_nm, source_wavelengths_nm, panel["emissivity"])
    panel["external_quantum_efficiency"] = np.interp(target_wavelengths_nm, source_wavelengths_nm, panel["external_quantum_efficiency"])
    panel["external_quantum_efficiency"] = np.minimum(panel["external_quantum_efficiency"], panel["absorptivity"])
    panel["wavelengths_nm"] = target_wavelengths_nm.copy()
    return panel


def load_panel_models(simulation):
    with open(panel_property_file, "rb") as input_file:
        saved = pickle.load(input_file)
    if saved.get("schema_version") != 3:
        raise ValueError("The property file is obsolete; rerun the revised panel_initialise.py.")
    source_wavelengths_nm = np.asarray(saved["wavelengths_nm"], dtype=float)
    target_wavelengths_nm = simulation["wavelengths_nm"]
    if source_wavelengths_nm.ndim != 1 or source_wavelengths_nm.size < 2 or np.any(np.diff(source_wavelengths_nm) <= 0.0):
        raise ValueError("The panel property file contains an invalid wavelength grid.")
    print(f"Loaded panel properties from {panel_property_file}.")
    print(f"Property wavelength source: {saved.get('wavelength_source', 'not recorded')}.")
    print(f"Simulation grid: {target_wavelengths_nm.size} bins from {target_wavelengths_nm[0]:g} to {target_wavelengths_nm[-1]:g} nm.")
    print(f"Property grid:   {source_wavelengths_nm.size} bins from {source_wavelengths_nm[0]:g} to {source_wavelengths_nm[-1]:g} nm.")

    models = {}
    for panel_type in np.unique(simulation["panel_types"]):
        panel_type = int(panel_type)
        if panel_type not in saved["panels"]:
            raise ValueError(f"Panel type {panel_type} has no material definition.")
        if panel_type not in coating_group_by_type:
            raise ValueError(f"Panel type {panel_type} has no coating group.")
        model = align_panel_wavelengths(saved["panels"][panel_type], source_wavelengths_nm, target_wavelengths_nm)
        if model["is_solar"]:
            if panel_type not in active_area_fraction_by_type or not 0.0 < active_area_fraction_by_type[panel_type] <= 1.0:
                raise ValueError(f"Panel type {panel_type} needs an active-area fraction in (0, 1].")
        models[panel_type] = model
    return models


def mpp_power_density_numeric(panel, current_density_a_m2):
    lookup = panel["iv_lookup"]
    current_grid = lookup["photogenerated_current_density_a_m2"]
    if current_density_a_m2 > current_grid[-1] * (1.0 + 1.0e-10):
        raise ValueError(f"{panel['panel_name']} exceeds its IV lookup; increase maximum_current_density_a_m2.")
    return float(np.interp(current_density_a_m2, current_grid, lookup["electrical_power_density_w_m2"]))


def create_mpp_interpolants(panel_models):
    interpolants = {}
    for panel_type, panel in panel_models.items():
        if panel["iv_lookup"] is None:
            continue
        lookup = panel["iv_lookup"]
        interpolants[panel_type] = ca.interpolant(f"mpp_type_{panel_type}", "linear", [lookup["photogenerated_current_density_a_m2"]], lookup["electrical_power_density_w_m2"])
    return interpolants


# ===================== PERFORMANCE CALCULATION ========================

def no_coating(group_names, wavelength_count):
    result = {}
    for group_name in group_names:
        result[group_name] = {"reflectivity": np.zeros(wavelength_count), "transmissivity": np.ones(wavelength_count), "absorptivity": np.zeros(wavelength_count), "emissivity": np.zeros(wavelength_count)}
    return result


def calculate_performance(simulation, panel_models, coatings):
    """Evaluate the same governing equations numerically for fixed coatings."""
    times_s = simulation["times_s"]
    wavelengths_m = simulation["wavelengths_nm"] * nm
    incident_power = simulation["incident_power_w_per_bin"]
    time_count, panel_count, _ = incident_power.shape
    names = ("incident", "reflected", "coating_absorbed", "panel_absorbed", "absorbed", "electrical", "heat", "transmitted")
    power = {name: np.zeros(time_count) for name in names}

    for time_index in range(time_count):
        for panel_index in range(panel_count):
            panel_type = int(simulation["panel_types"][panel_index])
            panel = panel_models[panel_type]
            coating = coatings[coating_group_by_type[panel_type]]
            bins = incident_power[time_index, panel_index]
            split = split_incident_power(bins, coating["reflectivity"], coating["transmissivity"], panel)
            incident_total = sum_spectrum(bins)
            reflected_total = sum_spectrum(split["reflected"])
            coating_absorbed_total = sum_spectrum(split["coating_absorbed"])
            panel_absorbed_total = sum_spectrum(split["panel_absorbed"])
            transmitted_total = sum_spectrum(split["transmitted"])
            absorbed_total = coating_absorbed_total + panel_absorbed_total
            electrical_total = 0.0

            if panel["is_solar"]:
                area_m2 = simulation["panel_areas_m2"][panel_index]
                current_density = photogenerated_current_density(bins, coating["transmissivity"], panel["external_quantum_efficiency"], wavelengths_m, area_m2)
                mpp_density = mpp_power_density_numeric(panel, current_density)
                electrical_total = electrical_power(mpp_density, area_m2, active_area_fraction_by_type[panel_type])
                if electrical_total > panel_absorbed_total + 1.0e-8 * max(1.0, panel_absorbed_total):
                    raise ValueError("Electrical power exceeds power absorbed by the solar-cell panel.")

            if not np.isclose(incident_total, reflected_total + absorbed_total + transmitted_total, rtol=1.0e-9, atol=1.0e-10):
                raise ValueError("The optical power split does not conserve energy.")
            power["incident"][time_index] += incident_total
            power["reflected"][time_index] += reflected_total
            power["coating_absorbed"][time_index] += coating_absorbed_total
            power["panel_absorbed"][time_index] += panel_absorbed_total
            power["absorbed"][time_index] += absorbed_total
            power["electrical"][time_index] += electrical_total
            power["transmitted"][time_index] += transmitted_total
        power["heat"][time_index] = heat_power(power["absorbed"][time_index], power["electrical"][time_index])

    energy = {name: float(integrate_energy(samples, times_s)) for name, samples in power.items()}
    if energy["heat"] < -1.0e-8 * max(1.0, energy["absorbed"]):
        raise ValueError("Total heat generation is negative.")
    return {"power_w": power, "energy_j": energy}


# ========================= OPTIMISER SETUP ============================

def setup_optimiser(simulation, panel_models, baseline):
    """Create decision variables, equations, bounds, cost and IPOPT settings."""
    if not 0.0 <= minimum_electrical_fraction <= 1.0 or no_coating_preference_weight < 0.0:
        raise ValueError("The electrical fraction or coating preference weight is invalid.")
    opti = ca.Opti()
    wavelengths_m = as_casadi_spectrum(simulation["wavelengths_nm"] * nm)
    incident_power = simulation["incident_power_w_per_bin"]
    wavelength_count = simulation["wavelengths_nm"].size
    panel_count = incident_power.shape[1]
    group_names = list(dict.fromkeys(coating_group_by_type[int(panel_type)] for panel_type in simulation["panel_types"]))
    mpp_interpolants = create_mpp_interpolants(panel_models)
    coatings = {}
    deviation = 0

    for group_name in group_names:
        reflectivity = opti.variable(wavelength_count)
        transmissivity = opti.variable(wavelength_count)
        add_coating_bounds(opti, reflectivity, transmissivity)
        opti.set_initial(reflectivity, 0.0)
        opti.set_initial(transmissivity, 1.0)
        coatings[group_name] = {"reflectivity": reflectivity, "transmissivity": transmissivity}
        deviation += sum_squares(reflectivity) + sum_squares(1.0 - transmissivity)

    deviation /= 2 * wavelength_count * len(group_names)

    electrical_power_by_time = []
    absorbed_power_by_time = []
    for time_index in range(simulation["times_s"].size):
        total_electrical = 0
        total_absorbed = 0
        for panel_index in range(panel_count):
            panel_type = int(simulation["panel_types"][panel_index])
            numeric_panel = panel_models[panel_type]
            panel = {"reflectivity": as_casadi_spectrum(numeric_panel["reflectivity"]), "absorptivity": as_casadi_spectrum(numeric_panel["absorptivity"]), "transmissivity": as_casadi_spectrum(numeric_panel["transmissivity"])}
            coating = coatings[coating_group_by_type[panel_type]]
            bins = as_casadi_spectrum(incident_power[time_index, panel_index])
            split = split_incident_power(bins, coating["reflectivity"], coating["transmissivity"], panel)
            panel_absorbed = sum_spectrum(split["panel_absorbed"])
            absorbed = sum_spectrum(split["absorbed"])
            panel_electrical = 0
            if numeric_panel["is_solar"]:
                area_m2 = simulation["panel_areas_m2"][panel_index]
                eqe = as_casadi_spectrum(numeric_panel["external_quantum_efficiency"])
                current_density = photogenerated_current_density(bins, coating["transmissivity"], eqe, wavelengths_m, area_m2)
                mpp_density = mpp_interpolants[panel_type](current_density)
                panel_electrical = electrical_power(mpp_density, area_m2, active_area_fraction_by_type[panel_type])
                opti.subject_to(panel_electrical <= panel_absorbed + 1.0e-10)
            total_electrical += panel_electrical
            total_absorbed += absorbed
        electrical_power_by_time.append(total_electrical)
        absorbed_power_by_time.append(total_absorbed)

    electrical_energy = integrate_energy(electrical_power_by_time, simulation["times_s"])
    absorbed_energy = integrate_energy(absorbed_power_by_time, simulation["times_s"])
    heat_energy = heat_power(absorbed_energy, electrical_energy)
    required_electrical_energy = minimum_electrical_energy(baseline["energy_j"]["electrical"])
    if required_electrical_energy > 0.0:
        opti.subject_to(electrical_energy >= required_electrical_energy)
    objective = objective_function(heat_energy, baseline["energy_j"]["incident"], deviation)
    opti.minimize(objective)
    options = {"max_iter": ipopt_maximum_iterations, "tol": ipopt_tolerance, "print_level": ipopt_print_level, "hessian_approximation": "limited-memory"}
    opti.solver("ipopt", {"expand": False}, options)
    return opti, coatings, objective


def solve_optimisation(simulation, panel_models):
    group_names = list(dict.fromkeys(coating_group_by_type[int(panel_type)] for panel_type in simulation["panel_types"]))
    baseline_coatings = no_coating(group_names, simulation["wavelengths_nm"].size)
    baseline = calculate_performance(simulation, panel_models, baseline_coatings)
    if baseline["energy_j"]["incident"] <= 0.0:
        raise ValueError("The total incident energy must be positive.")

    opti, coating_variables, objective = setup_optimiser(simulation, panel_models, baseline)
    solution = opti.solve()
    coatings = {}
    for group_name, variables in coating_variables.items():
        reflectivity = np.asarray(solution.value(variables["reflectivity"]), dtype=float).reshape(-1)
        transmissivity = np.asarray(solution.value(variables["transmissivity"]), dtype=float).reshape(-1)
        absorption = absorptivity(reflectivity, transmissivity)
        coatings[group_name] = {"reflectivity": reflectivity, "transmissivity": transmissivity, "absorptivity": absorption, "emissivity": absorption.copy()}

    optimised = calculate_performance(simulation, panel_models, coatings)
    required = minimum_electrical_energy(baseline["energy_j"]["electrical"])
    if optimised["energy_j"]["electrical"] < required - 1.0e-6 * max(1.0, required):
        raise RuntimeError("The returned solution violates the electrical-energy constraint.")
    return {"schema_version": 3, "wavelengths_nm": simulation["wavelengths_nm"].copy(), "panel_labels": simulation["panel_labels"], "panel_types": simulation["panel_types"].copy(), "panel_areas_m2": simulation["panel_areas_m2"].copy(), "simulation_data_kind": simulation["simulation_data_kind"], "coatings": coatings, "baseline": baseline, "optimised": optimised, "objective": float(solution.value(objective)), "minimum_electrical_fraction": minimum_electrical_fraction, "solver_status": solution.stats()["return_status"]}


# ========================= SAVE AND DISPLAY ===========================

def incident_energy_spectrum(simulation):
    spectral_power = np.zeros((simulation["times_s"].size, simulation["wavelengths_nm"].size))
    for time_index in range(spectral_power.shape[0]):
        for panel_index in range(simulation["incident_power_w_per_bin"].shape[1]):
            spectral_power[time_index] += simulation["incident_power_w_per_bin"][time_index, panel_index]
    return np.asarray(integrate_energy(spectral_power, simulation["times_s"]), dtype=float)


def plot_results(simulation, results):
    wavelengths_nm = simulation["wavelengths_nm"]
    figure, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(wavelengths_nm, incident_energy_spectrum(simulation))
    axes[0].set_ylabel("Incident energy (J/bin)")
    for group_name, coating in results["coatings"].items():
        axes[1].plot(wavelengths_nm, coating["reflectivity"], label=f"{group_name} R")
        axes[1].plot(wavelengths_nm, coating["absorptivity"], label=f"{group_name} A")
        axes[1].plot(wavelengths_nm, coating["transmissivity"], label=f"{group_name} T")
    axes[1].set_xlabel("Wavelength-bin centre (nm)")
    axes[1].set_ylabel("Coating fraction")
    axes[1].set_ylim(0.0, 1.05)
    axes[1].legend(ncol=3)
    for axis in axes:
        axis.set_xscale("log")
        axis.grid(True, alpha=0.3)
    figure.tight_layout()

    energy_figure, energy_axis = plt.subplots(figsize=(9, 4))
    categories = ["reflected", "absorbed", "electrical", "heat", "transmitted"]
    x = np.arange(len(categories))
    energy_axis.bar(x - 0.2, [results["baseline"]["energy_j"][name] for name in categories], width=0.4, label="No coating")
    energy_axis.bar(x + 0.2, [results["optimised"]["energy_j"][name] for name in categories], width=0.4, label="Optimised")
    energy_axis.set_xticks(x, categories)
    energy_axis.set_ylabel("Energy (J)")
    energy_axis.grid(True, axis="y", alpha=0.3)
    energy_axis.legend()
    energy_figure.tight_layout()

    if save_plots:
        results_directory.mkdir(parents=True, exist_ok=True)
        figure.savefig(results_directory / "coating_spectra.png", dpi=200)
        energy_figure.savefig(results_directory / "energy_comparison.png", dpi=200)
    if show_plots:
        plt.show()
    else:
        plt.close("all")


def save_results(results):
    results_directory.mkdir(parents=True, exist_ok=True)
    with open(results_pickle_file, "wb") as output_file:
        pickle.dump(results, output_file)
    print(f"Saved optimisation results to {results_pickle_file}")


def main():
    simulation = load_simulation_data()
    panel_models = load_panel_models(simulation)
    results = solve_optimisation(simulation, panel_models)
    print(f"Solver status: {results['solver_status']}")
    print(f"No-coating electrical energy: {results['baseline']['energy_j']['electrical']:.6g} J")
    print(f"Optimised electrical energy: {results['optimised']['energy_j']['electrical']:.6g} J")
    print(f"No-coating heat generation: {results['baseline']['energy_j']['heat']:.6g} J")
    print(f"Optimised heat generation: {results['optimised']['energy_j']['heat']:.6g} J")
    if save_results_to_file:
        save_results(results)
    plot_results(simulation, results)


if __name__ == "__main__":
    main()
