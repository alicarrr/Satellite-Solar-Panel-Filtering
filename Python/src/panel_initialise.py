"""Create the material-property file used by panel_optimise.py.

Model assumptions
-----------------
1. Wavelengths are bin centres in nm and exactly match the orbital simulation.
2. Every panel satisfies R + A + T = 1 at every wavelength.
3. R, A, T, emissivity and external quantum efficiency (EQE) lie in [0, 1].
4. EQE cannot exceed absorptivity because a collected photon must be absorbed.
5. Emissivity defaults to absorptivity (spectral Kirchhoff law). It is saved for
   later thermal work but is not used by panel_optimise.py.
6. The step model is deliberately a test model: its EQE is a literal step and
   every collected electron delivers energy q*V at one fixed operating voltage.
7. The Solcore model is a single-junction device at one fixed temperature,
   incidence angle and polarisation. Its MPP is represented by a lookup against
   photogenerated current density.
8. In the Solcore model, spectral shape affects the IV curve only through total
   photogenerated current density. This is the usual superposition assumption.
9. Solcore is evaluated only on the wavelength interval stated in the panel
   definition. Outside it, EQE is zero and the optical extension is either
   explicitly supplied or, for testing only, held at the nearest edge value.
10. The saved IV data are area-normalised. Cell interconnection, mismatch,
    degradation, wiring, converter and MPPT losses are not included.
11. Passive panel types have zero EQE and no IV lookup.

The two solar-cell functions create_step_panel_properties and
create_solcore_panel_properties return exactly the same dictionary structure.
The optimisation script therefore does not know which model created the data.
"""

import pickle
import warnings
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from constants import c, h, nm, q

# ============================= USER INPUTS =============================

# Choose "step" for a transparent test model or "solcore" for the device model.
eqe_model = "solcore"

# Read the exact wavelength array saved by the orbital simulation. This avoids
# independently recreating a grid with different endpoints or round-off.
use_simulation_wavelengths = True
simulation_pickle_file = Path(__file__).resolve().parent / "Simulation Data" / "simulation_data.pkl"

# These values are used only if use_simulation_wavelengths is False. The upper
# limit is excluded, matching np.arange(minimum, maximum, step).
minimum_wavelength_nm = 10.0
maximum_wavelength_nm = 50.0e3
wavelength_step_nm = 10.0

# Both models are evaluated at this fixed panel temperature.
panel_temperature_k = 300.0

# The IV lookup must cover every current density expected in the orbital data.
maximum_current_density_a_m2 = 700.0
minimum_positive_current_fraction = 1.0e-7
number_of_positive_current_points = 181

output_directory = Path(__file__).resolve().parent / "Solar Panel Properties"
output_pickle_file = output_directory / "panel_properties.pkl"
output_plot_file = output_directory / "panel_properties.png"
show_plot = True

# Each solar-panel type contains one definition for each EQE model. Add another
# integer key to define another solar panel type.
solar_panel_definitions = {
    0: {
        "panel_name": "Single-junction GaAs panel",
        "step": {
            "eqe_lower_nm": 400.0,
            "eqe_upper_nm": 900.0,
            "eqe_value": 0.85,
            "operating_voltage_v": 1.0,
            "reflectivity": 0.10,
            "transmissivity": 0.00,
            "emissivity": None,
        },
        "solcore": {
            "cell_material": "GaAs",
            "emitter_thickness_nm": 200.0,
            "emitter_doping_cm3": 3.0e18,
            "emitter_minority_diffusion_length_nm": 500.0,
            "base_thickness_nm": 3000.0,
            "base_doping_cm3": 1.0e17,
            "base_minority_diffusion_length_nm": 5000.0,
            "front_recombination_velocity_cm_s": 1.0e5,
            "back_recombination_velocity_cm_s": 1.0e5,
            "series_resistance_ohm_m2": 0.0,
            "shunt_resistance_ohm_m2": 1.0e14,
            "front_layers": [],
            "rear_layers": [{"material": "Al", "thickness_nm": 200.0, "coherent": True}],
            "minimum_wavelength_nm": 280.0,
            "maximum_wavelength_nm": 2500.0,
            # "edge" is convenient for testing but is not a material prediction.
            # Replace it with measured {"reflectivity": ..., "transmissivity": ...,
            # "emissivity": ...} data for thesis results outside the Solcore range.
            "shortwave_optics": "edge",
            "longwave_optics": "edge",
        },
    }
}

# Passive surfaces are saved beside solar panels so all material properties are
# defined in this file rather than partly inside the optimiser.
passive_panel_definitions = {
    1: {
        "panel_name": "Example aluminium surface",
        "reflectivity": 0.85,
        "transmissivity": 0.00,
        "emissivity": None,
    }
}

# Solcore calculation settings. The negative external voltage range matches the
# n-on-p orientation used by this GaAs definition; saved values are magnitudes.
reference_spectrum = "AM0"
incidence_angle_deg = 0.0
polarisation = "u"
position_step_nm = 1.0
optics_method = "TMM"
electrical_method = "DA"
da_mode = "green"
no_back_reflection = True
external_voltage_min_v = -1.50
external_voltage_max_v = 0.00
external_voltage_points = 3001
internal_voltage_min_v = -2.00
internal_voltage_max_v = 1.00
internal_voltage_points = 2001


# ======================== GOVERNING EQUATIONS =========================

def absorptivity(reflectivity, transmissivity):
    """A(lambda) = 1 - R(lambda) - T(lambda)."""
    return 1.0 - reflectivity - transmissivity


def step_eqe(wavelengths_nm, lower_nm, upper_nm, value):
    """EQE(lambda) = value inside [lower, upper], and zero outside it."""
    return value * ((wavelengths_nm >= lower_nm) & (wavelengths_nm <= upper_nm))


def electrical_power_density(current_density_a_m2, operating_voltage_v):
    """P_el/A = J_op V_op."""
    return current_density_a_m2 * operating_voltage_v


def photogenerated_current_density(eqe, photon_flux_per_m, wavelengths_m):
    """J_ph = q integral EQE(lambda) Phi(lambda) d lambda."""
    return q * np.trapezoid(eqe * photon_flux_per_m, wavelengths_m)


# ======================= INPUT CONVERSIONS ============================

def create_wavelengths():
    """Return the simulation grid, or an explicitly requested fallback grid."""
    if use_simulation_wavelengths:
        if not simulation_pickle_file.is_file():
            raise FileNotFoundError(f"Simulation data not found at {simulation_pickle_file}.")
        with open(simulation_pickle_file, "rb") as input_file:
            simulation_data = pickle.load(input_file)
        if "wavelengths" not in simulation_data:
            raise KeyError("The simulation pickle does not contain 'wavelengths'.")
        wavelengths_nm = np.asarray(simulation_data["wavelengths"], dtype=float)
        source = simulation_pickle_file
    else:
        if wavelength_step_nm <= 0.0 or maximum_wavelength_nm <= minimum_wavelength_nm:
            raise ValueError("The fallback wavelength limits and step are invalid.")
        wavelengths_nm = np.arange(minimum_wavelength_nm, maximum_wavelength_nm, wavelength_step_nm, dtype=float)
        source = "the fallback np.arange settings"

    if wavelengths_nm.ndim != 1 or wavelengths_nm.size < 2 or np.any(~np.isfinite(wavelengths_nm)) or np.any(np.diff(wavelengths_nm) <= 0.0):
        raise ValueError("Wavelengths must contain at least two finite, increasing bin centres.")
    print(f"Using {wavelengths_nm.size} wavelengths from {wavelengths_nm[0]:g} to {wavelengths_nm[-1]:g} nm, read from {source}.")
    return wavelengths_nm.copy()


def evaluate_spectrum(value, wavelengths_nm, name):
    """Convert a scalar, array or function into one value per wavelength bin."""
    if callable(value):
        result = np.asarray(value(wavelengths_nm), dtype=float)
    elif np.isscalar(value):
        result = np.full(wavelengths_nm.shape, float(value))
    else:
        result = np.asarray(value, dtype=float)
    if result.shape != wavelengths_nm.shape or np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must provide one finite value per wavelength bin.")
    return result


def create_optical_properties(wavelengths_nm, definition, name):
    """Evaluate R and T, derive A, and apply the stated emissivity assumption."""
    reflectivity = evaluate_spectrum(definition["reflectivity"], wavelengths_nm, f"{name} reflectivity")
    transmissivity = evaluate_spectrum(definition["transmissivity"], wavelengths_nm, f"{name} transmissivity")
    panel_absorptivity = absorptivity(reflectivity, transmissivity)
    emissivity_value = panel_absorptivity if definition.get("emissivity") is None else definition["emissivity"]
    emissivity = evaluate_spectrum(emissivity_value, wavelengths_nm, f"{name} emissivity")
    return reflectivity, panel_absorptivity, transmissivity, emissivity


def create_current_grid():
    """Create a zero point followed by logarithmically spaced positive currents."""
    if maximum_current_density_a_m2 <= 0.0 or not 0.0 < minimum_positive_current_fraction < 1.0:
        raise ValueError("The current-density lookup limits are invalid.")
    if number_of_positive_current_points < 2:
        raise ValueError("number_of_positive_current_points must be at least two.")
    positive = np.geomspace(minimum_positive_current_fraction * maximum_current_density_a_m2, maximum_current_density_a_m2, number_of_positive_current_points)
    return np.concatenate(([0.0], positive))


def validate_panel(panel):
    """Reject non-physical spectra or an inconsistent common output structure."""
    names = ("reflectivity", "absorptivity", "transmissivity", "emissivity", "external_quantum_efficiency")
    arrays = [np.asarray(panel[name], dtype=float) for name in names]
    if any(array.shape != panel["wavelengths_nm"].shape for array in arrays):
        raise ValueError(f"{panel['panel_name']} contains inconsistent spectral dimensions.")
    if any(np.any(~np.isfinite(array)) or np.any(array < -1.0e-9) or np.any(array > 1.0 + 1.0e-9) for array in arrays):
        raise ValueError(f"{panel['panel_name']} contains a spectral value outside [0, 1].")
    if not np.allclose(panel["reflectivity"] + panel["absorptivity"] + panel["transmissivity"], 1.0, atol=1.0e-8):
        raise ValueError(f"{panel['panel_name']} does not satisfy R + A + T = 1.")
    if np.any(panel["external_quantum_efficiency"] > panel["absorptivity"] + 1.0e-8):
        raise ValueError(f"{panel['panel_name']} has EQE greater than absorptivity.")
    lookup = panel["iv_lookup"]
    if lookup is not None:
        current = np.asarray(lookup["photogenerated_current_density_a_m2"])
        voltage = np.asarray(lookup["operating_voltage_v"])
        operating_current = np.asarray(lookup["operating_current_density_a_m2"])
        power = np.asarray(lookup["electrical_power_density_w_m2"])
        if current.ndim != 1 or current[0] != 0.0 or np.any(np.diff(current) <= 0.0):
            raise ValueError(f"{panel['panel_name']} has an invalid current lookup grid.")
        if any(value.shape != current.shape or np.any(~np.isfinite(value)) for value in (voltage, operating_current, power)):
            raise ValueError(f"{panel['panel_name']} has inconsistent IV lookup arrays.")
        if np.any(voltage < 0.0) or np.any(operating_current < 0.0) or np.any(power < 0.0):
            raise ValueError(f"{panel['panel_name']} has a negative IV lookup value.")
        if np.any(operating_current > current + 1.0e-7 * np.maximum(1.0, current)):
            raise ValueError(f"{panel['panel_name']} delivers more current than it photogenerates.")
        if np.any(np.diff(power) < -1.0e-8):
            raise ValueError(f"{panel['panel_name']} has non-monotonic operating power.")
        if not np.allclose(power, electrical_power_density(operating_current, voltage), rtol=1.0e-7, atol=1.0e-10):
            raise ValueError(f"{panel['panel_name']} does not satisfy P = IV at its operating point.")
    return panel


# ========================= STEP BEHAVIOUR =============================

def create_step_panel_properties(panel_type, definition, wavelengths_nm, temperature_k, current_grid):
    """Create a literal step-EQE, fixed-voltage panel for testing."""
    settings = definition["step"]
    reflectivity, panel_absorptivity, transmissivity, emissivity = create_optical_properties(wavelengths_nm, settings, definition["panel_name"])
    eqe = step_eqe(wavelengths_nm, settings["eqe_lower_nm"], settings["eqe_upper_nm"], settings["eqe_value"])
    voltage = np.full(current_grid.shape, settings["operating_voltage_v"])
    voltage[0] = 0.0
    power = electrical_power_density(current_grid, voltage)
    panel = {
        "panel_type": panel_type,
        "panel_name": definition["panel_name"],
        "is_solar": True,
        "eqe_model": "step",
        "temperature_k": temperature_k,
        "wavelengths_nm": wavelengths_nm.copy(),
        "reflectivity": reflectivity,
        "absorptivity": panel_absorptivity,
        "transmissivity": transmissivity,
        "emissivity": emissivity,
        "external_quantum_efficiency": eqe,
        "iv_lookup": {
            "photogenerated_current_density_a_m2": current_grid.copy(),
            "operating_current_density_a_m2": current_grid.copy(),
            "operating_voltage_v": voltage,
            "electrical_power_density_w_m2": power,
        },
    }
    return validate_panel(panel)


# ======================== SOLCORE BEHAVIOUR ===========================

def build_solcore_cell(settings, temperature_k):
    """Convert one user panel definition into a Solcore single-junction cell."""
    from solcore import material
    from solcore.solar_cell import SolarCell
    from solcore.structure import Junction, Layer

    elements = []
    coherency = []
    for layer in settings["front_layers"]:
        layer_material = material(layer["material"])(T=temperature_k)
        elements.append(Layer(width=layer["thickness_nm"] * nm, material=layer_material, role="front layer"))
        coherency.append("c" if layer["coherent"] else "i")

    emitter = material(settings["cell_material"])(T=temperature_k, Nd=settings["emitter_doping_cm3"] * 1.0e6, hole_diffusion_length=settings["emitter_minority_diffusion_length_nm"] * nm)
    base = material(settings["cell_material"])(T=temperature_k, Na=settings["base_doping_cm3"] * 1.0e6, electron_diffusion_length=settings["base_minority_diffusion_length_nm"] * nm)
    junction = Junction([Layer(width=settings["emitter_thickness_nm"] * nm, material=emitter, role="emitter"), Layer(width=settings["base_thickness_nm"] * nm, material=base, role="base")], sn=settings["front_recombination_velocity_cm_s"] * 1.0e-2, sp=settings["back_recombination_velocity_cm_s"] * 1.0e-2, R_shunt=settings["shunt_resistance_ohm_m2"], T=temperature_k, kind=electrical_method)
    junction_index = len(elements)
    elements.append(junction)
    coherency.extend(["c", "i"])

    for layer in settings["rear_layers"]:
        layer_material = material(layer["material"])(T=temperature_k)
        elements.append(Layer(width=layer["thickness_nm"] * nm, material=layer_material, role="rear layer"))
        coherency.append("c" if layer["coherent"] else "i")

    cell = SolarCell(elements, T=temperature_k, R_series=settings["series_resistance_ohm_m2"])
    return cell, junction_index, coherency


def solve_solcore_optics(cell, junction_index, coherency, wavelengths_nm):
    """Run Solcore once and return R, A, T and EQE on its valid wavelength grid."""
    from solcore.solar_cell_solver import solar_cell_solver

    wavelengths_m = wavelengths_nm * nm
    options = {"wavelength": wavelengths_m, "optics_method": optics_method, "da_mode": da_mode, "theta": incidence_angle_deg, "pol": polarisation, "position": position_step_nm * nm, "no_back_reflection": no_back_reflection, "coherency_list": coherency}
    with redirect_stdout(StringIO()):
        solar_cell_solver(cell, "qe", user_options=options)

    # Remove only solver round-off at the physical bounds, then derive A so the
    # saved optical fractions conserve energy exactly.
    reflectivity = np.clip(np.asarray(cell.reflected, dtype=float), 0.0, 1.0)
    transmissivity = np.clip(np.asarray(cell.transmitted, dtype=float), 0.0, 1.0)
    transmissivity = np.minimum(transmissivity, 1.0 - reflectivity)
    panel_absorptivity = absorptivity(reflectivity, transmissivity)
    eqe = np.clip(np.asarray(cell[junction_index].eqe(wavelengths_m), dtype=float), 0.0, 1.0)
    return reflectivity, panel_absorptivity, transmissivity, eqe


def extension_optics(definition, wavelengths_nm, edge_values, name):
    """Evaluate an explicit out-of-range optical extension."""
    if definition == "edge":
        warnings.warn(f"{name} optics are held at the nearest Solcore edge for testing only.", stacklevel=2)
        reflectivity = np.full(wavelengths_nm.shape, edge_values[0])
        transmissivity = np.full(wavelengths_nm.shape, edge_values[2])
        emissivity = np.full(wavelengths_nm.shape, edge_values[1])
        return reflectivity, absorptivity(reflectivity, transmissivity), transmissivity, emissivity
    if not isinstance(definition, dict) or "reflectivity" not in definition or "transmissivity" not in definition:
        raise ValueError(f"{name} optics must be 'edge' or a dictionary containing reflectivity and transmissivity; a Boolean is not valid.")
    return create_optical_properties(wavelengths_nm, definition, name)


def extend_solcore_properties(wavelengths_nm, core_wavelengths_nm, core_properties, settings):
    """Copy Solcore values into the full grid and fill only the two outer regions."""
    reflectivity = np.empty(wavelengths_nm.shape)
    panel_absorptivity = np.empty(wavelengths_nm.shape)
    transmissivity = np.empty(wavelengths_nm.shape)
    emissivity = np.empty(wavelengths_nm.shape)
    eqe = np.zeros(wavelengths_nm.shape)
    core_mask = (wavelengths_nm >= core_wavelengths_nm[0]) & (wavelengths_nm <= core_wavelengths_nm[-1])
    reflectivity[core_mask], panel_absorptivity[core_mask], transmissivity[core_mask], eqe[core_mask] = core_properties
    emissivity[core_mask] = panel_absorptivity[core_mask]

    regions = [(wavelengths_nm < core_wavelengths_nm[0], settings["shortwave_optics"], (core_properties[0][0], core_properties[1][0], core_properties[2][0]), "shortwave"), (wavelengths_nm > core_wavelengths_nm[-1], settings["longwave_optics"], (core_properties[0][-1], core_properties[1][-1], core_properties[2][-1]), "longwave")]
    for mask, definition, edge_values, name in regions:
        if np.any(mask):
            extension = extension_optics(definition, wavelengths_nm[mask], edge_values, name)
            reflectivity[mask], panel_absorptivity[mask], transmissivity[mask], emissivity[mask] = extension
    return reflectivity, panel_absorptivity, transmissivity, emissivity, eqe


def reference_current_density(wavelengths_nm, eqe):
    """Calculate the AM0 current used only to scale Solcore IV illumination."""
    from solcore.light_source import LightSource

    wavelengths_m = wavelengths_nm * nm
    source = LightSource(source_type="standard", version=reference_spectrum, x=wavelengths_m, output_units="photon_flux_per_m")
    _, photon_flux_per_m = source.spectrum(x=wavelengths_m, output_units="photon_flux_per_m")
    return float(photogenerated_current_density(eqe, photon_flux_per_m, wavelengths_m))


def solve_solcore_operating_point(cell, wavelengths_nm, illumination_scale):
    """Solve one IV curve and return current, voltage and power at its MPP."""
    from solcore.light_source import LightSource
    from solcore.solar_cell_solver import solar_cell_solver

    if illumination_scale == 0.0:
        return 0.0, 0.0, 0.0
    wavelengths_m = wavelengths_nm * nm
    source = LightSource(source_type="standard", version=reference_spectrum, x=wavelengths_m, output_units="photon_flux_per_m", concentration=float(illumination_scale))
    voltages = np.linspace(external_voltage_min_v, external_voltage_max_v, external_voltage_points)
    internal_voltages = np.linspace(internal_voltage_min_v, internal_voltage_max_v, internal_voltage_points)
    options = {"wavelength": wavelengths_m, "light_source": source, "light_iv": True, "mpp": True, "voltages": voltages, "internal_voltages": internal_voltages, "optics_method": None, "da_mode": da_mode}
    with redirect_stdout(StringIO()):
        solar_cell_solver(cell, "iv", user_options=options)
    iv = cell.iv
    voltage_step = abs(external_voltage_max_v - external_voltage_min_v) / (external_voltage_points - 1)
    if abs(iv["Vmpp"] - external_voltage_min_v) <= voltage_step or abs(iv["Vmpp"] - external_voltage_max_v) <= voltage_step:
        raise ValueError("The Solcore MPP reached a voltage-sweep boundary; widen the voltage range.")
    operating_current = abs(float(iv["Impp"]))
    operating_voltage = abs(float(iv["Vmpp"]))
    power = max(float(iv["Pmpp"]), 0.0)
    if not np.isclose(power, electrical_power_density(operating_current, operating_voltage), rtol=1.0e-7, atol=1.0e-10):
        raise ValueError("Solcore returned an MPP that does not satisfy P = IV.")
    return operating_current, operating_voltage, power


def create_solcore_panel_properties(panel_type, definition, wavelengths_nm, temperature_k, current_grid):
    """Create a Solcore EQE/optics model and MPP lookup with the common output."""
    settings = definition["solcore"]
    core_mask = (wavelengths_nm >= settings["minimum_wavelength_nm"]) & (wavelengths_nm <= settings["maximum_wavelength_nm"])
    core_wavelengths_nm = wavelengths_nm[core_mask]
    if core_wavelengths_nm.size < 2:
        raise ValueError("The requested grid has fewer than two wavelengths in the Solcore range.")

    cell, junction_index, coherency = build_solcore_cell(settings, temperature_k)
    core_properties = solve_solcore_optics(cell, junction_index, coherency, core_wavelengths_nm)
    reflectivity, panel_absorptivity, transmissivity, emissivity, eqe = extend_solcore_properties(wavelengths_nm, core_wavelengths_nm, core_properties, settings)
    reference_current = reference_current_density(core_wavelengths_nm, core_properties[3])
    if reference_current <= 0.0:
        raise ValueError("Solcore produced a non-positive AM0 reference current.")

    operating_current = np.zeros(current_grid.shape)
    operating_voltage = np.zeros(current_grid.shape)
    power = np.zeros(current_grid.shape)
    print(f"Solving {current_grid.size - 1} non-zero Solcore IV points for panel type {panel_type}...")
    for index in range(1, current_grid.size):
        operating_current[index], operating_voltage[index], power[index] = solve_solcore_operating_point(cell, core_wavelengths_nm, current_grid[index] / reference_current)

    panel = {
        "panel_type": panel_type,
        "panel_name": definition["panel_name"],
        "is_solar": True,
        "eqe_model": "solcore",
        "temperature_k": temperature_k,
        "wavelengths_nm": wavelengths_nm.copy(),
        "reflectivity": reflectivity,
        "absorptivity": panel_absorptivity,
        "transmissivity": transmissivity,
        "emissivity": emissivity,
        "external_quantum_efficiency": eqe,
        "iv_lookup": {
            "photogenerated_current_density_a_m2": current_grid.copy(),
            "operating_current_density_a_m2": operating_current,
            "operating_voltage_v": operating_voltage,
            "electrical_power_density_w_m2": power,
        },
    }
    return validate_panel(panel)


def create_passive_panel_properties(panel_type, definition, wavelengths_nm, temperature_k):
    """Create the same spectral structure for a panel with no electrical output."""
    reflectivity, panel_absorptivity, transmissivity, emissivity = create_optical_properties(wavelengths_nm, definition, definition["panel_name"])
    panel = {"panel_type": panel_type, "panel_name": definition["panel_name"], "is_solar": False, "eqe_model": None, "temperature_k": temperature_k, "wavelengths_nm": wavelengths_nm.copy(), "reflectivity": reflectivity, "absorptivity": panel_absorptivity, "transmissivity": transmissivity, "emissivity": emissivity, "external_quantum_efficiency": np.zeros(wavelengths_nm.shape), "iv_lookup": None}
    return validate_panel(panel)


# ========================= SAVE AND DISPLAY ===========================

def create_panel_models():
    """Create every configured panel type using the selected solar-cell model."""
    if eqe_model not in ("step", "solcore"):
        raise ValueError('eqe_model must be either "step" or "solcore".')
    wavelengths_nm = create_wavelengths()
    current_grid = create_current_grid()
    panels = {}
    for panel_type, definition in solar_panel_definitions.items():
        if eqe_model == "step":
            panels[panel_type] = create_step_panel_properties(panel_type, definition, wavelengths_nm, panel_temperature_k, current_grid)
        else:
            panels[panel_type] = create_solcore_panel_properties(panel_type, definition, wavelengths_nm, panel_temperature_k, current_grid)
    for panel_type, definition in passive_panel_definitions.items():
        if panel_type in panels:
            raise ValueError(f"Panel type {panel_type} is defined as both solar and passive.")
        panels[panel_type] = create_passive_panel_properties(panel_type, definition, wavelengths_nm, panel_temperature_k)
    return {"schema_version": 3, "wavelengths_nm": wavelengths_nm, "panels": panels, "eqe_model": eqe_model, "temperature_k": panel_temperature_k, "wavelength_source": str(simulation_pickle_file) if use_simulation_wavelengths else "fallback settings"}


def save_panel_models(models):
    output_directory.mkdir(parents=True, exist_ok=True)
    with open(output_pickle_file, "wb") as output_file:
        pickle.dump(models, output_file)
    print(f"Saved panel properties to {output_pickle_file}")


def plot_panel_models(models):
    figure, axes = plt.subplots(2, 1, figsize=(9, 8))
    for panel_type, panel in models["panels"].items():
        wavelengths_nm = panel["wavelengths_nm"]
        axes[0].plot(wavelengths_nm, panel["reflectivity"], label=f"type {panel_type} R")
        axes[0].plot(wavelengths_nm, panel["absorptivity"], label=f"type {panel_type} A")
        axes[0].plot(wavelengths_nm, panel["transmissivity"], label=f"type {panel_type} T")
        if panel["is_solar"]:
            axes[0].plot(wavelengths_nm, panel["external_quantum_efficiency"], label=f"type {panel_type} EQE")
            lookup = panel["iv_lookup"]
            axes[1].plot(lookup["photogenerated_current_density_a_m2"], lookup["electrical_power_density_w_m2"], label=f"type {panel_type} electrical power")
    axes[0].set_xscale("log")
    axes[0].set_ylim(0.0, 1.05)
    axes[0].set_xlabel("Wavelength-bin centre (nm)")
    axes[0].set_ylabel("Fraction")
    axes[1].set_xlabel("Photogenerated current density (A/m²)")
    axes[1].set_ylabel("Operating electrical power density (W/m²)")
    for axis in axes:
        axis.grid(True, alpha=0.3)
        axis.legend()
    figure.tight_layout()
    output_directory.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_plot_file, dpi=200)
    if show_plot:
        plt.show()
    else:
        plt.close(figure)


def main():
    models = create_panel_models()
    save_panel_models(models)
    plot_panel_models(models)


if __name__ == "__main__":
    main()
