"""Build the spectral lookup consumed by optimise_filters.py.

Edit the settings below, then run this file (no command-line arguments).
Requires: numpy, scipy, pandas, PyYAML, matplotlib, solcore; tested with
Solcore 5.10.1. Keep the project's constants.py beside these scripts.

Model: a bare, homogeneous silicon n-on-p wafer in vacuum, with phase-averaged
optics. This represents the cell material, not coverglass, contacts, coatings
or a complete spacecraft panel. Geometry, doping, diffusion lengths and relative
permittivity follow Solcore Education's Example 1a: Simple Si solar cell.
Surface recombination and shunt resistance use the documented Solcore DA defaults.
n,k alone determine optics; they do NOT specify the junction's efficiency.

R/A/T use the supplied n,k at its stated temperature. Electrical properties
vary with temperature through the band gap, intrinsic carrier density and dark
current; the example's diffusion lengths are fixed. Optical temperature shifts
and radiation damage are omitted. Other existing numerical settings and spectral
approximations are retained from the supplied script.
Electrical/thermal fractions are an allocation at one-sun AM0 maximum power,
as expected by the existing optimiser; they are not universal material constants.
All absorbed power not exported electrically is assigned to heat; escaping
luminescence is neglected. Directional emissivity follows TOTAL absorption,
not heat. Normal-incidence absorption is the existing optimiser's approximation
to emissivity; this script does not change that optimiser.

Sources:
https://qpv-research-group.github.io/solcore-education/solar-cell-simulation/notebooks/1a-simple_cell.html
https://refractiveindex.info/?shelf=main&book=Si&page=Franta-300K
https://docs.solcore.solar/en/latest/Optics/tmm.html
https://docs.solcore.solar/en/latest/Solvers/depletion.html
"""

from io import StringIO
from importlib.metadata import version
from pathlib import Path
import pickle

import numpy as np
import yaml
from solcore import material
from solcore.absorption_calculator import OptiStack, calculate_rat
from solcore.analytic_solar_cells.depletion_approximation import iv_depletion, qe_depletion
from solcore.light_source import LightSource
from solcore.state import State
from solcore.structure import Junction, Layer

from constants import c, h, k, nm, q

import matplotlib.pyplot as plt
# Latex format in plots
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})



# Settings: use the SAME simulation pickle in this file and the optimiser.
here = Path(__file__).resolve().parent
nk_file = here / "Solar Panel Properties" / "Si-Properties-Franta-300K.yml"
simulation_pickle_file = here / "Simulation Data" / "simulation_data.pkl"
# In optimise_filters.py, point lookup_file at this output_file.
output_file = here / "Solar Panel Properties" / "silicon_single_junction_spectral_lookup.pkl"
plot_file = here / "Solar Panel Properties" / "silicon_spectral_properties.png"
make_plot = True
show_titles = True           # True: titles above each panel; False: no titles (e.g. use report captions)
separate_figures = False     # False: one 2x2 figure; True: four separate figures
separate_figure_size = (6.5, 4.5)   # inches per figure when separate_figures is True
# Separate figures are saved beside plot_file as <name>_optics.png, _partition, _temperature, _power.

material_name = "Si"          # MUST match the material in nk_file.
temperature_c = np.arange(-100.0, 300.0 + 0.1, 5.0)
# Solcore Education, Example 1a: 600 nm n-type emitter + 200 um p-type base.
wafer_thickness_m = 200.6e-6
emitter_thickness_m = 0.6e-6
emitter_donors_m3 = 1e27     # Tutorial value: 1e21 cm^-3, converted to SI.
base_acceptors_m3 = 1e22    # Tutorial value: 1e16 cm^-3, converted to SI.
emitter_hole_diffusion_length_m = 10e-6
base_electron_diffusion_length_m = 400e-6
relative_permittivity = 11.7
# Example 1a omits these losses; these are the documented DA solver defaults.
front_srv_m_s = 0.0
rear_srv_m_s = 0.0
shunt_resistance_ohm_m2 = 1e14
voltage_grid_size = 4001

# Below 200 nm the simple single-pair junction model is not used: absorbed
# power is counted as heat. Above the temperature-dependent band gap, one
# electron-hole pair per absorbed photon is assumed (no carrier multiplication).
electrical_minimum_nm = 200.0
bin_quadrature_points = 3    # Flat-spectrum bin averages on the input bin edges.

# User-selected simplification below the first YAML wavelength (31 nm here).
# These are arbitrary optical fractions, not measured data; A = 1 - R - T.
short_wave_reflectivity = 0.0
short_wave_transmissivity = 0.0


def load_wavelengths(path):
    """Copy the simulation's centres exactly; pickle inputs must be trusted."""
    with open(path, "rb") as stream:
        saved = pickle.load(stream)
    wavelengths = np.asarray(saved["wavelengths"], dtype=float)
    if wavelengths.ndim != 1 or wavelengths.size < 2 or not np.isfinite(wavelengths).all():
        raise ValueError("The simulation must contain a finite wavelength vector in nm.")
    if wavelengths[0] <= 0 or np.any(np.diff(wavelengths) <= 0):
        raise ValueError("Wavelengths must be positive and strictly increasing.")
    return wavelengths.copy()


def bin_sampling(wavelengths):
    """Infer edges as in the optimiser; three-point quadrature averages each bin."""
    edges = np.r_[wavelengths[0] - 0.5 * np.diff(wavelengths)[0],
                  0.5 * (wavelengths[:-1] + wavelengths[1:]),
                  wavelengths[-1] + 0.5 * np.diff(wavelengths)[-1]]
    if edges[0] <= 0:
        raise ValueError("The inferred first wavelength-bin edge must be positive.")
    points, weights = np.polynomial.legendre.leggauss(bin_quadrature_points)
    nodes = 0.5 * (edges[:-1, None] + edges[1:, None]) + 0.5 * np.diff(edges)[:, None] * points
    return nodes, weights / 2.0, edges


def read_nk_yaml(path):
    """refractiveindex.info tabulated nk uses micrometres, dimensionless n,k."""
    with open(path, encoding="utf-8") as stream:
        document = yaml.safe_load(stream)
    blocks = [block for block in document.get("DATA", []) if block.get("type") == "tabulated nk"]
    if len(blocks) != 1:
        raise ValueError("Choose a YAML file containing one 'tabulated nk' block, not a formula-only file.")
    values = np.loadtxt(StringIO(blocks[0]["data"]), ndmin=2)
    if values.shape[1] != 3 or values.shape[0] < 2 or not np.isfinite(values).all():
        raise ValueError("Expected finite wavelength, n, k columns.")
    wavelengths, n, extinction = values.T
    wavelengths = wavelengths * 1000.0
    if wavelengths[0] <= 0 or np.any(np.diff(wavelengths) <= 0) or np.any(n <= 0) or np.any(extinction < 0):
        raise ValueError("Expected increasing positive wavelengths, n > 0 and k >= 0.")
    return wavelengths, n, extinction, document


def sample_nk(nodes_nm, bin_edges_nm):
    wavelengths, n_data, k_data, document = read_nk_yaml(nk_file)
    flat = nodes_nm.ravel()
    if bin_edges_nm[-1] > wavelengths[-1]:
        raise ValueError(f"YAML ends at {wavelengths[-1]:g} nm; need {bin_edges_nm[-1]:g} nm. Supply wider data.")
    if bin_edges_nm[0] < wavelengths[0] and wavelengths[0] > 50:
        raise ValueError("Only a missing region below 50 nm may use the arbitrary short-wave extension.")
    n = np.interp(flat, wavelengths, n_data)
    extinction = np.interp(flat, wavelengths, k_data)
    short = flat < wavelengths[0]
    if np.any(n <= 0) or np.any(extinction < 0):
        raise ValueError("The combined optical model has invalid n or k.")
    metadata = {"nk_file": str(nk_file), "nk_reference": document.get("REFERENCES", ""),
                "nk_temperature_K": document.get("CONDITIONS", {}).get("temperature"),
                "yaml_range_nm": [float(wavelengths[0]), float(wavelengths[-1])],
                "short_wave_assumption": {"below_nm": float(wavelengths[0]), "R": short_wave_reflectivity,
                                          "T": short_wave_transmissivity, "A": 1 - short_wave_reflectivity - short_wave_transmissivity}}
    yaml_bins = (bin_edges_nm[:-1] >= wavelengths[0]) & (bin_edges_nm[1:] <= wavelengths[-1])
    return n, extinction, short, yaml_bins, metadata


def calculate_optics(wavelengths_nm, n, extinction):
    """Vacuum/wafer/vacuum; retain rear reflection and finite transmission."""
    stack = OptiStack([[wafer_thickness_m / nm, wavelengths_nm, n, extinction]], no_back_reflection=False)
    result = calculate_rat(stack, wavelengths_nm, angle=0, pol="u", coherent=False, coherency_list=["i"], no_back_reflection=False)
    return tuple(np.asarray(result[key], dtype=float) for key in ("R", "A", "T"))


def make_junction(temperature_k):
    """Example 1a n-on-p junction with fixed, example-supplied diffusion lengths."""
    reference = material(material_name)(T=300.0)
    states_scale = (temperature_k / 300.0) ** 1.5
    common = dict(T=temperature_k, Nc=reference.Nc * states_scale,
                  Nv=reference.Nv * states_scale, relative_permittivity=relative_permittivity)
    emitter = material(material_name)(Nd=emitter_donors_m3,
                                     hole_diffusion_length=emitter_hole_diffusion_length_m, **common)
    base = material(material_name)(Na=base_acceptors_m3,
                                  electron_diffusion_length=base_electron_diffusion_length_m, **common)
    junction = Junction([Layer(emitter_thickness_m, emitter, role="emitter"),
                         Layer(wafer_thickness_m - emitter_thickness_m, base, role="base")],
                        kind="DA", T=temperature_k, sn=front_srv_m_s, sp=rear_srv_m_s,
                        R_shunt=shunt_resistance_ohm_m2)
    return junction, base.band_gap


def absorption_profile(wavelengths_m, n, extinction, absorption, band_gap_j):
    """Phase-averaged forward/backward Beer-Lambert profile, normalised to TMM A.

    Its depth integral equals the slab absorption. The band-gap gate prevents
    IR/free-carrier/lattice absorption from being counted as photogeneration.
    A hard edge is an approximation for indirect-gap Si and its phonon tails.
    """
    coefficient = 4 * np.pi * extinction / wavelengths_m
    back_reflectivity = abs((n + 1j * extinction - 1) / (n + 1j * extinction + 1)) ** 2
    passed = np.exp(-coefficient * wafer_thickness_m)
    denominator = -np.expm1(-coefficient * wafer_thickness_m) * (1 + back_reflectivity * passed)
    scale = np.divide(absorption * coefficient, denominator, out=np.zeros_like(absorption), where=denominator > 0)
    scale *= (h * c / wavelengths_m >= band_gap_j)

    def profile(depth_m):
        z = np.asarray(depth_m)[..., None]
        return scale * (np.exp(-coefficient * z) + back_reflectivity * np.exp(-coefficient * (2 * wafer_thickness_m - z)))

    return profile


def mpp_from_dark_iv(junction, temperature_k, photocurrent):
    """DA superposition: J_out(V)=J_ph-J_dark(V), maximised for the AM0 J_ph.

    Use the same zero-bias collection approximation as Solcore DA.
    Series resistance and bias-dependent collection are omitted in this model.
    Solcore uses negative forward voltages and negative dark currents for this
    n-on-p junction. Convert them to a positive forward-voltage/current axis
    here; saved Vmpp, Jmpp and Voc remain positive generating magnitudes.
    """
    ni = junction[1].material.ni
    built_in_v = k * temperature_k / q * np.log(emitter_donors_m3 * base_acceptors_m3 / ni ** 2)
    voltage = np.linspace(0, built_in_v - 0.002, voltage_grid_size)
    iv_depletion(junction, State(T=temperature_k, internal_voltages=-voltage, light_iv=False))
    dark = -np.asarray(junction.current)
    if not np.isfinite(dark).all() or np.any(np.diff(dark) <= 0) or dark[-1] < photocurrent:
        raise ValueError("The dark IV does not bracket the photocurrent; check the junction parameters.")
    power = voltage * (photocurrent - dark)
    index = np.argmax(power)
    vmpp = voltage[index]
    jmpp = photocurrent - dark[index]
    voc = np.interp(photocurrent, dark, voltage)
    return vmpp, jmpp, voc


def build_lookup():
    if not 0 < emitter_thickness_m < wafer_thickness_m:
        raise ValueError("The emitter thickness must lie inside the wafer.")
    if np.any(np.diff(temperature_c) <= 0) or np.any(temperature_c <= -273.15):
        raise ValueError("Temperatures must increase and be above absolute zero.")
    wavelengths = load_wavelengths(simulation_pickle_file)
    nodes, weights, edges = bin_sampling(wavelengths)
    flat_nm = nodes.ravel()
    n, extinction, short, yaml_bins, metadata = sample_nk(nodes, edges)
    r_nodes, a_nodes, t_nodes = calculate_optics(flat_nm, n, extinction)
    # Override this unmeasured interval with the requested arbitrary values.
    r_nodes[short] = short_wave_reflectivity
    t_nodes[short] = short_wave_transmissivity
    a_nodes[short] = 1 - short_wave_reflectivity - short_wave_transmissivity

    def average(values):
        return np.asarray(values).reshape(nodes.shape) @ weights

    r, a, transmission = map(average, (r_nodes, a_nodes, t_nodes))
    # The sim supplies W per bin, not W/nm. Use matching bin powers for AM0
    # normalisation; within each bin our response assumes uniform spectral power.
    source = LightSource(source_type="standard", version="AM0")
    _, spectral_power = source.spectrum(x=wavelengths * nm, output_units="power_density_per_m")
    reference_power = spectral_power * np.diff(edges) * nm
    shape = (temperature_c.size, wavelengths.size)
    eqe = np.zeros(shape)
    current_response = np.zeros(shape)
    electrical = np.zeros(shape)
    metrics = {key: [] for key in ("vmpp_V", "jmpp_A_m2", "isc_A_m2", "voc_V", "pmpp_W_m2", "band_gap_eV")}
    smallest_gap = min(make_junction(t + 273.15)[1] for t in temperature_c)
    active = (~short) & (flat_nm >= electrical_minimum_nm) & (h * c / (flat_nm * nm) >= smallest_gap)
    if active.sum() < 2:
        raise ValueError("The wavelength grid does not cover the electrical response.")

    for row, temp in enumerate(temperature_c):
        kelvin = temp + 273.15
        junction, band_gap = make_junction(kelvin)
        junction.absorbed = absorption_profile(flat_nm[active] * nm, n[active], extinction[active], a_nodes[active], band_gap)
        qe_depletion(junction, State(T=kelvin, wavelength=flat_nm[active] * nm, da_mode="green"))
        eqe_nodes = np.zeros_like(flat_nm)
        eqe_nodes[active] = junction.eqe(flat_nm[active] * nm)
        if not np.isfinite(eqe_nodes).all() or np.min(eqe_nodes) < -1e-8 or np.max(eqe_nodes - a_nodes) > 1e-4:
            raise ValueError(f"Unphysical EQE at {temp:g} C; check generation profile and collection model.")
        # Remove numerical roundoff only; the check above prevents hiding a bad model.
        eqe_nodes = np.clip(eqe_nodes, 0, a_nodes)
        eqe[row] = average(eqe_nodes)
        current_response[row] = average(q * eqe_nodes * flat_nm * nm / (h * c))
        jph = float(np.dot(current_response[row], reference_power))
        if jph <= 0:
            raise ValueError("The calculated AM0 photocurrent must be positive.")
        vmpp, jmpp, voc = mpp_from_dark_iv(junction, kelvin, jph)
        pmpp = vmpp * jmpp
        # EQE*q*Vmpp/E alone gives Jph*Vmpp, NOT Jmpp*Vmpp. The effective
        # voltage Pmpp/Jph makes the spectral allocation integrate to Pmpp.
        electrical[row] = current_response[row] * pmpp / jph
        for key, value in zip(metrics, (vmpp, jmpp, jph, voc, pmpp, band_gap / q)):
            metrics[key].append(value)
        print(f"{temp:6.1f} C: Voc={metrics['voc_V'][-1]:.3f} V, AM0 Pmpp={pmpp:.2f} W/m²")

    optical_shape = lambda values: np.broadcast_to(values, shape).copy()
    metadata.update({"schema_version": 1, "material": material_name, "model": "Solcore TMM + DA bare n-on-p homojunction", "solcore_version": version("solcore"),
                     "cell_parameter_source": "https://qpv-research-group.github.io/solcore-education/solar-cell-simulation/notebooks/1a-simple_cell.html",
                     "loss_parameter_source": "https://docs.solcore.solar/en/latest/Solvers/depletion.html",
                     "junction_polarity": "n-on-p",
                     "wavelength_source": str(simulation_pickle_file), "wavelength_unit": "nm", "temperature_unit": "degC",
                     "binning": "Gauss-Legendre flat-power bin averages; exact simulation centres and inferred edges",
                     "optical_temperature_assumption": "n,k fixed at YAML temperature; phase-averaged finite slab",
                     "electrical_absorptivity_definition": "AM0 allocation: current_response_A_per_W * Pmpp/Jph",
                     "thermal_absorptivity_definition": "total_absorptivity - electrical_absorptivity_mpp",
                     "operating_point_assumption": "Spectral conversion fractions are calibrated to one-sun AM0 MPP; changes of MPP under a filtered spectrum are not recalculated by the existing optimiser.",
                     "electrical_metric_sign_convention": "Positive generating magnitudes; Solcore n-on-p forward voltage and dark current are negated for the MPP calculation.",
                     "wafer_thickness_m": wafer_thickness_m, "emitter_thickness_m": emitter_thickness_m,
                     "emitter_donors_m3": emitter_donors_m3, "base_acceptors_m3": base_acceptors_m3,
                     "emitter_hole_diffusion_length_m": emitter_hole_diffusion_length_m,
                     "base_electron_diffusion_length_m": base_electron_diffusion_length_m,
                     "relative_permittivity": relative_permittivity,
                     "front_srv_m_s": front_srv_m_s, "rear_srv_m_s": rear_srv_m_s,
                     "shunt_resistance_ohm_m2": shunt_resistance_ohm_m2,
                     "illumination_for_operating_point": "AM0, one sun, on the simulation bins",
                     "electrical_minimum_nm": electrical_minimum_nm,
                     "assumptions": "Bare wafer; normal incidence for incident power; fixed example diffusion lengths; no series resistance, coverglass, contacts, encapsulation, radiation damage or luminescent escape."})
    table = {"temperature_C": temperature_c.copy(), "wavelength_nm": wavelengths, "wavelength_bin_edges_nm": edges,
             "total_absorptivity": optical_shape(a), "reflectivity": optical_shape(r), "transmissivity": optical_shape(transmission),
             "electrical_absorptivity_mpp": electrical,
             "thermal_absorptivity_mpp": a[None, :] - electrical, "external_quantum_efficiency": eqe,
             "current_response_A_per_W": current_response, "valid_optics_mask": yaml_bins,
             "valid_electrical_mask": active.reshape(nodes.shape).any(axis=1),
             "reference_power_W_m2_per_bin": reference_power, "metadata": metadata,
             **{key: np.asarray(value) for key, value in metrics.items()}}
    validate_lookup(table)
    return table


def validate_lookup(table):
    r, a, transmission = (table[key] for key in ("reflectivity", "total_absorptivity", "transmissivity"))
    electrical, thermal = table["electrical_absorptivity_mpp"], table["thermal_absorptivity_mpp"]
    for values in (r, a, transmission, electrical, thermal):
        if not np.isfinite(values).all() or values.min() < -1e-8 or values.max() > 1 + 1e-8:
            raise ValueError("A spectral fraction is outside [0, 1] or is non-finite.")
    if not np.allclose(r + a + transmission, 1, atol=1e-8, rtol=0):
        raise ValueError("R + A + T does not equal one.")
    if not np.allclose(electrical + thermal, a, atol=1e-10, rtol=0):
        raise ValueError("The electrical/heat partition does not conserve energy.")
    allocated_power = electrical @ table["reference_power_W_m2_per_bin"]
    if not np.allclose(allocated_power, table["pmpp_W_m2"], rtol=1e-10, atol=1e-9):
        raise ValueError("The reference electrical spectrum does not integrate to Pmpp.")


def _plot_optics(ax, table, row):
    wavelengths = table["wavelength_nm"]
    for key, label in (("reflectivity", r"$\rho(\lambda)$"), ("total_absorptivity", r"$\alpha(\lambda)$"), ("transmissivity", r"$\tau(\lambda)$")):
        ax.semilogx(wavelengths, table[key][row], label=label)
    if show_titles:
        ax.set_title(r"Optical Properties of a Bare Silicon Wafer (Fixed $n$, $k$)")
    ax.set_xlabel(r"Wavelength, $\lambda$ (nm)")
    ax.set_ylabel("Fraction of incident power")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)


def _plot_partition(ax, table, row):
    wavelengths = table["wavelength_nm"]
    for key, label in (("total_absorptivity", r"Total, $\alpha(\lambda)$"), ("electrical_absorptivity_mpp", r"Electrical, $\alpha_{\mathrm{el}}(\lambda)$ (AM0 MPP)"),
                       ("thermal_absorptivity_mpp", r"Heat, $\alpha_{\mathrm{th}}(\lambda)$ (AM0 MPP)")):
        ax.semilogx(wavelengths, table[key][row], label=label)
    if show_titles:
        ax.set_title(rf"Electrical and Thermal Absorptivity at ${table['temperature_C'][row]:g}\,^{{\circ}}\mathrm{{C}}$")
    ax.set_xlabel(r"Wavelength, $\lambda$ (nm)")
    ax.set_ylabel("Fraction of incident power")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)


def _plot_temperature(ax, table, row):
    wavelengths = table["wavelength_nm"]
    for index in (0, row, len(table["temperature_C"]) - 1):
        ax.plot(wavelengths, table["electrical_absorptivity_mpp"][index], label=rf"$T={table['temperature_C'][index]:g}\,^{{\circ}}\mathrm{{C}}$")
    ax.set_xlim(200, 1300)
    if show_titles:
        ax.set_title(r"Electrical Absorptivity at Different Cell Temperatures")
    ax.set_xlabel(r"Wavelength, $\lambda$ (nm)")
    ax.set_ylabel(r"$\alpha_{\mathrm{el}}$")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)


def _plot_power(ax, table, row):
    ax.plot(table["temperature_C"], table["pmpp_W_m2"])
    ax.set(xlabel=r"Cell temperature, $T$ ($^{\circ}$C)", ylabel=r"AM0 maximum power, $P_{\mathrm{mpp}}$ ($\mathrm{W\,m^{-2}}$)")
    if show_titles:
        ax.set_title("Maximum Power Point vs. Cell Temperature (AM0)")
    ax.grid(alpha=0.25)


def plot_lookup(table, path):
    """Draw the four panels as one 2x2 figure, or as four separate figures."""
    row = np.argmin(abs(table["temperature_C"] - 25))
    panels = (("optics", _plot_optics), ("partition", _plot_partition),
              ("temperature", _plot_temperature), ("power", _plot_power))
    path.parent.mkdir(parents=True, exist_ok=True)
    if separate_figures:
        for name, draw in panels:
            fig, ax = plt.subplots(figsize=separate_figure_size)
            draw(ax, table, row)
            fig.tight_layout()
            fig.savefig(path.with_name(f"{path.stem}_{name}{path.suffix}"), dpi=160)
            plt.close(fig)
    else:
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        for ax, (_, draw) in zip(axes.ravel(), panels):
            draw(ax, table, row)
        fig.tight_layout()
        fig.savefig(path, dpi=160)
        plt.close(fig)


def main():
    table = build_lookup()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "wb") as stream:
        pickle.dump(table, stream, protocol=pickle.HIGHEST_PROTOCOL)
    if make_plot:
        plot_lookup(table, plot_file)
    print(f"Saved {output_file}: {len(table['wavelength_nm'])} bins, {len(table['temperature_C'])} temperatures.")


if __name__ == "__main__":
    main()