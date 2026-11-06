"""Define and plot the simplified wafer and its opaque backing.

Keep both test scripts in Simple Optimiser Test, below constants.py.
Run it directly to plot the properties without running the optimiser.
"""

from pathlib import Path
import sys

# Find the shared modules one folder above this test folder.
sys.path.insert(1, str(Path(__file__).resolve().parent.parent))

import numpy as np
from constants import minute
import matplotlib.pyplot as plt

# latex in plots
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
    "patch.linewidth": 2.0 # increases line width in the stair plotting
})

# Operating conditions. The wafer and backing form one thermal node.
area_m2 = 0.02
mass_kg = 0.06349
heat_capacity_j_per_kg_k = 869.0
initial_temperature_k = 293.15
target_temperature_k = 293.15
temperature_reference_k = 10.0
temperature_min_k = 200.0
temperature_max_k = 400.0
duration_s = 30 * minute
time_step_s = 10.0
figure_dir = Path(__file__).resolve().parent / "simple_validation"

# Interval limits in nm, number of bins and total incident power in watts.
# Power is uniform in wavelength within each interval and constant in time.
spectral_regions = [(300, 400, 6, 2), (400, 900, 24, 10), (900, 2500, 24, 8), (2500, 50000, 66, 0)]
wavelength_edges_nm, incident_power_w = [], []
for lower, upper, count, power in spectral_regions:
    edges = np.geomspace(lower, upper, count + 1)
    wavelength_edges_nm.extend(edges[:-1])
    incident_power_w.extend(power * np.diff(edges) / (upper - lower))
wavelength_edges_nm = np.array(wavelength_edges_nm + [50000.0])
wavelength_nm = np.sqrt(wavelength_edges_nm[:-1] * wavelength_edges_nm[1:])
bandwidth_nm = np.diff(wavelength_edges_nm)
incident_power_w = np.array(incident_power_w)
time_s = np.arange(0, duration_s + time_step_s / 2, time_step_s)
interval_s = np.r_[np.diff(time_s), 0.0]  # No interval follows the final sample.

# Electricity is produced only in this band. Wafer transmission elsewhere
# reaches an ideal absorbing backing and therefore becomes heat.
conversion_band = (wavelength_nm >= 400) & (wavelength_nm < 900)
wafer_absorptivity = np.where(conversion_band, 0.8, 0.2)
wafer_reflectivity = np.full(wavelength_nm.size, 0.2)
wafer_transmissivity = 1 - wafer_absorptivity - wafer_reflectivity
electrical_conversion = np.where(conversion_band, 0.25, 0.0)
wafer_thermal_absorptivity = wafer_absorptivity - electrical_conversion

# Independently prescribed properties of the equivalent opaque panel.
panel_absorptivity = np.full(wavelength_nm.size, 0.8)
panel_reflectivity = np.full(wavelength_nm.size, 0.2)
panel_transmissivity = np.zeros(wavelength_nm.size)


def plot_properties(properties=None):
    """Save the optical properties and the electrical/thermal partition as PDFs."""
    # A saved solution can supply its original properties when replotting.
    if properties is None:
        properties = {
            "wavelength_edges_nm": wavelength_edges_nm,
            "wafer_absorptivity": wafer_absorptivity,
            "wafer_reflectivity": wafer_reflectivity,
            "wafer_transmissivity": wafer_transmissivity,
            "electrical_conversion": electrical_conversion,
            "wafer_thermal_absorptivity": wafer_thermal_absorptivity,
            "panel_absorptivity": panel_absorptivity,
            "panel_reflectivity": panel_reflectivity,
            "panel_transmissivity": panel_transmissivity
        }
    figure_dir.mkdir(exist_ok=True)
    edges_nm = properties["wavelength_edges_nm"]
    fig, axes = plt.subplots(2, 1, figsize=(6.3, 6.5), sharex=True)
    # Match the report: r is the receiving surface, p includes the backing.
    for ax, alpha, rho, tau, subscript, title in (
        (axes[0], properties["wafer_absorptivity"], properties["wafer_reflectivity"], properties["wafer_transmissivity"], "r", "Receiving surface"),
        (axes[1], properties["panel_absorptivity"], properties["panel_reflectivity"], properties["panel_transmissivity"], "p", "Receiving surface with backing")
    ):
        ax.stairs(alpha, edges_nm, baseline=None, linestyle="-", zorder=2, label=rf"$\alpha_{{{subscript}}}(\lambda)$")
        ax.stairs(rho, edges_nm, baseline=None, linestyle="--", zorder=4, label=rf"$\rho_{{{subscript}}}(\lambda)$")
        ax.stairs(tau, edges_nm, baseline=None, linestyle="-.", zorder=3, label=rf"$\tau_{{{subscript}}}(\lambda)$")
        # ax.stairs(alpha, edges_nm, baseline=None, label=rf"$\alpha_{{{subscript}}}(\lambda)$")
        # ax.stairs(rho, edges_nm, baseline=None, label=rf"$\rho_{{{subscript}}}(\lambda)$", linestyle="--")
        # ax.stairs(tau, edges_nm, baseline=None, label=rf"$\tau_{{{subscript}}}(\lambda)$")
        ax.set(xscale="log", xlim=(edges_nm[0], edges_nm[-1]), ylim=(-0.02, 1.02), ylabel="Fraction of incident power", title=title)
        ax.legend(ncol=3)
        ax.grid(alpha=0.25)
    axes[1].set_xlabel(r"Wavelength, $\lambda$ (nm)")
    fig.tight_layout()
    fig.savefig(figure_dir / "simple_panel_optical_properties.pdf", dpi=200, bbox_inches="tight")

    fig, ax = plt.subplots(figsize=(6.3, 4))
    # These fractions are relative to radiation reaching the receiving surface.
    ax.stairs(properties["electrical_conversion"], edges_nm, baseline=None, linestyle=":", zorder=4, label=r"Electrical, $\beta_r(\lambda)$")
    ax.stairs(properties["wafer_thermal_absorptivity"], edges_nm, baseline=None, linestyle="-", zorder=2, label=r"Wafer heat, $\alpha_{r,\mathrm{therm}}(\lambda)$")
    ax.stairs(properties["wafer_transmissivity"], edges_nm, baseline=None, linestyle="-.", zorder=3, label=r"Backing heat, $\tau_r(\lambda)$")
    ax.stairs(properties["panel_absorptivity"], edges_nm, baseline=None, linestyle="--", zorder=2, label=r"Total absorbed, $\alpha_p(\lambda)$")

    # ax.stairs(properties["electrical_conversion"], edges_nm, baseline=None, label=r"Electrical, $\beta_r(\lambda)$")
    # ax.stairs(properties["wafer_thermal_absorptivity"], edges_nm, baseline=None, label=r"Wafer heat, $\alpha_{r,\mathrm{therm}}(\lambda)$")
    # ax.stairs(properties["wafer_transmissivity"], edges_nm, baseline=None, label=r"Backing heat, $\tau_r(\lambda)$")
    # ax.stairs(properties["panel_absorptivity"], edges_nm, baseline=None, linestyle="--", label=r"Total absorbed, $\alpha_p(\lambda)$")
    ax.set(xscale="log", xlim=(edges_nm[0], edges_nm[-1]), ylim=(-0.02, 1.02), xlabel=r"Wavelength, $\lambda$ (nm)", ylabel="Fraction of incident power")
    ax.legend(ncol=2)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_dir / "simple_panel_absorption_partition.pdf", dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    plot_properties()
    plt.show()
