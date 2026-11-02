"""Initialise the opaque aluminium panel spectrum for optimise_filters.py.

Run this file before the optimiser, using the same Simulation Data/simulation_data.pkl
as panel_initialise.py. The supplied CSV contains the 206 reflectance values
from Table 2 of Rakic (1995), exported from the supplied paper, and the digital
n,k values from the CC0 refractiveindex.info database. The paper's reflectance
column was checked independently against the normal-incidence Fresnel equation:

    R = ((n - 1)**2 + k**2) / ((n + 1)**2 + k**2)

Use the printed R values, interpolate R linearly in wavelength, then average each bin using
three-point Gauss-Legendre quadrature and a flat within-bin incident spectrum.
For an opaque passive panel, alpha = 1 - R, tau = 0, beta = 0 and epsilon = alpha.

This is an ideal smooth, unoxidised aluminium reference. It is not a measured
Al-6061 panel with its native oxide, surface finish or anodising. Normal-incidence
absorption is also the optimiser's approximation to hemispherical emissivity.
All optical properties are held independent of temperature. The 1 mm thermal
panel thickness is not used to fit an optical spectrum: opacity is assumed.

Dependencies: numpy; pandas to unpickle the existing simulation's DataFrames;
matplotlib only if make_plot is True. Solcore and network access are not needed.

Paper: https://doi.org/10.1364/AO.34.004755
Data: https://refractiveindex.info/database/data/main/Al/nk/Rakic.yml
Database: https://doi.org/10.1038/s41597-023-02898-2
"""

from __future__ import annotations

import csv
import hashlib
import pickle
from pathlib import Path

import numpy as np

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


# Edit these settings, then run this script. Paths are relative to this file.
here = Path(__file__).resolve().parent
simulation_pickle_file = here / "Simulation Data" / "simulation_data.pkl"
source_csv_file = here / "Aluminium Panel Properties" / "aluminium_rakic1995.csv"
output_file = here / "Aluminium Panel Properties" / "aluminium_spectral_lookup.pkl"
output_csv_file = here / "Aluminium Panel Properties" / "aluminium_spectral_properties.csv"
plot_file = here / "Aluminium Panel Properties" / "aluminium_spectral_properties.pdf"
make_plot = True
bin_quadrature_points = 3

PAPER_DOI = "10.1364/AO.34.004755"
DATABASE_DOI = "10.1038/s41597-023-02898-2"
DATA_URL = "https://refractiveindex.info/database/data/main/Al/nk/Rakic.yml"
SOURCE_YAML_SHA256 = "e93b31fe821551836e5b000949c67c633b8a1ca1f2c5549f4db0f1c03b24844d"


def infer_bin_edges(wavelengths_nm: np.ndarray) -> np.ndarray:
    """Use the same midpoint/extrapolated edges as the silicon initialiser."""
    wavelengths = np.asarray(wavelengths_nm, dtype=float)
    if (wavelengths.ndim != 1 or wavelengths.size < 2
            or not np.isfinite(wavelengths).all() or wavelengths[0] <= 0
            or np.any(np.diff(wavelengths) <= 0)):
        raise ValueError("Wavelength centres must be finite, positive and strictly increasing.")
    edges = np.r_[wavelengths[0] - 0.5 * (wavelengths[1] - wavelengths[0]),
                  0.5 * (wavelengths[:-1] + wavelengths[1:]),
                  wavelengths[-1] + 0.5 * (wavelengths[-1] - wavelengths[-2])]
    if edges[0] <= 0:
        raise ValueError("The inferred first wavelength-bin edge must be positive.")
    return edges


def read_aluminium_data(path: Path = source_csv_file) -> tuple[np.ndarray, np.ndarray]:
    """Read Table 2 R and cross-check it against independently tabulated n,k."""
    path = Path(path)
    data = np.genfromtxt(path, names=True, delimiter=",", dtype=float, encoding="utf-8")
    required = {"wavelength_nm", "n", "k", "reflectivity_table2", "reflectivity_from_nk"}
    if not required.issubset(data.dtype.names or ()) or data.ndim != 1 or data.size < 2:
        raise ValueError(f"{path} must contain at least two rows with columns {sorted(required)}.")
    wavelengths, n, extinction, reflected = (np.asarray(data[name], dtype=float)
                                           for name in ("wavelength_nm", "n", "k", "reflectivity_table2"))
    if (not np.isfinite([wavelengths, n, extinction, reflected]).all()
            or wavelengths[0] <= 0 or np.any(np.diff(wavelengths) <= 0)
            or np.any(n <= 0) or np.any(extinction < 0)):
        raise ValueError("Aluminium data require increasing positive wavelengths, n > 0 and k >= 0.")
    fresnel = ((n - 1.0) ** 2 + extinction ** 2) / ((n + 1.0) ** 2 + extinction ** 2)
    if not np.allclose(data["reflectivity_from_nk"], fresnel, rtol=0, atol=1e-12):
        raise ValueError("The CSV's Fresnel check does not agree with its n,k values.")
    # R and n,k were rounded independently in the paper. The maximum observed
    # discrepancy is 7.925e-6 in absolute reflectivity; 1e-5 admits that rounding.
    if not np.allclose(reflected, fresnel, rtol=0, atol=1e-5):
        raise ValueError("The Table 2 reflectance does not agree with the independent n,k check.")
    if np.any(reflected < 0) or np.any(reflected >= 1):
        raise ValueError("The aluminium model requires 0 <= reflectivity < 1.")
    return wavelengths, reflected


def average_reflectivity(edges_nm: np.ndarray, source_nm: np.ndarray,
                         source_reflectivity: np.ndarray, points: int = 3) -> np.ndarray:
    """Average linearly interpolated R inside each bin; never extrapolate."""
    edges = np.asarray(edges_nm, dtype=float)
    wavelengths = np.asarray(source_nm, dtype=float)
    reflected = np.asarray(source_reflectivity, dtype=float)
    if (edges.ndim != 1 or edges.size < 2 or not np.isfinite(edges).all()
            or edges[0] <= 0 or np.any(np.diff(edges) <= 0)):
        raise ValueError("Bin edges must be finite, positive and strictly increasing.")
    if (wavelengths.ndim != 1 or wavelengths.size < 2 or reflected.shape != wavelengths.shape
            or not np.isfinite([wavelengths, reflected]).all()
            or wavelengths[0] <= 0 or np.any(np.diff(wavelengths) <= 0)
            or np.any(reflected < 0) or np.any(reflected >= 1)):
        raise ValueError("Invalid source wavelength/reflectivity arrays.")
    if not isinstance(points, (int, np.integer)) or points < 1:
        raise ValueError("The quadrature-point count must be a positive integer.")
    if edges[0] < wavelengths[0] or edges[-1] > wavelengths[-1]:
        raise ValueError(f"Requested bin edges {edges[0]:g} .. {edges[-1]:g} nm exceed "
                         f"the aluminium data range {wavelengths[0]:g} .. {wavelengths[-1]:g} nm.")
    nodes, weights = np.polynomial.legendre.leggauss(points)
    queries = (0.5 * (edges[:-1, None] + edges[1:, None])
               + 0.5 * np.diff(edges)[:, None] * nodes)
    return np.interp(queries, wavelengths, reflected) @ (0.5 * weights)


def build_lookup(wavelengths_nm: np.ndarray, data_path: Path = source_csv_file) -> dict:
    """Build temperature-independent fractions on the exact simulation grid."""
    wavelengths = np.asarray(wavelengths_nm, dtype=float).copy()
    edges = infer_bin_edges(wavelengths)
    source_nm, source_reflectivity = read_aluminium_data(data_path)
    reflected = average_reflectivity(edges, source_nm, source_reflectivity, bin_quadrature_points)
    absorbed = 1.0 - reflected
    table = {
        "schema_version": 1,
        "material": "Aluminium",
        "wavelength_nm": wavelengths,
        "wavelength_bin_edges_nm": edges,
        "total_absorptivity": absorbed,
        "reflectivity": reflected,
        "transmissivity": np.zeros_like(absorbed),
        "emissivity": absorbed.copy(),
        "electrical_absorptivity_mpp": np.zeros_like(absorbed),
        "metadata": {
            "paper_doi": PAPER_DOI,
            "database_doi": DATABASE_DOI,
            "source_url": DATA_URL,
            "source_retrieved": "2026-11-02",
            "source_yaml_sha256": SOURCE_YAML_SHA256,
            "source_pdf_sha256": "27287fe55f5821774a97fa496ea3b3c8e3dc0eb00fa7a7d868978d4a84c98f4f",
            "source_csv_sha256": hashlib.sha256(Path(data_path).read_bytes()).hexdigest(),
            "source_points": int(source_nm.size),
            "source_range_nm": [float(source_nm[0]), float(source_nm[-1])],
            "wavelength_unit": "nm",
            "fractions_unit": "dimensionless",
            "reflectivity_definition": "Printed normal-incidence R from Table 2, checked against Fresnel R from tabulated n,k",
            "table_location": "Rakic (1995), Table 2, printed pp. 4762-4765",
            "interpolation": "Linear interpolation of R in wavelength, not interpolation of n,k",
            "binning": f"{bin_quadrature_points}-point Gauss-Legendre flat-power bin averages",
            "assumptions": "Smooth unoxidised opaque aluminium; beta=0; optics fixed with temperature; normal-incidence absorptivity approximates hemispherical emissivity; no extrapolation",
            "license": "Source n,k data: CC0 1.0, refractiveindex.info database",
        },
    }
    validate_lookup(table)
    return table


def validate_lookup(table: dict, expected_wavelengths_nm: np.ndarray | None = None) -> None:
    """Reject an incompatible grid, missing data or nonphysical fractions."""
    required = {"schema_version", "material", "wavelength_nm", "wavelength_bin_edges_nm",
                "total_absorptivity", "reflectivity", "transmissivity", "emissivity",
                "electrical_absorptivity_mpp", "metadata"}
    missing = required.difference(table)
    if missing:
        raise ValueError(f"Aluminium lookup is missing: {sorted(missing)}")
    if table["schema_version"] != 1 or table["material"] != "Aluminium":
        raise ValueError("Unsupported aluminium lookup schema or material.")
    wavelengths = np.asarray(table["wavelength_nm"], dtype=float)
    edges = infer_bin_edges(wavelengths)
    if not np.array_equal(np.asarray(table["wavelength_bin_edges_nm"], dtype=float), edges):
        raise ValueError("Aluminium lookup bin edges do not match its wavelength centres.")
    if expected_wavelengths_nm is not None and not np.array_equal(wavelengths, np.asarray(expected_wavelengths_nm)):
        raise ValueError("The aluminium lookup grid does not match the simulation; rerun aluminium_initialise.py.")
    keys = ("total_absorptivity", "reflectivity", "transmissivity", "emissivity", "electrical_absorptivity_mpp")
    fractions = [np.asarray(table[key], dtype=float) for key in keys]
    if any(values.shape != wavelengths.shape or not np.isfinite(values).all()
           or values.min() < 0 or values.max() > 1 for values in fractions):
        raise ValueError("Aluminium fractions must be finite vectors in [0, 1] on the wavelength grid.")
    a, r, transmission, emissivity, electrical = fractions
    if not np.allclose(a + r + transmission, 1.0, rtol=0, atol=1e-12):
        raise ValueError("Aluminium A + R + T does not equal one.")
    if np.any(r >= 1) or np.any(transmission != 0) or np.any(electrical != 0):
        raise ValueError("The aluminium panel must be opaque and passive, with reflectivity < 1.")
    if not np.array_equal(emissivity, a):
        raise ValueError("Aluminium emissivity must equal total absorptivity in this model.")


def load_aluminium_lookup(path: Path, expected_wavelengths_nm: np.ndarray | None = None) -> dict:
    """Load a locally generated, trusted pickle; the optimiser calls this."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing aluminium lookup: {path}. Run aluminium_initialise.py first.")
    with path.open("rb") as stream:
        table = pickle.load(stream)
    validate_lookup(table, expected_wavelengths_nm)
    return table


def load_simulation_wavelengths(path: Path) -> np.ndarray:
    """Read only wavelengths from the existing trusted simulation pickle."""
    with Path(path).open("rb") as stream:
        saved = pickle.load(stream)
    wavelengths = np.asarray(saved["wavelengths"], dtype=float)
    infer_bin_edges(wavelengths)
    return wavelengths.copy()


def save_lookup(table: dict, pickle_path: Path, csv_path: Path) -> None:
    """Save a lookup for the optimiser and a readable CSV of the same bins."""
    validate_lookup(table)
    pickle_path, csv_path = Path(pickle_path), Path(csv_path)
    pickle_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with pickle_path.open("wb") as stream:
        pickle.dump(table, stream, protocol=pickle.HIGHEST_PROTOCOL)
    edges = table["wavelength_bin_edges_nm"]
    keys = ("total_absorptivity", "reflectivity", "transmissivity", "emissivity", "electrical_absorptivity_mpp")
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("wavelength_nm", "lower_edge_nm", "upper_edge_nm", "absorptivity", "reflectivity", "transmissivity", "emissivity", "electrical_conversion"))
        writer.writerows(zip(table["wavelength_nm"], edges[:-1], edges[1:], *(table[key] for key in keys)))


def plot_lookup(table: dict, path: Path) -> None:
    """Save a vector plot"""

    fig, ax = plt.subplots(figsize=(7, 4))
    for key, label in (("reflectivity", r"$\rho(\lambda)$"), ("total_absorptivity", r"$\alpha(\lambda)=\varepsilon(\lambda)$"), ("transmissivity", r"$\tau(\lambda)=0$")):
        ax.semilogx(table["wavelength_nm"], table[key], linewidth=1.0, label=label)
    ax.set(xlabel=r"Wavelength, $\lambda$ (nm)", ylabel="Fraction", ylim=(-0.02, 1.02))
    ax.set_xlim(table["wavelength_bin_edges_nm"][[0, -1]])
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    wavelengths = load_simulation_wavelengths(simulation_pickle_file)
    table = build_lookup(wavelengths)
    save_lookup(table, output_file, output_csv_file)
    if make_plot:
        plot_lookup(table, plot_file)
    lo, hi = table["wavelength_bin_edges_nm"][[0, -1]]
    print(f"Saved aluminium lookup: {wavelengths.size} bins covering {lo:g} .. {hi:g} nm.")
    print(f"Source: Rakic (1995), DOI {PAPER_DOI}; no wavelength extrapolation.")
    print(f"Lookup: {output_file}")
    print(f"CSV: {output_csv_file}")


if __name__ == "__main__":
    main()
