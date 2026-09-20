#!/usr/bin/env python3
"""Generate a temperature/wavelength lookup for a single-junction GaAs cell.

The wavelength grid is read verbatim from the simulation pickle.  Solcore is
used in two distinct ways:

1. A transfer-matrix calculation supplies R, A and T for a bare GaAs wafer
   stack at normal incidence in vacuum.  Solcore's bundled optical constants
   stop at 2.48 um, so beyond that R, A and T are calculated from tabulated
   far-infrared n,k with an incoherent thick-slab model.  Every column of the
   lookup is therefore physically calculated rather than extrapolated.
2. Solcore's depletion-approximation (DA) junction model supplies EQE and the
   AM0 maximum-power-point voltage at each cell temperature.

EQE is only an intermediate quantity.  The exported electrical power fraction
for incident monochromatic power is

    electrical_absorptivity_mpp = EQE * q * V_mpp / (h*c/lambda).

The residual absorbed fraction is reported as thermal absorptivity:

    thermal_absorptivity_mpp = total_absorptivity
                               - electrical_absorptivity_mpp.

This partition assumes operation at the AM0 maximum-power point and assigns
all absorbed power not exported electrically to heat.  Luminescent escape is
not modelled.
"""

from __future__ import annotations

import contextlib
import json
import os
import pickle
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml
from scipy.constants import c, e, h, k as boltzmann_constant
from scipy.interpolate import RegularGridInterpolator


# ---------------------------------------------------------------------------
# Model limits and baseline assumptions
# ---------------------------------------------------------------------------
# Solcore's bundled optical-constant tables only cover this interval, so its
# transfer-matrix result is used here and tabulated far-infrared n,k take over
# beyond it (see far_infrared_nk_file below).  valid_optics_mask marks the
# Solcore interval, but every column of the lookup is now physically
# calculated, not extrapolated.
solcore_optical_min_nm = 206.64028666666666
solcore_optical_max_nm = 2479.68344

# There is no useful GaAs photovoltaic response beyond the band edge.  A
# slightly wider solve interval is retained so that the calculated response
# naturally reaches zero rather than being truncated exactly at the edge.
electrical_solve_max_nm = 1205.0

# Solcore's bundled n and k data are not temperature-resolved.  Consequently,
# the optical calculation is performed once at this reference temperature and
# repeated along the temperature axis.  The electrical calculation is still
# repeated at every requested cell temperature.
reference_temperature_c = 25.0

# Baseline geometry and recombination assumptions for the generic bare cell.
# These values can be replaced in the user configuration below when measured
# or manufacturer-specific parameters are available.
default_wafer_thickness_um = 175.0
default_surface_recombination_velocity_m_s = 100.0
default_electron_minority_lifetime_s = 3e-6
default_hole_minority_lifetime_s = 2.5e-7


# ---------------------------------------------------------------------------
# Run configuration
#
# Edit these lower-case variables to change an input, model assumption, or
# output location.
# ---------------------------------------------------------------------------
# Resolve paths relative to this file, so the script can be launched from any
# working directory without changing where it searches for input data.
here = Path(__file__).resolve().parent

# Input pickle produced by the simulation-data listener.
pickle_file = here / "Simulation Data" / "simulation_data.pkl"

# Output path without an extension. This is where the pickle file will be written
output_base = here / "gaas_single_junction_spectral_lookup"

# Tabulated far-infrared optical constants, used beyond Solcore's range.
# This is a refractiveindex.info YAML file with a "tabulated nk" block:
#   Franta et al. 2026, undoped GaAs, n,k 0.895-400 um at 300 K
#   main/GaAs/nk/Franta-300K.yml
# Download once from
#   https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database
#     /main/database/data/main/GaAs/nk/Franta-300K.yml
# Companion files Franta-370K.yml and Franta-440K.yml exist if the optical
# constants are ever wanted temperature-resolved as well.
far_infrared_nk_file = here / "GaAs_Franta_300K_nk.yml"

# The summary figure is saved beside the lookup files.  Set show_plot to True
# when running on a desktop if the figure should also open interactively after
# it has been saved. Keeping it False is safer for headless/remote runs.
plot_file = here / "gaas_single_junction_spectral_properties.png"
show_plot = True
plot_dpi = 200

# Inclusive temperature grid used for the lookup table.
temperature_min_c = -50.0
temperature_max_c = 100.0
temperature_step_c = 5.0

# Generic device assumptions. These should be replaced by measured values
# when the exact cell construction is known.
wafer_thickness_um = default_wafer_thickness_um
surface_recombination_velocity_m_s = (
    default_surface_recombination_velocity_m_s
)
electron_minority_lifetime_s = default_electron_minority_lifetime_s
hole_minority_lifetime_s = default_hole_minority_lifetime_s

# Each temperature is independent, so several temperatures can be solved in
# parallel. Limiting the count avoids excessive memory use on large machines.
jobs = min(6, os.cpu_count() or 1)

# False keeps normal output to the two final saved-file paths. Set this to
# True only when Solcore diagnostics and per-temperature progress are needed.
verbose_solcore = False


# ---------------------------------------------------------------------------
# Small result container used by each temperature worker
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TemperatureResult:
    """Electrical results returned by one independent temperature solve."""

    temperature_c: float
    vmpp_v: float
    isc_a_m2: float
    voc_v: float
    pmpp_w_m2: float
    eqe: np.ndarray


# ---------------------------------------------------------------------------
# Input and grid construction
# ---------------------------------------------------------------------------
def load_simulation_wavelengths(simulation_pickle: str | Path) -> np.ndarray:
    """Load and validate the exact wavelength vector saved by the listener.

    Only the wavelength vector is read here. The irradiance arrays remain in
    the original simulation pickle and are not changed by this generator.
    """
    # pandas is used because it matches the writer used by the listener.
    data = pd.read_pickle(Path(simulation_pickle))

    # Fail early with a clear message if a different pickle schema is supplied.
    if not isinstance(data, dict) or "wavelengths" not in data:
        raise ValueError("Simulation pickle must contain a 'wavelengths' entry.")

    # Convert to a predictable numeric representation before validation.
    wavelengths_nm = np.asarray(data["wavelengths"], dtype=float)

    # Interpolation requires a finite, one-dimensional, strictly increasing
    # coordinate vector.
    if wavelengths_nm.ndim != 1 or wavelengths_nm.size < 2:
        raise ValueError("Simulation wavelengths must be a one-dimensional vector.")
    if not np.isfinite(wavelengths_nm).all():
        raise ValueError("Simulation wavelengths contain non-finite values.")
    if not np.all(np.diff(wavelengths_nm) > 0):
        raise ValueError("Simulation wavelengths must be strictly increasing.")
    return wavelengths_nm


def make_temperature_grid(
    minimum_c: float = -50.0,
    maximum_c: float = 100.0,
    step_c: float = 5.0,
) -> np.ndarray:
    """Return an inclusive, uniformly spaced Celsius temperature vector."""
    if step_c <= 0 or maximum_c < minimum_c:
        raise ValueError("Require step_c > 0 and maximum_c >= minimum_c.")

    # Constructing the grid by an integer count avoids the usual floating-point
    # endpoint surprises associated with np.arange(start, stop, float_step).
    count = int(round((maximum_c - minimum_c) / step_c))
    temperatures_c = minimum_c + step_c * np.arange(count + 1, dtype=float)
    if not np.isclose(temperatures_c[-1], maximum_c, rtol=0, atol=1e-10):
        raise ValueError("Temperature range must be an integer number of steps.")
    temperatures_c[-1] = maximum_c
    return temperatures_c


# ---------------------------------------------------------------------------
# Console-output control
# ---------------------------------------------------------------------------
@contextlib.contextmanager
def suppress_process_output(enabled: bool = True) -> Iterable[None]:
    """Temporarily redirect a process's stdout and stderr to the null device.

    Solcore initializes optional solvers when it is imported, and some of its
    numerical routines write directly to the operating-system file descriptors.
    Redirecting the descriptors, rather than only Python's ``sys.stdout``, also
    contains output produced by compiled Fortran/C extensions. Exceptions are
    not swallowed: the descriptors are restored and the exception propagates.
    """
    if not enabled:
        yield
        return

    # Flush buffered text before changing the underlying descriptors so that
    # output written before this context is not accidentally redirected.
    sys.stdout.flush()
    sys.stderr.flush()

    # Duplicate the original descriptors so they can always be restored in the
    # finally block, including when Solcore raises an exception.
    saved_stdout = os.dup(1)
    saved_stderr = os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        # File descriptors 1 and 2 are stdout and stderr respectively.
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        # Restore normal output before allowing control to leave the context.
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(saved_stdout, 1)
        os.dup2(saved_stderr, 2)
        os.close(saved_stdout)
        os.close(saved_stderr)
        os.close(devnull)


# ---------------------------------------------------------------------------
# Electrical device construction and solution
# ---------------------------------------------------------------------------
def _build_da_cell(
    temperature_k: float,
    surface_recombination_velocity_m_s: float,
    electron_minority_lifetime_s: float,
    hole_minority_lifetime_s: float,
):
    """Build a depletion-approximation version of Solcore's GaAs device."""
    # Imports are local so that Solcore is only initialized inside the quiet
    # context established by the calling function.
    from solcore import material
    from solcore.solar_cell import SolarCell, default_GaAs

    # Use the geometry and doping of Solcore's documented default GaAs cell.
    substrate = material("GaAs")(T=temperature_k)
    junction = default_GaAs(temperature_k)

    # Select the depletion-approximation electrical solver and apply the same
    # effective recombination velocity to the front and rear interfaces.
    junction.kind = "DA"
    junction.sn = float(surface_recombination_velocity_m_s)
    junction.sp = float(surface_recombination_velocity_m_s)

    # The DA model requires minority-carrier diffusion lengths. The emitter's
    # minority carriers are electrons and the base's minority carriers are
    # holes. Calculate both from the Einstein relation D = mu*k*T/q followed
    # by L = sqrt(D*tau). The configured lifetimes match Solcore's built-in
    # GaAs PDD fallback values unless the user changes them above.
    emitter = junction[1].material
    base = junction[2].material
    emitter.electron_diffusion_length = np.sqrt(
        emitter.electron_mobility
        * boltzmann_constant
        * temperature_k
        / e
        * electron_minority_lifetime_s
    )
    base.hole_diffusion_length = np.sqrt(
        base.hole_mobility
        * boltzmann_constant
        * temperature_k
        / e
        * hole_minority_lifetime_s
    )

    # Series resistance is set to zero because no specific contact or busbar
    # resistance was supplied. The returned object contains one junction and
    # an optically/electrically compatible GaAs substrate.
    return SolarCell(
        [junction],
        T=temperature_k,
        R_series=0,
        substrate=substrate,
    )


def _solve_one_temperature(
    temperature_c: float,
    electrical_wavelengths_nm: np.ndarray,
    surface_recombination_velocity_m_s: float,
    electron_minority_lifetime_s: float,
    hole_minority_lifetime_s: float,
    quiet: bool = True,
) -> TemperatureResult:
    """Run independent QE and illuminated-IV solves at one temperature."""
    # Solcore expects kelvin and wavelengths in metres for the solar-cell
    # solver, whereas the external lookup axes use degrees Celsius and nm.
    temperature_k = float(temperature_c + 273.15)
    wavelengths_m = np.asarray(electrical_wavelengths_nm, dtype=float) * 1e-9

    # Include the imports and object construction in this context because
    # Solcore emits optional-solver notices while it initializes its registry.
    with suppress_process_output(quiet):
        from solcore.light_source import LightSource
        from solcore.solar_cell_solver import solar_cell_solver

        # The operating voltage is spectrum dependent. It is evaluated under
        # the standard one-sun AM0 spectrum appropriate to an unfiltered cell
        # in space, separately at each cell temperature.
        light_source = LightSource(
            source_type="standard",
            version="AM0",
            x=wavelengths_m,
            output_units="photon_flux_per_m",
            concentration=1,
        )
        # TMM supplies the optical generation profile used by both the QE and
        # IV calculations. A 2 mV external voltage grid resolves the MPP.
        solver_options = {
            "T_ambient": temperature_k,
            "wavelength": wavelengths_m,
            "light_source": light_source,
            "optics_method": "TMM",
            "no_back_reflection": True,
            "position": 2e-9,
            "light_iv": True,
            "mpp": True,
            "voltages": np.linspace(0.0, 1.4, 701),
            "internal_voltages": np.linspace(-1.5, 1.5, 301),
        }

        # Independent objects keep the two solver tasks free of cached optical
        # or electrical state.
        cell_qe = _build_da_cell(
            temperature_k,
            surface_recombination_velocity_m_s,
            electron_minority_lifetime_s,
            hole_minority_lifetime_s,
        )
        cell_iv = _build_da_cell(
            temperature_k,
            surface_recombination_velocity_m_s,
            electron_minority_lifetime_s,
            hole_minority_lifetime_s,
        )
        solar_cell_solver(cell_iv, "iv", user_options=solver_options)
        solar_cell_solver(cell_qe, "qe", user_options=solver_options)

        # EQE is retained only long enough to convert incident photon energy
        # into exported electrical energy. It is not written to the lookup.
        eqe = np.asarray(cell_qe[0].eqe(wavelengths_m), dtype=float)
        vmpp_v = float(cell_iv.iv["Vmpp"])
        isc_a_m2 = float(cell_iv.iv["Isc"])
        voc_v = float(cell_iv.iv["Voc"])
        pmpp_w_m2 = float(cell_iv.iv["Pmpp"])

    # Numerical roundoff can produce values infinitesimally outside the
    # physical 0-to-1 interval, so normalize them before further use.
    eqe = np.clip(np.nan_to_num(eqe, nan=0.0), 0.0, 1.0)
    return TemperatureResult(
        temperature_c=float(temperature_c),
        vmpp_v=vmpp_v,
        isc_a_m2=isc_a_m2,
        voc_v=voc_v,
        pmpp_w_m2=pmpp_w_m2,
        eqe=eqe,
    )


# ---------------------------------------------------------------------------
# Far-infrared optical constants and the thick-slab optics that use them
# ---------------------------------------------------------------------------
def load_tabulated_nk(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a refractiveindex.info YAML file into (wavelength_um, n, k).

    Only the first "tabulated nk" block is used. That format is three
    whitespace-separated columns per line: wavelength in micrometres, n, k.
    """
    with open(path) as stream:
        document = yaml.safe_load(stream)

    # A single YAML file can hold several dispersion blocks of different kinds.
    for entry in document.get("DATA", []):
        if entry.get("type", "").strip() == "tabulated nk":
            rows = [line.split() for line in entry["data"].strip().splitlines()]
            table = np.array([[float(value) for value in row[:3]] for row in rows])
            return table[:, 0], table[:, 1], table[:, 2]

    raise ValueError(f"{Path(path).name} contains no 'tabulated nk' data block.")


def incoherent_slab_rat(
    n: np.ndarray,
    k: np.ndarray,
    wavelength_um: np.ndarray,
    thickness_um: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """R/A/T of a thick slab with two identical surfaces, ignoring interference.

    Beyond the band gap the wafer is hundreds of wavelengths thick and its real
    thickness varies by more than a wavelength across the die, so interference
    fringes are not physical and the incoherent sum of multiple internal passes
    is the right model. Summing that geometric series gives closed forms for R
    and T, and A is whatever is left, which makes R + A + T = 1 exact.
    """
    # Normal-incidence Fresnel reflectance of one vacuum/material interface.
    surface_r = ((n - 1.0) ** 2 + k ** 2) / ((n + 1.0) ** 2 + k ** 2)

    # Single-pass internal transmission through the slab (Beer-Lambert).
    internal_t = np.exp(-4.0 * np.pi * k * thickness_um / wavelength_um)

    # Geometric series over all internal round trips.
    denominator = 1.0 - (surface_r * internal_t) ** 2
    reflectivity = surface_r * (
        1.0 + (1.0 - surface_r) ** 2 * internal_t ** 2 / denominator
    )
    transmissivity = (1.0 - surface_r) ** 2 * internal_t / denominator

    # Clip the two computed quantities, then take absorption as the remainder so
    # that conservation holds exactly rather than to within rounding.
    reflectivity = np.clip(reflectivity, 0.0, 1.0)
    transmissivity = np.clip(transmissivity, 0.0, 1.0 - reflectivity)
    absorptivity = 1.0 - reflectivity - transmissivity
    return reflectivity, absorptivity, transmissivity


# ---------------------------------------------------------------------------
# Optical reflection, absorption, and transmission
# ---------------------------------------------------------------------------
def calculate_bare_cell_rat(
    wavelengths_nm: np.ndarray,
    wafer_thickness_um: float = default_wafer_thickness_um,
    temperature_c: float = reference_temperature_c,
    quiet: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Calculate R/A/T for a bare, freestanding GaAs wafer cell.

    The stack is Solcore's default GaAs junction followed by enough undoped GaAs
    substrate to give the requested total wafer thickness. The thick substrate
    is treated incoherently. No coverglass, adhesive, ARC, grid shadowing or
    rear metallisation is included.

    The three spectra are returned on the complete simulation grid. Solcore's
    transfer matrix covers its own validated interval; beyond it the spectra
    come from tabulated far-infrared n,k through an incoherent slab model. The
    returned mask marks the Solcore interval.
    """
    # Keep the complete simulation grid, even where optical data are missing.
    # A Boolean mask marks exactly where Solcore's bundled constants are valid.
    wavelengths_nm = np.asarray(wavelengths_nm, dtype=float)
    valid = (
        (wavelengths_nm >= solcore_optical_min_nm)
        & (wavelengths_nm <= solcore_optical_max_nm)
    )

    # Unsupported entries begin as NaN so that a failure to fill them later is
    # loud rather than silent. They are replaced by the edge-hold below before
    # this function returns.
    reflectivity = np.full(wavelengths_nm.shape, np.nan, dtype=float)
    absorptivity = np.full(wavelengths_nm.shape, np.nan, dtype=float)
    transmissivity = np.full(wavelengths_nm.shape, np.nan, dtype=float)
    if not valid.any():
        return reflectivity, absorptivity, transmissivity, valid

    with suppress_process_output(quiet):
        # Local imports keep Solcore's optional-solver initialization messages
        # inside this controlled section.
        from solcore import material
        from solcore.absorption_calculator import calculate_rat
        from solcore.solar_cell import default_GaAs
        from solcore.structure import Layer, Structure

        temperature_k = temperature_c + 273.15
        junction = default_GaAs(temperature_k)

        # The four active layers occupy 3.38 um. Add undoped GaAs behind them
        # until the requested total wafer thickness is reached.
        active_thickness_m = sum(float(layer.width) for layer in junction)
        substrate_thickness_m = wafer_thickness_um * 1e-6 - active_thickness_m
        if substrate_thickness_m <= 0:
            raise ValueError(
                "wafer_thickness_um must exceed the 3.38 um active-stack thickness."
            )

        substrate = Layer(
            width=substrate_thickness_m,
            material=material("GaAs")(T=temperature_k),
            role="substrate",
        )

        # The thin active layers retain interference phase. The thick wafer is
        # incoherent because its optical thickness and real thickness variation
        # make coherent fringes physically inappropriate.
        layers = list(junction) + [substrate]
        stack = Structure(layers)
        result = calculate_rat(
            stack,
            wavelength=wavelengths_nm[valid],
            angle=0,
            coherent=False,
            coherency_list=["c", "c", "c", "c", "i"],
            no_back_reflection=False,
        )

    # Solcore returns dimensionless incident-power fractions. Small numerical
    # excursions are clipped to the physical interval without extrapolation.
    reflectivity[valid] = np.clip(result["R"], 0.0, 1.0)
    absorptivity[valid] = np.clip(result["A"], 0.0, 1.0)
    transmissivity[valid] = np.clip(result["T"], 0.0, 1.0)

    # Beyond Solcore's range, calculate R/A/T from tabulated far-infrared n,k
    # instead of extrapolating. The wafer is optically thick and incoherent
    # there, so the slab formula above is the appropriate model and needs no
    # transfer matrix. Below Solcore's range the first calculated value is held:
    # that region carries a negligible share of both the incident spectrum and
    # the thermal emission, and no tabulated data covers it.
    first_index = int(np.flatnonzero(valid)[0])
    last_index = int(np.flatnonzero(valid)[-1])

    reflectivity[:first_index] = reflectivity[first_index]
    absorptivity[:first_index] = absorptivity[first_index]
    transmissivity[:first_index] = transmissivity[first_index]

    far_infrared = np.zeros(wavelengths_nm.shape, dtype=bool)
    far_infrared[last_index + 1:] = True
    if far_infrared.any():
        nk_wavelength_um, nk_n, nk_k = load_tabulated_nk(far_infrared_nk_file)
        query_um = wavelengths_nm[far_infrared] * 1e-3
        n_values = np.interp(query_um, nk_wavelength_um, nk_n)
        k_values = np.interp(query_um, nk_wavelength_um, nk_k)
        (
            reflectivity[far_infrared],
            absorptivity[far_infrared],
            transmissivity[far_infrared],
        ) = incoherent_slab_rat(n_values, k_values, query_um, wafer_thickness_um)

    return reflectivity, absorptivity, transmissivity, valid


# ---------------------------------------------------------------------------
# Temperature sweep orchestration
# ---------------------------------------------------------------------------
def _solve_temperatures(
    temperatures_c: np.ndarray,
    electrical_wavelengths_nm: np.ndarray,
    surface_recombination_velocity_m_s: float,
    electron_minority_lifetime_s: float,
    hole_minority_lifetime_s: float,
    jobs: int,
    quiet: bool,
) -> list[TemperatureResult]:
    """Solve each temperature serially or in independent worker processes."""
    # A single-worker path is useful for debugging and avoids multiprocessing
    # overhead when only one or two temperatures are requested.
    if jobs <= 1:
        results = []
        for index, temperature_c in enumerate(temperatures_c, start=1):
            result = _solve_one_temperature(
                float(temperature_c),
                electrical_wavelengths_nm,
                surface_recombination_velocity_m_s,
                electron_minority_lifetime_s,
                hole_minority_lifetime_s,
                quiet=quiet,
            )
            results.append(result)
            if not quiet:
                print(
                    f"[{index:02d}/{temperatures_c.size:02d}] "
                    f"T={temperature_c:6.1f} C, Vmp={result.vmpp_v:.4f} V",
                    flush=True,
                )
        return results

    # Each future owns a complete temperature solve. No mutable Solcore cell
    # objects are shared between processes.
    futures = {}
    results = []
    with ProcessPoolExecutor(max_workers=jobs) as executor:
        for temperature_c in temperatures_c:
            future = executor.submit(
                _solve_one_temperature,
                float(temperature_c),
                electrical_wavelengths_nm,
                surface_recombination_velocity_m_s,
                electron_minority_lifetime_s,
                hole_minority_lifetime_s,
                quiet,
            )
            futures[future] = float(temperature_c)

        # Futures finish out of order. Progress is optional and the final list
        # is sorted below before it is inserted into the lookup table.
        completed = 0
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            completed += 1
            if not quiet:
                print(
                    f"[{completed:02d}/{temperatures_c.size:02d}] "
                    f"T={result.temperature_c:6.1f} C, Vmp={result.vmpp_v:.4f} V",
                    flush=True,
                )

    return sorted(results, key=lambda item: item.temperature_c)


# ---------------------------------------------------------------------------
# Assemble the full temperature-by-wavelength lookup
# ---------------------------------------------------------------------------
def build_lookup(
    simulation_pickle: str | Path,
    minimum_temperature_c: float = -50.0,
    maximum_temperature_c: float = 100.0,
    temperature_step_c: float = 5.0,
    wafer_thickness_um: float = default_wafer_thickness_um,
    surface_recombination_velocity_m_s: float = (
        default_surface_recombination_velocity_m_s
    ),
    electron_minority_lifetime_s: float = default_electron_minority_lifetime_s,
    hole_minority_lifetime_s: float = default_hole_minority_lifetime_s,
    jobs: int = 1,
    quiet: bool = True,
) -> dict[str, Any]:
    """Build, validate, and return the complete lookup dictionary."""
    # The wavelength axis is copied exactly from the simulation input so no
    # downstream interpolation is needed merely to align array columns.
    wavelengths_nm = load_simulation_wavelengths(simulation_pickle)
    temperatures_c = make_temperature_grid(
        minimum_temperature_c,
        maximum_temperature_c,
        temperature_step_c,
    )

    # Restrict expensive QE/IV optics to wavelengths for which the bundled
    # GaAs material data are supported and an electrical response could occur.
    electrical_valid = (
        (wavelengths_nm >= solcore_optical_min_nm)
        & (wavelengths_nm <= electrical_solve_max_nm)
    )
    electrical_wavelengths_nm = wavelengths_nm[electrical_valid]
    if electrical_wavelengths_nm.size == 0:
        raise ValueError("No simulation wavelengths overlap the GaAs electrical range.")

    # The optical constants are not temperature-resolved, so calculate one
    # reference R/A/T spectrum and repeat it later across all temperatures.
    reflectivity_1d, total_absorptivity_1d, transmissivity_1d, optical_valid = (
        calculate_bare_cell_rat(
            wavelengths_nm,
            wafer_thickness_um=wafer_thickness_um,
            temperature_c=reference_temperature_c,
            quiet=quiet,
        )
    )

    # Check optical power conservation before using the spectrum in the
    # electrical/thermal partition. The edge-hold populates every column, so
    # this now covers the whole grid rather than only the calculated interval.
    rat_sum = (
        reflectivity_1d
        + total_absorptivity_1d
        + transmissivity_1d
    )
    if not np.allclose(rat_sum, 1.0, rtol=0, atol=5e-10):
        raise RuntimeError("Optical calculation failed the R + A + T = 1 check.")

    # Obtain EQE and the AM0 maximum-power operating point at each temperature.
    results = _solve_temperatures(
        temperatures_c,
        electrical_wavelengths_nm,
        surface_recombination_velocity_m_s,
        electron_minority_lifetime_s,
        hole_minority_lifetime_s,
        jobs=max(1, int(jobs)),
        quiet=quiet,
    )

    # Allocate complete two-dimensional arrays. EQE remains an internal work
    # array and is intentionally omitted from the returned table.
    n_temperature = temperatures_c.size
    n_wavelength = wavelengths_nm.size
    eqe_internal = np.zeros((n_temperature, n_wavelength), dtype=float)
    vmpp_v = np.empty(n_temperature, dtype=float)
    isc_a_m2 = np.empty(n_temperature, dtype=float)
    voc_v = np.empty(n_temperature, dtype=float)
    pmpp_w_m2 = np.empty(n_temperature, dtype=float)

    # Results have already been sorted into ascending temperature order.
    for row, result in enumerate(results):
        eqe_internal[row, electrical_valid] = result.eqe
        vmpp_v[row] = result.vmpp_v
        isc_a_m2[row] = result.isc_a_m2
        voc_v[row] = result.voc_v
        pmpp_w_m2[row] = result.pmpp_w_m2

    # Convert wavelength to photon energy in eV. Numerically, q*V joules per
    # collected electron divided by photon energy is equivalent to V/E[eV].
    photon_energy_ev = h * c / (wavelengths_nm * 1e-9) / e

    # This is the fraction of incident monochromatic power exported at the AM0
    # MPP, not the electron-per-photon EQE itself.
    electrical_absorptivity = eqe_internal * (
        vmpp_v[:, np.newaxis] / photon_energy_ev[np.newaxis, :]
    )

    # The electrical fraction is explicitly zero where no electrical solution
    # was performed. The accompanying mask tells consumers which region was
    # actually evaluated.
    electrical_absorptivity[:, ~electrical_valid] = 0.0
    electrical_absorptivity = np.clip(electrical_absorptivity, 0.0, 1.0)

    # Repeat the reference optical spectra along the temperature axis. copy()
    # creates writable arrays instead of read-only broadcast views.
    total_absorptivity = np.broadcast_to(
        total_absorptivity_1d, (n_temperature, n_wavelength)
    ).copy()
    reflectivity = np.broadcast_to(
        reflectivity_1d, (n_temperature, n_wavelength)
    ).copy()
    transmissivity = np.broadcast_to(
        transmissivity_1d, (n_temperature, n_wavelength)
    ).copy()

    # Prevent any tiny mismatch between the independent QE and TMM calculations
    # from assigning more electrical power than the cell optically absorbs.
    electrical_absorptivity[:, optical_valid] = np.minimum(
        electrical_absorptivity[:, optical_valid],
        total_absorptivity[:, optical_valid],
    )
    # Under the stated energy-partition assumption, everything absorbed but not
    # exported electrically becomes heat inside the cell.
    thermal_absorptivity = total_absorptivity - electrical_absorptivity

    # Store enough context for another program to interpret the numbers without
    # relying on this source file or the README.
    metadata = {
        "schema_version": 1,
        "model": "Solcore 5.10.1 default_GaAs DA plus TMM",
        "illumination_for_operating_point": "AM0, one sun",
        "incidence": "normal incidence from vacuum",
        "temperature_unit": "degC",
        "wavelength_unit": "nm",
        "fractions_are_dimensionless": True,
        "temperature_range_C": [
            float(temperatures_c[0]),
            float(temperatures_c[-1]),
        ],
        "temperature_step_C": float(temperature_step_c),
        "source_wavelength_count": int(wavelengths_nm.size),
        "source_wavelength_range_nm": [
            float(wavelengths_nm[0]),
            float(wavelengths_nm[-1]),
        ],
        "validated_optical_range_nm": [
            solcore_optical_min_nm,
            solcore_optical_max_nm,
        ],
        "electrical_solve_range_nm": [
            float(electrical_wavelengths_nm[0]),
            float(electrical_wavelengths_nm[-1]),
        ],
        "wafer_thickness_um": float(wafer_thickness_um),
        "effective_surface_recombination_velocity_m_s": float(
            surface_recombination_velocity_m_s
        ),
        "electron_minority_lifetime_s": float(electron_minority_lifetime_s),
        "hole_minority_lifetime_s": float(hole_minority_lifetime_s),
        "optical_stack": [
            "30 nm Al0.8Ga0.2As window",
            "150 nm p-GaAs emitter",
            "3000 nm n-GaAs base",
            "200 nm n-GaAs BSF",
            "GaAs substrate to requested total wafer thickness",
        ],
        "excluded_hardware": [
            "coverglass",
            "adhesive",
            "anti-reflection coating",
            "front grid shadowing",
            "rear metallisation",
            "panel substrate/backsheet",
        ],
        "electrical_absorptivity_definition": (
            "EQE(lambda,T) * q*Vmpp(T) / (h*c/lambda); EQE is internal only"
        ),
        "thermal_absorptivity_definition": (
            "total_absorptivity - electrical_absorptivity_mpp"
        ),
        "thermal_partition_assumption": (
            "All absorbed power not exported at MPP becomes heat; luminescent "
            "escape is neglected."
        ),
        "optical_temperature_assumption": (
            "R/A/T are held fixed at 25 C because Solcore's bundled GaAs and "
            "AlGaAs n,k tables are not temperature-resolved."
        ),
        "electrical_model_note": (
            "The DA solver is used because the bundled PDD solver failed a "
            "continuous-temperature convergence check below room temperature. "
            "Diffusion lengths are recalculated from mobility, temperature and "
            "the stated minority-carrier lifetimes."
        ),
        "optics_beyond_solcore_range": (
            "Above the Solcore interval R/A/T come from tabulated n,k through "
            "an incoherent thick-slab model; below it the first Solcore value "
            "is held, a region carrying a negligible share of both the incident "
            "spectrum and the thermal emission. valid_optics_mask marks the "
            "Solcore interval. Note that undoped GaAs is genuinely transparent "
            "from about 3 to 15 um, so absorptivity, and by Kirchhoff "
            "emissivity, is near zero across the peak of a 300 K blackbody."
        ),
        "far_infrared_nk_source": (
            "Franta et al. 2026, undoped GaAs n,k 0.895-400 um at 300 K, "
            "refractiveindex.info main/GaAs/nk/Franta-300K.yml"
        ),
    }

    # External field names are deliberately preserved for compatibility with
    # the previously supplied lookup. They include conventional capital unit
    # symbols such as C and V, but every Python variable remains lower case.
    table = {
        "temperature_C": temperatures_c,
        "wavelength_nm": wavelengths_nm,
        "total_absorptivity": total_absorptivity,
        "electrical_absorptivity_mpp": electrical_absorptivity,
        "thermal_absorptivity_mpp": thermal_absorptivity,
        "reflectivity": reflectivity,
        "transmissivity": transmissivity,
        "valid_optics_mask": optical_valid,
        "valid_electrical_mask": electrical_valid,
        "vmpp_V": vmpp_v,
        "isc_A_m2": isc_a_m2,
        "voc_V": voc_v,
        "pmpp_W_m2": pmpp_w_m2,
        "metadata": metadata,
    }
    # Validation occurs before the caller is allowed to save or interpolate the
    # lookup, so malformed or non-conservative data cannot silently propagate.
    validate_lookup(table)
    return table


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def validate_lookup(table: dict[str, Any]) -> None:
    """Raise if the lookup violates dimensional or energy-balance checks."""
    # Recover the two independent coordinates and therefore the required shape
    # of every spectral property array.
    temperatures_c = np.asarray(table["temperature_C"])
    wavelengths_nm = np.asarray(table["wavelength_nm"])
    shape = (temperatures_c.size, wavelengths_nm.size)

    # Every optical column is now populated: calculated inside Solcore's range
    # and held at the nearest calculated value outside it. The checks below
    # therefore cover the whole grid instead of only the calculated interval,
    # which is what catches a failure of the edge-hold to fill its columns.
    valid = np.ones(wavelengths_nm.size, dtype=bool)

    # All five output spectra must share [temperature, wavelength] indexing.
    for name in (
        "total_absorptivity",
        "electrical_absorptivity_mpp",
        "thermal_absorptivity_mpp",
        "reflectivity",
        "transmissivity",
    ):
        if np.asarray(table[name]).shape != shape:
            raise ValueError(f"{name} does not have shape {shape}.")

    # Optical fractions must be finite and physical at every wavelength. A NaN
    # surviving here means the edge-hold did not populate a column.
    for name in ("total_absorptivity", "reflectivity", "transmissivity"):
        values = np.asarray(table[name])[:, valid]
        if not np.isfinite(values).all() or np.any((values < 0) | (values > 1)):
            raise ValueError(f"{name} is invalid within the supported optical range.")

    # Verify the requested split of absorbed power into useful electricity and
    # residual heat.
    electrical = np.asarray(table["electrical_absorptivity_mpp"])[:, valid]
    thermal = np.asarray(table["thermal_absorptivity_mpp"])[:, valid]
    total = np.asarray(table["total_absorptivity"])[:, valid]
    if not np.isfinite(electrical).all() or not np.isfinite(thermal).all():
        raise ValueError("Electrical/thermal partition contains non-finite values.")
    if np.any(electrical < -1e-12) or np.any(thermal < -1e-12):
        raise ValueError("Electrical/thermal absorptivity became negative.")
    if not np.allclose(electrical + thermal, total, rtol=0, atol=2e-12):
        raise ValueError("Electrical plus thermal absorptivity does not equal total.")

    # Finally verify conservation of incident optical power.
    rat = (
        np.asarray(table["reflectivity"])[:, valid]
        + total
        + np.asarray(table["transmissivity"])[:, valid]
    )
    if not np.allclose(rat, 1.0, rtol=0, atol=5e-10):
        raise ValueError("R + A + T does not equal one.")


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------
def save_lookup(table: dict[str, Any], output_base: str | Path) -> tuple[Path, Path]:
    """Save Python-pickle lookup files."""
    # Treat output_base as a stem and create its parent directory if needed.
    output_base = Path(output_base)
    output_base.parent.mkdir(parents=True, exist_ok=True)
    pickle_path = output_base.with_suffix(".pkl")

    # The pickle contains the same dictionary structure and is supplied for
    # direct compatibility with the existing simulation workflow. Only load a
    # pickle when its source is trusted.
    with pickle_path.open("wb") as stream:
        pickle.dump(table, stream, protocol=pickle.HIGHEST_PROTOCOL)
    return pickle_path

# ---------------------------------------------------------------------------
# Optional lookup interpolation helper
# ---------------------------------------------------------------------------
def interpolate_properties(
    table: dict[str, Any],
    temperature_c: float | np.ndarray,
    wavelength_nm: float | np.ndarray,
    properties: tuple[str, ...],
    pickle_path: Path = (
        "total_absorptivity",
        "electrical_absorptivity_mpp",
        "thermal_absorptivity_mpp",
        "reflectivity",
        "transmissivity",
    ),
) -> dict[str, np.ndarray]:
    """Linearly interpolate lookup properties at paired T/lambda query points.

    NumPy broadcasting is applied to temperature_c and wavelength_nm. Optical
    queries outside Solcore's validated range return NaN by construction.
    """
    # Broadcasting permits scalar/scalar, scalar/vector, or paired array input.
    temperature_query, wavelength_query = np.broadcast_arrays(
        electrical_absorptivity_mpp=table["electrical_absorptivity_mpp"],
        thermal_absorptivity_mpp=table["thermal_absorptivity_mpp"],
        reflectivity=table["reflectivity"],
        transmissivity=table["transmissivity"],
        valid_optics_mask=table["valid_optics_mask"],
        valid_electrical_mask=table["valid_electrical_mask"],
        vmpp_V=table["vmpp_V"],
        isc_A_m2=table["isc_A_m2"],
        voc_V=table["voc_V"],
        pmpp_W_m2=table["pmpp_W_m2"],
        metadata_json=np.asarray(json.dumps(table["metadata"], indent=2)),
    )
    # The pickle contains the same dictionary structure and is supplied for
    # direct compatibility with the existing simulation workflow. Only load a
    # pickle when its source is trusted.
    with pickle_path.open("wb") as stream:
        pickle.dump(table, stream, protocol=pickle.HIGHEST_PROTOCOL)
    return pickle_path

# ---------------------------------------------------------------------------
# Optional lookup interpolation helper
# ---------------------------------------------------------------------------
def interpolate_properties(
    table: dict[str, Any],
    temperature_c: float | np.ndarray,
    wavelength_nm: float | np.ndarray,
    properties: tuple[str, ...] = (
        "total_absorptivity",
        "electrical_absorptivity_mpp",
        "thermal_absorptivity_mpp",
        "reflectivity",
        "transmissivity",
    ),
) -> dict[str, np.ndarray]:
    """Linearly interpolate lookup properties at paired T/lambda query points.

    NumPy broadcasting is applied to temperature_c and wavelength_nm. Optical
    queries outside Solcore's validated range return NaN by construction.
    """
    # Broadcasting permits scalar/scalar, scalar/vector, or paired array input.
    temperature_query, wavelength_query = np.broadcast_arrays(
        np.asarray(temperature_c, dtype=float),
        np.asarray(wavelength_nm, dtype=float),
    )
    # RegularGridInterpolator expects one [temperature, wavelength] pair per
    # row, so flatten broadcast inputs and restore their shape after evaluation.
    points = np.column_stack(
        (temperature_query.ravel(), wavelength_query.ravel())
    )
    output = {}
    for property_name in properties:
        # A separate interpolator is constructed for each requested property.
        # NaNs in unsupported optical columns naturally remain unavailable.
        interpolator = RegularGridInterpolator(
            (table["temperature_C"], table["wavelength_nm"]),
            table[property_name],
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        output[property_name] = interpolator(points).reshape(
            temperature_query.shape
        )
    return output


# ---------------------------------------------------------------------------
# Lookup visualization
# ---------------------------------------------------------------------------
def plot_lookup(
    table: dict[str, Any],
    output_path: str | Path,
    show_plot: bool = False,
    dpi: int = 200,
) -> Path:
    """Plot the spectral properties and electrical temperature dependence.

    The upper-left panel shows the optical R/A/T balance. The upper-right
    panel partitions total absorption into exported electrical power and heat
    at the available temperature nearest 25 C. The lower panel shows the
    complete electrical-absorptivity lookup as a temperature/wavelength map.

    The whole simulation grid is plotted on a logarithmic wavelength axis, and
    the calculated interval is shaded. Everything outside that shading is the
    held edge value, not a calculated one.
    """
    # Import matplotlib only when plotting is requested. This keeps the core
    # lookup functions importable in workflows that do not create figures.
    import matplotlib.pyplot as plt

    # Validate again so this function is also safe when called independently
    # with a lookup loaded from disk rather than directly from build_lookup().
    validate_lookup(table)

    temperatures_c = np.asarray(table["temperature_C"], dtype=float)
    wavelengths_nm = np.asarray(table["wavelength_nm"], dtype=float)
    valid_optics = np.asarray(table["valid_optics_mask"], dtype=bool)
    valid_electrical = np.asarray(
        table["valid_electrical_mask"], dtype=bool
    )

    # The optical constants are evaluated at 25 C. If a custom temperature
    # grid does not contain 25 C exactly, use and label its nearest row.
    reference_index = int(
        np.argmin(np.abs(temperatures_c - reference_temperature_c))
    )
    selected_temperature_c = temperatures_c[reference_index]

    total_absorptivity = np.asarray(table["total_absorptivity"], dtype=float)
    electrical_absorptivity = np.asarray(
        table["electrical_absorptivity_mpp"], dtype=float
    )
    thermal_absorptivity = np.asarray(
        table["thermal_absorptivity_mpp"], dtype=float
    )
    reflectivity = np.asarray(table["reflectivity"], dtype=float)
    transmissivity = np.asarray(table["transmissivity"], dtype=float)

    # A two-row layout leaves enough horizontal space for readable wavelength
    # axes while allowing the heatmap to span the full width of the figure.
    figure = plt.figure(figsize=(13.0, 9.0), constrained_layout=True)
    grid = figure.add_gridspec(2, 1, height_ratios=(1.0, 1.15))
    optical_axis = figure.add_subplot(grid[0, 0])
    partition_axis = figure.add_subplot(grid[1, 0])

    # Plot the complete grid. The shading added below marks which part of it
    # Solcore actually calculated, so the held region is never mistaken for it.
    optical_wavelengths_nm = wavelengths_nm
    calculated_min_nm = float(wavelengths_nm[valid_optics][0])
    calculated_max_nm = float(wavelengths_nm[valid_optics][-1])

    # Total absorptivity, reflectivity and transmissivity are the three direct
    # outputs of the transfer-matrix calculation and must sum to one.
    optical_axis.plot(
        optical_wavelengths_nm,
        total_absorptivity[reference_index, :],
        label="Total absorptivity",
        linewidth=1.8,
    )
    optical_axis.plot(
        optical_wavelengths_nm,
        reflectivity[reference_index, :],
        label="Reflectivity",
        linewidth=1.8,
    )
    optical_axis.plot(
        optical_wavelengths_nm,
        transmissivity[reference_index, :],
        label="Transmissivity",
        linewidth=1.8,
    )
    optical_axis.set_title("Optical spectral properties")
    optical_axis.set_xlabel("Wavelength (nm)")
    optical_axis.set_ylabel("Fraction of incident power")
    optical_axis.set_ylim(0.0, 1.05)
    # A log axis is needed because the simulation grid spans 15 nm to 50 um.
    optical_axis.set_xscale("log")
    optical_axis.axvspan(
        calculated_min_nm,
        calculated_max_nm,
        color="0.85",
        zorder=0,
        label="Calculated by Solcore (elsewhere: edge value held)",
    )
    optical_axis.grid(alpha=0.25)
    optical_axis.legend()

    # This panel makes the energy accounting explicit: electrical and thermal
    # absorptivity add to total absorptivity at every supported wavelength.
    partition_axis.plot(
        optical_wavelengths_nm,
        total_absorptivity[reference_index, :],
        label="Total absorptivity",
        linewidth=1.8,
    )
    partition_axis.plot(
        optical_wavelengths_nm,
        electrical_absorptivity[reference_index, :],
        label="Electrical absorptivity at MPP",
        linewidth=1.8,
    )
    partition_axis.plot(
        optical_wavelengths_nm,
        thermal_absorptivity[reference_index, :],
        label="Thermal absorptivity at MPP",
        linewidth=1.8,
    )
    partition_axis.set_title(
        f"Absorbed-power partition at {selected_temperature_c:g} C"
    )
    partition_axis.set_xlabel("Wavelength (nm)")
    partition_axis.set_ylabel("Fraction of incident power")
    partition_axis.set_ylim(0.0, 1.05)
    # A log axis is needed because the simulation grid spans 15 nm to 50 um.
    partition_axis.set_xscale("log")
    partition_axis.axvspan(
        calculated_min_nm,
        calculated_max_nm,
        color="0.85",
        zorder=0,
        label="Calculated by Solcore (elsewhere: edge value held)",
    )
    partition_axis.grid(alpha=0.25)
    partition_axis.legend()

    # Saving happens before display, so closing an interactive window cannot
    # prevent the PNG from being written. Non-interactive runs close the
    # figure immediately to release memory.
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=dpi)


    # plot Vmpp, Isc, Voc, and Pmpp vs temperature in a separate figure
    electrical_axis = plt.figure(figsize=(10.0, 6.0), constrained_layout=True).add_subplot(1, 1, 1)
    electrical_axis.plot(
        temperatures_c,
        table["vmpp_V"],
        label="Vmpp (V)",
        linewidth=1.8,
    )
    electrical_axis.set_title("Electrical properties vs temperature")
    electrical_axis.set_xlabel("Temperature (C)")
    electrical_axis.set_ylabel("Electrical properties")
    electrical_axis.grid(alpha=0.25)
    electrical_axis.legend()

    if show_plot:
        plt.show()
    plt.close(figure)


    return output_path


# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------
def main() -> None:
    """Generate and save the lookup using the configuration near the top."""
    # Give a direct path-specific error rather than allowing pandas to fail with
    # a less informative traceback later.
    if not pickle_file.is_file():
        raise FileNotFoundError(
            f"Could not find the simulation pickle at: {pickle_file}\n"
            "Place simulation_data.pkl in the 'Simulation Data' folder beside "
            "this script, or edit pickle_file in the configuration block."
        )

    # Build performs all optical/electrical calculations and validates the
    # resulting table before this function writes anything.
    table = build_lookup(
        simulation_pickle=pickle_file,
        minimum_temperature_c=temperature_min_c,
        maximum_temperature_c=temperature_max_c,
        temperature_step_c=temperature_step_c,
        wafer_thickness_um=wafer_thickness_um,
        surface_recombination_velocity_m_s=(
            surface_recombination_velocity_m_s
        ),
        electron_minority_lifetime_s=electron_minority_lifetime_s,
        hole_minority_lifetime_s=hole_minority_lifetime_s,
        jobs=jobs,
        quiet=not verbose_solcore,
    )

    # Save both formats only after the in-memory table has passed validation.
    pickle_path = save_lookup(table, output_base)

    # Plot only after the validated lookup has been saved. The figure contains
    # the optical spectra, electrical/thermal split, and full temperature map.
    saved_plot_path = plot_lookup(
        table,
        output_path=plot_file,
        show_plot=show_plot,
        dpi=plot_dpi,
    )
    print(f"Saved {pickle_path}")
    print(f"Saved {saved_plot_path}")


if __name__ == "__main__":
    main()