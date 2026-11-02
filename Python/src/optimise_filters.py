"""
optimise_filters.py
===================

Spectral-filter design for the rotating-satellite simulation, solved with CasADi + IPOPT.

The objective balances electrical energy against temperature tracking, with
separate weights. Tracking penalises squared deviations above or below t_set.
The unfiltered case is retained for comparison and initial guesses, using the
same incident power, initial temperatures and thermal model. Its electrical
power is not a constraint. Temperature also affects electrical conversion
through the lookup. Hard operating temperature bounds apply to the filtered
design. The unfiltered reference can exceed those bounds, but must stay at
positive temperature and within the electrical lookup's actual data range.

The script is laid out so that every block maps one-to-one onto a numbered equation in ``Solar Panel Math''.  Summations in the maths are written as explicit ``for'' loops so the code can be read side-by-side with the LaTeX and edited in the same places.

    Section 1   run configuration and model inputs  (edit these)
    Section 2   panel spectral properties (silicon and aluminium lookups)
    Section 3   data loading + wavelength-grid reduction
    Section 4   blackbody spectral emissive power M_lambda(lambda, T)
    Section 5   the NLP, equation by equation
    Section 6   initial guesses (nominal + randomised, for multi-start)
    Section 7   solve / multi-start driver
    Section 8   results

There are no command-line arguments: everything is set in Section 1.0 and the script is run as-is.
Set solve_new_filter = False to plot output_file without loading simulation or
lookup data, building the NLP or running the solver.

Decision variables (per filter type j, per wavelength bin k):
    alpha_j(lambda), rho_j(lambda), tau_j(lambda)
    eps_j(lambda) = alpha_j(lambda)                       (Kirchhoff, eliminated)
Auxiliary variables:
    T_i(t)                                                (panel temperatures)

Equation map
------------
    eq (1)  heat flows / spectral sum          ->  Section 5.4
    eq (2)  dQ_i = Qdot_i dt                   ->  Section 5.5
    eq (3)  T_i(t+dt) = T_i(t) + dQ_i/(mc)     ->  Section 5.5
    eq (4)  cost J                            ->  Section 5.6
    eq (5)  smoothness  |dp/dlambda| <= D      ->  Section 5.3
    eq (6)  alpha + rho + tau = 1              ->  Section 5.2
    Kirchhoff  eps_j = alpha_j                 ->  Section 5.1 (by elimination)

The filter and panel share one temperature. The filter has the same optical
properties on both sides. Incoherent reflections between them are summed.
External thermal emission uses the total absorptivity of this combined stack,
including absorption in the filter. Internal radiation does not remove energy
from the combined thermal node, but emission to space does.
Radiation transmitted through the silicon wafer is absorbed by an opaque
backing at the same temperature. It cannot escape through the rear of the panel.

Remaining approximations
    - Normal-incidence optical properties stand in for angular averages.
    - Only the outward face radiates. The rear of the panel is closed.
    - The backing absorbs all wafer transmission without returning light to the
      wafer. Its heat capacity is included in the chosen panel mass and heat capacity.
    - The bare-wafer lookup does not include coverglass or contacts.
    - Electrical conversion is the supplied one-sun AM0 MPP allocation. A new
      MPP is not calculated for the filtered illumination or illumination level.
    - Panels have no thermal coupling or filter heat capacity. A prescribed
      constant internal heat load can be distributed to selected panel types.
    - Radiation outside the simulation wavelength range is omitted.
    - The result covers the supplied transient, not a periodic orbital equilibrium.
    - Free spectral coefficients and smoothness bounds do not specify a
      manufacturable coating. This is an idealised spectral design problem.

The comments in Section 5.4 give the stack equations. The physical basis is
Kirchhoff's law for reciprocal isothermal structures and an incoherent sum of
repeated reflections. See https://arxiv.org/abs/2306.06151 and
https://docs.solcore.solar/en/latest/Optics/tmm.html .
"""

from __future__ import annotations

import pickle
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import casadi as ca
import numpy as np
import pandas as pd

from constants import *
from aluminium_initialise import load_aluminium_lookup

_here = Path(__file__).resolve().parent

# Set plotting to be with latex formatting
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

# =============================================================================
# 1.  RUN CONFIGURATION AND MODEL INPUTS
#     -------------------------------------------------------------------
#     The pickle supplies:  Qdot_in(t, lambda) per panel [W per bin], the
#     wavelength and time grids, and each panel's material type and area.
#     Shared physical constants are read from constants.py.  The remaining
#     values in this section are configuration choices or model inputs.
#     Values flagged  <<CHECK>>  are estimates, not measured data.
# =============================================================================

# ---- 1.0  Run configuration --------------------------------------------------
# All paths are relative to this script.
simulation_pickle_file = _here / "Simulation Data" / "simulation_data.pkl"
lookup_file = _here / "Solar Panel Properties" / "silicon_single_junction_spectral_lookup.pkl"
aluminium_lookup_file = _here / "Aluminium Panel Properties" / "aluminium_spectral_lookup.pkl"
output_file = _here / "Optimisation Results" / "filter_solution.pkl"  # saved by solving, read by plot-only mode
plot_file = _here / "Optimisation Results" / "filter_solution.png"

solve_new_filter = True         # True: optimise; False: plot the existing output_file
n_starts = 1                    # multi-start points; 1 = single solve
random_seed = 0                 # seed for the randomised starts
save_solution = True            # write output_file
make_plot = True                # plot after solving; plot-only mode always plots
verbose_ipopt = False           # True prints the full IPOPT iteration log

# ---- 1.2  Material type labels (as used by the C++ sim) ----------------------
solar = 0                       # M_i = 0  ->  solar panel
aluminium = 1                   # M_i = 1  ->  aluminium

# NOTE ON AREA:  the pickle's panel powers are already absolute powers [W]
# (i.e. area-integrated), so do not multiply q_in by A_i again. Area enters
# emission and, optionally, the internal heat distribution. If the sim ever
# sends fluxes [W m^-2] instead, multiply q_in by A_i in Section 3.

# ---- 1.4  Per-material thermal constants -------------------------------------
# Per panel, including the backing: 200.6 um Si on 1 mm Al for each
# 0.02 m^2 PV face; each 0.01 m^2 passive face is a 1 mm Al plate.
# rho_Si = 2330 kg/m^3, cp_Si = 714 J/(kg K); rho_Al = 2707 kg/m^3,
# cp_Al = 896 J/(kg K). Effective PV cp is the mass-weighted value.
panel_mass = {solar: 0.06348796, aluminium: 0.02707}          # m  [kg]
panel_heat_capacity = {solar: 869.20234, aluminium: 896.0}    # c  [J/kg/K]
panel_t0 = {solar: 293.15, aluminium: 293.15}                # T_i(0) [K], prescribed

# Per-panel overrides, keyed by panel label; keys: mass, heat_capacity, T0.
panel_overrides: dict[str, dict[str, float]] = {
    # "+x": {"mass": 0.15, "T0": 300.0},
}

# Constant TOTAL thermal power generated by internal components, not per panel.
# The full load is shared only among enabled panel types present in the data.
# Both type toggles False (or no matching panels) intentionally applies zero heat.
internal_heat_load_w = 10.0       # total internal heat load [W]        <<CHECK>>
internal_heat_to_metal = True   # aluminium panels receive internal heat
internal_heat_to_solar = True   # solar panels receive internal heat
internal_heat_by_area = True    # True: area weighted; False: equal per enabled panel

# ---- 1.5  Where each material's PANEL spectral properties come from ----------
# The material UNDER the filter.  Fixed data, not decision variables.
#   "lookup"     -> the temperature/wavelength table in lookup_file (1.6)
#   "aluminium_lookup" -> the opaque passive spectrum from aluminium_initialise.py
#   "parametric" -> the two-plateau model in panel_spectral below
panel_spectral_source = {solar: "lookup", aluminium: "aluminium_lookup"}
lookup_material = solar            # which material type the lookup describes

# Used only where the source is "parametric".  Kirchhoff on the panel as well,
# so eps_{lambda,i} = alpha_{lambda,i} and one curve serves both.  It is a
# smooth two-plateau blend in log(lambda):
#     alpha(lambda) = alpha_solar  ->  alpha_IR,  crossing over at lam_c_nm
panel_spectral = {
    solar: {                       # only used if the source is switched back
        "alpha_solar": 0.92, "alpha_IR": 0.85,
        "lam_c_nm": 3000.0, "blend_decades": 0.25,
    },
    aluminium: {                   # legacy assumed curve, disabled by default
        "alpha_solar": 0.15,       # absorptivity shortward of the crossover
        "alpha_IR": 0.05,          # absorptivity longward of the crossover
        "lam_c_nm": 3000.0,        # crossover wavelength                   [nm]
        "blend_decades": 0.25,     # crossover width, in decades of lambda
    },
}

# ---- 1.6  Spectral lookup table ----------------------------------------------
# lookup_file supplies, on the simulation's own wavelength grid:
#     total_absorptivity(lambda)           -> absorption in the wafer
#     transmissivity(lambda)              -> absorption in the opaque backing
#     electrical_absorptivity_mpp(T,lambda)-> alpha_{lambda,elec,i}(T)
# The backed panel has alpha_panel = alpha_wafer + tau_wafer = 1 - rho_wafer.
# Its emissivity follows this combined absorption by Kirchhoff. Electrical
# conversion still comes only from the wafer, so backing absorption becomes heat.
# Only the electrical part is temperature dependent; the generator holds R/A/T
# fixed at 25 C because its n,k tables are not temperature resolved.
#
# T interpolation: "bspline" is C2-continuous and passes through every knot,
# which suits an interior-point solver; "linear" reproduces the generator's own
# RegularGridInterpolator exactly but kinks at each knot.
t_interp_method = "bspline"        # "bspline" or "linear"
clamp_t_to_lookup = True           # override t_min/t_max with the lookup's range

# ---- 1.7  Design targets and cost weights ------------------------------------
t_set = 273.15                  # temperature set-point and emission band reference [K]
t_ref_scale = 10.0              # dT_ref, normalises the tracking term [K]
w_elec = 100.0                    # weight on electrical energy (maximised)
w_temp = 1.0                    # weight on temperature tracking (minimised)
w_smooth = 1.0e-3               # tie-breaker on filter roughness (see 5.6)

t_min, t_max = 0.0, 450.0     # hard bounds on T_i(t)                [K]

# ---- 1.8  Smoothness limits --------------------------------------------------
# eq. (5):   |dp_j / dlambda|      <= D_{j,p}     units: property per nm.
#            D = 1e-2 lets a property swing 0 -> 1 over 100 nm.
d_smooth = {
    solar:     {"alpha": 1.0e-2, "rho": 1.0e-2, "tau": 1.0e-2},
    aluminium: {"alpha": 1.0e-2, "rho": 1.0e-2, "tau": 1.0e-2},
}

# eq. (5'):  |dp_j / dln(lambda)|  <= D'_{j,p}    dimensionless.
# A per-nm limit is ~1e5 times looser at 50 um than at 500 nm, so on its own
# eq. (5) leaves the thermal IR unconstrained and the solution chatters
# bin-to-bin there. A real interference filter
# transitions over a fixed *fraction* of lambda, which is what this limits.
# D' = 5 allows a 0 -> 1 swing over roughly a factor of 1.2 in wavelength.
# Set the whole dict, or any single property, to None to use eq. (5) alone.
d_smooth_log = {
    solar:     {"alpha": 5.0, "rho": 5.0, "tau": 5.0},
    aluminium: {"alpha": 5.0, "rho": 5.0, "tau": 5.0},
}

# ---- 1.9  Wavelength grid reduction ------------------------------------------
# The sim ships 4999 bins.  A manufacturable filter cannot have 4999
# independent values, and the NLP is far better conditioned on a coarser grid.
# n_filter_bands = None  ->  use the native grid unchanged (slow: ~30k vars).
# At 120 log-spaced bands the emission quadrature error is ~2.5e-4.
n_filter_bands = 150
band_spacing = "log"            # "log" or "linear"

# ---- 1.10  Solver options ----------------------------------------------------
ipopt_options = {
    "ipopt.max_iter": 3000,
    "ipopt.tol": 1e-8,
    "ipopt.acceptable_tol": 1e-6,
    "ipopt.constr_viol_tol": 1e-8,
    "ipopt.acceptable_constr_viol_tol": 1e-8,
    "ipopt.mu_strategy": "adaptive",
    "ipopt.hessian_approximation": "exact",
    "ipopt.bound_relax_factor": 0.0,   # keep 0 <= alpha, rho, tau <= 1 exactly
    "ipopt.print_level": 5,
    "print_time": False,
    # "expand" would flatten the vector graph below into one scalar SX node per
    # (panel, time, wavelength).  That costs about 2 GB of RAM and 70 s on the
    # first solve here, for no change in the answer, so it stays off.
    "expand": False,
}
ipopt_options_quiet = {**ipopt_options, "ipopt.print_level": 0, "ipopt.sb": "yes"}


# =============================================================================
# 2.  MATERIAL MODELS
#     -------------------------------------------------------------------
#     Turns Sections 1.5/1.6 into three arrays per material, on the WORKING
#     wavelength grid:
#
#       alpha_abs[j][k]       alpha_{lambda,i}, for the two absorption terms
#       eps_emit[j][k]        eps_{lambda,i},   for the emission term
#       alpha_elec[j](T)[k]   alpha_{lambda,elec,i}(T), a CasADi function of T
#
#     alpha_abs and eps_emit come from the SAME underlying curve (Kirchhoff),
#     but are band-averaged with different weights -- incident power for what
#     is absorbed, the Planck spectrum for what is emitted.  That is the usual
#     reason a material's solar absorptance differs from its IR emittance, and
#     with these weights the band sums reproduce the native-grid totals exactly.
# =============================================================================

def _logistic(z: np.ndarray) -> np.ndarray:
    """Overflow-safe logistic, for building the fixed material curves."""
    return 0.5 * (1.0 + np.tanh(0.5 * np.clip(z, -60.0, 60.0)))


def parametric_absorptivity(lam_nm: np.ndarray, material: int) -> np.ndarray:
    """alpha_{lambda,i} from the two-plateau panel_spectral entry."""
    try:
        p = panel_spectral[material]
    except KeyError:
        raise ValueError(f"no panel_spectral entry for material type {material}") from None
    s = _logistic((np.log10(lam_nm) - np.log10(p["lam_c_nm"])) / p["blend_decades"])
    return p["alpha_solar"] + (p["alpha_IR"] - p["alpha_solar"]) * s


@dataclass
class PanelProperties:
    """Panel (not filter) spectral data on the working grid."""
    alpha_abs: dict[int, np.ndarray]            # [n_lam], incident-power weighted
    eps_emit: dict[int, np.ndarray]             # [n_lam], Planck weighted
    rho_abs: dict[int, np.ndarray]              # panel reflectivity, incident-power weighted
    rho_emit: dict[int, np.ndarray]             # panel reflectivity, Planck weighted
    alpha_elec: dict[int, ca.Function | None]   # T [K] -> [n_lam]; None if no PV
    source: dict[int, str]
    t_bounds: tuple[float, float]               # (t_min, t_max) actually used
    lookup_info: dict


def load_spectral_lookup(path: Path = lookup_file) -> dict:
    """Read the generator's pickle and check it against the simulation grid."""
    with open(path, "rb") as fh:
        table = pickle.load(fh)
    for key in ("temperature_C", "wavelength_nm", "total_absorptivity", "reflectivity", "transmissivity",
                "electrical_absorptivity_mpp", "valid_optics_mask"):
        if key not in table:
            raise ValueError(f"{path.name} is missing '{key}'")
    return table


def _band_average(native: np.ndarray, idx: np.ndarray, weight: np.ndarray,
                  n_lam: int, lam_native: np.ndarray | None = None,
                  lam_band: np.ndarray | None = None) -> np.ndarray:
    """Weighted mean of a native-grid curve within each working band.

    Three cases, in order of preference:
      1. weight sums to > 0 in the band  -> the weighted mean.  With the
         incident spectrum as the weight this makes the band sums reproduce
         the native-grid totals exactly.
      2. the band holds native bins but zero weight  -> the plain mean (e.g.
         the Planck weight underflows in the visible, where emission is nil).
      3. the band holds no native bins at all  -> interpolate at the band
         centre.  Happens at the short-wavelength end, where a log-spaced
         working grid is finer than the simulation's own 10 nm bins.
    """
    num = np.zeros(n_lam)
    den = np.zeros(n_lam)
    cnt = np.zeros(n_lam)
    flat = np.zeros(n_lam)
    np.add.at(num, idx, native * weight)
    np.add.at(den, idx, weight)
    np.add.at(flat, idx, native)
    np.add.at(cnt, idx, 1.0)

    out = np.zeros(n_lam)
    weighted = den > 0.0
    out[weighted] = num[weighted] / den[weighted]
    plain = (~weighted) & (cnt > 0)
    out[plain] = flat[plain] / cnt[plain]
    empty = cnt == 0
    if empty.any():
        if lam_native is None or lam_band is None:
            out[empty] = 0.0
        else:
            out[empty] = np.interp(lam_band[empty], lam_native, native)
    return out


def build_panel_properties(data: "SimData", verbose: bool = True) -> PanelProperties:
    """Assemble alpha_abs, eps_emit and alpha_elec(T) for every material present."""
    lam_native, idx, n_lam = data.native_lam_nm, data.native_idx, data.n_lam
    materials = sorted(set(int(m) for m in data.material))

    # --- band-averaging weights, on the native grid --------------------------
    # Absorption: incident energy summed over panels and time intervals.
    w_abs = data.native_q_spectrum.copy()
    # Emission: the Planck spectrum at the reference temperature. Work in metres so the
    # coefficients follow directly from h, c and k_B in constants.py.
    lam_m = lam_native * nm
    dlam_m = np.gradient(lam_native) * nm
    x = np.clip(h * c / (lam_m * k * t_set), 0.0, 700.0)
    e = np.exp(-x)
    w_emit = (2.0 * np.pi * h * c ** 2 / lam_m ** 5
              * e / (1.0 - e) * dlam_m)

    alpha_abs, eps_emit, rho_abs, rho_emit, alpha_elec, source = {}, {}, {}, {}, {}, {}
    t_lo, t_hi = t_min, t_max
    info: dict = {}
    aluminium_info = None

    table = None
    if any(panel_spectral_source.get(m) == "lookup" for m in materials):
        table = load_spectral_lookup()
        if not np.array_equal(np.asarray(table["wavelength_nm"]), lam_native):
            raise ValueError("The lookup's wavelength grid does not match the "
                             "simulation's; regenerate it from this pickle.")

    for j in materials:
        src = panel_spectral_source.get(j, "parametric")
        source[j] = src

        if src == "parametric":
            a_native = parametric_absorptivity(lam_native, j)
            r_native = 1.0 - a_native  # The parametric aluminium surface is opaque.
            alpha_elec[j] = None

        elif src == "aluminium_lookup":
            if j != aluminium:
                raise ValueError("The aluminium lookup can only describe the aluminium material type.")
            al_table = load_aluminium_lookup(aluminium_lookup_file, lam_native)
            a_native = np.asarray(al_table["total_absorptivity"], dtype=float).copy()
            r_native = np.asarray(al_table["reflectivity"], dtype=float).copy()
            alpha_elec[j] = None  # Passive metal: no electrical conversion.
            source[j] = "Rakic (1995) Table 2 reflectance + opaque aluminium"
            aluminium_info = al_table["metadata"].copy()

        elif src == "lookup":
            if j != lookup_material:
                raise ValueError(f"material {j} is set to 'lookup' but the lookup "
                                 f"describes material {lookup_material}")
            # The lookup describes the wafer, not the backing behind it.
            a_wafer = np.asarray(table["total_absorptivity"])[0].copy()
            r_native = np.asarray(table["reflectivity"])[0].copy()
            t_native = np.asarray(table["transmissivity"])[0]
            if not np.isfinite([a_wafer, r_native, t_native]).all() or not np.allclose(a_wafer + r_native + t_native, 1.0, atol=1e-8, rtol=0):
                raise ValueError("Lookup optical fractions must be finite and satisfy A + R + T = 1.")
            if min(a_wafer.min(), r_native.min(), t_native.min()) < -1e-8:
                raise ValueError("Lookup optical fractions cannot be negative.")

            # All light leaving the rear of the wafer is absorbed in the backing.
            # Both are in one thermal node. The panel has no rear transmission.
            # Keep the wafer's reflection and electrical response, since the
            # ideal absorbing backing returns no light to the wafer.
            a_native = a_wafer + t_native
            source[j] = "lookup + absorbing backing"

            # ---- alpha_elec(T, lambda) -> a CasADi function of T ------------
            # The Celsius-to-kelvin zero-point offset is definition-based and
            # is not currently provided by constants.py.
            t_grid_k = np.asarray(table["temperature_C"], dtype=float) + 273.15
            e_native = np.asarray(table["electrical_absorptivity_mpp"], dtype=float)
            if not np.isfinite(e_native).all() or e_native.min() < -1e-8 or np.any(e_native > a_wafer[None, :] + 1e-8):
                raise ValueError("Electrical absorption must be finite and lie between zero and wafer absorption.")
            # Band-average each temperature row, then build one interpolant that
            # maps T -> the whole n_lam vector.  CasADi wants output-fastest order.
            e_bands = np.stack([_band_average(row, idx, w_abs, n_lam,
                                              lam_native, data.lam_nm)
                                for row in e_native])        # [n_T, n_lam]
            alpha_elec[j] = ca.interpolant(f"alpha_elec_{j}", t_interp_method,
                                           [t_grid_k], e_bands.ravel(order="C"))
            if clamp_t_to_lookup:
                t_lo, t_hi = max(t_lo, t_grid_k[0]), min(t_hi, t_grid_k[-1])

            valid = np.asarray(table["valid_optics_mask"], dtype=bool)
            info = {
                "solcore_nm": (float(lam_native[valid][0]), float(lam_native[valid][-1])),
                "T_range_K": (float(t_grid_k[0]), float(t_grid_k[-1])),
                "eps_planck_weighted": float(
                    np.sum(a_native * w_emit) / np.sum(w_emit)),
            }
        else:
            raise ValueError(f"panel_spectral_source[{j}] must be 'lookup', "
                             f"'aluminium_lookup' or 'parametric', got {src!r}")

        alpha_abs[j] = _band_average(a_native, idx, w_abs, n_lam,
                                     lam_native, data.lam_nm)
        eps_emit[j] = _band_average(a_native, idx, w_emit, n_lam,
                                    lam_native, data.lam_nm)
        rho_abs[j] = _band_average(r_native, idx, w_abs, n_lam, lam_native, data.lam_nm)
        rho_emit[j] = _band_average(r_native, idx, w_emit, n_lam, lam_native, data.lam_nm)
        if min(r_native.min(), a_native.min()) < -1e-8 or r_native.max() >= 1.0:
            raise ValueError("This stack model requires non-negative panel absorption and 0 <= panel reflectivity < 1.")

    if aluminium_info is not None:
        info["aluminium"] = aluminium_info

    if verbose:
        print("  panel properties :")
        for j in materials:
            print(f"      M = {j}: {source[j]}"
                  + ("  (no PV conversion)" if alpha_elec[j] is None else
                     f"  alpha_elec(T) via {t_interp_method}"))
        if "solcore_nm" in info:
            lo, hi = info["solcore_nm"]
            print(f"      wafer optics from lookup, with an opaque absorbing backing "
                  f"(Solcore covered {lo:.0f} - {hi:.0f} nm)")
            print(f"      Planck-weighted eps at t_set : "
                  f"{info['eps_planck_weighted']:.4f}")
        if aluminium_info is not None:
            print("      aluminium optics: normal-incidence Rakic (1995) data, "
                  "opaque passive panel, no extrapolation")
        n_empty = int(np.sum(np.bincount(idx, minlength=n_lam) == 0))
        if n_empty:
            print(f"      note: {n_empty} of {n_lam} bands hold no simulation bin "
                  f"(below ~{data.lam_nm[n_empty]:.0f} nm the log grid out-resolves "
                  f"the sim's 10 nm bins); those are interpolated and carry no energy")
        print(f"      T bounds used    : {t_lo:.2f} .. {t_hi:.2f} K")

    return PanelProperties(alpha_abs=alpha_abs, eps_emit=eps_emit, rho_abs=rho_abs, rho_emit=rho_emit,
                           alpha_elec=alpha_elec, source=source,
                           t_bounds=(t_lo, t_hi), lookup_info=info)


# =============================================================================
# 3.  DATA LOADING + WAVELENGTH-GRID REDUCTION
# =============================================================================

@dataclass
class SimData:
    """Everything the NLP needs from the simulation, on the WORKING grid."""
    lam_nm: np.ndarray          # band centres (geometric)             [n_lam]
    dlam_nm: np.ndarray         # band widths                          [n_lam]
    lam_edges_nm: np.ndarray    # band edges                           [n_lam + 1]
    t: np.ndarray               # simulation times                     [n_t]
    dt: np.ndarray              # per-sample time weights              [n_t]
    q_in: np.ndarray            # Qdot_in(i, t, lambda)   [W]     [n_p, n_t, n_lam]
    material: np.ndarray        # M_i                                  [n_p]
    labels: list[str]           # panel labels
    mass: np.ndarray            # m_i     [kg]                         [n_p]
    heat_cap: np.ndarray        # c_i     [J/kg/K]                     [n_p]
    area: np.ndarray            # A_i     [m^2]                        [n_p]
    T0: np.ndarray              # T_i(0)  [K]                          [n_p]
    # --- the native (simulation) grid, kept so Section 2 can band-average ----
    native_lam_nm: np.ndarray       # the sim's own wavelengths     [n_native]
    native_idx: np.ndarray          # native bin -> working band    [n_native]
    native_q_spectrum: np.ndarray   # incident energy summed over i and t [n_native]

    @property
    def n_p(self) -> int: return len(self.labels)
    @property
    def n_t(self) -> int: return len(self.t)
    @property
    def n_lam(self) -> int: return len(self.lam_nm)


def _band_edges(lam_native: np.ndarray, n_bands: int | None, spacing: str) -> np.ndarray:
    """Edges of the working grid, covering the native grid exactly."""
    half = 0.5 * np.diff(lam_native)
    native_edges = np.concatenate(([lam_native[0] - half[0]],
                                   lam_native[:-1] + half,
                                   [lam_native[-1] + half[-1]]))
    if n_bands is None or n_bands >= len(lam_native):
        return native_edges
    lo, hi = native_edges[0], native_edges[-1]
    if spacing == "log":
        return np.geomspace(max(lo, 1e-3), hi, n_bands + 1)
    return np.linspace(lo, hi, n_bands + 1)


def load_sim_data(pickle_path: Path = simulation_pickle_file,
                  n_bands: int | None = n_filter_bands,
                  spacing: str = band_spacing) -> SimData:
    """Read the pickle and aggregate the incident power onto the working grid."""
    raw = pd.read_pickle(pickle_path)

    labels: list[str] = list(raw["panelLabels"])
    n_p = len(labels)
    t = np.asarray(raw["simulation_times"], dtype=float)
    lam_native = np.asarray(raw["wavelengths"], dtype=float)
    if t.ndim != 1 or t.size < 2 or not np.isfinite(t).all() or np.any(np.diff(t) <= 0):
        raise ValueError("Simulation times must be finite and strictly increasing.")

    # --- panel types and areas, in the same order as panelLabels ------------
    material = np.asarray(raw["panelTypes"], dtype=int)
    area = np.asarray(raw["panelAreas"], dtype=float)
    print(f"  panel types : {material.tolist()} (from pickle)")

    # Qdot_in(i, t, lambda), W per native bin
    Q_native = np.stack([np.asarray(raw["panelPowers"][lab], dtype=float)
                         for lab in labels])

    # --- aggregate native bins into working bands (powers add) ---------------
    edges = _band_edges(lam_native, n_bands, spacing)
    idx = np.clip(np.searchsorted(edges, lam_native, side="right") - 1, 0, len(edges) - 2)
    n_lam = len(edges) - 1
    q_in = np.zeros((*Q_native.shape[:2], n_lam))
    np.add.at(q_in.transpose(2, 0, 1), idx, Q_native.transpose(2, 0, 1))

    lam = np.sqrt(edges[:-1] * edges[1:]) if spacing == "log" else 0.5 * (edges[:-1] + edges[1:])
    dlam = np.diff(edges)

    # --- per-sample time weights: dt_n for the forward-Euler step n -> n+1 ----
    dt = np.diff(t)
    dt = np.concatenate((dt, [0.0]))        # The final sample has no following interval.

    # --- expand per-material constants to per-panel --------------------------
    def per_panel(table: dict[int, float], key: str) -> np.ndarray:
        vals = np.array([float(table[m]) for m in material])
        for p, over in panel_overrides.items():
            if key in over:
                vals[labels.index(p)] = float(over[key])
        return vals

    return SimData(
        lam_nm=lam, dlam_nm=dlam, lam_edges_nm=edges, t=t, dt=dt,
        q_in=q_in, material=material, labels=labels, area=area,
        mass=per_panel(panel_mass, "mass"),
        heat_cap=per_panel(panel_heat_capacity, "heat_capacity"),
        T0=per_panel(panel_t0, "T0"),
        native_lam_nm=lam_native, native_idx=idx,
        native_q_spectrum=np.sum(Q_native * dt[None, :, None], axis=(0, 1)),
    )


# =============================================================================
# 4.  BLACKBODY SPECTRAL EMISSIVE POWER
#     -------------------------------------------------------------------
#     Planck's law as hemispherical spectral emissive power, so the emitted
#     power is a total over all directions and no intensity enters:
#
#         M_lambda(lambda, T)
#             = 2 pi h c^2 / ( lambda^5 ( exp(h c / (lambda k_B T)) - 1 ) )
#
#     in [W m^-3] with lambda in metres.  The band power is M_lambda dlambda,
#     and summing over a complete spectrum recovers sigma T^4.  The function
#     below accepts wavelengths in nm, converts them using nm, and returns
#     the vector M_lambda(lambda_k, T) dlambda_k [W m^-2].
#
#     Written in the form exp(-x) / (1 - exp(-x)),
#     x = h c / (lambda k_B T), which underflows harmlessly to zero in the far
#     UV instead of overflowing.  No separate radiation constants are stored.
#     T is kept strictly positive by the t_min bound in Section 5.1.
# =============================================================================

def planck_band_power(T, lam_nm: ca.DM, dlam_nm: ca.DM):
    """M_lambda(lambda_k, T) * dlambda_k  [W m^-2], one entry per band.

    Vectorised: this is the emission spectrum itself, evaluated once per
    (panel, time), and then indexed bin-by-bin in the sum of eq. (1).
    """
    lam_m = lam_nm * nm
    dlam_m = dlam_nm * nm
    x = h * c / (lam_m * k * T)
    e = ca.exp(-x)
    M = 2.0 * np.pi * h * c ** 2 / (lam_m ** 5) * e / (1.0 - e)
    return M * dlam_m


def in_band_emission_fraction(T: float, lam_nm: np.ndarray, dlam_nm: np.ndarray) -> float:
    """Share of sigma T^4 that falls inside the modelled band, for reporting."""
    lam_m, dlam_m = lam_nm * nm, dlam_nm * nm
    x = h * c / (lam_m * k * T)
    e = np.exp(-np.clip(x, 0.0, 700.0))
    M = 2.0 * np.pi * h * c ** 2 / lam_m ** 5 * e / (1.0 - e)
    return float(np.sum(M * dlam_m) / (sigma_sb * T ** 4))


def internal_heat_per_panel(data: SimData) -> np.ndarray:
    """Distribute the total internal heat load as constant powers [W per panel]."""
    if not np.isfinite(internal_heat_load_w) or internal_heat_load_w < 0.0:
        raise ValueError("internal_heat_load_w must be finite and non-negative.")

    enabled = (((data.material == aluminium) & internal_heat_to_metal)
               | ((data.material == solar) & internal_heat_to_solar))
    heat = np.zeros(data.n_p)
    if internal_heat_load_w == 0.0 or not np.any(enabled):
        return heat

    weights = data.area[enabled] if internal_heat_by_area else np.ones(np.count_nonzero(enabled))
    if not np.isfinite(weights).all() or np.any(weights <= 0.0):
        raise ValueError("Area-weighted internal heating requires positive finite areas for enabled panels.")
    heat[enabled] = internal_heat_load_w * weights / weights.sum()
    return heat


def simulate_unfiltered(data: SimData, props: PanelProperties) -> tuple[np.ndarray, np.ndarray]:
    """Calculate baseline temperatures [K] and electrical powers [W], shaped [panel, time]."""
    temperature = np.zeros((data.n_p, data.n_t))
    power = np.zeros_like(temperature)
    temperature[:, 0] = data.T0
    lam_nm, dlam_nm = ca.DM(data.lam_nm), ca.DM(data.dlam_nm)
    internal_heat = internal_heat_per_panel(data)

    for i, material in enumerate(data.material):
        mc = data.mass[i] * data.heat_cap[i]
        electrical = props.alpha_elec[material]
        # Operating limits constrain the designed filter, not its comparison
        # reference. A passive aluminium spectrum has no temperature lookup.
        # For PV, never extrapolate the supplied electrical-conversion table.
        lookup_bounds = props.lookup_info.get("T_range_K", props.t_bounds)
        for n in range(data.n_t):
            ti = temperature[i, n]
            if not np.isfinite(ti) or ti <= 0:
                raise ValueError(f"Unfiltered panel {data.labels[i]} reaches invalid temperature {ti:.3f} K at {data.t[n]:.3f} s; check the loads and time step.")
            if electrical is not None and not lookup_bounds[0] <= ti <= lookup_bounds[1]:
                raise ValueError(f"Unfiltered PV panel {data.labels[i]} reaches {ti:.3f} K at {data.t[n]:.3f} s, outside the electrical lookup range {lookup_bounds[0]:.2f} .. {lookup_bounds[1]:.2f} K. Extend and regenerate the lookup; do not extrapolate it.")
            a_elec = np.asarray(electrical(ti)).ravel() if electrical is not None else np.zeros(data.n_lam)
            power[i, n] = np.dot(a_elec, data.q_in[i, n])
            # A clear filter has alpha = rho = 0 and tau = 1.
            q_absorbed = np.dot(props.alpha_abs[material] - a_elec, data.q_in[i, n])
            m_dlam = np.asarray(planck_band_power(ti, lam_nm, dlam_nm)).ravel()
            q_emitted = data.area[i] * np.dot(props.eps_emit[material], m_dlam)
            if n < data.n_t - 1:
                temperature[i, n + 1] = ti + (q_absorbed + internal_heat[i] - q_emitted) * data.dt[n] / mc

    return temperature, power


def stack_optics(alpha_f, rho_f, tau_f, alpha_p, rho_p):
    """Return panel illumination and total absorption for a symmetric filter over a panel.

    Reflections are incoherent. The filter and panel have a common temperature.
    alpha_p includes absorption in the opaque backing, so no rear power escapes.
    For thermal emission, Kirchhoff gives epsilon_stack = alpha_stack.
    """
    panel_illumination = tau_f / (1.0 - rho_f * rho_p)
    filter_absorption = alpha_f * (1.0 + rho_p * panel_illumination)
    stack_absorption = filter_absorption + alpha_p * panel_illumination
    return panel_illumination, stack_absorption


# =============================================================================
# 5.  THE NLP
# =============================================================================

@dataclass
class FilterProblem:
    """The built Opti problem plus handles needed for multi-start and readout."""
    opti: ca.Opti
    data: SimData
    filter_types: list[int]
    alpha_f: dict[int, ca.MX]          # whole-spectrum vectors (set_initial / readout)
    rho_f: dict[int, ca.MX]
    tau_f: dict[int, ca.MX]
    temperature: dict[int, ca.MX]
    J: ca.MX
    e_elec: ca.MX
    e_ref: float
    t_track: ca.MX                    # mean squared normalised temperature error
    p_elec: ca.MX                    # filtered electrical power [panel, time]
    unfiltered_temperature: np.ndarray
    unfiltered_power: np.ndarray
    e_unfiltered: float
    internal_heat_w: np.ndarray       # constant internal heat load [W per panel]
    temperature_bounds: tuple[float, float] = (0.0, float("inf"))


def build_problem(data: SimData, props: PanelProperties, verbose: bool = True,
                  solver_options: dict | None = None) -> FilterProblem:
    """Assemble the NLP.

    The Opti problem AND the IPOPT solver are built exactly once here.  Every
    later solve is just `set_initial` + `solve`, which is what makes
    multi-start cheap: only the first solve pays the (large) one-off cost of
    constructing the Jacobian/Hessian functions.
    """
    if not np.isfinite(t_ref_scale) or t_ref_scale <= 0.0:
        raise ValueError("t_ref_scale must be finite and positive.")
    opti = ca.Opti()

    n_p, n_t, n_lam = data.n_p, data.n_t, data.n_lam
    filter_types = sorted(set(int(m) for m in data.material))

    lam_nm = ca.DM(data.lam_nm)
    dlam_nm = ca.DM(data.dlam_nm)

    # --- fixed material data on the working grid (Section 2) -----------------
    alpha_elec_fn = props.alpha_elec   # T -> alpha_{lambda,elec,i}(T), or None
    t_lo, t_hi = props.t_bounds
    internal_heat = internal_heat_per_panel(data)
    unfiltered_temperature, unfiltered_power = simulate_unfiltered(data, props)
    e_unfiltered = float(np.sum(unfiltered_power * data.dt[None, :]))

    # -------------------------------------------------------------------------
    # 5.1  Decision variables
    #      alpha_j(lambda), rho_j(lambda), tau_j(lambda) for each filter j,
    #      and the auxiliary panel temperatures T_i(t).
    #
    #      Each is declared as one vector over the spectrum and then split into
    #      per-bin scalars, so the equations below can be written with the
    #      index k exactly as the maths does.
    # -------------------------------------------------------------------------
    alpha_f_vec = {j: opti.variable(n_lam) for j in filter_types}
    rho_f_vec = {j: opti.variable(n_lam) for j in filter_types}
    tau_f_vec = {j: opti.variable(n_lam) for j in filter_types}
    temperature_vec = {i: opti.variable(n_t) for i in range(n_p)}

    alpha_f = {j: ca.vertsplit(v) for j, v in alpha_f_vec.items()}     # alpha_f[j][k]
    rho_f = {j: ca.vertsplit(v) for j, v in rho_f_vec.items()}
    tau_f = {j: ca.vertsplit(v) for j, v in tau_f_vec.items()}
    temperature = {i: ca.vertsplit(v)                    # temperature[i][n]
                   for i, v in temperature_vec.items()}

    # Filter emission is included through the combined stack absorptivity below.
    # Kirchhoff gives eps_f = alpha_f for the symmetric reciprocal filter.
    eps_f = alpha_f  # noqa: F841

    for j in filter_types:
        opti.subject_to(opti.bounded(0.0, alpha_f_vec[j], 1.0))
        opti.subject_to(opti.bounded(0.0, rho_f_vec[j], 1.0))
        opti.subject_to(opti.bounded(0.0, tau_f_vec[j], 1.0))
    for i in range(n_p):
        opti.subject_to(opti.bounded(t_lo, temperature_vec[i], t_hi))

    # -------------------------------------------------------------------------
    # 5.2  Conservation of energy, eq. (6)
    #          alpha_j(lambda) + rho_j(lambda) + tau_j(lambda) = 1,  for all lambda
    # -------------------------------------------------------------------------
    for j in filter_types:
        for k in range(n_lam):
            opti.subject_to(alpha_f[j][k] + rho_f[j][k] + tau_f[j][k] == 1.0)

    # -------------------------------------------------------------------------
    # 5.3  Smoothness, eq. (5) and its log-wavelength companion eq. (5')
    #
    #          ( p_j(lambda + dlambda) - p_j(lambda) ) / dlambda       <=  D_{j,p}
    #          ( p_j(lambda + dlambda) - p_j(lambda) ) / dln(lambda)   <=  D'_{j,p}
    #
    #      both applied two-sided, since the intent is a magnitude limit.
    # -------------------------------------------------------------------------
    ln_lam = np.log(data.lam_nm)
    for j in filter_types:
        for name, p in (("alpha", alpha_f[j]), ("rho", rho_f[j]), ("tau", tau_f[j])):
            D = (d_smooth or {}).get(j, {}).get(name)
            D_log = (d_smooth_log or {}).get(j, {}).get(name)
            for k in range(n_lam - 1):
                if D is not None:
                    dlam_k = data.lam_nm[k + 1] - data.lam_nm[k]
                    opti.subject_to(opti.bounded(-D, (p[k + 1] - p[k]) / dlam_k, D))
                if D_log is not None:
                    dln_k = ln_lam[k + 1] - ln_lam[k]
                    opti.subject_to(opti.bounded(-D_log, (p[k + 1] - p[k]) / dln_k, D_log))

    # -------------------------------------------------------------------------
    # 5.4  Heat flows, eq. (1), and the electrical power
    #
    #      g = tau_f / (1 - rho_f rho_p), summing repeated reflections.
    #      alpha_stack = alpha_f (1 + rho_p g) + alpha_p g
    #      P_elec = sum_lambda g alpha_elec(T) Q_in
    #      Q_heat = sum_lambda alpha_stack Q_in - P_elec
    #      Q_emit = area * sum_lambda epsilon_stack M_lambda(T) dlambda
    #      Q_dot = Q_heat + Q_internal_i - Q_emit
    #      Q_internal_i is the same constant load in filtered and unfiltered cases.
    #
    #      epsilon_stack follows TOTAL absorption, not thermal absorption.
    #      Incident and emission band averages use their respective weights.
    # -------------------------------------------------------------------------
    q_dot = {}      # q_dot[i][n]   net thermal power into panel i at time n  [W]
    p_elec = {}     # p_elec[i][n]  electrical power off panel i at time n    [W]

    for i in range(n_p):
        m = int(data.material[i])
        q_dot[i], p_elec[i] = [], []
        alpha_p_col = ca.DM(props.alpha_abs[m])
        eps_p_col = ca.DM(props.eps_emit[m])
        illumination, stack_absorption = stack_optics(alpha_f_vec[m], rho_f_vec[m], tau_f_vec[m], alpha_p_col, ca.DM(props.rho_abs[m]))
        _, stack_emissivity = stack_optics(alpha_f_vec[m], rho_f_vec[m], tau_f_vec[m], eps_p_col, ca.DM(props.rho_emit[m]))

        for n in range(n_t):
            ti = temperature[i][n]
            m_dlam = planck_band_power(ti, lam_nm, dlam_nm)
            a_elec = alpha_elec_fn[m](ti) if alpha_elec_fn[m] is not None else ca.DM.zeros(n_lam)
            q_in_col = ca.DM(data.q_in[i, n])
            electrical_power = ca.sum1(illumination * a_elec * q_in_col)
            absorbed_heat = ca.sum1(stack_absorption * q_in_col) - electrical_power
            emitted_power = ca.sum1(stack_emissivity * m_dlam) * data.area[i]
            q_dot[i].append(absorbed_heat + internal_heat[i] - emitted_power)
            p_elec[i].append(electrical_power)

    # -------------------------------------------------------------------------
    # 5.5  Temperature dynamics, eqs. (2) and (3)
    #          dQ_i(t)  = Qdot_i(t) dt
    #          T_i(t + dt) = T_i(t) + dQ_i(t) / (m c)
    #      with T_i(0) = T_{i,0}.
    # -------------------------------------------------------------------------
    for i in range(n_p):
        opti.subject_to(temperature[i][0] == data.T0[i])
        mc = data.mass[i] * data.heat_cap[i]
        for n in range(n_t - 1):
            dQ = q_dot[i][n] * data.dt[n]
            opti.subject_to(temperature[i][n + 1] == temperature[i][n] + dQ / mc)

    # -------------------------------------------------------------------------
    # 5.6  Cost function, eq. (4)
    #
    #      J = -W_elec * e_elec / e_ref
    #          + W_temp * mean_i,n [((T_i(t_n) - T_set) / dT_ref)^2]
    #          + W_smooth * roughness
    #
    #      The first term maximises electrical energy. The second penalises
    #      temperature deviations on either side of the set-point, equally
    #      weighted across all panels and time samples, as in the original cost.
    #      There is no minimum electrical power relative to the unfiltered case.
    #      Hard temperature bounds still apply. Set w_smooth = 0 for two terms only.
    # -------------------------------------------------------------------------
    power_matrix = ca.vertcat(*[ca.horzcat(*p_elec[i]) for i in range(n_p)])

    e_terms, t_terms = [], []
    for i in range(n_p):
        for n in range(n_t):
            e_terms.append(p_elec[i][n] * data.dt[n])
            t_terms.append(((temperature[i][n] - t_set) / t_ref_scale) ** 2)

    s_terms = []
    for j in filter_types:
        for p in (alpha_f[j], rho_f[j], tau_f[j]):
            for k in range(n_lam - 1):
                s_terms.append((p[k + 1] - p[k]) ** 2)

    e_elec = ca.sum1(ca.vertcat(*e_terms))
    t_track = ca.sum1(ca.vertcat(*t_terms)) / (n_p * n_t)
    roughness = ca.sum1(ca.vertcat(*s_terms)) / (len(filter_types) * n_lam)

    e_ref = float(np.sum(data.q_in * data.dt[None, :, None]))
    if e_ref <= 0.0:
        raise ValueError("The incident energy must be positive to normalise the objective.")
    J = -w_elec * e_elec / e_ref + w_temp * t_track + w_smooth * roughness

    opti.minimize(J)
    opti.solver("ipopt", solver_options if solver_options is not None else ipopt_options_quiet)

    if verbose:
        frac = in_band_emission_fraction(t_set, data.lam_nm, data.dlam_nm)
        print(f"  panels           : {n_p}  {data.labels}")
        print(f"  material types   : {[int(m) for m in data.material]} "
              f"-> filters for {filter_types}")
        print(f"  areas       [m^2]: {[round(float(a), 5) for a in data.area]}")
        print(f"  internal heat [W]: {internal_heat.tolist()}  (per panel, total {internal_heat.sum():.4f} W)")
        print(f"  time samples     : {n_t}   span {data.t[0]:.2f} .. {data.t[-1]:.2f} s")
        print(f"  wavelength bands : {n_lam}  ({data.lam_nm[0]:.1f} .. {data.lam_nm[-1]:.0f} nm)")
        print(f"  incident energy  : e_ref = {e_ref:.4f} J")
        print(f"  unfiltered energy: {e_unfiltered:.4f} J")
        if np.any((unfiltered_temperature < t_lo) | (unfiltered_temperature > t_hi)):
            print(f"  reference T range: {unfiltered_temperature.min():.2f} .. "
                  f"{unfiltered_temperature.max():.2f} K "
                  "(reference exceeds filtered operating bounds)")
        print(f"  temperature cost : mean ((T - {t_set:.2f} K) / {t_ref_scale:g} K)^2")
        print(f"  cost weights     : electrical {w_elec:g}, temperature {w_temp:g}, smoothness {w_smooth:g}")
        print("  power comparison : diagnostic only, no minimum relative to unfiltered")
        print("  radiation model  : isothermal filter and panel, with repeated reflections")
        print("  electrical model : AM0-calibrated conversion fractions, not a new MPP solve")
        no_pv = [j for j in filter_types if alpha_elec_fn[j] is None]
        if no_pv:
            print(f"  note             : filters {no_pv} have no electrical objective or thermal link to PV panels")
        print(f"  in-band emission : {100 * frac:.2f} % of sigma T^4 at {t_set:.2f} K "
              f"(the rest lies beyond {data.lam_edges_nm[-1] / 1000:.0f} um)")

    return FilterProblem(opti=opti, data=data, filter_types=filter_types,
                         alpha_f=alpha_f_vec, rho_f=rho_f_vec, tau_f=tau_f_vec,
                         temperature=temperature_vec, J=J, e_elec=e_elec, e_ref=e_ref, t_track=t_track,
                         p_elec=power_matrix, unfiltered_temperature=unfiltered_temperature,
                         unfiltered_power=unfiltered_power, e_unfiltered=e_unfiltered,
                         internal_heat_w=internal_heat, temperature_bounds=(t_lo, t_hi))


# =============================================================================
# 6.  INITIAL GUESSES
#     -------------------------------------------------------------------
#     Filter guesses satisfy optical conservation and property bounds.
#     Temperatures are clipped to operating bounds; initial guesses need
#     not satisfy the thermal dynamics. `rng=None` is the nominal start.
# =============================================================================

def make_initial_guess(prob: FilterProblem, rng: np.random.Generator | None = None) -> dict:
    """Return {"filters": {j: (alpha, rho, tau)}, "T": {i: trajectory}}."""
    data = prob.data
    n_lam, n_t = data.n_lam, data.n_t

    guess = {"filters": {}, "T": {}}
    for j in prob.filter_types:
        if rng is None:
            # Start with a transparent filter and reference-based temperatures.
            alpha = np.zeros(n_lam)
            tau = np.ones(n_lam)
        else:
            # Random but smooth: a few low-frequency modes in log(lambda),
            # squashed into the simplex.  Respects eq. (5) far better than
            # white noise, so multi-start samples plausible filters.
            u = np.log10(data.lam_nm)
            u = (u - u.min()) / np.ptp(u)
            base = np.zeros((2, n_lam))
            for mode in range(1, 5):
                amp = rng.normal(scale=1.0 / mode, size=2)
                phase = rng.uniform(0, 2 * np.pi, size=2)
                base += amp[:, None] * np.sin(2 * np.pi * mode * u[None, :] + phase[:, None])
            e = np.exp(np.vstack((base, np.zeros(n_lam))))
            simplex = e / e.sum(axis=0)
            alpha, tau = simplex[0], simplex[1]
        rho = np.clip(1.0 - alpha - tau, 0.0, 1.0)
        alpha = np.clip(alpha, 0.0, 1.0)
        tau = np.clip(1.0 - alpha - rho, 0.0, 1.0)
        guess["filters"][j] = (alpha, rho, tau)

    for i in range(data.n_p):
        spread = 0.0 if rng is None else rng.normal(scale=5.0)
        guess["T"][i] = np.clip(prob.unfiltered_temperature[i] + spread,
                                *prob.temperature_bounds)
        guess["T"][i][0] = data.T0[i]
    return guess


def apply_initial_guess(prob: FilterProblem, guess: dict) -> None:
    for j in prob.filter_types:
        alpha, rho, tau = guess["filters"][j]
        prob.opti.set_initial(prob.alpha_f[j], alpha)
        prob.opti.set_initial(prob.rho_f[j], rho)
        prob.opti.set_initial(prob.tau_f[j], tau)
    for i, traj in guess["T"].items():
        prob.opti.set_initial(prob.temperature[i], traj)


# =============================================================================
# 7.  SOLVE / MULTI-START DRIVER
# =============================================================================

@dataclass
class Solution:
    success: bool
    J: float
    e_elec: float
    efficiency: float                  # e_elec / e_ref
    t_track: float                     # mean squared normalised temperature error
    power: np.ndarray                  # filtered electrical power [panel, time]
    min_power_margin_w: float          # min_t (total filtered - total unfiltered), diagnostic only
    filters: dict[int, dict[str, np.ndarray]]
    temperature: dict[str, np.ndarray]
    start: int = 0
    message: str = ""


def _extract(prob: FilterProblem, val: Callable, ok: bool, start: int, msg: str) -> Solution:
    filters = {j: {"alpha": np.asarray(val(prob.alpha_f[j])).ravel(),
                   "rho": np.asarray(val(prob.rho_f[j])).ravel(),
                   "tau": np.asarray(val(prob.tau_f[j])).ravel(),
                   "eps": np.asarray(val(prob.alpha_f[j])).ravel()}   # Kirchhoff
               for j in prob.filter_types}
    temperature = {lab: np.asarray(val(prob.temperature[i])).ravel()
         for i, lab in enumerate(prob.data.labels)}
    E = float(val(prob.e_elec))
    power = np.asarray(val(prob.p_elec)).reshape(prob.data.n_p, prob.data.n_t)
    min_margin = float(np.min(power.sum(axis=0) - prob.unfiltered_power.sum(axis=0)))
    return Solution(success=ok, J=float(val(prob.J)), e_elec=E,
                    efficiency=E / prob.e_ref, t_track=float(val(prob.t_track)),
                    power=power, min_power_margin_w=min_margin,
                    filters=filters, temperature=temperature, start=start, message=msg)


def solve_once(prob: FilterProblem, guess: dict, start: int = 0) -> Solution:
    """One IPOPT solve from `guess`.  Never raises; failures come back flagged.

    The solver was constructed in `build_problem`, so this is only
    `set_initial` + `solve` -- cheap enough to call hundreds of times.
    """
    apply_initial_guess(prob, guess)
    try:
        sol = prob.opti.solve()
        return _extract(prob, sol.value, True, start, "solved")
    except RuntimeError as exc:
        # Recover the last iterate so an "acceptable" point is not thrown away.
        return _extract(prob, prob.opti.debug.value, False, start, str(exc).splitlines()[0])


def multi_start(prob: FilterProblem, starts: int = n_starts, seed: int = random_seed
                ) -> tuple[Solution, list[Solution]]:
    """Solve from `n_starts` initial guesses and keep the best feasible one.

    Start 0 is the deterministic nominal guess; the rest are randomised by
    `make_initial_guess`.  Only the first solve pays the solver-construction
    cost, so additional starts are comparatively cheap.
    """
    rng = np.random.default_rng(seed)
    results: list[Solution] = []
    for s in range(starts):
        guess = make_initial_guess(prob, None if s == 0 else rng)
        t0 = time.perf_counter()
        res = solve_once(prob, guess, start=s)
        results.append(res)
        print(f"  start {s:>3d}: {'ok ' if res.success else 'FAIL'}  "
              f"J = {res.J: .6f}   eta = {100 * res.efficiency:6.3f} %   "
              f"T-track = {res.t_track:.4f}   [{time.perf_counter() - t0:5.1f} s]",
              flush=True)

    feasible = [r for r in results if r.success]
    if not feasible:
        raise RuntimeError("No solve satisfied the constraints. No filter solution has been saved.")
    best = min(feasible, key=lambda r: r.J)
    return best, results


# =============================================================================
# 8.  RESULTS
# =============================================================================

material_names = {solar: "solar-panel filter", aluminium: "aluminium filter"}


def report(best: Solution, prob: FilterProblem) -> None:
    data = prob.data
    print("\n" + "=" * 70)
    print(f"BEST SOLUTION  (start {best.start}, {best.message})")
    print("=" * 70)
    print(f"  J                      {best.J: .6f}")
    print(f"  electrical energy      {best.e_elec:.4f} J  "
          f"over {data.t[-1] - data.t[0]:.1f} s")
    print(f"  conversion efficiency  {100 * best.efficiency:.3f} %  of incident energy")
    print(f"  mean (dT/dT_ref)^2      {best.t_track:.5f}")
    print(f"  unfiltered energy      {prob.e_unfiltered:.4f} J")
    print(f"  electrical gain        {best.e_elec - prob.e_unfiltered:.6f} J")
    print(f"  minimum power difference {best.min_power_margin_w:.3e} W  (filtered - unfiltered, diagnostic only)")

    print("\n  panel   type   A [m^2]   T start -> end [K]    unfiltered end [K]")
    for i, (lab, traj) in enumerate(best.temperature.items()):
        m = int(data.material[i])
        print(f"    {lab:>3s}      {m}     {data.area[i]:.4f}   "
              f"{traj[0]:8.3f} -> {traj[-1]:8.3f}      {prob.unfiltered_temperature[i, -1]:8.3f}")

    print("\n  filter band averages")
    for j, f in best.filters.items():
        print(f"    {material_names.get(j, f'filter {j}')} (M = {j}):")
        for lo, hi, tag in ((300, 1150, "PV band   "),
                            (1150, 3000, "solar IR  "),
                            (3000, 50000, "thermal IR")):
            m = (data.lam_nm >= lo) & (data.lam_nm < hi)
            if not m.any():
                continue
            print(f"      {tag} {lo:>5.0f}-{hi:<6.0f} nm : "
                  f"alpha {f['alpha'][m].mean():.3f}   "
                  f"rho {f['rho'][m].mean():.3f}   "
                  f"tau {f['tau'][m].mean():.3f}")


def solution_payload(best: Solution, results: list[Solution], prob: FilterProblem,
                     props: PanelProperties) -> dict:
    """Collect the solved data for saving or plotting without rebuilding the model."""
    return {
        "lam_nm": prob.data.lam_nm,
        "dlam_nm": prob.data.dlam_nm,
        "lam_edges_nm": prob.data.lam_edges_nm,
        "t": prob.data.t,
        "material": prob.data.material,
        "panel_labels": prob.data.labels,
        "area": prob.data.area,
        "mass_kg": prob.data.mass,
        "specific_heat_J_kg_K": prob.data.heat_cap,
        "thermal_capacitance_J_K": prob.data.mass * prob.data.heat_cap,
        "filters": best.filters,
        "panel_alpha_abs": props.alpha_abs,
        "panel_eps_emit": props.eps_emit,
        "panel_alpha_elec": {j: np.stack([np.asarray(fn(Tq)).ravel() for Tq in props.t_bounds])
                             for j, fn in props.alpha_elec.items() if fn is not None},
        "panel_rho_abs": props.rho_abs,
        "panel_rho_emit": props.rho_emit,
        "panel_source": props.source,
        "t_bounds": props.t_bounds,
        "lookup_info": props.lookup_info,
        "model_assumptions": {
            "stack": "Symmetric reciprocal filter and panel at one temperature, incoherent reflections",
            "emission": "Combined stack absorptivity, normal-incidence approximation, outward face only",
            "electrical": "Supplied one-sun AM0 MPP allocation, temperature dependent",
            "thermal_coupling": "No conduction between panels or to a spacecraft body",
            "internal_heat": "Prescribed constant heat load shared among enabled panel types, identical in filtered and unfiltered cases",
            "transmission": "All wafer transmission is absorbed in an opaque backing at the panel temperature",
            "backing": "Perfect absorber, no return reflection, included in the panel thermal mass",
            "spectral_extent": "No thermal radiation outside the simulation wavelength range",
            "time_integration": "Forward Euler intervals with zero weight on the final sample",
        },
        "T": best.temperature,
        "J": best.J,
        "e_elec": best.e_elec,
        "t_track": best.t_track,
        "electrical_power_W": best.power,
        "unfiltered_electrical_power_W": prob.unfiltered_power,
        "unfiltered_temperature_K": prob.unfiltered_temperature,
        "unfiltered_electrical_energy_J": prob.e_unfiltered,
        "internal_heat_per_panel_W": prob.internal_heat_w,
        "min_power_margin_W": best.min_power_margin_w,
        "efficiency": best.efficiency,
        "all_starts": [(r.start, r.success, r.J, r.efficiency) for r in results],
        "settings": {"t_set": t_set, "t_ref_scale": t_ref_scale,
                     "w_elec": w_elec, "w_temp": w_temp,
                     "internal_heat_load_w": internal_heat_load_w,
                     "internal_heat_to_metal": internal_heat_to_metal,
                     "internal_heat_to_solar": internal_heat_to_solar,
                     "internal_heat_by_area": internal_heat_by_area,
                     "power_constraint": None,
                     "w_smooth": w_smooth, "d_smooth": d_smooth,
                     "d_smooth_log": d_smooth_log,
                     "n_filter_bands": n_filter_bands, "band_spacing": band_spacing},
    }


def save(best: Solution, results: list[Solution], prob: FilterProblem,
         props: PanelProperties, path: Path) -> None:
    payload = solution_payload(best, results, prob, props)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(payload, path)
    print(f"\n  saved -> {path}")


def plot(payload: dict, path: Path) -> None:
    """Plot saved result arrays. Older pickles can omit electrical-property curves."""
    lam_nm, t = payload["lam_nm"], payload["t"]
    filters = payload["filters"]
    n = len(filters)
    has_power = "electrical_power_W" in payload
    n_axes = n + 2 + int(has_power)
    fig, axes = plt.subplots(1, n_axes, figsize=(5.8 * n_axes, 4.2))
    axes = np.atleast_1d(axes)

    for ax, (j, f) in zip(axes[:n], filters.items()):
        for key, style in (("alpha", "-"), ("rho", "--"), ("tau", "-.")):
            if key == "alpha":
                ax.plot(lam_nm, f[key], style, label=rf"$\{key}_{{{j}}} = \epsilon_{{{j}}}$")
            else:
                ax.plot(lam_nm, f[key], style, label=rf"$\{key}_{{{j}}}$")
        ax.set_xscale("log")
        ax.set_xlabel("wavelength [nm]")
        ax.set_ylabel("property")
        ax.set_title(f"{material_names.get(j, f'filter {j}')}  (M = {j})")
        ax.set_ylim(-0.05, 1.05)
        ax.grid(alpha=0.3)
        ax.legend()

    # --- the panel properties the filters are designed against ---------------
    ax = axes[n]
    for j in sorted(payload["panel_alpha_abs"]):
        tag = f"M={j} ({payload['panel_source'][j]})"
        ax.plot(lam_nm, payload["panel_alpha_abs"][j], "-", label=rf"$\alpha_{{p,{j}}}$ {tag}")
        ax.plot(lam_nm, payload["panel_eps_emit"][j], ":", label=rf"$\epsilon_{{p,{j}}}$ {tag}")
        if j in payload.get("panel_alpha_elec", {}):
            for Tq, curve, style in zip(payload["t_bounds"], payload["panel_alpha_elec"][j], ("--", "-.")):
                ax.plot(lam_nm, curve, style,
                        label=rf"$\alpha_{{elec,{j}}}$ @ {Tq:.0f} K")
    ax.set_xscale("log")
    ax.set_xlabel("wavelength [nm]")
    ax.set_ylabel("property")
    ax.set_title("Panel properties (under the filter)")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)

    ax = axes[n + 1]
    for i, lab in enumerate(payload["panel_labels"]):
        line, = ax.plot(t, payload["T"][lab], label=f"{lab} filtered")
        if "unfiltered_temperature_K" in payload:
            ax.plot(t, payload["unfiltered_temperature_K"][i], "--", color=line.get_color(), alpha=0.7, label=f"{lab} unfiltered")
    if "t_set" in payload.get("settings", {}):
        ax.axhline(payload["settings"]["t_set"], color="k", ls=":", label=r"$T_\mathrm{set}$")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("panel temperature [K]")
    ax.set_title("Temperature trajectories")
    ax.grid(alpha=0.3)
    ax.legend(ncol=2, fontsize=7)

    if has_power:
        ax = axes[n + 2]
        ax.plot(t, np.asarray(payload["electrical_power_W"]).sum(axis=0), label="Filtered")
        if "unfiltered_electrical_power_W" in payload:
            ax.plot(t, np.asarray(payload["unfiltered_electrical_power_W"]).sum(axis=0), "--", label="Unfiltered")
        ax.set_xlabel("time [s]")
        ax.set_ylabel("total electrical power [W]")
        ax.set_title("Electrical power comparison")
        ax.grid(alpha=0.3)
        ax.legend()

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"  saved -> {path}")


# =============================================================================
#  MAIN  -- behaviour is controlled entirely by Section 1.0
# =============================================================================

def main() -> None:
    if not solve_new_filter:
        print(f"Loading saved optimisation {output_file} ...")
        plot(pd.read_pickle(output_file), plot_file)
        return

    print(f"Loading {simulation_pickle_file.name} and {lookup_file.name} ...")
    data = load_sim_data()
    props = build_panel_properties(data)

    print("Building NLP ...")
    t0 = time.perf_counter()
    prob = build_problem(data, props, solver_options=ipopt_options if verbose_ipopt
                         else ipopt_options_quiet)
    print(f"built in {time.perf_counter() - t0:.1f} s")

    print(f"Solving ({n_starts} start{'s' if n_starts > 1 else ''}) ... "
          f"the first solve also builds the derivative functions")
    best, results = multi_start(prob, starts=n_starts, seed=random_seed)

    report(best, prob)
    if save_solution:
        save(best, results, prob, props, output_file)
    if make_plot:
        plot(solution_payload(best, results, prob, props), plot_file)


if __name__ == "__main__":
    main()
