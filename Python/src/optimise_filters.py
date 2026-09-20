"""
optimise_filters.py
===================

Spectral-filter design for the rotating-satellite simulation, solved with
CasADi + IPOPT.

The script is laid out so that every block maps one-to-one onto a numbered
equation in ``Solar Panel Math``.  Summations in the maths are written as
explicit ``for`` loops so the code can be read side-by-side with the LaTeX and
edited in the same places.

    Section 1   run configuration and constants NOT in the pickle  (edit these)
    Section 2   material shape functions (numbers live in 1.5 and 1.6)
    Section 3   data loading + wavelength-grid reduction
    Section 4   blackbody spectral emissive power  M_lambda(lambda, T)
    Section 5   the NLP, equation by equation
    Section 6   initial guesses (nominal + randomised, for multi-start)
    Section 7   solve / multi-start driver
    Section 8   results

There are no command-line arguments: everything is set in Section 1.0 and the
script is run as-is.

Decision variables (per filter type j, per wavelength bin k):
    alpha_j(lambda),  rho_j(lambda),  tau_j(lambda)
    eps_j(lambda) = alpha_j(lambda)                       (Kirchhoff, eliminated)
Auxiliary variables:
    T_i(t)                                                (panel temperatures)

Equation map
------------
    eq (1)  heat flows / spectral sum          ->  Section 5.4
    eq (2)  dQ_i = Qdot_i dt                   ->  Section 5.5
    eq (3)  T_i(t+dt) = T_i(t) + dQ_i/(mc)     ->  Section 5.5
    eq (4)  cost J                             ->  Section 5.6
    eq (5)  smoothness  |dp/dlambda| <= D      ->  Section 5.3
    eq (6)  alpha + rho + tau = 1              ->  Section 5.2
    Kirchhoff  eps_j = alpha_j                 ->  Section 5.1 (by elimination)

Emission (Section 4).  The panel emits the blackbody spectral emissive power
scaled by its own emissivity and the filter transmissivity:

    Qdot_out(t, lambda) = -eps_{lambda,i} tau_{lambda,f_i} M_lambda(lambda, T_i) dlambda A_i

with M_lambda the hemispherical spectral emissive power [W m^-2 um^-1], so
direction and intensity do not enter and

    sum_{all lambda} M_lambda dlambda  ->  sigma T^4

over a complete spectrum.  The simulation's band (15 nm - 50 um) holds about
96.6 % of that at 293 K; the report prints the exact figure.

One departure from the literal spec remains, flagged where it occurs: a
companion smoothness limit in ln(lambda), eq (5'), and a small roughness
tie-breaker in J.  Both exist only to pick a sensible member of a degenerate
optimal set; see Sections 1.8 and 5.6.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import casadi as ca
import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent

# =============================================================================
# 1.  RUN CONFIGURATION AND CONSTANTS NOT PRESENT IN THE PICKLE
#     -------------------------------------------------------------------
#     The pickle supplies:  Qdot_in(t, lambda) per panel [W per bin], the
#     wavelength grid and the time grid.  Everything else is set here.
#     Values flagged  <<CHECK>>  are estimates, not measured data.
# =============================================================================

# ---- 1.0  Run configuration --------------------------------------------------
PICKLE_FILE = _HERE / "Simulation Data" / "simulation_data.pkl"
OUTPUT_FILE = _HERE / "filter_solution.pkl"
PLOT_FILE = _HERE / "filter_solution.png"

N_STARTS = 5                    # multi-start points; 1 = single solve
RANDOM_SEED = 0                 # seed for the randomised starts
SAVE_SOLUTION = True            # write OUTPUT_FILE
MAKE_PLOT = True                # write PLOT_FILE (needs matplotlib)
VERBOSE_IPOPT = False           # True prints the full IPOPT iteration log

# ---- 1.1  Physical constants -------------------------------------------------
SIGMA = 5.670374419e-8          # Stefan-Boltzmann constant            [W m^-2 K^-4]
C1_UM4 = 3.741771852e8          # 1st radiation const, 2 pi h c^2  [W um^4 m^-2]
C2_UM_K = 1.4387769e4           # 2nd radiation constant               [um K]

# ---- 1.2  Material type labels (as used by the C++ sim) ----------------------
SOLAR = 0                       # M_i = 0  ->  solar panel
ALUMINIUM = 1                   # M_i = 1  ->  aluminium

# ---- 1.3  Per-panel data missing from the pickle -----------------------------
# Both are listed in the pickle's panelLabels order: ["+x","-x","+y","-y","+z","-z"].
# Set either to None to fall back (types -> the pickle's panelTypes,
# areas -> PANEL_AREA_DEFAULT).
PANEL_TYPES = [SOLAR, SOLAR, SOLAR, SOLAR, ALUMINIUM, ALUMINIUM]     # <<CHECK>>
PANEL_AREAS = [0.02, 0.02, 0.02, 0.02, 0.01, 0.01]      # A_i  [m^2]  <<CHECK>>
PANEL_AREA_DEFAULT = 0.01                               #             [m^2]

# NOTE ON AREA:  the pickle's panel powers are already absolute powers [W]
# (i.e. area-integrated), so A_i enters ONLY through the emission term.  If the
# sim ever sends fluxes [W m^-2] instead, multiply Q_in by A_i in Section 3.

# ---- 1.4  Per-material thermal constants -------------------------------------
PANEL_MASS = {SOLAR: 0.120, ALUMINIUM: 0.090}            # m  [kg]        <<CHECK>>
PANEL_HEAT_CAPACITY = {SOLAR: 700.0, ALUMINIUM: 900.0}   # c  [J/kg/K]    <<CHECK>>
PANEL_T0 = {SOLAR: 293.15, ALUMINIUM: 293.15}            # T_i(0)  [K]    <<CHECK>>

# Per-panel overrides, keyed by panel label; keys: mass, heat_capacity, T0.
PANEL_OVERRIDES: dict[str, dict[str, float]] = {
    # "+x": {"mass": 0.15, "T0": 300.0},
}

# ---- 1.5  Per-material PANEL spectral properties -----------------------------
# The material UNDER the filter.  Fixed data, not decision variables.
# Kirchhoff on the panel as well, so eps_{lambda,i} = alpha_{lambda,i} and one
# curve serves both.  Each is a smooth two-plateau blend in log(lambda):
#
#     alpha(lambda) = alpha_solar  ->  alpha_IR,  crossing over at lam_c_nm
#
# Swap in a measured (lambda, alpha) table by editing panel_absorptivity() in
# Section 2; it is a two-line change and everything downstream is unaffected.
PANEL_SPECTRAL = {
    SOLAR: {                       # solar cell + coverglass            <<CHECK>>
        "alpha_solar": 0.92,       # absorptivity shortward of the crossover
        "alpha_IR": 0.85,          # absorptivity longward of the crossover
        "lam_c_nm": 3000.0,        # crossover wavelength                   [nm]
        "blend_decades": 0.25,     # crossover width, in decades of lambda
    },
    ALUMINIUM: {                   # polished / anodised aluminium      <<CHECK>>
        "alpha_solar": 0.15,
        "alpha_IR": 0.05,
        "lam_c_nm": 3000.0,
        "blend_decades": 0.25,
    },
}

# ---- 1.6  Electrical absorptivity,  alpha_elec(lambda, T) --------------------
#     alpha_elec(lambda, T) = alpha_elec_ref(lambda) * g(T)
#     alpha_therm(lambda, T) = alpha(lambda) - alpha_elec(lambda, T)
#
# alpha_elec_ref is a silicon-like EQE shape: zero in the UV, a plateau through
# the visible/NIR, cut off at the band gap.  Set a material to None for no PV.
ALPHA_ELEC_SHAPE = {
    SOLAR: {                       # silicon-like                       <<CHECK>>
        "peak": 0.45,              # peak fraction of incident power -> electricity
        "lam_rise_nm": 340.0,      # UV roll-on centre                      [nm]
        "rise_width_nm": 40.0,     # UV roll-on width                       [nm]
        "lam_cut_nm": 1110.0,      # band-gap cut-off centre                [nm]
        "cut_width_nm": 25.0,      # band-gap cut-off width                 [nm]
    },
    ALUMINIUM: None,               # no PV conversion
}

# g(T) is the "polynomial-type function that is pre-calculated" in the spec.
# Coefficients are on (T - T_ELEC_REF); c1 is the usual silicon power
# temperature coefficient of about -0.45 %/K.
T_ELEC_REF = 298.15                                    # [K]
ALPHA_ELEC_TEMP_POLY = {                               # [c0, c1, c2, ...]
    SOLAR:     [1.0, -4.5e-3, 0.0],
    ALUMINIUM: [0.0, 0.0, 0.0],                        # no PV conversion
}

# ---- 1.7  Design targets and cost weights ------------------------------------
T_SET = 293.15                  # T_set, temperature set-point         [K]
T_REF_SCALE = 10.0              # dT_ref, normalises the tracking term [K]
W_ELEC = 1.0                    # weight on electrical energy (maximised)
W_TEMP = 1.0                    # weight on temperature tracking (minimised)
W_SMOOTH = 1.0e-3               # tie-breaker on filter roughness (see 5.6)

T_MIN, T_MAX = 150.0, 450.0     # hard bounds on T_i(t)                [K]

# ---- 1.8  Smoothness limits --------------------------------------------------
# eq. (5):   |dp_j / dlambda|      <= D_{j,p}     units: property per nm.
#            D = 1e-2 lets a property swing 0 -> 1 over 100 nm.
D_SMOOTH = {
    SOLAR:     {"alpha": 1.0e-2, "rho": 1.0e-2, "tau": 1.0e-2},
    ALUMINIUM: {"alpha": 1.0e-2, "rho": 1.0e-2, "tau": 1.0e-2},
}

# eq. (5'):  |dp_j / dln(lambda)|  <= D'_{j,p}    dimensionless.
# A per-nm limit is ~1e5 times looser at 50 um than at 500 nm, so on its own
# eq. (5) leaves the thermal IR unconstrained and the solution chatters
# bin-to-bin there (the emission depends only on the band-average eps*tau, so
# that part of the problem is degenerate).  A real interference filter
# transitions over a fixed *fraction* of lambda, which is what this limits.
# D' = 5 allows a 0 -> 1 swing over roughly a factor of 1.2 in wavelength.
# Set the whole dict, or any single property, to None to use eq. (5) alone.
D_SMOOTH_LOG = {
    SOLAR:     {"alpha": 5.0, "rho": 5.0, "tau": 5.0},
    ALUMINIUM: {"alpha": 5.0, "rho": 5.0, "tau": 5.0},
}

# ---- 1.9  Wavelength grid reduction ------------------------------------------
# The sim ships 4999 bins.  A manufacturable filter cannot have 4999
# independent values, and the NLP is far better conditioned on a coarser grid.
# N_FILTER_BANDS = None  ->  use the native grid unchanged (slow: ~30k vars).
# At 120 log-spaced bands the emission quadrature error is ~2.5e-4.
N_FILTER_BANDS = 120
BAND_SPACING = "log"            # "log" or "linear"

# ---- 1.10  Solver options ----------------------------------------------------
IPOPT_OPTIONS = {
    "ipopt.max_iter": 3000,
    "ipopt.tol": 1e-8,
    "ipopt.acceptable_tol": 1e-6,
    "ipopt.mu_strategy": "adaptive",
    "ipopt.hessian_approximation": "exact",
    "ipopt.bound_relax_factor": 0.0,   # keep 0 <= alpha, rho, tau <= 1 exactly
    "ipopt.print_level": 5,
    "print_time": False,
    # Expand the MX graph to SX once, up front.  Costs a few seconds on the
    # first solve and makes every later multi-start solve ~5x faster.
    "expand": True,
}
IPOPT_OPTIONS_QUIET = {**IPOPT_OPTIONS, "ipopt.print_level": 0, "ipopt.sb": "yes"}


# =============================================================================
# 2.  MATERIAL MODELS
#     -------------------------------------------------------------------
#     Shape functions only.  Every number lives in Sections 1.5 and 1.6; these
#     just evaluate the tables onto the working wavelength grid, in nm.
# =============================================================================

def _logistic(z: np.ndarray) -> np.ndarray:
    """Overflow-safe logistic, for building the fixed material curves."""
    return 0.5 * (1.0 + np.tanh(0.5 * np.clip(z, -60.0, 60.0)))


def panel_absorptivity(lam_nm: np.ndarray, material: int) -> np.ndarray:
    """alpha_{lambda, i}: total (thermal + electrical) absorptivity of the panel.

    Two plateaus blended in log(lambda), per the PANEL_SPECTRAL entry.  To use
    measured data instead, replace the body with an interpolation, e.g.
        return np.interp(lam_nm, MEASURED_LAM_NM[material], MEASURED_ALPHA[material])
    """
    try:
        p = PANEL_SPECTRAL[material]
    except KeyError:
        raise ValueError(f"no PANEL_SPECTRAL entry for material type {material}") from None
    s = _logistic((np.log10(lam_nm) - np.log10(p["lam_c_nm"])) / p["blend_decades"])
    return p["alpha_solar"] + (p["alpha_IR"] - p["alpha_solar"]) * s


def panel_emissivity(lam_nm: np.ndarray, material: int) -> np.ndarray:
    """eps_{lambda, i}.  Kirchhoff on the panel as well: eps = alpha."""
    return panel_absorptivity(lam_nm, material)


def alpha_elec_reference(lam_nm: np.ndarray, material: int) -> np.ndarray:
    """alpha_elec_ref(lambda): spectral conversion efficiency at T = T_ELEC_REF.

    UV roll-on times a band-gap cut-off, per the ALPHA_ELEC_SHAPE entry.  To
    use a measured EQE curve, replace the body with an interpolation.
    """
    p = ALPHA_ELEC_SHAPE.get(material)
    if p is None:                                      # not a PV material
        return np.zeros_like(lam_nm)
    rise = _logistic((lam_nm - p["lam_rise_nm"]) / p["rise_width_nm"])
    fall = _logistic(-(lam_nm - p["lam_cut_nm"]) / p["cut_width_nm"])
    return p["peak"] * rise * fall


def g_elec(T, material: int):
    """Temperature factor of the electrical absorptivity (CasADi-safe)."""
    coeffs = ALPHA_ELEC_TEMP_POLY[material]
    dT = T - T_ELEC_REF
    out = 0.0
    for power, c in enumerate(coeffs):
        out = out + c * dT ** power
    return out


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
    Q_in: np.ndarray            # Qdot_in(i, t, lambda)   [W]     [n_p, n_t, n_lam]
    material: np.ndarray        # M_i                                  [n_p]
    labels: list[str]           # panel labels
    mass: np.ndarray            # m_i     [kg]                         [n_p]
    heat_cap: np.ndarray        # c_i     [J/kg/K]                     [n_p]
    area: np.ndarray            # A_i     [m^2]                        [n_p]
    T0: np.ndarray              # T_i(0)  [K]                          [n_p]

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


def load_sim_data(pickle_path: Path = PICKLE_FILE,
                  n_bands: int | None = N_FILTER_BANDS,
                  spacing: str = BAND_SPACING) -> SimData:
    """Read the pickle and aggregate the incident power onto the working grid."""
    raw = pd.read_pickle(pickle_path)

    labels: list[str] = list(raw["panelLabels"])
    n_p = len(labels)
    t = np.asarray(raw["simulation_times"], dtype=float)
    lam_native = np.asarray(raw["wavelengths"], dtype=float)

    # --- panel types: Section 1.3 wins over the pickle -----------------------
    if PANEL_TYPES is None:
        material = np.asarray(raw["panelTypes"], dtype=int)
    else:
        if len(PANEL_TYPES) != n_p:
            raise ValueError(f"PANEL_TYPES has {len(PANEL_TYPES)} entries, "
                             f"but the pickle has {n_p} panels {labels}")
        material = np.asarray(PANEL_TYPES, dtype=int)

    # --- panel areas ---------------------------------------------------------
    if PANEL_AREAS is None:
        area = np.full(n_p, float(PANEL_AREA_DEFAULT))
    else:
        if len(PANEL_AREAS) != n_p:
            raise ValueError(f"PANEL_AREAS has {len(PANEL_AREAS)} entries, "
                             f"but the pickle has {n_p} panels {labels}")
        area = np.asarray(PANEL_AREAS, dtype=float)

    # Qdot_in(i, t, lambda), W per native bin
    Q_native = np.stack([np.asarray(raw["panelPowers"][lab], dtype=float)
                         for lab in labels])

    # --- aggregate native bins into working bands (powers add) ---------------
    edges = _band_edges(lam_native, n_bands, spacing)
    idx = np.clip(np.searchsorted(edges, lam_native, side="right") - 1, 0, len(edges) - 2)
    n_lam = len(edges) - 1
    Q_in = np.zeros((*Q_native.shape[:2], n_lam))
    np.add.at(Q_in.transpose(2, 0, 1), idx, Q_native.transpose(2, 0, 1))

    lam = np.sqrt(edges[:-1] * edges[1:]) if spacing == "log" else 0.5 * (edges[:-1] + edges[1:])
    dlam = np.diff(edges)

    # --- per-sample time weights: dt_n for the forward-Euler step n -> n+1 ----
    dt = np.diff(t)
    dt = np.concatenate((dt, [dt[-1]]))     # hold the last step for the energy sum

    # --- expand per-material constants to per-panel --------------------------
    def per_panel(table: dict[int, float], key: str) -> np.ndarray:
        vals = np.array([float(table[m]) for m in material])
        for p, over in PANEL_OVERRIDES.items():
            if key in over:
                vals[labels.index(p)] = float(over[key])
        return vals

    return SimData(
        lam_nm=lam, dlam_nm=dlam, lam_edges_nm=edges, t=t, dt=dt,
        Q_in=Q_in, material=material, labels=labels, area=area,
        mass=per_panel(PANEL_MASS, "mass"),
        heat_cap=per_panel(PANEL_HEAT_CAPACITY, "heat_capacity"),
        T0=per_panel(PANEL_T0, "T0"),
    )


# =============================================================================
# 4.  BLACKBODY SPECTRAL EMISSIVE POWER
#     -------------------------------------------------------------------
#     Planck's law as hemispherical spectral emissive power, so the emitted
#     power is a total over all directions and no intensity enters:
#
#         M_lambda(lambda, T) = C1 / ( lambda^5 ( exp(C2 / (lambda T)) - 1 ) )
#
#     in [W m^-2 um^-1] with lambda in um.  The band power is M_lambda dlambda,
#     and summing over a complete spectrum recovers sigma T^4.  The function
#     below returns the whole vector  M_lambda(lambda_k, T) dlambda_k  [W m^-2].
#
#     Written in the form  exp(-x) / (1 - exp(-x)),  x = C2 / (lambda T), which
#     underflows harmlessly to zero in the far UV instead of overflowing.
#     T is kept strictly positive by the T_MIN bound in Section 5.1.
# =============================================================================

def planck_band_power(T, lam_um: ca.DM, dlam_um: ca.DM):
    """M_lambda(lambda_k, T) * dlambda_k  [W m^-2], one entry per band.

    Vectorised: this is the emission spectrum itself, evaluated once per
    (panel, time), and then indexed bin-by-bin in the sum of eq. (1).
    """
    x = C2_UM_K / (lam_um * T)
    e = ca.exp(-x)
    M = C1_UM4 / (lam_um ** 5) * e / (1.0 - e)
    return M * dlam_um


def in_band_emission_fraction(T: float, lam_nm: np.ndarray, dlam_nm: np.ndarray) -> float:
    """Share of sigma T^4 that falls inside the modelled band, for reporting."""
    lam_um, dlam_um = lam_nm * 1e-3, dlam_nm * 1e-3
    x = C2_UM_K / (lam_um * T)
    e = np.exp(-np.clip(x, 0.0, 700.0))
    M = C1_UM4 / lam_um ** 5 * e / (1.0 - e)
    return float(np.sum(M * dlam_um) / (SIGMA * T ** 4))


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
    T: dict[int, ca.MX]
    J: ca.MX
    E_elec: ca.MX
    E_ref: float
    T_track: ca.MX


def build_problem(data: SimData, verbose: bool = True,
                  ipopt_options: dict | None = None) -> FilterProblem:
    """Assemble the NLP.

    The Opti problem AND the IPOPT solver are built exactly once here.  Every
    later solve is just `set_initial` + `solve`, which is what makes
    multi-start cheap: only the first solve pays the (large) one-off cost of
    constructing the Jacobian/Hessian functions.
    """
    opti = ca.Opti()

    n_p, n_t, n_lam = data.n_p, data.n_t, data.n_lam
    filter_types = sorted(set(int(m) for m in data.material))

    lam_um = ca.DM(data.lam_nm * 1e-3)
    dlam_um = ca.DM(data.dlam_nm * 1e-3)

    # --- fixed material data on the working grid -----------------------------
    alpha_p = {j: panel_absorptivity(data.lam_nm, j) for j in filter_types}
    eps_p = {j: panel_emissivity(data.lam_nm, j) for j in filter_types}
    alpha_el = {j: alpha_elec_reference(data.lam_nm, j) for j in filter_types}

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
    T_vec = {i: opti.variable(n_t) for i in range(n_p)}

    alpha_f = {j: ca.vertsplit(v) for j, v in alpha_f_vec.items()}     # alpha_f[j][k]
    rho_f = {j: ca.vertsplit(v) for j, v in rho_f_vec.items()}
    tau_f = {j: ca.vertsplit(v) for j, v in tau_f_vec.items()}
    T = {i: ca.vertsplit(v) for i, v in T_vec.items()}                 # T[i][n]

    # Kirchhoff on the filter, eps_j(lambda) = alpha_j(lambda), is imposed by
    # elimination: there is no separate eps variable.  It does not appear in
    # the heat balance because the filter's own emission cancels (spec rule 3);
    # it is reported in Section 8 as eps = alpha.
    eps_f = alpha_f  # noqa: F841

    for j in filter_types:
        opti.subject_to(opti.bounded(0.0, alpha_f_vec[j], 1.0))
        opti.subject_to(opti.bounded(0.0, rho_f_vec[j], 1.0))
        opti.subject_to(opti.bounded(0.0, tau_f_vec[j], 1.0))
    for i in range(n_p):
        opti.subject_to(opti.bounded(T_MIN, T_vec[i], T_MAX))

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
            D = (D_SMOOTH or {}).get(j, {}).get(name)
            D_log = (D_SMOOTH_LOG or {}).get(j, {}).get(name)
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
    #      Qdot_f_in (t, lam) =  alpha_{lam, f_i} Qdot_in(t, lam)
    #      Qdot_p_in (t, lam) =  tau_{lam, f_i} alpha_{lam, therm, i}(T) Qdot_in(t, lam)
    #      Qdot_out  (t, lam) = -eps_{lam, i} tau_{lam, f_i} M_lam(lam, T_i) dlam A_i
    #      Qdot_i(t)          =  sum_lam ( Qdot_f_in + Qdot_p_in + Qdot_out )
    #
    #      with  alpha_{lam, therm, i}(T) = alpha_{lam, i} - alpha_{lam, elec, i}(T)
    #      and   P_elec_i(t)   =  sum_lam tau_{lam, f_i} alpha_{lam, elec, i}(T) Qdot_in
    # -------------------------------------------------------------------------
    Qdot = {}       # Qdot[i][n]   net thermal power into panel i at time n  [W]
    P_elec = {}     # P_elec[i][n] electrical power off panel i at time n    [W]

    for i in range(n_p):
        M = int(data.material[i])
        Qdot[i], P_elec[i] = [], []

        for n in range(n_t):
            Ti = T[i][n]
            M_dlam = ca.vertsplit(planck_band_power(Ti, lam_um, dlam_um))  # [W m^-2]
            g = g_elec(Ti, M)                                   # alpha_elec temp. factor

            Q_terms, P_terms = [], []
            for k in range(n_lam):
                Q_in_k = data.Q_in[i, n, k]

                a_elec = alpha_el[M][k] * g                 # alpha_{lam, elec}(T)
                a_therm = alpha_p[M][k] - a_elec            # alpha_{lam, therm}(T)

                Q_f_in = alpha_f[M][k] * Q_in_k
                Q_p_in = tau_f[M][k] * a_therm * Q_in_k
                Q_out = -eps_p[M][k] * tau_f[M][k] * M_dlam[k] * data.area[i]

                Q_terms.append(Q_f_in + Q_p_in + Q_out)
                P_terms.append(tau_f[M][k] * a_elec * Q_in_k)

            # sum over all lambda
            Qdot[i].append(ca.sum1(ca.vertcat(*Q_terms)))
            P_elec[i].append(ca.sum1(ca.vertcat(*P_terms)))

    # -------------------------------------------------------------------------
    # 5.5  Temperature dynamics, eqs. (2) and (3)
    #          dQ_i(t)  = Qdot_i(t) dt
    #          T_i(t + dt) = T_i(t) + dQ_i(t) / (m c)
    #      with T_i(0) = T_{i,0}.
    # -------------------------------------------------------------------------
    for i in range(n_p):
        opti.subject_to(T[i][0] == data.T0[i])
        mc = data.mass[i] * data.heat_cap[i]
        for n in range(n_t - 1):
            dQ = Qdot[i][n] * data.dt[n]
            opti.subject_to(T[i][n + 1] == T[i][n] + dQ / mc)

    # -------------------------------------------------------------------------
    # 5.6  Cost function, eq. (4)
    #
    #          J = -W_elec * E_elec / E_ref
    #              + W_temp   / (n_p n_t) * sum_i sum_t ( (T_i(t) - T_set) / dT_ref )^2
    #              + W_smooth / (n_j n_lam) * sum_j sum_lam sum_p ( p_j(l+dl) - p_j(l) )^2
    #
    #      Term 1 maximises the electrical energy harvested over the dataset,
    #      term 2 holds the panels at the set-point.  Both are O(1), so the
    #      weights trade off directly.  E_ref is the total incident energy, so
    #      term 1 is minus the overall conversion efficiency.
    #
    #      Term 3 is a tie-breaker, not physics.  In the thermal IR the cost
    #      depends only on the Planck-weighted average of eps*tau, so a whole
    #      family of filters is exactly optimal and IPOPT returns an arbitrary
    #      (ragged) member of it.  A small roughness penalty selects the
    #      smoothest member.  Keep W_SMOOTH small enough that it does not move
    #      E_elec (check: it should change the efficiency in the 4th decimal);
    #      set W_SMOOTH = 0 for the pure two-term cost.
    # -------------------------------------------------------------------------
    E_terms, T_terms = [], []
    for i in range(n_p):
        for n in range(n_t):
            E_terms.append(P_elec[i][n] * data.dt[n])
            T_terms.append(((T[i][n] - T_SET) / T_REF_SCALE) ** 2)

    S_terms = []
    for j in filter_types:
        for p in (alpha_f[j], rho_f[j], tau_f[j]):
            for k in range(n_lam - 1):
                S_terms.append((p[k + 1] - p[k]) ** 2)

    E_elec = ca.sum1(ca.vertcat(*E_terms))
    T_track = ca.sum1(ca.vertcat(*T_terms)) / (n_p * n_t)
    roughness = ca.sum1(ca.vertcat(*S_terms)) / (len(filter_types) * n_lam)

    E_ref = float(np.sum(data.Q_in * data.dt[None, :, None]))
    J = -W_ELEC * E_elec / E_ref + W_TEMP * T_track + W_SMOOTH * roughness

    opti.minimize(J)
    opti.solver("ipopt", ipopt_options if ipopt_options is not None else IPOPT_OPTIONS_QUIET)

    if verbose:
        frac = in_band_emission_fraction(T_SET, data.lam_nm, data.dlam_nm)
        print(f"  panels           : {n_p}  {data.labels}")
        print(f"  material types   : {[int(m) for m in data.material]} "
              f"-> filters for {filter_types}")
        print(f"  areas       [m^2]: {[round(float(a), 5) for a in data.area]}")
        print(f"  time samples     : {n_t}   span {data.t[0]:.2f} .. {data.t[-1]:.2f} s")
        print(f"  wavelength bands : {n_lam}  ({data.lam_nm[0]:.1f} .. {data.lam_nm[-1]:.0f} nm)")
        print(f"  incident energy  : E_ref = {E_ref:.4f} J")
        print(f"  in-band emission : {100 * frac:.2f} % of sigma T^4 at T_set "
              f"(the rest lies beyond {data.lam_edges_nm[-1] / 1000:.0f} um)")

    return FilterProblem(opti=opti, data=data, filter_types=filter_types,
                         alpha_f=alpha_f_vec, rho_f=rho_f_vec, tau_f=tau_f_vec,
                         T=T_vec, J=J, E_elec=E_elec, E_ref=E_ref, T_track=T_track)


# =============================================================================
# 6.  INITIAL GUESSES
#     -------------------------------------------------------------------
#     Every guess satisfies alpha + rho + tau = 1 and the bounds, which keeps
#     IPOPT out of trouble.  `rng=None` gives the deterministic nominal start.
# =============================================================================

def make_initial_guess(prob: FilterProblem, rng: np.random.Generator | None = None) -> dict:
    """Return {"filters": {j: (alpha, rho, tau)}, "T": {i: trajectory}}."""
    data = prob.data
    n_lam, n_t = data.n_lam, data.n_t

    guess = {"filters": {}, "T": {}}
    for j in prob.filter_types:
        if rng is None:
            # Nominal: a clear filter (tau = 1) -- neutral, feasible, unbiased.
            alpha = np.full(n_lam, 0.02)
            tau = np.full(n_lam, 0.96)
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
        guess["T"][i] = np.full(n_t, data.T0[i] + spread)
        guess["T"][i][0] = data.T0[i]
    return guess


def apply_initial_guess(prob: FilterProblem, guess: dict) -> None:
    for j in prob.filter_types:
        alpha, rho, tau = guess["filters"][j]
        prob.opti.set_initial(prob.alpha_f[j], alpha)
        prob.opti.set_initial(prob.rho_f[j], rho)
        prob.opti.set_initial(prob.tau_f[j], tau)
    for i, traj in guess["T"].items():
        prob.opti.set_initial(prob.T[i], traj)


# =============================================================================
# 7.  SOLVE / MULTI-START DRIVER
# =============================================================================

@dataclass
class Solution:
    success: bool
    J: float
    E_elec: float
    efficiency: float                  # E_elec / E_ref
    T_track: float
    filters: dict[int, dict[str, np.ndarray]]
    T: dict[str, np.ndarray]
    start: int = 0
    message: str = ""


def _extract(prob: FilterProblem, val: Callable, ok: bool, start: int, msg: str) -> Solution:
    filters = {j: {"alpha": np.asarray(val(prob.alpha_f[j])).ravel(),
                   "rho": np.asarray(val(prob.rho_f[j])).ravel(),
                   "tau": np.asarray(val(prob.tau_f[j])).ravel(),
                   "eps": np.asarray(val(prob.alpha_f[j])).ravel()}   # Kirchhoff
               for j in prob.filter_types}
    T = {lab: np.asarray(val(prob.T[i])).ravel()
         for i, lab in enumerate(prob.data.labels)}
    E = float(val(prob.E_elec))
    return Solution(success=ok, J=float(val(prob.J)), E_elec=E,
                    efficiency=E / prob.E_ref, T_track=float(val(prob.T_track)),
                    filters=filters, T=T, start=start, message=msg)


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


def multi_start(prob: FilterProblem, n_starts: int = N_STARTS, seed: int = RANDOM_SEED
                ) -> tuple[Solution, list[Solution]]:
    """Solve from `n_starts` initial guesses and keep the best feasible one.

    Start 0 is the deterministic nominal guess; the rest are randomised by
    `make_initial_guess`.  Only the first solve pays the solver-construction
    cost, so additional starts are comparatively cheap.
    """
    rng = np.random.default_rng(seed)
    results: list[Solution] = []
    for s in range(n_starts):
        guess = make_initial_guess(prob, None if s == 0 else rng)
        t0 = time.perf_counter()
        res = solve_once(prob, guess, start=s)
        results.append(res)
        print(f"  start {s:>3d}: {'ok ' if res.success else 'FAIL'}  "
              f"J = {res.J: .6f}   eta = {100 * res.efficiency:6.3f} %   "
              f"T-track = {res.T_track:.4f}   [{time.perf_counter() - t0:5.1f} s]",
              flush=True)

    feasible = [r for r in results if r.success] or results
    best = min(feasible, key=lambda r: r.J)
    return best, results


# =============================================================================
# 8.  RESULTS
# =============================================================================

MATERIAL_NAMES = {SOLAR: "solar-panel filter", ALUMINIUM: "aluminium filter"}


def report(best: Solution, prob: FilterProblem) -> None:
    data = prob.data
    print("\n" + "=" * 70)
    print(f"BEST SOLUTION  (start {best.start}, {best.message})")
    print("=" * 70)
    print(f"  J                      {best.J: .6f}")
    print(f"  electrical energy      {best.E_elec:.4f} J  "
          f"over {data.t[-1] - data.t[0]:.1f} s")
    print(f"  conversion efficiency  {100 * best.efficiency:.3f} %  of incident energy")
    print(f"  mean (dT/dT_ref)^2     {best.T_track:.5f}")

    print("\n  panel   type   A [m^2]   T start -> end [K]    max |T - T_set|")
    for i, (lab, traj) in enumerate(best.T.items()):
        m = int(data.material[i])
        print(f"    {lab:>3s}      {m}     {data.area[i]:.4f}   "
              f"{traj[0]:8.3f} -> {traj[-1]:8.3f}      {np.max(np.abs(traj - T_SET)):8.3f}")

    print("\n  filter band averages")
    for j, f in best.filters.items():
        print(f"    {MATERIAL_NAMES.get(j, f'filter {j}')} (M = {j}):")
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


def save(best: Solution, results: list[Solution], prob: FilterProblem, path: Path) -> None:
    payload = {
        "lam_nm": prob.data.lam_nm,
        "dlam_nm": prob.data.dlam_nm,
        "lam_edges_nm": prob.data.lam_edges_nm,
        "t": prob.data.t,
        "material": prob.data.material,
        "area": prob.data.area,
        "filters": best.filters,
        "T": best.T,
        "J": best.J,
        "E_elec": best.E_elec,
        "efficiency": best.efficiency,
        "all_starts": [(r.start, r.success, r.J, r.efficiency) for r in results],
        "settings": {"T_SET": T_SET, "W_ELEC": W_ELEC, "W_TEMP": W_TEMP,
                     "W_SMOOTH": W_SMOOTH, "D_SMOOTH": D_SMOOTH,
                     "D_SMOOTH_LOG": D_SMOOTH_LOG,
                     "N_FILTER_BANDS": N_FILTER_BANDS, "BAND_SPACING": BAND_SPACING},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(payload, path)
    print(f"\n  saved -> {path}")


def plot(best: Solution, prob: FilterProblem, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = prob.data
    n = len(best.filters)
    fig, axes = plt.subplots(1, n + 1, figsize=(5.8 * (n + 1), 4.2))
    axes = np.atleast_1d(axes)

    for ax, (j, f) in zip(axes, best.filters.items()):
        for key, style in (("alpha", "-"), ("rho", "--"), ("tau", "-.")):
            ax.plot(data.lam_nm, f[key], style, label=rf"$\{key}_{{{j}}}$")
        ax.set_xscale("log")
        ax.set_xlabel("wavelength [nm]")
        ax.set_ylabel("property")
        ax.set_title(f"{MATERIAL_NAMES.get(j, f'filter {j}')}  (M = {j})")
        ax.set_ylim(-0.05, 1.05)
        ax.grid(alpha=0.3)
        ax.legend()

    ax = axes[-1]
    for i, (lab, traj) in enumerate(best.T.items()):
        ax.plot(data.t, traj, ls="-" if data.material[i] == SOLAR else "--",
                label=f"{lab} (M={int(data.material[i])})")
    ax.axhline(T_SET, color="k", ls=":", label=r"$T_\mathrm{set}$")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("panel temperature [K]")
    ax.set_title("Temperature trajectories")
    ax.grid(alpha=0.3)
    ax.legend(ncol=2, fontsize=8)

    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"  saved -> {path}")


# =============================================================================
#  MAIN  -- behaviour is controlled entirely by Section 1.0
# =============================================================================

def main() -> None:
    print(f"Loading {PICKLE_FILE} ...")
    data = load_sim_data()

    print("Building NLP ...")
    t0 = time.perf_counter()
    prob = build_problem(data, ipopt_options=IPOPT_OPTIONS if VERBOSE_IPOPT
                         else IPOPT_OPTIONS_QUIET)
    print(f"  built in {time.perf_counter() - t0:.1f} s")

    print(f"Solving ({N_STARTS} start{'s' if N_STARTS > 1 else ''}) ... "
          f"the first solve also builds the derivative functions")
    best, results = multi_start(prob, n_starts=N_STARTS, seed=RANDOM_SEED)

    report(best, prob)
    if SAVE_SOLUTION:
        save(best, results, prob, OUTPUT_FILE)
    if MAKE_PLOT:
        plot(best, prob, PLOT_FILE)


if __name__ == "__main__":
    main()
