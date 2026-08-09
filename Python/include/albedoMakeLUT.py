#!/usr/bin/env python3

"""
Generate a libRadtran top-of-atmosphere Earth-reflected spectral-radiance LUT.

The LUT is orbit-independent.  At runtime, the C/C++ simulator computes the local
geometry of each visible Earth patch:

    solar zenith angle
    viewing zenith angle
    relative Sun/view azimuth
    local Lambertian surface albedo

and interpolates this LUT.

The data stored are TOA spectral radiance:

    W / (m^2 nm sr)

Binary data order:

    radiance[
        albedo_index,
        solar_zenith_index,
        view_zenith_index,
        relative_azimuth_index,
        wavelength_index
    ]

All axes and their sample counts are stored in the binary file before the data.
"""

from pathlib import Path
import math
import shutil
import struct
import subprocess
import tempfile

import numpy as np


# =============================================================================
# USER CONFIGURATION
# =============================================================================

# ---- Sampling ----------------------------------------------------------------

# Number of samples used for EACH angular coordinate:
#   solar zenith,
#   viewing zenith,
#   relative azimuth.
#
# Therefore the angular part of the LUT contains N_POSITION_SAMPLES^3 points
# for each albedo and wavelength.
N_POSITION_SAMPLES = 20

# Number of output wavelengths in the LUT.
N_WAVELENGTH_SAMPLES = 1000

# Surface albedo can also be interpolated.
# Set this to 1 for one fixed global albedo.
N_ALBEDO_SAMPLES = 1


# ---- Angular domain ----------------------------------------------------------

SOLAR_ZENITH_MIN_DEG = 0.0
SOLAR_ZENITH_MAX_DEG = 89.9

VIEW_ZENITH_MIN_DEG = 0.0
VIEW_ZENITH_MAX_DEG = 89.9

# For libRadtran's one-dimensional atmosphere only relative azimuth matters.
# 0..180 degrees spans all unique relative geometries.
RELATIVE_AZIMUTH_MIN_DEG = 0.0
RELATIVE_AZIMUTH_MAX_DEG = 180.0


# ---- Wavelength domain -------------------------------------------------------

WAVELENGTH_MIN_NM = 300.0
WAVELENGTH_MAX_NM = 2500.0


# ---- Surface albedo domain ---------------------------------------------------

# If N_ALBEDO_SAMPLES == 1, FIXED_SURFACE_ALBEDO is used.
FIXED_SURFACE_ALBEDO = 0.30

# These are used only when N_ALBEDO_SAMPLES > 1.
SURFACE_ALBEDO_MIN = 0.0
SURFACE_ALBEDO_MAX = 1.0


# ---- libRadtran settings -----------------------------------------------------

NUMBER_OF_STREAMS = 16

# REPTRAN is libRadtran's default spectral parameterization. "fine" requests
# its finer spectral representation.
MOLECULAR_ABSORPTION = "reptran fine"

# Use the full-resolution Kurucz spectrum already present in this project.
USE_PROJECT_KURUCZ_SPECTRUM = True

# libRadtran's solar spectral units are mW/(m^2 nm), and "uu" therefore has
# units mW/(m^2 nm sr). Store SI units in the LUT.
LIBRADTRAN_RADIANCE_TO_SI = 1.0e-3


# ---- Output ------------------------------------------------------------------

OUTPUT_FILENAME = "albedo_lut.bin"


# =============================================================================
# PROJECT PATHS
# =============================================================================

# This file is expected at:
#   <project>/Python/include/albedoMakeLUT.py
PROJECT_ROOT = Path(__file__).resolve().parents[2]

LIBRADTRAN_DIR = PROJECT_ROOT / "C" / "lib" / "libRadtran-2.0.6"
UVSPEC_PATH = LIBRADTRAN_DIR / "bin" / "uvspec"

LIBRADTRAN_DATA_DIR = LIBRADTRAN_DIR / "data"
ATMOSPHERE_FILE = LIBRADTRAN_DATA_DIR / "atmmod" / "afglus.dat"

PROJECT_KURUCZ_FILE = PROJECT_ROOT / "C" / "kurudz_full.dat"

OUTPUT_PATH = PROJECT_ROOT / "C" / OUTPUT_FILENAME


# =============================================================================
# LUT BINARY FORMAT
# =============================================================================

LUT_MAGIC = b"ALBLUT02"
LUT_VERSION = 2

# Fixed-size metadata header.  Axis arrays immediately follow this header.
HEADER_SIZE_BYTES = 256

# Data are stored as little-endian IEEE-754 float64.
DATA_DTYPE = np.dtype("<f8")
DATA_SCALAR_BYTES = DATA_DTYPE.itemsize


# Header field layout, in order:
#
#   byte  type      name
#   ----  --------  ----------------------------------------------------------
#      0  char[8]   magic = "ALBLUT02"
#      8  uint32    version
#     12  uint32    header_size_bytes
#     16  uint32    scalar_bytes
#     20  uint32    n_position_samples
#     24  uint32    n_wavelength_samples
#     28  uint32    n_albedo_samples
#     32  uint32    n_solar_zenith_samples
#     36  uint32    n_view_zenith_samples
#     40  uint32    n_relative_azimuth_samples
#     44  uint32    number_of_streams
#     48  uint64    axis_data_offset_bytes
#     56  uint64    radiance_data_offset_bytes
#     64  uint64    radiance_value_count
#     72  float64   wavelength_min_nm
#     80  float64   wavelength_max_nm
#     88  float64   solar_zenith_min_deg
#     96  float64   solar_zenith_max_deg
#    104  float64   view_zenith_min_deg
#    112  float64   view_zenith_max_deg
#    120  float64   relative_azimuth_min_deg
#    128  float64   relative_azimuth_max_deg
#    136  float64   surface_albedo_min
#    144  float64   surface_albedo_max
#    152  float64   fixed_surface_albedo
#
# Remaining bytes up to HEADER_SIZE_BYTES are reserved and set to zero.
#
# Following the header are the float64 axis arrays, in this exact order:
#
#   surface_albedo[n_albedo_samples]
#   solar_zenith_deg[n_solar_zenith_samples]
#   view_zenith_deg[n_view_zenith_samples]
#   relative_azimuth_deg[n_relative_azimuth_samples]
#   wavelength_nm[n_wavelength_samples]
#
# Following the axis arrays is the float64 radiance array:
#
#   [albedo][solar_zenith][view_zenith][relative_azimuth][wavelength]


# =============================================================================
# VALIDATION AND GRIDS
# =============================================================================

def validate_configuration() -> None:
    if N_POSITION_SAMPLES < 2:
        raise ValueError("N_POSITION_SAMPLES must be at least 2.")

    if N_WAVELENGTH_SAMPLES < 2:
        raise ValueError("N_WAVELENGTH_SAMPLES must be at least 2.")

    if N_ALBEDO_SAMPLES < 1:
        raise ValueError("N_ALBEDO_SAMPLES must be at least 1.")

    if not (0.0 <= SOLAR_ZENITH_MIN_DEG < SOLAR_ZENITH_MAX_DEG < 90.0):
        raise ValueError(
            "Solar zenith limits must satisfy "
            "0 <= min < max < 90 degrees."
        )

    if not (0.0 <= VIEW_ZENITH_MIN_DEG < VIEW_ZENITH_MAX_DEG < 90.0):
        raise ValueError(
            "Viewing zenith limits must satisfy "
            "0 <= min < max < 90 degrees."
        )

    if not (
        0.0 <= RELATIVE_AZIMUTH_MIN_DEG
        < RELATIVE_AZIMUTH_MAX_DEG
        <= 180.0
    ):
        raise ValueError(
            "Relative azimuth limits must satisfy "
            "0 <= min < max <= 180 degrees."
        )

    if WAVELENGTH_MIN_NM >= WAVELENGTH_MAX_NM:
        raise ValueError("WAVELENGTH_MIN_NM must be less than WAVELENGTH_MAX_NM.")

    if NUMBER_OF_STREAMS <= 0 or NUMBER_OF_STREAMS % 2 != 0:
        raise ValueError(
            "NUMBER_OF_STREAMS must be a positive even integer for DISORT."
        )

    if not (0.0 <= FIXED_SURFACE_ALBEDO <= 1.0):
        raise ValueError("FIXED_SURFACE_ALBEDO must lie in [0, 1].")

    if not (0.0 <= SURFACE_ALBEDO_MIN <= SURFACE_ALBEDO_MAX <= 1.0):
        raise ValueError("Surface albedo limits must lie in [0, 1].")


def make_grids() -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    solar_zenith_deg = np.linspace(
        SOLAR_ZENITH_MIN_DEG,
        SOLAR_ZENITH_MAX_DEG,
        N_POSITION_SAMPLES,
        dtype=np.float64,
    )

    view_zenith_deg = np.linspace(
        VIEW_ZENITH_MIN_DEG,
        VIEW_ZENITH_MAX_DEG,
        N_POSITION_SAMPLES,
        dtype=np.float64,
    )

    relative_azimuth_deg = np.linspace(
        RELATIVE_AZIMUTH_MIN_DEG,
        RELATIVE_AZIMUTH_MAX_DEG,
        N_POSITION_SAMPLES,
        dtype=np.float64,
    )

    wavelength_nm = np.linspace(
        WAVELENGTH_MIN_NM,
        WAVELENGTH_MAX_NM,
        N_WAVELENGTH_SAMPLES,
        dtype=np.float64,
    )

    if N_ALBEDO_SAMPLES == 1:
        surface_albedo = np.array(
            [FIXED_SURFACE_ALBEDO],
            dtype=np.float64,
        )
    else:
        surface_albedo = np.linspace(
            SURFACE_ALBEDO_MIN,
            SURFACE_ALBEDO_MAX,
            N_ALBEDO_SAMPLES,
            dtype=np.float64,
        )

    return (
        surface_albedo,
        solar_zenith_deg,
        view_zenith_deg,
        relative_azimuth_deg,
        wavelength_nm,
    )


# =============================================================================
# libRadtran
# =============================================================================

def find_uvspec() -> Path:
    if UVSPEC_PATH.is_file():
        return UVSPEC_PATH

    system_uvspec = shutil.which("uvspec")

    if system_uvspec is not None:
        return Path(system_uvspec)

    raise RuntimeError(
        "Could not find libRadtran uvspec.\n"
        f"Expected:\n  {UVSPEC_PATH}\n"
        "Build libRadtran first with ./configure and make."
    )


def write_wavelength_file(
    path: Path,
    wavelength_nm: np.ndarray,
) -> None:
    np.savetxt(
        path,
        wavelength_nm,
        fmt="%.15f",
    )


def build_uvspec_input(
    wavelength_file: Path,
    solar_zenith_deg: float,
    surface_albedo: float,
    view_zenith_deg: np.ndarray,
    relative_azimuth_deg: np.ndarray,
) -> str:
    # libRadtran requires umu values in INCREASING order.
    #
    # Our natural view grid is 0 deg -> horizon, whose cosine goes
    # 1 -> 0. Therefore we reverse the cosine array before passing it
    # to libRadtran. run_uvspec() reverses the view dimension back after
    # parsing the result.
    view_mu = np.cos(np.deg2rad(view_zenith_deg))
    libradtran_umu = view_mu[::-1]

    umu_text = " ".join(
        f"{float(value):.15f}"
        for value in libradtran_umu
    )

    phi_text = " ".join(
        f"{float(value):.15f}"
        for value in relative_azimuth_deg
    )

    source_line = "source solar"

    if USE_PROJECT_KURUCZ_SPECTRUM:
        source_line = f"source solar {PROJECT_KURUCZ_FILE}"

    return f"""\
data_files_path {LIBRADTRAN_DATA_DIR}/
atmosphere_file {ATMOSPHERE_FILE}

{source_line}
mol_abs_param {MOLECULAR_ABSORPTION}

wavelength {WAVELENGTH_MIN_NM:.15f} {WAVELENGTH_MAX_NM:.15f}
spline_file {wavelength_file}
output_process per_nm

sza {solar_zenith_deg:.15f}
phi0 0.0

umu {umu_text}
phi {phi_text}

zout TOA

rte_solver disort
number_of_streams {NUMBER_OF_STREAMS}

albedo {surface_albedo:.15f}

output_user lambda uu
quiet
"""


def run_uvspec(
    uvspec_path: Path,
    input_text: str,
    n_view: int,
    n_azimuth: int,
    n_wavelength: int,
) -> np.ndarray:
    process = subprocess.run(
        [str(uvspec_path)],
        input=input_text,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )

    if process.returncode != 0:
        raise RuntimeError(
            "uvspec failed.\n\n"
            f"Input:\n{input_text}\n\n"
            f"stderr:\n{process.stderr}"
        )

    numeric_output = np.fromstring(
        process.stdout,
        sep=" ",
        dtype=np.float64,
    )

    values_per_wavelength = 1 + n_view * n_azimuth
    expected_value_count = n_wavelength * values_per_wavelength

    if numeric_output.size != expected_value_count:
        raise RuntimeError(
            "Unexpected uvspec output size.\n"
            f"Expected {expected_value_count} numeric values but got "
            f"{numeric_output.size}.\n"
            f"Expected {n_wavelength} wavelength rows with "
            f"{values_per_wavelength} values per row.\n\n"
            f"First 2000 characters of stdout:\n{process.stdout[:2000]}"
        )

    rows = numeric_output.reshape(
        n_wavelength,
        values_per_wavelength,
    )

    returned_wavelengths = rows[:, 0]

    if returned_wavelengths.size != n_wavelength:
        raise RuntimeError("Unexpected number of returned wavelengths.")

    radiance = rows[:, 1:].reshape(
        n_wavelength,
        n_view,
        n_azimuth,
    )

    # libRadtran returned view directions in increasing umu order, which
    # corresponds to view zenith from near-horizon -> nadir.
    # Reverse that axis so stored LUT order matches our grid:
    # view zenith 0 deg -> horizon.
    radiance = radiance[:, ::-1, :]

    # Store as:
    # [view][azimuth][wavelength]
    radiance = np.transpose(
        radiance,
        (1, 2, 0),
    )

    return radiance * LIBRADTRAN_RADIANCE_TO_SI


# =============================================================================
# BINARY LUT CREATION
# =============================================================================

def calculate_axis_bytes(
    n_albedo: int,
    n_solar: int,
    n_view: int,
    n_azimuth: int,
    n_wavelength: int,
) -> int:
    axis_value_count = (
        n_albedo
        + n_solar
        + n_view
        + n_azimuth
        + n_wavelength
    )

    return axis_value_count * DATA_SCALAR_BYTES


def calculate_radiance_value_count(
    n_albedo: int,
    n_solar: int,
    n_view: int,
    n_azimuth: int,
    n_wavelength: int,
) -> int:
    return (
        n_albedo
        * n_solar
        * n_view
        * n_azimuth
        * n_wavelength
    )


def write_binary_header_and_axes(
    output_path: Path,
    surface_albedo: np.ndarray,
    solar_zenith_deg: np.ndarray,
    view_zenith_deg: np.ndarray,
    relative_azimuth_deg: np.ndarray,
    wavelength_nm: np.ndarray,
) -> tuple[int, int]:
    n_albedo = surface_albedo.size
    n_solar = solar_zenith_deg.size
    n_view = view_zenith_deg.size
    n_azimuth = relative_azimuth_deg.size
    n_wavelength = wavelength_nm.size

    axis_data_offset = HEADER_SIZE_BYTES

    axis_bytes = calculate_axis_bytes(
        n_albedo,
        n_solar,
        n_view,
        n_azimuth,
        n_wavelength,
    )

    radiance_data_offset = axis_data_offset + axis_bytes

    radiance_value_count = calculate_radiance_value_count(
        n_albedo,
        n_solar,
        n_view,
        n_azimuth,
        n_wavelength,
    )

    radiance_bytes = radiance_value_count * DATA_SCALAR_BYTES
    final_file_size = radiance_data_offset + radiance_bytes

    header = bytearray(HEADER_SIZE_BYTES)

    struct.pack_into("<8s", header, 0, LUT_MAGIC)
    struct.pack_into("<I", header, 8, LUT_VERSION)
    struct.pack_into("<I", header, 12, HEADER_SIZE_BYTES)
    struct.pack_into("<I", header, 16, DATA_SCALAR_BYTES)

    struct.pack_into("<I", header, 20, N_POSITION_SAMPLES)
    struct.pack_into("<I", header, 24, N_WAVELENGTH_SAMPLES)

    struct.pack_into("<I", header, 28, n_albedo)
    struct.pack_into("<I", header, 32, n_solar)
    struct.pack_into("<I", header, 36, n_view)
    struct.pack_into("<I", header, 40, n_azimuth)
    struct.pack_into("<I", header, 44, NUMBER_OF_STREAMS)

    struct.pack_into("<Q", header, 48, axis_data_offset)
    struct.pack_into("<Q", header, 56, radiance_data_offset)
    struct.pack_into("<Q", header, 64, radiance_value_count)

    struct.pack_into("<d", header, 72, WAVELENGTH_MIN_NM)
    struct.pack_into("<d", header, 80, WAVELENGTH_MAX_NM)

    struct.pack_into("<d", header, 88, SOLAR_ZENITH_MIN_DEG)
    struct.pack_into("<d", header, 96, SOLAR_ZENITH_MAX_DEG)

    struct.pack_into("<d", header, 104, VIEW_ZENITH_MIN_DEG)
    struct.pack_into("<d", header, 112, VIEW_ZENITH_MAX_DEG)

    struct.pack_into("<d", header, 120, RELATIVE_AZIMUTH_MIN_DEG)
    struct.pack_into("<d", header, 128, RELATIVE_AZIMUTH_MAX_DEG)

    struct.pack_into("<d", header, 136, float(surface_albedo[0]))
    struct.pack_into("<d", header, 144, float(surface_albedo[-1]))
    struct.pack_into("<d", header, 152, FIXED_SURFACE_ALBEDO)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open("wb") as file:
        file.write(header)

        np.asarray(
            surface_albedo,
            dtype=DATA_DTYPE,
        ).tofile(file)

        np.asarray(
            solar_zenith_deg,
            dtype=DATA_DTYPE,
        ).tofile(file)

        np.asarray(
            view_zenith_deg,
            dtype=DATA_DTYPE,
        ).tofile(file)

        np.asarray(
            relative_azimuth_deg,
            dtype=DATA_DTYPE,
        ).tofile(file)

        np.asarray(
            wavelength_nm,
            dtype=DATA_DTYPE,
        ).tofile(file)

        file.truncate(final_file_size)

    return radiance_data_offset, radiance_value_count


# =============================================================================
# GENERATION
# =============================================================================

def generate_lut() -> None:
    validate_configuration()

    uvspec_path = find_uvspec()

    (
        surface_albedo,
        solar_zenith_deg,
        view_zenith_deg,
        relative_azimuth_deg,
        wavelength_nm,
    ) = make_grids()

    n_albedo = surface_albedo.size
    n_solar = solar_zenith_deg.size
    n_view = view_zenith_deg.size
    n_azimuth = relative_azimuth_deg.size
    n_wavelength = wavelength_nm.size

    radiance_value_count = calculate_radiance_value_count(
        n_albedo,
        n_solar,
        n_view,
        n_azimuth,
        n_wavelength,
    )

    radiance_bytes = radiance_value_count * DATA_SCALAR_BYTES
    approximate_gib = radiance_bytes / (1024.0 ** 3)

    print(f"uvspec: {uvspec_path}")
    print()
    print("Configured LUT:")
    print(f"  position samples per angular axis : {N_POSITION_SAMPLES}")
    print(f"  wavelength samples                : {N_WAVELENGTH_SAMPLES}")
    print(f"  albedo samples                    : {N_ALBEDO_SAMPLES}")
    print()
    print("Data shape:")
    print(
        f"  [{n_albedo}] albedo x "
        f"[{n_solar}] solar zenith x "
        f"[{n_view}] view zenith x "
        f"[{n_azimuth}] relative azimuth x "
        f"[{n_wavelength}] wavelength"
    )
    print()
    print(f"Radiance values: {radiance_value_count:,}")
    print(f"Approx. radiance-data size: {approximate_gib:.3f} GiB")
    print()
    print(f"Output: {OUTPUT_PATH}")
    print()

    radiance_data_offset, stored_value_count = write_binary_header_and_axes(
        OUTPUT_PATH,
        surface_albedo,
        solar_zenith_deg,
        view_zenith_deg,
        relative_azimuth_deg,
        wavelength_nm,
    )

    if stored_value_count != radiance_value_count:
        raise RuntimeError("Internal LUT size mismatch.")

    radiance_lut = np.memmap(
        OUTPUT_PATH,
        mode="r+",
        dtype=DATA_DTYPE,
        offset=radiance_data_offset,
        shape=(
            n_albedo,
            n_solar,
            n_view,
            n_azimuth,
            n_wavelength,
        ),
        order="C",
    )

    total_uvspec_runs = n_albedo * n_solar
    current_run = 0

    with tempfile.TemporaryDirectory(
        prefix="libradtran_albedo_lut_"
    ) as temporary_directory:
        wavelength_file = (
            Path(temporary_directory)
            / "output_wavelengths.dat"
        )

        write_wavelength_file(
            wavelength_file,
            wavelength_nm,
        )

        for albedo_index in range(n_albedo):
            albedo = float(surface_albedo[albedo_index])

            for solar_index in range(n_solar):
                solar_zenith = float(
                    solar_zenith_deg[solar_index]
                )

                current_run += 1

                print(
                    f"[{current_run}/{total_uvspec_runs}] "
                    f"albedo={albedo:.6f}, "
                    f"solar zenith={solar_zenith:.6f} deg"
                )

                input_text = build_uvspec_input(
                    wavelength_file=wavelength_file,
                    solar_zenith_deg=solar_zenith,
                    surface_albedo=albedo,
                    view_zenith_deg=view_zenith_deg,
                    relative_azimuth_deg=relative_azimuth_deg,
                )

                radiance_lut[
                    albedo_index,
                    solar_index,
                    :,
                    :,
                    :,
                ] = run_uvspec(
                    uvspec_path=uvspec_path,
                    input_text=input_text,
                    n_view=n_view,
                    n_azimuth=n_azimuth,
                    n_wavelength=n_wavelength,
                )

                radiance_lut.flush()

    del radiance_lut

    final_size_gib = OUTPUT_PATH.stat().st_size / (1024.0 ** 3)

    print()
    print("Finished.")
    print(f"LUT written to: {OUTPUT_PATH}")
    print(f"Final file size: {final_size_gib:.3f} GiB")
    print()
    print("Stored radiance units:")
    print("  W / (m^2 nm sr)")
    print()
    print("The binary header stores the configured position/wavelength sample")
    print("counts, all individual axis counts and bounds, and the offsets needed")
    print("to read the coordinate arrays and radiance data directly from C/C++.")


if __name__ == "__main__":
    generate_lut()
