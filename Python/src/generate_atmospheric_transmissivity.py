#!/usr/bin/env python3
"""Generate a clear-sky atmospheric transmissivity spectrum with libRadtran.

The script calls libRadtran's ``uvspec`` executable, then writes a CSV with the
wavelength and three dimensionless quantities:

  * direct_beam_transmissivity
  * diffuse_downward_transmissivity
  * total_downward_transmissivity (direct + diffuse)

All three are normalised by the extraterrestrial irradiance projected onto a
horizontal plane, E0*cos(SZA).  Consequently, the direct column is the
atmospheric transmission along the solar ray without an extra cosine factor.

Molecular absorption is selected with MOLECULAR_ABSORPTION.  REPTRAN cannot
cover the whole range on its own: its thermal tables are tagged
``source thermal`` and uvspec refuses to pair them with ``source solar``, while
``source thermal`` provides no direct beam to transmit and disallows
``output_quantity reflectivity``.  LOWTRAN carries no source tag and spans
shortwave and longwave in one calculation, so it is the default.  The hybrid
mode uses REPTRAN in the shortwave, where it is more accurate and its band width
is selectable, and LOWTRAN above the join.  ``source solar`` with no spectrum
file supplies a unit collimated beam, so the result is atmospheric transmission
with no emission by the atmosphere or surface.

AIR_MASS uses the simple plane-parallel relation AM = sec(SZA).  AIR_MASS = 0
is handled specially as no atmosphere (AM0), for which direct = total = 1 and
diffuse = 0.  Otherwise AIR_MASS must be at least 1.

Set GENERATE_DATA to False to skip libRadtran entirely and plot an existing CSV.
Generating data needs only the standard library; plotting needs matplotlib.
Edit the constants immediately below to select the libRadtran path, wavelength
grid, air mass, and atmosphere.  Run with ``--check`` first for a fast
end-to-end smoke test.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from decimal import Decimal
from pathlib import Path


# ---------------------------------------------------------------------------
# User-editable constants
# ---------------------------------------------------------------------------

# True  -> run libRadtran and write OUTPUT_CSV.
# False -> read OUTPUT_CSV and draw OUTPUT_PLOT with matplotlib.  No libRadtran
#          call is made, so the plot reflects whatever produced the CSV.
GENERATE_DATA = False

SCRIPT_DIRECTORY = Path(__file__).resolve().parent

# Expected layout:
#   parent_directory/
#     your_project/generate_atmospheric_transmissivity.py
#     libRadtran-2.0.6/
# Edit this if your libRadtran directory has a different name or location.
LIBRADTRAN_DIRECTORY = (
    SCRIPT_DIRECTORY / "../../C/lib/libRadtran-2.0.6"
).resolve()

OUTPUT_CSV = SCRIPT_DIRECTORY / "earth_atmosphere_transmissivity.csv"
OUTPUT_PLOT = SCRIPT_DIRECTORY / "earth_atmosphere_transmissivity.png"

# LOWTRAN supports roughly 250 nm to just below 100000 nm (100 micrometres).
# The hybrid mode reaches down to 120 nm through the REPTRAN solar table.
# The default below produces 280 nm to 50 um.
MIN_WAVELENGTH_NM = 100.0
MAX_WAVELENGTH_NM = 100000.0 - 1
WAVELENGTH_INTERVAL_NM = 1.0

# 0.0 means AM0 (no atmosphere).  Common values are 1.0 and 1.5.
AIR_MASS = 1.5

# U.S. Standard atmosphere: representative, cloud-free atmosphere profile.
# Other supplied profiles are in libRadtran's data/atmmod directory.
ATMOSPHERE_FILE = Path("data/atmmod/afglus.dat")

# A simple clear blue-sky model.  If False, only molecular absorption and
# Rayleigh scattering are included.  If True, the Shettle rural aerosol model
# is also included.  No cloud is added in either case.
INCLUDE_DEFAULT_AEROSOL = True
AEROSOL_VISIBILITY_KM = 50.0

# Molecular absorption parameterisation.
#
#   "lowtran" - the LOWTRAN/SBDART band model across the whole grid.  One
#               source-agnostic calculation with no join, but its bands are
#               fixed at 20 cm^-1, so REPTRAN_RESOLUTION has no effect.
#   "hybrid"  - REPTRAN up to REPTRAN_SOLAR_THERMAL_SPLIT_NM, LOWTRAN above it.
#               REPTRAN is more accurate in the shortwave and its band width is
#               selectable, at the cost of a step in the spectrum at the join.
MOLECULAR_ABSORPTION = "lowtran"

# REPTRAN band width.  Used only when MOLECULAR_ABSORPTION is "hybrid".
#   "coarse" - 15 cm^-1, included in the normal libRadtran download
#   "medium" -  5 cm^-1, needs the separate REPTRAN archive from libradtran.org
#   "fine"   -  1 cm^-1, needs the separate REPTRAN archive from libradtran.org
REPTRAN_RESOLUTION = "coarse"

# Wavelengths sent to uvspec per call.  Progress is printed after each chunk.
# Set to 0 to run each spectral segment in a single call with no progress
# output.  Chunking is not entirely free: pseudo-spectral parameterisations
# build their internal grid from the wavelengths they are handed, so values
# very near a chunk boundary can differ slightly from a single-call run.
# Larger chunks reduce that; 0 removes it.
PROGRESS_CHUNK_SIZE = 500

RTE_SOLVER = "disort"
NUMBER_OF_STREAMS = 20

# LOWTRAN band-model limits.  The underlying band model resolves 20 cm^-1, which
# is about 5 nm in the visible but widens as lambda^2: roughly 200 nm at 10 um
# and a few micrometres at 50 um.  Sampling more finely than that costs solver
# time without adding spectral information.
LOWTRAN_MIN_WAVELENGTH_NM = 100.0
LOWTRAN_MAX_WAVELENGTH_NM_EXCLUSIVE = 100000.0

# REPTRAN solar table limits, and the wavelength at which the hybrid mode hands
# over to LOWTRAN.  The solar table runs out shortly above the split.
REPTRAN_MIN_WAVELENGTH_NM =100.0
REPTRAN_SOLAR_THERMAL_SPLIT_NM = 5000.0


def validate_absorption_settings() -> None:
    """Reject unusable parameterisation settings before any work is done."""

    if MOLECULAR_ABSORPTION not in {"lowtran", "hybrid"}:
        raise ValueError("MOLECULAR_ABSORPTION must be 'lowtran' or 'hybrid'")
    if MOLECULAR_ABSORPTION == "hybrid" and (
        REPTRAN_RESOLUTION.lower() not in {"coarse", "medium", "fine"}
    ):
        raise ValueError(
            "REPTRAN_RESOLUTION must be 'coarse', 'medium', or 'fine'"
        )


def supported_wavelength_bounds() -> tuple[float, float]:
    """Return the inclusive minimum and exclusive maximum for the mode."""

    if MOLECULAR_ABSORPTION == "hybrid":
        minimum = REPTRAN_MIN_WAVELENGTH_NM
    else:
        minimum = LOWTRAN_MIN_WAVELENGTH_NM
    return minimum, LOWTRAN_MAX_WAVELENGTH_NM_EXCLUSIVE


def decimal_wavelength_grid() -> list[Decimal]:
    """Return an inclusive, exactly equidistant wavelength grid."""

    minimum = Decimal(str(MIN_WAVELENGTH_NM))
    maximum = Decimal(str(MAX_WAVELENGTH_NM))
    interval = Decimal(str(WAVELENGTH_INTERVAL_NM))

    if minimum <= 0:
        raise ValueError("MIN_WAVELENGTH_NM must be greater than zero")
    if maximum < minimum:
        raise ValueError(
            "MAX_WAVELENGTH_NM must be greater than or equal to "
            "MIN_WAVELENGTH_NM"
        )
    if interval <= 0:
        raise ValueError("WAVELENGTH_INTERVAL_NM must be greater than zero")

    minimum_supported_nm, maximum_exclusive_nm = supported_wavelength_bounds()
    minimum_supported = Decimal(str(minimum_supported_nm))
    maximum_exclusive = Decimal(str(maximum_exclusive_nm))
    if minimum < minimum_supported or maximum >= maximum_exclusive:
        raise ValueError(
            f"MOLECULAR_ABSORPTION='{MOLECULAR_ABSORPTION}' supports "
            f"wavelengths from {minimum_supported_nm:g} nm up to, but not "
            f"including, {maximum_exclusive_nm:g} nm."
        )

    step_count = (maximum - minimum) / interval
    integral_step_count = step_count.to_integral_value()
    if step_count != integral_step_count:
        raise ValueError(
            "The wavelength interval must divide "
            "MAX_WAVELENGTH_NM - MIN_WAVELENGTH_NM exactly"
        )

    return [
        minimum + Decimal(index) * interval
        for index in range(int(integral_step_count) + 1)
    ]


def split_wavelengths_by_band_model(
    wavelengths: list[Decimal],
) -> list[tuple[str, list[Decimal]]]:
    """Split an ordered grid into the segments each band model will handle."""

    if not wavelengths:
        return []
    if MOLECULAR_ABSORPTION == "lowtran":
        return [("lowtran", wavelengths)]

    split = Decimal(str(REPTRAN_SOLAR_THERMAL_SPLIT_NM))
    reptran_wavelengths = [
        wavelength for wavelength in wavelengths if wavelength <= split
    ]
    lowtran_wavelengths = [
        wavelength for wavelength in wavelengths if wavelength > split
    ]

    segments: list[tuple[str, list[Decimal]]] = []
    if reptran_wavelengths:
        segments.append(("reptran", reptran_wavelengths))
    if lowtran_wavelengths:
        segments.append(("lowtran", lowtran_wavelengths))
    return segments


def chunk_wavelengths(wavelengths: list[Decimal]) -> list[list[Decimal]]:
    """Split a segment into the batches sent to individual uvspec calls."""

    if not wavelengths:
        return []
    if PROGRESS_CHUNK_SIZE <= 0:
        return [wavelengths]
    return [
        wavelengths[start:start + PROGRESS_CHUNK_SIZE]
        for start in range(0, len(wavelengths), PROGRESS_CHUNK_SIZE)
    ]


def solar_zenith_angle_deg() -> float:
    """Convert plane-parallel air mass to solar zenith angle."""

    if AIR_MASS == 0.0:
        return 0.0
    if AIR_MASS < 1.0:
        raise ValueError("AIR_MASS must be 0 (AM0) or at least 1")

    return math.degrees(math.acos(1.0 / AIR_MASS))


def write_am0_csv(wavelengths: list[Decimal]) -> None:
    """Write the exact no-atmosphere result without calling libRadtran."""

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(
            [
                "wavelength_nm",
                "direct_beam_transmissivity",
                "diffuse_downward_transmissivity",
                "total_downward_transmissivity",
            ]
        )
        for wavelength in wavelengths:
            writer.writerow([format(wavelength, "f"), "1", "0", "1"])


def check_reptran_data_files(wavelengths: list[Decimal]) -> None:
    """Confirm the REPTRAN solar tables exist, if this run will use them."""

    band_models = {
        band_model
        for band_model, _ in split_wavelengths_by_band_model(wavelengths)
    }
    if "reptran" not in band_models:
        return

    resolution = REPTRAN_RESOLUTION.lower()
    reptran_directory = LIBRADTRAN_DIRECTORY / "data/correlated_k/reptran"
    reptran_species = (
        "CH4",
        "CO",
        "CO2",
        "H2O",
        "N2",
        "N2O",
        "O2",
        "O3",
    )
    table_base = f"reptran_solar_{resolution}"
    required_reptran_files = [reptran_directory / f"{table_base}.cdf"]
    required_reptran_files.extend(
        reptran_directory / f"{table_base}.lookup.{species}.cdf"
        for species in reptran_species
    )

    missing_reptran_files = [
        path for path in required_reptran_files if not path.is_file()
    ]
    if missing_reptran_files:
        missing_names = "\n  ".join(
            path.name for path in missing_reptran_files
        )
        archive_hint = (
            "The normal libRadtran distribution includes REPTRAN coarse."
            if resolution == "coarse"
            else (
                "REPTRAN medium/fine lookup tables require the separate "
                "REPTRAN data archive from libradtran.org."
            )
        )
        raise FileNotFoundError(
            f"Missing REPTRAN {resolution} data files:\n  {missing_names}\n"
            f"{archive_hint} Extract the data into "
            f"{LIBRADTRAN_DIRECTORY / 'data'}."
        )


def check_libradtran_installation(
    wavelengths: list[Decimal],
) -> tuple[Path, Path]:
    """Resolve and validate the files needed for this calculation."""

    uvspec_executable = LIBRADTRAN_DIRECTORY / "bin/uvspec"
    atmosphere_path = LIBRADTRAN_DIRECTORY / ATMOSPHERE_FILE

    if not LIBRADTRAN_DIRECTORY.is_dir():
        raise FileNotFoundError(
            f"libRadtran directory not found: {LIBRADTRAN_DIRECTORY}\n"
            "Edit LIBRADTRAN_DIRECTORY near the top of this script."
        )
    if not uvspec_executable.is_file():
        raise FileNotFoundError(
            f"uvspec executable not found: {uvspec_executable}\n"
            "Build libRadtran with ./configure && make first."
        )
    if not os.access(uvspec_executable, os.X_OK):
        raise FileNotFoundError(
            f"uvspec is not executable: {uvspec_executable}\n"
            f"Run: chmod +x '{uvspec_executable}'"
        )
    if not atmosphere_path.is_file():
        raise FileNotFoundError(
            f"Atmosphere profile not found: {atmosphere_path}"
        )

    # Give a useful error before launching a long calculation if a dynamically
    # linked library (for example libgsl.so.27) is missing.  A static build, or
    # a system without ldd, simply skips this optional diagnostic.
    ldd_executable = shutil.which("ldd")
    if ldd_executable is not None:
        dependency_check = subprocess.run(
            [ldd_executable, str(uvspec_executable)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        missing_dependencies = [
            line.strip()
            for line in dependency_check.stdout.splitlines()
            if "not found" in line
        ]
        if missing_dependencies:
            raise RuntimeError(
                "uvspec has missing shared-library dependencies:\n  "
                + "\n  ".join(missing_dependencies)
                + "\nOn Ubuntu, install/reinstall the required runtime and "
                "development libraries with:\n"
                "  sudo apt update\n"
                "  sudo apt install libgsl-dev libnetcdf-dev"
            )

    check_reptran_data_files(wavelengths)
    return uvspec_executable, atmosphere_path


def make_uvspec_input(
    grid_file: Path,
    solar_zenith_deg: float,
    band_model: str,
) -> str:
    """Create a libRadtran uvspec input document."""

    if band_model not in {"lowtran", "reptran"}:
        raise ValueError("band_model must be 'lowtran' or 'reptran'")
    if AEROSOL_VISIBILITY_KM <= 0:
        raise ValueError("AEROSOL_VISIBILITY_KM must be greater than zero")
    if NUMBER_OF_STREAMS <= 0 or NUMBER_OF_STREAMS % 2 != 0:
        raise ValueError("NUMBER_OF_STREAMS must be a positive even integer")

    # output_quantity reflectivity is intentionally used here. In libRadtran,
    # that option normalises irradiances by E0*cos(SZA), whereas the option named
    # transmittance normalises by E0 and retains the geometric cosine factor.
    # Both require `source solar`; neither is available with `source thermal`.
    #
    # `source solar` with no spectrum file gives a unit extraterrestrial beam,
    # so the longwave part of the grid carries no Planck emission.
    lines = [
        "data_files_path data/",
        f"atmosphere_file {ATMOSPHERE_FILE.as_posix()}",
        "source solar",
        # The grid file alone defines the requested range and sampling; a
        # simultaneous `wavelength` directive is rejected.
        f"wavelength_grid_file {grid_file.as_posix()}",
    ]

    if band_model == "reptran":
        # Only the solar table is ever named here. The thermal table carries a
        # `source thermal` attribute that uvspec refuses to accept alongside
        # `source solar`, which is why LOWTRAN handles the longwave instead.
        resolution = REPTRAN_RESOLUTION.lower()
        lines.extend(
            [
                f"mol_abs_param reptran {resolution}",
                "reptran_file data/correlated_k/reptran/"
                f"reptran_solar_{resolution}",
            ]
        )
    else:
        lines.append("mol_abs_param lowtran")

    lines.extend(
        [
            f"sza {solar_zenith_deg:.12g}",
            "albedo 0.0",
            f"rte_solver {RTE_SOLVER}",
            f"number_of_streams {NUMBER_OF_STREAMS}",
            "output_quantity reflectivity",
            "output_user lambda edir edn eglo",
        ]
    )

    if INCLUDE_DEFAULT_AEROSOL:
        lines.extend(
            [
                "aerosol_default",
                f"aerosol_visibility {AEROSOL_VISIBILITY_KM:.12g}",
            ]
        )

    lines.append("quiet")
    return "\n".join(lines) + "\n"


def run_uvspec(
    uvspec_executable: Path,
    uvspec_input: str,
) -> str:
    """Run uvspec and return its standard output."""

    completed = subprocess.run(
        [str(uvspec_executable)],
        input=uvspec_input,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=LIBRADTRAN_DIRECTORY,
        check=False,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "uvspec failed with exit code "
            f"{completed.returncode}.\n\n"
            "Generated uvspec input:\n"
            f"{uvspec_input}\n"
            "uvspec diagnostics:\n"
            f"{completed.stderr.strip()}"
        )

    return completed.stdout


def run_calculation(
    uvspec_executable: Path,
    wavelengths: list[Decimal],
    solar_zenith_deg: float,
    band_model: str,
) -> list[tuple[float, float, float, float]]:
    """Run uvspec for a supplied output grid and validate the result."""

    with tempfile.TemporaryDirectory(
        prefix=f"libradtran_transmissivity_{band_model}_"
    ) as temp_dir:
        grid_file = Path(temp_dir) / "wavelength_grid_nm.dat"
        grid_file.write_text(
            "".join(
                f"{format(wavelength, 'f')}\n" for wavelength in wavelengths
            ),
            encoding="ascii",
        )

        uvspec_input = make_uvspec_input(
            grid_file,
            solar_zenith_deg,
            band_model,
        )
        output_text = run_uvspec(uvspec_executable, uvspec_input)

    return parse_uvspec_output(output_text, wavelengths)


def format_duration(seconds: float) -> str:
    """Format a duration in seconds as H:MM:SS."""

    whole_seconds = int(max(0.0, seconds))
    hours, remainder = divmod(whole_seconds, 3600)
    minutes, remaining_seconds = divmod(remainder, 60)
    return f"{hours:d}:{minutes:02d}:{remaining_seconds:02d}"


def report_progress(
    completed_count: int,
    total_count: int,
    started_at: float,
    band_model: str,
) -> None:
    """Overwrite a single stderr line with the current progress."""

    elapsed = time.monotonic() - started_at
    fraction = completed_count / total_count if total_count else 1.0
    remaining = (
        elapsed * (total_count - completed_count) / completed_count
        if completed_count
        else 0.0
    )
    sys.stderr.write(
        f"\r  {fraction * 100.0:6.2f}%  "
        f"{completed_count}/{total_count} wavelengths  "
        f"[{band_model}]  "
        f"elapsed {format_duration(elapsed)}  "
        f"remaining ~{format_duration(remaining)}    "
    )
    sys.stderr.flush()


def run_all_segments(
    uvspec_executable: Path,
    wavelengths: list[Decimal],
    solar_zenith_deg: float,
) -> list[tuple[float, float, float, float]]:
    """Calculate and join every spectral segment, reporting progress."""

    show_progress = PROGRESS_CHUNK_SIZE > 0
    total_count = len(wavelengths)
    completed_count = 0
    started_at = time.monotonic()

    rows: list[tuple[float, float, float, float]] = []
    for band_model, segment_wavelengths in (
        split_wavelengths_by_band_model(wavelengths)
    ):
        for chunk in chunk_wavelengths(segment_wavelengths):
            rows.extend(
                run_calculation(
                    uvspec_executable,
                    chunk,
                    solar_zenith_deg,
                    band_model,
                )
            )
            completed_count += len(chunk)
            if show_progress:
                report_progress(
                    completed_count,
                    total_count,
                    started_at,
                    band_model,
                )

    if show_progress:
        sys.stderr.write("\n")
        sys.stderr.flush()
    return rows


def clamp_transmissivity(value: float, wavelength_nm: float, name: str) -> float:
    """Clamp tiny numerical overshoots while rejecting nonphysical output."""

    tolerance = 5.0e-4
    if not math.isfinite(value):
        raise RuntimeError(
            f"uvspec returned non-finite {name} at {wavelength_nm:g} nm"
        )
    if value < -tolerance or value > 1.0 + tolerance:
        raise RuntimeError(
            f"uvspec returned {name}={value:g} at {wavelength_nm:g} nm, "
            "outside the expected range [0, 1]"
        )
    return min(1.0, max(0.0, value))


def parse_uvspec_output(
    output_text: str,
    wavelengths: list[Decimal],
) -> list[tuple[float, float, float, float]]:
    """Parse wavelength, direct, diffuse, and global output columns."""

    rows: list[tuple[float, float, float, float]] = []
    for line_number, line in enumerate(output_text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        fields = stripped.split()
        if len(fields) != 4:
            raise RuntimeError(
                "Unexpected uvspec output on line "
                f"{line_number}: expected 4 columns, got {len(fields)}\n{line}"
            )

        try:
            wavelength_nm, direct, diffuse, total = map(float, fields)
        except ValueError as error:
            raise RuntimeError(
                f"Could not parse uvspec output line {line_number}: {line}"
            ) from error

        direct = clamp_transmissivity(
            direct, wavelength_nm, "direct_beam_transmissivity"
        )
        diffuse = clamp_transmissivity(
            diffuse, wavelength_nm, "diffuse_downward_transmissivity"
        )
        total = clamp_transmissivity(
            total, wavelength_nm, "total_downward_transmissivity"
        )
        rows.append((wavelength_nm, direct, diffuse, total))

    if len(rows) != len(wavelengths):
        raise RuntimeError(
            "uvspec returned "
            f"{len(rows)} wavelengths, but {len(wavelengths)} were requested. "
            "Check the requested range against the selected band model."
        )

    wavelength_tolerance_nm = max(
        1.0e-5,
        float(Decimal(str(WAVELENGTH_INTERVAL_NM))) * 1.0e-5,
    )
    for expected, actual_row in zip(wavelengths, rows):
        actual = actual_row[0]
        if abs(actual - float(expected)) > wavelength_tolerance_nm:
            raise RuntimeError(
                "uvspec output grid does not match the requested grid: "
                f"expected {expected} nm, received {actual:g} nm"
            )

    return rows


def write_libradtran_csv(
    rows: list[tuple[float, float, float, float]],
) -> None:
    """Write parsed and validated libRadtran results to CSV."""

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(
            [
                "wavelength_nm",
                "direct_beam_transmissivity",
                "diffuse_downward_transmissivity",
                "total_downward_transmissivity",
            ]
        )
        for wavelength_nm, direct, diffuse, total in rows:
            writer.writerow(
                [
                    f"{wavelength_nm:.10g}",
                    f"{direct:.10g}",
                    f"{diffuse:.10g}",
                    f"{total:.10g}",
                ]
            )


def read_transmissivity_csv() -> dict[str, list[float]]:
    """Read a previously generated CSV into parallel float columns."""

    if not OUTPUT_CSV.is_file():
        raise FileNotFoundError(
            f"No data file to plot: {OUTPUT_CSV}\n"
            "Set GENERATE_DATA = True and run the script once first."
        )

    column_names = [
        "wavelength_nm",
        "direct_beam_transmissivity",
        "diffuse_downward_transmissivity",
        "total_downward_transmissivity",
    ]
    columns: dict[str, list[float]] = {name: [] for name in column_names}

    with OUTPUT_CSV.open(newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        missing_columns = [
            name
            for name in column_names
            if name not in (reader.fieldnames or [])
        ]
        if missing_columns:
            raise RuntimeError(
                f"{OUTPUT_CSV} is missing expected columns: "
                + ", ".join(missing_columns)
            )
        for row_number, row in enumerate(reader, start=2):
            try:
                for name in column_names:
                    columns[name].append(float(row[name]))
            except (TypeError, ValueError) as error:
                raise RuntimeError(
                    f"Could not parse {OUTPUT_CSV} line {row_number}"
                ) from error

    if not columns["wavelength_nm"]:
        raise RuntimeError(f"{OUTPUT_CSV} contains no data rows")
    return columns


def plot_transmissivity() -> None:
    """Draw the spectrum from OUTPUT_CSV and save it to OUTPUT_PLOT."""

    try:
        import matplotlib.pyplot as plt
        from matplotlib.ticker import ScalarFormatter
    except ImportError as error:
        raise RuntimeError(
            "matplotlib is needed to plot. Install it with:\n"
            "  python3 -m pip install matplotlib"
        ) from error

    columns = read_transmissivity_csv()
    wavelengths_um = [value / 1000.0 for value in columns["wavelength_nm"]]

    figure, axes = plt.subplots(figsize=(11.0, 6.0))
    axes.plot(
        wavelengths_um,
        columns["total_downward_transmissivity"],
        linewidth=0.7,
        label="total downward",
    )
    axes.plot(
        wavelengths_um,
        columns["direct_beam_transmissivity"],
        linewidth=0.7,
        label="direct beam",
    )
    axes.plot(
        wavelengths_um,
        columns["diffuse_downward_transmissivity"],
        linewidth=0.7,
        label="diffuse downward",
    )

    axes.set_xscale("log")
    axes.set_xlim(min(wavelengths_um), max(wavelengths_um))
    axes.set_ylim(0.0, 1.02)
    axes.xaxis.set_major_formatter(ScalarFormatter())
    axes.xaxis.set_minor_formatter(ScalarFormatter())
    axes.set_xlabel("Wavelength (\u00b5m)")
    axes.set_ylabel("Transmissivity")
    axes.set_title(
        f"Clear-sky atmospheric transmissivity, AM={AIR_MASS:g}, "
        f"{MOLECULAR_ABSORPTION}"
    )
    axes.grid(True, which="major", linewidth=0.4, alpha=0.5)
    axes.grid(True, which="minor", linewidth=0.2, alpha=0.3)
    axes.legend(loc="upper right")
    figure.tight_layout()

    OUTPUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PLOT, dpi=200)
    print(
        f"Plotted {len(wavelengths_um)} rows from {OUTPUT_CSV}\n"
        f"Wrote {OUTPUT_PLOT}"
    )
    plt.show()


def main(check_only: bool = False) -> None:
    if not GENERATE_DATA:
        if check_only:
            print("--check only applies when GENERATE_DATA is True.")
            return
        plot_transmissivity()
        return

    validate_absorption_settings()
    wavelengths = decimal_wavelength_grid()
    solar_zenith_deg = solar_zenith_angle_deg()

    if AIR_MASS == 0.0:
        if check_only:
            print("Configuration check passed. AM0 does not call uvspec.")
            return
        write_am0_csv(wavelengths)
        print(f"Wrote AM0 spectrum ({len(wavelengths)} rows) to {OUTPUT_CSV}")
        return

    uvspec_executable, _ = check_libradtran_installation(wavelengths)

    if check_only:
        # Exercise every band model this run will use, as well as the
        # executable, atmosphere, solver, and parser, before the full spectrum.
        tested_wavelength_count = 0
        tested_band_models: list[str] = []
        for band_model, segment_wavelengths in (
            split_wavelengths_by_band_model(wavelengths)
        ):
            sample_indexes = sorted(
                {
                    0,
                    len(segment_wavelengths) // 2,
                    len(segment_wavelengths) - 1,
                }
            )
            smoke_test_wavelengths = [
                segment_wavelengths[index] for index in sample_indexes
            ]
            run_calculation(
                uvspec_executable,
                smoke_test_wavelengths,
                solar_zenith_deg,
                band_model,
            )
            tested_wavelength_count += len(smoke_test_wavelengths)
            tested_band_models.append(band_model)

        print(
            "Configuration check passed: uvspec successfully calculated "
            f"{tested_wavelength_count} sample wavelengths using "
            f"{' and '.join(tested_band_models)}."
        )
        return

    print(
        f"Calculating {len(wavelengths)} wavelengths from "
        f"{MIN_WAVELENGTH_NM:g} to {MAX_WAVELENGTH_NM:g} nm..."
    )
    rows = run_all_segments(
        uvspec_executable,
        wavelengths,
        solar_zenith_deg,
    )
    write_libradtran_csv(rows)

    aerosol_description = (
        f"default rural aerosol, {AEROSOL_VISIBILITY_KM:g} km visibility"
        if INCLUDE_DEFAULT_AEROSOL
        else "no aerosol"
    )
    band_model_description = (
        "mol_abs_param lowtran"
        if MOLECULAR_ABSORPTION == "lowtran"
        else (
            f"mol_abs_param reptran {REPTRAN_RESOLUTION.lower()} below "
            f"{REPTRAN_SOLAR_THERMAL_SPLIT_NM:g} nm + lowtran above"
        )
    )
    print(
        f"Wrote {len(rows)} rows to {OUTPUT_CSV}\n"
        f"AM={AIR_MASS:g}, SZA={solar_zenith_deg:.6f} deg, "
        f"{band_model_description}, {aerosol_description}"
    )


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser(
        description=(
            "Generate a clear-sky atmospheric transmissivity spectrum with "
            "libRadtran, or plot a previously generated one."
        )
    )
    argument_parser.add_argument(
        "--check",
        action="store_true",
        help="run a fast installation and band-model smoke test",
    )
    arguments = argument_parser.parse_args()

    try:
        main(check_only=arguments.check)
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1) from error