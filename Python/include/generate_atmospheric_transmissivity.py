#!/usr/bin/env python3
"""Generate a clear-sky atmospheric transmissivity spectrum with libRadtran.

The script calls libRadtran's ``uvspec`` executable once or twice, then writes a
CSV with the wavelength and three dimensionless quantities:

  * direct_beam_transmissivity
  * diffuse_downward_transmissivity
  * total_downward_transmissivity (direct + diffuse)

All three are normalised by the extraterrestrial irradiance projected onto a
horizontal plane, E0*cos(SZA).  Consequently, the direct column is the
atmospheric transmission along the solar ray without an extra cosine factor.

libRadtran has separate REPTRAN band tables for the solar and thermal spectral
ranges.  This script uses the solar table through 5000 nm and the thermal table
above 5000 nm, joining both results into one ordered CSV.  Both calculations use
``source solar`` as a unit external collimated beam.  This is deliberate: it
calculates atmospheric transmission without adding thermal emission by the
atmosphere or surface.

AIR_MASS uses the simple plane-parallel relation AM = sec(SZA).  AIR_MASS = 0
is handled specially as no atmosphere (AM0), for which direct = total = 1 and
diffuse = 0.  Otherwise AIR_MASS must be at least 1.

Only Python's standard library is required.  Edit the constants immediately
below to select the libRadtran path, wavelength grid, air mass, and atmosphere.
Run with ``--check`` first for a fast end-to-end smoke test of every required
REPTRAN table.
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
from decimal import Decimal
from pathlib import Path


# ---------------------------------------------------------------------------
# User-editable constants
# ---------------------------------------------------------------------------

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

# The combined solar and thermal REPTRAN tables support 120 nm to just below
# 100000 nm (100 micrometres). The default below produces 280 nm to 50 um.
MIN_WAVELENGTH_NM = 280.0
MAX_WAVELENGTH_NM = 50000.0
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

# "coarse" is included in the normal libRadtran download.  "medium" and
# "fine" require the additional REPTRAN data archive from libradtran.org.
REPTRAN_RESOLUTION = "coarse"  # "coarse", "medium", or "fine"

RTE_SOLVER = "disort"
NUMBER_OF_STREAMS = 16

# Internal table-selection constants. The tables overlap from 2500-5000 nm;
# the solar parameterisation is used through 5000 nm and thermal above it.
REPTRAN_MIN_WAVELENGTH_NM = 120.0
REPTRAN_SOLAR_THERMAL_SPLIT_NM = 5000.0
REPTRAN_MAX_WAVELENGTH_NM_EXCLUSIVE = 100000.0


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
    minimum_supported = Decimal(str(REPTRAN_MIN_WAVELENGTH_NM))
    maximum_exclusive = Decimal(
        str(REPTRAN_MAX_WAVELENGTH_NM_EXCLUSIVE)
    )
    if minimum < minimum_supported or maximum >= maximum_exclusive:
        raise ValueError(
            "The combined REPTRAN configuration requires wavelengths from "
            f"{REPTRAN_MIN_WAVELENGTH_NM:g} nm up to, but not including, "
            f"{REPTRAN_MAX_WAVELENGTH_NM_EXCLUSIVE:g} nm."
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


def split_wavelengths_by_reptran_table(
    wavelengths: list[Decimal],
) -> list[tuple[str, list[Decimal]]]:
    """Split an ordered grid into solar- and thermal-REPTRAN segments."""

    split = Decimal(str(REPTRAN_SOLAR_THERMAL_SPLIT_NM))
    solar_wavelengths = [
        wavelength for wavelength in wavelengths if wavelength <= split
    ]
    thermal_wavelengths = [
        wavelength for wavelength in wavelengths if wavelength > split
    ]

    segments: list[tuple[str, list[Decimal]]] = []
    if solar_wavelengths:
        segments.append(("solar", solar_wavelengths))
    if thermal_wavelengths:
        segments.append(("thermal", thermal_wavelengths))
    return segments


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

    resolution = REPTRAN_RESOLUTION.lower()
    if resolution not in {"coarse", "medium", "fine"}:
        raise ValueError(
            "REPTRAN_RESOLUTION must be 'coarse', 'medium', or 'fine'"
        )

    reptran_directory = (
        LIBRADTRAN_DIRECTORY / "data/correlated_k/reptran"
    )
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
    required_reptran_files: list[Path] = []
    for table_name, _ in split_wavelengths_by_reptran_table(wavelengths):
        table_base = f"reptran_{table_name}_{resolution}"
        required_reptran_files.append(
            reptran_directory / f"{table_base}.cdf"
        )
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

    return uvspec_executable, atmosphere_path


def make_uvspec_input(
    grid_file: Path,
    solar_zenith_deg: float,
    reptran_table: str,
) -> str:
    """Create a libRadtran uvspec input document."""

    resolution = REPTRAN_RESOLUTION.lower()
    if reptran_table not in {"solar", "thermal"}:
        raise ValueError("reptran_table must be 'solar' or 'thermal'")
    if AEROSOL_VISIBILITY_KM <= 0:
        raise ValueError("AEROSOL_VISIBILITY_KM must be greater than zero")
    if NUMBER_OF_STREAMS <= 0 or NUMBER_OF_STREAMS % 2 != 0:
        raise ValueError("NUMBER_OF_STREAMS must be a positive even integer")

    # output_quantity reflectivity is intentionally used here. In libRadtran,
    # that option normalises irradiances by E0*cos(SZA), whereas the option named
    # transmittance normalises by E0 and retains the geometric cosine factor.
    #
    # `source solar` is also intentional for the longwave segment. It supplies
    # a unit collimated external beam without atmospheric/surface emission.
    # `reptran_file` independently selects the appropriate absorption table.
    reptran_file_base = (
        "data/correlated_k/reptran/"
        f"reptran_{reptran_table}_{resolution}"
    )
    lines = [
        "data_files_path data/",
        f"atmosphere_file {ATMOSPHERE_FILE.as_posix()}",
        "source solar",
        # With REPTRAN, libRadtran 2.0.6 rejects a simultaneous `wavelength`
        # directive and `wavelength_grid_file`.  The grid file alone defines
        # both the requested range and sampling.
        f"wavelength_grid_file {grid_file.as_posix()}",
        f"mol_abs_param reptran {resolution}",
        f"reptran_file {reptran_file_base}",
        f"sza {solar_zenith_deg:.12g}",
        "albedo 0.0",
        f"rte_solver {RTE_SOLVER}",
        f"number_of_streams {NUMBER_OF_STREAMS}",
        "output_quantity reflectivity",
        "output_user lambda edir edn eglo",
    ]

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
    reptran_table: str,
) -> list[tuple[float, float, float, float]]:
    """Run uvspec for a supplied output grid and validate the result."""

    with tempfile.TemporaryDirectory(
        prefix=f"libradtran_transmissivity_{reptran_table}_"
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
            reptran_table,
        )
        output_text = run_uvspec(uvspec_executable, uvspec_input)

    return parse_uvspec_output(output_text, wavelengths)


def run_all_reptran_segments(
    uvspec_executable: Path,
    wavelengths: list[Decimal],
    solar_zenith_deg: float,
) -> list[tuple[float, float, float, float]]:
    """Calculate and join all required REPTRAN spectral segments."""

    rows: list[tuple[float, float, float, float]] = []
    for reptran_table, segment_wavelengths in (
        split_wavelengths_by_reptran_table(wavelengths)
    ):
        rows.extend(
            run_calculation(
                uvspec_executable,
                segment_wavelengths,
                solar_zenith_deg,
                reptran_table,
            )
        )
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
            "Check the requested range and installed REPTRAN data."
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


def main(check_only: bool = False) -> None:
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
        # Exercise every required REPTRAN table as well as the executable,
        # atmosphere, solver, and parser before running the full spectrum.
        tested_wavelength_count = 0
        tested_tables: list[str] = []
        for reptran_table, segment_wavelengths in (
            split_wavelengths_by_reptran_table(wavelengths)
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
                reptran_table,
            )
            tested_wavelength_count += len(smoke_test_wavelengths)
            tested_tables.append(reptran_table)

        print(
            "Configuration check passed: uvspec successfully calculated "
            f"{tested_wavelength_count} sample wavelengths using the "
            f"{' and '.join(tested_tables)} REPTRAN table(s)."
        )
        return

    rows = run_all_reptran_segments(
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
    table_names = [
        table_name
        for table_name, _ in split_wavelengths_by_reptran_table(wavelengths)
    ]
    print(
        f"Wrote {len(rows)} rows to {OUTPUT_CSV}\n"
        f"AM={AIR_MASS:g}, SZA={solar_zenith_deg:.6f} deg, "
        f"REPTRAN={REPTRAN_RESOLUTION} "
        f"({' + '.join(table_names)} tables), {aerosol_description}"
    )


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser(
        description=(
            "Generate a clear-sky atmospheric transmissivity spectrum with "
            "libRadtran."
        )
    )
    argument_parser.add_argument(
        "--check",
        action="store_true",
        help="run a fast installation and REPTRAN-table smoke test",
    )
    arguments = argument_parser.parse_args()

    try:
        main(check_only=arguments.check)
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1) from error
