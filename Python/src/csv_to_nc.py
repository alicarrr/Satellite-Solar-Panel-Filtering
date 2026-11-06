import argparse
import pandas as pd
import xarray as xr
from pathlib import Path

SCRIPT_DIRECTORY = Path(__file__).resolve().parent


def csv_to_netcdf(csv_file: str, nc_file: str) -> None:
    """
    Convert a CSV file to a NetCDF file.

    Each CSV column becomes a NetCDF variable, and all variables share
    a common 'row' dimension.

    Parameters
    ----------
    csv_file : str
        Path to the input CSV file.
    nc_file : str
        Path to the output NetCDF file.
    """

    # Read CSV
    dataframe = pd.read_csv(csv_file)

    # Convert columns into an xarray Dataset
    dataset = xr.Dataset(
        {
            column: ("row", dataframe[column].to_numpy())
            for column in dataframe.columns
        }
    )

    # Save as NetCDF
    dataset.to_netcdf(nc_file)

    print(f"Converted '{csv_file}' -> '{nc_file}'")


if __name__ == "__main__":
    # parser = argparse.ArgumentParser(
    #     description="Convert a CSV file to NetCDF."
    # )

    # parser.add_argument(
    #     "csv_file",
    #     help="Input CSV file"
    # )

    # parser.add_argument(
    #     "nc_file",
    #     help="Output NetCDF file"
    # )

    # args = parser.parse_args()

    OUTPUT_DIRECTORY = (SCRIPT_DIRECTORY / "../../dat files").resolve()
    csv_file = OUTPUT_DIRECTORY / "./earth_atmosphere_transmissivity.csv"
    nc_file = OUTPUT_DIRECTORY / "./earth_atmosphere_transmissivity.nc"

    csv_to_netcdf(csv_file, nc_file)