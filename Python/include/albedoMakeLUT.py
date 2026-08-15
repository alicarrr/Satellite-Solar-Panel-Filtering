from pathlib import Path

import netCDF4 as nc
import numpy as np


HAMSTER_DIR = (
    Path(__file__).resolve().parents[2]
    / "C"
    / "lib"
    / "HAMSTER"
)


def get_albedo_spectrum(
    latitude_deg: float,
    longitude_deg: float,
    day_of_year: int,
) -> tuple[np.ndarray, np.ndarray]:

    filename = HAMSTER_DIR / f"DOY{day_of_year:03d}.nc"

    with nc.Dataset(filename, "r") as data:

        latitude = np.asarray(data.variables["lat"][:])
        longitude = np.asarray(data.variables["lon"][:])
        wavelength_nm = np.asarray(data.variables["wvl"][:])

        latitude_index = int(
            np.argmin(np.abs(latitude - latitude_deg))
        )

        longitude_index = int(
            np.argmin(np.abs(longitude - longitude_deg))
        )

        albedo = np.asarray(
            data.variables["Black_Sky_Albedo"][
                latitude_index,
                longitude_index,
                :
            ]
        )

    return wavelength_nm, albedo