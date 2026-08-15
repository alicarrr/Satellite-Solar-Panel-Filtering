import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


def get_albedo_spectrum(
    latitude_deg: float,
    longitude_deg: float,
    day_of_year: int = 1,
) -> tuple[np.ndarray, np.ndarray]:

    hamster_dir = "../../C/lib/HAMSTER"

    filename = f"{hamster_dir}/DOY{day_of_year:03d}.nc"

    with nc.Dataset(filename, "r") as data:

        latitude = np.asarray(data.variables["Latitude"][:])
        longitude = np.asarray(data.variables["Longitude"][:])
        wavelength = np.asarray(data.variables["Wavelength"][:])

        albedo = data.variables["Black_Sky_Albedo"]

        # Find nearest HAMSTER grid point
        latitude_index = int(
            np.argmin(np.abs(latitude - latitude_deg))
        )

        longitude_index = int(
            np.argmin(np.abs(longitude - longitude_deg))
        )

        # Black_Sky_Albedo has dimensions:
        # [latitude][longitude][wavelength]
        albedo_spectrum = np.asarray(
            albedo[
                latitude_index,
                longitude_index,
                :
            ]
        )

        print(
            f"Requested position: "
            f"{latitude_deg:.3f}°, {longitude_deg:.3f}°"
        )

        print(
            f"HAMSTER position:   "
            f"{latitude[latitude_index]:.3f}°, "
            f"{longitude[longitude_index]:.3f}°"
        )

    return wavelength, albedo_spectrum


if __name__ == "__main__":

    # Test location
    latitude_deg = -33.8688
    longitude_deg = 151.2093

    day_of_year = 1

    wavelength, albedo = get_albedo_spectrum(
        latitude_deg,
        longitude_deg,
        day_of_year,
    )

    plt.plot(wavelength, albedo)

    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Black-sky albedo")
    plt.title(
        f"HAMSTER spectral albedo\n"
        f"Lat: {latitude_deg:.3f}°, "
        f"Lon: {longitude_deg:.3f}°"
    )

    plt.grid()
    plt.tight_layout()
    plt.show()
    plt.savefig("albedo_spectrum.png", dpi=300)