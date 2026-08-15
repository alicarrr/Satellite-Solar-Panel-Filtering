import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# latex formatting for plots
plt.rcParams.update({"text.usetex": True})
plt.rcParams.update({"font.size": 14, "figure.figsize": (8, 6)})


def get_albedo_spectrum_at_lat_and_lon(latitude_deg: float, longitude_deg: float, day_of_year: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Get the HAMSTER albedo spectrum at a given latitude and longitude for a specific day of the year. This assumes that the HAMSTER data file for that day of the year is available in the specified directory.

    Args:
        latitude_deg (float): latitude in degrees, must be between -90 and 90.
        longitude_deg (float): longitude in degrees, should be between -180 and 180.
        day_of_year (int, optional): day of the year (1-365). Defaults to 1.

    Raises:
        ValueError: If the latitude is not within the valid range.
        ValueError: If the longitude is not within the valid range.
        ValueError: If the day of the year is not within the valid range.
        ValueError: If the HAMSTER data file for the specified day of the year is not found.

    Returns:
        tuple[np.ndarray, np.ndarray]: The wavelength (in 10nm steps) and albedo arrays at the specified latitude and longitude for the given day of the year.
    """
    # Check if the latitude is within valid range
    if latitude_deg < -90 or latitude_deg > 90:
        raise ValueError(f"Latitude must be between -90 and 90 degrees. Given: {latitude_deg}")

    # Check if valid day of year is given (must be between 1 and 365)
    if day_of_year < 1 or day_of_year > 365:
        raise ValueError(f"Day of year must be between 1 and 365. Given: {day_of_year}")

    # Wrap longitude to the range [-180, 180]
    longitude_deg = ((longitude_deg + 180) % 360) - 180

    # Find data for the given day of year
    hamster_dir = "../../C/lib/HAMSTER"
    filename = f"{hamster_dir}/DOY{day_of_year:03d}.nc"

    # Check if the file exists
    try:
        with open(filename, "r"):
            pass
    except FileNotFoundError:
        raise ValueError(f"HAMSTER data file for day {day_of_year} not found at {filename}. Check the file is downloaded")

    # Open the netCDF file and read the data
    with nc.Dataset(filename, "r") as data:
        # lat is from -90 to 90
        latitude = np.asarray(data.variables["Latitude"][:])

        # lon is from -180 to 180
        longitude = np.asarray(data.variables["Longitude"][:])

        # wavelength is from 400 to 2500 nm, in 10 nm increments
        wavelength = np.asarray(data.variables["Wavelength"][:])

        # Find nearest HAMSTER grid point
        latitude_index = int(
            np.argmin(np.abs(latitude - latitude_deg))
        )
        longitude_index = int(
            np.argmin(np.abs(longitude - longitude_deg))
        )

        # Invalid albedo values are automatically masked by netCDF4 using the variable's valid range [0.0, 1.0].
        albedo = data.variables["Black_Sky_Albedo"]
        albedo_spectrum = np.ma.filled(albedo[latitude_index, longitude_index, :], fill_value=np.nan)

        print(
            f"Requested position: \n {latitude_deg:.3f} deg, {longitude_deg:.3f} deg"
        )
        print(
            f"HAMSTER position: \n{latitude[latitude_index]:.3f} deg, "
            f"{longitude[longitude_index]:.3f} deg"
        )

    return wavelength, albedo_spectrum


if __name__ == "__main__":

    # Test location
    latitude_deg = 23
    longitude_deg = 13

    day_of_year = 1

    try:
        wavelength, albedo = get_albedo_spectrum_at_lat_and_lon(latitude_deg, longitude_deg, day_of_year)
    except ValueError as e:
        print(f"Error: {e}")
        exit(1)

    plt.plot(wavelength, albedo)

    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Black-sky albedo")
    plt.title(
        f"HAMSTER spectral albedo\n"
        f"Lat: {latitude_deg:.3f}$^\\circ$, "
        f"Lon: {longitude_deg:.3f}$^\\circ$"
    )

    plt.grid()
    plt.tight_layout()
    plt.show()
    plt.savefig("albedo_spectrum.png", dpi=300)