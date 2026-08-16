#pragma once

#include <eigen3/Eigen/Dense>


// Directory containing the HAMSTER DOYxxx.nc files.
// The executable is run from the project root.
const char hamster_dir[] = "C/lib/HAMSTER";

// HAMSTER day of year to load.
const unsigned int hamster_day_of_year = 265;

// Maximum memory used to store albedo data after loading.
// The spatial sampling interval is chosen automatically from the file size.
const unsigned int hamster_max_albedo_memory_mb = 512;


/**
 * @brief In-memory HAMSTER spectral albedo data.
 *
 * The source NetCDF file is read once by read_albedo_data() and then closed.
 * The complete wavelength axis is retained, while latitude and longitude are spatially subsampled as required to keep the albedo matrix within the configured memory limit.
 *
 * The albedo matrix is arranged as:
 *
 *     row = latitude_index * n_lon + longitude_index
 *     column = wavelength_index
 */
typedef struct
{
    unsigned int n_lat;
    unsigned int n_lon;
    unsigned int n_wavelength;

    Eigen::VectorXd latitude;
    Eigen::VectorXd longitude;
    Eigen::VectorXd wavelength;

    Eigen::Matrix<
        double,
        Eigen::Dynamic,
        Eigen::Dynamic,
        Eigen::RowMajor
    > albedo;

} AlbedoData;


/**
 * @brief Load the selected HAMSTER file into memory.
 *
 * The latitude, longitude, wavelength, validity limits, and albedo data are read directly from the NetCDF file. Spatial data are automatically subsampled if required to keep memory use below hamster_max_albedo_memory_mb.
 *
 * The NetCDF file is closed before this function returns.
 *
 * @return Loaded HAMSTER albedo data.
 *
 * @throws std::invalid_argument if hamster_day_of_year is outside [1, 365].
 * @throws std::runtime_error if the NetCDF file cannot be read or does not contain the expected HAMSTER variables.
 */
AlbedoData read_albedo_data();


/**
 * @brief Get the HAMSTER albedo spectrum at the nearest stored pixel.
 *
 * No interpolation and no file access are performed. The closest stored latitude and longitude values are found independently and the corresponding spectral albedo vector is returned.
 *
 * Longitude is treated as periodic, so values around -180 and +180 degrees are handled correctly.
 *
 * @param albedo_data Loaded HAMSTER data returned by read_albedo_data().
 * @param latitude_deg Latitude in degrees, in the range [-90, 90].
 * @param longitude_deg Longitude in degrees.
 *
 * @return Spectral albedo at the nearest stored HAMSTER pixel.
 *
 * @throws std::invalid_argument if latitude_deg is outside [-90, 90].
 */
Eigen::VectorXd get_albedo_spectrum(
    const AlbedoData& albedo_data,
    double latitude_deg,
    double longitude_deg
);