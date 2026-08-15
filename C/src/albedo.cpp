#include "albedo.hpp"
#include <netcdf.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>


/**
 * @brief Check the return value from a NetCDF function.
 *
 * @param status NetCDF status code.
 */
void check_nc(int status)
{
    if (status != NC_NOERR)
    {
        throw std::runtime_error(nc_strerror(status));
    }
}


/**
 * @brief Read a one-dimensional coordinate variable from the NetCDF file.
 *
 * The size of the coordinate axis is read directly from the file.
 *
 * @param file_id Open NetCDF file ID.
 * @param variable_name Name of the coordinate variable.
 * @param values Eigen vector that receives the coordinate values.
 *
 * @return Number of values read.
 */
int read_axis(
    int file_id,
    const char* variable_name,
    Eigen::VectorXf& values
)
{
    int variable_id;
    int dimension_id;

    // NetCDF uses size_t for dimension lengths.
    size_t n_values;

    check_nc(nc_inq_varid(file_id, variable_name, &variable_id));
    check_nc(nc_inq_vardimid(file_id, variable_id, &dimension_id));
    check_nc(nc_inq_dimlen(file_id, dimension_id, &n_values));

    values.resize((int)n_values);

    check_nc(
        nc_get_var_float(
            file_id,
            variable_id,
            values.data()
        )
    );

    return (int)n_values;
}


/**
 * @brief Find the first coordinate index surrounding a requested value.
 *
 * Works for both increasing and decreasing coordinate arrays.
 * The returned index i satisfies that the requested value lies between
 * values(i) and values(i + 1).
 *
 * @param values Coordinate array.
 * @param value Requested coordinate.
 *
 * @return Lower array index of the surrounding pair.
 */
int find_interval(
    const Eigen::VectorXf& values,
    float value
)
{
    int low = 0;
    int high = values.size() - 1;
    bool increasing = values(1) > values(0);

    while (high - low > 1)
    {
        int middle = (low + high) / 2;

        if (
            (increasing && values(middle) <= value) ||
            (!increasing && values(middle) >= value)
        )
        {
            low = middle;
        }
        else
        {
            high = middle;
        }
    }

    return low;
}


/**
 * @brief Read one complete albedo spectrum from one HAMSTER pixel.
 *
 * Only one latitude/longitude pixel is read from disk. Invalid HAMSTER
 * values are replaced with NaN using valid_min and valid_max from the file.
 *
 * @param albedo_data Loaded HAMSTER file information.
 * @param latitude_index Latitude pixel index.
 * @param longitude_index Longitude pixel index.
 *
 * @return Albedo spectrum at the selected pixel.
 */
Eigen::VectorXf read_pixel_spectrum(
    const AlbedoData& albedo_data,
    int latitude_index,
    int longitude_index
)
{
    Eigen::VectorXf spectrum(albedo_data.n_wavelength);

    // nc_get_vara_float requires size_t arrays for the start and count.
    size_t start[3] = {
        (size_t)latitude_index,
        (size_t)longitude_index,
        0
    };

    size_t count[3] = {
        1,
        1,
        (size_t)albedo_data.n_wavelength
    };

    check_nc(
        nc_get_vara_float(
            albedo_data.file_id,
            albedo_data.albedo_var_id,
            start,
            count,
            spectrum.data()
        )
    );

    for (int i = 0; i < albedo_data.n_wavelength; i++)
    {
        if (
            spectrum(i) < albedo_data.valid_min ||
            spectrum(i) > albedo_data.valid_max
        )
        {
            spectrum(i) = std::numeric_limits<float>::quiet_NaN();
        }
    }

    return spectrum;
}


/**
 * @brief Open the selected HAMSTER file and read its metadata.
 *
 * The latitude, longitude and wavelength axes are read directly from the
 * NetCDF file. The full albedo dataset remains on disk.
 *
 * @return AlbedoData containing the open NetCDF file and coordinate axes.
 */
AlbedoData read_albedo_data()
{
    AlbedoData albedo_data;

    char filename[256];

    std::snprintf(
        filename,
        sizeof(filename),
        "%s/DOY%03d.nc",
        hamster_dir,
        hamster_day_of_year
    );

    check_nc(
        nc_open(
            filename,
            NC_NOWRITE,
            &albedo_data.file_id
        )
    );

    albedo_data.n_lat = read_axis(
        albedo_data.file_id,
        "Latitude",
        albedo_data.latitude
    );

    albedo_data.n_lon = read_axis(
        albedo_data.file_id,
        "Longitude",
        albedo_data.longitude
    );

    albedo_data.n_wavelength = read_axis(
        albedo_data.file_id,
        "Wavelength",
        albedo_data.wavelength
    );

    check_nc(
        nc_inq_varid(
            albedo_data.file_id,
            "Black_Sky_Albedo",
            &albedo_data.albedo_var_id
        )
    );

    check_nc(
        nc_get_att_float(
            albedo_data.file_id,
            albedo_data.albedo_var_id,
            "valid_min",
            &albedo_data.valid_min
        )
    );

    check_nc(
        nc_get_att_float(
            albedo_data.file_id,
            albedo_data.albedo_var_id,
            "valid_max",
            &albedo_data.valid_max
        )
    );

    std::cout
        << "Loaded HAMSTER DOY "
        << hamster_day_of_year
        << ": "
        << albedo_data.n_lat
        << " latitude pixels x "
        << albedo_data.n_lon
        << " longitude pixels."
        << std::endl;

    return albedo_data;
}


/**
 * @brief Get the spectral albedo at a requested latitude and longitude.
 *
 * The four surrounding HAMSTER pixels are read from disk and bilinearly
 * interpolated. Longitude wraps across the -180/180 degree boundary.
 * Latitude is clamped to the first or last HAMSTER pixel at the poles.
 *
 * @param albedo_data Loaded HAMSTER data returned by read_albedo_data().
 * @param latitude_deg Latitude in degrees, in the range [-90, 90].
 * @param longitude_deg Longitude in degrees. Values outside [-180, 180]
 * are wrapped automatically.
 *
 * @return Interpolated albedo spectrum.
 */
Eigen::VectorXf get_albedo_spectrum(
    const AlbedoData& albedo_data,
    float latitude_deg,
    float longitude_deg
)
{
    if (latitude_deg < -90.0f || latitude_deg > 90.0f)
    {
        throw std::invalid_argument(
            "Latitude must be between -90 and 90 degrees."
        );
    }

    // Wrap longitude to [-180, 180).
    longitude_deg = std::fmod(longitude_deg + 180.0f, 360.0f);

    if (longitude_deg < 0.0f)
    {
        longitude_deg += 360.0f;
    }

    longitude_deg -= 180.0f;


    // Find the two latitude pixels surrounding the requested latitude.
    int lat_0;
    int lat_1;
    float lat_fraction;

    if (latitude_deg >= albedo_data.latitude(0))
    {
        lat_0 = 0;
        lat_1 = 0;
        lat_fraction = 0.0f;
    }
    else if (latitude_deg <= albedo_data.latitude(albedo_data.n_lat - 1))
    {
        lat_0 = albedo_data.n_lat - 1;
        lat_1 = lat_0;
        lat_fraction = 0.0f;
    }
    else
    {
        lat_0 = find_interval(albedo_data.latitude, latitude_deg);
        lat_1 = lat_0 + 1;

        lat_fraction =
            (latitude_deg - albedo_data.latitude(lat_0)) /
            (albedo_data.latitude(lat_1) - albedo_data.latitude(lat_0));
    }


    // Find the two longitude pixels surrounding the requested longitude.
    int lon_0;
    int lon_1;
    float lon_fraction;

    if (
        longitude_deg < albedo_data.longitude(0) ||
        longitude_deg > albedo_data.longitude(albedo_data.n_lon - 1)
    )
    {
        // The requested position lies between the last and first longitude
        // pixels across the -180/180 degree boundary.
        lon_0 = albedo_data.n_lon - 1;
        lon_1 = 0;

        float lon = longitude_deg;

        if (lon < albedo_data.longitude(0))
        {
            lon += 360.0f;
        }

        lon_fraction =
            (lon - albedo_data.longitude(lon_0)) /
            (albedo_data.longitude(lon_1) + 360.0f -
             albedo_data.longitude(lon_0));
    }
    else
    {
        lon_0 = find_interval(albedo_data.longitude, longitude_deg);
        lon_1 = lon_0 + 1;

        lon_fraction =
            (longitude_deg - albedo_data.longitude(lon_0)) /
            (albedo_data.longitude(lon_1) - albedo_data.longitude(lon_0));
    }


    // Read only the four surrounding spectra.
    Eigen::VectorXf albedo_00 =
        read_pixel_spectrum(albedo_data, lat_0, lon_0);

    Eigen::VectorXf albedo_01 =
        read_pixel_spectrum(albedo_data, lat_0, lon_1);

    Eigen::VectorXf albedo_10 =
        read_pixel_spectrum(albedo_data, lat_1, lon_0);

    Eigen::VectorXf albedo_11 =
        read_pixel_spectrum(albedo_data, lat_1, lon_1);


    // Bilinear interpolation. NaN values naturally propagate through
    // the Eigen arithmetic if any required HAMSTER pixel is invalid.
    Eigen::VectorXf albedo_lat_0 =
        albedo_00 + lon_fraction * (albedo_01 - albedo_00);

    Eigen::VectorXf albedo_lat_1 =
        albedo_10 + lon_fraction * (albedo_11 - albedo_10);

    return
        albedo_lat_0 +
        lat_fraction * (albedo_lat_1 - albedo_lat_0);
}


/**
 * @brief Close the HAMSTER NetCDF file.
 *
 * @param albedo_data HAMSTER data structure containing the open file ID.
 */
void close_albedo_data(AlbedoData& albedo_data)
{
    check_nc(nc_close(albedo_data.file_id));
    albedo_data.file_id = -1;
}