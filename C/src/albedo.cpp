#include "albedo.hpp"
#include <netcdf.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>


/**
 * @brief Check the return status of a NetCDF function.
 *
 * @param status NetCDF return code.
 *
 * @return Nothing.
 *
 * @throws std::runtime_error if the NetCDF call failed.
 */
static void check_nc(int status)
{
    if (status != NC_NOERR)
    {
        throw std::runtime_error(nc_strerror(status));
    }
}


/**
 * @brief Find the closest value in a coordinate vector.
 *
 * @param values Coordinate values.
 * @param requested_value Requested coordinate.
 *
 * @return Index of the closest coordinate value.
 */
static unsigned int find_closest_index(
    const Eigen::VectorXd& values,
    double requested_value
)
{
    unsigned int closest_index = 0;
    unsigned int n_values = (unsigned int)values.size();

    double closest_error = std::fabs(values(0) - requested_value);

    for (unsigned int i = 1; i < n_values; i++)
    {
        double error = std::fabs(values(i) - requested_value);

        if (error < closest_error)
        {
            closest_error = error;
            closest_index = i;
        }
    }

    return closest_index;
}


/**
 * @brief Find the closest longitude in a periodic longitude vector.
 *
 * Longitude values separated by 360 degrees are treated as identical.
 *
 * @param longitude Stored longitude coordinates.
 * @param requested_longitude Requested longitude in degrees.
 *
 * @return Index of the closest longitude value.
 */
static unsigned int find_closest_longitude_index(
    const Eigen::VectorXd& longitude,
    double requested_longitude
)
{
    unsigned int closest_index = 0;
    unsigned int n_values = (unsigned int)longitude.size();

    double difference = std::fabs(longitude(0) - requested_longitude);
    difference = std::fmod(difference, 360.0);

    if (difference > 180.0)
    {
        difference = 360.0 - difference;
    }

    double closest_error = difference;

    for (unsigned int i = 1; i < n_values; i++)
    {
        difference = std::fabs(longitude(i) - requested_longitude);
        difference = std::fmod(difference, 360.0);

        if (difference > 180.0)
        {
            difference = 360.0 - difference;
        }

        if (difference < closest_error)
        {
            closest_error = difference;
            closest_index = i;
        }
    }

    return closest_index;
}


/**
 * @brief Load the configured HAMSTER file into a bounded in-memory structure.
 *
 * All coordinate information is read from the NetCDF file. The complete
 * wavelength axis is retained. Latitude and longitude are spatially sampled
 * using one automatically selected stride so that the stored albedo matrix
 * remains below hamster_max_albedo_memory_mb.
 *
 * @return Loaded HAMSTER albedo data.
 *
 * @throws std::invalid_argument if hamster_day_of_year is outside [1, 365].
 * @throws std::runtime_error if the file cannot be read or has an unexpected
 * structure.
 */
AlbedoData read_albedo_data()
{
    if (hamster_day_of_year < 1 || hamster_day_of_year > 365)
    {
        throw std::invalid_argument(
            "HAMSTER day of year must be between 1 and 365."
        );
    }

    // Construct the selected HAMSTER filename.
    char filename[256];

    std::snprintf(
        filename,
        sizeof(filename),
        "%s/DOY%03u.nc",
        hamster_dir,
        hamster_day_of_year
    );

    int file_id;
    check_nc(nc_open(filename, NC_NOWRITE, &file_id));

    // Find the coordinate and albedo variables.
    int latitude_var_id;
    int longitude_var_id;
    int wavelength_var_id;
    int albedo_var_id;

    check_nc(nc_inq_varid(file_id, "Latitude", &latitude_var_id));
    check_nc(nc_inq_varid(file_id, "Longitude", &longitude_var_id));
    check_nc(nc_inq_varid(file_id, "Wavelength", &wavelength_var_id));
    check_nc(nc_inq_varid(file_id, "Black_Sky_Albedo", &albedo_var_id));

    // Read the dimension IDs used by the coordinate variables.
    int latitude_dim_id;
    int longitude_dim_id;
    int wavelength_dim_id;

    check_nc(nc_inq_vardimid(file_id, latitude_var_id, &latitude_dim_id));
    check_nc(nc_inq_vardimid(file_id, longitude_var_id, &longitude_dim_id));
    check_nc(nc_inq_vardimid(file_id, wavelength_var_id, &wavelength_dim_id));

    // NetCDF requires size_t for dimension lengths.
    size_t file_n_lat_nc;
    size_t file_n_lon_nc;
    size_t file_n_wavelength_nc;

    check_nc(nc_inq_dimlen(file_id, latitude_dim_id, &file_n_lat_nc));
    check_nc(nc_inq_dimlen(file_id, longitude_dim_id, &file_n_lon_nc));
    check_nc(nc_inq_dimlen(file_id, wavelength_dim_id, &file_n_wavelength_nc));

    unsigned int file_n_lat = (unsigned int)file_n_lat_nc;
    unsigned int file_n_lon = (unsigned int)file_n_lon_nc;
    unsigned int file_n_wavelength = (unsigned int)file_n_wavelength_nc;

    if (
        file_n_lat < 1 ||
        file_n_lon < 1 ||
        file_n_wavelength < 1
    )
    {
        nc_close(file_id);
        throw std::runtime_error(
            "HAMSTER file contains an empty coordinate dimension."
        );
    }

    // Confirm that Black_Sky_Albedo is ordered as [lat, lon, wavelength].
    int n_albedo_dimensions;
    int albedo_dim_ids[3];

    check_nc(
        nc_inq_varndims(
            file_id,
            albedo_var_id,
            &n_albedo_dimensions
        )
    );

    if (n_albedo_dimensions != 3)
    {
        nc_close(file_id);
        throw std::runtime_error(
            "Black_Sky_Albedo must have three dimensions."
        );
    }

    check_nc(
        nc_inq_vardimid(
            file_id,
            albedo_var_id,
            albedo_dim_ids
        )
    );

    if (
        albedo_dim_ids[0] != latitude_dim_id ||
        albedo_dim_ids[1] != longitude_dim_id ||
        albedo_dim_ids[2] != wavelength_dim_id
    )
    {
        nc_close(file_id);
        throw std::runtime_error(
            "Black_Sky_Albedo dimensions must be [lat, lon, wavelength]."
        );
    }

    // Read the physical validity range directly from the file.
    double valid_min;
    double valid_max;

    check_nc(
        nc_get_att_double(
            file_id,
            albedo_var_id,
            "valid_min",
            &valid_min
        )
    );

    check_nc(
        nc_get_att_double(
            file_id,
            albedo_var_id,
            "valid_max",
            &valid_max
        )
    );

    // Determine the spatial stride required to keep the stored data bounded.
    // Wavelength is never subsampled.
    double full_data_bytes =
        (double)file_n_lat *
        (double)file_n_lon *
        (double)file_n_wavelength *
        (double)sizeof(double);

    double maximum_data_bytes =
        (double)hamster_max_albedo_memory_mb *
        1024.0 *
        1024.0;

    unsigned int spatial_stride = 1;

    if (full_data_bytes > maximum_data_bytes)
    {
        spatial_stride = (unsigned int)std::ceil(
            std::sqrt(full_data_bytes / maximum_data_bytes)
        );
    }

    // Number of latitude and longitude samples retained after subsampling.
    unsigned int n_lat =
        (file_n_lat - 1) / spatial_stride + 1;

    unsigned int n_lon =
        (file_n_lon - 1) / spatial_stride + 1;

    AlbedoData albedo_data;

    albedo_data.n_lat = n_lat;
    albedo_data.n_lon = n_lon;
    albedo_data.n_wavelength = file_n_wavelength;

    albedo_data.latitude.resize(n_lat);
    albedo_data.longitude.resize(n_lon);
    albedo_data.wavelength.resize(file_n_wavelength);

    albedo_data.albedo.resize(
        n_lat * n_lon,
        file_n_wavelength
    );

    // NetCDF requires size_t and ptrdiff_t for strided reads.
    // These types are confined to the NetCDF API boundary.
    size_t axis_start[1] = {0};

    size_t latitude_count[1] = {(size_t)n_lat};
    size_t longitude_count[1] = {(size_t)n_lon};

    ptrdiff_t spatial_axis_stride[1] = {
        (ptrdiff_t)spatial_stride
    };

    // Read the stored latitude and longitude coordinates using the same stride
    // that will be used for the albedo data.
    check_nc(
        nc_get_vars_double(
            file_id,
            latitude_var_id,
            axis_start,
            latitude_count,
            spatial_axis_stride,
            albedo_data.latitude.data()
        )
    );

    check_nc(
        nc_get_vars_double(
            file_id,
            longitude_var_id,
            axis_start,
            longitude_count,
            spatial_axis_stride,
            albedo_data.longitude.data()
        )
    );

    // Keep the complete wavelength axis.
    check_nc(
        nc_get_var_double(
            file_id,
            wavelength_var_id,
            albedo_data.wavelength.data()
        )
    );

    // Read the spatially subsampled albedo cube in one operation.
    // NetCDF stores [lat][lon][wavelength] with wavelength varying fastest.
    // The row-major Eigen matrix therefore receives the data directly as:
    //
    //     row = lat_index * n_lon + lon_index
    //     column = wavelength_index
    size_t albedo_start[3] = {0, 0, 0};

    size_t albedo_count[3] = {
        (size_t)n_lat,
        (size_t)n_lon,
        (size_t)file_n_wavelength
    };

    ptrdiff_t albedo_stride[3] = {
        (ptrdiff_t)spatial_stride,
        (ptrdiff_t)spatial_stride,
        1
    };

    check_nc(
        nc_get_vars_double(
            file_id,
            albedo_var_id,
            albedo_start,
            albedo_count,
            albedo_stride,
            albedo_data.albedo.data()
        )
    );

    check_nc(nc_close(file_id));

    // Replace invalid HAMSTER values with NaN once during loading.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    unsigned int n_values = (unsigned int)albedo_data.albedo.size();

    for (unsigned int i = 0; i < n_values; i++)
    {
        if (
            albedo_data.albedo.data()[i] < valid_min ||
            albedo_data.albedo.data()[i] > valid_max
        )
        {
            albedo_data.albedo.data()[i] = nan;
        }
    }

    std::cout
        << "Loaded HAMSTER DOY "
        << hamster_day_of_year
        << ": "
        << albedo_data.n_lat
        << " latitude pixels x "
        << albedo_data.n_lon
        << " longitude pixels x "
        << albedo_data.n_wavelength
        << " wavelengths."
        << std::endl;

    return albedo_data;
}


/**
 * @brief Get the spectral albedo at the nearest stored latitude/longitude.
 *
 * This function performs no interpolation and no NetCDF file access.
 *
 * @param albedo_data Loaded HAMSTER data.
 * @param latitude_deg Latitude in degrees, in the range [-90, 90].
 * @param longitude_deg Longitude in degrees.
 *
 * @return Spectral albedo at the nearest stored geographic pixel.
 *
 * @throws std::invalid_argument if latitude_deg is outside [-90, 90].
 */
Eigen::VectorXd get_albedo_spectrum(
    const AlbedoData& albedo_data,
    double latitude_deg,
    double longitude_deg
)
{
    if (latitude_deg < -90.0 || latitude_deg > 90.0)
    {
        throw std::invalid_argument(
            "Latitude must be between -90 and 90 degrees."
        );
    }

    unsigned int latitude_index = find_closest_index(
        albedo_data.latitude,
        latitude_deg
    );

    unsigned int longitude_index = find_closest_longitude_index(
        albedo_data.longitude,
        longitude_deg
    );

    unsigned int row =
        latitude_index * albedo_data.n_lon +
        longitude_index;

    return albedo_data.albedo.row(row).transpose();
}