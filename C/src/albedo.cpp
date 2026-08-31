#include "albedo.hpp"
#include <netcdf.h>
#include <netcdf>

#include "transforms.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <fstream>
#include <sstream>

double d_lat_file;
double d_lon_file;

double lat_min_file;
double lat_max_file;
double lon_min_file;
double lon_max_file;

int albedo_spectra_file_len;

Eigen::VectorXd latitudes_file;
Eigen::VectorXd longitudes_file;
std::vector<double> wavelengths_file_nm;

netCDF::NcFile albedo_file;

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



void read_albedo_data(unsigned int doy)
{
    int width = 3;

    std::ostringstream stream;
    stream << hamster_dir << "/" << "DOY" << std::setfill('0') << std::setw(width) << doy << ".nc";
    std::string file_name = stream.str();

    std::cout << file_name << std::endl;

    albedo_file.open(file_name, netCDF::NcFile::read);
}



void printVariableInformation(const std::string& variable_name, const netCDF::NcFile& file)
{
    netCDF::NcVar variable = file.getVar(variable_name);

    if (variable.isNull())
    {
        throw std::runtime_error(
            "Variable not found: " + variable_name);
    }

    std::cout << "Variable: " << variable.getName() << '\n';
    std::cout << "Stored type: "
              << variable.getType().getName() << '\n';

    std::vector<netCDF::NcDim> dimensions = variable.getDims();

    std::cout << "Number of dimensions: "
              << dimensions.size() << '\n';

    for (int i = 0;
         i < static_cast<int>(dimensions.size());
         ++i)
    {
        std::cout
            << "  Dimension " << i
            << ": " << dimensions[i].getName()
            << ", size = " << dimensions[i].getSize()
            << '\n';
    }
}


double getVisibilityConeAngle(double distance, double r)
{
    double theta = acos(r/distance);

    return theta;
}

void openAlbedoFile()
{
    std::cout << "loading in Albedo data" << std::endl;
    read_albedo_data(hamster_day_of_year);
}

void getAlbedoFileParams()
{
    if (albedo_file.isNull())
    {
        openAlbedoFile();
    }


    #if useLargeHamster == true
        // Use for larger HAMSTER files
        netCDF::NcVar latitudeVar = albedo_file.getVar("Latitude");
        netCDF::NcVar longitudeVar = albedo_file.getVar("Longitude");
        netCDF::NcVar spectrumVar = albedo_file.getVar("Wavelength");
    #else
        // Use fpr smaller HAMSTER files
        netCDF::NcVar latitudeVar = albedo_file.getVar("lat");
        netCDF::NcVar longitudeVar = albedo_file.getVar("lon");
        netCDF::NcVar spectrumVar = albedo_file.getVar("spectral");
    #endif

    

    

    const int dimensionLat = latitudeVar.getDim(0).getSize();
    const int dimensionLon = longitudeVar.getDim(0).getSize();
    albedo_spectra_file_len = spectrumVar.getDim(0).getSize();


    std::vector<double> lats(2);
    std::vector<double> lons(2);

    std::vector<double> latitudes(dimensionLat);
    std::vector<double> longitudes(dimensionLon);


    std::vector<std::size_t> start(1);
    start[0] = 0;
    
    std::vector<std::size_t> count(1);
    count[0] = 2;

    std::vector<std::size_t> countLat(1);
    countLat[0] = dimensionLat;

    std::vector<std::size_t> countLon(1);
    countLon[0] = dimensionLon;

    std::vector<std::size_t> countWav(1);
    countWav[0] = albedo_spectra_file_len;

    latitudeVar.getVar(start, count, lats.data());
    longitudeVar.getVar(start, count, lons.data());

    latitudeVar.getVar(start, countLat, latitudes.data());
    longitudeVar.getVar(start, countLon, longitudes.data());

    wavelengths_file_nm = std::vector<double>(albedo_spectra_file_len);
    spectrumVar.getVar(start, countWav, wavelengths_file_nm.data());

    latitudes_file = Eigen::Map<Eigen::VectorXd>(latitudes.data(), latitudes.size());
    longitudes_file = Eigen::Map<Eigen::VectorXd>(longitudes.data(), longitudes.size());
    

    d_lat_file = std::abs(std::round((lats[1] - lats[0]) * 1000) / 1000);
    d_lon_file = std::abs(std::round((lons[1] - lons[0]) * 1000) / 1000);

    lat_min_file = latitudes_file.minCoeff();
    lat_max_file = latitudes_file.maxCoeff();
    lon_min_file = longitudes_file.minCoeff();
    lon_max_file = longitudes_file.maxCoeff();

}


int getIndexOfEntry(const Eigen::VectorXd& v, double val)
{
    Eigen::Index index;

    (v.array() - val).abs().minCoeff(&index);

    return (int)index;
}


double wrapValue(double val, double min_val, double max_val) {
    double range = max_val - min_val;
    if (range == 0.0) return min_val;
    
    // std::fmod handles remainder, but can return negative numbers if 'val - min_val' is negative
    double result = min_val + std::fmod(val - min_val, range);
    
    // If result is negative, shift it back into the positive range [min_val, max_val)
    if (result < min_val) {
        result += range;
    }
    
    return result;
}


std::vector<int> getLatLonIndices(Eigen::Vector3d LLH)
{
    double lat = LLH(0);
    double lon = LLH(1);

    lat = wrapValue(lat, lat_min_file, lat_max_file);
    lon = wrapValue(lon, lon_min_file, lon_max_file);

    std::vector<int> indices(2);

    // LARGE FILES: The indices in the file are ordered lat then lon
    // SMALL FILES: The indices in the file are ordered lon then lat

    #if useLargeHamster == true
        indices[0] = getIndexOfEntry(latitudes_file, lat);
        indices[1] = getIndexOfEntry(longitudes_file, lon);
    #else
        indices[1] = getIndexOfEntry(latitudes_file, lat);
        indices[0] = getIndexOfEntry(longitudes_file, lon);
    #endif

    return indices;
}


double fillAlbedoData(double lambda_nm)
{
    return 0;
}


SpectraVector_t getAlbedoSpectra(Eigen::Vector3d LLH)
{
    std::vector<int> indices = getLatLonIndices(LLH);

    std::vector<double> spectra(albedo_spectra_file_len);
    

    #if useLargeHamster == true
        netCDF::NcVar albedoData = albedo_file.getVar("Black_Sky_Albedo");
    #else
        netCDF::NcVar albedoData = albedo_file.getVar("albedo");
    #endif

    std::vector<double> albedoSpectra(albedo_spectra_file_len);

    std::vector<std::size_t> start(3);
    start[0] = indices[0];
    start[1] = indices[1];
    start[2] = 0;
    
    std::vector<std::size_t> count(3);
    count[0] = 1;
    count[1] = 1;
    count[2] = albedo_spectra_file_len;

    albedoData.getVar(start, count, albedoSpectra.data());

    // Need to rebin the data from the File
    std::vector<double> albedoSpectraExtended = resampleAndExtendSpectrum(wavelengths_file_nm, albedoSpectra, min_wavelength_nm, max_wavelength_nm, d_wavelength_nm, fillAlbedoData);

    SpectraVector_t albedoSpectraVec = Eigen::Map<SpectraVector_t>(albedoSpectraExtended.data(), albedoSpectraExtended.size());
    

    return albedoSpectraVec;
}


Eigen::Vector3d getNormalToEarthEcef(double lat, double lon)
{
    double e_earth = sqrt(1 - pow((b_earth/a_earth), 2));
    double b_over_a_earth = b_earth/a_earth;
    Eigen::Vector3d n_hat_ecef = Eigen::Vector3d({b_over_a_earth * cos(lat)*cos(lon), b_over_a_earth * cos(lat)*sin(lon), sin(lat)})/sqrt(1 - pow(e_earth * cos(lat), 2));

    return n_hat_ecef;
}

SpectraVector_t getPowerOnAlbedoSurface(Eigen::Vector3d LLH, Eigen::Vector3d r_sol_ECI, double t, const SimParams_t& simParams)
{
    

    double lat = LLH[0];
    double lon = LLH[1];

    Eigen::Vector3d r_patch_ecef = llh_to_ecef(LLH);

    /** Normal vector for oblate sphereoid */
    Eigen::Vector3d n_hat_ecef = getNormalToEarthEcef(lat, lon);

    Eigen::Vector3d r_sol_ecef = simParams.ephemeris->eciToEcef(r_sol_ECI, t);

    Eigen::Vector3d r_hat_patch_to_sun = (r_sol_ecef - r_patch_ecef).normalized();

    double r = (r_sol_ecef - r_patch_ecef).norm();

    /** Approximately correct, uses differential area of oblate sphereoid */
    double A_patch = a_earth * cos(lat) * sqrt(pow(a_earth * sin(lat), 2) + pow(b_earth * cos(lat), 2)) * d_lat * d_lon;

    double zenith_cosine = n_hat_ecef.dot(r_hat_patch_to_sun);
    if (zenith_cosine < 0)
    {
        zenith_cosine = 0;
    }
    double zenith_theta = acos(zenith_cosine);

    double AM_frac = AMKastenYoung(zenith_theta);

    SpectraVector_t transmission = atmosphere_transmittance.array().pow(AM_frac);
    

    // Flux is (transmittance^AM) * total power / (4pi r^2), power is flux * patch area * cos(zenith) (where negative zenith = 0)
    SpectraVector_t incidentPower = transmission.cwiseProduct(sun_spectra / (4.0 * M_PI * pow(r, 2))) * A_patch * zenith_cosine;

    return incidentPower;
}


SpectraVector_t getEarthIRPowerOnAlbedoSurface(Eigen::Vector3d LLH, double t)
{
    /** TODO: Fix with actual IR spectrum */
    return SpectraVector_t::Zero();
}


SpectraVector_t reBinAlbedoSpectra(unsigned int lat_i, unsigned int lon_i)
{

}


std::vector<SpectraVector_t> getIncidentAlbedoAndIRPowerSpectra(Eigen::Vector3d LLH, Eigen::Vector3d r_ECI, Eigen::Matrix3d R_eci_from_B, Eigen::Vector3d r_sol_ECI, double t, const SimParams_t& simParams)
{
    Eigen::Vector3d LLH_0 = LLH;
    LLH_0[2] = 0; // Just double check that the height is zero.

    double lat = LLH[0];
    double lon = LLH[1];

    SpectraVector_t incidentPowerPatch = getPowerOnAlbedoSurface(LLH_0, r_sol_ECI, t, simParams);

    SpectraVector_t emittedIRPower = getEarthIRPowerOnAlbedoSurface(LLH_0, t);

    Eigen::Vector3d r_patch_ecef = llh_to_ecef(LLH_0);
    Eigen::Vector3d r_sat_ecef = simParams.ephemeris->eciToEcef(r_ECI, t);

    Eigen::Vector3d r_patch_to_sat_ecef = r_sat_ecef - r_patch_ecef;
    Eigen::Vector3d r_hat_patch_to_sat_ecef = r_patch_to_sat_ecef.normalized();

    Eigen::Vector3d r_hat_patch_to_sat_eci = simParams.ephemeris->ecefToEci(r_hat_patch_to_sat_ecef, t);

    Eigen::Vector3d n_hat_ecef = getNormalToEarthEcef(lat, lon);

    double r = r_patch_to_sat_ecef.norm();

    double zenith_cosine = n_hat_ecef.dot(r_hat_patch_to_sat_ecef);
    if (zenith_cosine < 0)
    {
        zenith_cosine = 0;
    }
    double zenith_theta = acos(zenith_cosine);

    double AM_frac = AMKastenYoung(zenith_theta);

    SpectraVector_t transmission = atmosphere_transmittance.array().pow(AM_frac);

    SpectraVector_t albedo = getAlbedoSpectra(LLH);

    SpectraVector_t powerFromPatch = albedo.cwiseProduct(incidentPowerPatch) + emittedIRPower;

    int n_panels = cubesat.panels.size();
    std::vector<SpectraVector_t> panelIncidentPower(n_panels);

    for (int i = 0; i < n_panels; i++)
    {
        Eigen::Vector3d n_hat_panel_eci = R_eci_from_B * cubesat.panels[i].normal_vector;

        double incident_cosine = -n_hat_panel_eci.dot(r_hat_patch_to_sat_eci);
        if (incident_cosine < 0)
        {
            incident_cosine = 0;
        }

        panelIncidentPower[i] = transmission * (powerFromPatch / (M_PI * pow(r, 2))) * zenith_cosine * incident_cosine * cubesat.panels[i].size;
    }

    return panelIncidentPower;
}