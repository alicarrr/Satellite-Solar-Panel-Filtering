#include "albedo.hpp"
#include <netcdf.h>
#include <netcdf>

#include "transforms.hpp"

#include <atomic>
#include <cmath>
#include <cstdio>
#include <exception>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <thread>
#include <utility>

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

std::vector<float> albedo_tile_spectra;
unsigned int albedo_tile_dimension_0_size = 0;
unsigned int albedo_tile_dimension_1_size = 0;
unsigned int albedo_tile_wavelength_count = 0;

static unsigned int loaded_albedo_day_of_year = 0;

typedef struct
{
    Eigen::Vector3d LLH;
} AlbedoPatchTask_t;

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
static unsigned int find_closest_index(const Eigen::VectorXd& values, double requested_value)
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
static unsigned int find_closest_longitude_index(const Eigen::VectorXd& longitude, double requested_longitude)
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
 * @brief Load the complete HAMSTER albedo cube into contiguous float storage.
 */
void loadAlbedoTileSpectraIntoMemory()
{
    if (albedo_file.isNull())
    {
        throw std::runtime_error("Cannot load HAMSTER albedo data before opening the file.");
    }

    #if useLargeHamster == true
        netCDF::NcVar albedoData = albedo_file.getVar("Black_Sky_Albedo");
    #else
        netCDF::NcGroup maps_group = albedo_file.getGroup("maps");
        netCDF::NcGroup lambertian_group = maps_group.getGroup("lambertian");
        netCDF::NcVar albedoData = lambertian_group.getVar("albedo");
    #endif

    if (albedoData.isNull())
    {
        throw std::runtime_error("Could not find the HAMSTER albedo variable.");
    }

    std::vector<netCDF::NcDim> dimensions = albedoData.getDims();
    if (dimensions.size() != 3)
    {
        throw std::runtime_error("The HAMSTER albedo variable must have three dimensions.");
    }

    #if useLargeHamster == true
        if (dimensions[0].getName() != "Latitude" || dimensions[1].getName() != "Longitude" || dimensions[2].getName() != "Wavelength")
        {
            throw std::runtime_error("Black_Sky_Albedo must be ordered as [Latitude, Longitude, Wavelength].");
        }
    #else
        if (dimensions[0].getName() != "lon" || dimensions[1].getName() != "lat" || dimensions[2].getName() != "spectral")
        {
            throw std::runtime_error("Small-file albedo must be ordered as [lon, lat, spectral].");
        }
    #endif

    if (
        dimensions[0].getSize() > std::numeric_limits<unsigned int>::max() || dimensions[1].getSize() > std::numeric_limits<unsigned int>::max() || dimensions[2].getSize() > std::numeric_limits<unsigned int>::max()
    )
    {
        throw std::runtime_error("A HAMSTER albedo dimension is too large for this build.");
    }

    albedo_tile_dimension_0_size = static_cast<unsigned int>(dimensions[0].getSize());
    albedo_tile_dimension_1_size = static_cast<unsigned int>(dimensions[1].getSize());
    albedo_tile_wavelength_count = static_cast<unsigned int>(dimensions[2].getSize());

    unsigned long long total_value_count = static_cast<unsigned long long>(albedo_tile_dimension_0_size) * static_cast<unsigned long long>(albedo_tile_dimension_1_size) * static_cast<unsigned long long>(albedo_tile_wavelength_count);

    if (
        total_value_count > static_cast<unsigned long long>(albedo_tile_spectra.max_size())
    )
    {
        throw std::length_error( "The HAMSTER albedo cube is too large for std::vector<float>.");
    }

    albedo_tile_spectra.resize(static_cast<std::size_t>(total_value_count));

    std::vector<std::size_t> start(3, 0);
    std::vector<std::size_t> count(3);
    count[0] = albedo_tile_dimension_0_size;
    count[1] = albedo_tile_dimension_1_size;
    count[2] = albedo_tile_wavelength_count;

    // This is the only full-cube NetCDF read. The vector retains the file's
    // native dimension ordering and wavelength remains the fastest index.
    albedoData.getVar(start, count, albedo_tile_spectra.data());

    albedo_spectra_file_len = static_cast<int>(albedo_tile_wavelength_count);

    double memory_mebibytes = static_cast<double>(albedo_tile_spectra.size()) * static_cast<double>(sizeof(float)) / (1024.0 * 1024.0);

    std::cout
        << "Loaded HAMSTER albedo cube into RAM: " << albedo_tile_dimension_0_size << " x " << albedo_tile_dimension_1_size << " x " << albedo_tile_wavelength_count << " floats (" << memory_mebibytes << " MiB)." << std::endl;
}



void read_albedo_data(unsigned int doy)
{
    int width = 3;

    std::ostringstream stream;
    stream << hamster_dir << "/" << "DOY" << std::setfill('0') << std::setw(width) << doy << ".nc";
    std::string file_name = stream.str();

    if (
        !albedo_tile_spectra.empty() && loaded_albedo_day_of_year == doy
    )
    {
        // The startup sequence asks for this day twice. Avoid reloading the
        // approximately 835 MiB cube, but reopen the file if another caller
        // closed it and still needs the coordinate metadata.
        if (albedo_file.isNull())
        {
            albedo_file.open(file_name, netCDF::NcFile::read);
        }

        return;
    }

    std::cout << file_name << std::endl;

    if (!albedo_file.isNull())
    {
        albedo_file.close();
    }

    albedo_file.open(file_name, netCDF::NcFile::read);
    loadAlbedoTileSpectraIntoMemory();

    loaded_albedo_day_of_year = doy;
}


// TODO: could delete this function
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

    for (int i = 0; i < static_cast<int>(dimensions.size()); i++)
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
        // Use for smaller HAMSTER files
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


double wrapValue(double val, double min_val, double max_val) 
{
    double range = max_val - min_val;
    if (range == 0.0) return min_val;
    
    // std::fmod handles remainder, but can return negative numbers if 'val - min_val' is negative
    double result = min_val + std::fmod(val - min_val, range);
    
    // If result is negative, shift it back into the positive range [min_val, max_val)
    if (result < min_val) 
    {
        result += range;
    }
    
    return result;
}


// std::vector<int> getLatLonIndices(Eigen::Vector3d LLH)
// {   
//     // in LLH, lat and lon are in radians, so convert to degrees for Hamster
//     double lat = LLH(0)*RAD2DEG;
//     double lon = LLH(1)*RAD2DEG;

//     lat = wrapValue(lat, lat_min_file, lat_max_file);
//     lon = wrapValue(lon, lon_min_file, lon_max_file);

//     std::vector<int> indices(2);

//     // LARGE FILES: The indices in the file are ordered lat then lon
//     // SMALL FILES: The indices in the file are ordered lon then lat

//     #if useLargeHamster == true
//         indices[0] = getIndexOfEntry(latitudes_file, lat);
//         indices[1] = getIndexOfEntry(longitudes_file, lon);
//     #else
//         indices[1] = getIndexOfEntry(latitudes_file, lat);
//         indices[0] = getIndexOfEntry(longitudes_file, lon);
//     #endif

//     return indices;
// }

std::vector<int> getLatLonIndices(Eigen::Vector3d LLH)
{
    // Patch coordinates are radians; HAMSTER coordinates are degrees.
    double latitude_deg = LLH(0) * RAD2DEG;
    double longitude_deg = LLH(1) * RAD2DEG;

    int latitude_index = getIndexOfEntry(latitudes_file, latitude_deg);

    int longitude_index =
        static_cast<int>(find_closest_longitude_index(longitudes_file, longitude_deg));

    std::vector<int> indices(2);

#if useLargeHamster == true
    // Large file ordering: [latitude, longitude, wavelength]
    indices[0] = latitude_index;
    indices[1] = longitude_index;
#else
    // Small file ordering: [longitude, latitude, wavelength]
    indices[0] = longitude_index;
    indices[1] = latitude_index;
#endif

    return indices;
}


double fillAlbedoData(double lambda_nm)
{
    return 0;
}


/**
 * @brief Copy one raw HAMSTER spectrum from the in-memory float cube.
 *
 * The returned vector converts the stored floats to doubles for the existing
 * spectral resampling function. The global cube is read-only after startup,
 * so this function can safely be called concurrently by worker threads.
 */
static std::vector<double> readRawAlbedoSpectrumFromMemory(const Eigen::Vector3d& LLH)
{
    if (albedo_tile_spectra.empty())
    {
        throw std::runtime_error("HAMSTER albedo data have not been loaded into memory.");
    }

    std::vector<int> indices = getLatLonIndices(LLH);

    if (
        indices[0] < 0 || indices[1] < 0 || static_cast<unsigned int>(indices[0]) >= albedo_tile_dimension_0_size || static_cast<unsigned int>(indices[1]) >= albedo_tile_dimension_1_size
    )
    {
        throw std::out_of_range("HAMSTER latitude/longitude index is outside the RAM cube.");
    }

    unsigned long long spectrum_start =
        (static_cast<unsigned long long>(indices[0]) * static_cast<unsigned long long>(albedo_tile_dimension_1_size) + static_cast<unsigned long long>(indices[1])) * static_cast<unsigned long long>(albedo_tile_wavelength_count);

    unsigned long long spectrum_end = spectrum_start + static_cast<unsigned long long>(albedo_tile_wavelength_count);

    if (spectrum_end > static_cast<unsigned long long>(albedo_tile_spectra.size()))
    {
        throw std::out_of_range("HAMSTER spectrum extends beyond the RAM cube.");
    }

    std::vector<double> raw_albedo_spectrum(albedo_tile_wavelength_count);

    for (unsigned int wavelength_index = 0; wavelength_index < albedo_tile_wavelength_count; wavelength_index++)
    {
        unsigned long long source_index = spectrum_start + static_cast<unsigned long long>(wavelength_index);

        raw_albedo_spectrum[wavelength_index] = static_cast<double>(albedo_tile_spectra[static_cast<std::size_t>(source_index)]);
    }

    return raw_albedo_spectrum;
}


SpectraVector_t getAlbedoSpectra(Eigen::Vector3d LLH)
{
    std::vector<double> raw_albedo_spectrum = readRawAlbedoSpectrumFromMemory(LLH);

    // Need to rebin the data from the File
    std::vector<double> albedoSpectraExtended = resampleAndExtendSpectrum(wavelengths_file_nm, raw_albedo_spectrum, min_wavelength_nm, max_wavelength_nm, d_wavelength_nm, fillAlbedoData);

    SpectraVector_t albedoSpectraVec = Eigen::Map<SpectraVector_t>(albedoSpectraExtended.data(), albedoSpectraExtended.size());

    return albedoSpectraVec;
}


Eigen::Vector3d getNormalToEarthEcef(double lat, double lon)
{
    double e_earth = sqrt(1 - pow((b_earth/a_earth), 2));
    double b_over_a_earth = b_earth/a_earth;
    // Eigen::Vector3d n_hat_ecef = Eigen::Vector3d({b_over_a_earth * cos(lat)*cos(lon), b_over_a_earth * cos(lat)*sin(lon), sin(lat)})/sqrt(1 - pow(e_earth * cos(lat), 2));
    Eigen::Vector3d n_hat_ecef = Eigen::Vector3d({cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)});

    return n_hat_ecef;
}


double patchArea(Eigen::Vector3d LLH)
{
    double lat = LLH[0];
    double lon = LLH[1];

    // double A_patch = a_earth * cos(lat) * sqrt(pow(a_earth * sin(lat), 2) + pow(b_earth * cos(lat), 2)) * (DEG2RAD * d_lat) * (DEG2RAD * d_lon);
    double e2 = 1.0 - pow(b_earth / a_earth, 2);
    double q = 1.0 - e2 * pow(sin(lat), 2);

    double N = a_earth / sqrt(q);
    double M = a_earth * (1.0 - e2) / pow(q, 1.5);

    double A_patch = M * N * cos(lat) * (DEG2RAD * d_lat) * (DEG2RAD * d_lon);

    return A_patch;
}


/**
 * @brief Test the same patch-to-satellite visibility condition used by the
 *        full patch calculation.
 */
static bool isAlbedoPatchVisible(const Eigen::Vector3d& LLH, const Eigen::Vector3d& r_sat_ecef)
{
    Eigen::Vector3d LLH_0 = LLH;
    LLH_0[2] = 0;

    Eigen::Vector3d r_patch_ecef = llh_to_ecef(LLH_0);
    Eigen::Vector3d r_hat_patch_to_sat_ecef = (r_sat_ecef - r_patch_ecef).normalized();

    Eigen::Vector3d n_hat_ecef = getNormalToEarthEcef(LLH_0[0], LLH_0[1]);

    double zenith_cosine = n_hat_ecef.dot(r_hat_patch_to_sat_ecef);

    return zenith_cosine >= 0.0;
}


/**
 * @brief Calculate one patch using data that were loaded before worker threads
 *        started.
 *
 * This helper contains no NetCDF, HDF5, or CSPICE calls. All inputs are either
 * immutable global spectra or values owned by the current worker.
 */
static std::vector<SpectraVector_t> getPatchIncidentAlbedoAndIRPowerSpectraFromPreloadedData(const Eigen::Vector3d& LLH, const Eigen::Vector3d& r_sat_ecef, const Eigen::Vector3d& r_sol_ecef, const Eigen::Matrix3d& R_eci_from_ecef, const Eigen::Matrix3d& R_eci_from_B, const SpectraVector_t& earth_ir_power_per_square_metre, const SpectraVector_t& albedo)
{
    Eigen::Vector3d LLH_0 = LLH;
    LLH_0[2] = 0; // height = 0, i.e. at the Earth's surface.

    // patch lat, lon, height
    double lat = LLH_0[0];
    double lon = LLH_0[1];

    int n_panels = cubesat.panels.size();
    std::vector<SpectraVector_t> panelIncidentPower(n_panels, SpectraVector_t::Zero());

    // use ECEF as the global frame until the patch-to-satellite vector is transformed into the ECI frame for the panel incidence calculation
    Eigen::Vector3d r_patch_ecef = llh_to_ecef(LLH_0);
    Eigen::Vector3d n_hat_ecef = getNormalToEarthEcef(lat, lon);

    Eigen::Vector3d r_patch_to_sat_ecef = r_sat_ecef - r_patch_ecef;

    double r_patch_to_sat = r_patch_to_sat_ecef.norm();
    Eigen::Vector3d r_hat_patch_to_sat_ecef = r_patch_to_sat_ecef.normalized();

    double satellite_zenith_cosine = n_hat_ecef.dot(r_hat_patch_to_sat_ecef);

    // Check if the patch is visible to the satellite. If not, return zero power for all panels.
    if (satellite_zenith_cosine < 0.0)
    {
        return panelIncidentPower;
    }

    // area of the patch on Earth surface 
    double A_patch = patchArea(LLH_0);

    SpectraVector_t incidentPowerPatch = SpectraVector_t::Zero();

    Eigen::Vector3d r_patch_to_sun_ecef = r_sol_ecef - r_patch_ecef;
    double r_patch_to_sun = r_patch_to_sun_ecef.norm();
    Eigen::Vector3d r_hat_patch_to_sun_ecef = r_patch_to_sun_ecef.normalized();

    double solar_zenith_cosine = n_hat_ecef.dot(r_hat_patch_to_sun_ecef);

    // check if the patch is illuminated by the sun. If not, incident power from the sun is zero.
    if (solar_zenith_cosine >= 0.0)
    {
        // if visible
        double solar_zenith_theta = acos(std::clamp(solar_zenith_cosine, -1.0, 1.0));     // angle between the sun and the patch normal (in rads)
        double solar_AM_frac = AMKastenYoung(solar_zenith_theta);   // air mass fraction for incident solar radiation

        SpectraVector_t solar_transmission = atmosphere_transmittance;
        // scale the transmissivity by the air mass fraction to account for the increased path length through the atmosphere at varied zenith angles
        solar_transmission = solar_transmission.array().pow(solar_AM_frac);
        
        // incident power on the patch from the sun, accounting for the distance to the sun, the area of the patch, and the angle of incidence of sunlight on the patch
        incidentPowerPatch = solar_transmission.cwiseProduct(sun_spectra / (4.0 * M_PI * pow(r_patch_to_sun, 2))) * A_patch * solar_zenith_cosine;
    }

    // // calculate the emitted IR power from the patch, which is independent of the sun's position and only depends on the patch's temperature and area
    SpectraVector_t emittedIRPower = earth_ir_power_per_square_metre * A_patch;

    // angle between the patch normal and the satellite panel normal.
    double satellite_zenith_theta = acos(std::clamp(satellite_zenith_cosine, -1.0, 1.0));
    double satellite_AM_frac = AMKastenYoung(satellite_zenith_theta); // air mass fraction for radiation from the patch to the satellite

    SpectraVector_t satellite_transmission = atmosphere_transmittance;
    satellite_transmission = satellite_transmission.array().pow(satellite_AM_frac);

    SpectraVector_t powerFromPatch = albedo.cwiseProduct(incidentPowerPatch) + emittedIRPower;

    Eigen::Vector3d r_hat_patch_to_sat_eci = R_eci_from_ecef * r_hat_patch_to_sat_ecef;

    for (int i = 0; i < n_panels; i++)
    {
        Eigen::Vector3d n_hat_panel_eci = R_eci_from_B * cubesat.panels[i].normal_vector;

        double incident_cosine = -n_hat_panel_eci.dot(r_hat_patch_to_sat_eci);

        if (incident_cosine <= 0.0)
        {
            continue;
        }

        panelIncidentPower[i] = satellite_transmission.cwiseProduct(powerFromPatch / (M_PI * pow(r_patch_to_sat, 2))) * satellite_zenith_cosine * incident_cosine * cubesat.panels[i].size;
    }
    

    return panelIncidentPower;
}


Eigen::Vector3d getLatLonOfSat(Eigen::Vector3d r_ECI, double t, const SimParams_t* simParams)
{
    Eigen::Vector3d r_ecef = simParams->ephemeris->eciToEcef(r_ECI, t);

    Eigen::Vector3d LLH = ecef_to_llh(r_ecef);
    Eigen::Vector3d LLH_0 = LLH;
    LLH_0[2] = 0;

    return LLH_0;
}


std::vector<SpectraVector_t> getTotalIncidentAlbedoAndIRPowerSpectra(Eigen::Vector3d r_ECI, Eigen::Matrix3d R_eci_from_B, Eigen::Vector3d r_sol_ECI, double t, const SimParams_t* simParams)
{
    int n_panels = cubesat.panels.size();
    std::vector<SpectraVector_t> panelIncidentPower(n_panels, SpectraVector_t::Zero());

    // This is the difference in lat/lon either side of the direct projected LLH
    double theta_vis = getVisibilityConeAngle(r_ECI.norm(), b_earth); // Angle formed on Earth by cone
    // std::cout << "theta_vis: " << RAD2DEG * theta_vis << std::endl;

    // This refreshes the ephemeris frame cache for t. All remaining ephemeris
    // conversions are also completed before any worker is started.
    Eigen::Vector3d LLH_sat = getLatLonOfSat(r_ECI, t, simParams);

    int N_lat = 180.0 / d_lat;//theta_vis / (DEG2RAD * d_lat);
    int N_lon = 360.0 / d_lon;//theta_vis / (DEG2RAD * d_lon);

    // std::cout << "N_lat: " << N_lat << ", N_lon: " << N_lon << std::endl;

    // int patch_offset_count = N_lat * N_lon;
    // int candidate_task_count = 4 * patch_offset_count;

    int candidate_task_count = N_lat * N_lon;

    if (candidate_task_count <= 0)
    {
        return panelIncidentPower;
    }

    // Resolve every library-dependent quantity on the calling thread. No
    // worker below accesses NetCDF, HDF5, or CSPICE.
    Eigen::Vector3d r_sat_ecef = simParams->ephemeris->eciToEcef(r_ECI, t);
    Eigen::Vector3d r_sol_ecef = simParams->ephemeris->eciToEcef(r_sol_ECI, t);

    Eigen::Matrix3d R_eci_from_ecef;
    R_eci_from_ecef.col(0) = simParams->ephemeris->ecefToEci(Eigen::Vector3d::UnitX(), t);
    R_eci_from_ecef.col(1) = simParams->ephemeris->ecefToEci(Eigen::Vector3d::UnitY(), t);
    R_eci_from_ecef.col(2) = simParams->ephemeris->ecefToEci(Eigen::Vector3d::UnitZ(), t);

    // This wavelength-dependent part of Earth IR is identical for every
    // patch. Workers only multiply it by their patch area.
    SpectraVector_t earth_ir_power_per_square_metre = SpectraVector_t::Zero();
    for (int i = 0; i < earth_ir_power_per_square_metre.size(); i++)
    {
        earth_ir_power_per_square_metre(i) = blackBodySpectrum(T_earth, spectra_wavelengths_nm(i)) * d_wavelength_nm;
    }

    if (albedo_tile_spectra.empty())
    {
        throw std::runtime_error("HAMSTER albedo data must be loaded before calculating albedo.");
    }

    // Tasks now store only their patch coordinate. Workers copy the required
    // 211 float values from the immutable RAM cube and convert them to double.
    // No per-patch NetCDF reads or raw-spectrum preload are required here.
    std::vector<AlbedoPatchTask_t> patch_tasks;
    patch_tasks.reserve(candidate_task_count);

    for (int candidate_task_index = 0; candidate_task_index < candidate_task_count; candidate_task_index++)
    {
        int i = candidate_task_index / N_lon;
        int j = candidate_task_index % N_lon;

        Eigen::Vector3d LLH_ij(
            (-90.0 + (i + 0.5) * d_lat) * DEG2RAD,
            (-180.0 + (j + 0.5) * d_lon) * DEG2RAD,
            0.0
        );

        if (!isAlbedoPatchVisible(LLH_ij, r_sat_ecef))
        {
            continue;
        }

        AlbedoPatchTask_t patch_task;
        patch_task.LLH = LLH_ij;

        patch_tasks.push_back(std::move(patch_task));
    }

    int task_count = static_cast<int>(patch_tasks.size());
    if (task_count == 0)
    {
        return panelIncidentPower;
    }

    unsigned int worker_count = std::thread::hardware_concurrency();
    if (worker_count == 0)
    {
        worker_count = 1;
    }
    if (worker_count > static_cast<unsigned int>(task_count))
    {
        worker_count = static_cast<unsigned int>(task_count);
    }

    // Each worker owns its accumulator, so the expensive vector additions do
    // not need a mutex and cannot race with another worker.
    std::vector<std::vector<SpectraVector_t> > worker_sums;
    worker_sums.reserve(worker_count);
    for (unsigned int worker_index = 0; worker_index < worker_count; worker_index++)
    {
        worker_sums.push_back(
            std::vector<SpectraVector_t>(n_panels, SpectraVector_t::Zero())
        );
    }

    // Workers dynamically claim the next visible, RAM-backed patch. This
    // balances differences in the spectral maths between patches.
    std::atomic<int> next_task(0);
    std::atomic<bool> stop_workers(false);

    std::exception_ptr worker_exception;
    std::mutex worker_exception_mutex;

    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    try
    {
        for (unsigned int worker_index = 0; worker_index < worker_count; worker_index++)
        {
            workers.push_back(std::thread([&, worker_index]()
                    {
                        try
                        {
                            while (!stop_workers.load(std::memory_order_relaxed))
                            {
                                int task_index = next_task.fetch_add(1, std::memory_order_relaxed);

                                if (task_index >= task_count)
                                {
                                    break;
                                }

                                const AlbedoPatchTask_t& patch_task = patch_tasks[task_index];

                                std::vector<double> raw_albedo_spectrum = readRawAlbedoSpectrumFromMemory(patch_task.LLH);

                                std::vector<double> albedo_spectrum_extended = resampleAndExtendSpectrum(wavelengths_file_nm, raw_albedo_spectrum, min_wavelength_nm, max_wavelength_nm, d_wavelength_nm, fillAlbedoData);

                                SpectraVector_t albedo = Eigen::Map<SpectraVector_t>(albedo_spectrum_extended.data(), albedo_spectrum_extended.size());

                                std::vector<SpectraVector_t> panelIncidentPower_ij = getPatchIncidentAlbedoAndIRPowerSpectraFromPreloadedData(patch_task.LLH, r_sat_ecef, r_sol_ecef, R_eci_from_ecef, R_eci_from_B, earth_ir_power_per_square_metre, albedo);

                                for (int k = 0; k < n_panels; k++)
                                {
                                    worker_sums[worker_index][k] += panelIncidentPower_ij[k];
                                }
                            }
                        }
                        catch (...)
                        {
                            stop_workers.store(true, std::memory_order_relaxed);

                            std::lock_guard<std::mutex> exception_lock(worker_exception_mutex);

                            if (!worker_exception)
                            {
                                worker_exception = std::current_exception();
                            }
                        }
                    }
                )
            );
        }
    }
    catch (...)
    {
        stop_workers.store(true, std::memory_order_relaxed);

        for (unsigned int worker_index = 0; worker_index < static_cast<unsigned int>(workers.size()); worker_index++)
        {
            if (workers[worker_index].joinable())
            {
                workers[worker_index].join();
            }
        }

        throw;
    }

    for (unsigned int worker_index = 0; worker_index < static_cast<unsigned int>(workers.size()); worker_index++)
    {
        workers[worker_index].join();
    }

    if (worker_exception)
    {
        std::rethrow_exception(worker_exception);
    }

    // No workers are running now, so this final reduction also needs no mutex.
    for (unsigned int worker_index = 0; worker_index < worker_count; worker_index++)
    {
        for (int k = 0; k < n_panels; k++)
        {
            panelIncidentPower[k] += worker_sums[worker_index][k];
        }
    }

    return panelIncidentPower;
}



// OLD FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------------------


// SpectraVector_t getEarthIRPowerOnAlbedoSurface(Eigen::Vector3d LLH, double t, const SimParams_t* simParams)
// {
//     SpectraVector_t IRPatchPower = SpectraVector_t::Zero();
//     // Calculate blackbody power at each wavelength in W/m^2/nm
//     for (int i = 0; i < IRPatchPower.size(); i++)
//     {
//         IRPatchPower(i) = blackBodySpectrum(T_earth, spectra_wavelengths_nm(i));
//     }

//     // Multiply by d_wavelength_nm to get power in W/m^2
//     IRPatchPower *= d_wavelength_nm;

//     // Multiply by patch area to get power in W
//     double patchArea_m2 = patchArea(LLH);
//     IRPatchPower *= patchArea_m2;

//     return IRPatchPower;
// }




// SpectraVector_t getPowerOnAlbedoSurface(Eigen::Vector3d LLH, Eigen::Vector3d r_sol_ECI, double t, const SimParams_t* simParams)
// {
//     double lat = LLH[0];
//     double lon = LLH[1];

//     SpectraVector_t incidentPower = SpectraVector_t::Zero();

//     Eigen::Vector3d r_patch_ecef = llh_to_ecef(LLH);

//     /** Normal vector for oblate sphereoid */
//     Eigen::Vector3d n_hat_ecef = getNormalToEarthEcef(lat, lon);

//     Eigen::Vector3d r_sol_ecef = simParams->ephemeris->eciToEcef(r_sol_ECI, t);

//     Eigen::Vector3d r_hat_patch_to_sun = (r_sol_ecef - r_patch_ecef).normalized();

//     double r = (r_sol_ecef - r_patch_ecef).norm();

//     /** Approximately correct, uses differential area of oblate sphereoid */
//     double A_patch = patchArea(LLH);

//     double zenith_cosine = n_hat_ecef.dot(r_hat_patch_to_sun);
//     if (zenith_cosine < 0)
//     {
//         zenith_cosine = 0;

//         // No need to do the rest of the calculations, this is not in sunlight
//         return incidentPower;
//     }
//     double zenith_theta = acos(zenith_cosine);

//     double AM_frac = AMKastenYoung(zenith_theta);

//     SpectraVector_t transmission = atmosphere_transmittance;
//     transmission = transmission.array().pow(AM_frac);

//     // Flux is (transmittance^AM) * total power / (4pi r^2), power is flux * patch area * cos(zenith) (where negative zenith = 0)
//     incidentPower = transmission.cwiseProduct(sun_spectra / (4.0 * M_PI * pow(r, 2))) * A_patch * zenith_cosine;

//     return incidentPower;
// }




// std::vector<SpectraVector_t> getPatchIncidentAlbedoAndIRPowerSpectra(Eigen::Vector3d LLH, Eigen::Vector3d r_ECI, Eigen::Matrix3d R_eci_from_B, Eigen::Vector3d r_sol_ECI, double t, const SimParams_t* simParams)
// {
//     Eigen::Vector3d LLH_0 = LLH;
//     LLH_0[2] = 0; // Just double check that the height is zero.

//     double lat = LLH[0];
//     double lon = LLH[1];


//     int n_panels = cubesat.panels.size();
//     std::vector<SpectraVector_t> panelIncidentPower(n_panels, SpectraVector_t::Zero());
    

//     SpectraVector_t incidentPowerPatch = getPowerOnAlbedoSurface(LLH_0, r_sol_ECI, t, simParams);

//     SpectraVector_t emittedIRPower = getEarthIRPowerOnAlbedoSurface(LLH_0, t, simParams);

//     Eigen::Vector3d r_patch_ecef = llh_to_ecef(LLH_0);
//     Eigen::Vector3d r_sat_ecef = simParams->ephemeris->eciToEcef(r_ECI, t);

//     Eigen::Vector3d r_patch_to_sat_ecef = r_sat_ecef - r_patch_ecef;
//     Eigen::Vector3d r_hat_patch_to_sat_ecef = r_patch_to_sat_ecef.normalized();

//     Eigen::Vector3d r_hat_patch_to_sat_eci = simParams->ephemeris->ecefToEci(r_hat_patch_to_sat_ecef, t);

//     Eigen::Vector3d n_hat_ecef = getNormalToEarthEcef(lat, lon);

//     double r = r_patch_to_sat_ecef.norm();

//     double zenith_cosine = n_hat_ecef.dot(r_hat_patch_to_sat_ecef);
//     if (zenith_cosine < 0)
//     {
//         zenith_cosine = 0;

//         // No need to do the rest of the calculations, this is not visible
//         return panelIncidentPower;
//     }
//     double zenith_theta = acos(zenith_cosine);

//     double AM_frac = AMKastenYoung(zenith_theta);

//     SpectraVector_t transmission = atmosphere_transmittance;
//     transmission = transmission.array().pow(AM_frac);

//     SpectraVector_t albedo = getAlbedoSpectra(LLH);

//     SpectraVector_t powerFromPatch = albedo.cwiseProduct(incidentPowerPatch) + emittedIRPower;

    

//     for (int i = 0; i < n_panels; i++)
//     {
//         Eigen::Vector3d n_hat_panel_eci = R_eci_from_B * cubesat.panels[i].normal_vector;

//         double incident_cosine = -n_hat_panel_eci.dot(r_hat_patch_to_sat_eci);
//         if (incident_cosine <= 0)
//         {
//             // If this panel has no flux on it then just skip the maths
//             continue;
//         }

//         panelIncidentPower[i] = transmission.cwiseProduct(powerFromPatch / (M_PI * pow(r, 2))) * zenith_cosine * incident_cosine * cubesat.panels[i].size;
//     }

//     return panelIncidentPower;
// }