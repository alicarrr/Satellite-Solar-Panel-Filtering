#pragma once

#include <eigen3/Eigen/Dense>
#include "constants.hpp"
#include "spectra.hpp"
#include "simParams.hpp"

#include <vector>
#include <netcdf>

#include <iomanip>
#include <sstream>
#include <string>

#define useLargeHamster false

// Directory containing the HAMSTER DOYxxx.nc files.
// The executable is run from the project root.

#if useLargeHamster == true
    const char hamster_dir[] = "C/lib/HAMSTER";
#else
    const char hamster_dir[] = "C/lib/HAMSTER/hamster_1140_720/data/albedo/hamster";
#endif



// HAMSTER day of year to load.
const unsigned int hamster_day_of_year = 244;

// Maximum memory used to store albedo data after loading.
// The spatial sampling interval is chosen automatically from the file size.
const unsigned int hamster_max_albedo_memory_mb = 512;

extern unsigned int n_spectra_albedo; // Number of wavelengths in the albedo data

inline constexpr double d_lambda = 10 * nm; // metres
inline constexpr double d_lat = 0.25; // degrees
inline constexpr double d_lon = 0.25; // degrees

extern double d_lat_file;
extern double d_lon_file;

extern double lat_min_file;
extern double lat_max_file;
extern double lon_min_file;
extern double lon_max_file;

extern int albedo_spectra_file_len;

extern Eigen::VectorXd latitudes_file;
extern Eigen::VectorXd longitudes_file;
extern std::vector<double> wavelengths_file_nm;

extern netCDF::NcFile albedo_file;

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

void read_albedo_data(unsigned int doy);
void printVariableInformation(const std::string& variable_name, const netCDF::NcFile& file);

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


/**
 * @brief The angle either side of the lat/lon that is directly below the object that is visible.
 * 
 * @param distance The distance from the centre of the object
 * @param r The radius of the object (assuming a uniform sphere)
 * 
 * @returns The angle
 */
double getVisibilityConeAngle(double distance, double r);

void getAlbedoFileParams();

void openAlbedoFile();

int getIndexOfEntry(const Eigen::VectorXd& v, double val);

std::vector<int> getLatLonIndices(Eigen::Vector3d LLH);

SpectraVector_t getAlbedoSpectra(Eigen::Vector3d LLH);

SpectraVector_t getPowerOnAlbedoSurface(Eigen::Vector3d LLH, Eigen::Vector3d r_sol_ECI, double t, const SimParams_t* simParams);

SpectraVector_t getEarthIRPowerOnAlbedoSurface(Eigen::Vector3d LLH, double t, const SimParams_t* simParams);

std::vector<SpectraVector_t> getPatchIncidentAlbedoAndIRPowerSpectra(Eigen::Vector3d LLH, Eigen::Vector3d r_ECI, Eigen::Matrix3d R_eci_from_B, Eigen::Vector3d r_sol_ECI, double t, const SimParams_t* simParams);

double fillAlbedoData(double lambda_nm);

Eigen::Vector3d getLatLonOfSat(Eigen::Vector3d r_ECI, double t, const SimParams_t* simParams);

std::vector<SpectraVector_t> getTotalIncidentAlbedoAndIRPowerSpectra(Eigen::Vector3d r_ECI, Eigen::Matrix3d R_eci_from_B, Eigen::Vector3d r_sol_ECI, double t, const SimParams_t* simParams);