#pragma once

#include "constants.hpp"
#include "satParams.hpp"
#include <eigen3/Eigen/Dense>
#include <string>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <netcdf>


inline constexpr double min_wavelength_nm = 10;
inline constexpr double max_wavelength_nm = 50e3;
inline constexpr double d_wavelength_nm = 10;

// Number of frequencies with data for input and output spectra
inline constexpr int N_samples_spectra = ceil((max_wavelength_nm - min_wavelength_nm)/d_wavelength_nm);

typedef Eigen::Vector<double, N_samples_spectra> SpectraVector_t;

extern SpectraVector_t spectra_wavelengths_nm;
extern SpectraVector_t sun_spectra; // In watts, binned 
extern SpectraVector_t atmosphere_transmittance; // Dimensionless, put to power of AM
extern SpectraVector_t earth_ir_spectra; // In watts, binned

inline const std::string atmosphere_file_path = "dat files/earth_atmosphere_transmissivity.nc";
inline const std::string AM_spectra_file_name = "dat files/astmg173/ASTMG173.csv";
inline const std::string sun_spectrum_path = "dat files/kurudz_full.dat";

extern netCDF::NcFile atmosphere_transmissivity_file;

/**
 * @brief Spectrum with corresponding wavelength values.
 */
typedef struct {
    Eigen::VectorXd wavelengths;
    Eigen::VectorXd values;
} Spectrum_t;


typedef struct {
    std::vector<double> wavelengths_nm;
    std::vector<double> spectra;
} FileSpectrum_t;

/**
 * @brief Struct to hold the final results of the simulation for each panel at each time step.
 */
typedef struct {
    // store time step and which panel this is for
    double time; // Time step in seconds
    int panel_index; // Index of the panel in the cubesat panels array
    Spectrum_t incident_sun_spectrum; // Incident sun spectrum on the panel at each time step
    Spectrum_t Earth_IR_spectrum; // Earth IR spectrum on the panel at each time step
    Spectrum_t Earth_Albedo_spectrum; // Earth Albedo spectrum on the panel at each time step
} PanelResults_t;

/**
 * @brief Struct to hold the final results at each time step for the entire satellite.
 */
typedef struct {
    double time; // Time step in seconds
    CubesatState_t cubesat_state; // Cubesat state at each time step
    std::vector<PanelResults_t> panel_results; // Vector of incident fluxes on each panel for each time step
} SatelliteResults_t;


/**
 * @brief Read a two-column spectrum .dat file.
 *
 * @param file_name Name of spectrum file.
 *
 * @return Spectrum containing wavelength and spectrum values.
 */
Spectrum_t readSpectrumDat(const std::string& file_name);


/**
 * @brief Bin a spectrum to a uniformly spaced wavelength range.
 *
 * The input wavelengths do not need to be uniformly spaced. For a bin the values are interpolated, accounting for wavelength spacing. The output spectrum is uniformly spaced with the specified step size. Each output bin contains the weighted mean value over that wavelength interval.
 *
 * @param input Input spectrum and corresponding wavelengths.
 * @param min_wavelength Minimum output wavelength in nm.
 * @param max_wavelength Maximum output wavelength in nm.
 * @param wavelength_step_size Output bin width in nm.
 *
 * @return Binned spectrum values and wavelength array
 */
Spectrum_t binSpectrum(const Spectrum_t& input, double min_wavelength, double max_wavelength, double wavelength_step_size);

/**
 * @brief rotate a vector from the body frame to the global frame using the rotation matrix.
 * 
 * @param vector_B Vector in body frame.
 * @param R_G_from_B Rotation matrix from body frame to global frame.
 * 
 * @return Vector in global frame.
 */
Eigen::Vector3d rotateVectorFromBToG(Eigen::Vector3d vector_B, Eigen::Matrix3d R_G_from_B);


/**
 * @brief Convert panel position (in body frame relative to GC) to global frame using the rotation matrix.
 * 
 * @param r_panel_B             Panel position in body frame relative to the satellite center of mass.
 * @param R_G_from_B            Rotation matrix from body frame to global frame.
 * @param S_sat_pos_GC_G        Satellite geometric center absolute position in global frame.
 * 
 * @return Panel position in global frame.
 */
Eigen::Vector3d panelPosInG(Eigen::Vector3d r_panel_B, Eigen::Matrix3d R_G_from_B, Eigen::Vector3d S_sat_pos_GC_G);


/**
 * @brief Calculate incident sun spectrum at a surface (in space) given the sun spectrum, and the relative Sun and satellite positions in the global frame. The panel normal vector in the global frame is used.
 * 
 * @param input Input sun spectrum and corresponding wavelengths.
 * @param sun_pos_G Distance from the sun in AU.
 * @param panel_pos_G Panel position in global frame.
 * @param panel_normal_vect_B Satellite panel normal vector in body frame.
 * 
 * @return Scaled spectrum values and wavelength array as well as angle between the sun ray and the satellite panel normal vector.
 */
Spectrum_t calculateIncidentSunSpectrumOnPanel(const Spectrum_t& input, Eigen::Vector3d sun_pos_G, Eigen::Vector3d panel_pos_G, Eigen::Vector3d panel_normal_vect_B);


/**
 * @brief Calculate panel normal vector direction in the global frame given the satellite position and rotation matrix.
 * 
 * @param S_sat_pos_GC_G Satellite position in global frame.
 * @param R_G_from_B Rotation matrix from body frame to global frame.
 * @param panel_normal_vect_B Satellite panel normal vector in body frame.
 * 
 * @return Panel normal vector in global frame.
 */
Eigen::Vector3d PanelNormVectInG(Eigen::Vector3d S_sat_pos_GC_G, Eigen::Matrix3d R_G_from_B, Eigen::Vector3d panel_normal_vect_B);

/**
 * @brief Calculate the optical thickness which scales to the power of m where m is the air-mass fraction that the light passes through.
 * @param AM0_spectrum_binned Binned solar spectrum at the top of the atmosphere (AM0).
 * @param AM1_5_spectrum_binned Binned solar spectrum at the Earth's surface (AM1.5).
 * @return Optical thickness spectrum.
 */
Spectrum_t opticalTransmittanceAtmosphereFromData(const Spectrum_t& AM0_spectrum_binned, const Spectrum_t& AM1_5_spectrum_binned);

/**
 * @brief An approximate surface-reflected direct-beam contribution to directional TOA spectral radiance, using Beer–Lambert air-mass scaling and Lambertian reflection.
 * 
 * Calculate the top of atmosphere (TOA) spectrum of incident light from an albedo patch on Earth. This accounts for an approximate atmospheric attenuation based on decay of the solar spectrum from AM0 to AM1.5. The albedo spectrum is then scaled by the cosine of the angle between the sun and the panel normal vector.
 * 
 * @param AM0_spectrum_binned Binned solar spectrum at the top of the atmosphere (AM0).
 * @param AM1_5_spectrum_binned Binned solar spectrum at the Earth's surface (AM1.5).
 * @param albedo_spectrum_at_loc Albedo spectrum at the location of the albedo patch on Earth.
 * @param albedo_patch_pos_G Position of the albedo patch in global frame.
 * @param albedo_patch_norm_vector Normal vector of the albedo patch in global frame.
 * @param panel_pos_G Position of the satellite panel in global frame.
 * @param sun_pos_G Position of the sun in global frame.
 */
Spectrum_t albedoRadianceTOA(const Spectrum_t& transmission_ratio_AM1_5, const Spectrum_t& AM0, Spectrum_t albedo_spectrum_at_loc, Eigen::Vector3d albedo_patch_pos_G, Eigen::Vector3d albedo_patch_norm_vector, Eigen::Vector3d panel_pos_G, Eigen::Vector3d sun_pos_G);


Spectrum_t readASTMGColumn(
    const std::string& file_name,
    const int spectrum_column);

Spectrum_t readASTMGFileForAM0(const std::string& file_name);

Spectrum_t readASTMGFileForAM1_5(const std::string& file_name);

Spectrum_t atmosphereTransissionRatio(const Spectrum_t& AM0_spectrum_binned, const Spectrum_t& AM1_5_spectrum_binned);

Spectrum_t interpolateSpectrum(Spectrum_t input, Eigen::VectorXd new_wavelengths);

double cosineOfAngleBetweenVectors(Eigen::Vector3d vector1, Eigen::Vector3d vector2);

double AMKastenYoung(double zenith_angle_rad);

Spectrum_t scaledTransmissivityForAM(double air_mass, const Spectrum_t& atmosphereTransmittance);

// Note this returns in W/m^2/nm
Spectrum_t blackBodySpectrum(double temperature_K, const Spectrum_t& wavelength_spectrum);

Spectrum_t EarthIRRadianceTOA(const Spectrum_t& transmission_ratio_AM1_5, const Spectrum_t& AM0, Eigen::Vector3d patch_pos_G, Eigen::Vector3d panel_pos_G);



/**
 * @brief Conservatively rebins and extends a point-sampled spectral
 * irradiance dataset.
 *
 * The input spectral irradiances are treated as point values connected by
 * piecewise-linear segments. Each output value is calculated by integrating
 * this reconstructed spectrum over the output bin and dividing by the bin
 * width. Consequently, the returned values are bin-average spectral
 * irradiances rather than point samples.
 *
 * The requested minimum and maximum wavelengths are treated as the outer
 * edges of the output domain. The output wavelength values are the centres
 * of bins having width step_wavelength_nm. The output range must therefore
 * contain an integer number of bins.
 *
 * Source data outside the requested output domain are discarded without
 * renormalisation. Regions inside the requested output domain but outside
 * the source dataset are filled using domain_extension_function. Passing
 * nullptr uses zero outside the source domain.
 *
 * @param wavelengths_nm Strictly increasing source wavelengths in nm.
 *
 * @param spectral_irradiances Source spectral irradiance point values. These
 * must have the same length as wavelengths_nm and must be expressed per nm.
 *
 * @param min_wavelength_nm Lower outer edge of the output domain in nm.
 *
 * @param max_wavelength_nm Upper outer edge of the output domain in nm.
 *
 * @param step_wavelength_nm Width of every output wavelength bin in nm.
 *
 * @param domain_extension_function Pointer to a function that takes a
 * wavelength in nm and returns the extended spectral irradiance at that
 * wavelength. The returned value must use the same units as the source
 * spectral irradiances. A nullptr selects zero extension.
 *
 * @return A std::vector containing the average
 * spectral irradiance of every output bin.
 *
 * @throws std::invalid_argument If the input vectors are invalid, the
 * wavelengths are not strictly increasing, the output parameters are
 * invalid, or the output range is not divisible by the requested step.
 *
 * @throws std::domain_error If the extension function returns a non-finite
 * value.
 */
std::vector<double> rebinSpectrumConservingFlux(
    const std::vector<double>& wavelengths_nm,
    const std::vector<double>& spectral_irradiances,
    double min_wavelength_nm,
    double max_wavelength_nm,
    double step_wavelength_nm,
    double (*domain_extension_function)(double));



/**
 * @brief Resamples and extends spectral data without conserving integrated
 * area.
 *
 * The minimum and maximum wavelengths are treated as the outer edges of the
 * output domain. One output value is calculated for each uniformly spaced
 * bin.
 *
 * If one or more source data points lie inside an output bin, their values
 * are combined using an arithmetic mean. If no source points lie inside the
 * bin, the value is calculated by linearly interpolating the source data at
 * the bin centre.
 *
 * If the bin centre lies outside the source wavelength domain,
 * domain_extension_function is evaluated at the bin centre. Passing nullptr
 * uses zero outside the source domain.
 *
 * This function does not preserve the integrated area of the input spectrum.
 * It is suitable for dimensionless spectral properties such as
 * transmissivity, reflectivity, absorptivity or quantum efficiency.
 *
 * @param wavelengths_nm Strictly increasing source wavelengths in nm.
 *
 * @param spectral_values Source spectral values corresponding to each
 * wavelength.
 *
 * @param min_wavelength_nm Lower outer edge of the output domain in nm.
 *
 * @param max_wavelength_nm Upper outer edge of the output domain in nm.
 *
 * @param step_wavelength_nm Width of each output wavelength bin in nm.
 *
 * @param domain_extension_function Pointer to a function that takes a
 * wavelength in nm and returns the extended spectral value at that
 * wavelength. A nullptr selects zero extension.
 *
 * @return A vector containing one resampled spectral value for each output
 * wavelength bin.
 */
std::vector<double> resampleAndExtendSpectrum(
    const std::vector<double>& wavelengths_nm,
    const std::vector<double>& spectral_values,
    double min_wavelength_nm,
    double max_wavelength_nm,
    double step_wavelength_nm,
    double (*domain_extension_function)(double));


/**
 * @brief Calculates the centre wavelength of each uniformly spaced bin.
 *
 * The minimum and maximum wavelengths are treated as the outer edges of the
 * complete wavelength domain. The range must contain an integer number of
 * bins having width step_wavelength_nm.
 *
 * @param min_wavelength_nm Lower outer edge of the wavelength domain in nm.
 *
 * @param max_wavelength_nm Upper outer edge of the wavelength domain in nm.
 *
 * @param step_wavelength_nm Width of each wavelength bin in nm.
 *
 * @return A vector containing the centre wavelength of each bin in nm.
 */
std::vector<double> calculateBinCentreWavelengths(double min_wavelength_nm, double max_wavelength_nm, double step_wavelength_nm);


/** 
 * @brief Calculates the sun spectrum if the dataset is not defined.
 * 
 * @param wavelength_nm The wavelength at which to evaluate the spectral irradiance (in W/m^2/nm), measured in nm
 * 
 * @return The spectral irradiance in W/m^2/nm
*/
double sun_spectra_extend(double wavelength_nm);


/** 
 * @brief Calculates the Earth IR spectrum if the dataset is not defined.
 * 
 * @param wavelength_nm The wavelength at which to evaluate the spectral irradiance (in W/m^2/nm), measured in nm
 * 
 * @return The spectral irradiance in W/m^2/nm
*/
double earth_ir_spectra_extend(double wavelength_nm);


/** 
 * @brief Calculates the AM0 spectrum if the dataset is not defined.
 * 
 * @param wavelength_nm The wavelength at which to evaluate the spectral irradiance (in W/m^2/nm), measured in nm
 * 
 * @return The spectral irradiance in W/m^2/nm
*/
double am_0_spectra_extend(double wavelength_nm);


/** 
 * @brief Calculates the AM1.5 spectrum if the dataset is not defined.
 * 
 * @param wavelength_nm The wavelength at which to evaluate the spectral irradiance (in W/m^2/nm), measured in nm
 * 
 * @return The spectral irradiance in W/m^2/nm
*/
double am_1_5_spectra_extend(double wavelength_nm);


/** 
 * @brief Calculates the Earth atmosphere transmissivity if the dataset is not defined
 * 
 * @param wavelength_nm The wavelength at which to evaluate the transmissivity, measured in nm
 * 
 * @return The transmissivity
*/
double transmissivity_extend(double wavelength_nm);


void setWavelengthArray();

void getSolarSpectrumFromFile(const std::string& fileName);

void getTransmittanceFromFile(const std::string& fileName);

void getEarthIRSpectrumFromFile(const std::string& fileName);
