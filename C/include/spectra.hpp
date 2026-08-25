#pragma once

#include "constants.hpp"
#include "satParams.hpp"
#include <eigen3/Eigen/Dense>
#include <string>

// all spectrum vector types have the same length
#define VectorSpectrum_t Eigen::Vector<double, N_samples_spectra>

/**
 * @brief Spectrum with corresponding wavelength values.
 */
typedef struct {
    Eigen::VectorXd wavelengths;
    Eigen::VectorXd values;
} Spectrum_t;

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