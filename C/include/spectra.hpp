#pragma once

#include "constants.hpp"
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