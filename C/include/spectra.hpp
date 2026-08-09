#pragma once

#include <eigen3/Eigen/Dense>
#include "constants.hpp"
#include <string>

// all spectrum vector types have the same length
#define VectorSpectrum_t Eigen::Vector<double, N_samples_spectra>

/**
 * @brief Get sun spectrum from csv. This will have the number of steps given by the constant N_samples_spectra, and given wavelength range.
 * @param csv_name Name of file with sun spectra
 * @param max_wavelength Maximum wavelength to consider. Default is 2500.0 nm
 * @param min_wavelength Minimum wavelength to consider. Default is 300.0 nm
 * @return vector of sun intensity over wavelength
 */
VectorSpectrum_t getSunSpectrum(std::string csv_name,
                                 double max_wavelength = 2500.0,
                                 double min_wavelength = 300.0);

/**
 * 
 */