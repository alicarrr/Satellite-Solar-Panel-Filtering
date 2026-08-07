#pragma once

#include <eigen3/Eigen/Dense>
#include "constants.hpp"
#include <string>

// all spectrum vector types have the same length
#define VectorSpectrum_t Eigen::Vector<double, N_samples_spectra>

/**
 * @brief Get sun spectrum from csv
 * @param csv_name Name of file with sun spectra
 * @return vector of sun intensity over wavelength
 */
VectorSpectrum_t getSunSpectrum(std::string csv_name);
