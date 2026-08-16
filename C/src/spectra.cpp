#include "spectra.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>


Spectrum_t readSpectrumDat(const std::string& file_name)
{
    std::ifstream file(file_name);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + file_name);
    }

    std::string line;
    int num_samples = 0;

    // First pass: count data rows.
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        num_samples++;
    }

    if (num_samples == 0)
    {
        throw std::runtime_error("No spectrum data found.");
    }

    // Allocate the Eigen vectors once.
    Spectrum_t spectrum;

    spectrum.wavelengths.resize(num_samples);
    spectrum.values.resize(num_samples);

    // Reset the file to the beginning.
    file.clear();
    file.seekg(0);

    // Skip the header comments.
    while (file >> std::ws && file.peek() == '#')
    {
        file.ignore(10000, '\n');
    }

    // Second pass: read directly into Eigen.
    for (int i = 0; i < num_samples; ++i)
    {
        if (!(file >> spectrum.wavelengths(i) >> spectrum.values(i)))
        {
            throw std::runtime_error("Error reading spectrum data.");
        }
    }

    return spectrum;
}


Spectrum_t binSpectrum(const Spectrum_t& input, double min_wavelength, double max_wavelength, double wavelength_step_size)
{
    // Calculate the number of bins
    int num_bins;
    num_bins = static_cast<int>((max_wavelength - min_wavelength) / wavelength_step_size) + 1;

    // Create output spectrum
    Spectrum_t binned_spectrum;
    binned_spectrum.wavelengths.resize(num_bins);
    binned_spectrum.values.resize(num_bins);

    // Fill the wavelengths for the binned spectrum
    for (int i = 0; i < num_bins; ++i)
    {
        binned_spectrum.wavelengths(i) = min_wavelength + i * wavelength_step_size;
    }

    // Bin the input spectrum
    for (int i = 0; i < num_bins; ++i)
    {
        double bin_start = binned_spectrum.wavelengths(i);
        double bin_end = bin_start + wavelength_step_size;

        // Find indices of input wavelengths that fall within the current bin
        std::vector<int> indices;
        for (int j = 0; j < input.wavelengths.size(); ++j)
        {
            if (input.wavelengths(j) >= bin_start && input.wavelengths(j) < bin_end)
            {
                indices.push_back(j);
            }
        }

        // If there are values in this bin, calculate the mean value
        if (!indices.empty())
        {
            double sum = 0.0;
            for (int idx : indices)
            {
                sum += input.values(idx);
            }
            binned_spectrum.values(i) = sum / indices.size();
        }
        else
        {
            binned_spectrum.values(i) = 0.0; // or some other default value
        }
    }

    return binned_spectrum;
}