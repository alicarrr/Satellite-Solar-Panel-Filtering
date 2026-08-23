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

// from sim get satellite pos (GC) and sun pos in ECI (i.e. Earth at origin). Need to find
/**
 * - Panel normal vector direction in global frame for each panel
 * - Panel position in global frame for each panel
 * - Incident sun spectrum on each panel at each time step
 */


Eigen::Vector3d rotateVectorFromBToG(Eigen::Vector3d vector_B, Eigen::Matrix3d R_G_from_B)
{
    return R_G_from_B * vector_B;
}


Eigen::Vector3d panelPosInG(Eigen::Vector3d r_panel_B, Eigen::Matrix3d R_G_from_B, Eigen::Vector3d S_sat_pos_GC_G)
{
    // Transform the panel position from body frame to global frame
    Eigen::Vector3d panel_pos_G = R_G_from_B * r_panel_B + S_sat_pos_GC_G;
    return panel_pos_G;
}


Spectrum_t calculateIncidentSunSpectrumOnPanel(const Spectrum_t& input, Eigen::Vector3d sun_pos_G, Eigen::Vector3d panel_pos_G, Eigen::Vector3d panel_normal_vect_B)
{
    // Calculate the distance from the sun to the panel in AU. Panel position is in meters in global frame, so we need to convert it to AU.
    double distance_AU = (sun_pos_G - panel_pos_G).norm() / AU;

    // Scale the input spectrum based on the inverse square law
    Spectrum_t panel_solar_incident_spectrum;
    panel_solar_incident_spectrum.wavelengths = input.wavelengths;
    panel_solar_incident_spectrum.values = input.values / (distance_AU * distance_AU);

    // Calculate the angle between the sun ray and the satellite panel normal vector
    Eigen::Vector3d vect_panel_to_sun = (sun_pos_G - panel_pos_G).normalized();
    double cos_theta = vect_panel_to_sun.dot(panel_normal_vect_B.normalized()); // should be negative if the sun is behind the panel

    panel_solar_incident_spectrum.values *= std::max(0.0, cos_theta); // Only consider positive angles (sun facing the panel)

    return panel_solar_incident_spectrum;
}


Eigen::Vector3d PanelNormVectInG(Eigen::Vector3d S_sat_pos_GC_G, Eigen::Matrix3d R_G_from_B, Eigen::Vector3d panel_normal_vect_B)
{
    // Transform the panel normal vector from body frame to global frame
    Eigen::Vector3d panel_normal_vect_G = R_G_from_B * panel_normal_vect_B;
    return panel_normal_vect_G.normalized(); // Return the normalized vector
}



SatelliteResults_t computeOrbitalSpectralLoads(CubesatState_t cubesat_state, const Spectrum_t& binned_solar_spectrum)
{
    // at a timestep compute the satellite orbital spectral loads on each panel. This includes solar radiation, Earth IR, and Earth Albedo. 
    SatelliteResults_t satellite_results;

    /** TODO: change later */
    Eigen::Vector3d sun_pos_G = Eigen::Vector3d(0, 0, 0); 
    

    // loop over and calculate each spectral load for each panel
    for (int i = 0; i < cubesat.panels.size(); i++)
    {
        PanelResults_t panel_results;
        panel_results.time = cubesat_state.time; // Using the time parameter
        panel_results.panel_index = i;

        // Current panel
        Panel_t panel = cubesat.panels[i];

        // Calculate the panel normal vector in global frame
        Eigen::Vector3d panel_normal_vect_G = PanelNormVectInG(cubesat_state.position_GC_ECI, cubesat_state.rotation_matrix_ECI_from_B, panel.normal_vector);

        // Calculate the panel position in global frame
        Eigen::Vector3d panel_pos_G = panelPosInG(panel.position_rel_GC_in_B, cubesat_state.rotation_matrix_ECI_from_B, cubesat_state.position_GC_ECI);

        // Incident spectrum from solar radiation on the panel
        panel_results.incident_sun_spectrum = calculateIncidentSunSpectrumOnPanel(binned_solar_spectrum, sun_pos_G, panel_pos_G, panel_normal_vect_G);


        // Incident spectrum from Earth IR on the panel
        /** TODO: */
        

        // Incident spectrum from Earth Albedo on the panel
        /** TODO: */


    }
    return satellite_results;
}
