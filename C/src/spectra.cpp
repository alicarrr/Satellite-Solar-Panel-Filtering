#include "spectra.hpp"
#include "constants.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <cstdio>

#include <sstream>



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



Spectrum_t readASTMGColumn(
    const std::string& file_name,
    const int spectrum_column)
{
    if (spectrum_column < 1 || spectrum_column > 3)
    {
        throw std::invalid_argument(
            "ASTM G173 spectrum column must be between 1 and 3.");
    }

    std::ifstream file(file_name);
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + file_name);
    }

    std::string line;

    // Skip the ASTM G173 description and column-heading lines.
    if (!std::getline(file, line) || !std::getline(file, line))
    {
        throw std::runtime_error(
            "ASTM G173 file does not contain the expected header.");
    }

    std::vector<double> wavelengths;
    std::vector<double> values;

    int line_number = 2;

    while (std::getline(file, line))
    {
        ++line_number;

        // Ignore completely empty or whitespace-only lines.
        if (line.find_first_not_of(" \t\r\n") == std::string::npos)
        {
            continue;
        }

        std::stringstream row(line);
        std::string cell;
        std::string wavelength_cell;
        std::string spectrum_cell;

        // Column 0 is wavelength. Columns 1, 2 and 3 contain
        // AM0, AM1.5 global and AM1.5 direct respectively.
        for (int column = 0; column <= spectrum_column; ++column)
        {
            if (!std::getline(row, cell, ','))
            {
                throw std::runtime_error(
                    "Malformed ASTM G173 data at line " +
                    std::to_string(line_number) + ".");
            }

            if (column == 0)
            {
                wavelength_cell = cell;
            }

            if (column == spectrum_column)
            {
                spectrum_cell = cell;
            }
        }

        try
        {
            wavelengths.push_back(std::stod(wavelength_cell));
            values.push_back(std::stod(spectrum_cell));
        }
        catch (const std::exception&)
        {
            throw std::runtime_error(
                "Could not parse numerical data at line " +
                std::to_string(line_number) + ".");
        }
    }

    if (wavelengths.empty())
    {
        throw std::runtime_error(
            "No spectral data were found in file: " + file_name);
    }

    Spectrum_t spectrum;

    const int number_of_samples =
        static_cast<int>(wavelengths.size());

    spectrum.wavelengths.resize(number_of_samples);
    spectrum.values.resize(number_of_samples);

    for (int i = 0; i < number_of_samples; ++i)
    {
        spectrum.wavelengths(i) = wavelengths[i];
        spectrum.values(i) = values[i];
    }

    return spectrum;
}

Spectrum_t readASTMGFileForAM0(const std::string& file_name)
{
    // Column 1: Extraterrestrial spectrum, conventionally called AM0.
    return readASTMGColumn(file_name, 1);
}

Spectrum_t readASTMGFileForAM1_5(const std::string& file_name)
{
    // Column 3: AM1.5 direct plus circumsolar spectrum.
    return readASTMGColumn(file_name, 3);
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


Spectrum_t opticalTransmittanceAtmosphereFromData(const Spectrum_t& AM0_spectrum_binned, const Spectrum_t& AM1_5_spectrum_binned)
{
    Spectrum_t opticalThickness;
    opticalThickness.wavelengths = AM0_spectrum_binned.wavelengths;

    // Calculate fraction
    opticalThickness.values = AM1_5_spectrum_binned.values.cwiseProduct(AM0_spectrum_binned.values.cwiseInverse());

    // Raise fraction to 2/3
    for (int i = 0; i < opticalThickness.values.size(); i++)
    {
        opticalThickness.values(i) = pow(opticalThickness.values(i), 2.0/3.0);
    }

    return opticalThickness;

}

Spectrum_t albedoRadianceTOA(const Spectrum_t& AM0_spectrum_binned, const Spectrum_t& AM1_5_spectrum_binned, Spectrum_t albedo_spectrum_at_loc, Eigen::Vector3d albedo_patch_pos_G, Eigen::Vector3d albedo_patch_norm_vector, Eigen::Vector3d panel_pos_G, Eigen::Vector3d sun_pos_G)
{   
    std::cout << "Hi" << std::endl;
    // AM1.5/AM0 is the transmission through 1.5 air masses. Keeps the ratio between 0 and 1, and avoids divide-by-zero errors.
    Eigen::ArrayXd transmission_ratio_AM1_5 = (AM1_5_spectrum_binned.values.array() / AM0_spectrum_binned.values.array().max(1.0e-12)).max(0.0).min(1.0);

    // Initialise the output spectrum with zeros and the same wavelengths as the albedo spectrum at the location.
    Spectrum_t albedo_toa_spectrum;
    albedo_toa_spectrum.values = Eigen::VectorXd::Zero(albedo_spectrum_at_loc.values.size());
    albedo_toa_spectrum.wavelengths = albedo_spectrum_at_loc.wavelengths;

    Eigen::Vector3d patch_to_sun_G = sun_pos_G - albedo_patch_pos_G;
    Eigen::Vector3d patch_to_panel_G = panel_pos_G - albedo_patch_pos_G;

    double patch_normal_magnitude = albedo_patch_norm_vector.norm();
    double sun_distance = patch_to_sun_G.norm(); // m
    double patch_to_panel_distance = patch_to_panel_G.norm(); // m

    if (patch_normal_magnitude <= 0.0 || sun_distance <= 0.0 || patch_to_panel_distance <= 0.0)
    {
        return albedo_toa_spectrum;
    }

    Eigen::Vector3d patch_norm_dir_vector_G = albedo_patch_norm_vector.normalized();
    Eigen::Vector3d sun_dir_G = patch_to_sun_G.normalized();
    Eigen::Vector3d patch_to_panel_dir_G = patch_to_panel_G.normalized();

    double cos_zenith = patch_norm_dir_vector_G.dot(sun_dir_G);
    double cos_exit = patch_norm_dir_vector_G.dot(patch_to_panel_dir_G);

    // Return zero if the patch is unilluminated or not visible from the satellite.
    if (cos_zenith <= 0.0 || cos_exit <= 0.0)
    {
        // std::cout << "Patch not visible" << std::endl;
        return albedo_toa_spectrum;
    }

    // Prevent acos() producing NaN from numerical round-off.
    cos_zenith = std::min(1.0, cos_zenith);
    cos_exit = std::min(1.0, cos_exit);

    double zenith_angle_rad = std::acos(cos_zenith);
    double zenith_angle_deg = zenith_angle_rad * RAD2DEG;

    double exit_angle_rad = std::acos(cos_exit);
    double exit_angle_deg = exit_angle_rad * RAD2DEG;

    // Kasten-Young air mass, including the near-horizon curvature correction.
    // Simple calc would be m=1/cos(zenith_angle_rad), but this is inaccurate for large zenith angles.
    double m_z = 1.0 /(cos_zenith + 0.50572 * std::pow(96.07995 - zenith_angle_deg, -1.6364));
    double m_z_exit = 1.0 /(cos_exit + 0.50572 * std::pow(96.07995 - exit_angle_deg, -1.6364));

    // Total air mass scale is the sum of the air mass for the incoming solar radiation and the outgoing reflected radiation.
    double total_AM = m_z + m_z_exit;

    Spectrum_t atmosphereTransmittance = opticalTransmittanceAtmosphereFromData(AM0_spectrum_binned, AM1_5_spectrum_binned);

    std::cout << "Total AM: " << total_AM << std::endl;

    Eigen::VectorXd test_t_coefficient = pow(atmosphereTransmittance.values.array(), total_AM) * AM0_spectrum_binned.values.array();

    // print T coefficient
    // std::cout << "Transmission coefficient:" << std::endl << transmission_ratio_AM1_5.pow(total_AM_scale) << std::endl;


    // plot original spectrum and attenuated spectrum using gnuplot
    FILE* gnuplot = popen("gnuplot -persistent", "w");

    std::fprintf(gnuplot,
                 "set title 'Albedo TOA spectrum'\n"
                 "set xlabel \"Wavelength\"\n"
                 "set ylabel \"Spectral value\"\n"
                 "set grid\n"
                 "plot '-' using 1:2 with lines title 'AM0', "
                 "'-' using 1:2 with lines title 'AM1.5', "
                 "'-' using 1:2 with lines title 'Attenuated AM0'\n");
    for (Eigen::Index i = 0; i < AM0_spectrum_binned.wavelengths.size(); ++i)
        std::fprintf(gnuplot, "%.17g %.17g\n", AM0_spectrum_binned.wavelengths(i), AM0_spectrum_binned.values(i));
    std::fprintf(gnuplot, "e\n");

    for (Eigen::Index i = 0; i < AM1_5_spectrum_binned.wavelengths.size(); ++i)
        std::fprintf(gnuplot, "%.17g %.17g\n", AM1_5_spectrum_binned.wavelengths(i), AM1_5_spectrum_binned.values(i));
    std::fprintf(gnuplot, "e\n");

    // now plot test_t_coefficient NOT albedo spectrum
    for (Eigen::Index i = 0; i < AM1_5_spectrum_binned.wavelengths.size(); ++i)
    {
        std::fprintf(gnuplot, "%.17g %.17g\n",
                     AM1_5_spectrum_binned.wavelengths(i),
                     test_t_coefficient(i));
    }
    std::fprintf(gnuplot, "e\n");
    std::fflush(gnuplot);
    pclose(gnuplot);

    // very big bug - fix binning

    // albedo_toa_spectrum.values = (AM0_spectrum_binned.values.array() * transmission_ratio_AM1_5.pow(total_AM_scale) * albedo_spectrum_at_loc.values.array() * cos_zenith * (AU * AU / patch_to_sun_G.squaredNorm())/M_PI).matrix();

    return albedo_toa_spectrum;
}

/** TODO: make function for Earth IR that has a similar attenuation but upward, and with the original spectrum as a blackbody */



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
