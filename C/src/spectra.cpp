#include "spectra.hpp"
#include "constants.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <cstdio>

#include <cmath>
#include <limits>

#include <sstream>

SpectraVector_t spectra_wavelengths_nm;
SpectraVector_t sun_spectra; // In watts, binned 
SpectraVector_t atmosphere_transmittance; // Dimensionless, put to power of AM
SpectraVector_t earth_ir_spectra; // In watts, binned

netCDF::NcFile atmosphere_transmissivity_file;


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
    Spectrum_t opticalTransmittance;
    opticalTransmittance.wavelengths = AM0_spectrum_binned.wavelengths;

    // Calculate fraction
    opticalTransmittance.values = AM1_5_spectrum_binned.values;
    for (int i = 0; i < AM0_spectrum_binned.wavelengths.size(); i++)
    {
        if (AM0_spectrum_binned.values(i) == 0)
        {
            opticalTransmittance.values(i) = 0;
        }
        else
        {
            opticalTransmittance.values(i) /= AM0_spectrum_binned.values(i);
        }
    }

    // Raise fraction to 2/3 (i.e. 1/1.5)
    for (int i = 0; i < opticalTransmittance.values.size(); i++)
    {
        opticalTransmittance.values(i) = pow(opticalTransmittance.values(i), 2.0/3.0);
    }

    return opticalTransmittance;
}


Spectrum_t atmosphereTransissionRatio(const Spectrum_t& AM0_spectrum_binned, const Spectrum_t& AM1_5_spectrum_binned)
{   
    // Transmission ratio is the ratio of the AM1.5 spectrum to the AM0 spectrum, which gives the fraction of light that is transmitted through the atmosphere at each wavelength.
    Spectrum_t transmission_ratio_AM1_5;
    transmission_ratio_AM1_5.values = (AM1_5_spectrum_binned.values.array() / AM0_spectrum_binned.values.array());
    transmission_ratio_AM1_5.wavelengths = AM0_spectrum_binned.wavelengths;

    return transmission_ratio_AM1_5;
}


Spectrum_t interpolateSpectrum(Spectrum_t input, Eigen::VectorXd new_wavelengths)
{
    Spectrum_t interpolated_spectrum;
    interpolated_spectrum.wavelengths = new_wavelengths;
    interpolated_spectrum.values.resize(new_wavelengths.size());

    for (int i = 0; i < new_wavelengths.size(); i++)
    {
        double wavelength = new_wavelengths(i);

        // Find the two nearest wavelengths in the input spectrum
        auto it = std::lower_bound(input.wavelengths.data(), input.wavelengths.data() + input.wavelengths.size(), wavelength);
        if (it == input.wavelengths.data())
        {
            // Wavelength is below the range of the input spectrum
            interpolated_spectrum.values(i) = input.values(0);
        }
        else if (it == input.wavelengths.data() + input.wavelengths.size())
        {
            // Wavelength is above the range of the input spectrum
            interpolated_spectrum.values(i) = input.values(input.values.size() - 1);
        }
        else
        {
            // Perform linear interpolation
            int idx_upper = it - input.wavelengths.data();
            int idx_lower = idx_upper - 1;

            double x0 = input.wavelengths(idx_lower);
            double x1 = input.wavelengths(idx_upper);
            double y0 = input.values(idx_lower);
            double y1 = input.values(idx_upper);

            double slope = (y1 - y0) / (x1 - x0);
            interpolated_spectrum.values(i) = y0 + slope * (wavelength - x0);
        }
    }

    return interpolated_spectrum;
}


double cosineOfAngleBetweenVectors(Eigen::Vector3d vector1, Eigen::Vector3d vector2)
{
    double dot_product = vector1.dot(vector2);
    double magnitude_product = vector1.norm() * vector2.norm();

    if (magnitude_product == 0.0)
    {
        throw std::runtime_error("One or both vectors have zero magnitude.");
    }

    double cos_angle = dot_product / magnitude_product;
    return cos_angle;
}


double AMKastenYoung(double zenith_angle_rad)
{
    double zenith_angle_deg = zenith_angle_rad * RAD2DEG;

    // Kasten-Young air mass, including the near-horizon curvature correction.
    double m_z = 1.0 /(std::cos(zenith_angle_rad) + 0.50572 * std::pow(96.07995 - zenith_angle_deg, -1.6364));

    return m_z;
}

Spectrum_t scaledTransmissivityForAM(double air_mass, const Spectrum_t& atmosphereTransmittance)
{
    Spectrum_t transmissivity_spectrum =
    {   
        .wavelengths = atmosphereTransmittance.wavelengths,
        .values = pow(atmosphereTransmittance.values.array(), air_mass)        
    };
    
    return transmissivity_spectrum;
}


Spectrum_t albedoRadianceTOA(const Spectrum_t& transmission_ratio_AM1_5, const Spectrum_t& AM0, Spectrum_t albedo_spectrum_at_loc, Eigen::Vector3d albedo_patch_pos_G, Eigen::Vector3d albedo_patch_norm_vector, Eigen::Vector3d panel_pos_G, Eigen::Vector3d sun_pos_G)
{   
    // Interpolate albedo spectrum to match the wavelengths of the transmission ratio spectrum. 
    albedo_spectrum_at_loc = interpolateSpectrum(albedo_spectrum_at_loc, transmission_ratio_AM1_5.wavelengths);

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
    // Use cosineOfAngleBetweenVectors to calculate the needed zenith and exit angles
    double cos_zenith = cosineOfAngleBetweenVectors(albedo_patch_norm_vector, patch_to_sun_G);
    double cos_exit = cosineOfAngleBetweenVectors(albedo_patch_norm_vector, patch_to_panel_G);

    // Eigen::Vector3d patch_norm_dir_vector_G = albedo_patch_norm_vector.normalized();
    // Eigen::Vector3d sun_dir_G = patch_to_sun_G.normalized();
    // Eigen::Vector3d patch_to_panel_dir_G = patch_to_panel_G.normalized();

    // double cos_zenith = patch_norm_dir_vector_G.dot(sun_dir_G);
    // double cos_exit = patch_norm_dir_vector_G.dot(patch_to_panel_dir_G);

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
    double m_z = AMKastenYoung(zenith_angle_rad);
    double m_z_exit = AMKastenYoung(exit_angle_rad);

    // Total air mass scale is the sum of the air mass for the incoming solar radiation and the outgoing reflected radiation.
    double total_AM = m_z + m_z_exit;

    std::cout << "Total AM: " << total_AM << std::endl;

    // Scaled AM transmittance
    Spectrum_t scaled_transmittance = scaledTransmissivityForAM(total_AM, transmission_ratio_AM1_5);
    
    Spectrum_t transmissionSpectrumAMTotal =
    {
        .wavelengths = scaled_transmittance.wavelengths,
        .values = scaled_transmittance.values.array() * AM0.values.array()
    };

    // radiance in W/m^2/sr/nm at the top of the atmosphere (TOA)

    albedo_toa_spectrum.values = transmissionSpectrumAMTotal.values.array() * albedo_spectrum_at_loc.values.array() / M_PI;

    // NOTE: to get irradiance at panel multiply by cos(angle between panel and patch), patch area, and divide by distance squared. 

    return albedo_toa_spectrum;
}

/** TODO: make function for Earth IR that has a similar attenuation but upward, and with the original spectrum as a blackbody */

Spectrum_t blackBodySpectrum(double temperature_K, const Spectrum_t& wavelength_spectrum)
{
    // Planck's law for blackbody radiation
    Spectrum_t blackbody_spectrum;
    blackbody_spectrum.wavelengths = wavelength_spectrum.wavelengths;
    blackbody_spectrum.values.resize(wavelength_spectrum.wavelengths.size());

    for (int i = 0; i < wavelength_spectrum.wavelengths.size(); i++)
    {
        double wavelength_m = wavelength_spectrum.wavelengths(i) * nm; // Convert nm to m
        double exponent = (h * c) / (wavelength_m * k * temperature_K);
        blackbody_spectrum.values(i) = (2.0 * h * c * c) / (std::pow(wavelength_m, 5) * (std::exp(exponent) - 1.0)); 
    }

    return blackbody_spectrum;
}



// function to compute Earth IR spectrum at the top of the atmosphere (TOA) based on a blackbody spectrum at 288K, scaled by the cosine of the angle between the sun and the panel normal vector, and attenuated by the atmosphere. This is a simplified model and does not account for all atmospheric effects.
Spectrum_t EarthIRRadianceTOA(const Spectrum_t& transmission_ratio_AM1_5, const Spectrum_t& AM0, Eigen::Vector3d patch_pos_G, Eigen::Vector3d panel_pos_G)
{
    // for earth IR the distance calc for equivalent air mass is just the path travelled from the earth emitting the IR to the satellite panel, which is just the distance from the earth to the satellite panel. The zenith angle is the angle between the panel normal vector and the vector from the panel to the earth center (which is at the origin in ECI). 
    Spectrum_t earth_ir_toa_spectrum;
    earth_ir_toa_spectrum.wavelengths = transmission_ratio_AM1_5.wavelengths;
    earth_ir_toa_spectrum.values = Eigen::VectorXd::Zero(transmission_ratio_AM1_5.values.size());

    Eigen::Vector3d patch_to_panel_G = panel_pos_G - patch_pos_G;

    double patch_normal_magnitude = patch_pos_G.norm();
    double patch_to_panel_distance = patch_to_panel_G.norm(); // m

    // Use cosineOfAngleBetweenVectors to calculate the needed zenith and exit angles
    double cos_out = cosineOfAngleBetweenVectors(patch_pos_G, patch_to_panel_G);

    // check if the patch is visible from the satellite panel
    if (cos_out <= 0.0)
    {
        return earth_ir_toa_spectrum;
    }

    cos_out = std::min(1.0, cos_out);
    
    double exit_angle_rad = std::acos(cos_out);
    double exit_angle_deg = exit_angle_rad * RAD2DEG;

    // Kasten-Young air mass, including the near-horizon curvature correction.
    double m_z_exit = AMKastenYoung(exit_angle_rad);

    // Scaled AM transmittance
    Spectrum_t scaled_transmittance = scaledTransmissivityForAM(m_z_exit, transmission_ratio_AM1_5);

    Spectrum_t Earth_blackbody_spectrum = blackBodySpectrum(288.0, transmission_ratio_AM1_5);

    earth_ir_toa_spectrum.values = (scaled_transmittance.values.array() * Earth_blackbody_spectrum.values.array() * cos_out * (AU * AU / patch_to_panel_G.squaredNorm())/M_PI).matrix();

    return earth_ir_toa_spectrum;

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



////////////////////////////////////////////////////////////////////////
/////////////               SPECTRAL RE-BINNING          ///////////////
////////////////////////////////////////////////////////////////////////

static double evaluateExtensionFunction(double (*extension_function)(double), double wavelength_nm)
{
    const double value = extension_function(wavelength_nm);

    if (!std::isfinite(value))
    {
        throw std::domain_error(
            "The domain extension function returned a non-finite value.");
    }

    return value;
}


static double simpsonEstimate(double lower_wavelength_nm, double upper_wavelength_nm, double lower_value, double midpoint_value, double upper_value)
{
    return (upper_wavelength_nm - lower_wavelength_nm)
        * (lower_value + 4.0 * midpoint_value + upper_value)
        / 6.0;
}


static double adaptiveSimpsonRecursive(
    double (*extension_function)(double),
    double lower_wavelength_nm,
    double upper_wavelength_nm,
    double lower_value,
    double midpoint_value,
    double upper_value,
    double previous_estimate,
    double tolerance,
    unsigned int remaining_depth)
{
    const double midpoint_nm =
        lower_wavelength_nm
        + 0.5 * (upper_wavelength_nm - lower_wavelength_nm);

    const double left_midpoint_nm =
        lower_wavelength_nm
        + 0.5 * (midpoint_nm - lower_wavelength_nm);

    const double right_midpoint_nm =
        midpoint_nm
        + 0.5 * (upper_wavelength_nm - midpoint_nm);

    const double left_midpoint_value =
        evaluateExtensionFunction(
            extension_function,
            left_midpoint_nm);

    const double right_midpoint_value =
        evaluateExtensionFunction(
            extension_function,
            right_midpoint_nm);

    const double left_estimate = simpsonEstimate(
        lower_wavelength_nm,
        midpoint_nm,
        lower_value,
        left_midpoint_value,
        midpoint_value);

    const double right_estimate = simpsonEstimate(
        midpoint_nm,
        upper_wavelength_nm,
        midpoint_value,
        right_midpoint_value,
        upper_value);

    const double refined_estimate =
        left_estimate + right_estimate;

    const double estimated_error =
        refined_estimate - previous_estimate;

    if (remaining_depth == 0
        || std::fabs(estimated_error) <= 15.0 * tolerance)
    {
        // Richardson correction for Simpson integration.
        return refined_estimate + estimated_error / 15.0;
    }

    return adaptiveSimpsonRecursive(
               extension_function,
               lower_wavelength_nm,
               midpoint_nm,
               lower_value,
               left_midpoint_value,
               midpoint_value,
               left_estimate,
               0.5 * tolerance,
               remaining_depth - 1)
        + adaptiveSimpsonRecursive(
               extension_function,
               midpoint_nm,
               upper_wavelength_nm,
               midpoint_value,
               right_midpoint_value,
               upper_value,
               right_estimate,
               0.5 * tolerance,
               remaining_depth - 1);
}


static double integrateExtensionFunction(
    double (*extension_function)(double),
    double lower_wavelength_nm,
    double upper_wavelength_nm)
{
    if (extension_function == nullptr
        || upper_wavelength_nm <= lower_wavelength_nm)
    {
        return 0.0;
    }

    const double midpoint_nm =
        lower_wavelength_nm
        + 0.5 * (upper_wavelength_nm - lower_wavelength_nm);

    const double lower_value =
        evaluateExtensionFunction(
            extension_function,
            lower_wavelength_nm);

    const double midpoint_value =
        evaluateExtensionFunction(
            extension_function,
            midpoint_nm);

    const double upper_value =
        evaluateExtensionFunction(
            extension_function,
            upper_wavelength_nm);

    const double initial_estimate = simpsonEstimate(
        lower_wavelength_nm,
        upper_wavelength_nm,
        lower_value,
        midpoint_value,
        upper_value);

    const double absolute_tolerance = 1.0e-12;
    const double relative_tolerance = 1.0e-10;

    const double tolerance =
        absolute_tolerance
        + relative_tolerance * std::fabs(initial_estimate);

    const unsigned int maximum_recursion_depth = 20;

    return adaptiveSimpsonRecursive(
        extension_function,
        lower_wavelength_nm,
        upper_wavelength_nm,
        lower_value,
        midpoint_value,
        upper_value,
        initial_estimate,
        tolerance,
        maximum_recursion_depth);
}

static double evaluateCumulativeSourceFlux(const std::vector<double>& wavelengths_nm, const std::vector<double>& spectral_irradiances, const std::vector<double>& cumulative_flux, double wavelength_nm)
{
    if (wavelength_nm <= wavelengths_nm.front())
    {
        return 0.0;
    }

    if (wavelength_nm >= wavelengths_nm.back())
    {
        return cumulative_flux.back();
    }

    const std::vector<double>::const_iterator upper_iterator = std::upper_bound(wavelengths_nm.begin(), wavelengths_nm.end(), wavelength_nm);

    const unsigned int upper_index = static_cast<unsigned int>(upper_iterator - wavelengths_nm.begin());

    const unsigned int lower_index = upper_index - 1;

    const double interval_width_nm = wavelengths_nm[upper_index] - wavelengths_nm[lower_index];

    const double slope = (spectral_irradiances[upper_index] - spectral_irradiances[lower_index]) / interval_width_nm;

    const double offset_nm = wavelength_nm - wavelengths_nm[lower_index];

    return cumulative_flux[lower_index] + spectral_irradiances[lower_index] * offset_nm + 0.5 * slope * offset_nm * offset_nm;
}



std::vector<double> rebinSpectrumConservingFlux(
    const std::vector<double>& wavelengths_nm,
    const std::vector<double>& spectral_irradiances,
    double min_wavelength_nm,
    double max_wavelength_nm,
    double step_wavelength_nm,
    double (*domain_extension_function)(double))
{
    if (wavelengths_nm.size() != spectral_irradiances.size())
    {
        throw std::invalid_argument(
            "The wavelength and spectral irradiance vectors "
            "must have the same length.");
    }

    if (wavelengths_nm.size() < 2)
    {
        throw std::invalid_argument(
            "At least two source spectrum points are required.");
    }

    if (!std::isfinite(min_wavelength_nm)
        || !std::isfinite(max_wavelength_nm)
        || !std::isfinite(step_wavelength_nm))
    {
        throw std::invalid_argument(
            "The output wavelength parameters must be finite.");
    }

    if (max_wavelength_nm <= min_wavelength_nm)
    {
        throw std::invalid_argument(
            "The maximum wavelength must exceed the minimum wavelength.");
    }

    if (step_wavelength_nm <= 0.0)
    {
        throw std::invalid_argument(
            "The wavelength step must be positive.");
    }

    if (wavelengths_nm.size()
        > std::numeric_limits<unsigned int>::max())
    {
        throw std::length_error(
            "The input spectrum contains too many points.");
    }

    const unsigned int number_of_input_points =
        static_cast<unsigned int>(wavelengths_nm.size());

    for (unsigned int i = 0; i < number_of_input_points; ++i)
    {
        if (!std::isfinite(wavelengths_nm[i])
            || !std::isfinite(spectral_irradiances[i]))
        {
            throw std::invalid_argument(
                "All source wavelengths and irradiances must be finite.");
        }

        if (i > 0 && wavelengths_nm[i] <= wavelengths_nm[i - 1])
        {
            throw std::invalid_argument(
                "The source wavelengths must be strictly increasing.");
        }
    }

    const double unrounded_number_of_bins =
        (max_wavelength_nm - min_wavelength_nm)
        / step_wavelength_nm;

    const double rounded_number_of_bins =
        std::round(unrounded_number_of_bins);

    const double bin_count_tolerance =
        1.0e-10
        * std::max(1.0, std::fabs(unrounded_number_of_bins));

    if (std::fabs(
            unrounded_number_of_bins - rounded_number_of_bins)
        > bin_count_tolerance)
    {
        throw std::invalid_argument(
            "The output wavelength range must be divisible "
            "by the wavelength step.");
    }

    if (rounded_number_of_bins < 1.0
        || rounded_number_of_bins
            > static_cast<double>(
                std::numeric_limits<unsigned int>::max()))
    {
        throw std::length_error(
            "The requested number of output bins is invalid.");
    }

    const unsigned int number_of_output_bins =
        static_cast<unsigned int>(rounded_number_of_bins);

    /*
     * cumulative_flux[i] is the integrated piecewise-linear source
     * spectrum between wavelengths_nm[0] and wavelengths_nm[i].
     */
    std::vector<double> cumulative_flux(
        number_of_input_points,
        0.0);

    for (unsigned int i = 0;
         i + 1 < number_of_input_points;
         ++i)
    {
        const double interval_width_nm =
            wavelengths_nm[i + 1] - wavelengths_nm[i];

        const double interval_flux =
            0.5
            * (spectral_irradiances[i]
               + spectral_irradiances[i + 1])
            * interval_width_nm;

        cumulative_flux[i + 1] =
            cumulative_flux[i] + interval_flux;

        if (!std::isfinite(cumulative_flux[i + 1]))
        {
            throw std::overflow_error(
                "The integrated source flux overflowed.");
        }
    }

    const double source_min_wavelength_nm =
        wavelengths_nm.front();

    const double source_max_wavelength_nm =
        wavelengths_nm.back();

    std::vector<double> output_spectrum;
    std::vector<double> output_wavelengths_nm;

    output_wavelengths_nm.reserve(number_of_output_bins);
    output_spectrum.reserve(number_of_output_bins);

    for (unsigned int j = 0;
         j < number_of_output_bins;
         ++j)
    {
        const double bin_lower_nm =
            min_wavelength_nm
            + static_cast<double>(j) * step_wavelength_nm;

        /*
         * Set the final edge directly from max_wavelength_nm to avoid
         * accumulation of floating-point error.
         */
        const double bin_upper_nm =
            (j + 1 == number_of_output_bins)
            ? max_wavelength_nm
            : min_wavelength_nm
                + static_cast<double>(j + 1)
                    * step_wavelength_nm;

        const double bin_width_nm =
            bin_upper_nm - bin_lower_nm;

        double bin_flux = 0.0;

        /*
         * Integrate the portion covered by the source dataset.
         */
        const double known_lower_nm =
            std::max(bin_lower_nm, source_min_wavelength_nm);

        const double known_upper_nm =
            std::min(bin_upper_nm, source_max_wavelength_nm);

        if (known_upper_nm > known_lower_nm)
        {
            bin_flux += evaluateCumulativeSourceFlux(
                            wavelengths_nm,
                            spectral_irradiances,
                            cumulative_flux,
                            known_upper_nm)
                - evaluateCumulativeSourceFlux(
                            wavelengths_nm,
                            spectral_irradiances,
                            cumulative_flux,
                            known_lower_nm);
        }

        /*
         * Integrate any part of the bin below the source domain.
         */
        const double lower_extension_upper_nm =
            std::min(bin_upper_nm, source_min_wavelength_nm);

        if (lower_extension_upper_nm > bin_lower_nm)
        {
            bin_flux += integrateExtensionFunction(
                domain_extension_function,
                bin_lower_nm,
                lower_extension_upper_nm);
        }

        /*
         * Integrate any part of the bin above the source domain.
         */
        const double upper_extension_lower_nm =
            std::max(bin_lower_nm, source_max_wavelength_nm);

        if (bin_upper_nm > upper_extension_lower_nm)
        {
            bin_flux += integrateExtensionFunction(
                domain_extension_function,
                upper_extension_lower_nm,
                bin_upper_nm);
        }

        const double bin_centre_nm =
            bin_lower_nm + 0.5 * bin_width_nm;

        const double bin_average_irradiance =
            bin_flux / bin_width_nm;

        output_wavelengths_nm.push_back(bin_centre_nm);
        output_spectrum.push_back(
            bin_average_irradiance);
    }

    return output_spectrum;
}



std::vector<double> calculateBinCentreWavelengths(double min_wavelength_nm, double max_wavelength_nm, double step_wavelength_nm)
{
    if (!std::isfinite(min_wavelength_nm)
        || !std::isfinite(max_wavelength_nm)
        || !std::isfinite(step_wavelength_nm))
    {
        throw std::invalid_argument(
            "Wavelength parameters must be finite.");
    }

    if (max_wavelength_nm <= min_wavelength_nm)
    {
        throw std::invalid_argument(
            "Maximum wavelength must exceed minimum wavelength.");
    }

    if (step_wavelength_nm <= 0.0)
    {
        throw std::invalid_argument(
            "Wavelength step must be positive.");
    }

    const double unrounded_number_of_bins =
        (max_wavelength_nm - min_wavelength_nm)
        / step_wavelength_nm;

    const double rounded_number_of_bins =
        std::round(unrounded_number_of_bins);

    const double tolerance =
        1.0e-10
        * std::max(1.0, std::fabs(unrounded_number_of_bins));

    if (std::fabs(
            unrounded_number_of_bins - rounded_number_of_bins)
        > tolerance)
    {
        throw std::invalid_argument(
            "Wavelength range must be divisible by the wavelength step.");
    }

    if (rounded_number_of_bins
        > static_cast<double>(
            std::numeric_limits<unsigned int>::max()))
    {
        throw std::length_error(
            "Requested number of wavelength bins is too large.");
    }

    const unsigned int number_of_bins =
        static_cast<unsigned int>(rounded_number_of_bins);

    std::vector<double> bin_centres(number_of_bins);

    for (unsigned int i = 0; i < number_of_bins; ++i)
    {
        bin_centres[i] =
            min_wavelength_nm
            + (static_cast<double>(i) + 0.5)
                * step_wavelength_nm;
    }

    return bin_centres;
}




static double interpolateSpectralValue(
    const std::vector<double>& wavelengths_nm,
    const std::vector<double>& spectral_values,
    double wavelength_nm)
{
    if (wavelength_nm <= wavelengths_nm.front())
    {
        return spectral_values.front();
    }

    if (wavelength_nm >= wavelengths_nm.back())
    {
        return spectral_values.back();
    }

    const std::vector<double>::const_iterator upper_iterator =
        std::lower_bound(
            wavelengths_nm.begin(),
            wavelengths_nm.end(),
            wavelength_nm);

    const unsigned int upper_index =
        static_cast<unsigned int>(
            upper_iterator - wavelengths_nm.begin());

    /*
     * Return the original value directly if the requested wavelength
     * exactly matches a source wavelength.
     */
    if (wavelengths_nm[upper_index] == wavelength_nm)
    {
        return spectral_values[upper_index];
    }

    const unsigned int lower_index = upper_index - 1;

    const double interpolation_fraction =
        (wavelength_nm - wavelengths_nm[lower_index])
        / (wavelengths_nm[upper_index]
           - wavelengths_nm[lower_index]);

    return spectral_values[lower_index]
        + interpolation_fraction
            * (spectral_values[upper_index]
               - spectral_values[lower_index]);
}


std::vector<double> resampleAndExtendSpectrum(
    const std::vector<double>& wavelengths_nm,
    const std::vector<double>& spectral_values,
    double min_wavelength_nm,
    double max_wavelength_nm,
    double step_wavelength_nm,
    double (*domain_extension_function)(double))
{
    if (wavelengths_nm.size() != spectral_values.size())
    {
        throw std::invalid_argument(
            "Wavelength and spectral value vectors must have equal lengths.");
    }

    if (wavelengths_nm.size() < 2)
    {
        throw std::invalid_argument(
            "At least two source spectrum points are required.");
    }

    if (!std::isfinite(min_wavelength_nm)
        || !std::isfinite(max_wavelength_nm)
        || !std::isfinite(step_wavelength_nm))
    {
        throw std::invalid_argument(
            "Output wavelength parameters must be finite.");
    }

    if (max_wavelength_nm <= min_wavelength_nm)
    {
        throw std::invalid_argument(
            "Maximum wavelength must exceed minimum wavelength.");
    }

    if (step_wavelength_nm <= 0.0)
    {
        throw std::invalid_argument(
            "Wavelength step must be positive.");
    }

    if (wavelengths_nm.size()
        > std::numeric_limits<unsigned int>::max())
    {
        throw std::length_error(
            "The source spectrum contains too many points.");
    }

    const unsigned int number_of_source_points =
        static_cast<unsigned int>(wavelengths_nm.size());

    for (unsigned int i = 0; i < number_of_source_points; ++i)
    {
        if (!std::isfinite(wavelengths_nm[i])
            || !std::isfinite(spectral_values[i]))
        {
            throw std::invalid_argument(
                "All source wavelengths and values must be finite.");
        }

        if (i > 0 && wavelengths_nm[i] <= wavelengths_nm[i - 1])
        {
            throw std::invalid_argument(
                "Source wavelengths must be strictly increasing.");
        }
    }

    const double unrounded_number_of_bins =
        (max_wavelength_nm - min_wavelength_nm)
        / step_wavelength_nm;

    const double rounded_number_of_bins =
        std::round(unrounded_number_of_bins);

    const double tolerance =
        1.0e-10
        * std::max(1.0, std::fabs(unrounded_number_of_bins));

    if (std::fabs(
            unrounded_number_of_bins - rounded_number_of_bins)
        > tolerance)
    {
        throw std::invalid_argument(
            "The wavelength range must be divisible by the wavelength step.");
    }

    if (rounded_number_of_bins
        > static_cast<double>(
            std::numeric_limits<unsigned int>::max()))
    {
        throw std::length_error(
            "The requested number of output bins is too large.");
    }

    const unsigned int number_of_output_bins =
        static_cast<unsigned int>(rounded_number_of_bins);

    std::vector<double> output_values;
    output_values.reserve(number_of_output_bins);

    for (unsigned int i = 0; i < number_of_output_bins; ++i)
    {
        const double bin_lower_nm =
            min_wavelength_nm
            + static_cast<double>(i) * step_wavelength_nm;

        const double bin_upper_nm =
            (i + 1 == number_of_output_bins)
            ? max_wavelength_nm
            : min_wavelength_nm
                + static_cast<double>(i + 1)
                    * step_wavelength_nm;

        const double bin_centre_nm =
            bin_lower_nm + 0.5 * (bin_upper_nm - bin_lower_nm);

        /*
         * Use the extension function whenever the bin centre lies outside
         * the source dataset's domain.
         */
        if (bin_centre_nm < wavelengths_nm.front()
            || bin_centre_nm > wavelengths_nm.back())
        {
            double extended_value = 0.0;

            if (domain_extension_function != nullptr)
            {
                extended_value =
                    domain_extension_function(bin_centre_nm);

                if (!std::isfinite(extended_value))
                {
                    throw std::domain_error(
                        "The domain extension function returned "
                        "a non-finite value.");
                }
            }

            output_values.push_back(extended_value);
            continue;
        }

        /*
         * Use half-open bins [lower, upper), except for the final bin,
         * which includes its upper boundary.
         */
        const std::vector<double>::const_iterator first_iterator =
            std::lower_bound(
                wavelengths_nm.begin(),
                wavelengths_nm.end(),
                bin_lower_nm);

        std::vector<double>::const_iterator end_iterator;

        if (i + 1 == number_of_output_bins)
        {
            end_iterator = std::upper_bound(
                wavelengths_nm.begin(),
                wavelengths_nm.end(),
                bin_upper_nm);
        }
        else
        {
            end_iterator = std::lower_bound(
                wavelengths_nm.begin(),
                wavelengths_nm.end(),
                bin_upper_nm);
        }

        const unsigned int first_index =
            static_cast<unsigned int>(
                first_iterator - wavelengths_nm.begin());

        const unsigned int end_index =
            static_cast<unsigned int>(
                end_iterator - wavelengths_nm.begin());

        const unsigned int number_of_points_in_bin =
            end_index - first_index;

        if (number_of_points_in_bin > 0)
        {
            double value_sum = 0.0;

            for (unsigned int j = first_index;
                 j < end_index;
                 ++j)
            {
                value_sum += spectral_values[j];
            }

            output_values.push_back(
                value_sum
                / static_cast<double>(number_of_points_in_bin));
        }
        else
        {
            output_values.push_back(
                interpolateSpectralValue(
                    wavelengths_nm,
                    spectral_values,
                    bin_centre_nm));
        }
    }

    return output_values;
}



double sun_spectra_extend(double wavelength_nm)
{
    return 0;
}

double earth_ir_spectra_extend(double wavelength_nm)
{
    return 0;
}

double am_0_spectra_extend(double wavelength_nm)
{
    return 0;
}

double am_1_5_spectra_extend(double wavelength_nm)
{
    return 0;
}


void getTransmittanceFromFile(const std::string& fileName)
{
    atmosphere_transmissivity_file.open(fileName, netCDF::NcFile::read);

    netCDF::NcVar wavelength_var = atmosphere_transmissivity_file.getVar("wavelength_nm");
    netCDF::NcVar transmissivity_var = atmosphere_transmissivity_file.getVar("total_downward_transmissivity");

    int dimension = wavelength_var.getDim(0).getSize();

    std::vector<double> wavelengths_nm(dimension);
    std::vector<double> transmissivity_vals(dimension);

    std::vector<std::size_t> start(1);
    start[0] = 0;

    std::vector<std::size_t> count(1);
    count[0] = dimension;

    wavelength_var.getVar(start, count, wavelengths_nm.data());

    std::cout << wavelengths_nm[101] << std::endl;

}