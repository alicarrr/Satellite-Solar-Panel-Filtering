// Your main file doesn't have to be called main, but it must contain the main function.
#include "stuff.hpp"
#include "spectra.hpp"
#include "albedo.hpp"

#include <lcm/lcm-cpp.hpp>
#include <vector>
#include "simMessages/sim_data_t.hpp"
#include <iostream>
#include <chrono>
#include <thread>
#include <iostream>


/*
This function is your executable function, it is the function that gets called when you run your code.
*/
int main()
{
    // std::this_thread::sleep_for(std::chrono::milliseconds(3000));

    // lcm::LCM lcm;

    // if(!lcm.good())
    // {
    //     return 1;
    // }
    
    // // std::cout<<getSunSpectrum("hi")<<std::endl;

    // for (int i=0; i<3; i++)
    // {
    //     simMessages::sim_data_t test;
    //     test.t = 0;
    //     test.pos[0] = 0;
    //     test.pos[1] = 1;
    //     test.pos[2] = -2;

    //     lcm.publish("TEST", &test);
    //     std::cout<<"Published"<<std::endl;
    // }

    // // Call albedo data reading function to test it
    // AlbedoData albedo_data = read_albedo_data();

    // double latitude_deg = 23.0f; // Example latitude
    // double longitude_deg = 30.0f; // Example longitude

    // Eigen::VectorXd albedo_spectrum = get_albedo_spectrum(albedo_data, latitude_deg, longitude_deg);

    // std::cout << "Albedo Spectrum:" << std::endl;
    // for (int i = 0; i < albedo_spectrum.size(); ++i) {
    //     std::cout << albedo_spectrum(i) << " ";
    // }
    // std::cout << std::endl;


    // FILE* gnuplot = popen("gnuplot -persistent", "w");

    // if (gnuplot == nullptr)
    // {
    //     std::cerr << "Could not open gnuplot." << std::endl;
    //     return 1;
    // }

    // fprintf(gnuplot, "set title 'HAMSTER spectral albedo'\n");
    // fprintf(gnuplot, "set xlabel 'Wavelength [nm]'\n");
    // fprintf(gnuplot, "set ylabel 'Black-sky albedo'\n");
    // fprintf(gnuplot, "set grid\n");

    // fprintf(
    //     gnuplot,
    //     "plot '-' with lines title 'Lat %.2f, Lon %.2f'\n",
    //     latitude_deg,
    //     longitude_deg
    // );

    // for (unsigned int i = 0; i < albedo_data.n_wavelength; i++)
    // {
    //     fprintf(
    //         gnuplot,
    //         "%f %f\n",
    //         albedo_data.wavelength(i),
    //         albedo_spectrum(i)
    //     );
    // }

    // fprintf(gnuplot, "e\n");

    // pclose(gnuplot);

    //Test spectrum reading from file and binning
    Spectrum_t spectrum = readSpectrumDat("dat files/kurudz_full.dat");

    // Test binning the spectrum to a new wavelength range and step size
    double min_wavelength = 300.0;
    double max_wavelength = 2500.0;
    double wavelength_step_size = 10.0;

    Spectrum_t binned_spectrum = binSpectrum(spectrum, min_wavelength, max_wavelength, wavelength_step_size);

    
    
    // check that the binned spectrum has the same integrated value as the original spectrum
    double unbinned_power = 0.0;

    for (Eigen::Index i = 1; i < spectrum.values.size() - 1; ++i)
    {
        if (spectrum.wavelengths(i) < min_wavelength ||
            spectrum.wavelengths(i) > max_wavelength)
        {
            continue;
        }

        const double wavelength_width =
            0.5 * (spectrum.wavelengths(i + 1) - spectrum.wavelengths(i - 1));

        unbinned_power += spectrum.values(i) * wavelength_width;
    }

    const double binned_power =
        binned_spectrum.values.sum() * wavelength_step_size;

    std::cout << "Unbinned power: " << unbinned_power << std::endl;
    std::cout << "Binned power:   " << binned_power << std::endl;

    // Make a plot of the unbinned spectrum, and the binned spectrum overlaid to compare
    FILE* gnuplot = popen("gnuplot -persistent", "w");

    if (gnuplot == nullptr)
    {
        std::cerr << "Could not open gnuplot." << std::endl;
        return 1;
    }
    // Set x axis to be from min_wavelength to max_wavelength
    fprintf(gnuplot, "set xrange [%f:%f]\n", min_wavelength, max_wavelength);

    fprintf(gnuplot, "set title 'Spectrum Binning Comparison'\n");
    fprintf(gnuplot, "set xlabel 'Wavelength [nm]'\n");
    fprintf(gnuplot, "set ylabel 'Value'\n");
    fprintf(gnuplot, "set grid\n");

    fprintf(gnuplot, "plot '-' with dots title 'Unbinned', '-' with lines title 'Binned'\n");

    // Plot the unbinned spectrum
    for (unsigned int i = 0; i < spectrum.wavelengths.size(); i++)
    {
        fprintf(gnuplot, "%f %f\n", spectrum.wavelengths(i), spectrum.values(i));
    }
    fprintf(gnuplot, "e\n");

    // Plot the binned spectrum
    for (unsigned int i = 0; i < binned_spectrum.wavelengths.size(); i++)
    {
        fprintf(gnuplot, "%f %f\n", binned_spectrum.wavelengths(i), binned_spectrum.values(i));
    }
    fprintf(gnuplot, "e\n");

    pclose(gnuplot);

    // close_albedo_data(albedo_data);

    // leave at bottom of main - shows no errors
    return 0;
    
}