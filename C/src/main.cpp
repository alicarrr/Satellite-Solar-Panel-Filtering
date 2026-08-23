// Your main file doesn't have to be called main, but it must contain the main function.
#include "stuff.hpp"
#include "spectra.hpp"
#include "albedo.hpp"
#include "satParams.hpp"

#include <lcm/lcm-cpp.hpp>
#include <vector>
#include "simMessages/sim_data_t.hpp"
#include <iostream>
#include <chrono>
#include <thread>
#include <iostream>
#include <string>


#include "simParams.hpp"



SimParams_t simParams;

const double dt = 1.0; // Sim Time Step
const std::string date_0 = "2025-01-09T01:01:10.050000";




void setSimParams()
{
    try {
        std::string lskPath = "lib/spice/kernels/naif0012.tls";
        std::string spkPath = "lib/spice/kernels/de442s.bsp";
        std::string earthHighPrecPath = "lib/spice/kernels/earth_000101_260601_260305.bpc";

        std::string earthGravPath = "lib/geographiclib/share/GeographicLib/gravity";
        std::string earthMagPath = "lib/geographiclib/share/GeographicLib/magnetic";

        simParams.ephemeris = new SolarSystemEphemeris(lskPath,
                                                       spkPath,
                                                       date_0);

        simParams.ephemeris->loadEarthFieldModels(earthHighPrecPath,
                                                  "egm2008", earthGravPath, 10, 10,//2106 max
                                                  "wmm2025", earthMagPath, 11, 11 // 11 max
                                                ); 
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }

    simParams.time_step = dt;
    simParams.initialUtc = date_0;
    simParams.et_0 = simParams.ephemeris->etAtSeconds(0);
    simParams.tYear_0 = decimalYearFromEtUtc(date_0, 0.0);
}


void clearMem()
{
    delete simParams.ephemeris;
}




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

    // leave at bottom of main - shows no errors
    return 0;
    
}