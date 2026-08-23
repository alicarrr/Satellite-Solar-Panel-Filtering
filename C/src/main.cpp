// Your main file doesn't have to be called main, but it must contain the main function.
#include "stuff.hpp"
#include "spectra.hpp"
#include "albedo.hpp"
#include "satParams.hpp"
#include "vectorMaths.hpp"

#include <lcm/lcm-cpp.hpp>
#include <vector>
#include "simMessages/sim_data_t.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <string>
#include <unistd.h>
#include <cstdlib>


#include "simParams.hpp"

#include "satDynamics.hpp"
#include "LSODA.hpp"



SimParams_t simParams;

SimState_t simState;
CubesatState_t cubesatState;

const double dt_data = 1.0; 
const double dt = 1e-1; // Sim Time Step
const std::string date_0 = "2025-01-09T01:01:10.050000";

const double t_final_sim = 1 * minute;

double t_last_data = -2 * dt_data;


void setSimParams()
{
    try {
        std::string lskPath = "C/lib/spice/kernels/naif0012.tls";
        std::string spkPath = "C/lib/spice/kernels/de442s.bsp";
        std::string earthHighPrecPath = "C/lib/spice/kernels/earth_000101_260601_260305.bpc";

        std::string earthGravPath = "C/lib/geographiclib/share/GeographicLib/gravity";
        std::string earthMagPath = "C/lib/geographiclib/share/GeographicLib/magnetic";

        simParams.ephemeris = new SolarSystemEphemeris(lskPath,
                                                       spkPath,
                                                       date_0);

        simParams.ephemeris->loadEarthFieldModels(earthHighPrecPath,
                                                    "egm2008", earthGravPath, 2, 2,//2106 max
                                                    "wmm2025", earthMagPath, 11, 11 // 11 max
                                                 ); 
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }

    simParams.time_step = dt;
    simParams.t_final = t_final_sim;
    simParams.initialUtc = date_0;
    simParams.et_0 = simParams.ephemeris->etAtSeconds(0);
    simParams.tYear_0 = decimalYearFromEtUtc(date_0, 0.0);
}


void clearMem()
{
    delete simParams.ephemeris;
}


void setInitialOrbitState()
{
    double period = 95.2 * minute; // seconds

    double mean_motion = (2.0 * M_PI) / period;

    double e = 0.002; // Eccentricity
    double M = 0.0 * M_PI; // Mean anomaly at epoch

    double argp = 30 * DEG2RAD;
    double raan = -90 * DEG2RAD;
    double inc = 97.5 * DEG2RAD;


    tuple<Eigen::Vector3d, Eigen::Vector3d> r_v_eci = angular_parameters_to_eci(M, e, mean_motion, inc, raan, argp, mu_earth);

    Eigen::Vector3d s_com_in_G = std::get<0>(r_v_eci); // m
    Eigen::Vector3d v_com_in_G = std::get<1>(r_v_eci); // m/s

    std::cout << "Initial position in G frame (km): " << s_com_in_G.transpose()/km << std::endl;
    std::cout << "Initial velocity in G frame (km/s): " << v_com_in_G.transpose()/km << std::endl;

    Eigen::Vector3d thetas(-1.0, 0.2, 1.0); // rad -- Euler Roll Pitch Yaw angles 
    Eigen::Vector3d theta_dots(0*DEG2RAD, 0*DEG2RAD, 0*DEG2RAD); // rad/s

    Eigen::Matrix3d R_G_from_B = euler_to_rotation_matrix(thetas);
    Eigen::Matrix3d J_omega = euler_rates_omega_jac(thetas);

    Eigen::Vector3d omega_B_in_G = J_omega * theta_dots;


    simState.s_COM_in_G = s_com_in_G;
    simState.v_COM_in_G = v_com_in_G;
    simState.R_G_from_B = R_G_from_B;
    simState.omega_B_in_G = omega_B_in_G;

    cubesatState = getSatState(simState);
}


void waitForPythonReadyMessage()
{
    /** TODO: Loop LCM receiving until Python sends a ready message*/
}


void doSimLoop()
{
    LSODA solver;

    using clock = std::chrono::steady_clock;

    auto wall_start_sim = clock::now();

    std::cout << "Starting simulation loop..." << std::endl;

    Eigen::Vector3d s_sol_ECI;

    double t = 0;
    while (t < t_final_sim)
    {
        auto wall_start = clock::now();
        
        if (t - t_last_data >= dt_data)
        {
            t_last_data = t;

            /** TODO: Calculate incident flux stuff here. */

            /** TODO: Send flux data to Python Here */
        }

        // Do simulation step
        ControlInput_t u;
        u.inputForce_in_G = Eigen::Vector3d::Zero();
        u.inputTorque_in_G = Eigen::Vector3d::Zero();

        simState = doSimStep(&solver, t, dt, simState, u, simParams);
        cubesatState = getSatState(simState);

        s_sol_ECI = get_sun_location_eci(t, simParams.ephemeris);
        

        
        // Log time and completion percentage
        auto wall_now = clock::now();
        double wall_elapsed = std::chrono::duration<double>(wall_now - wall_start).count();

        double speed_of_sim = simParams.time_step / wall_elapsed;

        

        std::cout << '\r' << "Sim percent: " << std::fixed << std::setprecision(2) << (t/simParams.t_final)*100 << "%, Speed of sim: " << std::fixed  << std::setprecision(2) << speed_of_sim << "x real-time     " << std::flush;

        t += dt;
    }
    std::cout << std::endl;
    std::cout << s_sol_ECI.transpose()/AU << ", " << s_sol_ECI.norm()/AU << std::endl;
}


/*
This function is your executable function, it is the function that gets called when you run your code.
*/
int main()
{
    setSimParams();
    setInitialOrbitState();

    waitForPythonReadyMessage();

    doSimLoop();

    std::cout << "Simulation Finished" << std::endl;

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


    clearMem();
    // leave at bottom of main - shows no errors
    return 0;
    
}