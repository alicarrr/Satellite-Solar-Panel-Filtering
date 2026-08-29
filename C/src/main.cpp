// Your main file doesn't have to be called main, but it must contain the main function.
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

#include <netcdf>


#include "simParams.hpp"

#include "satDynamics.hpp"
#include "LSODA.hpp"




SimParams_t simParams;

SimState_t simState;
CubesatState_t cubesatState;

const double dt_data = 1.0; 
const double dt = 1e-1; // Sim Time Step
const std::string date_0 = "2025-01-09T01:01:10.050000";

const double t_final_sim = 0.5 * second;

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
    albedo_file.close();
}


void setInitialOrbitState()
{
    double period = 95.2 * minute; // seconds

    double mean_motion = (2.0 * M_PI) / period;

    double e = 0.002; // Eccentricity
    double M = 0.15 * M_PI; // Mean anomaly at epoch

    double argp = 30 * DEG2RAD;
    double raan = -90 * DEG2RAD;
    double inc = 60 * DEG2RAD;


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


void preLoadData()
{
    std::cout << "Reading AM0 and AM1.5 spectra from file: " << AM_spectra_file_name << std::endl;
    Spectrum_t AM0 = readASTMGFileForAM0(AM_spectra_file_name);
    Spectrum_t AM1_5 = readASTMGFileForAM1_5(AM_spectra_file_name);

    std::vector<double> AM0_wavelengths_nm = to_std_vector(AM0.wavelengths);
    std::vector<double> AM1_5_wavelengths_nm = to_std_vector(AM1_5.wavelengths);
    std::vector<double> AM0_irradiance = to_std_vector(AM0.values);
    std::vector<double> AM1_5_irradiance = to_std_vector(AM1_5.values);

    std::vector<double> wavelengths_nm_vec = calculateBinCentreWavelengths(min_wavelength_nm, max_wavelength_nm, d_wavelength_nm);

    SpectraVector_t wavelengths_nm = Eigen::Map<SpectraVector_t>(wavelengths_nm_vec.data(), wavelengths_nm_vec.size());

    SpectraVector_t AM0_irradiances_binned = Eigen::Map<SpectraVector_t>(rebinSpectrumConservingFlux(AM0_wavelengths_nm, AM0_irradiance, min_wavelength_nm, max_wavelength_nm, d_wavelength_nm, am_0_spectra_extend).data(), wavelengths_nm_vec.size());

    SpectraVector_t AM1_5_irradiances_binned = Eigen::Map<SpectraVector_t>(rebinSpectrumConservingFlux(AM1_5_wavelengths_nm, AM1_5_irradiance, min_wavelength_nm, max_wavelength_nm, d_wavelength_nm, am_1_5_spectra_extend).data(), wavelengths_nm_vec.size());

    AM0.wavelengths = wavelengths_nm;
    AM0.values = AM0_irradiances_binned;

    AM1_5.wavelengths = wavelengths_nm;
    AM1_5.values = AM1_5_irradiances_binned;

    std::cout << "Arrays are of size: " << AM0.wavelengths.size() << std::endl;
    std::cout << "Arrays end with " << AM0.values(N_samples_spectra - 2) << std::endl;

    Spectrum_t atmosphereTransmittance = opticalTransmittanceAtmosphereFromData(AM0, AM1_5);

    std::cout << "Transmittance ends with " << AM0.values(N_samples_spectra - 2) << std::endl;


    // get hamster albedo spectrum at a given location on Earth
    std::cout << "Reading albedo data from file: " << std::endl;
    read_albedo_data(hamster_day_of_year);

    std::cout << "Reading Atmosphere Data" << std::endl;
    getTransmittanceFromFile(atmosphere_file_path);

    // printVariableInformation("total_downward_transmissivity", atmosphere_transmissivity_file);

    

    // Calculate transmission ratio from AM0 and AM1.5 spectra
    Spectrum_t transmission_ratio_AM1_5 = atmosphereTransissionRatio(AM0, AM1_5);


    // Testing plotting
    FILE* gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot == nullptr)
    {
        std::cerr << "Could not open gnuplot." << std::endl;
        return;
    }
    std::fprintf(gnuplot,
                    "set title 'Albedo and reflected radiance spectra'\n"
                    "set xlabel 'Wavelength'\n"
                    "set ylabel 'Spectral value'\n"
                    "set grid\n"
                    "plot '-' using 1:2 with lines title 'AM0', "
                    "'-' using 1:2 with lines title 'AM1.5', "
                    "'-' using 1:2 with lines title 'Atmosphere Transmittance', "
                    "'-' using 1:2 with lines title 'Reflected radiance'\n");

    for (Eigen::Index i = 0; i < AM0.wavelengths.size(); ++i)
        std::fprintf(gnuplot, "%.17g %.17g\n", AM0.wavelengths(i), AM0.values(i));
    std::fprintf(gnuplot, "e\n");

    for (Eigen::Index i = 0; i < AM1_5.wavelengths.size(); ++i)
        std::fprintf(gnuplot, "%.17g %.17g\n", AM1_5.wavelengths(i), AM1_5.values(i));
    std::fprintf(gnuplot, "e\n");

    std::fflush(gnuplot);
    pclose(gnuplot);


    // Testing plotting
    gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot == nullptr)
    {
        std::cerr << "Could not open gnuplot." << std::endl;
        return;
    }
    std::fprintf(gnuplot,
                    "set title 'Albedo and reflected radiance spectra'\n"
                    "set xlabel 'Wavelength'\n"
                    "set ylabel 'Spectral value'\n"
                    "set grid\n"
                    "plot '-' using 1:2 with lines title 'Atmosphere Transmittance'\n");

    for (Eigen::Index i = 0; i < atmosphereTransmittance.wavelengths.size(); ++i)
        std::fprintf(gnuplot, "%.17g %.17g\n", atmosphereTransmittance.wavelengths(i), atmosphereTransmittance.values(i));
    std::fprintf(gnuplot, "e\n");

    std::fflush(gnuplot);
    pclose(gnuplot);
}


void doSimLoop()
{
    preLoadData();

    LSODA solver;

    using clock = std::chrono::steady_clock;

    auto wall_start_sim = clock::now();

    std::cout << "Starting simulation loop..." << std::endl;

    Eigen::Vector3d s_sol_ECI;
    // Get initial sun position in ECI frame at t=0
    s_sol_ECI = get_sun_location_eci(0, simParams.ephemeris);

    double t = 0;
    while (t < t_final_sim)
    {
        auto wall_start = clock::now();
        
        if (t - t_last_data >= dt_data)
        {
            t_last_data = t;

            /** TODO: Calculate incident flux stuff here. */

            // test albedo toa calculation with one position in orbit and one visible albedo patch. This is a test to see if the function works and returns a reasonable spectrum.
            

            

            // AlbedoData albedo_data = read_albedo_data();
            // take current llh of satellite and get albedo spectrum at that location
            Eigen::Vector3d sat_pos_ECEF = simParams.ephemeris->eciToEcef(cubesatState.position_GC_ECI, t);
            Eigen::Vector3d sat_llh = ecef_to_llh(sat_pos_ECEF);
            // Eigen::VectorXd albedo_values = get_albedo_spectrum(albedo_data, sat_llh(0), sat_llh(1));
            
            // Spectrum_t albedo_spectrum_at_loc = 
            // {
            //     .wavelengths = albedo_data.wavelength,
            //     .values = albedo_values
            // };

            // test the interpolation function by interpolating the albedo spectrum to the wavelengths of the AM1.5 spectrum
            // Spectrum_t albedo_spectrum_at_loc_interp = interpolateSpectrum(albedo_spectrum_at_loc, AM1_5.wavelengths);

            Eigen::Vector3d albedo_patch_pos_llh = Eigen::Vector3d(sat_llh(0), sat_llh(1), 0); // lat, lon, alt= 6371000 m (Earth radius)

            Eigen::Vector3d albedo_patch_pos_ecef = llh_to_ecef(albedo_patch_pos_llh); 
            Eigen::Vector3d albedo_patch_pos_eci = simParams.ephemeris->ecefToEci(albedo_patch_pos_ecef, t);
            // convert to ECI at t=0
            Eigen::Vector3d albedo_patch_norm_vector = albedo_patch_pos_eci.cwiseProduct(1.0/albedo_patch_pos_eci.norm() * Eigen::Vector3d::Ones()); // normal vector pointing up from the surface of the Earth

            Eigen::Vector3d panel_pos_G = panelPosInG(cubesat.panels[0].position_rel_GC_in_B, cubesatState.rotation_matrix_ECI_from_B, cubesatState.position_GC_ECI);
            Eigen::Vector3d sun_pos_G = s_sol_ECI;
            
            

            std::cout << "Entering function" << std::endl;
            // Spectrum_t test_spectra = albedoRadianceTOA(transmission_ratio_AM1_5, AM0, albedo_spectrum_at_loc, albedo_patch_pos_eci, albedo_patch_norm_vector, panel_pos_G, sun_pos_G);

            // Plot the albedo and reflected spectra using inline gnuplot data.
            // FILE* gnuplot = popen("gnuplot -persistent", "w");
            // if (gnuplot == nullptr)
            // {
            //     std::cerr << "Could not open gnuplot." << std::endl;
            //     return;
            // }

            

            // // calculate earth ir spectrum at the patch and plot it
            // Spectrum_t earth_ir_spectrum = EarthIRRadianceTOA(transmission_ratio_AM1_5, AM0, albedo_patch_pos_eci, panel_pos_G);

            // // get black body spectrum at earth temp
            // double earth_temp_K = 288.0; // K
            // Spectrum_t earth_blackbody_spectrum = blackBodySpectrum(earth_temp_K, AM0);

            // gnuplot = popen("gnuplot -persistent", "w");
            // if (gnuplot == nullptr)
            // {
            //     std::cerr << "Could not open gnuplot." << std::endl;
            //     return;
            // }
            // std::fprintf(gnuplot,
            //              "set title 'Earth IR radiance spectrum'\n"
            //              "set xlabel 'Wavelength'\n"
            //              "set ylabel 'Spectral value'\n"
            //              "set grid\n"
            //              "plot '-' using 1:2 with lines title 'Earth IR radiance'\n"
            //              "'-' using 1:2 with lines title 'Earth black body'\n");
            // for (Eigen::Index i = 0; i < earth_ir_spectrum.wavelengths.size(); i++)
            // {
            //     std::fprintf(gnuplot, "%.17g %.17g\n",
            //                  earth_ir_spectrum.wavelengths(i),
            //                  earth_ir_spectrum.values(i));
            // }
            // // plot black body spectrum for comparison
            // for (Eigen::Index i = 0; i < earth_blackbody_spectrum.wavelengths.size(); i++)
            // {
            //     std::fprintf(gnuplot, "%.17g %.17g\n",
            //                  earth_blackbody_spectrum.wavelengths(i),
            //                  earth_blackbody_spectrum.values(i));
            // }
            // std::fprintf(gnuplot, "e\n");
            // std::fflush(gnuplot);
            // pclose(gnuplot);
            

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
    double max_wavelength = 50000.0;
    double wavelength_step_size = 10.0;

    Spectrum_t binned_spectrum = binSpectrum(spectrum, min_wavelength, max_wavelength, wavelength_step_size);


    clearMem();
    // leave at bottom of main - shows no errors
    return 0;
    
}