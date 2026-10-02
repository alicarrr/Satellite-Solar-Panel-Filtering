// Your main file doesn't have to be called main, but it must contain the main function.
#include "spectra.hpp"
#include "albedo.hpp"
#include "spectra.hpp"
#include "satParams.hpp"
#include "vectorMaths.hpp"

#include <lcm/lcm-cpp.hpp>
#include "simMessages/sim_data_t.hpp"
#include "simMessages/sim_ready_t.hpp"
#include "simMessages/spectrum_metadata_t.hpp"
#include "simMessages/spectrum_data_t.hpp"

#include <vector>
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


bool ready = false;

inline constexpr bool sendDynamics = false;

Eigen::Vector3d s_sol_ECI;

SimParams_t simParams;

SimState_t simState;
CubesatState_t cubesatState;

std::vector<SpectraVector_t> panelPowerState(cubesat.panels.size(), SpectraVector_t::Zero());

const double dt_data = 10.0; 
const double dt = 1e-1; // Sim Time Step
const std::string date_0 = "2025-01-09T01:01:10.050000";

const double t_final_sim = 96*60*60; // seconds

double t_last_data = -2 * dt_data;

double t = 0;


class ReadyHandler
{
    public:
        ~ReadyHandler() {}

        void handleMessage(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const simMessages::sim_ready_t* msg)
        {
            if (msg->is_ready)
            {
                ready = true;
                std::cout << "Python Listener is Ready" << std::endl;
            }
        }
};

void publishSimDone(lcm::LCM* lcm)
{
    simMessages::sim_ready_t simDoneMsg;

    simDoneMsg.is_ready = true;

    lcm->publish("SimStatus", &simDoneMsg);
}

void publishSimData(lcm::LCM* lcm)
{
    if constexpr (sendDynamics)
    {
        simMessages::sim_data_t dynData;
        dynData.t = t;

        for (int i = 0; i < 3; i++)
        {
            dynData.pos_sat_eci[i] = simState.s_COM_in_G(i);
            dynData.pos_sun_eci[i] = s_sol_ECI(i);
            for (int j = 0; j < 3; j++)
            {
                dynData.R_eci_from_B[i][j] = simState.R_G_from_B(i, j);
            }
        }

        lcm->publish("DynamicsData", &dynData);
        
    }

    simMessages::spectrum_data_t simData;

    simData.t = t;
    simData.N_panel = cubesat.panels.size();
    simData.N_spectra = N_samples_spectra;

    simData.panelPowers.resize(simData.N_panel);

    for (int i = 0; i < simData.N_panel; i++)
    {
        simData.panelPowers[i].resize(N_samples_spectra);
    }

    for (int i = 0; i < simData.N_panel; i++)
    {
        for (int j = 0; j < N_samples_spectra; j++)
        {
            simData.panelPowers[i][j] = panelPowerState[i][j];
        }
    }

    lcm->publish("SpectraData", &simData);

}


void publishSimMetadata(lcm::LCM* lcm)
{
    simMessages::spectrum_metadata_t spectraData;

    spectraData.N_spectra = N_samples_spectra;
    spectraData.N_panel = cubesat.panels.size();
    spectraData.wavelengths.resize(N_samples_spectra);
    spectraData.panelAreas.resize(cubesat.panels.size());
    spectraData.panelTypes.resize(cubesat.panels.size());

    for (int i = 0; i < N_samples_spectra; i++)
    {
        spectraData.wavelengths[i] = spectra_wavelengths_nm(i);
    }

    for (int i = 0; i < spectraData.N_panel; i++)
    {
        spectraData.panelTypes[i] = cubesat.panels[i].panel_type;
        spectraData.panelAreas[i] = cubesat.panels[i].size;
    }

    lcm->publish("SpectraMetaData", &spectraData);


}


void calculateSpectralPowersOnSatellite()
{
    std::vector<SpectraVector_t> albedoPowerPerPanel = getTotalIncidentAlbedoAndIRPowerSpectra(simState.s_COM_in_G, simState.R_G_from_B, s_sol_ECI, t, &simParams);

    std::vector<SpectraVector_t> solarPowersPerPanel = solarPowerPerPanel(simState.s_COM_in_G, simState.R_G_from_B, s_sol_ECI);

    // double totalPower = 0.0;
    double totalAlbedoPower = 0.0;
    double totalSolarPower = 0.0;

    for (int i = 0; i < cubesat.panels.size(); i++)
    {
        panelPowerState[i] = albedoPowerPerPanel[i] + solarPowersPerPanel[i];
        
        // totalPower += panelPowerState[i].sum();
        totalAlbedoPower += albedoPowerPerPanel[i].sum();
        totalSolarPower += solarPowersPerPanel[i].sum();
    }
    
    // std::cout << "Total power on satellite at t=" << t << "s: " << totalPower*1000 << " mW" << std::endl;
    // std::cout << "Total IR and albedo power: " << totalAlbedoPower*1000 << " mW" << std::endl;
    // std::cout << "Total solar power: " << totalSolarPower*1000 << " mW" << std::endl;
}


void setSimParams()
{
    try {
        std::string lskPath = "C/lib/spice/kernels/naif0012.tls";
        std::string spkPath = "C/lib/spice/kernels/de442s.bsp";
        std::string earthHighPrecPath = "C/lib/spice/kernels/earth_000101_260601_260305.bpc";

        std::string earthGravPath = "C/lib/geographiclib/share/GeographicLib/gravity";
        std::string earthMagPath = "C/lib/geographiclib/share/GeographicLib/magnetic";

        simParams.ephemeris = new SolarSystemEphemeris(lskPath, spkPath, date_0);

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
    atmosphere_transmissivity_file.close();
}


void setInitialOrbitState()
{
    double period = 95.2 * minute; // seconds

    double mean_motion = (2.0 * M_PI) / period;

    double e = 0.002; // Eccentricity
    double M = -0.5 * M_PI; // Mean anomaly at epoch

    double argp = 30 * DEG2RAD;
    double raan = -90 * DEG2RAD;
    double inc = 10 * DEG2RAD;


    tuple<Eigen::Vector3d, Eigen::Vector3d> r_v_eci = angular_parameters_to_eci(M, e, mean_motion, inc, raan, argp, mu_earth);

    Eigen::Vector3d s_com_in_G = std::get<0>(r_v_eci); // m
    Eigen::Vector3d v_com_in_G = std::get<1>(r_v_eci); // m/s

    std::cout << "Initial position in G frame (km): " << s_com_in_G.transpose()/km << std::endl;
    std::cout << "Initial velocity in G frame (km/s): " << v_com_in_G.transpose()/km << std::endl;

    Eigen::Vector3d thetas(-1.0, 0.2, 1.0); // rad -- Euler Roll Pitch Yaw angles 
    Eigen::Vector3d theta_dots(3*DEG2RAD, 2*DEG2RAD, 0*DEG2RAD); // rad/s

    Eigen::Matrix3d R_G_from_B = euler_to_rotation_matrix(thetas);
    Eigen::Matrix3d J_omega = euler_rates_omega_jac(thetas);

    Eigen::Vector3d omega_B_in_G = J_omega * theta_dots;

    simState.s_COM_in_G = s_com_in_G;
    simState.v_COM_in_G = v_com_in_G;
    simState.R_G_from_B = R_G_from_B;
    simState.omega_B_in_G = omega_B_in_G;

    cubesatState = getSatState(simState);
}


void waitForPythonReadyMessage(lcm::LCM* lcm)
{
    ReadyHandler handler;
    lcm->subscribe("ListenerReady", &ReadyHandler::handleMessage, &handler);
    while (!ready)
    {
        lcm->handle();
    }
}


void preLoadData()
{
    setWavelengthArray();

    getAlbedoFileParams();

    // Get initial sun position in ECI frame at t=0
    s_sol_ECI = get_sun_location_eci(0, simParams.ephemeris);

    std::cout << d_lat_file << ", " << d_lon_file << std::endl;

    std::cout << latitudes_file.minCoeff() << ", " << latitudes_file.maxCoeff() << std::endl;
    std::cout << longitudes_file.minCoeff() << ", " << longitudes_file.maxCoeff() << std::endl;


    std::cout << "Reading Atmosphere Data" << std::endl;
    getTransmittanceFromFile(atmosphere_file_path);

    // get hamster albedo spectrum at a given location on Earth
    std::cout << "Reading albedo data from file: " << std::endl;
    read_albedo_data(hamster_day_of_year);

    std::cout << "Getting sun spectrum" << std::endl;
    getSolarSpectrumFromFile(sun_spectrum_path);

    double totalSunPower = sun_spectra.sum();
    std::cout << "Total Solar Power is " << totalSunPower << "W" << std::endl;
}


void doSimLoop(lcm::LCM* lcm)
{

    LSODA solver;

    using clock = std::chrono::steady_clock;

    auto wall_start_sim = clock::now();

    std::cout << "Starting simulation loop..." << std::endl;

    // Get initial sun position in ECI frame at t=0
    s_sol_ECI = get_sun_location_eci(0, simParams.ephemeris);

    
    while (t < t_final_sim)
    {
        auto wall_start = clock::now();
        
        if (t - t_last_data >= dt_data)
        {
            t_last_data = t;

            /** Calculate incident flux stuff here. */
            calculateSpectralPowersOnSatellite();
            
            /** Send flux data to Python Here */
            publishSimData(lcm);
            
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
        double wall_elapsed_sim = std::chrono::duration<double>(wall_now - wall_start_sim).count();

        double speed_of_sim = simParams.time_step / wall_elapsed;
        double total_speed_of_sim = t / wall_elapsed_sim;

        std::cout << '\r' << "Sim percent: " << std::fixed << std::setprecision(2) << (t/simParams.t_final)*100 << "%, Current speed of sim: " << std::fixed  << std::setprecision(2) << speed_of_sim << "x real-time     " << "Total speed of sim: " << std::fixed  << std::setprecision(2) << total_speed_of_sim << "x real-time     " << std::flush;

        t += dt;
    }
}


/*
This function is your executable function, it is the function that gets called when you run your code.
*/
int main()
{
    lcm::LCM lcm;

    if(!lcm.good())
    {
        return 1;
    }


    waitForPythonReadyMessage(&lcm);

    setSimParams();
    setInitialOrbitState();

    preLoadData();

    publishSimMetadata(&lcm);

    doSimLoop(&lcm);

    std::cout << "Simulation Finished" << std::endl;

    publishSimDone(&lcm);


    clearMem();
    // leave at bottom of main - shows no errors
    return 0;
    
}