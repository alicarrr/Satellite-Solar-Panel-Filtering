#pragma once

#include <eigen3/Eigen/Dense>
#include "lib/spice/cspice/include/SpiceUsr.h"
#include "ephem.hpp"
#include "transforms.hpp"
#include "simParams.hpp"

#include <tuple>

inline constexpr int n_state = 15;
inline constexpr int n_orbit_state = 6;
inline constexpr int n_rotate_state = 9;

typedef Eigen::Vector<double, n_state> StateVector;

typedef Eigen::Vector<double, n_orbit_state> OrbitStateVector;
typedef Eigen::Vector<double, n_rotate_state> RigidBodyStateVector;


// Will use rotation matrices for cubesat movement at an orbital position

// Can take a TLE or orbital elements and convert that to an initial position and velocity (Space 1/2 Code)
// Will integrate with RK8 the gravity equation plus perturbations for orbital position and velocity

typedef struct {
    Eigen::Vector3d s_COM_G; // COM position in global frame
    Eigen::Vector3d s_GC_G; // GC position in global frame
    Eigen::Matrix3d R_G_from_B; // Rotation matrix from body frame to global frame
} satPosState_t;


typedef struct {
    Eigen::Vector3d s_COM_in_G;
    Eigen::Vector3d x_hat_B_in_G;
    Eigen::Vector3d y_hat_B_in_G;
    Eigen::Vector3d v_COM_in_G;
    Eigen::Vector3d omega_B_in_G;

} SimState_t;

typedef struct {
    Eigen::Vector3d s_COM_in_G;
    Eigen::Vector3d v_COM_in_G;

} OrbitSimState_t;

typedef struct {
    Eigen::Vector3d x_hat_B_in_G;
    Eigen::Vector3d y_hat_B_in_G;
    Eigen::Vector3d omega_B_in_G;

} RigidBodySimState_t;

/**
 * @brief Convert SimState_t struct to a state vector for the simulation
 * 
 * @param simState The sim state as a SimState_t struct
 * 
 * @return The sim state as a StateVector (Eigen Vector)
 */
StateVector SimState_to_StateVector(SimState_t simState);

/**
 * @brief Convert SimState_t struct to a orbit state vector for the simulation
 * 
 * @param simState The sim state as a SimState_t struct
 * 
 * @return The orbit sim state as a OrbitStateVector (Eigen Vector)
 */
OrbitStateVector SimState_to_OrbitStateVector(SimState_t simState);

/**
 * @brief Convert OrbitSimState_t struct to a orbit state vector for the simulation
 * 
 * @param simState The orbit sim state as a OrbitSimState_t struct
 * 
 * @return The orbit sim state as a OrbitStateVector (Eigen Vector)
 */
OrbitStateVector SimState_to_OrbitStateVector(OrbitSimState_t simState);

/**
 * @brief Convert SimState_t struct to a rigid body rotation state vector for the simulation
 * 
 * @param simState The sim state as a SimState_t struct
 * 
 * @return The sim state as a RigidBodyStateVector (Eigen Vector)
 */
RigidBodyStateVector SimState_to_RigidBodyStateVector(SimState_t simState);


/**
 * @brief Convert RigidBodySimState_t struct to a rigid body rotation state vector for the simulation
 * 
 * @param simState The rigid body sim state as a RigidBodySimState_t struct
 * 
 * @return The sim state as a RigidBodyStateVector (Eigen Vector)
 */
RigidBodyStateVector SimState_to_RigidBodyStateVector(RigidBodySimState_t simState);


/**
 * @brief Convert sim state vector to a simState_t Struct
 * 
 * @param simState The simulation state as a StateVector
 * 
 * @return The sim state as a SimState_t struct
 */
SimState_t StateVector_to_SimState(StateVector simState);

/**
 * @brief Convert orbit sim state vector to a OrbitSimState_t Struct
 * 
 * @param simState The orbit simulation state as a OrbitStateVector
 * 
 * @return The sim state as a OrbitSimState_t struct
 */
OrbitSimState_t StateVector_to_SimState(OrbitStateVector simState);

/**
 * @brief Convert rigid body rotation state vector to a RigidBodySimState_t Struct
 * 
 * @param simState The simulation state as a RigidBodyStateVector
 * 
 * @return The rigid body rotation sim state as a RigidBodySimState_t struct
 */
RigidBodySimState_t StateVector_to_SimState(RigidBodyStateVector simState);

/**
 * @brief Combine the two simulation states into one state object
 * 
 * @param orbitState The orbit state as an OrbitSimState_t struct
 * @param rigidBodyState The rigid body state as a RigidBodySimState_t struct
 * 
 * @return simulation state as a SimState_t struct
 */
SimState_t combineSimulationStates(OrbitSimState_t orbitState, RigidBodySimState_t rigidBodyState);


/**
 * @brief Calculate the derivative of the orbit simulation state
 * 
 * @param state The current simulation state to evaluate the derivative at
 * @param simParams The simulation parameters as a SimParams_t struct
 * 
 * @return The state derivative as a vector
 */
OrbitStateVector orbitDynamics(OrbitStateVector state, const SimParams_t& simParams);

/**
 * @brief Calculate the derivative of the rigid body simulation state
 * 
 * @param state The current simulation state to evaluate the derivative at
 * @param simParams The simulation parameters as a SimParams_t struct
 * 
 * @return The state derivative as a vector
 */
RigidBodyStateVector rotationDynamics(RigidBodyStateVector state, const SimParams_t& simParams);

