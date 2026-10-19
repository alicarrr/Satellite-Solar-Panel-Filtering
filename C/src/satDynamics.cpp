#include "satDynamics.hpp"
#include "satParams.hpp"
#include "LSODA.hpp"
#include "RK78.hpp"
#include "vectorMaths.hpp"

#include <vector>
#include <iostream>
#include <stdexcept>

StateVector SimState_to_StateVector(SimState_t simState)
{
    StateVector state;

    state.segment(0, 3) = simState.s_COM_in_G;
    state.segment(3, 3) = simState.R_G_from_B.col(0);
    state.segment(6, 3) = simState.R_G_from_B.col(1);
    state.segment(9, 3) = simState.v_COM_in_G;
    state.segment(12, 3) = simState.omega_B_in_G;

    return state;
}


OrbitStateVector SimState_to_OrbitStateVector(SimState_t simState)
{
    OrbitStateVector state;

    state.segment(0, 3) = simState.s_COM_in_G;
    state.segment(3, 3) = simState.v_COM_in_G;

    return state;
}


OrbitStateVector SimState_to_OrbitStateVector(OrbitSimState_t simState)
{
    OrbitStateVector state;

    state.segment(0, 3) = simState.s_COM_in_G;
    state.segment(3, 3) = simState.v_COM_in_G;

    return state;
}


RigidBodyStateVector SimState_to_RigidBodyStateVector(SimState_t simState)
{
    RigidBodyStateVector state;

    state.segment(0, 3) = simState.R_G_from_B.col(0);
    state.segment(3, 3) = simState.R_G_from_B.col(1);
    state.segment(6, 3) = simState.omega_B_in_G;

    return state;
}


RigidBodyStateVector SimState_to_RigidBodyStateVector(RigidBodySimState_t simState)
{
    RigidBodyStateVector state;

    state.segment(0, 3) = simState.x_hat_B_in_G;
    state.segment(3, 3) = simState.y_hat_B_in_G;
    state.segment(6, 3) = simState.omega_B_in_G;

    return state;
}


SimState_t StateVector_to_SimState(StateVector simState)
{
    SimState_t state;

    Eigen::Vector3d x_hat_B_in_G = simState.segment(3, 3);
    Eigen::Vector3d y_hat_B_in_G = simState.segment(6, 3);

    state.s_COM_in_G = simState.segment(0, 3);
    state.R_G_from_B.col(0) = x_hat_B_in_G;
    state.R_G_from_B.col(1) = y_hat_B_in_G;
    state.R_G_from_B.col(2) = x_hat_B_in_G.cross(y_hat_B_in_G);
    state.v_COM_in_G = simState.segment(9, 3);
    state.omega_B_in_G = simState.segment(12, 3);

    return state;
}


OrbitSimState_t StateVector_to_SimState(OrbitStateVector simState)
{
    OrbitSimState_t state;

    state.s_COM_in_G = simState.segment(0, 3);
    state.v_COM_in_G = simState.segment(3, 3);

    return state;
    
}


RigidBodySimState_t StateVector_to_SimState(RigidBodyStateVector simState)
{
    RigidBodySimState_t state;

    state.x_hat_B_in_G = simState.segment(0, 3);
    state.y_hat_B_in_G = simState.segment(3, 3);
    state.omega_B_in_G = simState.segment(6, 3);

    return state;
}

SimState_t combineSimulationStates(OrbitSimState_t orbitState, RigidBodySimState_t rigidBodyState)
{
    SimState_t state;

    state.s_COM_in_G = orbitState.s_COM_in_G;
    state.R_G_from_B.col(0) = rigidBodyState.x_hat_B_in_G;
    state.R_G_from_B.col(1) = rigidBodyState.y_hat_B_in_G;
    state.R_G_from_B.col(2) = rigidBodyState.x_hat_B_in_G.cross(rigidBodyState.y_hat_B_in_G);
    state.v_COM_in_G = orbitState.v_COM_in_G;
    state.omega_B_in_G = rigidBodyState.omega_B_in_G;

    return state;
}


ControlInput_t inputVector_to_ControlInput(InputVector u)
{
    ControlInput_t input;

    input.inputForce_in_G = u.segment(0, 3);
    input.inputTorque_in_G = u.segment(3, 3);

    return input;
}


InputVector ControlInput_to_InputVector(ControlInput_t u)
{
    InputVector input;

    input.segment(0, 3) = u.inputForce_in_G;
    input.segment(3, 3) = u.inputTorque_in_G;

    return input;
}


OrbitStateVector orbitDynamics(OrbitStateVector state, Eigen::Vector3d inputForce_in_G, double t, const SimParams_t& simParams)
{
    OrbitSimState_t simState = StateVector_to_SimState(state);

    Eigen::Vector3d g_in_G = simParams.ephemeris->earthGravityFieldJ2000Mps2(simState.s_COM_in_G, simParams.initialUtc, t);

    OrbitSimState_t simStateDerivative;
    simStateDerivative.s_COM_in_G = simState.v_COM_in_G;
    simStateDerivative.v_COM_in_G = g_in_G + inputForce_in_G / cubesat.mass;

    OrbitStateVector stateDerivative = SimState_to_OrbitStateVector(simStateDerivative);

    return stateDerivative;
}


RigidBodyStateVector rotationDynamics(RigidBodyStateVector state, Eigen::Vector3d inputTorque_in_G, double t, const SimParams_t& simParams)
{
    RigidBodySimState_t simState = StateVector_to_SimState(state);

    RigidBodySimState_t simStateDerivative;

    Eigen::Matrix3d R_G_from_B;
    R_G_from_B.col(0) = simState.x_hat_B_in_G;
    R_G_from_B.col(1) = simState.y_hat_B_in_G;
    R_G_from_B.col(2) = simState.x_hat_B_in_G.cross(simState.y_hat_B_in_G);

    Eigen::Matrix3d J_B_in_B = cubesat.inertia_tensor;

    Eigen::Vector3d omega_B_in_B = R_G_from_B.transpose() * simState.omega_B_in_G;

    Eigen::Vector3d H_B_in_B = J_B_in_B * omega_B_in_B;

    double k_p = 1e-6; // Constant to pull the Rotation matrix back to a proper rotation matrix from numerical error build up

    simStateDerivative.x_hat_B_in_G = simState.omega_B_in_G.cross(simState.x_hat_B_in_G) + k_p * (1 - simState.x_hat_B_in_G.norm()) * simState.x_hat_B_in_G;
    simStateDerivative.y_hat_B_in_G = simState.omega_B_in_G.cross(simState.y_hat_B_in_G) + k_p * (1 - simState.y_hat_B_in_G.norm()) * simState.y_hat_B_in_G;

    simStateDerivative.omega_B_in_G = R_G_from_B * J_B_in_B.inverse() * ((R_G_from_B.transpose() * inputTorque_in_G) - omega_B_in_B.cross(H_B_in_B));

    RigidBodyStateVector stateDerivative = SimState_to_RigidBodyStateVector(simStateDerivative);

    return stateDerivative;
}


void getOrbitSimDerivatives(double t, double* x, double* x_dot, void* params)
{
    OrbitSimParams_t* orbitSimParams = static_cast<OrbitSimParams_t*>(params);

    Eigen::Vector3d inputForce_in_G = orbitSimParams->inputForce_in_G;
    SimParams_t simParams = orbitSimParams->simParams;

    OrbitStateVector state;
    for (int i = 0; i < n_orbit_state; i++)
    {
        state(i) = x[i];
    }

    OrbitStateVector stateDerivative = orbitDynamics(state, inputForce_in_G, t, simParams);

    for (int i = 0; i < n_orbit_state; i++)
    {
        x_dot[i] = stateDerivative(i);
    }
}


void getRotationSimDerivatives(double t, double* x, double* x_dot, void* params)
{
    RotationSimParams_t* rotationSimParams = static_cast<RotationSimParams_t*>(params);

    Eigen::Vector3d inputTorque_in_G = rotationSimParams->inputTorque_in_G;
    SimParams_t simParams = rotationSimParams->simParams;

    RigidBodyStateVector state;
    for (int i = 0; i < n_rotate_state; i++)
    {
        state(i) = x[i];
    }

    RigidBodyStateVector stateDerivative = rotationDynamics(state, inputTorque_in_G, t, simParams);

    for (int i = 0; i < n_rotate_state; i++)
    {
        x_dot[i] = stateDerivative(i);
    }
}


double rk_error_estimate(Eigen::VectorXd atol, double rtol, Eigen::VectorXd x, Eigen::VectorXd error, double h)
{
    Eigen::VectorXd scale = atol + rtol * x.cwiseAbs();
    double error_norm = (error.cwiseQuotient(scale)).norm() / sqrt(x.size());
    return error_norm;
}


SimState_t doSimStep(LSODA* solver, double t, double dt, SimState_t state, ControlInput_t u, const SimParams_t& simParams)
{
    double t0 = t;
    double t0_rot = t;
    double tf = t + dt;

    OrbitStateVector orbitStateVector = SimState_to_OrbitStateVector(state);
    RigidBodyStateVector rotationStateVector = SimState_to_RigidBodyStateVector(state);

    OrbitSimParams_t orbitSimParams;
    RotationSimParams_t rotationSimParams;

    orbitSimParams.inputForce_in_G = u.inputForce_in_G;
    orbitSimParams.simParams = simParams;

    rotationSimParams.inputTorque_in_G = u.inputTorque_in_G;
    rotationSimParams.simParams = simParams;

    double x_orbit_0[n_orbit_state];
    for (int i = 0; i < n_orbit_state; i++)
    {
        x_orbit_0[i] = orbitStateVector(i);
    }

    std::vector<double> x_rotation = to_std_vector(rotationStateVector);

    double x_orbit_next[n_orbit_state];
    std::vector<double> y_rotation_next;

    // Do Rotation Step with LSODA stiff integrator
    int is_valid = 1;

    solver->lsoda_update(getRotationSimDerivatives, n_rotate_state, x_rotation, y_rotation_next, &t0_rot, tf, &is_valid, &rotationSimParams, rtol_lsoda, atol_lsoda);

    if (is_valid <= 0)
    {
        std::cerr << "LSODA integration failed at time " << t0 << "s, State was " << std::endl;
        std::cerr << SimState_to_StateVector(state).transpose() << std::endl;

        throw std::runtime_error("LSODA failed");
    }

    std::vector<double> x_rotation_next(n_rotate_state);
    // for (size_t i = 0; i < y_rotation_next.size(); i++)
    for (size_t i = 0; i < n_rotate_state; i++)
    {
        x_rotation_next[i] = y_rotation_next[i+1];
    }

    RigidBodyStateVector rotationStateVectorNext = RigidBodyStateVector(x_rotation_next.data());

    // Now do orbit sim step with RK78, loop to reduce error
    // Advance the orbit to tf using accepted RK78 substeps.
    double dt_current = dt;
    while (t0 < tf)
    {
        // Limit the final substep to the remaining interval.
        if (dt_current > tf - t0)
        {
            dt_current = tf - t0;
        }

        double x_orbit_err[n_orbit_state];

        OrbitStateVector orbitStateErr = OrbitStateVector::Zero();

        int status = rkf78_step(t0, dt_current, x_orbit_0, n_orbit_state, &orbitSimParams, getOrbitSimDerivatives, x_orbit_next, x_orbit_err);

        for (int i = 0; i < n_orbit_state; i++)
        {
            orbitStateErr(i) = x_orbit_err[i];
        }

        double err = rk_error_estimate(atol_rk_vec, rtol_rk, orbitStateVector, orbitStateErr, dt_current);
        if (err < 1.0)
        {
            // Continue from the accepted state and time.
            for (int i = 0; i < n_orbit_state; i++)
            {
                x_orbit_0[i] = x_orbit_next[i];
                orbitStateVector(i) = x_orbit_next[i];
            }
            t0 += dt_current;
        }
        else
        {
            dt_current = dt_current * 0.8 * pow(1.0 / err, 1.0 / 8.0);
            if (dt_current < 1e-6)
            {
                throw std::runtime_error("RK78 failed to converge, dt_current is too small.");
            }
        }
    }

    OrbitStateVector orbitStateVectorNext = OrbitStateVector(x_orbit_next);

    RigidBodySimState_t rotationStateNext = StateVector_to_SimState(rotationStateVectorNext);
    OrbitSimState_t orbitStateNext = StateVector_to_SimState(orbitStateVectorNext);

    SimState_t stateNext = combineSimulationStates(orbitStateNext, rotationStateNext);

    return stateNext;
}


CubesatState_t getSatState(SimState_t simState)
{
    CubesatState_t cubesatState;

    cubesatState.position_COM_ECI = simState.s_COM_in_G;
    cubesatState.position_GC_ECI = simState.s_COM_in_G - simState.R_G_from_B * cubesat.position_COM_rel_GC;
    cubesatState.rotation_matrix_ECI_from_B = simState.R_G_from_B;
    cubesatState.velocity_COM_ECI = simState.v_COM_in_G;
    cubesatState.angular_velocity_ECI = simState.omega_B_in_G;

    return cubesatState;
}