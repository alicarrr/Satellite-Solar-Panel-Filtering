#include "satDynamics.hpp"

StateVector SimState_to_StateVector(SimState_t simState)
{
    StateVector state;

    state.segment(0, 3) = simState.s_COM_in_G;
    state.segment(3, 3) = simState.x_hat_B_in_G;
    state.segment(6, 3) = simState.y_hat_B_in_G;
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

    state.segment(0, 3) = simState.x_hat_B_in_G;
    state.segment(3, 3) = simState.y_hat_B_in_G;
    state.segment(6, 3) = simState.omega_B_in_G;

    return state;
}