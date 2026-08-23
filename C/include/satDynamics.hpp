#pragma once

#include <eigen3/Eigen/Dense>

// Will use rotation matrices for cubesat movement at an orbital position

// Can take a TLE or orbital elements and convert that to an initial position and velocity (Space 1/2 Code)
// Will integrate with RK8 the gravity equation plus perturbations for orbital position and velocity

typedef struct {
    Eigen::Vector3d s_COM_G; // COM position in global frame
    Eigen::Vector3d s_GC_G; // GC position in global frame
    Eigen::Matrix3d R_G_from_B; // Rotation matrix from body frame to global frame
} satPosState_t;
