#pragma once

// Will use rotation matrices for cubesat movement at an orbital position

// Can take a TLE or orbital elements and convert that to an initial position and velocity (Space 1/2 Code)
// Will integrate with RK8 the gravity equation plus perturbations for orbital position and velocity