// Headers have header guards which prevent infinitely including code
#pragma once

// Include Eigen Vectors and Matrices
#include <eigen3/Eigen/Dense>

// Header files contain function prototypes. They tell you the name, inputs and output of the function but not the code for the function.

// Comments like the one below show up when you hover a function, they are Doxygen format like docstrings in python.

/**
* @brief This function converts from polar coordinates to cartesian where the input
*
* @param r_polar The vector of (r, theta)
*
* @return The cartesian coordinates (x, y)
**/
Eigen::Vector2d polarToCartesian(Eigen::Vector2d r_polar);