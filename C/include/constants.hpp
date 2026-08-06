// Headers have header guards which prevent infinitely including code
#pragma once

// We want to use Pi so we include the standard C math header
#include <math.h>
// Include Eigen Vectors and Matrices
#include <eigen3/Eigen/Dense>

// There are multiple ways of defining something as constant

// a #define replaces all instances of the word following (usually in all caps) everywhere in your code
#define RAD2DEG 180.0f / M_PI // We divide by 180.0f to ensure that we don't accidentally round down when dividing by an integer

// Sometimes a better way to define a constant is with the const keyword.
// When a variable is a const its value cannot be changed after the variable is created.
// The inline keyword means that every file that uses DEG2RAD will get a new copy of that constant value.
// The inline keyword is a bit less memory efficient but prevents multiple pieces of code from accidentally clashing over the same memory.
inline const double DEG2RAD = M_PI / 180.0f; // A comment here appears when you hover the variable

// If there is no comment after a variable then the comments above it will appear when you hover it.

// An even better way for single constants is with constexpr, 
// it's effectively the same but instead of an unchangable variable it creates a constant value in your program.
// This means it behaves truly like a constant rather than an unchangable variable
inline constexpr double HOUR2MINUTE = 60;

// Vectors/arrays cannot be constexpr values, because they are lists in memory. 
// When you want to define a vector/list as constant you can define it as an inline const
inline const Eigen::Vector3d x_hat = Eigen::Vector3d(1, 0, 0);