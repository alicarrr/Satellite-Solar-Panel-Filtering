// Headers have header guards which prevent infinitely including code
#pragma once

// We want to use Pi so we include the standard C math header
#include <math.h>
// Include Eigen Vectors and Matrices
#include <eigen3/Eigen/Dense>


// Multiply degrees by this to convert to radians
inline constexpr double DEG2RAD = M_PI / 180.0f; 

// Multiply radians by this to convert to degrees
inline constexpr double RAD2DEG = 180.0f / M_PI;

// If there is no comment after a variable then the comments above it will appear when you hover it.

// An even better way for single constants is with constexpr, 
// it's effectively the same but instead of an unchangable variable it creates a constant value in your program.
// This means it behaves truly like a constant rather than an unchangable variable
inline constexpr double HOUR2MINUTE = 60;

// Vectors/arrays cannot be constexpr values, because they are lists in memory. 
// When you want to define a vector/list as constant you can define it as an inline const
inline const Eigen::Vector3d x_hat = Eigen::Vector3d(1, 0, 0);





inline constexpr double m = 1.0; // m
inline constexpr double nm = 1e-9; // m
inline constexpr double km = 1000.0; // m
inline constexpr double AU = 149597870700.0; // m


inline constexpr double c = 299792458.0; // m/s


inline constexpr double kg = 1.0; // kg
inline constexpr double K = 1.0; // Kelvin




inline constexpr double second = 1.0; // seconds
inline constexpr double minute = 60.0; // seconds
inline constexpr double hour = 60.0 * minute; // seconds
inline constexpr double day = 24.0 * hour; // seconds
inline constexpr double year = 365.24218967 * day; // seconds

inline constexpr double siderial_day = 23 * hour + 56 * minute + 4 * second; // seconds
inline constexpr double solar_day = 24 * hour; // seconds

inline constexpr double N = kg * m / (second * second); // Newton
inline constexpr double J = N * m; // Joule
inline constexpr double W = J / second; // Watt

inline constexpr double amp = 1.0; // A
inline constexpr double ohm = kg * m * m / (second * second * second * amp * amp); // Ohm
inline constexpr double volt = ohm * amp; // Volt
inline constexpr double henry = ohm * second; // Henry
inline constexpr double tesla = N / (amp * m); // Tesla
inline constexpr double farad = amp * second / volt; // Farad




inline constexpr double G = 6.67430e-11 * m*m*m / (kg * second*second); // Gravitational constant (m^3 kg^-1 s^-2)

inline constexpr double sigma_sb = 5.670374419e-8 * W / (m * m * (K*K*K*K)); // Stefan-Boltzmann constant (W m^-2 K^-4)

inline constexpr double M_sun = 1988500e24 * kg; // Mass of the Sun (kg)
inline constexpr double mu_sun = G * M_sun; // Gravitational parameter of the Sun (m^3 s^-2)
inline constexpr double L_sol = 3.828e26 * W; // Solar luminosity (W)
inline constexpr double R_sol = 695700e3 * m; // Radius of the Sun (m)

inline constexpr double M_earth = 5.97217e24 * kg; // Mass of Earth (kg)
inline constexpr double mu_earth = G * M_earth; // Gravitational parameter of Earth (m^3 s^-2)
inline constexpr double R_earth = 6371e3 * m; // Radius of Earth (m)
inline constexpr double albedo_earth = 0.306; // Average albedo of Earth
inline constexpr double emissivity_IR_earth = 0.95; // Emissivity of Earth in the infrared spectrum
inline constexpr double T_earth = 288.0 * K; // Average temperature of Earth (K)
inline constexpr double a_earth = 6378.137 * km;
inline constexpr double b_earth = 6356.7523142 * km;
inline constexpr double i_earth = 23.45 * DEG2RAD; // Inclination of Earth's axis
inline constexpr double angular_vel_earth = 2.0 * M_PI / siderial_day;

inline constexpr double lat_geomagnetic_north = 80.7 * DEG2RAD; // Geomagnetic north pole latitude
inline constexpr double lon_geomagnetic_north = -72.7 * DEG2RAD; // Geomagnetic north pole longitude
inline constexpr double magnetic_dipole_earth = 7.71e22 * amp * m * m; // Earth's magnetic dipole moment (A m^2)

inline constexpr double mu_0 = 4.0 * M_PI * 1e-7 * N / (amp * amp); // Vacuum permeability (N/A^2)

inline constexpr double K_sun = 93.0 / W; // Luminous Efficacy of the Sun (lm/W)
inline constexpr double lux_to_footcandela = 1.0/10.7639; // Conversion factor from lux to foot-candela (lm/m^2 to lm/ft^2)

inline constexpr double TT_delta = 69.184 * second; // Difference between Terrestrial Time and Universal Time (seconds)

inline constexpr double J2000_JD = 2451545;

inline constexpr double h = 6.62607015e-34 * J * second; // Planck's constant (J s)
inline constexpr double k = 1.380649e-23 * J / K; // Boltzmann's constant (J/K)