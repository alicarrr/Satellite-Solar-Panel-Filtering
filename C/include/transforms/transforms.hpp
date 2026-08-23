#pragma once

#include <eigen3/Eigen/Dense>
#include <tuple>

#include "time_conv.hpp"



Eigen::Matrix3d eci_to_ecliptic_matrix();

Eigen::Vector3d eci_to_ecliptic(const Eigen::Vector3d& r_eci);

std::tuple<Eigen::Vector3d, Eigen::Vector3d> eci_to_ecliptic(const Eigen::Vector3d& r_eci, const Eigen::Vector3d& v_eci);

Eigen::Vector3d ecliptic_to_eci(const Eigen::Vector3d& r_ecliptic);

std::tuple<Eigen::Vector3d, Eigen::Vector3d> ecliptic_to_eci(const Eigen::Vector3d& r_ecliptic, const Eigen::Vector3d& v_ecliptic);

Eigen::Vector3d ecef_to_llh(const Eigen::Vector3d& r_ecef);

Eigen::Vector3d llh_to_ecef(const Eigen::Vector3d& llh);

Eigen::Vector3d eci_to_ecef(const Eigen::Vector3d& r_eci, std::string date_utc);

std::tuple<Eigen::Vector3d, Eigen::Vector3d> eci_to_ecef(const Eigen::Vector3d& r_eci, const Eigen::Vector3d& v_eci, std::string date_utc);

Eigen::Vector3d ecef_to_eci(const Eigen::Vector3d& r_ecef, std::string date_utc);

std::tuple<Eigen::Vector3d, Eigen::Vector3d> ecef_to_eci(const Eigen::Vector3d& r_ecef, const Eigen::Vector3d& v_ecef, std::string date_utc);


double kepler_equation_ellipse(double E, double e);

double kepler_rate_ellipse(double E, double e);

double kepler_equation_hyperbola(double H, double e);

double kepler_rate_hyperbola(double H, double e);

double newton_method_kepler_ellipse(double M, double e);

double newton_method_kepler_hyperbola(double M, double e);

double mean_anomaly_from_true_anomaly(double theta, double e);

double semi_major_axis_from_n(double n, double mu);

double specific_angular_momentum(double n, double e, double mu);

double compute_true_anomaly(double t, double M_0, double e, double n);

double r_from_true_anomaly(double theta, double e, double h, double mu);


std::tuple<Eigen::Vector3d, Eigen::Vector3d>  angular_parameters_to_perifocal(double theta, double e, double n, double mu);

Eigen::Matrix3d R_perifocal_to_eci(double i, double Omega, double omega);

std::tuple<Eigen::Vector3d, Eigen::Vector3d> angular_parameters_to_eci(double theta, double e, double n, double i, double Omega, double omega, double mu);
