#pragma once


#include <eigen3/Eigen/Dense>


double wrapNumber(double N, double lower_bound, double upper_bound);

void matmul(double A[3][3], double B[3][3], double C[3][3]);

Eigen::Matrix3d matmul(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B);

void matmul(double A[3][3], double B[3], double C[3]);

Eigen::Vector3d matmul(const Eigen::Matrix3d& A, const Eigen::Vector3d& B);

void Rx(double angle, double R[3][3]);

void Ry(double angle, double R[3][3]);

void Rz(double angle, double R[3][3]);

Eigen::Matrix3d Rx(double angle);

Eigen::Matrix3d Ry(double angle);

Eigen::Matrix3d Rz(double angle);


void transpose(double A[3][3], double A_T[3][3]);

void cross_product(double a[3], double b[3], double result[3]);

void cross_product_matrix(double a[3], double A[3][3]);

Eigen::Matrix3d cross_product_matrix(const Eigen::Vector3d& a);

void cross_product_matrix_inv(double A[3][3], double a[3]);

Eigen::Vector3d cross_product_matrix_inv(const Eigen::Matrix3d& A);

void euler_to_rotation_matrix(double euler_angles[3], double R[3][3]);

Eigen::Matrix3d euler_to_rotation_matrix(const Eigen::Vector3d& euler_angles);

void rotation_matrix_to_euler(double R[3][3], double euler_angles[3]);

Eigen::Vector3d rotation_matrix_to_euler(const Eigen::Matrix3d& R);

void euler_rates_omega_jac(double euler_angles[3], double J[3][3]);

Eigen::Matrix3d euler_rates_omega_jac(const Eigen::Vector3d& euler_angles);

void euler_rates_omega_jac_dot(double euler_angles[3], double euler_rates[3], double J_dot[3][3]);

Eigen::Matrix3d euler_rates_omega_jac_dot(const Eigen::Vector3d& euler_angles, const Eigen::Vector3d& euler_rates);

void euler_accels_from_alpha(double euler_angles[3], double euler_rates[3], double alpha_B_in_G[3], double euler_accels[3]);

Eigen::Vector3d euler_accels_from_alpha(const Eigen::Vector3d& euler_angles, const Eigen::Vector3d& euler_rates, const Eigen::Vector3d& alpha_B_in_G);

void unit_vector(double a[3], double result[3]);

Eigen::Vector3d unit_vector(const Eigen::Vector3d& a);

void vecmul(double a[3], double b, double result[3]);

void inverse(double A[3][3], double A_inv[3][3]);