#include <math.h>
#include <stdio.h>
#include <iostream>

// #include <Eigen/Dense>
#include "vectorMaths.hpp"


double wrapNumber(double N, double lower_bound, double upper_bound)
{
    return lower_bound + fmod((N - lower_bound), (upper_bound - lower_bound));
}


void matmul(double A[3][3], double B[3][3], double C[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            C[i][j] = 0.0;
            for (int k = 0; k < 3; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

Eigen::Matrix3d matmul(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B)
{
    return A * B;
}


void matmul(double A[3][3], double B[3], double C[3])
{
    for (int i = 0; i < 3; i++)
    {
        C[i] = 0.0;
        for (int k = 0; k < 3; k++)
        {
            C[i] += A[i][k] * B[k];
        }
    }
}

Eigen::Vector3d matmul(const Eigen::Matrix3d& A, const Eigen::Vector3d& B)
{
    return A * B;
}


void inverse(double A[3][3], double A_inv[3][3])
{
    double det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) -
                 A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) +
                 A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);

    if (det == 0) {
        // Handle singular matrix case (not invertible)
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                A_inv[i][j] = 0.0;
        return;
    }

    double inv_det = 1.0 / det;

    A_inv[0][0] = inv_det * (A[1][1]*A[2][2] - A[1][2]*A[2][1]);
    A_inv[0][1] = inv_det * (A[0][2]*A[2][1] - A[0][1]*A[2][2]);
    A_inv[0][2] = inv_det * (A[0][1]*A[1][2] - A[0][2]*A[1][1]);

    A_inv[1][0] = inv_det * (A[1][2]*A[2][0] - A[1][0]*A[2][2]);
    A_inv[1][1] = inv_det * (A[0][0]*A[2][2] - A[0][2]*A[2][0]);
    A_inv[1][2] = inv_det * (A[0][2]*A[1][0] - A[0][0]*A[1][2]);

    A_inv[2][0] = inv_det * (A[1][0]*A[2][1] - A[1][1]*A[2][0]);
    A_inv[2][1] = inv_det * (A[0][1]*A[2][0] - A[0][0]*A[2][1]);
    A_inv[2][2] = inv_det * (A[0][0]*A[1][1] - A[0][1]*A[1][0]);

}


void Rx(double angle, double R[3][3])
{
    double c = cos(angle);
    double s = sin(angle);

    R[0][0] = 1.0; R[0][1] = 0.0; R[0][2] = 0.0;
    R[1][0] = 0.0; R[1][1] = c;   R[1][2] = -s;
    R[2][0] = 0.0; R[2][1] = s;   R[2][2] = c;
}

void Ry(double angle, double R[3][3])
{
    double c = cos(angle);
    double s = sin(angle);

    R[0][0] = c;   R[0][1] = 0.0; R[0][2] = s;
    R[1][0] = 0.0; R[1][1] = 1.0; R[1][2] = 0.0;
    R[2][0] = -s;  R[2][1] = 0.0; R[2][2] = c;
}

void Rz(double angle, double R[3][3])
{
    double c = cos(angle);
    double s = sin(angle);

    R[0][0] = c;   R[0][1] = -s;  R[0][2] = 0.0;
    R[1][0] = s;   R[1][1] = c;   R[1][2] = 0.0;
    R[2][0] = 0.0; R[2][1] = 0.0; R[2][2] = 1.0;
}



Eigen::Matrix3d Rx(double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    Eigen::Matrix3d R;
    R << 1.0, 0.0, 0.0,
         0.0, c,   -s,
         0.0, s,   c;
    return R;
}


Eigen::Matrix3d Ry(double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    Eigen::Matrix3d R;
    R << c,   0.0, s,
         0.0, 1.0, 0.0,
         -s,  0.0, c;
    return R;
}


Eigen::Matrix3d Rz(double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    Eigen::Matrix3d R;
    R << c,   -s,  0.0,
         s,   c,   0.0,
         0.0, 0.0, 1.0;
    return R;
}



void transpose(double A[3][3], double A_T[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            A_T[j][i] = A[i][j];
        }
    }
}

void cross_product(double a[3], double b[3], double result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}


void cross_product_matrix(double a[3], double A[3][3])
{
    A[0][0] = 0.0;   A[0][1] = -a[2]; A[0][2] = a[1];
    A[1][0] = a[2];  A[1][1] = 0.0;   A[1][2] = -a[0];
    A[2][0] = -a[1]; A[2][1] = a[0];  A[2][2] = 0.0;
}

Eigen::Matrix3d cross_product_matrix(const Eigen::Vector3d& a)
{
    Eigen::Matrix3d A;
    A << 0.0, -a(2), a(1),
         a(2), 0.0, -a(0),
         -a(1), a(0), 0.0;
    return A;
}



void cross_product_matrix_inv(double A[3][3], double a[3])
{
    a[0] = A[2][1];
    a[1] = A[0][2];
    a[2] = A[1][0];
}

Eigen::Vector3d cross_product_matrix_inv(const Eigen::Matrix3d& A)
{
    Eigen::Vector3d a;
    a(0) = A(2, 1);
    a(1) = A(0, 2);
    a(2) = A(1, 0);
    return a;
}


void euler_to_rotation_matrix(double euler_angles[3], double R[3][3])
{
    double precess, nutate, spin;
    spin = euler_angles[0];
    nutate = euler_angles[1];
    precess = euler_angles[2];


    double R_1[3][3], R_2[3][3], R_3[3][3];
    Rz(spin, R_1);
    Ry(nutate, R_2);
    Rz(precess, R_3);

    double R_temp[3][3];
    matmul(R_3, R_2, R_temp);
    matmul(R_temp, R_1, R);
}


Eigen::Matrix3d euler_to_rotation_matrix(const Eigen::Vector3d& euler_angles)
{
    double spin = euler_angles(0);
    double nutate = euler_angles(1);
    double precess = euler_angles(2);

    Eigen::Matrix3d R_1 = Rz(spin);
    Eigen::Matrix3d R_2 = Ry(nutate);
    Eigen::Matrix3d R_3 = Rz(precess);

    return R_3 * R_2 * R_1;
}


void rotation_matrix_to_euler(double R[3][3], double euler_angles[3])
{
    double precess, nutate, spin;
    // double Rt[3][3];
    precess = atan2(R[1][2], R[0][2]);
    nutate = atan2(sqrt(1 - R[2][2]*R[2][2]), R[2][2]);
    spin = atan2(R[2][1], -R[2][0]);

    euler_angles[0] = spin;
    euler_angles[1] = nutate;
    euler_angles[2] = precess;
}


Eigen::Vector3d rotation_matrix_to_euler(const Eigen::Matrix3d& R)
{
    double precess = atan2(R(1, 2), R(0, 2));
    double nutate = atan2(sqrt(1 - R(2, 2)*R(2, 2)), R(2, 2));
    double spin = atan2(R(2, 1), -R(2, 0));

    Eigen::Vector3d euler_angles;
    euler_angles << spin, nutate, precess;
    return euler_angles;
}



void euler_rates_omega_jac(double euler_angles[3], double J[3][3])
{
    double phi = euler_angles[1];
    double psi = euler_angles[2];

    double c_phi = cos(phi);
    double s_phi = sin(phi);
    double c_psi= cos(psi);
    double s_psi = sin(psi);

    J[0][0] = c_psi*s_phi; J[0][1] = -s_psi; J[0][2] = 0.0;
    J[1][0] = s_psi*s_phi; J[1][1] = c_psi;  J[1][2] = 0.0;
    J[2][0] = c_phi;      J[2][1] = 0.0;    J[2][2] = 1.0;
}

Eigen::Matrix3d euler_rates_omega_jac(const Eigen::Vector3d& euler_angles)
{
    double phi = euler_angles(1);
    double psi = euler_angles(2);

    double c_phi = cos(phi);
    double s_phi = sin(phi);
    double c_psi= cos(psi);
    double s_psi = sin(psi);

    Eigen::Matrix3d J;
    J << c_psi*s_phi, -s_psi, 0.0,
         s_psi*s_phi, c_psi, 0.0,
         c_phi,      0.0,   1.0;
    return J;
}


void euler_rates_omega_jac_dot(double euler_angles[3], double euler_rates[3], double J_dot[3][3])
{
    double phi = euler_angles[1];
    double psi = euler_angles[2];

    double phi_dot = euler_rates[1];
    double psi_dot = euler_rates[2];

    double c_phi = cos(phi);
    double s_phi = sin(phi);
    double c_psi= cos(psi);
    double s_psi = sin(psi);

    J_dot[0][0] = -psi_dot*s_psi*s_phi + phi_dot * c_psi*c_phi;
    J_dot[0][1] = -psi_dot*c_psi;
    J_dot[0][2] = 0.0;

    J_dot[1][0] = psi_dot*c_psi*s_phi + phi_dot * s_psi*c_phi;
    J_dot[1][1] = -psi_dot*s_psi;
    J_dot[1][2] = 0.0;

    J_dot[2][0] = -phi_dot * s_phi;
    J_dot[2][1] = 0.0;
    J_dot[2][2] = 0.0;
}


Eigen::Matrix3d euler_rates_omega_jac_dot(const Eigen::Vector3d& euler_angles, const Eigen::Vector3d& euler_rates)
{
    double phi = euler_angles(1);
    double psi = euler_angles(2);

    double phi_dot = euler_rates(1);
    double psi_dot = euler_rates(2);

    double c_phi = cos(phi);
    double s_phi = sin(phi);
    double c_psi= cos(psi);
    double s_psi = sin(psi);

    Eigen::Matrix3d J_dot;
    J_dot << -psi_dot*s_psi*s_phi + phi_dot * c_psi*c_phi, -psi_dot*c_psi, 0.0,
              psi_dot*c_psi*s_phi + phi_dot * s_psi*c_phi, -psi_dot*s_psi, 0.0,
              -phi_dot * s_phi, 0.0, 0.0;
    return J_dot;
}


void euler_accels_from_alpha(double euler_angles[3], double euler_rates[3], double alpha_B_in_G[3], double euler_accels[3])
{
    double J_omega[3][3];
    double J_omega_inv[3][3];
    double J_omega_dot[3][3];

    euler_rates_omega_jac(euler_angles, J_omega);
    euler_rates_omega_jac_dot(euler_angles, euler_rates, J_omega_dot);


    inverse(J_omega, J_omega_inv);

    

    double Cor_vec[3];
    matmul(J_omega_dot, euler_rates, Cor_vec);

    double rhs[3];
    for (int i = 0; i < 3; i++)
    {
        rhs[i] = alpha_B_in_G[i] - Cor_vec[i];
    }

    matmul(J_omega_inv, rhs, euler_accels);

}


Eigen::Vector3d euler_accels_from_alpha(const Eigen::Vector3d& euler_angles, const Eigen::Vector3d& euler_rates, const Eigen::Vector3d& alpha_B_in_G)
{
    Eigen::Matrix3d J_omega = euler_rates_omega_jac(euler_angles);
    Eigen::Matrix3d J_omega_dot = euler_rates_omega_jac_dot(euler_angles, euler_rates);
    Eigen::Matrix3d J_omega_inv = J_omega.inverse();

    Eigen::Vector3d Cor_vec = J_omega_dot * euler_rates;
    Eigen::Vector3d rhs = alpha_B_in_G - Cor_vec;

    return J_omega_inv * rhs;
}


void unit_vector(double a[3], double result[3])
{
    double norm = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    if (norm > 0) {
        result[0] = a[0] / norm;
        result[1] = a[1] / norm;
        result[2] = a[2] / norm;
    } else {
        result[0] = 0.0;
        result[1] = 0.0;
        result[2] = 0.0;
    }
}


Eigen::Vector3d unit_vector(const Eigen::Vector3d& a)
{
    double norm = a.norm();
    if (norm > 0) {
        return a / norm;
    } else {
        return Eigen::Vector3d::Zero();
    }
}


void vecmul(double a[3], double b, double result[3])
{
    result[0] = a[0] * b;
    result[1] = a[1] * b;
    result[2] = a[2] * b;
}


std::vector<double> to_std_vector(const Eigen::VectorXd& x) {
    return std::vector<double>(x.data(), x.data() + x.size());
}

