#include "transforms.hpp"

#include <eigen3/Eigen/Dense>
#include <tuple>

#include "time_conv.hpp"
#include "constants.hpp"
#include "vectorMaths.hpp"



Eigen::Matrix3d eci_to_ecliptic_matrix()
{
    double epsilon = -i_earth; // Obliquity of the ecliptic in radians

    return Rx(epsilon);
}

Eigen::Vector3d eci_to_ecliptic(const Eigen::Vector3d& r_eci)
{
    Eigen::Matrix3d R = eci_to_ecliptic_matrix();
    return R * r_eci;
}


std::tuple<Eigen::Vector3d, Eigen::Vector3d> eci_to_ecliptic(const Eigen::Vector3d& r_eci, const Eigen::Vector3d& v_eci)
{
    Eigen::Matrix3d R = eci_to_ecliptic_matrix();
    
    return std::make_tuple(R * r_eci, R * v_eci);
}


Eigen::Vector3d ecliptic_to_eci(const Eigen::Vector3d& r_ecliptic)
{
    Eigen::Matrix3d R = eci_to_ecliptic_matrix();
    return R.transpose() * r_ecliptic;
}

std::tuple<Eigen::Vector3d, Eigen::Vector3d> ecliptic_to_eci(const Eigen::Vector3d& r_ecliptic, const Eigen::Vector3d& v_ecliptic)
{
    Eigen::Matrix3d R = eci_to_ecliptic_matrix();
    
    return std::make_tuple(R.transpose() * r_ecliptic, R.transpose() * v_ecliptic);
}


Eigen::Vector3d ecef_to_llh(const Eigen::Vector3d& r_ecef)
{
    double x = r_ecef(0);
    double y = r_ecef(1);
    double z = r_ecef(2);

    double a = a_earth; // Semi-major axis of the Earth
    double b = b_earth; // Semi-minor axis of the Earth
    double e_E = sqrt((a*a - b*b) / (a*a)); // Eccentricity of the Earth
    double e_E_prime = sqrt((a*a - b*b) / (b*b)); // Second eccentricity of the Earth

    double p = sqrt(x*x + y*y);
    double F = 54*b*b * z*z;
    double G = p*p + (1 - e_E*e_E)*z*z - e_E*e_E*(a*a - b*b);
    double c = e_E*e_E*e_E*e_E * F * p*p / (G*G*G);
    double s = cbrt(1 + c + sqrt(c*c + 2*c));
    double k = s + 1 + 1/s;
    double P = F / (3 * k*k * G*G);
    double Q = sqrt(1 + 2 * e_E*e_E*e_E*e_E * P);
    double r0 = (-P * e_E*e_E * p)/(1 + Q) + sqrt(0.5*a*a * (1 + 1/Q) - (P * (1 - e_E*e_E) * z*z)/(Q * (1 + Q)) - 0.5 * P * p*p);
    double U = sqrt((p - e_E*e_E * r0)*(p - e_E*e_E * r0) + z*z);
    double V = sqrt((p - e_E*e_E * r0)*(p - e_E*e_E * r0) + (1 - e_E*e_E)*z*z);
    double z0 = (b*b * z) / (a*V);

    double latitude = atan((z + e_E_prime*e_E_prime * z0) / p);
    double longitude = atan2(y, x);
    double h = U * (1 - b*b/(a*V));

    return Eigen::Vector3d(latitude, longitude, h);
}


Eigen::Vector3d llh_to_ecef(const Eigen::Vector3d& llh)
{
    double latitude = llh(0);
    double longitude = llh(1);
    double h = llh(2);

    double a = a_earth; // Semi-major axis of the Earth
    double b = b_earth; // Semi-minor axis of the Earth

    double N = a*a / sqrt(a*a*cos(latitude)*cos(latitude) + b*b*sin(latitude)*sin(latitude));

    double x = (N + h) * cos(latitude) * cos(longitude);
    double y = (N + h) * cos(latitude) * sin(longitude);
    double z = (N * (b*b/(a*a)) + h) * sin(latitude);

    return Eigen::Vector3d(x, y, z);
}


Eigen::Vector3d eci_to_ecef(const Eigen::Vector3d& r_eci, std::string date_utc)
{
    double GMST = epoch_to_GMST(date_utc);

    Eigen::Matrix3d R = Rz(GMST * angular_vel_earth);

    return R * r_eci;
}


std::tuple<Eigen::Vector3d, Eigen::Vector3d> eci_to_ecef(const Eigen::Vector3d& r_eci, const Eigen::Vector3d& v_eci, std::string date_utc)
{
    double GMST = epoch_to_GMST(date_utc);

    Eigen::Vector3d omega_earth(0, 0, angular_vel_earth);

    Eigen::Matrix3d R = Rz(GMST * angular_vel_earth);
    Eigen::Matrix3d R_dot = cross_product_matrix(omega_earth) * R;

    Eigen::Vector3d r_ecef = R * r_eci;
    Eigen::Vector3d v_ecef = R * v_eci + R_dot * r_eci;

    return std::make_tuple(r_ecef, v_ecef);
}


Eigen::Vector3d ecef_to_eci(const Eigen::Vector3d& r_ecef, std::string date_utc)
{
    double GMST = epoch_to_GMST(date_utc);

    Eigen::Matrix3d R = Rz(GMST * angular_vel_earth);

    return R.transpose() * r_ecef;
}


std::tuple<Eigen::Vector3d, Eigen::Vector3d> ecef_to_eci(const Eigen::Vector3d& r_ecef, const Eigen::Vector3d& v_ecef, std::string date_utc)
{
    double GMST = epoch_to_GMST(date_utc);

    Eigen::Vector3d omega_earth(0, 0, angular_vel_earth);

    Eigen::Matrix3d R = Rz(GMST * angular_vel_earth);
    Eigen::Matrix3d R_dot = cross_product_matrix(omega_earth) * R;

    Eigen::Vector3d r_eci = R.transpose() * r_ecef;
    Eigen::Vector3d v_eci = R.transpose() * (v_ecef - R_dot * r_eci);

    return std::make_tuple(r_eci, v_eci);
}


double kepler_equation_ellipse(double E, double e)
{
    return E - e * sin(E);
}


double kepler_rate_ellipse(double E, double e)
{
    return 1.0 - e * cos(E);
}


double kepler_equation_hyperbola(double H, double e)
{
    return e * sinh(H) - H;
}


double kepler_rate_hyperbola(double H, double e)
{
    return e * cosh(H) - 1.0;
}


double newton_method_kepler_ellipse(double M, double e)
{
    double E = M; // Initial guess
    for (int i = 0; i < 1000; i++) // Max iterations
    {
        double f = kepler_equation_ellipse(E, e) - M;
        double f_prime = kepler_rate_ellipse(E, e);
        double E_next = E - f / f_prime;
        if (fabs(E_next - E) < 1e-12) // Convergence check
            break;
        E = E_next;
    }
    return E;
}


double newton_method_kepler_hyperbola(double M, double e)
{
    double H = log(2*M/e + 1); // Initial guess
    for (int i = 0; i < 1000; i++) // Max iterations
    {
        double f = kepler_equation_hyperbola(H, e) - M;
        double f_prime = kepler_rate_hyperbola(H, e);
        double H_next = H - f / f_prime;
        if (fabs(H_next - H) < 1e-12) // Convergence check
            break;
        H = H_next;
    }
    return H;
}


double mean_anomaly_from_true_anomaly(double theta, double e)
{
    double M;
    if (e < 1.0) // Elliptical orbit
    {
        M = 2*atan(sqrt((1.0-e)/(1.0+e))*tan(theta/2.0)) - (e * sin(theta) * sqrt(1.0 - e*e))/(1.0 + e*cos(theta));
    }
    else if(e == 1) // Parabolic
    {
        M = M = (1.0/2.0)*tan(theta/2.0) + (1.0/6.0)*(pow(tan(theta/2.0), 3));
    }
    else
    {
        // Hyperbolic orbit
        M = (e * sqrt(e*e - 1) * sin(theta))/(1.0 + e*cos(theta)) - log((sqrt(e + 1.0) + sqrt(e - 1.0) * tan(theta/2.0))/(sqrt(e + 1.0) - sqrt(e - 1.0) * tan(theta/2.0)));
    }

    return M;
}


double semi_major_axis_from_n(double n, double mu)
{
    return cbrt(mu / (n*n));
}


double specific_angular_momentum(double n, double e, double mu)
{
    double a = semi_major_axis_from_n(n, mu);
    return sqrt(a * abs(1.0 - e*e) * mu);
}


double compute_true_anomaly(double t, double M_0, double e, double n)
{
    double M = M_0 + n * t;

    double theta;
    if (e < 1.0) // Elliptical orbit
    {
        double E = newton_method_kepler_ellipse(M, e);
        theta = 2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(E/2.0));
    }
    else if(e == 1) // Parabolic
    {
        double z = cbrt(3*M + sqrt(1 + 9*M*M));
        theta = 2*atan(z - 1.0/z);
    }
    else
    {
        double H = newton_method_kepler_hyperbola(M, e);
        theta = 2.0*atan(sqrt((e+1.0)/(e-1.0))*tanh(H/2.0));
    }

    theta = wrapNumber(theta, -M_PI, M_PI);

    return theta;
}


double r_from_true_anomaly(double theta, double e, double h, double mu)
{
    return (h*h/mu) / (1.0 + e*cos(theta));
}


std::tuple<Eigen::Vector3d, Eigen::Vector3d>  angular_parameters_to_perifocal(double theta, double e, double n, double mu)
{
    double h = specific_angular_momentum(n, e, mu);
    double r = r_from_true_anomaly(theta, e, h, mu);

    double p = r * cos(theta);
    double q = r * sin(theta);
    
    double p_dot = -(mu/h) * sin(theta);
    double q_dot = (mu/h) * (e + cos(theta));


    return std::make_tuple(Eigen::Vector3d(p, q, 0), Eigen::Vector3d(p_dot, q_dot, 0));
}


Eigen::Matrix3d R_perifocal_to_eci(double i, double Omega, double omega)
{
    Eigen::Matrix3d R;
    R = Rz(Omega) * Rx(i) * Rz(omega);
    return R;
}


std::tuple<Eigen::Vector3d, Eigen::Vector3d> angular_parameters_to_eci(double theta, double e, double n, double i, double Omega, double omega, double mu)
{
    std::tuple<Eigen::Vector3d, Eigen::Vector3d> perifocal = angular_parameters_to_perifocal(theta, e, n, mu);
    Eigen::Vector3d r_perifocal = std::get<0>(perifocal);
    Eigen::Vector3d v_perifocal = std::get<1>(perifocal);

    Eigen::Matrix3d R = R_perifocal_to_eci(i, Omega, omega);
    return std::make_tuple(R * r_perifocal, R * v_perifocal);
}
