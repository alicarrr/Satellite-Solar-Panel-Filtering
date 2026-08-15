#pragma once

#include <eigen3/Eigen/Dense>

// HAMSTER data selection.
// The executable is run from the project root, which is the parent of C/.
const int hamster_day_of_year = 1;
const char hamster_dir[] = "C/lib/HAMSTER";


typedef struct
{
    int file_id;
    int albedo_var_id;

    int n_lat;
    int n_lon;
    int n_wavelength;

    float valid_min;
    float valid_max;

    Eigen::VectorXf latitude;
    Eigen::VectorXf longitude;
    Eigen::VectorXf wavelength;

} AlbedoData;


AlbedoData read_albedo_data();

Eigen::VectorXf get_albedo_spectrum(
    const AlbedoData& albedo_data,
    float latitude_deg,
    float longitude_deg
);

void close_albedo_data(AlbedoData& albedo_data);