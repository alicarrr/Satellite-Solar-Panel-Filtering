#pragma once

#include <eigen3/Eigen/Dense>
#include <string>
#include <vector>

typedef enum {
    SOLAR_PANEL,
    METAL
} PanelType_t;


/**
 * What to include:
 * Side Panel struct
 * - NOTE: all solar panel properties will be done in python. inc. absorptivity, emissivity, transmissivity, and reflectivity of the solar panel. These are all wavelength dependent, so we need to know the spectrum of the solar panel.
 * - panel size
 * - normal vector in body frame (of cubesat) of the solar panel
 * - Position in body frame (of cubesat) of the solar panel
 * - Panel type (is it a solar panel or no)
 */

/**
 * Cubesat struct
 * - position of center of mass in body frame of the cubesat
 * - inertia tensor
 * - mass of the cubesat
 * - array of panel structs
 */

/**
 * cubesat state struct:
 * - position in ECI frame of center of mass of the cubesat
 * - rotation matrix from body to ECI frame
 * - velocity of center of mass in ECI frame
 * - angular velocity in ECI frame
 */

// Panel struct
typedef struct {
    double size; // Size of the panel in m^2
    Eigen::Vector3d normal_vector; // Normal vector of the panel in body frame
    Eigen::Vector3d position_rel_GC_in_B; // Position of the panel in body frame
    PanelType_t panel_type=SOLAR_PANEL; // Type of the panel (solar panel or metal). Default is solar panel.
} Panel_t;

// Cubesat struct (made of panels)
typedef struct {
    Eigen::Vector3d position_COM_rel_GC; // Position of the center of mass in body frame
    Eigen::Matrix3d inertia_tensor; // Inertia tensor of the cubesat
    double mass; // Mass of the cubesat in kg
    std::vector<Panel_t> panels; // Array of panel structs
} Cubesat_t;

// Cubesat state struct
typedef struct {
    Eigen::Vector3d position_COM_ECI; // Position of the center of mass in ECI frame
    Eigen::Matrix3d rotation_matrix_ECI_from_B; // Rotation matrix from body frame to ECI frame
    Eigen::Vector3d velocity_COM_ECI; // Velocity of the center of mass in ECI frame
    Eigen::Vector3d angular_velocity_ECI; // Angular velocity in ECI frame
} CubesatState_t;


// Initialise cubesat and panels
inline const Panel_t solar_panel_px = {
    .size = 20e-2 * 10e-2, // 20 cm x 10 cm
    .normal_vector = Eigen::Vector3d(1.0, 0.0, 0.0),
    .position_rel_GC_in_B = Eigen::Vector3d(5e-2, 0.0, 0.0),
    .panel_type = SOLAR_PANEL 
};
inline const Panel_t solar_panel_nx = {
    .size = 20e-2 * 10e-2, // 20 cm x 10 cm
    .normal_vector = Eigen::Vector3d(-1.0, 0.0, 0.0),
    .position_rel_GC_in_B = Eigen::Vector3d(-5e-2, 0.0, 0.0),
    .panel_type = SOLAR_PANEL
};
inline const Panel_t solar_panel_py = {
    .size = 20e-2 * 10e-2, // 20 cm x 10 cm
    .normal_vector = Eigen::Vector3d(0.0, 1.0, 0.0),
    .position_rel_GC_in_B = Eigen::Vector3d(0.0, 5e-2, 0.0),
    .panel_type = SOLAR_PANEL
};
inline const Panel_t solar_panel_ny = {
    .size = 20e-2 * 10e-2, // 20 cm x 10 cm
    .normal_vector = Eigen::Vector3d(0.0, -1.0, 0.0),
    .position_rel_GC_in_B = Eigen::Vector3d(0.0, -5e-2, 0.0),
    .panel_type = SOLAR_PANEL
};
inline const Panel_t panel_pz = {
    .size = 10e-2 * 10e-2, // 10 cm x 10 cm
    .normal_vector = Eigen::Vector3d(0.0, 0.0, 1.0),
    .position_rel_GC_in_B = Eigen::Vector3d(0.0, 0.0, 10e-2),
    .panel_type = METAL
};
inline const Panel_t panel_nz = {
    .size = 10e-2 * 10e-2, // 10 cm x 10 cm
    .normal_vector = Eigen::Vector3d(0.0, 0.0, -1.0),
    .position_rel_GC_in_B = Eigen::Vector3d(0.0, 0.0, -10e-2),
    .panel_type = METAL
};


static const Eigen::Matrix3d J_sat_in_B = Eigen::Matrix3d({{1.058e7, -44930.65, 1.471e5}, {-44930.65, 1.035e7, -2.626e5}, {1.471e5, -2.626e5, 4.722e6}}) * 1e-9; // g mm^2 -> kg m^2


inline const Cubesat_t cubesat = {
    .position_COM_rel_GC = Eigen::Vector3d(0.0, 0.0, 0.0),
    .inertia_tensor = J_sat_in_B,
    .mass = 2.0, // in kg
    .panels = {
        solar_panel_px,
        solar_panel_nx,
        solar_panel_py,
        solar_panel_ny,
        panel_pz,
        panel_nz
    }
};

