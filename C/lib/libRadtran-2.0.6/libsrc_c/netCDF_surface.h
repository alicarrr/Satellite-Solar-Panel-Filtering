/************************************************************************
 * $Id$
 *
 * MYSTIC - Monte Carlo code for the physically correct tracing of
 *          photons in cloudy atmospheres.
 *
 * Copyright (c) 2000-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * Correspondence: bernhard.mayer@lmu.de
 *
 ************************************************************************/

/** @file 
 * The files `netCDF_surface.h`, `netCDF_surface.c` and parts of `mystic.c` contain
 * code related to loading surface data from netCDF files.
 * @author Cedric Wehrum
 * @date Nov 1 2023
 */

#ifndef __ncdf_surface_h
#define __ncdf_surface_h

#if defined(__cplusplus)
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "errors.h"
#include "mystic.h"

// Error codes associated with netCDF_surface.c
#define NCDF_NOERR 0
#define NCDF_ERROR -1
#define NCDF_NETCDF_ERROR -2
#define NCDF_UNSUPPORTED_MAP_TYPE -3
#define NCDF_TYPE_NOT_ACCEPTABLE -4
#define NCDF_UNSUPPORTED_VAR_TYPE -5
#define NCDF_ERROR_STANDARD_NAME -6
#define NCDF_ERROR_AXIS -7
#define NCDF_INVALID_DIMENSION -8
#define NCDF_INVALID_PROJECTION -9
#define NCDF_UNSUPPORTED_DIM_TYPE -10
#define NCDF_UNKNOWN_UNIT -11
#define NCDF_FLAG_ERROR -12
#define NCDF_PACK_ERROR -13
#define NCDF_INVALID_MODEL_PARAM -14
#define NCDF_INTERPOLATION_ERROR -15
#define NCDF_PRINT_ERROR -16
#define NCDF_MODELS_DOUBLE_DEFINITION -17
#define NCDF_PARAMETER_DOUBLE_INIT -18
#define NCDF_NO_NETCDF4 -19
#define NCDF_INTERNAL_ERROR -20
#define NCDF_MISSING_PARAMETER -21
#define NCDF_LAST_ERROR -22
typedef int ncdf_status_t;

/* // Type of projections, currently unused (TODO) */
/* #define MC_PROJ_UNDEFINED  0 */
/* #define MC_PROJ_XY         1 */
/* #define MC_PROJ_LON_LAT    2 */
/* typedef unsigned char ncdf_projection_t; */

/// surface data is saved contiguously in a flat one-dimensional array since
/// NetCDF functions return the same format and it can be more memory efficient.
/// To comfortably access this array this struct keeps track of the
/// offsets (strides) from one column to another
/* typedef struct { */
/*   size_t param;     ///< in case an array contains multiple parameters */
/*   size_t wl;        ///< wavelength */
/*   size_t y;         ///< y coordinate or latitude */
/*   size_t x;         ///< x coordinate or longitude */
/*   /// Some functions also need the total lenght of the array */
/*   size_t length; */
/* } mc_ncdf_array_stride; */

/* /// An array that holds surface data of several parameters */
/* typedef struct { */
/*   /// The offsets between different columns */
/*   mc_ncdf_array_stride stride; */
/*   /// The actual data */
/*   /// can have up to 4 dimensions in any order, the offsets in the stride member */
/*   /// indirectly define the actual order */
/*   float* data; */
/*   /// Given one of the parameter constants (MCALB_<model>_PARAM_<param-name>) */
/*   /// in `mystic_definitions.h` this identifies the corresponding column */
/*   size_t param_indices[MCALB_NUM_PARAMS]; */
/*   /// Given an albedo method from `mystic_definitions.h` (MCALB_<model>) this */
/*   /// identifies how many parameters are saved for each model */
/*   size_t param_count[MCALB_NUM_METHODS]; */
/*   /// The actual count of parameter columns (equivalent to max(param_count)) */
/*   size_t max_param_count; */
/* } mc_ncdf_array; */

/* /// This struct contains all relevant information regarding a surface */
/* typedef struct { */
/*   /// array of all longitudes */
/*   double* x;          //TODO check: equidistant and increasing */
/*   /// array of all latitudes */
/*   double* y;          //TODO check: equidistant and increasing */
/*   /// array of all wavelengths */
/*   double* spectrum;   //TODO check: equidistant and increasing */
/*   double* int_spectrum;    ///< array of all interpolated wavelengths on an */
/*                            ///< internal wavelength grid */
/*   size_t nx;               ///< number of longitudes */
/*   size_t ny;               ///< number of latitudes */
/*   size_t nspectrum;        ///< number of wavelengths */
/*   size_t int_nspectrum;    ///< number of wavelengths on the internal grid */
/*   ncdf_projection_t projection; ///< type of projection (currently unused) */

/*   mc_surface_method* models;                  ///< 2D array that associates */
/*   mc_ncdf_array_stride models_stride;             ///< a model for each pixel */

/*   mc_ncdf_array parameters;                       ///< constant parameters */
/*   mc_ncdf_array spatial_parameters;               ///< 2D spatial parameters */

/*   mc_ncdf_array spectral_parameters;              ///< 1D spectral parameters */
/*   mc_ncdf_array spectral_spatial_parameters;      ///< 3D spatial/spectral */
/*                                                   ///< parameters */

/*   mc_ncdf_array int_spectral_parameters;          ///< same, but interpolated */
/*   mc_ncdf_array int_spectral_spatial_parameters;  ///< to the internal grid */

/*   mc_ncdf_array* param_link[MCALB_NUM_PARAMS];    ///< given an albedo method */
/*   mc_ncdf_array* int_param_link[MCALB_NUM_PARAMS];///< this links to the correct */
/*                                                   ///< array */
/* } mc_ncdf_surface; */

/**
 * @brief Loads a NedtCDF surface file
 * @author Cedric Wehrum
 * @date Nov 1 2023
 *
 * @param[in]  filenames   NetCDF file names
 * @param[in]  n_filenames Number of files
 * @param[out] surface     The surface struct that will be initialized and
 *                         loaded
 * @param[in]  quiet       boolean value that specifies if debug level messages
 *                         shall be printed to `stderr`
 * @return                 NCDF_NOERR (i.e. 0) on succes, error code on failure
 */
ncdf_status_t mc_ncdf_load_n (char** filenames, size_t n_filenames, mc_ncdf_surface* surface, int quiet);

/**
 * @brief Loads a NedtCDF surface file
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * Exactly like `mc_ncdf_load_n`, but `filenames` is null terminated. I.e. there
 * is a trailing null pointer that terminates the array. Also directories are
 * ignored.
 *
 * @param[in]  filenames   NetCDF file names
 * @param[out] surface     The surface struct that will be initialized and
 *                         loaded
 * @param[in]  quiet       boolean value that specifies if debug level messages
 *                         shall be printed to `stderr`
 * @return                 NCDF_NOERR (i.e. 0) on succes, error code on failure
 */
ncdf_status_t mc_ncdf_load (char** filenames, mc_ncdf_surface* data, int quiet);

/**
 * @brief Interpolates values along a wavelength grid
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * This function creates an evenly spaced wavelength grid and interpolates
 * parameters that were loaded in `mc_ncdf_load`
 *
 * @param[in,out] surface The surface whose parameters shall be interpolated
 * @param[in]     wl_min  Minimum wavelength of the grid
 * @param[in]     wl_max  Maximum wavelength of the grid
 * @param[in]     count   Number of wavelengths on the grid
 * @param[in]     quiet   boolean value that specifies if debug level messages
 *                        shall be printed to `stderr`
 * @return                NCDF_NOERR (i.e. 0) on succes, error code on failure
 */
ncdf_status_t mc_ncdf_int_wl (mc_ncdf_surface* surface, double wl_min, double wl_max, size_t count, int quiet);

/**
 * @brief Interpolates values along a wavelength grid
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * This function uses an given wavelength grid and interpolates
 * parameters that were loaded in `mc_ncdf_load`
 *
 * @param[in,out] surface  The surface whose parameters shall be interpolated
 * @param[in]     spectrum Array of wavelengths
 * @param[in]     count    Number of wavelengths on the grid
 * @param[in]     quiet    Boolean value that specifies if debug level messages
 *                         shall be printed to `stderr`
 * @return                 NCDF_NOERR (i.e. 0) on succes, error code on failure
 */
ncdf_status_t mc_ncdf_int_wl_fromspec (mc_ncdf_surface* surface, float* spectrum, size_t count, int quiet);

/**
 * @brief Frees all resources in a surface struct
 * @author Cedric Wehrum
 * @date Nov 1 2023
 *
 * @param[in,out] surface  The surface to be freed
 */
void mc_ncdf_free (mc_ncdf_surface* surface);

/**
 * @brief Retrieves the albedo method for a specific pixel
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * Given a pixel coordinate this will return the requested albedo method.
 * The coordinates `i` and `j` refer to the members `x` and `y` of the given
 * surface.
 *
 * @param[in] surface The surface whose model will be returned
 * @param[in] i       Longitude or x coordinate index
 * @param[in] j       Latitude or y coordinate index
 * @return            One of the albedo methods specified in
 *                    `mystic_definitions.h`
 */
mc_surface_method mc_ncdf_get_model (mc_ncdf_surface* surface, size_t i, size_t j);

/**
 * @brief Retrieves a parameter for a pixel at a specific wavelength
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * Given a pixel coordinate and a wavelength this will return the requested
 * parameter from the uninterpolated values. That means that the wavelength is
 * an index along the spectral axis as it was defined in the NetCDF file.
 * For the interpolated values see `mc_ncdf_get_param_int`. The coordinates
 * `i`, `j` and `wl` refer to the members `x`, `y` and `spectrum` of the given
 * surface.
 *
 * @param[in] surface The surface whose parameter will be interpolated
 * @param[in] i       Longitude or x coordinate index
 * @param[in] j       Latitude or y coordinate index
 * @param[in] wl      Wavelength index
 * @param[in] param   Requested parameter; possible values for
 *                    `mc_surface_parameter` are defined in
 *                    `mystic_definitions.h`
 * @return            Value of the requested parameter
 */
float mc_ncdf_get_param (mc_ncdf_surface* surface, size_t i, size_t j, size_t wl, mc_surface_parameter param);

/**
 * @brief Retrieves a parameter for a pixel at a specific wavelength
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * Given a pixel coordinate and a wavelength this will return the requested
 * parameter from the interpolated values. That means that the wavelength is
 * an index along the spectral axis of the internal wavelength grid.
 * For the uninterpolated values see `mc_ncdf_get_param`. The coordinates
 * `i`, `j` and `wl` refer to the members `x`, `y` and `int_spectrum` of the
 * given surface.
 *
 * @param[in] surface The surface whose parameter will be returned
 * @param[in] i       Longitude or x coordinate index
 * @param[in] j       Latitude or y coordinate index
 * @param[in] wl      Wavelength index
 * @param[in] param   Requested parameter; possible values for
 *                    `mc_surface_parameter` are defined in
 *                    `mystic_definitions.h`
 * @return            Value of the requested parameter
 */
float mc_ncdf_get_param_int (mc_ncdf_surface* surface, size_t i, size_t j, long wl, mc_surface_parameter param);

/**
 * @brief Prints model and all parameters at a specific pixel and wavelength
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * This is a convenience and debugging function that will print the albedo
 * method and all its uninterpolated parameters to `stderr`. That means that the
 * wavelength is an index along the spectral axis as it was defined in the
 * NetCDF file. For the interpolated values see `mc_ncdf_print_params_int`. The
 * coordinates `i`, `j` and `wl` refer to the members `x`, `y` and `spectrum`
 * of the given surface.
 *
 * @param[in] surface The surface whose parameters will be printed
 * @param[in] i       Longitude or x coordinate index
 * @param[in] j       Latitude or y coordinate index
 * @param[in] wl      Wavelength index
 * @return            NCDF_NOERR (i.e. 0) on succes, error code on failure
 */
ncdf_status_t mc_ncdf_print_params (mc_ncdf_surface* surface, size_t i, size_t j, size_t wl);

/**
 * @brief Prints model and all parameters at a specific pixel and wavelength
 * @author Cedric Wehrum
 * @date Nov 1 2023
 * This is a convenience and debugging function that will print the albedo
 * method and all its interpolated parameters to `stderr`. That means that the
 * wavelength is an index along the spectral axis of the internal wavelength
 * grid. For the uninterpolated values see `mc_ncdf_print_params`. The
 * coordinates `i`, `j` and `wl` refer to the members `x`, `y` and
 * `int_spectrum` of the given surface.
 *
 * @param[in] surface The surface whose parameters will be printed
 * @param[in] i       Longitude or x coordinate index
 * @param[in] j       Latitude or y coordinate index
 * @param[in] wl      Wavelength index
 * @return            NCDF_NOERR (i.e. 0) on succes, error code on failure
 */
ncdf_status_t mc_ncdf_print_params_int (mc_ncdf_surface* surface, size_t i, size_t j, long wl);

int reflection_matrix_mixed_surface_properties (albedo_struct* albedo,
                                                photon_struct* p,
                                                int            nstokes,
                                                int            backward,
                                                double         mu_inc,
                                                double         phi_inc,
                                                double         mu_esc,
                                                double         phi_esc,
                                                int            ia,
                                                int            ja,
                                                int            il,
                                                int            iv,
                                                /* Output */
                                                double*** refl_mat);

#if defined(__cplusplus)
}
#endif

#endif  //__ncdf_surface_h
