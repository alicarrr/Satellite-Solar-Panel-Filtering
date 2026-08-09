/*-------------------------------------------------------------------
 * This file is part of libRadtran.
 * Copyright (c) 1997-2018 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras,
 *			   
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
//#include <f2c_pprts.h>
#include "common_math.h"
#include "netCDF_functions.h"
#include "locate.h"
#include "equation.h"
#include "twostrebe.h"
#include "dynamic_tenstream.h"
#include "solver.h"
#include "ascii.h"
//#include "f2c_solver_ids.h"
#include "errors.h"
#include "allocnd.h"

/* --------------------------------------------------------------------

Look-Up Table related functions

----------------------------------------------------------------------*/

/* struct that contains all LUT relevant parameters */
typedef struct LUT {

  // LUT type
  char* type;

  // LUT dimension sizes
  int N_tauz;
  int N_w0;
  int N_aspect_zx;
  int N_g;
  int N_phi;
  int N_theta;

  // LUT sampling points
  float* tauz;
  float* w0;
  float* aspect_zx;
  float* g;
  float* phi;
  float* theta;

  // LUT
  float* LUT;

  // LUT offsets
  size_t* offsets;

} lut;

/*
Function that determines an offset vector that can be used to transform an
n-dimensional index vector to an 1D index value that is needed to address
the entries of flattened arrays such as the look-up tables. 

Input:
    - N_coeff:
        number of coefficients in the LUT that are stored for every individual 
        index value (e.g. 9 for the LUT of direct to direct coefficients)
    - LUT_dims:
        number of look-up table dimensions (4 for diffuse, 6 for direct LUTs)
    - LUT_dimsizes:
        array containing the number of sampling points for every dimension
        of the LUT (these numbers are stored in an initialized LUT struct)

Output:
    - LUT_offsets:
        array containing offsets so that the dot product of LUT_offsets and
        an index array of size LUT_dims returns the 1D index needed to address
        values in the LUTs.
*/
int calculate_LUT_offsets (size_t N_coeff, int LUT_dims, int* LUT_dimsizes, size_t* LUT_offsets) {

  LUT_offsets[0] = N_coeff;
  for (int i = 1; i < LUT_dims; i++) {
    LUT_offsets[i] = LUT_offsets[i - 1] * LUT_dimsizes[i - 1];
  }

  return 0;
}

/*
Function that opens a netCDF file containing a TenStream look-up table to initialize 
all the variables in the struct LUT_data

Input:
    - type:
        type of the look-up table given in filename. Can be either "dir2dir", 
        "dir2diff" or "diff2diff"
    - file:
        path to the file containing the look-up table whose data we want to get 
    - N_direct:
        number of direct streams
    - N_diffuse:
        number of diffuse streams

Output:
    - LUT_data: 
        initialized struct that contains the LUT, its offsets, as well as the sampling points 
        and number of sampling points for every dimension of the LUT
*/
int get_LUT_data (char* type, char* file, int N_direct, int N_diffuse, int quiet, lut* LUT_data) {

  if (strcmp (type, "dir2dir") == 0 || strcmp (type, "dir2diff") == 0 || strcmp (type, "diff2diff") == 0) {

#define MAXDIM 7
    if (!quiet) {
      fprintf (stderr, "FILEPATH: %s \n", file);
    }

    // Load the netCDF file
    int ncid, ierr;
    ierr = open_netcdf (file, &ncid);
    if (ierr != 0) {
      fprintf (stderr,
               "*******************************************************************************\n"
               "Error! Did not find TenStream lookup tables in data/tenstream_lookup.\n"
               "Make sure that you have them stored there in order to run dynamic_tenstream.\n"
               "*******************************************************************************\n");
    }
    CHKERR (ierr);

    int    N_dims;
    size_t dimsizes[MAXDIM] = {
      0,
    };

    // Get the sampling points and number of sampling points for every dimension of the LUT
    ierr = read_netcdf_float (ncid, "pspace.tau.values", &(LUT_data->tauz), &N_dims, dimsizes);
    CHKERR (ierr);
    LUT_data->N_tauz = dimsizes[0];

    //LUT_data->w0 = NULL; //if LUT_data is initialized with malloc instead of calloc
    ierr = read_netcdf_float (ncid, "pspace.w0.values", &(LUT_data->w0), &N_dims, dimsizes);
    CHKERR (ierr);
    LUT_data->N_w0 = dimsizes[0];

    ierr = read_netcdf_float (ncid, "pspace.aspect_zx.values", &(LUT_data->aspect_zx), &N_dims, dimsizes);
    CHKERR (ierr);
    LUT_data->N_aspect_zx = dimsizes[0];

    ierr = read_netcdf_float (ncid, "pspace.g.values", &(LUT_data->g), &N_dims, dimsizes);
    CHKERR (ierr);
    LUT_data->N_g = dimsizes[0];

    if (strcmp (type, "dir2dir") == 0 || strcmp (type, "dir2diff") == 0) {
      ierr = read_netcdf_float (ncid, "pspace.phi.values", &(LUT_data->phi), &N_dims, dimsizes);
      CHKERR (ierr);
      LUT_data->N_phi = dimsizes[0];

      ierr = read_netcdf_float (ncid, "pspace.theta.values", &(LUT_data->theta), &N_dims, dimsizes);
      CHKERR (ierr);
      LUT_data->N_theta = dimsizes[0];
    } else {
      LUT_data->N_phi   = 0;
      LUT_data->N_theta = 0;
    }

    /*
    Depending on the LUT type, get the LUT data and store them in the variable LUT and
    calculate its offsets
    */
    if (strcmp (type, "diff2diff") == 0) {
      LUT_data->type = "diff2diff";

      ierr = read_netcdf_float (ncid, "diffuse.S", &(LUT_data->LUT), &N_dims, dimsizes);
      CHKERR (ierr);

      const int N_diff_dim           = 4;
      int       diff_LUT_dimsizes[4] = {LUT_data->N_tauz, LUT_data->N_w0, LUT_data->N_aspect_zx, LUT_data->N_g};
      size_t*   a_LUT_offsets        = calloc_size_t_1D (N_diff_dim, "a_LUT_offsets");
      ierr                           = calculate_LUT_offsets (N_diffuse * N_diffuse, N_diff_dim, diff_LUT_dimsizes, a_LUT_offsets);
      CHKERR (ierr);
      LUT_data->offsets = a_LUT_offsets;
    } else if ((strcmp (type, "dir2diff") == 0) || (strcmp (type, "dir2dir") == 0)) {
      const int N_dir_dim = 6;
      int       dir_LUT_dimsizes[6] =
        {LUT_data->N_tauz, LUT_data->N_w0, LUT_data->N_aspect_zx, LUT_data->N_g, LUT_data->N_phi, LUT_data->N_theta};

      if (strcmp (type, "dir2diff") == 0) {
        LUT_data->type = "dir2diff";

        ierr = read_netcdf_float (ncid, "direct.S", &(LUT_data->LUT), &N_dims, dimsizes);
        CHKERR (ierr);

        size_t* b_LUT_offsets = calloc_size_t_1D (N_dir_dim, "b_LUT_offsets");
        ierr                  = calculate_LUT_offsets (N_diffuse * N_direct, N_dir_dim, dir_LUT_dimsizes, b_LUT_offsets);
        CHKERR (ierr);
        LUT_data->offsets = b_LUT_offsets;
      } else {
        LUT_data->type = "dir2dir";

        ierr = read_netcdf_float (ncid, "direct.T", &(LUT_data->LUT), &N_dims, dimsizes);
        CHKERR (ierr);

        size_t* c_LUT_offsets = calloc_size_t_1D (N_dir_dim, "c_LUT_offsets");
        ierr                  = calculate_LUT_offsets (N_direct * N_direct, N_dir_dim, dir_LUT_dimsizes, c_LUT_offsets);
        CHKERR (ierr);
        LUT_data->offsets = c_LUT_offsets;
      }
    }

    ierr = close_netcdf (ncid);
    CHKERR (ierr);

    return 0;

  } else {
    fprintf (stderr, "Invalid look-up table type (%s). Must be either 'dir2dir', 'dir2diff' or 'diff2diff'. \n", type);
    return 1;
  }
}

/* 
Function that returns the floating point index of a value within an 1D array.
The floating point index is set to 0 if the value is smaller than the smallest 
value in the array and to N_sampling_points-1 in case it is larger than the 
largest value in the array.

Input:
    - value:
        floating point number whose floating point index within the 1D array 
        sampling_points we want to determine
    - sampling_points:
        1D array for which we want to determine the floating point index
    - N_sampling_points:
        number of elements in the array sampling_points

Output:
    - floating_point_index:
        (pointer to the) floating point index of value in sampling_points
*/
int get_floating_point_index (float value, float* sampling_points, size_t N_sampling_points, float* floating_point_index) {

  if (value < sampling_points[0]) {
    *floating_point_index = 0;
    return 1;
  } else if (value > sampling_points[N_sampling_points - 1]) {
    *floating_point_index = N_sampling_points - 1;
    return 2;
  } else {
    size_t lower_int_index = flocate (sampling_points, N_sampling_points, value);
    *floating_point_index  = lower_int_index + (value - sampling_points[lower_int_index]) /
                                                (sampling_points[lower_int_index + 1] - sampling_points[lower_int_index]);
    return 0;
  }
}

/* Function determining solar angle dependent quantities
    - l_switch_east and l_switch_east:  values that determine in which quadrant of a 2D coordinate
                                        system you are to account for the fact that TenStream LUT's
                                        are only calculated for values between 0 and 90 degrees azimuth
    - phi_deg and theta_deg:            converted solar azimuth and zenith angles with a value between
                                        0 and PI/2 needed for calling the TenStream LUT's
 */
int determine_angular_LUT_quantities (double mu0,
                                      double phi0,
                                      float* theta_deg,
                                      float* phi_deg,
                                      int*   lswitch_east,
                                      int*   lswitch_north) {

  if (phi0 > 90.0 && phi0 <= 270.0) {
    /* y-direction of solar incident angle is switched between 90 and 270 degrees solar azimuth angle */
    *lswitch_north = 1;
  }
  if (phi0 > 180.0) {
    /* x-direction of solar incident angle is switched for solar azimuth angles larger than 180 degrees */
    *lswitch_east = 1;
  }

  /* convert the solar azimuth angle to a value between 0 and PI/2 as it is done in
     external/gitsrc_tenstream/src/pprts.F90 and convert it to degrees afterwards */
  *phi_deg = asin (sin (acos (cos (phi0 * M_PI / 180.0)))) * 180.0 / M_PI;

  /* convert solar zenith angle to degrees for TenStream routines */
  *theta_deg = acos (mu0) * 180.0 / M_PI;

  return 0;
}

/*
Function that determines the floating point index for the coefficients in a TenStream
LUT for a given set of optical properties. This index is required to interpolate the 
coefficients stored in the LUTs to the set of optical properties handed to this function.

Input:
    - tauz:
        optical thickness of the grid box
    - w0:
        single scattering albedo of the grid box
    - aspect_zx:
        aspect ratio of the grid box
    - g:
        asymmetry parameter of the grid box
    - phi_deg:
        solar azimuth angle
    - theta_deg:
        solar zenith angle
    - LUT_data:
        struct containing all LUT related data (has to be initialized with the
        "get_LUT_data" function first!)

Output:
    - interpolation_index:
        4- or 6-dimensional index specifying the coefficients of the grid box 
        with the given set of optical properties in the TenStream LUT
*/
int get_interpolation_index (float  tauz,
                             float  w0,
                             float  aspect_zx,
                             float  g,
                             float  phi_deg,
                             float  theta_deg,
                             lut*   LUT_data,
                             float* interpolation_index) {

  int ierr;

  ierr = get_floating_point_index (tauz, LUT_data->tauz, LUT_data->N_tauz, &interpolation_index[0]);
  if (ierr != 0) {
    fprintf (
      stderr,
      "CAUTION! The value of tau = %g is outside the range of the LUT (MIN: %g, MAX: %g). It was thus set to tau_dir = %g.\n",
      tauz,
      LUT_data->tauz[0],
      LUT_data->tauz[LUT_data->N_tauz - 1],
      LUT_data->tauz[(int) interpolation_index[0]]);
  }

  ierr = get_floating_point_index (w0, LUT_data->w0, LUT_data->N_w0, &interpolation_index[1]);
  if (ierr != 0) {
    fprintf (stderr,
             "CAUTION! The value of w0 = %g is outside the range of the LUT (MIN: %g, MAX: %g). It was thus set to w0 = %g.\n",
             w0,
             LUT_data->w0[0],
             LUT_data->w0[LUT_data->N_w0 - 1],
             LUT_data->w0[(int) interpolation_index[1]]);
  }

  ierr = get_floating_point_index (aspect_zx, LUT_data->aspect_zx, LUT_data->N_aspect_zx, &interpolation_index[2]);
  if (ierr != 0) {
    fprintf (stderr,
             "CAUTION! The value of aspect_zx = %g is outside the range of the LUT (MIN: %g, MAX: %g). It was thus set to "
             "aspect_zx = %g.\n",
             aspect_zx,
             LUT_data->aspect_zx[0],
             LUT_data->aspect_zx[LUT_data->N_aspect_zx - 1],
             LUT_data->aspect_zx[(int) interpolation_index[2]]);
  }

  ierr = get_floating_point_index (g, LUT_data->g, LUT_data->N_g, &interpolation_index[3]);
  if (ierr != 0) {
    fprintf (stderr,
             "CAUTION! The value of g = %g is outside the range of the LUT (MIN: %g, MAX: %g). It was thus set to g = %g.\n",
             g,
             LUT_data->g[0],
             LUT_data->g[LUT_data->N_g - 1],
             LUT_data->g[(int) interpolation_index[3]]);
  }

  if (strcmp (LUT_data->type, "dir2dir") == 0 || strcmp (LUT_data->type, "dir2diff") == 0) {
    ierr = get_floating_point_index (phi_deg, LUT_data->phi, LUT_data->N_phi, &interpolation_index[4]);
    if (ierr != 0) {
      fprintf (
        stderr,
        "CAUTION! The value of phi_deg = %g is outside the range of the LUT (MIN: %g, MAX: %g). It was thus set to phi_deg = %g.\n",
        phi_deg,
        LUT_data->phi[0],
        LUT_data->phi[LUT_data->N_phi - 1],
        LUT_data->phi[(int) interpolation_index[4]]);
    }

    ierr = get_floating_point_index (theta_deg, LUT_data->theta, LUT_data->N_theta, &interpolation_index[5]);
    if (ierr != 0) {
      fprintf (stderr,
               "CAUTION! The value of theta_deg = %g is outside the range of the LUT (MIN: %g, MAX: %g). It was thus set to "
               "theta_deg = %g.\n",
               theta_deg,
               LUT_data->theta[0],
               LUT_data->theta[LUT_data->N_theta - 1],
               LUT_data->theta[(int) interpolation_index[5]]);
    }
  }

  return 0;
}

/*
Recursive function that linearily interpolates the discrete look-up table data given in
LUT_data to a specific floating point index given in interpolation_index

Input:
    - N_dim:    
        integer that specifies the number of dimensions that still need interpolation; 
        starts with the size of the interpolation_index array for the very first call 
        of the function
    - interpolation_index:
        floating point array that specifies the index to which the data in the 
        LUTs should be interpolated. Its size should equal the number of dimensions
        in the LUTs. Given that the LUT is e.g. 2-dimensional, an interpolation index
        of [1.2, 2.7] would indicate that the data should be interpolated onto a
        point that is in between indices 1 and 2 in the first and in between indices
        2 and 3 in the second dimension
    - N_coeff:
        number of coefficients that we want to interpolate to the point given in
        interpolation_index. It must be equal to the number of LUTs stored in
        LUT_data
    - LUT_data:
        flattened array that contains the flattened n-dimensional LUT for every 
        coefficient that should be interpolated onto the specific interpolation_index. 
        For interpolating direct to direct radiation coefficients, for example, the 
        size of LUT_data would be 3 * 3 * N_LUT_tau * N_LUT_w0 * N_LUT_zx * N_LUT_g *
        N_LUT_phi * N_LUT_theta, as we have to interpolate all 9 transport 
        coefficients to the optical properties given by the interpolation_index of
        a certain grid box. 
    - LUT_offsets:
        1D integer array containing the offsets that map an n-dimensional index such as
        the interpolation_index onto the flattened LUTs in LUT_data in a such a way 
        that the dot product of LUT_offsets and an n-dimensional index array should 
        result in the corresponding flattened index.
        For N_coeff = 2 three-dimensional LUTs of shape (Nx, Ny, Nz), for example, 
        LUT_offsets would be (N_coeff, N_coeff*Nx, N_coeff*Nx*Ny), so that the dot 
        product of an index (ix, iy, iz) with LUT_offsets returns the index 
        N_coeff*ix + N_coeff*Nx*iy + N_coeff*Nx*Ny*iz of the flattened array
    - offset: 
        current offset that specifies the flattened index pointing to the beginning of
        a sequence of values of lower dimensions that we are currently looking at. 
        It should be initialized with 0 at the very first call of the function. To 
        demonstrate its values in further recursive calls of the function, we can take
        a look at an example where we want to interpolate a 2D LUT (N_coeff = 1) to 
        the point [1.2, 2.7]. Therefore, we would first like to interpolate the LUT 
        between points iy = 2 and iy = 3 at ix = 1.2. Therefore, we have to determine 
        the LUT value between points ix = 1 and ix = 2 at both iy = 2 and iy = 3. 
        Therefore, we call the recursive function again with offsets of 0 + 2*N_coeff*Nx
        and 0 + 3*N_coeff*Nx pointing to the beginning of the x-values that we have to
        consider for this interpolation one dimension down.
    - weight:
        current weight of the data point or data point range (in case N_dim != 0) that
        this call of the recursive function is pointing to. It should be initialized 
        with 1. 

Output:
    - interpolated_values:
        1D floating point array that contains the LUT values interpolated onto the
        index given by interpolation_index for all the coefficients 
*/
int interpolate_LUT (size_t  N_dim,
                     float*  interpolation_index,
                     size_t  N_coeff,
                     float*  LUT_data,
                     size_t* LUT_offsets,
                     size_t  offset,
                     float   weight,
                     float*  interpolated_values) {

  int ierr;

  /*
    GENERAL IDEA OF THE ALGORITHM:
    In dimension N_dim, the interpolation_index[N_dim-1] tells us the two points
    between which we have to interpolate the LUT. Given the distance to the closer 
    LUT point int(interpolation_index[N_dim-1]), one of these points gets more weight
    than the other. We repeat this procedure for both of these two points 
    one dimension down and calculate resulting weights for now four points in total,
    and continue with this pattern until we end up with weights for all LUT data 
    points in the hypercube surrounding the location given by the interpolation_index.

    Looking for example at a 2D interpolation, we first determine the distances to
    the two data points that are closest to the y component of the interpolation_index.
    These distances determine how much weight will be given to the top and bottom
    LUT values in this example. From thereon, we repeat this procedure and determine 
    the distances and corresponding resulting weights for the closest two points in 
    x direction for both the upper and lower y values. At this point, we have no more 
    dimension left and all the weights specify the exact amount that the four LUT data 
    points closest to the interpolation point contribute to the final result. Adding up
    these weights multiplied by their respective values results in the linarily 
    interpolated value at the interpolation location.
    */

  /*
    When we have reached a single data point in the LUT, we can add its contribution
    (= weight * value) to the respective final interpolated values for all the
    coefficients
    */

  if (N_dim == 0) {
    for (size_t i = 0; i < N_coeff; i++) {
      /*
            printf("prev value for coeff %lu: %8.5f, adding %8.5f (weight %8.5f and datapoint %5.3f) \n", 
                    i, 
                    interpolated_values[i], 
                    weight * LUT_data[offset + i], 
                    weight,
                    LUT_data[offset + i]);
            */
      interpolated_values[i] += weight * LUT_data[offset + i];
    }
    return 0;
  }

  /* 
    If we have not yet reached a single data point, we have to perform a linear
    interpolation in dimension N_dim. But before we do that, we have to check
    whether maybe the interpolation_index in this dimension is basically directly
    on a LUT data point.
    */
  int         dim_needs_interp;
  const float interp_threshold      = 1E-3;
  const float interp_index_this_dim = interpolation_index[N_dim - 1];
  if (interp_index_this_dim - (size_t) (interp_index_this_dim) <= interp_threshold) {
    dim_needs_interp = 0;
  } else if (interp_index_this_dim - (size_t) (interp_index_this_dim) >= 1.0 - interp_threshold) {
    dim_needs_interp = 0;
  } else {
    dim_needs_interp = 1;
  }

  if (dim_needs_interp) {
    /*
        We have to linearily interpolate the LUT in dimension N_dim between
        indices A = int(interpolation_index[N_dim-1]) and B = A + 1. The weights
        of points A and B are given by 1 - weight_B and weight_B, respectively 
        (B gets less weight the closer our interpolation point is to index A).
        */
    const size_t LUT_index_A = (size_t) (interpolation_index[N_dim - 1]);
    const float  weight_B    = interpolation_index[N_dim - 1] - (float) (LUT_index_A);

    /*
        We then call the function recursively again to perform the interpolation
        around points A and B one dimension down. Point A refers to a (series) of
        datapoints (series in case N_dim != 0) starting at an offset of
        LUT_offsets[N_dim-1] * LUT_index_A further down the flattened LUT with an
        additional weight of 1.0 - weight_B. Point B is pointing to a series of 
        datapoints even further down the flattened LUT with an additional weight 
        of weight_B.
        */
    ierr = interpolate_LUT (N_dim - 1,
                            interpolation_index,
                            N_coeff,
                            LUT_data,
                            LUT_offsets,
                            offset + (LUT_offsets[N_dim - 1] * LUT_index_A),
                            weight * (1.0 - weight_B),
                            interpolated_values);
    CHKERR (ierr);

    ierr = interpolate_LUT (N_dim - 1,
                            interpolation_index,
                            N_coeff,
                            LUT_data,
                            LUT_offsets,
                            offset + (LUT_offsets[N_dim - 1] * (LUT_index_A + 1)),
                            weight * weight_B,
                            interpolated_values);
    CHKERR (ierr);
  } else {
    /*
        If the interpolation index in dimension N_dim is basically directly
        on a LUT point, we simply use this data point with a weight of 1.
        */
    ierr = interpolate_LUT (N_dim - 1,
                            interpolation_index,
                            N_coeff,
                            LUT_data,
                            LUT_offsets,
                            offset + (LUT_offsets[N_dim - 1] * (size_t) (round (interpolation_index[N_dim - 1]))),
                            weight * 1.0,
                            interpolated_values);
    CHKERR (ierr);
  }
  return 0;
}

/*
Function that modifies the direct to diffuse look-up table for azimuth angles
between 90 and 360 degrees.

Input:
    - lswitch_east:
        integer value that is set to 1 for azimuth angles larger than 180 
        degrees (Sun in the East)
    - lswitch_north:
        integer value that is set to 1 for azimuth angles larger between 90
        and 270 degrees (Sun in the North)
    - N_diffuse:
        integer specifying the number of diffuse TenStream streams
    - N_direct:
        integer specifying the number of direct TenStream streams

In-/Output:
    - b:
        direct to diffuse TenStream coefficients for a specific grid box; 
        takes the address of the first element for gridbox [iz][ix][iy] 
        in b[iz][ix][iy][idiff][idir] as input
*/
int apply_dir2diff_coeff_symmetry (int lswitch_east, int lswitch_north, int N_diffuse, int N_direct, float* b) {

  float b_old[N_direct * N_diffuse];

  if (lswitch_east) {
    for (int i = 2 * N_direct; i < 6 * N_direct; i++) {
      b_old[i] = b[i];
    }
    for (int i = 0; i < N_direct; i++) {
      b[2 * N_direct + i] = b_old[3 * N_direct + i];
      b[3 * N_direct + i] = b_old[2 * N_direct + i];
      b[4 * N_direct + i] = b_old[5 * N_direct + i];
      b[5 * N_direct + i] = b_old[4 * N_direct + i];
    }
  }

  if (lswitch_north) {
    for (int i = 6 * N_direct; i < 10 * N_direct; i++) {
      b_old[i] = b[i];
    }
    for (int i = 0; i < N_direct; i++) {
      b[6 * N_direct + i] = b_old[7 * N_direct + i];
      b[7 * N_direct + i] = b_old[6 * N_direct + i];
      b[8 * N_direct + i] = b_old[9 * N_direct + i];
      b[9 * N_direct + i] = b_old[8 * N_direct + i];
    }
  }

  return 0;
}

/*
Function that assigns the interpolated diffuse to diffuse, direct to diffuse and
direct to direct TenStream coefficients for a grid box with given optical 
properties to the corresponding arrays a, b and c.

Input:
    - tau_dir and tau_diff:
        direct and diffuse optical thickness of the grid box
    - w0_dir and w0_diff:
        direct and diffuse single scattering albedo of the grid box
    - g_dir and g_diff:
        direct and diffuse asymmetry parameter of the grid box
    - aspect_zx:
        z-x-aspect ratio of the grid box
    - phi_deg:
        solar azimuth angle in degrees (must be converted to a value between 
        0 and 90 degrees using the function determine_angular_LUT_quantities)
    - theta_deg:
        solar zenith angle in degrees (must be properly calculated using the
        function determine_angular_LUT_quantities before)
    - a_LUT:
        lut struct containing all data for the diffuse to diffuse look-up tables
    - b_LUT:
        lut struct containing all data for the direct to diffuse look-up tables
    - c_LUT:
        lut struct containing all data for the direct to direct look-up tables
    - lswitch_east and lswitch_north:
        integer flags that determine in which quadrant of a 2D coordinate system 
        we are to account for the fact that TenStream LUT's are only calculated 
        for values between 0 and 90 degrees azimuth
    - N_diffuse and N_direct:
        integers quantifying the number of diffuse and direct TenStream streams

Output:
    - a, b, c:
        Interpolated TenStream coefficients for the given grid box. These variables
        should be assigned with the first element for the corresponding grid box
        in the full 5D arrays a, b or c, i.e. e.g. with &(a[iz][ix][iy][0][0]).
*/
int get_tenstream_coefficients (float  tauz_dir,
                                float  tauz_diff,
                                float  w0_dir,
                                float  w0_diff,
                                float  g_dir,
                                float  g_diff,
                                double aspect_zx,
                                float  phi_deg,
                                float  theta_deg,
                                lut*   a_LUT,
                                lut*   b_LUT,
                                lut*   c_LUT,
                                int    lswitch_east,
                                int    lswitch_north,
                                int    N_diffuse,
                                int    N_direct,
                                float* a,
                                float* b,
                                float* c) {

  int          ierr;
  const size_t N_diff_dim = 4;
  const size_t N_dir_dim  = 6;

  /*
    I: Diffuse to diffuse coefficients ("a"-coefficients)
    */
  float diff_interp_index[4];
  ierr = get_interpolation_index (tauz_diff, w0_diff, aspect_zx, g_diff, phi_deg, theta_deg, a_LUT, diff_interp_index);
  CHKERR (ierr);

  ierr = interpolate_LUT (N_diff_dim, diff_interp_index, (size_t) (N_diffuse * N_diffuse), a_LUT->LUT, a_LUT->offsets, 0, 1.0, a);
  CHKERR (ierr);

  /*
    II: Direct to diffuse coefficients ("b"-coefficients)
    */
  float dir_interp_index[6];
  ierr = get_interpolation_index (tauz_dir, w0_dir, aspect_zx, g_dir, phi_deg, theta_deg, b_LUT, dir_interp_index);
  CHKERR (ierr);

  ierr = interpolate_LUT (N_dir_dim, dir_interp_index, (size_t) (N_diffuse * N_direct), b_LUT->LUT, b_LUT->offsets, 0, 1.0, b);
  CHKERR (ierr);
  ierr = apply_dir2diff_coeff_symmetry (lswitch_east, lswitch_north, N_diffuse, N_direct, b);
  CHKERR (ierr);

  /*
    III: Direct to direct coefficients ("c"-coefficients)
    */
  ierr = interpolate_LUT (N_dir_dim, dir_interp_index, (size_t) (N_direct * N_direct), c_LUT->LUT, c_LUT->offsets, 0, 1.0, c);
  CHKERR (ierr);

  return 0;
}

/*
 Function retrieving diffuse to diffuse (a), direct to diffuse (b) and direct to direct (c)
 TenStream coefficients for a grid box with given optical properties
 */
/*
int get_tenstream_LUT_coefficients(void*       opp,
                                   int         ierr,
                                   float       tau_dir,
                                   float       tau_diff,
                                   float       omega0_dir,
                                   float       omega0_diff,
                                   float       g_dir,
                                   float       g_diff,
                                   double      aspect_zx,
                                   float       phi_deg,
                                   float       theta_deg,
                                   int         lswitch_east,
                                   int         lswitch_north,
                                   int         N_diffuse,
                                   int         N_direct,
                                   float*      a,               // address of the first of the N_diffuse * N_diffuse elements for gridbox a[iz][ix][iy] in a
                                   float*      b,               // address of the first of the N_diffuse * N_direct elements for gridbox a[iz][ix][iy] in b 
                                   float*      c) {             // address of the first of the N_direct * N_direct elements for gridbox a[iz][ix][iy] in c 
    
    // diffuse to diffuse coefficients ("a"-coefficients)
    int imode = 3;
    pprts_f2c_opp_get_coeff(
        opp,
        tau_diff,
        omega0_diff,
        g_diff,
        aspect_zx,
        phi_deg,
        theta_deg,
        imode,
        lswitch_east,
        lswitch_north,
        N_diffuse * N_diffuse,
        a, 
        &ierr);
    CHKERR(ierr);
    
    // direct to diffuse coefficients ("b"-coefficients)
    imode = 2;
    pprts_f2c_opp_get_coeff(
        opp,
        tau_dir,
        omega0_dir,
        g_dir,
        aspect_zx,
        phi_deg,
        theta_deg,
        imode,
        lswitch_east,
        lswitch_north,
        N_diffuse * N_direct,
        b, 
        &ierr);
    CHKERR(ierr);
    
    // direct to direct coefficients ("c"-coefficients)
    imode = 1;
    pprts_f2c_opp_get_coeff(
        opp,
        tau_dir,
        omega0_dir,
        g_dir,
        aspect_zx,
        phi_deg,
        theta_deg,
        imode,
        lswitch_east,
        lswitch_north,
        N_direct * N_direct,
        c, 
        &ierr);
    CHKERR(ierr);
    
    return 0;
}
*/

/* ---------------------------------------------------------------------

Other routines

---------------------------------------------------------------------- */

int delta_scale (double tau, double omega0, double g, double g_max, double* tau_scale, double* omega0_scale, double* g_scale) {

  double f = g * g;

  /*
    if (g >= g_max) {
	    // choose f so that g = g_max
    	f = (g - g_max) / (1.0 - g_max);
    }
    */

  *tau_scale    = (1.0 - f * omega0) * tau;
  *omega0_scale = (1.0 - f) * omega0 / (1.0 - f * omega0);
  *g_scale      = (g - f) / (1.0 - f);

  return 0;
}

/*
 Function that reads the values from files named
    E_direct_iv_%IV%_ib_%IB%_nlev_%NLEV%_nx_%NX%_ny_%NY%.dat
    E_diffuse_iv_%IV%_ib_%IB%_nlev_%NLEV%_nx_%NX%_ny_%NY%.dat
 where
    %IV%    is the wavelength index
    %IB%    is the subband index
    %NLEV%  is the number of vertical levels used for the current calculation
    %NX%    is the number of gridboxes in x direction used for the current calculation
    %NY%    is the number of gridboxes in y direction used for the current calculation
 and stores them as initial values in E_direct and E_diffuse
 */
int read_initial_irradiances (int       iv,
                              int       ib,
                              int       nlev,
                              int       nx,
                              int       ny,
                              int       N_direct,
                              int       N_diffuse,
                              int*      success,
                              float**** E_direct,
                              float**** E_diffuse) {

  char E_direct_filename[100], E_diffuse_filename[100];
  sprintf (E_direct_filename, "E_direct_iv_%i_ib_%i_nlev_%i_nx_%i_ny_%i.dat", iv, ib, nlev, nx, ny);
  sprintf (E_diffuse_filename, "E_diffuse_iv_%i_ib_%i_nlev_%i_nx_%i_ny_%i.dat", iv, ib, nlev, nx, ny);

  FILE* f_E_direct  = fopen (E_direct_filename, "r");
  FILE* f_E_diffuse = fopen (E_diffuse_filename, "r");

  int read_elements;
  if ((f_E_direct != NULL) && (f_E_direct != NULL)) {
    fprintf (stderr,
             "CAUTION! Using precalculated irradiances as initial values. If this is NOT what you want, please delete all E_*.dat "
             "files in your folder before starting a new calculation. \n");
    read_elements = fread (&(E_direct[0][0][0][0]), sizeof (float), nlev * nx * ny * N_direct, f_E_direct);
    read_elements = fread (&(E_diffuse[0][0][0][0]), sizeof (float), nlev * nx * ny * N_diffuse, f_E_diffuse);
    fclose (f_E_direct);
    fclose (f_E_diffuse);
    *success = 1;
  } else {
    fprintf (stderr, "Found no previously calculated file to read in. Starting dynamic tenstream calculation from scratch. \n");
  }

  return 0;
}

/*
 Function that writes the values stored in E_direct and E_diffuse to files named
    E_direct_iv_%IV%_ib_%IB%_nlev_%NLEV%_nx_%NX%_ny_%NY%.dat
    E_diffuse_iv_%IV%_ib_%IB%_nlev_%NLEV%_nx_%NX%_ny_%NY%.dat
 where
    %IV%    is the wavelength index
    %IB%    is the subband index
    %NLEV%  is the number of vertical levels used for the current calculation
    %NX%    is the number of gridboxes in x direction used for the current calculation
    %NY%    is the number of gridboxes in y direction used for the current calculation
 */
int dump_irradiances (int       iv,
                      int       ib,
                      int       nlev,
                      int       nx,
                      int       ny,
                      int       N_direct,
                      int       N_diffuse,
                      float**** E_direct,
                      float**** E_diffuse) {

  char E_direct_filename[100], E_diffuse_filename[100];
  sprintf (E_direct_filename, "E_direct_iv_%i_ib_%i_nlev_%i_nx_%i_ny_%i.dat", iv, ib, nlev, nx, ny);
  sprintf (E_diffuse_filename, "E_diffuse_iv_%i_ib_%i_nlev_%i_nx_%i_ny_%i.dat", iv, ib, nlev, nx, ny);

  FILE* f_E_direct  = fopen (E_direct_filename, "w");
  FILE* f_E_diffuse = fopen (E_diffuse_filename, "w");

  if ((f_E_direct != NULL) && (f_E_direct != NULL)) {
    fwrite (&(E_direct[0][0][0][0]), sizeof (float), nlev * nx * ny * N_direct, f_E_direct);
    fwrite (&(E_diffuse[0][0][0][0]), sizeof (float), nlev * nx * ny * N_diffuse, f_E_diffuse);
    fclose (f_E_direct);
    fclose (f_E_diffuse);
  } else {
    fprintf (stderr, "Could not open files to dump calculated direct and diffuse irradiance fields. \n");
  }

  return 0;
}

/*
 Functions calculating effective thermal irradiance for each of the ten streams
 for a given model grid box
 */
int thermal_emission (float B_top, float B_bottom, float tau, int number_of_angles, float* B_eff) {

  double B_eff_top = 0.0, B_eff_bottom = 0.0;

  for (int i = 0; i < number_of_angles; i++) {

    double mu = 1.0 / (2.0 * number_of_angles) + i * 1.0 / number_of_angles;

    if (tau / mu < 1E-6) {
      B_eff_top += 2 * M_PI * (B_top + B_bottom) * 0.5 * mu * (1.0 / number_of_angles);
      B_eff_bottom += 2 * M_PI * (B_top + B_bottom) * 0.5 * mu * (1.0 / number_of_angles);
    } else {
      B_eff_top += 2 * M_PI * ((B_bottom * exp (-tau / mu) - B_top) / (expm1 (-tau / mu)) + (B_bottom - B_top) / (tau / mu)) * mu *
                   (1.0 / number_of_angles);
      B_eff_bottom += 2 * M_PI * ((B_top * exp (-tau / mu) - B_bottom) / (expm1 (-tau / mu)) + (B_top - B_bottom) / (tau / mu)) *
                      mu * (1.0 / number_of_angles);
    }
  }

  B_eff[0] = B_eff_top;
  B_eff[1] = B_eff_bottom;
  B_eff[2] = B_eff_bottom / 2.0;
  B_eff[3] = B_eff_bottom / 2.0;
  B_eff[4] = B_eff_top / 2.0;
  B_eff[5] = B_eff_top / 2.0;
  B_eff[6] = B_eff_bottom / 2.0;
  B_eff[7] = B_eff_bottom / 2.0;
  B_eff[8] = B_eff_top / 2.0;
  B_eff[9] = B_eff_top / 2.0;

  return 0;
}

static void thermal_emission_2 (double B0, double B1, float dtau_ext, float g_Parameter, float omega_0, float e, float* B_eff) {

  double alpha1;
  double alpha2;
  if (omega_0 == 0) {
    alpha1 = 1.75;
    alpha2 = -0.25;
  } else {
    alpha1 = (1.0 - omega_0) + 0.75 * (1.0 - omega_0 * g_Parameter);
    alpha2 = -(1.0 - omega_0) + 0.75 * (1.0 - omega_0 * g_Parameter);
  }

  double lambda = sqrt (alpha1 * alpha1 - alpha2 * alpha2);

  double A;
  double a11;
  double a12;
  if (alpha2 != 0) {
    A   = 1.0 / (alpha2 / (alpha1 - lambda) * exp (lambda * dtau_ext) - alpha2 / (alpha1 + lambda) * exp (-lambda * dtau_ext));
    a11 = A * 2.0 * lambda / alpha2;
    a12 = A * (exp (lambda * dtau_ext) - exp (-lambda * dtau_ext));
  } else {  // special solution for alpha2 == 0 which causes A = 0/0 == nan
    a11 = exp (-lambda * dtau_ext);
    a12 = 0;
  }

  double kappa  = 2.0 * M_PI * (1.0 - omega_0);
  double alpha7 = kappa * (B0 / (alpha1 - alpha2) + B1 / lambda / lambda);
  double alpha9 = kappa * (B0 / (alpha1 - alpha2) - B1 / lambda / lambda);

  double alpha8  = kappa * B1 / (alpha1 - alpha2);
  double alpha10 = alpha8;

  float B_eff_top    = -a11 * (alpha7 + alpha8 * dtau_ext) - a12 * (alpha9) + alpha7;
  float B_eff_bottom = -a12 * (alpha7 + alpha8 * dtau_ext) - a11 * (alpha9) + alpha9 + alpha10 * dtau_ext;

  B_eff[0] = B_eff_top / e;
  B_eff[1] = B_eff_bottom / e;
  B_eff[2] = B_eff_bottom / (e * 2.0);
  B_eff[3] = B_eff_bottom / (e * 2.0);
  B_eff[4] = B_eff_top / (e * 2.0);
  B_eff[5] = B_eff_top / (e * 2.0);
  B_eff[6] = B_eff_bottom / (e * 2.0);
  B_eff[7] = B_eff_bottom / (e * 2.0);
  B_eff[8] = B_eff_top / (e * 2.0);
  B_eff[9] = B_eff_top / (e * 2.0);
}

/* Function adjusting loop indices to account for the sun direction and whether you are in thermal 
   spectral range or not
 */
int adjust_loop_indices (int ix0, int iy0, int iz0, int phi0, int planck, int nx, int ny, int nlyr, int* ix, int* iy, int* iz) {

  /* if the sun is in the East, iterate from the East to the West rather than from the West to
       to the East */
  if (phi0 > 180.0) {
    *ix = nx - 1 - ix0;
  } else {
    *ix = ix0;
  }

  /* if the sun is in the North, iterate from North to South rather than from the South to the 
       North */
  if (phi0 > 90.0 && phi0 <= 270.0) {
    *iy = ny - 1 - iy0;
  } else {
    *iy = iy0;
  }

  /* if we are in the thermal spectral range, iterate from bottom to top rather than from top
       to bottom */
  if (planck) {
    *iz = iz0;
  } else {
    *iz = nlyr - 1 - iz0;
  }

  return 0;
}

/* Function determining loop indices considering the number of iterations performed, the sun
   direction in the solar spectral range and whether you are in thermal spectral range or not
 */
int determine_loop_indices (int           nit,
                            int           phi0,
                            int           planck,
                            int           nx,
                            int           ny,
                            int           nlyr,
                            unsigned int* ix_start,
                            unsigned int* ix_n_it,
                            unsigned int* ix_increment,
                            unsigned int* iy_start,
                            unsigned int* iy_n_it,
                            unsigned int* iy_increment,
                            unsigned int* iz_start,
                            unsigned int* iz_n_it,
                            unsigned int* iz_increment) {

  /* if the sun is in the East, start to iterate from the East to the West rather than from
       the West to to the East and alternate that iteration direction from there on */
  if ((phi0 > 180.0 && nit % 2 == 0) || (phi0 <= 180.0 && nit % 2 == 1)) {
    *ix_start     = nx - 1;
    *ix_increment = -1;
  } else {
    *ix_start     = 0;
    *ix_increment = 1;
  }
  *ix_n_it = nx;

  /* if the sun is in the North, start to iterate from the North to South rather than from 
       the South to the North and alternate that iteration direction from there on */
  if ((phi0 > 90.0 && phi0 <= 270.0 && nit % 2 == 0) || ((phi0 <= 90.0 || phi0 > 270.0) && nit % 2 == 1)) {
    *iy_start     = ny - 1;
    *iy_increment = -1;
  } else {
    *iy_start     = 0;
    *iy_increment = 1;
  }
  *iy_n_it = ny;

  /* if we are in the thermal spectral range, start to iterate from bottom to top rather than
       from top to bottom and alternate that iteration direction from there on */
  if ((planck && nit % 2 == 0) || (!planck && nit % 2 == 1)) {
    *iz_start     = 0;
    *iz_increment = 1;
  } else {
    *iz_start     = nlyr - 1;
    *iz_increment = -1;
  }
  *iz_n_it = nlyr;

  return 0;
}

/* Function returning ecRad Eddington coefficients as implemented in eddington.F90 in the
   TenStream solver
 */
int eddington_ec (double  dtau,
                  double  g,
                  double  omega0,
                  double  mu0,
                  double* t,
                  double* r,
                  double* rdir,
                  double* sdir,
                  double* tdir) {

  double f  = 0;
  double g1 = 0, g2 = 0, g3 = 0, g4 = 0;
  double alpha1 = 0, alpha2 = 0, A = 0, beta = 0;
  double e0 = 0, e = 0, e2 = 0;
  double k_mu0 = 0, k_g3 = 0, k_g4 = 0, k_2_e = 0;
  double dtau_slant = 0;

  f  = 0.75 * g;
  g1 = 2.0 - (1.25 + f) * omega0;
  g2 = (0.75 - f) * omega0;
  g3 = 0.5 - mu0 * f;

  dtau_slant = MAX (dtau / MAX (sqrt (DBL_MIN), mu0), 0.0);

  if (dtau_slant > 1E-6) {
    g4     = 1.0 - g3;
    alpha1 = g1 * g4 + g2 * g3;
    alpha2 = g1 * g3 + g2 * g4;

    A     = sqrt (MAX ((g1 - g2) * (g1 + g2), 1E-12));
    k_mu0 = A * mu0;
    k_g3  = A * g3;
    k_g4  = A * g4;

    e0    = exp (-dtau_slant);
    *tdir = e0;

    e     = exp (-A * dtau);
    e2    = e * e;
    k_2_e = 2.0 * A * e;
    if (fabs (k_mu0 - 1.0) < DBL_EPSILON) {
      k_mu0 = 1.0 - 10.0 * DBL_EPSILON;
    }
    beta = 1.0 / (A + g1 + (A - g1) * e2);

    *r = g2 * (1.0 - e2) * beta;
    *t = k_2_e * beta;

    beta = omega0 * beta / (1.0 - k_mu0 * k_mu0);

    *sdir = beta * (k_2_e * (g4 + alpha1 * mu0) - e0 * ((1 + k_mu0) * (alpha1 + k_g4) - (1.0 - k_mu0) * (alpha1 - k_g4) * e2));
    *rdir = beta * ((1.0 - k_mu0) * (alpha2 + k_g3) - (1 + k_mu0) * (alpha2 - k_g3) * e2 - k_2_e * (g3 - alpha2 * mu0) * e0);
  } else {
    *t    = (1.0 - g1 * dtau);
    *r    = g2 * dtau;
    *sdir = (1.0 - g3) * (omega0 * dtau);
    *rdir = g3 * (omega0 * dtau);
    *tdir = 1.0 - dtau_slant;
  }
  return 0;
}

/*
 Function calculating the emissivity of a certain grid box for each of the ten streams
 */
int calculate_emissivities (int N_diffuse, float** a, float* e) {

  int opposite_index[10] = {1, 0, 5, 4, 3, 2, 9, 8, 7, 6};

  for (int idst = 0; idst < N_diffuse; idst++) {
    e[idst] = 1.0;
    for (int isrc = 0; isrc < N_diffuse; isrc++) {
      e[idst] -= a[isrc][opposite_index[idst]];
    }
  }

  return 0;
}

/*
 Function determining effective Planck values as well as the TenStream coefficients a, b, c and e for the
 whole domain using either Eddington or TenStream LUT coefficients (depending on zx aspect ratio)
 */
int determine_tenstream_coefficients_and_effective_B (int                N_direct,
                                                      int                N_diffuse,
                                                      int                nlyr,
                                                      int                nx,
                                                      int                ny,
                                                      double             mu0,
                                                      double             phi0,
                                                      int                delta,
                                                      int                planck,
                                                      float*             zd,
                                                      float*             temper,
                                                      float              wvnmlo,
                                                      float              wvnmhi,
                                                      float*             dtau_org,
                                                      float*             omega0_org,
                                                      float*             g_org,
                                                      int                found3D_ic,
                                                      int                found3D_wc,
                                                      caoth3d_out_struct ic3D,
                                                      caoth3d_out_struct wc3D,
                                                      char*              datapath,
                                                      int                quiet,
                                                      float*****         a,
                                                      float*****         b,
                                                      float*****         c,
                                                      float****          e,
                                                      float****          B_eff) {

  int ierr;

  float theta_deg = 0, phi_deg = 0; /* solar zenith and azimuth angle for TenStream Lookup Tables */

  int lswitch_east = 0, lswitch_north = 0; /* variables that are not 0 if the TenStream solar azimuth angle is
                                                 initially larger than 90 degrees to account for the fact that the
                                                 LUT's are only calculated for angles between 0 and 90 degrees azimuth */
  ierr = determine_angular_LUT_quantities (mu0, phi0, &theta_deg, &phi_deg, &lswitch_east, &lswitch_north);
  CHKERR (ierr);

  /*

    // VERSION 1: TenStream Interface
    const int mpi_comm  = 0;
    const int solver_id = SOLVER_ID_PPRTS_3_10;

    static void* opp = NULL;  // this is pointer to the tenstream optical properties object

    static float diff_tauz_range[2];
    static float diff_w0_range[2];
    static float diff_g_range[2];
    static float diff_aspect_zx_range[2];
    static float dir_tauz_range[2];
    static float dir_w0_range[2];
    static float dir_g_range[2];
    static float dir_aspect_zx_range[2];
    static float dir_phi_range[2];
    static float dir_theta_range[2];

    if (!opp) {
        pprts_f2c_opp_init(mpi_comm, solver_id, &opp, &ierr);
        CHKERR(ierr);

        pprts_f2c_opp_get_info(opp,
                               &N_direct,
                               &N_diffuse,
                               diff_tauz_range,
                               diff_w0_range,
                               diff_g_range,
                               diff_aspect_zx_range,
                               dir_tauz_range,
                               dir_w0_range,
                               dir_g_range,
                               dir_aspect_zx_range,
                               dir_phi_range,
                               dir_theta_range,
                               &ierr);
        CHKERR(ierr);
    }
    */
  /*
    fprintf(stderr, "LUT has %d direct and %d diffuse streams\n",         N_direct,               N_diffuse);
    fprintf(stderr, "LUT range %20s :: %12.6g %12.6g\n", "diff_tauz",     diff_tauz_range[0],     diff_tauz_range[1]);
    fprintf(stderr, "LUT range %20s :: %12.6g %12.6g\n", "diff_w0",       diff_w0_range[0],       diff_w0_range[1]);
    fprintf(stderr, "LUT range %20s :: %12.6g %12.6g\n", "diff_g",        diff_g_range[0],        diff_g_range[1]);
    fprintf(stderr, "LUT range %20s :: %12.6g %12.6g\n", "diff_aspect_zx",diff_aspect_zx_range[0],diff_aspect_zx_range[1]);
    fprintf(stderr, "LUT range %20s :: %12.6g %12.6g\n", "dir_tauz",      dir_tauz_range[0],      dir_tauz_range[1]);
    fprintf(stderr, "LUT range %20s :: %12.6g %12.6g\n", "dir_w0",        dir_w0_range[0],        dir_w0_range[1]);
    fprintf(stderr, "LUT range %20s :: %12.6g %12.6g\n", "dir_g",         dir_g_range[0],         dir_g_range[1]);
    fprintf(stderr, "LUT range %20s :: %12.6g %12.6g\n", "dir_aspect_zx", dir_aspect_zx_range[0], dir_aspect_zx_range[1]);
    fprintf(stderr, "LUT range %20s :: %12.6g %12.6g\n", "dir_phi",       dir_phi_range[0],       dir_phi_range[1]);
    fprintf(stderr, "LUT range %20s :: %12.6g %12.6g\n", "dir_theta",     dir_theta_range[0],     dir_theta_range[1]);
    */

  // **********************************************
  // Version 2: TenStream Independent
  // **********************************************

  static lut* a_LUT = NULL;
  static lut* b_LUT = NULL;
  static lut* c_LUT = NULL;
  char        fileloc[1000];

  if (!a_LUT) {
    if (!quiet) {
      fprintf (stderr, "Reading diff2diff LUT \n");
    }
    strcpy (fileloc, datapath);
    strcat (fileloc, "tenstream_lookup/LUT_diffuse_10.tau31.w020.aspect_zx23.g6.ds1000.nc");
    a_LUT = calloc (1, sizeof (lut));
    ierr  = get_LUT_data ("diff2diff", fileloc, N_direct, N_diffuse, quiet, a_LUT);
    CHKERR (ierr);
  }
  if (!b_LUT) {
    if (!quiet) {
      fprintf (stderr, "Reading dir2diff LUT \n");
    }
    strcpy (fileloc, datapath);
    strcat (fileloc, "tenstream_lookup/LUT_direct_3_10.tau31.w020.aspect_zx23.g6.phi19.theta19.ds1000.nc");
    b_LUT = calloc (1, sizeof (lut));
    ierr  = get_LUT_data ("dir2diff", fileloc, N_direct, N_diffuse, quiet, b_LUT);
    CHKERR (ierr);
  }
  if (!c_LUT) {
    if (!quiet) {
      fprintf (stderr, "Reading dir2dir LUT \n");
    }
    strcpy (fileloc, datapath);
    strcat (fileloc, "tenstream_lookup/LUT_direct_3_10.tau31.w020.aspect_zx23.g6.phi19.theta19.ds1000.nc");
    c_LUT = calloc (1, sizeof (lut));
    ierr  = get_LUT_data ("dir2dir", fileloc, N_direct, N_diffuse, quiet, c_LUT);
    CHKERR (ierr);
  }

  // Determine effective Planck value and TenStream coefficients a, b, c, e for every grid box of the domain
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nlyr; iz++) {

        /* z-index in 1D, since 1D profiles are sorted from top to bottom,
                 while 3D fields are sorted from bottom to top */
        int iz1D = nlyr - 1 - iz;

        /* vertical gridbox thickness in meters */
        double delta_z = (zd[iz1D] - zd[iz1D + 1]) * 1000.0;

        /* maximum allowed asymmetry parameter g_max according to LUT */
        //double g_max = MIN(dir_g_range[1], diff_g_range[1]);
        double g_max = MIN (a_LUT->g[0], c_LUT->g[0]);

        // Calculation of total optical thickness, single scattering albedo and asymmetry parameter
        double tau_wc = 0, omega0_wc = 0, g_wc = 0;
        if (found3D_wc) {
          if (wc3D.threed[iz]) {

            tau_wc    = wc3D.ext[iz][ix][iy] * delta_z;
            omega0_wc = wc3D.ssa[iz][ix][iy];
            g_wc      = wc3D.g1[iz][ix][iy];

            if (delta) {
              delta_scale (tau_wc, omega0_wc, g_wc, g_max, &tau_wc, &omega0_wc, &g_wc);
            }
          }
        }

        double tau_ic = 0, omega0_ic = 0, g_ic = 0;
        if (found3D_ic) {
          if (ic3D.threed[iz]) {

            tau_ic    = ic3D.ext[iz][ix][iy] * delta_z;
            omega0_ic = ic3D.ssa[iz][ix][iy];
            g_ic      = ic3D.g1[iz][ix][iy];

            if (delta) {
              delta_scale (tau_ic, omega0_ic, g_ic, g_max, &tau_ic, &omega0_ic, &g_ic);
            }
          }
        }

        float tau    = dtau_org[iz1D] + tau_wc + tau_ic;
        float omega0 = (omega0_org[iz1D] * dtau_org[iz1D] + omega0_wc * tau_wc + omega0_ic * tau_ic) / tau;
        float g      = (omega0_org[iz1D] * dtau_org[iz1D] * g_org[iz1D] + omega0_wc * tau_wc * g_wc + omega0_ic * tau_ic * g_ic) /
                  (omega0 * tau);

        // TenStream coefficients
        double aspect_zx = 0;
        if (wc3D.delX != 0) {
          aspect_zx = delta_z / wc3D.delX;
        } else if (ic3D.delX != 0) {
          aspect_zx = delta_z / ic3D.delX;
        }

        if (aspect_zx >= 2.0) {

          /* use 1D Eddington coefficients for grid boxes with an z-x-aspect ratio of larger than 2 */

          const double taumax = 100.0;  // maximum layer optical thickness allowed to avoid problems with the Eddington code
          double       r, t, r_dir, t_dir, s_dir;

          tau = MIN (taumax, tau);
          eddington_ec (tau, g, omega0, mu0, &t, &r, &r_dir, &s_dir, &t_dir);

          a[iz][ix][iy][0][0] = t;
          a[iz][ix][iy][0][1] = r;
          a[iz][ix][iy][1][0] = r;
          a[iz][ix][iy][1][1] = t;
          b[iz][ix][iy][0][0] = r_dir;
          b[iz][ix][iy][1][0] = s_dir;
          c[iz][ix][iy][0][0] = t_dir;

        } else {
          /*
                    if (tau > diff_tauz_range[1]) {
                        delta_scale(tau, omega0, g, 0.0, &tau, &omega0, &g);
                    }
                    */

          /* Set the optical properties of the grid box to that of the tenstream LUT's borders if necessary */
          /*
                    float g_dir = MAX(dir_g_range[0], MIN(g, dir_g_range[1]));
                    float g_diff = MAX(diff_g_range[0], MIN(g, diff_g_range[1]));
                    float omega0_dir = MAX(dir_w0_range[0], MIN(omega0, dir_w0_range[1]));
                    float omega0_diff = MAX(diff_w0_range[0], MIN(omega0, diff_w0_range[1]));
                    float tau_dir = MAX(dir_tauz_range[0], MIN(tau, dir_tauz_range[1]));
                    float tau_diff = MAX(diff_tauz_range[0], MIN(tau, diff_tauz_range[1]));
      
                    ierr = get_tenstream_LUT_coefficients(opp,
                                                          ierr,
                                                          tau_dir,
                                                          tau_diff,
                                                          omega0_dir,
                                                          omega0_diff,
                                                          g_dir,
                                                          g_diff,
                                                          aspect_zx,
                                                          phi_deg,
                                                          theta_deg,
                                                          lswitch_east,
                                                          lswitch_north,
                                                          N_diffuse,
                                                          N_direct,
                                                          &(a[iz][ix][iy][0][0]),
                                                          &(b[iz][ix][iy][0][0]),
                                                          &(c[iz][ix][iy][0][0]));
                    CHKERR(ierr);
                    */

          // **********************************************
          // Version 2: TenStream Independent
          // **********************************************

          float tau_dir     = MAX (c_LUT->tauz[0], MIN (tau, c_LUT->tauz[c_LUT->N_tauz - 1]));
          float tau_diff    = MAX (a_LUT->tauz[0], MIN (tau, a_LUT->tauz[a_LUT->N_tauz - 1]));
          float omega0_dir  = MAX (c_LUT->w0[0], MIN (omega0, c_LUT->w0[c_LUT->N_w0 - 1]));
          float omega0_diff = MAX (a_LUT->w0[0], MIN (omega0, a_LUT->w0[a_LUT->N_w0 - 1]));
          float g_dir       = MAX (c_LUT->g[0], MIN (g, c_LUT->g[c_LUT->N_g - 1]));
          float g_diff      = MAX (a_LUT->g[0], MIN (g, a_LUT->g[a_LUT->N_g - 1]));

          ierr = get_tenstream_coefficients (tau_dir,
                                             tau_diff,
                                             omega0_dir,
                                             omega0_diff,
                                             g_dir,
                                             g_diff,
                                             aspect_zx,
                                             phi_deg,
                                             theta_deg,
                                             a_LUT,
                                             b_LUT,
                                             c_LUT,
                                             lswitch_east,
                                             lswitch_north,
                                             N_diffuse,
                                             N_direct,
                                             &(a[iz][ix][iy][0][0]),
                                             &(b[iz][ix][iy][0][0]),
                                             &(c[iz][ix][iy][0][0]));
          CHKERR (ierr);

          /*
                    fprintf(stderr, "phi = %g, theta = %g \n", phi_deg, theta_deg);
                    for (int first = 0; first < N_diffuse; first++) {
                        for (int second = 0; second < N_diffuse; second++) {
                            fprintf(stderr, "a[%i][%i] = %g \n", first, second, a[iz][ix][iy][first][second]);
                        }
                    }
                    for (int first = 0; first < N_direct; first++) {
                        for (int second = 0; second < N_diffuse; second++) {
                            fprintf(stderr, "b[%i][%i] = %g \n", first, second, b[iz][ix][iy][first][second]);
                        }
                    }
                    for (int first = 0; first < N_direct; first++) {
                        for (int second = 0; second < N_direct; second++) {
                            fprintf(stderr, "c[%i][%i] = %g \n", first, second, c[iz][ix][iy][first][second]);
                        }
                    }
                    */
        }

        /* Calculation of thermal emissivities using Kirchhoff's law: emissivity in a certain direction is absorptivity in the opposite
                 direction, which in turn is one minus the transmittance in that direction */
        if (planck) {
          ierr = calculate_emissivities (N_diffuse, a[iz][ix][iy], e[iz][ix][iy]);
          CHKERR (ierr);
        }

        /* Calculation of effective Planck values */
        float plkavg_up = 0.0, plkavg_dn = 0.0;
        // if (planck && ix == 16 && iy == 16 && iz == 125) {
        if (planck) {
          plkavg_up = c_planck_func1 (wvnmlo, wvnmhi, temper[iz1D]);
          plkavg_dn = c_planck_func1 (wvnmlo, wvnmhi, temper[iz1D + 1]);

          int number_of_angles = 1;
          ierr                 = thermal_emission (plkavg_up, plkavg_dn, tau, number_of_angles, B_eff[iz][ix][iy]);
          CHKERR (ierr);
        }
      }
    }
  }
  //pprts_f2c_opp_destroy(opp, &ierr);
  //CHKERR(ierr);

  return 0;
}

/*
 Function initializing the borders of the updated direct and diffuse irradiance fields
 (_new fields) for JACOBI iterations where iteration goes from TOP to BOTTOM
 (starting at z = nlev - 1), WEST to EAST (starting at x = 0) and SOUTH to NORTH
 (starting at y = 0) with the corresponding values of the direct and diffuse irradiance
 fields of the previous time step and the solar irradiance at the top
 */
int initialize_borders (int       nlev,
                        int       nx,
                        int       ny,
                        double    S0,
                        double    mu0,
                        float**** E_direct,
                        float**** E_diffuse,
                        float**** E_direct_new,
                        float**** E_diffuse_new) {

  for (int iz = 0; iz < nlev; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {

        // Boundary Condition at the upper border of the domain
        /*
                // point source
                if (iz == nlev - 1 && ix == (int)(nx/2) && iy == (int)(ny/2)) {
                    E_direct_new[iz][ix][iy][0] = S0 * mu0;
                }
                */

        if (iz == nlev - 1) {
          E_direct_new[iz][ix][iy][0] = S0 * mu0;
        }

        if (ix == 0 && iz != nlev - 1) {
          E_diffuse_new[iz][ix][iy][3] = E_diffuse[iz][ix][iy][3];
          E_diffuse_new[iz][ix][iy][5] = E_diffuse[iz][ix][iy][5];

          E_diffuse_new[iz][ix][iy][2] = E_diffuse[iz][ix][iy][2];
          E_diffuse_new[iz][ix][iy][4] = E_diffuse[iz][ix][iy][4];

          E_direct_new[iz][ix][iy][1] = E_direct[iz][ix][iy][1];
        }

        if (iy == 0 && iz != nlev - 1) {
          E_diffuse_new[iz][ix][iy][7] = E_diffuse[iz][ix][iy][7];
          E_diffuse_new[iz][ix][iy][9] = E_diffuse[iz][ix][iy][9];

          E_diffuse_new[iz][ix][iy][6] = E_diffuse[iz][ix][iy][6];
          E_diffuse_new[iz][ix][iy][8] = E_diffuse[iz][ix][iy][8];

          E_direct_new[iz][ix][iy][2] = E_direct[iz][ix][ix][2];
        }
      }
    }
  }
  return 0;
}

/* Function determining the ten diffuse and three direct incoming RADIATIVE FLUXES
 for a certain grid box (iz, ix, iy) and assigning them to an array E_in of length 13 */
int determine_E_in (int       ix,
                    int       iy,
                    int       iz,
                    int       nx,
                    int       ny,
                    double    dx,
                    double    dy,
                    double    dz,
                    double    phi0,
                    float**** E_direct,
                    float**** E_diffuse,
                    float*    E_in) {

  /*  ingoing diffuse fluxes */
  E_in[0]  = E_diffuse[iz][ix % nx][iy % ny][0] * dx * dy;
  E_in[1]  = E_diffuse[iz + 1][ix % nx][iy % ny][1] * dx * dy;
  E_in[2]  = E_diffuse[iz][(ix + 1) % nx][iy % ny][2] * dy * dz;
  E_in[3]  = E_diffuse[iz][ix % nx][iy % ny][3] * dy * dz;
  E_in[4]  = E_diffuse[iz][(ix + 1) % nx][iy % ny][4] * dy * dz;
  E_in[5]  = E_diffuse[iz][ix % nx][iy % ny][5] * dy * dz;
  E_in[6]  = E_diffuse[iz][ix % nx][(iy + 1) % ny][6] * dx * dz;
  E_in[7]  = E_diffuse[iz][ix % nx][iy % ny][7] * dx * dz;
  E_in[8]  = E_diffuse[iz][ix % nx][(iy + 1) % ny][8] * dx * dz;
  E_in[9]  = E_diffuse[iz][ix % nx][iy % ny][9] * dx * dz;
  E_in[10] = E_direct[iz + 1][ix % nx][iy % ny][0] * dx * dy;

  /* ingoing direct fluxes */
  if (phi0 <= 90.0) {
    E_in[11] = E_direct[iz][ix % nx][iy % ny][1] * dy * dz;
    E_in[12] = E_direct[iz][ix % nx][iy % ny][2] * dx * dz;
  } else if (phi0 > 90.0 && phi0 <= 180.0) {
    E_in[11] = E_direct[iz][ix % nx][iy % ny][1] * dy * dz;
    E_in[12] = E_direct[iz][ix % nx][(iy + 1) % ny][2] * dx * dz;
  } else if (phi0 > 180.0 && phi0 <= 270.0) {
    E_in[11] = E_direct[iz][(ix + 1) % nx][iy % ny][1] * dy * dz;
    E_in[12] = E_direct[iz][ix % nx][(iy + 1) % ny][2] * dx * dz;
  } else {
    E_in[11] = E_direct[iz][(ix + 1) % nx][iy % ny][1] * dy * dz;
    E_in[12] = E_direct[iz][ix % nx][iy % ny][2] * dx * dz;
  }

  return 0;
}

/* function applying a Jacobian interation to update the direct and diffuse
 irradiance fields */
int one_jacobian_iteration (int        nx,
                            int        ny,
                            int        nlyr,
                            int        N_direct,
                            int        N_diffuse,
                            double     S0,
                            double     mu0,
                            double     phi0,
                            double     Ag,
                            float      Bg,
                            float****  B_eff,
                            double     dx,
                            double     dy,
                            float*     zd,
                            float***** a,
                            float***** b,
                            float***** c,
                            float****  e,
                            float****  E_direct_new,
                            float****  E_diffuse_new,
                            float****  E_direct,
                            float****  E_diffuse) {

  int ierr;
  int nlev = nlyr + 1;

  // Boundary Conditions at iz == nlev - 1, ix == 0 and iy == 0
  ierr = initialize_borders (nlev, nx, ny, S0, mu0, E_direct, E_diffuse, E_direct_new, E_diffuse_new);
  CHKERR (ierr);

  // ONE JACOBIAN ITERATION FORWARD in the dynamic tenstream algorithm
  for (int iz = nlyr - 1; iz >= 0; iz--) {
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {

        // z-index in 1D, since 1D profiles are sorted from top to bottom,
        // while 3D fields are sorted from bottom to top
        int iz1D = nlyr - 1 - iz;

        // Dynamic Radiative Transfer: For every grid box, new output irradiances fluxes are calculated using the input irradiances of the previous timestep
        double dz = (zd[iz1D] - zd[iz1D + 1]) * 1.0E3;

        float E_in[N_diffuse + N_direct];
        ierr = determine_E_in (ix, iy, iz, nx, ny, dx, dy, dz, phi0, E_direct, E_diffuse, E_in);
        CHKERR (ierr);

        float E_out[N_diffuse + N_direct];
        for (int i = 0; i < (N_diffuse + N_direct); i++) {
          E_out[i] = 0;
        }

        for (int idiff = 0; idiff < N_diffuse; idiff++) {
          for (int jdiff = 0; jdiff < N_diffuse; jdiff++) {
            E_out[idiff] += a[iz][ix][iy][idiff][jdiff] * E_in[jdiff];
          }
          for (int jdir = 0; jdir < N_direct; jdir++) {
            E_out[idiff] += b[iz][ix][iy][idiff][jdir] * E_in[N_diffuse + jdir];
          }
          if (idiff == 0 || idiff == 1) {
            E_out[idiff] += e[iz][ix][iy][idiff] * B_eff[iz][ix][iy][idiff] * dx * dy;
          } else if (idiff >= 2 && idiff <= 5) {
            E_out[idiff] += e[iz][ix][iy][idiff] * B_eff[iz][ix][iy][idiff] * dy * dz;
          } else {
            E_out[idiff] += e[iz][ix][iy][idiff] * B_eff[iz][ix][iy][idiff] * dx * dz;
          }
        }
        for (int idir = 0; idir < N_direct; idir++) {
          for (int jdir = 0; jdir < N_direct; jdir++) {
            E_out[N_diffuse + idir] += c[iz][ix][iy][idir][jdir] * E_in[N_diffuse + jdir];
          }
        }

        E_diffuse_new[iz + 1][ix % nx][iy % ny][0]   = E_out[0] / (dx * dy);
        E_diffuse_new[iz][ix % nx][iy % ny][1]       = E_out[1] / (dx * dy);
        E_diffuse_new[iz][ix % nx][iy % ny][2]       = E_out[2] / (dy * dz);
        E_diffuse_new[iz][(ix + 1) % nx][iy % ny][3] = E_out[3] / (dy * dz);
        E_diffuse_new[iz][ix % nx][iy % ny][4]       = E_out[4] / (dy * dz);
        E_diffuse_new[iz][(ix + 1) % nx][iy % ny][5] = E_out[5] / (dy * dz);
        E_diffuse_new[iz][ix % nx][iy % ny][6]       = E_out[6] / (dx * dz);
        E_diffuse_new[iz][ix % nx][(iy + 1) % ny][7] = E_out[7] / (dx * dz);
        E_diffuse_new[iz][ix % nx][iy % ny][8]       = E_out[8] / (dx * dz);
        E_diffuse_new[iz][ix % nx][(iy + 1) % ny][9] = E_out[9] / (dx * dz);
        E_direct_new[iz][ix % nx][iy % ny][0]        = E_out[10] / (dx * dy);

        if (phi0 <= 90.0) {
          E_direct_new[iz][(ix + 1) % nx][iy % ny][1] = E_out[11] / (dy * dz);
          E_direct_new[iz][ix % nx][(iy + 1) % ny][2] = E_out[12] / (dx * dz);
        } else if (phi0 > 90.0 && phi0 <= 180.0) {
          E_direct_new[iz][(ix + 1) % nx][iy % ny][1] = E_out[11] / (dy * dz);
          E_direct_new[iz][ix % nx][iy % ny][2]       = E_out[12] / (dx * dz);
        } else if (phi0 > 180.0 && phi0 <= 270.0) {
          E_direct_new[iz][ix % nx][iy % ny][1] = E_out[11] / (dy * dz);
          E_direct_new[iz][ix % nx][iy % ny][2] = E_out[12] / (dx * dz);
        } else {
          E_direct_new[iz][ix % nx][iy % ny][1]       = E_out[11] / (dy * dz);
          E_direct_new[iz][ix % nx][(iy + 1) % ny][2] = E_out[12] / (dx * dz);
        }
      }
    }
  }

  // Boundary Conditions II: iz == 0
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      E_diffuse_new[0][ix][iy][0] = Ag * (E_diffuse[0][ix][iy][1] + E_direct[0][ix][iy][0]) + (1.0 - Ag) * M_PI * Bg;
    }
  }

  // Copying to old fields
  for (int iz = nlev - 1; iz >= 0; iz--) {
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
        for (int idir = 0; idir < N_direct; idir++) {
          E_direct[iz][ix][iy][idir] = E_direct_new[iz][ix][iy][idir];
        }
        for (int idiff = 0; idiff < N_diffuse; idiff++) {
          E_diffuse[iz][ix][iy][idiff] = E_diffuse_new[iz][ix][iy][idiff];
        }
      }
    }
  }

  return 0;
}

/* function applying a Gauss-Seidel interation to update the direct and diffuse
 irradiance fields */
int one_gauss_seidel_iteration (int        nit,
                                int        nx,
                                int        ny,
                                int        nlyr,
                                int        N_direct,
                                int        N_diffuse,
                                int        planck,
                                double     S0,
                                double     mu0,
                                double     phi0,
                                double     Ag,
                                float      Bg,
                                float****  B_eff,
                                double     dx,
                                double     dy,
                                float*     zd,
                                float***** a,
                                float***** b,
                                float***** c,
                                float****  e,
                                float****  E_direct,
                                float****  E_diffuse) {

  int ierr;
  int nlev = nlyr + 1;

  const double inv_dx_dy = 1.0 / (dx * dy);

  // Boundary Conditions at iz == nlev - 1, ix == 0 and iy == 0
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      E_direct[nlev - 1][ix][iy][0] = S0 * mu0;
    }
  }

  /* Determine loop indices depending on sun direction in the solar spectral
       range and whether you are in the thermal spectral range 
       */
  unsigned int ix_start, ix_n_it, ix_increment;
  unsigned int iy_start, iy_n_it, iy_increment;
  unsigned int iz_start, iz_n_it, iz_increment;

  ierr = determine_loop_indices (nit,
                                 phi0,
                                 planck,
                                 nx,
                                 ny,
                                 nlyr,
                                 &ix_start,
                                 &ix_n_it,
                                 &ix_increment,
                                 &iy_start,
                                 &iy_n_it,
                                 &iy_increment,
                                 &iz_start,
                                 &iz_n_it,
                                 &iz_increment);
  CHKERR (ierr);

  // ONE GAUSS-SEIDEL ITERATION FORWARD in the dynamic tenstream algorithm
  for (unsigned int iz = iz_start; iz < iz_n_it; iz += iz_increment) {

    /* z-index in 1D, since 1D profiles are sorted from top to bottom,
           while 3D fields are sorted from bottom to top */
    int    iz1D = nlyr - 1 - iz;
    double dz   = (zd[iz1D] - zd[iz1D + 1]) * 1.0E3;

    const double inv_dx_dz = 1.0 / (dx * dz);
    const double inv_dy_dz = 1.0 / (dy * dz);

    for (unsigned int ix = ix_start; ix < ix_n_it; ix += ix_increment) {
      for (unsigned int iy = iy_start; iy < iy_n_it; iy += iy_increment) {

        /* Dynamic Radiative Transfer: For every grid box, new output irradiances 
                   are calculated using the input irradiances of the previous timestep */
        float E_in[N_diffuse + N_direct];
        ierr = determine_E_in (ix, iy, iz, nx, ny, dx, dy, dz, phi0, E_direct, E_diffuse, E_in);
        CHKERR (ierr);

        float E_out[N_diffuse + N_direct];
        for (int i = 0; i < (N_diffuse + N_direct); i++) {
          E_out[i] = 0;
        }

        for (int idiff = 0; idiff < N_diffuse; idiff++) {
          for (int jdiff = 0; jdiff < N_diffuse; jdiff++) {
            E_out[idiff] += a[iz][ix][iy][idiff][jdiff] * E_in[jdiff];
          }
          for (int jdir = 0; jdir < N_direct; jdir++) {
            // if (ix == 16 && iy == 16 && iz == 125) {
            E_out[idiff] += b[iz][ix][iy][idiff][jdir] * E_in[N_diffuse + jdir];
          }
          if (idiff == 0 || idiff == 1) {
            E_out[idiff] += e[iz][ix][iy][idiff] * B_eff[iz][ix][iy][idiff] * dx * dy;
          } else if (idiff >= 2 && idiff <= 5) {
            E_out[idiff] += e[iz][ix][iy][idiff] * B_eff[iz][ix][iy][idiff] * dy * dz;
          } else {
            E_out[idiff] += e[iz][ix][iy][idiff] * B_eff[iz][ix][iy][idiff] * dx * dz;
          }
        }
        for (int idir = 0; idir < N_direct; idir++) {
          for (int jdir = 0; jdir < N_direct; jdir++) {
            E_out[N_diffuse + idir] += c[iz][ix][iy][idir][jdir] * E_in[N_diffuse + jdir];
          }
        }

        E_diffuse[iz + 1][ix % nx][iy % ny][0]   = E_out[0] * inv_dx_dy;
        E_diffuse[iz][ix % nx][iy % ny][1]       = E_out[1] * inv_dx_dy;
        E_diffuse[iz][ix % nx][iy % ny][2]       = E_out[2] * inv_dy_dz;
        E_diffuse[iz][(ix + 1) % nx][iy % ny][3] = E_out[3] * inv_dy_dz;
        E_diffuse[iz][ix % nx][iy % ny][4]       = E_out[4] * inv_dy_dz;
        E_diffuse[iz][(ix + 1) % nx][iy % ny][5] = E_out[5] * inv_dy_dz;
        E_diffuse[iz][ix % nx][iy % ny][6]       = E_out[6] * inv_dx_dz;
        E_diffuse[iz][ix % nx][(iy + 1) % ny][7] = E_out[7] * inv_dx_dz;
        E_diffuse[iz][ix % nx][iy % ny][8]       = E_out[8] * inv_dx_dz;
        E_diffuse[iz][ix % nx][(iy + 1) % ny][9] = E_out[9] * inv_dx_dz;

        E_direct[iz][ix % nx][iy % ny][0] = E_out[10] * inv_dx_dy;
        if (phi0 <= 90.0) {
          E_direct[iz][(ix + 1) % nx][iy % ny][1] = E_out[11] * inv_dy_dz;
          E_direct[iz][ix % nx][(iy + 1) % ny][2] = E_out[12] * inv_dx_dz;
        } else if (phi0 > 90.0 && phi0 <= 180.0) {
          E_direct[iz][(ix + 1) % nx][iy % ny][1] = E_out[11] * inv_dy_dz;
          E_direct[iz][ix % nx][iy % ny][2]       = E_out[12] * inv_dx_dz;
        } else if (phi0 > 180.0 && phi0 <= 270.0) {
          E_direct[iz][ix % nx][iy % ny][1] = E_out[11] * inv_dy_dz;
          E_direct[iz][ix % nx][iy % ny][2] = E_out[12] * inv_dx_dz;
        } else {
          E_direct[iz][ix % nx][iy % ny][1]       = E_out[11] * inv_dy_dz;
          E_direct[iz][ix % nx][(iy + 1) % ny][2] = E_out[12] * inv_dx_dz;
        }
      }
    }
  }

  // Boundary Conditions II: iz == 0
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      E_diffuse[0][ix][iy][0] = Ag * (E_diffuse[0][ix][iy][1] + E_direct[0][ix][iy][0]) + (1.0 - Ag) * M_PI * Bg;
    }
  }
  return 0;
}

/* function calculating irradiances */
int calculate_irradiances (int                n_iterations,
                           int                nx,
                           int                ny,
                           int                nlyr,
                           int                N_direct,
                           int                N_diffuse,
                           int                delta,
                           int                planck,
                           double             S0,
                           double             mu0,
                           double             phi0,
                           double             dx,
                           double             dy,
                           float*             zd,
                           float*             temper,
                           float              wvnmlo,
                           float              wvnmhi,
                           float*             dtau_org,
                           float*             omega0_org,
                           float*             g_org,
                           int                found3D_ic,
                           int                found3D_wc,
                           caoth3d_out_struct ic3D,
                           caoth3d_out_struct wc3D,
                           char*              datapath,
                           int                quiet,
                           float              Ag,
                           float              Bg,
                           float****          B_eff,
                           float*****         a,
                           float*****         b,
                           float*****         c,
                           float****          e,
                           float****          E_direct,
                           float****          E_diffuse) {

  int ierr;

  // Setup of Planck Emissivity and TenStream Coefficients
  ierr = determine_tenstream_coefficients_and_effective_B (N_direct,
                                                           N_diffuse,
                                                           nlyr,
                                                           nx,
                                                           ny,
                                                           mu0,
                                                           phi0,
                                                           delta,
                                                           planck,
                                                           zd,
                                                           temper,
                                                           wvnmlo,
                                                           wvnmhi,
                                                           dtau_org,
                                                           omega0_org,
                                                           g_org,
                                                           found3D_ic,
                                                           found3D_wc,
                                                           ic3D,
                                                           wc3D,
                                                           datapath,
                                                           quiet,
                                                           a,
                                                           b,
                                                           c,
                                                           e,
                                                           B_eff);
  CHKERR (ierr);

  // Memory allocation of arrays needed inside the dynamic tenstream loops

  /* Updated Radiative Fields calculated in time step i */
  float**** E_direct_new  = calloc_float_4D (nlyr + 1, nx, ny, N_direct, "E_direct_new");
  float**** E_diffuse_new = calloc_float_4D (nlyr + 1, nx, ny, N_diffuse, "E_diffuse_new");

  for (int i = 0; i < n_iterations; i++) {

    /*
        ierr = one_jacobian_iteration(nx,
                                      ny,
                                      nlyr,
                                      N_direct,
                                      N_diffuse,
                                      S0,
                                      mu0,
                                      phi0,
                                      Ag,
                                      cld3D,
                                      Bg,
                                      B_eff,
                                      zd,
                                      a,
                                      b,
                                      c,
                                      e,
                                      E_direct_new,
                                      E_diffuse_new,
                                      E_direct,
                                      E_diffuse);
         CHKERR(ierr);
         */

    ierr = one_gauss_seidel_iteration (i,
                                       nx,
                                       ny,
                                       nlyr,
                                       N_direct,
                                       N_diffuse,
                                       planck,
                                       S0,
                                       mu0,
                                       phi0,
                                       Ag,
                                       Bg,
                                       B_eff,
                                       dx,
                                       dy,
                                       zd,
                                       a,
                                       b,
                                       c,
                                       e,
                                       E_direct,
                                       E_diffuse);
    CHKERR (ierr);
  }

  free_float_4D (E_direct_new);
  free_float_4D (E_diffuse_new);

  return 0;
}

/* function calculating heating rates */
int calculate_heating_rates (int        nx,
                             int        ny,
                             int        nlyr,
                             int        N_direct,
                             int        N_diffuse,
                             double     phi0,
                             float****  B_eff,
                             double     dx,
                             double     dy,
                             float*     zd,
                             float***** a,
                             float***** b,
                             float***** c,
                             float****  e,
                             float****  E_direct,
                             float****  E_diffuse,
                             float***   divE) {

  int ierr;

  for (int iz = 0; iz < nlyr; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {

        int iz1D = nlyr - 1 - iz;

        divE[iz][ix][iy] = 0.0;

        double dz = (zd[iz1D] - zd[iz1D + 1]) * 1.0E3;

        float E_in[N_diffuse + N_direct];
        ierr = determine_E_in (ix, iy, iz, nx, ny, dx, dy, dz, phi0, E_direct, E_diffuse, E_in);
        CHKERR (ierr);

        double divE_z = 0.0, divE_x = 0.0, divE_y = 0.0;

        divE_z += E_in[10];
        for (int idiff = 0; idiff <= 1; idiff++) {
          divE_z += E_in[idiff] - e[iz][ix][iy][idiff] * B_eff[iz][ix][iy][idiff] * dx * dy;
          for (int jdiff = 0; jdiff < N_diffuse; jdiff++) {
            divE_z -= a[iz][ix][iy][idiff][jdiff] * E_in[jdiff];
          }
        }
        for (int idir = 0; idir < N_direct; idir++) {
          divE_z -= c[iz][ix][iy][0][idir] * E_in[N_diffuse + idir];
          for (int jdiff = 0; jdiff <= 1; jdiff++) {
            divE_z -= b[iz][ix][iy][jdiff][idir] * E_in[N_diffuse + idir];
          }
        }

        divE_x += E_in[11];
        for (int idiff = 2; idiff <= 5; idiff++) {
          divE_x += E_in[idiff] - e[iz][ix][iy][idiff] * B_eff[iz][ix][iy][idiff] * dy * dz;
          for (int jdiff = 0; jdiff < N_diffuse; jdiff++) {
            divE_x -= a[iz][ix][iy][idiff][jdiff] * E_in[jdiff];
          }
        }
        for (int idir = 0; idir < N_direct; idir++) {
          divE_x -= c[iz][ix][iy][1][idir] * E_in[N_diffuse + idir];
          for (int jdiff = 2; jdiff <= 5; jdiff++) {
            divE_x -= b[iz][ix][iy][jdiff][idir] * E_in[N_diffuse + idir];
          }
        }

        divE_y += E_in[12];
        for (int idiff = 6; idiff <= 9; idiff++) {
          divE_y += E_in[idiff] - e[iz][ix][iy][idiff] * B_eff[iz][ix][iy][idiff] * dx * dz;
          for (int jdiff = 0; jdiff < N_diffuse; jdiff++) {
            divE_y -= a[iz][ix][iy][idiff][jdiff] * E_in[jdiff];
          }
        }
        for (int idir = 0; idir < N_direct; idir++) {
          divE_y -= c[iz][ix][iy][2][idir] * E_in[N_diffuse + idir];
          for (int jdiff = 6; jdiff <= 9; jdiff++) {
            divE_y -= b[iz][ix][iy][jdiff][idir] * E_in[N_diffuse + idir];
          }
        }

        divE[iz][ix][iy] = (divE_z + divE_x + divE_y) / (dx * dy);
      }
    }
  }

  return 0;
}

int clearsky_spinup (int                spinup_iterations,
                     int                nx,
                     int                ny,
                     int                nlyr,
                     int                N_direct,
                     int                N_diffuse,
                     int                delta,
                     int                planck,
                     double             S0,
                     double             mu0,
                     double             phi0,
                     double             dx,
                     double             dy,
                     float*             zd,
                     float*             temper,
                     float              wvnmlo,
                     float              wvnmhi,
                     float*             dtau_org,
                     float*             omega0_org,
                     float*             g_org,
                     caoth3d_out_struct ic3D,
                     caoth3d_out_struct wc3D,
                     char*              datapath,
                     int                quiet,
                     float              Ag,
                     float              Bg,
                     float****          E_direct,
                     float****          E_diffuse) {

  int ierr;

  /* Define a 1D clearsky column for the spin-up */
  int nx_1D      = 1;
  int ny_1D      = 1;
  int found3D_wc = 0;
  int found3D_ic = 0;

  float**** B_eff_1D     = calloc_float_4D (nlyr, nx_1D, ny_1D, N_diffuse, "B_eff");
  float**** E_direct_1D  = calloc_float_4D (nlyr + 1, nx_1D, ny_1D, N_direct, "E_direct");
  float**** E_diffuse_1D = calloc_float_4D (nlyr + 1, nx_1D, ny_1D, N_diffuse, "E_diffuse");

  float***** a_1D = calloc_float_5D (nlyr, nx_1D, ny_1D, N_diffuse, N_diffuse, "a");
  float***** b_1D = calloc_float_5D (nlyr, nx_1D, ny_1D, N_diffuse, N_direct, "b");
  float***** c_1D = calloc_float_5D (nlyr, nx_1D, ny_1D, N_direct, N_direct, "c");
  float****  e_1D = calloc_float_4D (nlyr, nx_1D, ny_1D, N_diffuse, "e");

  ierr = calculate_irradiances (spinup_iterations,
                                nx_1D,
                                ny_1D,
                                nlyr,
                                N_direct,
                                N_diffuse,
                                delta,
                                planck,
                                S0,
                                mu0,
                                phi0,
                                dx,
                                dy,
                                zd,
                                temper,
                                wvnmlo,
                                wvnmhi,
                                dtau_org,
                                omega0_org,
                                g_org,
                                found3D_ic,
                                found3D_wc,
                                ic3D,
                                wc3D,
                                datapath,
                                quiet,
                                Ag,
                                Bg,
                                B_eff_1D,
                                a_1D,
                                b_1D,
                                c_1D,
                                e_1D,
                                E_direct_1D,
                                E_diffuse_1D);
  CHKERR (ierr);

  // Copy 1D clearsky results to 3D fields
  for (int iz = 0; iz < nlyr + 1; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
        for (int idir = 0; idir < N_direct; idir++) {
          E_direct[iz][ix][iy][idir] = E_direct_1D[iz][0][0][idir];
        }
        for (int idiff = 0; idiff < N_diffuse; idiff++) {
          E_diffuse[iz][ix][iy][idiff] = E_diffuse_1D[iz][0][0][idiff];
        }
      }
    }
  }

  return 0;
}

int dynamic_tenstream (int                 n_iterations,
                       int                 history,
                       caoth3d_out_struct* caoth3d,
                       int                 n_caoth,
                       float*              dtau_org,
                       float*              omega0_org,
                       float*              g_org,
                       int                 nlev,
                       double              S0,
                       double              mu0,
                       double              phi0,
                       double              Ag,
                       int                 planck,
                       int                 delta,
                       int                 nzout,
                       float*              zd,
                       float*              temper,
                       float               btemp,
                       float               wvnmlo,
                       float               wvnmhi,
                       float*              zout,
                       int                 iv,
                       int                 ib,
                       char*               datapath,
                       int                 quiet,
                       float*              fldn,
                       float*              flup,
                       float*              fldir,
                       float*              uavg,
                       float***            rfldir3d,
                       float***            rfldn3d,
                       float***            flup3d,
                       float***            abs3d) {

  int                ierr;
  int                found3D_wc = 0, found3D_ic = 0;
  caoth3d_out_struct wc3D, ic3D;
  int                nlyr = nlev - 1;

  if (!quiet) {
    fprintf (stderr, " ... *** Richard Maier's dynamic tenstream solver with %d iterations\n", n_iterations);
    fprintf (stderr, " ... *** history = %d\n", history);
  }

  /* output quantities at user levels */
  for (int lu = 0; lu < nzout; lu++) {
    fldn[lu]  = 0.0;
    flup[lu]  = 0.0;
    fldir[lu] = 0.0;
    uavg[lu]  = 0.0;
  }

  for (int isp = MCCAOTH_FIR; isp <= n_caoth; isp++) {
    int ispo = isp - MCCAOTH_FIR;

    /* transfer indices of water and ice from coath3d to atmos structure
         for subsequent lookups since caoth arrangement can be arbitrary */
    if (strncasecmp (caoth3d[ispo].name, "wc", 2) == 0) {

      found3D_wc = 1;

      /* only 3D clouds with more than zero 3D layers are considered */
      if (caoth3d[ispo].nlyr == 0 || caoth3d[ispo].nthreed == 0) {
        found3D_wc = 0;
      } else { /* assign water cloud caoth to wc3D */
        wc3D = caoth3d[ispo];
      }
    }
    if (strncasecmp (caoth3d[ispo].name, "ic", 2) == 0) {

      found3D_ic = 1;

      /* only 3D clouds with more than zero 3D layers are considered */
      if (caoth3d[ispo].nlyr == 0 || caoth3d[ispo].nthreed == 0) {
        found3D_ic = 0;
      } else { /* assign ice cloud caoth to ic3D */
        ic3D = caoth3d[ispo];
      }
    }
  }

  /* if 3D cloud was found, found3D is set to 1 and the 3D cloud is accessible in caoth3d_out_struct cld3D */

  if (!quiet) {
    if (found3D_wc == 0 && found3D_ic == 0) {
      fprintf (stderr, " ... no 3D cloud found\n");
    } else {
      if (found3D_wc == 1) {
        fprintf (stderr, " ... found 3D water cloud %d x %d with %d 3D layers\n", wc3D.Nx, wc3D.Ny, wc3D.nthreed);
      }
      if (found3D_ic == 1) {
        fprintf (stderr, " ... found 3D ice cloud %d x %d with %d 3D layers\n", ic3D.Nx, ic3D.Ny, ic3D.nthreed);
      }
    }
  }

  if (found3D_wc == 1 || found3D_ic == 1) {

    /* check whether dimenesion of atmopshere and cloud fields are consistent */
    if (found3D_wc) {
      if (wc3D.nlyr + 1 != nlev) {
        fprintf (stderr,
                 "Fatal error, different number of layers in atmosphere (%d) and 3D water cloud profile (%d)\n",
                 nlev - 1,
                 wc3D.nlyr);
        return -1;
      }
    }
    if (found3D_ic) {
      if (ic3D.nlyr + 1 != nlev) {
        fprintf (stderr,
                 "Fatal error, different number of layers in atmosphere (%d) and 3D ice cloud profile (%d)\n",
                 nlev - 1,
                 ic3D.nlyr);
        return -1;
      }
    }
    if (found3D_wc == 1 && found3D_ic == 1) {
      if (wc3D.Nx != ic3D.Nx || wc3D.Ny != ic3D.Ny) {
        fprintf (stderr,
                 "Fatal error, different horizontal grid shape for 3D water (%d x %d) and ice cloud (%d x %d) domains \n",
                 wc3D.Nx,
                 wc3D.Ny,
                 ic3D.Nx,
                 ic3D.Ny);
        return -1;
      }
    }

    /* create an array that combines wc3D.threed and ic3D.threed */
    int* threed = calloc (nlev, sizeof (int));
    if (found3D_wc) {
      for (int ilyr = 0; ilyr < nlyr; ilyr++) {
        threed[ilyr] = wc3D.threed[ilyr];
      }
    }
    if (found3D_ic) {
      for (int ilyr = 0; ilyr < nlyr; ilyr++) {
        if (ic3D.threed[ilyr]) {
          threed[ilyr] = ic3D.threed[ilyr];
        }
      }
    }

    // TenStream Setup: Number of diffuse and direct streams
    int N_diffuse = 10;
    int N_direct  = 3;

    // Horizontal domain size
    int nx, ny;
    if (found3D_wc) {
      nx = wc3D.Nx, ny = wc3D.Ny;
    } else {
      nx = ic3D.Nx, ny = ic3D.Ny;
    }

    // Horizontal domain spacing
    double dx, dy;
    if (found3D_wc) {
      dx = wc3D.delX, dy = wc3D.delY;
    } else {
      dx = ic3D.delX, dy = ic3D.delY;
    }

    // Setup of Radiative Quantities
    /* allocate these 4D arrays as 1D arrays using the calloc functions in allocnd.c in order to be able
         to read and write them easily to files */
    float**** E_direct  = calloc_float_4D (nlev, nx, ny, N_direct, "E_direct");
    float**** E_diffuse = calloc_float_4D (nlev, nx, ny, N_diffuse, "E_diffuse");

    /* if available, read in previously calculated E_direct and E_diffuse fields for the current subband */
    int success = 0;
    if (history) {
      ierr = read_initial_irradiances (iv, ib, nlev, nx, ny, N_direct, N_diffuse, &success, E_direct, E_diffuse);
      CHKERR (ierr);
    }

    // Planck Setup
    float**** B_eff = calloc_float_4D (nlyr, nx, ny, N_diffuse, "B_eff");
    float     Bg    = 0.0;  // black body radiance by the ground within the given wavelength spectrum
    if (planck) {
      Bg = c_planck_func1 (wvnmlo, wvnmhi, btemp);
    }

    // Setup of Transport Coefficient Quantities
    /* TenStream coefficients for the whole domain and quantities needed for the retrieval of these coefficients from the LUT's */
    float***** a =
      calloc_float_5D (nlyr, nx, ny, N_diffuse, N_diffuse, "a");  // TenStream coefficients for diffuse to diffuse streams
    float***** b =
      calloc_float_5D (nlyr, nx, ny, N_diffuse, N_direct, "b");              // TenStream coefficients for diffuse to direct streams
    float***** c = calloc_float_5D (nlyr, nx, ny, N_direct, N_direct, "c");  // TenStream coefficients for direct to direct streams
    float****  e = calloc_float_4D (nlyr, nx, ny, N_diffuse, "e");           // thermal emissivities

    // Calculation of Irradiances
    if (!quiet) {
      fprintf (stderr, "CURRENT WAVELENGTH BAND: IV = %i, IB = %i \n", iv, ib);
    }

    /* if no previously calculated E_direct and E_diffuse fields could be retrieved, we start with a 
           clearsky spin-up */
    if (!success) {

      int spinup_iterations = 2000;
      if (!quiet) {
        fprintf (stderr, "Performing clearsky spin-up with %5i spin-up iterations. \n", spinup_iterations);
      }

      ierr = clearsky_spinup (spinup_iterations,
                              nx,
                              ny,
                              nlyr,
                              N_direct,
                              N_diffuse,
                              delta,
                              planck,
                              S0,
                              mu0,
                              phi0,
                              dx,
                              dy,
                              zd,
                              temper,
                              wvnmlo,
                              wvnmhi,
                              dtau_org,
                              omega0_org,
                              g_org,
                              ic3D,
                              wc3D,
                              datapath,
                              quiet,
                              Ag,
                              Bg,
                              E_direct,
                              E_diffuse);
      CHKERR (ierr);
    }

    if (!quiet) {
      fprintf (stderr, "Performing %5i full dynamic TenStream iterations. \n", n_iterations);
    }
    ierr = calculate_irradiances (n_iterations,
                                  nx,
                                  ny,
                                  nlyr,
                                  N_direct,
                                  N_diffuse,
                                  delta,
                                  planck,
                                  S0,
                                  mu0,
                                  phi0,
                                  dx,
                                  dy,
                                  zd,
                                  temper,
                                  wvnmlo,
                                  wvnmhi,
                                  dtau_org,
                                  omega0_org,
                                  g_org,
                                  found3D_ic,
                                  found3D_wc,
                                  ic3D,
                                  wc3D,
                                  datapath,
                                  quiet,
                                  Ag,
                                  Bg,
                                  B_eff,
                                  a,
                                  b,
                                  c,
                                  e,
                                  E_direct,
                                  E_diffuse);
    CHKERR (ierr);

    // Calculation of Heating Rates
    float*** divE = calloc_float_3D (nlyr, nx, ny, "divE");

    ierr =
      calculate_heating_rates (nx, ny, nlyr, N_direct, N_diffuse, phi0, B_eff, dx, dy, zd, a, b, c, e, E_direct, E_diffuse, divE);
    CHKERR (ierr);

    // Dump direct and diffuse irradiance fields for new dynamic_tenstream runs starting from these fields
    if (history) {
      ierr = dump_irradiances (iv, ib, nlev, nx, ny, N_direct, N_diffuse, E_direct, E_diffuse);
      CHKERR (ierr);
    }

    // Copy Results to Final Fields
    for (int iz = 0; iz < nlev; iz++) {
      for (int jz = 0; jz < nzout; jz++) {
        if (zout[jz] == zd[iz]) {
          int iz3D = nlev - 1 - iz;
          for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
              rfldir3d[jz][ix][iy] = E_direct[iz3D][ix][iy][0];
              rfldn3d[jz][ix][iy]  = E_diffuse[iz3D][ix][iy][1];
              flup3d[jz][ix][iy]   = E_diffuse[iz3D][ix][iy][0];

              fldir[jz] += E_direct[iz3D][ix][iy][0];
              fldn[jz] += E_diffuse[iz3D][ix][iy][1];
              flup[jz] += E_diffuse[iz3D][ix][iy][0];
            }
          }
          fldir[jz] /= (nx * ny);
          fldn[jz] /= (nx * ny);
          flup[jz] /= (nx * ny);
          uavg[jz] = NAN;
        }
      }
    }
    for (int ilyr = 0; ilyr < nlyr; ilyr++) {
      if (threed[ilyr]) {
        for (int ix = 0; ix < nx; ix++) {
          for (int iy = 0; iy < ny; iy++) {
            abs3d[ilyr][ix][iy] = divE[ilyr][ix][iy];
          }
        }
      }
    }

    // Free allocated arrays
    free (threed);

    free_float_4D (E_diffuse);
    free_float_4D (E_direct);

    free_float_3D (divE);

    free_float_4D (B_eff);

    free_float_5D (a);
    free_float_5D (b);
    free_float_5D (c);
    free_float_4D (e);

    /* !!! atmospheric profile top to bottom, 3D caoth bottom to top !!! */
    //    for (int ilev=0; ilev<nlev; ilev++)
    //      fprintf (stderr, "%d %f %f\n", ilev, cld3D.zd[nlev-1-ilev], zd[ilev]);

    /*
        fprintf (stderr, "# z[km]    tau_clr    LWC reff  tau_cld   omega0     asym\n");

        for (int ilyr = 0; ilyr < nlyr; ilyr++) {
            fprintf (stderr, "%7.2f %10.6f  ", zd[ilyr + 1], dtau_org[ilyr]);
            int icld = nlyr - 1 - ilyr;
            if (wc3D.threed[icld]) {
                double tau = cld3D.ext[icld][0][0] * (zd[ilyr] - zd[ilyr + 1]) * 1000.0;
                fprintf (stderr,
                         "%5.3f %4.1f  %7.2f %.6f %.6f\n",
                         wc3D.lwc[icld][0][0],
                         wc3D.reff[icld][0][0],
                         tau,
                         wc3D.ssa[icld][0][0],
                         wc3D.g1[icld][0][0]);
            }
            else {
                fprintf (stderr, "%5.3f %4.1f  %7.2f %.6f %.6f\n", 0.0, 0.0, 0.0, 0.0, 0.0);
            }
        }
        
        */

    /*
        // example fluxes at zout's
        for (int iz = 0; iz < nzout; iz++) {
            rfldir3d[iz][0][0]                       = (float)(iz + 1) / 2.0;
            rfldir3d[iz][wc3D.Nx - 1][wc3D.Ny - 1] = (float)(iz + 1) / 2.0;
        }
        // example absorption in all 3d model layers
        for (int ilyr = 0; ilyr < nlyr; ilyr++) {
            if (wc3D.threed[ilyr]) {
                abs3d[ilyr][0][0] = 2.0;
                //  abs3d[ilyr][wc3D.Nx-1][wc3D.Ny-1] = 1.0;
            }
        }
        */
  }

  return 0;
}
