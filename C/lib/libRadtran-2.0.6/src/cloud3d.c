/*--------------------------------------------------------------------
 * $Id$
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
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
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "uvspec.h"
#include "ckdfu.h"
#include "cloud3d.h"
#include "ipa.h"
#include "wcloud3d.h"
#include "elevation2d.h"
#include "temperature2d.h"
#include "common_math.h"
#include "f77-uscore.h"
#include "ascii.h"
#include "errors.h"

#if HAVE_VROOM
#include "vroom.h"
#endif
#if HAVE_LIDAR
#include "lidar.h"
#endif
#if HAVE_TIPA
#include "tipa.h"
#endif

/* if uncommented, calculate area for elevation numerically in */
/* addition to the analytical calculation; for debugging only; */
/* output is then additionally written to area2D.dat           */

/* #define AREA_DEBUG_N 100000 */

static inline int float_equal (float a, float b);

/* prototypes of internal functions */
/* static int compare_levels_df (double *z1, int n1, float  *z2, int n2); */

#if HAVE_MYSTIC
static int setup_sample_grid (float*         zd,
                              int            nlev,
                              int            nxout,
                              int            nyout,
                              double         dxout,
                              double         dyout,
                              int            nx_caoth,
                              int            ny_caoth,
                              double         dx_caoth,
                              double         dy_caoth,
                              float*         zout,
                              int            nzout,
                              int            radpat_Nr,
                              int            radpat_Nt,
                              float          radpat_dt,
                              float*         umu,
                              int            numu,
                              float          alpha,
                              float*         phi,
                              int            nphi,
                              int            aziconv,
                              sample_struct* sample,
                              int            escape,
                              int            vroom,
                              int            calcstd,
                              double         std_value,
                              int            backward,
                              int            output,
                              int            backward_islower,
                              int            backward_jslower,
                              int            backward_isupper,
                              int            backward_jsupper,
                              int            backward_isstep,
                              int            backward_jsstep,
                              int            backward_writeallpixels,
                              int            backward_writeback,
                              int            refraction,
                              int            spectral_is,
                              int            boxairmass,
                              int            boxairmass3D,
                              int            spherical,
                              int            spherical3D,
                              int            spherical3D_scene,
                              double         spherical3D_scene_lon_min,
                              double         spherical3D_scene_lon_max,
                              double         spherical3D_scene_lat_min,
                              double         spherical3D_scene_lat_max,
                              int            DoLE,
                              int            ixmin,
                              int            ixmax,
                              int            iymin,
                              int            iymax,
                              int            ipa3d,
                              int            tipa,
                              int            polarisation,
                              int            polarisation_state,
                              int            nstokes,
                              int            bcond,
                              int            panorama,
                              int            panorama_forward,
                              char*          panorama_filename,
                              char*          backward_ray_filename,
                              int            maxscatters,
                              int            minscatters,
                              int            ncirc,
                              double         dcirc,
                              double*        blitz_position,
                              int            locest,
                              int            jacobian,
                              int            jacobian3D,
                              int            jacobian_std,
                              int            delta_scaling_start,
                              char*          lidarfilename,
                              int            reference_to_NN,
                              int            coherent_backscatter,
                              int            coherent_backscatter_lmode,
                              int            coherent_backscatter_off,
                              double         ris_factor,
                              double         ris_optical_depth,
                              int            sample_cldprp,
                              int            quiet);
#endif

static int read_caoth3D_file (char*            filename,
                              input_struct     input,
                              caoth_inp_struct input_caoth,
                              output_struct*   output,
                              /* Output */
                              caoth_out_struct*   output_caoth,
                              caoth3d_out_struct* caoth3d,
                              int*                nmom);

static float interpolate_f (phase_function_table* p, int i1, int i2, double delta_r);
static float interpolate_dscale (phase_function_table* p, int i1, int i2, double delta_r);

double calc_elevation_area (elevation_struct* elev, double x1, double y1, double x2, double y2);

/***********************************************************************************/
/* Function: equal                                                        @62_30i@ */
/* Description:                                                                    */
/*  Two numbers equal?                                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

#if HAVE_MYSTIC
static int equal (double x, double y) {
  if ((float) x != (float) y)
    return 0;

  return 1;
}
#endif

/***********************************************************************************/
/* Function: setup_sample_grid                                            @62_30i@ */
/* Description:                                                                    */
/*  Initialize sample_struct.                                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/*  float *zout:  array of altitudes to be sampled                                 */
/*  int   nzout:  number of altitudes to be sampled                                */
/*                                                                                 */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

#if HAVE_MYSTIC
static int setup_sample_grid (float*         zd,
                              int            nlev,
                              int            nxout,
                              int            nyout,
                              double         dxout,
                              double         dyout,
                              int            nx_caoth,
                              int            ny_caoth,
                              double         dx_caoth,
                              double         dy_caoth,
                              float*         zout,
                              int            nzout,
                              int            radpat_Nr,
                              int            radpat_Nt,
                              float          radpat_dt,
                              float*         umu,
                              int            numu,
                              float          alpha,
                              float*         phi,
                              int            nphi,
                              int            aziconv,
                              sample_struct* sample,
                              int            escape,
                              int            vroom,
                              int            calcstd,
                              double         std_value,
                              int            backward,
                              int            output,
                              int            backward_islower,
                              int            backward_jslower,
                              int            backward_isupper,
                              int            backward_jsupper,
                              int            backward_isstep,
                              int            backward_jsstep,
                              int            backward_writeallpixels,
                              int            backward_writeback,
                              int            refraction,
                              int            spectral_is,
                              int            boxairmass,
                              int            boxairmass3D,
                              int            spherical,
                              int            spherical3D,
                              int            spherical3D_scene,
                              double         spherical3D_scene_lon_min,
                              double         spherical3D_scene_lon_max,
                              double         spherical3D_scene_lat_min,
                              double         spherical3D_scene_lat_max,
                              int            DoLE,
                              int            ixmin,
                              int            ixmax,
                              int            iymin,
                              int            iymax,
                              int            ipa3d,
                              int            tipa,
                              int            polarisation,
                              int            polarisation_state,
                              int            nstokes,
                              int            bcond,
                              int            panorama,
                              int            panorama_forward,
                              char*          panorama_filename,
                              char*          backward_ray_filename,
                              int            maxscatters,
                              int            minscatters,
                              int            ncirc,
                              double         dcirc,
                              double*        blitz_position,
                              int            locest,
                              int            jacobian,
                              int            jacobian3D,
                              int            jacobian_std,
                              int            delta_scaling_start,
                              char*          lidarfilename,
                              int            reference_to_NN,
                              int            coherent_backscatter,
                              int            coherent_backscatter_lmode,
                              int            coherent_backscatter_off,
                              double         ris_factor,
                              double         ris_optical_depth,
                              int            sample_cldprp,
                              int            quiet) {
  int     status = 0;
  int     lu = 0, iu = 0, j = 0, id = 0, kc = 0;
  double* z_hist = NULL;
  int     Nd     = 0;

  double xr = 0, yr = 0;

  double *stheta = NULL, *sphi = NULL, *salpha = NULL;

  FILE* ftmp = NULL;

#if HAVE_LIDAR
  double cosd_temp;
  int    ili = 0, i = 0, io = 0, ia = 0;
  double alphaact = 0.0, scr = 0.0;
#endif

  if (!quiet)
    fprintf (stderr, " ... reading sampling information\n");

  /* test if a sample2D.dat is in the current directory */
  if ((ftmp = fopen ("./sample2D.dat", "r")) != NULL) {

    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!! ATTENTION: Found an old sample file sample2D.dat. Before October 25, 2008 !!!\n");
    fprintf (stderr, "!!! MYSTIC got information about sample grid, radiance directions, etc. from  !!!\n");
    fprintf (stderr, "!!! this file. Since October 25, 2008 it is ignored because all the           !!!\n");
    fprintf (stderr, "!!! information can now be directly specified in the uvspec input file.       !!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

    fclose (ftmp);
  }

  if (panorama) {
    escape = 1;
    numu   = 1;
    nphi   = 1;
  }

  if (panorama_forward) {
    escape = 0;
    numu   = 1; /* ??? check ??? */
    nphi   = 1; /* ??? check ??? */
  }

  /* read angles from mc_panorama_file; we need it here because it determines sample->Nx and sample->Ny */
  if (strlen (panorama_filename) > 0) {
    status =
      read_2c_file (panorama_filename, &(sample->pan_phi_from_file), &(sample->pan_theta_from_file), &(sample->pan_n_from_file));
    if (status != 0) {
      fprintf (stderr, "Error reading %s\n", panorama_filename);
      CHKERR (-1);
    }
    // check if theta between 0 and 180; don't care about phi
    for (int itheta = 0; itheta < sample->pan_n_from_file; itheta++)
      if (sample->pan_theta_from_file[itheta] < 0 || sample->pan_theta_from_file[itheta] > 180.0) {
        fprintf (stderr, "Error, theta smaller than 0 or larger than 180 degree\n");
        CHKERR (-1);
      }

    switch (aziconv) {
    case MCAZIMUTH_NEW:
      for (int iphi = 0; iphi < sample->pan_n_from_file; iphi++)
        sample->pan_phi_from_file[iphi] += 180.0;
      break;
    case MCAZIMUTH_OLD:
      break;
    default:
      fprintf (stderr, "Error, unknown azimuth convention %d\n", aziconv);
      CHKERR (-1);
    }

    /* Nx = 1; Ny = pan_n_from_file */
    nxout = 1;
    nyout = sample->pan_n_from_file;
  }

  /* read angles from mc_panorama_ray_file; we need it here because it determines sample->Nx and sample->Ny */
  if (strlen (backward_ray_filename) > 0) {
    int mincol, maxcol;
    status = ASCII_file2double_nrows (backward_ray_filename,
                                      -1,
                                      &(sample->backward_n_from_ray_file),
                                      &mincol,
                                      &maxcol,
                                      &(sample->backward_ray_from_file));
    CHKERR (status);
    CHKERROUT (mincol - 6, "expected 6 columns in panorama ray file");
    CHKERROUT (maxcol - 6, "expected 6 columns in panorama ray file");

    for (size_t i = 0; i < sample->backward_n_from_ray_file; ++i) {
      if (!quiet) {
        if (i < 20 || i > sample->backward_n_from_ray_file - 20) {
          fprintf (stderr,
                   "Backward ray %lu origin [ %8.4g %8.4g %8.4g ] direction [ %8.4g %8.4g %8.4g ]\n",
                   i,
                   sample->backward_ray_from_file[i][0],
                   sample->backward_ray_from_file[i][1],
                   sample->backward_ray_from_file[i][2],
                   sample->backward_ray_from_file[i][3],
                   sample->backward_ray_from_file[i][4],
                   sample->backward_ray_from_file[i][5]);
        } else if (i == 20) {
          fprintf (stderr, "[...]\n");
        }
      }
    }

    /* Nx = 1; Ny = pan_n_from_file */
    nxout = 1;
    nyout = sample->backward_n_from_ray_file;
  } else {
    sample->backward_n_from_ray_file = 0;
    sample->backward_ray_from_file   = NULL;
  }

  /* if ipa3d calculation, copy sample grid from cloud grid */
  if (ipa3d) {
    nxout = nx_caoth;
    nyout = ny_caoth;

    dxout = dx_caoth / 1000.0;
    dyout = dy_caoth / 1000.0;

    if (!quiet)
      fprintf (stderr, " ... setting sample grid to %d x %d for ipa3d calculcation\n", nxout, nyout);
  }

  /* if no user-defined sample grid */
  if (nxout <= 0 || nyout <= 0 || dxout <= 0 || dyout <= 0) {
    if (!quiet)
      fprintf (stderr, " ... no user-defined sample grid ");

    if (nxout <= 0 || nyout <= 0) {
      nxout = 1;
      nyout = 1;
    }

    if (dx_caoth > 0 && dy_caoth > 0) {
      dxout = dx_caoth * (double) nx_caoth / (double) nxout / 1000.0;
      dyout = dy_caoth * (double) ny_caoth / (double) nyout / 1000.0;

      if (!quiet)
        fprintf (stderr, "-> copying domain size from cloud grid\n");
    } else {

      if (panorama || panorama_forward) {
        dxout = 1 / (double) nxout;
        dyout = 1 / (double) nyout;
        if (!quiet)
          fprintf (stderr, "-> in panorama mode, we take 1km by default\n");
      } else {
        nxout = 1;
        nyout = 1;

        if (spherical) {
          /*CE Numerical problems may occur if sample domain is too large in spherical geometry ... why??? */
          dxout = 100000;
          dyout = 100000;
        } else {
          /* use a domain size DBL_MAX / 10.0 - the factor 1 / 10 allows some photon movement */
          /* as the photon is allowed to leave the domain before it is forced back by the     */
          /* periodic boundary conditions and we don't want an overflow in this case.         */
          dxout = DBL_MAX / 10.0 / 1000.0; /* convert to km */
          dyout = DBL_MAX / 10.0 / 1000.0; /* convert to km */
        }

        if (!quiet) {
          fprintf (stderr, "-> defining an \"infinitely\" large domain\n");
          fprintf (stderr, " ...   of length %gm with one single sample pixel!\n", dxout);
        }
      }
    }
  }

  /* allocate memory for sample flag */
  sample->sample = calloc ((size_t) nlev, sizeof (int));

  sample->delX = dxout * 1000.0; /* convert to m */
  sample->delY = dyout * 1000.0; /* convert to m */

  sample->Nx = nxout;
  sample->Ny = nyout;

  sample->ixmin = ixmin;
  sample->ixmax = ixmax;
  sample->iymin = iymin;
  sample->iymax = iymax;

  /* This must be set to 1 else the mystic loop is not started for non-Lidar */
  sample->Nli = 1;

  if (!quiet)
    fprintf (stderr, " ... sample grid %d x %d pixels with area %gm x %gm\n", sample->Nx, sample->Ny, sample->delX, sample->delY);

  if (!quiet && spherical) {
    fprintf (stderr, " ... ATTENTION: A sample grid is defined for 1D spherical geometry. \n");
    fprintf (stderr, " ...      Be aware, that for 1D spherical geometry, absorbing boundary \n");
    fprintf (stderr, " ...      conditions are used and make sure that the domain is sufficiently large. \n");
  }

  /* process user altitudes */
  if (nzout > 0) {

    /* first check if zout is sorted in ascending order */
    for (lu = 1; lu < nzout; lu++)
      if (zout[lu] - zout[lu - 1] < 0) {
        fprintf (stderr, "Error, zout not sorted in ascending order\n");
        status--;
      }

    /* ??? special treatment of -999, surface ??? */

    if (!quiet)
      fprintf (stderr, " ... using %d altitudes from the uvspec input file:\n", nzout);

    for (lu = 0; lu < nzout; lu++) {
      if (zout[lu] >= 0) {
        for (kc = 0; kc < nlev; kc++) {
          if (float_equal (zd[nlev - 1 - kc], zout[lu])) {
            sample->sample[kc] = 1;
            break;
          }
        }
        if (kc == nlev) {
          fprintf (stderr, "\nFATAL error: altitude grid does not contain level %g\n", zout[lu]);
          fprintf (stderr, "which has been specified as output altitude\n");
          CHKERR (-1);
        }
      }
    }
  }

  if (!quiet) {
    fprintf (stderr, " ... sampling result at the following altitudes:\n");

    for (kc = 0; kc < nlev; kc++)
      if (sample->sample[kc])
        fprintf (stderr, " ...  %8.3f km\n", zd[nlev - 1 - kc]);

    fprintf (stderr, "\n");
  }

  sample->DoLE = DoLE;

  if (zout[0] < 0)
    sample->surface = 1;

  /* process radiance directions */
  if (numu > 0 && nphi > 0) { /* use data passed as function arguments */

    if (!quiet && !panorama && !panorama_forward)
      fprintf (stderr, " ... using %d radiance angles from the uvspec input file:\n", numu * nphi);

    /* allocate memory for sample radiances */
    Nd         = numu * nphi;
    sample->Nd = Nd;
    if (Nd > 0)
      sample->rad = calloc (Nd + sample->DoLE, sizeof (radang));

    for (iu = 0; iu < numu; iu++) {
      for (j = 0; j < nphi; j++) {

        if (!panorama && !panorama_forward) {
          /* copy uvspec radiance to sample array */
          sample->rad[id].theta = acosd (umu[iu]);
          sample->rad[id].phi   = phi[j];
          sample->rad[id].alpha = alpha;
        } else {
          sample->rad[id].theta = 0.0;
          sample->rad[id].phi   = 0.0;
          sample->rad[id].alpha = alpha;
        }

        if (!quiet && !panorama && !panorama_forward)
          fprintf (stderr, " ...  (%g, %g)\n", sample->rad[id].theta, sample->rad[id].phi);

        id++;
      }
    }
  }

  /****************************************************/
  /* Changed MYSTIC azimuth convention March 1, 2004! */
  /* From then on, the MYSTIC azimuth convention      */
  /* follows the DISORT convention, that is, as seen  */
  /* from the photon rather than from the detector.   */
  /* azimuth 0 = coming from the direction of the sun */
  /* (looking towards the sun) while azimuth 180 is   */
  /* going into the direction of the sun (or for the  */
  /* observer, sun in the back. If azigconv is set to  */
  /* MCAZIMUTH_OLD then the old convention is used.   */
  /*                                                  */
  /* Inside MYSTIC nothing is changed; azimuths       */
  /* still refer to the observing instrument.         */
  /*                                                  */
  /* To achieve this, we need to reverse the azimuths */
  /* at this location and store the original values   */
  /* for the output.                                  */
  /****************************************************/

  for (id = 0; id < sample->Nd; id++) {
    sample->rad[id].externalphi = sample->rad[id].phi;

    switch (aziconv) {
    case MCAZIMUTH_NEW:
      sample->rad[id].phi = 180.0 + sample->rad[id].phi;
      break;
    case MCAZIMUTH_OLD:
      sample->rad[id].phi = sample->rad[id].phi;
      break;
    default:
      fprintf (stderr, "Error, unknown azimuth convention %d\n", aziconv);
      CHKERR (-1);
    }

    /* shift to interval [0 ... 360] */
    while (sample->rad[id].phi < 0.0)
      sample->rad[id].phi += 360.0;

    while (sample->rad[id].phi > 360.0)
      sample->rad[id].phi -= 360.0;
  }

  sample->delta_scaling = delta_scaling_start;

  if (stheta != NULL)
    free (stheta);
  if (sphi != NULL)
    free (sphi);
  if (salpha != NULL)
    free (salpha);

  if (sample->Nd > 0) {
    sample->escape = escape;

    if (!quiet && !panorama && !panorama_forward)
      fprintf (stderr, " ... calculating radiances for %d directions\n", sample->Nd);
  }

  sample->ncirc = ncirc;
  sample->dcirc = dcirc;

#if HAVE_LIDAR
  if (jacobian && locest) {
    sample->LLE_jacobian = jacobian;
    if (jacobian_std)
      sample->LLE_jacobian_std = jacobian_std;
    sample->jac_eps = 1e-1;
  }
#endif

  if ((jacobian || jacobian3D) && escape) {
    if (backward) {
      sample->jacobian   = jacobian;
      sample->jacobian3D = jacobian3D;
      if (jacobian_std)
        sample->jacobian_std = jacobian_std;
      sample->jac_eps = 1e-3;
    } else {
      CHKERROUT (-1, "molecular jacobian can only be calculated for backward montecarlo");
    }
  }

  if (calcstd) {
    sample->std       = 1;
    sample->std_value = std_value;
  } else {
    sample->std       = 0;
    sample->std_value = 0;
  }

  /* radial pathlength sampling */
  sample->Nr = radpat_Nr;
  sample->Nt = radpat_Nt;
  sample->dt = radpat_dt;

  /* fit circle into domain */
  xr = 0.5 * (double) sample->Nx * sample->delX;
  yr = 0.5 * (double) sample->Ny * sample->delY;

  sample->dr = (xr < yr ? xr : yr) / (double) sample->Nr;

  if (!quiet && sample->Nr > 0)
    fprintf (stderr, " ... radial pathlength sampling with intervals dr = %.6e and dt = %.6e\n", sample->dr, sample->dt);

  free (z_hist);

  /* number of levels to be sampled */
  for (kc = 0; kc < nlev; kc++)
    sample->nzout += sample->sample[kc];

  sample->backward_writeallpixels = backward_writeallpixels;
  sample->backward_writeback      = backward_writeback;

  if (backward)
    sample->backward_Nz = sample->nzout - 1;
  else
    sample->backward_Nz = 0;

  /* backward sampling area - may not be needed but is initialized anyway */
  sample->backward_islower = backward_islower;
  sample->backward_jslower = backward_jslower;
  sample->backward_isupper = backward_isupper;
  sample->backward_jsupper = backward_jsupper;
  sample->backward_isstep  = backward_isstep;
  sample->backward_jsstep  = backward_jsstep;

  if (sample->backward_islower == -999 && sample->backward_jslower == -999 && sample->backward_isupper == -999 &&
      sample->backward_jsupper == -999) {
    /* special case: no user-defined backward pixels - calculating all sample pixels */
  } else {
    if (sample->backward_islower < 0)
      sample->backward_islower = 0;
    if (sample->backward_jslower < 0)
      sample->backward_jslower = 0;
    if (sample->backward_isupper < 0)
      sample->backward_isupper = 0;
    if (sample->backward_jsupper < 0)
      sample->backward_jsupper = 0;

    if (sample->backward_islower >= sample->Nx)
      sample->backward_islower = sample->Nx - 1;
    if (sample->backward_jslower >= sample->Ny)
      sample->backward_jslower = sample->Ny - 1;
    if (sample->backward_isupper >= sample->Nx)
      sample->backward_isupper = sample->Nx - 1;
    if (sample->backward_jsupper >= sample->Ny)
      sample->backward_jsupper = sample->Ny - 1;
  }

  sample->reference_to_NN = 0;

  /* calculate zenith and azimuth angles */
  if (backward) {
    if (Nd > 0) {
      sample->backward        = MCBACKWARD_RADIANCE;
      sample->reference_to_NN = reference_to_NN;
    } else
      sample->backward = output; /* irradiance or actinic flux - only one at*/
  }

  if (sample->backward && !sample->escape) {
    sample->escape = 1;

    if (!quiet)
      fprintf (stderr, " ... switching on local estimate (default for backward calculations)\n");
  }

  if (sample->backward && (sample->Nr > 0 && sample->dt > 0)) {
    CHKERROUT (-1, "radial/pathlength and backward cannot be combined");
  }

  /* boundary conditions */
  sample->bcond = bcond;

#ifdef CLDPRP
  /* activate sampling of cloud properties (reff,tau,lwc gradients) by photons */
  if (sample_cldprp)
    sample->cldprp = 1;
#endif

  sample->coherent_backscatter       = coherent_backscatter;
  sample->coherent_backscatter_lmode = coherent_backscatter_lmode;
  sample->ris_factor                 = ris_factor;
  sample->ris_optical_depth          = ris_optical_depth;

  /* spherical geometry - always use backwards and absorbing boundary conditions! */
  sample->spherical = spherical;
  if (spherical) {
    if (Nd > 0)
      sample->backward = MCBACKWARD_RADIANCE;
    else
      sample->backward = output;

    sample->escape                  = 1;
    sample->backward_writeallpixels = 0;

    if (!panorama) {

      if (!quiet) {
        if (sample->backward_islower != 0 || sample->backward_jslower != 0 || sample->backward_isupper != 0 ||
            sample->backward_jsupper != 0) {
          fprintf (stderr, " ... backward sampling area ignored in spherical geometry!\n");
          fprintf (stderr, " ... (photons are always started at the center of the domain)\n");
        }
      }
      sample->backward_islower = 0;
      sample->backward_jslower = 0;
      sample->backward_isupper = 0;
      sample->backward_jsupper = 0;
      sample->backward_isstep  = 1;
      sample->backward_jsstep  = 1;
    }

    if (sample->bcond != MCBCOND_ABSORB && !quiet)
      fprintf (stderr, " ... using absorbing boundary conditions for spherical\n");

    sample->bcond = MCBCOND_ABSORB;
  }
  sample->refraction = refraction;

  /* High spectral resolution calculation using importance sampling method */
  sample->spectral_is = spectral_is;

  /* Box airmass factors */
  sample->boxairmass = boxairmass;

  /* Box airmass factors in 3D */
  sample->boxairmass3D = boxairmass3D;

  /* boxairmass factors can only be calculated in backward mode */
  if (!spherical && (boxairmass || boxairmass3D)) {
    if (Nd > 0)
      sample->backward = MCBACKWARD_RADIANCE;
    else
      sample->backward = output;

    sample->escape                  = 1;
    sample->backward_writeallpixels = 0;
  }

  /* spherical geometry in 3D! */
  sample->spherical3D = spherical3D;

  sample->spherical3D_scene         = spherical3D_scene;
  sample->spherical3D_scene_lon_min = spherical3D_scene_lon_min;
  sample->spherical3D_scene_lon_max = spherical3D_scene_lon_max;
  sample->spherical3D_scene_lat_min = spherical3D_scene_lat_min;
  sample->spherical3D_scene_lat_max = spherical3D_scene_lat_max;

  /* define if using tipa-dir */
  if (tipa == TIPA_DIR3D)
    sample->tipadir = tipa;
  else
    sample->tipadir = 0;

  sample->polarisation       = polarisation;
  sample->polarisation_state = polarisation_state;
  sample->nstokes            = nstokes;
  if (polarisation)
    sample->nphamat = NPHAMAT;
  else
    sample->nphamat = 1;

  /* if non-zero, photon is stopped as soon as it reaches maxscatters scatterings */
  sample->maxscatters = maxscatters;
  sample->minscatters = minscatters;

  if (blitz_position != NULL)
    sample->blitz_position = blitz_position;

  /* Switch to use get_phase_phase_matrix_total in scattering in mystic.c  */
  /* instead of get_phase_matrix_caoth. Using the latter can be dangerous! */
  /* Always use get_phase_matrix_total when using lidar or radar!          */
  sample->use_phase_matrix_total = 1;

#if HAVE_VROOM
  /* cases where vroom should be off */
  if (sample->delta_scaling != -1) {
    if (vroom == 1) {
      CHKERROUT (-1, "You are trying to use vroom AND delta-scaling at the same time, this does not make sense");
    }
    vroom = 0;
  }

  if (!sample->escape)
    vroom = 0;

  if (vroom == 1 && sample->Nd > 1) {
    CHKERROUT (-1, "You can not use mc_vroom when you are sampling more than one direction");
  }

  sample->vroom = vroom;

  //  status = set_vroom_settings (vroom, sample, quiet);
  //
  //  if (status!=0) {
  //    fprintf (stderr, "Error %d returned from set_vroom_settings\n", status);
  //    return status;
  //  }

#endif

#if HAVE_LIDAR
  /*two experimental switches added by Christian Pause*/
  /*planar switch:*/
  /*0 = standard lidar/radar*/
  /*1 = transmitter uses planar waves*/
  /*2 = detector uses planar waves*/
  /*3 = both transmitter and detector planar waves*/
  if (coherent_backscatter_lmode == 2)
    sample->planar_switch = 3;
  else
    sample->planar_switch = 0;
  /*switch off CB SV when using CB, resulting SV is only average of both paths*/
  sample->nocb_switch = 0;
  if (coherent_backscatter_off == 1)
    sample->nocb_switch = 1;

  sample->LLE_eps_ddis_upf = 0.0;
  sample->LLE_eps_ddis_uda = 0.0;

  /* Read laser and detector data for local estimator */
  if (locest) {

    status = set_lidar_settings (locest, sample);
    if (status != 0) {
      fprintf (stderr, "Error %d returned by set_lidar_settings \n", status);
      return status;
    }

    /* check if lidar file exists and can be opened for reading */
    ftmp = fopen (lidarfilename, "r");
    if (ftmp == NULL) {
      CHKERROUT (-1, "Either define the lidar in the input file or provide lidar.dat");
    } else {
      fclose (ftmp);

      /* read lidar.dat */
      status = read_lidar (lidarfilename, sample, quiet);

      if (status != 0) {
        fprintf (stderr, "Error %d reading %s\n", status, lidarfilename);
        return status;
      }
    }

    /* define some variables */
    if (sample->LLE_Nt > 0)
      sample->maxpathlength = sample->lidar_t[sample->LLE_Nt - 1];

    for (ili = 0; ili < sample->Nli; ili++) {

      sample->laser[ili].cosalpha  = cos (sample->laser[ili].alpha);
      sample->laser[ili].omega     = 2.0 * M_PI * (1.0 - sample->laser[ili].cosalpha);
      sample->laser[ili].dir.dx[0] = 0.;
      sample->laser[ili].dir.dx[1] = 0.;
      sample->laser[ili].dir.dx[2] = -1.;
      sample->laser[ili].area      = M_PI * sample->laser[ili].radius * sample->laser[ili].radius;

      /* -90.: conformity with phi0 definition */
      cosd_temp = cosd (sample->laser[ili].theta);

      /* XXX quick fix for numerical error */
      if (sample->LidarLocEst == MCRADAR) {
        if (cosd_temp >= -1.e-10 && cosd_temp <= 1.e-10)
          cosd_temp = 0.;
      }

      new_direction (cosd_temp, sample->laser[ili].phi - 90., &(sample->laser[ili].dir), 0.);

      sample->lidar[ili].cosalpha = calloc (sample->LLE_No, sizeof (double));
      sample->lidar[ili].omega    = calloc (sample->LLE_No, sizeof (double));
      sample->lidar[ili].azphi    = calloc (sample->LLE_Na + 1, sizeof (double));

      if (sample->LLE_islinear == 2) {
        /* cartesian geometry */
        sample->lidar[ili].cosalpha[0] = cos (sample->lidar[ili].alpha);
        sample->LLE_afac               = 2.0 * sin (sample->lidar[ili].alpha) / (double) sample->LLE_No;
        for (io = 0; io < sample->LLE_No + 1; io++) {
          /* this is the geometry */
          sample->lidar[ili].azphi[io] = -sin (sample->lidar[ili].alpha) + (double) io * sample->LLE_afac;
          sample->lidar[ili].omega[io] = sample->lidar[ili].azphi[io];
        }
      } else {
        alphaact = sample->lidar[ili].alpha;
        for (io = 0; io < sample->LLE_No; io++) {
          sample->lidar[ili].cosalpha[io] = cos (alphaact);
          sample->lidar[ili].omega[io]    = 2.0 * M_PI * (1.0 - sample->lidar[ili].cosalpha[io]);
          if (sample->LLE_islinear)
            alphaact += sample->LLE_afac * sample->lidar[ili].alpha;
          else
            alphaact *= sample->LLE_afac;
        }
        for (ia = 0; ia < sample->LLE_Na + 1; ia++)
          sample->lidar[ili].azphi[ia] = 360. * ((double) ia) / ((double) sample->LLE_Na);
      }

      sample->lidar[ili].dir.dx[0] = 0.;
      sample->lidar[ili].dir.dx[1] = 0.;
      sample->lidar[ili].dir.dx[2] = -1.;
      sample->lidar[ili].area      = M_PI * sample->lidar[ili].radius * sample->lidar[ili].radius;

      /* -90.: conformity with phi0 definition */
      new_direction (cosd (sample->lidar[ili].theta), sample->lidar[ili].phi - 90., &(sample->lidar[ili].dir), 0.);

      /* perp1 always looking downwards, in line with dir.dx when seen from above */
      /* if dir.dx is vertical, perp1= (-1,0,0) */
      sample->lidar[ili].dx_perp1[2] = -sqrt (sample->lidar[ili].dir.dx[0] * sample->lidar[ili].dir.dx[0] +
                                              sample->lidar[ili].dir.dx[1] * sample->lidar[ili].dir.dx[1]);
      if (sample->lidar[ili].dx_perp1[2] == 0.0) {
        sample->lidar[ili].dx_perp1[0] = -1.0;
        sample->lidar[ili].dx_perp1[1] = 0.0;
      } else {
        scr                            = -sample->lidar[ili].dir.dx[2] / sample->lidar[ili].dx_perp1[2];
        sample->lidar[ili].dx_perp1[0] = sample->lidar[ili].dir.dx[0] * scr;
        sample->lidar[ili].dx_perp1[1] = sample->lidar[ili].dir.dx[1] * scr;
      }

      v_cross_product (sample->lidar[ili].dir.dx, sample->lidar[ili].dx_perp1, sample->lidar[ili].dx_perp2);

      sample->lidar[ili].s_det = sin (sample->lidar[ili].alpha);
      sample->lidar[ili].c_det = cos (sample->lidar[ili].alpha);
      sample->lidar[ili].t_det = tan (sample->lidar[ili].alpha);
      sample->lidar[ili].z_det = -sample->lidar[ili].radius / sample->lidar[ili].t_det;
      for (i = 0; i < 3; i++)
        sample->lidar[ili].x_cc[i] = sample->lidar[ili].x[i] + sample->lidar[ili].dir.dx[i] * sample->lidar[ili].z_det;

      if (!quiet) {
        fprintf (stderr,
                 "laser %d %e %e %e\n",
                 ili,
                 sample->laser[ili].dir.dx[0],
                 sample->laser[ili].dir.dx[1],
                 sample->laser[ili].dir.dx[2]);
        fprintf (stderr,
                 "lidar %d %e %e %e\n",
                 ili,
                 sample->lidar[ili].dir.dx[0],
                 sample->lidar[ili].dir.dx[1],
                 sample->lidar[ili].dir.dx[2]);
      }
    }
  } else
    /* This must be set to 1 else the mystic loop is not started for non-Lidar */
    sample->Nli = 1;
#else
  if (locest) {
    CHKERROUT (-1, "This version of MYSTIC does not contain a lidar/radar simulator");
  }
#endif

  return 0;
}
#endif

/***********************************************************************************/
/* Function:                                                              @62_30i@ */
/* Description:                                                                    */
/*  Setup a 3D atmosphere structure, to be used as input to MYSTIC.                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int setup_caoth3D (input_struct input, output_struct* output) {
  int  i = 0, isp = 0, kc = 0, status = 0, hori_agree = 0;
  char caoth3d_filename[FILENAME_MAX];

  int    Nx = 0, Ny = 0;
  double delX = 0.0, delY = 0.0;

  /***********************/
  /* Read 3D caoth files */
  /***********************/

  output->caoth3d = calloc ((size_t) input.n_caoth, sizeof (caoth3d_out_struct));
  if (output->caoth3d == NULL)
    return mem_err_out ("output->caoth3d", ERROR_POSITION);

  for (isp = 0; isp < input.n_caoth; isp++) {
    /* reset counters */
    output->caoth3d[isp].nthreed = 0;

    /* copy no_scattering and no_absorption information */
    output->caoth3d[isp].no_scattering = input.caoth[isp].no_scattering;
    output->caoth3d[isp].no_absorption = input.caoth[isp].no_absorption;

    /* copy oriented particle information */
    output->caoth3d[isp].oriented_crystal_fraction      = input.caoth[isp].oriented_crystal_fraction;
    output->caoth3d[isp].oriented_crystal_angdist_width = input.caoth[isp].oriented_crystal_angdist_width;
    if (input.caoth[isp].oriented_crystal_fraction > 0 && !input.quiet)
      fprintf (stderr,
               "\n!!!! Warning !!!! \nThere can be spikes in the results due to the very narrow phase function for scattering on "
               "oriented ice crystals,\nespecially if the oriented_crystal_fraction is too large or the "
               "oriented_crystal_angdist_width is too small.\n\n ");

    /* copy cldproperties/nonHGcld information to 3D structure */
    /* even if no 3D caoth but only a 1D caoth is defined      */
    output->caoth3d[isp].name = calloc (strlen (input.caoth[isp].name) + 1, sizeof (char));
    strcpy (output->caoth3d[isp].name, input.caoth[isp].name);
    output->caoth3d[isp].fullname = calloc (strlen (input.caoth[isp].fullname) + 1, sizeof (char));
    strcpy (output->caoth3d[isp].fullname, input.caoth[isp].fullname);
    output->caoth3d[isp].nonHG         = output->caoth[isp].nonHG;
    output->caoth3d[isp].cldproperties = output->caoth[isp].cldproperties;

    if (input.caoth[isp].source == CAOTH_FROM_3D) {
      strcpy (caoth3d_filename, input.caoth[isp].filename);
      if (!input.quiet)
        fprintf (stderr, " ... reading 3D data from %s for %s\n", caoth3d_filename, input.caoth[isp].fullname);

      /* read 3D caoth data */
      status = read_caoth3D_file (caoth3d_filename,
                                  input,
                                  input.caoth[isp],
                                  output,
                                  &(output->caoth[isp]),
                                  &(output->caoth3d[isp]),
                                  &(output->atm.nmom));
      if (status != 0) {
        fprintf (stderr,
                 "Error %d reading 3D %s file %s in (line %d, function '%s', file '%s')\n",
                 status,
                 input.caoth[isp].fullname,
                 caoth3d_filename,
                 ERROR_POSITION);
        return status;
      }

#if HAVE_TIPA
      /* ulrike: implementation of tica, where the "tilted caoth matrix"
	 is used as input for the calculation of the direct AND the diffuse radiation
	 This is possible with 1D-solvers, such as disort2 and rodents,
	 where tipa dirdiff must be specified in the input-file
	 as well as with mystic where mc_tipa dirdiff must be specified */

      if (input.tipa == TIPA_DIRDIFF || input.rte.mc.tipa == TIPA_DIRDIFF)
        status = tipa_dirdiff (&(output->caoth3d[isp]),
                               &(output->caoth[isp].tipa),
                               input.tipa,
                               input.atm.sza,
                               input.atm.phi0,
                               output->atm.zout_sur,
                               input.rte.mc.tipa);
      if (status)
        return fct_err_out (status, "tipa_dirdiff", ERROR_POSITION);
#endif

      /* add 3D caoth levels to atmospheric z-grid */

      /* sort by descending altitude */
      {
        float zd_temp[output->caoth3d[isp].nlyr + 1];
        for (kc = 0; kc <= output->caoth3d[isp].nlyr; kc++)
          zd_temp[kc] = output->caoth3d[isp].zd[output->caoth3d[isp].nlyr - kc];

        set_common_z_grid (zd_temp, output->caoth3d[isp].nlyr + 1, &output->atm.zd_common, &output->atm.nlev_common);
      }

      if (!input.quiet)
        fprintf (stderr,
                 " ... adding %d 3D %s levels to the common grid\n",
                 output->caoth3d[isp].nlyr + 1,
                 input.caoth[isp].fullname);
    } else { /* no 3D caoth data */
      output->caoth3d[isp].nthreed = 0;
      output->caoth3d[isp].nlyr    = 0;
    }

    if (output->caoth[isp].n_raytracing_prop > 0) {
      output->caoth3d[isp].raytracing_prop   = calloc (output->caoth[isp].n_raytracing_prop, sizeof (crystal_prop_struct));
      output->caoth3d[isp].n_raytracing_prop = output->caoth[isp].n_raytracing_prop;

      for (i = 0; i < output->caoth3d[isp].n_raytracing_prop; i++) {
        output->caoth3d[isp].raytracing_prop[i].name =
          calloc (strlen (output->caoth[isp].raytracing_prop[i].name) + 1, sizeof (char));
        strcpy (output->caoth3d[isp].raytracing_prop[i].name, output->caoth[isp].raytracing_prop[i].name);

        output->caoth3d[isp].raytracing_prop[i].fraction          = output->caoth[isp].raytracing_prop[i].fraction;
        output->caoth3d[isp].raytracing_prop[i].oriented_fraction = output->caoth[isp].raytracing_prop[i].oriented_fraction;
        output->caoth3d[isp].raytracing_prop[i].angdist_width     = output->caoth[isp].raytracing_prop[i].angdist_width;
        output->caoth3d[isp].raytracing_prop[i].orientation_dof   = output->caoth[isp].raytracing_prop[i].orientation_dof;
      }
    }

    /* we require that caoth are defined */
    /* on identical horizontal grids       */
    if (output->caoth3d[isp].nthreed > 0) {
      delX = output->caoth3d[isp].delX;
      delY = output->caoth3d[isp].delY;
      Nx   = output->caoth3d[isp].Nx;
      Ny   = output->caoth3d[isp].Ny;
    }
  } /* end loop isp */

  for (isp = 0; isp < input.n_caoth; isp++) {
    if (output->caoth3d[isp].nthreed > 0) {
      if (output->caoth3d[isp].delX != delX || output->caoth3d[isp].delY != delY || output->caoth3d[isp].Nx != Nx ||
          output->caoth3d[isp].Ny != Ny) {
        CHKERROUT (-1, "3D input must be defined on identical horizontal grids");
      } else
        hori_agree = 1;
    } else {
      output->caoth3d[isp].delX = delX;
      output->caoth3d[isp].delY = delY;
      output->caoth3d[isp].Nx   = Nx;
      output->caoth3d[isp].Ny   = Ny;
    }
  }

  /* ulrike */
  //  if (input.rte.solver != SOLVER_MONTECARLO && hori_agree) { /* hori_agree means that threed exists!*/
  // BM, 15.7.2020: Shouldn't that be rather ipa3d than solver != MONTECARLO */
  if (input.ipa3d && hori_agree) { /* hori_agree means that threed exists!*/
    output->niipa = Nx;
    output->njipa = Ny;
  }

  // initialize output->islower etc and sample grid size - it is set by mystic() anyway
  // but we require it for non-mystic applications like dynamic_tenstream()

  output->islower = 0;
  output->jslower = 0;
  output->isupper = Nx - 1;
  output->jsupper = Ny - 1;
  output->isstep  = 1;
  output->jsstep  = 1;

  output->mc.sample.Nx = Nx;
  output->mc.sample.Ny = Ny;

  output->atm.Nxcld = Nx;
  output->atm.Nycld = Ny;

  if (!input.quiet && hori_agree)
    fprintf (stderr, " ... horizontal grids of 3D fields agree\n");

  return 0;
}

static int set_sensor_position (const input_struct input,
                                const int          sample_Nx,
                                const int          sample_Ny,
                                const double       sample_delX,
                                const double       sample_delY,
                                const double       domain_height, /* in [m] */
                                int*               sensorposition,
                                double             sensorposition_x[3],
                                int*               sensordirection,
                                double             sensordirection_dx[3],
                                int*               sample_panorama,
                                int*               sample_panorama_forward,
                                int*               sample_pan_no_pixel,
                                int*               sample_pan_quicklook,
                                int*               sample_pan_distr_photons_over_pixel,
                                int*               sample_pan_with_direct_rad,
                                int*               sample_pan_weight_with_cos,
                                int*               sample_pan_circumsolar_var_red,
                                int*               sample_pan_transmittance,
                                int*               sample_pan_alignment,
                                double*            sample_pan_alpha,
                                double*            sample_pan_phi_min,
                                double*            sample_pan_phi_max,
                                double*            sample_pan_theta_min,
                                double*            sample_pan_theta_max,
                                double*            sample_pan_umu_min,
                                double*            sample_pan_umu_max,
                                double*            sample_pan_phi_from_file,
                                double*            sample_pan_theta_from_file,
                                int                sample_pan_n_from_file,
                                double*            sample_align_theta,
                                double*            sample_align_phi,
                                double*            sample_pan_roll,
                                double*            sample_pan_pitch,
                                double*            sample_pan_yaw,
                                double*            sample_add_phi_when_aligning) {

  /* user-defined sensor position */
  switch (input.rte.mc.sensorposition) {
  case MC_SENSORPOSITION_NONE:
    break;
  case MC_SENSORPOSITION_CARTESIAN:
    *sensorposition     = 1;
    sensorposition_x[0] = input.rte.mc.sensorposition_x;
    sensorposition_x[1] = input.rte.mc.sensorposition_y;
    sensorposition_x[2] = input.rte.mc.sensorposition_z;

    if (sensorposition_x[2] > domain_height) {
      fprintf (stderr,
               "Cannot have a camera position (%f %f %f km) above the domain top (%f km)!\n"
               "Increase your atmosphere file to contain the elevation of the camera position\n",
               1e-3 * sensorposition_x[0],
               1e-3 * sensorposition_x[1],
               1e-3 * sensorposition_x[2],
               1e-3 * domain_height);
      CHKERROUT (-1, "Camera cannot be above the domain top");
    }

    /* check if sensorposition is within domain; if not, shift periodically */
    int shifted = 0;
    while (sensorposition_x[0] < 0) {
      sensorposition_x[0] += sample_Nx * sample_delX;
      shifted++;
    }

    while (sensorposition_x[0] >= sample_Nx * sample_delX) {
      sensorposition_x[0] -= sample_Nx * sample_delX;
      shifted++;
    }

    while (sensorposition_x[1] < 0) {
      sensorposition_x[1] += sample_Ny * sample_delY;
      shifted++;
    }

    while (sensorposition_x[1] >= sample_Ny * sample_delY) {
      sensorposition_x[1] -= sample_Ny * sample_delY;
      shifted++;
    }

    if (!input.quiet && shifted)
      fprintf (stderr, " ... shifting sensorposition perioducally to %f %f\n", sensorposition_x[0], sensorposition_x[1]);

    break;
  case MC_SENSORPOSITION_SPHERICAL:
    *sensorposition     = 1;
    const double temp   = 1e3 * input.r_earth + input.rte.mc.sensorposition_alt;
    sensorposition_x[0] = temp * cosd (input.rte.mc.sensorposition_lat) * cosd (input.rte.mc.sensorposition_lon);
    sensorposition_x[1] = temp * cosd (input.rte.mc.sensorposition_lat) * sind (input.rte.mc.sensorposition_lon);
    sensorposition_x[2] = temp * sind (input.rte.mc.sensorposition_lat);
    break;
  default:
    fprintf (stderr, "Error, unknown sensorposition %d\n", input.rte.mc.sensorposition);
    CHKERR (-1);
  }

  if (!input.quiet && sensorposition) {
    fprintf (stderr, " ... user-defined sensor position\n");
    fprintf (stderr, "     (%f, %f, %f)\n", sensorposition_x[0], sensorposition_x[1], sensorposition_x[2]);
  }

  /* panorama option */

  switch (input.rte.mc.panorama) {
  case PAN_MODE_NONE:
    break;
  case PAN_MODE_CAMERA:
    *sample_pan_phi_min = input.rte.mc.pan_angles[0];
    *sample_pan_phi_max = input.rte.mc.pan_angles[1];

    if (*sample_pan_phi_min > *sample_pan_phi_max) {
      fprintf (stderr,
               "Error! panorama phi_min is larger than phi_max! phi_min %e phi_max %e\n",
               *sample_pan_phi_min,
               *sample_pan_phi_max);
      CHKERR (-1);
    }

    *sample_pan_theta_min = input.rte.mc.pan_angles[2];
    *sample_pan_theta_max = input.rte.mc.pan_angles[3];

    if (*sample_pan_theta_min > *sample_pan_theta_max) {
      fprintf (stderr,
               "Error! panorama theta_min is larger than theta_max! theta_min %e theta_max %e\n",
               *sample_pan_theta_min,
               *sample_pan_theta_max);
      CHKERR (-1);
    }

    if (*sample_pan_theta_min < 0.0) {
      fprintf (stderr, "Error! panorama theta_min is smaller than 0.0! theta_min %e\n", *sample_pan_theta_min);
      CHKERR (-1);
    }

    if (*sample_pan_theta_max > 180.0) {
      fprintf (stderr, "Error! panorama theta_max is larger than 180.0! theta_max %e\n", *sample_pan_theta_max);
      CHKERR (-1);
    }

    *sample_pan_umu_min = cosd (*sample_pan_theta_min);
    *sample_pan_umu_max = cosd (*sample_pan_theta_max);

    /* new definition of panorama azimuth */
    if (!(*sample_pan_umu_min < -2.0)) {
      switch (input.rte.mc.azimuth) {
      case MCAZIMUTH_NEW:
        *sample_pan_phi_min += 180.0;
        *sample_pan_phi_max += 180.0;
        break;
      case MCAZIMUTH_OLD:
        break;
      default:
        fprintf (stderr, "Error, unknown azimuth convention %d\n", input.rte.mc.azimuth);
        CHKERR (-1);
      }
      /* make sure that phi_min is between 0 and 360 degrees;
	 then phi_max is between 0 and 720 degrees */
      while (*sample_pan_phi_min < 0.0) {
        *sample_pan_phi_min += 360.0;
        *sample_pan_phi_max += 360.0;
      }
      while (*sample_pan_phi_min >= 360.0) {
        *sample_pan_phi_min -= 360.0;
        *sample_pan_phi_max -= 360.0;
      }
    }

    /* from here on it's identical for PAN_MODE_FILE */
  case PAN_MODE_FILE:
    *sample_panorama = 1;

    switch (input.rte.mc.pan_alignment) {
    case (PAN_ALIGNMENT_SUN):
      *sample_pan_alignment = MCPAN_ALIGNMENT_SUN;
      *sample_align_theta   = 180.0 - input.atm.sza;
      *sample_align_phi     = 180.0 + input.atm.phi0;
      break;
    case (PAN_ALIGNMENT_MU):
      *sample_pan_alignment = MCPAN_ALIGNMENT_MU;
      *sample_align_theta   = acosd (input.rte.umu[0]);
      *sample_align_phi     = input.rte.phi[0];
      break;
    case (PAN_ALIGNMENT_AIRCRAFT):
      *sample_pan_alignment = MCPAN_ALIGNMENT_AIRCRAFT;
      // pitch, roll, and yaw implemented by Lea Volkmer
      // default: flying towards North, looking straight down
      // align coordinate system such that z-axis points horizontally towards the North
      // and such that default (roll = pitch = yaw = 0) is flying towards the North
      *sample_align_theta = 90.0;
      *sample_align_phi   = 180.;

      // write roll, pitch and yaw angles to sample_struct
      *sample_pan_roll  = input.rte.mc.pan_roll;
      *sample_pan_pitch = input.rte.mc.pan_pitch;
      *sample_pan_yaw   = input.rte.mc.pan_yaw;

      break;

    case (PAN_ALIGNMENT_ZENITH):
      if (!input.rte.mc.spherical[DIM_3D]) {
        /* ZENITH means NONE in absence of spherical3D */
        *sample_pan_alignment = MCPAN_ALIGNMENT_NONE;
        break;
      }
      *sample_pan_alignment = MCPAN_ALIGNMENT_MU;

      *sample_align_theta =
        180.0 -
        atand (sqrt (sensorposition_x[0] * sensorposition_x[0] + sensorposition_x[1] * sensorposition_x[1]) / sensorposition_x[2]);
      *sample_align_phi = 90.0 - atand (sensorposition_x[1] / sensorposition_x[0]);

      if (sensorposition_x[0] < 0.0)
        *sample_align_phi += 180;

      /* or align_theta < 0 */
      if (sensorposition_x[2] < 0.0)
        *sample_align_theta = 180.0 + *sample_align_theta;

      break;
    case (PAN_ALIGNMENT_NONE):
      *sample_pan_alignment = MCPAN_ALIGNMENT_NONE;
      break;
    default:
      fprintf (stderr, "Error, unknown panorama alignment %d\n", input.rte.mc.pan_alignment);
      CHKERR (-1);
    }

    break;

  case PAN_MODE_BACKWARD_FILE:
    *sample_panorama = 1;

    break;

  case PAN_MODE_SATELLITE:
    *sample_panorama = 1;

    *sample_pan_phi_min = 180. + input.rte.mc.pan_angles[1];
    *sample_pan_phi_max = 180. + input.rte.mc.pan_angles[0];

    *sample_pan_theta_min = 90.0 - input.rte.mc.pan_angles[2];
    *sample_pan_theta_max = 90.0 - input.rte.mc.pan_angles[3];

    *sample_pan_umu_min = cosd (*sample_pan_theta_min);
    *sample_pan_umu_max = cosd (*sample_pan_theta_max);

    /* satellite automatically means pan_alignment mu */
    *sample_pan_alignment = MCPAN_ALIGNMENT_MU;

    if (sensorposition_x[0] == 0.0 && sensorposition_x[1] == 0.0) {
      /* special case satellite above pole */
      *sample_align_theta = 90.0;
      if (sensorposition_x[2] < 0.0)
        *sample_align_phi = 270.0;
      else
        *sample_align_phi = 90.0;
    } else {
      *sample_align_theta =
        atand (sensorposition_x[2] / sqrt (sensorposition_x[0] * sensorposition_x[0] + sensorposition_x[1] * sensorposition_x[1]));
      if (sensorposition_x[0] != 0) {
        *sample_align_phi = 90.0 - atand (sensorposition_x[1] / sensorposition_x[0]);
        if (sensorposition_x[0] < 0.0)
          *sample_align_phi += 180.0;
      } else if (sensorposition_x[1] > 0)
        *sample_align_phi = 0.0;
      else
        *sample_align_phi = 180.0;
    }

    /* in case satellite south of equator, azimuth must be turned by 180 degrees */
    if (sensorposition_x[2] < 0.0) {
      *sample_pan_phi_min += 180.0;
      *sample_pan_phi_max += 180.0;
    }
    if (sensorposition_x[2] == 0.0)
      *sample_add_phi_when_aligning = 180.0 + *sample_align_phi;

    /* make sure that phi_min is between 0 and 360 degrees;
       then phi_max is between 0 and 720 degrees */
    while (*sample_pan_phi_min < 0.0) {
      *sample_pan_phi_min += 360.0;
      *sample_pan_phi_max += 360.0;
    }
    while (*sample_pan_phi_min >= 360.0) {
      *sample_pan_phi_min -= 360.0;
      *sample_pan_phi_max -= 360.0;
    }

    if (*sample_pan_phi_max - *sample_pan_phi_min > 360.0) {
      fprintf (stderr,
               "Error! satellite has more than 360 degrees azimuthal viewing angle! lon_max-lon_min = %e - %e > 360\n",
               *sample_pan_phi_max,
               *sample_pan_phi_min);
      CHKERR (-1);
    }

    if (*sample_pan_theta_max < 0.0 || *sample_pan_theta_min > 180.0) {
      fprintf (stderr,
               "Error! satellite theta angles must be within -90 and +90 degrees! lat_max, lat_min = %e, %e\n",
               90.0 - *sample_pan_theta_min,
               90.0 - *sample_pan_theta_max);
      CHKERR (-1);
    }

    break;
  default:
    fprintf (stderr, "Error, unknown panorama mode %d\n", input.rte.mc.panorama);
    CHKERR (-1);
  }

  /* strict forward tracing - count photons into angular bins covering full 4pi */
  if (input.rte.mc.panorama_forward) {
    *sample_panorama_forward = 1;

    *sample_pan_phi_min   = input.rte.mc.pan_angles[0];
    *sample_pan_phi_max   = input.rte.mc.pan_angles[1];
    *sample_pan_theta_min = input.rte.mc.pan_angles[2];
    *sample_pan_theta_max = input.rte.mc.pan_angles[3];

    *sample_pan_umu_min = cosd (*sample_pan_theta_min);
    *sample_pan_umu_max = cosd (*sample_pan_theta_max);

    /* new definition of panorama azimuth */
    switch (input.rte.mc.azimuth) {
    case MCAZIMUTH_NEW:
      *sample_pan_phi_min += 180.0;
      *sample_pan_phi_max += 180.0;
      break;
    case MCAZIMUTH_OLD:
      break;
    default:
      fprintf (stderr, "Error, unknown azimuth convention %d\n", input.rte.mc.azimuth);
      CHKERR (-1);
    }
  }

  *sample_pan_no_pixel                 = input.rte.mc.pan[PAN_FLAG_NO_PIXEL];
  *sample_pan_quicklook                = input.rte.mc.pan[PAN_FLAG_QUICKLOOK];
  *sample_pan_distr_photons_over_pixel = input.rte.mc.pan[PAN_FLAG_DISTR_PHOTONS_OVER_PIXEL];
  *sample_pan_with_direct_rad          = input.rte.mc.pan[PAN_FLAG_WITH_DIRECT_RAD];
  *sample_pan_weight_with_cos          = input.rte.mc.pan[PAN_FLAG_WEIGHT_WITH_COS];
  *sample_pan_circumsolar_var_red      = input.rte.mc.pan[PAN_FLAG_CIRCUMSOLAR_VAR_RED];
  *sample_pan_transmittance            = input.rte.mc.pan[PAN_FLAG_TRANSMITTANCE];

  *sample_pan_alpha = input.rte.mc.pan_alpha;

  /* user-defined surface normal */
  *sensordirection = input.rte.mc.sensordirection;
  if (*sensordirection) {
    const double temp     = sqrt (input.rte.mc.sensordirection_dx * input.rte.mc.sensordirection_dx +
                              input.rte.mc.sensordirection_dy * input.rte.mc.sensordirection_dy +
                              input.rte.mc.sensordirection_dz * input.rte.mc.sensordirection_dz);
    sensordirection_dx[0] = input.rte.mc.sensordirection_dx / temp;
    sensordirection_dx[1] = input.rte.mc.sensordirection_dy / temp;
    sensordirection_dx[2] = input.rte.mc.sensordirection_dz / temp;

    fprintf (stderr, " ... user-defined surface normal for surface irradiance calculations:\n");
    fprintf (stderr, "     (%f, %f, %f)\n", sensordirection_dx[0], sensordirection_dx[1], sensordirection_dx[2]);
  }
  return 0;
}

/***********************************************************************************/
/* Function:                                                              @62_30i@ */
/* Description:                                                                    */
/*  Setup a 3D atmosphere structure, to be used as input to MYSTIC.                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int setup_atmosphere3D (input_struct input, output_struct* output) {

#if HAVE_MYSTIC
  int    kc = 0, status = 0;
  int    is = 0, js = 0, isp = 0;
  int    nx_caoth = 0, ny_caoth = 0;
  double dx_caoth = 0.0, dy_caoth = 0.0;

  int counter = 0;

#ifdef AREA_DEBUG_N
  double        n0 = 0, n1 = 0, n2 = 0;
  double        area_numeric = 0, tmpstd = 0;
  photon_struct ptmp;
  int           i = 0, ie = 0, je = 0;
#endif

  FILE* f_area = NULL;

  output->atm.nthreed = 0;

  // initialize output->atm.Nzcld after redistribution;
  // we require it for non-mystic applications like dynamic_tenstream()
  output->atm.Nzcld = output->atm.nlyr;

  // BM, 22.1.2021: Moved that part here because the threed[] information is also required for dynamic_tenstream()

  /* determine the combined caoth 3D levels */
  output->atm.nthreed = 0;
  output->atm.threed  = calloc (output->atm.nlyr, sizeof (int));

  for (kc = 0; kc < output->atm.nlyr; kc++) {
    output->atm.threed[kc] = 0;

    for (isp = 0; isp < input.n_caoth; isp++)
      if (output->caoth3d[isp].nthreed > 0)
        output->atm.threed[kc] += output->caoth3d[isp].threed[kc];

    if (output->atm.threed[kc] > 0)
      output->atm.nthreed++;
  }

  if (input.make_all_layers_3d) {
    if (!input.quiet)
      fprintf (stderr, " ... setting all layers to be 3D, because option make_all_layers_3d\n");
    for (size_t k = 0; k < output->atm.nlyr; k++) {
      output->atm.threed[k]++;
    }
    output->atm.nthreed = vec_count_true (output->atm.nlyr, output->atm.threed);
  }

  /* for 3d boxairmass factors and 3D Jacobians, all layers should be treated 3d, even when there are no 3D clouds */
  if (input.rte.mc.boxairmass[DIM_3D] || input.rte.mc.jacobian[DIM_3D])
    for (kc = 0; kc < output->atm.nlyr; kc++) {
      output->atm.threed[kc] = 1;
      output->atm.nthreed    = output->atm.nlyr;
    }

  if (!input.quiet)
    fprintf (stderr,
             " ... found %d 3D layers (cloud inputs combined) on a %d x %d x %d grid\n",
             output->atm.nthreed,
             output->atm.Nxcld,
             output->atm.Nycld,
             output->atm.Nzcld);

  /* for ipa3D heating rates we need to calculate irradiances above and below each 3D layer;          */
  /* careful: for this application we require ONLY levels above and below 3D layers and no extra ones */
  if (input.atm.zout_source == OUTLEVEL_3D_LEVELS) {

    /* first count required levels */
    output->atm.nzout = 0;
    if (output->atm.threed[0]) /* bottom layer 3D? */
      output->atm.nzout = 2;

    for (int lc = 1; lc < output->atm.nlyr; lc++) {
      if (output->atm.threed[lc]) {
        if (output->atm.threed[lc - 1])
          output->atm.nzout++;
        else
          output->atm.nzout += 2;
      }
    }

    output->atm.zout_sea = calloc (output->atm.nzout, sizeof (float));
    output->atm.zout_sur = calloc (output->atm.nzout, sizeof (float));

    int lu = 0;
    if (output->atm.threed[0]) { /* bottom layer 3D? */
      output->atm.zout_sea[lu++] = output->atm.zd[output->atm.nlyr];
      output->atm.zout_sea[lu++] = output->atm.zd[output->atm.nlyr - 1];
    }

    for (int lc = 1; lc < output->atm.nlyr; lc++) {
      if (output->atm.threed[lc]) {
        if (output->atm.threed[lc - 1])
          output->atm.zout_sea[lu++] = output->atm.zd[output->atm.nlyr - lc - 1]; /* only top */
        else {
          output->atm.zout_sea[lu++] = output->atm.zd[output->atm.nlyr - lc];     /* bottom */
          output->atm.zout_sea[lu++] = output->atm.zd[output->atm.nlyr - lc - 1]; /* top    */
        }
      }
    }

    if (lu != output->atm.nzout) {
      CHKERROUT (-1, "Fatal error determining number of 3D levels");
    }

    for (lu = 0; lu < output->atm.nzout; lu++)
      output->atm.zout_sur[lu] = output->atm.zout_sea[lu] - output->alt.altitude;
  }

  // BM, 22.1.2021  CHECK IF WE REALLY NEED THE REST FOR IPA3D
  /* don't need the rest if no MYSTIC calculation */
  if (input.rte.solver != SOLVER_MONTECARLO && input.ipa3d == 0) /*ulrike 03.05.2010*/
    return 0;

  /****************************************/
  /* convert caoth altitudes from km to m */
  /****************************************/

  /* ulrike 3.5.2010: I use setup_atmosphere3D for ipa3d, but later on I do not need (?)
     the height levels zd in m. Thus, conversion from km->m should only be done
     in case of montecarlo */

  if (input.ipa3d != 1) {
    for (isp = 0; isp < input.n_caoth; isp++)
      if (output->caoth3d[isp].nthreed > 0)
        for (kc = 0; kc <= output->caoth3d[isp].nlyr; kc++)
          output->caoth3d[isp].zd[kc] *= 1000.0;
  }

  // Does it make sense to loop though isp if all are the same size anyway? Probably not.
  for (isp = 0; isp < input.n_caoth; isp++) {
    nx_caoth = output->caoth3d[isp].Nx;
    ny_caoth = output->caoth3d[isp].Ny;

    dx_caoth = output->caoth3d[isp].delX;
    dy_caoth = output->caoth3d[isp].delY;
  }

  /* load triangular surface data */
  if (strlen (input.rte.mc.filename[FN_MC_TRIANGULAR_SURFACE]) > 0) {
    int ierr;
    ierr =
      setup_triangular_surface (input.rte.mc.filename[FN_MC_TRIANGULAR_SURFACE], &(output->mc.triangular_surface), input.quiet);
    CHKERR (ierr);

    const t_triangular_surface* srfc = &(output->mc.triangular_surface);

    if (dx_caoth <= 0 && dy_caoth <= 0) {  // find vertices bounding box and set it as domain boundaries
      double vbox_min[3], vbox_max[3];
      ierr = get_global_vertices_bounding_box (srfc, vbox_min, vbox_max);
      CHKERR (ierr);
      dx_caoth = (vbox_max[0] + vbox_max[0] * FLT_EPSILON * 2);
      dy_caoth = (vbox_max[1] + vbox_max[1] * FLT_EPSILON * 2);

      nx_caoth = 1;
      ny_caoth = 1;
    }

    unsigned have_critical_triangles = 0;
    {  // check if we have axis aligned triangles
      for (size_t i = 0; i < srfc->N_triangles; ++i) {
        const size_t verts[3] = {srfc->triangles[i][0], srfc->triangles[i][1], srfc->triangles[i][2]};

        const double dimlimits[6] = {0, dx_caoth, 0, dy_caoth, 0, output->atm.zd[0] * 1e3};
        for (unsigned m = 0; m < 2; ++m) {    // check min or max
          for (unsigned k = 0; k < 3; ++k) {  // check for alignement along dimension k
            const int aligned[3] = {
              approx_rel (srfc->vertices[verts[0]][k], dimlimits[k * 2 + m], FLT_EPSILON),
              approx_rel (srfc->vertices[verts[1]][k], dimlimits[k * 2 + m], FLT_EPSILON),
              approx_rel (srfc->vertices[verts[2]][k], dimlimits[k * 2 + m], FLT_EPSILON),
            };
            if (vec_count_true (3, aligned) == 3) {
              fprintf (
                stderr,
                "Warning, triangle %lu is very close and axis aligned (%i %i %i) to the domain boundary in dimension %u (%e)!\n"
                " ... vertex 0: %e %e %e\n"
                " ... vertex 1: %e %e %e\n"
                " ... vertex 2: %e %e %e\n",
                i,
                aligned[0],
                aligned[1],
                aligned[2],
                k,
                dimlimits[k * 2 + m],
                srfc->vertices[verts[0]][0],
                srfc->vertices[verts[0]][1],
                srfc->vertices[verts[0]][2],
                srfc->vertices[verts[1]][0],
                srfc->vertices[verts[1]][1],
                srfc->vertices[verts[1]][2],
                srfc->vertices[verts[2]][0],
                srfc->vertices[verts[2]][1],
                srfc->vertices[verts[2]][2]);
              ++have_critical_triangles;
            }
          }
        }
      }
    }
    if (have_critical_triangles) {
      fprintf (stderr,
               " ... found %u axis aligned triangles at the border of the domain\n"
               " ... you have to be careful!\n"
               " ... Because we are jumping back and forth (cyclic boundaries), all the time,\n"
               " ... triangles in the vicinity of the domain boundary might be ignored!\n"
               " ... Try to move them a little bit.\n",
               have_critical_triangles);
    }
  } else {
    output->mc.triangular_surface.N_vertices              = 0;
    output->mc.triangular_surface.N_triangles             = 0;
    output->mc.triangular_surface.N_materials             = 0;
    output->mc.triangular_surface.vertices                = NULL;
    output->mc.triangular_surface.triangles               = NULL;
    output->mc.triangular_surface.material_of_triangle    = NULL;
    output->mc.triangular_surface.material_albedo         = NULL;
    output->mc.triangular_surface.temperature_of_triangle = NULL;
    output->mc.triangular_surface.Bplanck_of_triangle     = NULL;
    output->mc.triangular_surface.star_engine_tracer      = NULL;
    output->mc.triangular_surface.native_tracer           = NULL;
  }

  /********************************/
  /* first, setup the sample grid */
  /********************************/

  status = setup_sample_grid (output->atm.zd,
                              output->atm.nlev,
                              input.rte.mc.Nx_sample,
                              input.rte.mc.Ny_sample,
                              input.rte.mc.dx_sample,
                              input.rte.mc.dy_sample,
                              nx_caoth,
                              ny_caoth,
                              dx_caoth,
                              dy_caoth,
                              output->atm.zout_sur,
                              output->atm.nzout,
                              input.rte.mc.Nr,
                              input.rte.mc.Nt,
                              input.rte.mc.dt,
                              input.rte.umu,
                              input.rte.numu,
                              input.rte.mc.alpha,
                              input.rte.phi,
                              input.rte.nphi,
                              input.rte.mc.azimuth,
                              &(output->mc.sample),
                              input.rte.mc.escape,
                              input.rte.mc.vroom,
                              input.rte.mc.std,
                              input.rte.mc.std_value,
                              input.rte.mc.backward.yes,
                              input.rte.mc.backward.output,
                              input.rte.mc.backward.islower,
                              input.rte.mc.backward.jslower,
                              input.rte.mc.backward.isupper,
                              input.rte.mc.backward.jsupper,
                              input.rte.mc.backward.isstep,
                              input.rte.mc.backward.jsstep,
                              input.rte.mc.backward.writeallpixels,
                              input.rte.mc.backward.writeback,
                              input.rte.mc.refraction,
                              input.rte.mc.spectral_is,
                              input.rte.mc.boxairmass[DIM_1D],
                              input.rte.mc.boxairmass[DIM_3D],
                              input.rte.mc.spherical[DIM_1D],
                              input.rte.mc.spherical[DIM_3D],
                              input.rte.mc.spherical3D_scene,
                              input.rte.mc.spherical3D_scene_lon_min,
                              input.rte.mc.spherical3D_scene_lon_max,
                              input.rte.mc.spherical3D_scene_lat_min,
                              input.rte.mc.spherical3D_scene_lat_max,
                              input.rte.mc.DoLE,
                              input.rte.mc.ixmin,
                              input.rte.mc.ixmax,
                              input.rte.mc.iymin,
                              input.rte.mc.iymax,
                              input.ipa3d,
                              input.rte.mc.tipa,
                              input.rte.mc.polarisation,
                              input.rte.mc.polarisation_state,
                              input.rte.mc.nstokes,
                              input.rte.mc.bcond,
                              input.rte.mc.panorama,
                              input.rte.mc.panorama_forward,
                              input.rte.mc.filename[FN_MC_PANORAMA_ANGLES],
                              input.rte.mc.filename[FN_MC_BACKWARD_FILE],
                              input.rte.mc.maxscatters,
                              input.rte.mc.minscatters,
                              input.rte.mc.ncirc,
                              input.rte.mc.dcirc,
                              input.rte.mc.blitz_position,
                              input.rte.mc.locest,
                              input.rte.mc.jacobian[DIM_1D],
                              input.rte.mc.jacobian[DIM_3D],
                              input.rte.mc.jacobian_std,
                              input.rte.mc.delta_scaling_start,
                              input.rte.mc.filename[FN_MC_LIDAR],
                              input.rte.mc.reference_to_NN,
                              input.rte.mc.coherent_backscatter,
                              input.rte.mc.coherent_backscatter_lmode,
                              input.rte.mc.coherent_backscatter_off,
                              input.rte.mc.ris[MC_RIS_FACTOR],
                              input.rte.mc.ris[MC_RIS_OPTICAL_DEPTH],
                              input.rte.mc.sample_cldprp,
                              input.quiet);
  CHKERR (status);

  /* for memory allocation of the MYSTIC output passed */
  /* to uvspec we need to know which pixels are        */
  /* actually calculated because we allocate memory    */
  /* only for these                                    */
  if (input.rte.mc.backward.writeallpixels || (input.rte.mc.backward.islower == -999 && input.rte.mc.backward.isupper == -999 &&
                                               input.rte.mc.backward.jslower == -999 && input.rte.mc.backward.jsupper == -999)) {

    output->islower = 0;
    output->jslower = 0;
    output->isupper = output->mc.sample.Nx - 1;
    output->jsupper = output->mc.sample.Ny - 1;
    output->isstep  = 1;
    output->jsstep  = 1;
  } else {
    output->islower = input.rte.mc.backward.islower;
    output->isupper = input.rte.mc.backward.isupper;
    output->jslower = input.rte.mc.backward.jslower;
    output->jsupper = input.rte.mc.backward.jsupper;
    output->isstep  = input.rte.mc.backward.isstep;
    output->jsstep  = input.rte.mc.backward.jsstep;

    if (output->isupper >= output->mc.sample.Nx)
      output->isupper = output->mc.sample.Nx - 1;

    if (output->jsupper >= output->mc.sample.Ny)
      output->jsupper = output->mc.sample.Ny - 1;
  }

  if (output->islower < 0 || output->isupper < 0 || output->jslower < 0 || output->jsupper < 0) {
    fprintf (stderr, "Error with backward pixel definition in (line %d, function '%s', file '%s'):\n", ERROR_POSITION);
    fprintf (stderr,
             "islower = %d, isupper = %d, jslower = %d, jsupper = %d.\n",
             output->islower,
             output->isupper,
             output->jslower,
             output->jsupper);
    fprintf (stderr, "Please tell the programmer!\n");
    CHKERR (-1);
  }

  if (input.verbose)
    fprintf (stderr,
             " ... allocating memory for 3D fields, is = %d ... %d,  js = %d ... %d\n",
             output->islower,
             output->isupper,
             output->jslower,
             output->jsupper);

  status = set_sensor_position (input,
                                output->mc.sample.Nx,
                                output->mc.sample.Ny,
                                output->mc.sample.delX,
                                output->mc.sample.delY,
                                output->atm.zd[0] * 1e3,
                                &output->mc.sample.sensorposition,
                                output->mc.sample.sensorposition_x,
                                &output->mc.sample.sensordirection,
                                output->mc.sample.sensordirection_dx,
                                &output->mc.sample.panorama,
                                &output->mc.sample.panorama_forward,
                                &output->mc.sample.pan_no_pixel,
                                &output->mc.sample.pan_quicklook,
                                &output->mc.sample.pan_distr_photons_over_pixel,
                                &output->mc.sample.pan_with_direct_rad,
                                &output->mc.sample.pan_weight_with_cos,
                                &output->mc.sample.pan_circumsolar_var_red,
                                &output->mc.sample.pan_transmittance,
                                &output->mc.sample.pan_alignment,
                                &output->mc.sample.pan_alpha,
                                &output->mc.sample.pan_phi_min,
                                &output->mc.sample.pan_phi_max,
                                &output->mc.sample.pan_theta_min,
                                &output->mc.sample.pan_theta_max,
                                &output->mc.sample.pan_umu_min,
                                &output->mc.sample.pan_umu_max,
                                output->mc.sample.pan_phi_from_file,
                                output->mc.sample.pan_theta_from_file,
                                output->mc.sample.pan_n_from_file,
                                &output->mc.sample.align_theta,
                                &output->mc.sample.align_phi,
                                &output->mc.sample.pan_roll,
                                &output->mc.sample.pan_pitch,
                                &output->mc.sample.pan_yaw,
                                &output->mc.sample.add_phi_when_aligning);
  CHKERR (status);

  /* read 2D surface elevation */
  /* read 2D elevation data */

  if (strlen (input.rte.mc.filename[FN_MC_ELEVATION]) > 0) {

    /* read elevation data from file */
    status = setup_elevation2D (input.rte.mc.filename[FN_MC_ELEVATION],
                                output->mc.sample.Nx,
                                output->mc.sample.Ny,
                                output->mc.sample.delX,
                                output->mc.sample.delY,
                                &(output->mc.elev),
                                input.rte.mc.visualize,
                                input.quiet);
    CHKERR (status);

    output->mc.elev.elev2D = 1;

    /* surface-parallel or horizontal irradiance */
    output->mc.sample.surfaceparallel = input.rte.mc.surfaceparallel;

    if (!input.quiet) {
      if (output->mc.sample.surfaceparallel)
        fprintf (stderr, " ... calculating surface parallel irradiance\n");
      else
        fprintf (stderr, " ... calculating horizontal irradiance\n");
    }
  }

  /* calculate surface area for each sample pixel in case of 2D elevation and surface parallel irradiance */
  if (output->mc.elev.elev2D && output->mc.sample.surfaceparallel) {

    output->mc.sample.surface_area         = calloc (output->mc.sample.Nx, sizeof (double*));
    output->mc.sample.surface_area_std     = calloc (output->mc.sample.Nx, sizeof (double*));
    output->mc.sample.surface_area_counter = calloc (output->mc.sample.Nx, sizeof (int*));
    for (is = 0; is < output->mc.sample.Nx; is++) {
      output->mc.sample.surface_area[is]         = calloc (output->mc.sample.Ny, sizeof (double));
      output->mc.sample.surface_area_std[is]     = calloc (output->mc.sample.Ny, sizeof (double));
      output->mc.sample.surface_area_counter[is] = calloc (output->mc.sample.Ny, sizeof (int));
    }

    f_area = fopen ("area2D.dat", "w");
    if (!input.quiet) {
      fprintf (stderr, " ... calculating surface area of sample pixels for surface-parallel irradiance  ");
    }

    /* calculate surface area of each sample pixel */
    for (is = 0; is < output->mc.sample.Nx; is++)
      for (js = 0; js < output->mc.sample.Ny; js++) {

        output->mc.sample.surface_area_std[is][js]     = 0;
        output->mc.sample.surface_area_counter[is][js] = 0;

        counter++;

#ifdef AREA_DEBUG_N

        /* visualize progress to proof that mystic is still working */
        if (counter % 4 == 0)
          fprintf (stderr, "\b|");
        else if (counter % 4 == 1)
          fprintf (stderr, "\b/");
        else if (counter % 4 == 2)
          fprintf (stderr, "\b-");
        else if (counter % 4 == 3)
          fprintf (stderr, "\b\\");

        tmpstd = 0;

        /* numerical calculation via Monte Carlo; this is not used anymore because we now have the */
        /* analytical solution; but it is left there because it might be useful for future testing */
        for (i = 0; i < AREA_DEBUG_N; i++) {
          ptmp.x[0] = ((double) is + uvspec_random()) * output->mc.sample.delX;
          ptmp.x[1] = ((double) js + uvspec_random()) * output->mc.sample.delY;

          /* elevation coordinates */
          elev_coord (&ptmp, &(output->mc.elev), &ie, &je);

          /* upward normal to the slant surface */
          n0 = -output->mc.elev.surf[ie][je].a - output->mc.elev.surf[ie][je].c * (ptmp.x[1] - (double) je * output->mc.elev.delY);
          n1 = -output->mc.elev.surf[ie][je].b - output->mc.elev.surf[ie][je].c * (ptmp.x[0] - (double) ie * output->mc.elev.delX);
          n2 = 1;

          /* normalize n */
          temp = sqrt (n0 * n0 + n1 * n1 + n2 * n2);
          n2 /= temp;

          output->mc.sample.surface_area[is][js] += 1.0 / fabs (n2);
          output->mc.sample.surface_area_std[is][js] += 1.0 / n2 / n2;
        }

        output->mc.sample.surface_area[is][js] /= (double) AREA_DEBUG_N;
        output->mc.sample.surface_area_std[is][js] /= (double) AREA_DEBUG_N;
        output->mc.sample.surface_area_counter[is][js] = AREA_DEBUG_N;

        output->mc.sample.surface_area_std[is][js] =
          sqrt (output->mc.sample.surface_area_std[is][js] -
                output->mc.sample.surface_area[is][js] * output->mc.sample.surface_area[is][js]) /
          sqrt ((double) AREA_DEBUG_N);

        /* store numeric area */
        area_numeric = output->mc.sample.surface_area[is][js];

#endif

        /* calculate surface area analytically */
        output->mc.sample.surface_area[is][js] = calc_elevation_area (&(output->mc.elev),
                                                                      (double) is * output->mc.sample.delX,
                                                                      (double) js * output->mc.sample.delY,
                                                                      (double) (is + 1) * output->mc.sample.delX,
                                                                      (double) (js + 1) * output->mc.sample.delY);

#ifdef AREA_DEBUG_N
        fprintf (f_area,
                 "%3d %3d %.6e %.6e  %.6e\n",
                 is,
                 js,
                 output->mc.sample.surface_area[is][js],
                 area_numeric,
                 output->mc.sample.surface_area_std[is][js]);
#else
        fprintf (f_area, "%3d %3d %.6e\n", is, js, output->mc.sample.surface_area[is][js]);
#endif
      }

    fclose (f_area);
    if (!input.quiet) {
      fprintf (stderr, "\b\b\n");
    }
  }

  /* surface temperature: either defined using surface_temperature */
  /* or with a 2D temperature file. Otherwise the atmospheric      */
  /* temperature at the respective altitude is used by mystic()    */
  output->mc.surftemp.btemp = output->surface_temperature;

  /* flag if user defined surface_temperature */
  output->mc.surftemp.user_defined_btemp = 0;
  if (input.surface_temperature >= 0)
    output->mc.surftemp.user_defined_btemp = 1;

  if (strlen (input.rte.mc.filename[FN_MC_TEMPERATURE]) > 0) {

    /* read 2D surface temperature data from file */
    status = setup_temperature2D (input.rte.mc.filename[FN_MC_TEMPERATURE],
                                  output->mc.sample.Nx,
                                  output->mc.sample.Ny,
                                  output->mc.sample.delX,
                                  output->mc.sample.delY,
                                  &(output->mc.surftemp),
                                  input.rte.mc.ixmin,
                                  input.rte.mc.ixmax,
                                  input.rte.mc.iymin,
                                  input.rte.mc.iymax,
                                  input.quiet);
    CHKERR (status);

    output->mc.surftemp.surf2D = 1;
  }

  /* copy 3D caoth structures to atmosphere structure */
  for (isp = 0; isp < input.n_caoth; isp++) {
    if (output->caoth3d[isp].nlyr > 0) {
      output->atm.Nxcld = output->caoth3d[isp].Nx;
      output->atm.Nycld = output->caoth3d[isp].Ny;
      output->atm.Nzcld = output->caoth3d[isp].nlyr;
      output->atm.dxcld = output->caoth3d[isp].delX;
      output->atm.dycld = output->caoth3d[isp].delY;
    }
  }

  { /* if no 3D caoth defined, copy sampling area size to caoth area size */
    int any_3d_caoth = 0;
    for (isp = 0; isp < input.n_caoth; isp++)
      if (output->caoth3d[isp].nthreed > 0)
        any_3d_caoth++;

    if (!any_3d_caoth) {
      output->atm.Nxcld = 1;
      output->atm.Nycld = 1;
      output->atm.dxcld = (double) (output->mc.sample.Nx) * (output->mc.sample.delX);
      output->atm.dycld = (double) (output->mc.sample.Ny) * (output->mc.sample.delY);
      for (isp = 0; isp < input.n_caoth; isp++) {
        output->caoth3d[isp].Nx   = 1;
        output->caoth3d[isp].Ny   = 1;
        output->caoth3d[isp].delX = (double) (output->mc.sample.Nx) * (output->mc.sample.delX);
        output->caoth3d[isp].delY = (double) (output->mc.sample.Ny) * (output->mc.sample.delY);
      }
    }
  }

  /* check if caoth grid and sample grids are compatible */
  for (isp = 0; isp < input.n_caoth; isp++) {
    if (output->caoth3d[isp].nthreed > 0) {
      if (!equal ((double) (output->mc.sample.Nx) * (output->mc.sample.delX),
                  (double) output->caoth3d[isp].Nx * output->caoth3d[isp].delX)) {
        fprintf (stderr,
                 "Error, x-size of sample domain (%.10e) does not equal x-size of %s domain (%.10e)\n",
                 (double) (output->mc.sample.Nx) * (output->mc.sample.delX),
                 input.caoth[isp].fullname,
                 (double) output->caoth3d[isp].Nx * output->caoth3d[isp].delX);
        CHKERR (-1);
      }

      if (!equal ((double) (output->mc.sample.Ny) * (output->mc.sample.delY),
                  (double) output->caoth3d[isp].Ny * output->caoth3d[isp].delY)) {
        fprintf (stderr,
                 "Error, y-size of sample domain (%.10e) does not equal y-size of %s domain (%.10e)\n",
                 (double) (output->mc.sample.Ny) * (output->mc.sample.delY),
                 input.caoth[isp].fullname,
                 (double) output->caoth3d[isp].Ny * output->caoth3d[isp].delY);
        fprintf (stderr,
                 "Ny_sample = %d delY_sample = %f Ny_caoth = %d delY_caoth = %f\n",
                 (output->mc.sample.Ny),
                 (output->mc.sample.delY),
                 output->caoth3d[isp].Ny,
                 output->caoth3d[isp].delY);
        CHKERR (-1);
      }
    }
  }

  // Changing layers to be 3D layers for all layers where a surface triangle intersects.
  // All those layers have to be 3D because of the way we check for ray <--> triangle intersections
  // which is tricky in case of cyclic boundary conditions. If the layers are 3D,
  // the photons are actually registered when crossing the horizontal boundary,
  // which is needed for correct intersection distances.
  // Conceptually we determine if a 1D layers bounding box intersects any triangle,
  // if so we flag this layer as 3D.
  if (output->mc.triangular_surface.triangles) {
    t_triangular_surface* srfc = &(output->mc.triangular_surface);

    float heightlvls_rev[output->atm.nlyr + 1];
    for (size_t i = 0; i < output->atm.nlyr + 1; ++i)
      heightlvls_rev[i] = output->atm.zd[i];
    int ierr = reverse_1d_f4 (output->atm.nlyr + 1, heightlvls_rev);
    CHKERR (ierr);

    for (kc = 0; kc < output->atm.nlyr; kc++) {
      if (!output->atm.threed[kc]) {  // if not yet 3D

        double boxcenter[3];
        double boxhalfsize[3];
        int    ierr = cellid2boundingbox (kc,
                                       -1,
                                       -1,
                                       output->atm.Nxcld,
                                       output->atm.Nycld,
                                       0,
                                       output->atm.dxcld,
                                       output->atm.dycld,
                                       heightlvls_rev,
                                       boxcenter,
                                       boxhalfsize);
        CHKERR (ierr);

        for (size_t iface = 0; iface < srfc->N_triangles; ++iface) {
          int has_intersection = check_triangle_boundingbox_overlap (iface, srfc, boxcenter, boxhalfsize);
          if (has_intersection) {
            output->atm.threed[kc]++;
            output->atm.nthreed++;
            break;
          }
        }
      }
    }
    // Nzcld has to be consistent with nlyr if there shall be any 3D layers
    output->atm.Nzcld = output->atm.nlyr;

    if (0) {  // For testing purposes: also set neighbouring layers to 3D
      int tmp[output->atm.nlyr];
      for (size_t k = 0; k < output->atm.nlyr; k++)
        tmp[k] = output->atm.threed[k];

      for (size_t k = 0; k < output->atm.nlyr; k++) {
        if (tmp[k]) {
          output->atm.threed[MAX (0, k - 1)]++;
          output->atm.threed[MIN (output->atm.nlyr - 1, k + 1)]++;
        }
      }
      output->atm.nthreed = vec_count_true (output->atm.nlyr, output->atm.threed);
      for (size_t k = 0; k < output->atm.nlyr; k++) {
        fprintf (stderr,
                 "layer %lu is 3D? %d -- height [%f %f]\n",
                 k,
                 output->atm.threed[k],
                 heightlvls_rev[k],
                 heightlvls_rev[k + 1]);
      }
    }

    // in case that there is no 3D cloud provided,
    // lets split the atmosphere so that on average,
    // there are only 100 triangles in each voxel
#if !HAVE_STAR_ENGINE
    if (1) {
      if (output->atm.Nxcld == 1 && output->atm.Nycld == 1) {
        const int Nx = MIN (100, 1 + srfc->N_triangles / (2 * 100 * output->atm.nthreed + 1));
        const int Ny = MIN (100, 1 + srfc->N_triangles / (2 * 100 * output->atm.nthreed + 1));

        output->atm.dxcld = output->atm.dxcld * output->atm.Nxcld / Nx;
        output->atm.dycld = output->atm.dycld * output->atm.Nycld / Ny;
        output->atm.Nxcld = Nx;
        output->atm.Nycld = Ny;
      }
    }
#endif

    // Then build acceleration structures with the enlarged threed flags.
    // ATTENTION: Note that at this point the mesh which is used in the solver
    //            has to be the same as it is stated here!
    //            Otherwise the acceleration structure is probably not working correctly anymore.
    //            You however, destroy it and build it anew...
    ierr = init_native_tracer (output->atm.Nxcld,
                               output->atm.Nycld,
                               output->atm.nlyr,
                               output->atm.dxcld,
                               output->atm.dycld,
                               output->atm.zd,
                               output->atm.threed,
                               srfc,
                               input.quiet);
    CHKERR (ierr);
  }  // end of triangle setup

  if (!input.quiet)
    fprintf (stderr,
             " ... found %d 3D layers (all input combined) on a %d x %d x %d grid\n",
             output->atm.nthreed,
             output->atm.Nxcld,
             output->atm.Nycld,
             output->atm.nlyr);

#endif /* just do nothing if MYSTIC not included */
  return 0;
}

/***********************************************************************************/
/* Function: compare_levels_df                                            @62_30i@ */
/* Description:                                                                    */
/*  Compare two altitude grids. First is expected to be double, second float       */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

/*
  static int compare_levels_df (double *z1, int n1,
  float *z2, int n2)
  {
  int kc=0;

  if (n1!=n2) {
  fprintf (stderr, "Error, number of layers differing\n");
  return -1;
  }

  for (kc=0; kc<=n1; kc++)
  if (fabs((float) z1[kc] - z2[kc]) > MC_EPSILON * z1[kc]) {
  fprintf (stderr, "Error, difference at level %d, altitude %g vs. %g\n",	
  kc, z1[kc], z2[kc]);
  return -1;
  }

  return 0;
  }
*/

/* Read 3D caoth data from file */
static int read_caoth3D_file (char*            filename,
                              input_struct     input,
                              caoth_inp_struct input_caoth,
                              output_struct*   output,
                              /* Output */
                              caoth_out_struct*   output_caoth,
                              caoth3d_out_struct* caoth3d,
                              int*                nmom) {
  int kc = 0, status = 0;
  int ix = 0, iy = 0;

  void F77_FUNC (wcloud, WCLOUD) (float* lambda_r,
                                  int*   newsiz,
                                  int*   nlyr,
                                  char*  filepath,
                                  int*   nstring,
                                  float* lwc,
                                  float* wceffr,
                                  float* tmp_wc_dtau,
                                  float* tmp_wc_gg,
                                  float* tmp_wc_ssa,
                                  float* zd,
                                  int*   wclyr);

  /* read caoth_file */
  status = read_3D_caoth (filename,
                          &(caoth3d->Nx),
                          &(caoth3d->Ny),
                          &(caoth3d->nlyr),
                          &(caoth3d->delX),
                          &(caoth3d->delY),
                          &(caoth3d->zd),
                          &(caoth3d->lwc),
                          &(caoth3d->reff),
                          &(caoth3d->cf),
                          &(caoth3d->ext),
                          &(caoth3d->g1),
                          &(caoth3d->g2),
                          &(caoth3d->ff),
                          &(caoth3d->f),
                          &(caoth3d->ssa),
                          &(caoth3d->dscale),
                          &(caoth3d->reffmin),
                          &(caoth3d->reffmax),
                          &(caoth3d->cldproperties),
                          &(caoth3d->threed),
                          input.rte.mc.ixmin,
                          input.rte.mc.ixmax,
                          input.rte.mc.iymin,
                          input.rte.mc.iymax,
                          input.quiet);
  if (status)
    return fct_err_out (status, "read_3D_caoth", ERROR_POSITION);

  if (!input.quiet)
    fprintf (stderr, " ... read %d x %d x %d %s data set\n", caoth3d->Nx, caoth3d->Ny, caoth3d->nlyr, caoth3d->fullname);

  /* we got the 3D caoth structure in caoth3d = output->wc.threed now; */
  /* this structure is passed to mystic()                              */

  for (kc = 0; kc < caoth3d->nlyr; kc++)
    if (caoth3d->threed[kc])
      caoth3d->nthreed++;

  if (!input.quiet) {
    switch (caoth3d->cldproperties) {
    case CLD_OPTPROP:
      fprintf (stderr, " ... read extinction coefficient, asymmetry parameter,\n");
      fprintf (stderr, " ... and single scattering albedo 3D for %sfile\n", caoth3d->fullname);
      break;
    case CLD_EXTREFF:
      fprintf (stderr, " ... read extinction coefficient and effective droplet radius from 3D file for %s\n", caoth3d->fullname);
      break;
    case CLD_LWCREFF:
      fprintf (stderr, " ... read liquid water content and effective droplet radius from 3D file for %s\n", caoth3d->fullname);
      break;
    case CLD_LWCREFFCF:
      fprintf (stderr,
               " ... read liquid water content, effective droplet radius, and cloud fraction from 3D file for %s\n",
               caoth3d->fullname);
      break;
    default:
      fprintf (stderr, "Error, unknown properties %d for %s\n", caoth3d->cldproperties, caoth3d->fullname);
      CHKERR (-1);
    }
  }

  if (!input.quiet)
    fprintf (stderr, " ... found %d 3D layers in %s\n", caoth3d->nthreed, caoth3d->fullname);

  if (!input.quiet) {
    for (kc = 0; kc < caoth3d->nlyr; kc++)
      if (caoth3d->threed[kc])
        fprintf (stderr,
                 " ... found 3D %s at level %d, altitude range %g - %g km\n",
                 caoth3d->fullname,
                 kc,
                 caoth3d->zd[kc],
                 caoth3d->zd[kc + 1]);
  }

  /* no need to load properties, there is no caoth present; probably
     the user wanted a dummy 3D grid */
  if (caoth3d->nthreed == 0)
    return 0;

  /* default: Henyey-Greenstein - flag needed by MYSTIC */
  caoth3d->nonHG = 0;

  switch (input_caoth.properties) {
  case PROP_HU:
  case PROP_ECHAM4:
  case PROP_YANG:
  case PROP_KEY:
    break;

  case PROP_MIE:
  case PROP_FILE:
  case PROP_IC_MIE:
  case PROP_BAUM:
  case PROP_BAUM_V36:
  case PROP_HEY:
  case PROP_YANG2013:

    if (caoth3d->cldproperties != CLD_OPTPROP) { /* no need to read optical properies if HG anyway */

      /* need to read caoth property file only if not yet already done */
      if (!(output_caoth->ssprop.alloc && output_caoth->ssprop.alloc_moments)) {
        /* Read optical properties file */
        status = read_caoth_prop (input_caoth.properties_filename,
                                  output->wl.lambda_r,
                                  output->wl.nlambda_r,
                                  output->wl.iv_lower,
                                  output->wl.iv_upper,
                                  input_caoth.interpolate,
                                  0,
                                  output_caoth,
                                  input_caoth.properties,
                                  NULL,
                                  NULL,
                                  0,
                                  nmom,
                                  input.rte.pol_nstokes,
                                  input.rte.nstr,
                                  input.rte.solver,
                                  input.rte.disort_icm,
                                  input.verbose,
                                  input.quiet);
        if (status != 0) {
          fprintf (stderr,
                   "Error %d reading property file %s for %s\n",
                   status,
                   input_caoth.properties_filename,
                   input_caoth.fullname);
          return status;
        }

        /* set flag for non-Henyey-Greenstein, that is, tabulated phase function */
        caoth3d->nonHG = 1;
      } else {
        fprintf (stderr,
                 " ... no need to read property file %s again for %s\n",
                 input_caoth.properties_filename,
                 input_caoth.fullname);
      }
    } else {
      fprintf (stderr,
               " ... ignoring property file %s for %s, using Henyey-Greenstein\n",
               input_caoth.properties_filename,
               input_caoth.fullname);
    }
    break;

  case PROP_FU:

    /* check if users wants Fu or Yang definition of effective radius */

    if (input_caoth.fu2yang) {
      if (!input.quiet)
        fprintf (stderr, " ... converting from Fu (1996/98) to Key et al. (2002) definition of effective radius \n");
      if (caoth3d->nthreed > 0) /* if at least one 3D layer */
        for (kc = 0; kc < caoth3d->nlyr; kc++)
          if (caoth3d->threed[kc]) /* 3D layer */
            for (ix = 0; ix < caoth3d->Nx; ix++)
              for (iy = 0; iy < caoth3d->Ny; iy++) {
                if (caoth3d->reff[kc][ix][iy] > 0.0) {
                  if (input.verbose)
                    fprintf (stderr, "fu2yang: scaling effective radius from %f to", caoth3d->reff[kc][ix][iy]);
                  caoth3d->reff[kc][ix][iy] /= (3.0 * sqrt (3.0) / 4.0);
                  if (input.verbose)
                    fprintf (stderr, " %f\n", caoth3d->reff[kc][ix][iy]);
                }
              }

      caoth3d->reffmin /= (3.0 * sqrt (3.0) / 4.0);
      caoth3d->reffmax /= (3.0 * sqrt (3.0) / 4.0);
    }

    break;

  case PROP_RAYTRACING:
    if (strncasecmp (output_caoth->name, "ic", 2) == 0) {
      status = read_raytracing_file (input.filename[FN_RAYTRACING],
                                     &(output_caoth->raytracing_prop),
                                     &(output_caoth->n_raytracing_prop),
                                     input.quiet);

      if (status != 0) {
        fprintf (stderr, "Error %d reading %s\n", status, input.filename[FN_RAYTRACING]);
        return status;
      }
    }
    break;

  default:
    fprintf (stderr, "Error, unknown property %d for %s. This is not\n", input_caoth.properties, input_caoth.fullname);
    fprintf (stderr, "supposed to happen and indicates a coding error.\n");

    CHKERR (-1);
    break;
  }

  /* scale tau RPB */
  if (input_caoth.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE] >= 0.0)
    if (caoth3d->nthreed > 0) /* if at least one 3D layer */
      for (kc = 0; kc < caoth3d->nlyr; kc++)
        if (caoth3d->threed[kc]) /* 3D layer */
          for (ix = 0; ix < caoth3d->Nx; ix++)
            for (iy = 0; iy < caoth3d->Ny; iy++)
              caoth3d->lwc[kc][ix][iy] *= input_caoth.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE];

  return 0;
}

/********************************************************************/
/* Convert 1D and 3D caoth data to MYSTIC optical properties        */
/* for a given wavelength index iv. In particular, need to          */
/* extinction, single scattering albedo, and delta-scaling factor   */
/* to the effective radius of each 3D grid box. Also, need to       */
/* calculate the cumulative scattering function table for           */
/* the respective wavelength.                                       */
/********************************************************************/

int convert_caoth_mystic (input_struct        input,
                          caoth_inp_struct    input_caoth,
                          output_struct*      output,
                          caoth_out_struct*   output_caoth,
                          caoth3d_out_struct* caoth3d,
                          int                 iv) {
  int status  = 0;
  int nstring = (int) strlen (input.filename[FN_PATH]);

  float   lay_lwc[2]  = {0, 0};
  float   lay_reff[2] = {0, 0};
  float   lay_dtau[2] = {0, 0};
  float   lay_g1[2]   = {0, 0};
  float   lay_g2[2]   = {0, 0};
  float   lay_f[2]    = {0, 0};
  float   lay_ff[2]   = {0, 0};
  float   lay_ssa[2]  = {0, 0};
  float   lay_zd[2]   = {1, 0};
  float** wvl_dtau    = NULL;
  float** wvl_g1      = NULL;
  float** wvl_ssa     = NULL;

  float rmin = 0, rmax = 0;
  int   nlyr = 1, nlev = 2;
  int   newsiz   = 1;
  int   lyr_flag = 1;
  int   ix = 0, iy = 0, iz = 0;

  int newkey      = 0;
  int caoth_habit = 0;

  int    i1     = 0;
  double extinc = 0, scattr = 0, omega = 0, dscale = 0, delta_r = 0, f = 0, gg = 0;
  double ssa_us = 0.0;

  char  path[FILENAME_MAX] = "";
  float wavelength         = 0;

  float rho_ice = 0.917;

  void F77_FUNC (wcloud, WCLOUD) (float* lambda_r,
                                  int*   newsiz,
                                  int*   nlyr,
                                  char*  filepath,
                                  int*   nstring,
                                  float* lwc,
                                  float* wceffr,
                                  float* tmp_wc_dtau,
                                  float* tmp_wc_gg,
                                  float* tmp_wc_ssa,
                                  float* zd,
                                  int*   wclyr);

  status = ASCII_calloc_float (&wvl_dtau, 18, 2);
  CHKERR (status);

  status = ASCII_calloc_float (&wvl_g1, 18, 2);
  CHKERR (status);

  status = ASCII_calloc_float (&wvl_ssa, 18, 2);
  CHKERR (status);

  switch (caoth3d->cldproperties) {
  case CLD_OPTPROP:
    break;

  case CLD_EXTREFF:
  case CLD_LWCREFF:
  case CLD_LWCREFFCF:

    /* first create phase function table */
    switch (input_caoth.properties) {
    case PROP_HU:
    case PROP_ECHAM4:
    case PROP_FU:
    case PROP_RAYTRACING:
      break;

    case PROP_KEY:
      newkey      = 0;
      caoth_habit = input_caoth.habit;
      break;

    case PROP_YANG:
      newkey      = 1;
      caoth_habit = input_caoth.habit;
      break;

    case PROP_MIE:
    case PROP_FILE:
    case PROP_IC_MIE:
    case PROP_BAUM:
    case PROP_BAUM_V36:
    case PROP_HEY:
    case PROP_YANG2013:

      /* combine 1D and 3D rmin/rmax */
      rmax = (caoth3d->reffmax > output_caoth->microphys.effrmax ? caoth3d->reffmax : output_caoth->microphys.effrmax);

      /* take the smaller of both which is larger than zero */
      if (caoth3d->reffmin > 0.0) {
        if (output_caoth->microphys.effrmin > 0.0)
          rmin = (caoth3d->reffmin < output_caoth->microphys.effrmin ? caoth3d->reffmin : output_caoth->microphys.effrmin);
        else
          rmin = caoth3d->reffmin;
      } else
        rmin = output_caoth->microphys.effrmin;

      if (rmin <= 0 || rmin == FLT_MAX) /* no need to do anything */
        break;

      if (!input.quiet)
        fprintf (stderr, " ... generating phase function tables between rmin=%f and rmax=%f\n", rmin, rmax);

      /* generate phase function table */
      status = setup_Legendre_table (&(caoth3d->phase),
                                     &(output_caoth->ssprop),
                                     iv,
                                     rmin,
                                     rmax,
                                     input.rte.mc.delta_scaling_mucut,
                                     input.rte.mc.truncate,
                                     input.quiet);
      CHKERR (status);
      break;

    default:
      fprintf (stderr, "Error, unknown property %d for %s. This is not\n", input_caoth.properties, input_caoth.fullname);
      fprintf (stderr, "supposed to happen and indicates a coding error.\n");

      CHKERR (-1);
      break;
    }

    for (iz = 0; iz < caoth3d->nlyr; iz++) {
      if (caoth3d->nthreed > 0) {
        if (caoth3d->threed[iz]) { /* 3D layer */

          for (ix = 0; ix < caoth3d->Nx; ix++) {
            for (iy = 0; iy < caoth3d->Ny; iy++) {

              if (caoth3d->ext[iz][ix][iy] > 0.0 || caoth3d->lwc[iz][ix][iy] > 0.0) {

                /* use dummy LWC - it doesn't matter anyway */
                if (caoth3d->cldproperties == CLD_EXTREFF) {
                  lay_lwc[0] = 1.0;
                  lay_lwc[1] = 1.0;
                }

                if (caoth3d->cldproperties == CLD_LWCREFF || caoth3d->cldproperties == CLD_LWCREFFCF) {
                  lay_lwc[0] = caoth3d->lwc[iz][ix][iy];
                  lay_lwc[1] = caoth3d->lwc[iz][ix][iy];
                }

                lay_reff[0] = caoth3d->reff[iz][ix][iy];
                lay_reff[1] = caoth3d->reff[iz][ix][iy];

                switch (input_caoth.properties) {
                case PROP_HU:

                  switch (input.ck_scheme) {
                  case CK_FU:

                    /* call Fu and Liou code by Fred Rose */
                    status = ckdfucld (lay_zd, lay_reff, lay_lwc, nlev, iv, wvl_dtau, wvl_g1, wvl_ssa);
                    CHKERR (status);

                    /* write results to final destination */

                    if (caoth3d->cldproperties == CLD_LWCREFF || caoth3d->cldproperties == CLD_LWCREFFCF)
                      caoth3d->ext[iz][ix][iy] = wvl_dtau[iv][0] / 1000.0;

                    caoth3d->ssa[iz][ix][iy] = wvl_ssa[iv][0];
                    caoth3d->g1[iz][ix][iy]  = wvl_g1[iv][0];
                    caoth3d->g2[iz][ix][iy]  = 0;
                    caoth3d->ff[iz][ix][iy]  = 1.0;

                    break;
                  default: /* old uvspec wcloud.f for Hu and Stamnes */

                    /* wcloud calculates lwc_layer and reff_layer internally */
                    F77_FUNC (wcloud, WCLOUD)
                    (&output->wl.lambda_r[iv],
                     &newsiz,
                     &nlyr,
                     input.filename[FN_PATH],
                     &nstring,
                     lay_lwc,
                     lay_reff,
                     lay_dtau,
                     lay_g1,
                     lay_ssa,
                     lay_zd,
                     &lyr_flag);

                    /* write results to final destination */

                    /* need to convert again from 1/km to 1/m */
                    /* because wcloud() thinks in km rather   */
                    /* than m                                 */

                    if (caoth3d->cldproperties == CLD_LWCREFF || caoth3d->cldproperties == CLD_LWCREFFCF)
                      caoth3d->ext[iz][ix][iy] = lay_dtau[0] / 1000.0;

                    caoth3d->ssa[iz][ix][iy] = lay_ssa[0];
                    caoth3d->g1[iz][ix][iy]  = lay_g1[0];
                    caoth3d->g2[iz][ix][iy]  = 0;
                    caoth3d->ff[iz][ix][iy]  = 1.0;

                    /*
		      fprintf (stderr, "WCLOUD3D %d %d %d  %f %f %f %f\n",
		      iz, ix, iy,
		      caoth3d->ssa[iz][ix][iy],
		      caoth3d->g1 [iz][ix][iy],
		      caoth3d->g2 [iz][ix][iy],
		      caoth3d->ff [iz][ix][iy]);
		    */

                    break;
                  }

                  break;

                case PROP_FU:

                  /* calculate caoth optical properties */
                  switch (input.ck_scheme) {
                  case CK_FU:

                    /* call Fu and Liou code by Fred Rose */
                    status = ckdfuice (lay_zd, lay_reff, lay_lwc, nlev, iv, wvl_dtau, wvl_g1, wvl_ssa, input_caoth.unscaled);
                    CHKERR (status);

                    /* write results to final destination */
                    if (caoth3d->cldproperties == CLD_LWCREFF || caoth3d->cldproperties == CLD_LWCREFFCF)
                      caoth3d->ext[iz][ix][iy] = wvl_dtau[iv][0] / 1000.;

                    caoth3d->ssa[iz][ix][iy] = wvl_ssa[iv][0];
                    caoth3d->g1[iz][ix][iy]  = wvl_g1[iv][0];
                    caoth3d->g2[iz][ix][iy]  = 0;
                    caoth3d->ff[iz][ix][iy]  = 1.0;
                    break;

                  default:
                    if (output->wl.lambda_r[iv] < 4000.0) {

                      status = ic_fu96 (output->wl.lambda_r[iv],
                                        nlev - 1,
                                        input.filename[FN_PATH],
                                        lay_lwc,
                                        lay_reff,
                                        lay_dtau,
                                        lay_g1,
                                        lay_ssa,
                                        lay_f,
                                        lay_zd,
                                        lyr_flag,
                                        input_caoth.unscaled);
                      CHKERR (status);
                    } else {
                      status = ic_fu98 (output->wl.lambda_r[iv],
                                        nlev - 1,
                                        input.filename[FN_PATH],
                                        lay_lwc,
                                        lay_reff,
                                        lay_dtau,
                                        lay_g1,
                                        lay_ssa,
                                        lay_zd,
                                        lyr_flag);
                      CHKERR (status);
                    }

                    /* write results to final destination */
                    if (caoth3d->cldproperties == CLD_LWCREFF || caoth3d->cldproperties == CLD_LWCREFFCF)
                      caoth3d->ext[iz][ix][iy] = lay_dtau[0] / 1000.;

                    caoth3d->ssa[iz][ix][iy] = lay_ssa[0];
                    caoth3d->g1[iz][ix][iy]  = lay_g1[0];
                    caoth3d->g2[iz][ix][iy]  = 0;
                    caoth3d->ff[iz][ix][iy]  = 1.0;
                    caoth3d->f[iz][ix][iy]   = lay_f[0];

                  } /* end switch (input.ck_scheme) */

                  break;

                case PROP_KEY:
                case PROP_YANG:

                  /* call yang routine to obtain optical properties */
                  lay_reff[0] = caoth3d->reff[iz][ix][iy];
                  strcpy (path, input.filename[FN_PATH]);
                  strcat (path, "/ic/yang56/");
                  wavelength = output->wl.lambda_r[iv];

                  status =
                    yang56 (wavelength / 1000.0, lay_reff[0], caoth_habit, path, lay_dtau, lay_ssa, lay_g1, lay_g2, lay_ff, newkey);
                  CHKERR (status);

                  if (caoth3d->cldproperties == CLD_LWCREFF || caoth3d->cldproperties == CLD_LWCREFFCF)
                    caoth3d->ext[iz][ix][iy] = caoth3d->lwc[iz][ix][iy] * lay_dtau[0] / 1000.0;

                  caoth3d->ssa[iz][ix][iy] = lay_ssa[0];
                  caoth3d->g1[iz][ix][iy]  = lay_g1[0];
                  caoth3d->g2[iz][ix][iy]  = lay_g2[0];
                  caoth3d->ff[iz][ix][iy]  = lay_ff[0];

                  break;

                case PROP_ECHAM4:

                  if (strcasecmp (input_caoth.name, "wc") == 0) {
                    status =
                      wc_echam4 (output->wl.lambda_r[iv], nlev - 1, lay_lwc, lay_reff, lay_dtau, lay_g1, lay_ssa, lay_zd, lyr_flag);
                    CHKERR (status);
                  } else if (strcasecmp (input_caoth.name, "ic") == 0) {
                    status =
                      ic_echam4 (output->wl.lambda_r[iv], nlev, lay_lwc, lay_reff, lay_dtau, lay_g1, lay_ssa, lay_zd, lyr_flag);
                    CHKERR (status);
                  } else {
                    fprintf (stderr, "Error! You can not use echam4 for profile %s\n", input_caoth.name);
                    CHKERR (-1);
                  }

                  /* write results to final destination */

                  /* need to convert again from 1/km to 1/m */
                  /* because wcloud() thinks in km rather   */
                  /* than m                                 */
                  if (caoth3d->cldproperties == CLD_LWCREFF || caoth3d->cldproperties == CLD_LWCREFFCF)
                    caoth3d->ext[iz][ix][iy] = lay_dtau[0] / 1000.0;

                  caoth3d->ssa[iz][ix][iy] = lay_ssa[0];
                  caoth3d->g1[iz][ix][iy]  = lay_g1[0];
                  caoth3d->g2[iz][ix][iy]  = 0;
                  caoth3d->ff[iz][ix][iy]  = 1.0;
                  if (strcasecmp (input_caoth.name, "ic") == 0)
                    caoth3d->f[iz][ix][iy] = 0.0;

                  break;

                case PROP_MIE:
                case PROP_FILE:
                case PROP_IC_MIE:
                case PROP_BAUM:
                case PROP_BAUM_V36:
                case PROP_HEY:
                case PROP_YANG2013:

                  /* interpolate extinction and absorption coefficient */
                  /* of each box to given effective radius             */
                  i1 = (int) ((caoth3d->reff[iz][ix][iy] - caoth3d->phase.r0) / caoth3d->phase.dr);

                  /* careful with rounding errors at the array boundaries */
                  i1 = (i1 > caoth3d->phase.n - 1 ? caoth3d->phase.n - 1 : i1);
                  i1 = (i1 < 0 ? 0 : i1);

                  /* BM, changed interpolation: now interpolating scattering coefficient */
                  /* instead of single scattering albedo                                 */

                  if (i1 == caoth3d->phase.n - 1) {
                    extinc = caoth3d->phase.ext[i1];
                    gg     = caoth3d->phase.gg[i1];
                    scattr = caoth3d->phase.sca[i1];
                    dscale = caoth3d->phase.dscale[i1];
                    f      = caoth3d->phase.f[i1];
                  } else {
                    delta_r = caoth3d->reff[iz][ix][iy] - (caoth3d->phase.r0 + (double) i1 * caoth3d->phase.dr);

                    extinc = caoth3d->phase.ext[i1] + caoth3d->phase.dext[i1] * delta_r;
                    gg     = caoth3d->phase.gg[i1] + caoth3d->phase.dgg[i1] * delta_r;
                    scattr = caoth3d->phase.sca[i1] + caoth3d->phase.dsca[i1] * delta_r;

                    /* ??? DDDDD interpolate the delta scaling factor f   DDDDD ??? */
                    if (caoth3d->phase.f[i1] != 0 || caoth3d->phase.f[i1 + 1] != 0)
                      f = interpolate_f (&(caoth3d->phase), i1, i1 + 1, delta_r);
                    else
                      f = 0;

                    if (f < 0 || f > 1) {
                      fprintf (stderr, "Error in %s, line %d, function %s:\n", __FILE__, __LINE__, __func__);
                      fprintf (stderr, "delta-scaling factor %f out of range!\n", f);
                      fprintf (stderr,
                               "Tried to interpolate between radii %f and %f,\n",
                               caoth3d->phase.r0 + (double) i1 * caoth3d->phase.dr,
                               caoth3d->phase.r0 + (double) (i1 + 1) * caoth3d->phase.dr);
                      fprintf (stderr, "f values %f and %f\n", caoth3d->phase.f[i1], caoth3d->phase.f[i1 + 1]);
                      CHKERR (-1);
                    }

                    /* ??? DDDDD interpolate the delta scaling factor dscale DDDDD ??? */
                    /* dscale  = caoth3d->phase.dscale[i1] + caoth3d->phase.ddscale[i1] * delta_r; */
                    if (caoth3d->phase.dscale[i1] != 0 || caoth3d->phase.dscale[i1 + 1] != 0)
                      dscale = interpolate_dscale (&(caoth3d->phase), i1, i1 + 1, delta_r);
                    else
                      dscale = 0;

                    if (dscale < 0 || dscale > 1) {
                      fprintf (stderr, "Error in %s, line %d, function %s:\n", __FILE__, __LINE__, __func__);
                      fprintf (stderr, "delta-scaling factor %f out of range!\n", dscale);
                      fprintf (stderr,
                               "Tried to interpolate between radii %f and %f,\n",
                               caoth3d->phase.r0 + (double) i1 * caoth3d->phase.dr,
                               caoth3d->phase.r0 + (double) (i1 + 1) * caoth3d->phase.dr);
                      fprintf (stderr, "f values %f and %f\n", caoth3d->phase.dscale[i1], caoth3d->phase.dscale[i1 + 1]);
                      CHKERR (-1);
                    }
                  }

                  /* single scattering albedo */
                  omega = scattr / extinc;

                  /* extinction coefficient */
                  if (caoth3d->cldproperties == CLD_LWCREFF || caoth3d->cldproperties == CLD_LWCREFFCF)
                    caoth3d->ext[iz][ix][iy] = caoth3d->lwc[iz][ix][iy] * extinc / 1000.0;

                  /* single scattering albedo */
                  caoth3d->ssa[iz][ix][iy] = omega;

                  /* asymmetry parameter */
                  caoth3d->g1[iz][ix][iy] = gg;

                  /* external delta-scaling factor CE */
                  caoth3d->f[iz][ix][iy] = f;

                  /* internal delta-scaling factor TZ */
                  caoth3d->dscale[iz][ix][iy] = dscale;

                  break;

                case PROP_RAYTRACING:
                  /* BM, AW, 4.4.2024: For geometrical optics we need extinction efficiency Qext=2 to be consistent with */
                  /* all other cloud parameterizations; the fact that geometrical optics without diffraction actually    */
                  /* expects Qext=1 is taken care of in mystic.c:scattering() where only 50% of the photons are          */
                  /* actually scattered.                                                                                 */
                  caoth3d->ext[iz][ix][iy] = 3.0 * caoth3d->lwc[iz][ix][iy] / 2.0 / rho_ice / caoth3d->reff[iz][ix][iy];

                  /* single scattering albedo */
                  caoth3d->ssa[iz][ix][iy] = 1.0;

                  /* probably not used: */
                  caoth3d->f[iz][ix][iy]      = 0;
                  caoth3d->dscale[iz][ix][iy] = 0;
                  caoth3d->g1[iz][ix][iy]     = 0;
                  caoth3d->g2[iz][ix][iy]     = 0;
                  caoth3d->ff[iz][ix][iy]     = 0;

                  break;
                default:
                  fprintf (stderr,
                           "Error, unknown property %d for %s. This is not\n",
                           input_caoth.properties,
                           input_caoth.fullname);
                  fprintf (stderr, "supposed to happen and indicates a coding error.\n");
                  CHKERR (-1);
                  break;
                } /* end switch (properties) */

                /* CE: Apply delta scaling to extinction. Here we
                   assume that the specified input files correspond to
                   unscaled extinction */
                /* RPB: except for ic, f is normally 0, and this does
		   not do anything, this is good! */
                if (caoth3d->cldproperties != CLD_LWCREFF || caoth3d->cldproperties == CLD_LWCREFFCF) {
                  ssa_us = caoth3d->ssa[iz][ix][iy] / (1.0 + caoth3d->f[iz][ix][iy] * (caoth3d->ssa[iz][ix][iy] - 1.0));
                  caoth3d->ext[iz][ix][iy] *= (1.0 - ssa_us * caoth3d->f[iz][ix][iy]);
                }

              } /* end if ext or lwc > 0 */
            }   /* end for iy */
          }     /* end for ix */
        }       /* end if (caoth3d->threed[iz]) */
      }         /* end if (caoth3d->nthreed>0) */
    }           /* for (iz=0; iz<caoth3d->nlyr; iz++) */

    /* Switch of scattering or/and absorption if no_scattering/no_absorption are specified */
    if (caoth3d->no_scattering) {
      if (!input.quiet)
        fprintf (stderr, " ... Switch off scattering for 3D profile %s \n", caoth3d->name);

      for (iz = 0; iz < caoth3d->nlyr; iz++) {
        if (caoth3d->threed[iz]) {
          for (ix = 0; ix < caoth3d->Nx; ix++) {
            for (iy = 0; iy < caoth3d->Ny; iy++) {
              caoth3d->ext[iz][ix][iy] = (1. - caoth3d->ssa[iz][ix][iy]) * caoth3d->ext[iz][ix][iy];
              caoth3d->ssa[iz][ix][iy] = 0.0;
            }
          }
        }
      }
    }

    if (caoth3d->no_absorption) {
      if (!input.quiet)
        fprintf (stderr, " ... Switch off absorption for 3D profile %s \n", caoth3d->name);
      for (iz = 0; iz < caoth3d->nlyr; iz++) {
        if (caoth3d->threed[iz]) {
          for (ix = 0; ix < caoth3d->Nx; ix++) {
            for (iy = 0; iy < caoth3d->Ny; iy++) {
              caoth3d->ext[iz][ix][iy] = caoth3d->ssa[iz][ix][iy] * caoth3d->ext[iz][ix][iy];
              caoth3d->ssa[iz][ix][iy] = 1.0;
            }
          }
        }
      }
    }

    for (iz = 0; iz < output_caoth->nlev - 1; iz++) {
      /* for a 1D caoth we only need to interpolate the delta scaling factor;  */
      /* this is a bit unfortunate but we need to do that here because         */
      /* all other properties are already interpolated in setup_caoth()       */
      /* where dscale is not yet known                                         */

      if (output_caoth->microphys.lwc_layer != NULL) { /* should check on 1d */

        if (output_caoth->microphys.lwc_layer[iz] > 0) {

          switch (input_caoth.properties) {
          case PROP_HU:
          case PROP_FU:
          case PROP_ECHAM4:
          case PROP_KEY:
          case PROP_YANG:
          case PROP_RAYTRACING:

            /* no need to do anything because delta scaling only for */
            /* non-HG phase functions                                */

            break;

          case PROP_MIE:
          case PROP_FILE:
          case PROP_IC_MIE:
          case PROP_BAUM:
          case PROP_BAUM_V36:
          case PROP_HEY:
          case PROP_YANG2013:

            /* determine index of effective radius */
            i1 = (int) ((output_caoth->microphys.effr_layer[iz] - caoth3d->phase.r0) / caoth3d->phase.dr);

            /* careful with rounding errors at the array boundaries */
            i1 = (i1 > caoth3d->phase.n - 1 ? caoth3d->phase.n - 1 : i1);
            i1 = (i1 < 0 ? 0 : i1);

            /* ??? DDDDD interpolate delta-scaling factor dscale DDDDD ??? */
            if (caoth3d->phase.dscale[i1] != 0 || caoth3d->phase.dscale[i1 + 1] != 0) {
              delta_r = output_caoth->microphys.effr_layer[iz] - (caoth3d->phase.r0 + (double) i1 * caoth3d->phase.dr);
              dscale  = interpolate_dscale (&(caoth3d->phase), i1, i1 + 1, delta_r);
            } else
              dscale = 0;

            if (dscale < 0 || dscale > 1) {
              fprintf (stderr, "Error in %s, line %d, function %s:\n", __FILE__, __LINE__, __func__);
              fprintf (stderr, "delta-scaling factor %f out of range!\n", dscale);
              fprintf (stderr,
                       "Tried to interpolate between radii %f and %f,\n",
                       caoth3d->phase.r0 + (double) i1 * caoth3d->phase.dr,
                       caoth3d->phase.r0 + (double) (i1 + 1) * caoth3d->phase.dr);
              fprintf (stderr, "f values %f and %f\n", caoth3d->phase.dscale[i1], caoth3d->phase.dscale[i1 + 1]);
              CHKERR (-1);
            }

            output_caoth->optprop.dscale[iv][iz] = dscale;

            /*
	      output_caoth->optprop.dscale[iv][iz] = caoth3d->phase.dscale[i1] + caoth3d->phase.ddscale[i1] * 
	      (output_caoth->microphys.effr_layer[iz] - (caoth3d->phase.r0 + (double) i1 * caoth3d->phase.dr));
	    */

            break;
          default:
            fprintf (stderr, "Error, unknown property %d for %s. This is not\n", input_caoth.properties, input_caoth.fullname);
            fprintf (stderr, "supposed to happen and indicates a coding error.\n");
            CHKERR (-1);
          }
        } /* end if lwc_layer */
      }   /* end if exists lwc_layer */

    } /* for (iz=0; iz<caoth3d->nlyr; iz++) */
    break;

  default:
    fprintf (stderr, "Error, unknown cldproperty %d for %s\n", caoth3d->cldproperties, input_caoth.fullname);
    CHKERR (-1);
  }

  /* ulrike: for mystic, calculate tausol3d for caoth */
#if HAVE_MYSTIC
#if HAVE_TIPA
  if (input.rte.mc.tipa == TIPA_DIR) {
    fprintf (stderr, " ... calling mc_tipa_dir for %s\n", input_caoth.fullname);
    status = mc_tipa_dir (output_caoth, caoth3d, &(output_caoth->tipa), iv, input, &(output->wl));
    CHKERR (status);
    /*give adress of caoth3d, since tausol will be allocated/written
      in mc_tipa_dir otherwise only a copy would be generated */
  }
#endif
#endif

  status = ASCII_free_float (wvl_dtau, 18);
  CHKERR (status);

  status = ASCII_free_float (wvl_g1, 18);
  CHKERR (status);

  status = ASCII_free_float (wvl_ssa, 18);
  CHKERR (status);

  return 0;
}

/* Free 3D caoth optical properties */

int free_caoth_mystic (int properties, caoth3d_out_struct* caoth3d) {
  switch (caoth3d->cldproperties) {
  case CLD_OPTPROP:
    break;

  case CLD_EXTREFF:
  case CLD_LWCREFF:
  case CLD_LWCREFFCF:

    /* free phase function table */
    switch (properties) {
    case PROP_HU:
    case PROP_ECHAM4:
    case PROP_KEY:
    case PROP_YANG:
    case PROP_FU:
      break;

    case PROP_MIE:
    case PROP_FILE:
    case PROP_IC_MIE:
    case PROP_BAUM:
    case PROP_BAUM_V36:
    case PROP_HEY:
    case PROP_YANG2013:
      free_phase_function_table (&(caoth3d->phase));
      break;

    default:
      fprintf (stderr, "Error, unknown property %d for %s in free_caoth_mystic. This is not\n", properties, caoth3d->fullname);
      fprintf (stderr, " This is not supposed to happen and indicates a coding error.\n");
      CHKERR (-1);
      break;
    }
    break;
  default:
    fprintf (stderr, "Error, unknown cldproperty %d for %s in free_caoth_mystic.\n", caoth3d->cldproperties, caoth3d->fullname);
    fprintf (stderr, " This is not supposed to happen and indicates a coding error.\n");
    CHKERR (-1);
    break;
  }

  return 0;
}

/* check if two floats are equal - still not optmical, but how to improve? */

static inline int float_equal (float a, float b) {

  if (a == 0 || b == 0) {
    if (fabs (a - b) < MC_EPSILON)
      return 1;
  } else { /* relative difference smaller than MC_EPSILON */
    if (fabs (a - b) < MC_EPSILON * fabs (a))
      return 1;
  }

  return 0;
}

/* ??? DDDDD interpolate the delta scaling factor f     DDDDD ??? */
/* Interpolate the external delta scaling factor; the idea is to interpolate the      */
/* scaled and unscaled extinction and scattering coefficients and to calculate the    */
/* delta scaling factor from those - this seems more appropriate than a simple linear */
/* interpolation of the delta scaling factor. In the phase function table we          */
/* store the scaled optical properties (with respect to the external delta scaling    */
/* factor) and need to calculate the unscaled ones.                                   */

static float interpolate_f (phase_function_table* p, int i1, int i2, double delta_r) {
  double f1 = 0, f2 = 0, f = 0;
  double om1 = 0, om2 = 0, om1_unscaled = 0, om2_unscaled = 0, om_unscaled = 0;
  double ext1 = 0, ext2 = 0, ext = 0, ext1_unscaled = 0, ext2_unscaled = 0, ext_unscaled = 0;
  double sca_unscaled = 0;

  om1 = p->sca[i1] / p->ext[i1];
  om2 = p->sca[i2] / p->ext[i2];

  f1 = p->f[i1];
  f2 = p->f[i2];

  /* compare linear and improved interpolation for debugging purposes */
  /*  fprintf (stderr, " ... f linear %.6f ", f1 + delta_r / p->dr * (f2 - f1)); */

  om1_unscaled = om1 / (1.0 - f1 + f1 * om1);
  om2_unscaled = om2 / (1.0 - f2 + f2 * om2);

  ext1 = p->ext[i1];
  ext2 = p->ext[i2];

  ext1_unscaled = ext1 / (1.0 - f1 * om1_unscaled);
  ext2_unscaled = ext2 / (1.0 - f2 * om2_unscaled);

  /* interpolate scaled and unscaled extinction as well as unscaled single scattering albedo */
  ext          = ext1 + delta_r / p->dr * (ext2 - ext1);
  ext_unscaled = ext1_unscaled + delta_r / p->dr * (ext2_unscaled - ext1_unscaled);
  sca_unscaled = ext1_unscaled * om1_unscaled + delta_r / p->dr * (ext2_unscaled * om2_unscaled - ext1_unscaled * om2_unscaled);
  om_unscaled  = sca_unscaled / ext_unscaled;

  f = (1.0 - ext / ext_unscaled) / om_unscaled;
  /*  fprintf (stderr, "vs. improved %.6f interpolation\n", f); */

  return f;
}

/* ??? DDDDD interpolate the delta scaling factor dscale   DDDDD ??? */
/* Interpolate the internal delta scaling factor; the idea is to interpolate the      */
/* scaled and unscaled extinction and scattering coefficients and to calculate the    */
/* delta scaling factor from those - this seems more appropriate than a simple linear */
/* interpolation of the delta scaling factor. In the phase function table we          */
/* store the unscaled optical properties (with respect to the internal delta scaling  */
/* factor) and need to calculate the scaled ones.                                     */

static float interpolate_dscale (phase_function_table* p, int i1, int i2, double delta_r) {
  double ext1 = 0, ext2 = 0, ext = 0, ext1_scaled = 0, ext2_scaled = 0, ext_scaled = 0;
  double om1 = 0, om2 = 0, om = 0;
  double dscale1 = 0, dscale2 = 0;
  double sca    = 0;
  double dscale = 0;

  ext1 = p->ext[i1];
  ext2 = p->ext[i2];

  om1 = p->sca[i1] / p->ext[i1];
  om2 = p->sca[i2] / p->ext[i2];

  dscale1 = p->dscale[i1];
  dscale2 = p->dscale[i2];

  /* compare linear and improved interpolation for debugging purposes */
  /* fprintf (stderr, " ... dscale linear %.6f ", dscale1 + delta_r / p->dr * (dscale2 - dscale1)); */

  ext1_scaled = (1.0 - dscale1 * om1) * ext1;
  ext2_scaled = (1.0 - dscale2 * om2) * ext2;

  ext        = ext1 + delta_r / p->dr * (ext2 - ext1);
  ext_scaled = ext1_scaled + delta_r / p->dr * (ext2_scaled - ext1_scaled);
  sca        = ext1 * om1 + delta_r / p->dr * (ext2 * om2 - ext1 * om2);
  om         = sca / ext;

  dscale = (1.0 - ext_scaled / ext) / om;
  /* fprintf (stderr, "vs. improved %.6f interpolation\n", dscale); */

  return dscale;
}

/* function used by analytic_area() */
double analytic_area_stammfunk (double a, double b, double c, double x1, double x2, double y) {
  double T3 = b + c * x1;
  double T4 = b + c * x2;
  double T5 = a + c * y;
  double T1 = sqrt (1.0 + T3 * T3 + T5 * T5);
  double T2 = sqrt (1.0 + T4 * T4 + T5 * T5);

  return (2.0 * T5 * (T4 * T2 - T3 * T1) + 2.0 * atan ((T3 * T5) / T1) - 2.0 * atan ((T4 * T5) / T2) +
          T5 * (3.0 + T5 * T5) * (log (b + c * x2 + T2) - log (b + c * x1 + T1)) + T4 * (T4 * T4 + 3.0) * log (a + c * y + T2) -
          T3 * (T3 * T3 + 3.0) * log (a + c * y + T1)) /
         (6.0 * c * c);
}

/***************************************************************************************************/
/* analytic_area calculates the surface area of an elevation pixel between limits x1,y1 and x2,y2; */
/* The integration was done by Petra Hausmann in March 2009; http://www.integrals.com provides a   */
/* solution in principle, but the formula covers several pages; the result was simplified and      */
/* rearranged to obtain a numerically stable solution as the original solution proved to be        */
/* highly unstable for small c's. For very small c the solution might still diverge but the        */
/* threshold set now (fabs(c)*(x2-x1)*(y2-y1)<MC_EPSILON) seems to be appropriate and safe;        */
/* the accuracy of the formula can still be tested by comparison with the numerical solution,      */
/* to be activated by uncommenting AREA_DEBUG_N in which case the numerical result is written to   */
/* area2D.dat in addition to the analytical one.                                                   */
/***************************************************************************************************/

double analytic_area (double a, double b, double c, double x1, double x2, double y1, double y2) {
  /* in the limit of c->0 we get 0/0 which is numerically unstable; */
  /* therefore we choose an upper limit which is somehow related    */
  /* to the area size; this is not fully tested!                    */
  if (fabs (c) * (x2 - x1) * (y2 - y1) < MC_EPSILON)
    return sqrt (1.0 + a * a + b * b) * (x2 - x1) * (y2 - y1);
  else
    return analytic_area_stammfunk (a, b, c, x1, x2, y2) - analytic_area_stammfunk (a, b, c, x1, x2, y1);
}

/*************************************************************************************/
/* integrate surface area in a 2D elevation grid between limits (x1,y1) and (x2,y2); */
/* uses the analytical formula provided by Petra Hausmann in March 2009; see         */
/* description of analytic_area().                                                   */
/*************************************************************************************/

#if HAVE_MYSTIC
double calc_elevation_area (elevation_struct* elev, double x1, double y1, double x2, double y2) {
  photon_struct ptmp;
  int           ie = 0, je = 0, ie1 = 0, je1 = 0, ie2 = 0, je2 = 0;
  double        area  = 0;
  double        xleft = 0, xright = 0, yleft = 0, yright = 0;
  double        xe1 = 0, ye1 = 0; /*xe2=0, ye2=0*/

  /* first determine which elevation pixels we need */
  /* elevation coordinates */
  ptmp.x[0] = x1;
  ptmp.x[1] = y1;
  elev_coord (&ptmp, elev, &ie1, &je1);

  /* elevation coordinates */
  ptmp.x[0] = x2;
  ptmp.x[1] = y2;
  elev_coord (&ptmp, elev, &ie2, &je2);

  /* loop over all elevation pixels in the sample pixel under consideration;     */
  /* we also include some pixels with area 0 but that doesn't cause any problems */
  /* except a very small computational overhead.                                 */
  for (ie = ie1; ie <= ie2; ie++)
    for (je = je1; je <= je2; je++) {
      xe1 = (double) ie * elev->delX;
      /* commented out by RB, because was not used but caused compiler warnings */
      /* xe2 = xe1 + elev->delX; */

      ye1 = (double) je * elev->delY;
      /* commented out by RB, because was not used but caused compiler warnings */
      /* ye2 = ye1 + elev->delY; */

      if (ie == ie1)
        xleft = x1 - xe1;
      else
        xleft = 0.0;

      if (ie == ie2)
        xright = x2 - xe1;
      else
        xright = elev->delX;

      if (je == je1)
        yleft = y1 - ye1;
      else
        yleft = 0.0;

      if (je == je2)
        yright = y2 - ye1;
      else
        yright = elev->delY;

      area += analytic_area (elev->surf[ie][je].a, elev->surf[ie][je].b, elev->surf[ie][je].c, xleft, xright, yleft, yright);
    }

  area /= ((x2 - x1) * (y2 - y1));

  return area;
}
#endif
