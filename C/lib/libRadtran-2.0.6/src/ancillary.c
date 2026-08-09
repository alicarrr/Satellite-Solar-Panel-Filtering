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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocnd.h"
#include "ascii.h"
#include "cdisort.h"
#include "ckdfu.h"
#include "common_math.h"
#include "errors.h"
#include "numeric.h"
#include "solver.h"
#include "uvspec.h"
#include "options_stash.h"

#if HAVE_NETCDF4
#include <netcdf.h>
#include "netCDF_functions.h"
#endif

#define EPSILON 1E-6
#define s2day 3600.0 * 24.0 /* seconds to day */
#define mW2W 1.E-3;         /* mW to W */

#define SWAP(a, b, t)                                                                                                              \
  {                                                                                                                                \
    t = a;                                                                                                                         \
    a = b;                                                                                                                         \
    b = t;                                                                                                                         \
  }

/* prototypes of internal functions */
static int cnvlv (double* x_spec, float* y_spec, int n_spec, double* x_slit, double* y_slit, int n_slit, int std);

static float radiance2bt (float rad, float* wvnmlo, float* wvnmhi, float* filter, int n, int processing, int ivi);

static int output2bt (input_struct input, output_struct* output, int iv, int is_3d);

static int read_photon_file (char* filename, float* lambda_r, int nlamdba_r, float** fraction);

static int set_raman_wl_grid (wl_inp_struct* wl_inp, wl_out_struct* wl_out, char* filename, int quiet);
static int find_raman_closest_wl_in_solar_file (float* wls, float* wle, char* filename, int quiet);
static int set_transmittance_wl_grid (wl_inp_struct* wl_inp, wl_out_struct* wl_out, int quiet);
static int set_transmittance_wl_grid_lowtran (wl_inp_struct* wl_inp, wl_out_struct* wl_out, int quiet);

static int
set_transmittance_wl_grid_reptran (input_struct input, float** lambda_lower, float** lambda_upper, wl_out_struct* wl_out);

static int set_rte_wl_grid_reptran (input_struct input, output_struct* output);

static int write_spectral3D (input_struct input, output_struct* output);
static int write_spectral3D_nc (input_struct input, output_struct* output);

static int        spec2rgb (input_struct input, output_struct* output);
static int        spec2rgb3D (input_struct input, output_struct* output);
static int        raman_spec2spec (input_struct input, output_struct* output);
static inline int float_equal (float a, float b);
double            linpol (double x1, double x2, double y1, double y2, double x);

double nm_to_inv_cm (double wavelength_nm);
double inv_cm_to_nm (double wavenumber_inv_cm);
double polarizability_anisotropy_N2 (double nu);
double polarizability_anisotropy_O2 (double nu);
int    pfraction_reptran (wl_out_struct* wl_out);

// determine the model layer corresponding to a given zout layer
// for a given sample height, e.g. from output->atm.zout_sur,
// returns the index k and k-1 of the atm.zd that is to be considered
// also return the corresponding heights above and below that layer
// output parameters are not set if NULL
static int mc_backward_absorption_find_vert_index (const atm_out_struct* atm,
                                                   const float           sample_z,
                                                   int*                  kbot,
                                                   int*                  ktop,
                                                   float*                zbot,
                                                   float*                ztop) {
  int lc;
  for (lc = 0; lc < atm->nlev; lc++) {
    if (approx_rel (atm->zd[lc], sample_z, 10 * FLT_EPSILON))
      break;
  }
  lc = MIN (atm->nlev - 1, lc);
  if (kbot)
    *kbot = lc;
  if (zbot)
    *zbot = atm->zd[lc];
  if (lc == 0) {
    CHKERROUT (-1, "found sample height to be at TOA, beware that there is no layer defined above");
  }
  if (ktop)
    *ktop = lc - 1;
  if (ztop)
    *ztop = atm->zd[lc - 1];
  return 0;
}

/******************************************************************************/
/* select the required wavelength range from an array of wavelengths (lambda) */
/******************************************************************************/

static int select_wavelength_indices (const float* lambda,
                                      const int    nlambda,
                                      int          start_index,
                                      int          end_index,
                                      int          quiet,
                                      int          raman,
                                      int*         lower,
                                      int*         upper,
                                      float*       lambda_lower,
                                      float*       lambda_upper) {
  int iv = 0;
  if (start_index > 0 && end_index > 0) {
    /* if wavelength indices are defined */

    /* check if start_index and end_index sorted; if not, sort */
    if (start_index > end_index) {
      int tmp;
      SWAP (start_index, end_index, tmp);
    }

    /* check if we are out of range */
    if (end_index > nlambda) {
      fprintf (stderr, "end index %d > nlambda %d\n", end_index, nlambda);
      CHKERR (-1);
    }

    /* the selected range is valid */
    *lower = start_index - 1;
    *upper = end_index - 1;

    *lambda_lower = lambda[*lower];
    *lambda_upper = lambda[*upper];

  } else {
    /* if wavelength is specified search for suitable indices */

    while (lambda[iv] <= *lambda_lower)
      if (++iv == nlambda)
        break;

    if (iv > 0)
      iv -= 1;

    *lower = iv;

    iv = 0;
    while (lambda[iv] < *lambda_upper) {
      if (++iv == nlambda)
        break;
    }

    /* Need to include larger wavelength range if Raman scattering is included */
    if (raman) {
      iv = 0;

      while (lambda[iv] <= *lambda_lower)
        if (++iv == nlambda)
          break;
      if (iv > 0)
        iv -= 1;
      *lower = iv;

      iv = 0;
      while (lambda[iv] < *lambda_upper) {
        if (++iv == nlambda)
          break;
      }

      if (iv == nlambda)
        iv -= 1;

      *upper = iv;
    }

    if (iv == nlambda)
      iv -= 1;

    *upper = iv;
  }

  return 0;
}

/********************************************/
/* Setup the transmittance wavelength grid. */
/********************************************/
int setup_wlgrid (input_struct input, output_struct* output) {
  int    iv = 0, status = 0, monochromatic = 0;
  float *lambda_lower = NULL, *lambda_upper = NULL;
  float  test_lower = 0, test_upper = 0;

  if (!input.quiet) {
    fprintf (stderr,
             " ... setup_wlgrid:"
             " inp.wl.start %f inp.wl.end %f"
             " inp.wl.start_idx %d inp.wl.end_idx %d\n",
             input.wl.start,
             input.wl.end,
             input.wl.start_index,
             input.wl.end_index);
  }

  /* Check whether representative wavelengths will be used. */
  if (input.ck_scheme == CK_REPTRAN || input.ck_scheme == CK_REPTRAN_CHANNEL)
    output->wl.use_reptran = 1;
  else
    output->wl.use_reptran = 0;

  /* in case of RGB conversion we need an internal wavelength grid */
  if (input.processing == PROCESS_RGB || input.processing == PROCESS_RGBRAW) {
    /* reset wavelength range */
    input.wl.start = -999;
    input.wl.end   = -999;

    switch (input.source) {
    case SRC_SOLAR:
    case SRC_BLITZ: /* BCA */
    case SRC_LIDAR: /* BCA */
      /* and replace transmittance grid */
      strcpy (input.filename[FN_WLTRANS], input.filename[FN_PATH]);
      strcat (input.filename[FN_WLTRANS], "solar_flux/rgb");

      if (!input.quiet) {
        fprintf (stderr, " ... ignoring user-defined wavelength range and using\n");
        fprintf (stderr, " ...   internal wavelength grid from %s\n", input.filename[FN_WLTRANS]);
      }

      break;

    case SRC_THERMAL:
      /* and replace thermal bands file */
      strcpy (input.filename[FN_WLBANDS], input.filename[FN_PATH]);
      strcat (input.filename[FN_WLBANDS], "solar_flux/rgb_bands");

      if (!input.quiet) {
        fprintf (stderr, " ... ignoring user-defined wavelength range and using\n");
        fprintf (stderr, " ...   internal thermal bands from %s\n", input.filename[FN_WLBANDS]);
      }

      break;

    default:
      fprintf (stderr, "Error, unsupported source %d\n", input.source);
      CHKERR (-1);
    }
  }

  /* set default unit for band width */
  output->bandwidth_unit = input.bandwidth_unit;

  if (output->bandwidth_unit == UNIT_NOT_DEFINED) {
    switch (input.source) {
    case SRC_SOLAR:
    case SRC_BLITZ: /* BCA */
    case SRC_LIDAR: /* BCA */
      output->bandwidth_unit = UNIT_PER_NM;
      break;
    case SRC_THERMAL:
      output->bandwidth_unit = UNIT_PER_CM_1;
      break;
    default:
      fprintf (stderr, "Error, unsupported source %d\n", input.source);
      CHKERR (-1);
    }
  }

  if (input.raman) {
    /* Internally the lower and upper wavelengths for Raman scattering are different from */
    /* what the user specifies. This is setup in setup_wlgrid (ancillary.c).              */
    /* For Raman scattering only include wavelengths that the user asked.                 */
    /* Internally we have to include more wavelengths to account for                      */
    /* Raman scattered radiation                                                          */
    /* The values 196.8269 and 194.3015 are the maximum shifts (in cm-1) for N2 and O2    */
    /* as given in report: Accounting for Raman Scattering in DOAS, J.F. de Haan,         */
    /* SN-OMIE-KNMI-409, Version 1.0,  May 18, 2003. See also functions: crs_raman_N2 and */
    /* crs_raman_O2.                                                                      */

    output->wl.delta_wvl_raman_extra = 1.0;
    output->wl.delta_wvl_raman_lower = input.wl.start - inv_cm_to_nm (nm_to_inv_cm (input.wl.start) + 196.8269);
    input.wl.start                   = inv_cm_to_nm (nm_to_inv_cm (input.wl.start) + 196.8269) - output->wl.delta_wvl_raman_extra;
    output->wl.delta_wvl_raman_upper = inv_cm_to_nm (nm_to_inv_cm (input.wl.end) - 194.3015) - input.wl.end;
    input.wl.end                     = inv_cm_to_nm (nm_to_inv_cm (input.wl.end) - 194.3015) + output->wl.delta_wvl_raman_extra;
    status = find_raman_closest_wl_in_solar_file (&input.wl.start, &input.wl.end, input.filename[FN_EXTRATERRESTRIAL], input.quiet);
    CHKERR (status);
  }

  if (input.wl.start > 0 && input.wl.end > 0) {
    output->wl.start = input.wl.start;
    output->wl.end   = input.wl.end;
  }

  output->wl.type = WLGRID_NONE;

  if (strlen (input.filename[FN_FILTERFUNCTION]) > 0 && input.ck_scheme == CK_REPTRAN_CHANNEL)
    fprintf (stderr, "Error: Combining options 'mol_abs_param reptran_channel' and 'filter_function_file' is not allowed.");

  if (output->wl.use_reptran || input.ck_scheme == CK_CRS || input.ck_scheme == CK_LOWTRAN ||
      input.ck_scheme == CK_RAMAN) { /* no real correlated-k */

    if (strlen (input.filename[FN_WLBANDS]) > 0 && input.source == SRC_THERMAL) { /* thermal_bands_file */

      output->wl.type              = WLGRID_BANDS;
      output->wl.ignore_solar_file = 1; /* ignore solar file if FN_WLBANDS is used */

      if (!input.quiet)
        fprintf (stderr, " ... reading thermal_bands_file from %s\n", input.filename[FN_WLBANDS]);

      /* read center wavelength and band limits [wavenumbers] from file */
      status = read_3c_file_float (input.filename[FN_WLBANDS],
                                   &(output->wl.lambda_t),
                                   &(lambda_lower),
                                   &(lambda_upper),
                                   &output->wl.nlambda_t);

      if (status != 0) {
        fprintf (stderr, "Error %d opening %s\n", status, input.filename[FN_WLBANDS]);
        CHKERR (status);
      }

    }

    else if (strlen (input.filename[FN_WLTRANS]) > 0) { /* transmittance_wl_file */

      output->wl.type = WLGRID_USER;

      /* read internal wavelength grid from transmittance file */
      status = read_1c_file_float (input.filename[FN_WLTRANS], &(output->wl.lambda_t), &output->wl.nlambda_t);

      if (status != 0) {
        fprintf (stderr, "Error %d opening %s\n", status, input.filename[FN_WLTRANS]);
        CHKERR (status);
      }

      if (output->wl.use_reptran) {

        status = set_transmittance_wl_grid_reptran (input, &lambda_lower, &lambda_upper, &output->wl);
        CHKERR (status);
      }

    }

    else if (strlen (input.filename[FN_MOL_TAU_ABS]) > 0) { /* moltau_file */

      output->wl.type = WLGRID_MOLABS;

      if (!input.quiet) {
        fprintf (stderr, " ... molecular_tau_file specified but computational wavelength grid\n");
        fprintf (stderr, "     not explicitely defined; reading the wavelength grid\n");
        fprintf (stderr, "     from molecular_tau_file %s\n", input.filename[FN_MOL_TAU_ABS]);
      }
      status = read_molecular_absorption_lambda (input.filename[FN_MOL_TAU_ABS],
                                                 input.quiet,
                                                 &output->wl.lambda_t,
                                                 &output->wl.nlambda_t,
                                                 &monochromatic);
      if (status != 0) {
        fprintf (stderr, "Error %d reading wavelength grid from %s\n", status, input.filename[FN_MOL_TAU_ABS]);
        CHKERR (-1);
      }
    }

    if (output->wl.type == WLGRID_MOLABS && monochromatic == 1) {
      output->wl.lambda_t[0] = input.wl.start;
    } else if (output->wl.type == WLGRID_NONE ||
               (output->wl.type == WLGRID_MOLABS && monochromatic == 0 && input.ck_scheme == CK_RAMAN)) {

      output->wl.type = WLGRID_UVSPEC;

      if (output->wl.use_reptran) {

        /* transmittance wavelength grid consists of band centers */
        status = set_transmittance_wl_grid_reptran (input, &lambda_lower, &lambda_upper, &output->wl);
        CHKERR (status);

      } else if (input.ck_scheme == CK_CRS) {

        /* set up a reasonable wavelength grid for the radiative transfer calculation */
        status = set_transmittance_wl_grid (&input.wl, &output->wl, input.quiet);
        CHKERR (status);

      } else if (input.ck_scheme == CK_RAMAN) {

        /* Set the internal radiative transfer grid equal to the grid specified in */
        /* the extraterrestrial spectrum, because we use the absolute value of     */
        /* the solar source in all calculations.                                   */
        status = set_raman_wl_grid (&input.wl, &output->wl, input.filename[FN_EXTRATERRESTRIAL], input.quiet);
        CHKERR (status);

        /* For Raman scattering only include wavelengths that the user asked. */
        /* Internally we have to include more wavelengths because we need     */
        /* these cross sections to account for Raman scattered radiation.     */
        test_lower = output->wl.lambda_t[0] + output->wl.delta_wvl_raman_lower + output->wl.delta_wvl_raman_extra;
        test_upper =
          output->wl.lambda_t[output->wl.nlambda_t - 1] - output->wl.delta_wvl_raman_upper - output->wl.delta_wvl_raman_extra;
        for (iv = 0; iv < output->wl.nlambda_t; iv++) {
          if (output->wl.lambda_t[iv] < test_lower)
            output->wl.raman_start_id = iv + 1;
          if (output->wl.lambda_t[output->wl.nlambda_t - 1 - iv] > test_upper)
            output->wl.raman_end_id = output->wl.nlambda_t - iv - 2;
        }

      } else {
        status = set_transmittance_wl_grid_lowtran (&input.wl, &output->wl, input.quiet);
        CHKERR (status);
      }
    }

    /* here we need to calculate the full wavelength range */
    output->wl.iv_lower = 0;
    output->wl.iv_upper = output->wl.nlambda_t - 1;

  } else { /* correlated-k */

    output->wl.type = WLGRID_CK;

    /* read information about wavelength grid and quadrature points */
    switch (input.ck_scheme) {
    case CK_KATO:
    case CK_KATO2:
    case CK_KATO2_96:
    case CK_KATO2ANDWANDJI:

      /* read Kato et al. [1999] tables */
      status = kato_readtables (input.ck_scheme, &(output->ck), input.filename[FN_PATH], input.rte.mc.filename[FN_MC_PHOTONS]);
      CHKERR (status);
      break;

    case CK_FU:

      /* read Fu and Liou [1992/93] tables */
      status = fu_readtables (&(output->ck), input.filename[FN_PATH], input.rte.mc.filename[FN_MC_PHOTONS]);
      CHKERR (status);
      break;

    case CK_AVHRR_KRATZ:

      /* read Kratz [1999] tables */
      status = avhrr_kratz_readtables (&(output->ck), input.filename[FN_PATH], input.rte.mc.filename[FN_MC_PHOTONS]);
      CHKERR (status);
      break;

    case CK_FILE:

      /* read generic tables in CDF format */
      status = ck_generic_readtables (&(output->ck), input.ck_scheme_filename, input.rte.mc.filename[FN_MC_PHOTONS]);
      CHKERR (status);
      break;

    default:
      CHKERROUT (-1, "unsupported correlated-k scheme");
    }

    /* copy center wavelengths to transmittance grid  */
    /* and set wavenumber intervals                   */

    output->wl.nlambda_t = output->ck.n_wvl;
    output->wl.lambda_t  = calloc (output->wl.nlambda_t, sizeof (float));

    for (iv = 0; iv < output->wl.nlambda_t; iv++) {
      output->wl.lambda_t[iv] = output->ck.wvlc[iv + 1];
    }

  } /* end correlated-k */

  /* setup array containing the band limits */
  output->wl.wvnmlo_t = calloc (output->wl.nlambda_t, sizeof (float));
  output->wl.wvnmhi_t = calloc (output->wl.nlambda_t, sizeof (float));

  if (output->wl.type == WLGRID_CK) {
    for (iv = 0; iv < output->wl.nlambda_t; iv++) {
      output->wl.wvnmlo_t[iv] = output->ck.wvnlo[iv + 1];
      output->wl.wvnmhi_t[iv] = output->ck.wvnhi[iv + 1];
    }
  } else if (output->wl.type == WLGRID_BANDS && input.source == SRC_THERMAL) {
    /* band limits from thermal_bands_file */
    for (iv = 0; iv < output->wl.nlambda_t; iv++) {
      output->wl.wvnmlo_t[iv] = 1.0E7 / lambda_upper[iv];
      output->wl.wvnmhi_t[iv] = 1.0E7 / lambda_lower[iv];
    }
  } else {
    /* the default bandwidth is 1cm-1 to get the emittance per cm-1; */
    /* input.bandwidth can be set with thermal_bandwidth             */
    if (output->bandwidth_unit == UNIT_PER_CM_1) {
      for (iv = 0; iv < output->wl.nlambda_t; iv++) {
        output->wl.wvnmlo_t[iv] = 1.0E7 / output->wl.lambda_t[iv] - input.bandwidth / 2.0;
        output->wl.wvnmhi_t[iv] = output->wl.wvnmlo_t[iv] + input.bandwidth;
      }
    } else if (output->bandwidth_unit == UNIT_PER_NM) {
      for (iv = 0; iv < output->wl.nlambda_t; iv++) {
        output->wl.wvnmlo_t[iv] = 1.0E7 / (output->wl.lambda_t[iv] + input.bandwidth / 2.0);
        output->wl.wvnmhi_t[iv] = 1.0E7 / (output->wl.lambda_t[iv] - input.bandwidth / 2.0);
      }
    } else {
      fprintf (stderr, "Error, unsupported bandwidth_unit %d\n", output->bandwidth_unit);
      CHKERR (-1);
    }
  }

  /* free temporary wavelengths arrays */
  if (lambda_lower != NULL)
    free (lambda_lower);
  if (lambda_upper != NULL)
    free (lambda_upper);

  if (output->wl.type == WLGRID_CK || output->wl.type == WLGRID_BANDS || output->wl.type == WLGRID_USER ||
      (output->wl.type == WLGRID_MOLABS && monochromatic == 0)) {

    /* if start and end wavelength have not been set, use entire wavelength range */
    if (input.wl.start < 0 || input.wl.end < 0) {
      output->wl.start = output->wl.lambda_t[0];
      output->wl.end   = output->wl.lambda_t[output->wl.nlambda_t - 1];
      if (!input.quiet)
        fprintf (stderr, " ... setting wavelength range to %f - %fnm\n", output->wl.start, output->wl.end);
    }

    /* select wavelength range */
    status = select_wavelength_indices (output->wl.lambda_t,
                                        output->wl.nlambda_t,
                                        input.wl.start_index,
                                        input.wl.end_index,
                                        input.quiet,
                                        input.raman,
                                        &(output->wl.iv_lower),
                                        &(output->wl.iv_upper),
                                        &(output->wl.start),
                                        &(output->wl.end));
    CHKERR (status);

    /* check if the end wavelength is larger than 850nm in which case */
    /* we require user-selected molecular absorption properties       */
    /*
    if (output->wl.lambda_t[output->wl.iv_upper]>850 &&
	  (input.ck_scheme==CK_CRS && strlen(input.filename[FN_MOL_TAU_ABS])==0))  {
      fprintf (stderr, "Error, you want to do a spectral calculation for wavelengths larger than 850 nm. While uvspec\n");
      fprintf (stderr, "    treats ozone absorption correctly, molecular absorption is NOT considered in monochromatic\n");
      fprintf (stderr, "    uvspec calculations, as absorption cross-section are highly variable with wavelength.\n");
      fprintf (stderr, "    To consider molecular absorption other than ozone you have two choices \n");
      fprintf (stderr, "    with uvspec:\n");
      fprintf (stderr, "     (1) Do a line-by-line calculation using 'molecular_tau_file' to specify\n");
      fprintf (stderr, "         the wavelength-dependent absorption profile; to calculate the\n");
      fprintf (stderr, "         latter, you need something like David Edwards' genln2.\n");
      fprintf (stderr, "         ATTENTION: line-by-line calculations are very time-consuming!\n");
      fprintf (stderr, "     (2) Use the correlated-k approximation which is the most accurate\n");
      fprintf (stderr, "         solution after the line-by-line calculation; use either the\n");
      fprintf (stderr, "         pre-defined parameterization that come with libRadtran or provide\n");
      fprintf (stderr, "         your own; both options are selected with 'mol_abs_param ...'\n");
      fprintf (stderr, "\n");
      return -1;
    }
    */
  }

  /* now setup the ck structure for the LOWTRAN/SBDART table */
  if (input.ck_scheme == CK_LOWTRAN) {

    /* read LOWTRAN/SBDART tables */
    status = sbdart_readtables (&(output->ck),
                                output->wl.nlambda_t,
                                input.filename[FN_PATH],
                                input.rte.mc.filename[FN_MC_PHOTONS],
                                input.quiet);
    CHKERR (status);
  }

  if (output->wl.type != WLGRID_CK && output->wl.use_reptran == 0) {
    /* check if the internal grid covers the required wavelength range */
    if (output->wl.lambda_t[0] > output->wl.start || output->wl.lambda_t[output->wl.nlambda_t - 1] < output->wl.end) {
      fprintf (stderr,
               "Error, internal wavelength grid (%f - %f nm)\n",
               output->wl.lambda_t[0],
               output->wl.lambda_t[output->wl.nlambda_t - 1]);
      fprintf (stderr, "does not cover the user-defined wavelength range (%f - %f nm)\n", output->wl.start, output->wl.end);
      CHKERR (-1);
    }
  }

  /* inform the user about wavelength selection */
  if (!input.quiet) {
    if (output->wl.iv_upper != output->wl.nlambda_t - 1 || output->wl.iv_lower != 0) {
      if (!approx_rel (input.wl.start, NOT_DEFINED_FLOAT, 0)) {
        fprintf (stderr, "     user wavelength range: %f - %f nm\n", input.wl.start, input.wl.end);
      }
      fprintf (stderr, "     number of all wavelengths %d\n", output->wl.nlambda_t);
      fprintf (stderr, "     selected wavelength indices %d - %d\n", output->wl.iv_lower, output->wl.iv_upper);
      fprintf (stderr,
               "     wavelength bands boundaries: %f - %f nm\n",
               output->wl.lambda_t[output->wl.iv_lower],
               output->wl.lambda_t[output->wl.iv_upper]);
    }
  }

  return 0;
}

/**********************************************************************/
/* Setup the wavelength grid for the radiative transfer calculations. */
/**********************************************************************/
int setup_rte_wlgrid (input_struct input, output_struct* output) {

  int i, iv, status;

  if (output->wl.use_reptran) {

    /* checking compatibility of representative wavelengths with selected wavelength grid */
    if (output->wl.type == WLGRID_BANDS)
      CHKERROUT (-1, "combining representative wavelengths and thermal_bands_file is not allowed");

    if (output->wl.type == WLGRID_MOLABS)
      CHKERROUT (-1, "combining representative wavelengths and molecular_tau_file is not allowed");

    if (output->wl.type == WLGRID_UVSPEC || output->wl.type == WLGRID_USER) {

      status = set_rte_wl_grid_reptran (input, output);
      CHKERR (status);

      output->wl.iv_lower = 0;
      output->wl.iv_upper = output->wl.nlambda_r - 1;

      if (input.verbose) {
        fprintf (stderr, "      transmittance wavelength | radiative transfer wavelength | weight\n");
        for (iv = 0; iv < output->wl.nlambda_t; iv++)
          for (i = 0; i < output->wl.nlambda_in_reptran_band[output->wl.reptran_band_t[iv]]; i++)
            fprintf (stderr,
                     "               %12.6f nm |               %12.6f nm | %f  \n",
                     output->wl.lambda_t[iv],
                     output->wl.lambda_r[output->wl.reptran_band[output->wl.reptran_band_t[iv]][i]],
                     output->wl.weight_reptran_band[output->wl.reptran_band_t[iv]][i]);
      }

    } else {
      CHKERROUT (-1, "incompatible wavelength grid type");
    }

  } else {

    /* If representative wavelengths are not used, the wavelength grids for */
    /* radiative transfer and transmission are identical                    */
    output->wl.nlambda_r = output->wl.nlambda_t;
    output->wl.lambda_r  = output->wl.lambda_t;
    output->wl.wvnmlo_r  = output->wl.wvnmlo_t;
    output->wl.wvnmhi_r  = output->wl.wvnmhi_t;
  }

  /* need to read the photons file */
  if (output->wl.use_reptran || input.ck_scheme == CK_CRS || input.ck_scheme == CK_RAMAN) {

    if (strlen (input.rte.mc.filename[FN_MC_PHOTONS]) > 0) {
      status =
        read_photon_file (input.rte.mc.filename[FN_MC_PHOTONS], output->wl.lambda_r, output->wl.nlambda_r, &(output->wl.pfraction));
      if (status != 0) {
        fprintf (stderr, "Error %d reading photon file name %s\n", status, input.rte.mc.filename[FN_MC_PHOTONS]);
        CHKERR (status);
      }
    } else {
      output->wl.pfraction = calloc (output->wl.nlambda_r, sizeof (float));
      if (output->wl.use_reptran) {
        status = pfraction_reptran (&(output->wl));
        CHKERR (status);
      } else {
        for (iv = output->wl.iv_lower; iv <= output->wl.iv_upper; iv++) {
          output->wl.pfraction[iv] = 1.0 / (float) (output->wl.iv_upper - output->wl.iv_lower + 1);
        }
      }
    }
  }

  if (input.rte.mc.spectral_is) {
    if (!input.quiet)
      fprintf (stderr, " ... using %d wavelengths for spectral importance sampling \n", output->wl.nlambda_r);
    output->mc.alis.nlambda_abs = output->wl.nlambda_r;
    output->mc.alis.lambda      = calloc (output->mc.alis.nlambda_abs, sizeof (float));
    for (iv = 0; iv < output->mc.alis.nlambda_abs; iv++) {
      output->mc.alis.lambda[iv] = output->wl.lambda_r[iv];
    }

    /* Initialization for number of concentrations */
    output->mc.alis.Nc = 1;
  }

  if (input.rte.mc.concentration_is)
    output->mc.alis.nlambda_abs = 1;

  return 0;
}

double linpol (double x1, double x2, double y1, double y2, double x) {
  /* Linearly interpolate between two points to get wanted y value for x. */
  double y;
  double a = 0, b = 0;

  if (x2 > 0 && x1 > 0 && fabs (x2 - x1) < 0.0000001) {
    y = y1;
  } else {
    a = (y2 - y1) / (x2 - x1);
    b = y1 - a * x1;
    y = a * x + b;
  }
  return y;
}

/***********************************************************************/
/* Interpolate old_y (on the old_x grid) to new_y (on the new_x grid), */
/* either with                                                         */
/*     natural cubic splines (linear=0),                               */
/*     linear (linear=1), or                                           */
/*     logarithmic (linear=2)                                          */
/* interpolation methods.                                              */
/* The new_y must already be allocated with the right size (n_new_x)!! */
/* If the input is sorted in descending order,                         */
/*                   'descend' needs to be set to 1.                   */
/***********************************************************************/

int arb_wvn (const int    n_old_x,
             const float* old_x,
             const float* old_y,
             const int    n_new_x,
             const float* new_x,
             float*       new_y,
             const int    linear,
             const int    descend) {
  int     i = 0, j = 0, status = 0;
  double *a0 = NULL, *a1 = NULL, *a2 = NULL, *a3 = NULL;
  double *x = NULL, *y = NULL;
  double* tmp_y = NULL;
  double  ynew  = 0;
  double  tst1 = 0, tst2 = 0;
  double  sum = 0;

  x = calloc (n_old_x, sizeof (double));
  y = calloc (n_old_x, sizeof (double));

  if (!descend) {
    for (i = 0; i < n_old_x; i++) {
      x[i] = (double) old_x[i];
      y[i] = (double) old_y[i];
    }
  } else {
    for (i = 0; i < n_old_x; i++) {
      x[i] = (double) old_x[n_old_x - 1 - i];
      y[i] = (double) old_y[n_old_x - 1 - i];
    }
  }

  /* check if the array is now sorted in ascending order */
  for (i = 0; i < n_old_x - 1; i++)
    if (x[i] >= x[i + 1]) {
      fprintf (stderr, "Error, x not sorted in ascending order\n");
      fprintf (stderr, "x[%d] =  %f, x[%d] = %f\n", i, x[i], i + 1, x[i + 1]);
      CHKERR (-1);
    }

  /* check that range of x_new is <= range of x_old (if not checked -> segfault) */
  /* check minimum */
  for (i = 0; i < n_new_x; i++)
    if (new_x[i] < x[0]) {
      fprintf (stderr, "Error, x_new[%d] = %f < min(x_old) = %f, arb_wvn() (in ancillary.c) \n", i, new_x[i], x[0]);
      CHKERR (-2);
    }
  /* check maximum */
  for (i = 0; i < n_new_x; i++)
    if (new_x[i] > x[n_old_x - 1]) {
      fprintf (stderr, "Error, x_new[%d] = %f > max(x_old) = %f, arb_wvn() (in ancillary.c) \n", i, new_x[i], x[n_old_x - 1]);
      CHKERR (-3);
    }

  /* special check for log_spline interpolation */
  if (linear == 4) {
    sum = 0.0;
    for (i = 0; i < n_old_x; i++) {
      sum += y[i];
      if (y[i] < 0.0) {
        fprintf (stderr, "Error, cannot use log_spline interpolation for profile with negative values\n");
        CHKERR (-4);
      }
    }
    if (sum == 0.0) { /* x is zero in each layer, but grid must be changed anyway */
      free (x);
      free (y);
      for (i = 0; i < n_new_x; i++)
        new_y[i] = 0.0;
      return 0;
    } else {
      for (i = 0; i < n_old_x; i++) {
        if (y[i] > 0.0)
          y[i] = log (y[i]);
        else
          y[i] = -10.E+99; /* this is a little bit cheating, but it works */
      }
    }
  }

  switch (linear) {
  case 0: /* spline */
  case 4: /* log spline */

    status = spline_coeffc (x, y, n_old_x, &a0, &a1, &a2, &a3);
    CHKERR (status);
    break;

  case 1: /* linear */
    status = linear_coeffc (x, y, n_old_x, &a0, &a1, &a2, &a3);
    CHKERR (status);
    break;

  case 2: /* log */
    tmp_y = calloc (n_new_x, sizeof (double));
    for (i = 0; i < n_new_x; i++) {

      j = 0;
      while (new_x[i] > x[j])
        j++;

      if (j > 0)
        j--;

      /* logarithmic interpolation (if reasonable) */
      tst1 = fabs (y[j + 1] - y[j]);
      tst2 = (y[j + 1] < y[j] ? y[j + 1] : y[j]);
      if (tst1 <= 0.001 * y[j] || tst2 <= 0) /* linear */
        tmp_y[i] = y[j] + (new_x[i] - x[j]) / (x[j + 1] - x[j]) * (y[j + 1] - y[j]);
      else /* logarithmic */
        tmp_y[i] = exp (log (y[j]) + (new_x[i] - x[j]) / (x[j + 1] - x[j]) * (log (y[j + 1]) - log (y[j])));
    }

    /* first interpolate, then copy because then source */
    /* and target may be one and the same               */
    for (i = 0; i < n_new_x; i++)
      new_y[i] = (float) tmp_y[i];

    free (tmp_y);
    break;

  case 3: /* linear mixing ratio */
    fprintf (stderr, "Error, linear mixing interpolation not possible with arb_wvn.\n");
    fprintf (stderr, "Please use interpolate_density instead. \n");
    CHKERR (-5);
    break;

  default:
    fprintf (stderr, "Error, unknown interpolation type %d\n", linear);
    CHKERR (-6);
  }

  switch (linear) {
  case 0: /* spline */
  case 1: /* linear */
  case 4: /* log spline */

    for (i = 0; i < n_new_x; i++) {

      if (linear == 1) {
        status = calc_linear_value ((double) new_x[i], &ynew, x, n_old_x, a0, a1);
      } else {
        status = calc_splined_value ((double) new_x[i], &ynew, x, n_old_x, a0, a1, a2, a3);
      }

      /* this check is new in 0.99-alpha-6; so far, no error */
      /* was reported if a value could not be interpolated   */

      if (status != 0) {
        fprintf (stderr, "Error %d returned by calc_splined_value (%g)\n", status, new_x[i]);
        fprintf (stderr, " %d data points, x[0] = %g, x[%d] = %g\n", n_old_x, x[0], n_old_x - 1, x[n_old_x - 1]);
        CHKERR (status);
      }

      new_y[i] = (float) ynew;
    }

    free (a0);
    free (a1);
    free (a2);
    free (a3);
    break;

  case 2:
    /* no need to do anything because interpolation has already been done above */
    break;

  case 3: /* linear mixing ratio */
    fprintf (stderr, "Error, linear mixing interpolation not possible with arb_wvn.\n");
    fprintf (stderr, "Please use interpolate_density instead. \n");
    CHKERR (-1);
    break;

  default:
    fprintf (stderr, "Error, unknown interpolation type %d\n", linear);
    CHKERR (-7);
  }

  if (linear == 4)
    for (i = 0; i < n_new_x; i++)
      new_y[i] = exp (new_y[i]);

  free (x);
  free (y);

  return 0;
}

int arb_wvn_double (int n_old_x, double* old_x, double* old_y, int n_new_x, double* new_x, double* new_y, int linear, int descend) {
  int     i = 0, j = 0, status = 0;
  double *a0 = NULL, *a1 = NULL, *a2 = NULL, *a3 = NULL;
  double *x = NULL, *y = NULL;
  double* tmp_y = NULL;
  double  ynew  = 0;
  double  tst1 = 0, tst2 = 0;
  double  sum = 0;

  x = calloc (n_old_x, sizeof (double));
  y = calloc (n_old_x, sizeof (double));

  if (!descend) {
    for (i = 0; i < n_old_x; i++) {
      x[i] = (double) old_x[i];
      y[i] = (double) old_y[i];
    }
  } else {
    for (i = 0; i < n_old_x; i++) {
      x[i] = (double) old_x[n_old_x - 1 - i];
      y[i] = (double) old_y[n_old_x - 1 - i];
    }
  }

  /* check if the array is now sorted in ascending order */
  for (i = 0; i < n_old_x - 1; i++)
    if (x[i] >= x[i + 1]) {
      fprintf (stderr, "x[%d] =  %f, x[%d] = %f\n", i, x[i], i + 1, x[i + 1]);
      CHKERROUT (-1, "x not sorted in ascending order");
    }

  /* check that range of x_new is <= range of x_old (if not checked -> segfault) */
  /* check minimum */
  for (i = 0; i < n_new_x; i++)
    if (new_x[i] < x[0]) {
      fprintf (stderr, "Error, x_new[%d] = %f < min(x_old) = %f, arb_wvn() (in ancillary.c) \n", i, new_x[i], x[0]);
      CHKERR (-2);
    }
  /* check maximum */
  for (i = 0; i < n_new_x; i++)
    if (new_x[i] > x[n_old_x - 1]) {
      fprintf (stderr, "Error, x_new[%d] = %f > max(x_old) = %f, arb_wvn() (in ancillary.c) \n", i, new_x[i], x[n_old_x - 1]);
      CHKERR (-3);
    }

  /* special check for log_spline interpolation */
  if (linear == 4) {
    sum = 0.0;
    for (i = 0; i < n_old_x; i++) {
      sum += y[i];
      if (y[i] < 0.0) {
        CHKERROUT (-4, "cannot use log_spline interpolation for profile with negative values");
      }
    }
    if (sum == 0.0) { /* x is zero in each layer, but grid must be changed anyway */
      free (x);
      free (y);
      for (i = 0; i < n_new_x; i++)
        new_y[i] = 0.0;
      return 0;
    } else {
      for (i = 0; i < n_old_x; i++) {
        if (y[i] > 0.0)
          y[i] = log (y[i]);
        else
          y[i] = -10.E+99; /* this is a little bit cheating, but it works */
      }
    }
  }

  switch (linear) {
  case 0: /* spline */
  case 4: /* log spline */

    status = spline_coeffc (x, y, n_old_x, &a0, &a1, &a2, &a3);
    CHKERR (status);
    break;

  case 1: /* linear */
    status = linear_coeffc (x, y, n_old_x, &a0, &a1, &a2, &a3);
    CHKERR (status);
    break;

  case 2: /* log */
    tmp_y = calloc (n_new_x, sizeof (double));
    for (i = 0; i < n_new_x; i++) {

      j = 0;
      while (new_x[i] > x[j])
        j++;

      if (j > 0)
        j--;

      /* logarithmic interpolation (if reasonable) */
      tst1 = fabs (y[j + 1] - y[j]);
      tst2 = (y[j + 1] < y[j] ? y[j + 1] : y[j]);
      if (tst1 <= 0.001 * y[j] || tst2 <= 0) /* linear */
        tmp_y[i] = y[j] + (new_x[i] - x[j]) / (x[j + 1] - x[j]) * (y[j + 1] - y[j]);
      else /* logarithmic */
        tmp_y[i] = exp (log (y[j]) + (new_x[i] - x[j]) / (x[j + 1] - x[j]) * (log (y[j + 1]) - log (y[j])));
    }

    /* first interpolate, then copy because then source */
    /* and target may be one and the same               */
    for (i = 0; i < n_new_x; i++)
      new_y[i] = (float) tmp_y[i];

    free (tmp_y);
    break;

  case 3: /* linear mixing ratio */
    CHKERROUT (-5, "linear mixing interpolation not possible with arb_wvn. Please use interpolate_density instead.");
    break;

  default:
    fprintf (stderr, "Error, unknown interpolation type %d\n", linear);
    CHKERROUT (-6, "unknown interpolation type");
  }

  switch (linear) {
  case 0: /* spline */
  case 1: /* linear */
  case 4: /* log spline */

    for (i = 0; i < n_new_x; i++) {

      if (linear == 1) {
        status = calc_linear_value ((double) new_x[i], &ynew, x, n_old_x, a0, a1);
        CHKERR (status);
      } else {
        status = calc_splined_value ((double) new_x[i], &ynew, x, n_old_x, a0, a1, a2, a3);
        CHKERR (status);
      }

      /* this check is new in 0.99-alpha-6; so far, no error */
      /* was reported if a value could not be interpolated   */

      if (status != 0) {
        fprintf (stderr, "Error %d returned by calc_splined_value (%g)\n", status, new_x[i]);
        fprintf (stderr, " %d data points, x[0] = %g, x[%d] = %g\n", n_old_x, x[0], n_old_x - 1, x[n_old_x - 1]);
        CHKERR (status);
      }

      new_y[i] = (float) ynew;
    }

    free (a0);
    free (a1);
    free (a2);
    free (a3);
    break;

  case 2:
    /* no need to do anything because interpolation has already been done above */
    break;

  case 3: /* linear mixing ratio */
    CHKERROUT (-1, "linear mixing interpolation not possible with arb_wvn. Please use interpolate_density instead.");
    break;

  default:
    fprintf (stderr, "Error, unknown interpolation type %d\n", linear);
    CHKERROUT (-7, "unknown interpolation type");
  }

  if (linear == 4)
    for (i = 0; i < n_new_x; i++)
      new_y[i] = exp (new_y[i]);

  free (x);
  free (y);

  return 0;
}

/***************************************************************/
/* Interpolate from one wavelength grid to another, either     */
/* with natural cubic splines (linear=0) or linear (linear=1); */
/* if the input is sorted in descending order, descend needs   */
/* to be set to 1. In contrast to arb_wvn(), values that       */
/* cannot be interpolated will be set to 0.                    */
/***************************************************************/

int arb_wvn_zero (int n_old_x, float* old_x, float* old_y, int n_new_x, float* new_x, float* new_y, int linear, int descend) {
  int     i = 0, status = 0;
  double *a0 = NULL, *a1 = NULL, *a2 = NULL, *a3 = NULL;
  double *x = NULL, *y = NULL;
  double  ynew = 0;

  x = calloc (n_old_x, sizeof (double));
  y = calloc (n_old_x, sizeof (double));

  if (!descend) {
    for (i = 0; i < n_old_x; i++) {
      x[i] = (double) old_x[i];
      y[i] = (double) old_y[i];
    }
  } else {
    for (i = 0; i < n_old_x; i++) {
      x[i] = (double) old_x[n_old_x - 1 - i];
      y[i] = (double) old_y[n_old_x - 1 - i];
    }
  }

  switch (linear) {
  case 0:
    status = spline_coeffc (x, y, n_old_x, &a0, &a1, &a2, &a3);
    CHKERR (status);
    break;

  case 1:
    status = linear_coeffc (x, y, n_old_x, &a0, &a1, &a2, &a3);
    CHKERR (status);
    break;

  case 2:
  default:
    fprintf (stderr, "Error, unknown interpolation type %d\n", linear);
    return -1;
  }

  switch (linear) {
  case 0:
  case 1:
    for (i = 0; i < n_new_x; i++) {

      if (linear == 1)
        status = calc_linear_value ((double) new_x[i], &ynew, x, n_old_x, a0, a1);
      else
        status = calc_splined_value ((double) new_x[i], &ynew, x, n_old_x, a0, a1, a2, a3);

      if (status == 0)
        new_y[i] = (float) ynew;
      else
        new_y[i] = 0;
    }

    free (a0);
    free (a1);
    free (a2);
    free (a3);

    break;

  case 2:
  default:
    fprintf (stderr, "Error, unknown interpolation type %d\n", linear);
    return -1;
  }

  free (x);
  free (y);

  return 0;
}

/***************************************************************/
/* Interpolate from one wavelength grid to another, either     */
/* with natural cubic splines (linear=0) or linear (linear=1); */
/* interpolate only to a selected range of the output grid     */
/* (including n_new_x_lower and n_new_x_upper)                 */
/***************************************************************/

int arb_wvn2 (int    n_old_x,
              float* old_x,
              float* old_y,
              int    n_new_x_lower,
              int    n_new_x_upper,
              float* new_x,
              float* new_y,
              int    linear) {
  int     i = 0, status = 0;
  double *a0 = NULL, *a1 = NULL, *a2 = NULL, *a3 = NULL;
  double *x = NULL, *y = NULL;
  double  ynew = 0;

  x = calloc (n_old_x, sizeof (double));
  y = calloc (n_old_x, sizeof (double));

  for (i = 0; i < n_old_x; i++) {
    x[i] = (double) old_x[i];
    y[i] = (double) old_y[i];
  }

  if (linear) {
    status = linear_coeffc (x, y, n_old_x, &a0, &a1, &a2, &a3);
    CHKERR (status);
  } else {
    status = spline_coeffc (x, y, n_old_x, &a0, &a1, &a2, &a3);
    CHKERR (status);
  }

  for (i = n_new_x_lower; i <= n_new_x_upper; i++) {

    if (linear)
      status = calc_linear_value ((double) new_x[i], &ynew, x, n_old_x, a0, a1);
    else
      status = calc_splined_value ((double) new_x[i], &ynew, x, n_old_x, a0, a1, a2, a3);

    /* this check is new in 0.99-alpha-6; so far, no error */
    /* was reported if a value could not be interpolated   */

    if (status != 0) {
      fprintf (stderr, "Error %d returned by calc_splined_value (%g)\n", status, new_x[i]);
      fprintf (stderr, " %d data points, x[0] = %g, x[%d] = %g\n", n_old_x, x[0], n_old_x - 1, x[n_old_x - 1]);
      CHKERR (status);
    }

    new_y[i] = (float) ynew;
  }

  free (a0);
  free (a1);
  free (a2);
  free (a3);
  free (x);
  free (y);

  return 0;
}

int arb_wvn2_double (int     n_old_x,
                     double* old_x,
                     double* old_y,
                     int     n_new_x_lower,
                     int     n_new_x_upper,
                     double* new_x,
                     double* new_y,
                     int     linear) {
  int     i = 0, status = 0;
  double *a0 = NULL, *a1 = NULL, *a2 = NULL, *a3 = NULL;
  double *x = NULL, *y = NULL;
  double  ynew = 0;

  x = calloc (n_old_x, sizeof (double));
  y = calloc (n_old_x, sizeof (double));

  for (i = 0; i < n_old_x; i++) {
    x[i] = (double) old_x[i];
    y[i] = (double) old_y[i];
  }

  if (linear) {
    status = linear_coeffc (x, y, n_old_x, &a0, &a1, &a2, &a3);
    CHKERR (status);
  } else {
    status = spline_coeffc (x, y, n_old_x, &a0, &a1, &a2, &a3);
    CHKERR (status);
  }

  for (i = n_new_x_lower; i <= n_new_x_upper; i++) {

    if (linear)
      status = calc_linear_value ((double) new_x[i], &ynew, x, n_old_x, a0, a1);
    else
      status = calc_splined_value ((double) new_x[i], &ynew, x, n_old_x, a0, a1, a2, a3);

    /* this check is new in 0.99-alpha-6; so far, no error */
    /* was reported if a value could not be interpolated   */

    if (status != 0) {
      fprintf (stderr, "Error %d returned by calc_splined_value (%g)\n", status, new_x[i]);
      fprintf (stderr, " %d data points, x[0] = %g, x[%d] = %g\n", n_old_x, x[0], n_old_x - 1, x[n_old_x - 1]);
      CHKERR (status);
    }

    new_y[i] = (double) ynew;
  }

  free (a0);
  free (a1);
  free (a2);
  free (a3);
  free (x);
  free (y);

  return 0;
}

/***************************************************/
/* Distribute the photons according to the weights */
/*   when represenative wavelengths are used       */
/***************************************************/
int pfraction_reptran (wl_out_struct* wl_out) {
  int   i, i_t, i_band;
  float weight_sum;

  /* calculate the sum of the weights over all bands */
  weight_sum = 0;
  for (i_t = 0; i_t < wl_out->nlambda_t; i_t++) {

    i_band = wl_out->reptran_band_t[i_t];

    /* Test if previous _t wavelength was in the same band. */
    /* If yes, go to next _t wavelength since we do not     */
    /* want to add the weights for a band multiple times.   */
    if (i_t > 0 && (i_band == wl_out->reptran_band_t[i_t - 1]))
      continue;

    for (i = 0; i < wl_out->nlambda_in_reptran_band[i_band]; i++)
      weight_sum += fabs (wl_out->weight_reptran_band[i_band][i]) * wl_out->extra_reptran_r[wl_out->reptran_band[i_band][i]];
  }

  /* Calculate the pfraction */
  for (i_t = 0; i_t < wl_out->nlambda_t; i_t++) {

    i_band = wl_out->reptran_band_t[i_t];

    for (i = 0; i < wl_out->nlambda_in_reptran_band[i_band]; i++)
      wl_out->pfraction[wl_out->reptran_band[i_band][i]] =
        fabs (wl_out->weight_reptran_band[i_band][i]) * wl_out->extra_reptran_r[wl_out->reptran_band[i_band][i]] / weight_sum;
  }

  return 0;
}

/**************************************************************************************************************/
/* Calculate results at transmission grid (results_t) from results at representative wavelengths (results_r). */
/* In case of solar source, the results are weighted with the weights (weight_reptran_band) and the            */
/*   extraterrestrial spectrum (extra_reptran_r) in accordance with the approach taken during finding          */
/*   the representative wavelengths.                                                                          */
/* In case of thermal source, only the weights (weight_reptran_band) are relevant                              */
/*   (extra_reptran_r was set to 1, it is neutral).                                                            */
/* Set is_variance=1 for weighting variances, and is_variance=0 for weighting other quantities                */
/**************************************************************************************************************/
int weighting_reptran (const wl_out_struct* wl_out, int is_variance, float* results_r, float* results_t, int reptran_option) {
  int   i, i_t, i_band;
  float weight;
  float sum;
  float weight_sum;

  for (i_t = 0; i_t < wl_out->nlambda_t; i_t++) {

    i_band = wl_out->reptran_band_t[i_t];

    sum        = 0;
    weight_sum = 0;

    for (i = 0; i < wl_out->nlambda_in_reptran_band[i_band]; i++) {

      weight = wl_out->weight_reptran_band[i_band][i] * wl_out->extra_reptran_r[wl_out->reptran_band[i_band][i]];

      if (is_variance)
        sum += results_r[wl_out->reptran_band[i_band][i]] * weight * weight;
      else {
        sum += results_r[wl_out->reptran_band[i_band][i]] * weight;
      }

      weight_sum += weight;
    }

    /* Integral must not be devided by weight_sum for "reptran int" because weights are not normalized. */
    if (reptran_option == REPTRAN_OPTION_INT)
      weight_sum = 1.0;

    if (is_variance)
      results_t[i_t] = sum / (weight_sum * weight_sum);
    else
      results_t[i_t] = sum / weight_sum;
  }

  return 0;
}

/************************************************************************************************/
/* This function returns the name of the file with the representative wavelengths (if i_mol<=0) */
/*   or the name of the corresponding absorption lookup table files (if i_mol>0).               */
/************************************************************************************************/
int reptran_filename (input_struct input, int i_mol, char* filename) {

  int   len;
  char* gas = NULL;

  if (strlen (input.filename[FN_REPTRAN]) > 0)

    strcpy (filename, input.filename[FN_REPTRAN]);

  else {

    strcpy (filename, input.filename[FN_PATH]);

    if (input.ck_reptran_option == REPTRAN_OPTION_INT) {
      strcat (filename, "correlated_k/repint/");
      switch (input.ck_repint_option) {
      case REPINT_DEMOURGUES10:
        strcat (filename, "demourgues10/");
        break;
      case REPINT_DEMOURGUES30:
        strcat (filename, "demourgues30/");
        break;
      case REPINT_DEMOURGUES50:
        strcat (filename, "demourgues50/");
        break;
      case REPINT_DEMOURGUES100:
        strcat (filename, "demourgues100/");
        break;
      case REPINT_MEIER50:
        strcat (filename, "meier50/");
        break;
      case REPINT_MEIER100:
        strcat (filename, "meier100/");
        break;
      default:
        fprintf (stderr, "Error, unknown REPINT option\n");
        return -1;
      }
    } else
      strcat (filename, "correlated_k/reptran/");

    strcat (filename, "reptran_");

    if (input.source == SRC_THERMAL)
      strcat (filename, "thermal_");
    else if (input.source == SRC_SOLAR)
      strcat (filename, "solar_");
    else
      return err_out ("Error: Unsupported source in reptran_filename().\n", -1);

    if (input.ck_scheme == CK_REPTRAN) {
      if (input.ck_reptran_option == REPTRAN_OPTION_FINE)
        strcat (filename, "fine");
      else if (input.ck_reptran_option == REPTRAN_OPTION_MEDIUM)
        strcat (filename, "medium");
      else if (input.ck_reptran_option == REPTRAN_OPTION_INT)
        strcat (filename, "int");
      else if (input.ck_reptran_option == REPTRAN_OPTION_COARSE || input.ck_reptran_option == REPTRAN_OPTION_NONE)
        strcat (filename, "coarse");
      else
        return err_out ("Error: Unknown ck_reptran_option in function reptran_filename().\n", -1);
    } else if (input.ck_scheme ==
               CK_REPTRAN_CHANNEL) {  // reptran filename contains only first part of channel name until first underscore;
      len = strlen (input.ck_reptran_channel) - strlen (strchr (input.ck_reptran_channel, '_'));
      if (
        isdigit (
          input.ck_reptran_channel
            [len -
             1]))  // reptran filename does not contain the last or the last two characters if there are digits (e.g. sentinel3 --> sentinel or sentinel2a --> sentinel)
        len = len - 1;
      else if (isdigit (input.ck_reptran_channel[len - 2]))
        len = len - 2;
      strncat (filename, input.ck_reptran_channel, len);
    }
  }

  if (i_mol > 0) {

    strcat (filename, ".lookup.");
    gas = gas_number2string (i_mol);
    strcat (filename, gas);
    if (filename[strlen (filename) - 1] == ' ')
      filename[strlen (filename) - 1] = '\0'; /* remove upto two space characters from the species names */
    if (filename[strlen (filename) - 1] == ' ')
      filename[strlen (filename) - 1] = '\0';

    free (gas);
  }

  strcat (filename, ".cdf");

  return 0;
}

/**************************************************************/
/* Interpolate a given profile (x,y) to a new grid xnew;      */
/* data are written to original array y and memory is         */
/* automatically reallocated                                  */
/* automatically check for negative values                    */
/*                                                            */
/* Ulrich Hamann                                              */
/**************************************************************/

int interpolate_profile (float* x, float** y, int n, float* xnew, int nnew, int interpol_method, int quiet) {
  int    status = 0;
  float* tmp    = NULL;
  int    lc = 0, lc2 = 0, first = 0;

  tmp = calloc (nnew, sizeof (float));

  switch (interpol_method) {
  case INTERP_METHOD_SPLINE:
  case INTERP_METHOD_LINEAR:
  case INTERP_METHOD_LOG:
  case INTERP_METHOD_LOG_SPLINE:
    status = arb_wvn (n, x, *y, nnew, xnew, tmp, interpol_method, 1);
    break;

  case INTERP_METHOD_LINMIX:
    CHKERROUT (-1, "linear mixing ratio not possible with interpolate_profile Please use interpolate_density instead!");
    break;

  default:
    fprintf (stderr, "Error, unknown interpolation method in interpolate_profile\n");
    return -1;
  }

  /* testing for negativ density values */

  switch (interpol_method) {
  case INTERP_METHOD_SPLINE:
    for (lc = 0; lc < nnew; lc++) {
      if (tmp[lc] < 0.0) {
        if (!quiet) {
          if (first == 0)
            fprintf (stderr, "*** Warning: In interpolate_density, automatic correction of negative density to 0.0\n");
          fprintf (stderr, "*** Warning: z[%d]= %5.1f, dens[%d]= %12.7e -> dens[%d]= 0.0\n", lc, x[lc], lc, tmp[lc], lc);
          tmp[lc] = 0.0;
          first   = 1;
        }
      }
    }
    break;

  case INTERP_METHOD_LINEAR:
  case INTERP_METHOD_LOG:
  case INTERP_METHOD_LINMIX:
  case INTERP_METHOD_LOG_SPLINE:
    for (lc = 0; lc < nnew; lc++) {
      if (tmp[lc] < 0.0) {
        if (fabs (tmp[lc]) < EPSILON) {
          if (!quiet) {
            fprintf (stderr, "*** Warning, small negative density detected during interpolate_profile!\n");
            fprintf (stderr, "*** in layer = %d, dens = %e.\n", lc, tmp[lc]);
            fprintf (stderr, "*** This may happen, when dens profiles contain 0.0 values.\n");
            fprintf (stderr, "*** Setting to 0.0 automatically.\n");
          }
          tmp[lc] = 0.0;
        } else {
          fprintf (stderr, "Error, negative density detected in interpolate_profile! lc=%d \n", lc);
          for (lc2 = 0; lc2 < nnew; lc2++)
            fprintf (stderr, " ### dens[%d]=%e\n", lc2, tmp[lc2]);
          return -1;
        }
      }
    }
    break;
  default:
    fprintf (stderr, "Error, unknown interpolation method in interpolate_profile\n");
    return -1;
  }

  CHKERR (status);

  free (*y);
  *y = tmp;

  return 0;
}

/**************************************************************/
/* Interpolate a given profile (x,y) to a new grid xnew;      */
/* data are written to original array y and memory is         */
/* automatically reallocated                                  */
/* automatically check for negative values                    */
/*                                                            */
/* Ulrich Hamann                                              */
/**************************************************************/

int interpolate_density (float*  x,
                         float** y,
                         int     n,
                         float*  xnew,
                         int     nnew,
                         int     interpol_method,
                         float*  dens_air_old,
                         float*  dens_air,
                         int     quiet) {
  int    status;
  float* tmp = calloc (nnew, sizeof (float));

  switch (interpol_method) {
  case INTERP_METHOD_SPLINE:
  case INTERP_METHOD_LINEAR:
  case INTERP_METHOD_LOG:
  case INTERP_METHOD_LOG_SPLINE:
    status = arb_wvn (n, x, *y, nnew, xnew, tmp, interpol_method, 1);
    CHKERR (status);
    break;

  case INTERP_METHOD_LINMIX:
    /* convert from density to mixing ratio */
    for (int lc = 0; lc < n; lc++) {
      if (dens_air_old[lc] <= 0.0) {
        CHKERROUT (-1, "cannot use linmix interpolation, when air density is <= 0.0");
      } else {
        (*y)[lc] = (*y)[lc] / dens_air_old[lc];
      }
    }
    /*linear interpolation of the mixing ratio*/
    status = arb_wvn (n, x, *y, nnew, xnew, tmp, INTERP_METHOD_LINEAR, 1);
    CHKERR (status);

    /* convert back to number density */
    for (int lc = 0; lc < nnew; lc++)
      tmp[lc] = tmp[lc] * dens_air[lc];
    break;
  default:
    CHKERROUT (-1, "unknown interpolation method in interpolate_profile");
  }

  /* testing for negativ density values */

  int first = 0;
  switch (interpol_method) {
  case INTERP_METHOD_SPLINE:
    for (int lc = 0; lc < nnew; lc++) {
      if (tmp[lc] < 0.0) {
        if (!quiet) {
          if (first == 0)
            fprintf (stderr, "***  Warning: In interpolate_density, automatic correction of negative density to 0.0\n");
          fprintf (stderr, "***  Warning: z[%3d]= %5.1f, dens[%3d]= %12.7e -> dens[%3d]= 0.0\n", lc, x[lc], lc, tmp[lc], lc);
          tmp[lc] = 0.0;
          first   = 1;
        }
      }
    }
    break;

  case INTERP_METHOD_LINEAR:
  case INTERP_METHOD_LOG:
  case INTERP_METHOD_LINMIX:
  case INTERP_METHOD_LOG_SPLINE:
    for (int lc = 0; lc < nnew; lc++) {
      if (tmp[lc] < 0.0) {
        if (fabs (tmp[lc]) < EPSILON) {
          if (!quiet) {
            fprintf (stderr, "*** Warning, small negative density detected during interpolate_density!\n");
            fprintf (stderr, "*** in layer = %d, dens = %e, abs(dens)= %e\n", lc, tmp[lc], fabs (tmp[lc]));
            fprintf (stderr, "*** This may happen, when dens profiles contain 0.0 values.\n");
            fprintf (stderr, "*** Setting to 0.0 automatically.\n");
          }
          tmp[lc] = 0.0;
        } else {
          fprintf (stderr, "Error, negative density detected in interpolate_density!, lc=%3d \n", lc);
          for (int lc2 = 0; lc2 < nnew; lc2++)
            fprintf (stderr, " ### z[%3d] = %7.2f, dens[%3d] = %e\n", lc2, xnew[lc2], lc2, tmp[lc2]);
          CHKERR (-1);
        }
      }
    }
    break;
  default:
    CHKERROUT (-1, "unknown interpolation method in interpolate_profile");
  }

  CHKERR (status);

  free (*y);
  *y = tmp;

  return 0;
}

/*******************************************************************/
/* Interpolate all given atmospheric profiles to a new z-grid      */
/* profiles are written to original arrays and memory is           */
/* automatically reallocated                                       */
/*******************************************************************/

int interpolate_atmosphere (float*     z,
                            float****  p_p,
                            float****  p_T,
                            float***** p_dens,
                            float****  p_Tavg,
                            float***** p_densavg,
                            int        n,
                            float*     znew,
                            int        nnew,
                            int        interpol_method_press,
                            int        interpol_method_temper,
                            int*       interpol_method_gas,
                            int        quiet)

{
  /* n is old number of levels, nnew is new number of levels */

  /* "p_xxx"  here means "pointer to xxx" */

  int    status;
  int    lc = -999, i = -999;
  float* dens_air_old = NULL;

  /* Pressure */
  status = interpolate_profile (z, p_p[0][0], n, znew, nnew, interpol_method_press, quiet);
  CHKERR (status);

  /* Temperature */
  status = interpolate_profile (z, p_T[0][0], n, znew, nnew, interpol_method_temper, quiet);
  CHKERR (status);

  /* copy old air number density */
  dens_air_old = calloc (n, sizeof (float));
  for (lc = 0; lc < n; lc++)
    dens_air_old[lc] = (*p_dens)[MOL_AIR][0][0][lc];

  /* Air */
  /* this is not 100% consistent with interpolation of p and T */
  status = interpolate_profile (z, &((*p_dens)[MOL_AIR][0][0]), n, znew, nnew, interpol_method_gas[MOL_AIR], quiet);
  CHKERR (status);

  /*   /\* determine air number density from pressure and temperature *\/ */
  /*   free((*p_dens)[MOL_AIR]); */
  /*   (*p_dens)[MOL_AIR] = calloc (nnew, sizeof(float)); */
  /*   for (lc=0; lc<nnew; lc++) */
  /*     (*p_dens)[MOL_AIR][lc] =  (*p_p)[lc] / (BOLTZMANN * (*p_T)[lc]) * 100.0 / 1e6; */

  /* Ozone, O2, water vapour, CO2, NO2, BRO, OClO, HCHO, O4 */
  for (i = 0; i < MOL_NN; i++) {
    if (i != MOL_AIR) {
      status = interpolate_density (z,
                                    &((*p_dens)[i][0][0]),
                                    n,
                                    znew,
                                    nnew,
                                    interpol_method_gas[i],
                                    dens_air_old,
                                    (*p_dens)[MOL_AIR][0][0],
                                    quiet);
      CHKERR (status);
    }
  }

  /* recalculate layer average temperature and densities   */
  /* they are needed to convert heating rates to K_per_day */

  if (*p_Tavg != NULL)
    ASCII_free_float_3D ((*p_Tavg), 1, 1);
  ASCII_calloc_float_3D (&(*p_Tavg), 1, 1, nnew);

  status = average_dens ((*p_T)[0][0], (*p_dens)[MOL_AIR][0][0], znew, nnew, interpol_method_temper, &((*p_Tavg)[0][0]), NO);
  CHKERR (status);

  /* We need this check because otherwise the cloud overlap examples crash! */
  /* It seems that memory for the average densities is not allocated        */
  /* correctly in that case.                                                */

  if (*p_densavg != NULL) {
    for (i = 0; i < MOL_NN; i++) {
      if (*p_densavg != NULL)
        free ((*p_densavg)[i][0][0]);

      status = average_dens ((*p_dens)[i][0][0],
                             (*p_dens)[MOL_AIR][0][0],
                             znew,
                             nnew,
                             interpol_method_gas[i],
                             &(*p_densavg)[i][0][0],
                             YES);
      CHKERR (status);
    }
  }

  free (dens_air_old);
  return 0;
}

/******************************************************************/
/* Define the internal wavelength grid for the radiative transfer */
/* calculation.                                                   */
/******************************************************************/

static int set_transmittance_wl_grid (wl_inp_struct* wl_inp, wl_out_struct* wl_out, int quiet) {
  int        iv     = 0;
  float      lambda = 0.0, lambdanew = 0.0;
  static int firstsr = 1, firsthz = 1;

  /* Determine number of wavelengths needed for transmittance calculation;  */
  /* careful - if something is changed, it needs to be changed twice in the */
  /* following code!                                                        */
  lambda = wl_inp->start;
  iv     = 0;

  if (wl_inp->start != wl_inp->end) { /* non-monochromatic calculation */
    while (lambda < wl_inp->end) {

      lambdanew = lambda;

      if (lambda < 121.0) {
        lambdanew += 1.0;
        if (lambdanew > 121.0) {
          lambda    = 121.0;
          lambdanew = 121.0;
        }
      }

      /* Lyman alpha */
      if (lambda >= 121.0 && lambda < 122.0)
        lambdanew += 0.01;

      /* Schumann-Runge continuum */
      if (lambda >= 122.0 && lambda < 130) {
        lambdanew += 0.1;
        if (lambdanew >= 1.0E7 / 57000.0)
          lambda += 0.1;
      }

      if (lambda >= 130.0 && lambda < 1.0E7 / 57000.0) {
        lambdanew += 0.5;
        if (lambdanew >= 1.0E7 / 57000.0)
          lambda += 0.5;
      }

      /* Schumann-Runge bands */
      if (lambda <= 1.0E7 / 49000.5 && lambda >= 1.0E7 / 57000.0) {
        if (firstsr) {
          if (iv > 0)
            lambdanew = 1.0E7 / 57000.0;
          else
            lambdanew = 1.0E7 / (floor (2.0 * 1.0E7 / lambda) / 2.0);

          firstsr = 0;
        } else
          lambdanew = 1.0E7 / (1.0E7 / lambda - 0.5);
      }

      /* Herzberg continuum, Hartley-Huggins bands */
      if (lambda > 1.0E7 / 49000.5 && lambda < 350.0) {
        if (firsthz) {
          lambdanew = ceil (2.0 * lambda) / 2.0;
          if (lambdanew == lambda)
            lambdanew += 0.5;
          firsthz = 0;
        } else
          lambdanew += 0.5;
      }

      if (lambda >= 350.0)
        lambdanew += 1.0;

      lambda = lambdanew;

      iv++;
    }
  }

  firstsr           = 1;
  firsthz           = 1;
  wl_out->nlambda_t = iv + 1;
  wl_out->lambda_t  = calloc (wl_out->nlambda_t, sizeof (float));

  /* Set wavelengths for radiative transfer calculation */
  lambda = wl_inp->start;

  for (iv = 0; iv < wl_out->nlambda_t; iv++) {
    wl_out->lambda_t[iv] = lambda;
    lambdanew            = lambda;

    if (lambda < 121.0) {
      lambdanew += 1.0;
      if (lambdanew > 121.0) {
        lambda    = 121.0;
        lambdanew = 121.0;
      }
    }

    /* Lyman alpha */
    if (lambda >= 121.0 && lambda < 122.0)
      lambdanew += 0.01;

    /* Schumann-Runge continuum */
    if (lambda >= 122.0 && lambda < 130) {
      lambdanew += 0.1;
      if (lambdanew >= 1.0E7 / 57000.0)
        lambda += 0.1;
    }

    if (lambda >= 130.0 && lambda < 1.0E7 / 57000.0) {
      lambdanew += 0.5;
      if (lambdanew >= 1.0E7 / 57000.0)
        lambda += 0.5;
    }

    /* Schumann-Runge bands */
    if (lambda <= 1.0E7 / 49000.5 && lambda >= 1.0E7 / 57000.0) {
      if (firstsr) {
        if (iv > 0)
          lambdanew = 1.0E7 / 57000.0;
        else { /* start of spectrum */
          lambdanew = 1.0E7 / (floor (2.0 * 1.0E7 / lambda) / 2.0);
          if (lambdanew == lambda) /* avoid that we use the same wavelength twice */
            lambdanew = 1.0E7 / (1.0E7 / lambda - 0.5);
        }

        firstsr = 0;
      } else
        lambdanew = 1.0E7 / (1.0E7 / lambda - 0.5);
    }

    /* Herzberg continuum, Hartley-Huggins bands */
    if (lambda > 1.0E7 / 49000.5 && lambda < 350.0) {
      if (firsthz) {
        lambdanew = ceil (2.0 * lambda) / 2.0;
        if (lambdanew == lambda)
          lambdanew += 0.5;
        firsthz = 0;
      } else
        lambdanew += 0.5;
    }

    if (lambda >= 350.0)
      lambdanew += 1.0;

    lambda = lambdanew;
  }

  return 0;
}

static int set_raman_wl_grid (wl_inp_struct* wl_inp, wl_out_struct* wl_out, char* filename, int quiet) {
  int    iv = 0, i = 0;
  int    nlambda = 0, status = 0;
  int    ivs = 0, ive = 0; /* Start and end indices */
  float *tmp_lambda = NULL, *tmp_fbeam = NULL;

  /* read extraterrestrial irradiance */
  status = read_2c_file_float (filename, &tmp_lambda, &tmp_fbeam, &nlambda);
  CHKERR (status);

  wl_out->nlambda_t = 0;
  for (iv = 0; iv < nlambda; iv++) {
    if (tmp_lambda[iv] < wl_inp->start)
      ivs = iv + 1;
    if (tmp_lambda[iv] < wl_inp->end)
      ive = iv + 1;
    if (tmp_lambda[iv] >= wl_inp->start && tmp_lambda[iv] <= wl_inp->end) {
      wl_out->nlambda_t++;
    }
  }
  wl_out->nlambda_t = ive - ivs + 1;
  wl_out->lambda_t  = calloc (wl_out->nlambda_t, sizeof (float));
  i                 = 0;
  for (iv = ivs; iv <= ive; iv++) {
    wl_out->lambda_t[i] = tmp_lambda[iv];
    i++;
  }

  free (tmp_lambda);
  free (tmp_fbeam);

  return 0;
}

static int find_raman_closest_wl_in_solar_file (float* wls, float* wle, char* filename, int quiet) {
  int    iv      = 0;
  int    nlambda = 0, status = 0;
  int    ivs = 0, ive = 0; /* Start and end indices */
  float *tmp_lambda = NULL, *tmp_fbeam = NULL;

  /* read extraterrestrial irradiance */
  status = read_2c_file_float (filename, &tmp_lambda, &tmp_fbeam, &nlambda);
  CHKERR (status);

  for (iv = 0; iv < nlambda; iv++) {
    if (tmp_lambda[iv] < *wls)
      ivs = iv;
    if (tmp_lambda[iv] < *wle)
      ive = iv + 1;
  }

  *wls = tmp_lambda[ivs];
  *wle = tmp_lambda[ive];

  free (tmp_lambda);
  free (tmp_fbeam);

  return 0;
}

/******************************************************************/
/* Define the internal wavelength grid for the radiative transfer */
/* calculation. LOWTRAN requires much less detail, in             */
/* particular below 300nm                                         */
/******************************************************************/

static int set_transmittance_wl_grid_lowtran (wl_inp_struct* wl_inp, wl_out_struct* wl_out, int quiet) {
  int   iv     = 0;
  float lambda = 0.0, wvn_step_t = 0.0;

  /* Determine number of wavelengths needed for radiative transfer calculation */
  lambda = wl_inp->start;
  iv     = 1;
  while (lambda < wl_inp->end) {

    if (lambda < 350.0)
      wvn_step_t = 0.5;
    else
      wvn_step_t = 1.0;

    lambda += wvn_step_t;
    iv++;
  }

  wl_out->nlambda_t = iv;
  wl_out->lambda_t  = calloc (wl_out->nlambda_t, sizeof (float));

  /* Set wavelengths for radiative transfer calculation */
  lambda = wl_inp->start;

  for (iv = 0; iv < wl_out->nlambda_t; iv++) {
    wl_out->lambda_t[iv] = lambda;

    if (lambda < 350.0)
      wvn_step_t = 0.5;
    else
      wvn_step_t = 1.0;

    lambda += wvn_step_t;
  }

  return 0;
}

/***********************************************************************************/
/* Define the transmittance wavelength grid when using representative wavelengths. */
/* The transmittance wavelengths are the center wavelengths of the bands.          */
/***********************************************************************************/
static int
set_transmittance_wl_grid_reptran (input_struct input, float** lambda_lower, float** lambda_upper, wl_out_struct* wl_out) {

#if HAVE_NETCDF4
  int status            = 0;
  int nbands            = 0;
  int max_len_band_name = 0;

  double* wvlmin;
  double* wvlmax;

  char** band_name;

  char cdf_filename[FILENAME_MAX] = "";
  int  ncid                       = 0;

  int idd_nbands            = 0;
  int idd_max_len_band_name = 0;

  int id_wvlmin    = 0;
  int id_wvlmax    = 0;
  int id_band_name = 0;

  size_t dimlen = 0;

  int i_band;
  int i_band_start = -1;
  int i_band_end   = -1;
  int it;

  size_t start[] = {0, 0};
  size_t count[] = {1, 1};

  /* determine filename with the parameterization */
  status = reptran_filename (input, -1, cdf_filename);
  CHKERR (status);

  if (input.verbose)
    fprintf (stderr, "     reading bands from %s.\n", cdf_filename);

  /* open netcdf file and read the minimum and maximum wavelengths of the bands */
  status = nc_open (cdf_filename, NC_NOWRITE, &ncid);
  if (status == NC_NOERR) {

    status = NC_NOERR;
    status += nc_inq_dimid (ncid, "nbands", &idd_nbands);
    status += nc_inq_dimid (ncid, "max_len_band_name", &idd_max_len_band_name);

    status += nc_inq_dimlen (ncid, idd_nbands, &dimlen);
    nbands = dimlen;
    status += nc_inq_dimlen (ncid, idd_max_len_band_name, &dimlen);
    max_len_band_name = dimlen;

    status += nc_inq_varid (ncid, "wvlmin", &id_wvlmin);
    status += nc_inq_varid (ncid, "wvlmax", &id_wvlmax);
    CHKERROUT (status - NC_NOERR, "failed reading the representative wavelengths file");

    wvlmin = calloc (nbands, sizeof (double));
    wvlmax = calloc (nbands, sizeof (double));

    status += nc_get_var_double (ncid, id_wvlmin, wvlmin);
    status += nc_get_var_double (ncid, id_wvlmax, wvlmax);

    status += nc_inq_varid (ncid, "band_name", &id_band_name);

    ASCII_calloc_char (&band_name, nbands, max_len_band_name);

    count[1] = max_len_band_name;
    for (i_band = 0; i_band < nbands; i_band++) {
      start[0] = i_band;
      status += nc_get_vara_text (ncid, id_band_name, start, count, band_name[i_band]);
    }
    CHKERROUT (status - NC_NOERR, "failed reading the representative wavelengths file");

    nc_close (ncid);

  } else {
    if (status) {
      fprintf (stderr, "**********************************************************************************\n");
      fprintf (stderr, "*Error: Data files for REPTRAN not found in directory data/correlated_k/reptran. *'\n");
      fprintf (stderr, "*       Please check whether you have downloaded the required REPTRAN data files *\n");
      fprintf (stderr, "*       from http://www.libradtran.org/doku.php?id=download and unzipped the data*\n");
      fprintf (stderr, "*       in the libRadtran folder.                                                *\n");
      fprintf (stderr, "**********************************************************************************\n");
    }
    CHKERROUT (status - NC_NOERR, "failed opening the representative wavelengths file");
  }

  if (wl_out->type == WLGRID_UVSPEC) {

    if (input.ck_scheme == CK_REPTRAN_CHANNEL) { /* channel is searched here in the netcdf file */

      for (i_band = 0; i_band < nbands; i_band++)

        if (strncasecmp (input.ck_reptran_channel, band_name[i_band], strlen (input.ck_reptran_channel)) == 0) {

          i_band_start = i_band;
          i_band_end   = i_band;

          /* setting the wavelength range to the wavelength range required for selected channel */
          wl_out->start = wvlmin[i_band];
          wl_out->end   = wvlmax[i_band];
        }

      if (i_band_start < 0 || i_band_end < 0)
        CHKERROUT (-1, "Error: Channel not found in reptran_file");

    } else if (input.ck_reptran_option == REPTRAN_OPTION_INT) {
      i_band       = 0;
      i_band_start = i_band;
      i_band_end   = i_band;

      /* setting the wavelength range to the wavelength range required for selected channel */
      wl_out->start = wvlmin[i_band];
      wl_out->end   = wvlmax[i_band];
    }

    else if (input.wl.start > 0) {

      /* select bands required for user-specified wavelength range */
      i_band_start = nbands;
      for (i_band = nbands - 1; i_band >= 0; i_band--)
        if (input.wl.start < (float) wvlmax[i_band])
          i_band_start = i_band;

      i_band_end = -1;
      for (i_band = 0; i_band < nbands; i_band++)
        if (input.wl.end > (float) wvlmin[i_band])
          i_band_end = i_band;

      if (input.wl.start == input.wl.end &&
          i_band_start >
            i_band_end) /* special case when a) only a single wavelength was specified and b) this wavelength is a band boundary */
        i_band_end = i_band_start;

      if (i_band_start == nbands || i_band_end == -1 || input.wl.end > (float) wvlmax[nbands - 1] ||
          input.wl.start < (float) wvlmin[0]) {
        fprintf (stderr, "*****************************************************************\n");
        fprintf (stderr, "Error: User-specified wavelength range not covered by REPTRAN.\n");
        fprintf (stderr, "       Wavelength range covered by REPTRAN is from %f nm to %f nm\n", wvlmin[0], wvlmax[nbands - 1]);
        fprintf (stderr, "*****************************************************************\n");
        return err_out ("Error: User-specified wavelength range not covered by the representative wavelengths parameterization.\n",
                        -1);
      }
    } else if (input.wl.start_index > 0) {

      i_band_start = input.wl.start_index - 1;
      i_band_end   = input.wl.end_index - 1;

      if (i_band_end > nbands - 1)
        return err_out ("Error: Wavelength index too large for activated representative wavelengths parameterization.\n", -1);

      wl_out->start = wvlmin[i_band_start];
      wl_out->end   = wvlmax[i_band_end];

    } else {

      fprintf (stderr, "***************************************************************************************\n");
      fprintf (stderr, "Error: No wavelength range selected. Please use options wavelength or wavelength_index.\n");
      fprintf (stderr, "***************************************************************************************\n");
      return err_out ("Error: No wavelength range selected.\n", -1);

      /*
	i_band_start=0;
	i_band_end=nbands-1;
	wl_out->start=wvlmin[i_band_start];
	wl_out->end=wvlmax[i_band_end];
      */
    }

    if (wl_out->start < 1.0e7 / 49000.5) {
      fprintf (stderr, "**************************************************************************************\n");
      fprintf (stderr, "Warning: The wavelength resolution of REPTRAN at wavelength < 204.1nm might be too low\n");
      fprintf (stderr, "         for fully resolving the spectral absorption features.\n");
      fprintf (stderr, "         Use 'mol_abs_param crs' if you need higher spectral resolution.\n");
      fprintf (stderr, "**************************************************************************************\n");
    }

    wl_out->nlambda_t = i_band_end - i_band_start + 1;
    wl_out->lambda_t  = calloc (wl_out->nlambda_t, sizeof (float));

    *lambda_lower = calloc (wl_out->nlambda_t, sizeof (float));
    *lambda_upper = calloc (wl_out->nlambda_t, sizeof (float));

    wl_out->reptran_band_t = calloc (wl_out->nlambda_t, sizeof (int));

    for (i_band = i_band_start; i_band <= i_band_end; i_band++) {

      wl_out->lambda_t[i_band - i_band_start] = (wvlmin[i_band] + wvlmax[i_band]) / 2;

      (*lambda_lower)[i_band - i_band_start]        = wvlmin[i_band];
      (*lambda_upper)[i_band - i_band_start]        = wvlmax[i_band];
      wl_out->reptran_band_t[i_band - i_band_start] = i_band;
    }

  } else if (wl_out->type == WLGRID_USER) {

    if (input.wl.start > 0)
      return err_out ("Error: Combination of 'mol_abs_param reptran' and 'wavelength_grid_file' with 'wavelength' not supported.\n",
                      -1);

    /* only need to find suitable bands for the transmission wavelengths given by the user */
    wl_out->reptran_band_t = calloc (wl_out->nlambda_t, sizeof (int));
    for (it = 0; it < wl_out->nlambda_t; it++) {

      wl_out->reptran_band_t[it] = -1;

      for (i_band = 0; i_band < nbands; i_band++)
        if (wl_out->lambda_t[it] >= wvlmin[i_band] && wl_out->lambda_t[it] < wvlmax[i_band])
          wl_out->reptran_band_t[it] = i_band;

      if (wl_out->reptran_band_t[it] == -1)
        return err_out (
          "Error: Wavelengths in wavelength_grid_file not covered by range of representative wavelengths parameterization.\n",
          -1);
    }
  }

  if (input.verbose)
    fprintf (stderr, "       %d wavelengths set by set_transmittance_wl_grid_reptran().\n", wl_out->nlambda_t);

  ASCII_free_char (band_name, nbands);
  free (wvlmin);
  free (wvlmax);

  return 0;
#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
  fprintf (stderr, " ***********************************************************************\n");
  return -1;
#endif
}

/*****************************************************************************/
/* Define the rte wavelength grid when using the representative wavelengths. */
/*****************************************************************************/
static int set_rte_wl_grid_reptran (input_struct input, output_struct* output) {
#if HAVE_NETCDF4

  int status           = 0;
  int nbands           = 0;
  int nwvl             = 0;
  int max_nwvl_in_band = 0;

  char cdf_filename[FILENAME_MAX] = "";
  int  ncid                       = 0;

  int idd_nwvl             = 0;
  int idd_nbands           = 0;
  int idd_max_nwvl_in_band = 0;
  int id_wvl               = 0;
  int id_iwvl              = 0;
  int id_weight_wvl        = 0;
  int id_nwvl_in_band      = 0;
  int id_extra             = 0;
  int id_wvl_integral      = 0;

  size_t dimlen = 0;

  double* wvl;
  int*    wvl_active;
  int     n_active;
  int     i_wvl;
  int*    index_in_lambda_r;
  double* extra;

  int no_extra;

  int i_band;
  int i;

  int*    i_wvl_tmp;
  double* weight_wvl_tmp;

  int isp;
  int add_tau_wvl_to_lambda_r;

  double wvl_min_tmp = NOT_A_NUMBER;
  double wvl_max_tmp = NOT_A_NUMBER;
  double tau_wvl;

  /* determine file with the parameterization */
  status = reptran_filename (input, -1, cdf_filename);
  CHKERR (status);

  /* open this file for reading */
  status = nc_open (cdf_filename, NC_NOWRITE, &ncid);
  CHKERR (status);

  /* read dimensions */
  status = 0;
  status += nc_inq_dimid (ncid, "nwvl", &idd_nwvl);
  status += nc_inq_dimid (ncid, "nbands", &idd_nbands);
  status += nc_inq_dimid (ncid, "max_nwvl_in_band", &idd_max_nwvl_in_band);

  status += nc_inq_dimlen (ncid, idd_nwvl, &dimlen);
  nwvl = dimlen;
  status += nc_inq_dimlen (ncid, idd_nbands, &dimlen);
  nbands = dimlen;
  status += nc_inq_dimlen (ncid, idd_max_nwvl_in_band, &dimlen);
  max_nwvl_in_band = dimlen;
  if (status != 0) {
    fprintf (stderr, "Error %d reading dimensions from %s\n", status, cdf_filename);
    CHKERR (status);
  }

  status = nc_inq_varid (ncid, "wvl", &id_wvl);
  CHKERR (status);
  wvl    = calloc (nwvl, sizeof (double));
  status = nc_get_var_double (ncid, id_wvl, wvl);
  CHKERR (status);

  status = nc_inq_varid (ncid, "nwvl_in_band", &id_nwvl_in_band);
  CHKERR (status);
  output->wl.nlambda_in_reptran_band = calloc (nbands, sizeof (int));
  status                             = nc_get_var_int (ncid, id_nwvl_in_band, (output->wl.nlambda_in_reptran_band));
  CHKERR (status);

  status = nc_inq_varid (ncid, "iwvl", &id_iwvl);
  CHKERR (status);
  status = nc_inq_varid (ncid, "iwvl_weight", &id_weight_wvl);
  CHKERR (status);
  status = ASCII_calloc_int (&(output->wl.reptran_band), nbands, max_nwvl_in_band);
  CHKERR (status);
  status = ASCII_calloc_double (&(output->wl.weight_reptran_band), nbands, max_nwvl_in_band);
  CHKERR (status);
  i_wvl_tmp      = calloc (nbands * max_nwvl_in_band, sizeof (int));
  weight_wvl_tmp = calloc (nbands * max_nwvl_in_band, sizeof (double));
  status         = nc_get_var_int (ncid, id_iwvl, i_wvl_tmp);
  CHKERR (status);
  status = nc_get_var_double (ncid, id_weight_wvl, weight_wvl_tmp);
  CHKERR (status);
  for (i_band = 0; i_band < nbands; i_band++) {
    for (i_wvl = 0; i_wvl < max_nwvl_in_band; i_wvl++) {
      output->wl.reptran_band[i_band][i_wvl]        = i_wvl_tmp[(i_band + i_wvl * nbands)];
      output->wl.weight_reptran_band[i_band][i_wvl] = weight_wvl_tmp[(i_band + i_wvl * nbands)];
    }
  }
  free (i_wvl_tmp);
  free (weight_wvl_tmp);

  status = nc_inq_varid (ncid, "extra", &id_extra);
  if (status != 0) {
    no_extra = 1; /* representative wavelengths for thermal calculations */
    if (input.source != SRC_THERMAL) {
      fprintf (stderr, "Error: Representative wavelengths file %s was created for\n", cdf_filename);
      fprintf (stderr, "  thermal source, but in the input file 'source' is not set to 'thermal'.");
      CHKERR (-4);
    }
  } else {
    no_extra = 0; /* representative wavelengths for solar calculations  */
    if (input.source != SRC_SOLAR) {
      fprintf (stderr, "Error: Representative wavelengths file %s was created for\n", cdf_filename);
      fprintf (stderr, "  solar source, but in the input file 'source' is not set to 'solar'.");
      CHKERR (-5);
    }
  }
  extra = calloc (nwvl, sizeof (double));

  if (no_extra) {
    status = 0;
    for (i_wvl = 0; i_wvl < nwvl; i_wvl++)
      extra[i_wvl] = 1; /* In case of thermal calculations extra is set to 1 */
  } else
    status = nc_get_var_double (ncid, id_extra, extra);

  status += nc_inq_varid (ncid, "wvl_integral", &id_wvl_integral);
  output->wl.width_of_reptran_band = calloc (nbands, sizeof (double));
  status += nc_get_var_double (ncid, id_wvl_integral, output->wl.width_of_reptran_band);

  if (status != 0) {
    fprintf (stderr, "Error %d reading %s\n", status, cdf_filename);
    CHKERR (status);
  }

  nc_close (ncid);

  /* find required ('active') representative wavelengths */
  wvl_active = calloc (nwvl, sizeof (int));
  for (i = 0; i < nwvl; i++)
    wvl_active[i] = 0;

  /* go through transmission wavelength grid */
  for (i_wvl = 0; i_wvl < output->wl.nlambda_t; i_wvl++)
    for (i = 0; i < output->wl.nlambda_in_reptran_band[output->wl.reptran_band_t[i_wvl]]; i++)
      wvl_active[output->wl.reptran_band[output->wl.reptran_band_t[i_wvl]][i] - 1] = 1;

  /* count active wavelengths */
  n_active = 0;
  for (i_wvl = 0; i_wvl < nwvl; i_wvl++)
    if (wvl_active[i_wvl] == 1)
      n_active++;

  /* find out wavelength for tau scaling*/
  tau_wvl = -1;

  if (input.aer.modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET] >= 0) {
    if (input.aer.tau_wvl_lambda > 0)
      tau_wvl = input.aer.tau_wvl_lambda;
    else
      tau_wvl = 550;
  }

  for (isp = 0; isp < input.n_caoth; isp++) {
    if (input.caoth[isp].modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET] >= 0)
      tau_wvl = 550;
  }

  /* find out if tau_wvl needs to be added to wavelength grid */
  add_tau_wvl_to_lambda_r = 0;
  if (tau_wvl > 0) {

    /* find out wavelength range */
    wvl_min_tmp = 1e12;
    wvl_max_tmp = 0;
    for (i_wvl = 0; i_wvl < nwvl; i_wvl++) {
      if (wvl_active[i_wvl] == 1) {
        if (wvl_min_tmp > wvl[i_wvl])
          wvl_min_tmp = wvl[i_wvl];
        if (wvl_max_tmp < wvl[i_wvl])
          wvl_max_tmp = wvl[i_wvl];
      }
    }

    /* check if tau_wvl is outside the wavelength range */
    if (wvl_max_tmp < tau_wvl || wvl_min_tmp > tau_wvl) {
      add_tau_wvl_to_lambda_r = 1;
      n_active += 1;
    }
  }

  /* allocate required arrays */
  output->wl.nlambda_r              = n_active;
  output->wl.lambda_r               = calloc (output->wl.nlambda_r, sizeof (float));
  output->wl.iwvl_in_reptran_file_r = calloc (output->wl.nlambda_r, sizeof (int));
  output->wl.extra_reptran_r        = calloc (output->wl.nlambda_r, sizeof (float));
  output->wl.wvnmlo_r               = calloc (output->wl.nlambda_r, sizeof (float));
  output->wl.wvnmhi_r               = calloc (output->wl.nlambda_r, sizeof (float));
  index_in_lambda_r                 = calloc (nwvl, sizeof (int));

  /* reset counter */
  n_active = 0;

  /* if tau_wvl needs to be added before the representative wavelengths */
  if (add_tau_wvl_to_lambda_r == 1 && wvl_min_tmp > tau_wvl) {
    output->wl.lambda_r[n_active]               = tau_wvl;
    output->wl.iwvl_in_reptran_file_r[n_active] = -1;
    output->wl.extra_reptran_r[n_active]        = 0;
    n_active++;
    if (input.verbose)
      fprintf (stderr, "       Added the hidden wavelength %f nm at which the aerosol/cloud amount is scaled.\n", tau_wvl);
  }

  /* add the required representative wavelengths */
  for (i_wvl = 0; i_wvl < nwvl; i_wvl++) {
    if (wvl_active[i_wvl] == 1) {
      output->wl.lambda_r[n_active]               = wvl[i_wvl];
      output->wl.iwvl_in_reptran_file_r[n_active] = i_wvl;
      output->wl.extra_reptran_r[n_active]        = extra[i_wvl];
      index_in_lambda_r[i_wvl]                    = n_active;
      n_active++;
    } else
      index_in_lambda_r[i_wvl] = -1;
  }

  /* if tau_wvl needs to be added after the representative wavelengths */
  if (add_tau_wvl_to_lambda_r == 1 && wvl_max_tmp < tau_wvl) {
    output->wl.lambda_r[n_active]               = tau_wvl;
    output->wl.iwvl_in_reptran_file_r[n_active] = -1;
    output->wl.extra_reptran_r[n_active]        = 0;
    n_active++;
    if (input.verbose)
      fprintf (stderr, "       Added the hidden wavelength %f nm at which the aerosol amount is scaled.\n", tau_wvl);
  }

  /* calculate wavenumbers according to given bandwidth */
  if (output->bandwidth_unit == UNIT_PER_CM_1) {
    for (i_wvl = 0; i_wvl < n_active; i_wvl++) {
      output->wl.wvnmlo_r[i_wvl] = 1.0E7 / output->wl.lambda_r[i_wvl] - input.bandwidth / 2.0;
      output->wl.wvnmhi_r[i_wvl] = output->wl.wvnmlo_r[i_wvl] + input.bandwidth;
    }
  } else if (output->bandwidth_unit == UNIT_PER_NM) {
    for (i_wvl = 0; i_wvl < n_active; i_wvl++) {
      output->wl.wvnmlo_r[i_wvl] = 1.0E7 / (output->wl.lambda_r[i_wvl] + input.bandwidth / 2.0);
      output->wl.wvnmhi_r[i_wvl] = 1.0E7 / (output->wl.lambda_r[i_wvl] - input.bandwidth / 2.0);
    }
  } else {
    fprintf (stderr, "Error, unsupported bandwidth_unit %d in while setting up rte_wavelength_grid.\n", output->bandwidth_unit);
    CHKERR (-1);
  }

  /* change reptran_band from pointing to wvl (netcdf variable) to pointing to lambda_r (internal grid) */
  for (i_band = 0; i_band < nbands; i_band++)
    for (i = 0; i < output->wl.nlambda_in_reptran_band[i_band]; i++)
      if (output->wl.reptran_band[i_band][i] > -1)
        output->wl.reptran_band[i_band][i] = index_in_lambda_r[output->wl.reptran_band[i_band][i] - 1];

  if (input.verbose)
    fprintf (stderr, "       %d wavelengths set by set_rte_wl_grid_reptran().\n", output->wl.nlambda_r);

  free (wvl);
  free (wvl_active);
  free (extra);
  free (index_in_lambda_r);

  return 0;
#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
  fprintf (stderr, " ***********************************************************************\n");
  return -1;
#endif
}

/**************************************************************/
/* Convolve the spectrum with a user-defined slit function    */
/**************************************************************/

int convolve (input_struct input, output_struct* output) {
  int     status;
  double* x_data = NULL;

  int     n_slit = 0;
  double *x_slit = NULL, *y_slit = NULL;
  int     nstokes = 0;

  nstokes = input.rte.mc.nstokes;

  /* read slit function from file */
  status = read_2c_file (input.filename[FN_SLITFUNCTION], &x_slit, &y_slit, &n_slit);
  CHKERR (status);

  x_data = calloc (output->wl.nlambda_h, sizeof (double));

  for (int iv = 0; iv < output->wl.nlambda_h; iv++)
    x_data[iv] = (double) output->wl.lambda_h[iv];

  /* convolve albmed and trnmed */

  /* albmed */
  for (int iu = 0; iu < input.rte.numu; iu++) {
    status = cnvlv (x_data, output->albmed[iu], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
    CHKERR (status);
  }

  /* trnmed */
  for (int iu = 0; iu < input.rte.numu; iu++) {
    status = cnvlv (x_data, output->trnmed[iu], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
    CHKERR (status);
  }

  for (int lev = 0; lev < output->atm.nzout; lev++) {

    /* Convolve rfldir, rfldn, flup and uavg for all cases */

    /* rfldir */
    status = cnvlv (x_data, output->rfldir[lev], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
    CHKERR (status);

    /* rfldn */
    status = cnvlv (x_data, output->rfldn[lev], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
    CHKERR (status);

    /* flup */
    status = cnvlv (x_data, output->flup[lev], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
    CHKERR (status);

    /* uavg */
    status = cnvlv (x_data, output->uavg[lev], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
    CHKERR (status);

    /* 3D fields */
    if (input.output3D)
      for (int is = output->islower; is <= output->isupper; is += output->isstep)
        for (int js = output->jslower; js <= output->jsupper; js += output->jsstep) {

          status = cnvlv (x_data, output->rfldir3d[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0); /* rfldir3d */
          CHKERR (status);

          status = cnvlv (x_data, output->rfldn3d[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0); /* rfldn3d */
          CHKERR (status);

          status = cnvlv (x_data, output->flup3d[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0); /* flup3d */
          CHKERR (status);

          status = cnvlv (x_data, output->uavgso3d[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0); /* uavgso3d */
          CHKERR (status);

          status = cnvlv (x_data, output->uavgdn3d[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0); /* uavgdn3d */
          CHKERR (status);

          status = cnvlv (x_data, output->uavgup3d[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0); /* uavgup3d */
          CHKERR (status);

          for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
            status =
              cnvlv (x_data, output->fl3d_is[lev][is][js][ic], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0); /* fl3d_is */
            CHKERR (status);
          }

          for (int ip = 0; ip < nstokes; ip++) {
            for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
              status = cnvlv (x_data,
                              output->radiance3d[lev][is][js][ip][ic],
                              output->wl.nlambda_h,
                              x_slit,
                              y_slit,
                              n_slit,
                              0); /* radiance3d */
              CHKERR (status);
            }
          }

          if (input.rte.mc.backward.absorption) {
            status =
              cnvlv (x_data, output->absback3d[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0); /* absback3d */
            CHKERR (status);
          }

          /* variances */
          if (input.rte.mc.std) {

            status =
              cnvlv (x_data, output->rfldir3d_var[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 1); /* rfldir3d_var */
            CHKERR (status);

            status =
              cnvlv (x_data, output->rfldn3d_var[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 1); /* rfldn3d_var */
            CHKERR (status);

            status =
              cnvlv (x_data, output->flup3d_var[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 1); /* flup3d_var */
            CHKERR (status);

            status =
              cnvlv (x_data, output->uavgso3d_var[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 1); /* uavgso3d_var */
            CHKERR (status);

            status =
              cnvlv (x_data, output->uavgdn3d_var[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 1); /* uavgdn3d_var */
            CHKERR (status);

            status =
              cnvlv (x_data, output->uavgup3d_var[lev][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 1); /* uavgup3d_var */
            CHKERR (status);

            for (int ip = 0; ip < nstokes; ip++) {
              status = cnvlv (x_data,
                              output->radiance3d_var[lev][is][js][ip],
                              output->wl.nlambda_h,
                              x_slit,
                              y_slit,
                              n_slit,
                              1); /* radiance3d_var */
              CHKERR (status);
            }

            if (input.rte.mc.backward.absorption) {
              status = cnvlv (x_data,
                              output->absback3d_var[lev][is][js],
                              output->wl.nlambda_h,
                              x_slit,
                              y_slit,
                              n_slit,
                              1); /* absback3d_var */
              CHKERR (status);
            }
          }
        }

    if (input.rte.solver != SOLVER_FTWOSTR) {

      /* Convolve uavgso, uavgdn, uavgup */
      status = cnvlv (x_data, output->uavgso[lev], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
      CHKERR (status);
      status = cnvlv (x_data, output->uavgdn[lev], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
      CHKERR (status);
      status = cnvlv (x_data, output->uavgup[lev], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
      CHKERR (status);

      /* Convolve u0u */
      for (int iu = 0; iu < input.rte.numu; iu++) {
        status = cnvlv (x_data, output->u0u[lev][iu], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
        CHKERR (status);
      }

      if (output->print_phi > 0) {
        /* Convolve uu */
        for (int j = 0; j < input.rte.nphi; j++) {
          for (int iu = 0; iu < input.rte.numu; iu++) {
            status = cnvlv (x_data, output->uu[lev][j][iu], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
            CHKERR (status);
          }
        }
      }
    }
  }

  /* 3D absorption */
  if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
      input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
    const size_t N3d_layers = vec_count_true (output->atm.Nzcld, output->atm.threed);
    for (size_t k = 0; k < N3d_layers; k++) {
      for (int is = 0; is < output->atm.Nxcld; is++) {
        for (int js = 0; js < output->atm.Nycld; js++) {

          status = cnvlv (x_data, output->abs3d[k][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
          CHKERR (status);

          if (input.rte.mc.std) {
            status = cnvlv (x_data, output->abs3d_var[k][is][js], output->wl.nlambda_h, x_slit, y_slit, n_slit, 0);
            CHKERR (status);
          }
        }
      }
    }
  }

  free (x_slit);
  free (y_slit);
  free (x_data);

  return 0;
}

/**************************************************************/
/* Convolve a spectrum with a slit function.                  */
/* Internal function, used by convolve().                     */
/**************************************************************/

static int cnvlv (double* x_spec, float* y_spec, int n_spec, double* x_slit, double* y_slit, int n_slit, int std) {
  int     iv = 0, status = 0;
  double* y_data = calloc (n_spec, sizeof (double));
  double* tmp    = NULL;

  for (iv = 0; iv < n_spec; iv++)
    y_data[iv] = (double) y_spec[iv];

  status = int_convolute (x_spec, y_data, n_spec, x_slit, y_slit, n_slit, &tmp, std);
  CHKERR (status);

  for (iv = 0; iv < n_spec; iv++)
    y_spec[iv] = (float) tmp[iv];

  free (y_data);
  free (tmp);

  return 0;
}

/**************************************************************/
/* Interpolate the high resolution spectrum  *_h              */
/* to a user-defined wavelength grid.        *_s              */
/**************************************************************/

int spline_interpolate (input_struct input, output_struct* output) {
  int     status = 0;
  double* x_user = NULL;
  int     linear = 0;

  if (strlen (input.filename[FN_SPLINE]) > 0) {

    /* read file with user x values */
    status = read_1c_file (input.filename[FN_SPLINE], &x_user, &output->wl.nlambda_s);
    CHKERR (status);

    output->wl.lambda_s = calloc (output->wl.nlambda_s, sizeof (float));

    for (int iv = 0; iv < output->wl.nlambda_s; iv++)
      output->wl.lambda_s[iv] = (float) x_user[iv];

    /* Check wavelength ranges */
    status = 0;
    if (output->wl.start > output->wl.lambda_s[0]) {
      fprintf (stderr, "Error, spline wavelength %f is smaller than wvn %f\n", output->wl.lambda_s[0], output->wl.start);
      status--;
    }

    if (output->wl.end < output->wl.lambda_s[output->wl.nlambda_s - 1]) {
      fprintf (stderr,
               "Error, spline wavelength %f is greater than wvn %f\n",
               output->wl.lambda_s[output->wl.nlambda_s - 1],
               output->wl.end);
      status--;
    }
    CHKERR (status);

  } else {
    output->wl.nlambda_s = (int) ((input.spline_lambda_1 - input.spline_lambda_0) / input.spline_lambda_step + 1);

    output->wl.lambda_s = calloc (output->wl.nlambda_s, sizeof (float));

    for (int iv = 0; iv < output->wl.nlambda_s; iv++)
      output->wl.lambda_s[iv] = input.spline_lambda_0 + (float) iv * input.spline_lambda_step;
  }

  /* Allocation of the solar zenith angle */
  output->sza_s = calloc (output->wl.nlambda_s, sizeof (float));
  CHKPOINTEROUT (output->sza_s, "failed allocation of sza_s in spline_interpolate");

  /* Interpolation of the solar zenith angle to output wavelength grid */
  status += arb_wvn (output->wl.nlambda_h,
                     output->wl.lambda_h,
                     output->sza_h,
                     output->wl.nlambda_s,
                     output->wl.lambda_s,
                     output->sza_s,
                     linear,
                     0);
  CHKERR (status);

  /* copy to sza_h for print_output */
  free (output->sza_h);
  output->sza_h = calloc (output->wl.nlambda_s, sizeof (float));
  for (int iv = 0; iv < output->wl.nlambda_s; iv++)
    output->sza_h[iv] = (float) output->sza_s[iv];

  /* Interpolation to output wavelength grid */

  /* albmed */
  for (int iu = 0; iu < input.rte.numu; iu++) {
    status = arb_wvn (output->wl.nlambda_h,
                      output->wl.lambda_h,
                      output->albmed[iu],
                      output->wl.nlambda_s,
                      output->wl.lambda_s,
                      output->albmed[iu],
                      linear,
                      0);
    CHKERR (status);
  }

  /* trnmed */
  for (int iu = 0; iu < input.rte.numu; iu++) {
    status = arb_wvn (output->wl.nlambda_h,
                      output->wl.lambda_h,
                      output->trnmed[iu],
                      output->wl.nlambda_s,
                      output->wl.lambda_s,
                      output->trnmed[iu],
                      linear,
                      0);
    CHKERR (status);
  }

  for (int lev = 0; lev < output->atm.nzout; lev++) {

    /* Spline interpolate rfldir, rfldn, flup, uavg and heat for all cases */

    /* rfldir */
    status = arb_wvn (output->wl.nlambda_h,
                      output->wl.lambda_h,
                      output->rfldir[lev],
                      output->wl.nlambda_s,
                      output->wl.lambda_s,
                      output->rfldir[lev],
                      linear,
                      0);
    CHKERR (status);

    /* rfldn */
    status = arb_wvn (output->wl.nlambda_h,
                      output->wl.lambda_h,
                      output->rfldn[lev],
                      output->wl.nlambda_s,
                      output->wl.lambda_s,
                      output->rfldn[lev],
                      linear,
                      0);
    CHKERR (status);

    /* flup */
    status = arb_wvn (output->wl.nlambda_h,
                      output->wl.lambda_h,
                      output->flup[lev],
                      output->wl.nlambda_s,
                      output->wl.lambda_s,
                      output->flup[lev],
                      linear,
                      0);
    CHKERR (status);

    /* uavg */
    status = arb_wvn (output->wl.nlambda_h,
                      output->wl.lambda_h,
                      output->uavg[lev],
                      output->wl.nlambda_s,
                      output->wl.lambda_s,
                      output->uavg[lev],
                      linear,
                      0);
    CHKERR (status);

    /* heating rates */
    if (input.heating != HEAT_NONE) {
      status = arb_wvn (output->wl.nlambda_h,
                        output->wl.lambda_h,
                        output->heat[lev],
                        output->wl.nlambda_s,
                        output->wl.lambda_s,
                        output->heat[lev],
                        linear,
                        0);
      CHKERR (status);

      status = arb_wvn (output->wl.nlambda_h,
                        output->wl.lambda_h,
                        output->emis[lev],
                        output->wl.nlambda_s,
                        output->wl.lambda_s,
                        output->emis[lev],
                        linear,
                        0);
      CHKERR (status);

      status = arb_wvn (output->wl.nlambda_h,
                        output->wl.lambda_h,
                        output->w_zout[lev],
                        output->wl.nlambda_s,
                        output->wl.lambda_s,
                        output->w_zout[lev],
                        linear,
                        0);
      CHKERR (status);
    }

    /* 3D fields */
    if (input.output3D)
      for (int is = output->islower; is <= output->isupper; is += output->isstep)
        for (int js = output->jslower; js <= output->jsupper; js += output->jsstep) {

          /* rfldir3d */
          status = arb_wvn (output->wl.nlambda_h,
                            output->wl.lambda_h,
                            output->rfldir3d[lev][is][js],
                            output->wl.nlambda_s,
                            output->wl.lambda_s,
                            output->rfldir3d[lev][is][js],
                            linear,
                            0);
          CHKERR (status);

          /* rfldn3d */
          status = arb_wvn (output->wl.nlambda_h,
                            output->wl.lambda_h,
                            output->rfldn3d[lev][is][js],
                            output->wl.nlambda_s,
                            output->wl.lambda_s,
                            output->rfldn3d[lev][is][js],
                            linear,
                            0);
          CHKERR (status);

          /* rflup3d */
          status = arb_wvn (output->wl.nlambda_h,
                            output->wl.lambda_h,
                            output->flup3d[lev][is][js],
                            output->wl.nlambda_s,
                            output->wl.lambda_s,
                            output->flup3d[lev][is][js],
                            linear,
                            0);
          CHKERR (status);

          /* uavgso3d */
          status = arb_wvn (output->wl.nlambda_h,
                            output->wl.lambda_h,
                            output->uavgso3d[lev][is][js],
                            output->wl.nlambda_s,
                            output->wl.lambda_s,
                            output->uavgso3d[lev][is][js],
                            linear,
                            0);
          CHKERR (status);

          /* uavgdn3d */
          status = arb_wvn (output->wl.nlambda_h,
                            output->wl.lambda_h,
                            output->uavgdn3d[lev][is][js],
                            output->wl.nlambda_s,
                            output->wl.lambda_s,
                            output->uavgdn3d[lev][is][js],
                            linear,
                            0);
          CHKERR (status);

          /* uavgup3d */
          status = arb_wvn (output->wl.nlambda_h,
                            output->wl.lambda_h,
                            output->uavgup3d[lev][is][js],
                            output->wl.nlambda_s,
                            output->wl.lambda_s,
                            output->uavgup3d[lev][is][js],
                            linear,
                            0);
          CHKERR (status);

          if (input.rte.mc.concentration_is || input.rte.mc.spectral_is) {
            for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
              status = arb_wvn (output->wl.nlambda_h,
                                output->wl.lambda_h,
                                output->fl3d_is[lev][is][js][ic],
                                output->wl.nlambda_s,
                                output->wl.lambda_s,
                                output->fl3d_is[lev][is][js][ic],
                                linear,
                                0);
              CHKERR (status);
            }
          }

          for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
            for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
              /* radiance3d */
              status = arb_wvn (output->wl.nlambda_h,
                                output->wl.lambda_h,
                                output->radiance3d[lev][is][js][ip][ic],
                                output->wl.nlambda_s,
                                output->wl.lambda_s,
                                output->radiance3d[lev][is][js][ip][ic],
                                linear,
                                0);
              CHKERR (status);
            }
          }

          /* absback3d */
          if (input.rte.mc.backward.absorption) {
            status = arb_wvn (output->wl.nlambda_h,
                              output->wl.lambda_h,
                              output->absback3d[lev][is][js],
                              output->wl.nlambda_s,
                              output->wl.lambda_s,
                              output->absback3d[lev][is][js],
                              linear,
                              0);
            CHKERR (status);
          }

          /* variances */
          if (input.rte.mc.std) {

            /* rfldir3d_var */
            status = arb_wvn (output->wl.nlambda_h,
                              output->wl.lambda_h,
                              output->rfldir3d_var[lev][is][js],
                              output->wl.nlambda_s,
                              output->wl.lambda_s,
                              output->rfldir3d_var[lev][is][js],
                              linear,
                              0);
            CHKERR (status);

            /* rfldn3d_var */
            status = arb_wvn (output->wl.nlambda_h,
                              output->wl.lambda_h,
                              output->rfldn3d_var[lev][is][js],
                              output->wl.nlambda_s,
                              output->wl.lambda_s,
                              output->rfldn3d_var[lev][is][js],
                              linear,
                              0);
            CHKERR (status);

            /* rflup3d_var */
            status = arb_wvn (output->wl.nlambda_h,
                              output->wl.lambda_h,
                              output->flup3d_var[lev][is][js],
                              output->wl.nlambda_s,
                              output->wl.lambda_s,
                              output->flup3d_var[lev][is][js],
                              linear,
                              0);
            CHKERR (status);

            /* uavgso3d_var */
            status = arb_wvn (output->wl.nlambda_h,
                              output->wl.lambda_h,
                              output->uavgso3d_var[lev][is][js],
                              output->wl.nlambda_s,
                              output->wl.lambda_s,
                              output->uavgso3d_var[lev][is][js],
                              linear,
                              0);
            CHKERR (status);

            /* uavgdn3d_var */
            status = arb_wvn (output->wl.nlambda_h,
                              output->wl.lambda_h,
                              output->uavgdn3d_var[lev][is][js],
                              output->wl.nlambda_s,
                              output->wl.lambda_s,
                              output->uavgdn3d_var[lev][is][js],
                              linear,
                              0);
            CHKERR (status);

            /* uavgup3d_var */
            status = arb_wvn (output->wl.nlambda_h,
                              output->wl.lambda_h,
                              output->uavgup3d_var[lev][is][js],
                              output->wl.nlambda_s,
                              output->wl.lambda_s,
                              output->uavgup3d_var[lev][is][js],
                              linear,
                              0);
            CHKERR (status);

            for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {

              /* radiance3d_var */
              status = arb_wvn (output->wl.nlambda_h,
                                output->wl.lambda_h,
                                output->radiance3d_var[lev][is][js][ip],
                                output->wl.nlambda_s,
                                output->wl.lambda_s,
                                output->radiance3d_var[lev][is][js][ip],
                                linear,
                                0);
              CHKERR (status);
            }

            /* absback3d_var */
            if (input.rte.mc.backward.absorption) {
              status = arb_wvn (output->wl.nlambda_h,
                                output->wl.lambda_h,
                                output->absback3d_var[lev][is][js],
                                output->wl.nlambda_s,
                                output->wl.lambda_s,
                                output->absback3d_var[lev][is][js],
                                linear,
                                0);
              CHKERR (status);
            }
          }
        }

    if (input.rte.solver != SOLVER_FTWOSTR) {

      /* Spline interpolate uavgso, uavgdn, uavgup */

      /* uavgso */
      status = arb_wvn (output->wl.nlambda_h,
                        output->wl.lambda_h,
                        output->uavgso[lev],
                        output->wl.nlambda_s,
                        output->wl.lambda_s,
                        output->uavgso[lev],
                        linear,
                        0);
      CHKERR (status);

      /* uavgdn */
      status = arb_wvn (output->wl.nlambda_h,
                        output->wl.lambda_h,
                        output->uavgdn[lev],
                        output->wl.nlambda_s,
                        output->wl.lambda_s,
                        output->uavgdn[lev],
                        linear,
                        0);
      CHKERR (status);

      /* uavgup */
      status = arb_wvn (output->wl.nlambda_h,
                        output->wl.lambda_h,
                        output->uavgup[lev],
                        output->wl.nlambda_s,
                        output->wl.lambda_s,
                        output->uavgup[lev],
                        linear,
                        0);
      CHKERR (status);

      /* Spline interpolate u0u */
      for (int iu = 0; iu < input.rte.numu; iu++) {
        status = arb_wvn (output->wl.nlambda_h,
                          output->wl.lambda_h,
                          output->u0u[lev][iu],
                          output->wl.nlambda_s,
                          output->wl.lambda_s,
                          output->u0u[lev][iu],
                          linear,
                          0);
        CHKERR (status);
      }

      if (output->print_phi > 0) {

        /* Spline interpolate uu */
        for (int j = 0; j < input.rte.nphi; j++) {
          for (int iu = 0; iu < input.rte.numu; iu++) {
            status = arb_wvn (output->wl.nlambda_h,
                              output->wl.lambda_h,
                              output->uu[lev][j][iu],
                              output->wl.nlambda_s,
                              output->wl.lambda_s,
                              output->uu[lev][j][iu],
                              linear,
                              0);
            CHKERR (status);
          }
        }
      }
    }
  }

  if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
      input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
    const size_t N3d_layers = vec_count_true (output->atm.Nzcld, output->atm.threed);
    for (size_t k = 0; k < N3d_layers; k++) {
      for (int is = 0; is < output->atm.Nxcld; is++) {
        for (int js = 0; js < output->atm.Nycld; js++) {

          status = arb_wvn (output->wl.nlambda_h,
                            output->wl.lambda_h,
                            output->abs3d[k][is][js],
                            output->wl.nlambda_s,
                            output->wl.lambda_s,
                            output->abs3d[k][is][js],
                            linear,
                            0);
          CHKERR (status);

          if (input.rte.mc.std) {
            status = arb_wvn (output->wl.nlambda_h,
                              output->wl.lambda_h,
                              output->abs3d_var[k][is][js],
                              output->wl.nlambda_s,
                              output->wl.lambda_s,
                              output->abs3d_var[k][is][js],
                              linear,
                              0);
            CHKERR (status);
          }
        }
      }
    }
  }

  free (output->wl.lambda_h);

  output->wl.nlambda_h = output->wl.nlambda_s;

  output->wl.lambda_h = calloc (output->wl.nlambda_s, sizeof (float));

  for (int iv = 0; iv < output->wl.nlambda_s; iv++)
    output->wl.lambda_h[iv] = (float) output->wl.lambda_s[iv];

  return status;
}

/****************************************************************************************************************/
/* Transfer radiative transfer results from radiative transfer wavelength grid to transmittance wavelength grid */
/****************************************************************************************************************/
int internal_to_transmittance_grid (input_struct input, output_struct* output) {
  int status;

  double ffactor = 0.0, rfactor = 0.0, hfactor = 0.0;
  double unit_factor;

  /* Pointers to results are just copied, if no representative wavelengths are used */
  if (output->wl.use_reptran == 0) {
    CHKERROUT (output->wl.nlambda_t - output->wl.nlambda_r,
               "if we dont copy values here, shouldnt _t and _r grid be the same size?");

    output->albmed_t = output->albmed_r;
    output->trnmed_t = output->trnmed_r;

    output->flup_t            = output->flup_r;
    output->rfldir_t          = output->rfldir_r;
    output->rfldn_t           = output->rfldn_r;
    output->uavg_t            = output->uavg_r;
    output->uavgdn_t          = output->uavgdn_r;
    output->uavgso_t          = output->uavgso_r;
    output->uavgup_t          = output->uavgup_r;
    output->sslidar_nphot_t   = output->sslidar_nphot_r;
    output->sslidar_nphot_q_t = output->sslidar_nphot_q_r;
    output->sslidar_ratio_t   = output->sslidar_ratio_r;
    if (input.heating != HEAT_NONE) {
      output->heat_t   = output->heat_r;
      output->emis_t   = output->emis_r;
      output->w_zout_t = output->w_zout_r;
    }

    /* radiances */
    output->u0u_t = output->u0u_r;
    output->uu_t  = output->uu_r;

    /* 3D fields */
    if (input.output3D) {
      output->rfldir3d_t   = output->rfldir3d_r;
      output->rfldn3d_t    = output->rfldn3d_r;
      output->flup3d_t     = output->flup3d_r;
      output->uavgso3d_t   = output->uavgso3d_r;
      output->uavgdn3d_t   = output->uavgdn3d_r;
      output->uavgup3d_t   = output->uavgup3d_r;
      output->radiance3d_t = output->radiance3d_r;

      if (input.rte.mc.concentration_is || input.rte.mc.spectral_is)
        output->fl3d_is_t = output->fl3d_is_r;

      if (input.rte.mc.jacobian[DIM_1D])
        output->jacobian_t = output->jacobian_r;

      if (input.rte.mc.backward.absorption)
        output->absback3d_t = output->absback3d_r;

      /* variances */
      if (input.rte.mc.std) {

        output->rfldir3d_var_t   = output->rfldir3d_var_r;
        output->rfldn3d_var_t    = output->rfldn3d_var_r;
        output->flup3d_var_t     = output->flup3d_var_r;
        output->uavgso3d_var_t   = output->uavgso3d_var_r;
        output->uavgdn3d_var_t   = output->uavgdn3d_var_r;
        output->uavgup3d_var_t   = output->uavgup3d_var_r;
        output->radiance3d_var_t = output->radiance3d_var_r;
        if (input.rte.mc.backward.absorption)
          output->absback3d_var_t = output->absback3d_var_r;
      }
    }

    /* absorption is defined on the 3D caoth grid, not on the user-defined grid */
    if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
        input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
      output->abs3d_t = output->abs3d_r;
      if (input.rte.mc.std) {
        output->abs3d_var_t = output->abs3d_var_r;
      }
    }

    /* copy solar zenith angle */
    output->atm.sza_t = output->atm.sza_r;

    output->triangle_results_t = output->triangle_results_r;
  } else { /* if representative wavelengths are used, the results at the radiative transfer wavelength grid need to be weighted  */

    output->atm.sza_t = calloc (output->wl.nlambda_t, sizeof (float));
    status            = weighting_reptran (&output->wl, 0, output->atm.sza_r, output->atm.sza_t, input.ck_reptran_option);
    CHKERR (status);

    output->albmed_t = calloc (input.rte.numu, sizeof (float*));
    output->trnmed_t = calloc (input.rte.numu, sizeof (float*));
    for (int iu = 0; iu < input.rte.numu; iu++) {

      output->albmed_t[iu] = calloc (output->wl.nlambda_t, sizeof (float));
      output->trnmed_t[iu] = calloc (output->wl.nlambda_t, sizeof (float));

      status = weighting_reptran (&output->wl, 0, output->albmed_r[iu], output->albmed_t[iu], input.ck_reptran_option);
      CHKERR (status);
      status = weighting_reptran (&output->wl, 0, output->trnmed_r[iu], output->trnmed_t[iu], input.ck_reptran_option);
      CHKERR (status);
    }

    if (input.heating != HEAT_NONE) {
      output->heat_t   = calloc (output->atm.nzout, sizeof (float*));
      output->emis_t   = calloc (output->atm.nzout, sizeof (float*));
      output->w_zout_t = calloc (output->atm.nzout, sizeof (float*));
    }

    output->flup_t            = calloc (output->atm.nzout, sizeof (float*));
    output->rfldir_t          = calloc (output->atm.nzout, sizeof (float*));
    output->rfldn_t           = calloc (output->atm.nzout, sizeof (float*));
    output->uavg_t            = calloc (output->atm.nzout, sizeof (float*));
    output->uavgdn_t          = calloc (output->atm.nzout, sizeof (float*));
    output->uavgso_t          = calloc (output->atm.nzout, sizeof (float*));
    output->uavgup_t          = calloc (output->atm.nzout, sizeof (float*));
    output->sslidar_nphot_t   = calloc (output->atm.nzout, sizeof (float*));
    output->sslidar_nphot_q_t = calloc (output->atm.nzout, sizeof (float*));
    output->sslidar_ratio_t   = calloc (output->atm.nzout, sizeof (float*));

    if (input.output3D) { /** Allocate memory for 3D results output->*_t */
      const size_t Nz      = output->atm.nzout;
      const size_t Nx      = output->isupper - output->islower + 1;
      const size_t Ny      = output->jsupper - output->jslower + 1;
      const size_t Nlambda = output->wl.nlambda_t;
      const size_t Nalis   = output->mc.alis.Nc;
      const size_t Nstokes = input.rte.mc.nstokes;

      {  // alloc 4d vars like flxs ans absorption
        const size_t dims[4]    = {Nz, Nx, Ny, Nlambda};
        const int    offsets[4] = {0, -output->islower, -output->jslower, 0};

        if (!input.quiet)
          fprintf (stderr,
                   " ... allocating Nz %lu x Nx %lu x Ny %lu x Nwvl %lu = %lu entries (%f Gbytes) for 3D fluxes on "
                   "transmittance grid\n",
                   dims[0],
                   dims[1],
                   dims[2],
                   dims[3],
                   dims[0] * dims[1] * dims[2] * dims[3],
                   dims[0] * dims[1] * dims[2] * dims[3] * sizeof (float) / 1024. / 1024. / 1024.);

        status = alloc_array_with_offsets_4d (dims, offsets, "rfldir3d_t", &output->rfldir3d_t);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "rfldn3d_t", &output->rfldn3d_t);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "flup3d_t", &output->flup3d_t);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "uavgso3d_t", &output->uavgso3d_t);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "uavgdn3d_t", &output->uavgdn3d_t);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "uavgup3d_t", &output->uavgup3d_t);
        CHKERR (status);

        if (input.rte.mc.backward.absorption) {
          status = alloc_array_with_offsets_4d (dims, offsets, "absback3d_t", &output->absback3d_t);
          CHKERR (status);
        }

        if (input.rte.mc.std) {
          status = alloc_array_with_offsets_4d (dims, offsets, "rfldir3d_var_t", &output->rfldir3d_var_t);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "rfldn3d_var_t", &output->rfldn3d_var_t);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "flup3d_var_t", &output->flup3d_var_t);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "uavgso3d_var_t", &output->uavgso3d_var_t);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "uavgdn3d_var_t", &output->uavgdn3d_var_t);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "uavgup3d_var_t", &output->uavgup3d_var_t);
          CHKERR (status);

          if (input.rte.mc.backward.absorption) {
            status = alloc_array_with_offsets_4d (dims, offsets, "absback3d_var_t", &output->absback3d_var_t);
            CHKERR (status);
          }
        }
      }

      /* 3d absorption */
      if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
          input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
        const size_t N3d_layers = vec_count_true (output->atm.Nzcld, output->atm.threed);
        const size_t dims[4]    = {N3d_layers, output->atm.Nxcld, output->atm.Nycld, Nlambda};
        const int    offsets[4] = {0, 0, 0, 0};
        if (!input.quiet)
          fprintf (stderr,
                   " ... allocating Nz3D %lu x Nx %lu x Ny %lu x Nwvl %lu = %lu entries (%f Gbytes) for 3D absorption\n",
                   dims[0],
                   dims[1],
                   dims[2],
                   dims[3],
                   dims[0] * dims[1] * dims[2] * dims[3],
                   dims[0] * dims[1] * dims[2] * dims[3] * sizeof (float) / 1024. / 1024. / 1024.);

        status = alloc_array_with_offsets_4d (dims, offsets, "abs3d_t", &output->abs3d_t);
        CHKERR (status);

        if (input.rte.mc.std) {
          status = alloc_array_with_offsets_4d (dims, offsets, "abs3d_var_t", &output->abs3d_var_t);
          CHKERR (status);
        }
      }

      if (output->radiance3d_r) {
        const size_t dims[6]    = {Nz, Nx, Ny, Nstokes, Nalis, Nlambda};
        const int    offsets[6] = {0, -output->islower, -output->jslower, 0, 0, -output->wl.iv_lower};

        status = alloc_array_with_offsets_6d (dims, offsets, "radiance3d_t", &output->radiance3d_t);
        CHKERR (status);
      }

      if (output->radiance3d_var_r) {
        if (input.rte.mc.std) {
          const size_t dims[5]    = {Nz, Nx, Ny, Nstokes, Nlambda};
          const int    offsets[5] = {0, -output->islower, -output->jslower, 0, -output->wl.iv_lower};

          status = alloc_array_with_offsets_5d (dims, offsets, "radiance3d_var_t", &output->radiance3d_var_t);
          CHKERR (status);
        }
      }

      if (input.rte.mc.spectral_is || input.rte.mc.concentration_is) {
        const size_t dims[5]    = {Nz, Nx, Ny, Nalis, Nlambda};
        const int    offsets[5] = {0, -output->islower, -output->jslower, 0, -output->wl.iv_lower};

        status = alloc_array_with_offsets_5d (dims, offsets, "fl3d_is_t", &output->fl3d_is_t);
        CHKERR (status);
      }

      /* So far we allow only 1D output for postprocessing */
      if (input.rte.mc.jacobian[DIM_1D]) {
        const size_t dims[7]    = {Nz, 1, 1, input.n_caoth + 2, 2, output->atm.nlev - 1, Nlambda};
        const int    offsets[7] = {0, 0, 0, 0, -output->wl.iv_lower};

        status = alloc_array_with_offsets_7d (dims, offsets, "fl3d_is_t", &output->jacobian_t);
        CHKERR (status);
      }
    } /** End allocate memory for 3D results output->*_t */

    output->u0u_t = calloc (output->atm.nzout, sizeof (float**));
    output->uu_t  = calloc (output->atm.nzout, sizeof (float***));

    for (int lu = 0; lu < output->atm.nzout; lu++) {

      if (input.heating != HEAT_NONE) {

        output->heat_t[lu]   = calloc (output->wl.nlambda_t, sizeof (float));
        output->emis_t[lu]   = calloc (output->wl.nlambda_t, sizeof (float));
        output->w_zout_t[lu] = calloc (output->wl.nlambda_t, sizeof (float));
        status               = weighting_reptran (&output->wl, 0, output->heat_r[lu], output->heat_t[lu], input.ck_reptran_option);
        CHKERR (status);
        status = weighting_reptran (&output->wl, 0, output->emis_r[lu], output->emis_t[lu], input.ck_reptran_option);
        CHKERR (status);
        status = weighting_reptran (&output->wl, 0, output->w_zout_r[lu], output->w_zout_t[lu], input.ck_reptran_option);
        CHKERR (status);
      }

      output->flup_t[lu]            = calloc (output->wl.nlambda_t, sizeof (float));
      output->rfldir_t[lu]          = calloc (output->wl.nlambda_t, sizeof (float));
      output->rfldn_t[lu]           = calloc (output->wl.nlambda_t, sizeof (float));
      output->uavg_t[lu]            = calloc (output->wl.nlambda_t, sizeof (float));
      output->uavgdn_t[lu]          = calloc (output->wl.nlambda_t, sizeof (float));
      output->uavgso_t[lu]          = calloc (output->wl.nlambda_t, sizeof (float));
      output->uavgup_t[lu]          = calloc (output->wl.nlambda_t, sizeof (float));
      output->sslidar_nphot_t[lu]   = calloc (output->wl.nlambda_t, sizeof (float));
      output->sslidar_nphot_q_t[lu] = calloc (output->wl.nlambda_t, sizeof (float));
      output->sslidar_ratio_t[lu]   = calloc (output->wl.nlambda_t, sizeof (float));
      status = weighting_reptran (&output->wl, 0, output->flup_r[lu], output->flup_t[lu], input.ck_reptran_option);
      CHKERR (status);
      status = weighting_reptran (&output->wl, 0, output->rfldir_r[lu], output->rfldir_t[lu], input.ck_reptran_option);
      CHKERR (status);
      status = weighting_reptran (&output->wl, 0, output->rfldn_r[lu], output->rfldn_t[lu], input.ck_reptran_option);
      CHKERR (status);
      status = weighting_reptran (&output->wl, 0, output->uavg_r[lu], output->uavg_t[lu], input.ck_reptran_option);
      CHKERR (status);
      status = weighting_reptran (&output->wl, 0, output->uavgdn_r[lu], output->uavgdn_t[lu], input.ck_reptran_option);
      CHKERR (status);
      status = weighting_reptran (&output->wl, 0, output->uavgso_r[lu], output->uavgso_t[lu], input.ck_reptran_option);
      CHKERR (status);
      status = weighting_reptran (&output->wl, 0, output->uavgup_r[lu], output->uavgup_t[lu], input.ck_reptran_option);
      CHKERR (status);
      status =
        weighting_reptran (&output->wl, 0, output->sslidar_nphot_r[lu], output->sslidar_nphot_t[lu], input.ck_reptran_option);
      CHKERR (status);
      status =
        weighting_reptran (&output->wl, 0, output->sslidar_nphot_q_r[lu], output->sslidar_nphot_q_t[lu], input.ck_reptran_option);
      CHKERR (status);
      status =
        weighting_reptran (&output->wl, 0, output->sslidar_ratio_r[lu], output->sslidar_ratio_t[lu], input.ck_reptran_option);
      CHKERR (status);

      /* 3D fields */
      if (input.output3D) {
        for (int is = output->islower; is <= output->isupper; is += output->isstep) {
          for (int js = output->jslower; js <= output->jsupper; js += output->jsstep) {
            status = weighting_reptran (&output->wl,
                                        0,
                                        output->rfldir3d_r[lu][is][js],
                                        output->rfldir3d_t[lu][is][js],
                                        input.ck_reptran_option);
            CHKERR (status);
            status = weighting_reptran (&output->wl,
                                        0,
                                        output->rfldn3d_r[lu][is][js],
                                        output->rfldn3d_t[lu][is][js],
                                        input.ck_reptran_option);
            CHKERR (status);
            status = weighting_reptran (&output->wl,
                                        0,
                                        output->flup3d_r[lu][is][js],
                                        output->flup3d_t[lu][is][js],
                                        input.ck_reptran_option);
            CHKERR (status);
            status = weighting_reptran (&output->wl,
                                        0,
                                        output->uavgso3d_r[lu][is][js],
                                        output->uavgso3d_t[lu][is][js],
                                        input.ck_reptran_option);
            CHKERR (status);
            status = weighting_reptran (&output->wl,
                                        0,
                                        output->uavgdn3d_r[lu][is][js],
                                        output->uavgdn3d_t[lu][is][js],
                                        input.ck_reptran_option);
            CHKERR (status);
            status = weighting_reptran (&output->wl,
                                        0,
                                        output->uavgup3d_r[lu][is][js],
                                        output->uavgup3d_t[lu][is][js],
                                        input.ck_reptran_option);
            CHKERR (status);

            if (output->radiance3d_r) {
              if (input.rte.mc.concentration_is || input.rte.mc.spectral_is) {
                for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
                  status = weighting_reptran (&output->wl,
                                              0,
                                              output->fl3d_is_r[lu][is][js][ic],
                                              output->fl3d_is_t[lu][is][js][ic],
                                              input.ck_reptran_option);
                  CHKERR (status);
                }
              }

              for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
                for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
                  status = weighting_reptran (&output->wl,
                                              0,
                                              output->radiance3d_r[lu][is][js][ip][ic],
                                              output->radiance3d_t[lu][is][js][ip][ic],
                                              input.ck_reptran_option);
                  CHKERR (status);
                }
              }
            }

            if (input.rte.mc.jacobian[DIM_1D]) {

              for (int isp = 0; isp < input.n_caoth + 2; isp++) {

                for (int i = 0; i < 2; i++) {

                  for (int lc = 0; lc < output->atm.nlev - 1; lc++) {

                    status = weighting_reptran (&output->wl,
                                                0,
                                                output->jacobian_r[lu][is][js][isp][i][lc],
                                                output->jacobian_t[lu][is][js][isp][i][lc],
                                                input.ck_reptran_option);
                    CHKERR (status);
                  }
                }
              }
            }

            if (input.rte.mc.backward.absorption) {
              status = weighting_reptran (&output->wl,
                                          0,
                                          output->absback3d_r[lu][is][js],
                                          output->absback3d_t[lu][is][js],
                                          input.ck_reptran_option);
              CHKERR (status);
            }

            /* variances */
            if (input.rte.mc.std) {
              status = weighting_reptran (&output->wl,
                                          1,
                                          output->rfldir3d_var_r[lu][is][js],
                                          output->rfldir3d_var_t[lu][is][js],
                                          input.ck_reptran_option);
              CHKERR (status);
              status = weighting_reptran (&output->wl,
                                          1,
                                          output->rfldn3d_var_r[lu][is][js],
                                          output->rfldn3d_var_t[lu][is][js],
                                          input.ck_reptran_option);
              CHKERR (status);
              status = weighting_reptran (&output->wl,
                                          1,
                                          output->flup3d_var_r[lu][is][js],
                                          output->flup3d_var_t[lu][is][js],
                                          input.ck_reptran_option);
              CHKERR (status);
              status = weighting_reptran (&output->wl,
                                          1,
                                          output->uavgso3d_var_r[lu][is][js],
                                          output->uavgso3d_var_t[lu][is][js],
                                          input.ck_reptran_option);
              CHKERR (status);
              status = weighting_reptran (&output->wl,
                                          1,
                                          output->uavgdn3d_var_r[lu][is][js],
                                          output->uavgdn3d_var_t[lu][is][js],
                                          input.ck_reptran_option);
              CHKERR (status);
              status = weighting_reptran (&output->wl,
                                          1,
                                          output->uavgup3d_var_r[lu][is][js],
                                          output->uavgup3d_var_t[lu][is][js],
                                          input.ck_reptran_option);
              CHKERR (status);

              if (output->radiance3d_var_r) {
                for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
                  status = weighting_reptran (&output->wl,
                                              1,
                                              output->radiance3d_var_r[lu][is][js][ip],
                                              output->radiance3d_var_t[lu][is][js][ip],
                                              input.ck_reptran_option);
                  CHKERR (status);
                }
              }

              if (input.rte.mc.backward.absorption) {
                status = weighting_reptran (&output->wl,
                                            1,
                                            output->absback3d_var_r[lu][is][js],
                                            output->absback3d_var_t[lu][is][js],
                                            input.ck_reptran_option);
                CHKERR (status);
              }
            }
          }
        }
      }

      output->u0u_t[lu] = calloc (input.rte.numu, sizeof (float*));
      for (int iu = 0; iu < input.rte.numu; iu++) {
        output->u0u_t[lu][iu] = calloc (output->wl.nlambda_t, sizeof (float));

        status = weighting_reptran (&output->wl, 0, output->u0u_r[lu][iu], output->u0u_t[lu][iu], input.ck_reptran_option);
        CHKERR (status);
      }

      output->uu_t[lu] = calloc (input.rte.nphi, sizeof (float**));
      for (int j = 0; j < input.rte.nphi; j++) {
        output->uu_t[lu][j] = calloc (input.rte.numu, sizeof (float*));
        for (int iu = 0; iu < input.rte.numu; iu++) {
          output->uu_t[lu][j][iu] = calloc (output->wl.nlambda_t, sizeof (float));

          status = weighting_reptran (&output->wl, 0, output->uu_r[lu][j][iu], output->uu_t[lu][j][iu], input.ck_reptran_option);
          CHKERR (status);
        }
      }
    } /* end of loop over layers */

    if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
        input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
      const size_t N3d_layers = vec_count_true (output->atm.Nzcld, output->atm.threed);
      for (size_t k = 0; k < N3d_layers; k++) {
        for (int is = 0; is < output->atm.Nxcld; is++) {
          for (int js = 0; js < output->atm.Nycld; js++) {

            status =
              weighting_reptran (&output->wl, 0, output->abs3d_r[k][is][js], output->abs3d_t[k][is][js], input.ck_reptran_option);
            CHKERR (status);

            if (input.rte.mc.std) { /* **CK added for forward mc_std */
              status = weighting_reptran (&output->wl,
                                          1,
                                          output->abs3d_var_r[k][is][js],
                                          output->abs3d_var_t[k][is][js],
                                          input.ck_reptran_option);
              CHKERR (status);
            }
          }
        }
      }
    }

    {  // copy (and weight) triangle_result_r -> triangle_results_t
      int ierr;
      ierr = init_spectral_triangular_surface_result_struct (output->wl.nlambda_t,
                                                             output->mc.triangular_surface.N_triangles,
                                                             &(output->triangle_results_t));
      CHKERR (ierr);

      // because the triangle_result output does not have wavelength as leading dimension
      // we have this wrapper here to copy and weight triangle_result_r -> triangle_results_t
      for (size_t id = 0; id < output->mc.triangular_surface.N_triangles; ++id) {
        const size_t Nvar = 6;
        float        var_r[Nvar][output->wl.nlambda_r];
        float        var_t[Nvar][output->wl.nlambda_t];
        for (size_t iv_r = 0; iv_r < output->wl.nlambda_r; ++iv_r) {
          var_r[0][iv_r] = output->triangle_results_r[iv_r]->ndir[id];
          var_r[1][iv_r] = output->triangle_results_r[iv_r]->ndn[id];
          var_r[2][iv_r] = output->triangle_results_r[iv_r]->nup[id];
          var_r[3][iv_r] = output->triangle_results_r[iv_r]->edir[id];
          var_r[4][iv_r] = output->triangle_results_r[iv_r]->edn[id];
          var_r[5][iv_r] = output->triangle_results_r[iv_r]->eup[id];
        }
        for (size_t ivar = 0; ivar < Nvar; ++ivar) {
          const int ierr = weighting_reptran (&output->wl, 0, var_r[ivar], var_t[ivar], input.ck_reptran_option);
          CHKERR (ierr);
        }
        for (size_t iv_t = 0; iv_t < output->wl.nlambda_t; ++iv_t) {
          output->triangle_results_t[iv_t]->ndir[id] = var_t[0][iv_t];
          output->triangle_results_t[iv_t]->ndn[id]  = var_t[1][iv_t];
          output->triangle_results_t[iv_t]->nup[id]  = var_t[2][iv_t];
          output->triangle_results_t[iv_t]->edir[id] = var_t[3][iv_t];
          output->triangle_results_t[iv_t]->edn[id]  = var_t[4][iv_t];
          output->triangle_results_t[iv_t]->eup[id]  = var_t[5][iv_t];
        }
      }
    }

    /* scale output to user-requested output unit
         * (when no representative wavelengths are used, this was already done in solve_rte())
         * */
    for (int iv = 0; iv < output->wl.nlambda_t; iv++) {

      if (input.source == SRC_THERMAL && output->spectrum_unit == UNIT_PER_CM_1) {

        switch (input.output_unit) {
        case UNIT_PER_CM_1:
          unit_factor = 1;
          break;
        case UNIT_PER_NM:
          unit_factor = 1.0e+7 / (output->wl.lambda_t[iv] * output->wl.lambda_t[iv]);
          break;
        case UNIT_PER_BAND:
          unit_factor = output->wl.width_of_reptran_band[output->wl.reptran_band_t[iv]];
          break;
        case UNIT_NOT_DEFINED:
          unit_factor = 1.0;
          break;
        default:
          fprintf (stderr, "Error: Program bug, unsupported output unit %d\n", input.output_unit);
          CHKERR (-1);
        }

      } else if (input.source == SRC_SOLAR && output->spectrum_unit == UNIT_PER_NM) {
        switch (input.output_unit) {
        case UNIT_PER_CM_1:
          unit_factor = (output->wl.lambda_t[iv] * output->wl.lambda_t[iv]) / 1.0e+7;
          break;
        case UNIT_PER_NM:
          unit_factor = 1;
          break;
        case UNIT_PER_BAND:
          unit_factor = output->wl.width_of_reptran_band[output->wl.reptran_band_t[iv]];
          break;
        case UNIT_NOT_DEFINED:
          unit_factor = 1.0;
          break;
        default:
          fprintf (stderr, "Error: Program bug, unsupported output unit %d\n", input.output_unit);
          CHKERR (-1);
        }

      } else {

        fprintf (stderr,
                 "Error: Program bug, unsupported comination of spectrum unit %d with source %d\n",
                 input.spectrum_unit,
                 input.source);
        CHKERR (-1);
      }

      ffactor = unit_factor;
      rfactor = unit_factor;
      hfactor = unit_factor;
      status  = scale_output (input,
                             &(output->rfldir_t),
                             &(output->rfldn_t),
                             &(output->flup_t),
                             &(output->albmed_t),
                             &(output->trnmed_t),
                             &(output->uavgso_t),
                             &(output->uavgdn_t),
                             &(output->uavgup_t),
                             &(output->uavg_t),
                             &(output->u0u_t),
                             &(output->uu_t),
                             &(output->heat_t),
                             &(output->emis_t),
                             &(output->w_zout_t),
                             &(output->rfldir3d_t),
                             &(output->rfldn3d_t),
                             &(output->flup3d_t),
                             &(output->fl3d_is_t),
                             &(output->uavgso3d_t),
                             &(output->uavgdn3d_t),
                             &(output->uavgup3d_t),
                             &(output->radiance3d_t),
                             &(output->jacobian_t),
                             &(output->absback3d_t),
                             &(output->rfldir3d_var_t),
                             &(output->rfldn3d_var_t),
                             &(output->flup3d_var_t),
                             &(output->uavgso3d_var_t),
                             &(output->uavgdn3d_var_t),
                             &(output->uavgup3d_var_t),
                             &(output->radiance3d_var_t),
                             &(output->abs3d_var_t),
                             &(output->absback3d_var_t),
                             output->atm.nzout,
                             output->atm.Nxcld,
                             output->atm.Nycld,
                             output->atm.Nzcld,
                             output->mc.alis.Nc,
                             output->atm.nlev - 1,
                             output->atm.threed,
                             output->islower,
                             output->isupper,
                             output->jslower,
                             output->jsupper,
                             output->isstep,
                             output->jsstep,
                             &(output->abs3d_t),
                             output->triangle_results_t,
                             ffactor,
                             rfactor,
                             hfactor,
                             iv); /* in ancillary.c */
      CHKERR (status);
    }  // iv
  }

  return 0;
}

int closest (float lambda, float* lambda_raw, int n_crs) {
  int iv = 0;

  float min_delta = FLT_MAX;

  int result = -1;

  for (iv = 0; iv < n_crs; iv++) {
    if (fabs (lambda - lambda_raw[iv]) < min_delta) {
      min_delta = fabs (lambda - lambda_raw[iv]);
      result    = iv;
    }
  }
  return result;
}

/*****************************************************************************/
/* Interpolate the transmittance from the transmission wavelength grid to    */
/* high resolution grid.                                                     */
/*****************************************************************************/

int interpolate_transmittance (input_struct input, output_struct* output) {
  int linear = 1, status;

  /* If only one wavelength is required, or in correlated-k mode, */
  /* arrays are not interpolated but copied                       */
  /* ??? need to add the condition that the extraterrestrial  ??? */
  /* ??? spectrum was not read from file but set to 1         ??? */

  if (output->wl.nlambda_t == 1 || output->wl.use_reptran == 1 ||
      (output->wl.use_reptran == 0 && input.ck_scheme != CK_CRS && input.ck_scheme != CK_RAMAN && input.ck_scheme != CK_LOWTRAN)) {

    for (int ivh = 0; ivh < output->wl.nlambda_h; ivh++) {

      int iv;
      if (output->wl.nlambda_t == 1)
        iv = ivh;
      else if (output->wl.use_reptran == 1)
        iv = closest (output->wl.lambda_h[ivh], output->wl.lambda_t, output->wl.nlambda_t);
      else
        iv = output->wl.map_e2h[ivh];

      for (int iu = 0; iu < input.rte.numu; iu++) {
        output->albmed[iu][ivh] = output->albmed_t[iu][iv];
        output->trnmed[iu][ivh] = output->trnmed_t[iu][iv];
      }

      for (int lu = 0; lu < output->atm.nzout; lu++) {
        output->flup[lu][ivh]            = output->flup_t[lu][iv];
        output->rfldir[lu][ivh]          = output->rfldir_t[lu][iv];
        output->rfldn[lu][ivh]           = output->rfldn_t[lu][iv];
        output->uavg[lu][ivh]            = output->uavg_t[lu][iv];
        output->uavgdn[lu][ivh]          = output->uavgdn_t[lu][iv];
        output->uavgso[lu][ivh]          = output->uavgso_t[lu][iv];
        output->uavgup[lu][ivh]          = output->uavgup_t[lu][iv];
        output->sslidar_nphot[lu][ivh]   = output->sslidar_nphot_t[lu][iv];
        output->sslidar_nphot_q[lu][ivh] = output->sslidar_nphot_q_t[lu][iv];
        output->sslidar_ratio[lu][ivh]   = output->sslidar_ratio_t[lu][iv];
        if (input.heating != HEAT_NONE) {
          output->heat[lu][ivh]   = output->heat_t[lu][iv];
          output->emis[lu][ivh]   = output->emis_t[lu][iv];
          output->w_zout[lu][ivh] = output->w_zout_t[lu][iv];
        }

        /* radiances */
        for (int iu = 0; iu < input.rte.numu; iu++) {
          output->u0u[lu][iu][ivh] = output->u0u_t[lu][iu][iv];

          for (int j = 0; j < input.rte.nphi; j++)
            output->uu[lu][j][iu][ivh] = output->uu_t[lu][j][iu][iv];
        }

        /* 3D fields */
        if (input.output3D) {
          for (int is = output->islower; is <= output->isupper; is += output->isstep) {
            for (int js = output->jslower; js <= output->jsupper; js += output->jsstep) {
              output->rfldir3d[lu][is][js][ivh] = output->rfldir3d_t[lu][is][js][iv];
              output->rfldn3d[lu][is][js][ivh]  = output->rfldn3d_t[lu][is][js][iv];
              output->flup3d[lu][is][js][ivh]   = output->flup3d_t[lu][is][js][iv];
              output->uavgso3d[lu][is][js][ivh] = output->uavgso3d_t[lu][is][js][iv];
              output->uavgdn3d[lu][is][js][ivh] = output->uavgdn3d_t[lu][is][js][iv];
              output->uavgup3d[lu][is][js][ivh] = output->uavgup3d_t[lu][is][js][iv];
              if (output->radiance3d_t) {
                for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
                  for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
                    output->radiance3d[lu][is][js][ip][ic][ivh] = output->radiance3d_t[lu][is][js][ip][ic][iv];
                  }
                }
              }

              if (input.rte.mc.concentration_is || input.rte.mc.spectral_is) {
                for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
                  output->fl3d_is[lu][is][js][ic][ivh] = output->fl3d_is_t[lu][is][js][ic][iv];
                }
              }

              if (input.rte.mc.jacobian[DIM_1D]) {
                for (int isp = 0; isp < input.n_caoth + 2; isp++) {
                  for (int i = 0; i < 2; i++) {
                    for (int lc = 0; lc < output->atm.nlev - 1; lc++) {
                      output->jacobian[lu][is][js][isp][i][lc][ivh] = output->jacobian_t[lu][is][js][isp][i][lc][iv];
                    }
                  }
                }
              }

              if (input.rte.mc.backward.absorption) {
                output->absback3d[lu][is][js][ivh] = output->absback3d_t[lu][is][js][iv];
              }

              /* variances */
              if (input.rte.mc.std) {

                output->rfldir3d_var[lu][is][js][ivh] = output->rfldir3d_var_t[lu][is][js][iv];
                output->rfldn3d_var[lu][is][js][ivh]  = output->rfldn3d_var_t[lu][is][js][iv];
                output->flup3d_var[lu][is][js][ivh]   = output->flup3d_var_t[lu][is][js][iv];
                output->uavgso3d_var[lu][is][js][ivh] = output->uavgso3d_var_t[lu][is][js][iv];
                output->uavgdn3d_var[lu][is][js][ivh] = output->uavgdn3d_var_t[lu][is][js][iv];
                output->uavgup3d_var[lu][is][js][ivh] = output->uavgup3d_var_t[lu][is][js][iv];
                if (output->radiance3d_var_t) {
                  for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
                    output->radiance3d_var[lu][is][js][ip][ivh] = output->radiance3d_var_t[lu][is][js][ip][iv];
                  }
                }

                if (input.rte.mc.backward.absorption) {
                  output->absback3d_var[lu][is][js][ivh] = output->absback3d_var_t[lu][is][js][iv];
                }
              }
            }
          }
        }
      }

      /* absorption is defined on the 3D caoth grid, not on the user-defined grid */
      if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
          input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
        const size_t N3d_layers = vec_count_true (output->atm.Nzcld, output->atm.threed);
        for (size_t k = 0; k < N3d_layers; k++) {
          for (int is = 0; is < output->atm.Nxcld; is++) {
            for (int js = 0; js < output->atm.Nycld; js++) {
              output->abs3d[k][is][js][ivh] = output->abs3d_t[k][is][js][iv];
              if (input.rte.mc.std) {
                output->abs3d_var[k][is][js][ivh] = output->abs3d_var_t[k][is][js][iv];
              }
            }
          }
        }
      }

      /* copy solar zenith angle */
      output->sza_h[ivh] = output->atm.sza_t[iv];

      if (output->triangle_results_t) {  // link triangle_results_t -> triangle_results_h
        if (!output->triangle_results_o)
          output->triangle_results_o = malloc (output->wl.nlambda_h * sizeof (t_triangle_radiation_field*));
        output->triangle_results_o[ivh] = output->triangle_results_t[iv];
      }
    }

  } else { /* interpolation to different wavelengths required */
    // wavelength offset because output->rfldir3d_t etc. may be defined from output->wl.iv_lower to output->wl.iv_upper
    const size_t iv_off = output->wl.iv_lower;
    const size_t Niv_t  = output->wl.iv_upper - output->wl.iv_lower + 1;
    const size_t Niv_h  = output->wl.nlambda_h;
    float*       iv_t   = &output->wl.lambda_t[iv_off];
    float*       iv_h   = output->wl.lambda_h;

    /* interpolate solar zenith angle */
    status = arb_wvn (Niv_t, iv_t, &output->atm.sza_t[iv_off], Niv_h, iv_h, output->sza_h, linear, 0);
    CHKERR (status);

    for (int iu = 0; iu < input.rte.numu; iu++) {
      status = arb_wvn (Niv_t, iv_t, &output->albmed_t[iu][iv_off], Niv_h, iv_h, output->albmed[iu], linear, 0);
      CHKERR (status);
      status = arb_wvn (Niv_t, iv_t, &output->trnmed_t[iu][iv_off], Niv_h, iv_h, output->trnmed[iu], linear, 0);
      CHKERR (status);
    }

    for (int lu = 0; lu < output->atm.nzout; lu++) {

      status = arb_wvn (Niv_t, iv_t, &output->flup_t[lu][iv_off], Niv_h, iv_h, output->flup[lu], linear, 0);
      CHKERR (status);
      status = arb_wvn (Niv_t, iv_t, &output->rfldir_t[lu][iv_off], Niv_h, iv_h, output->rfldir[lu], linear, 0);
      CHKERR (status);
      status = arb_wvn (Niv_t, iv_t, &output->rfldn_t[lu][iv_off], Niv_h, iv_h, output->rfldn[lu], linear, 0);
      CHKERR (status);
      status = arb_wvn (Niv_t, iv_t, &output->uavg_t[lu][iv_off], Niv_h, iv_h, output->uavg[lu], linear, 0);
      CHKERR (status);
      status = arb_wvn (Niv_t, iv_t, &output->uavgdn_t[lu][iv_off], Niv_h, iv_h, output->uavgdn[lu], linear, 0);
      CHKERR (status);
      status = arb_wvn (Niv_t, iv_t, &output->uavgso_t[lu][iv_off], Niv_h, iv_h, output->uavgso[lu], linear, 0);
      CHKERR (status);
      status = arb_wvn (Niv_t, iv_t, &output->uavgup_t[lu][iv_off], Niv_h, iv_h, output->uavgup[lu], linear, 0);
      CHKERR (status);
      status = arb_wvn (Niv_t, iv_t, &output->sslidar_nphot_t[lu][iv_off], Niv_h, iv_h, output->sslidar_nphot[lu], linear, 0);
      CHKERR (status);
      status = arb_wvn (Niv_t, iv_t, &output->sslidar_nphot_q_t[lu][iv_off], Niv_h, iv_h, output->sslidar_nphot_q[lu], linear, 0);
      CHKERR (status);
      status = arb_wvn (Niv_t, iv_t, &output->sslidar_ratio_t[lu][iv_off], Niv_h, iv_h, output->sslidar_ratio[lu], linear, 0);
      CHKERR (status);

      if (input.heating != HEAT_NONE) {
        status = arb_wvn (Niv_t, iv_t, &output->heat_t[lu][iv_off], Niv_h, iv_h, output->heat[lu], linear, 0);
        CHKERR (status);
        status = arb_wvn (Niv_t, iv_t, &output->emis_t[lu][iv_off], Niv_h, iv_h, output->emis[lu], linear, 0);
        CHKERR (status);
        status = arb_wvn (Niv_t, iv_t, &output->w_zout_t[lu][iv_off], Niv_h, iv_h, output->w_zout[lu], linear, 0);
        CHKERR (status);
      }

      /* 3D fields */
      if (input.output3D) {
        for (int is = output->islower; is <= output->isupper; is += output->isstep) {
          for (int js = output->jslower; js <= output->jsupper; js += output->jsstep) {
            status =
              arb_wvn (Niv_t, iv_t, &output->rfldir3d_t[lu][is][js][iv_off], Niv_h, iv_h, output->rfldir3d[lu][is][js], linear, 0);
            CHKERR (status);

            status =
              arb_wvn (Niv_t, iv_t, &output->rfldn3d_t[lu][is][js][iv_off], Niv_h, iv_h, output->rfldn3d[lu][is][js], linear, 0);
            CHKERR (status);

            status =
              arb_wvn (Niv_t, iv_t, &output->flup3d_t[lu][is][js][iv_off], Niv_h, iv_h, output->flup3d[lu][is][js], linear, 0);
            CHKERR (status);

            status =
              arb_wvn (Niv_t, iv_t, &output->uavgso3d_t[lu][is][js][iv_off], Niv_h, iv_h, output->uavgso3d[lu][is][js], linear, 0);
            CHKERR (status);

            status =
              arb_wvn (Niv_t, iv_t, &output->uavgdn3d_t[lu][is][js][iv_off], Niv_h, iv_h, output->uavgdn3d[lu][is][js], linear, 0);
            CHKERR (status);

            status =
              arb_wvn (Niv_t, iv_t, &output->uavgup3d_t[lu][is][js][iv_off], Niv_h, iv_h, output->uavgup3d[lu][is][js], linear, 0);
            CHKERR (status);

            if (output->radiance3d_t) {
              for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
                for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
                  status = arb_wvn (Niv_t,
                                    iv_t,
                                    &output->radiance3d_t[lu][is][js][ip][ic][iv_off],
                                    Niv_h,
                                    iv_h,
                                    output->radiance3d[lu][is][js][ip][ic],
                                    linear,
                                    0);
                  CHKERR (status);
                }
              }
            }

            if (input.rte.mc.concentration_is || input.rte.mc.spectral_is) {
              for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
                status = arb_wvn (Niv_t,
                                  iv_t,
                                  &output->fl3d_is_t[lu][is][js][ic][iv_off],
                                  Niv_h,
                                  iv_h,
                                  output->fl3d_is[lu][is][js][ic],
                                  linear,
                                  0);
                CHKERR (status);
              }
            }

            if (input.rte.mc.jacobian[DIM_1D])
              for (int isp = 0; isp < input.n_caoth + 2; isp++)
                for (int i = 0; i < 2; i++)
                  for (int lc = 0; lc < output->atm.nlev - 1; lc++) {
                    status = arb_wvn (Niv_t,
                                      iv_t,
                                      &output->jacobian_t[lu][is][js][isp][i][lc][iv_off],
                                      Niv_h,
                                      iv_h,
                                      output->jacobian[lu][is][js][isp][i][lc],
                                      linear,
                                      0);
                    CHKERR (status);
                  }

            if (input.rte.mc.backward.absorption) {
              status = arb_wvn (Niv_t,
                                iv_t,
                                &output->absback3d_t[lu][is][js][iv_off],
                                Niv_h,
                                iv_h,
                                output->absback3d[lu][is][js],
                                linear,
                                0);
              CHKERR (status);
            }

            /* variances */
            if (input.rte.mc.std) {

              status = arb_wvn (Niv_t,
                                iv_t,
                                &output->rfldir3d_var_t[lu][is][js][iv_off],
                                Niv_h,
                                iv_h,
                                output->rfldir3d_var[lu][is][js],
                                linear,
                                0);
              CHKERR (status);

              status = arb_wvn (Niv_t,
                                iv_t,
                                &output->rfldn3d_var_t[lu][is][js][iv_off],
                                Niv_h,
                                iv_h,
                                output->rfldn3d_var[lu][is][js],
                                linear,
                                0);
              CHKERR (status);

              status = arb_wvn (Niv_t,
                                iv_t,
                                &output->flup3d_var_t[lu][is][js][iv_off],
                                Niv_h,
                                iv_h,
                                output->flup3d_var[lu][is][js],
                                linear,
                                0);
              CHKERR (status);

              status = arb_wvn (Niv_t,
                                iv_t,
                                &output->uavgso3d_var_t[lu][is][js][iv_off],
                                Niv_h,
                                iv_h,
                                output->uavgso3d_var[lu][is][js],
                                linear,
                                0);
              CHKERR (status);

              status = arb_wvn (Niv_t,
                                iv_t,
                                &output->uavgdn3d_var_t[lu][is][js][iv_off],
                                Niv_h,
                                iv_h,
                                output->uavgdn3d_var[lu][is][js],
                                linear,
                                0);
              CHKERR (status);

              status = arb_wvn (Niv_t,
                                iv_t,
                                &output->uavgup3d_var_t[lu][is][js][iv_off],
                                Niv_h,
                                iv_h,
                                output->uavgup3d_var[lu][is][js],
                                linear,
                                0);
              CHKERR (status);

              if (output->radiance3d_var_t) {
                for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
                  status = arb_wvn (Niv_t,
                                    iv_t,
                                    &output->radiance3d_var_t[lu][is][js][ip][iv_off],
                                    Niv_h,
                                    iv_h,
                                    output->radiance3d_var[lu][is][js][ip],
                                    linear,
                                    0);
                  CHKERR (status);
                }
              }

              if (input.rte.mc.backward.absorption) {
                status = arb_wvn (Niv_t,
                                  iv_t,
                                  &output->absback3d_var_t[lu][is][js][iv_off],
                                  Niv_h,
                                  iv_h,
                                  output->absback3d_var[lu][is][js],
                                  linear,
                                  0);
                CHKERR (status);
              }
            }
          }
        }
      }

      for (int iu = 0; iu < input.rte.numu; iu++) {
        status = arb_wvn (Niv_t, iv_t, &output->u0u_t[lu][iu][iv_off], Niv_h, iv_h, output->u0u[lu][iu], linear, 0);
        CHKERR (status);
      }

      for (int j = 0; j < input.rte.nphi; j++) {
        for (int iu = 0; iu < input.rte.numu; iu++) {

          status = arb_wvn (Niv_t, iv_t, &output->uu_t[lu][j][iu][iv_off], Niv_h, iv_h, output->uu[lu][j][iu], linear, 0);
          CHKERR (status);
        }
      }
    }

    if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
        input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
      const size_t N3d_layers = vec_count_true (output->atm.Nzcld, output->atm.threed);
      for (size_t k = 0; k < N3d_layers; k++) {
        for (int is = 0; is < output->atm.Nxcld; is++) {
          for (int js = 0; js < output->atm.Nycld; js++) {

            status = arb_wvn (Niv_t, iv_t, &output->abs3d_t[k][is][js][iv_off], Niv_h, iv_h, output->abs3d[k][is][js], linear, 0);
            CHKERR (status);

            if (input.rte.mc.std) { /* **CK added for forward mc_std */
              status = arb_wvn (Niv_t,
                                iv_t,
                                &output->abs3d_var_t[k][is][js][iv_off],
                                Niv_h,
                                iv_h,
                                output->abs3d_var[k][is][js],
                                linear,
                                0);
              CHKERR (status);
            }
          }
        }
      }
    }
  }

  return 0;
}

/**************************************************************/
/* Multiply the transmittance with the extraterrestrial       */
/* irradiance.                                                */
/**************************************************************/

int multiply_extraterrestrial (input_struct input, output_struct* output) {
  int    iv = 0, status = 0;
  double ffactor = 0, rfactor = 0, hfactor = 0;
  double watt_factor = 1.0;
  double hfactor2    = 0.0;

  /* Figure out if we have to convert from mW to W. This should       */
  /* cover the input stuff that comes with uvspec. Now, if the user   */
  /* decides to change the input solar flux it may get interesting... */

  /* THIS mW to W CONVERTION IS USED FOR HEATING RATES, SHOULD IT BECOME DEFAULT? */
  switch (input.ck_scheme) {
  case CK_FU:
  case CK_KATO:
  case CK_KATO2:
  case CK_KATO2_96:
  case CK_KATO2ANDWANDJI:
  case CK_AVHRR_KRATZ:
    watt_factor = 1.0;
    break;
  case CK_LOWTRAN:
  case CK_FILE:
  case CK_CRS:
  case CK_REPTRAN:
  case CK_REPTRAN_CHANNEL:
  case CK_RAMAN:
    switch (input.source) {
    case SRC_SOLAR:
    case SRC_BLITZ: /* BCA */
    case SRC_LIDAR: /* BCA */
      watt_factor = mW2W;
      break;
    case SRC_THERMAL:
      watt_factor = 1.0;
      break;
    default:
      fprintf (stderr, "Error, unknown source %d\n", input.source);
      return -1;
    }
    break;
  default:
    fprintf (stderr, "Error: unsupported correlated-k scheme %d\n", input.ck_scheme);
    return -1;
    break;
  }

  if (input.ck_repint_option == REPINT_MEIER50)
    watt_factor = 1.0;

  hfactor  = watt_factor * s2day;
  hfactor2 = hfactor;

  /* multiply with extraterrestrial irradiance */
  for (iv = 0; iv < output->wl.nlambda_h; iv++) {

    switch (input.source) {
    case SRC_THERMAL:
      ffactor = output->wl.filter[iv];
      rfactor = output->wl.filter[iv];
      break;

    case SRC_SOLAR:
    case SRC_BLITZ: /* BCA */
    case SRC_LIDAR: /* BCA */
      switch (input.processing) {
      case PROCESS_INT:
      case PROCESS_SUM:
      case PROCESS_RGB:
      case PROCESS_RGBRAW:
        ffactor = output->wl.fbeam[iv] * output->sunshine_fraction * output->wl.filter[iv];
        rfactor = output->wl.fbeam[iv] * output->sunshine_fraction * output->wl.filter[iv];
        break;

      case PROCESS_NONE:
      case PROCESS_RAMAN:
        switch (input.calibration) {
        case OUTCAL_ABSOLUTE:
          ffactor = output->wl.fbeam[iv] * output->sunshine_fraction * output->wl.filter[iv];
          rfactor = output->wl.fbeam[iv] * output->sunshine_fraction * output->wl.filter[iv];
          break;

        case OUTCAL_TRANSMITTANCE:
          ffactor = 1.0 * output->wl.filter[iv];
          rfactor = 1.0 * output->wl.filter[iv];
          break;

        case OUTCAL_REFLECTIVITY:
          ffactor = 1.0 / cos (output->sza_h[iv] * M_PI / 180.0) * output->wl.filter[iv];
          rfactor = M_PI / cos (output->sza_h[iv] * M_PI / 180.0) * output->wl.filter[iv];
          break;

        default:
          fprintf (stderr, "Error, unknown output calibration %d\n", input.calibration);
          return -1;
        }

        break;

      default:
        fprintf (stderr, "Error, unknown output processing %d\n", input.processing);
        return -1;
      }

      break;

    default:
      fprintf (stderr, "Error, unknown source %d\n", input.source);
      return -1;
    }

    switch (input.source) {
    case SRC_SOLAR:
    case SRC_LIDAR: /* BCA */
    case SRC_BLITZ: /* BCA */
      hfactor = output->wl.fbeam[iv] * output->sunshine_fraction * hfactor2;
      break;
    case SRC_THERMAL:
      break;
    default:
      fprintf (stderr, "Error, unknown source %d\n", input.source);
      return -1;
    }

    /*****************************************************************************************/
    /* now scale irradiances with ffactor, radiances with rfactor, heating rate with hfactor */
    /*****************************************************************************************/
    status = scale_output (input,
                           &(output->rfldir),
                           &(output->rfldn),
                           &(output->flup),
                           &(output->albmed),
                           &(output->trnmed),
                           &(output->uavgso),
                           &(output->uavgdn),
                           &(output->uavgup),
                           &(output->uavg),
                           &(output->u0u),
                           &(output->uu),
                           &(output->heat),
                           &(output->emis),
                           &(output->w_zout),
                           &(output->rfldir3d),
                           &(output->rfldn3d),
                           &(output->flup3d),
                           &(output->fl3d_is),
                           &(output->uavgso3d),
                           &(output->uavgdn3d),
                           &(output->uavgup3d),
                           &(output->radiance3d),
                           &(output->jacobian),
                           &(output->absback3d),
                           &(output->rfldir3d_var),
                           &(output->rfldn3d_var),
                           &(output->flup3d_var),
                           &(output->uavgso3d_var),
                           &(output->uavgdn3d_var),
                           &(output->uavgup3d_var),
                           &(output->radiance3d_var),
                           &(output->abs3d_var),
                           &(output->absback3d_var),
                           output->atm.nzout,
                           output->atm.Nxcld,
                           output->atm.Nycld,
                           output->atm.Nzcld,
                           output->mc.alis.Nc,
                           output->atm.nlyr,
                           output->atm.threed,
                           output->islower,
                           output->isupper,
                           output->jslower,
                           output->jsupper,
                           output->isstep,
                           output->jsstep,
                           &(output->abs3d),
                           output->triangle_results_o,
                           ffactor,
                           rfactor,
                           hfactor,
                           iv);
    CHKERR (status);
  }

  return 0;
}

// helper function to convert a caoth3d to consecutive memory and write it to a netcdf file
static int dump_caoth3d_var (int ncid, int isp, char* caothname, caoth3d_out_struct caoth3d, char* varname, float*** caoth_data) {
  const size_t Nx = caoth3d.Nx;
  const size_t Ny = caoth3d.Ny;
  const size_t Nz = caoth3d.nlyr;

  if (Nz == 0)
    return 0;  // nothing to do

  char varprefix[FILENAME_MAX];
  char ncvarname[FILENAME_MAX];
  char xlabel[FILENAME_MAX];
  char ylabel[FILENAME_MAX];
  char zlabel[FILENAME_MAX];
  char isp_char[11];

  snprintf (isp_char, sizeof (isp_char), "%d", isp);
  strcpy (varprefix, "caoth3d_");
  strcat (varprefix, isp_char);
  strcat (varprefix, "_");
  strcat (varprefix, caothname);

  strcpy (ncvarname, varprefix);
  strcat (ncvarname, varname);
  strcat (strcpy (ylabel, varprefix), "_Nx");
  strcat (strcpy (zlabel, varprefix), "_Ny");
  strcpy (xlabel, "nlyr");

  float*** data = calloc_float_3D (Nz, Nx, Ny, ncvarname);
  for (size_t k = 0; k < caoth3d.nlyr; ++k) {
    if (caoth3d.threed[k]) {
      for (size_t i = 0; i < Nx; ++i) {
        for (size_t j = 0; j < Ny; ++j) {
          data[k][i][j] = caoth_data[k][i][j];
        }
      }
    }
  }

  const int NDIMS       = 3;
  size_t    dimsizes[3] = {Nz, Nx, Ny};
  char**    dimnames    = malloc_char_2D (3, NC_MAX_NAME, "dim names for netcdf");
  strcpy (dimnames[0], xlabel);
  strcpy (dimnames[1], ylabel);
  strcpy (dimnames[2], zlabel);
  int status = write_netcdf_contiguous_ND_float (ncid, &(data[0][0][0]), ncvarname, NDIMS, dimsizes, dimnames);
  CHKERR (status);
  free_char_2D (dimnames);
  free_float_3D (data);
  return 0;
}

/***********************************************************************************/
/* Function: optical_properties                                           @61_30i@ */
/* Description:                                                                    */
/*  Calculates optical depth, single scattering albedo and                         */
/*  phase function from absorption and scattering cross sections                   */
/*  corresponding gas and particulate matter concentrations for                    */
/*  one wavelength.                                                                */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i61_30@ */
/***********************************************************************************/

int optical_properties (input_struct   input,
                        output_struct* output,
                        double         wvl,
                        int            ir,
                        int            iv1,
                        int            iv2,
                        int            iq,
                        int            verbose,
                        int            skip_optical_properties) {

  static int first = 1;  //TODO: Is this good, static int? in function interpolate_profile same variable first
  int        k = 0, lc = 0, ip = 0, ips = 0, iv = 0, iv1r = 0, iv2r = 0, isp = 0, ispo = 0;
  int        status = 0;

  /* vectors for ALL caoths; 0 is reserved for MOL, 1 for AER */
  double *babs_s = NULL, babs_tot = 0.0;
  double *bsca_s = NULL, bsca_tot = 0.0;

  double babs_mol_md = 0.0, babs_tot_md = 0.0;
  double bext_tot_md = 0.0; /* The absorption of gases minus the one specified with a
			     matrix using the dens_file command (for airmass factor calculations). */

  double *ssa_s_unsc = NULL, *bsca_s_unsc = NULL, *babs_s_unsc = NULL;
  double *bsca_s_unsc_int = NULL, *babs_s_unsc_int = NULL;
  double* gg_s_unsc = NULL;

  double* babs_s_int = NULL;
  double* bsca_s_int = NULL;

  double tau_babs_sum = 0.0; /* ulrike: to add up the dtau's of babsa,... for all layers above the level considered*/

  double** mom_s    = NULL;
  double   bext_tot = 0.0;
  int      nphamat  = 0;

  /* Output Extinction; needed for ARLEM */
  char  extfilename[FILENAME_MAX] = "";
  FILE* extfile                   = NULL;

  int nlev = output->atm.nlev_common;

  int nlyr          = nlev - 1;
  int n_caoth       = input.n_caoth + 2; /* mol and aer included */
  int n_caoth_alloc = 0;
  int i_wc = 0, i_ic = 0;
  int i_wcn = 0, i_wck = 0;
  int i_icn = 0, i_ick = 0;

  double *mom_s_g1 = NULL, *mom_s_g2 = NULL;

  double rayleigh_mom2 = 0;

  double  rayleigh_depol = 0.0, momk = 0.0;
  double *ssa_s = NULL, *dscale_s = NULL, *f_s = NULL, *ff_s = NULL, *g1_s = NULL, *g2_s = NULL, *dtau_s = NULL, *mom0_s = NULL;
  double  wvl1 = 0, wvl2 = 0;

  int*     tmp_ntheta_in = 0;
  float ** tmp_theta_in = NULL, *tmp_theta_new = NULL;
  double **tmp_mu_in = NULL, *tmp_mu_new = NULL;
  float ** tmp_phase_in = NULL, *tmp_phase_new = NULL;
  int      n_in = 0, n_tot = 0, i = 0, ntheta_new = 0;
  double*  interp_optprop_weight = NULL;
  int*     phase_calloced        = NULL;

  double one = 1.0;

  /* phase matrix element indices for polarized RT */
  int   ip_simp[6] = {0, 0, 0, 0, 0, 0};
  int   ip_full[6] = {0, 1, 2, 3, 4, 5};
  int   ip_red[6]  = {0, 1, 2, 3, 0, 2};
  int** ip_act     = NULL;

  n_caoth_alloc = n_caoth;
  i_wc          = input.i_wc + 2;
  i_ic          = input.i_ic + 2;
  /* in case either wc or ic is not existent, allocate a dummy  */
  /* caoth, which is zero, and can be used for verbose output */
  if (i_ic == 1 || i_wc == 1)
    n_caoth_alloc++;
  /* set index of wc/ic to dummy caoth */
  if (i_wc == 1)
    i_wc = n_caoth_alloc - 1;
  if (i_ic == 1)
    i_ic = n_caoth_alloc - 1;

  /* allocate caoth arrays */
  babs_s     = calloc (n_caoth_alloc, sizeof (double));
  bsca_s     = calloc (n_caoth_alloc, sizeof (double));
  ssa_s      = calloc (n_caoth_alloc, sizeof (double));
  babs_s_int = calloc (n_caoth_alloc, sizeof (double));
  bsca_s_int = calloc (n_caoth_alloc, sizeof (double));

  mom_s    = calloc (n_caoth_alloc, sizeof (double*));
  mom_s_g1 = calloc (n_caoth_alloc, sizeof (double));
  mom_s_g2 = calloc (n_caoth_alloc, sizeof (double));
  dscale_s = calloc (n_caoth_alloc, sizeof (double));
  f_s      = calloc (n_caoth_alloc, sizeof (double));
  ff_s     = calloc (n_caoth_alloc, sizeof (double));
  g1_s     = calloc (n_caoth_alloc, sizeof (double));
  g2_s     = calloc (n_caoth_alloc, sizeof (double));
  dtau_s   = calloc (n_caoth_alloc, sizeof (double));
  mom0_s   = calloc (n_caoth_alloc, sizeof (double));

  babs_s_unsc     = calloc (n_caoth_alloc, sizeof (double));
  bsca_s_unsc     = calloc (n_caoth_alloc, sizeof (double));
  ssa_s_unsc      = calloc (n_caoth_alloc, sizeof (double));
  babs_s_unsc_int = calloc (n_caoth_alloc, sizeof (double));
  bsca_s_unsc_int = calloc (n_caoth_alloc, sizeof (double));
  gg_s_unsc       = calloc (n_caoth_alloc, sizeof (double));

  ip_act = calloc (n_caoth_alloc, sizeof (int*));

  /* */

  if (input.rte.pol_nstokes == 1) {
    nphamat = 1;
  } else if (input.rte.pol_nstokes > 1 && input.rte.pol_nstokes < 5) {
    nphamat = 6;
  } else {
    fprintf (stderr, "Input variable pol_nstokes is wrong! %d\n", input.rte.pol_nstokes);
    CHKERR (-1);
  }

  output->nphamat = nphamat;

  if (nphamat == 6) {
    /* phase matrix element indices for last two elements in
       case of spherical symmetric particles */
    if (output->aer.optprop.nphamat == 4)
      ip_act[CAOTH_AER] = ip_red;
    else
      ip_act[CAOTH_AER] = ip_full;

    for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
      ispo = isp - CAOTH_FIR;
      /* phase matrix element indices for last two elements in
	 case of spherical symmetric particles */
      if (output->caoth[ispo].optprop.nphamat == 4)
        ip_act[isp] = ip_red;
      else
        ip_act[isp] = ip_full;
    }

  } /* end if nphamat == 6 */
  else
    for (isp = CAOTH_AER; isp < n_caoth; isp++)
      ip_act[isp] = ip_simp;

  for (isp = 0; isp < n_caoth; isp++)
    mom_s[isp] = calloc (nphamat, sizeof (double));

  iv = iv1; /* iv needed both for Raman and not Raman when calculating phase function moments */
  if (input.raman) {
    iv1r = iv1 + output->wl.iv_lower;
    iv2r = iv2 + output->wl.iv_lower;
  }
  if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1)) {
    wvl1           = output->wl.lambda_r[iv1];
    wvl2           = output->wl.lambda_r[iv2];
    rayleigh_depol = linpol (wvl1, wvl2, output->rayleigh_depol[iv1r], output->rayleigh_depol[iv2r], wvl);
  } else {
    rayleigh_depol = output->rayleigh_depol[iv];
  }

  /* 2nd moment for Rayleigh scattering, including depolarization */
  rayleigh_mom2 = 0.2 * (1.0 - rayleigh_depol) / (2.0 + rayleigh_depol);

  /* Deltas in Eq. 2.16, Hansen and Travis, Space Science Rev., 16, 527-610, 1974.*/
  //Delta  = (1.0 - rayleigh_depol) / (1.0 + 0.5 * rayleigh_depol);
  //Deltap = (1.0 + 2.0 * rayleigh_depol) / (1.0 - rayleigh_depol);

  if (verbose) {
    fprintf (stderr, "... second moment for Rayleigh scattering: %f\n", rayleigh_mom2);
    fprintf (stderr, "*** optical_properties()\n");
    if (output->cf.nlev == 0 || input.rte.solver == SOLVER_TWOMAXRND || input.rte.solver == SOLVER_TWOMAXRND3C ||
        input.rte.solver == SOLVER_DYNAMIC_TWOSTREAM || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) {
      fprintf (stderr,
               " ------------------------------------------------------------------------------------------------------------------"
               "------------------------------------------------\n");
      fprintf (stderr,
               "   lc |   z[km]  |     Rayleigh |      Aerosol              |       Water cloud             |        Ice cloud     "
               "                                  |   Molecular \n");
      fprintf (stderr,
               "      |          |       dtau   |  scatter.      abs.  asy. |    scatter.        abs.  asy. |    scatter.        "
               "abs.  asy.   ff   g1     g2     f   |  absorption \n");
      fprintf (stderr,
               " ------------------------------------------------------------------------------------------------------------------"
               "------------------------------------------------\n");
    } else {
      fprintf (stderr,
               " ------------------------------------------------------------------------------------------------------------------"
               "---------------------------\n");
      fprintf (stderr,
               "   lc |   z[km]  |     Rayleigh |      Aerosol              |       effective cloud (water and ice)                "
               "  |   Molecular | Cloud   \n");
      fprintf (stderr,
               "      |          |       dtau   |  scatter.      abs.  asy. |    scatter.        abs.  asy.   ff    g1   g2      f "
               "  |  absorption | fraction\n");
      fprintf (stderr,
               " ------------------------------------------------------------------------------------------------------------------"
               "---------------------------\n");
    }
  }

  if (first != 0) {
    first = 0;

    /* allocate memory for profiles */

    output->dtauc    = calloc (nlyr, sizeof (float));
    output->dtauc_md = calloc (nlyr, sizeof (float));
    output->ssalb    = calloc (nlyr, sizeof (float));

    output->dtauc_clr = calloc (nlyr, sizeof (float));
    output->ssalb_clr = calloc (nlyr, sizeof (float));

    output->dtauc_cldk = calloc (nlyr, sizeof (float));
    output->ssalb_cldk = calloc (nlyr, sizeof (float));

    output->dtauc_cldn = calloc (nlyr, sizeof (float));
    output->ssalb_cldn = calloc (nlyr, sizeof (float));

    if (input.tipa == TIPA_DIR)
      output->tausol = calloc (nlyr, sizeof (float)); /* ulrike: for tipa dir */
    if (input.raman) {                                //TODO: Delete this, never used
      output->ssalbR  = calloc (nlyr, sizeof (float));
      output->ssalbRL = calloc (nlyr, sizeof (float));
    }

    if (input.rte.solver == SOLVER_MONTECARLO) {
      output->mc.z = calloc (nlev, sizeof (float));

      if (!input.atmosphere3d)
        ASCII_calloc_float_3D (&output->mc.temper,
                               1,
                               1,
                               nlyr + 1); /* CE: temper is defined on levels, so cahnged nlyr -> nlyr+1 ??*/

      /* alternative: output->mc.caoth = calloc (n_caoth, sizeof(mc_caoth_struct)); */
      output->mc.dt = calloc (n_caoth, sizeof (float*));
      output->mc.om = calloc (n_caoth, sizeof (float*));
      output->mc.g1 = calloc (n_caoth, sizeof (float*));
      output->mc.g2 = calloc (n_caoth, sizeof (float*));
      output->mc.ff = calloc (n_caoth, sizeof (float*));
      output->mc.ds = calloc (n_caoth, sizeof (float*));
      output->mc.re = calloc (n_caoth, sizeof (float*));

      for (isp = 0; isp < n_caoth; isp++) {
        output->mc.dt[isp] = calloc (nlev, sizeof (float));
        output->mc.om[isp] = calloc (nlev, sizeof (float));
        output->mc.g1[isp] = calloc (nlev, sizeof (float));
        output->mc.g2[isp] = calloc (nlev, sizeof (float));
        output->mc.ff[isp] = calloc (nlev, sizeof (float));
        output->mc.ds[isp] = calloc (nlev, sizeof (float));
        output->mc.re[isp] = calloc (nlev, sizeof (float));
      }

      output->mc.refind = calloc (nlev, sizeof (float));

      /* SBCA clean this later */
      output->mc.nmomaer = calloc (nlev, sizeof (int));
      output->mc.momaer  = calloc (nlev, sizeof (float**));

      //TODO: size: nlyr , nphamat, ntheta, Speicherzugriffsfehler, wenn nlev
      output->mc.nthetaaer = calloc (nlev, sizeof (int*));
      output->mc.thetaaer  = calloc (nlev, sizeof (float**));
      output->mc.muaer     = calloc (nlev, sizeof (double**));
      output->mc.phaseaer  = calloc (nlev, sizeof (float**));

      for (isp = 0; isp < n_caoth; isp++) {
        for (lc = 0; lc < nlev; lc++) {
          output->mc.ff[isp][lc] = 1.0;
        }
      }
    }

    /* output->atm.nmom+1 stores the maximum number of moments */
    /* for all wavelengths and layers, minus 1                 */

    status = ASCII_calloc_float_3D (&output->pmom, nlyr, nphamat, output->atm.nmom + 1);
    CHKERROUT (status, "failed allocating memory for output->pmom");

    output->pmom01_clr  = calloc (nlyr, sizeof (float));
    output->pmom01_cldk = calloc (nlyr, sizeof (float));
    output->pmom01_cldn = calloc (nlyr, sizeof (float));

    output->ntheta = calloc (nlyr, sizeof (int*));
    output->theta  = calloc (nlyr, sizeof (float**));
    output->mu     = calloc (nlyr, sizeof (double**));
    output->phase  = calloc (nlyr, sizeof (float**));
  } /* end if first */
  else {
    if (!skip_optical_properties) {
      /* need to free mu theta phase */
      for (lc = 0; lc < nlyr; lc++) {
        if (output->ntheta[lc] != NULL) {
          for (ip = 0; ip < nphamat; ip++) {
            if (output->theta[lc][ip] != NULL)
              free (output->theta[lc][ip]);
            if (output->mu[lc][ip] != NULL)
              free (output->mu[lc][ip]);
            if (output->phase[lc][ip] != NULL)
              free (output->phase[lc][ip]);
          }
          free (output->ntheta[lc]);
          free (output->theta[lc]);
          free (output->mu[lc]);
          free (output->phase[lc]);
        }
      }
    }
  }

  if (!skip_optical_properties) {

    /* identify indices of "wcn" and "wck" CAOTHs (thiN and thicK clouds) for twomaxrnd3c */
    if (input.rte.solver == SOLVER_TWOMAXRND3C) {
      i_wcn = -1;
      i_wck = -1;
      i_icn = -1;
      i_ick = -1;

      for (isp = 0; isp < input.n_caoth; isp++) {
        if (strcasecmp (output->caoth[isp].name, "wcn") == 0)
          i_wcn = isp;

        if (strcasecmp (output->caoth[isp].name, "wck") == 0)
          i_wck = isp;

        if (strcasecmp (output->caoth[isp].name, "ick") == 0)
          i_ick = isp;

        if (strcasecmp (output->caoth[isp].name, "icn") == 0)
          i_icn = isp;
      }

      if (i_wcn < 0 || i_wck < 0) {
        fprintf (stderr, "Error. Please define wck and wcn profiles for solver twomaxrnd3c!\n");
        return -1;
      }

      if ((i_icn >= 0 && i_ick < 0) || (i_ick >= 0 && i_icn < 0)) {
        fprintf (stderr, "Error, please define both thin and thick ice clouds, icn and ick!\n");
        return -1;
      }

      if (!input.quiet)
        if (i_icn >= 0 && i_ick >= 0)
          fprintf (stderr, " ... ice clouds for twomaxrnd3c found!\n");

      /* add offset CAOTH_FIR because elements 0 and 1 are reserved for MOL and AER */
      i_wcn += CAOTH_FIR;
      i_wck += CAOTH_FIR;

      if (i_icn > 0)
        i_icn += CAOTH_FIR;

      if (i_ick > 0)
        i_ick += CAOTH_FIR;

      if (input.verbose)
        fprintf (stderr, " ... twomaxrnd3c(): i_wcn=%d, i_wck=%d\n", i_wcn, i_wck);

      if (!input.quiet)
        if (i_icn > 0 && i_ick > 0)
          fprintf (stderr, " ... twomaxrnd3c(): i_icn=%d, i_ick=%d\n", i_icn, i_ick);
    }

    for (isp = 0; isp < n_caoth; isp++) {
      babs_s_int[isp] = 0.0;
      bsca_s_int[isp] = 0.0;
    }

    if (input.write_ext_to_file) {
      strcpy (extfilename, input.rte.mc.filename[FN_MC_BASENAME]);
      strcat (extfilename, ".ext_r");
      if ((extfile = fopen (extfilename, "w")) == NULL)
        return -1;
    }

    for (lc = 0; lc < nlyr; lc++) {

      switch (input.rte.solver) {
      case SOLVER_SSSI:
        /* special treatment for SOLVER_SSSI: caoth single scattering albedo is */
        /* set to 0 because caoth scattering is treated explicitely      */
        /* through the tabulated caoth reflectivity; caoth then only     */
        /* reduce radiance through Lambert-Beer                          */

        for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
          ispo                                    = isp - CAOTH_FIR;
          output->caoth[ispo].optprop.ssa[iv][lc] = 0;
        }

        break;

      default:
        break;
      }

      /* only in the Fu and Liou case, the Rayleigh scattering */
      /* cross section depends on the subband                  */
      switch (input.ck_scheme) {
      case CK_FU:
        bsca_s[CAOTH_MOL] = output->atm.optprop.tau_rayleigh_r[0][0][lc][iv][iq];
        break;

      case CK_KATO:
      case CK_KATO2:
      case CK_KATO2_96:
      case CK_KATO2ANDWANDJI:
      case CK_AVHRR_KRATZ:
      case CK_FILE:
      case CK_LOWTRAN:
      case CK_CRS:
      case CK_REPTRAN:
      case CK_REPTRAN_CHANNEL:
        bsca_s[CAOTH_MOL] = output->atm.optprop.tau_rayleigh_r[0][0][lc][iv][0];
        break;

      case CK_RAMAN:
        if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1))
          bsca_s[CAOTH_MOL] = linpol (wvl1,
                                      wvl2,
                                      output->atm.optprop.tau_rayleigh_r[0][0][lc][iv1r][0],
                                      output->atm.optprop.tau_rayleigh_r[0][0][lc][iv2r][0],
                                      wvl);
        else
          bsca_s[CAOTH_MOL] = output->atm.optprop.tau_rayleigh_r[0][0][lc][iq][0];
        break;

      default:
        fprintf (stderr, "Error: unsupported correlated-k scheme %d\n", input.ck_scheme);
        return -1;

        break;
      }
      switch (output->atm.molabs) {
      case MOLABS_CALC:
      case MOLABS_LOOKUP:
        if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1))
          babs_s[CAOTH_MOL] = linpol (wvl1,
                                      wvl2,
                                      output->atm.optprop.tau_molabs_r[0][0][lc][iv1r][0],
                                      output->atm.optprop.tau_molabs_r[0][0][lc][iv2r][0],
                                      wvl);
        else if (input.raman_fast && ir == 0)
          babs_s[CAOTH_MOL] = output->atm.optprop.tau_molabs_r[0][0][lc][iv][0];
        else
          babs_s[CAOTH_MOL] = output->atm.optprop.tau_molabs_r[0][0][lc][iv][iq];
        if (input.rte.solver == SOLVER_SDISORT)
          babs_mol_md = output->atm.optprop.tau_molabs_md_r[0][0][lc][iv][iq];
        else
          babs_mol_md = 0.0;
        break;

      case MOLABS_FILE_MONO:
        if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1))
          babs_s[CAOTH_MOL] = linpol (wvl1,
                                      wvl2,
                                      output->atm.optprop.tau_molabs_r[0][0][lc][iv1r][0],
                                      output->atm.optprop.tau_molabs_r[0][0][lc][iv2r][0],
                                      wvl);
        else if (input.raman_fast && ir == 0)
          babs_s[CAOTH_MOL] = output->atm.optprop.tau_molabs_r[0][0][lc][iv][0];
        else
          /* babs_s[CAOTH_MOL] = output->atm.optprop.tau_molabs_user[lc]; */
          babs_s[CAOTH_MOL] = output->atm.optprop.tau_molabs_r[0][0][lc][iv][iq];

        babs_mol_md = 0;
        break;

      case MOLABS_FILE_SPEC:
        if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1))
          babs_s[CAOTH_MOL] = linpol (wvl1,
                                      wvl2,
                                      output->atm.optprop.tau_molabs_r[0][0][lc][iv1r][0],
                                      output->atm.optprop.tau_molabs_r[0][0][lc][iv2r][0],
                                      wvl);
        else if (input.raman_fast && ir == 0)
          babs_s[CAOTH_MOL] = output->atm.optprop.tau_molabs_r[0][0][lc][iq][0];
        else
          babs_s[CAOTH_MOL] = output->atm.optprop.tau_molabs_r[0][0][lc][iv][iq];

        babs_mol_md = 0;
        break;

      case MOLABS_NONE:
        babs_s[CAOTH_MOL] = 0;
        babs_mol_md       = 0;
        break;

      default:
        fprintf (stderr, "Error, unknown molecular absorption option %d\n", output->atm.molabs);
        return -1;
      }

      if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1)) {
        ssa_s[CAOTH_AER]  = linpol (wvl1, wvl2, output->aer.optprop.ssa[iv1r][lc], output->aer.optprop.ssa[iv2r][lc], wvl);
        f_s[CAOTH_AER]    = linpol (wvl1, wvl2, output->aer.optprop.f[iv1r][lc], output->aer.optprop.f[iv2r][lc], wvl);
        ff_s[CAOTH_AER]   = linpol (wvl1, wvl2, output->aer.optprop.ff[iv1r][lc], output->aer.optprop.ff[iv2r][lc], wvl);
        g1_s[CAOTH_AER]   = linpol (wvl1, wvl2, output->aer.optprop.g1[iv1r][lc], output->aer.optprop.g1[iv2r][lc], wvl);
        g2_s[CAOTH_AER]   = linpol (wvl1, wvl2, output->aer.optprop.g2[iv1r][lc], output->aer.optprop.g2[iv2r][lc], wvl);
        dtau_s[CAOTH_AER] = linpol (wvl1, wvl2, output->aer.optprop.dtau[iv1r][lc], output->aer.optprop.dtau[iv2r][lc], wvl);

        for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
          ispo = isp - CAOTH_FIR;
          ssa_s[isp] =
            linpol (wvl1, wvl2, output->caoth[ispo].optprop.ssa[iv1r][lc], output->caoth[ispo].optprop.ssa[iv2r][lc], wvl);
          f_s[isp]  = linpol (wvl1, wvl2, output->caoth[ispo].optprop.f[iv1r][lc], output->caoth[ispo].optprop.f[iv2r][lc], wvl);
          ff_s[isp] = linpol (wvl1, wvl2, output->caoth[ispo].optprop.ff[iv1r][lc], output->caoth[ispo].optprop.ff[iv2r][lc], wvl);
          g1_s[isp] = linpol (wvl1, wvl2, output->caoth[ispo].optprop.g1[iv1r][lc], output->caoth[ispo].optprop.g1[iv2r][lc], wvl);
          g2_s[isp] = linpol (wvl1, wvl2, output->caoth[ispo].optprop.g2[iv1r][lc], output->caoth[ispo].optprop.g2[iv2r][lc], wvl);
          dtau_s[isp] =
            linpol (wvl1, wvl2, output->caoth[ispo].optprop.dtau[iv1r][lc], output->caoth[ispo].optprop.dtau[iv2r][lc], wvl);
          dscale_s[isp] =
            linpol (wvl1, wvl2, output->caoth[ispo].optprop.dscale[iv1r][lc], output->caoth[ispo].optprop.dscale[iv2r][lc], wvl);
        }

        //20120816ak refind is not used, commented
        //	refind   = linpol( wvl1, wvl2, output->atm.microphys.refind[iv1r][lc],
        //			   output->atm.microphys.refind[iv2r][lc], wvl);
      } else {
        ssa_s[CAOTH_AER]  = output->aer.optprop.ssa[iv][lc];
        f_s[CAOTH_AER]    = output->aer.optprop.f[iv][lc];
        ff_s[CAOTH_AER]   = output->aer.optprop.ff[iv][lc];
        g1_s[CAOTH_AER]   = output->aer.optprop.g1[iv][lc];
        g2_s[CAOTH_AER]   = output->aer.optprop.g2[iv][lc];
        dtau_s[CAOTH_AER] = output->aer.optprop.dtau[iv][lc];

        for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
          ispo = isp - CAOTH_FIR;
          if (output->caoth[ispo].optprop.ssa != NULL) {
            ssa_s[isp]    = output->caoth[ispo].optprop.ssa[iv][lc];
            f_s[isp]      = output->caoth[ispo].optprop.f[iv][lc];
            ff_s[isp]     = output->caoth[ispo].optprop.ff[iv][lc];
            g1_s[isp]     = output->caoth[ispo].optprop.g1[iv][lc];
            g2_s[isp]     = output->caoth[ispo].optprop.g2[iv][lc];
            dtau_s[isp]   = output->caoth[ispo].optprop.dtau[iv][lc];
            dscale_s[isp] = output->caoth[ispo].optprop.dscale[iv][lc];
            /* ?????????? */
          }
        }

        //20120816ak refind is not used, commented
        //	refind   = output->atm.microphys.refind[iv][lc];
      }

      for (isp = CAOTH_FIR; isp < n_caoth; isp++)
        ssa_s_unsc[isp] = ssa_s[isp] / (1.0 + f_s[isp] * (ssa_s[isp] - 1.0));

      /* absorption by aerosols and caoths */
      if (input.absorption) {
        for (isp = CAOTH_AER; isp < n_caoth; isp++)
          babs_s[isp] = (1.0 - ssa_s[isp]) * dtau_s[isp];
        for (isp = CAOTH_FIR; isp < n_caoth; isp++)
          babs_s_unsc[isp] = (1.0 - ssa_s_unsc[isp]) * dtau_s[isp] / (1.0 - ssa_s_unsc[isp] * f_s[isp]);

        babs_tot = 0.0;
        for (isp = 0; isp < n_caoth; isp++)
          babs_tot += babs_s[isp];

        babs_tot_md = babs_tot - babs_s[CAOTH_MOL] + babs_mol_md;

        if (babs_tot > FLT_MAX)
          babs_tot = FLT_MAX;
        if (babs_tot_md > FLT_MAX)
          babs_tot_md = FLT_MAX;
      } else {
        for (isp = 0; isp < n_caoth; isp++)
          babs_s[isp] = 0.0;

        babs_tot    = 0.0;
        babs_tot_md = 0.0;
      }

      /* scattering by non-molecules */
      for (isp = CAOTH_AER; isp < n_caoth; isp++)
        bsca_s[isp] = ssa_s[isp] * dtau_s[isp];

      for (isp = CAOTH_FIR; isp < n_caoth; isp++)
        bsca_s_unsc[isp] = ssa_s_unsc[isp] * dtau_s[isp] / (1.0 - ssa_s_unsc[isp] * f_s[isp]);

      /* switch scattering off, if user wants so */
      if (input.aer.no_scattering) {
        bsca_s[CAOTH_AER]      = 0.0;
        bsca_s_unsc[CAOTH_AER] = 0.0;
        ssa_s[CAOTH_AER]       = 0.0;
      }

      for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
        ispo = isp - CAOTH_FIR;
        if (input.caoth[ispo].no_scattering) {
          bsca_s[isp]      = 0.0;
          bsca_s_unsc[isp] = 0.0;
          ssa_s[isp]       = 0.0;
        }
      }

      /* Scattering coefficient */
      if (input.scattering) {
        bsca_tot = 0.0;
        for (isp = 0; isp < n_caoth; isp++)
          bsca_tot += bsca_s[isp];
        if (bsca_tot > FLT_MAX)
          bsca_tot = FLT_MAX;
        if (bsca_tot < FLT_MIN)
          bsca_tot = FLT_MIN;
      } else {
        for (isp = 0; isp < n_caoth; isp++) {
          bsca_s[isp]      = 0.0;
          bsca_s_unsc[isp] = 0.0;
          ssa_s[isp]       = 0.0;
        }
        bsca_tot = FLT_MIN;
      }

      /* Extinction coefficient */
      bext_tot = babs_tot + bsca_tot;
      if (bext_tot > FLT_MAX)
        bext_tot = FLT_MAX;
      if (bext_tot < FLT_MIN)
        bext_tot = FLT_MIN;

      bext_tot_md = babs_tot_md + bsca_tot;
      if (bext_tot_md > FLT_MAX)
        bext_tot_md = FLT_MAX;
      if (bext_tot_md < FLT_MIN)
        bext_tot_md = FLT_MIN;

      for (isp = CAOTH_FIR; isp < n_caoth; isp++)
        gg_s_unsc[isp] = g1_s[isp] * (1.0 - f_s[isp]) + f_s[isp];

      if (input.write_ext_to_file) {
        if (verbose && lc == 0)
          printf ("...saving Extinction Data\n");

        /* ??? do we need this here: */
        /*        if (output->cf.nlev == 0) { */

        /*	fprintf (extfile, "%5d %8.4f %9.6f %9.6f %9.6f %5.3f %11.6f %11.6f %5.3f %11.6f %11.6f %5.3f %5.3f %5.3f %6.3f %5.3f %11.6f\n", */
        fprintf (extfile,
                 "%5d %8.4f %11.6e %11.6e %11.6e %5.3f %11.6e %11.6e %5.3f %11.6e %11.6e %5.3f %5.3f %5.3f %6.3f %5.3f %11.6e\n",
                 lc,
                 output->atm.zd[lc + 1] + output->alt.altitude,
                 /* Rayleigh */
                 bsca_s[CAOTH_MOL],
                 /* Aerosol, averaging of delta scaling factors for mixture of aerosol types not yet implemented !!!*/
                 bsca_s[CAOTH_AER],
                 babs_s[CAOTH_AER],
                 g1_s[CAOTH_AER] * ff_s[CAOTH_AER] + (1.0 - ff_s[CAOTH_AER]) * g2_s[CAOTH_AER],
                 /* Water cloud */
                 bsca_s_unsc[i_wc],
                 babs_s_unsc[i_wc],
                 gg_s_unsc[i_wc] * ff_s[i_wc] + (1.0 - ff_s[i_wc]) * g2_s[i_wc],
                 /* ice cloud */
                 bsca_s_unsc[i_ic],
                 babs_s_unsc[i_ic],
                 gg_s_unsc[i_ic] * ff_s[i_ic] + (1.0 - ff_s[i_ic]) * g2_s[i_ic],
                 ff_s[i_ic],
                 gg_s_unsc[i_ic],
                 g2_s[i_ic],
                 f_s[i_ic],
                 /* Molecules */
                 babs_s[CAOTH_MOL]);
      }

      /* CE: Modified verbose output. Un-deltascaled optical properties are printed.*/
      if (verbose) {
        if (output->cf.nlev == 0 || input.rte.solver == SOLVER_TWOMAXRND || input.rte.solver == SOLVER_TWOMAXRND3C ||
            input.rte.solver == SOLVER_DYNAMIC_TWOSTREAM || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) {
          fprintf (stderr,
                   "%5d | %8.4f | %9.6e | %9.6f %9.6f %5.3f | %11.6f %11.6f %5.3f | %11.6f %11.6f %5.3f %5.3f %5.3f %6.3f %5.3f | "
                   "%11.6e\n",
                   lc,
                   output->atm.zd[lc + 1] + output->alt.altitude,
                   /* Rayleigh */
                   bsca_s[CAOTH_MOL],
                   /* Aerosol, averaging of delta scaling factors for mixture of aerosol types not yet implemented !!!*/
                   bsca_s[CAOTH_AER],
                   babs_s[CAOTH_AER],
                   g1_s[CAOTH_AER] * ff_s[CAOTH_AER] + (1.0 - ff_s[CAOTH_AER]) * g2_s[CAOTH_AER],
                   /* Water cloud */
                   bsca_s_unsc[i_wc],
                   babs_s_unsc[i_wc],
                   gg_s_unsc[i_wc] * ff_s[i_wc] + (1.0 - ff_s[i_wc]) * g2_s[i_wc],
                   /* ice cloud */
                   bsca_s_unsc[i_ic],
                   babs_s_unsc[i_ic],
                   gg_s_unsc[i_ic] * ff_s[i_ic] + (1.0 - ff_s[i_ic]) * g2_s[i_ic],
                   ff_s[i_ic],
                   gg_s_unsc[i_ic],
                   g2_s[i_ic],
                   f_s[i_ic],
                   /* Molecules */
                   babs_s[CAOTH_MOL]);
        } else {
          fprintf (stderr,
                   "%5d | %8.4f | %9.6e | %9.6f %9.6f %5.3f | %11.6f %11.6f %5.3f %5.3f %5.3f %6.3f %5.3f | %11.6e | %.3f\n",
                   lc,
                   output->atm.zd[lc + 1] + output->alt.altitude,
                   /* Rayleigh */
                   bsca_s[CAOTH_MOL],
                   /* Aerosol, averaging of delta scaling factors for mixture of aerosol types not yet implemented !!!*/
                   bsca_s[CAOTH_AER],
                   babs_s[CAOTH_AER],
                   g1_s[CAOTH_AER] * ff_s[CAOTH_AER] + (1.0 - ff_s[CAOTH_AER]) * g2_s[CAOTH_AER],
                   /* effective cloud, both  */
                   bsca_s_unsc[i_wc],
                   babs_s_unsc[i_wc],
                   gg_s_unsc[i_wc] * ff_s[i_wc] + (1.0 - ff_s[i_wc]) * g2_s[i_wc],
                   ff_s[i_wc],
                   gg_s_unsc[i_ic],
                   g2_s[i_wc],
                   f_s[i_ic], /* question to CE: does this make sense??? mixing wc and ic properties! */
                              /* Molecules, Cloudfraction */
                   babs_s[CAOTH_MOL],
                   output->cf.cf[lc]);
        }
      }

      for (isp = 0; isp < n_caoth; isp++) {
        mom_s_g1[isp] = bsca_s[isp];
        mom_s_g2[isp] = bsca_s[isp];
      }

      output->pmom[lc][0][0] = 1.0;

      /* check that raman interpolation can work */
      if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1)) {
        if (output->aer.optprop.nmom[iv1r][lc] != output->aer.optprop.nmom[iv2r][lc]) {
          CHKERROUT (-1, "Number of aerosol moments must be equal for all wavelengths when Raman scattering is included");
        }
        for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
          ispo = isp - CAOTH_FIR;
          if (output->caoth[ispo].optprop.nmom[iv1][lc] != output->caoth[ispo].optprop.nmom[iv2][lc]) {
            status = -1;
            fprintf (stderr, "Number of %s moments must be equal for all wavelengths\n", output->caoth[ispo].fullname);
            fprintf (stderr, "when Raman scattering is included.\n");
            CHKERR (status);
          }
        }
      }

      /* zeroth moment needed for normalization */
      if (output->aer.optprop.nmom[iv][lc] > 0) {
        if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1))
          mom0_s[CAOTH_AER] =
            linpol (wvl1, wvl2, output->aer.optprop.moment[iv1][lc][0][0], output->aer.optprop.moment[iv2][lc][0][0], wvl);
        else
          mom0_s[CAOTH_AER] = output->aer.optprop.moment[iv][lc][0][0];
      }

      for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
        ispo = isp - CAOTH_FIR;
        if (output->caoth[ispo].optprop.nmom != NULL) { /* rather ask if !montecarlo,
							   or if mom0_s needed, SBCA */
          if (output->caoth[ispo].optprop.nmom[iv][lc] > 0) {
            if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1))
              mom0_s[isp] = linpol (wvl1,
                                    wvl2,
                                    output->caoth[ispo].optprop.moment[iv1][lc][0][0],
                                    output->caoth[ispo].optprop.moment[iv2][lc][0][0],
                                    wvl);
            else
              mom0_s[isp] = output->caoth[ispo].optprop.moment[iv][lc][0][0];
          }
        }
      }

      for (k = 1; k <= output->atm.nmom; k++) {
        /* aerosol */
        if (output->aer.optprop.nmom[iv][lc] > 0) {
          if (k < output->aer.optprop.nmom[iv][lc]) {
            if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1))
              momk =
                linpol (wvl1, wvl2, output->aer.optprop.moment[iv1r][lc][0][k], output->aer.optprop.moment[iv2r][lc][0][k], wvl);
            else
              momk = output->aer.optprop.moment[iv][lc][0][k];

            mom_s[CAOTH_AER][0] = bsca_s[CAOTH_AER] * momk / mom0_s[CAOTH_AER];
          } else {
            mom_s[CAOTH_AER][0] = 0;
          }
        } else { /* double Henyey-Greenstein */
          mom_s_g1[CAOTH_AER] *= g1_s[CAOTH_AER];
          mom_s_g2[CAOTH_AER] *= g2_s[CAOTH_AER];

          mom_s[CAOTH_AER][0] = (mom_s_g1[CAOTH_AER] * ff_s[CAOTH_AER] + (1.0 - ff_s[CAOTH_AER]) * mom_s_g2[CAOTH_AER]);
        }

        /* caoth */
        for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
          ispo = isp - CAOTH_FIR;
          if (output->caoth[ispo].optprop.nmom != NULL) { /* dito, SBCA */
            if (output->caoth[ispo].optprop.nmom[iv][lc] > 0) {
              if (k < output->caoth[ispo].optprop.nmom[iv][lc]) {
                if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1))
                  momk = linpol (wvl1,
                                 wvl2,
                                 output->caoth[ispo].optprop.moment[iv1r][lc][0][k],
                                 output->caoth[ispo].optprop.moment[iv2r][lc][0][k],
                                 wvl);
                else
                  momk = output->caoth[ispo].optprop.moment[iv][lc][0][k];

                mom_s[isp][0] = bsca_s[isp] * momk / mom0_s[isp];
              } else {
                mom_s[isp][0] = 0;
              }
            } else { /* double Henyey-Greenstein */
              mom_s_g1[isp] *= g1_s[isp];
              mom_s_g2[isp] *= g2_s[isp];

              mom_s[isp][0] = (mom_s_g1[isp] * ff_s[isp] + (1.0 - ff_s[isp]) * mom_s_g2[isp]);
            }
          }
        }

        if (k == 2) /* Rayleigh scattering, including depolarization */
          mom_s[CAOTH_MOL][0] = bsca_s[CAOTH_MOL] * rayleigh_mom2;
        else
          mom_s[CAOTH_MOL][0] = 0.0;

        /* sum all moments, and normalize with bsca_tot */
        output->pmom[lc][0][k] = 0.0;
        if (k == 1) {
          output->pmom01_clr[lc]  = 0.0;
          output->pmom01_cldk[lc] = 0.0;
          output->pmom01_cldn[lc] = 0.0;
        }

        for (isp = 0; isp < n_caoth; isp++) {
          output->pmom[lc][0][k] += mom_s[isp][0];

          if (input.rte.solver == SOLVER_TWOMAXRND || input.rte.solver == SOLVER_DYNAMIC_TWOSTREAM ||
              input.rte.solver == SOLVER_DYNAMIC_TENSTREAM)
            if (isp != i_wc && isp != i_ic && k == 1)
              output->pmom01_clr[lc] += mom_s[isp][0];

          if (input.rte.solver == SOLVER_TWOMAXRND3C) {

            // clear := everything except thin and thick water and ice clouds
            if (k == 1)
              if (isp != i_wcn && isp != i_wck && isp != i_icn && isp != i_ick)
                output->pmom01_clr[lc] += mom_s[isp][0];

            // wck := everything except wcn and icn
            if (k == 1)
              if (isp != i_wcn && isp != i_icn)
                output->pmom01_cldk[lc] += mom_s[isp][0];

            // wcn := everything except wck and ick
            if (k == 1)
              if (isp != i_wck && isp != i_ick)
                output->pmom01_cldn[lc] += mom_s[isp][0];
          }
        }

        output->pmom[lc][0][k] /= bsca_tot;

        if (input.rte.solver == SOLVER_TWOMAXRND || input.rte.solver == SOLVER_DYNAMIC_TWOSTREAM ||
            input.rte.solver == SOLVER_DYNAMIC_TENSTREAM)
          if (k == 1) {
            if (bsca_tot - bsca_s[i_wc] - bsca_s[i_ic] > 0)
              output->pmom01_clr[lc] /= (bsca_tot - bsca_s[i_wc] - bsca_s[i_ic]);
            else
              output->pmom01_clr[lc] = 0.0;
          }

        if (input.rte.solver == SOLVER_TWOMAXRND3C)
          if (k == 1) {
            double bsca_clr  = bsca_tot - bsca_s[i_wck] - bsca_s[i_wcn];
            double bsca_cldk = bsca_tot - bsca_s[i_wcn];
            double bsca_cldn = bsca_tot - bsca_s[i_wck];

            if (i_icn > 0 && i_ick > 0) {
              bsca_clr -= (bsca_s[i_ick] - bsca_s[i_icn]);
              bsca_cldn -= bsca_s[i_ick];
              bsca_cldk -= bsca_s[i_icn];
            }

            if (bsca_clr > 0)
              output->pmom01_clr[lc] /= bsca_clr;
            else
              output->pmom01_clr[lc] = 0.0;

            if (bsca_cldk > 0)
              output->pmom01_cldk[lc] /= bsca_cldk;
            else
              output->pmom01_cldk[lc] = 0.0;

            if (bsca_cldn > 0)
              output->pmom01_cldn[lc] /= bsca_cldn;
            else
              output->pmom01_cldn[lc] = 0.0;
          }

      } /* end loop k */

      /* phase functions */
      /* this should NOT be performed in case of MYSTIC solver!
	 Basically, this should ONLY be done in case of
	 SOLVER_FDISORT2/SOLVER_DISORT with new ICM, or with SOLVER_SSLIDAR */
      if (((input.rte.solver == SOLVER_FDISORT2 || input.rte.solver == SOLVER_DISORT) &&
           input.rte.disort_icm == DISORT_ICM_PHASE) ||
          input.rte.solver == SOLVER_SSLIDAR) {
        /* the following has strong resemblance with the function
	   sort_and_add_weighted_phase in phasetable.c */

        n_in = 99; /* maximal possible number of n_ins to be defined */

        tmp_ntheta_in         = calloc (n_in, sizeof (int));
        tmp_theta_in          = calloc (n_in, sizeof (float*));
        tmp_mu_in             = calloc (n_in, sizeof (double*));
        tmp_phase_in          = calloc (n_in, sizeof (float*));
        phase_calloced        = calloc (n_in, sizeof (int));
        interp_optprop_weight = calloc (n_in, sizeof (double));

        output->ntheta[lc] = calloc (nphamat, sizeof (int));
        output->theta[lc]  = calloc (nphamat, sizeof (float*));
        output->mu[lc]     = calloc (nphamat, sizeof (double*));
        output->phase[lc]  = calloc (nphamat, sizeof (float*));

        /* for sslidar, we want this not only for P_11, but also for
	   P_12 and P_22 */
        for (ip = 0; ip < nphamat; ip++) {

          /* find out number of phase functions already existent */
          n_in  = 0;
          n_tot = 0;
          if (bsca_s[CAOTH_AER] > 0.0) {
            ips = ip_act[CAOTH_AER][ip];

            if (output->aer.optprop.ntheta[iv][lc][ips] > 0) {
              if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1)) {
                /* for raman, FIXME, THIS MAY NEEDED TO BE IMPROVED AKY20111605 */
                tmp_ntheta_in[n_in]         = output->aer.optprop.ntheta[iv1r][lc][ips];
                tmp_theta_in[n_in]          = output->aer.optprop.theta[iv1r][lc][ips];
                tmp_mu_in[n_in]             = output->aer.optprop.mu[iv1r][lc][ips];
                tmp_phase_in[n_in]          = output->aer.optprop.phase[iv1r][lc][ips];
                interp_optprop_weight[n_in] = bsca_s[CAOTH_AER] / bsca_tot * (wvl2 - wvl) / (wvl2 - wvl1);
                n_tot += tmp_ntheta_in[n_in];
                n_in++;
              } else {
                tmp_ntheta_in[n_in]         = output->aer.optprop.ntheta[iv][lc][ips];
                tmp_theta_in[n_in]          = output->aer.optprop.theta[iv][lc][ips];
                tmp_mu_in[n_in]             = output->aer.optprop.mu[iv][lc][ips];
                tmp_phase_in[n_in]          = output->aer.optprop.phase[iv][lc][ips];
                interp_optprop_weight[n_in] = bsca_s[CAOTH_AER] / bsca_tot;
                n_tot += tmp_ntheta_in[n_in];
                n_in++;
              }
            } else {
              if (output->aer.optprop.nmom[iv][lc] > 0) {
                /* only moments defined */
                fprintf (stderr,
                         "Error, you need to specify 'disort_intcor moments' in order to use these aerosol phase functions !\n");
                return -1;
              } else {
                if (ip > 0) {
                  fprintf (stderr,
                           "Error, you are trying to combine aerosol Henyey-Greenstein with polarisation! This does not work!\n");
                  return -1;
                }

                /* HG */
                status = create_phase_from_HG (g1_s[CAOTH_AER],
                                               &(tmp_ntheta_in[n_in]),
                                               &(tmp_theta_in[n_in]),
                                               &(tmp_mu_in[n_in]),
                                               &(tmp_phase_in[n_in]),
                                               input.rte.solver == SOLVER_SSLIDAR);
                CHKERR (status);

                phase_calloced[n_in] = 1;

                interp_optprop_weight[n_in] = bsca_s[CAOTH_AER] / bsca_tot * ff_s[CAOTH_AER];
                n_tot += tmp_ntheta_in[n_in];
                n_in++;

                /* double HG */
                if (ff_s[CAOTH_AER] < 1.) {
                  status = create_phase_from_HG (g2_s[CAOTH_AER],
                                                 &(tmp_ntheta_in[n_in]),
                                                 &(tmp_theta_in[n_in]),
                                                 &(tmp_mu_in[n_in]),
                                                 &(tmp_phase_in[n_in]),
                                                 input.rte.solver == SOLVER_SSLIDAR);
                  CHKERR (status);

                  phase_calloced[n_in] = 1;

                  interp_optprop_weight[n_in] = bsca_s[CAOTH_AER] / bsca_tot * (1. - ff_s[CAOTH_AER]);
                  n_tot += tmp_ntheta_in[n_in];
                  n_in++;
                }
              }
            }
          } /* if bsca_s[CAOTH_AER] > 0 */

          for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
            ispo = isp - CAOTH_FIR;
            if (bsca_s[isp] > 0.0) {
              ips = ip_act[isp][ip];
              if (output->caoth[ispo].optprop.ntheta[iv][lc][ips] > 0) {
                if ((input.raman && !input.raman_fast) || (input.raman_fast && ir == 1)) {
                  /* for raman, FIXME, THIS MAY NEEDED TO BE IMPROVED AKY20111605 */
                  tmp_ntheta_in[n_in]         = output->caoth[ispo].optprop.ntheta[iv1r][lc][ips];
                  tmp_theta_in[n_in]          = output->caoth[ispo].optprop.theta[iv1r][lc][ips];
                  tmp_mu_in[n_in]             = output->caoth[ispo].optprop.mu[iv1r][lc][ips];
                  tmp_phase_in[n_in]          = output->caoth[ispo].optprop.phase[iv1r][lc][ips];
                  interp_optprop_weight[n_in] = bsca_s[isp] / bsca_tot * (wvl2 - wvl) / (wvl2 - wvl1);
                  n_tot += tmp_ntheta_in[n_in];
                  n_in++;
                } else {
                  tmp_ntheta_in[n_in]         = output->caoth[ispo].optprop.ntheta[iv][lc][ips];
                  tmp_theta_in[n_in]          = output->caoth[ispo].optprop.theta[iv][lc][ips];
                  tmp_mu_in[n_in]             = output->caoth[ispo].optprop.mu[iv][lc][ips];
                  tmp_phase_in[n_in]          = output->caoth[ispo].optprop.phase[iv][lc][ips];
                  interp_optprop_weight[n_in] = bsca_s[isp] / bsca_tot;
                  n_tot += tmp_ntheta_in[n_in];
                  n_in++;
                }
              } else {
                if (output->caoth[ispo].optprop.nmom[iv][lc] > 0) {
                  /* only moments defined */
                  fprintf (stderr,
                           "Error, you need to specify 'disort_intcor moments' in order to use these %s phase functions !\n",
                           output->caoth[ispo].fullname);
                  return -1;
                } else {
                  if (ip > 0) {
                    fprintf (stderr,
                             "Error, you are trying to combine %s Henyey-Greenstein with polarisation! This does not work!\n",
                             output->caoth[ispo].fullname);
                    return -1;
                  }

                  /* HG */
                  status = create_phase_from_HG (g1_s[isp],
                                                 &(tmp_ntheta_in[n_in]),
                                                 &(tmp_theta_in[n_in]),
                                                 &(tmp_mu_in[n_in]),
                                                 &(tmp_phase_in[n_in]),
                                                 input.rte.solver == SOLVER_SSLIDAR);
                  CHKERR (status);

                  phase_calloced[n_in] = 1;

                  interp_optprop_weight[n_in] = bsca_s[isp] / bsca_tot * ff_s[isp];
                  n_tot += tmp_ntheta_in[n_in];
                  n_in++;

                  /* double HG */
                  if (ff_s[isp] < 1.) {
                    status = create_phase_from_HG (g2_s[isp],
                                                   &(tmp_ntheta_in[n_in]),
                                                   &(tmp_theta_in[n_in]),
                                                   &(tmp_mu_in[n_in]),
                                                   &(tmp_phase_in[n_in]),
                                                   input.rte.solver == SOLVER_SSLIDAR);
                    CHKERR (status);

                    phase_calloced[n_in] = 1;

                    interp_optprop_weight[n_in] = bsca_s[isp] / bsca_tot * (1. - ff_s[isp]);
                    n_tot += tmp_ntheta_in[n_in];
                    n_in++;
                  }
                }
              }
            } /* end if (bsca_s[isp] > 0.0) */
          }   /* end for isp */

          /* rayleigh */
          if (bsca_s[CAOTH_MOL] > 0.0) {
            status = create_phase_from_Rayleigh (rayleigh_depol,
                                                 &(tmp_ntheta_in[n_in]),
                                                 &(tmp_theta_in[n_in]),
                                                 &(tmp_mu_in[n_in]),
                                                 &(tmp_phase_in[n_in]),
                                                 input.rte.solver == SOLVER_SSLIDAR,
                                                 ip);
            CHKERR (status);
            phase_calloced[n_in] = 1;

            interp_optprop_weight[n_in] = bsca_s[CAOTH_MOL] / bsca_tot;
            n_tot += tmp_ntheta_in[n_in];
            n_in++;
          }

          if (input.rte.solver == SOLVER_SSLIDAR) {

            /* allocate theta dimension */
            output->theta[lc][ip] = calloc (1, sizeof (float));
            output->mu[lc][ip]    = calloc (1, sizeof (double));
            output->phase[lc][ip] = calloc (1, sizeof (float));

            /* set trivial values for "phase function" */
            output->ntheta[lc][ip]   = 1;
            output->theta[lc][ip][0] = 0.0;
            output->mu[lc][ip][0]    = -1.0;

            /* interpolate backscatter direction */
            output->phase[lc][ip][0] = 0.0;
            for (i = 0; i < n_in; i++)
              output->phase[lc][ip][0] += interp_optprop_weight[i] * tmp_phase_in[i][0];

          } else {
            tmp_theta_new = calloc (n_tot, sizeof (float));
            tmp_mu_new    = calloc (n_tot, sizeof (double));

            output->ntheta[lc][ip] = sort_theta_and_mu (n_in, tmp_ntheta_in, tmp_theta_in, tmp_theta_new, tmp_mu_in, tmp_mu_new);

            if (output->ntheta[lc][ip] == -1)
              return fct_err_out (-1, "sort_theta_and_mu", ERROR_POSITION);

            /* allocate theta dimension */
            output->theta[lc][ip] = calloc (output->ntheta[lc][ip], sizeof (float));
            output->mu[lc][ip]    = calloc (output->ntheta[lc][ip], sizeof (double));
            output->phase[lc][ip] = calloc (output->ntheta[lc][ip], sizeof (float));

            /* copy new theta grid into target */
            for (i = 0; i < output->ntheta[lc][ip]; i++) {
              output->theta[lc][ip][i] = tmp_theta_new[i];
              output->mu[lc][ip][i]    = tmp_mu_new[i];
            }

            /* interpolate phase */
            status = interpolate_phase_weighted (n_in,
                                                 tmp_ntheta_in,
                                                 tmp_mu_in,
                                                 tmp_phase_in,
                                                 interp_optprop_weight,
                                                 output->ntheta[lc][ip],
                                                 output->mu[lc][ip],
                                                 output->phase[lc][ip],
                                                 0);
            CHKERR (status);

            free (tmp_theta_new);
            free (tmp_mu_new);

            for (i = 0; i < n_in; i++) {
              if (phase_calloced[i]) {
                free (tmp_theta_in[i]);
                free (tmp_mu_in[i]);
                free (tmp_phase_in[i]);
              }
            }
          } /* end else (SOLVER_SSLIDAR) */
        }   /* end for ip */

        free (tmp_ntheta_in);
        free (tmp_theta_in);
        free (tmp_mu_in);
        free (tmp_phase_in);
        free (interp_optprop_weight);
        free (phase_calloced);

        if (input.rte.solver != SOLVER_SSLIDAR)
          normalize_phase (output->mu[lc], output->phase[lc], NULL, output->ntheta[lc], 1, 0);

      } /* end (disort && disort_icm phase) || sslidar */

      /* total column for verbose output */
      /* CE: print un-deltascaled optical thickness.
         For Aerosol the averaging of the delta scaling
         factor is not yet implemented !!! */
      for (isp = 0; isp < n_caoth; isp++) {
        babs_s_int[isp] += babs_s[isp];
        bsca_s_int[isp] += bsca_s[isp];
      }

      for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
        babs_s_unsc_int[isp] += babs_s_unsc[isp];
        bsca_s_unsc_int[isp] += bsca_s_unsc[isp];
      }

      /* Optical properties for the mc model */
      if (input.rte.solver == SOLVER_MONTECARLO) {
        for (isp = 0; isp < n_caoth; isp++)
          output->mc.dt[isp][nlyr - 1 - lc] = bsca_s[isp] + babs_s[isp];

        if (bsca_s[CAOTH_MOL] + babs_s[CAOTH_MOL] > 0.0)
          output->mc.om[CAOTH_MOL][nlyr - 1 - lc] = bsca_s[CAOTH_MOL] / (bsca_s[CAOTH_MOL] + babs_s[CAOTH_MOL]);
        else
          output->mc.om[CAOTH_MOL][nlyr - 1 - lc] = 0.0;

        /* in case absorption is turned off, molecular absorption is
	   turned off somewhere else. Look for input.molabs */
        if (input.absorption)
          for (isp = CAOTH_AER; isp < n_caoth; isp++)
            output->mc.om[isp][nlyr - 1 - lc] = ssa_s[isp];
        else
          for (isp = CAOTH_AER; isp < n_caoth; isp++)
            output->mc.om[isp][nlyr - 1 - lc] = 1.0;

        for (isp = CAOTH_AER; isp < n_caoth; isp++) {
          output->mc.g1[isp][nlyr - 1 - lc] = g1_s[isp];
          output->mc.g2[isp][nlyr - 1 - lc] = g2_s[isp];
          output->mc.ff[isp][nlyr - 1 - lc] = ff_s[isp];
        }

        for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
          ispo = isp - CAOTH_FIR;
          /* effective droplet radius */
          output->mc.ds[isp][nlyr - 1 - lc] = dscale_s[isp];
          if (output->caoth[ispo].microphys.effr_layer != NULL)
            output->mc.re[isp][nlyr - 1 - lc] = output->caoth[ispo].microphys.effr_layer[lc];
        }

        /* ??? consequently, one should apply input.atm.interpol_method_refind; for reasons */
        /* ??? of lazyness we simply average the refractive index at the adjacent levels    */
        /* ??? to get the layer property for MYSTIC                                         */
        output->mc.refind[nlyr - 1 - lc] = 0.5 * (output->atm.microphys.refind[iv][lc] + output->atm.microphys.refind[iv][lc + 1]);

        if (input.rte.mc.spectral_is) {
          for (isp = 0; isp < n_caoth; isp++) {
            output->mc.alis.dt[iv][isp][nlyr - 1 - lc] = output->mc.dt[isp][nlyr - 1 - lc];
            output->mc.alis.om[iv][isp][nlyr - 1 - lc] = output->mc.om[isp][nlyr - 1 - lc];
          }
        }
      }

      /* Set optical depth and single scattering albedo */
      output->dtauc[lc]    = bext_tot;
      output->dtauc_md[lc] = bext_tot_md;
      output->ssalb[lc]    = bsca_tot / bext_tot;

      if (input.rte.solver == SOLVER_TWOMAXRND || input.rte.solver == SOLVER_DYNAMIC_TWOSTREAM ||
          input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) {
        output->dtauc_clr[lc] = bext_tot - bsca_s[i_wc] - babs_s[i_wc] - bsca_s[i_ic] - babs_s[i_ic];

        if (output->dtauc_clr[lc] > 0)
          output->ssalb_clr[lc] = (bsca_tot - bsca_s[i_wc] - bsca_s[i_ic]) / output->dtauc_clr[lc];
        else
          output->ssalb_clr[lc] = 1.0;
      }

      if (input.rte.solver == SOLVER_TWOMAXRND3C) {
        double bsca_clr  = 0;
        double bsca_cldn = 0;
        double bsca_cldk = 0;

        // clear := everything except thin and thick water and ice clouds
        output->dtauc_clr[lc] = bext_tot - bsca_s[i_wcn] - babs_s[i_wcn] - bsca_s[i_wck] - babs_s[i_wck];
        if (i_icn > 0 && i_ick > 0)
          output->dtauc_clr[lc] -= (bsca_s[i_icn] + babs_s[i_icn] + bsca_s[i_ick] + babs_s[i_ick]);

        bsca_clr = bsca_tot - bsca_s[i_wcn] - bsca_s[i_wck];
        if (i_icn > 0 && i_ick > 0)
          bsca_clr -= (bsca_s[i_icn] + bsca_s[i_ick]);

        if (output->dtauc_clr[lc] > 0)
          output->ssalb_clr[lc] = bsca_clr / output->dtauc_clr[lc];
        else
          output->ssalb_clr[lc] = 1.0;

        // wck := everything except thin water and ice clouds
        output->dtauc_cldk[lc] = bext_tot - bsca_s[i_wcn] - babs_s[i_wcn];
        if (i_icn > 0 && i_ick > 0)
          output->dtauc_cldk[lc] -= (bsca_s[i_icn] - babs_s[i_icn]);

        bsca_cldk = bsca_tot - bsca_s[i_wcn];
        if (i_icn > 0 && i_ick > 0)
          bsca_cldk -= bsca_s[i_icn];

        if (output->dtauc_cldk[lc] > 0)
          output->ssalb_cldk[lc] = bsca_cldk / output->dtauc_cldk[lc];
        else
          output->ssalb_cldk[lc] = 1.0;

        // wcn := everything except thick water and ice clouds
        output->dtauc_cldn[lc] = bext_tot - bsca_s[i_wck] - babs_s[i_wck];
        if (i_icn > 0 && i_ick > 0)
          output->dtauc_cldn[lc] -= (bsca_s[i_ick] + babs_s[i_ick]);

        bsca_cldn = bsca_tot - bsca_s[i_wck];
        if (i_icn > 0 && i_ick > 0)
          bsca_cldn -= bsca_s[i_ick];

        if (output->dtauc_cldn[lc] > 0)
          output->ssalb_cldn[lc] = bsca_cldn / output->dtauc_cldn[lc];
        else
          output->ssalb_cldn[lc] = 1.0;
      }

      /* ulrike: add dtau_mol, dtau_aer and tau_wc and tau_ic (the
                 latter two are obtained from tipa (dir), see
                 solve_rte.c (tipa_calcdtau); With tausol[lc] the
                 direct radiation is then calculated and used for the
                 calculation of the diffuse radiation */
      /* note: tausol is the SUM over all dtau's in layers above the
	 level (lc) considered!!!*/
      if (input.tipa == TIPA_DIR) {
        if (!input.quiet && lc == 0)
          fprintf (stderr, " ... (tipa dir) calculate total tau(sol) for every level\n");

        for (isp = CAOTH_FIR; isp < n_caoth; isp++)
          tau_babs_sum += babs_s[isp];

        output->tausol[lc] = tau_babs_sum;

        /* add dtau for caoth at the levels where caoth
	   (due to tipa dir) contribute */
        /* output->caoth.tipa.taudircld is the tau of caoth along
	   the beam, i.e. a sum over all layers up to 'caoth top' */
        for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
          ispo = isp - CAOTH_FIR;
          if ((nlyr - lc) <= output->caoth[ispo].tipa.nztilt && output->caoth[ispo].tipa.nztilt > 0)
            output->tausol[lc] += output->caoth[ispo].tipa.taudircld[iv][nlyr - lc - 1];
        }
      }

      if (1.0 - output->ssalb[lc] < FLT_MIN)
        output->ssalb[lc] = 1.0 - FLT_MIN;

    } /* endfor (lc=0; lc<nlyr; lc++) */

    /* end writing extinction file for ARLEM */
    if (input.write_ext_to_file)
      fclose (extfile);

    if ((input.rte.solver == SOLVER_FDISORT2 || input.rte.solver == SOLVER_DISORT) && input.rte.disort_icm == DISORT_ICM_PHASE) {
      /* we need to generalize the mu grid for the disort_icm phase,
         the explicit phase function is needed.  It must be defined on
         a mu grid which has to be identical for all atmospheric
         layers.  Before this point, the mu grid is not identical,
         now, we define a common mu grid. */

      tmp_ntheta_in = calloc ((unsigned int) nlyr, sizeof (int));
      tmp_theta_in  = calloc ((unsigned int) nlyr, sizeof (float*));
      tmp_mu_in     = calloc ((unsigned int) nlyr, sizeof (double*));

      /* point to mu grid of each layer which contains phase function */
      n_tot = 0;
      n_in  = 0;
      for (lc = 0; lc < nlyr; lc++) {
        if (output->ntheta[lc][0] > 0) {
          n_tot += output->ntheta[lc][0];
          tmp_ntheta_in[n_in] = output->ntheta[lc][0];
          tmp_theta_in[n_in]  = output->theta[lc][0];
          tmp_mu_in[n_in]     = output->mu[lc][0];
          n_in++;
        }
      }

      tmp_theta_new = calloc (n_tot, sizeof (float));
      tmp_mu_new    = calloc (n_tot, sizeof (double));

      /* find common mu grid */
      ntheta_new = sort_theta_and_mu (n_in, tmp_ntheta_in, tmp_theta_in, tmp_theta_new, tmp_mu_in, tmp_mu_new);
      if (ntheta_new == -1)
        return fct_err_out (-1, "sort_theta_and_mu", ERROR_POSITION);

      if (verbose)
        fprintf (stderr, "The phase function for SOLVER_FDISORT2/cdisort is using %d grid points ...\n", output->ntheta[0][0]);

      for (lc = 0; lc < nlyr; lc++) {
        /* interpolate phase */
        tmp_phase_new = calloc (ntheta_new, sizeof (float));

        status = interpolate_phase_weighted (1,
                                             &(output->ntheta[lc][0]),
                                             &(output->mu[lc][0]),
                                             &(output->phase[lc][0]),
                                             &one,
                                             ntheta_new,
                                             tmp_mu_new,
                                             tmp_phase_new,
                                             0);
        CHKERR (status);

        free (output->phase[lc][0]);
        output->phase[lc][0] = tmp_phase_new;
      }

      /* actually, only lc=0 is used from now on, but to avoid
	 programming errors, we redefine all ntheta's, theta's, and
	 mu's */
      for (lc = 0; lc < nlyr; lc++) {
        output->ntheta[lc][0] = ntheta_new;

        free (output->theta[lc][0]);
        free (output->mu[lc][0]);

        /* allocate theta dimension */
        output->theta[lc][0] = calloc (output->ntheta[0][0], sizeof (float));
        output->mu[lc][0]    = calloc (output->ntheta[0][0], sizeof (double));

        /* copy new theta grid into target */
        for (i = 0; i < output->ntheta[0][0]; i++) {
          output->theta[lc][0][i] = tmp_theta_new[i];
          output->mu[lc][0][i]    = tmp_mu_new[i];
        }
      }

      free (tmp_ntheta_in);
      free (tmp_theta_in);
      free (tmp_mu_in);
      free (tmp_theta_new);
      free (tmp_mu_new);

      /* phase function needs to be renormalized */
      for (lc = 0; lc < nlyr; lc++)
        normalize_phase (output->mu[0], output->phase[lc], NULL, output->ntheta[0], 1, 0);
    } /* end if disort && disort_icm phase */

    if (input.rte.solver == SOLVER_MONTECARLO) {

      for (isp = 0; isp < n_caoth; isp++) {
        output->mc.dt[isp][nlyr] = 0.0;
        output->mc.om[isp][nlyr] = 0.0;
      }

      /* loop includes nlyr because altitude and temperature */
      /* are defined per level                               */

      for (lc = 0; lc <= nlyr; lc++) {
        output->mc.z[nlyr - lc] = output->atm.zd[lc];
        if (!input.atmosphere3d)
          output->mc.temper[0][0][nlyr - lc] = output->atm.microphys.temper[0][0][lc];
      }

      /* loop stops at nlyr-1 because aerosol properties */
      /* are defined per layer*/
      output->mc.nphamataer = output->aer.optprop.nphamat;

      for (lc = 0; lc < nlyr; lc++) {

        output->mc.nmomaer[nlyr - lc - 1] = output->aer.optprop.nmom[iv][lc];
        output->mc.momaer[nlyr - lc - 1]  = calloc (output->aer.optprop.nphamat, sizeof (float*));

        for (ip = 0; ip < output->aer.optprop.nphamat; ip++) {

          output->mc.momaer[nlyr - lc - 1][ip] = calloc (output->mc.nmomaer[nlyr - lc - 1], sizeof (float));

          /* ??? should free this memory later ??? */

          for (k = 0; k < output->mc.nmomaer[nlyr - lc - 1]; k++) {
            output->mc.momaer[nlyr - lc - 1][ip][k] =
              output->aer.optprop.moment[iv][lc][ip][k] / output->aer.optprop.moment[iv][lc][0][0];
          }
        }
      }

      /* we need to copy phase functions too */
      if (output->aer.optprop.ntheta[iv] != NULL) {
        for (lc = 0; lc < nlyr; lc++) {
          output->mc.nthetaaer[nlyr - lc - 1] = output->aer.optprop.ntheta[iv][lc];
          output->mc.thetaaer[nlyr - lc - 1]  = output->aer.optprop.theta[iv][lc];  //TODO: Wird nicht benutzt?
          output->mc.muaer[nlyr - lc - 1]     = output->aer.optprop.mu[iv][lc];
          output->mc.phaseaer[nlyr - lc - 1]  = output->aer.optprop.phase[iv][lc];
        }
      }
    }

#if HAVE_SSSI
    if (input.rte.solver == SOLVER_SSSI) {
      /* Total caoth optical thickness and top caoth level */
      /* for the SOLVER_SSSI approximation                        */
      output->sssi.tautot = 0.0;
      for (isp = 0; isp < n_caoth; isp++)
        output->sssi.tautot += bsca_s_int[isp];
      for (lc = 0; lc < nlyr; lc++) {
        for (isp = CAOTH_FIR; isp < n_caoth; isp++) {
          ispo = isp - CAOTH_FIR;
          if (output->caoth[ispo].optprop.dtau[iv][lc] > 0.0)
            break;
        }
        if (isp != n_caoth)
          break;
      }

      output->sssi.lctop = lc;

      /* the uppermost layer determines cloud phase */
      output->sssi.type =
        (output->caoth[i_wc].optprop.dtau[iv][lc] > output->caoth[i_ic].optprop.dtau[iv][lc] ? ISCCP_WATER : ISCCP_ICE);
    }
#endif

    /* output vertical (total) sum */
    if (verbose) {
      if (output->cf.nlev == 0 || input.rte.solver == SOLVER_TWOMAXRND || input.rte.solver == SOLVER_TWOMAXRND3C ||
          input.rte.solver == SOLVER_DYNAMIC_TWOSTREAM || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) {
        fprintf (stderr,
                 " ----------------------------------------------------------------------------------------------------------------"
                 "-------------------------------------------------\n");
        fprintf (stderr,
                 "%5s |  %7.3f | %12.6e | %9.6f %9.6f %5.3f | %11.6f %11.6f %5.3f | %11.6f %11.6f %5.3f %5.3f %5.3f %6.3f %5.3f | "
                 "%11.6f\n",
                 "sum",
                 0.0 / 0.0,
                 bsca_s_int[CAOTH_MOL],
                 bsca_s_int[CAOTH_AER],
                 babs_s_int[CAOTH_AER],
                 0.0 / 0.0,
                 bsca_s_unsc_int[i_wc],
                 babs_s_unsc_int[i_wc],
                 0.0 / 0.0,
                 bsca_s_unsc_int[i_ic],
                 babs_s_unsc_int[i_ic],
                 0.0 / 0.0,
                 0.0 / 0.0,
                 0.0 / 0.0,
                 0.0 / 0.0,
                 0.0 / 0.0,
                 babs_s_int[CAOTH_MOL]);
        fprintf (stderr,
                 " ----------------------------------------------------------------------------------------------------------------"
                 "-------------------------------------------------\n");
      } else {
        fprintf (stderr,
                 " ----------------------------------------------------------------------------------------------------------------"
                 "----------------------\n");
        fprintf (stderr,
                 "%5s |  %7.3f | %12.6e | %9.6f %9.6f %5.3f | %11.6f %11.6f %5.3f %5.3f %5.3f %6.3f %5.3f | %11.6f\n",
                 "sum",
                 0.0 / 0.0,
                 bsca_s_int[CAOTH_MOL],
                 bsca_s_int[CAOTH_AER],
                 babs_s_int[CAOTH_AER],
                 0.0 / 0.0,
                 bsca_s_unsc_int[i_wc],
                 babs_s_unsc_int[i_wc],
                 0.0 / 0.0,
                 0.0 / 0.0,
                 0.0 / 0.0,
                 0.0 / 0.0,
                 0.0 / 0.0,
                 babs_s_int[CAOTH_MOL]);
        fprintf (stderr,
                 " ----------------------------------------------------------------------------------------------------------------"
                 "---------------------\n");
      }
    }
  } /*  if (!skip_optical_properties) { */
  else {
    if (verbose)
      fprintf (stderr, " *** skip calculation of optical properties! iv = %4d, iq = %4d \n", iv, iq);
  }

  /* print only integrated optical properties to stderr */
  /* fprintf (stderr, "%.0f  %7.3f  %12.6e  %9.6f %9.6f %5.3f  %11.6f %11.6f %5.3f  %11.6f %11.6f %5.3f %5.3f %5.3f %6.3f %5.3f  %11.6f\n",  */
  /*          0.0, 0.0/0.0, */
  /*          bsca_s_int[CAOTH_MOL], */
  /*          bsca_s_int[CAOTH_AER], babs_s_int[CAOTH_AER], 0.0/0.0, */
  /*          bsca_s_unsc_int[i_wc], babs_s_unsc_int[i_wc], 0.0/0.0, */
  /*          bsca_s_unsc_int[i_ic], babs_s_unsc_int[i_ic], 0.0/0.0, */
  /*          0.0/0.0, 0.0/0.0, 0.0/0.0, 0.0/0.0, */
  /*           babs_s_int[CAOTH_MOL]); */

  for (isp = 0; isp < n_caoth; isp++)
    free (mom_s[isp]);

  free (babs_s);
  free (bsca_s);
  free (ssa_s);
  free (babs_s_int);
  free (bsca_s_int);

  free (mom_s);
  free (mom_s_g1);
  free (mom_s_g2);
  free (dscale_s);
  free (f_s);
  free (ff_s);
  free (g1_s);
  free (g2_s);
  free (dtau_s);
  free (mom0_s);

  free (babs_s_unsc);
  free (bsca_s_unsc);
  free (ssa_s_unsc);
  free (babs_s_unsc_int);
  free (bsca_s_unsc_int);
  free (gg_s_unsc);

  free (ip_act);
#if HAVE_NETCDF4

  // write optical properties to file for test suite //
  if (input.write_optical_properties) {
    int  ncid;
    char optpropfilename[FILENAME_MAX + 1] = "";
    int  nwvl                              = 1;
    nwvl                                   = output->wl.nlambda_r;

    if (nwvl == 1)
      snprintf (optpropfilename, sizeof (optpropfilename), "optical_properties.nc");
    else
      snprintf (optpropfilename, sizeof (optpropfilename), "optical_properties.iv%d.nc", iv);

    /**********  Create netcdf file **********/
    status = create_new_netcdf (optpropfilename, &ncid);
    CHKERR (status);

    status = write_netcdf_3D_float (ncid,
                                    output->pmom,
                                    nlyr,
                                    nphamat,
                                    output->atm.nmom + 1,
                                    "output_pmom",
                                    "nlyr",
                                    "nphamat",
                                    "nmom+1",
                                    NULL);
    CHKERR (status);
    if (((input.rte.solver == SOLVER_FDISORT2 || input.rte.solver == SOLVER_DISORT) && input.rte.disort_icm == DISORT_ICM_PHASE) ||
        input.rte.solver == SOLVER_SSLIDAR) {
      status = write_netcdf_2D_int (ncid, output->ntheta, nlyr, nphamat, "output_ntheta", "nlyr", "nphamat", NULL);
      CHKERR (status);

      status = write_netcdf_3Dirr_float (ncid,
                                         output->phase,
                                         nlyr,
                                         nphamat,
                                         output->ntheta,
                                         "output_phase",
                                         "nlyr",
                                         "nphamat",
                                         "output_ntheta");
      CHKERR (status);
      status = write_netcdf_3Dirr_float (ncid,
                                         output->theta,
                                         nlyr,
                                         nphamat,
                                         output->ntheta,
                                         "output_theta",
                                         "nlyr",
                                         "nphamat",
                                         "output_ntheta");
      CHKERR (status);
      status = write_netcdf_3Dirr_double (ncid,
                                          output->mu,
                                          nlyr,
                                          nphamat,
                                          output->ntheta,
                                          "output_mu",
                                          "nlyr",
                                          "nphamat",
                                          "output_ntheta");
      CHKERR (status);
      status = write_netcdf_1D_float (ncid, output->dtauc, nlyr, "output_dtauc", "nlyr", NULL);
      CHKERR (status);
      status = write_netcdf_1D_float (ncid, output->dtauc_md, nlyr, "output_dtauc_md", "nlyr", NULL);
      CHKERR (status);
      status = write_netcdf_1D_float (ncid, output->ssalb, nlyr, "output_ssalb", "nlyr", NULL);
      CHKERR (status);
    }
    if (input.tipa == TIPA_DIR) {
      status = write_netcdf_1D_float (ncid, output->tausol, nlyr, "output_tausol", "nlyr", NULL);
      CHKERR (status);
    }
#if HAVE_SSSI
    if (input.rte.solver == SOLVER_SSSI) {
      status = write_netcdf_0D_float (ncid, output->sssi.tautot, "output_sssi.tautot");
      CHKERR (status);
      status = write_netcdf_0D_int (ncid, output->sssi.lctop, "output_sssi.lctop");
      CHKERR (status);
      status = write_netcdf_0D_int (ncid, output->sssi.type, "output_sssi.type");
      CHKERR (status);
    }
#endif
    if (input.rte.solver == SOLVER_MONTECARLO) {
      if (output->mc.alis.dt != NULL) {
        status =
          write_netcdf_2D_double (ncid, output->mc.alis.dt[iv], n_caoth, nlev, "output_mc.alis.dt[iv]", "ncaoth", "nlev", NULL);
        CHKERR (status);
      }
      if (output->mc.alis.om != NULL) {
        status =
          write_netcdf_2D_double (ncid, output->mc.alis.om[iv], n_caoth, nlev, "output_mc.alis.om[iv]", "ncaoth", "nlev", NULL);
        CHKERR (status);
      }
      status = write_netcdf_2D_float (ncid, output->mc.dt, n_caoth, nlev, "output_dt", "ncaoth", "nlev", NULL);
      CHKERR (status);
      status = write_netcdf_2D_float (ncid, output->mc.om, n_caoth, nlev, "output_om", "ncaoth", "nlev", NULL);
      CHKERR (status);
      status = write_netcdf_2D_float (ncid, output->mc.g1, n_caoth, nlev, "output_g1", "ncaoth", "nlev", NULL);
      CHKERR (status);
      status = write_netcdf_2D_float (ncid, output->mc.g2, n_caoth, nlev, "output_g2", "ncaoth", "nlev", NULL);
      CHKERR (status);
      status = write_netcdf_2D_float (ncid, output->mc.ff, n_caoth, nlev, "output_ff", "ncaoth", "nlev", NULL);
      CHKERR (status);
      status = write_netcdf_2D_float (ncid, output->mc.ds, n_caoth, nlev, "output_ds", "ncaoth", "nlev", NULL);
      CHKERR (status);
      status = write_netcdf_2D_float (ncid, output->mc.re, n_caoth, nlev, "output_re", "ncaoth", "nlev", NULL);
      CHKERR (status);
      status = write_netcdf_1D_float (ncid, output->mc.refind, nlev, "output_mc.refind", "nlev", NULL);
      CHKERR (status);
      status = write_netcdf_1D_float (ncid, output->mc.z, nlev, "output_mc.z", "nlev", NULL);
      CHKERR (status);
      status = write_netcdf_1D_float (ncid, output->mc.temper[0][0], nlev, "output_mc.temper", "nlev", NULL);
      CHKERR (status);
      status = write_netcdf_0D_int (ncid, output->mc.nphamataer, "output_mc.nphamataer");
      CHKERR (status);
      {
        size_t max_num = 0;
        for (size_t k = 0; k < nlev; k++) {
          max_num = MAX (max_num, output->mc.nmomaer[k]);
        }
        if (max_num) {
          status = write_netcdf_3Dirr_row_float (ncid,
                                                 output->mc.momaer,
                                                 nlev,
                                                 output->aer.optprop.nphamat,
                                                 output->mc.nmomaer,
                                                 "output_mc.momaer",
                                                 "nlev",
                                                 "output_aer.optprop.nphamat",
                                                 "output_mc.nmomaer");
          CHKERR (status);
        }
      }
      if (output->mc.nthetaaer != NULL) {
        status = write_netcdf_2D_int (ncid,
                                      output->mc.nthetaaer,
                                      nlyr,
                                      output->aer.optprop.nphamat,
                                      "output_mc.nthetaaer",
                                      "nlyr",
                                      "output_aer.optprop.nphamat",
                                      NULL);
        CHKERR (status);
        status = write_netcdf_3Dirr_float (ncid,
                                           output->mc.thetaaer,
                                           nlyr,
                                           output->aer.optprop.nphamat,
                                           output->mc.nthetaaer,
                                           "output_mc.thetaaer",
                                           "nlyr",
                                           "output_aer.optprop.nphamat",
                                           "output_mc.nthetaaer");
        CHKERR (status);
        status = write_netcdf_3Dirr_double (ncid,
                                            output->mc.muaer,
                                            nlyr,
                                            output->aer.optprop.nphamat,
                                            output->mc.nthetaaer,
                                            "output_mc.muaer",
                                            "nlyr",
                                            "output_aer.optprop.nphamat",
                                            "output_mc.nthetaaer");
        CHKERR (status);
        status = write_netcdf_3Dirr_float (ncid,
                                           output->mc.phaseaer,
                                           nlyr,
                                           output->aer.optprop.nphamat,
                                           output->mc.nthetaaer,
                                           "output_mc.phaseaer",
                                           "nlyr",
                                           "output_aer.optprop.nphamat",
                                           "output_mc.nthetaaer");
        CHKERR (status);
      }
    }

    // Dump 3D optical properties in caoths
    for (isp = 0; isp < input.n_caoth; isp++) {
      status = dump_caoth3d_var (ncid, isp, input.caoth[isp].name, output->caoth3d[isp], "_lwc", output->caoth3d[isp].lwc);
      CHKERR (status);
      status = dump_caoth3d_var (ncid, isp, input.caoth[isp].name, output->caoth3d[isp], "_reff", output->caoth3d[isp].reff);
      CHKERR (status);
      status = dump_caoth3d_var (ncid, isp, input.caoth[isp].name, output->caoth3d[isp], "_ext", output->caoth3d[isp].ext);
      CHKERR (status);
      status = dump_caoth3d_var (ncid, isp, input.caoth[isp].name, output->caoth3d[isp], "_ssa", output->caoth3d[isp].ssa);
      CHKERR (status);
      status = dump_caoth3d_var (ncid, isp, input.caoth[isp].name, output->caoth3d[isp], "_g1", output->caoth3d[isp].g1);
      CHKERR (status);
      status = dump_caoth3d_var (ncid, isp, input.caoth[isp].name, output->caoth3d[isp], "_g2", output->caoth3d[isp].g2);
      CHKERR (status);
      status = dump_caoth3d_var (ncid, isp, input.caoth[isp].name, output->caoth3d[isp], "_ff", output->caoth3d[isp].ff);
      CHKERR (status);
      status = dump_caoth3d_var (ncid, isp, input.caoth[isp].name, output->caoth3d[isp], "_f", output->caoth3d[isp].f);
      CHKERR (status);
    }  // end Dump 3D optical properties in caoths

    /********** Close the file. This frees up any internal netCDF resources
    * associated with the file, and flushes any buffers. **********/
    status = close_netcdf (ncid);
    CHKERR (status);
  }
  // finish writing optical properties to file for test suite //
  return 0;

#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
  fprintf (stderr, " ***********************************************************************\n");
  return -1;
#endif
}

/***************************************/
/* Post - processing of the output     */
/* sum or integration over wavelength, */
/* calculation of heating rates or     */
/* conversion to RGB                   */
/***************************************/

int processing1D (input_struct input, output_struct* output) {
  int status = 0, iv = 0;

  switch (input.processing) {

  case PROCESS_NONE:
    if (input.calibration == OUTCAL_BRIGHTNESS) {
      if (!input.quiet)
        fprintf (stderr, " ... converting radiances to brightness temperatures\n");

      /* convert irradiances / radiances to brightness temperatures */
      for (iv = 0; iv < output->wl.nlambda_h; iv++) {
        status = output2bt (input, output, iv, 0);
        CHKERR (status);
      }
    }
    break;

  case PROCESS_SUM:
    status = sum1D (input, output); /* multiply with extraterrestrial and sum up (function in ancillary.c) */
    CHKERR (status);
    break;

  case PROCESS_INT:
    status = integrate1D (input, output); /* multiply with extraterrestrial and integrate (function in ancillary.c) */
    CHKERR (status);
    break;

  case PROCESS_RGB:
  case PROCESS_RGBRAW:
    status = spec2rgb (input, output); /* convert spectrum to red,green,blue space (function in ancillary.c) */
    CHKERR (status);
    break;

  case PROCESS_RAMAN:
    status = raman_spec2spec (input, output); /* convert Raman spectrum with extra wavelengths to user spectrum */
    CHKERR (status);
    break;

  default:
    CHKERROUT (1000 + input.processing, "unknown processing scheme")
  }

  return 0;
}

/***********************************************************************************/
/* Sum the 1D results over wavelength and write data to the first array element    */
/***********************************************************************************/

int sum1D (input_struct input, output_struct* output) {
  int status = 0;

  /* calculate incident flux (output->incident) and sum extraterrestrial irradiance (wl.fbeam[0]) */
  for (size_t iv = 0; iv < output->wl.nlambda_h; iv++) {
    output->incident += output->wl.filter[iv] * output->wl.fbeam[iv] * cos (output->sza_h[iv] * M_PI / 180.0);
    if (iv != 0)
      output->wl.fbeam[0] += output->wl.filter[iv] * output->wl.fbeam[iv];
    else /* iv == 0 */
      output->wl.fbeam[0] = output->wl.filter[0] * output->wl.fbeam[0];
  }

  /* callocate integrated values */
  status = calloc_int_values (input, output);
  CHKERR (status);

  /* fluxes and radiances */
  for (size_t lev = 0; lev < output->atm.nzout; lev++) {
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++) {

      output->rfldir_int[lev] += (double) output->rfldir[lev][iv];
      output->rfldn_int[lev] += (double) output->rfldn[lev][iv];
      output->flup_int[lev] += (double) output->flup[lev][iv];
      output->uavg_int[lev] += (double) output->uavg[lev][iv];
      output->uavgso_int[lev] += (double) output->uavgso[lev][iv];
      output->uavgdn_int[lev] += (double) output->uavgdn[lev][iv];
      output->uavgup_int[lev] += (double) output->uavgup[lev][iv];
      if (input.heating != HEAT_NONE) {
        output->heat_int[lev] += (double) output->heat[lev][iv];
        output->emis_int[lev] += (double) output->emis[lev][iv];
        output->w_zout_int[lev] += (double) output->w_zout[lev][iv];
      }

      for (size_t iu = 0; iu < input.rte.numu; iu++) {
        output->u0u_int[lev][iu] += output->u0u[lev][iu][iv];

        for (size_t j = 0; j < input.rte.nphi; j++)
          output->uu_int[lev][j][iu] += output->uu[lev][j][iu][iv];
      }
    }
  }

  /* albedo, transmittance of total atmosphere */
  for (size_t iu = 0; iu < input.rte.numu; iu++) {
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++) {
      output->albmed_int[iu] += (double) output->albmed[iu][iv];
      output->trnmed_int[iu] += (double) output->trnmed[iu][iv];
    }
  }

  /* store integrated values in the first entry of the wavelength index, if any result fields are defined */
  if (output->wl.nlambda_h > 0)
    status += double2float_integrated_values (input, output);

  if (status != 0) {
    fprintf (stderr, "Error, conversion from double to float not possible");
    fprintf (stderr, "sum1D(), ancillary.c, status %d\n", status);
    CHKERR (status);
  }

  status += scaling_integrated_values (input, output);
  if (status != 0) {
    fprintf (stderr, "Error, scaling integrated values.\n");
    fprintf (stderr, "scaling_integrated_values(), in ancillary.c, status %d\n", status);
    CHKERR (status);
  }

  return status;
}

/**************************************************************************************/
/* Integrate the 1D results over wavelength and write data to the first array element */
/**************************************************************************************/

int integrate1D (input_struct input, output_struct* output) {
  int     status = 0, used_unit = 0;
  double *xint = NULL, *yint = NULL;
  double* cos_SZA = NULL;

  /* callocate integrated values */
  status = calloc_int_values (input, output);
  CHKERR (status);

  /* calculate incident flux */
  xint = calloc (output->wl.nlambda_h, sizeof (double));
  CHKPOINTER (xint);
  yint = calloc (output->wl.nlambda_h, sizeof (double));
  CHKPOINTER (yint);

  cos_SZA = calloc (output->wl.nlambda_h, sizeof (double));
  CHKPOINTEROUT (cos_SZA, "failed allocating memory for cos_SZA");

  switch (input.source) {
  case SRC_SOLAR:
  case SRC_LIDAR: /* BCA */
  case SRC_BLITZ: /* BCA */
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      cos_SZA[iv] = cos (output->sza_h[iv] * M_PI / 180.0);
    break;
  case SRC_THERMAL:
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      cos_SZA[iv] = 1.0;
    break;
  default:
    CHKERROUT (1000 + input.source, "unknown source");
  }

  used_unit = input.output_unit;
  if (used_unit == UNIT_NOT_DEFINED) /* if no output unit is specified, assume the same units as the input spectrum */
    used_unit = output->spectrum_unit;

  switch (used_unit) {
  case UNIT_PER_NM:
    if (input.verbose)
      fprintf (stderr, " *** integration in wavelength space \n");
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++) {
      xint[iv] = (double) output->wl.lambda_h[iv];
      yint[iv] = output->wl.filter[iv] * output->wl.fbeam[iv] * cos_SZA[iv];
    }
    break;
  case UNIT_PER_CM_1:
    /* integration over wavenumber k, minus in order to get ascending wavenumbers */
    if (input.verbose)
      fprintf (stderr, " *** integration in wavenumber space \n");
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++) {
      xint[iv] = -(double) 1.0e+7 / output->wl.lambda_h[iv]; /* k = 10**7 / lambda */ /* 10**7 == nm -> cm */
      yint[iv] = output->wl.filter[iv] * output->wl.fbeam[iv] * cos_SZA[iv];
    }
    break;
  case UNIT_PER_BAND:
    fprintf (stderr, "Error, combination 'output_process integrate' and 'output_process per_band',\n");
    fprintf (stderr, "    or combination 'output_process integrate' and an input spectrum defined in \n");
    fprintf (stderr, "       (wavelength or wavenumber) bands does not make sense \n");
    fprintf (stderr, "       please use 'output_process sum', when dealing with band parametrisations. \n\n");
    return -1;
    break;
  case UNIT_NOT_DEFINED:
    fprintf (stderr, "Error, in order to use 'output_process integrate' it is nessesary to specify the \n");
    fprintf (stderr, "       unit of the extraterrestrial spectrum with 'source solar filename unit' or\n");
    fprintf (stderr, "       unit of the output with 'output_process per_nm' or 'output_process per_cm-1' \n\n");
    return -1;
    break;
  default:
    fprintf (stderr,
             "Error: Program bug, unsupported unit of extraterrestial flux %d or output unit %d\n",
             output->spectrum_unit,
             input.output_unit);
    return -1;
  }

  output->incident = integrate (xint, yint, output->wl.nlambda_h);

  /* sum extraterrestrial irradiance */
  for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
    yint[iv] = output->wl.filter[iv] * output->wl.fbeam[iv];
  /* unit factor should alway be 1 here == output always in W/(m2 nm) */
  /* therefore it is not included here in the code */

  output->wl.fbeam[0] = (float) integrate (xint, yint, output->wl.nlambda_h);

  /* fluxes and radiances */
  for (size_t iu = 0; iu < input.rte.numu; iu++) {
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      yint[iv] = output->albmed[iu][iv];
    output->albmed_int[iu] = integrate (xint, yint, output->wl.nlambda_h);
  }
  for (size_t iu = 0; iu < input.rte.numu; iu++) {
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      yint[iv] = output->trnmed[iu][iv];
    output->trnmed_int[iu] = integrate (xint, yint, output->wl.nlambda_h);
  }

  for (size_t lev = 0; lev < output->atm.nzout; lev++) {

    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      yint[iv] = output->rfldir[lev][iv];
    output->rfldir_int[lev] = integrate (xint, yint, output->wl.nlambda_h);
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      yint[iv] = output->rfldn[lev][iv];
    output->rfldn_int[lev] = integrate (xint, yint, output->wl.nlambda_h);
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      yint[iv] = output->flup[lev][iv];
    output->flup_int[lev] = integrate (xint, yint, output->wl.nlambda_h);
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      yint[iv] = output->uavg[lev][iv];
    output->uavg_int[lev] = integrate (xint, yint, output->wl.nlambda_h);
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      yint[iv] = output->uavgso[lev][iv];
    output->uavgso_int[lev] = integrate (xint, yint, output->wl.nlambda_h);
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      yint[iv] = output->uavgdn[lev][iv];
    output->uavgdn_int[lev] = integrate (xint, yint, output->wl.nlambda_h);
    for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
      yint[iv] = output->uavgup[lev][iv];
    output->uavgup_int[lev] = integrate (xint, yint, output->wl.nlambda_h);
    if (input.heating != HEAT_NONE) {
      for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
        yint[iv] = output->heat[lev][iv];
      output->heat_int[lev] = integrate (xint, yint, output->wl.nlambda_h);
      for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
        yint[iv] = output->emis[lev][iv];
      output->emis_int[lev] = integrate (xint, yint, output->wl.nlambda_h);
      for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
        yint[iv] = output->w_zout[lev][iv];
      output->w_zout_int[lev] = integrate (xint, yint, output->wl.nlambda_h);
    }

    for (size_t iu = 0; iu < input.rte.numu; iu++) {
      for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
        yint[iv] = (double) output->u0u[lev][iu][iv];
      output->u0u_int[lev][iu] = integrate (xint, yint, output->wl.nlambda_h);

      for (size_t j = 0; j < input.rte.nphi; j++) {
        for (size_t iv = 0; iv < output->wl.nlambda_h; iv++)
          yint[iv] = (double) output->uu[lev][j][iu][iv];
        output->uu_int[lev][j][iu] = integrate (xint, yint, output->wl.nlambda_h);
      }
    }
  }

  free (xint);
  free (yint);

  /* store integrated values in the first entry of the wavelength index */
  if (output->wl.nlambda_h > 0)
    status += double2float_integrated_values (input, output);

  if (status != 0) {
    fprintf (stderr, "Error, conversion from double to float not possible");
    fprintf (stderr, "integrate1D(), ancillary.c, status %d\n", status);
    CHKERR (status);
  }

  status += scaling_integrated_values (input, output);
  if (status != 0) {
    fprintf (stderr, "Error, scaling integrated values.\n");
    fprintf (stderr, "integrate1D(), in ancillary.c, status %d\n", status);
    CHKERR (status);
  }

  return status;
}

static int write_triangular_surface_results (const input_struct input, const output_struct* output) {
  if (!output->triangle_results_o)
    return 0;
  if (!output->triangle_results_o[0]->N_triangles)
    return 0;
  const t_triangle_radiation_field* out = output->triangle_results_o[0];

  const int ldebug = 0;
  if (ldebug) {
    for (size_t i = 0; i < out->N_triangles; ++i) {
      fprintf (stderr, "Surface fluxes on triangle %lu => %f %f %f\n", i, out->edir[i], out->edn[i], out->eup[i]);
    }
  }

  {  // write netcdf
    int  ierr;
    int  ncid;
    char filename[FILENAME_MAX] = "";
    strcpy (filename, input.rte.mc.filename[FN_MC_BASENAME]);
    strcat (filename, ".triangle.nc");
    // Create the file. overwrite this file, if it already exists
    ierr = create_new_netcdf (filename, &ncid);
    CHKERR (ierr);

    ierr = write_netcdf_1D_size_t (ncid, out->ndir, out->N_triangles, "count_edir", "Ntriangles", NULL);
    CHKERR (ierr);
    ierr = write_netcdf_1D_size_t (ncid, out->ndn, out->N_triangles, "count_edn", "Ntriangles", NULL);
    CHKERR (ierr);
    ierr = write_netcdf_1D_size_t (ncid, out->nup, out->N_triangles, "count_eup", "Ntriangles", NULL);
    CHKERR (ierr);

    ierr = write_netcdf_1D_double (ncid, out->edir, out->N_triangles, "edir", "Ntriangles", NULL);
    CHKERR (ierr);
    ierr = write_netcdf_1D_double (ncid, out->edn, out->N_triangles, "edn", "Ntriangles", NULL);
    CHKERR (ierr);
    ierr = write_netcdf_1D_double (ncid, out->eup, out->N_triangles, "eup", "Ntriangles", NULL);
    CHKERR (ierr);

    ierr = close_netcdf (ncid);
    CHKERR (ierr);
  }
  return 0;
}

/***************************************/
/* Post - processing of the output     */
/* sum or integration over wavelength, */
/* calculation of heating rates or     */
/* conversion to RGB                   */
/***************************************/

int processing3D (input_struct input, output_struct* output) {
  int status;

  if (input.rte.mc.locest)
    return 0;

  output->wl.nlambda_h_print3D = output->wl.nlambda_h;

  switch (input.processing) {

  case PROCESS_NONE:
    if (input.calibration == OUTCAL_BRIGHTNESS) {
      if (!input.quiet)
        fprintf (stderr, " ... converting radiances to brightness temperatures\n");

      /* convert irradiances / radiances to brightness temperatures */
      for (int iv = 0; iv < output->wl.nlambda_h; iv++) {
        status = output2bt (input, output, iv, 1);
        CHKERR (status);
      }
    }
    break;

  case PROCESS_SUM:
    status = sum3D (input, output);
    CHKERR (status);
    output->wl.nlambda_h_print3D = 1;
    break;

  case PROCESS_INT:
    status = integrate3D (input, output);
    CHKERR (status);
    output->wl.nlambda_h_print3D = 1;
    break;

  case PROCESS_RGB:
  case PROCESS_RGBRAW:
    status = spec2rgb3D (input, output); /* convert spectrum to red,green,blue space (function in ancillary.c) */
    CHKERR (status);
    output->wl.nlambda_h_print3D = 3;
    break;

  default:
    fprintf (stderr, "Error, unknown processing scheme %d\n", input.processing);
    CHKERR (-1);
  }

  /* convert unit of absorbed irradiance */

  switch (input.rte.mc.abs_unit) {
  case MCABS_UNIT_W_PER_M2_AND_DZ:
    /* default -> no change */
    break;

  case MCABS_UNIT_W_PER_M3:
  case MCABS_UNIT_K_PER_DAY:

    if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
        input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
      size_t k3d = 0;  // count 3d layers, because output->abs3d_r has only entries for 3d layers
      for (int ks = 0; ks < output->atm.Nzcld; ks++) {
        if (output->atm.threed[ks]) { /* only for 3D layers, BM07122005 */

          float scale_factor_abs3d = 1.0;

          /* find model level corresponding to the user level - there must be a better way to do that! */
          if (input.rte.mc.abs_unit == MCABS_UNIT_K_PER_DAY) {

            const float rho_air = output->atm.microphys.dens_avg[MOL_AIR][0][0][output->atm.Nzcld - 1 - ks] * 1.e+6 * 1.e-3 *
                                  input.atm.mol_mass[MOL_AIR] / AVOGADRO;

            /* 1.e+6: convert from cm-3 to m-3; 1.e-3: convert g -> kg */
            float c_p;
            status = specific_heat_capacity_moist_air (output->atm.microphys.temper_avg[0][0][output->atm.Nzcld - 1 - ks],
                                                       output->atm.microphys.dens_avg[MOL_AIR][0][0][output->atm.Nzcld - 1 - ks],
                                                       output->atm.microphys.dens_avg[MOL_H2O][0][0][output->atm.Nzcld - 1 - ks],
                                                       &(c_p),
                                                       input.quiet);
            CHKERR (status);
            scale_factor_abs3d = scale_factor_abs3d * s2day / (c_p * rho_air);
          }

          const float dz =
            (output->atm.zd[output->atm.Nzcld - ks - 1] - output->atm.zd[output->atm.Nzcld - ks]) * 1000.0; /* 1000 == km -> m */
          scale_factor_abs3d = scale_factor_abs3d / dz;

          for (int is = 0; is < output->atm.Nxcld; is++) {
            for (int js = 0; js < output->atm.Nycld; js++) {
              for (int iv = 0; iv < output->wl.nlambda_h_print3D; iv++) {
                output->abs3d[k3d][is][js][iv] *= scale_factor_abs3d;
                if (input.rte.mc.std) {
                  output->abs3d_var[k3d][is][js][iv] *= (scale_factor_abs3d * scale_factor_abs3d);
                }
              }
            }
          }
          k3d++;
        }
      }
    }

    if (input.rte.mc.backward.absorption) {
      for (int ks = 0; ks < output->atm.nzout; ks++) {

        float scale_factor_abs3d = 1.0;

        float zbot, ztop;
        int   kbot, ktop;
        status = mc_backward_absorption_find_vert_index (&output->atm, output->atm.zout_sur[ks], &kbot, &ktop, &zbot, &ztop);
        CHKERR (status);

        if (input.rte.mc.abs_unit == MCABS_UNIT_K_PER_DAY) {
          const float rho_air =
            output->atm.microphys.dens_avg[MOL_AIR][0][0][ktop] * 1.e+6 * 1.e-3 * input.atm.mol_mass[MOL_AIR] / AVOGADRO;

          /* 1.e+6: convert from cm-3 to m-3; 1.e-3: convert g -> kg */
          float c_p;
          status = specific_heat_capacity_moist_air (output->atm.microphys.temper_avg[0][0][ktop],
                                                     output->atm.microphys.dens_avg[MOL_AIR][0][0][ktop],
                                                     output->atm.microphys.dens_avg[MOL_H2O][0][0][ktop],
                                                     &(c_p),
                                                     input.quiet);
          CHKERR (status);
          scale_factor_abs3d = scale_factor_abs3d * s2day / (c_p * rho_air);
        }

        const float dz     = (ztop - zbot) * 1000;
        scale_factor_abs3d = scale_factor_abs3d / dz;

        if (!input.quiet && input.ipa3d != 1) /*ulrike added && input.ipa3d!=1*/
          fprintf (stderr,
                   "converting to heating rate, level %d %.3f - %.3f km, dens=%g, temper=%.3f, scale_factor %f\n",
                   kbot,
                   zbot,
                   ztop,
                   output->atm.microphys.dens_avg[MOL_AIR][0][0][ktop],
                   output->atm.microphys.temper_avg[0][0][ktop],
                   scale_factor_abs3d);

        /* ulrike 4.5.2010: absback3d is used to save the heating rates in case of ipa3d and
	   absback3d is written into mc.abs.spc when mc_backward_output heat K_per_day
	   is specified in the input-file; however, no further scaling of absback3d is necessary
	   (scaling of absback3d was done in scale_output in ancillary.c) */
        if (input.ipa3d != 1)
          for (int is = output->islower; is <= output->isupper; is += output->isstep)
            for (int js = output->jslower; js <= output->jsupper; js += output->jsstep)
              for (int iv = 0; iv < output->wl.nlambda_h_print3D; iv++) {
                output->absback3d[ks][is][js][iv] *= scale_factor_abs3d;
                /* 27.02.2013 **CK **BM: for thermal backward heating rates std */
                if (input.rte.mc.std)
                  output->absback3d_var[ks][is][js][iv] *= (scale_factor_abs3d * scale_factor_abs3d);
              }
      } /* endfor (ks=0; ks<output->atm.nzout; ks++) */
    }   /* endif (input.rte.mc.backward.absorption) */

    break;
  default:
    fprintf (stderr, "Error, unknown abs_unit %d\n", input.rte.mc.abs_unit);
    CHKERR (-1);
  }

  /* write 3D data to files */
  if (input.write_output_as_netcdf) {
    status = write_spectral3D_nc (input, output);
    CHKERR (status);
  } else {
    status = write_spectral3D (input, output);
    CHKERR (status);
  }

  status = write_triangular_surface_results (input, output);
  CHKERR (status);

  return 0;
}

/***********************************************************************************/
/* Sum the 3D results over wavelength and print data to mc.sum                     */
/***********************************************************************************/

static int write_spectral3D (input_struct input, output_struct* output) {
  const time_t starttime = time (NULL);

  int doflx = 1; /* switch that turns off writing of flx files. Should be determined automatically */

  char flxfilename[FILENAME_MAX]   = "";
  char radfilename[FILENAME_MAX]   = "";
  char absfilename[FILENAME_MAX]   = "";
  char jacfilename[FILENAME_MAX]   = "";
  char flxisfilename[FILENAME_MAX] = "";

  char flxvarfilename[FILENAME_MAX] = "";
  char radvarfilename[FILENAME_MAX] = "";
  char absvarfilename[FILENAME_MAX] = ""; /* 27.02.2013 **CK **BM: add new variable for thermal backward heating rates std */

  FILE *fflx = NULL, *fflxvar = NULL, *fabs = NULL, *fabsvar = NULL, *frad = NULL, *fradvar = NULL, *fjac = NULL,
       *fflxis = NULL; /* 27.02.2013 **CK **BM: add *fabsvar=NULL for thermal backward heating rates std */

  /* generate output file names */
  strcpy (flxfilename, input.rte.mc.filename[FN_MC_BASENAME]);
  strcpy (absfilename, input.rte.mc.filename[FN_MC_BASENAME]);
  strcpy (radfilename, input.rte.mc.filename[FN_MC_BASENAME]);

  strcpy (flxvarfilename, input.rte.mc.filename[FN_MC_BASENAME]);
  strcpy (radvarfilename, input.rte.mc.filename[FN_MC_BASENAME]);
  strcpy (absvarfilename,
          input.rte.mc.filename[FN_MC_BASENAME]); /* 27.02.2013 **CK **BM: added for thermal backward heating rates std */

  strcat (flxfilename, ".flx.spc");
  strcat (radfilename, ".rad.spc");

  strcat (flxvarfilename, ".flx.std.spc");
  strcat (radvarfilename, ".rad.std.spc");

  if (input.rte.mc.jacobian[DIM_1D]) {
    strcpy (jacfilename, input.rte.mc.filename[FN_MC_BASENAME]);
    strcat (jacfilename, ".jac.spc");
  }

  if (input.rte.mc.concentration_is || input.rte.mc.spectral_is) {
    strcpy (flxisfilename, input.rte.mc.filename[FN_MC_BASENAME]);
    strcat (flxisfilename, ".flx.is.spc");
  }

  /* extension for absorption/heating/actinic etc. file */
  switch (input.rte.mc.absorption) {
  case MCFORWARD_ABS_ACTINIC:
    strcat (absfilename, ".act.spc");
    strcat (absvarfilename, ".act.std.spc"); /* 27.02.2013 **CK **BM: added for thermal backward heating rates std */
    break;

  case MCFORWARD_ABS_HEATING:
  case MCFORWARD_ABS_NONE:
    strcat (absfilename, ".abs.spc");
    strcat (absvarfilename, ".abs.std.spc"); /* 27.02.2013 **CK **BM: added for thermal backward heating rates std */
    break;

  default:
    fprintf (stderr, "Error, unknown absorption type %d\n", input.rte.mc.absorption);
    break;
  }

  if (!input.quiet) {
    int Nx, Ny;
    for (Nx = output->islower; Nx <= output->isupper; Nx += output->isstep)
      ;
    for (Ny = output->jslower; Ny <= output->jsupper; Ny += output->jsstep)
      ;
    const int Nz = output->atm.nzout;
    const int Nv = output->wl.nlambda_h_print3D;
    if (doflx == 1) {
      fprintf (stderr, " ... writing flx output to <%s> shape (Nx %d Ny %d Nz %d Nwvl %d)\n", flxfilename, Nx, Ny, Nz, Nv);
    }
    fprintf (stderr,
             " ... writing rad output to <%s> shape (Nx %d Ny %d Nz %d Nwvl %d Nstokes %d Nalis %d)\n",
             radfilename,
             Nx,
             Ny,
             Nz,
             Nv,
             input.rte.mc.nstokes,
             output->mc.alis.Nc);
  }

  if (doflx == 1) {
    if ((fflx = fopen (flxfilename, "w")) == NULL) {
      fprintf (stderr, "Error opening %s for writing\n", flxfilename);
      CHKPOINTER (fflx);
    }

    if (input.rte.mc.concentration_is || input.rte.mc.spectral_is) {
      if ((fflxis = fopen (flxisfilename, "w")) == NULL) {
        fprintf (stderr, "Error opening %s for writing\n", flxisfilename);
        CHKPOINTER (fflxis);
      }
    }
  }

  if ((frad = fopen (radfilename, "w")) == NULL) {
    fprintf (stderr, "Error opening %s for writing\n", radfilename);
    CHKPOINTER (frad);
  }

  if (input.rte.mc.jacobian[DIM_1D]) {
    if ((fjac = fopen (jacfilename, "w")) == NULL) {
      fprintf (stderr, "Error opening %s for writing\n", jacfilename);
      CHKPOINTER (fjac);
    }
  }

  /* variances */
  if (input.rte.mc.std) {
    if (doflx == 1) {
      if ((fflxvar = fopen (flxvarfilename, "w")) == NULL) {
        fprintf (stderr, "Error opening %s for writing\n", flxvarfilename);
        CHKPOINTER (fflxvar);
      }
    }

    if ((fradvar = fopen (radvarfilename, "w")) == NULL) {
      fprintf (stderr, "Error opening %s for writing\n", radvarfilename);
      CHKPOINTER (fradvar);
    }
  }

  for (int iv = 0; iv < output->wl.nlambda_h_print3D; iv++)
    for (int ks = 0; ks < output->atm.nzout; ks++)
      for (int is = output->islower; is <= output->isupper; is += output->isstep)
        for (int js = output->jslower; js <= output->jsupper; js += output->jsstep) {

          if (doflx == 1) {
            fprintf (fflx,
                     "%9.5f %4d %4d %4d %.8e %.8e %.8e %.8e %.8e %.8e\n",
                     output->wl.lambda_h[iv],
                     is,
                     js,
                     ks,
                     output->rfldir3d[ks][is][js][iv],
                     output->rfldn3d[ks][is][js][iv],
                     output->flup3d[ks][is][js][iv],
                     output->uavgso3d[ks][is][js][iv],
                     output->uavgdn3d[ks][is][js][iv],
                     output->uavgup3d[ks][is][js][iv]);

            if (input.rte.mc.concentration_is || input.rte.mc.spectral_is) {
              for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
                fprintf (fflxis,
                         "%9.5f %4d %4d %4d %4d %g\n",
                         output->wl.lambda_h[iv],
                         ic,
                         is,
                         js,
                         ks,
                         output->fl3d_is[ks][is][js][ic][iv]);
              }
            }
          }

          /* FIXCE reformat output for concentration_is ????  */
          if (output->radiance3d) {
            for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
              for (int ic = 0; ic < output->mc.alis.Nc; ic++) {
                fprintf (frad,
                         "%9.5f %4d %4d %4d %g\n",
                         output->wl.lambda_h[iv],
                         is,
                         js,
                         ks,
                         output->radiance3d[ks][is][js][ip][ic][iv]);
              }
            }
          }

          if (input.rte.mc.jacobian[DIM_1D]) {
            for (int kc = 0; kc < output->atm.nlev - 1; kc++) {
              fprintf (fjac, "%9.5f %4d %4d %4d %4d ", output->wl.lambda_h[iv], is, js, ks, kc);
              for (int isp = 0; isp < input.n_caoth + 2; isp++) {
                fprintf (fjac,
                         "%.6e %.6e ",
                         output->jacobian[ks][is][js][isp][0][kc][iv],
                         output->jacobian[ks][is][js][isp][1][kc][iv]);
              }
              fprintf (fjac, "\n");
            }
          }

          /* variances */
          if (input.rte.mc.std) {

            if (doflx == 1)
              fprintf (fflxvar,
                       "%9.5f %4d %4d %4d %.8e %.8e %.8e %.8e %.8e %.8e\n",
                       output->wl.lambda_h[iv],
                       is,
                       js,
                       ks,
                       sqrt (output->rfldir3d_var[ks][is][js][iv]),
                       sqrt (output->rfldn3d_var[ks][is][js][iv]),
                       sqrt (output->flup3d_var[ks][is][js][iv]),
                       sqrt (output->uavgso3d_var[ks][is][js][iv]),
                       sqrt (output->uavgdn3d_var[ks][is][js][iv]),
                       sqrt (output->uavgup3d_var[ks][is][js][iv]));

            if (output->radiance3d_var) {
              for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
                fprintf (fradvar,
                         "%9.5f %4d %4d %4d %g\n",
                         output->wl.lambda_h[iv],
                         is,
                         js,
                         ks,
                         sqrt (output->radiance3d_var[ks][is][js][ip][iv]));
              }
            }
          }
        }

  if (doflx == 1)
    fclose (fflx);
  fclose (frad);

  if (input.rte.mc.concentration_is)
    fclose (fflxis);

  if (input.rte.mc.jacobian[DIM_1D])
    fclose (fjac);

  if (input.rte.mc.std) {
    if (doflx == 1)
      fclose (fflxvar);
    fclose (fradvar);
  }

  if (input.output3D && ((input.rte.mc.absorption != MCFORWARD_ABS_NONE) || input.rte.mc.backward.absorption ||
                         (input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) || ((input.ipa3d > 0) && (input.heating != HEAT_NONE)))) {

    if ((fabs = fopen (absfilename, "w")) == NULL) {
      fprintf (stderr, "Error opening %s for writing\n", absfilename);
      CHKPOINTER (fabs);
    }

    /* 27.02.2013 **CK **BM: added for thermal backward heating rates std */
    if (input.rte.mc.std)
      if ((fabsvar = fopen (absvarfilename, "w")) == NULL) {
        fprintf (stderr, "Error opening %s for writing\n", radvarfilename);
        CHKERR (-1);
      }

    if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
        input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {

      if (!input.quiet) {
        fprintf (stderr,
                 " ... writing absorption to <%s> shape (Nx %d Ny %d Nz %d Nwvl %d)\n",
                 absfilename,
                 output->atm.Nxcld,
                 output->atm.Nycld,
                 output->atm.nthreed,
                 output->wl.nlambda_h_print3D);
      }

      size_t k3d = 0;  // count 3d layers, because output->abs3d_r has only entries for 3d layers
      for (int ks = 0; ks < output->atm.Nzcld; ks++) {
        if (output->atm.threed[ks]) { /* only for 3D layers, BM07122005 */
          for (int is = 0; is < output->atm.Nxcld; is++) {
            for (int js = 0; js < output->atm.Nycld; js++) {
              for (int iv = 0; iv < output->wl.nlambda_h_print3D; iv++) {
                fprintf (fabs, "%9.5f %4d %4d %4d %.8e\n", output->wl.lambda_h[iv], is, js, ks, output->abs3d[k3d][is][js][iv]);
                /* 27.02.2013 **CK **BM: added for thermal backward heating rates std */
                if (input.rte.mc.std) {
                  fprintf (fabsvar,
                           "%9.5f %4d %4d %4d %.8e\n",
                           output->wl.lambda_h[iv],
                           is,
                           js,
                           ks,
                           sqrt (output->abs3d_var[k3d][is][js][iv]));
                }
              }
            }
          }
          k3d++;
        }
      }
    }

    if (input.rte.mc.backward.absorption && !(input.ipa3d > 0)) {
      if (!input.quiet) {
        int Nx, Ny;
        for (Nx = output->islower; Nx <= output->isupper; Nx += output->isstep)
          ;
        for (Ny = output->jslower; Ny <= output->jsupper; Ny += output->jsstep)
          ;
        fprintf (stderr,
                 " ... writing absorption to <%s> shape (Nx %d Ny %d Nz %d Nwvl %d)\n",
                 absfilename,
                 Nx,
                 Ny,
                 output->atm.nzout,
                 output->wl.nlambda_h_print3D);
      }
      for (int iv = 0; iv < output->wl.nlambda_h_print3D; iv++) {
        for (int ks = 0; ks < output->atm.nzout; ks++) {
          for (int is = output->islower; is <= output->isupper; is += output->isstep) {
            for (int js = output->jslower; js <= output->jsupper; js += output->jsstep) {
              fprintf (fabs, "%9.5f %4d %4d %4d %.8e\n", output->wl.lambda_h[iv], is, js, ks, output->absback3d[ks][is][js][iv]);

              /* 27.02.2013 **CK **BM: added for thermal backward heating rates std */
              if (input.rte.mc.std) {
                fprintf (fabsvar,
                         "%9.5f %4d %4d %4d %.8e\n",
                         output->wl.lambda_h[iv],
                         is,
                         js,
                         ks,
                         sqrt (output->absback3d_var[ks][is][js][iv]));
              }
            }
          }
        }
      }
    }

    fclose (fabs);

    /* 27.02.2013 **CK **BM: added for thermal backward heating rates std */
    if (input.rte.mc.std)
      fclose (fabsvar);
  }

  if (!input.quiet) {
    const time_t endtime  = time (NULL);
    const int    timediff = difftime (endtime, starttime);
    fprintf (stderr, " ... used %d seconds for write_spectral3D\n", timediff);
  }

  return 0;
}

static int write_spectral3D_nc (input_struct input, output_struct* output) {
  const time_t starttime              = time (NULL);
  char         filename[FILENAME_MAX] = "";
  const char   extension[4]           = ".nc";
  int          status;

  strcpy (filename, input.rte.mc.filename[FN_MC_BASENAME]);

  strcat (filename, ".spc");

  strcat (filename, extension);

  int ncid;
  status = create_new_netcdf (filename, &ncid);
  CHKERR (status);

  {
    char       datestr[100];
    struct tm* t = localtime (&starttime);
    strftime (datestr, sizeof (datestr), "%Y-%m-%d_%H:%M", t);
    status = write_netcdf_attribute_text (ncid, NULL, "created_time", datestr);
    CHKERR (status);
  }

  {
    char* options_stash;
    status = options_stash_get (&options_stash);
    CHKERR (status);
    status = write_netcdf_attribute_text (ncid, NULL, "options", options_stash);
    CHKERR (status);
  }

  int Nx = 0;
  int Ny = 0;
  for (int i = output->islower; i <= output->isupper; i += output->isstep, ++Nx)
    ;
  for (int j = output->jslower; j <= output->jsupper; j += output->jsstep, ++Ny)
    ;
  const int Nz = output->atm.nzout;
  const int Nv = output->wl.nlambda_h_print3D;

  if (!input.quiet) {
    fprintf (stderr, " ... writing netcdf flx output to <%s>\n", filename);
    fprintf (stderr, " ... array shape Nx %d Ny %d Nz %d Nwvl %d\n", Nx, Ny, Nz, Nv);
  }

  {
    float zlev[Nz];
    for (size_t k = 0; k < Nz; ++k) {
      zlev[k] = output->atm.zout_sur[k] * 1e3;
    }
    status = write_netcdf_1D_float (ncid, zlev, Nz, "zlev", "zlev", NULL);
    CHKERR (status);
  }
  {
    status = write_netcdf_1D_float (ncid, output->wl.lambda_h, Nv, "wvl", "wvl", NULL);
    CHKERR (status);
  }
  {
    float xdim[Nx];
    for (size_t i = output->islower, cnt = 0; i < output->isupper + 1; i += output->isstep, ++cnt) {
      xdim[cnt] = (i + 0.5) * output->mc.sample.delX;
    }
    status = write_netcdf_1D_float (ncid, xdim, Nx, "x", "x", NULL);
    CHKERR (status);
  }
  {
    float ydim[Ny];
    for (size_t j = output->jslower, cnt = 0; j < output->jsupper + 1; j += output->jsstep, ++cnt) {
      ydim[cnt] = (j + 0.5) * output->mc.sample.delX;
    }
    status = write_netcdf_1D_float (ncid, ydim, Ny, "y", "y", NULL);
    CHKERR (status);
  }

  if (input.output3D) {
    size_t slice_z[3] = {0, Nz, 1};
    size_t slice_x[3] = {output->islower, output->isupper + 1, output->isstep};
    size_t slice_y[3] = {output->jslower, output->jsupper + 1, output->jsstep};
    size_t slice_v[3] = {0, Nv, 1};
    status            = write_netcdf_4D_float_sliced (ncid,
                                           output->rfldir3d,
                                           slice_z,
                                           slice_x,
                                           slice_y,
                                           slice_v,
                                           "edir",
                                           "zlev",
                                           "x",
                                           "y",
                                           "wvl",
                                           NULL);
    CHKERR (status);
    status = write_netcdf_4D_float_sliced (ncid,
                                           output->rfldn3d,
                                           slice_z,
                                           slice_x,
                                           slice_y,
                                           slice_v,
                                           "edn",
                                           "zlev",
                                           "x",
                                           "y",
                                           "wvl",
                                           NULL);
    CHKERR (status);
    status =
      write_netcdf_4D_float_sliced (ncid, output->flup3d, slice_z, slice_x, slice_y, slice_v, "eup", "zlev", "x", "y", "wvl", NULL);
    CHKERR (status);
    status = write_netcdf_4D_float_sliced (ncid,
                                           output->uavgso3d,
                                           slice_z,
                                           slice_x,
                                           slice_y,
                                           slice_v,
                                           "udir",
                                           "zlev",
                                           "x",
                                           "y",
                                           "wvl",
                                           NULL);
    CHKERR (status);
    status = write_netcdf_4D_float_sliced (ncid,
                                           output->uavgdn3d,
                                           slice_z,
                                           slice_x,
                                           slice_y,
                                           slice_v,
                                           "udn",
                                           "zlev",
                                           "x",
                                           "y",
                                           "wvl",
                                           NULL);
    CHKERR (status);
    status = write_netcdf_4D_float_sliced (ncid,
                                           output->uavgup3d,
                                           slice_z,
                                           slice_x,
                                           slice_y,
                                           slice_v,
                                           "uup",
                                           "zlev",
                                           "x",
                                           "y",
                                           "wvl",
                                           NULL);
    CHKERR (status);

    if (input.rte.mc.std) {
      status = write_netcdf_4D_float_sliced (ncid,
                                             output->rfldir3d_var,
                                             slice_z,
                                             slice_x,
                                             slice_y,
                                             slice_v,
                                             "edir_std",
                                             "zlev",
                                             "x",
                                             "y",
                                             "wvl",
                                             &sqrtf);
      CHKERR (status);
      status = write_netcdf_4D_float_sliced (ncid,
                                             output->rfldn3d_var,
                                             slice_z,
                                             slice_x,
                                             slice_y,
                                             slice_v,
                                             "edn_std",
                                             "zlev",
                                             "x",
                                             "y",
                                             "wvl",
                                             &sqrtf);
      CHKERR (status);
      status = write_netcdf_4D_float_sliced (ncid,
                                             output->flup3d_var,
                                             slice_z,
                                             slice_x,
                                             slice_y,
                                             slice_v,
                                             "eup_std",
                                             "zlev",
                                             "x",
                                             "y",
                                             "wvl",
                                             &sqrtf);
      CHKERR (status);
      status = write_netcdf_4D_float_sliced (ncid,
                                             output->uavgso3d_var,
                                             slice_z,
                                             slice_x,
                                             slice_y,
                                             slice_v,
                                             "udir_std",
                                             "zlev",
                                             "x",
                                             "y",
                                             "wvl",
                                             &sqrtf);
      CHKERR (status);
      status = write_netcdf_4D_float_sliced (ncid,
                                             output->uavgdn3d_var,
                                             slice_z,
                                             slice_x,
                                             slice_y,
                                             slice_v,
                                             "udn_std",
                                             "zlev",
                                             "x",
                                             "y",
                                             "wvl",
                                             &sqrtf);
      CHKERR (status);
      status = write_netcdf_4D_float_sliced (ncid,
                                             output->uavgup3d_var,
                                             slice_z,
                                             slice_x,
                                             slice_y,
                                             slice_v,
                                             "uup_std",
                                             "zlev",
                                             "x",
                                             "y",
                                             "wvl",
                                             &sqrtf);
      CHKERR (status);
    }

    if (output->radiance3d) {
      size_t slice_alis[3]   = {0, output->mc.alis.Nc, 1};
      size_t slice_stokes[3] = {0, input.rte.mc.nstokes, 1};
      status                 = write_netcdf_6D_float_sliced (ncid,
                                             output->radiance3d,
                                             slice_z,
                                             slice_x,
                                             slice_y,
                                             slice_stokes,
                                             slice_alis,
                                             slice_v,
                                             "radiance",
                                             "zlev",
                                             "x",
                                             "y",
                                             "Nstokes",
                                             "Nalis",
                                             "wvl",
                                             NULL);
      CHKERR (status);

      if (output->radiance3d_var) {
        status = write_netcdf_5D_float_sliced (ncid,
                                               output->radiance3d_var,
                                               slice_z,
                                               slice_x,
                                               slice_y,
                                               slice_stokes,
                                               slice_v,
                                               "radiance_std",
                                               "zlev",
                                               "x",
                                               "y",
                                               "Nstokes",
                                               "wvl",
                                               &sqrtf);
        CHKERR (status);
      }
    }
  }

  if (input.output3D && ((input.rte.mc.absorption != MCFORWARD_ABS_NONE) || input.rte.mc.backward.absorption ||
                         (input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) || ((input.ipa3d > 0) && (input.heating != HEAT_NONE)))) {

    if ((input.rte.solver == SOLVER_MONTECARLO && (input.rte.mc.absorption != MCFORWARD_ABS_NONE)) ||
        (input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) || (input.ipa3d > 0)) {
      {
        float  zlyr3D[output->atm.nthreed];  // layer mean heights for 3D layers (only have 3D absorption in 3D layers)
        size_t cnt = 0;
        for (size_t k = 0; k < output->atm.nlyr; ++k) {
          if (output->atm.threed[k]) {
            size_t k_rev  = output->atm.nlyr - 1 - k;
            zlyr3D[cnt++] = (output->atm.zd[k_rev] + output->atm.zd[k_rev + 1]) * .5 * 1e3;
          }
        }
        status = write_netcdf_1D_float (ncid, zlyr3D, cnt, "zlay3D", "zlay3D", NULL);
        CHKERR (status);
      }
      {
        float xdim_abs[output->atm.Nxcld];
        for (size_t i = 0; i < output->atm.Nxcld; ++i) {
          xdim_abs[i] = (i + 0.5) * output->atm.dxcld;
        }
        status = write_netcdf_1D_float (ncid, xdim_abs, output->atm.Nxcld, "x_abs", "x_abs", NULL);
        CHKERR (status);
      }
      {
        float ydim_abs[output->atm.Nycld];
        for (size_t j = 0; j < output->atm.Nycld; ++j) {
          ydim_abs[j] = (j + 0.5) * output->atm.dycld;
        }
        status = write_netcdf_1D_float (ncid, ydim_abs, output->atm.Nycld, "y_abs", "y_abs", NULL);
        CHKERR (status);
      }

      status = write_netcdf_4D_float (ncid,
                                      output->abs3d,
                                      output->atm.nthreed,
                                      output->atm.Nxcld,
                                      output->atm.Nycld,
                                      output->wl.nlambda_h_print3D,
                                      "absorption",
                                      "zlay3D",
                                      "x_abs",
                                      "y_abs",
                                      "wvl",
                                      NULL);
      CHKERR (status);

      if (input.rte.mc.std) {
        printf ("abso_std start\n");
        status = write_netcdf_4D_float (ncid,
                                        output->abs3d_var,
                                        output->atm.nthreed,
                                        output->atm.Nxcld,
                                        output->atm.Nycld,
                                        output->wl.nlambda_h_print3D,
                                        "absorption_std",
                                        "zlay3D",
                                        "x_abs",
                                        "y_abs",
                                        "wvl",
                                        &sqrtf);
        CHKERR (status);
        printf ("abso_std end\n");
      }
    }

    if (input.rte.mc.backward.absorption && !(input.ipa3d > 0)) {
      float zlyr3D[output->atm.nzout];
      for (size_t k = 0; k < output->atm.nzout; ++k) {
        float zbot, ztop;
        status    = mc_backward_absorption_find_vert_index (&output->atm, output->atm.zout_sur[k], NULL, NULL, &zbot, &ztop);
        zlyr3D[k] = (zbot + ztop) * .5 * 1e3;
      }
      status = write_netcdf_1D_float (ncid, zlyr3D, output->atm.nzout, "zlay3D", "zlay3D", NULL);
      CHKERR (status);
      //TODO: 06.10.2023, P.G.: x_abs, y_abs here as well? is Nx and Ny in this case appropriate?

      status = write_netcdf_4D_float (ncid,
                                      output->absback3d,
                                      output->atm.nzout,
                                      Nx,
                                      Ny,
                                      output->wl.nlambda_h_print3D,
                                      "absorption",
                                      "zlay3D",
                                      "x",
                                      "y",
                                      "wvl",
                                      NULL);
      CHKERR (status);
      if (input.rte.mc.std) {
        status = write_netcdf_4D_float (ncid,
                                        output->absback3d_var,
                                        output->atm.nzout,
                                        Nx,
                                        Ny,
                                        output->wl.nlambda_h_print3D,
                                        "absorption_std",
                                        "zlay3D",
                                        "x",
                                        "y",
                                        "wvl",
                                        &sqrtf);
        CHKERR (status);
      }
    }
  }

  status = close_netcdf (ncid);
  CHKERR (status);

  if (!input.quiet) {
    const time_t endtime  = time (NULL);
    const int    timediff = difftime (endtime, starttime);
    fprintf (stderr, " ... used %d seconds for write_spectral3D_nc\n", timediff);
  }

  return 0;
}

/***********************************************************************************/
/* Sum the 3D results over wavelength and write data to the first array element    */
/***********************************************************************************/

int sum3D (input_struct input, output_struct* output) {
  int    iv = 0, is = 0, js = 0, ks = 0, ip = 0, ic = 0, isp = 0, i = 0, lc = 0;
  double ffactor = 0, ffactor2 = 0, rfactor = 0, rfactor2 = 0;
  float  incident = 0, fbeam = 0;
  int    status = 0;

  /* calculate incident flux (incident) and sum extraterrestrial irradiance (fbeam) */
  for (iv = 0; iv < output->wl.nlambda_h; iv++) {
    incident += output->wl.filter[iv] * output->wl.fbeam[iv] * cos (output->sza_h[iv] * M_PI / 180.0);
    fbeam += output->wl.filter[iv] * output->wl.fbeam[iv];
  }

  /* default scaling factors                   */
  ffactor = 1.0; /* irradiance multiplicator */
  rfactor = 1.0; /* radiance multiplicator   */

  /* if transmittance, then divide by the extraterrestrial flux */
  if (input.source != SRC_THERMAL) {
    switch (input.calibration) {
    case OUTCAL_ABSOLUTE:
      break;

    case OUTCAL_TRANSMITTANCE:
    case OUTCAL_BRIGHTNESS:
      ffactor = 1.0 / fbeam; /* irradiance multiplicator */
      rfactor = 1.0 / fbeam; /* radiance multiplicator   */
      break;

    case OUTCAL_REFLECTIVITY:

      ffactor = 1.0 / incident;
      rfactor = M_PI / incident;
      break;

    default:
      fprintf (stderr, "Error, unknown output calibration %d\n", input.calibration);
      return -1;
    }
  }

  /* fluxes and radiances */
  for (ks = 0; ks < output->atm.nzout; ks++)
    for (is = output->islower; is <= output->isupper; is += output->isstep)
      for (js = output->jslower; js <= output->jsupper; js += output->jsstep) {

        for (iv = 1; iv < output->wl.nlambda_h; iv++) {
          output->rfldir3d[ks][is][js][0] += output->rfldir3d[ks][is][js][iv];
          output->rfldn3d[ks][is][js][0] += output->rfldn3d[ks][is][js][iv];
          output->flup3d[ks][is][js][0] += output->flup3d[ks][is][js][iv];
          output->uavgso3d[ks][is][js][0] += output->uavgso3d[ks][is][js][iv];
          output->uavgdn3d[ks][is][js][0] += output->uavgdn3d[ks][is][js][iv];
          output->uavgup3d[ks][is][js][0] += output->uavgup3d[ks][is][js][iv];

          if (output->radiance3d)
            for (ip = 0; ip < input.rte.mc.nstokes; ip++)
              for (ic = 0; ic < output->mc.alis.Nc; ic++)
                output->radiance3d[ks][is][js][ip][ic][0] += output->radiance3d[ks][is][js][ip][ic][iv];

          if (input.rte.mc.concentration_is || input.rte.mc.spectral_is)
            for (ic = 0; ic < output->mc.alis.Nc; ic++) {
              /* fprintf(stderr, "sum3D %f is %f \n", output->rfldn3d[ks][is][js][iv], output->fl3d_is [ks][is][js][ic][iv]); */
              output->fl3d_is[ks][is][js][ic][0] += output->fl3d_is[ks][is][js][ic][iv];
            }

          if (input.rte.mc.jacobian[DIM_1D])
            for (isp = 0; isp < input.n_caoth + 2; isp++)
              for (i = 0; i < 2; i++)
                for (lc = 0; lc < output->atm.nlev - 1; lc++) {
                  /* fprintf(stderr, "sum3D %f \n", output->jacobian [ks][is][js][isp][i][lc][iv]); */
                  output->jacobian[ks][is][js][isp][i][lc][0] += output->jacobian[ks][is][js][isp][i][lc][iv];
                }

          if (input.rte.mc.backward.absorption)
            output->absback3d[ks][is][js][0] += output->absback3d[ks][is][js][iv];

          /* variances */
          if (input.rte.mc.std) {

            output->rfldir3d_var[ks][is][js][0] += output->rfldir3d_var[ks][is][js][iv];
            output->rfldn3d_var[ks][is][js][0] += output->rfldn3d_var[ks][is][js][iv];
            output->flup3d_var[ks][is][js][0] += output->flup3d_var[ks][is][js][iv];
            output->uavgso3d_var[ks][is][js][0] += output->uavgso3d_var[ks][is][js][iv];
            output->uavgdn3d_var[ks][is][js][0] += output->uavgdn3d_var[ks][is][js][iv];
            output->uavgup3d_var[ks][is][js][0] += output->uavgup3d_var[ks][is][js][iv];

            if (output->radiance3d_var)
              for (ip = 0; ip < input.rte.mc.nstokes; ip++)
                output->radiance3d_var[ks][is][js][ip][0] += output->radiance3d_var[ks][is][js][ip][iv];

            if (input.rte.mc.backward.absorption) {
              output->absback3d_var[ks][is][js][0] += output->absback3d_var[ks][is][js][iv];
            }
          }
        }

        /* calibrate */
        if (output->wl.nlambda_h > 0) {
          output->rfldir3d[ks][is][js][0] *= ffactor;
          output->rfldn3d[ks][is][js][0] *= ffactor;
          output->flup3d[ks][is][js][0] *= ffactor;
          output->uavgso3d[ks][is][js][0] *= ffactor;
          output->uavgdn3d[ks][is][js][0] *= ffactor;
          output->uavgup3d[ks][is][js][0] *= ffactor;

          if (input.rte.mc.concentration_is || input.rte.mc.spectral_is)
            for (ic = 0; ic < output->mc.alis.Nc; ic++)
              output->fl3d_is[ks][is][js][ic][0] *= ffactor;

          if (output->radiance3d)
            for (ip = 0; ip < input.rte.mc.nstokes; ip++)
              for (ic = 0; ic < output->mc.alis.Nc; ic++)
                output->radiance3d[ks][is][js][ip][ic][0] *= rfactor;

          if (input.rte.mc.jacobian[DIM_1D]) {
            for (isp = 0; isp < input.n_caoth + 2; isp++)
              for (i = 0; i < 2; i++)
                for (lc = 0; lc < output->atm.nlev - 1; lc++)
                  if (input.rte.mc.backward.output == MCBACKWARD_EDN || input.rte.mc.backward.output == MCBACKWARD_EDNPV ||
                      input.rte.mc.backward.output == MCBACKWARD_EUP)
                    output->jacobian[ks][is][js][isp][i][lc][0] *= ffactor;
                  else  // radiance
                    output->jacobian[ks][is][js][isp][i][lc][0] *= rfactor;
          }

          if (input.rte.mc.backward.absorption)
            output->absback3d[ks][is][js][0] *= ffactor;

          /* variances */
          if (input.rte.mc.std) {

            ffactor2 = ffactor * ffactor;
            rfactor2 = rfactor * rfactor;

            output->rfldir3d_var[ks][is][js][0] *= ffactor2;
            output->rfldn3d_var[ks][is][js][0] *= ffactor2;
            output->flup3d_var[ks][is][js][0] *= ffactor2;
            output->uavgso3d_var[ks][is][js][0] *= ffactor2;
            output->uavgdn3d_var[ks][is][js][0] *= ffactor2;
            output->uavgup3d_var[ks][is][js][0] *= ffactor2;

            if (output->radiance3d_var)
              for (ip = 0; ip < input.rte.mc.nstokes; ip++)
                output->radiance3d_var[ks][is][js][ip][0] *= rfactor2;

            if (input.rte.mc.backward.absorption)
              output->absback3d_var[ks][is][js][0] *= ffactor2;
          }
        }
      }

  if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
      input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
    const size_t N3d_layers = vec_count_true (output->atm.Nzcld, output->atm.threed);
    for (size_t k = 0; k < N3d_layers; k++) {
      for (is = 0; is < output->atm.Nxcld; is++) {
        for (js = 0; js < output->atm.Nycld; js++) {

          for (iv = 1; iv < output->wl.nlambda_h; iv++) {
            output->abs3d[k][is][js][0] += output->abs3d[k][is][js][iv];
            if (input.rte.mc.std) {
              output->abs3d_var[k][is][js][0] += output->abs3d_var[k][is][js][iv];
            }
          }
          output->abs3d[k][is][js][0] *= ffactor;
          if (input.rte.mc.std) {
            output->abs3d_var[k][is][js][0] *= ffactor2;
          }
        }
      }
    }
  }

  /* finally, convert to brightness temperature if requested */
  if (input.source == SRC_THERMAL && input.calibration == OUTCAL_BRIGHTNESS) {

    if (!input.quiet)
      fprintf (stderr, " ... converting 3D radiances and fluxes to brightness temperatures\n");

    /* convert irradiances / radiances to brightness temperatures */
    iv     = 0;
    status = output2bt (input, output, iv, 1);
    if (status != 0) {
      fprintf (stderr, "Error %d returned by output2bt()\n", status);
      CHKERR (status);
    }
  }

  if (output->triangle_results_o) {  // sum up triangle_results
    for (iv = 1; iv < output->wl.nlambda_h; iv++) {
      status = add_triangular_surface_result (ffactor, output->triangle_results_o[iv], output->triangle_results_o[0]);
      CHKERR (status);
    }
  }
  return status;
}

/**************************************************************************************/
/* Integrate the 3D results over wavelength and write data to the first array element */
/**************************************************************************************/

int integrate3D (input_struct input, output_struct* output) {
  int     iv = 0, is = 0, js = 0, ks = 0, ip = 0, ic = 0, isp = 0, i = 0, lc = 0;
  int     status = 0, used_unit = 0;
  double* cos_SZA = NULL;
  double  ffactor = 0, ffactor2 = 0, rfactor = 0, rfactor2 = 0;
  float * xint = NULL, *yint = NULL;
  float   incident = 0, fbeam = 0;

  /* calculate integrated incident flux */
  xint = calloc (output->wl.nlambda_h, sizeof (float));
  CHKPOINTER (xint);
  yint = calloc (output->wl.nlambda_h, sizeof (float));
  CHKPOINTER (yint);

  cos_SZA = calloc (output->wl.nlambda_h, sizeof (double));
  CHKPOINTER (cos_SZA);

  switch (input.source) {
  case SRC_SOLAR:
  case SRC_LIDAR: /* BCA */
  case SRC_BLITZ: /* BCA */
    for (iv = 0; iv < output->wl.nlambda_h; iv++)
      cos_SZA[iv] = cos (output->sza_h[iv] * M_PI / 180.0);
    break;
  case SRC_THERMAL:
    for (iv = 0; iv < output->wl.nlambda_h; iv++)
      cos_SZA[iv] = 1.0;
    break;
  default:
    fprintf (stderr, "Error, unknown source %d\n", input.source);
    CHKERR (-1);
  }

  used_unit = input.output_unit;
  if (used_unit == UNIT_NOT_DEFINED) /* if no output unit is specified, assume the same units as the input spectrum */
    used_unit = output->spectrum_unit;

  switch (used_unit) {
  case UNIT_PER_NM:
    if (input.verbose)
      fprintf (stderr, " *** integration in wavelength space \n");
    for (iv = 0; iv < output->wl.nlambda_h; iv++) {
      xint[iv] = (double) output->wl.lambda_h[iv];
      yint[iv] = output->wl.filter[iv] * output->wl.fbeam[iv] * cos_SZA[iv];
    }
    break;
  case UNIT_PER_CM_1:
    /* integration over wavenumber k, minus in order to get ascending wavenumbers */
    if (input.verbose)
      fprintf (stderr, " *** integration in wavenumber space \n");
    for (iv = 0; iv < output->wl.nlambda_h; iv++) {
      xint[iv] = -(double) 1.0e+7 / output->wl.lambda_h[iv]; /* k = 10**7 / lambda */ /* 10**7 == nm -> cm */
      yint[iv] = output->wl.filter[iv] * output->wl.fbeam[iv] * cos_SZA[iv];
    }
    break;
  case UNIT_PER_BAND:
    fprintf (stderr, "Error, combination 'output_process integrate' and 'output_process per_band',\n");
    fprintf (stderr, "    or combination 'output_process integrate' and an input spectrum defined in \n");
    fprintf (stderr, "       (wavelength or wavenumber) bands does not make sense \n");
    fprintf (stderr, "       please use 'output_process sum', when dealing with band parametrisations. \n\n");
    CHKERR (-1);
    break;
  case UNIT_NOT_DEFINED:
    fprintf (stderr, "Error, in order to use 'output_process integrate' it is nessesary to specify the \n");
    fprintf (stderr, "       unit of the extraterrestrial spectrum with 'source solar filename unit' or\n");
    fprintf (stderr, "       unit of the output with 'output_process per_nm' or 'output_process per_cm-1' \n\n");
    CHKERR (-1);
    break;
  default:
    fprintf (stderr,
             "Error: Program bug, unsupported unit of extraterrestrial flux %d or output unit %d\n",
             output->spectrum_unit,
             input.output_unit);
    CHKERR (-1);
  }

  incident = integrate_float (xint, yint, output->wl.nlambda_h);

  /* calculate integrated extraterrestrial irradiance */
  for (iv = 0; iv < output->wl.nlambda_h; iv++)
    yint[iv] = output->wl.filter[iv] * output->wl.fbeam[iv];

  fbeam = integrate_float (xint, yint, output->wl.nlambda_h);

  /* default scaling factors                   */
  ffactor = 1.0; /* irradiance multiplicator */
  rfactor = 1.0; /* radiance multiplicator   */

  /* if transmittance, then divide by the extraterrestrial flux */
  if (input.source != SRC_THERMAL) {
    switch (input.calibration) {
    case OUTCAL_ABSOLUTE:
      break;

    case OUTCAL_TRANSMITTANCE:
    case OUTCAL_BRIGHTNESS:
      ffactor = 1.0 / fbeam; /* irradiance multiplicator */
      rfactor = 1.0 / fbeam; /* radiance multiplicator   */
      break;

    case OUTCAL_REFLECTIVITY:

      ffactor = 1.0 / incident;
      rfactor = M_PI / incident;
      break;

    default:
      fprintf (stderr, "Error, unknown output calibration %d\n", input.calibration);
      CHKERR (-1);
    }
  }

  /* fluxes and radiances */
  for (ks = 0; ks < output->atm.nzout; ks++) {
    for (is = output->islower; is <= output->isupper; is += output->isstep) {
      for (js = output->jslower; js <= output->jsupper; js += output->jsstep) {

        output->rfldir3d[ks][is][js][0] = integrate_float (xint, output->rfldir3d[ks][is][js], output->wl.nlambda_h);
        output->rfldn3d[ks][is][js][0]  = integrate_float (xint, output->rfldn3d[ks][is][js], output->wl.nlambda_h);
        output->flup3d[ks][is][js][0]   = integrate_float (xint, output->flup3d[ks][is][js], output->wl.nlambda_h);
        output->uavgso3d[ks][is][js][0] = integrate_float (xint, output->uavgso3d[ks][is][js], output->wl.nlambda_h);
        output->uavgdn3d[ks][is][js][0] = integrate_float (xint, output->uavgdn3d[ks][is][js], output->wl.nlambda_h);
        output->uavgup3d[ks][is][js][0] = integrate_float (xint, output->uavgup3d[ks][is][js], output->wl.nlambda_h);

        if (input.rte.mc.concentration_is || input.rte.mc.spectral_is)
          for (ic = 0; ic < output->mc.alis.Nc; ic++)
            output->fl3d_is[ks][is][js][ic][0] = integrate_float (xint, output->fl3d_is[ks][is][js][ic], output->wl.nlambda_h);

        if (output->radiance3d)
          for (ip = 0; ip < input.rte.mc.nstokes; ip++)
            for (ic = 0; ic < output->mc.alis.Nc; ic++)
              output->radiance3d[ks][is][js][ip][ic][0] =
                integrate_float (xint, output->radiance3d[ks][is][js][ip][ic], output->wl.nlambda_h);

        if (input.rte.mc.jacobian[DIM_1D])
          for (isp = 0; isp < input.n_caoth + 2; isp++)
            for (i = 0; i < 2; i++)
              for (lc = 0; lc < output->atm.nlev - 1; lc++)
                output->jacobian[ks][is][js][isp][i][lc][0] =
                  integrate_float (xint, output->jacobian[ks][is][js][isp][i][lc], output->wl.nlambda_h);

        if (input.rte.mc.backward.absorption)
          output->absback3d[ks][is][js][0] = integrate_float (xint, output->absback3d[ks][is][js], output->wl.nlambda_h);

        /* variances */
        if (input.rte.mc.std) {

          output->rfldir3d_var[ks][is][js][0] = integrate_variance_float (xint, output->rfldir3d_var[ks][is][js], output->wl.nlambda_h);
          output->rfldn3d_var[ks][is][js][0]  = integrate_variance_float (xint, output->rfldn3d_var[ks][is][js], output->wl.nlambda_h);
          output->flup3d_var[ks][is][js][0]   = integrate_variance_float (xint, output->flup3d_var[ks][is][js], output->wl.nlambda_h);
          output->uavgso3d_var[ks][is][js][0] = integrate_variance_float (xint, output->uavgso3d_var[ks][is][js], output->wl.nlambda_h);
          output->uavgdn3d_var[ks][is][js][0] = integrate_variance_float (xint, output->uavgdn3d_var[ks][is][js], output->wl.nlambda_h);
          output->uavgup3d_var[ks][is][js][0] = integrate_variance_float (xint, output->uavgup3d_var[ks][is][js], output->wl.nlambda_h);

          if (output->radiance3d_var)
            for (ip = 0; ip < input.rte.mc.nstokes; ip++)
              output->radiance3d_var[ks][is][js][ip][0] =
                integrate_variance_float (xint, output->radiance3d_var[ks][is][js][ip], output->wl.nlambda_h);

          if (input.rte.mc.backward.absorption)
            output->absback3d_var[ks][is][js][0] = integrate_variance_float (xint, output->absback3d_var[ks][is][js], output->wl.nlambda_h);
        }

        /* calibrate */
        if (output->wl.nlambda_h > 0) {
          output->rfldir3d[ks][is][js][0] *= ffactor;
          output->rfldn3d[ks][is][js][0]  *= ffactor;
          output->flup3d[ks][is][js][0]   *= ffactor;
          output->uavgso3d[ks][is][js][0] *= ffactor;
          output->uavgdn3d[ks][is][js][0] *= ffactor;
          output->uavgup3d[ks][is][js][0] *= ffactor;

          if (input.rte.mc.concentration_is || input.rte.mc.spectral_is)
            for (ic = 0; ic < output->mc.alis.Nc; ic++)
              output->fl3d_is[ks][is][js][ic][0] *= ffactor;

          if (output->radiance3d)
            for (ip = 0; ip < input.rte.mc.nstokes; ip++)
              for (ic = 0; ic < output->mc.alis.Nc; ic++)
                output->radiance3d[ks][is][js][ip][ic][0] *= rfactor;

          if (input.rte.mc.jacobian[DIM_1D]) {
            for (isp = 0; isp < input.n_caoth + 2; isp++)
              for (i = 0; i < 2; i++)
                for (lc = 0; lc < output->atm.nlev - 1; lc++)
                  if (input.rte.mc.backward.output == MCBACKWARD_EDN || input.rte.mc.backward.output == MCBACKWARD_EDNPV ||
                      input.rte.mc.backward.output == MCBACKWARD_EUP)
                    output->jacobian[ks][is][js][isp][i][lc][0] *= ffactor;
                  else
                    output->jacobian[ks][is][js][isp][i][lc][0] *= rfactor;
          }

          if (input.rte.mc.backward.absorption)
            output->absback3d[ks][is][js][0] *= ffactor;

          /* variances */
          if (input.rte.mc.std) {
            ffactor2 = ffactor * ffactor;
            rfactor2 = rfactor * rfactor;

            output->rfldir3d_var[ks][is][js][0] *= ffactor2;
            output->rfldn3d_var[ks][is][js][0]  *= ffactor2;
            output->flup3d_var[ks][is][js][0]   *= ffactor2;
            output->uavgso3d_var[ks][is][js][0] *= ffactor2;
            output->uavgdn3d_var[ks][is][js][0] *= ffactor2;
            output->uavgup3d_var[ks][is][js][0] *= ffactor2;

            if (output->radiance3d_var)
              for (ip = 0; ip < input.rte.mc.nstokes; ip++)
                output->radiance3d_var[ks][is][js][ip][0] *= rfactor2;

            if (input.rte.mc.backward.absorption)
              output->absback3d_var[ks][is][js][0] *= ffactor2;
          }
        }
      }
    }
  }

  if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
      input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
    const size_t N3d_layers = vec_count_true (output->atm.Nzcld, output->atm.threed);
    for (size_t k = 0; k < N3d_layers; k++) {
      for (is = 0; is < output->atm.Nxcld; is++) {
        for (js = 0; js < output->atm.Nycld; js++) {

          for (iv = 1; iv < output->wl.nlambda_h; iv++) {
            output->abs3d[k][is][js][0]     = integrate_float (xint, output->abs3d[k][is][js], output->wl.nlambda_h);
            output->abs3d_var[k][is][js][0] = integrate_variance_float (xint, output->abs3d_var[k][is][js], output->wl.nlambda_h);
          }

          if (output->wl.nlambda_h > 0) {
            output->abs3d[k][is][js][0] *= ffactor;
            output->abs3d_var[k][is][js][0] *= ffactor2;
          }
        }
      }
    }
  }

  /* finally, convert to brightness temperature if requested */
  if (input.source == SRC_THERMAL && input.calibration == OUTCAL_BRIGHTNESS) {

    if (!input.quiet)
      fprintf (stderr, " ... converting 3D radiances and fluxes to brightness temperatures\n");

    /* convert irradiances / radiances to brightness temperatures */
    status = output2bt (input, output, 0, 1);
    CHKERR (status);
  }

  if (output->triangle_results_o) {
    CHKERROUT (-1, "Integrate not yet implemented for triangles");
  }

  return status;
}

/*********************************/
/* used by integrate1D and sum1D */
/*********************************/
int calloc_int_values (input_struct input, output_struct* output) {
  int status = 0;

  /* callocate wavelength integrated fluxes and radiances */
  output->rfldir_int = calloc (output->atm.nzout, sizeof (double));
  CHKPOINTER (output->rfldir_int);
  output->rfldn_int = calloc (output->atm.nzout, sizeof (double));
  CHKPOINTER (output->rfldn_int);
  output->flup_int = calloc (output->atm.nzout, sizeof (double));
  CHKPOINTER (output->flup_int);
  output->uavg_int = calloc (output->atm.nzout, sizeof (double));
  CHKPOINTER (output->uavg_int);
  output->uavgso_int = calloc (output->atm.nzout, sizeof (double));
  CHKPOINTER (output->uavgso_int);
  output->uavgdn_int = calloc (output->atm.nzout, sizeof (double));
  CHKPOINTER (output->uavgdn_int);
  output->uavgup_int = calloc (output->atm.nzout, sizeof (double));
  CHKPOINTER (output->uavgup_int);
  if (input.heating != HEAT_NONE) {
    output->heat_int = calloc (output->atm.nzout, sizeof (double));
    CHKPOINTER (output->heat_int);
    output->emis_int = calloc (output->atm.nzout, sizeof (double));
    CHKPOINTER (output->emis_int);
    output->w_zout_int = calloc (output->atm.nzout, sizeof (double));
    CHKPOINTER (output->w_zout_int);
  }

  output->u0u_int = calloc (output->atm.nzout, sizeof (double*));
  CHKPOINTER (output->u0u_int);
  for (size_t lev = 0; lev < output->atm.nzout; lev++) {
    output->u0u_int[lev] = calloc (input.rte.numu, sizeof (double));
    CHKPOINTER (output->u0u_int[lev]);
  }

  output->uu_int = calloc (output->atm.nzout, sizeof (double**));
  CHKPOINTER (output->uu_int);
  for (size_t lev = 0; lev < output->atm.nzout; lev++) {
    output->uu_int[lev] = calloc (input.rte.nphi, sizeof (double*));
    CHKPOINTER (output->uu_int[lev]);
    for (size_t j = 0; j < input.rte.nphi; j++) {
      output->uu_int[lev][j] = calloc (input.rte.numu, sizeof (double));
      CHKPOINTER (output->uu_int[lev][j]);
    }
  }

  output->albmed_int = calloc (input.rte.numu, sizeof (double));
  CHKPOINTER (output->albmed_int);
  output->trnmed_int = calloc (input.rte.numu, sizeof (double));
  CHKPOINTER (output->trnmed_int);

  return status;
}

/***************************************************************************/
/* abuse first entry of the wavelength index to store the integrated value */
/* and to convert double (nessesary for heating rates) to floats           */
/* entry (n+1) or entry (-1) would be nicer                                */
/***************************************************************************/

int double2float_integrated_values (input_struct input, output_struct* output) {
  int status = 0;

  for (size_t iu = 0; iu < input.rte.numu; iu++) {

    output->albmed[iu][0] = (float) output->albmed_int[iu];
    output->trnmed[iu][0] = (float) output->trnmed_int[iu];
  }

  for (size_t lev = 0; lev < output->atm.nzout; lev++) {

    /* irradiance / actinic flux */
    output->rfldir[lev][0] = (float) output->rfldir_int[lev];
    output->rfldn[lev][0]  = (float) output->rfldn_int[lev];
    output->flup[lev][0]   = (float) output->flup_int[lev];
    output->uavg[lev][0]   = (float) output->uavg_int[lev];
    output->uavgso[lev][0] = (float) output->uavgso_int[lev];
    output->uavgdn[lev][0] = (float) output->uavgdn_int[lev];
    output->uavgup[lev][0] = (float) output->uavgup_int[lev];
    if (input.heating != HEAT_NONE) {
      output->heat[lev][0]   = (float) output->heat_int[lev];
      output->emis[lev][0]   = (float) output->emis_int[lev];
      output->w_zout[lev][0] = (float) output->w_zout_int[lev];
    }

    /* radiances   */
    for (size_t iu = 0; iu < input.rte.numu; iu++) {
      output->u0u[lev][iu][0] = (float) output->u0u_int[lev][iu];

      for (size_t j = 0; j < input.rte.nphi; j++)
        output->uu[lev][j][iu][0] = (float) output->uu_int[lev][j][iu];
    }
  }

  return status;
}

/*********************************/
/* used by integrate1D and sum1D */
/*********************************/

int scaling_integrated_values (input_struct input, output_struct* output) {
  double ffactor = 0, rfactor = 0, hfactor = 1;
  int    status = 0;

  /* if transmittance, then divide by the extraterrestrial flux */
  if (input.source != SRC_THERMAL) {
    switch (input.calibration) {
    case OUTCAL_ABSOLUTE:
      ffactor = 1.0; /* irradiance multiplicator */
      rfactor = 1.0; /* radiance multiplicator   */
      break;

    case OUTCAL_TRANSMITTANCE:
    case OUTCAL_BRIGHTNESS:
      ffactor = 1.0 / output->wl.fbeam[0]; /* irradiance multiplicator */
      rfactor = 1.0 / output->wl.fbeam[0]; /* radiance multiplicator   */
      break;

    case OUTCAL_REFLECTIVITY:
      ffactor = 1.0 / output->incident;
      rfactor = M_PI / output->incident;
      break;

    default:
      fprintf (stderr, "Error, unknown output calibration %d\n", input.calibration);
      CHKERR (-1);
    }

    /**************************************************************/
    /* now scale irradiances with ffactor, radiances with rfactor */
    /**************************************************************/

    if (output->wl.nlambda_h > 0) {
      /* in iv==0 the integrated values are stored */
      status = scale_output (input,
                             &(output->rfldir),
                             &(output->rfldn),
                             &(output->flup),
                             &(output->albmed),
                             &(output->trnmed),
                             &(output->uavgso),
                             &(output->uavgdn),
                             &(output->uavgup),
                             &(output->uavg),
                             &(output->u0u),
                             &(output->uu),
                             &(output->heat),
                             &(output->emis),
                             &(output->w_zout),
                             &(output->rfldir3d),
                             &(output->rfldn3d),
                             &(output->flup3d),
                             &(output->fl3d_is),
                             &(output->uavgso3d),
                             &(output->uavgdn3d),
                             &(output->uavgup3d),
                             &(output->radiance3d),
                             &(output->jacobian),
                             &(output->absback3d),
                             &(output->rfldir3d_var),
                             &(output->rfldn3d_var),
                             &(output->flup3d_var),
                             &(output->uavgso3d_var),
                             &(output->uavgdn3d_var),
                             &(output->uavgup3d_var),
                             &(output->radiance3d_var),
                             &(output->abs3d_var),
                             &(output->absback3d_var),
                             output->atm.nzout,
                             output->atm.Nxcld,
                             output->atm.Nycld,
                             output->atm.Nzcld,
                             output->mc.alis.Nc,
                             output->atm.nlyr,
                             output->atm.threed,
                             output->islower,
                             output->isupper,
                             output->jslower,
                             output->jsupper,
                             output->isstep,
                             output->jsstep,
                             &(output->abs3d),
                             output->triangle_results_o,
                             ffactor,
                             rfactor,
                             hfactor,
                             0);
      CHKERR (status);
    }
  } else {
    if (input.calibration == OUTCAL_BRIGHTNESS) {

      if (!input.quiet)
        fprintf (stderr, " ... converting radiances to brightness temperatures\n");

      /* convert irradiances / radiances to brightness temperatures */
      /* the first element (iv==0) stores the spectrally integrated quantities; */
      status = output2bt (input, output, 0, 0);
      CHKERR (status);
    }
  }

  /* abuse first entry to store the integrated value, this is not so nice programmed */
  if (output->wl.nlambda_h > 0)
    output->wl.nlambda_h = 1;

  return 0;
}

/***********************************************************************************/
/* Convert Raman spectrum to user spectum                                          */
/***********************************************************************************/

static int raman_spec2spec (input_struct input, output_struct* output) {
  int   lev = 0, iu = 0, j = 0, status = 0, iv = 0;
  int   start_id = 0, end_id = 0;
  float test_lower = 0, test_upper = 0;

  test_lower = output->wl.lambda_h[0] + output->wl.delta_wvl_raman_lower + output->wl.delta_wvl_raman_extra;
  test_upper = output->wl.lambda_h[output->wl.nlambda_h - 1] - output->wl.delta_wvl_raman_upper - output->wl.delta_wvl_raman_extra;
  for (iv = 0; iv < output->wl.nlambda_h; iv++) {
    if (output->wl.lambda_h[iv] < test_lower)
      start_id = iv + 1;
    if (output->wl.lambda_h[output->wl.nlambda_h - 1 - iv] > test_upper)
      end_id = output->wl.nlambda_h - iv - 2;
  }

  /* The old nlambda_h is always larger than what we are going to output. Furthermore, */
  /* start_id is always > 0 and end_id always < nlambda_h, so we just shift and      */
  /* arys decrease nlambda_h                                                           */

  output->wl.nlambda_h = end_id - start_id + 1;

  for (iv = 0; iv < output->wl.nlambda_h; iv++) {
    output->wl.lambda_h[iv] = output->wl.lambda_h[iv + start_id];
    for (lev = 0; lev < output->atm.nzout; lev++) {
      output->rfldir[lev][iv] = output->rfldir[lev][iv + start_id];
      output->rfldn[lev][iv]  = output->rfldn[lev][iv + start_id];
      output->flup[lev][iv]   = output->flup[lev][iv + start_id];
      output->uavgso[lev][iv] = output->uavgso[lev][iv + start_id];
      output->uavgdn[lev][iv] = output->uavgdn[lev][iv + start_id];
      output->uavgup[lev][iv] = output->uavgup[lev][iv + start_id];
      for (iu = 0; iu < input.rte.numu; iu++) {
        output->u0u[lev][iu][iv] = output->u0u[lev][iu][iv + start_id];
        for (j = 0; j < input.rte.nphi; j++)
          output->uu[lev][j][iu][iv] = output->uu[lev][j][iu][iv + start_id];
      }
    }
  }

  return status;
}

/***********************************************************************************/
/* Convert spectra to RGB and write data to the first three array elements         */
/***********************************************************************************/

static int spec2rgb (input_struct input, output_struct* output) {
  int status;

  if (output->wl.nlambda_h < 3) {
    CHKERROUT (-1, "Fatal error, need at least 3 wavelengths to store RGB!");
  }

  for (size_t lev = 0; lev < output->atm.nzout; lev++) {
    status = spectrum_to_rgb_overwrite (output->rfldir[lev], output->wl.nlambda_h);
    CHKERR (status);
    status = spectrum_to_rgb_overwrite (output->rfldn[lev], output->wl.nlambda_h);
    CHKERR (status);
    status = spectrum_to_rgb_overwrite (output->flup[lev], output->wl.nlambda_h);
    CHKERR (status);
    status = spectrum_to_rgb_overwrite (output->uavgso[lev], output->wl.nlambda_h);
    CHKERR (status);
    status = spectrum_to_rgb_overwrite (output->uavgdn[lev], output->wl.nlambda_h);
    CHKERR (status);
    status = spectrum_to_rgb_overwrite (output->uavgup[lev], output->wl.nlambda_h);
    CHKERR (status);

    for (size_t iu = 0; iu < input.rte.numu; iu++) {
      status = spectrum_to_rgb_overwrite (output->u0u[lev][iu], output->wl.nlambda_h);
      CHKERR (status);

      for (size_t j = 0; j < input.rte.nphi; j++) {
        status = spectrum_to_rgb_overwrite (output->uu[lev][j][iu], output->wl.nlambda_h);
        CHKERR (status);
      }
    }
  }

  /* normalized and gamma corrected */
  if (input.processing == PROCESS_RGB) {
    status = image_wide_rgb_correction1d (output->rfldir, output->atm.nzout);
    CHKERR (status);
    status = image_wide_rgb_correction1d (output->rfldn, output->atm.nzout);
    CHKERR (status);
    status = image_wide_rgb_correction1d (output->flup, output->atm.nzout);
    CHKERR (status);
    status = image_wide_rgb_correction1d (output->uavgso, output->atm.nzout);
    CHKERR (status);
    status = image_wide_rgb_correction1d (output->uavgdn, output->atm.nzout);
    CHKERR (status);
    status = image_wide_rgb_correction1d (output->uavgup, output->atm.nzout);
    CHKERR (status);
    status = image_wide_rgb_correction2d (output->u0u, output->atm.nzout, input.rte.numu - 1);
    CHKERR (status);
    status = image_wide_rgb_correction3d (output->uu, output->atm.nzout, 0, input.rte.nphi - 1, 1, 0, input.rte.numu - 1, 1);
    CHKERR (status);
  }

  if (input.verbose)
    fprintf (stderr, "LEAVING spec2rgb()\n");

  if (output->wl.nlambda_h > 0)
    output->wl.nlambda_h = 3;

  output->wl.lambda_h[0] = 600;
  output->wl.lambda_h[1] = 500;
  output->wl.lambda_h[2] = 400;

  return 0;
}

/***********************************************************************************/
/* Convert 3D spectra to RGB and write data to the first three array elements      */
/***********************************************************************************/

static int spec2rgb3D (input_struct input, output_struct* output) {
  int is = 0, js = 0, ks = 0, status = 0;

  if (output->wl.nlambda_h < 3) {
    CHKERROUT (-1, "Fatal error, need at least 3 wavelengths to store RGB!");
  }

  for (ks = 0; ks < output->atm.nzout; ks++) {
    for (is = output->islower; is <= output->isupper; is += output->isstep) {
      for (js = output->jslower; js <= output->jsupper; js += output->jsstep) {

        if (input.verbose)
          fprintf (stderr, "R3D = %f\n", output->radiance3d[ks][is][js][0][0][0]);
        status = spectrum_to_rgb_overwrite (output->rfldir3d[ks][is][js], output->wl.nlambda_h);
        CHKERR (status);
        status = spectrum_to_rgb_overwrite (output->rfldn3d[ks][is][js], output->wl.nlambda_h);
        CHKERR (status);
        status = spectrum_to_rgb_overwrite (output->flup3d[ks][is][js], output->wl.nlambda_h);
        CHKERR (status);
        status = spectrum_to_rgb_overwrite (output->uavgso3d[ks][is][js], output->wl.nlambda_h);
        CHKERR (status);
        status = spectrum_to_rgb_overwrite (output->uavgdn3d[ks][is][js], output->wl.nlambda_h);
        CHKERR (status);
        status = spectrum_to_rgb_overwrite (output->uavgup3d[ks][is][js], output->wl.nlambda_h);
        CHKERR (status);
        /* ignore polarization for rgb output */
        status = spectrum_to_rgb_overwrite (output->radiance3d[ks][is][js][0][0], output->wl.nlambda_h);
        CHKERR (status);
        /* ignore variances for rgb output */

        if (input.verbose) {
          fprintf (stderr, "SPEC2RGB %d %d %d  %d\n", ks, is, js, output->wl.nlambda_h);
          fprintf (stderr, "R3D = %f\n", output->radiance3d[ks][is][js][0][0][0]);
        }
      }
    }
  }
  /* normalized and gamma corrected */
  if (input.processing == PROCESS_RGB) {
    status = image_wide_rgb_correction3d (output->rfldir3d,
                                          output->atm.nzout,
                                          output->islower,
                                          output->isupper,
                                          output->isstep,
                                          output->jslower,
                                          output->jsupper,
                                          output->jsstep);
    CHKERR (status);
    status = image_wide_rgb_correction3d (output->rfldn3d,
                                          output->atm.nzout,
                                          output->islower,
                                          output->isupper,
                                          output->isstep,
                                          output->jslower,
                                          output->jsupper,
                                          output->jsstep);
    CHKERR (status);
    status = image_wide_rgb_correction3d (output->flup3d,
                                          output->atm.nzout,
                                          output->islower,
                                          output->isupper,
                                          output->isstep,
                                          output->jslower,
                                          output->jsupper,
                                          output->jsstep);
    CHKERR (status);
    status = image_wide_rgb_correction3d (output->uavgso3d,
                                          output->atm.nzout,
                                          output->islower,
                                          output->isupper,
                                          output->isstep,
                                          output->jslower,
                                          output->jsupper,
                                          output->jsstep);
    CHKERR (status);
    status = image_wide_rgb_correction3d (output->uavgdn3d,
                                          output->atm.nzout,
                                          output->islower,
                                          output->isupper,
                                          output->isstep,
                                          output->jslower,
                                          output->jsupper,
                                          output->jsstep);
    CHKERR (status);
    status = image_wide_rgb_correction3d (output->uavgup3d,
                                          output->atm.nzout,
                                          output->islower,
                                          output->isupper,
                                          output->isstep,
                                          output->jslower,
                                          output->jsupper,
                                          output->jsstep);
    CHKERR (status);
    status = image_wide_rgb_correction3d (output->uavgup3d,
                                          output->atm.nzout,
                                          output->islower,
                                          output->isupper,
                                          output->isstep,
                                          output->jslower,
                                          output->jsupper,
                                          output->jsstep);
    CHKERR (status);
    status = image_wide_rgb_correction3drad (output->radiance3d,
                                             output->atm.nzout,
                                             output->islower,
                                             output->isupper,
                                             output->isstep,
                                             output->jslower,
                                             output->jsupper,
                                             output->jsstep);
    CHKERR (status);
  }

  /* don't adjust lambda_h and nlambda_h here - we need the */
  /* original numbers for processing1D()!                   */
  /*
  if (output->wl.nlambda_h>0)
    output->wl.nlambda_h = 3;

  output->wl.lambda_h[0] = 600;
  output->wl.lambda_h[1] = 500;
  output->wl.lambda_h[2] = 400;
  */

  if (input.verbose)
    fprintf (stderr, "LEAVING spec2rgb3D()\n");

  if (output->triangle_results_o)
    CHKERROUT (-1, "Integrate not yet implemented for triangles");

  return 0;
}

/***********************************************************************************/
/* Scale uvspec output fields with a given factor; ffactor is the scaling factor   */
/* for irradiances/actinic fluxes and rfactor is the scaling factor for radiances  */
/***********************************************************************************/

int scale_output (input_struct                 input,
                  float***                     p_rfldir,
                  float***                     p_rfldn,
                  float***                     p_flup,
                  float***                     p_albmed,
                  float***                     p_trnmed,
                  float***                     p_uavgso,
                  float***                     p_uavgdn,
                  float***                     p_uavgup,
                  float***                     p_uavg,
                  float****                    p_u0u,
                  float*****                   p_uu,
                  float***                     p_heat,
                  float***                     p_emis,
                  float***                     p_w_zout,
                  float*****                   p_rfldir3d,
                  float*****                   p_rfldn3d,
                  float*****                   p_flup3d,
                  float******                  p_fl3d_is,
                  float*****                   p_uavgso3d,
                  float*****                   p_uavgdn3d,
                  float*****                   p_uavgup3d,
                  float*******                 p_radiance3d,
                  float********                p_jacobian,
                  float*****                   p_absback3d,
                  float*****                   p_rfldir3d_var,
                  float*****                   p_rfldn3d_var,
                  float*****                   p_flup3d_var,
                  float*****                   p_uavgso3d_var,
                  float*****                   p_uavgdn3d_var,
                  float*****                   p_uavgup3d_var,
                  float******                  p_radiance3d_var,
                  float*****                   p_abs3d_var,
                  float*****                   p_absback3d_var,
                  int                          nzout,
                  int                          Nx,
                  int                          Ny,
                  int                          Nz,
                  int                          Nc,
                  int                          Nlyr,
                  int*                         threed,
                  int                          islower,
                  int                          isupper,
                  int                          jslower,
                  int                          jsupper,
                  int                          isstep,
                  int                          jsstep,
                  float*****                   p_abs3d,
                  t_triangle_radiation_field** triangle_result,
                  double                       ffactor,
                  double                       rfactor,
                  double                       hfactor,
                  int                          iv)

/* changed output_struc to a bunch of pointers, as this function is also used in solve_rte() with different arguements, UH, 2006-02 */
{
  double ffactor2 = 0, rfactor2 = 0, hfactor2 = 0; /* 27.02.2013 **CK **BM: add hfactor2 for thermal backward heating rate std */

  /* Arve 20160121: spherical albedo and transmittance should not be scaled
     if disort_spherical_albedo is set. */
  if (!input.rte.ibcnd) {
    for (int iu = 0; iu < input.rte.numu; iu++) {
      (*p_albmed)[iu][iv] *= ffactor;
      (*p_trnmed)[iu][iv] *= ffactor;
    }
  }

  for (int lev = 0; lev < nzout; lev++) {

    /* irradiance / actinic flux / heating rate */
    (*p_rfldir)[lev][iv] *= ffactor;
    (*p_rfldn)[lev][iv] *= ffactor;
    (*p_flup)[lev][iv] *= ffactor;
    (*p_uavg)[lev][iv] *= ffactor;
    (*p_uavgso)[lev][iv] *= ffactor;
    (*p_uavgdn)[lev][iv] *= ffactor;
    (*p_uavgup)[lev][iv] *= ffactor;
    if (input.heating != HEAT_NONE) {
      (*p_heat)[lev][iv] *= hfactor;
      (*p_emis)[lev][iv] *= hfactor;
      (*p_w_zout)[lev][iv] *= hfactor;
    }

    /* radiances   */
    for (int iu = 0; iu < input.rte.numu; iu++) {
      (*p_u0u)[lev][iu][iv] *= rfactor;

      for (int j = 0; j < input.rte.nphi; j++)
        (*p_uu)[lev][j][iu][iv] *= rfactor;
    }

    /* 3D irradiances and radiances */
    if (input.output3D) {
      for (int is = islower; is <= isupper; is += isstep) {
        for (int js = jslower; js <= jsupper; js += jsstep) {
          (*p_rfldir3d)[lev][is][js][iv] *= ffactor;
          (*p_rfldn3d)[lev][is][js][iv] *= ffactor;
          (*p_flup3d)[lev][is][js][iv] *= ffactor;
          (*p_uavgso3d)[lev][is][js][iv] *= ffactor;
          (*p_uavgdn3d)[lev][is][js][iv] *= ffactor;
          (*p_uavgup3d)[lev][is][js][iv] *= ffactor;

          if (input.rte.mc.concentration_is || input.rte.mc.spectral_is) {
            for (int ic = 0; ic < Nc; ic++)
              (*p_fl3d_is)[lev][is][js][ic][iv] *= ffactor;
          }

          if (*p_radiance3d) {
            for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
              for (int ic = 0; ic < Nc; ic++) {
                (*p_radiance3d)[lev][is][js][ip][ic][iv] *= rfactor;
              }
            }
          }

          if (input.rte.mc.jacobian[DIM_1D]) {
            for (int isp = 0; isp < input.n_caoth + 2; isp++) {
              for (int i = 0; i < 2; i++) {
                for (int lc = 0; lc < Nlyr; lc++) {
                  if (input.rte.mc.backward.output == MCBACKWARD_EDN || input.rte.mc.backward.output == MCBACKWARD_EDNPV ||
                      input.rte.mc.backward.output == MCBACKWARD_EUP) {
                    (*p_jacobian)[lev][is][js][isp][i][lc][iv] *= ffactor;
                  } else {
                    (*p_jacobian)[lev][is][js][isp][i][lc][iv] *= rfactor;
                  }
                }
              }
            }
          }

          if (input.rte.mc.backward.absorption && input.ipa3d != 1) { /*ulrike: added && input.ipa3d!=1*/
            (*p_absback3d)[lev][is][js][iv] *= ffactor;
          }

          /* ulrike 3.5.2010: if we have ipa3d
	       then we use absback3d to save the heating rate; thus, p_absback3d should not be multiplied
	       with the ffactor, but with the hfactor, which is the heating rate factor*/
          if (input.rte.mc.backward.absorption && input.ipa3d) {
            (*p_absback3d)[lev][is][js][iv] *= hfactor;
          }

          /* variances */
          if (input.rte.mc.std) {

            ffactor2 = ffactor * ffactor;
            rfactor2 = rfactor * rfactor;
            hfactor2 = hfactor * hfactor; /* 27.02.2013 **CK **BM: add hfactor2 for thermal backward heating rate std */

            (*p_rfldir3d_var)[lev][is][js][iv] *= ffactor2;
            (*p_rfldn3d_var)[lev][is][js][iv] *= ffactor2;
            (*p_flup3d_var)[lev][is][js][iv] *= ffactor2;
            (*p_uavgso3d_var)[lev][is][js][iv] *= ffactor2;
            (*p_uavgdn3d_var)[lev][is][js][iv] *= ffactor2;
            (*p_uavgup3d_var)[lev][is][js][iv] *= ffactor2;

            if (*p_radiance3d_var) {
              for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
                (*p_radiance3d_var)[lev][is][js][ip][iv] *= rfactor2;
              }
            }

            /* 27.02.2013 **CK **BM: add "&&" for thermal backward heating rate std */
            if (input.rte.mc.backward.absorption && input.ipa3d != 1) {
              (*p_absback3d_var)[lev][is][js][iv] *= ffactor2;
            }

            /* 27.02.2013 **CK **BM: add hfactor2 for thermal backward heating rate std */
            if (input.rte.mc.backward.absorption && input.ipa3d) {
              (*p_absback3d_var)[lev][is][js][iv] *= hfactor2;
            }
          }
        }
      }
    }
  }

  /* 3D absorption fields; here the loop goes over all 3D boxes */
  if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
      input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d > 0) {
    const size_t N3d_layers = vec_count_true (Nz, threed);
    for (size_t k = 0; k < N3d_layers; k++) {
      for (int is = 0; is < Nx; is++) {
        for (int js = 0; js < Ny; js++) { /* **CK added bracket */
          (*p_abs3d)[k][is][js][iv] *= ffactor;

          if (input.rte.mc.std) {
            /** CK added for forward mc_std
                 *  ???????????????
                 *  CK shouldn't we also switch between ffactor und hfactor for ipa3d/not ipa3d as in backward mode?
                 **/
            (*p_abs3d_var)[k][is][js][iv] *= ffactor2;
          }
        }
      }
    }
  }

  if (triangle_result) {
    const int status = scale_output_triangle_surface (ffactor, triangle_result[iv]);
    CHKERR (status);
  }

  return 0;
}

int scale_output_triangle_surface (const double ffactor, t_triangle_radiation_field* result) {
  for (size_t i = 0; i < result->N_triangles; ++i) {
    result->edir[i] *= ffactor;
    result->edn[i] *= ffactor;
    result->eup[i] *= ffactor;
  }
  return 0;
}

/***********************************************************************************/
/* Convert radiance to brightness temperature for a given filter function.         */
/***********************************************************************************/

float radiance2bt (float rad, float* wvnmlo, float* wvnmhi, float* filter, int n, int processing, int ivi) {
  int     iv = 0, it = 0;
  double  t1 = 0, t2 = 0, plkrad1 = 0, plkrad2 = 0, r = 0, wvlmlo = 0, wvlmhi = 0;
  double *xint = NULL, *yint1 = NULL, *yint2 = NULL;

  /* careful: if accuracy is set too small, it may not be reachable because */
  /* the difference between two neighbouring floats may be larger!          */
  /* Should better tie this to the actual numerical precision               */
  double accur = 5e-5;
  double dt    = 0.001;

  /* special treatment of NaN */
  if (rad != rad)
    return rad;

  t1 = 273;
  t2 = t1 + dt;

  if (rad <= 0)
    return 0;

  plkrad1 = 0;
  plkrad2 = 0;

  if (processing == PROCESS_INT) {

    xint  = calloc (n, sizeof (double));
    yint1 = calloc (n, sizeof (double));
    yint2 = calloc (n, sizeof (double));

    for (iv = 0; iv < n; iv++) {

      wvlmlo = wvnmlo[iv];
      wvlmhi = wvnmhi[iv];

      // wavenumber is multiplied by -1 since the integrate()-function assumes increasing x-values
      xint[iv] = -(wvnmlo[iv] + wvnmhi[iv]) / 2.0;

      // the difference between wvlmlo and wvlmhi is assumed to be 1 cm^-1, otherwise the integrated value is probably wrong
      r         = c_planck_func1 (wvlmlo, wvlmhi, t1);
      yint1[iv] = filter[iv] * r;
      r         = c_planck_func1 (wvlmlo, wvlmhi, t2);
      yint2[iv] = filter[iv] * r;
    }

    plkrad1 = integrate (xint, yint1, n);
    plkrad2 = integrate (xint, yint2, n);

  } else if (processing == PROCESS_SUM) {
    for (iv = 0; iv < n; iv++) {

      wvlmlo = wvnmlo[iv];
      wvlmhi = wvnmhi[iv];
      r      = c_planck_func1 (wvlmlo, wvlmhi, t1);
      plkrad1 += (double) filter[iv] * r;
      r = c_planck_func1 (wvlmlo, wvlmhi, t2);
      plkrad2 += (double) filter[iv] * r;
    }
  } else {

    wvlmlo  = wvnmlo[ivi];
    wvlmhi  = wvnmhi[ivi];
    plkrad1 = c_planck_func1 (wvlmlo, wvlmhi, t1);
    plkrad2 = c_planck_func1 (wvlmlo, wvlmhi, t2);
  }

  it = 0;
  while (fabs ((plkrad1 - rad) / rad) > accur) {

    t1 = t1 + (rad - plkrad1) / (plkrad2 - plkrad1) * dt;
    t2 = t1 + dt;

    plkrad1 = 0;
    plkrad2 = 0;
    if (processing == PROCESS_INT) {

      for (iv = 0; iv < n; iv++) {

        wvlmlo = wvnmlo[iv];
        wvlmhi = wvnmhi[iv];

        r         = c_planck_func1 (wvlmlo, wvlmhi, t1);
        yint1[iv] = filter[iv] * r;
        r         = c_planck_func1 (wvlmlo, wvlmhi, t2);
        yint2[iv] = filter[iv] * r;
      }

      plkrad1 = integrate (xint, yint1, n);
      plkrad2 = integrate (xint, yint2, n);

    } else if (processing == PROCESS_SUM) {
      for (iv = 0; iv < n; iv++) {

        wvlmlo = wvnmlo[iv];
        wvlmhi = wvnmhi[iv];
        r      = c_planck_func1 (wvlmlo, wvlmhi, t1);
        plkrad1 += (double) filter[iv] * r;
        r = c_planck_func1 (wvlmlo, wvlmhi, t2);
        plkrad2 += (double) filter[iv] * r;
      }
    } else {

      wvlmlo  = wvnmlo[ivi];
      wvlmhi  = wvnmhi[ivi];
      plkrad1 = c_planck_func1 (wvlmlo, wvlmhi, t1);
      plkrad2 = c_planck_func1 (wvlmlo, wvlmhi, t2);
    }
    if (it > 999) {
      fprintf (stderr, "While loop in function %s, file %s, did not converge.\n", __func__, __FILE__);
      fprintf (stderr, "Continuing anyway, but be careful with results.\n");
      fprintf (stderr, "The relevant variables in %s have the following values:\n", __func__);
      fprintf (stderr, "t1      = %12.6f\n", t1);
      fprintf (stderr, "t2      = %12.6f\n", t2);
      fprintf (stderr, "plkrad1 = %12.6e\n", plkrad1);
      fprintf (stderr, "plkrad2 = %12.6e\n", plkrad2);
      fprintf (stderr, "rad     = %12.6e\n", rad);
      fprintf (stderr, "accur   = %12.6e\n", accur);
      fprintf (stderr, "fabs((plkrad1 - rad)/rad) = %12.6e\n", fabs ((plkrad1 - rad) / rad));
      break;
    }
    it++;
  }

  if (processing == PROCESS_INT) {

    free (xint);
    free (yint1);
    free (yint2);
  }

  return t1;
}

/***********************************************************************************/
/* Convert uvspec output to brightness temperatures; currently it is assumed that  */
/* the first element (iv==0) stores the spectrally integrated quantities;          */
/* only these will be considered.                                                  */
/***********************************************************************************/

static int output2bt (input_struct input, output_struct* output, int iv, int is_3d) {
  int lev = 0, iu = 0, j = 0, is = 0, ih = 0, js = 0, ic = 0, ip = 0;

  float *wvnmlo, *wvnmhi, *weight;
  int    nlambda;

  int processing;

  processing = input.processing;

  /* Use the representative wavelengths (_r-grid) only when no integration or summation was requested by the user */
  if (output->wl.use_reptran && input.processing == PROCESS_NONE) {

    processing = PROCESS_SUM;

    nlambda = output->wl.nlambda_in_reptran_band[output->wl.reptran_band_t[output->wl.map_e2h[iv]]];

    wvnmlo = calloc (nlambda, sizeof (float));
    wvnmhi = calloc (nlambda, sizeof (float));
    weight = calloc (nlambda, sizeof (float));

    for (ih = 0; ih < nlambda; ih++) {

      wvnmlo[ih] = output->wl.lambda_r[output->wl.reptran_band[output->wl.reptran_band_t[output->wl.map_e2h[iv]]][ih]];
      wvnmlo[ih] = 1.0E7 / wvnmlo[ih] - input.bandwidth / 2.0;
      wvnmhi[ih] = wvnmlo[ih] + input.bandwidth;
      weight[ih] = output->wl.weight_reptran_band[output->wl.reptran_band_t[output->wl.map_e2h[iv]]][ih];
    }

  } else {

    wvnmlo  = output->wl.wvnmlo_h;
    wvnmhi  = output->wl.wvnmhi_h;
    weight  = output->wl.filter;
    nlambda = output->wl.nlambda_h;
  }

  if (is_3d) { /* 3D processing */

    for (lev = 0; lev < output->atm.nzout; lev++) {

      /* 3D fields */
      if (input.output3D)
        for (is = output->islower; is <= output->isupper; is += output->isstep)
          for (js = output->jslower; js <= output->jsupper; js += output->jsstep) {

            output->rfldir3d[lev][is][js][iv] =
              radiance2bt (output->rfldir3d[lev][is][js][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

            output->rfldn3d[lev][is][js][iv] =
              radiance2bt (output->rfldn3d[lev][is][js][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

            output->flup3d[lev][is][js][iv] =
              radiance2bt (output->flup3d[lev][is][js][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

            output->uavgso3d[lev][is][js][iv] =
              radiance2bt (output->uavgso3d[lev][is][js][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

            output->uavgdn3d[lev][is][js][iv] =
              radiance2bt (output->uavgdn3d[lev][is][js][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

            output->uavgup3d[lev][is][js][iv] =
              radiance2bt (output->uavgup3d[lev][is][js][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

            if (input.rte.mc.concentration_is || input.rte.mc.spectral_is)
              for (ic = 0; ic < output->mc.alis.Nc; ic++)
                output->fl3d_is[lev][is][js][ic][iv] =
                  radiance2bt (output->fl3d_is[lev][is][js][ic][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

            for (ip = 0; ip < input.rte.mc.nstokes; ip++)
              for (ic = 0; ic < output->mc.alis.Nc; ic++)
                output->radiance3d[lev][is][js][ip][ic][iv] =
                  radiance2bt (output->radiance3d[lev][is][js][ip][ic][iv], wvnmlo, wvnmhi, weight, nlambda, processing, iv);

            /* it would certainly be nonsense to convert standard deviations to */
            /* brightness temperature; therefore set to NaN                     */

            /* variances */
            if (input.rte.mc.std) {

              output->rfldir3d_var[lev][is][js][iv] = 0.0 / 0.0;
              output->rfldn3d_var[lev][is][js][iv]  = 0.0 / 0.0;
              output->flup3d_var[lev][is][js][iv]   = 0.0 / 0.0;
              output->uavgso3d_var[lev][is][js][iv] = 0.0 / 0.0;
              output->uavgdn3d_var[lev][is][js][iv] = 0.0 / 0.0;
              output->uavgup3d_var[lev][is][js][iv] = 0.0 / 0.0;

              for (ip = 0; ip < input.rte.mc.nstokes; ip++)
                output->radiance3d_var[lev][is][js][ip][iv] = 0.0 / 0.0;

              if (input.rte.mc.backward.absorption)
                output->absback3d_var[lev][is][js][iv] = 0.0 / 0.0;
            }
          }
    }

    /* conversion of 3D absorption to BT is probably useless, hence we don't do it */

  } else { /* 1D processing */

    for (lev = 0; lev < output->atm.nzout; lev++) {

      /* irradiance / actinic flux */
      output->rfldir[lev][iv] = radiance2bt (output->rfldir[lev][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

      output->rfldn[lev][iv] = radiance2bt (output->rfldn[lev][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

      output->flup[lev][iv] = radiance2bt (output->flup[lev][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

      output->uavg[lev][iv] = radiance2bt (output->uavg[lev][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

      output->uavgso[lev][iv] = radiance2bt (output->uavgso[lev][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

      output->uavgdn[lev][iv] = radiance2bt (output->uavgdn[lev][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

      output->uavgup[lev][iv] = radiance2bt (output->uavgup[lev][iv] / M_PI, wvnmlo, wvnmhi, weight, nlambda, processing, iv);

      /* radiances */
      for (iu = 0; iu < input.rte.numu; iu++) {
        output->u0u[lev][iu][iv] = radiance2bt (output->u0u[lev][iu][iv], wvnmlo, wvnmhi, weight, nlambda, processing, iv);

        for (j = 0; j < input.rte.nphi; j++)
          output->uu[lev][j][iu][iv] = radiance2bt (output->uu[lev][j][iu][iv], wvnmlo, wvnmhi, weight, nlambda, processing, iv);
      }
    }
  }

  if (output->wl.use_reptran && input.processing == PROCESS_NONE) {
    free (wvnmlo);
    free (wvnmhi);
    free (weight);
  }

  return 0;
}

static int read_photon_file (char* filename, float* lambda_r, int nlambda_r, float** fraction) {
  float* wvl = NULL;
  int    iv = 0, n = 0, status = 0;

  status = read_2c_file_float (filename, &wvl, fraction, &n);
  if (status != 0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    CHKERR (status);
  }

  /* compare wavelength grids */
  if (nlambda_r != n) {
    fprintf (stderr, "Error, number of wavelengths in %s differing from\n", filename);
    fprintf (stderr, "internal wavelength grid\n");
    return -1;
  }

  for (iv = 0; iv < n; iv++)
    if (lambda_r[iv] != wvl[iv]) {
      fprintf (stderr,
               "Error, wavelengths in %s differ from internal wavelength grid: %f vs. %f\n",
               filename,
               lambda_r[iv],
               wvl[iv]);
      return -1;
    }

  return 0;
}

/**************************************************************************/
/*  average return the average density per layer in the array **y_average */
/*  **y_average will be automatically callocated                          */
/* OUTPUT: dens_avg */
/**************************************************************************/

int average_dens (float* dens, float* dens_air, float* zd, int nlev, int interpol_method, float** dens_avg, int allocate) {
  int    status = 0;
  int    n_layer;
  float* mix; /* mixing ratio*/
  int    lc;

  n_layer = nlev - 1;

  if (allocate)
    if (((*dens_avg) = calloc (n_layer, sizeof (float))) == NULL) {
      fprintf (stderr, "Error allocating memory for (*dens_avg) in average (ancillary.c)\n");
      return -1;
    }

  switch (interpol_method) {
  case INTERP_METHOD_SPLINE:
    status = spline_average (zd, dens, nlev, dens_avg);
    break;
  case INTERP_METHOD_LINEAR:
    for (lc = 0; lc < n_layer; lc++)
      (*dens_avg)[lc] = 0.5 * (dens[lc] + dens[lc + 1]);
    break;
  case INTERP_METHOD_LOG:
    for (lc = 0; lc < n_layer; lc++) {
      (*dens_avg)[lc] = log_average (dens[lc], dens[lc + 1]);
    }
    break;
  case INTERP_METHOD_LINMIX:     /* linear mixing ratio integration                */
  case INTERP_METHOD_LOG_SPLINE: /* no log_spline intergration implemented jet !!! */
    mix = calloc (nlev, sizeof (float));
    if (mix == NULL) {
      fprintf (stderr, "Error allocating memory for mixing_ratio in average (ancillary.c)\n");
      return -1;
    }
    for (lc = 0; lc < nlev; ++lc)
      mix[lc] = dens[lc] / dens_air[lc]; /* calculating to mixing ratio */

    for (lc = 0; lc < n_layer; lc++) {
      (*dens_avg)[lc] = linmix_average (mix[lc], mix[lc + 1], dens_air[lc], dens_air[lc + 1]);
    }
    free (mix);
    break;
  default:
    fprintf (stderr, "Error, unknown interpolation method input.");
    return -1;
  }

  if (status != 0) {
    fprintf (stderr, "Error %d calculating average concentration (average, in ancillary.c)\n", status);
    CHKERR (status);
  }

  return status;
}

/****************************************************************/
/*  log_average returns the logarithmic average density between */
/*  two layers with d[ens]1 und d[ens]2                         */
/*  assuming logarithmic variation with hight                   */
/****************************************************************/

float log_average (float d1, float d2)
/* small function to calculate the average density */
/* assuming a function d=exp(-kx) */
{
  float avg   = 0;
  float test1 = -666., test2 = -666.0;

  test1 = (d1 < d2 ? d1 : d2);
  test2 = fabs (d2 - d1);

  if (test1 <= 0 || test2 <= 0.001 * d1) /* in this case numerically unstable => */
    avg = 0.5 * (d1 + d2);               /* use linear interpolation instead     */
  else
    avg = (d2 - d1) / log (d2 / d1);

  return avg;
}

double dlog_average (double d1, double d2)
/* small function to calculate the average density */
/* assuming a function d=exp(-kx) */
{
  double avg   = 0;
  double test1 = -666., test2 = -666.0;

  test1 = (d1 < d2 ? d1 : d2);
  test2 = fabs (d2 - d1);

  if (test1 <= 0 || test2 <= 0.001 * d1) /* in this case numerically unstable => */
    avg = 0.5 * (d1 + d2);               /* use linear interpolation instead     */
  else
    avg = (d2 - d1) / log (d2 / d1);

  return avg;
}

/*************************************************************************/
/*  linmix_average returns the linear mixing ratio average density n_bar */
/*  between two layers with m[ixing ratio]1 and m2 and                   */
/*  air number density n1 and n2                                         */
/*  assuming linear variation of the mix ratio with height               */
/*  and logarithmic variation of the air number dens with height         */
/*************************************************************************/

float linmix_average (float mmr1, float mmr2, float n1, float n2) {
  float avg   = 0;
  float test1 = -666., test2 = -666.0;
  float log_n = 0;

  test1 = (n1 < n2 ? n1 : n2);
  test2 = fabs (n1 - n2);

  if (test1 <= 0 || test2 <= 0.001 * n1) /* in this case numerically unstable => */
    avg = 0.5 * (mmr1 * n1 + mmr2 * n2); /* use linear interpolation instead     */
  else {
    log_n = log (n2 / n1);
    avg   = 1 / (log_n * log_n) * ((n2 * mmr2 - n1 * mmr1) * log_n - (mmr2 - mmr1) * (n2 - n1));
  }

  return avg;
}

/***********************************************************************/
/*  mass_weighted_average returns the mass weighted average of x       */
/*  between two layers with properties x1,n1 and x2,n2                 */
/*  air number density n1 and n2                                       */
/*  assuming linear variation of the x with height                     */
/*  and logarithmic variation of the air number dens with height       */
/***********************************************************************/

float mass_weighted_average (float x1, float x2, float n1, float n2) {
  float avg   = 0;
  float test1 = -666., test2 = -666.0;

  test1 = (n1 < n2 ? n1 : n2);
  test2 = fabs (n1 - n2);

  if (test1 <= 0 || test2 <= 0.001 * n1)   /* in this case numerically unstable =>    */
    avg = (n1 * x1 + n2 * x2) / (n1 + n2); /* use simple mass weighte average instead */
  else {
    avg = (n2 * x2 - n1 * x1) / (n2 - n1) - (x2 - x1) / log (n2 / n1);
  }

  return avg;
}

/************************************************************/
/* small function to calculate the average density          */
/* assuming cubic spline variation                          */
/* returns an array instead of number as previous functions */
/************************************************************/

int spline_average (float* x, float* y, int n, float** y_average)

{
  double *a0 = NULL, *a1 = NULL, *a2 = NULL, *a3 = NULL;
  int     m;
  int     i;
  double  dx, dx2, dx3;
  int     status;
  int     descend = 0;
  double *x_sort = NULL, *y_sort = NULL;
  double* y_average_d;

  x_sort = calloc (n, sizeof (double));
  y_sort = calloc (n, sizeof (double));

  if (x[1] <= x[0])
    descend = 1;

  /* number of layers = number of level - 1 */
  m = n - 1;

  if (!descend) {
    for (i = 0; i < n; i++) {
      x_sort[i] = (double) x[i];
      y_sort[i] = (double) y[i];
    }
  } else {
    for (i = 0; i < n; i++) {
      x_sort[i] = (double) x[n - 1 - i];
      y_sort[i] = (double) y[n - 1 - i];
    }
  }

  status = spline_coeffc (x_sort, y_sort, n, &a0, &a1, &a2, &a3);
  CHKERR (status);

  y_average_d = calloc (m, sizeof (double));

  for (i = 0; i < m; i++) {
    dx             = (double) (x_sort[i + 1] - x_sort[i]);
    dx2            = (double) (dx * dx);
    dx3            = (double) (dx * dx2);
    y_average_d[i] = ((double) 0.25) * a3[i] * dx3 + ((double) 1. / 3.) * a2[i] * dx2 + ((double) 0.5) * a1[i] * dx + a0[i];
  }

  if (!descend)
    for (i = 0; i < m; i++)
      (*y_average)[i] = (float) y_average_d[i];
  else
    for (i = 0; i < m; i++)
      (*y_average)[i] = (float) y_average_d[m - 1 - i];

  return 0;
}

/**********************************************************/
/* get_grid_index                                         */
/*                                                        */
/* small function that searches the index                 */
/* where "number" is closest to an element of "array"     */
/* periodic takes care of 360 degree periodic boundaries  */
/*                                                        */
/* Ulrich Hamann, Mai 2007                                */
/**********************************************************/

int get_grid_index (float number, double* array, int n_array, int periodic) {

  int   index       = NOT_DEFINED_INTEGER;
  int   i           = 0;
  float min_dist    = NOT_DEFINED_FLOAT;
  float dist        = NOT_DEFINED_FLOAT;
  float grid_dist_1 = NOT_DEFINED_FLOAT;
  float grid_dist_2 = NOT_DEFINED_FLOAT;

  /* if we have only one entry than choose the first and only entry */
  if (n_array == 1) {
    index = 0;
    return index;
  }

  if (periodic) {
    /* initialisation with periodic boundary considering periodicity */
    min_dist = fabs (array[n_array - 1] - 360.0 - number);
    index    = n_array - 1;
    if (fabs (array[n_array - 1] + 360.0 - number) < min_dist)
      min_dist = fabs (array[n_array - 1] + 360.0 - number);
  } else {
    /* initialisation with boundary */
    min_dist = fabs (array[n_array - 1] - number);
    index    = n_array - 1;
  }

  /* find minimum distance between number and grid elements */
  for (i = 0; i < n_array; i++) {
    dist = fabs (array[i] - number);
    if (dist < min_dist) {
      min_dist = dist;
      index    = i;
    }
  }

  /* check 1: index must be defined */
  if (index == NOT_DEFINED_INTEGER) {
    fprintf (stderr, "Error, did NOT find pixel index for pixel %6.0f !!! \n", number);
    CHKERR (-2);
  }

  /* check 2: min_dist should be smaller than distance to next neighbour */
  if (index != 0) {
    grid_dist_1 = fabs (array[index] - array[index - 1]);
  }
  if (index != n_array - 1) {
    grid_dist_2 = fabs (array[index] - array[index + 1]);
    /* choose the largest grid space distance (more relaxed constraint) */
    if (grid_dist_2 > grid_dist_1)
      grid_dist_1 = grid_dist_2;
  }
  if (index == 0 && index == n_array - 1)
    grid_dist_1 = 0;

  if (min_dist > grid_dist_1) {
    fprintf (stderr, "Error while searching grid index, minimum distance (%f) between element %8.2f and grid \n", number, min_dist);
    fprintf (stderr, "       is larger than one grid spacing! (grid spacing = %f) \n", grid_dist_1);
    CHKERR (-1);
  }

  /* correct periodic behaviour at the eastern boundary */
  if (periodic) {
    if (fabs (array[0] + 360.0 - number) < min_dist || fabs (array[0] - 360.0 - number) < min_dist) {
      index = 0;
    }
  }

  return index;
}

/****************************************************************/
/* get_all_netCDF_indices                                       */
/*                                                              */
/* Purpose:                                                     */
/* get size and index for lat, lon, time                        */
/*                                                              */
/*  input:                                                      */
/*  ------                                                      */
/*  ncid       id_number of the open netCDF file                */
/*  lat        latitude                                         */
/*  lon        longitude                                        */
/*  UTC        universal time correlated                        */
/*  time_interpolate    switch for time interpolation           */
/*  verbose    additional verbose output                        */
/*  quiet      no verbose output                                */
/*                                                              */
/*  output:                                                     */
/*  -------                                                     */
/*  ilat       index for latitude                               */
/*  nlat       size of latitude in netCDF file                  */
/*  ilon       index for longitude                              */
/*  nlon       size of longitude in netCDF file                 */
/*  nt         number of time steps (2 if time interpolation)   */
/*  itime1     first  index for time                            */
/*  itime2     second index for time                            */
/*  dt         time step fraction where time                    */
/*             is in between time1 and time2                    */
/*                                                              */
/*  September 2007  by Ulrich Hamann                            */
/****************************************************************/

int get_all_netCDF_indices (char*     filename,
                            float     lat,
                            float     lon,
                            int*      ncid,
                            long*     ilat,
                            size_t*   nlat,
                            long*     ilon,
                            size_t*   nlon,
                            double**  lat_grid,
                            double**  lon_grid,
                            struct tm UTC,
                            int       time_interpolate,
                            int*      nt,
                            int*      itime1,
                            int*      itime2,
                            float*    dt,
                            int       verbose,
                            int       quiet) {

  DEPRECATE (2023);
#if HAVE_NETCDF4

  int status = 0;

  int i = 0, j = 0;

  float lon_tmp = NOT_DEFINED_FLOAT;
  float lat_min = NOT_DEFINED_FLOAT, lat_max = NOT_DEFINED_FLOAT;
  float lon_min = NOT_DEFINED_FLOAT, lon_max = NOT_DEFINED_FLOAT;

  if (lat < -90.0 || 90.0 < lat) {
    fprintf (stderr, "Error, latitude %f outside range, maybe not specified (if == -999)\n", lat);
    CHKERR (-1);
  }

  if (lon < -360.0 || 360.0 < lon) {
    fprintf (stderr, "Error, longitude %f outside range, maybe not specified (if == -999)\n", lon);
    CHKERR (-1);
  }

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, ncid);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d opening netCDF file %s\n", status, filename);
    CHKERR (status);
  }

  /* read latitude array */
  status = alloc_and_read_netCDF_1D_double ((*ncid), "lat", nlat, "lat", lat_grid);
  CHKERR (status);

  /* search correct latitude index */
  (*ilat) = get_grid_index (lat, *lat_grid, (*nlat), FALSE);
  if (*ilat < 0) {
    fprintf (stderr, "Error -1 finding index for lat=%5.2f in %s\n", lat, filename);
    CHKERR (-1);
  }

  /* get range of latitude */
  lat_min = (*lat_grid)[0];
  lat_max = (*lat_grid)[0];
  for (i = 1; i < (*nlat); i++) {
    if ((*lat_grid)[i] > lat_max)
      lat_max = (*lat_grid)[i];
    if ((*lat_grid)[i] < lat_min)
      lat_min = (*lat_grid)[i];
  }

  /* read longitude */
  status = alloc_and_read_netCDF_1D_double ((*ncid), "lon", nlon, "lon", lon_grid);
  CHKERR (status);

  /* get range of longitude */
  lon_min = (*lon_grid)[0];
  lon_max = (*lon_grid)[0];
  for (j = 1; j < (*nlon); j++) {
    if ((*lon_grid)[j] > lon_max)
      lon_max = (*lon_grid)[j];
    if ((*lon_grid)[j] < lon_min)
      lon_min = (*lon_grid)[j];
  }

  lon_tmp = lon;

  /* longitude periodicity */
  /* correct periodicity of latitude so, that longitude is inside the range of the map */
  if (lon_tmp < lon_min)
    lon_tmp += 360.0;
  if (lon_tmp > lon_max)
    lon_tmp -= 360.0;

  /* search correct longitude index */
  (*ilon) = get_grid_index (lon_tmp, *lon_grid, (*nlon), TRUE);
  if (*ilon < 0) {
    fprintf (stderr, "Error -2 finding index for lon=%5.2f in %s\n", lon, filename);
    CHKERR (-2);
  }

  if (verbose)
    fprintf (stderr, "     found %zd x %zd data points\n", (*nlat), (*nlon));

  if (verbose) {
    fprintf (stderr,
             "     map size = [%8.3f (South),%8.3f (North)] x [%8.3f (West),%8.3f (East)]\n",
             lat_min,
             lat_max,
             lon_min,
             lon_max);
  }

  /* get time index */
  status = get_time_index ((*ncid), UTC, time_interpolate, nt, itime1, itime2, dt, verbose, quiet);
  CHKERR (status);

  if (verbose) {
    if (UTC.tm_mday > 0 && (*nt) == 2)
      fprintf (stderr,
               "     read data at: lat =%7.2f (%7.2f), lon =%7.2f (%7.2f), time = %s",
               (*lat_grid)[(*ilat)],
               lat,
               (*lon_grid)[(*ilon)],
               lon_tmp,
               asctime (&UTC));
    else
      fprintf (stderr,
               "     read data at: lat =%7.2f (%7.2f), lon =%7.2f (%7.2f), \n",
               (*lat_grid)[(*ilat)],
               lat,
               (*lon_grid)[(*ilon)],
               lon_tmp);

    if ((*nt) == 1)
      fprintf (stderr, "     element: i=%6ld (lon), j=%6ld (lat), t=%6d (time)\n", (*ilon), (*ilat), (*itime1));
    if ((*nt) == 2)
      fprintf (stderr,
               "     element: i=%6ld (lon), j=%6ld (lat), (%4d/%4d) (time1/time2)\n",
               (*ilon),
               (*ilat),
               (*itime1),
               (*itime2));
  }

  return status;

#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
  fprintf (stderr, " ***********************************************************************\n");
  CHKERR (-1);
#endif
}

/****************************************************************/
/* read_p_T_z_from_ECMWF_file                                   */
/*                                                              */
/* Purpose:                                                     */
/* read pressure, temperature from ECMWF file                   */
/* interpolate data with respect to time if wanted              */
/* calculate z-grid for the layers                              */
/*                                                              */
/*  input:                                                      */
/*  ------                                                      */
/*  ncid       id_number of the open netCDF file                */
/*  ilat       index for latitude                               */
/*  nlat       size of latitude in netCDF file                  */
/*  ilon       index for longitude                              */
/*  nlon       size of longitude in netCDF file                 */
/*  nt         number of time steps (2 if time interpolation)   */
/*  itime1     first  index for time                            */
/*  itime2     second index for time                            */
/*  dt         time step fraction where time                    */
/*             is in between time1 and time2                    */
/*  time_interpolate    switch for time interpolation           */
/*  verbose    additional verbose output                        */
/*  quiet      no verbose output                                */
/*                                                              */
/*  output:                                                     */
/*  -------                                                     */
/*  nlay       number of layers + 1 for surface                 */
/*  p_layer    layer averaged pressure                          */
/*  T_layer    layer averaged temperature                       */
/*  z_layer    z-levels for layer midpoints                     */
/*                                                              */
/*  September 2007  by Ulrich Hamann                            */
/****************************************************************/

int read_p_T_z_from_ECMWF_file (int     ncid,
                                long    ilat,
                                size_t  nlat,
                                long    ilon,
                                size_t  nlon,
                                int     nt,
                                int     itime1,
                                int     itime2,
                                float   dt,
                                size_t* nlay,
                                size_t* nlev,
                                float   altitude,
                                float** p_layer,
                                float** T_layer,
                                float** z_layer,
                                int     verbose,
                                int     quiet) {
  DEPRECATE (2023);
  int status = 0;

  float** tmp_p_level = NULL;
  float** tmp_p_layer = NULL;
  float** tmp_T_layer = NULL;

  int   itime = 0;
  int   t     = NOT_DEFINED_INTEGER;
  float SP    = -999.0;
  float dx    = NOT_DEFINED_FLOAT;

  int lc = 0;

  /* alloc nt timesteps for pressure */
  tmp_p_level = calloc (nt, sizeof (float*));
  CHKPOINTER (tmp_p_level);

  /* alloc nt timesteps for temperature */
  tmp_T_layer = calloc (nt, sizeof (float*));
  CHKPOINTER (tmp_T_layer);

  /* read nt (= 1 or 2) time steps */
  for (t = 0; t <= nt - 1; t++) {

    if (t == 0)
      itime = itime1;
    if (t == 1)
      itime = itime2;

    /* alloc and read pressure */
    /* requires that SP/LNSP, hyai and hybi are in the ECMWF netCDF file */ /* in ancillary.c */
    status =
      alloc_and_read_ECMWF_netCDF_pressure (ncid, &(tmp_p_level[t]), nlev, &(tmp_p_layer[t]), nlay, itime, ilat, ilon, FALSE);
    CHKERR (status);

    tmp_p_layer[t] = realloc (tmp_p_layer[t], ((*nlay) + 1) * sizeof (float));
    CHKPOINTER (tmp_p_layer[t]);

#if HAVE_NETCDF4
    /* read surface pressure from ECMWF file */
    status = read_ECMWF_surface_pressure (ncid, itime, ilat, ilon, &(SP), FALSE); /* FALSE == no verbose */
    CHKERR (status);
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use ECMWF input options. Please get netcdf and rebuild.             *\n");
    fprintf (stderr, " ***********************************************************************\n");
    CHKERR (-1);
#endif

    tmp_p_layer[t][(*nlay)] = SP / 100.0; /* Pa -> hPa */

    /* allocate and read temperature, defined on layers */
    status = alloc_and_read_netCDF_column_float (ncid, "T", &tmp_T_layer[t], (*nlay), itime, ilat, ilon, FALSE);
    CHKERR (status);

    tmp_T_layer[t] = realloc (tmp_T_layer[t], ((*nlay) + 1) * sizeof (float));
    CHKPOINTER (tmp_T_layer[t]);

    /* extrapolate temperature from last layer midpoint to surface */
    dx = (log (tmp_p_level[0][(*nlev) - 1] / tmp_p_level[0][(*nlev) - 1]) -
          log (tmp_p_layer[0][(*nlay) - 1] / tmp_p_level[0][(*nlev) - 1])) /
         (log (tmp_p_layer[0][(*nlay) - 1] / tmp_p_level[0][(*nlev) - 1]) -
          log (tmp_p_layer[0][(*nlay) - 2] / tmp_p_level[0][(*nlev) - 1]));
    tmp_T_layer[t][(*nlay)] = (1.0 + dx) * tmp_T_layer[0][(*nlay) - 1] - dx * tmp_T_layer[0][(*nlay) - 2];
  }

  /* number of final layers == number of layers plus one for surface */
  (*nlay) = (*nlay) + 1;

  if (verbose)
    fprintf (stderr, "     found %4zd levels,  \n", (*nlev));

  if (nt == 1) {
    /* no time interpolation needed, just copy data */
  } else {
    /* time interpolation */
    for (lc = 0; lc < (*nlay); lc++) {
      tmp_p_layer[0][lc] = (1.0 - dt) * tmp_p_layer[0][lc] + dt * tmp_p_layer[1][lc];
      tmp_T_layer[0][lc] = (1.0 - dt) * tmp_T_layer[0][lc] + dt * tmp_T_layer[1][lc];
    }
  }

  /* calculate z from T and p using hydrostatic equation (in ancillary.c) */
  status = calculate_z_from_p_and_T (tmp_p_layer[0], tmp_T_layer[0], z_layer, (*nlay), altitude, FALSE);
  CHKERR (status);

  /* allocate space for final results */
  *p_layer = calloc ((*nlay), sizeof (float));
  CHKPOINTER (*p_layer);

  /* allocate space for final results */
  *T_layer = calloc ((*nlay), sizeof (float));
  CHKPOINTER (*T_layer);

  for (lc = 0; lc < (*nlay); lc++) {
    (*p_layer)[lc] = tmp_p_layer[0][lc];
    (*T_layer)[lc] = tmp_T_layer[0][lc];
  }

  return status;
}

/****************************************************************/
/* get_number_from_netCDF_map                                   */
/*                                                              */
/* Purpose:                                                     */
/* search one entry (float)                                     */
/* by latitude, longitude, time from one netCDF file            */
/*                                                              */
/*  input:                                                      */
/*  ------                                                      */
/*  lat        latitude                                         */
/*  lon        longitude                                        */
/*  UTC        universal time correlated                        */
/*  filename   name of the netCDF file including the map        */
/*                                                              */
/*  output:                                                     */
/*  -------                                                     */
/*  result     data read from the map                           */
/*  February 2007  by Ulrich Hamann                             */
/****************************************************************/

int get_number_from_netCDF_map (float     lat,
                                float     lon,
                                struct tm UTC,
                                int       time_interpolate,
                                char*     filename,
                                void*     data,
                                int       external_type,
                                char*     variable_name,
                                int       verbose,
                                int       quiet) {
  int status = 0;

#if HAVE_NETCDF4

  int     ncid    = 0;
  int     id_data = 0;
  nc_type netCDF_type;

  char lat_name[FILENAME_MAX] = "";
  char lon_name[FILENAME_MAX] = "";

  size_t nlat = 0;
  size_t nlon = 0;

  int ilat  = NOT_DEFINED_INTEGER;
  int ilon  = NOT_DEFINED_INTEGER;
  int dummy = NOT_DEFINED_INTEGER;

  double* lat_grid = NULL;
  double* lon_grid = NULL;

  int   i = 0, j = 0, t = 0;
  int   nt     = -1;
  int   itime1 = NOT_DEFINED_INTEGER, itime2 = NOT_DEFINED_INTEGER;
  float dt = NOT_DEFINED_FLOAT;

  int dimensions = NOT_DEFINED_INTEGER;

  size_t* index = NULL;

  float lon_tmp = NOT_DEFINED_FLOAT;
  float lat_min = NOT_DEFINED_FLOAT, lat_max = NOT_DEFINED_FLOAT;
  float lon_min = NOT_DEFINED_FLOAT, lon_max = NOT_DEFINED_FLOAT;

  unsigned char* data_uchar  = NULL;
  short*         data_short  = NULL;
  int*           data_int    = NULL;
  float*         data_float  = NULL;
  double*        data_double = NULL;

  /* int    data_unit_len = NOT_DEFINED_INTEGER; */
  /* nc_type unit_type; */
  char data_unit[50] = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";

  float missing_value = NOT_DEFINED_FLOAT;
  float unit_factor   = 1.0;
  float scale_factor  = 1.0;
  float add_offset    = 0.0;

  if (verbose) {
    fprintf (stderr, " ... read %s from map: \n", variable_name);
    fprintf (stderr, "     %s \n", filename);
  }

  lon_tmp = lon;

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s' opening netCDF file %s\n", nc_strerror (status), filename);
    fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    CHKERR (status);
  }

  /* check format */
  status = nc_inq_varid (ncid, "lat", &id_data);
  if (status == NC_NOERR) {
    strcpy (lat_name, "lat");
    strcpy (lon_name, "lon");
  } else {
    status = nc_inq_varid (ncid, "latitude", &id_data);
    if (status == NC_NOERR) {
      strcpy (lat_name, "latitude");
      strcpy (lon_name, "longitude");
    } else {
      fprintf (stderr, "Error, neither 'lat' nor 'latitude' in %s while reading '%s'\n", filename, variable_name);
      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      CHKERR (status);
    }
  }

  /* get id for "data to read" (variable_name) */
  status = nc_inq_varid (ncid, variable_name, &id_data);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s', while getting ncid for %s from %s\n", nc_strerror (status), variable_name, filename);
    fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    CHKERR (status);
  }

  /* get type of the variable */
  status = nc_inq_vartype (ncid, id_data, &(netCDF_type));
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s', while getting type for %s from %s\n", nc_strerror (status), variable_name, filename);
    fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    CHKERR (status);
  }

  if (verbose)
    fprintf (stderr,
             "     data type = %d (NC_BYTE %d, NC_CHAR %d, NC_SHORT %d, NC_INT %d, NC_FLOAT %d, NC_DOUBLE %d)\n",
             netCDF_type,
             NC_BYTE,
             NC_CHAR,
             NC_SHORT,
             NC_INT,
             NC_FLOAT,
             NC_DOUBLE);

  /* read attribute unit and interprete it as much as possible */
  status = nc_get_att_text (ncid, id_data, "units", data_unit);

  if (status != NC_NOERR) {
    /* no given units -> default unit_factor = 1.0 */
    unit_factor = 1.0;
    /* that's OK, put status to 0 */
    status = 0;
  } else {
    /* some information is given about the unit, try to interprete it */
    if (strncasecmp ("per cent", data_unit, 8) == 0)
      unit_factor = 100.0; /* % -> 1  */
    if (strcasecmp ("m", data_unit) == 0)
      unit_factor = 1000.0; /* m -> km */

    if (strncasecmp ("meter", data_unit, 5) == 0)
      unit_factor = 1000.0; /* m -> km */
    if (strcasecmp ("m**2 s**-2", data_unit) == 0)
      unit_factor = 9.80 * 1000.0; /* gpm -> km */
    if (verbose)
      fprintf (stderr, "     data unit = %s, scale data with %f\n", data_unit, unit_factor);
  }
  CHKERR (status);

  /* read attribute scale_factor */
  status = nc_get_att_float (ncid, id_data, "scale_factor", &scale_factor);
  if (status != NC_NOERR) {
    /* no given scale_factor -> default scale_factor = 1.0 */
    /* that's OK, put status to 0 */
    status = 0;
  } else {
    /* if (verbose) */
    /*  fprintf (stderr,"     scale_factor = %f \n", scale_factor); */
  }
  CHKERR (status);

  /* read attribute add_offset */
  status = nc_get_att_float (ncid, id_data, "add_offset", &add_offset);
  if (status != NC_NOERR) {
    /* no given offset -> default offset = 0.0 */
    /* that's OK, put status to 0 */
    status = 0;
  } else {
    /* if (verbose) */
    /*   fprintf (stderr,"     add_offset = %f \n", add_offset); */
  }
  CHKERR (status);

  /* read attribute _FillValue / missing_value */
  status = nc_get_att_float (ncid, id_data, "missing_value", &missing_value);
  if (status != NC_NOERR) {
    /* no given missing_value -> try to read _FillValue */
    status = nc_get_att_float (ncid, id_data, "_FillValue", &missing_value);
    if (status != NC_NOERR) {
      /* no given attribute, that's OK, we can' make checks, but we continue */
      status = 0;
    } else {
      /* if (verbose) */
      /*  fprintf (stderr,"     missing_value = %f \n", missing_value); */
    }
  }
  CHKERR (status);

  /* read latitude array */
  status = alloc_and_read_netCDF_1D_double (ncid, lat_name, &nlat, lat_name, &(lat_grid));
  CHKERR (status);

  /* read longitude */
  status = alloc_and_read_netCDF_1D_double (ncid, lon_name, &nlon, lon_name, &(lon_grid));
  CHKERR (status);

  if (verbose)
    fprintf (stderr, "     found %zd (south-north) x %zd (west-east) data points\n", nlat, nlon);

  /* get range of latitude */
  lat_min = lat_grid[0];
  lat_max = lat_grid[0];
  for (i = 1; i < nlat; i++) {
    if (lat_grid[i] > lat_max)
      lat_max = lat_grid[i];
    if (lat_grid[i] < lat_min)
      lat_min = lat_grid[i];
  }

  /* get range of longitude */
  lon_min = lon_grid[0];
  lon_max = lon_grid[0];

  for (j = 1; j < nlon; j++) {
    if (lon_grid[j] > lon_max)
      lon_max = lon_grid[j];
    if (lon_grid[j] < lon_min)
      lon_min = lon_grid[j];
  }

  /* longitude periodicity */
  /* correct periodicity of latitude so, that longitude is inside the range of the map */
  if (lon_tmp < lon_min)
    lon_tmp += 360.0;
  if (lon_tmp > lon_max)
    lon_tmp -= 360.0;

  if (verbose) {
    fprintf (stderr,
             "     map size = [%8.3f (South),%8.3f (North)] x [%8.3f (West),%8.3f (East)]\n",
             lat_min,
             lat_max,
             lon_min,
             lon_max);
  }

  /* read time in netCDF file (if present)  */
  status = nc_inq_dimid (ncid, "time", &dummy);
  if (status == NC_NOERR) {

    /* get time index */
    status = get_time_index (ncid, UTC, time_interpolate, &(nt), &(itime1), &(itime2), &(dt), verbose, quiet);
    CHKERR (status);

    /* data in the netCDF file is assumed to have 3 dimensions: lat/lon/time */
    dimensions = 3;
  } else {
    /* no time variable in the netCDF file */
    status     = 0; /* it's OK, even if there is no time information */
    nt         = 1; /* take the first entry */
    dimensions = 2; /* lon / lat */
  }

  /* netCDF index might have 2 (lat/lon) or 3 (time/lat/lon) elements */
  if ((index = calloc (dimensions, sizeof (size_t))) == NULL) {
    fprintf (stderr, "Error allocation of memory for 'index'\n");
    fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -1;
  }

  /* start_2D[0] as type size_t can not be negative, that why first index is used here */
  ilat = get_grid_index (lat, lat_grid, nlat, FALSE);
  if (ilat < 0) {
    fprintf (stderr, "Error %d, while searching index in lat_grid \n", ilat);
    fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return ilat;
  }

  /* start_2D[1] as type size_t can not be negative, that why first index is used here */
  ilon = get_grid_index (lon_tmp, lon_grid, nlon, TRUE);
  if (ilon < 0) {
    fprintf (stderr, "Error %d, while searching index in lon_grid \n", ilon);
    fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return ilon;
  }

  /* alloc nt timesteps for data */
  switch (netCDF_type) {
  case (NC_BYTE):
  case (NC_CHAR):
    if ((data_uchar = calloc (nt, sizeof (unsigned char))) == NULL) {
      fprintf (stderr, "Error, allocating memory for data_intern\n");
      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -10;
    }
    break;
  case (NC_SHORT):
    if ((data_short = calloc (nt, sizeof (short))) == NULL) {
      fprintf (stderr, "Error, allocating memory for data_intern\n");
      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -10;
    }
    break;
  case (NC_INT):
    if ((data_int = calloc (nt, sizeof (int))) == NULL) {
      fprintf (stderr, "Error, allocating memory for data_intern\n");
      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -10;
    }
    break;
  case (NC_FLOAT):
    if ((data_float = calloc (nt, sizeof (float))) == NULL) {
      fprintf (stderr, "Error, allocating memory for data_intern\n");
      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -10;
    }
    break;
  case (NC_DOUBLE):
    if ((data_double = calloc (nt, sizeof (double))) == NULL) {
      fprintf (stderr, "Error, allocating memory for data_intern\n");
      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -10;
    }
    break;
  default:
    fprintf (stderr, "Error, unknown type (short, float, double ...) of variable %d\n", netCDF_type);
    fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -1;
  }

  if (dimensions == 2) {

    index[0] = ilat;
    index[1] = ilon;

    if (verbose) {
      fprintf (stderr, "     read data at: lat =%7.2f (%7.2f), lon =%7.2f (%7.2f), ", lat_grid[ilat], lat, lon_grid[ilon], lon_tmp);
      fprintf (stderr, " element: (j=%6d) x (i=%6d) \n", ilat, ilon);
    }

    switch (netCDF_type) {
    case (NC_BYTE):
    case (NC_CHAR):
      status = nc_get_var1_uchar (ncid, id_data, index, &(data_uchar[0]));
      /* fprintf (stderr, "     read data value is %d (unscaled)\n", (data_uchar[0])); */
      if (data_uchar[0] == missing_value)
        status = ERROR_READ_MISSING_VALUE;
      break;
    case (NC_SHORT):
      status = nc_get_var1_short (ncid, id_data, index, &(data_short[0]));
      /* fprintf (stderr, "     read data value is %d (unscaled)\n", (data_short[0])); */
      if (data_short[0] == missing_value)
        status = ERROR_READ_MISSING_VALUE;
      break;
    case (NC_INT):
      status = nc_get_var1_int (ncid, id_data, index, &(data_int[0]));
      /* fprintf (stderr, "     read data value is %d (unscaled)\n", (data_int[0])); */
      if (data_int[0] == missing_value)
        status = ERROR_READ_MISSING_VALUE;
      break;
    case (NC_FLOAT):
      status = nc_get_var1_float (ncid, id_data, index, &(data_float[0]));
      /* fprintf (stderr, "     read data value is %f (unscaled)\n", (data_float[0])); */
      if (data_float[0] == missing_value)
        status = ERROR_READ_MISSING_VALUE;
      break;
    case (NC_DOUBLE):
      status = nc_get_var1_double (ncid, id_data, index, &(data_double[0]));
      /* fprintf (stderr, "     read data value is %lf (unscaled)\n", (data_double[0])); */
      if (data_double[0] == missing_value)
        status = ERROR_READ_MISSING_VALUE;
      break;
    default:
      fprintf (stderr,
               "Error, unknown variable type %d (short, float, ...) in data file %s (line %d, function %s in %s) \n",
               netCDF_type,
               filename,
               __LINE__,
               __func__,
               __FILE__);
      return -1;
    }

    if (status == ERROR_READ_MISSING_VALUE) {
      fprintf (stderr, " !!! Error 'Read missing_value' reading %s from %s\n", variable_name, filename);
      CHKERR (status);
    }

    if (status != NC_NOERR) {
      fprintf (stderr, "Error '%s'\n", nc_strerror (status));
      fprintf (stderr, "      reading %s from %s\n", variable_name, filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      CHKERR (status);
    }

    switch (netCDF_type) {
    case (NC_BYTE):
    case (NC_CHAR):
      switch (external_type) {
      case (TYPE_CHAR):
        (*(char*) data) = (char) (data_uchar[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_SHORT):
        (*(short*) data) = (short) (data_uchar[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_INT):
        (*(int*) data) = (int) (data_uchar[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_FLOAT):
        (*(float*) data) = (float) (data_uchar[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_DOUBLE):
        (*(double*) data) = (double) (data_uchar[0] * scale_factor + add_offset) / unit_factor;
        break;
      }
      break;
    case (NC_SHORT):
      switch (external_type) {
      case (TYPE_CHAR):
        (*(char*) data) = (char) (data_short[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_SHORT):
        (*(short*) data) = (short) (data_short[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_INT):
        (*(int*) data) = (int) (data_short[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_FLOAT):
        (*(float*) data) = (float) (data_short[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_DOUBLE):
        (*(double*) data) = (double) (data_short[0] * scale_factor + add_offset) / unit_factor;
        break;
      }
      break;
    case (NC_INT):
      switch (external_type) {
      case (TYPE_CHAR):
        (*(char*) data) = (char) (data_int[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_SHORT):
        (*(short*) data) = (short) (data_int[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_INT):
        (*(int*) data) = (int) (data_int[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_FLOAT):
        (*(float*) data) = (float) (data_int[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_DOUBLE):
        (*(double*) data) = (double) (data_int[0] * scale_factor + add_offset) / unit_factor;
        break;
      }
      break;
    case (NC_FLOAT):
      switch (external_type) {
      case (TYPE_CHAR):
        (*(char*) data) = (char) (data_float[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_SHORT):
        (*(short*) data) = (short) (data_float[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_INT):
        (*(int*) data) = (int) (data_float[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_FLOAT):
        (*(float*) data) = (float) (data_float[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_DOUBLE):
        (*(double*) data) = (double) (data_float[0] * scale_factor + add_offset) / unit_factor;
        break;
      }
      break;
    case (NC_DOUBLE):
      switch (external_type) {
      case (TYPE_CHAR):
        (*(char*) data) = (char) (data_double[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_SHORT):
        (*(short*) data) = (short) (data_double[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_INT):
        (*(int*) data) = (int) (data_double[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_FLOAT):
        (*(float*) data) = (float) (data_double[0] * scale_factor + add_offset) / unit_factor;
        break;
      case (TYPE_DOUBLE):
        (*(double*) data) = (double) (data_double[0] * scale_factor + add_offset) / unit_factor;
        break;
      }
      break;
    default:
      fprintf (stderr, "Error, type of variable %d in get_float_from_netCDF_map (ancillary.c) \n", netCDF_type);
      return -1;
    }

    switch (external_type) {
    case (TYPE_CHAR):
      if (verbose)
        fprintf (stderr, "     read data value is %d\n", (*(char*) data));
      break;
    case (TYPE_SHORT):
      if (verbose)
        fprintf (stderr, "     read data value is %d\n", (*(short*) data));
      break;
    case (TYPE_INT):
      if (verbose)
        fprintf (stderr, "     read data value is %d\n", (*(int*) data));
      break;
    case (TYPE_FLOAT):
      if (verbose)
        fprintf (stderr, "     read data value is %f\n", (*(float*) data));
      break;
    case (TYPE_DOUBLE):
      if (verbose)
        fprintf (stderr, "     read data value is %lf\n", (*(double*) data));
      break;
    }
  } else if (dimensions == 3) {

    index[1] = ilat;
    index[2] = ilon;

    if (verbose) {
      if (UTC.tm_mday > 0) {
        fprintf (stderr,
                 "     read data at: lat =%7.2f (%7.2f), lon =%7.2f (%7.2f), time = %s",
                 lat_grid[ilat],
                 lat,
                 lon_grid[ilon],
                 lon_tmp,
                 asctime (&UTC));
        if (nt == 1)
          fprintf (stderr, "     element: (j=%6d) x (i=%6d) x (t=%6d)\n", ilat, ilon, itime1);
        if (nt == 2)
          fprintf (stderr, "     element: (j=%6d) x (i=%6d) x (t1/t2)=(%4d/%4d)\n", ilat, ilon, itime1, itime2);
      } else {
        fprintf (stderr,
                 "     read data at: lat =%7.2f (%7.2f), lon =%7.2f (%7.2f)\n",
                 lat_grid[ilat],
                 lat,
                 lon_grid[ilon],
                 lon_tmp);
        fprintf (stderr, "     element: (j=%6d) x (i=%6d) x (t=%6d)\n", ilat, ilon, itime1);
      }
    }

    /* read nt (= 1 or 2) time steps */
    for (t = 0; t <= nt - 1; t++) {

      if (t == 0)
        index[0] = itime1;
      if (t == 1)
        index[0] = itime2;

      switch (netCDF_type) {
      case (NC_BYTE):
      case (NC_CHAR):
        status = nc_get_var1_uchar (ncid, id_data, index, &(data_uchar[t]));
        /* fprintf (stderr, "     read data value is %d (unscaled)\n", (data_uchar[t])); */
        if (data_uchar[t] == missing_value)
          status = ERROR_READ_MISSING_VALUE;
        break;
      case (NC_SHORT):
        status = nc_get_var1_short (ncid, id_data, index, &(data_short[t]));
        /* fprintf (stderr, "     read data value is %d (unscaled)\n", (data_short[t])); */
        if (data_short[t] == missing_value)
          status = ERROR_READ_MISSING_VALUE;
        break;
      case (NC_INT):
        status = nc_get_var1_int (ncid, id_data, index, &(data_int[t]));
        /* fprintf (stderr, "     read data value is %d (unscaled)\n", (data_int[t])); */
        if (data_int[t] == missing_value)
          status = ERROR_READ_MISSING_VALUE;
        break;
      case (NC_FLOAT):
        status = nc_get_var1_float (ncid, id_data, index, &(data_float[t]));
        /* fprintf (stderr, "     read data value is %f (unscaled)\n", (data_float[t])); */
        if (data_float[t] == missing_value)
          status = ERROR_READ_MISSING_VALUE;
        break;
      case (NC_DOUBLE):
        status = nc_get_var1_double (ncid, id_data, index, &(data_double[t]));
        /* fprintf (stderr, "     read data value is %lf (unscaled)\n", (data_double[t])); */
        if (data_double[t] == missing_value)
          status = ERROR_READ_MISSING_VALUE;
        break;
      default:
        fprintf (stderr,
                 "Error, unknown variable type %d (short, float, ...) in data file %s (line %d, function %s in %s)  \n",
                 netCDF_type,
                 filename,
                 __LINE__,
                 __func__,
                 __FILE__);
        return -1;
      }

      if (status == ERROR_READ_MISSING_VALUE) {
        fprintf (stderr, "Error 'Read missing_value' reading %s from %s\n", variable_name, filename);
        CHKERR (status);
      }

      if (status != NC_NOERR) {
        fprintf (stderr, "Error '%s'\n", nc_strerror (status));
        fprintf (stderr, "      reading %s from %s\n", variable_name, filename);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        CHKERR (status);
      }
    }

    if (nt == 1) {
      /* no time interpolation nessesary */
      switch (netCDF_type) {
      case (NC_BYTE):
      case (NC_CHAR):
        switch (external_type) {
        case (TYPE_CHAR):
          (*(char*) data) = (char) (data_uchar[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_SHORT):
          (*(short*) data) = (short) (data_uchar[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_INT):
          (*(int*) data) = (int) (data_uchar[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_FLOAT):
          (*(float*) data) = (float) (data_uchar[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_DOUBLE):
          (*(double*) data) = (double) (data_uchar[0] * scale_factor + add_offset) / unit_factor;
          break;
        }
        break;
      case (NC_SHORT):
        switch (external_type) {
        case (TYPE_CHAR):
          (*(char*) data) = (char) (data_short[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_SHORT):
          (*(short*) data) = (short) (data_short[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_INT):
          (*(int*) data) = (int) (data_short[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_FLOAT):
          (*(float*) data) = (float) (data_short[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_DOUBLE):
          (*(double*) data) = (double) (data_short[0] * scale_factor + add_offset) / unit_factor;
          break;
        }
        break;
      case (NC_INT):
        switch (external_type) {
        case (TYPE_CHAR):
          (*(char*) data) = (char) (data_int[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_SHORT):
          (*(short*) data) = (short) (data_int[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_INT):
          (*(int*) data) = (int) (data_int[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_FLOAT):
          (*(float*) data) = (float) (data_int[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_DOUBLE):
          (*(double*) data) = (double) (data_int[0] * scale_factor + add_offset) / unit_factor;
          break;
        }
        break;
      case (NC_FLOAT):
        switch (external_type) {
        case (TYPE_CHAR):
          (*(char*) data) = (char) (data_float[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_SHORT):
          (*(short*) data) = (short) (data_float[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_INT):
          (*(int*) data) = (int) (data_float[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_FLOAT):
          (*(float*) data) = (float) (data_float[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_DOUBLE):
          (*(double*) data) = (double) (data_float[0] * scale_factor + add_offset) / unit_factor;
          break;
        }
        break;
      case (NC_DOUBLE):
        switch (external_type) {
        case (TYPE_CHAR):
          (*(char*) data) = (char) (data_double[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_SHORT):
          (*(short*) data) = (short) (data_double[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_INT):
          (*(int*) data) = (int) (data_double[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_FLOAT):
          (*(float*) data) = (float) (data_double[0] * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_DOUBLE):
          (*(double*) data) = (double) (data_double[0] * scale_factor + add_offset) / unit_factor;
          break;
        }
        break;
      default:
        fprintf (stderr, "Error, type of variable %d in get_float_from_netCDF_map (ancillary.c) \n", netCDF_type);
        return -1;
      }
      if (verbose)
        fprintf (stderr, "     read data value is");
    } else {
      /* time interpolation */
      switch (netCDF_type) {
      case (NC_BYTE):
      case (NC_CHAR):
        switch (external_type) {
        case (TYPE_CHAR):
          (*(char*) data) = (char) (((1.0 - dt) * data_uchar[0] + dt * data_uchar[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_SHORT):
          (*(short*) data) = (short) (((1.0 - dt) * data_uchar[0] + dt * data_uchar[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_INT):
          (*(int*) data) = (int) (((1.0 - dt) * data_uchar[0] + dt * data_uchar[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_FLOAT):
          (*(float*) data) = (float) (((1.0 - dt) * data_uchar[0] + dt * data_uchar[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_DOUBLE):
          (*(double*) data) =
            (double) (((1.0 - dt) * data_uchar[0] + dt * data_uchar[1]) * scale_factor + add_offset) / unit_factor;
          break;
        }
        break;
      case (NC_SHORT):
        switch (external_type) {
        case (TYPE_CHAR):
          (*(char*) data) = (char) (((1.0 - dt) * data_short[0] + dt * data_short[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_SHORT):
          (*(short*) data) = (short) (((1.0 - dt) * data_short[0] + dt * data_short[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_INT):
          (*(int*) data) = (int) (((1.0 - dt) * data_short[0] + dt * data_short[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_FLOAT):
          (*(float*) data) = (float) (((1.0 - dt) * data_short[0] + dt * data_short[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_DOUBLE):
          (*(double*) data) =
            (double) (((1.0 - dt) * data_short[0] + dt * data_short[1]) * scale_factor + add_offset) / unit_factor;
          break;
        }
        break;
      case (NC_INT):
        switch (external_type) {
        case (TYPE_CHAR):
          (*(char*) data) = (char) (((1.0 - dt) * data_int[0] + dt * data_int[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_SHORT):
          (*(short*) data) = (short) (((1.0 - dt) * data_int[0] + dt * data_int[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_INT):
          (*(int*) data) = (int) (((1.0 - dt) * data_int[0] + dt * data_int[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_FLOAT):
          (*(float*) data) = (float) (((1.0 - dt) * data_int[0] + dt * data_int[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_DOUBLE):
          (*(double*) data) = (double) (((1.0 - dt) * data_int[0] + dt * data_int[1]) * scale_factor + add_offset) / unit_factor;
          break;
        }
        break;
      case (NC_FLOAT):
        switch (external_type) {
        case (TYPE_CHAR):
          (*(char*) data) = (char) (((1.0 - dt) * data_float[0] + dt * data_float[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_SHORT):
          (*(short*) data) = (short) (((1.0 - dt) * data_float[0] + dt * data_float[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_INT):
          (*(int*) data) = (int) (((1.0 - dt) * data_float[0] + dt * data_float[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_FLOAT):
          (*(float*) data) = (float) (((1.0 - dt) * data_float[0] + dt * data_float[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_DOUBLE):
          (*(double*) data) =
            (double) (((1.0 - dt) * data_float[0] + dt * data_float[1]) * scale_factor + add_offset) / unit_factor;
          break;
        }
        break;
      case (NC_DOUBLE):
        switch (external_type) {
        case (TYPE_CHAR):
          (*(char*) data) = (char) (((1.0 - dt) * data_double[0] + dt * data_double[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_SHORT):
          (*(short*) data) =
            (short) (((1.0 - dt) * data_double[0] + dt * data_double[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_INT):
          (*(int*) data) = (int) (((1.0 - dt) * data_double[0] + dt * data_double[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_FLOAT):
          (*(float*) data) =
            (float) (((1.0 - dt) * data_double[0] + dt * data_double[1]) * scale_factor + add_offset) / unit_factor;
          break;
        case (TYPE_DOUBLE):
          (*(double*) data) =
            (double) (((1.0 - dt) * data_double[0] + dt * data_double[1]) * scale_factor + add_offset) / unit_factor;
          break;
        }
        break;
      default:
        fprintf (stderr, "Error, type of variable %d in get_float_from_netCDF_map (ancillary.c) \n", netCDF_type);
        return -1;
      }

      if (verbose)
        fprintf (stderr, "     interpolated data value is");
    }

    if (verbose) {
      switch (external_type) {
      case (TYPE_CHAR):
        fprintf (stderr, " %d\n", (*(char*) data));
        break;
      case (TYPE_SHORT):
        fprintf (stderr, " %d\n", (*(short*) data));
        break;
      case (TYPE_INT):
        fprintf (stderr, " %d\n", (*(int*) data));
        break;
      case (TYPE_FLOAT):
        fprintf (stderr, " %f\n", (*(float*) data));
        break;
      case (TYPE_DOUBLE):
        fprintf (stderr, " %lf\n", (*(double*) data));
        break;
      }
    }
  }

  nc_close (ncid);

  free (lat_grid);
  free (lon_grid);
  free (index);

  switch (netCDF_type) {
  case (NC_BYTE):
  case (NC_CHAR):
    free (data_uchar);
    break;
  case (NC_SHORT):
    free (data_short);
    break;
  case (NC_INT):
    free (data_int);
    break;
  case (NC_FLOAT):
    free (data_float);
    break;
  case (NC_DOUBLE):
    free (data_double);
    break;
  default:
    fprintf (stderr,
             "Error, unknown variable type %d (short, float, ...) in data file %s (line %d, function %s in %s)  \n",
             netCDF_type,
             filename,
             __LINE__,
             __func__,
             __FILE__);
    return -1;
  }

#else
  fprintf (stderr, " ******************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot       *\n");
  fprintf (stderr, " * use the any map option. Please get netcdf and rebuild.         *\n");
  fprintf (stderr, " ******************************************************************\n");
  return -1;
#endif

  return status;
}

/*****************************************************************/
/* read satellite geometry                                       */
/*                                                               */
/* Purpose:                                                      */
/* read satellite geometry from netCDF file                      */
/*                                                               */
/*  input:                                                       */
/*  ------                                                       */
/*  sat_pixel_x      pixel index in x direction                  */
/*  sat_pixel_y      pixel index in y direction                  */
/*  UTC              simulated time (universal time correlated)  */
/*  filename         name of netCDF file (contain sat geometry)  */
/*                                                               */
/*  output:                                                      */
/*  -------                                                      */
/*  latitude         latitude                                    */
/*  longitude        longitude                                   */
/*  numu             number of cos(theta_sat) angles             */
/*  maxumu           number of cos(theta_sat) angles             */
/*  umu              cos(theta_sat) angles                       */
/*  nphi             number of azimith_sat angles                */
/*  maxphi           number of azimith_sat angles                */
/*  phi              azimith_sat angles in degrees               */
/*                                                               */
/*  Mai 2007    by Ulrich Hamann                                 */
/*****************************************************************/

int read_sat_geometry (int       pixel_x,
                       int       pixel_y,
                       struct tm UTC,
                       char*     filename,
                       float*    latitude,
                       float*    longitude,
                       int*      numu,
                       int*      maxumu,
                       float**   umu,
                       int*      nphi,
                       int*      maxphi,
                       float**   phi,
                       int       solver,
                       int       verbose,
                       int       quiet) {

#if HAVE_NETCDF4

  int     status           = 0;
  int     ncid             = NOT_DEFINED_INTEGER;
  int     id_lat           = NOT_DEFINED_INTEGER;
  int     id_lon           = NOT_DEFINED_INTEGER;
  int     id_vza           = NOT_DEFINED_INTEGER;
  int     id_vaa           = NOT_DEFINED_INTEGER;
  size_t  n_pixel_x        = NOT_DEFINED_INTEGER;
  size_t  n_pixel_y        = NOT_DEFINED_INTEGER;
  double* sat_pixel_x_grid = NULL;
  double* sat_pixel_y_grid = NULL;
  int     ix               = NOT_DEFINED_INTEGER;
  int     iy               = NOT_DEFINED_INTEGER;
  int     dummy            = NOT_DEFINED_INTEGER;
  size_t* index            = NULL;
  int     nt               = -1;
  int     itime1 = NOT_DEFINED_INTEGER, itime2 = NOT_DEFINED_INTEGER;
  float   dt         = NOT_DEFINED_FLOAT;
  int     dimensions = NOT_DEFINED_INTEGER;

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s' opening netCDF file %s\n", nc_strerror (status), filename);
    CHKERR (status);
  }

  /* get id for latitude */
  status = nc_inq_varid (ncid, "lat", &id_lat);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s', while getting id for 'latitude' from %s\n", nc_strerror (status), filename);
    CHKERR (status);
  }
  /* get id for longitude */
  status = nc_inq_varid (ncid, "lon", &id_lon);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s', while getting id for 'longitude' from %s\n", nc_strerror (status), filename);
    CHKERR (status);
  }
  /* get id for viewing zenith angle */
  status = nc_inq_varid (ncid, "vza", &id_vza);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s', while getting id for 'viewing zenith angle' from %s\n", nc_strerror (status), filename);
    CHKERR (status);
  }
  /* get id for viewing azimuth angle */
  status = nc_inq_varid (ncid, "vaa", &id_vaa);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s', while getting id for 'viewing azimuth angle' from %s\n", nc_strerror (status), filename);
    CHKERR (status);
  }

  /* read pixel_x array (must already be stored in the result file) */
  status = alloc_and_read_netCDF_1D_double (ncid, "pixel_x", &n_pixel_x, "pixel_x", &(sat_pixel_x_grid));
  CHKERR (status);

  /* read pixel_y array (must already be stored in the file) */
  status = alloc_and_read_netCDF_1D_double (ncid, "pixel_y", &n_pixel_y, "pixel_y", &(sat_pixel_y_grid));
  CHKERR (status);

  /* search index */
  ix = get_grid_index (pixel_x, sat_pixel_x_grid, n_pixel_x, FALSE);
  if (ix < 0) {
    fprintf (stderr, "Error %d, while searching index for pixel %6d in pixel-x_grid \n", ix, pixel_x);
    return ix;
  }

  iy = get_grid_index (pixel_y, sat_pixel_y_grid, n_pixel_y, FALSE);
  if (iy < 0) {
    fprintf (stderr, "Error %d, while searching index for pixel %6d in pixel-y_grid \n", iy, pixel_y);
    return iy;
  }

  /* read time in netCDF file (if present)  */
  status = nc_inq_dimid (ncid, "time", &dummy);
  if (status == NC_NOERR) {

    /* get time index */
    status = get_time_index (ncid, UTC, TIME_NEAREST_DATE, &(nt), &(itime1), &(itime2), &(dt), verbose, quiet);
    CHKERR (status);

    /* data in the netCDF file is assumed to have 3 dimensions: lat/lon/time */
    dimensions = 3;
  } else {
    /* no time variable in the netCDF file */
    nt         = 1;
    dimensions = 2;
  }

  /* netCDF index might have 2 (lat/lon) or 3 (time/lat/lon) elements */
  index = calloc (dimensions, sizeof (size_t));
  CHKPOINTER (index);

  (*numu)   = 1;
  (*maxumu) = 1;
  (*nphi)   = 1;

  switch (solver) {
  case SOLVER_SDISORT:
  case SOLVER_SPSDISORT:
  case SOLVER_FDISORT1:
  case SOLVER_FDISORT2:
    (*maxphi) = 3; /* Minimum number required by disort. */
    break;
  case SOLVER_DISORT:
  case SOLVER_FTWOSTR:
  case SOLVER_SOS:
  case SOLVER_MONTECARLO:
  case SOLVER_TZS:
  case SOLVER_SSS:
  case SOLVER_SSSI:
  case SOLVER_NULL:
  case SOLVER_RODENTS:
  case SOLVER_SCHWARZSCHILD:
  case SOLVER_TWOSTREBE:
  case SOLVER_TWOMAXRND:
  case SOLVER_TWOMAXRND3C:
  case SOLVER_DYNAMIC_TWOSTREAM:
  case SOLVER_DYNAMIC_TENSTREAM:
  case SOLVER_TWOSTR:
  case SOLVER_SSLIDAR:
    (*maxphi) = 1;
    break;
  default:
    fprintf (stderr, "Error, unknown rte_solver %d in read_sat_geometry\n", solver);
    return -1;
  }

  *umu = calloc ((*numu), sizeof (float));
  CHKPOINTER ((*umu));

  *phi = calloc ((*nphi), sizeof (float));
  CHKPOINTER ((*phi));

  if (dimensions == 2) {
    index[0] = iy;
    index[1] = ix;

    if (verbose) {
      fprintf (stderr,
               "*** read satellite geometry for pixel_x =%5.0f (%5d), pixel_y =%7.2f (%5d), ",
               sat_pixel_x_grid[ix],
               pixel_x,
               sat_pixel_y_grid[iy],
               pixel_y);
      fprintf (stderr, " element: %6d x %6d \n", ix + 1, iy + 1);
    }

  } else if (dimensions == 3) {

    index[0] = itime1;
    index[1] = iy;
    index[2] = ix;

    if (verbose) {
      if (UTC.tm_mday > 0) {
        fprintf (stderr,
                 " ... read satellite geometry for pixel_x =%6.0f (%5d), lon =%6.0f (%5d), time = %s",
                 sat_pixel_x_grid[ix],
                 pixel_x,
                 sat_pixel_y_grid[iy],
                 pixel_y,
                 asctime (&UTC));
        if (nt == 1)
          fprintf (stderr, "     element: %6d x %6d x %6d\n", ix + 1, iy + 1, itime1 + 1);
        if (nt == 2)
          fprintf (stderr, "     element: %6d x %6d x (%4d/%4d)\n", ix + 1, iy + 1, itime1 + 1, itime2 + 1);
      } else {
        fprintf (stderr,
                 " ... read satellite geometry for pixel_x =%6.0f (%5d), lon =%6.0f (%5d)\n",
                 sat_pixel_x_grid[ix],
                 pixel_x,
                 sat_pixel_y_grid[iy],
                 pixel_y);
        fprintf (stderr, "     element: %6d x %6d x %6d\n", ix + 1, iy + 1, itime1 + 1);
      }
    }

  } else {
    fprintf (stderr, "Error, wrong number for dimensions %d in read_sat_geometry (ancillary.c)\n", dimensions);
    return -1;
  }

  /* fprintf (stderr,"index= [%d, %d], status = %d/%d \n", index[0], index[1], status, NC_NOERR ); */

  status = nc_get_var1_float (ncid, id_lat, index, latitude);
  CHKERR (status);
  status = nc_get_var1_float (ncid, id_lon, index, longitude);
  CHKERR (status);
  status = nc_get_var1_float (ncid, id_vza, index, &((*umu)[0]));
  CHKERR (status);
  status = nc_get_var1_float (ncid, id_vaa, index, &((*phi)[0]));
  CHKERR (status);

  if (verbose)
    fprintf (stderr,
             "    latitude: %9.5f, longitude: %9.5f, theta_sat=%8.3f, phi_sat=%8.3f \n",
             (*latitude),
             (*longitude),
             (*umu)[0],
             (*phi)[0]);

  (*umu)[0] = cos (M_PI / 180 * (*umu)[0]);

  return status;

#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use the satellite_geometry option. Please get netcdf and rebuild.   *\n");
  fprintf (stderr, " ***********************************************************************\n");
  return -1;
#endif
}

/****************************************************************/
/* read_u10_from_map                                            */
/*                                                              */
/* Purpose:                                                     */
/* read wind in 10m height from an ECMWF file                   */
/*                                                              */
/*  input:                                                      */
/*  ------                                                      */
/*  lat        latitude                                         */
/*  lon        longitude                                        */
/*  UTC        universal time correlated                        */
/*  filename   name of the netCDF file including the wind map   */
/*                                                              */
/*  output:                                                     */
/*  -------                                                     */
/*  u10        wind velocity (vector norm) in m/s               */
/*                                                              */
/*  September 2007  by Ulrich Hamann                            */
/****************************************************************/

int read_u10_from_map (float     lat,
                       float     lon,
                       struct tm UTC,
                       int       time_interpolate,
                       char*     filename,
                       float*    u10,
                       int       verbose,
                       int       quiet) {

#if HAVE_NETCDF4

  int status = 0;

  float u_10 = 0.0;
  float v_10 = 0.0;

  float* u      = NULL;
  float* v      = NULL;
  float* w      = NULL;
  float* z_wind = NULL;
  size_t nlev   = 0;

  /* read u in 10m height */
  status = get_number_from_netCDF_map (lat, lon, UTC, time_interpolate, filename, &(u_10), TYPE_FLOAT, "U10", verbose, quiet);
  CHKERR (status);

  /* read v in 10m height */
  status = get_number_from_netCDF_map (lat, lon, UTC, time_interpolate, filename, &(v_10), TYPE_FLOAT, "V10", verbose, quiet);
  CHKERR (status);

  if (status == 0) {

    /* everything OK, calculate norm of the vector */
    (*u10) = sqrt (u_10 * u_10 + v_10 * v_10);
    CHKERR (status);

  } else {

    if (verbose) {
      fprintf (stderr, " ... didn't find u10 and v10 in netCDF file %s \n", filename);
      fprintf (stderr, "     try to read profiles u(z) and v(z) \n");
    }

    status = read_wind_from_ECMWF_file (lat,
                                        lon,
                                        UTC,
                                        time_interpolate,
                                        filename,
                                        0.0,
                                        &(u),
                                        &(v),
                                        &(w),
                                        &(z_wind),
                                        &(nlev),
                                        verbose,
                                        quiet);
    CHKERR (status);

    if (verbose)
      fprintf (stderr,
               "     take data from lowest level, lc=%3zd, z=%7.3f km, u=%7.2f m/s, v=%7.2f m/s \n",
               (nlev - 2),
               z_wind[nlev - 2],
               u[nlev - 2],
               v[nlev - 2]);

    /* (nlev-1) = surface, where u=0, v=0; (nlev-2)=last layer before surface, in ECMWF data approx 20m */
    (*u10) = sqrt (u[nlev - 2] * u[nlev - 2] + v[nlev - 2] * v[nlev - 2]);

    free (u);
    free (v);
    free (z_wind);
  }

  return status;

#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use the cox_and_munk_u10_map option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ***********************************************************************\n");
  CHKERR (-1);
#endif
}

/****************************************************************/
/* read_wind_from_ECMWF_file                                    */
/*                                                              */
/* Purpose:                                                     */
/* read wind components u and v from an ECMWF file              */
/*                                                              */
/*  input:                                                      */
/*  ------                                                      */
/*  lat        latitude                                         */
/*  lon        longitude                                        */
/*  UTC        universal time correlated                        */
/*  filename   name of the netCDF file including the wind map   */
/*                                                              */
/*  output:                                                     */
/*  -------                                                     */
/*  u-profil   wind (west-east component)                       */
/*  v-profil   wind (south north component)                     */
/*  z_wind     height scale for u and v                         */
/*  nlev_wind  number of levels for u, v, and z                 */
/*                                                              */
/*  September 2007  by Ulrich Hamann                            */
/****************************************************************/

int read_wind_from_ECMWF_file (float     lat,
                               float     lon,
                               struct tm UTC,
                               int       time_interpolate,
                               char*     filename,
                               float     altitude,
                               float**   u,
                               float**   v,
                               float**   w,
                               float**   z_wind,
                               size_t*   nlev,
                               int       verbose,
                               int       quiet) {

  int status = 0;

#if HAVE_NETCDF4

  int ncid = NOT_DEFINED_INTEGER;

  int  lc   = NOT_DEFINED_INTEGER;
  int  t    = NOT_DEFINED_INTEGER;
  int  nt   = -1;
  long ilat = NOT_DEFINED_INTEGER, ilon = NOT_DEFINED_INTEGER,
       itime = NOT_DEFINED_INTEGER; /* index for lat, lon, time in netCDF file */

  size_t nlat = 0;
  size_t nlon = 0;

  double* ECMWF_lat = NULL;
  double* ECMWF_lon = NULL;

  int   itime1 = -1, itime2 = -1;
  float dt = NOT_DEFINED_FLOAT;

  size_t tmp_nlev = 0;
  size_t nlay     = 0;

  float** p_level = NULL;
  float** p_layer = NULL;

  float** T_layer = NULL;
  float** u_layer = NULL;
  float** v_layer = NULL;
  float** w_layer = NULL;

  int  id_var_test = 0;
  char T_name[2]   = "";
  char U_name[2]   = "";
  char V_name[2]   = "";
  char W_name[2]   = "";

  float dx = NOT_DEFINED_FLOAT;

  float SP = -999.0;

  if (verbose)
    fprintf (stderr, " ... read wind profiles from netCDF file %s \n", filename);

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d opening netCDF file %s\n", status, filename);
    CHKERR (status);
  }

  /* read latitude array */
  status = alloc_and_read_netCDF_1D_double (ncid, "lat", &nlat, "lat", &(ECMWF_lat));
  CHKERR (status);

  /* search correct latitude index */
  ilat = get_grid_index (lat, ECMWF_lat, nlat, FALSE);
  if (ilat < 0) {
    fprintf (stderr, "Error -1 finding index for lat=%5.2f in %s\n", lat, filename);
    CHKERR (-1);
  }

  /* read longitude */
  status = alloc_and_read_netCDF_1D_double (ncid, "lon", &nlon, "lon", &(ECMWF_lon));
  CHKERR (status);

  /* search correct longitude index */
  ilon = get_grid_index (lon, ECMWF_lon, nlon, TRUE);
  if (ilon < 0) {
    fprintf (stderr, "Error -2 finding index for lon=%5.2f in %s\n", lon, filename);
    CHKERR (-2);
  }

  /* get time index */
  status = get_time_index (ncid, UTC, time_interpolate, &(nt), &(itime1), &(itime2), &(dt), verbose, quiet);
  CHKERR (status);

  /* alloc nt timesteps for pressure */
  p_level = calloc (nt, sizeof (float*));
  CHKPOINTER (p_level);

  /* alloc nt timesteps for pressure */
  p_layer = calloc (nt, sizeof (float*));
  CHKPOINTER (p_layer);

  /* alloc nt timesteps for temperature */
  T_layer = calloc (nt, sizeof (float*));
  CHKPOINTER (T_layer);

  /* alloc nt timesteps for west-east wind component u */
  u_layer = calloc (nt, sizeof (float*));
  CHKPOINTER (u_layer);

  /* alloc nt timesteps for south north wind component v */
  v_layer = calloc (nt, sizeof (float*));
  CHKPOINTER (v_layer);

  /* alloc nt timesteps for vertical wind component w */
  w_layer = calloc (nt, sizeof (float*));
  CHKPOINTER (w_layer);

  /* read nt (= 1 or 2) time steps */
  for (t = 0; t <= nt - 1; t++) {

    if (t == 0)
      itime = itime1;
    if (t == 1)
      itime = itime2;

    /* alloc and read pressure */
    /* requires that SP/LNSP, hyai and hybi are in the ECMWF netCDF file */ /* in ancillary.c */
    status =
      alloc_and_read_ECMWF_netCDF_pressure (ncid, &(p_level[t]), &(tmp_nlev), &(p_layer[t]), &(nlay), itime, ilat, ilon, verbose);
    CHKERR (status);

    p_layer[t] = realloc (p_layer[t], (nlay + 1) * sizeof (float));
    CHKPOINTER (p_layer[t]);

    /* read surface pressure from ECMWF file */
    status = read_ECMWF_surface_pressure (ncid, itime, ilat, ilon, &(SP), verbose);
    CHKERR (status);

    p_layer[t][nlay] = SP / 100.0;

    status = nc_inq_varid (ncid, "t", &id_var_test);
    if (status == NC_NOERR) {
      strcpy (T_name, "t");
      strcpy (U_name, "u");
      strcpy (V_name, "v");
      strcpy (W_name, "w");
    } else {
      status = nc_inq_varid (ncid, "T", &id_var_test);
      if (status == NC_NOERR) {
        strcpy (T_name, "T");
        strcpy (U_name, "U");
        strcpy (V_name, "V");
        strcpy (W_name, "W");
      } else {
        fprintf (stderr, "Error, unknown format of the ECMWF wind file %s\n", filename);
        CHKERR (status);
      }
    }

    /* allocate and read temperature, defined on layers */
    status = alloc_and_read_netCDF_column_float (ncid, T_name, &T_layer[t], nlay, itime, ilat, ilon, verbose);
    CHKERR (status);

    T_layer[t] = realloc (T_layer[t], (nlay + 1) * sizeof (float));
    CHKPOINTER (T_layer[t]);

    /* extrapolate temperature from last layer midpoint to surface */
    dx = (log (p_level[0][tmp_nlev - 1] / p_level[0][tmp_nlev - 1]) - log (p_layer[0][nlay - 1] / p_level[0][tmp_nlev - 1])) /
         (log (p_layer[0][nlay - 1] / p_level[0][tmp_nlev - 1]) - log (p_layer[0][nlay - 2] / p_level[0][tmp_nlev - 1]));
    T_layer[t][nlay] = (1.0 + dx) * T_layer[0][nlay - 1] - dx * T_layer[0][nlay - 2];

    /* allocate and read wind component u, defined on layers */
    status = alloc_and_read_netCDF_column_float (ncid, U_name, &u_layer[t], nlay, itime, ilat, ilon, verbose);
    CHKERR (status);

    u_layer[t] = realloc (u_layer[t], (nlay + 1) * sizeof (float));
    CHKPOINTER (u_layer[t]);

    u_layer[t][nlay] = 0.0;

    /* allocate and read wind component v, defined on layers */
    status = alloc_and_read_netCDF_column_float (ncid, V_name, &v_layer[t], nlay, itime, ilat, ilon, verbose);
    CHKERR (status);

    v_layer[t] = realloc (v_layer[t], (nlay + 1) * sizeof (float));
    CHKPOINTER (v_layer[t]);

    v_layer[t][nlay] = 0.0;

    /* allocate and read wind component w, defined on layers */
    status = alloc_and_read_netCDF_column_float (ncid, W_name, &w_layer[t], nlay, itime, ilat, ilon, verbose);
    CHKERR (status);

    w_layer[t] = realloc (w_layer[t], (nlay + 1) * sizeof (float));
    CHKPOINTER (w_layer[t]);

    w_layer[t][nlay] = 0.0;
  }

  if (verbose) {
    fprintf (stderr, "     found %zd x %zd x %zd (lat x lon x lev) data points\n", nlat, nlon, (*nlev));
    fprintf (stderr, "     reading pixel lat = %5.2f (%4ld), lon = %5.2f (%4ld)\n", lat, ilat, lon, ilon);
  }

  /* number of final levels == number of layers plus one, which was additional realloced above */
  (*nlev) = nlay + 1;

  if (nt == 1) {
    /* no time interpolation needed, do nothing */
  } else {
    /* write time interpolated data into the zero'th entry */
    for (lc = 0; lc < (*nlev); lc++) {
      p_layer[0][lc] = (1.0 - dt) * p_layer[0][lc] + dt * p_layer[1][lc];
      T_layer[0][lc] = (1.0 - dt) * T_layer[0][lc] + dt * T_layer[1][lc];
      u_layer[0][lc] = (1.0 - dt) * u_layer[0][lc] + dt * u_layer[1][lc];
      v_layer[0][lc] = (1.0 - dt) * v_layer[0][lc] + dt * v_layer[1][lc];
      w_layer[0][lc] = (1.0 - dt) * w_layer[0][lc] + dt * w_layer[1][lc];
    }
  }

  /* calculate z from T and p using hydrostatic equation (in ancillary.c) */
  status = calculate_z_from_p_and_T (p_layer[0], T_layer[0], z_wind, (*nlev), altitude, verbose);
  CHKERR (status);

  /* allocate space for final results */
  *u = calloc ((*nlev), sizeof (float));
  CHKPOINTER (*u);
  /* allocate space for final results */
  *v = calloc ((*nlev), sizeof (float));
  CHKPOINTER (*v);
  /* allocate space for final results */
  *w = calloc ((*nlev), sizeof (float));
  CHKPOINTER (*w);

  /* copy data */
  for (lc = 0; lc < (*nlev); lc++) {
    (*u)[lc] = u_layer[0][lc];
    (*v)[lc] = v_layer[0][lc];
    (*w)[lc] = w_layer[0][lc];
    /* /\* additional verbose output at layer mid-points *\/ */
    /* fprintf (stderr, " lc =%3d, z=%7.3f, u=%7.3f, v=%7.3f, w=%7.3f ( p_layer=%9.3f, T_layer=%7.3f ) \n", */
    /*                   lc, (*z_wind)[lc], (*u)[lc], (*v)[lc], (*w)[lc], p_layer[0][lc], T_layer[0][lc]); */
  }

  status = ASCII_free_float (p_layer, nt);
  CHKERR (status);
  status = ASCII_free_float (T_layer, nt);
  CHKERR (status);
  status = ASCII_free_float (u_layer, nt);
  CHKERR (status);
  status = ASCII_free_float (v_layer, nt);
  CHKERR (status);
  status = ASCII_free_float (w_layer, nt);
  CHKERR (status);

#else
  fprintf (stderr, " ******************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot       *\n");
  fprintf (stderr, " * use the ECMWF data file option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ******************************************************************\n");
  CHKERR (-1);
#endif

  return status;
}

/*****************************************************************/
/* get_time_index                                                */
/*                                                               */
/* Purpose:                                                      */
/* read time grid from netCDF file and search                    */
/* the suitable time index to input.UTC                          */
/* if only one time is specified in the netCDF file,             */
/* take that field                                               */
/*                                                               */
/*  input:                                                       */
/*  ------                                                       */
/*  ncid       id number of the netCDF file                      */
/*  UTC        universal time correlated                         */
/*                                                               */
/*  output:                                                      */
/*  -------                                                      */
/*  nt         number of time points for futher calculations     */
/*  t1         time index 1                                      */
/*  t2         time index 2                                      */
/*  dt         normalized fraction between time1 and time2       */
/*                                                               */
/*  February 2007  by Ulrich Hamann                              */
/*                                                               */
/*  changes:                                                     */
/*  - July 2007 by Ulrich Hamann                                 */
/*    interpretation of the attribute time:units                 */
/*****************************************************************/

int get_time_index (int ncid, struct tm UTC, int time_interpolate, int* nt, int* t1, int* t2, float* dt, int verbose, int quiet) {

  DEPRECATE (2023);

#if HAVE_NETCDF4

#define FORMAT_UNKNOWN 0
#define YYYYMMDD_FF 1 /*  day as %Y%m%d.%f                 */
#define MONTH_OF_YEAR 2
#define HOURS_SINCE 3 /*  hours since 1900-01-01 00:00:0.0 */
#define SECONDS_SINCE 4

  int status = 0;

  time_t time_in_s   = NOT_DEFINED_INTEGER;
  time_t min_delta_t = 999999999;

  int        id_time        = NOT_DEFINED_INTEGER;
  int        time_format    = NOT_DEFINED_INTEGER;
  double*    time_grid      = NULL;
  time_t*    time_grid_in_s = NULL;
  struct tm* date_grid      = NULL;

  char tmp_string[FILENAME_MAX] = "";
  int  i                        = 0;
  int  year = 0, month = 0, day = 0, hour = 0, min = 0;

  size_t ntime = NOT_DEFINED_INTEGER;
  int    t     = NOT_DEFINED_INTEGER;
  int    t_min = NOT_DEFINED_INTEGER;

  double tmp_time = NOT_DEFINED_FLOAT;
  char*  timestr  = NULL;

  char data_unit[50] = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";

  if (UTC.tm_mday < 0) {
    if (!quiet) {
      fprintf (stderr, " ******* WARNING >>>>>> no time specified\n");
      fprintf (stderr, "     take first entry in the netCDF file\n");
    }
    (*t1) = 0; /* read first entry in the netCDF file */
    (*t2) = 0;
    (*dt) = 0.0;
    (*nt) = 1; /* read one entry in the netCDF file */
    CHKERR (status);
  }

  /* read time in ECMWF file */
  status = alloc_and_read_netCDF_1D_double (ncid, "time", &(ntime), "time", &(time_grid));
  CHKERR (status);

  /* alloc ECMWF time in s - array */
  if ((time_grid_in_s = calloc (ntime, sizeof (time_t))) == NULL) {
    fprintf (stderr, "Error, allocating memory for time_grid_in_s (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__);
    return -1;
  }

  /* get id for time */
  status = nc_inq_varid (ncid, "time", &id_time);

  if (status == NC_NOERR) { /* time is specified in the netCDF file */

    /* try to read attribute 'unit' and interprete it as much as possible */
    status = nc_get_att_text (ncid, id_time, "units", data_unit);

    if (status != NC_NOERR) {
      /* no given units -> take first entry in netCDF file */

      /*assume that it is day as %Y%m%d.%f */
      if (!quiet)
        fprintf (stderr,
                 "*** Warning, no time unit specified in the netCDF file, assume that it is \'day as \%%Y\%%m\%%d.\%%f\'\n");

      time_format = YYYYMMDD_FF;

      /* we hope this is OK, so set status to OK again */
      status = 0;
    } else {

      /* some information is given about the unit, try to interprete it */
      if (strncasecmp ("day as %Y%m%d.%f", data_unit, 15) == 0)
        time_format = YYYYMMDD_FF;
      else if (strncasecmp ("hours since", data_unit, 10) == 0)
        time_format = HOURS_SINCE;
      else if (strncasecmp ("seconds since", data_unit, 12) == 0)
        time_format = SECONDS_SINCE;
      else if (strncasecmp ("month of year", data_unit, 12) == 0)
        time_format = MONTH_OF_YEAR;
      else {
        fprintf (stderr,
                 " ... Error, unknown time format %s (line %d, function %s in %s)\n",
                 data_unit,
                 __LINE__,
                 __func__,
                 __FILE__);
        return -1;
      }
    }
  } else {
    /* no time specified in netCDF file -> take first entry in netCDF file */
    status      = 0;
    time_format = FORMAT_UNKNOWN;
  }

  /* alloc date_grid ( structure - ) array */
  if ((date_grid = calloc (ntime, sizeof (struct tm))) == NULL) {
    fprintf (stderr, "Error, allocating memory for date_grid (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -10;
  }

  switch (time_format) {
  case FORMAT_UNKNOWN:
    /* do nothing, than automatically the first entry is specified (index=0) */
    break;

  case YYYYMMDD_FF:
  case MONTH_OF_YEAR:

    /* if 'month of year', then shift the month by 2 digits to the left, and add 15 for the middle of the month */
    if (time_format == MONTH_OF_YEAR)
      for (t = 0; t < ntime; t++)
        time_grid[t] = time_grid[t] * 100.0 + 15.;

    if (ntime > 1) {

      /* for year 0 AD we assume climatologic data, add current year for time selection */
      if (0101.0 <= time_grid[0] && time_grid[ntime - 1] < 1232.0)
        for (t = 0; t < ntime; t++)
          time_grid[t] += (UTC.tm_year + 1900) * 10000.0;

      /* time range check */
      /* dates after the year 1600 AD and before 3000 AD */
      if (time_grid[0] < 16000101.00 || 30000101.00 < time_grid[0]) {
        fprintf (stderr,
                 "Error, time [%f,%f] in netCDF file outside assumed range [1600AD, 3000AD] \n",
                 time_grid[0],
                 time_grid[ntime - 1]);
        return -1;
      }

      /* convert format "YYYYMMDD.FF" into "C-standard format" */
      for (t = 0; t < ntime; t++) {
        date_grid[t].tm_year = floor (time_grid[t]) / 10000 - 1900;
        tmp_time             = time_grid[t] - (date_grid[t].tm_year + 1900.0) * 10000.0;
        date_grid[t].tm_mon  = floor (tmp_time) / 100 - 1;
        tmp_time             = tmp_time - (date_grid[t].tm_mon + 1.0) * 100.0;
        date_grid[t].tm_mday = floor (tmp_time);
        tmp_time             = tmp_time - date_grid[t].tm_mday;
        date_grid[t].tm_hour = floor (tmp_time * 24.0);
        tmp_time             = tmp_time - date_grid[t].tm_hour / 24.0;
        date_grid[t].tm_min  = floor (tmp_time * 24.0 * 60.0);
        tmp_time             = tmp_time - date_grid[t].tm_min / (24.0 * 60.0);
        date_grid[t].tm_sec  = floor (tmp_time * (24.0 * 60.0 * 60.0));

        date_grid[t].tm_wday = weekday (date_grid[t].tm_year + 1900, date_grid[t].tm_mon + 1, date_grid[t].tm_mday);

        time_grid_in_s[t] = my_timegm (&(date_grid[t]));

        /* timestr = asctime( &(date_grid[t]) ); */
        /* fprintf (stderr,"time_grid[%d] = %10.5f, %ld , %s", t,  time_grid[t], time_grid_in_s[t], timestr); */
      }
    }
    /* else, do nothing, than automatically the first entry is specified (index=0) */

    break;

  case HOURS_SINCE:
  case SECONDS_SINCE:

    if (time_format == HOURS_SINCE)
      i = 0;
    if (time_format == SECONDS_SINCE)
      i = 2;

    year  = atoi (substr (tmp_string, data_unit, 12 + i, 4));
    month = atoi (substr (tmp_string, data_unit, 17 + i, 2));
    day   = atoi (substr (tmp_string, data_unit, 20 + i, 2));
    hour  = atoi (substr (tmp_string, data_unit, 23 + i, 2));
    min   = atoi (substr (tmp_string, data_unit, 26 + i, 2));
    /* there are not always seconds in the format */

    /* fprintf (stderr, " time scince: %d %d %d %d %d \n", year, month, day, hour, min ); */

    /*  /\* boundary check *\/ */
    /*  if ( time_grid[0] < 0.0 || 500.0 < time_grid[0] ) {  /\* assuming more than 0 and less than 500 forecast hours *\/ */
    /*    fprintf (stderr, "Error -1, wrong format for time %f in the netCDF file \n", time_grid[0] ); */
    /*    return -1;  */
    /*  } */

    /* convert format "hours since 1900-01-01 00:00:0.0" into "C-standard format" */
    for (t = 0; t < ntime; t++) {
      date_grid[t].tm_year = year - 1900;
      date_grid[t].tm_mon  = month - 1;
      date_grid[t].tm_mday = day;
      if (time_format == HOURS_SINCE)
        date_grid[t].tm_hour = hour + time_grid[t];
      else
        date_grid[t].tm_hour = hour;
      date_grid[t].tm_min = min;
      if (time_format == SECONDS_SINCE)
        date_grid[t].tm_sec = 0.0 + time_grid[t];
      else
        date_grid[t].tm_sec = 0.0;

      date_grid[t].tm_wday = weekday (date_grid[t].tm_year + 1900, date_grid[t].tm_mon + 1, date_grid[t].tm_mday);

      time_grid_in_s[t] = my_timegm (&(date_grid[t]));

      /* timestr = asctime (&(date_grid[t])); */
      /* fprintf (stderr,"time_grid[%d] = %10.5f, %ld , %s", t,  time_grid[t], time_grid_in_s[t], timestr); */
    }
    break;

  default:
    fprintf (stderr, "Error, unknown time_format = %d (line %d, function %s in %s)\n", time_format, __LINE__, __func__, __FILE__);
    return -1;
  }

  /* convert time to search for into seconds scince 01.01.1970 00:00:00 UTC */
  time_in_s = my_timegm (&(UTC));

  /* search time index */
  switch (time_interpolate) {
  case TIME_NEAREST_DATE:
    min_delta_t = labs (time_grid_in_s[0] - time_in_s); /* first guess for time difference */
    t_min       = 0;                                    /* first guess for index */
    for (t = 1; t < ntime; t++)
      if (labs (time_grid_in_s[t] - time_in_s) < min_delta_t) {
        min_delta_t = labs (time_grid_in_s[t] - time_in_s);
        t_min       = t;
      }
    (*t1) = t_min;
    (*t2) = NOT_DEFINED_INTEGER;
    (*dt) = (float) t_min;
    (*nt) = 1; /* read first and only time step from file */

    if (verbose) {
      timestr = asctime (&(UTC));
      fprintf (stderr, "     specified time:              %s", timestr);
      timestr = asctime (&(date_grid[(*t1)]));
      fprintf (stderr, "     nearest time in netCDF file: %s", timestr);
    }

    break;
  case TIME_INTERPOLATION:
    if (ntime == 1) {
      if ((time_grid_in_s[0] != time_in_s) && !quiet) {
        timestr = asctime (&(UTC));
        fprintf (stderr, "\n ******* WARNING >>>>>> specified time %s", timestr);
        timestr = asctime (&(date_grid[0]));
        fprintf (stderr, " ******* WARNING >>>>>> is NOT the time found in the ECMWF atmosphere data file %s\n", timestr);
      }
      (*t1) = 0;
      (*t2) = NOT_DEFINED_INTEGER;
      (*dt) = 0.0;
      (*nt) = 1; /* read first and only time step from file */
    } else {
      for (t = 1; t < ntime; t++) {
        if (time_grid_in_s[t - 1] <= time_in_s && time_in_s <= time_grid_in_s[t]) {
          /* fprintf (stderr,"E[t-1] = %ld, t_s = %ld, E[t] = %ld\n", time_grid_in_s[t-1], time_in_s, time_grid_in_s[t]); */
          (*t1) = t - 1;
          (*t2) = t;
          (*dt) = (float) (time_in_s - time_grid_in_s[t - 1]) / (float) (time_grid_in_s[t] - time_grid_in_s[t - 1]);
          (*nt) = 2; /* read two time steps from file */
          /* fprintf (stderr,"t1 = %3d, t2 = %3d , dt = %10.5f\n", t1, t2, dt); */
          break;
        }
      }
      if (*t1 < 0) {
        fprintf (stderr, "Error -1 finding index for time = %s in get_time_index (ancillary.c)\n", asctime (&(UTC)));
        CHKERR (-1);
      }

      if (verbose) {
        fprintf (stderr, "     time interpolation of data between \n");
        timestr = asctime (&(date_grid[(*t1)]));
        fprintf (stderr, "     %s", timestr);
        timestr = asctime (&(date_grid[(*t2)]));
        fprintf (stderr, "     %s", timestr);
      }
    }
    break;
  default:
    fprintf (stderr, "Error, unknown time_interpolate = %d in get_time_index (in ancillary.c)\n", time_interpolate);
    return -1;
  }

  free (date_grid);
  free (time_grid);
  free (time_grid_in_s);
  return status;

#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
  fprintf (stderr, " ***********************************************************************\n");
  return -1;
#endif
}

/*****************************************************************/
/* get_local_apparent_time_index                                 */
/* (almost the same as get_time_index)                           */
/*                                                               */
/* Purpose:                                                      */
/* read time grid from netCDF file and search                    */
/* the suitable time index to input.LAT                          */
/* if only one time is specified in the netCDF file,             */
/* take that field                                               */
/*                                                               */
/*  input:                                                       */
/*  ------                                                       */
/*  ncid       id number of the netCDF file                      */
/*  LAT        local apparent time                               */
/*                                                               */
/*  output:                                                      */
/*  -------                                                      */
/*  nt         number of time points for futher calculations     */
/*  t1         time index 1                                      */
/*  t2         time index 2                                      */
/*  dt         normalized fraction between time1 and time2       */
/*                                                               */
/*  April 2007  by Ulrich Hamann                                 */
/*****************************************************************/

int get_local_apparent_time_index (int       ncid,
                                   struct tm LAT,
                                   int       time_interpolate,
                                   int*      nt,
                                   int*      t1,
                                   int*      t2,
                                   float*    dt,
                                   int       verbose,
                                   int       quiet) {
  DEPRECATE (2023);
  int status = 0;

  time_t time_in_s   = NOT_DEFINED_INTEGER;
  time_t min_delta_t = 999999999;

  double*    time_grid      = NULL;
  time_t*    time_grid_in_s = NULL;
  struct tm* date_grid      = NULL;

  size_t ntime = NOT_DEFINED_INTEGER;
  int    t     = NOT_DEFINED_INTEGER;
  int    t_min = NOT_DEFINED_INTEGER;

  double tmp_time = NOT_DEFINED_FLOAT;
  char*  timestr  = NULL;

  if (LAT.tm_mday < 0) {
    fprintf (stderr, " ******* WARNING >>>>>> no time specified\n");
    fprintf (stderr, "     take first entry in the netCDF file\n");
    (*t1) = 0; /* read first entry in the netCDF file */
    (*t2) = 0;
    (*dt) = 0.0;
    (*nt) = 1; /* read one entry in the netCDF file */
    CHKERR (status);
  }

  /* time in seconds scince 01.01.1970 00:00:00 LAT */
  time_in_s = my_timegm (&(LAT));

  /* read time in netCDF file */
  status = alloc_and_read_netCDF_1D_double (ncid, "local_apparent_time", &(ntime), "local_apparent_time", &(time_grid));
  CHKERR (status);

  /* alloc netCDF date structure */
  date_grid = calloc (ntime, sizeof (struct tm));
  CHKPOINTER (date_grid);

  /* alloc netCDF time in s - array */
  time_grid_in_s = calloc (ntime, sizeof (time_t));
  CHKPOINTER (time_grid_in_s);

  /* convert format "YYYYMMDD.day_fraction" into "C-standard format" */
  for (t = 0; t < ntime; t++) {
    date_grid[t].tm_year = floor (time_grid[t]) / 10000 - 1900;
    tmp_time             = time_grid[t] - (date_grid[t].tm_year + 1900.0) * 10000.0;
    date_grid[t].tm_mon  = floor (tmp_time) / 100 - 1;
    tmp_time             = tmp_time - (date_grid[t].tm_mon + 1.0) * 100.0;
    date_grid[t].tm_mday = floor (tmp_time);
    tmp_time             = tmp_time - date_grid[t].tm_mday;
    date_grid[t].tm_hour = floor (tmp_time * 24.0);
    tmp_time             = tmp_time - date_grid[t].tm_hour / 24.0;
    date_grid[t].tm_min  = floor (tmp_time * 24.0 * 60.0);
    tmp_time             = tmp_time - date_grid[t].tm_min / (24.0 * 60.0);
    date_grid[t].tm_min  = floor (tmp_time * (24.0 * 60.0 * 60.0));

    date_grid[t].tm_wday = weekday (date_grid[t].tm_year + 1900, date_grid[t].tm_mon + 1, date_grid[t].tm_mday);

    time_grid_in_s[t] = my_timegm (&(date_grid[t]));

    /* timestr = asctime (&(date_grid[t])); */
    /* fprintf (stderr,"time_grid[%d] = %10.5f, %ld , %s", t,  time_grid[t], time_grid_in_s[t], timestr); */
  }

  /* search time index */

  switch (time_interpolate) {
  case TIME_NEAREST_DATE:
    min_delta_t = labs (time_grid_in_s[0] - time_in_s); /* first guess for time difference */
    t_min       = 0;                                    /* first guess for index */
    for (t = 1; t < ntime; t++)
      if (labs (time_grid_in_s[t] - time_in_s) < min_delta_t) {
        min_delta_t = labs (time_grid_in_s[t] - time_in_s);
        t_min       = t;
      }
    (*t1) = t_min;
    (*t2) = NOT_DEFINED_INTEGER;
    (*dt) = (float) t_min;
    (*nt) = 1; /* read first and only time step from file */

    if (verbose) {
      timestr = asctime (&(LAT));
      fprintf (stderr, "     specified time:              %s", timestr);
      timestr = asctime (&(date_grid[(*t1)]));
      fprintf (stderr, "     nearest time in netCDF file: %s", timestr);
    }

    break;
  case TIME_INTERPOLATION:
    if (ntime == 1) {
      if ((time_grid_in_s[0] != time_in_s) && !quiet) {
        timestr = asctime (&(LAT));
        fprintf (stderr, "\n ******* WARNING >>>>>> specified time %s", timestr);
        timestr = asctime (&(date_grid[0]));
        fprintf (stderr, " ******* WARNING >>>>>> is NOT the time found in the netCDF atmosphere data file %s\n", timestr);
      }
      (*t1) = 0;
      (*t2) = NOT_DEFINED_INTEGER;
      (*dt) = 0.0;
      (*nt) = 1; /* read first and only time step from file */
    } else {
      for (t = 1; t < ntime; t++) {
        if (time_grid_in_s[t - 1] <= time_in_s && time_in_s <= time_grid_in_s[t]) {
          /* fprintf (stderr,"E[t-1] = %ld, t_s = %ld, E[t] = %ld\n", time_grid_in_s[t-1], time_in_s, time_grid_in_s[t]); */
          (*t1) = t - 1;
          (*t2) = t;
          (*dt) = (float) (time_in_s - time_grid_in_s[t - 1]) / (float) (time_grid_in_s[t] - time_grid_in_s[t - 1]);
          (*nt) = 2; /* read two time steps from file */
          /* fprintf (stderr,"t1 = %3d, t2 = %3d , dt = %10.5f\n", t1, t2, dt); */
          break;
        }
      }
      if (*t1 < 0) {
        fprintf (stderr, "Error -1 finding index for time = %s\n", asctime (&(LAT)));
        CHKERR (-1);
      }

      if (verbose) {
        fprintf (stderr, "     time interpolation of data between \n");
        timestr = asctime (&(date_grid[(*t1)]));
        fprintf (stderr, "     %s", timestr);
        timestr = asctime (&(date_grid[(*t2)]));
        fprintf (stderr, "     %s", timestr);
      }
    }
    break;
  default:
    fprintf (stderr, "Error, unknown time_interpolate = %d\n", time_interpolate);
    CHKERR (-1);
  }

  return status;
}

/******************************************************************************/
/* alloc_and_read_ECMWF_netCDF_pressure                                       */
/* small function to read hybrid coefficients hyai, hybi and surface pressure */
/* from ECMWF netCDF file and convert those to a pressure array               */
/* requires that SP, hyai and hybi are in the ECMWF netCDF file               */
/******************************************************************************/

int alloc_and_read_ECMWF_netCDF_pressure (int     ncid,
                                          float** p_level,
                                          size_t* nlev,
                                          float** p_layer,
                                          size_t* nlay,
                                          size_t  itime,
                                          size_t  ilat,
                                          size_t  ilon,
                                          int     verbose) {
  DEPRECATE (2023);

#if HAVE_NETCDF4

  int status = 0;

  float SP = -999.0;

  double* hyai = NULL;
  double* hybi = NULL;
  double* hyam = NULL;
  double* hybm = NULL;
  int     lc   = 0;

  /* read surface pressure from ECMWF file */
  status = read_ECMWF_surface_pressure (ncid, itime, ilat, ilon, &(SP), verbose);
  if (status != NC_NOERR) {
    fprintf (stderr,
             "Error '%s' returned by read_ECMWF_surface_pressure (line %d, function '%s' in '%s')\n",
             nc_strerror (status),
             __LINE__,
             __func__,
             __FILE__);
    CHKERR (status);
  }

  /* read hybrid pressure coefficients on levels (layer boundaries) */
  status = alloc_and_read_netCDF_1D_double (ncid, "ilev", nlev, "hyai", &(hyai));
  CHKERR (status);
  status = alloc_and_read_netCDF_1D_double (ncid, "ilev", nlev, "hybi", &(hybi));
  CHKERR (status);

  /* read hybrid pressure coefficients at layers (layer midpoints) */
  status = alloc_and_read_netCDF_1D_double (ncid, "mlev", nlay, "hyam", &(hyam));
  CHKERR (status);
  status = alloc_and_read_netCDF_1D_double (ncid, "mlev", nlay, "hybm", &(hybm));
  CHKERR (status);

  /* allocate pressure */
  *p_level = calloc ((*nlev), sizeof (float));
  CHKPOINTEROUT (*p_level, "failed allocating memory for 'p_level'");

  /*   /\***************************************************************\/ */
  /*   /\* calculate level pressure (layer boundaries)                 *\/ */
  /*   /\* starting from 1, because we DO NOT like the uppermost p=0.0 *\/ */
  /*   /\* it is bad for merging with background atmosphere            *\/ */
  /*   /\***************************************************************\/ */

  /*   /\* allocate pressure *\/ */
  /*   if (((*p_level) = calloc((*nlev)-1, sizeof(float)))==NULL) { */
  /*     fprintf (stderr, "Error allocating memory for 'p_level' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__); */
  /*     return -1; */
  /*   }   */

  /*   /\* fprintf (stderr,"level pressures\n"); *\/ */
  /*   for (lc=1; lc<(*nlev); lc++) { */
  /*     (*p_level)[lc-1]= 0.01*(hyai[lc]+hybi[lc]*SP); /\* 0.01 == Pa -> hPa *\/ */
  /*     /\* fprintf (stderr,"lc=%3d a=%10.2f, b=%10.6f, p=%10.4f, SP=%10.4f, SP-p=%10.4f\n", *\/ */
  /*     /\*                  lc-1,hyai[lc],hybi[lc],(*p_level)[lc-1],SP/100.,SP/100.-(*p_level)[lc-1]); *\/ */
  /*   } */

  /*   (*nlev)=(*nlev)-1; */

  /* fprintf (stderr,"level pressures\n"); */
  for (lc = 0; lc < (*nlev); lc++) {
    (*p_level)[lc] = 0.01 * (hyai[lc] + hybi[lc] * SP); /* 0.01 == Pa -> hPa */
    /* fprintf (stderr,"lc=%3d a=%10.2f, b=%10.6f, p=%10.4f, SP=%10.4f, SP-p=%10.4f\n", */
    /*                 lc-1,hyai[lc],hybi[lc],(*p_level)[lc],SP/100.,SP/100.-(*p_level)[lc]); */
  }

  /* allocate pressure */
  if (((*p_layer) = calloc ((*nlay), sizeof (float))) == NULL) {
    fprintf (stderr, "Error allocating memory for 'p_layer' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -1;
  }

  /* calculate layer pressure (layer midpoints) */
  /* fprintf (stderr,"layer pressures\n"); */
  for (lc = 0; lc < (*nlay); lc++) {
    (*p_layer)[lc] = 0.01 * (hyam[lc] + hybm[lc] * SP); /* 0.01 == Pa -> hPa */
    /* fprintf (stderr,"lc=%3d a=%10.2f, b=%10.6f, p=%10.4f, SP=%10.4f, SP-p=%10.4f\n", */
    /*                  lc,hyam[lc],hybm[lc],(*p_layer)[lc],SP/100.,SP/100.-(*p_layer)[lc]); */
  }

  return status;

#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use any ECMWF input option. Please get netcdf and rebuild.          *\n");
  fprintf (stderr, " ***********************************************************************\n");
  return -1;
#endif
}

int read_ECMWF_surface_pressure (int ncid, size_t itime, size_t ilat, size_t ilon, float* SP, int verbose) {

#if HAVE_NETCDF4

  int status = 0;

  int id_data = 0;
  int log_SP  = FALSE;

  nc_type netCDF_type;
  float   scale_factor                = 1.0;
  float   add_offset                  = 0.0;
  char    variable_name[FILENAME_MAX] = "";

  size_t* index = NULL;

  /* get variable id for "SP" (surface pressure) */
  if ((status = nc_inq_varid (ncid, "SP", &id_data)) == NC_NOERR)
    strcpy (variable_name, "SP");
  else if ((status = nc_inq_varid (ncid, "sp", &id_data)) == NC_NOERR)
    strcpy (variable_name, "sp");
  else if ((status = nc_inq_varid (ncid, "LNSP", &id_data)) == NC_NOERR) {
    strcpy (variable_name, "LNSP");
    log_SP = TRUE;
  } else if ((status = nc_inq_varid (ncid, "lnsp", &id_data)) == NC_NOERR) {
    strcpy (variable_name, "lnsp");
    log_SP = TRUE;
  } else {
    fprintf (stderr, "Error '%s', while getting id for surface pressure\n", nc_strerror (status));
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    CHKERR (status);
  }

  /* get type of the variable */
  status = nc_inq_vartype (ncid, id_data, &(netCDF_type));
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s', while getting type (float, double ...) of '%s'\n", nc_strerror (status), variable_name);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    CHKERR (status);
  }

  /* read attribute scale_factor */
  status = nc_get_att_float (ncid, id_data, "scale_factor", &scale_factor);
  if (status != NC_NOERR) {
    /* no given scale_factor -> default scale_factor = 1.0 */
    /* that's OK, put status to 0 */
    status = 0;
  } else {
    /* if (verbose) */
    /*  fprintf (stderr,"     scale_factor = %f \n", scale_factor); */
  }
  CHKERR (status);

  /* read attribute add_offset */
  status = nc_get_att_float (ncid, id_data, "add_offset", &add_offset);
  if (status != NC_NOERR) {
    /* no given offset -> default offset = 0.0 */
    /* that's OK, put status to 0 */
    status = 0;
  } else {
    /* if (verbose) */
    /*   fprintf (stderr,"     add_offset = %f \n", add_offset); */
  }
  CHKERR (status);

  if (log_SP == FALSE) {

    if (((index) = calloc (3, sizeof (size_t))) == NULL) {
      fprintf (stderr, "Error, allocting memory for index\n");
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
      return -1;
    }

    /* index for surface pressure */
    index[0] = itime;
    index[1] = ilat;
    index[2] = ilon;

  } else { /* that means -> if (log_SP == TRUE) */

    if (((index) = calloc (4, sizeof (size_t))) == NULL) {
      fprintf (stderr, "Error, allocting memory for index\n");
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
      return -1;
    }

    /* index for log surface pressure */
    index[0] = itime;
    index[1] = 0; /* level == 0, there is only one level! */
    index[2] = ilat;
    index[3] = ilon;
  }

  /* read surface pressure */
  status = nc_get_var1_float (ncid, id_data, index, SP);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error '%s', while reading surface pressure\n", nc_strerror (status));
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    CHKERR (status);
  }

  (*SP) = (*SP) * scale_factor + add_offset;

  /* if (verbose) */
  /*   fprintf (stderr,"     %s = %f \n", variable_name, (*SP)); */

  if (log_SP == TRUE)
    (*SP) = exp ((*SP));

  free (index);

  return status;

#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use any ECMWF input option. Please get netcdf and rebuild.          *\n");
  fprintf (stderr, " ***********************************************************************\n");
  return -1;
#endif
}

/**********************************************************************/
/* alloc_and_read_netCDF_column_float                                 */
/* function extracts a column of size nlev of data from a netCDF file */
/* at the indeces itime, ilat, ilon                                   */
/* variable is the name (string), which is stored in the netCDF file  */
/**********************************************************************/

int alloc_and_read_netCDF_column_float (int     ncid,
                                        char*   variable_name,
                                        float** data,
                                        size_t  nlev,
                                        size_t  itime,
                                        size_t  ilat,
                                        size_t  ilon,
                                        int     verbose) {

#if HAVE_NETCDF4

  int     status  = 0;
  int     id_data = 0;
  nc_type netCDF_type;
  double  scale_factor = 1.0;
  double  add_offset   = 0.0;

  unsigned char data_uchar  = 0;
  short         data_short  = -999;
  int           data_int    = -999;
  float         data_float  = -999.0;
  double        data_double = -999.0;

  size_t index4D[4] = {0, 0, 0, 0}; /* time, lev, lat, lon */
  int    lc         = 0;

  /* get variable id for data */
  status = nc_inq_varid (ncid, variable_name, &id_data);
  if (status != NC_NOERR) {
    fprintf (stderr,
             "Error '%s', while getting id for '%s' (line %d, function '%s' in '%s')\n",
             nc_strerror (status),
             variable_name,
             __LINE__,
             __func__,
             __FILE__);
    CHKERR (status);
  }

  /* get type of the variable */
  status = nc_inq_vartype (ncid, id_data, &(netCDF_type));
  if (status != NC_NOERR) {
    fprintf (stderr,
             "Error '%s', while getting type for %s (line %d, function '%s' in '%s')\n",
             nc_strerror (status),
             variable_name,
             __LINE__,
             __func__,
             __FILE__);
    CHKERR (status);
  }

  if (verbose)
    fprintf (stderr,
             "     reading '%s' type %d (NC_BYTE %d, NC_CHAR %d, NC_SHORT %d, NC_INT %d, NC_FLOAT %d, NC_DOUBLE %d)\n",
             variable_name,
             netCDF_type,
             NC_BYTE,
             NC_CHAR,
             NC_SHORT,
             NC_INT,
             NC_FLOAT,
             NC_DOUBLE);

  /* read attribute scale_factor */
  status = nc_get_att_double (ncid, id_data, "scale_factor", &scale_factor);
  if (status != NC_NOERR) {
    /* no given scale_factor -> default scale_factor = 1.0 */
    /* that's OK, put status to 0 */
    status = 0;
  } else {
    /* if (verbose) */
    /*   fprintf (stderr,"     scale_factor = %f \n", scale_factor); */
  }
  CHKERR (status);

  /* read attribute add_offset */
  status = nc_get_att_double (ncid, id_data, "add_offset", &add_offset);
  if (status != NC_NOERR) {
    /* no given offset -> default offset = 0.0 */
    /* that's OK, put status to 0 */
    status = 0;
  } else {
    /* if (verbose) */
    /*   fprintf (stderr,"     add_offset = %f \n", add_offset); */
  }
  CHKERR (status);

  /* allocate temporary data array */
  *data = calloc (nlev, sizeof (float));
  CHKPOINTER (*data);

  /* read temperature (in netCDF = float *T(time, mlev, lat, lon);)*/
  index4D[0] = itime;
  /* index4D[1] see loops */
  index4D[2] = ilat;
  index4D[3] = ilon;

  switch (netCDF_type) {
  case (NC_BYTE):
  case (NC_CHAR):
    for (lc = 0; lc < nlev; lc++) {
      index4D[1]  = lc;
      status      = nc_get_var1_uchar (ncid, id_data, index4D, &(data_uchar));
      (*data)[lc] = (float) (((double) data_uchar) * scale_factor + add_offset);
      /* fprintf (stderr," lc = %3d, %s = %e \n", lc, variable_name, (*data)[lc]); */
    }
    if (status != NC_NOERR) {
      fprintf (stderr,
               "Error '%s', reading %s (uchar) \n  (line %d, function '%s' in '%s')\n",
               nc_strerror (status),
               variable_name,
               __LINE__,
               __func__,
               __FILE__);
      CHKERR (status);
    }
    break;
  case (NC_SHORT):
    for (lc = 0; lc < nlev; lc++) {
      index4D[1]  = lc;
      status      = nc_get_var1_short (ncid, id_data, index4D, &(data_short));
      (*data)[lc] = (float) (((double) data_short) * scale_factor + add_offset);
      if (fabs ((*data)[lc]) < scale_factor * 10E-6 && scale_factor != 1.0) /* THIS IS REALLY NOT NICE !!!! */
        (*data)[lc] = 0.0;                                                  /* THIS IS REALLY NOT NICE !!!! */
      /* fprintf (stderr," lc = %3d, %s = %e \n", lc, variable_name, (*data)[lc]); */
    }
    if (status != NC_NOERR) {
      fprintf (stderr,
               "Error '%s', reading %s (short) \n  (line %d, function '%s' in '%s')\n",
               nc_strerror (status),
               variable_name,
               __LINE__,
               __func__,
               __FILE__);
      CHKERR (status);
    }
    break;
  case (NC_INT):
    for (lc = 0; lc < nlev; lc++) {
      index4D[1]  = lc;
      status      = nc_get_var1_int (ncid, id_data, index4D, &(data_int));
      (*data)[lc] = (float) (((double) data_int) * scale_factor + add_offset);
      /* fprintf (stderr," lc = %3d, %s = %e \n", lc, variable_name, (*data)[lc]); */
    }
    if (status != NC_NOERR) {
      fprintf (stderr,
               "Error '%s', reading %s (integer) \n  (line %d, function '%s' in '%s')\n",
               nc_strerror (status),
               variable_name,
               __LINE__,
               __func__,
               __FILE__);
      CHKERR (status);
    }
    break;
  case (NC_FLOAT):
    for (lc = 0; lc < nlev; lc++) {
      index4D[1]  = lc;
      status      = nc_get_var1_float (ncid, id_data, index4D, &(data_float));
      (*data)[lc] = (float) (((double) data_float) * scale_factor + add_offset);
      /* if (status != NC_NOERR) fprintf (stderr,"Error '%s' reading '%s', index lc = %3d \n",
         nc_strerror(status), variable_name, lc); */
    }
    if (status != NC_NOERR) {
      fprintf (stderr,
               "Error '%s', reading %s (float) \n  (line %d, function '%s' in '%s')\n",
               nc_strerror (status),
               variable_name,
               __LINE__,
               __func__,
               __FILE__);
      CHKERR (status);
    }
    break;
  case (NC_DOUBLE):
    for (lc = 0; lc < nlev; lc++) {
      index4D[1]  = lc;
      status      = nc_get_var1_double (ncid, id_data, index4D, &(data_double));
      (*data)[lc] = (float) ((data_double) *scale_factor + add_offset);
      /* fprintf (stderr," lc = %3d, %s = %e \n", lc, variable_name, (*data)[lc]); */
    }
    if (status != NC_NOERR) {
      fprintf (stderr,
               "Error '%s', reading %s (double)\n  (line %d, function '%s' in '%s')\n",
               nc_strerror (status),
               variable_name,
               __LINE__,
               __func__,
               __FILE__);
      CHKERR (status);
    }
    break;
  default:
    fprintf (stderr,
             "Error, unknown type of variable %d \n  (line %d, function '%s' in '%s') \n",
             netCDF_type,
             __LINE__,
             __func__,
             __FILE__);
    return -1;
  }

  return status;

#else
  fprintf (stderr, " ******************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot       *\n");
  fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.          *\n");
  fprintf (stderr, " ******************************************************************\n");
  CHKERR (-1);
#endif
}

/***********************************************************************************/
/* Function: calculate_z_from_p_and_T                                              */
/* Description:                                                                    */
/*      calculates a z-grid using the hydrostatic equation                         */
/*      2 possible equation:                                                       */
/*      first:  z_{i-1} - z_{i} = R / g *0.5* bar{T} * log(p_{i}/p_{i-1}           */
/*      second: (longer derivation takes                                           */
/*                                                                                 */
/* Parameters:                                                                     */
/*      p         pressure in hPa (input)                                          */
/*      T         temperature in K (input)                                         */
/*      z         height above sea level in km (output)                            */
/*      n         number of height levels                                          */
/*      altitude  altitude of the ground in km (input)                             */
/*                                                                                 */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    ancillary.c                                                           */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    xx  200?   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int calculate_z_from_p_and_T (float* p, float* T, float** z, int n, float altitude, int verbose) {
  int status = 0;
  int lc     = NOT_DEFINED_INTEGER;
  //20120816ak g is not used, commented
  //  float   g  = NOT_DEFINED_FLOAT;

  /* Allocate z */
  *z = calloc (n, sizeof (float));
  CHKPOINTEROUT (*z, "failed allocating memory for z");

  if (verbose) {
    fprintf (stderr, " ... converting p and T to a z-grid\n");
    fprintf (stderr, "     assuming: R = %12.6f, g(z=0km) = %12.6f\n", R_AIR, G_SURFACE);
  }

  /* initialise lowest level */
  if (altitude == 0.0) {
    (*z)[n - 1] = 0.0;
    //20120816ak g is not used, commented
    //    g = G_SURFACE;
  } else {
    (*z)[n - 1] = altitude;
    //20120816ak g is not used, commented
    //    g = G_SURFACE * (R_EARTH * R_EARTH) / ((R_EARTH+(*z)[n-1]*1000.)*(R_EARTH+(*z)[n-1]*1000.));
  }

  /* calculate radiosonde z-grid (hydrostatic equation) */
  for (lc = n - 1; lc >= 1; lc--) {

    /* equation (1) for no variation of g inside one layer */
    /*  (*z)[lc-1]=(*z)[lc] + 0.001 * R/g*0.5*(T[lc]+T[lc-1])*log(p[lc]/p[lc-1]);     */
    /*                        /\* 0.001 <==> m -> km *\/                                          */
    /*  g = G_SURFACE * (R_EARTH * R_EARTH) / ((R_EARTH+1000.0*(*z)[lc-1])*(R_EARTH+1000.0*(*z)[lc-1])); */

    /* equation (2) taking care of variation of g inside one layer */
    (*z)[lc - 1] = -R_EARTH + 1 / (1 / (R_EARTH + 1000.0 * (*z)[lc]) -
                                   R_AIR / (G_SURFACE * R_EARTH * R_EARTH) * 0.5 * (T[lc] + T[lc - 1]) * log (p[lc] / p[lc - 1]));
    (*z)[lc - 1] *= 0.001; /* m -> km */

    /*     /\* addiotional verbose output *\/ */
    /*     fprintf (stderr,"     lc-1=%3d, dz(g=const)[km]=%7.4f, dz(g=g(z))[km]=%7.4f, z_sur[km]=%8.4f, T=%5.2f, p=%10.4f \n", */
    /*                      lc-1, 0.001 * R_AIR/g*0.5*(T[lc]+T[lc-1])*log(p[lc]/p[lc-1]), */
    /*                      (*z)[lc-1]-(*z)[lc],(*z)[lc-1],T[lc-1], p[lc-1]); */
  }

  return status;
}

/*  given month, day, year, returns day of week, eg. Monday = 0 etc. */
/*  tested for 1901 to 2099 (seems to work from 1800 on too)         */

int weekday (int year, int month, int day) {
  int ix = NOT_DEFINED_INTEGER, tx = NOT_DEFINED_INTEGER, vx = NOT_DEFINED_INTEGER;

  switch (month) {
  case 2:
  case 6:
    vx = 0;
    break;
  case 8:
    vx = 4;
    break;
  case 10:
    vx = 8;
    break;
  case 9:
  case 12:
    vx = 12;
    break;
  case 3:
  case 11:
    vx = 16;
    break;
  case 1:
  case 5:
    vx = 20;
    break;
  case 4:
  case 7:
    vx = 24;
    break;
  default:
    fprintf (stderr, "Error, determening day of week in function weekday (in ancillary.c) (month = %d)\n", month);
    return -1000;
    break;
  }

  if (year > 1900) /* 1900 was not a leap year */
    year -= 1900;
  ix = ((year - 21) % 28) + vx + (month > 2); /* take care of February  */
  tx = (ix + (ix / 4)) % 7 + day;             /* take care of leap year */
  return ((tx + 1) % 7);
}

char* strtrim (char* str, const char* trim) {
  return strltrim (strrtrim (str, trim), trim);
}

char* strrtrim (char* str, const char* trim) {
  char* end;

  if (!str)
    return NULL;

  if (!trim)
    trim = " \t\n\r";

  end = str + strlen (str);

  while (end-- > str) {
    if (!strchr (trim, *end))
      return str;
    *end = 0;
  }
  return str;
}

char* strltrim (char* str, const char* trim) {
  if (!str)
    return NULL;

  if (!trim)
    trim = " \t\r\n";

  while (*str) {
    if (!strchr (trim, *str))
      return str;
    ++str;
  }
  return str;
}

/***********************************************************************************/
/* Function: specific_heat_capacity_moist_air                                      */
/* Description:                                                                    */
/*  calculated the specific heat capacity of dry air and water vapour and          */
/*  returns a mass weighted average according to the specific humidity             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     float T_in_K       temperature in Kelvin                                    */
/*     float dens_air     number density of air moleculs in cm-3                   */
/*     float dens_wv      number density of water vapour moleculs in cm-3          */
/*     float mol_mass_air weight of one 'air' molecul in u                         */
/*     float mol_mass_wv  weight of one 'water vapour' molecul in u                */
/*     float *c_p         returned result                                          */
/*                          specific heat capacity in J/(kg K)                     */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    ancillary.c                                                           */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Feb 2008   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int specific_heat_capacity_moist_air (float T_in_K, float dens_air, float dens_wv, float* c_p, int quiet) {

  int status = 0;

  float mmr_h2o     = NOT_DEFINED_FLOAT;
  float c_p_dry_air = NOT_DEFINED_FLOAT;
  float c_p_wv      = NOT_DEFINED_FLOAT;

  mmr_h2o = (dens_wv * MOL_MASS_WV) / (dens_air * MOL_MASS_AIR);

  /* temperature dependent heat capacity of dry air */
  c_p_dry_air = specific_heat_capacity_dry_air (T_in_K, quiet);
  if (c_p_dry_air < 0.0) {
    fprintf (stderr, "Error during specific_heat_capacity_dry_air\n");
    CHKERR (-1);
  }

  /* temperature dependent heat capacity of water vapour */
  c_p_wv = specific_heat_capacity_water_vapour (T_in_K, quiet);
  if (c_p_wv < 0.0) {
    fprintf (stderr, "Error during specific_heat_capacity_water_vapour\n");
    CHKERR (-1);
  }

  /* mass weighted mean of specific heat capacity (moist air), e.g. Etling page 34 */
  *c_p = (1 - mmr_h2o) * c_p_dry_air + mmr_h2o * c_p_wv;

  return status;
}

/***********************************************************************************************/
/* function: specific_heat_capacity_dry_air                                                    */
/* interpolate a table to get the temperature dependent specific heat capacity of dry air      */
/* input     float T in K                                                                      */
/* output    float specific heat capacity of dry air in J/(kg K)                               */
/* ulrich hamann 2007-03-31                                                                    */
/***********************************************************************************************/

float specific_heat_capacity_dry_air (float T_in_K, int quiet) {

  const int n           = 14;
  int       i           = NOT_DEFINED_INTEGER;
  float     c_p_dry_air = NOT_DEFINED_FLOAT;

  float c_p_grid[14] =
    {1002.3, 1002.5, 1002.7, 1003.1, 1003.8, 1004.9, 1006.3, 1008.2, 1010.6, 1013.5, 1020.6, 1029.5, 1039.8, 1051.1};
  float T_grid[14] = {175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 450.0, 500.0, 550.0, 600.0};

  if (T_in_K < T_grid[0]) {
    if (!quiet) {
      fprintf (stderr, " *** Warning, getting specific heat capacity of dry air,\n");
      fprintf (stderr, "     temperature %7.2f is below of the tabled temperature region.\n", T_in_K);
    }
    return c_p_grid[0];
  } else if (T_in_K > T_grid[n - 1]) {
    if (!quiet) {
      fprintf (stderr, " *** Warning, getting specific heat capacity of dry air,\n");
      fprintf (stderr, "     temperature %7.2f is above of the tabled temperature region.\n", T_in_K);
    }
    return c_p_grid[n - 1];
  } else {

    for (i = 1; i < n; i++) {
      if (T_grid[i - 1] <= T_in_K && T_in_K <= T_grid[i]) {
        c_p_dry_air = c_p_grid[i - 1] + (c_p_grid[i] - c_p_grid[i - 1]) / (T_grid[i] - T_grid[i - 1]) * (T_in_K - T_grid[i - 1]);
        /* fprintf (stderr,"T_(i-1)=%5.2f, T=%5.2f, T_(i)=%5.2f, c_(i-1)=%6.2f, c=%6.2f, c_(i)=%6.2f (dry air) \n", */
        /*                  T_grid[i-1],T_in_K,T_grid[i],c_p_grid[i-1],c_p_dry_air,c_p_grid[i]);                    */
        break;
      }
    }

    if (c_p_dry_air < 0.0) {
      fprintf (stderr, "Error, determening specific heat of dry air, T = %6.2f K (in ancillary.c)\n", T_in_K);
      return -1.0;
    }

    return c_p_dry_air;
  }
}

/***********************************************************************************************/
/* function: specific_heat_capacity_water_vapour                                               */
/* interpolate a table to get the temperature dependent specific heat capacity of water vapour */
/* input     float T in K                                                                      */
/* output    float specific heat capacity of water vapour in J/(kg K)                          */
/* ulrich hamann: 2007-03-31                                                                   */
/***********************************************************************************************/

float specific_heat_capacity_water_vapour (float T_in_K, int quiet) {

  const int n      = 14;
  int       i      = NOT_DEFINED_INTEGER;
  float     c_p_wv = NOT_DEFINED_FLOAT;
  float     c_p_wv_grid[14] =
    {1850.0, 1851.0, 1852.0, 1855.0, 1859.0, 1864.0, 1871.0, 1880.0, 1890.0, 1901.0, 1926.0, 1954.0, 1984.0, 2015.0};
  float T_grid[14] = {175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 450.0, 500.0, 550.0, 600.0};

  if (T_in_K < T_grid[0]) {
    if (!quiet) {
      fprintf (stderr, " *** Warning, getting specific heat capacity of water vapour,\n");
      fprintf (stderr, "     temperature %7.2f is below of the tabled temperature region.\n", T_in_K);
    }
    return c_p_wv_grid[0];
  } else if (T_in_K > T_grid[n - 1]) {
    if (!quiet) {
      fprintf (stderr, " *** Warning, getting specific heat capacity of water vapour,\n");
      fprintf (stderr, "     temperature %7.2f is above of the tabled temperature region.\n", T_in_K);
    }
    return c_p_wv_grid[n - 1];
  } else {

    for (i = 1; i < n; i++) {
      if (T_grid[i - 1] <= T_in_K && T_in_K <= T_grid[i]) {
        c_p_wv =
          c_p_wv_grid[i - 1] + (c_p_wv_grid[i] - c_p_wv_grid[i - 1]) / (T_grid[i] - T_grid[i - 1]) * (T_in_K - T_grid[i - 1]);
        /* fprintf (stderr,"T_(i-1)=%5.2f, T=%5.2f, T_(i)=%5.2f, c_(i-1)=%6.2f, c=%6.2f, c_(i)=%6.2f (water vapour) \n", */
        /*                  T_grid[i-1],T_in_K,T_grid[i],c_p_wv_grid[i-1],c_p_wv,c_p_wv_grid[i]);                        */
        break;
      }
    }

    if (c_p_wv < 0.0) {
      fprintf (stderr, "Error, determening specific heat of water vapour, T = %6.2f K (in ancillary.c)\n", T_in_K);
      return -1.0;
    }

    return c_p_wv;
  }
}

/***********************************************************************************/
/* Function: my_timegm                                                             */
/* Description:                                                                    */
/*  converts a tm-struct (date as year, month, day, hour, min, ...)                */
/*  into time in s since 1.1.1970                                                  */
/*  (reverse function to gmtime)                                                   */
/*  (essentially the same as timegm, which is part of many GNU C libraries,        */
/*   but not all, and therefor not portable)                                       */
/*                                                                                 */
/* Parameters:                                                                     */
/*     struct tm *tm      input date                                               */
/*                                                                                 */
/* Return value:                                                                   */
/*     time_t             time                                                     */
/*                        == -1 if there was some error                            */
/*                                                                                 */
/* Example:                                                                        */
/* Files:      ancillary.c                                                         */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2003   Roger Dingledine                                                  */
/*               Delivered-to the SEUL-project under the Gnu Public Licence        */
/*               http://archives.seul.org/or/cvs/Oct-2003/msg00123.html            */
/*    Jan 2007   Ulrich Hamann                                                     */
/*               implemented into libRadtran                                       */
/*                                                                                 */
/***********************************************************************************/

time_t my_timegm (struct tm* tm) {
  time_t ret;
  char*  tz;

  tz = getenv ("TZ");
  setenv ("TZ", "", 1);
  tzset();
  ret = mktime (tm);
  if (tz)
    setenv ("TZ", tz, 1);
  else
    unsetenv ("TZ");
  tzset();
  return ret;
}

/* check if two floats are equal - still not optimal, but how to improve? */
/* THIS CODE IS IDENTICALLY IN CLOUD3D.C                                   */

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

double nm_to_inv_cm (double wavelength_nm) {
  double wavenumber_inv_cm;

  wavenumber_inv_cm = 1 / (wavelength_nm * 1E-07);

  return wavenumber_inv_cm;
}

double inv_cm_to_nm (double wavenumber_inv_cm) {
  double wavelength_nm;

  wavelength_nm = 1E+07 / wavenumber_inv_cm;

  return wavelength_nm;
}

double polarizability_anisotropy_N2 (double nu) {

  double gam = 0;

  gam = -6.01466e-25 + 2.38557e-14 / (1.86099e+10 - nu * nu);

  return gam;
}

double polarizability_anisotropy_O2 (double nu) {
  double gam;

  gam = 7.149e-26 + 4.59364e-15 / (4.82716e+9 - nu * nu);

  return gam;
}

struct E_rot_N2 {

  int    gJ;
  int    J;
  double E;
};

struct E_rot_O2 {

  int    N;
  int    J;
  double E;
};

struct E_rot_trans_N2 {

  int    J;
  int    Jp; /* J prime */
  int    gJ; /* Statistical weight */
  double E;
  double delta_E;
  double cpt; /* Plazcek-Teller coefficient */
};

struct E_rot_trans_O2 {

  int    N;
  int    J;
  int    Np; /* N prime */
  int    Jp; /* J prime */
  double E;
  double delta_E;
  double cpt; /* Plazcek-Teller coefficient */
};

int crs_raman_N2 (double lambda, int n_transitions, float temper, double*** crs, int* iv, int verbose) {

  /* Calculate Raman shifted wavelengths for input wavelength lambda in nm and */
  /* temperature temper. n_transitions are considered. Wavelength are returned */
  /* in first column of crs and cross section in the second column. Cross      */
  /* section is weighted by the volume mixing ratio.                           */

  double c  = 299792458;      /* Speed of light, m/s, from wikipedia */
  double h  = 6.62606896e-34; /* Planck constant, Js, from wikipedia */
  double hc = h * c * 100;    /* Converted to J*cm */

  double nm_to_cm = 1e-07; /* Convert from nm to cm */

  double raman_const = 256.0 * pow (M_PI, 5.0) / 27.0;

  /* All the numbers below for the Rotational energy states of O2 are */
  /* from the report: Accounting for Raman Scattering in DOAS,        */
  /* J.F. de Haan, SN-OMIE-KNMI-409, Version 1.0,  May 18, 2003       */

  /*
      N: Nuclear rotational angular momentum quantum number
      gJ: statistical weight factor
      E: Rotational energy of state in cm-1
      From Table A-1
   */
  int             N_E = 31;
  struct E_rot_N2 E_rot[N_E];
  E_rot[0].J   = 0;
  E_rot[0].gJ  = 6;
  E_rot[0].E   = 0.0000;
  E_rot[1].J   = 1;
  E_rot[1].gJ  = 3;
  E_rot[1].E   = 3.9791;
  E_rot[2].J   = 2;
  E_rot[2].gJ  = 6;
  E_rot[2].E   = 11.9373;
  E_rot[3].J   = 3;
  E_rot[3].gJ  = 3;
  E_rot[3].E   = 23.8741;
  E_rot[4].J   = 4;
  E_rot[4].gJ  = 6;
  E_rot[4].E   = 39.7892;
  E_rot[5].J   = 5;
  E_rot[5].gJ  = 3;
  E_rot[5].E   = 59.6821;
  E_rot[6].J   = 6;
  E_rot[6].gJ  = 6;
  E_rot[6].E   = 83.5521;
  E_rot[7].J   = 7;
  E_rot[7].gJ  = 3;
  E_rot[7].E   = 111.3983;
  E_rot[8].J   = 8;
  E_rot[8].gJ  = 6;
  E_rot[8].E   = 143.2197;
  E_rot[9].J   = 9;
  E_rot[9].gJ  = 3;
  E_rot[9].E   = 179.0154;
  E_rot[10].J  = 10;
  E_rot[10].gJ = 6;
  E_rot[10].E  = 218.7839;
  E_rot[11].J  = 11;
  E_rot[11].gJ = 3;
  E_rot[11].E  = 262.5240;
  E_rot[12].J  = 12;
  E_rot[12].gJ = 6;
  E_rot[12].E  = 310.2341;
  E_rot[13].J  = 13;
  E_rot[13].gJ = 3;
  E_rot[13].E  = 361.9126;
  E_rot[14].J  = 14;
  E_rot[14].gJ = 6;
  E_rot[14].E  = 417.5576;
  E_rot[15].J  = 15;
  E_rot[15].gJ = 3;
  E_rot[15].E  = 477.1673;
  E_rot[16].J  = 16;
  E_rot[16].gJ = 6;
  E_rot[16].E  = 540.7395;
  E_rot[17].J  = 17;
  E_rot[17].gJ = 3;
  E_rot[17].E  = 608.2722;
  E_rot[18].J  = 18;
  E_rot[18].gJ = 6;
  E_rot[18].E  = 679.7628;
  E_rot[19].J  = 19;
  E_rot[19].gJ = 3;
  E_rot[19].E  = 755.2090;
  E_rot[20].J  = 20;
  E_rot[20].gJ = 6;
  E_rot[20].E  = 834.6081;
  E_rot[21].J  = 21;
  E_rot[21].gJ = 3;
  E_rot[21].E  = 917.9574;
  E_rot[22].J  = 22;
  E_rot[22].gJ = 6;
  E_rot[22].E  = 1005.2540;
  E_rot[23].J  = 23;
  E_rot[23].gJ = 3;
  E_rot[23].E  = 1096.4948;
  E_rot[24].J  = 24;
  E_rot[24].gJ = 6;
  E_rot[24].E  = 1191.6766;
  E_rot[25].J  = 25;
  E_rot[25].gJ = 3;
  E_rot[25].E  = 1290.7963;
  E_rot[26].J  = 26;
  E_rot[26].gJ = 6;
  E_rot[26].E  = 1393.8503;
  E_rot[27].J  = 27;
  E_rot[27].gJ = 3;
  E_rot[27].E  = 1500.8350;
  E_rot[28].J  = 28;
  E_rot[28].gJ = 6;
  E_rot[28].E  = 1611.7467;
  E_rot[29].J  = 29;
  E_rot[29].gJ = 3;
  E_rot[29].E  = 1726.5816;
  E_rot[30].J  = 30;
  E_rot[30].gJ = 6;
  E_rot[30].E  = 1845.3358;

  /*
      N: Nuclear rotational angular momentum quantum number
      J: Total angular momentum quantum number
      E: Rotational energy of state in cm-1
      From Table A-2
   */
  int                   N_E_trans = 48;
  struct E_rot_trans_N2 Etr[N_E_trans];
  Etr[0].J        = 25;
  Etr[0].Jp       = 23;
  Etr[0].gJ       = 3;
  Etr[0].E        = 1290.7963;
  Etr[0].delta_E  = -194.3015;
  Etr[0].cpt      = 0.3601;
  Etr[1].J        = 24;
  Etr[1].Jp       = 22;
  Etr[1].gJ       = 6;
  Etr[1].E        = 1191.6766;
  Etr[1].delta_E  = -186.4226;
  Etr[1].cpt      = 0.3595;
  Etr[2].J        = 23;
  Etr[2].Jp       = 21;
  Etr[2].gJ       = 3;
  Etr[2].E        = 1096.4948;
  Etr[2].delta_E  = -178.5374;
  Etr[2].cpt      = 0.3589;
  Etr[3].J        = 22;
  Etr[3].Jp       = 20;
  Etr[3].gJ       = 6;
  Etr[3].E        = 1005.2540;
  Etr[3].delta_E  = -170.6459;
  Etr[3].cpt      = 0.3581;
  Etr[4].J        = 21;
  Etr[4].Jp       = 19;
  Etr[4].gJ       = 3;
  Etr[4].E        = 917.9574;
  Etr[4].delta_E  = -162.7484;
  Etr[4].cpt      = 0.3573;
  Etr[5].J        = 20;
  Etr[5].Jp       = 18;
  Etr[5].gJ       = 6;
  Etr[5].E        = 834.6081;
  Etr[5].delta_E  = -154.8453;
  Etr[5].cpt      = 0.3565;
  Etr[6].J        = 19;
  Etr[6].Jp       = 17;
  Etr[6].gJ       = 3;
  Etr[6].E        = 755.2090;
  Etr[6].delta_E  = -146.9368;
  Etr[6].cpt      = 0.3555;
  Etr[7].J        = 18;
  Etr[7].Jp       = 16;
  Etr[7].gJ       = 6;
  Etr[7].E        = 679.7628;
  Etr[7].delta_E  = -139.0233;
  Etr[7].cpt      = 0.3544;
  Etr[8].J        = 17;
  Etr[8].Jp       = 15;
  Etr[8].gJ       = 3;
  Etr[8].E        = 608.2722;
  Etr[8].delta_E  = -131.1049;
  Etr[8].cpt      = 0.3532;
  Etr[9].J        = 16;
  Etr[9].Jp       = 14;
  Etr[9].gJ       = 6;
  Etr[9].E        = 540.7395;
  Etr[9].delta_E  = -123.1819;
  Etr[9].cpt      = 0.3519;
  Etr[10].J       = 15;
  Etr[10].Jp      = 13;
  Etr[10].gJ      = 3;
  Etr[10].E       = 477.1673;
  Etr[10].delta_E = -115.2547;
  Etr[10].cpt     = 0.3504;
  Etr[11].J       = 14;
  Etr[11].Jp      = 12;
  Etr[11].gJ      = 6;
  Etr[11].E       = 417.5576;
  Etr[11].delta_E = -107.3235;
  Etr[11].cpt     = 0.3487;
  Etr[12].J       = 13;
  Etr[12].Jp      = 11;
  Etr[12].gJ      = 3;
  Etr[12].E       = 361.9126;
  Etr[12].delta_E = -99.3886;
  Etr[12].cpt     = 0.3467;
  Etr[13].J       = 12;
  Etr[13].Jp      = 10;
  Etr[13].gJ      = 6;
  Etr[13].E       = 310.2341;
  Etr[13].delta_E = -91.4502;
  Etr[13].cpt     = 0.3443;
  Etr[14].J       = 11;
  Etr[14].Jp      = 9;
  Etr[14].gJ      = 3;
  Etr[14].E       = 262.5240;
  Etr[14].delta_E = -83.5086;
  Etr[14].cpt     = 0.3416;
  Etr[15].J       = 10;
  Etr[15].Jp      = 8;
  Etr[15].gJ      = 6;
  Etr[15].E       = 218.7839;
  Etr[15].delta_E = -75.5642;
  Etr[15].cpt     = 0.3383;
  Etr[16].J       = 9;
  Etr[16].Jp      = 7;
  Etr[16].gJ      = 3;
  Etr[16].E       = 179.0154;
  Etr[16].delta_E = -67.6171;
  Etr[16].cpt     = 0.3344;
  Etr[17].J       = 8;
  Etr[17].Jp      = 6;
  Etr[17].gJ      = 6;
  Etr[17].E       = 143.2197;
  Etr[17].delta_E = -59.6676;
  Etr[17].cpt     = 0.3294;
  Etr[18].J       = 7;
  Etr[18].Jp      = 5;
  Etr[18].gJ      = 3;
  Etr[18].E       = 111.3983;
  Etr[18].delta_E = -51.7162;
  Etr[18].cpt     = 0.3231;
  Etr[19].J       = 6;
  Etr[19].Jp      = 4;
  Etr[19].gJ      = 6;
  Etr[19].E       = 83.5521;
  Etr[19].delta_E = -43.7629;
  Etr[19].cpt     = 0.3147;
  Etr[20].J       = 5;
  Etr[20].Jp      = 3;
  Etr[20].gJ      = 3;
  Etr[20].E       = 59.6821;
  Etr[20].delta_E = -35.8080;
  Etr[20].cpt     = 0.3030;
  Etr[21].J       = 4;
  Etr[21].Jp      = 2;
  Etr[21].gJ      = 6;
  Etr[21].E       = 39.7892;
  Etr[21].delta_E = -27.8519;
  Etr[21].cpt     = 0.2857;
  Etr[22].J       = 3;
  Etr[22].Jp      = 1;
  Etr[22].gJ      = 3;
  Etr[22].E       = 23.8741;
  Etr[22].delta_E = -19.8950;
  Etr[22].cpt     = 0.2571;
  Etr[23].J       = 2;
  Etr[23].Jp      = 0;
  Etr[23].gJ      = 6;
  Etr[23].E       = 11.9373;
  Etr[23].delta_E = -11.9373;
  Etr[23].cpt     = 0.2000;
  Etr[24].J       = 0;
  Etr[24].Jp      = 2;
  Etr[24].gJ      = 6;
  Etr[24].E       = 0.0000;
  Etr[24].delta_E = 11.9373;
  Etr[24].cpt     = 1.0000;
  Etr[25].J       = 1;
  Etr[25].Jp      = 3;
  Etr[25].gJ      = 3;
  Etr[25].E       = 3.9791;
  Etr[25].delta_E = 19.8950;
  Etr[25].cpt     = 0.6000;
  Etr[26].J       = 2;
  Etr[26].Jp      = 4;
  Etr[26].gJ      = 6;
  Etr[26].E       = 11.9373;
  Etr[26].delta_E = 27.8519;
  Etr[26].cpt     = 0.5143;
  Etr[27].J       = 3;
  Etr[27].Jp      = 5;
  Etr[27].gJ      = 3;
  Etr[27].E       = 23.8741;
  Etr[27].delta_E = 35.8080;
  Etr[27].cpt     = 0.4762;
  Etr[28].J       = 4;
  Etr[28].Jp      = 6;
  Etr[28].gJ      = 6;
  Etr[28].E       = 39.7892;
  Etr[28].delta_E = 43.7629;
  Etr[28].cpt     = 0.4545;
  Etr[29].J       = 5;
  Etr[29].Jp      = 7;
  Etr[29].gJ      = 3;
  Etr[29].E       = 59.6821;
  Etr[29].delta_E = 51.7162;
  Etr[29].cpt     = 0.4406;
  Etr[30].J       = 6;
  Etr[30].Jp      = 8;
  Etr[30].gJ      = 6;
  Etr[30].E       = 83.5521;
  Etr[30].delta_E = 59.6676;
  Etr[30].cpt     = 0.4308;
  Etr[31].J       = 7;
  Etr[31].Jp      = 9;
  Etr[31].gJ      = 3;
  Etr[31].E       = 111.3983;
  Etr[31].delta_E = 67.6171;
  Etr[31].cpt     = 0.4235;
  Etr[32].J       = 8;
  Etr[32].Jp      = 10;
  Etr[32].gJ      = 6;
  Etr[32].E       = 143.2197;
  Etr[32].delta_E = 75.5642;
  Etr[32].cpt     = 0.4180;
  Etr[33].J       = 9;
  Etr[33].Jp      = 11;
  Etr[33].gJ      = 3;
  Etr[33].E       = 179.0154;
  Etr[33].delta_E = 83.5086;
  Etr[33].cpt     = 0.4135;
  Etr[34].J       = 10;
  Etr[34].Jp      = 12;
  Etr[34].gJ      = 6;
  Etr[34].E       = 218.7839;
  Etr[34].delta_E = 91.4502;
  Etr[34].cpt     = 0.4099;
  Etr[35].J       = 11;
  Etr[35].Jp      = 13;
  Etr[35].gJ      = 3;
  Etr[35].E       = 262.5240;
  Etr[35].delta_E = 99.3886;
  Etr[35].cpt     = 0.4070;
  Etr[36].J       = 12;
  Etr[36].Jp      = 14;
  Etr[36].gJ      = 6;
  Etr[36].E       = 310.2341;
  Etr[36].delta_E = 107.3235;
  Etr[36].cpt     = 0.4044;
  Etr[37].J       = 13;
  Etr[37].Jp      = 15;
  Etr[37].gJ      = 3;
  Etr[37].E       = 361.9126;
  Etr[37].delta_E = 115.2547;
  Etr[37].cpt     = 0.4023;
  Etr[38].J       = 14;
  Etr[38].Jp      = 16;
  Etr[38].gJ      = 6;
  Etr[38].E       = 417.5576;
  Etr[38].delta_E = 123.1819;
  Etr[38].cpt     = 0.4004;
  Etr[39].J       = 15;
  Etr[39].Jp      = 17;
  Etr[39].gJ      = 3;
  Etr[39].E       = 477.1673;
  Etr[39].delta_E = 131.1049;
  Etr[39].cpt     = 0.3988;
  Etr[40].J       = 16;
  Etr[40].Jp      = 18;
  Etr[40].gJ      = 6;
  Etr[40].E       = 540.7395;
  Etr[40].delta_E = 139.0233;
  Etr[40].cpt     = 0.3974;
  Etr[41].J       = 17;
  Etr[41].Jp      = 19;
  Etr[41].gJ      = 3;
  Etr[41].E       = 608.2722;
  Etr[41].delta_E = 146.9368;
  Etr[41].cpt     = 0.3961;
  Etr[42].J       = 18;
  Etr[42].Jp      = 20;
  Etr[42].gJ      = 6;
  Etr[42].E       = 679.7628;
  Etr[42].delta_E = 154.8453;
  Etr[42].cpt     = 0.3950;
  Etr[43].J       = 19;
  Etr[43].Jp      = 21;
  Etr[43].gJ      = 3;
  Etr[43].E       = 755.2090;
  Etr[43].delta_E = 162.7484;
  Etr[43].cpt     = 0.3940;
  Etr[44].J       = 20;
  Etr[44].Jp      = 22;
  Etr[44].gJ      = 6;
  Etr[44].E       = 834.6081;
  Etr[44].delta_E = 170.6459;
  Etr[44].cpt     = 0.3931;
  Etr[45].J       = 21;
  Etr[45].Jp      = 23;
  Etr[45].gJ      = 3;
  Etr[45].E       = 917.9574;
  Etr[45].delta_E = 178.5374;
  Etr[45].cpt     = 0.3922;
  Etr[46].J       = 22;
  Etr[46].Jp      = 24;
  Etr[46].gJ      = 6;
  Etr[46].E       = 1005.2540;
  Etr[46].delta_E = 186.4226;
  Etr[46].cpt     = 0.3915;
  Etr[47].J       = 23;
  Etr[47].Jp      = 25;
  Etr[47].gJ      = 3;
  Etr[47].E       = 1096.4948;
  Etr[47].delta_E = 194.3015;
  Etr[47].cpt     = 0.3908;

  double E_J           = 0; /* Rotational energy of state J */
  double W_J           = 0; /* Fraction of molecules in the rotational state J at temperature T */
  double b_J           = 0; /* Placzek-Teller coefficient */
  double lambda_inv_cm = 0, lambda_cm = 0, lambda_shifted_inv_cm = 0, lambda_shifted_nm = 0, lambda_shifted_cm = 0;
  double delta_nu = 0;
  double crs_i = 0, crs_j = 0;
  double fNJ = 0, Z = 0; /* Eq 5 and 10 */
  double gam_i = 0, gam_j = 0;
  double volume_mixing_ratio = 0.7905;
  double sum                 = 0;

  int J   = 0; /* Rotational state J */
  int g_J = 0;

  int status = 0;
  int ivi = 0, is = 0;

  if (verbose) {
    fprintf (stderr, "Calculating Raman scattering wavelength shifts and cross sections ");
    fprintf (stderr, "for N2, wavelength %5.1f nm, number of transitions %d, T=%6.2f.\n", lambda, n_transitions, temper);
    fprintf (stderr,
             "             %3s %1s %7s %8s %7s %9s %10s %9s %17s %9s %12s %7s %12s %12s %11s %14s %12s\n",
             "ivi",
             "J",
             "E_J",
             "g_J",
             "W_J",
             "b_J",
             "W_J*b_J",
             "wvl",
             "wvl_shift_cm-1",
             "delta_nu",
             "wvl_shift_nm",
             "gam_i",
             "gam_j",
             "crs_i",
             "crs_j",
             "crs_i*vmr",
             "crs_j*vmr");
  }

  ivi = *iv;

  sum = 0;
  for (is = 0; is < N_E; is++) {
    g_J = E_rot[is].gJ;
    J   = E_rot[is].J;
    E_J = E_rot[is].E * hc;
    sum += g_J * (2 * J + 1) * exp (-E_J / (BOLTZMANN * temper));
  }
  Z = sum;

  for (is = 0; is < n_transitions; is++) {

    g_J      = Etr[is].gJ;
    J        = Etr[is].J;
    E_J      = Etr[is].E * hc;
    b_J      = Etr[is].cpt;
    delta_nu = -Etr[is].delta_E; /* Negative delta_E corresponds to a photon with a larger */
                                 /* energy, shorter wavelength than the incident photon    */
                                 /* delta_E is the change in rotational energy of the      */
                                 /* molecule. The minus puts the photon in the right       */
                                 /* wavelength.                                            */

    W_J                   = g_J * (2 * J + 1) * exp (-E_J / (BOLTZMANN * temper));
    fNJ                   = W_J / Z;
    lambda_inv_cm         = nm_to_inv_cm (lambda);
    lambda_shifted_inv_cm = lambda_inv_cm + delta_nu;
    lambda_shifted_nm     = inv_cm_to_nm (lambda_shifted_inv_cm);
    lambda_shifted_cm     = lambda_shifted_nm * nm_to_cm;
    lambda_cm             = lambda * nm_to_cm;
    gam_i                 = polarizability_anisotropy_N2 (lambda_shifted_inv_cm);
    gam_j                 = polarizability_anisotropy_N2 (lambda_inv_cm);
    crs_i                 = fNJ * gam_i * gam_i * raman_const * b_J / pow (lambda_shifted_cm, 4);
    crs_j                 = fNJ * gam_j * gam_j * raman_const * b_J / pow (lambda_cm, 4);
    (*crs)[ivi][0]        = lambda_shifted_nm;
    (*crs)[ivi][1]        = crs_i * volume_mixing_ratio;
    (*crs)[ivi][2]        = crs_j * volume_mixing_ratio;
    if (verbose)
      fprintf (stderr,
               "Raman_crs_N2 %2d %2d %12.6e %2d %12.6e %6.3f %12.6e %10.6f %10.6f %10.4f %10.6f %12.6e %12.6e %12.6e %12.6e %12.6e "
               "%12.6e\n",
               ivi,
               J,
               E_J,
               g_J,
               W_J,
               b_J,
               W_J * b_J,
               lambda,
               lambda_shifted_inv_cm,
               delta_nu,
               lambda_shifted_nm,
               gam_i,
               gam_j,
               crs_i,
               crs_j,
               crs_i * volume_mixing_ratio,
               crs_j * volume_mixing_ratio);
    ivi++;
  }

  *iv = ivi;

  return status;
}

int crs_raman_O2 (double lambda, int n_transitions, float temper, double*** crs, int* iv, int verbose) {

  /* Calculate Raman shifted wavelengths for input wavelength lambda in nm and */
  /* temperature temper. n_transitions are considered. Wavelength are returned */
  /* in first column of crs and cross section in the second column. Cross      */
  /* section is weighted by the volume mixing ratio.                           */

  double c  = 299792458;      /* Speed of light, m/s, from wikipedia */
  double h  = 6.62606896e-34; /* Planck constant, Js, from wikipedia */
  double hc = h * c * 100;    /* Converted to J*cm */

  double nm_to_cm = 1e-07; /* Convert from nm to cm */

  double raman_const = 256.0 * pow (M_PI, 5.0) / 27.0;

  /* All the numbers below for the Rotational energy states of O2 are */
  /* from the report: Accounting for Raman Scattering in DOAS,        */
  /* J.F. de Haan, SN-OMIE-KNMI-409, Version 1.0,  May 18, 2003       */

  /*
      N: Nuclear rotational angular momentum quantum number
      J: Total angular momentum quantum number
      E: Rotational energy of state in cm-1
      From Table A-1
   */
  int             N_E = 54;
  struct E_rot_O2 E_rot[N_E];
  E_rot[0].N  = 1;
  E_rot[0].J  = 0;
  E_rot[0].E  = 0.0000;
  E_rot[1].N  = 1;
  E_rot[1].J  = 2;
  E_rot[1].E  = 2.0843;
  E_rot[2].N  = 1;
  E_rot[2].J  = 1;
  E_rot[2].E  = 3.9611;
  E_rot[3].N  = 3;
  E_rot[3].J  = 2;
  E_rot[3].E  = 16.2529;
  E_rot[4].N  = 3;
  E_rot[4].J  = 4;
  E_rot[4].E  = 16.3876;
  E_rot[5].N  = 3;
  E_rot[5].J  = 3;
  E_rot[5].E  = 18.3372;
  E_rot[6].N  = 5;
  E_rot[6].J  = 4;
  E_rot[6].E  = 42.2001;
  E_rot[7].N  = 5;
  E_rot[7].J  = 6;
  E_rot[7].E  = 42.2240;
  E_rot[8].N  = 5;
  E_rot[8].J  = 5;
  E_rot[8].E  = 44.2117;
  E_rot[9].N  = 7;
  E_rot[9].J  = 8;
  E_rot[9].E  = 79.5646;
  E_rot[10].N = 7;
  E_rot[10].J = 6;
  E_rot[10].E = 79.6070;
  E_rot[11].N = 7;
  E_rot[11].J = 7;
  E_rot[11].E = 81.5805;
  E_rot[12].N = 9;
  E_rot[12].J = 10;
  E_rot[12].E = 128.3978;
  E_rot[13].N = 9;
  E_rot[13].J = 8;
  E_rot[13].E = 128.4921;
  E_rot[14].N = 9;
  E_rot[14].J = 9;
  E_rot[14].E = 130.4376;
  E_rot[15].N = 11;
  E_rot[15].J = 12;
  E_rot[15].E = 188.7135;
  E_rot[16].N = 11;
  E_rot[16].J = 10;
  E_rot[16].E = 188.8532;
  E_rot[17].N = 11;
  E_rot[17].J = 11;
  E_rot[17].E = 190.7749;
  E_rot[18].N = 13;
  E_rot[18].J = 14;
  E_rot[18].E = 260.5011;
  E_rot[19].N = 13;
  E_rot[19].J = 12;
  E_rot[19].E = 260.6826;
  E_rot[20].N = 13;
  E_rot[20].J = 13;
  E_rot[20].E = 262.5829;
  E_rot[21].N = 15;
  E_rot[21].J = 16;
  E_rot[21].E = 343.7484;
  E_rot[22].N = 15;
  E_rot[22].J = 14;
  E_rot[22].E = 343.9697;
  E_rot[23].N = 15;
  E_rot[23].J = 15;
  E_rot[23].E = 345.8500;
  E_rot[24].N = 17;
  E_rot[24].J = 18;
  E_rot[24].E = 438.4418;
  E_rot[25].N = 17;
  E_rot[25].J = 16;
  E_rot[25].E = 438.7015;
  E_rot[26].N = 17;
  E_rot[26].J = 17;
  E_rot[26].E = 440.5620;
  E_rot[27].N = 19;
  E_rot[27].J = 20;
  E_rot[27].E = 544.5658;
  E_rot[28].N = 19;
  E_rot[28].J = 18;
  E_rot[28].E = 544.8628;
  E_rot[29].N = 19;
  E_rot[29].J = 19;
  E_rot[29].E = 546.7050;
  E_rot[30].N = 21;
  E_rot[30].J = 22;
  E_rot[30].E = 662.1030;
  E_rot[31].N = 21;
  E_rot[31].J = 20;
  E_rot[31].E = 662.4368;
  E_rot[32].N = 21;
  E_rot[32].J = 21;
  E_rot[32].E = 664.2610;
  E_rot[33].N = 23;
  E_rot[33].J = 24;
  E_rot[33].E = 791.0344;
  E_rot[34].N = 23;
  E_rot[34].J = 22;
  E_rot[34].E = 791.4045;
  E_rot[35].N = 23;
  E_rot[35].J = 23;
  E_rot[35].E = 793.2100;
  E_rot[36].N = 25;
  E_rot[36].J = 26;
  E_rot[36].E = 931.3390;
  E_rot[37].N = 25;
  E_rot[37].J = 24;
  E_rot[37].E = 931.7450;
  E_rot[38].N = 25;
  E_rot[38].J = 25;
  E_rot[38].E = 933.5330;
  E_rot[39].N = 27;
  E_rot[39].J = 28;
  E_rot[39].E = 1082.9941;
  E_rot[40].N = 27;
  E_rot[40].J = 26;
  E_rot[40].E = 1083.4356;
  E_rot[41].N = 27;
  E_rot[41].J = 27;
  E_rot[41].E = 1085.2060;
  E_rot[42].N = 29;
  E_rot[42].J = 30;
  E_rot[42].E = 1245.9750;
  E_rot[43].N = 29;
  E_rot[43].J = 28;
  E_rot[43].E = 1246.4518;
  E_rot[44].N = 29;
  E_rot[44].J = 29;
  E_rot[44].E = 1248.2040;
  E_rot[45].N = 31;
  E_rot[45].J = 32;
  E_rot[45].E = 1420.2552;
  E_rot[46].N = 31;
  E_rot[46].J = 30;
  E_rot[46].E = 1420.7672;
  E_rot[47].N = 31;
  E_rot[47].J = 31;
  E_rot[47].E = 1422.5020;
  E_rot[48].N = 33;
  E_rot[48].J = 34;
  E_rot[48].E = 1605.8064;
  E_rot[49].N = 33;
  E_rot[49].J = 32;
  E_rot[49].E = 1606.3533;
  E_rot[50].N = 33;
  E_rot[50].J = 33;
  E_rot[50].E = 1608.0710;
  E_rot[51].N = 35;
  E_rot[51].J = 36;
  E_rot[51].E = 1802.5983;
  E_rot[52].N = 35;
  E_rot[52].J = 34;
  E_rot[52].E = 1803.1802;
  E_rot[53].N = 35;
  E_rot[53].J = 35;
  E_rot[53].E = 1804.8810;

  /*
      N: Nuclear rotational angular momentum quantum number
      J: Total angular momentum quantum number
      E: Rotational energy of state in cm-1
      From Table A-3
   */
  int                   N_E_trans = 185;
  struct E_rot_trans_O2 Etr[N_E_trans];
  Etr[0].N         = 33;
  Etr[0].J         = 32;
  Etr[0].Np        = 31;
  Etr[0].Jp        = 30;
  Etr[0].E         = 1606.3533;
  Etr[0].delta_E   = -185.5861;
  Etr[0].cpt       = 0.3630;
  Etr[1].N         = 33;
  Etr[1].J         = 33;
  Etr[1].Np        = 31;
  Etr[1].Jp        = 31;
  Etr[1].E         = 1608.0710;
  Etr[1].delta_E   = -185.5690;
  Etr[1].cpt       = 0.3630;
  Etr[2].N         = 33;
  Etr[2].J         = 34;
  Etr[2].Np        = 31;
  Etr[2].Jp        = 32;
  Etr[2].E         = 1605.8064;
  Etr[2].delta_E   = -185.5512;
  Etr[2].cpt       = 0.3637;
  Etr[3].N         = 31;
  Etr[3].J         = 30;
  Etr[3].Np        = 29;
  Etr[3].Jp        = 28;
  Etr[3].E         = 1420.7672;
  Etr[3].delta_E   = -174.3154;
  Etr[3].cpt       = 0.3622;
  Etr[4].N         = 31;
  Etr[4].J         = 31;
  Etr[4].Np        = 29;
  Etr[4].Jp        = 29;
  Etr[4].E         = 1422.5020;
  Etr[4].delta_E   = -174.2980;
  Etr[4].cpt       = 0.3622;
  Etr[5].N         = 31;
  Etr[5].J         = 32;
  Etr[5].Np        = 29;
  Etr[5].Jp        = 30;
  Etr[5].E         = 1420.2552;
  Etr[5].delta_E   = -174.2802;
  Etr[5].cpt       = 0.3630;
  Etr[6].N         = 29;
  Etr[6].J         = 28;
  Etr[6].Np        = 27;
  Etr[6].Jp        = 26;
  Etr[6].E         = 1246.4518;
  Etr[6].delta_E   = -163.0162;
  Etr[6].cpt       = 0.3613;
  Etr[7].N         = 29;
  Etr[7].J         = 29;
  Etr[7].Np        = 27;
  Etr[7].Jp        = 27;
  Etr[7].E         = 1248.2040;
  Etr[7].delta_E   = -162.9980;
  Etr[7].cpt       = 0.3613;
  Etr[8].N         = 29;
  Etr[8].J         = 30;
  Etr[8].Np        = 27;
  Etr[8].Jp        = 28;
  Etr[8].E         = 1245.9750;
  Etr[8].delta_E   = -162.9809;
  Etr[8].cpt       = 0.3622;
  Etr[9].N         = 27;
  Etr[9].J         = 26;
  Etr[9].Np        = 25;
  Etr[9].Jp        = 24;
  Etr[9].E         = 1083.4355;
  Etr[9].delta_E   = -151.6906;
  Etr[9].cpt       = 0.3602;
  Etr[10].N        = 27;
  Etr[10].J        = 27;
  Etr[10].Np       = 25;
  Etr[10].Jp       = 25;
  Etr[10].E        = 1085.2061;
  Etr[10].delta_E  = -151.6730;
  Etr[10].cpt      = 0.3602;
  Etr[11].N        = 27;
  Etr[11].J        = 28;
  Etr[11].Np       = 25;
  Etr[11].Jp       = 26;
  Etr[11].E        = 1082.9941;
  Etr[11].delta_E  = -151.6551;
  Etr[11].cpt      = 0.3612;
  Etr[12].N        = 25;
  Etr[12].J        = 24;
  Etr[12].Np       = 23;
  Etr[12].Jp       = 22;
  Etr[12].E        = 931.7450;
  Etr[12].delta_E  = -140.3405;
  Etr[12].cpt      = 0.3589;
  Etr[13].N        = 25;
  Etr[13].J        = 25;
  Etr[13].Np       = 23;
  Etr[13].Jp       = 23;
  Etr[13].E        = 933.5330;
  Etr[13].delta_E  = -140.3230;
  Etr[13].cpt      = 0.3589;
  Etr[14].N        = 25;
  Etr[14].J        = 26;
  Etr[14].Np       = 23;
  Etr[14].Jp       = 24;
  Etr[14].E        = 931.3390;
  Etr[14].delta_E  = -140.3046;
  Etr[14].cpt      = 0.3601;
  Etr[15].N        = 23;
  Etr[15].J        = 22;
  Etr[15].Np       = 21;
  Etr[15].Jp       = 20;
  Etr[15].E        = 791.4045;
  Etr[15].delta_E  = -128.9677;
  Etr[15].cpt      = 0.3574;
  Etr[16].N        = 23;
  Etr[16].J        = 23;
  Etr[16].Np       = 21;
  Etr[16].Jp       = 21;
  Etr[16].E        = 793.2100;
  Etr[16].delta_E  = -128.9490;
  Etr[16].cpt      = 0.3574;
  Etr[17].N        = 23;
  Etr[17].J        = 24;
  Etr[17].Np       = 21;
  Etr[17].Jp       = 22;
  Etr[17].E        = 791.0344;
  Etr[17].delta_E  = -128.9314;
  Etr[17].cpt      = 0.3589;
  Etr[18].N        = 21;
  Etr[18].J        = 20;
  Etr[18].Np       = 19;
  Etr[18].Jp       = 18;
  Etr[18].E        = 662.4368;
  Etr[18].delta_E  = -117.5740;
  Etr[18].cpt      = 0.3556;
  Etr[19].N        = 21;
  Etr[19].J        = 21;
  Etr[19].Np       = 19;
  Etr[19].Jp       = 19;
  Etr[19].E        = 664.2610;
  Etr[19].delta_E  = -117.5560;
  Etr[19].cpt      = 0.3556;
  Etr[20].N        = 21;
  Etr[20].J        = 22;
  Etr[20].Np       = 19;
  Etr[20].Jp       = 20;
  Etr[20].E        = 662.1030;
  Etr[20].delta_E  = -117.5372;
  Etr[20].cpt      = 0.3573;
  Etr[21].N        = 19;
  Etr[21].J        = 19;
  Etr[21].Np       = 17;
  Etr[21].Jp       = 18;
  Etr[21].E        = 546.7050;
  Etr[21].delta_E  = -108.2632;
  Etr[21].cpt      = 0.0021;
  Etr[22].N        = 19;
  Etr[22].J        = 18;
  Etr[22].Np       = 17;
  Etr[22].Jp       = 16;
  Etr[22].E        = 544.8628;
  Etr[22].delta_E  = -106.1613;
  Etr[22].cpt      = 0.3533;
  Etr[23].N        = 19;
  Etr[23].J        = 19;
  Etr[23].Np       = 17;
  Etr[23].Jp       = 17;
  Etr[23].E        = 546.7050;
  Etr[23].delta_E  = -106.1430;
  Etr[23].cpt      = 0.3534;
  Etr[24].N        = 19;
  Etr[24].J        = 20;
  Etr[24].Np       = 17;
  Etr[24].Jp       = 18;
  Etr[24].E        = 544.5658;
  Etr[24].delta_E  = -106.1240;
  Etr[24].cpt      = 0.3555;
  Etr[25].N        = 19;
  Etr[25].J        = 18;
  Etr[25].Np       = 17;
  Etr[25].Jp       = 17;
  Etr[25].E        = 544.8628;
  Etr[25].delta_E  = -104.3008;
  Etr[25].cpt      = 0.0022;
  Etr[26].N        = 17;
  Etr[26].J        = 17;
  Etr[26].Np       = 15;
  Etr[26].Jp       = 16;
  Etr[26].E        = 440.5620;
  Etr[26].delta_E  = -96.8136;
  Etr[26].cpt      = 0.0026;
  Etr[27].N        = 17;
  Etr[27].J        = 16;
  Etr[27].Np       = 15;
  Etr[27].Jp       = 14;
  Etr[27].E        = 438.7015;
  Etr[27].delta_E  = -94.7318;
  Etr[27].cpt      = 0.3505;
  Etr[28].N        = 17;
  Etr[28].J        = 17;
  Etr[28].Np       = 15;
  Etr[28].Jp       = 15;
  Etr[28].E        = 440.5620;
  Etr[28].delta_E  = -94.7120;
  Etr[28].cpt      = 0.3506;
  Etr[29].N        = 17;
  Etr[29].J        = 18;
  Etr[29].Np       = 15;
  Etr[29].Jp       = 16;
  Etr[29].E        = 438.4418;
  Etr[29].delta_E  = -94.6934;
  Etr[29].cpt      = 0.3532;
  Etr[30].N        = 17;
  Etr[30].J        = 16;
  Etr[30].Np       = 15;
  Etr[30].Jp       = 15;
  Etr[30].E        = 438.7015;
  Etr[30].delta_E  = -92.8515;
  Etr[30].cpt      = 0.0028;
  Etr[31].N        = 15;
  Etr[31].J        = 15;
  Etr[31].Np       = 13;
  Etr[31].Jp       = 14;
  Etr[31].E        = 345.8500;
  Etr[31].delta_E  = -85.3489;
  Etr[31].cpt      = 0.0033;
  Etr[32].N        = 15;
  Etr[32].J        = 14;
  Etr[32].Np       = 13;
  Etr[32].Jp       = 12;
  Etr[32].E        = 343.9697;
  Etr[32].delta_E  = -83.2871;
  Etr[32].cpt      = 0.3468;
  Etr[33].N        = 15;
  Etr[33].J        = 15;
  Etr[33].Np       = 13;
  Etr[33].Jp       = 13;
  Etr[33].E        = 345.8500;
  Etr[33].delta_E  = -83.2671;
  Etr[33].cpt      = 0.3471;
  Etr[34].N        = 15;
  Etr[34].J        = 16;
  Etr[34].Np       = 13;
  Etr[34].Jp       = 14;
  Etr[34].E        = 343.7484;
  Etr[34].delta_E  = -83.2473;
  Etr[34].cpt      = 0.3504;
  Etr[35].N        = 15;
  Etr[35].J        = 14;
  Etr[35].Np       = 13;
  Etr[35].Jp       = 13;
  Etr[35].E        = 343.9697;
  Etr[35].delta_E  = -81.3868;
  Etr[35].cpt      = 0.0036;
  Etr[36].N        = 13;
  Etr[36].J        = 13;
  Etr[36].Np       = 11;
  Etr[36].Jp       = 12;
  Etr[36].E        = 262.5829;
  Etr[36].delta_E  = -73.8694;
  Etr[36].cpt      = 0.0044;
  Etr[37].N        = 13;
  Etr[37].J        = 12;
  Etr[37].Np       = 11;
  Etr[37].Jp       = 10;
  Etr[37].E        = 260.6826;
  Etr[37].delta_E  = -71.8294;
  Etr[37].cpt      = 0.3418;
  Etr[38].N        = 13;
  Etr[38].J        = 13;
  Etr[38].Np       = 11;
  Etr[38].Jp       = 11;
  Etr[38].E        = 262.5829;
  Etr[38].delta_E  = -71.8080;
  Etr[38].cpt      = 0.3422;
  Etr[39].N        = 13;
  Etr[39].J        = 14;
  Etr[39].Np       = 11;
  Etr[39].Jp       = 12;
  Etr[39].E        = 260.5011;
  Etr[39].delta_E  = -71.7876;
  Etr[39].cpt      = 0.3467;
  Etr[40].N        = 13;
  Etr[40].J        = 12;
  Etr[40].Np       = 11;
  Etr[40].Jp       = 11;
  Etr[40].E        = 260.6826;
  Etr[40].delta_E  = -69.9077;
  Etr[40].cpt      = 0.0048;
  Etr[41].N        = 11;
  Etr[41].J        = 11;
  Etr[41].Np       = 9;
  Etr[41].Jp       = 10;
  Etr[41].E        = 190.7749;
  Etr[41].delta_E  = -62.3771;
  Etr[41].cpt      = 0.0062;
  Etr[42].N        = 11;
  Etr[42].J        = 10;
  Etr[42].Np       = 9;
  Etr[42].Jp       = 9;
  Etr[42].E        = 188.8532;
  Etr[42].delta_E  = -60.3611;
  Etr[42].cpt      = 0.3348;
  Etr[43].N        = 11;
  Etr[43].J        = 11;
  Etr[43].Np       = 9;
  Etr[43].Jp       = 9;
  Etr[43].E        = 190.7749;
  Etr[43].delta_E  = -60.3373;
  Etr[43].cpt      = 0.3354;
  Etr[44].N        = 11;
  Etr[44].J        = 12;
  Etr[44].Np       = 9;
  Etr[44].Jp       = 10;
  Etr[44].E        = 188.7135;
  Etr[44].delta_E  = -60.3157;
  Etr[44].cpt      = 0.3416;
  Etr[45].N        = 11;
  Etr[45].J        = 10;
  Etr[45].Np       = 9;
  Etr[45].Jp       = 9;
  Etr[45].E        = 188.8532;
  Etr[45].delta_E  = -58.4156;
  Etr[45].cpt      = 0.0068;
  Etr[46].N        = 9;
  Etr[46].J        = 9;
  Etr[46].Np       = 7;
  Etr[46].Jp       = 8;
  Etr[46].E        = 130.4376;
  Etr[46].delta_E  = -50.8730;
  Etr[46].cpt      = 0.0093;
  Etr[47].N        = 9;
  Etr[47].J        = 8;
  Etr[47].Np       = 7;
  Etr[47].Jp       = 6;
  Etr[47].E        = 128.4921;
  Etr[47].delta_E  = -48.8851;
  Etr[47].cpt      = 0.3220;
  Etr[48].N        = 9;
  Etr[48].J        = 9;
  Etr[48].Np       = 7;
  Etr[48].Jp       = 7;
  Etr[48].E        = 130.4376;
  Etr[48].delta_E  = -48.8571;
  Etr[48].cpt      = 0.3251;
  Etr[49].N        = 9;
  Etr[49].J        = 10;
  Etr[49].Np       = 7;
  Etr[49].Jp       = 8;
  Etr[49].E        = 128.3978;
  Etr[49].delta_E  = -48.8332;
  Etr[49].cpt      = 0.3344;
  Etr[50].N        = 9;
  Etr[50].J        = 8;
  Etr[50].Np       = 7;
  Etr[50].Jp       = 7;
  Etr[50].E        = 128.4921;
  Etr[50].delta_E  = -46.9116;
  Etr[50].cpt      = 0.0113;
  Etr[51].N        = 5;
  Etr[51].J        = 4;
  Etr[51].Np       = 1;
  Etr[51].Jp       = 2;
  Etr[51].E        = 42.2001;
  Etr[51].delta_E  = -40.1158;
  Etr[51].cpt      = 0.0011;
  Etr[52].N        = 7;
  Etr[52].J        = 7;
  Etr[52].Np       = 5;
  Etr[52].Jp       = 6;
  Etr[52].E        = 81.5805;
  Etr[52].delta_E  = -39.3565;
  Etr[52].cpt      = 0.0139;
  Etr[53].N        = 7;
  Etr[53].J        = 6;
  Etr[53].Np       = 5;
  Etr[53].Jp       = 4;
  Etr[53].E        = 79.6070;
  Etr[53].delta_E  = -37.4069;
  Etr[53].cpt      = 0.3013;
  Etr[54].N        = 7;
  Etr[54].J        = 7;
  Etr[54].Np       = 5;
  Etr[54].Jp       = 5;
  Etr[54].E        = 81.5805;
  Etr[54].delta_E  = -37.3688;
  Etr[54].cpt      = 0.3077;
  Etr[55].N        = 7;
  Etr[55].J        = 8;
  Etr[55].Np       = 5;
  Etr[55].Jp       = 6;
  Etr[55].E        = 79.5646;
  Etr[55].delta_E  = -37.3406;
  Etr[55].cpt      = 0.3223;
  Etr[56].N        = 7;
  Etr[56].J        = 6;
  Etr[56].Np       = 5;
  Etr[56].Jp       = 5;
  Etr[56].E        = 79.6070;
  Etr[56].delta_E  = -35.3953;
  Etr[56].cpt      = 0.0198;
  Etr[57].N        = 5;
  Etr[57].J        = 5;
  Etr[57].Np       = 3;
  Etr[57].Jp       = 4;
  Etr[57].E        = 44.2117;
  Etr[57].delta_E  = -27.8241;
  Etr[57].cpt      = 0.0261;
  Etr[58].N        = 5;
  Etr[58].J        = 4;
  Etr[58].Np       = 3;
  Etr[58].Jp       = 2;
  Etr[58].E        = 42.2001;
  Etr[58].delta_E  = -25.9472;
  Etr[58].cpt      = 0.2544;
  Etr[59].N        = 5;
  Etr[59].J        = 5;
  Etr[59].Np       = 3;
  Etr[59].Jp       = 3;
  Etr[59].E        = 44.2117;
  Etr[59].delta_E  = -25.8745;
  Etr[59].cpt      = 0.2727;
  Etr[60].N        = 5;
  Etr[60].J        = 6;
  Etr[60].Np       = 3;
  Etr[60].Jp       = 4;
  Etr[60].E        = 42.2240;
  Etr[60].delta_E  = -25.8364;
  Etr[60].cpt      = 0.3020;
  Etr[61].N        = 5;
  Etr[61].J        = 4;
  Etr[61].Np       = 3;
  Etr[61].Jp       = 4;
  Etr[61].E        = 42.2001;
  Etr[61].delta_E  = -25.8125;
  Etr[61].cpt      = 0.0015;
  Etr[62].N        = 5;
  Etr[62].J        = 4;
  Etr[62].Np       = 3;
  Etr[62].Jp       = 3;
  Etr[62].E        = 42.2001;
  Etr[62].delta_E  = -23.8629;
  Etr[62].cpt      = 0.0434;
  Etr[63].N        = 3;
  Etr[63].J        = 2;
  Etr[63].Np       = 1;
  Etr[63].Jp       = 0;
  Etr[63].E        = 16.2529;
  Etr[63].delta_E  = -16.2529;
  Etr[63].cpt      = 0.0923;
  Etr[64].N        = 3;
  Etr[64].J        = 3;
  Etr[64].Np       = 1;
  Etr[64].Jp       = 2;
  Etr[64].E        = 18.3372;
  Etr[64].delta_E  = -16.2529;
  Etr[64].cpt      = 0.0660;
  Etr[65].N        = 3;
  Etr[65].J        = 3;
  Etr[65].Np       = 1;
  Etr[65].Jp       = 1;
  Etr[65].E        = 18.3372;
  Etr[65].delta_E  = -14.3761;
  Etr[65].cpt      = 0.1714;
  Etr[66].N        = 3;
  Etr[66].J        = 4;
  Etr[66].Np       = 1;
  Etr[66].Jp       = 2;
  Etr[66].E        = 16.3876;
  Etr[66].delta_E  = -14.3033;
  Etr[66].cpt      = 0.2571;
  Etr[67].N        = 3;
  Etr[67].J        = 2;
  Etr[67].Np       = 1;
  Etr[67].Jp       = 2;
  Etr[67].E        = 16.2529;
  Etr[67].delta_E  = -14.1686;
  Etr[67].cpt      = 0.0184;
  Etr[68].N        = 3;
  Etr[68].J        = 2;
  Etr[68].Np       = 1;
  Etr[68].Jp       = 1;
  Etr[68].E        = 16.2529;
  Etr[68].delta_E  = -12.2918;
  Etr[68].cpt      = 0.1615;
  Etr[69].N        = 19;
  Etr[69].J        = 19;
  Etr[69].Np       = 19;
  Etr[69].Jp       = 20;
  Etr[69].E        = 546.7050;
  Etr[69].delta_E  = -2.1392;
  Etr[69].cpt      = 0.0020;
  Etr[70].N        = 17;
  Etr[70].J        = 17;
  Etr[70].Np       = 17;
  Etr[70].Jp       = 18;
  Etr[70].E        = 440.5620;
  Etr[70].delta_E  = -2.1202;
  Etr[70].cpt      = 0.0024;
  Etr[71].N        = 15;
  Etr[71].J        = 15;
  Etr[71].Np       = 15;
  Etr[71].Jp       = 16;
  Etr[71].E        = 345.8500;
  Etr[71].delta_E  = -2.1016;
  Etr[71].cpt      = 0.0031;
  Etr[72].N        = 1;
  Etr[72].J        = 2;
  Etr[72].Np       = 1;
  Etr[72].Jp       = 0;
  Etr[72].E        = 2.0843;
  Etr[72].delta_E  = -2.0843;
  Etr[72].cpt      = 0.1077;
  Etr[73].N        = 3;
  Etr[73].J        = 3;
  Etr[73].Np       = 3;
  Etr[73].Jp       = 2;
  Etr[73].E        = 18.3372;
  Etr[73].delta_E  = -2.0843;
  Etr[73].cpt      = 0.0769;
  Etr[74].N        = 13;
  Etr[74].J        = 13;
  Etr[74].Np       = 13;
  Etr[74].Jp       = 14;
  Etr[74].E        = 262.5829;
  Etr[74].delta_E  = -2.0818;
  Etr[74].cpt      = 0.0041;
  Etr[75].N        = 11;
  Etr[75].J        = 11;
  Etr[75].Np       = 11;
  Etr[75].Jp       = 12;
  Etr[75].E        = 190.7749;
  Etr[75].delta_E  = -2.0614;
  Etr[75].cpt      = 0.0057;
  Etr[76].N        = 9;
  Etr[76].J        = 9;
  Etr[76].Np       = 9;
  Etr[76].Jp       = 10;
  Etr[76].E        = 130.4376;
  Etr[76].delta_E  = -2.0398;
  Etr[76].cpt      = 0.0083;
  Etr[77].N        = 7;
  Etr[77].J        = 7;
  Etr[77].Np       = 7;
  Etr[77].Jp       = 8;
  Etr[77].E        = 81.5805;
  Etr[77].delta_E  = -2.0159;
  Etr[77].cpt      = 0.0122;
  Etr[78].N        = 5;
  Etr[78].J        = 5;
  Etr[78].Np       = 5;
  Etr[78].Jp       = 4;
  Etr[78].E        = 44.2117;
  Etr[78].delta_E  = -2.0116;
  Etr[78].cpt      = 0.0284;
  Etr[79].N        = 5;
  Etr[79].J        = 5;
  Etr[79].Np       = 5;
  Etr[79].Jp       = 6;
  Etr[79].E        = 44.2117;
  Etr[79].delta_E  = -1.9877;
  Etr[79].cpt      = 0.0221;
  Etr[80].N        = 7;
  Etr[80].J        = 7;
  Etr[80].Np       = 7;
  Etr[80].Jp       = 6;
  Etr[80].E        = 81.5805;
  Etr[80].delta_E  = -1.9735;
  Etr[80].cpt      = 0.0147;
  Etr[81].N        = 3;
  Etr[81].J        = 3;
  Etr[81].Np       = 3;
  Etr[81].Jp       = 4;
  Etr[81].E        = 18.3372;
  Etr[81].delta_E  = -1.9496;
  Etr[81].cpt      = 0.0513;
  Etr[82].N        = 9;
  Etr[82].J        = 9;
  Etr[82].Np       = 9;
  Etr[82].Jp       = 8;
  Etr[82].E        = 130.4376;
  Etr[82].delta_E  = -1.9455;
  Etr[82].cpt      = 0.0083;
  Etr[83].N        = 11;
  Etr[83].J        = 11;
  Etr[83].Np       = 11;
  Etr[83].Jp       = 10;
  Etr[83].E        = 190.7749;
  Etr[83].delta_E  = -1.9217;
  Etr[83].cpt      = 0.0056;
  Etr[84].N        = 13;
  Etr[84].J        = 13;
  Etr[84].Np       = 13;
  Etr[84].Jp       = 12;
  Etr[84].E        = 262.5829;
  Etr[84].delta_E  = -1.9003;
  Etr[84].cpt      = 0.0041;
  Etr[85].N        = 15;
  Etr[85].J        = 15;
  Etr[85].Np       = 15;
  Etr[85].Jp       = 14;
  Etr[85].E        = 345.8500;
  Etr[85].delta_E  = -1.8803;
  Etr[85].cpt      = 0.0031;
  Etr[86].N        = 1;
  Etr[86].J        = 1;
  Etr[86].Np       = 1;
  Etr[86].Jp       = 2;
  Etr[86].E        = 3.9611;
  Etr[86].delta_E  = -1.8768;
  Etr[86].cpt      = 0.2308;
  Etr[87].N        = 17;
  Etr[87].J        = 17;
  Etr[87].Np       = 17;
  Etr[87].Jp       = 16;
  Etr[87].E        = 440.5620;
  Etr[87].delta_E  = -1.8605;
  Etr[87].cpt      = 0.0024;
  Etr[88].N        = 19;
  Etr[88].J        = 19;
  Etr[88].Np       = 19;
  Etr[88].Jp       = 18;
  Etr[88].E        = 546.7050;
  Etr[88].delta_E  = -1.8422;
  Etr[88].cpt      = 0.0020;
  Etr[89].N        = 3;
  Etr[89].J        = 4;
  Etr[89].Np       = 3;
  Etr[89].Jp       = 2;
  Etr[89].E        = 16.3876;
  Etr[89].delta_E  = -0.1347;
  Etr[89].cpt      = 0.0021;
  Etr[90].N        = 3;
  Etr[90].J        = 2;
  Etr[90].Np       = 3;
  Etr[90].Jp       = 4;
  Etr[90].E        = 16.2529;
  Etr[90].delta_E  = 0.1347;
  Etr[90].cpt      = 0.0038;
  Etr[91].N        = 19;
  Etr[91].J        = 18;
  Etr[91].Np       = 19;
  Etr[91].Jp       = 19;
  Etr[91].E        = 544.8628;
  Etr[91].delta_E  = 1.8422;
  Etr[91].cpt      = 0.0021;
  Etr[92].N        = 17;
  Etr[92].J        = 16;
  Etr[92].Np       = 17;
  Etr[92].Jp       = 17;
  Etr[92].E        = 438.7015;
  Etr[92].delta_E  = 1.8605;
  Etr[92].cpt      = 0.0026;
  Etr[93].N        = 1;
  Etr[93].J        = 2;
  Etr[93].Np       = 1;
  Etr[93].Jp       = 1;
  Etr[93].E        = 2.0843;
  Etr[93].delta_E  = 1.8768;
  Etr[93].cpt      = 0.1385;
  Etr[94].N        = 15;
  Etr[94].J        = 14;
  Etr[94].Np       = 15;
  Etr[94].Jp       = 15;
  Etr[94].E        = 343.9697;
  Etr[94].delta_E  = 1.8803;
  Etr[94].cpt      = 0.0033;
  Etr[95].N        = 13;
  Etr[95].J        = 12;
  Etr[95].Np       = 13;
  Etr[95].Jp       = 13;
  Etr[95].E        = 260.6826;
  Etr[95].delta_E  = 1.9003;
  Etr[95].cpt      = 0.0044;
  Etr[96].N        = 11;
  Etr[96].J        = 10;
  Etr[96].Np       = 11;
  Etr[96].Jp       = 11;
  Etr[96].E        = 188.8532;
  Etr[96].delta_E  = 1.9217;
  Etr[96].cpt      = 0.0062;
  Etr[97].N        = 9;
  Etr[97].J        = 8;
  Etr[97].Np       = 9;
  Etr[97].Jp       = 9;
  Etr[97].E        = 128.4921;
  Etr[97].delta_E  = 1.9455;
  Etr[97].cpt      = 0.0092;
  Etr[98].N        = 3;
  Etr[98].J        = 4;
  Etr[98].Np       = 3;
  Etr[98].Jp       = 3;
  Etr[98].E        = 16.3876;
  Etr[98].delta_E  = 1.9496;
  Etr[98].cpt      = 0.0399;
  Etr[99].N        = 7;
  Etr[99].J        = 6;
  Etr[99].Np       = 7;
  Etr[99].Jp       = 7;
  Etr[99].E        = 79.6070;
  Etr[99].delta_E  = 1.9735;
  Etr[99].cpt      = 0.0170;
  Etr[100].N       = 5;
  Etr[100].J       = 6;
  Etr[100].Np      = 5;
  Etr[100].Jp      = 5;
  Etr[100].E       = 42.2240;
  Etr[100].delta_E = 1.9877;
  Etr[100].cpt     = 0.0187;
  Etr[101].N       = 5;
  Etr[101].J       = 4;
  Etr[101].Np      = 5;
  Etr[101].Jp      = 5;
  Etr[101].E       = 42.2001;
  Etr[101].delta_E = 2.0116;
  Etr[101].cpt     = 0.0347;
  Etr[102].N       = 7;
  Etr[102].J       = 8;
  Etr[102].Np      = 7;
  Etr[102].Jp      = 7;
  Etr[102].E       = 79.5646;
  Etr[102].delta_E = 2.0159;
  Etr[102].cpt     = 0.0108;
  Etr[103].N       = 9;
  Etr[103].J       = 10;
  Etr[103].Np      = 9;
  Etr[103].Jp      = 9;
  Etr[103].E       = 128.3978;
  Etr[103].delta_E = 2.0398;
  Etr[103].cpt     = 0.0075;
  Etr[104].N       = 11;
  Etr[104].J       = 12;
  Etr[104].Np      = 11;
  Etr[104].Jp      = 11;
  Etr[104].E       = 188.7135;
  Etr[104].delta_E = 2.0614;
  Etr[104].cpt     = 0.0052;
  Etr[105].N       = 13;
  Etr[105].J       = 14;
  Etr[105].Np      = 13;
  Etr[105].Jp      = 13;
  Etr[105].E       = 260.5011;
  Etr[105].delta_E = 2.0818;
  Etr[105].cpt     = 0.0038;
  Etr[106].N       = 1;
  Etr[106].J       = 0;
  Etr[106].Np      = 1;
  Etr[106].Jp      = 2;
  Etr[106].E       = 0.0000;
  Etr[106].delta_E = 2.0843;
  Etr[106].cpt     = 0.5383;
  Etr[107].N       = 3;
  Etr[107].J       = 2;
  Etr[107].Np      = 3;
  Etr[107].Jp      = 3;
  Etr[107].E       = 16.2529;
  Etr[107].delta_E = 2.0843;
  Etr[107].cpt     = 0.1077;
  Etr[108].N       = 15;
  Etr[108].J       = 16;
  Etr[108].Np      = 15;
  Etr[108].Jp      = 15;
  Etr[108].E       = 343.7484;
  Etr[108].delta_E = 2.1016;
  Etr[108].cpt     = 0.0029;
  Etr[109].N       = 17;
  Etr[109].J       = 18;
  Etr[109].Np      = 17;
  Etr[109].Jp      = 17;
  Etr[109].E       = 438.4418;
  Etr[109].delta_E = 2.1202;
  Etr[109].cpt     = 0.0023;
  Etr[110].N       = 19;
  Etr[110].J       = 20;
  Etr[110].Np      = 19;
  Etr[110].Jp      = 19;
  Etr[110].E       = 544.5658;
  Etr[110].delta_E = 2.1392;
  Etr[110].cpt     = 0.0019;
  Etr[111].N       = 1;
  Etr[111].J       = 1;
  Etr[111].Np      = 3;
  Etr[111].Jp      = 2;
  Etr[111].E       = 3.9611;
  Etr[111].delta_E = 12.2918;
  Etr[111].cpt     = 0.2692;
  Etr[112].N       = 1;
  Etr[112].J       = 2;
  Etr[112].Np      = 3;
  Etr[112].Jp      = 2;
  Etr[112].E       = 2.0843;
  Etr[112].delta_E = 14.1686;
  Etr[112].cpt     = 0.0184;
  Etr[113].N       = 1;
  Etr[113].J       = 2;
  Etr[113].Np      = 3;
  Etr[113].Jp      = 4;
  Etr[113].E       = 2.0843;
  Etr[113].delta_E = 14.3033;
  Etr[113].cpt     = 0.4628;
  Etr[114].N       = 1;
  Etr[114].J       = 1;
  Etr[114].Np      = 3;
  Etr[114].Jp      = 3;
  Etr[114].E       = 3.9611;
  Etr[114].delta_E = 14.3761;
  Etr[114].cpt     = 0.4000;
  Etr[115].N       = 1;
  Etr[115].J       = 0;
  Etr[115].Np      = 3;
  Etr[115].Jp      = 2;
  Etr[115].E       = 0.0000;
  Etr[115].delta_E = 16.2529;
  Etr[115].cpt     = 0.4617;
  Etr[116].N       = 1;
  Etr[116].J       = 2;
  Etr[116].Np      = 3;
  Etr[116].Jp      = 3;
  Etr[116].E       = 2.0843;
  Etr[116].delta_E = 16.2529;
  Etr[116].cpt     = 0.0923;
  Etr[117].N       = 3;
  Etr[117].J       = 3;
  Etr[117].Np      = 5;
  Etr[117].Jp      = 4;
  Etr[117].E       = 18.3372;
  Etr[117].delta_E = 23.8629;
  Etr[117].cpt     = 0.0558;
  Etr[118].N       = 3;
  Etr[118].J       = 4;
  Etr[118].Np      = 5;
  Etr[118].Jp      = 4;
  Etr[118].E       = 16.3876;
  Etr[118].delta_E = 25.8125;
  Etr[118].cpt     = 0.0015;
  Etr[119].N       = 3;
  Etr[119].J       = 4;
  Etr[119].Np      = 5;
  Etr[119].Jp      = 6;
  Etr[119].E       = 16.3876;
  Etr[119].delta_E = 25.8364;
  Etr[119].cpt     = 0.4362;
  Etr[120].N       = 3;
  Etr[120].J       = 3;
  Etr[120].Np      = 5;
  Etr[120].Jp      = 5;
  Etr[120].E       = 18.3372;
  Etr[120].delta_E = 25.8745;
  Etr[120].cpt     = 0.4286;
  Etr[121].N       = 3;
  Etr[121].J       = 2;
  Etr[121].Np      = 5;
  Etr[121].Jp      = 4;
  Etr[121].E       = 16.2529;
  Etr[121].delta_E = 25.9472;
  Etr[121].cpt     = 0.4579;
  Etr[122].N       = 3;
  Etr[122].J       = 4;
  Etr[122].Np      = 5;
  Etr[122].Jp      = 5;
  Etr[122].E       = 16.3876;
  Etr[122].delta_E = 27.8241;
  Etr[122].cpt     = 0.0319;
  Etr[123].N       = 5;
  Etr[123].J       = 5;
  Etr[123].Np      = 7;
  Etr[123].Jp      = 6;
  Etr[123].E       = 44.2117;
  Etr[123].delta_E = 35.3953;
  Etr[123].cpt     = 0.0234;
  Etr[124].N       = 5;
  Etr[124].J       = 6;
  Etr[124].Np      = 7;
  Etr[124].Jp      = 8;
  Etr[124].E       = 42.2240;
  Etr[124].delta_E = 37.3406;
  Etr[124].cpt     = 0.4214;
  Etr[125].N       = 5;
  Etr[125].J       = 5;
  Etr[125].Np      = 7;
  Etr[125].Jp      = 7;
  Etr[125].E       = 44.2117;
  Etr[125].delta_E = 37.3688;
  Etr[125].cpt     = 0.4196;
  Etr[126].N       = 5;
  Etr[126].J       = 4;
  Etr[126].Np      = 7;
  Etr[126].Jp      = 6;
  Etr[126].E       = 42.2001;
  Etr[126].delta_E = 37.4069;
  Etr[126].cpt     = 0.4352;
  Etr[127].N       = 5;
  Etr[127].J       = 6;
  Etr[127].Np      = 7;
  Etr[127].Jp      = 7;
  Etr[127].E       = 42.2240;
  Etr[127].delta_E = 39.3565;
  Etr[127].cpt     = 0.0160;
  Etr[128].N       = 1;
  Etr[128].J       = 2;
  Etr[128].Np      = 5;
  Etr[128].Jp      = 4;
  Etr[128].E       = 2.0843;
  Etr[128].delta_E = 40.1158;
  Etr[128].cpt     = 0.0019;
  Etr[129].N       = 7;
  Etr[129].J       = 7;
  Etr[129].Np      = 9;
  Etr[129].Jp      = 8;
  Etr[129].E       = 81.5805;
  Etr[129].delta_E = 46.9116;
  Etr[129].cpt     = 0.0128;
  Etr[130].N       = 7;
  Etr[130].J       = 8;
  Etr[130].Np      = 9;
  Etr[130].Jp      = 10;
  Etr[130].E       = 79.5646;
  Etr[130].delta_E = 48.8332;
  Etr[130].cpt     = 0.4130;
  Etr[131].N       = 7;
  Etr[131].J       = 7;
  Etr[131].Np      = 9;
  Etr[131].Jp      = 9;
  Etr[131].E       = 81.5805;
  Etr[131].delta_E = 48.8571;
  Etr[131].cpt     = 0.4118;
  Etr[132].N       = 7;
  Etr[132].J       = 6;
  Etr[132].Np      = 9;
  Etr[132].Jp      = 8;
  Etr[132].E       = 79.6070;
  Etr[132].delta_E = 48.8851;
  Etr[132].cpt     = 0.4210;
  Etr[133].N       = 7;
  Etr[133].J       = 8;
  Etr[133].Np      = 9;
  Etr[133].Jp      = 9;
  Etr[133].E       = 79.5646;
  Etr[133].delta_E = 50.8730;
  Etr[133].cpt     = 0.0104;
  Etr[134].N       = 9;
  Etr[134].J       = 9;
  Etr[134].Np      = 11;
  Etr[134].Jp      = 10;
  Etr[134].E       = 130.4376;
  Etr[134].delta_E = 58.4156;
  Etr[134].cpt     = 0.0075;
  Etr[135].N       = 9;
  Etr[135].J       = 10;
  Etr[135].Np      = 11;
  Etr[135].Jp      = 12;
  Etr[135].E       = 128.3978;
  Etr[135].delta_E = 60.3157;
  Etr[135].cpt     = 0.4067;
  Etr[136].N       = 9;
  Etr[136].J       = 9;
  Etr[136].Np      = 11;
  Etr[136].Jp      = 11;
  Etr[136].E       = 130.4376;
  Etr[136].delta_E = 60.3373;
  Etr[136].cpt     = 0.4060;
  Etr[137].N       = 9;
  Etr[137].J       = 8;
  Etr[137].Np      = 11;
  Etr[137].Jp      = 10;
  Etr[137].E       = 128.4921;
  Etr[137].delta_E = 60.3611;
  Etr[137].cpt     = 0.4135;
  Etr[138].N       = 9;
  Etr[138].J       = 10;
  Etr[138].Np      = 11;
  Etr[138].Jp      = 11;
  Etr[138].E       = 128.3978;
  Etr[138].delta_E = 62.3771;
  Etr[138].cpt     = 0.0068;
  Etr[139].N       = 11;
  Etr[139].J       = 11;
  Etr[139].Np      = 13;
  Etr[139].Jp      = 12;
  Etr[139].E       = 190.7749;
  Etr[139].delta_E = 69.9077;
  Etr[139].cpt     = 0.0052;
  Etr[140].N       = 11;
  Etr[140].J       = 12;
  Etr[140].Np      = 13;
  Etr[140].Jp      = 14;
  Etr[140].E       = 188.7135;
  Etr[140].delta_E = 71.7876;
  Etr[140].cpt     = 0.4021;
  Etr[141].N       = 11;
  Etr[141].J       = 11;
  Etr[141].Np      = 13;
  Etr[141].Jp      = 13;
  Etr[141].E       = 190.7749;
  Etr[141].delta_E = 71.8080;
  Etr[141].cpt     = 0.4017;
  Etr[142].N       = 11;
  Etr[142].J       = 10;
  Etr[142].Np      = 13;
  Etr[142].Jp      = 12;
  Etr[142].E       = 188.8532;
  Etr[142].delta_E = 71.8294;
  Etr[142].cpt     = 0.4070;
  Etr[143].N       = 11;
  Etr[143].J       = 12;
  Etr[143].Np      = 13;
  Etr[143].Jp      = 13;
  Etr[143].E       = 188.7135;
  Etr[143].delta_E = 73.8694;
  Etr[143].cpt     = 0.0048;
  Etr[144].N       = 13;
  Etr[144].J       = 13;
  Etr[144].Np      = 15;
  Etr[144].Jp      = 14;
  Etr[144].E       = 262.5829;
  Etr[144].delta_E = 81.3868;
  Etr[144].cpt     = 0.0038;
  Etr[145].N       = 13;
  Etr[145].J       = 14;
  Etr[145].Np      = 15;
  Etr[145].Jp      = 16;
  Etr[145].E       = 260.5011;
  Etr[145].delta_E = 83.2473;
  Etr[145].cpt     = 0.3987;
  Etr[146].N       = 13;
  Etr[146].J       = 13;
  Etr[146].Np      = 15;
  Etr[146].Jp      = 15;
  Etr[146].E       = 262.5829;
  Etr[146].delta_E = 83.2671;
  Etr[146].cpt     = 0.3985;
  Etr[147].N       = 13;
  Etr[147].J       = 12;
  Etr[147].Np      = 15;
  Etr[147].Jp      = 14;
  Etr[147].E       = 260.6826;
  Etr[147].delta_E = 83.2871;
  Etr[147].cpt     = 0.4023;
  Etr[148].N       = 13;
  Etr[148].J       = 14;
  Etr[148].Np      = 15;
  Etr[148].Jp      = 15;
  Etr[148].E       = 260.5011;
  Etr[148].delta_E = 85.3489;
  Etr[148].cpt     = 0.0036;
  Etr[149].N       = 15;
  Etr[149].J       = 15;
  Etr[149].Np      = 17;
  Etr[149].Jp      = 16;
  Etr[149].E       = 345.8500;
  Etr[149].delta_E = 92.8515;
  Etr[149].cpt     = 0.0029;
  Etr[150].N       = 15;
  Etr[150].J       = 16;
  Etr[150].Np      = 17;
  Etr[150].Jp      = 18;
  Etr[150].E       = 343.7484;
  Etr[150].delta_E = 94.6934;
  Etr[150].cpt     = 0.3961;
  Etr[151].N       = 15;
  Etr[151].J       = 15;
  Etr[151].Np      = 17;
  Etr[151].Jp      = 17;
  Etr[151].E       = 345.8500;
  Etr[151].delta_E = 94.7120;
  Etr[151].cpt     = 0.3959;
  Etr[152].N       = 15;
  Etr[152].J       = 14;
  Etr[152].Np      = 17;
  Etr[152].Jp      = 16;
  Etr[152].E       = 343.9697;
  Etr[152].delta_E = 94.7318;
  Etr[152].cpt     = 0.3988;
  Etr[153].N       = 15;
  Etr[153].J       = 16;
  Etr[153].Np      = 17;
  Etr[153].Jp      = 17;
  Etr[153].E       = 343.7484;
  Etr[153].delta_E = 96.8136;
  Etr[153].cpt     = 0.0028;
  Etr[154].N       = 17;
  Etr[154].J       = 17;
  Etr[154].Np      = 19;
  Etr[154].Jp      = 18;
  Etr[154].E       = 440.5620;
  Etr[154].delta_E = 104.3008;
  Etr[154].cpt     = 0.0023;
  Etr[155].N       = 17;
  Etr[155].J       = 18;
  Etr[155].Np      = 19;
  Etr[155].Jp      = 20;
  Etr[155].E       = 438.4418;
  Etr[155].delta_E = 106.1240;
  Etr[155].cpt     = 0.3939;
  Etr[156].N       = 17;
  Etr[156].J       = 17;
  Etr[156].Np      = 19;
  Etr[156].Jp      = 19;
  Etr[156].E       = 440.5620;
  Etr[156].delta_E = 106.1430;
  Etr[156].cpt     = 0.3938;
  Etr[157].N       = 17;
  Etr[157].J       = 16;
  Etr[157].Np      = 19;
  Etr[157].Jp      = 18;
  Etr[157].E       = 438.7015;
  Etr[157].delta_E = 106.1613;
  Etr[157].cpt     = 0.3961;
  Etr[158].N       = 17;
  Etr[158].J       = 18;
  Etr[158].Np      = 19;
  Etr[158].Jp      = 19;
  Etr[158].E       = 438.4418;
  Etr[158].delta_E = 108.2632;
  Etr[158].cpt     = 0.0022;
  Etr[159].N       = 19;
  Etr[159].J       = 19;
  Etr[159].Np      = 21;
  Etr[159].Jp      = 20;
  Etr[159].E       = 546.7050;
  Etr[159].delta_E = 115.7318;
  Etr[159].cpt     = 0.0019;
  Etr[160].N       = 19;
  Etr[160].J       = 20;
  Etr[160].Np      = 21;
  Etr[160].Jp      = 22;
  Etr[160].E       = 544.5658;
  Etr[160].delta_E = 117.5372;
  Etr[160].cpt     = 0.3922;
  Etr[161].N       = 19;
  Etr[161].J       = 19;
  Etr[161].Np      = 21;
  Etr[161].Jp      = 21;
  Etr[161].E       = 546.7050;
  Etr[161].delta_E = 117.5560;
  Etr[161].cpt     = 0.3921;
  Etr[162].N       = 19;
  Etr[162].J       = 18;
  Etr[162].Np      = 21;
  Etr[162].Jp      = 20;
  Etr[162].E       = 544.8628;
  Etr[162].delta_E = 117.5740;
  Etr[162].cpt     = 0.3940;
  Etr[163].N       = 19;
  Etr[163].J       = 20;
  Etr[163].Np      = 21;
  Etr[163].Jp      = 21;
  Etr[163].E       = 544.5658;
  Etr[163].delta_E = 119.6952;
  Etr[163].cpt     = 0.0018;
  Etr[164].N       = 21;
  Etr[164].J       = 22;
  Etr[164].Np      = 23;
  Etr[164].Jp      = 24;
  Etr[164].E       = 662.1030;
  Etr[164].delta_E = 128.9314;
  Etr[164].cpt     = 0.3908;
  Etr[165].N       = 21;
  Etr[165].J       = 21;
  Etr[165].Np      = 23;
  Etr[165].Jp      = 23;
  Etr[165].E       = 664.2610;
  Etr[165].delta_E = 128.9490;
  Etr[165].cpt     = 0.3907;
  Etr[166].N       = 21;
  Etr[166].J       = 20;
  Etr[166].Np      = 23;
  Etr[166].Jp      = 22;
  Etr[166].E       = 662.4368;
  Etr[166].delta_E = 128.9677;
  Etr[166].cpt     = 0.3922;
  Etr[167].N       = 23;
  Etr[167].J       = 24;
  Etr[167].Np      = 25;
  Etr[167].Jp      = 26;
  Etr[167].E       = 791.0344;
  Etr[167].delta_E = 140.3046;
  Etr[167].cpt     = 0.3895;
  Etr[168].N       = 23;
  Etr[168].J       = 23;
  Etr[168].Np      = 25;
  Etr[168].Jp      = 25;
  Etr[168].E       = 793.2100;
  Etr[168].delta_E = 140.3230;
  Etr[168].cpt     = 0.3895;
  Etr[169].N       = 23;
  Etr[169].J       = 22;
  Etr[169].Np      = 25;
  Etr[169].Jp      = 24;
  Etr[169].E       = 791.4045;
  Etr[169].delta_E = 140.3405;
  Etr[169].cpt     = 0.3908;
  Etr[170].N       = 25;
  Etr[170].J       = 26;
  Etr[170].Np      = 27;
  Etr[170].Jp      = 28;
  Etr[170].E       = 931.3390;
  Etr[170].delta_E = 151.6551;
  Etr[170].cpt     = 0.3885;
  Etr[171].N       = 25;
  Etr[171].J       = 25;
  Etr[171].Np      = 27;
  Etr[171].Jp      = 27;
  Etr[171].E       = 933.5330;
  Etr[171].delta_E = 151.6730;
  Etr[171].cpt     = 0.3885;
  Etr[172].N       = 25;
  Etr[172].J       = 24;
  Etr[172].Np      = 27;
  Etr[172].Jp      = 26;
  Etr[172].E       = 931.7450;
  Etr[172].delta_E = 151.6906;
  Etr[172].cpt     = 0.3896;
  Etr[173].N       = 27;
  Etr[173].J       = 28;
  Etr[173].Np      = 29;
  Etr[173].Jp      = 30;
  Etr[173].E       = 1082.9941;
  Etr[173].delta_E = 162.9809;
  Etr[173].cpt     = 0.3876;
  Etr[174].N       = 27;
  Etr[174].J       = 27;
  Etr[174].Np      = 29;
  Etr[174].Jp      = 29;
  Etr[174].E       = 1085.2061;
  Etr[174].delta_E = 162.9980;
  Etr[174].cpt     = 0.3876;
  Etr[175].N       = 27;
  Etr[175].J       = 26;
  Etr[175].Np      = 29;
  Etr[175].Jp      = 28;
  Etr[175].E       = 1083.4355;
  Etr[175].delta_E = 163.0162;
  Etr[175].cpt     = 0.3885;
  Etr[176].N       = 29;
  Etr[176].J       = 30;
  Etr[176].Np      = 31;
  Etr[176].Jp      = 32;
  Etr[176].E       = 1245.9750;
  Etr[176].delta_E = 174.2802;
  Etr[176].cpt     = 0.3868;
  Etr[177].N       = 29;
  Etr[177].J       = 29;
  Etr[177].Np      = 31;
  Etr[177].Jp      = 31;
  Etr[177].E       = 1248.2040;
  Etr[177].delta_E = 174.2980;
  Etr[177].cpt     = 0.3868;
  Etr[178].N       = 29;
  Etr[178].J       = 28;
  Etr[178].Np      = 31;
  Etr[178].Jp      = 30;
  Etr[178].E       = 1246.4518;
  Etr[178].delta_E = 174.3154;
  Etr[178].cpt     = 0.3876;
  Etr[179].N       = 31;
  Etr[179].J       = 32;
  Etr[179].Np      = 33;
  Etr[179].Jp      = 34;
  Etr[179].E       = 1420.2552;
  Etr[179].delta_E = 185.5512;
  Etr[179].cpt     = 0.3861;
  Etr[180].N       = 31;
  Etr[180].J       = 31;
  Etr[180].Np      = 33;
  Etr[180].Jp      = 33;
  Etr[180].E       = 1422.5020;
  Etr[180].delta_E = 185.5690;
  Etr[180].cpt     = 0.3861;
  Etr[181].N       = 31;
  Etr[181].J       = 30;
  Etr[181].Np      = 33;
  Etr[181].Jp      = 32;
  Etr[181].E       = 1420.7672;
  Etr[181].delta_E = 185.5861;
  Etr[181].cpt     = 0.3868;
  Etr[182].N       = 33;
  Etr[182].J       = 34;
  Etr[182].Np      = 35;
  Etr[182].Jp      = 36;
  Etr[182].E       = 1605.8064;
  Etr[182].delta_E = 196.7919;
  Etr[182].cpt     = 0.3855;
  Etr[183].N       = 33;
  Etr[183].J       = 33;
  Etr[183].Np      = 35;
  Etr[183].Jp      = 35;
  Etr[183].E       = 1608.0710;
  Etr[183].delta_E = 196.8100;
  Etr[183].cpt     = 0.3855;
  Etr[184].N       = 33;
  Etr[184].J       = 32;
  Etr[184].Np      = 35;
  Etr[184].Jp      = 34;
  Etr[184].E       = 1606.3533;
  Etr[184].delta_E = 196.8269;
  Etr[184].cpt     = 0.3861;

  double E_J           = 0; /* Rotational energy of state J */
  double W_J           = 0; /* Fraction of molecules in the rotational state J at temperature T */
  double b_J           = 0; /* Placzek-Teller coefficient */
  double lambda_inv_cm = 0, lambda_cm = 0, lambda_shifted_inv_cm = 0, lambda_shifted_nm = 0, lambda_shifted_cm = 0;
  double delta_nu = 0;
  double crs_i = 0, crs_j = 0;
  double fNJ = 0, Z = 0; /* Eq 5 and 10 */
  double gam_i = 0, gam_j = 0;
  double volume_mixing_ratio = 0.2095;
  double sum                 = 0;

  int J   = 0; /* Rotational state J */
  int g_J = 0;
  int is  = 0;

  int status = 0;
  int ivi    = 0;

  if (verbose) {
    fprintf (stderr, "Calculating Raman scattering wavelength shifts and cross sections ");
    fprintf (stderr,
             "for O2, wavelength %5.1f nm, number of transitions %d, T=%6.2f, iv=%d.\n",
             lambda,
             n_transitions,
             temper,
             *iv);
    fprintf (stderr,
             "             %3s %1s %7s %8s %7s %9s %10s %9s %17s %9s %12s %7s %12s %12s %11s %14s %12s\n",
             "ivi",
             "J",
             "E_J",
             "g_J",
             "W_J",
             "b_J",
             "W_J*b_J",
             "wvl",
             "wvl_shift_cm-1",
             "delta_nu",
             "wvl_shift_nm",
             "gam_i",
             "gam_j",
             "crs_i",
             "crs_j",
             "crs_i*vmr",
             "crs_j*vmr");
  }

  ivi = *iv;

  sum = 0;
  for (is = 0; is < N_E; is++) {
    g_J = 1; /* Only states for odd J are included in the summation */
    J   = E_rot[is].J;
    E_J = E_rot[is].E * hc;
    sum += g_J * (2 * J + 1) * exp (-E_J / (BOLTZMANN * temper));
  }
  Z = sum;
  /*   fprintf(stderr," %13.6e", Z); */

  for (is = 0; is < n_transitions; is++) {
    g_J = 1; /* Only states for odd J are included in the summation */
    J   = Etr[is].J;
    E_J = Etr[is].E * hc;
    W_J = g_J * (2 * J + 1) * exp (-E_J / (BOLTZMANN * temper));
    fNJ = W_J / Z;
    b_J = Etr[is].cpt;

    delta_nu = -Etr[is].delta_E; /* Negative delta_E corresponds to a photon with a larger */
                                 /* energy, shorter wavelength than the incident photon    */
                                 /* delta_E is the change in rotational energy of the      */
                                 /* molecule. The minus puts the photon in the right       */
                                 /* wavelength.                                            */

    lambda_inv_cm         = nm_to_inv_cm (lambda);
    lambda_shifted_inv_cm = lambda_inv_cm + delta_nu;
    lambda_shifted_nm     = inv_cm_to_nm (lambda_shifted_inv_cm);
    lambda_shifted_cm     = lambda_shifted_nm * nm_to_cm;
    lambda_cm             = lambda * nm_to_cm;
    gam_i                 = polarizability_anisotropy_O2 (lambda_shifted_inv_cm);
    gam_j                 = polarizability_anisotropy_O2 (lambda_inv_cm);
    crs_i                 = fNJ * gam_i * gam_i * raman_const * b_J / pow (lambda_shifted_cm, 4);
    crs_j                 = fNJ * gam_j * gam_j * raman_const * b_J / pow (lambda_cm, 4);
    (*crs)[ivi][0]        = lambda_shifted_nm;
    (*crs)[ivi][1]        = crs_i * volume_mixing_ratio;
    (*crs)[ivi][2]        = crs_j * volume_mixing_ratio;
    if (verbose)
      fprintf (stderr,
               "Raman_crs_O2 %2d %2d %12.6e %2d %12.6e %6.3f %12.6e %10.6f %10.6f %10.4f %10.6f %12.6e %12.6e %12.6e %12.6e %12.6e "
               "%12.6e\n",
               ivi,
               J,
               E_J,
               g_J,
               W_J,
               b_J,
               W_J * b_J,
               lambda,
               lambda_shifted_inv_cm,
               delta_nu,
               lambda_shifted_nm,
               gam_i,
               gam_j,
               crs_i,
               crs_j,
               crs_i * volume_mixing_ratio,
               crs_j * volume_mixing_ratio);
    ivi++;
  }

  *iv = ivi;

  return status;
}

/***********************************************************************************/
/* Function: dewpoint                                                              */
/* Description:                                                                    */
/*  Calculates the dewpoint temperature in K                                       */
/*                                                                                 */
/* Parameters:                                                                     */
/*     float press_h2o    partial pressure of water vapour in hPa                  */
/* Return value:                                                                   */
/*     float dewpoint     dewpoint temperature in K                                */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    ancillary.c                                                           */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Jan 2009   U. Hamann     converted from Fortran to C (see wvapour.f)         */
/*                                                                                 */
/***********************************************************************************/

float dewpoint (float press_h2o) {
  float dewpoint = NOT_DEFINED_FLOAT;

  if (press_h2o != 0.0) {

    press_h2o = log (press_h2o);
    dewpoint  = 273.15 + (243.5 * press_h2o - 440.8) / (19.48 - press_h2o);
  } else {
    dewpoint = 0.0; /* 0 Kelvin */
  }
  return dewpoint;
}

/***********************************************************************************/
/* Function: Tcon                                                                  */
/* Description:                                                                    */
/*  THIS FUNCTION RETURNS THE TEMPERATURE TCON (CELSIUS) AT                        */
/*   THE LIFTING CONDENSATION LEVEL, GIVEN THE TEMPERATURE T (CELSIUS)             */
/*   AND THE DEW POINT D (CELSIUS).                                                */
/*                                                                                 */
/* Parameters:                                                                     */
/*     T          - REAL TEMPERATURE (K)                                           */
/*     TD         - REAL DEWPOINT TEMPERATURE (K)                                  */
/* Return value:                                                                   */
/*     float Tcon   temperature at the lifting condensation level (CELSIUS)        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    ancillary.c                                                           */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    May 1982  D. Baker, T. Schlatter   original version                          */
/*    Jan 2009  U. Hamann  converted from Fortran to C (see wvapour.f)             */
/*                         change Celsius to Kelvin                                */
/*                                                                                 */
/***********************************************************************************/

float Tcon (float T, float TD) {
  float S    = NOT_DEFINED_FLOAT;
  float dT   = NOT_DEFINED_FLOAT;
  float Tcon = NOT_DEFINED_FLOAT;

  /* compute the dew point depression S. */
  S = T - TD;
  /* the approximation below, a third order polynomial in S and T, */
  /* is due to herman wobus. the source of data for fitting the    */
  /* polynomial is unknown.                                        */

  dT   = S * (1.2185 + 1.278e-3 * (T - 273.15) + S * (-2.19e-3 + 1.173e-5 * S - 5.2e-6 * (T - 273.15)));
  Tcon = T - dT;
  return Tcon;
}

/***********************************************************************************/
/* Function: EPT                                                                   */
/* Description:                                                                    */
/*    This function returns the equivalent potential temperature EPT               */
/*    (Celsius) for a parcel of air initially at temperature t (celsius),          */
/*    dew point Td (celsius) and pressure p (millibars). The formula used          */
/*    is eq.(43) in Bolton, David, 1980: "the computation of equivalent            */
/*    potential temperature," Monthly Weather Review, vol. 108, no. 7              */
/*    (july), pp. 1046-1053. the maximum error in ept in 0.3c.  in most            */
/*    cases the error is less than 0.1c.                                           */
/*                                                                                 */
/* Parameters:                                                                     */
/*     T        - REAL TEMPERATURE (K)                                             */
/*     TD       - REAL DEWPOINT TEMPERATURE (K)                                    */
/*     p        - REAL PRESSURE (hPa)                                              */
/* Return value:                                                                   */
/*     float EPT   equivalent potential temperature (CELSIUS)                      */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    ancillary.c                                                           */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    May 1982  T. Schlatter  original version                                     */
/*    Jan 2009  U. Hamann     converted from Fortran to C (see wvapour.f)          */
/*                            replaced function WMR by real calculaton             */
/*                            change Celsius to Kelvin                             */
/*                                                                                 */
/***********************************************************************************/

float EPT (float T, float TD, float p, float N_AIR, float N_H2O) {

  float kappa = (C_P_DRY_STD - C_V_DRY_STD) / C_P_DRY_STD;
  float W     = NOT_DEFINED_FLOAT;
  float TL    = NOT_DEFINED_FLOAT;
  float PT    = NOT_DEFINED_FLOAT;
  float EPT   = NOT_DEFINED_FLOAT;

  /*  replaced the wmr (compute the water mixing ratio) function as */
  /*  it is not valid for all wanted temperatures and pressures !!! */
  /*  1000 == kg(water)/kg(dry air) -> g(water)/kg(dry air)         */
  W = N_H2O * MOL_MASS_WV / (N_AIR * MOL_MASS_AIR) * 1000.;

  /* compute the temperature (celsius) at the lifting condensation level. */
  TL  = Tcon (T, TD);
  PT  = T * pow (1000. / p, kappa * (1. - 0.00028 * W));
  EPT = PT * exp ((3.376 / TL - 0.00254) * W * (1. + 0.00081 * W));

  return EPT;
}
