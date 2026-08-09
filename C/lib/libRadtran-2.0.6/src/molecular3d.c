/*--------------------------------------------------------------------
 * $Id$ 
 *
 * This file is part of libRadtran.
 * Copyright (c) 1997-2017 by Arve Kylling, Bernhard Mayer,
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
#include <string.h>

#include "allocnd.h"
#include "common_math.h"
#include "ancillary.h"
#include "ascii.h"
#include "errors.h"
#include "netCDF_functions.h"
#include "redistribute.h"
#include "uvspec.h"

#if HAVE_NETCDF4
#include <netcdf.h>
#define RUN_NC(command, msg)                                                                                                       \
  do {                                                                                                                             \
    int status = (command);                                                                                                        \
    if (status != NC_NOERR) {                                                                                                      \
      fprintf (stderr, "netCDF error %d: \"%s\" %s from %s\n", status, nc_strerror (status), (msg), filename);                     \
      CHKERR (status);                                                                                                             \
    }                                                                                                                              \
  } while (0)
#endif

#ifndef MAX_LENGTH_OF_LINE
#define MAX_LENGTH_OF_LINE 65536
#endif

/************************************/
/* prototypes of internal functions */
/************************************/
static int read_atmosphere_3d (char*           filename,
                               atm_out_struct* out,
                               float*          mol_mass,
                               int             interpol_method_temper,
                               int*            interpol_method_gas,
                               int             quiet,
                               int             verbose);

static int write_dummy_profile_file (char* filename, int quiet, int verbose);

static int
redistribute_3D_temp (float**** data, int nx, int ny, float* z_old, int nlyr_old, float* z_new, int nlyr_new, double* temper_1D);

int setup_molecular3d (input_struct input, output_struct* output) {
  if (!input.atmosphere3d)
    return 0;

  /* Setup 1D calculations for trace gases that are not 3D */
  /* must be set to 0, because first 1D atmospheric profiles are set up */
  output->molecular3d = 0;

  /* set dens column of included species to 0 for 1D calculations, should be specified in input file, first consider only H2O,
     use 1e-8 because 0.0 produces nan in Fu band 18.*/
  if (input.ck_scheme == CK_FU)
    input.atm.column[MOL_H2O] = 1.e-8;
  else
    input.atm.column[MOL_H2O] = 0.0;

  /* write dummy profile file */
  int status = write_dummy_profile_file (input.atmosphere3d_filename, input.quiet, input.verbose);
  CHKERR (status);

  return 0;
}

int setup_optprop_molecular3d (input_struct input, output_struct* output) {

  int       status;
  double ***tmp_molabs = NULL, ***tmp_rayleigh = NULL;
  float *   tmp_zd_common = NULL, *tmp_zd = NULL;
  int       tmp_nlev_common = 0, tmp_nlev = 0;
  double*   tmp_temper = NULL;

  if (!input.atmosphere3d)
    return 0;

  output->molecular3d = 1;

  /* We want to use the same functions for the calculation of
       molecular absorption for 1D and 3d, therefore we need to save
       the variables output->atm.optprop.tau_molabs_r and
       output->atm.optprop.tau_rayleigh_r. These are reserved for 1D
       optical properties and temporarily also used for 3D optical
       properties. */

  /* copy 1D molecular optical properties to temporary arrays */
  status = ASCII_calloc_double_3D_arylen_restricted (&tmp_molabs,
                                                     output->atm.nlev_common - 1,
                                                     output->wl.nlambda_r,
                                                     output->wl.iv_lower,
                                                     output->wl.iv_upper,
                                                     output->atm.nq_r);
  CHKERR (status);
  status = ASCII_calloc_double_3D_arylen_restricted (&tmp_rayleigh,
                                                     output->atm.nlev_common - 1,
                                                     output->wl.nlambda_r,
                                                     output->wl.iv_lower,
                                                     output->wl.iv_upper,
                                                     output->atm.nq_r);
  CHKERR (status);

  for (int lc = 0; lc < output->atm.nlev_common - 1; lc++)
    for (int iv = output->wl.iv_lower; iv <= output->wl.iv_upper; iv++)
      for (int iq = 0; iq < output->atm.nq_r[iv]; iq++) {
        tmp_molabs[lc][iv][iq]   = output->atm.optprop.tau_molabs_r[0][0][lc][iv][iq];
        tmp_rayleigh[lc][iv][iq] = output->atm.optprop.tau_rayleigh_r[0][0][lc][iv][iq];
      }

  /* We assume that z-grid of the 3d molecular atmosphere is the
       same as that of the "profile_files", therefore the common
       z-grid should not be changed after including 3D molecular
       properties. Variables are changed in functions below, so we
       need to store common z-grid temporarily. */

  // BM, 2.4.2022. Why not allow different z-grid for 3D molecular atmosphere?
  // Would need to be combined with 1D profile: 3D should ovewrite where available

  tmp_nlev = output->atm.nlev;
  tmp_zd   = calloc (output->atm.nlev, sizeof (float));
  for (int lc = 0; lc < output->atm.nlev; lc++)
    tmp_zd[lc] = output->atm.zd[lc];

  tmp_nlev_common = output->atm.nlev_common;
  tmp_zd_common   = calloc (output->atm.nlev_common, sizeof (float));
  tmp_temper      = calloc (output->atm.nlev_common, sizeof (double));
  for (int lc = 0; lc < output->atm.nlev_common; lc++)
    tmp_zd_common[lc] = output->atm.zd_common[lc];

  for (int lc = 0; lc < output->atm.nlev_common; lc++)
    tmp_temper[lc] = output->atm.microphys.temper[0][0][lc];

  /* Calculate 3D molecular optical properties */
  status = read_atmosphere_3d (input.atmosphere3d_filename,
                               &(output->atm),
                               input.atm.mol_mass,
                               input.atm.interpol_method_temper,
                               input.atm.interpol_method_gas,
                               input.quiet,
                               input.verbose);
  CHKERR (status);

  /**** Read cross section files, like O3, BrO, OclO, NO2, etc. ****/
  pmesg (" ... calling setup_crs(), reading cross section files\n", input.verbose);
  status = setup_crs (input, output); /* in molecular.c */
  CHKERR (status);

  /**** Calculate Rayleigh scattering cross section ****/
  pmesg (" ... calling setup_rayleigh(), calculating Rayleigh scattering for 3D atmosphere \n", input.verbose);
  status = setup_rayleigh (input, output); /* in molecular.c */
  CHKERR (status);

  /****  Optical properties of trace gases ****/
  pmesg (" ... calling setup_gases(), generating optical properties of trace gases for 3D atmosphere \n", input.verbose);
  status = setup_gases (input, output); /* in atmosphere.c */
  CHKERR (status);

  /**** Redistribute 3D molecular profiles to required vertical resolution ****/
  pmesg (" ... redistributing 3D molecular optical properties\n", input.verbose);

  /* common z_grid */
  free (output->atm.zd_common);
  output->atm.zd_common   = calloc (tmp_nlev_common, sizeof (float));
  output->atm.nlev_common = tmp_nlev_common;
  output->atm.nlyr_common = tmp_nlev_common - 1;

  for (int lc = 0; lc < output->atm.nlev_common; lc++)
    output->atm.zd_common[lc] = tmp_zd_common[lc];

  status = redistribute_molecular (&(output->atm.optprop),
                                   output->wl.nlambda_r,
                                   output->wl.iv_lower,
                                   output->wl.iv_upper,
                                   output->atm.nq_r,
                                   output->nipa,
                                   output->atm.zd,
                                   output->atm.nlev - 1,
                                   output->atm.Nxatm,
                                   output->atm.Nyatm,
                                   output->atm.zd_common,
                                   output->atm.nlev_common - 1);
  CHKERR (status);

  /* redistribute temperature */
  status = redistribute_3D_temp (&(output->atm.microphys.temper),
                                 output->atm.Nxatm,
                                 output->atm.Nyatm,
                                 output->atm.zd,
                                 output->atm.nlev,
                                 output->atm.zd_common,
                                 output->atm.nlev_common,
                                 tmp_temper);
  CHKERR (status);

  ASCII_free_float_3D (output->atm.microphys.temper_avg, output->atm.Nxatm, output->atm.Nyatm);
  ASCII_calloc_float_3D (&output->atm.microphys.temper_avg, output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev_common - 1);

  for (int ix = 0; ix < output->atm.Nxatm; ix++)
    for (int iy = 0; iy < output->atm.Nyatm; iy++)
      for (int lc = 0; lc < output->atm.nlev_common - 1; lc++)
        output->atm.microphys.temper_avg[ix][iy][lc] =
          0.5 * (output->atm.microphys.temper[ix][iy][lc] + output->atm.microphys.temper[ix][iy][lc + 1]);

  /* copy optical properties of 3D species to mc structure. */
  output->mc.kabs3D = output->atm.optprop.tau_molabs_r;
  output->mc.ksca3D = output->atm.optprop.tau_rayleigh_r;

  for (int ix = 0; ix < output->atm.Nxatm; ix++)
    for (int iy = 0; iy < output->atm.Nyatm; iy++)
      for (int lc = 0; lc < output->atm.nlev_common - 1; lc++)
        for (int iv = output->wl.iv_lower; iv <= output->wl.iv_upper; iv++)
          for (int iq = 0; iq < output->atm.nq_r[iv]; iq++) {
            // divide by layer thickness to convert optical thickness to extinction coefficient
            output->mc.kabs3D[ix][iy][lc][iv][iq] /= (output->atm.zd_common[lc] - output->atm.zd_common[lc + 1]);
            output->mc.ksca3D[ix][iy][lc][iv][iq] /= (output->atm.zd_common[lc] - output->atm.zd_common[lc + 1]);
          }

  // free memory
  switch (input.ck_scheme) {
  case CK_FU:
  case CK_FILE:        /* ??? need to check if that is appropriate here, probably not ??? */
  case CK_AVHRR_KRATZ: /* ??? need to check if this is appropriate here, probably not ??? */
  case CK_LOWTRAN:     /* ??? need to check if this is appropriate here, probably not ??? */
    break;

  case CK_KATO:
  case CK_KATO2:
  case CK_KATO2_96:
  case CK_KATO2ANDWANDJI:
  default:

    ASCII_free_float_5D (output->crs.crs, output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev, MOL_NN);
    break;
  }

  /* copy back 1D molecular optical properties. */
  output->atm.optprop.tau_molabs_r = calloc_float_5D_restricted (1,
                                                                 1,
                                                                 (size_t) output->atm.nlev_common - 1,
                                                                 (size_t) output->wl.iv_lower,
                                                                 (size_t) output->wl.iv_upper,
                                                                 output->atm.nq_r,
                                                                 "output->atm.optprop.tau_molabs_r");
  CHKPOINTER (output->atm.optprop.tau_molabs_r);

  output->atm.optprop.tau_rayleigh_r = calloc_float_5D_restricted (1,
                                                                   1,
                                                                   (size_t) output->atm.nlev_common - 1,
                                                                   (size_t) output->wl.iv_lower,
                                                                   (size_t) output->wl.iv_upper,
                                                                   output->atm.nq_r,
                                                                   "output->atm.optprop.tau_rayleigh_r");
  CHKPOINTER (output->atm.optprop.tau_rayleigh_r);

  free (output->atm.zd);
  output->atm.zd   = calloc (tmp_nlev, sizeof (float));
  output->atm.nlev = tmp_nlev;
  output->atm.nlyr = tmp_nlev - 1;
  for (int lc = 0; lc < output->atm.nlev; lc++)
    output->atm.zd[lc] = tmp_zd[lc];

  for (int lc = 0; lc < output->atm.nlev_common - 1; lc++)
    for (int iv = output->wl.iv_lower; iv <= output->wl.iv_upper; iv++)
      for (int iq = 0; iq < output->atm.nq_r[iv]; iq++) {
        output->atm.optprop.tau_molabs_r[0][0][lc][iv][iq]   = tmp_molabs[lc][iv][iq];
        output->atm.optprop.tau_rayleigh_r[0][0][lc][iv][iq] = tmp_rayleigh[lc][iv][iq];
      }

  /* ??? BM: quick and dirty hack: set 1D Rayleigh scattering to 0; need to change that for cases where 1D and 3D atmospheres are combined ??? */
  if (!input.quiet)
    fprintf (stderr, "... using 3D Rayleigh optical thickness\n");
  for (int lc = 0; lc < output->atm.nlev_common - 1; lc++)
    for (int iv = output->wl.iv_lower; iv <= output->wl.iv_upper; iv++)
      for (int iq = 0; iq < output->atm.nq_r[iv]; iq++)
        output->atm.optprop.tau_rayleigh_r[0][0][lc][iv][iq] = 0;

  /* free temporary variables */
  ASCII_free_double_3D (tmp_molabs, output->atm.nlev_common - 1, output->wl.nlambda_r);
  ASCII_free_double_3D (tmp_rayleigh, output->atm.nlev_common - 1, output->wl.nlambda_r);

  free (tmp_zd_common);
  free (tmp_zd);
  free (tmp_temper);

  return 0;
}

/***********************************************************************************/
/* Function: optical_properties_atmosphere3D                                       */
/* Description:                                                                    */
/*  Calculates 3D Rayleigh scattering coefficient and 3D gas absorption            */
/*  coefficient. Only molecules are included in this function, other atmospheric   */
/*  constituents should be included with options ic_file 3D, wc_file 3D,           */
/*  and profile_file 3D.                                                           */
/*                                                                                 */
/* Author: Claudia Emde                                                            */
/* Date:   May 2017                                                                */
/*                                                                                 */
/***********************************************************************************/

int optical_properties_molecular3d (input_struct input, output_struct* output, caoth3d_out_struct* caoth3d, int iv, int iq) {
  int start_lc = 0;
  int isp = 0, isp_mol = 0;
  int nlev = output->atm.nlev_common;
  int nlyr = nlev - 1;
  int Nx   = output->atm.Nxatm;
  int Ny   = output->atm.Nyatm;

  //output->mc.z       = (float *) calloc (nlev, sizeof(float));
  ASCII_calloc_float_3D (&output->mc.temper, Nx, Ny, nlyr + 1);

  /* ??? refind should also be 3D */
  //output->mc.refind  = (float *) calloc (nlev, sizeof(float));

  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      for (int lc = 0; lc < nlyr + 1; lc++) {
        if (output->atm.molabs == MOLABS_NONE)
          output->mc.kabs3D[ix][iy][lc][iv][iq] = 0.0;

        output->mc.temper[ix][iy][nlyr - lc] = output->atm.microphys.temper[ix][iy][lc];
      }
    }
  }

  // ext and ssa in caoth3d are only allocated for 3D layers, thus we need different indices here

  for (isp = 0; isp < input.n_caoth; isp++)
    if (strcmp (output->caoth3d[isp].name, "molecular_3d") == 0)
      isp_mol = isp;

  start_lc = nlyr - output->caoth3d[isp_mol].nthreed;

  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      for (int lc = start_lc; lc < nlyr; lc++) {

        switch (input.ck_scheme) {
        case CK_FU:
          /* only in the Fu and Liou case, the Rayleigh scattering */
          /* cross section depends on the subband                  */
          output->caoth3d[isp_mol].ext[nlyr - lc - 1][ix][iy] =
            0.001 * (output->mc.ksca3D[ix][iy][lc][iv][iq] + output->mc.kabs3D[ix][iy][lc][iv][iq]);
          if (output->caoth3d[isp_mol].ext[nlyr - lc - 1][ix][iy] > 0)
            output->caoth3d[isp_mol].ssa[nlyr - lc - 1][ix][iy] =
              0.001 * output->mc.ksca3D[ix][iy][lc][iv][iq] / output->caoth3d[isp_mol].ext[nlyr - lc - 1][ix][iy];
          else /* avoid dividing 0/0 in case of no extinction */
            output->caoth3d[isp_mol].ssa[nlyr - lc - 1][ix][iy] = 0;
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
          output->caoth3d[isp_mol].ext[nlyr - lc - 1][ix][iy] =
            0.001 * (output->mc.ksca3D[ix][iy][lc][iv][0] + output->mc.kabs3D[ix][iy][lc][iv][iq]);
          if (output->caoth3d[isp_mol].ext[nlyr - lc - 1][ix][iy] > 0)
            output->caoth3d[isp_mol].ssa[nlyr - lc - 1][ix][iy] =
              0.001 * output->mc.ksca3D[ix][iy][lc][iv][0] / output->caoth3d[isp_mol].ext[nlyr - lc - 1][ix][iy];
          else /* avoid dividing 0/0 in case of no extinction */
            output->caoth3d[isp_mol].ssa[nlyr - lc - 1][ix][iy] = 0;
          break;
        default:
          fprintf (stderr, "Error: unsupported correlated-k scheme %d\n", input.ck_scheme);
          return -1;

          break;
        }
      }
    }
  }

  return 0;
}

/*********************************************************/
/* Function: read_atmosphere_3d                          */
/* Description:                                          */
/* Read 3D atmosphere file and write to atm_out_struct.  */
/*                                                       */
/* Autor:                                                */
/* March 2017    C. Emde   Created                       */
/*********************************************************/
static int read_atmosphere_3d_netcdf (char*   filename,
                                      int     quiet,
                                      int     verbose,
                                      int     load_data,
                                      size_t* Nx,
                                      size_t* Ny,
                                      size_t* Nlev,
                                      double* dx,
                                      double* dy,
                                      float** zlev,
                                      float** press,
                                      float** temp,
                                      float** h2o) {
#if HAVE_NETCDF4
  int status = is_valid_netcdf (filename);
  if (status) {
    if (!quiet) {
      fprintf (stderr, " ... skipping netcdf read of 3D atmosphere data for %s ... seems it is not a valid netcdf\n", filename);
    }
    return status;
  } else {
    if (!quiet) {
      fprintf (stderr, " ... reading 3D atmosphere data from netCDF file %s\n", filename);
    }
  }

  int ncid;

  status = open_netcdf (filename, &ncid);
  CHKERROUT (status, "could not open 3D atmosphere file with netcdf reader");

  int    Ndims;
  size_t dimsizes[3];
  {
    float* x_coord = NULL;
    status         = read_netcdf_float (ncid, "x", &x_coord, &Ndims, dimsizes);
    CHKERR (status);
    *Nx = dimsizes[0];
    if (!quiet) {
      fprintf (stderr, " ... reading 3D atmosphere x coordinate N=%lu dx=%f\n", dimsizes[0], x_coord[1] - x_coord[0]);
    }
    UVSPEC_FREE (x_coord);
  }
  {
    float* y_coord = NULL;
    status         = read_netcdf_float (ncid, "y", &y_coord, &Ndims, dimsizes);
    CHKERR (status);
    *Ny = dimsizes[0];
    if (!quiet) {
      fprintf (stderr, " ... reading 3D atmosphere y coordinate N=%lu dy=%f\n", dimsizes[0], y_coord[1] - y_coord[0]);
    }
    UVSPEC_FREE (y_coord);
  }
  {
    status = read_netcdf_float (ncid, "zlev", zlev, &Ndims, dimsizes);
    CHKERR (status);
    *Nlev = dimsizes[0];
    if (!quiet) {
      fprintf (stderr, " ... reading 3D atmosphere z coordinate Nlev=%lu\n", dimsizes[0]);
    }
  }
  if (load_data) {
    status = read_netcdf_float (ncid, "press", press, &Ndims, dimsizes);
    CHKERR (status);
    if (!quiet)
      fprintf (stderr, " ... reading 3D atmosphere pressure (Nx %lu Ny %lu Nz %lu) \n", dimsizes[0], dimsizes[1], dimsizes[2]);
    status = read_netcdf_float (ncid, "temp", temp, &Ndims, dimsizes);
    CHKERR (status);
    if (!quiet)
      fprintf (stderr, " ... reading 3D atmosphere temperature (Nx %lu Ny %lu Nz %lu) \n", dimsizes[0], dimsizes[1], dimsizes[2]);
    status = read_netcdf_float (ncid, "qH2O", h2o, &Ndims, dimsizes);
    CHKERR (status);
    if (!quiet)
      fprintf (stderr, " ... reading 3D atmosphere qH2O (Nx %lu Ny %lu Nz %lu) \n", dimsizes[0], dimsizes[1], dimsizes[2]);
  }

  {
    double* attr_val = NULL;
    status           = read_netcdf_attribute_double (ncid, NULL, "dx", &attr_val);
    CHKERR (status);
    *dx = attr_val[0];

    status = read_netcdf_attribute_double (ncid, NULL, "dy", &attr_val);
    CHKERR (status);
    *dy = attr_val[0];
    UVSPEC_FREE (attr_val);

    if (!quiet) {
      fprintf (stderr, " ... reading 3D atmosphere attributes dx=%g dy=%g\n", *dx, *dy);
    }
  }
#else
  if (!quiet)
    fprintf (stderr, "... not trying to read <%s> with netcdf because not compiled with netcdf support\n", filename);
  return -1;
#endif
  return 0;
}

static int read_atmosphere_3d_ascii (char*   filename,
                                     int     quiet,
                                     int     verbose,
                                     int     load_data,
                                     size_t* Nx,
                                     size_t* Ny,
                                     size_t* Nlev,
                                     double* dx,
                                     double* dy,
                                     float** zlev,
                                     float** pres,
                                     float** temp,
                                     float** h2o) {

  int rows;
  int min_columns;
  int max_columns;
  int max_length;

  /* check input file */
  int status = ASCII_checkfile (filename, &rows, &min_columns, &max_columns, &max_length);
  CHKERR (status);

  if (!quiet)
    fprintf (stderr, " ... reading 3D atmosphere data from ascii file %s\n", filename);

  if (rows < 2) {
    fprintf (stderr, "Error: found less than two rows in %s\n", filename);
    CHKERR (LESS_THAN_TWO_COLUMNS);
  }

  if (min_columns < 3) {
    fprintf (stderr, "Error: found less than three columns in %s\n", filename);
    CHKERR (LESS_THAN_THREE_COLUMNS);
  }

  /* reset string to the beginning of the memory block */
  char  line[MAX_LENGTH_OF_LINE + 1] = "";
  char* string                       = line;

  /* open file */
  FILE* file;
  if ((file = fopen (filename, "r")) == NULL)
    return ASCIIFILE_NOT_FOUND;

  char** array = NULL;
  char*  dummy;
  /* first line */
  int number = 0;
  while (number == 0) {
    dummy  = fgets (string, MAX_LENGTH_OF_LINE, file);
    status = ASCII_parsestring (string, &array, &number);
    CHKERR (status);
  }

  if (number < 3) {
    fprintf (stderr, "Error: expected at least 3 values in 1st line, found only %d\n", number);
    CHKERR (-1);
  }

  /* first 3 values are number of cells */
  *Nx   = strtol (array[0], &dummy, 0);
  *Ny   = strtol (array[1], NULL, 0);
  *Nlev = strtol (array[2], NULL, 0);

  free (array);

  /* 2nd line */
  number = 0;
  while (number == 0) {
    dummy  = fgets (string, MAX_LENGTH_OF_LINE, file);
    status = ASCII_parsestring (string, &array, &number);
    CHKERROUT (status, "failed reading 1st line");
  }

  if (number != (*Nlev) + 2) {
    fprintf (stderr, "Error: found %d z levels, expected %lu\n", number - 2, *Nlev);
    CHKERR (-1);
  }

  /* horizontal grid sizes, in km */
  *dx = strtod (array[0], NULL);
  *dy = strtod (array[1], NULL);

  /* rest of line 2 */
  /* set altitude levels */
  *zlev = malloc ((*Nlev) * sizeof (float));
  for (int i = 2; i < (*Nlev) + 2; i++) {
    /*altitudes in out->zd are sorted in descending order,*/
    /* 3D input files in ascending order */
    (*zlev)[i - 2] = (float) strtod (array[(*Nlev) + 3 - i], NULL);
  }
  free (array);

  if (load_data) {

    const size_t Ndim        = 3;
    const size_t dimsizes[3] = {(*Nx), (*Ny), (*Nlev)};
    size_t       offsets[3];
    status = nd_offsets (3, dimsizes, offsets);
    CHKERR (status);

    *pres            = calloc_float_1D ((*Nx) * (*Ny) * (*Nlev), "pressure from 3D atmosphere ascii file");
    *temp            = calloc_float_1D ((*Nx) * (*Ny) * (*Nlev), "temperature from 3D atmosphere ascii file");
    *h2o             = calloc_float_1D ((*Nx) * (*Ny) * (*Nlev), "density from 3D atmosphere ascii file");
    char* have_entry = calloc_char_1D ((*Nx) * (*Ny) * (*Nlev), "3D atmosphere file has to have entries");

    for (int i = 2; i < rows; i++) {
      number = 0;
      while (number == 0) {
        dummy  = fgets (string, MAX_LENGTH_OF_LINE, file);
        status = ASCII_parsestring (string, &array, &number);

        if (status != 0) {
          fprintf (stderr, "Error %d reading line nr %i of %s\n", status, i, filename);
          CHKERR (status);
        }
      }

      /* cell indices */
      const int ix = strtol (array[0], NULL, 0) - 1;
      const int iy = strtol (array[1], NULL, 0) - 1;
      const int iz = strtol (array[2], NULL, 0) - 1;

      if (ix < 0 || ix >= *Nx) {
        fprintf (stderr, "Error, ix = %d out of bounds in line %d\n", ix + 1, i + 1);
        CHKERR (-1);
      }

      if (iy < 0 || iy >= *Ny) {
        fprintf (stderr, "Error, iy = %d out of bounds in line %d\n", iy + 1, i + 1);
        CHKERR (-1);
      }

      if (iz < 0 || iz >= *Nlev) {
        fprintf (stderr, "Error, iz = %d out of bounds in line %d\n", iz + 1, i + 1);
        CHKERR (-1);
      }

      const size_t idx3d[3] = {ix, iy, (*Nlev) - 1 - iz};  // note: vert coordinate reversed
      const size_t idx      = idx1d (Ndim, idx3d, offsets);
      (*pres)[idx]          = strtod (array[3], NULL);
      (*temp)[idx]          = strtod (array[4], NULL);
      (*h2o)[idx]           = strtod (array[5], NULL);
      have_entry[idx]       = 1;

      /*other trace gases can follow in further columns, to be implemented. */

      free (array);
    }  // for row

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", rows - 2, filename);

    status     = 0;
    size_t idx = 0;
    for (size_t i = 0; i < (*Nx); ++i) {
      for (size_t j = 0; j < (*Ny); ++j) {
        for (size_t k = 0; k < (*Nlev); ++k, ++idx) {
          if (!have_entry[idx]) {
            fprintf (stderr, "3D atmosphere has not set values for level: i=%lu j=%lu k=%lu\n", i + 1, j + 1, (*Nlev) - k);
            status++;
          }
        }
      }
    }
    UVSPEC_FREE (have_entry);
    CHKERROUT (status, "3D atmosphere has to have values for all levels");
  }

  fclose (file);

  return 0;
}

static int read_atmosphere_3d (char*           filename,
                               atm_out_struct* out,
                               float*          mol_mass,
                               int             interpol_method_temper,
                               int*            interpol_method_gas,
                               int             quiet,
                               int             verbose) {
  int status;

  const int load_data = 1;
  size_t    Nx, Ny, Nlev;
  double    dx, dy;
  float*    zlev  = NULL;
  float*    press = NULL;
  float*    temp  = NULL;
  float*    h2o   = NULL;

  status = read_atmosphere_3d_netcdf (filename, quiet, verbose, load_data, &Nx, &Ny, &Nlev, &dx, &dy, &zlev, &press, &temp, &h2o);

  if (status) {
    if (!quiet) {
      fprintf (stderr, " ... reading 3D atmosphere data %s with netCDF reader failed... trying with ASCII reader\n", filename);
    }

    status = read_atmosphere_3d_ascii (filename, quiet, verbose, load_data, &Nx, &Ny, &Nlev, &dx, &dy, &zlev, &press, &temp, &h2o);
    CHKERROUT (status, "3D atmosphere file could not be read with netcdf nor with ascii routine");
  }

  if (zlev[0] <= zlev[Nlev - 1])
    CHKERROUT (1, "Wrong ordering of height levels, needs to be reversed");
  if (out->nlev != Nlev) {
    for (size_t k = 0; k < MIN (out->nlev, Nlev); ++k) {
      if (!approx_abs (out->zd[k], zlev[k], 1e-7))
        fprintf (stderr, "out->zd and zlev from 3D atmosphere file differ at k=%lu : %g vs %g\n", k, out->zd[k], zlev[k]);
    }
    CHKERROUT (1, "Currently, a 3D atmosphere file has to have the same number of levels");
  }

  out->Nxatm       = Nx;
  out->Nyatm       = Ny;
  out->nlev        = Nlev;
  out->nlyr        = out->nlev - 1;
  out->nlev_common = Nlev;
  out->nlyr_common = out->nlev_common - 1;
  out->dxatm       = dx;
  out->dyatm       = dy;

  /* free and reallocate memory for altitude levels */
  free (out->zd);
  free (out->zd_common);
  out->zd        = calloc (out->nlev, sizeof (float));
  out->zd_common = calloc (out->nlev, sizeof (float));

  /* 3DAbs for now we assume the same 3D grid for the full atmosphere, so we do not need */
  /* the variables z_3d_layer and threed */
  /* out->z_3d_layer = calloc ((size_t) out->Nlyr_3d, sizeof(float)); */
  /* out->threed = calloc ((size_t) out->Nlyr_3d, sizeof(int)); */

  /* free and reallocate memory for pressure, temperature, water vapor ... */
  ASCII_free_float_3D (out->microphys.press, 1, 1);
  ASCII_calloc_float_3D (&out->microphys.press, out->Nxatm, out->Nyatm, out->nlev);
  ASCII_free_float_3D (out->microphys.temper, 1, 1);
  ASCII_calloc_float_3D (&out->microphys.temper, out->Nxatm, out->Nyatm, out->nlev);
  ASCII_free_float_4D (out->microphys.dens, MOL_NN, 1, 1);
  ASCII_calloc_float_4D (&out->microphys.dens, MOL_NN, out->Nxatm, out->Nyatm, out->nlev);

  //3DAbs avg may be allocated later
  ASCII_calloc_float_3D (&out->microphys.temper_avg, out->Nxatm, out->Nyatm, out->nlev - 1);
  ASCII_calloc_float_4D (&out->microphys.dens_avg, MOL_NN, out->Nxatm, out->Nyatm, out->nlev - 1);

  /* copy to out struct */
  for (size_t k = 0; k < out->nlev; ++k) {
    out->zd[k]        = zlev[k];
    out->zd_common[k] = zlev[k];
  }
  UVSPEC_FREE (zlev);

  {
    const size_t Ndim        = 3;
    const size_t dimsizes[3] = {Nx, Ny, Nlev};
    size_t       offsets[3];
    status = nd_offsets (3, dimsizes, offsets);
    CHKERR (status);
    for (int ix = 0; ix < out->Nxatm; ix++) {
      for (int iy = 0; iy < out->Nyatm; iy++) {
        for (int iz = 0; iz < out->nlev; iz++) {
          const size_t idx3d[3]                    = {ix, iy, iz};
          const size_t idx                         = idx1d (Ndim, idx3d, offsets);
          out->microphys.press[ix][iy][iz]         = press[idx];
          out->microphys.temper[ix][iy][iz]        = temp[idx];
          out->microphys.dens[MOL_H2O][ix][iy][iz] = h2o[idx];
        }
      }
    }
    UVSPEC_FREE (press);
    UVSPEC_FREE (temp);
    UVSPEC_FREE (h2o);
  }

  /* calculate air number density from pressure and temperature, as well as layer average quantities */
  for (int ix = 0; ix < out->Nxatm; ix++) {
    for (int iy = 0; iy < out->Nyatm; iy++) {
      for (int iz = 0; iz < out->nlev; iz++) {
        out->microphys.dens[MOL_AIR][ix][iy][iz] =
          out->microphys.press[ix][iy][iz] * 1E-04 / (BOLTZMANN * out->microphys.temper[ix][iy][iz]);

        /* convert from kg/kg (specific humidity) to concentration 1/cm^3  */
        out->microphys.dens[MOL_H2O][ix][iy][iz] = mol_mass[MOL_AIR] / mol_mass[MOL_H2O] *
                                                   out->microphys.dens[MOL_H2O][ix][iy][iz] *
                                                   out->microphys.dens[MOL_AIR][ix][iy][iz];
      }

      /* layer average temperature */
      status = average_dens (out->microphys.temper[ix][iy],
                             out->microphys.dens[MOL_AIR][ix][iy],
                             out->zd_common,
                             out->nlyr + 1,
                             interpol_method_temper,
                             &(out->microphys.temper_avg[ix][iy]),
                             NO);
      CHKERR (status);

      /* layer average air density  */
      status = average_dens (out->microphys.dens[MOL_AIR][ix][iy],
                             out->microphys.dens[MOL_AIR][ix][iy],
                             out->zd_common,
                             out->nlyr + 1,
                             interpol_method_gas[MOL_H2O],
                             &(out->microphys.dens_avg[MOL_AIR][ix][iy]),
                             NO);
      CHKERR (status);

      /* layer average water vapor */
      status = average_dens (out->microphys.dens[MOL_H2O][ix][iy],
                             out->microphys.dens[MOL_AIR][ix][iy],
                             out->zd_common,
                             out->nlyr + 1,
                             interpol_method_gas[MOL_H2O],
                             &(out->microphys.dens_avg[MOL_H2O][ix][iy]),
                             NO);
      CHKERR (status);
    }
  }

  return 0;
}

/*********************************************************/
/* Function: write_dummy_profile_file                    */
/* Description:                                          */
/* Generate dummy profile_file from 3d molecular file    */
/*                                                       */
/* Author:                                                */
/* August 2017    C. Emde   Created                      */
/*********************************************************/
static int write_dummy_profile_file (char* filename, int quiet, int verbose) {

  int status;

  if (!quiet)
    fprintf (stderr, " ... generate dummy profile file from 3d atmosphere file %s\n", filename);

  size_t Nx, Ny, Nlev;
  double dx, dy;
  float* zlev = NULL;

  {
    const int load_data = 0;
    status = read_atmosphere_3d_netcdf (filename, quiet, verbose, load_data, &Nx, &Ny, &Nlev, &dx, &dy, &zlev, NULL, NULL, NULL);

    if (status) {
      if (!quiet) {
        fprintf (
          stderr,
          " ... generate dummy profile file from 3d atmosphere file %s with netCDF reader failed... trying with ASCII reader\n",
          filename);
      }

      status = read_atmosphere_3d_ascii (filename, quiet, verbose, load_data, &Nx, &Ny, &Nlev, &dx, &dy, &zlev, NULL, NULL, NULL);
      CHKERROUT (status, "3D atmosphere file could not be read with netcdf nor with ascii routine");
    }
  }

  // ??? check disabled 12.9.2018, BM. At present, all layers are assumed 3D anyway
  /* BM, would need to read rest of atm file - we need to know which layers have 3D absorption */
  // Need to reactivate later!
  int threed[Nlev - 1];
  for (int k = 0; k < Nlev - 1; k++)
    threed[k] = 1;

  if (!quiet) {
    fprintf (stderr, "Attention! It is not checked if there are any non-3D layers in the 3D atmosphere file\n");
    fprintf (stderr, "since this is not possible at the moment anyway; however, once this is possible,\n");
    fprintf (stderr, "the check should be activated again in molecular3d.c\n");
  }

  for (int k = 0; k < Nlev - 1; k++)
    if (threed[k] && !quiet)
      fprintf (stderr, "FOUND absorber in layer %d\n", k);

  FILE* file_dummy = NULL;
  /* write dummy profile file */
  if ((file_dummy = fopen ("profile_mol3d_dummy.dat", "w")) == NULL)
    return ASCIIFILE_NOT_FOUND;

  fprintf (file_dummy, "%lu %lu %lu 3\n", Nx, Ny, Nlev - 1);
  fprintf (file_dummy, "%f %f ", dx, dy);
  for (int i = 0; i < Nlev; i++)
    fprintf (file_dummy, "%f ", zlev[Nlev - i - 1]);
  fprintf (file_dummy, "\n");

  /* BM, write all layers with 3D absorption */
  for (int k = 0; k < Nlev - 1; k++)
    if (threed[k])
      fprintf (file_dummy, "1 1 %d 0.0 10 \n", k + 1);

  fclose (file_dummy);

  UVSPEC_FREE (zlev);

  return 0;
}

/****************************************************************/
/* Regrid a 3D float array [ix, iy, iz] from a given input grid */
/* (z_old) to a given output grid.                              */
/****************************************************************/

static int
redistribute_3D_temp (float**** data, int nx, int ny, float* z_old, int nlev_old, float* z_new, int nlev_new, double* temper_1D) {
  float ***tmp = NULL, *z_tmp = NULL, *data_tmp = NULL;

  ASCII_calloc_float_3D (&tmp, nx, ny, nlev_old);
  /* copy original stuff to temporary array */

  /* common grid range, assume that 3D grid is defined inside lower part of 1D standard atmosphere */
  int start_lc = 0;
  while (z_new[start_lc] > z_old[0])
    start_lc++;

  const int nlev_tmp = nlev_new - start_lc;
  z_tmp              = calloc (nlev_tmp, sizeof (float));

  for (int lc = 0; lc < nlev_new - start_lc; lc++)
    z_tmp[lc] = z_new[lc + start_lc];

  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      for (int lc = 0; lc < nlev_old; lc++)
        tmp[ix][iy][lc] = (*data)[ix][iy][lc];

  data_tmp = calloc (nlev_tmp, sizeof (float));
  ASCII_free_float_3D (*data, nx, ny);
  ASCII_calloc_float_3D (data, nx, ny, nlev_new);

  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {

      int status = arb_wvn (nlev_old, z_old, tmp[ix][iy], nlev_tmp, z_tmp, data_tmp, 1, 1);

      if (status != 0) {
        fprintf (stderr, "Error %d returned by arb_wvn() in redistribute_temp (molecular_3d.c)\n", status);
        CHKERR (status);
      }

      for (int lc = 0; lc < start_lc; lc++)
        (*data)[ix][iy][lc] = temper_1D[lc];

      for (int lc = start_lc; lc < nlev_new; lc++)
        (*data)[ix][iy][lc] = data_tmp[lc - start_lc];
    }
  }

  ASCII_free_float_3D (tmp, nx, ny);
  free (data_tmp);
  free (z_tmp);

  return 0;
}
