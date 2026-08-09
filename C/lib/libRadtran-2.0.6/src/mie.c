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

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "ascii.h"
#include "errors.h"
#include "mie.h"
#include "numeric.h"

#if HAVE_NETCDF4
#include <netcdf.h>
#endif

#define MAXLINE 255

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

int aerosol_input_data (/* Input */
                        char* data_files_path,
                        int   aerosol_type,
                        int   verbose,
                        /* Output */
                        aerdist_struct* aerdistr);

int aerosol_refractive_index (/* Input */
                              double* lambda,
                              int     nlambda,
                              char*   data_files_path,
                              int     aerosol_type,
                              int     n_hum,
                              int     verbose,
                              /* Output */
                              double*** ref_real,
                              double*** ref_imag);

int aerosol_size_distribution (/* Input */
                               double rmin,
                               double rmax,
                               double rmod,
                               double sigma,
                               double lambda,
                               int    verbose,
                               /* Output */
                               int*     ndens,
                               double** radius,
                               double** number_dens);

int calc_distribution_radius_grid (/* Input */
                                   int     distribution,
                                   char*   sd_filename,
                                   double  r_eff_min,
                                   double  r_eff_max,
                                   double* reff,
                                   int     nreff,
                                   double  alpha,
                                   int     nalpha,
                                   double  step,
                                   double  lambda, /*wavelength in micrometer*/
                                   double  dx_max,
                                   double  n_r_max,
                                   int     verbose,
                                   /* Output */
                                   dist_struct** dist);

int calc_distribution_number_dens (/* Input */
                                   int     distribution,
                                   double  alpha,
                                   double  r_eff,
                                   double  rho_medium,
                                   int     n_dens,
                                   double* radius,
                                   int     verbose,
                                   /* Output */
                                   double* number_dens);

int calc_phase_function (/* Input */
                         int      nmom,
                         double   lambda,
                         double   reff,
                         int      nstokes,
                         int      nthetamax,
                         int      phase_from_amplitudes,
                         double** phase_in,
                         double*  xmu_in,
                         int      ntheta_in,
                         int      ntheta_user,
                         double   accuracy,
                         /* Output */
                         float*** pmom,
                         float*** phase,
                         float*** theta,
                         int**    ntheta,
                         int*     warnings_ntheta);

double double_min (double x, double y);

double double_max (double x, double y);

double gamma_ln (double xx);

int integrate_sizedist (/* Input */
                        float*   qext_array,
                        float*   qsca_array,
                        float*   gsca_array,
                        float*   qback_array,
                        float*** pmom_array,
                        int      nstokes,
                        int      nmom,
                        int      numang,
                        double*  x_size,
                        double*  y_size,
                        int      n_size,
                        int      verbose,
                        /* Output */
                        double*  beta,
                        double*  omega,
                        double*  g,
                        double*  back,
                        float*** pmom);

int mie_calc_all_radii (/* Input */
                        int            program,
                        int            medium,
                        int            water_properties,
                        mie_complex    crefin,
                        double         wavelength,
                        float          temperature,
                        int            nstokes,
                        int            numang,
                        dist_struct**  dist,
                        int            ia,
                        int            ir,
                        int            nalpha,
                        int            n_r_eff,
                        mie_inp_struct input,
                        int            ncid_phase,
                        int            varid,
                        int            dimids[3],
                        int            write_all_phasefunctions,
                        int            verbose,
                        /* Output */
                        mie_out_struct* output,
                        mie_complex*    ref,
                        float**         qext_array,
                        float**         qsca_array,
                        float**         gsca_array,
                        float**         qback_array,
                        float****       pmom_array,
                        double****      phase,
                        double**        phase_normalization);

int read_refractive_index (/* Input */
                           char* filename,
                           int   verbose,
                           /* Output */
                           int*      nlambda,
                           double**  lambda,
                           double*** ri_real,
                           double*** ri_imag);

int setup_alpha (/* Output */
                 double** alpha,
                 int*     nalpha,
                 /* Input */
                 char*  alpha_filename,
                 double input_alpha,
                 int    verbose);

int setup_reff_grid (/* Input */
                     int     distribution,
                     double  r_eff_step,
                     double  r_eff_min,
                     double  r_eff_max,
                     int     r_eff_log,
                     int     n_hum, /* nhum is used for aerosol*/
                     double* hum,
                     int     verbose,
                     /* Output */
                     int*     n_r_eff,
                     double** r_eff);

int setup_wavelength_grid (/* Input */
                           double wl_start,
                           double wl_end,
                           double wl_step,
                           char*  data_files_path,
                           char*  filename_ck_generic,
                           int    ck_scheme,
                           int    lambda_i_start,
                           int    lambda_i_end,
                           int    verbose,
                           /* Output */
                           int*     nlambda,
                           double** lambda);

static int write_output (input_struct input, output_struct* output, int ir, int iv);
static int sortasc (const void* x1, const void* x2);

#if HAVE_NETCDF4

int create_netcdf_file (/* Input */
                        char*   basename,
                        int     medium,
                        int     aerosol_type,
                        int     distribution,
                        double  alpha,
                        double  step,
                        int     nlam,
                        double* lambda,
                        int     nreff,
                        double* reff,
                        int     nmommax,
                        int     nstokes,
                        int     nthetamax,
                        int     nrho,
                        int     complete_moments,
                        double  n_r_max,
                        double  dx_max,
                        char*   sd_filename,
                        /* Output */
                        netcdf_ids* ncid,
                        char**      outputfile);

int create_phase_netcdf_file (/* Input */
                              char*   phasefilenamecdf,
                              int     nreff,
                              int     ntheta,
                              int     nstokes,
                              double* reff,
                              double* theta,
                              /* Output */
                              netcdf_ids* ncid);

static int write_output_netcdf (netcdf_ids ncid, input_struct input, output_struct* output, int iv, int ir, int ia);

#endif

/**
 * mie
 * 
 * Main function of mie tool.
 *
 *
 * @author Bernhard Mayer, Ulrich Hamann, Claudia Emde
 * @date  ????? 
 *        2010-02-11 Restructured and cleaned (CE)        
 */
int mie (input_struct input, output_struct* output) {

  int ip = 0, im = 0, ir = 0, iv = 0, ia = 0, status = 0;
  int firstwarn = 1;

  double wavelength = 0;

  float*   qext_array;
  float*   qsca_array;
  float*   gsca_array;
  float*   qback_array;
  float*** pmom_array; /* Phase function moments */
  double   norm = 0.0;

  double** ref_real;
  double** ref_imag;

  double ext, vol;
  int    lambda_i_start = 0;
  int    lambda_i_end   = 0;

  mie_complex* crefin = NULL;

  output->ref.re = 0.0;
  output->ref.im = 0.0;

  int varid_r, varid_theta, varid_phamat;
  int dimids[3];
  int ncid_phase;

#if HAVE_NETCDF4
  netcdf_ids ncid;
  int        nmommax;
  int        complete_moments = 0;
#endif

  /***********************************************************/
  /* Setup                                                   */
  /***********************************************************/

  /* Read input data and specify wavelengths grid*/
  if (strlen (input.data_files_path) == 0)
    input.data_files_path = "../data/";

  /* Check whether refractive index file is defined and if yes get
     wavelength grid from refractive index file */
  if (strlen (input.ri_filename) > 0) {
    status =
      read_refractive_index (input.ri_filename, input.verbose, &(output->wl.nlambda), &(output->wl.lambda), &ref_real, &ref_imag);
    if (status != 0) {
      fprintf (stderr, "Error %d returned by read_refractive_index()\n", status);
      return status;
    }
  }
  /* Else set up the wavelenth grid from input file specification */
  else {
    status = setup_wavelength_grid (input.wl.wl_start,
                                    input.wl.wl_end,
                                    input.wl.wl_step,
                                    input.data_files_path,
                                    input.filename_ck_generic,
                                    input.ck_scheme,
                                    input.lambda_i_start,
                                    input.lambda_i_end,
                                    input.verbose,
                                    &(output->wl.nlambda),
                                    &(output->wl.lambda));
    if (status != 0) {
      fprintf (stderr, "Error %d returned by setup_wavelengths_grid()\n", status);
      return status;
    }
  }
  if (input.lambda_i_start != -999) {
    lambda_i_start = input.lambda_i_start;
    lambda_i_end   = input.lambda_i_end;
  } else {
    lambda_i_start = 0;
    lambda_i_end   = output->wl.nlambda - 1;
  }

  /* Read aerosol distribution parameters from OPAC*/
  if (input.distribution == DIST_AER) {
    status = aerosol_input_data (input.data_files_path, input.aerosol_type, input.verbose, &(output->aerdistr));
    if (status != 0) {
      fprintf (stderr, "Error %d returned by aerosol_input_data()\n", status);
      return status;
    }
    status = aerosol_refractive_index (output->wl.lambda,
                                       output->wl.nlambda,
                                       input.data_files_path,
                                       input.aerosol_type,
                                       output->aerdistr.n_hum,
                                       input.verbose,
                                       &ref_real,
                                       &ref_imag);
    if (status != 0) {
      fprintf (stderr, "Error %d returned by aerosol_refractive_index()\n", status);
      return status;
    }
  }

  /* Width of size distribution, parameter alpha */
  status = setup_alpha (&(output->alpha), &(output->nalpha), input.alpha_filename, input.alpha, input.verbose);

  if (output->nalpha > 1 && output->wl.nlambda > 1) {
    fprintf (stderr, "Error, option alpha_file only implemented for calculations at one wavelength. \n");
    return -1;
  }

  /* Set up the effective radius grid */
  status = setup_reff_grid (input.distribution,
                            input.r_eff_step,
                            input.r_eff_min,
                            input.r_eff_max,
                            input.r_eff_log,
                            output->aerdistr.n_hum,
                            output->aerdistr.rel_hum,
                            input.verbose,
                            &(output->n_r_eff),
                            &(output->r_eff));
  if (status != 0) {
    fprintf (stderr, "Error %d returned by setup_reff_grid()\n", status);
    return status;
  }

  /**********************************************************/
  /* Memory allocation                                      */
  /**********************************************************/

  /* refractive index may be size dependent (aerosol) */
  crefin = calloc (output->n_r_eff, sizeof (mie_complex));

  /* Legendre polynomials of all phase matrix elements */
  if (input.mie.nmom > 0) {
    /* Output of Mie program for one size */
    output->mie.pmom = calloc (input.nstokes, sizeof (float*));
    for (ip = 0; ip < input.nstokes; ip++)
      output->mie.pmom[ip] = calloc (input.mie.nmom + 1, sizeof (float));

    /* Integrated, normalized and converted moments */
    output->pmom = calloc (input.nstokes, sizeof (float*));
    for (ip = 0; ip < input.nstokes; ip++)
      output->pmom[ip] = calloc (input.mie.nmom + 1, sizeof (float));
  }

  /* Memory for scattering amplitudes s1 and s2 */
  if (input.ntheta_user > 0)
    output->numang = input.ntheta_user;
  else {
    if (input.output_format == OUTPUT_FORMAT_NETCDF && input.phase_from_amplitudes)
      output->numang = 180001;
    else
      output->numang = input.mie.numang;
  }

  int dim = 0;
  if (input.program == BHMIE)
    dim = 2 * MXNANG - 1;  // do we really need that or would input.mie.numang ne sufficient?
                           // dimension in bhmie.f is 2*MXNANG-1
  else
    dim = output->numang;

  /* Output of Mie program for one size; directly passed to MIEV0/BHMIE */
  output->mie.xmu = calloc (dim, sizeof (double));
  output->mie.s1  = calloc (dim, sizeof (mie_complex));
  output->mie.s2  = calloc (dim, sizeof (mie_complex));

  if (output->numang == 1)
    output->mie.xmu[0] = 0;
  else if (output->numang > 1) {  // MIEV0 expects angles in xmu[]
    fprintf (stderr, "allocating %d angles for xmu\n", output->numang);

    double dang      = M_PI / (double) (output->numang - 1);
    double theta_off = 0;
    if (input.ntheta_user > 0) {
      dang      = (input.thetamax_user - input.thetamin_user) / 180.0 * M_PI / (double) (output->numang - 1);
      theta_off = input.thetamin_user / 180.0 * M_PI;
    }

    for (int iang = 0; iang < output->numang; iang++) {
      double theta          = theta_off + iang * dang;  // angles MUST be increasing in MIEV0
      output->mie.xmu[iang] = cos (theta);
    }
  }

  /* Warnings */
  output->warnings_nmom   = calloc (output->wl.nlambda, sizeof (int*));
  output->warnings_ntheta = calloc (output->wl.nlambda, sizeof (int*));
  for (iv = 0; iv < output->wl.nlambda; iv++) {
    output->warnings_nmom[iv]   = calloc (output->n_r_eff, sizeof (int));
    output->warnings_ntheta[iv] = calloc (output->n_r_eff, sizeof (int));
  }

  if (input.output_format == OUTPUT_FORMAT_NETCDF) {

    output->outputfile = calloc (FILENAME_MAX, sizeof (char));

    /* allocate memory for phase function*/
    if (input.nthetamax > 0) {
      output->ntheta = calloc (input.nstokes, sizeof (int));
      output->theta  = calloc (input.nstokes, sizeof (float*));
      output->phase  = calloc (input.nstokes, sizeof (float*));
      for (ip = 0; ip < input.nstokes; ip++) {
        output->theta[ip] = calloc (input.nthetamax, sizeof (float));
        output->phase[ip] = calloc (input.nthetamax, sizeof (float));
      }
    }
  }

  /**************************************************************/
  /* Start calculation                                          */
  /**************************************************************/
  /* loop over wavelengths */
  for (iv = lambda_i_start; iv <= lambda_i_end; iv++) {

    if (input.verbose)
      fprintf (stderr, " ... start Mie-calculation for lambda = %10.3f nm  \n", output->wl.lambda[iv]);

    /* Convert from nm to micro meter */
    wavelength = (output->wl.lambda[iv]) * 1.e-03;

    /* Aerosol treatment different from clouds for consistency reasons
       with OPAC, size distribution calculated below in loop over
       relative humudities. For each relative humidity a new Mie
       calculation is required, because the refractive index changes. */

    // allocate memory for nalpha x n_r_eff size distributions */
    output->dist = calloc (output->nalpha, sizeof (dist_struct*));
    for (ia = 0; ia < output->nalpha; ia++)
      output->dist[ia] = calloc (output->n_r_eff, sizeof (dist_struct));

    if (input.distribution != DIST_AER) {

      /* generate the radius grid of the size distribution */
      /* memory for radius grid and number density is      */
      /* allocated in this function                        */
      status = calc_distribution_radius_grid (input.distribution,
                                              input.sd_filename,
                                              input.r_eff_min,
                                              input.r_eff_max,
                                              output->r_eff,
                                              output->n_r_eff,
                                              input.alpha,
                                              output->nalpha,
                                              input.step,
                                              wavelength,
                                              input.dx_max,
                                              input.n_r_max,
                                              input.verbose,
                                              output->dist);
      if (status != 0) {
        fprintf (stderr, "Error %d returned by calc_distribution_radius_grid()\n", status);
        return status;
      }
    }

    /* NEXT STEP: move generation of size distribution out of loop */
    /* CAREFUL: for aerosol the number of points in the size distribution is different for each relative humidity (ir) */

    for (ir = 0; ir < output->n_r_eff; ir++) { /* Loop over effective radii (cloud) / relative humidities (aerosol) */

      /* refractive index */
      if (strlen (input.ri_filename) > 0) {
        crefin[ir].re = ref_real[0][iv];
        crefin[ir].im = ref_imag[0][iv];
      } else if (input.distribution == DIST_AER) {
        crefin[ir].re = ref_real[ir][iv];
        crefin[ir].im = ref_imag[ir][iv];
        fprintf (stderr, " ... aerosol refractive index %g %g \n", crefin[ir].re, crefin[ir].im);
      } else {
        crefin[ir].re = input.crefin.re;
        crefin[ir].im = input.crefin.im;
      }
    }

    /* Allocate memory for scattering matrix elements */
    output->phase_amp = calloc (input.nstokes, sizeof (double***));

    for (ip = 0; ip < input.nstokes; ip++) {
      output->phase_amp[ip] = calloc (output->nalpha, sizeof (double**));

      for (ia = 0; ia < output->nalpha; ia++) {
        output->phase_amp[ip][ia] = calloc (output->n_r_eff, sizeof (double*));

        for (ir = 0; ir < output->n_r_eff; ir++)
          output->phase_amp[ip][ia][ir] = calloc (output->numang, sizeof (double));
      }
    }

    /* Allocate memory for normalization factor for phase matrix */
    output->phase_normalization = calloc (output->nalpha, sizeof (double*));
    for (ia = 0; ia < output->nalpha; ia++)
      output->phase_normalization[ia] = calloc (output->n_r_eff, sizeof (double));

    // pre-calculate all size distributions because we need all of them in the first call to mie_calc_all_radii()
    for (ia = 0; ia < output->nalpha; ia++) {    /* loop over width of size distribution if specified in alpha_filename */
      for (ir = 0; ir < output->n_r_eff; ir++) { /* Loop over effective radii (cloud) / relative humidities (aerosol) */

        if (input.distribution != DIST_AER) {
          /* number density distribution per radius intervall*/
          status = calc_distribution_number_dens (input.distribution,
                                                  output->alpha[ia],
                                                  output->r_eff[ir],
                                                  output->rho_medium,
                                                  output->dist[ia][ir].n_dens,
                                                  output->dist[ia][ir].radius,
                                                  input.verbose,
                                                  output->dist[ia][ir].number_dens);
          if (status != 0) {
            fprintf (stderr, "Error %d returned by calc_distribution_number_dens()\n", status);
            return status;
          }
        } else {
          /* Set up aerosol size distribution radius grid and calculate
	     number densities.*/
          status = aerosol_size_distribution (output->aerdistr.r_min[ir],
                                              output->aerdistr.r_max[ir],
                                              output->aerdistr.rmod[ir],
                                              output->aerdistr.sigma,
                                              output->wl.lambda[iv],
                                              input.verbose,
                                              &output->dist[ia][ir].n_dens,
                                              &(output->dist[ia][ir].radius),
                                              &(output->dist[ia][ir].number_dens));
          if (status != 0) {
            fprintf (stderr, "Error %d returned by aerosol_size_distribution()\n", status);
            return status;
          }
        }
      }
    }

    float* r     = calloc (output->dist[0][0].n_dens, sizeof (float));
    float* theta = calloc (output->numang, sizeof (float));

    for (int iang = 0; iang < output->numang; iang++)
      theta[iang] = acos (output->mie.xmu[iang]) * 180 / M_PI;

    for (int is = 0; is < output->dist[0][0].n_dens; is++)
      r[is] = output->dist[0][0].radius[is];

    // Create NetCDF file to store all phase functions
    if (input.write_all_phasefunctions && input.phase_from_amplitudes) {

      if ((status = nc_create ("phase.cdf", NC_CLOBBER, &ncid_phase)) != NC_NOERR) {
        fprintf (stderr, "Error creating NetCDF file: %s\n", nc_strerror (status));
        return status;
      }

      // Define dimensions
      if ((status = nc_def_dim (ncid_phase, "nr", output->dist[0][0].n_dens, &dimids[0])) != NC_NOERR ||
          (status = nc_def_dim (ncid_phase, "nstokes", input.nstokes, &dimids[1])) != NC_NOERR ||
          (status = nc_def_dim (ncid_phase, "ntheta", output->numang, &dimids[2])) != NC_NOERR) {
        fprintf (stderr, "Error defining dimensions: %s\n", nc_strerror (status));
        return status;
      }

      // Define variable r
      if ((status = nc_def_var (ncid_phase, "r", NC_FLOAT, 1, &dimids[0], &varid_r)) != NC_NOERR) {
        fprintf (stderr, "Error defining variable: %s\n", nc_strerror (status));
        return status;
      }

      // Define variable theta
      if ((status = nc_def_var (ncid_phase, "theta", NC_FLOAT, 1, &dimids[2], &varid_theta)) != NC_NOERR) {
        fprintf (stderr, "Error defining variable: %s\n", nc_strerror (status));
        return status;
      }

      // Define variable phamat
      if ((status = nc_def_var (ncid_phase, "phamat", NC_FLOAT, 3, dimids, &varid_phamat)) != NC_NOERR) {
        fprintf (stderr, "Error defining variable: %s\n", nc_strerror (status));
        return status;
      }

      // End definition
      if ((status = nc_enddef (ncid_phase)) != NC_NOERR) {
        fprintf (stderr, "Error ending definition: %s\n", nc_strerror (status));
        return status;
      }

      if ((status = nc_put_var_float (ncid_phase, varid_theta, theta)) != NC_NOERR) {
        fprintf (stderr, "Error writing data: %s\n", nc_strerror (status));
        return status;
      }

      if ((status = nc_put_var_float (ncid_phase, varid_r, r)) != NC_NOERR) {
        fprintf (stderr, "Error writing data: %s\n", nc_strerror (status));
        return status;
      }
    }  // if (input.write_all_phasefunctions && input.phase_from_amplitudes) {

    // main loop
    for (ia = 0; ia < output->nalpha; ia++) {    /* loop over width of size distribution if specified in alpha_filename */
      for (ir = 0; ir < output->n_r_eff; ir++) { /* Loop over effective radii (cloud) / relative humidities (aerosol) */

        /* copy density of medium for given relative humidity to rho_medium which is used in write_output */
        if (input.distribution == DIST_AER)
          output->rho_medium = output->aerdistr.rho[ir];

        /* Allocate memory for optical properties for each size in the size distribution; */
        /* only once for cloud (because n_dens is independent of radius) and each time for aerosol*/
        if ((ir == 0 && ia == 0) || input.distribution == DIST_AER) {
          qext_array  = calloc (output->dist[ia][ir].n_dens, sizeof (float));
          qsca_array  = calloc (output->dist[ia][ir].n_dens, sizeof (float));
          gsca_array  = calloc (output->dist[ia][ir].n_dens, sizeof (float));
          qback_array = calloc (output->dist[ia][ir].n_dens, sizeof (float));

          pmom_array = calloc (input.nstokes, sizeof (float**));
          for (ip = 0; ip < input.nstokes; ip++) {
            pmom_array[ip] = calloc (input.mie.nmom + 1, sizeof (float*));
            for (im = 0; im < input.mie.nmom + 1; im++)
              pmom_array[ip][im] = calloc (output->dist[ia][ir].n_dens, sizeof (float));
          }
        }

        /***************************************************************/
        /* call Mie programs                                           */
        /***************************************************************/

        switch (input.distribution) {
        case DIST_NONE:
          /* one single radius */
          status = mie_calc (input.mie,
                             &(output->mie),
                             input.program,
                             input.medium,
                             input.water_properties,
                             crefin[ir],
                             wavelength,
                             output->r_eff[ir],
                             input.temperature,
                             input.nstokes,
                             output->numang,
                             &(output->ref));
          if (status != 0) {
            fprintf (stderr, "Error %d returned by mie_calc()\n", status);
            return status;
          }

          /* calculate beta and omega */
          ext           = M_PI * output->r_eff[ir] * output->r_eff[ir] * output->mie.qext * 1e-12;
          vol           = (4.0 / 3.0) * M_PI * pow ((output->r_eff[ir] * 1E-6), 3.0);
          output->beta  = ext / vol / 1000.0;
          output->omega = output->mie.qsca / output->mie.qext;
          output->back  = M_PI * output->r_eff[ir] * output->r_eff[ir] * output->mie.qback * 1e-12 / vol / 1000.0;

          /* Moments of phase matrix */
          if (input.mie.nmom > 0) {
            if (input.nstokes == 4) {
              /* Copy moments */
              for (ip = 0; ip < input.nstokes; ip++)
                for (im = 0; im <= input.mie.nmom; im++)
                  output->pmom[ip][im] = (output->mie.pmom)[ip][im];

              /*  Convert from Wiscombe definition to Evans definition of Mueller matrix */
              norm = output->mie.pmom[0][0] + output->mie.pmom[1][0];
              for (im = 0; im <= input.mie.nmom; im++) {
                output->pmom[0][im] = (output->mie.pmom[0][im] + output->mie.pmom[1][im]) / norm;
                output->pmom[1][im] = (output->mie.pmom[1][im] - output->mie.pmom[0][im]) / norm;
                output->pmom[2][im] = 2. * output->mie.pmom[2][im] / norm;
                output->pmom[3][im] = 2. * output->mie.pmom[3][im] / norm;
              }
            } else if (input.nstokes == 1) {
              norm = output->mie.pmom[0][0];
              for (im = 0; im <= input.mie.nmom; im++)
                output->pmom[0][im] = output->mie.pmom[0][im] / norm;
            } else {
              fprintf (stderr, "Error, nstokes must be either 1 or 4\n");
              return -1;
            }
          }

          break;
        case DIST_FILE:
        case DIST_GAMMA:
        case DIST_LOGNORMAL:
        case DIST_TRIANGULAR:
        case DIST_AER:

          /* droplet radius distribution */
          if ((ir == 0 && ia == 0) || input.distribution == DIST_AER) {
            status = mie_calc_all_radii (input.program,
                                         input.medium,
                                         input.water_properties,
                                         crefin[ir],
                                         wavelength,
                                         input.temperature,
                                         input.nstokes,
                                         output->numang,
                                         output->dist,
                                         ia,
                                         ir,
                                         output->nalpha,
                                         output->n_r_eff,
                                         input.mie,
                                         ncid_phase,
                                         varid_phamat,
                                         dimids,
                                         input.write_all_phasefunctions && input.phase_from_amplitudes,
                                         input.verbose,
                                         &(output->mie),
                                         &(output->ref),
                                         &qext_array,
                                         &qsca_array,
                                         &gsca_array,
                                         &qback_array,
                                         &pmom_array,
                                         output->phase_amp,
                                         output->phase_normalization);

            if (status != 0) {
              fprintf (stderr, "Error %d returned by mie_calc_all_radii()\n", status);
              return status;
            }
          }

          if (input.verbose) {
            fprintf (stderr, "Integration over size distribution i_reff %d i_alpha %d\n", ir, ia);
            fflush (stderr);
          }

          status = integrate_sizedist (qext_array,
                                       qsca_array,
                                       gsca_array,
                                       qback_array,
                                       pmom_array,
                                       input.nstokes,
                                       input.mie.nmom,
                                       input.mie.numang,
                                       output->dist[ia][ir].radius,
                                       output->dist[ia][ir].number_dens,
                                       output->dist[ia][ir].n_dens,
                                       input.verbose,
                                       &(output->beta),
                                       &(output->omega),
                                       &(output->g),
                                       &(output->back),
                                       &(output->pmom));

          if (status != 0) {
            fprintf (stderr, "Error %d returned by integrate_sizedist()\n", status);
            return status;
          }
          break;
        default:
          fprintf (stderr, "Error, distribution number in mie.c %d\n", input.distribution);
          return -1;

          break;
        }

        if (input.output_format != OUTPUT_FORMAT_NETCDF) {
          /* write output */
          status = write_output (input, output, ir, iv);
          if (status != 0) {
            fprintf (stderr, "Error %d returned by write_output()\n", status);
            return status;
          }
        }
#if HAVE_NETCDF4
        else {

          /* Find number of polynomials that are non-zero, this is done for P11.*/
          /*   For P12, P34 the */
          /* polynomials might oscillate around 0 so it is here not a good idea to cut the Legendre series */
          /* when the values are small */
          if (input.mie.nmom > 0) {
            output->nmommax = 0;
            for (int im = 0; im <= input.mie.nmom; im++) {
              if (fabs ((2 * im + 1) * (output->pmom)[0][im]) < 1E-5) {
                if (im > (output->nmommax))
                  (output->nmommax) = im;
                break;
              }
            }

            if (output->nmommax != 0)
              for (int ip = 0; ip < input.nstokes; ip++) {
                /* set values above nmommax to 0 */
                for (im = output->nmommax + 1; im <= input.mie.nmom; im++)
                  (output->pmom)[ip][im] = FLT_MAX;
              }

            /* remember to write a warning if there could be more
	       Legendre coefficients, only phase function is checked, the
	       other phase matrix elements not */
            if (output->nmommax == 0) {
              output->nmommax                 = input.mie.nmom;
              (output->warnings_nmom)[iv][ir] = 1;
            }
          }

          if (input.nthetamax > 0) {

            // phase function from amplitudes
            int      ntheta_tmp = output->numang;
            double** phase_tmp  = calloc (input.nstokes, sizeof (double*));
            for (int ip = 0; ip < input.nstokes; ip++)
              phase_tmp[ip] = calloc (ntheta_tmp, sizeof (double));

            double* xmu_tmp = calloc (ntheta_tmp, sizeof (double));

            for (int iang = 0; iang < ntheta_tmp; iang++) {
              xmu_tmp[iang] = output->mie.xmu[ntheta_tmp - 1 - iang];
              for (int ip = 0; ip < input.nstokes; ip++)
                phase_tmp[ip][iang] = output->phase_amp[ip][ia][ir][ntheta_tmp - 1 - iang];
            }

            status = calc_phase_function (  //input.mie.nmom, CE this should be nmommax??
              output->nmommax,
              wavelength,
              output->r_eff[ir],
              input.nstokes,
              input.nthetamax,
              input.phase_from_amplitudes,
              phase_tmp,
              xmu_tmp,
              ntheta_tmp,
              input.ntheta_user,
              input.accuracy,
              &(output->pmom),
              &(output->phase),
              &(output->theta),
              &(output->ntheta),
              &(output->warnings_ntheta)[iv][ir]);
            if (status != 0) {
              fprintf (stderr, "Error %d returned by calc_phase_function()\n", status);
              return status;
            }
          }

          // normalize phase function if calculated from amplitudes
          /* BM 30.4.2024: No need for that anymore because phase function is already normalized in mie_calc_all_radii() */
          /* If that is confirmed, make output->phase_normalization a local variable within mie_calc_all_radii() ??????? */
          /*
	  if (input.phase_from_amplitudes) {
	    double norm = 0;
	    for (int iang=0; iang<output->ntheta[0]-1; iang++)
	      norm += 0.5*(output->phase[0][iang] + output->phase[0][iang+1])*
		(cos(output->theta[0][iang+1]/180.0*M_PI)-cos(output->theta[0][iang]/180.0*M_PI));
	    
	    fprintf (stderr, "NORMALIZING phase function by %f\n", norm);
	    for (int ip=0; ip<input.nstokes; ip++)
	      for (int iang=0; iang<output->ntheta[ip]; iang++)
		output->phase[ip][iang] /= (norm/2.0);
	  }
	  */

          if (iv == 0 && ir == 0) {
            if (input.nmom_netcdf == 0) {
              nmommax          = input.mie.nmom;
              complete_moments = 1;
            } else {
              nmommax          = input.nmom_netcdf;
              complete_moments = 0;
            }

            status = create_netcdf_file (input.basename,
                                         input.medium,
                                         input.aerosol_type,
                                         input.distribution,
                                         output->alpha[ia],
                                         input.step,
                                         output->wl.nlambda,
                                         output->wl.lambda,
                                         output->n_r_eff,
                                         output->r_eff,
                                         nmommax,
                                         input.nstokes,
                                         input.nthetamax,
                                         output->aerdistr.n_hum,
                                         complete_moments,
                                         input.n_r_max,
                                         input.dx_max,
                                         input.sd_filename,
                                         &ncid,
                                         &(output->outputfile));
            if (status != 0) {
              fprintf (stderr, "Error %d returned by create_netcdf_file()\n", status);
              return status;
            }
          }

          write_output_netcdf (ncid, input, output, iv, ir, ia);
        }
#else
        fprintf (stderr, " ***********************************************************************\n");
        fprintf (stderr, " * You have built mie without libnetcdf and hence cannot               *\n");
        fprintf (stderr, " * use the option *output_user netcdf*. Please get netcdf and rebuild. *\n");
        fprintf (stderr, " ***********************************************************************\n");
        return -1;
#endif

        /* Free variables, for aerosol the number of points in the size
	   distribution is different for each relative humidity (ir),
	   therefore these variables are allocated for each ir.*/
        if ((ir == output->n_r_eff - 1 && ia == output->nalpha - 1) || input.distribution == DIST_AER) {
          for (ip = 0; ip < input.nstokes; ip++) {
            for (im = 0; im < input.mie.nmom + 1; im++)
              free (pmom_array[ip][im]);
            free (pmom_array[ip]);
          }
          free (pmom_array);

          free (qext_array);
          free (qsca_array);
          free (gsca_array);
        }

      }  // for (ir=0; ir<output->n_r_eff; ir++)
    }    // for (ia = 0; ia < output->nalpha; ia++)

    // Close the NetCDF phase functions file
    if (input.write_all_phasefunctions && input.phase_from_amplitudes) {
      if ((status = nc_close (ncid_phase)) != NC_NOERR) {
        fprintf (stderr, "Error closing NetCDF file: %s\n", nc_strerror (status));
        return status;
      }
    }

    for (ia = 0; ia < output->nalpha; ia++) {
      for (ir = 0; ir < output->n_r_eff; ir++) {
        free (output->dist[ia][ir].radius);
        free (output->dist[ia][ir].number_dens);
      }
      free (output->dist[ia]);
    }
    free (output->dist);

    for (ip = 0; ip < input.nstokes; ip++) {
      for (ia = 0; ia < output->nalpha; ia++) {
        for (ir = 0; ir < output->n_r_eff; ir++)
          free (output->phase_amp[ip][ia][ir]);
        free (output->phase_amp[ip][ia]);
      }
      free (output->phase_amp[ip]);
    }
    free (output->phase_amp);
  } /* loop over iv */

  /* write warnings*/
  if (input.output_format == OUTPUT_FORMAT_CLOUDPROP || input.output_format == OUTPUT_FORMAT_AERPROP ||
      input.output_format == OUTPUT_FORMAT_NETCDF) {

    for (iv = 0; iv < output->wl.nlambda; iv++) {
      for (ir = 0; ir < output->n_r_eff; ir++) {
        if (output->warnings_nmom[iv][ir] == 1) {
          if (firstwarn == 1) { /* true */
            fprintf (stderr, "*** Warning: maxmimum accuracy of the phase function not reached\n");
            fprintf (stderr, "*** with %4d Legendre coefficients.\n", input.mie.nmom);
            fprintf (stderr, "*** In order to get more accurate results\n");
            fprintf (stderr, "*** please increase 'nmom' in the mie input file!\n");
            firstwarn = 0;
          }
          if (input.medium != AEROSOL)
            fprintf (stderr, "          lambda = %10.3f nm, r_eff = %8.3f micro meter\n", output->wl.lambda[iv], output->r_eff[ir]);
          else
            fprintf (stderr,
                     "          lambda = %10.3f nm, rel_hum = %2f percent\n",
                     output->wl.lambda[iv],
                     output->aerdistr.rel_hum[ir]);
        }
      }
    }

    firstwarn = 1;
    for (iv = 0; iv < output->wl.nlambda; iv++) {
      for (ir = 0; ir < output->n_r_eff; ir++) {
        if (output->warnings_ntheta[iv][ir] == 1) {
          if (firstwarn == 1) {
            fprintf (stderr, "*** Warning: maxmimum accuracy of the phase function not reached\n");
            fprintf (stderr, "*** with %4d scattering angles.\n", input.nthetamax);
            fprintf (stderr, "*** In order to get more accurate results\n");
            fprintf (stderr, "*** please increase 'nthetamax' in the mie input file!\n");
            firstwarn = 0;
          }
          if (input.medium != AEROSOL)
            fprintf (stderr, "          lambda = %10.3f nm, r_eff = %8.3f micro meter\n", output->wl.lambda[iv], output->r_eff[ir]);
          else
            fprintf (stderr,
                     "          lambda = %10.3f nm, rel_hum = %2f percent\n",
                     output->wl.lambda[iv],
                     output->aerdistr.rel_hum[ir]);
        }
      }
    }
  }

  /*free memory */
  for (iv = 0; iv < output->wl.nlambda; iv++) {
    free (output->warnings_nmom[iv]);
    free (output->warnings_ntheta[iv]);
  }
  free (output->warnings_nmom);
  free (output->warnings_ntheta);

  if (input.output_format == OUTPUT_FORMAT_NETCDF) {

    free (output->outputfile);
    if (input.nthetamax > 0) {
      for (ip = 0; ip < input.nstokes; ip++) {
        free (output->phase[ip]);
        free (output->theta[ip]);
      }
      free (output->phase);
      free (output->theta);
      free (output->ntheta);
    }
  }

  free (output->mie.s1);
  free (output->mie.s2);
  free (output->mie.xmu);

  return 0;
}

/**
 * setup_reff_grid
 * 
 * Set up the effective radius grid for the Mie calculation from user
 * input parameters. 
 * 
 * @author      ..., Claudia Emde
 * @date        2010-02-11 put into subroutine by C. Emde
 */
int setup_reff_grid (/* Input */
                     int     distribution,
                     double  r_eff_step,
                     double  r_eff_min,
                     double  r_eff_max,
                     int     r_eff_log,
                     int     n_hum, /* nhum is used for aerosol*/
                     double* hum,
                     int     verbose,
                     /* Output */
                     int*     n_r_eff,
                     double** r_eff) {
  int status = 0;

  int ir = 0;

  /* determine the number of r_eff to be simulated */
  switch (distribution) {

  case DIST_NONE:
  case DIST_GAMMA:
  case DIST_LOGNORMAL:
  case DIST_TRIANGULAR:
    /* determine number of r_eff to be calculated */
    if (r_eff_step != 0.0) {
      if (!r_eff_log) /* 1e-7 avoids rouding errors */
        *n_r_eff = (int) floor ((r_eff_max - r_eff_min + 1e-7) / r_eff_step) + 1;
      else {
        if (r_eff_step <= 1) {
          fprintf (stderr, "Error, effective radius increment needs to be larger\n");
          fprintf (stderr, "than 1 for log steps.\n");
          return -1;
        }
        *n_r_eff = (int) (ceil (log (r_eff_max / r_eff_min) / log (r_eff_step))) + 1;
      }
    } else
      *n_r_eff = 1;

    /* setup r_eff grid */
    *r_eff = calloc (*n_r_eff, sizeof (double));
    if (verbose)
      fprintf (stderr, "... effective radii: \n");

    (*r_eff)[0] = r_eff_min;
    for (ir = 1; ir < *n_r_eff; ir++) {
      if (!r_eff_log)
        (*r_eff)[ir] = r_eff_min + (double) ir * r_eff_step;
      else
        (*r_eff)[ir] = (*r_eff)[ir - 1] * r_eff_step;

      if (verbose)
        fprintf (stderr, "%g ", (*r_eff)[ir]);
    }

    if (verbose)
      fprintf (stderr, "\n");

    break;

  case DIST_FILE:
    /* size distribution file means only one r_eff */
    *n_r_eff    = 1;
    *r_eff      = calloc (*n_r_eff, sizeof (double));
    (*r_eff)[0] = -999.0; /* effective radius value is not required */
    break;

  case DIST_AER:
    /* Effective radius is "replaced" by relative humidity for aerosols */
    *n_r_eff = n_hum;
    *r_eff   = calloc (*n_r_eff, sizeof (double));
    for (ir = 0; ir < n_hum; ir++)
      (*r_eff)[ir] = hum[ir];
    fprintf (stderr, "aerosol nhum %d\n", *n_r_eff);

    break;
  default:
    fprintf (stderr, "Error, distribution number in mie.c %d\n", distribution);
    return -3;
    break;
  }
  return status;
}

/**
 * read_refractive_index
 *
 * Read the refractive index from a file
 *
 * @author ...., Claudia Emde
 * @date   2010-02-11 restructured by C. Emde 
 */
int read_refractive_index (/* Input */
                           char* filename,
                           int   verbose,
                           /* Output */
                           int*      nlambda,
                           double**  lambda,
                           double*** ref_real,
                           double*** ref_imag) {
  int status = 0;

  int iv = 0;

  /* Here we have only one refractive index for each wavelength, this
     is different for aerosols. */
  *ref_real = calloc (1, sizeof (double*));
  *ref_imag = calloc (1, sizeof (double*));
  status    = read_3c_file (filename, &(*lambda), &(*ref_real)[0], &(*ref_imag)[0], &(*nlambda));

  if (status != 0) {
    fprintf (stderr, "error %d reading %s\n", status, filename);
    return status;
  }

  fprintf (stderr, " ... read %d data points from %s\n", *nlambda, filename);

  /* check if data points are sorted in ascending order */
  for (iv = 0; iv < *nlambda - 1; iv++) {
    if ((*lambda)[iv + 1] < (*lambda)[iv]) {
      fprintf (stderr, "Error: %s not sorted in ascending order\n", filename);
      return -1;
    }
  }

  return 0;
}

/**
 * setup_wavelength_grid
 *
 * Setup the wavelength grid for the Mie calculation from user
 * specified input parameters, read wavelengths from file, or specify
 * wavelengths for k-distributions.
 *
 * @author ...., Claudia Emde
 * @date   2010-02-11 restructured by C. Emde
 * 
 */

int setup_wavelength_grid (/* Input */
                           double wl_start,
                           double wl_end,
                           double wl_step,
                           char*  data_files_path,
                           char*  filename_ck_generic,
                           int    ck_scheme,
                           int    lambda_i_start,
                           int    lambda_i_end,
                           int    verbose,
                           /* Output */
                           int*     nlambda,
                           double** lambda) {
  int   status       = 0;
  int   iv           = 0;
  char* wvl_filename = NULL;

  int min_columns = 0;
  int max_columns = 0;
  int max_length  = 0;

  char*** string = NULL;
  char*   dummy  = NULL;

  int unit_factor = 1; /* specifies if wvl contain micro meter or nm */

  char function_name[] = "setup_wavelength_grid";
  char file_name[]     = "mie.c";

  if (ck_scheme == CK_CRS) {
    if (wl_step != 0.0)
      *nlambda = (int) (wl_end - wl_start) / wl_step + 1;
    else
      *nlambda = 1;

    *lambda = calloc (*nlambda, sizeof (double));

    if (verbose)
      fprintf (stderr, "... wavelengths:\n");
    for (iv = 0; iv < *nlambda; iv++) {
      (*lambda)[iv] = wl_start + wl_step * iv;
      if (verbose)
        fprintf (stderr, "%g ", (*lambda)[iv]);
    }
    if (verbose)
      fprintf (stderr, "\n");
  }

  else {

    wvl_filename = calloc (FILENAME_MAX, sizeof (char));

    strcpy (wvl_filename, data_files_path);

    /* Note: For the correlated-k schemes the optical properties are
       calculated for the center wavelengths. Evans' cloudprop does a
       weighting over the Planck functions (Sun, Earth). */
    switch (ck_scheme) {
    case CK_KATO:
      strcat (wvl_filename, "correlated_k/kato/wvl.dat");
      fprintf (stderr, " ... Kato wavelength bands\n");
      unit_factor = 1; /* Kato wvl.dat contains nano meter */
      break;
    case CK_KATO2:
      strcat (wvl_filename, "correlated_k/kato2/wvl.dat");
      fprintf (stderr, " ... Kato2 wavelength bands\n");
      unit_factor = 1; /* Kato2 wvl.dat contains nano meter */
      break;
    case CK_KATO2ANDWANDJI:
      strcat (wvl_filename, "correlated_k/kato2andwandji/wvl.dat");
      fprintf (stderr, " ... Kato2andwandji wavelength bands\n");
      unit_factor = 1; /* Kato2 wvl.dat contains nano meter */
      break;
    case CK_KATO2_96:
      strcat (wvl_filename, "correlated_k/kato2.hitran96/wvl.dat");
      fprintf (stderr, " ... Kato2_96 wavelength bands\n");
      unit_factor = 1; /* Kato2_96 wvl.dat contains nano meter */
      break;
    case CK_FU:
      strcat (wvl_filename, "correlated_k/fu/wvl.dat");
      fprintf (stderr, " ... Fu wavelength bands\n");
      unit_factor = 1000; /* Fu wvl.dat contains micro meter */
      break;
    case CK_AVHRR_KRATZ:
      strcat (wvl_filename, "correlated_k/kratz/avhrr_wvl.dat");
      fprintf (stderr, " ... AVHRR-Kratz wavelength bands\n");
      unit_factor = 1000; /* AVHRR-Kratz wvl.dat contains micro meter */
      break;
    case CK_GENERIC:
      strcpy (wvl_filename, filename_ck_generic);
      fprintf (stderr, " ... generic correlated_k wavelength bands\n");
      break;
    default:
      fprintf (stderr, "Error: unsupported correlated-k scheme, nr %2d\n", ck_scheme);
      return -1;
    }

    status = ASCII_checkfile (wvl_filename, &(*nlambda), &min_columns, &max_columns, &max_length);
    if (status != 0) {
      fprintf (stderr, "Error %d reading wavelength_file %s in %s (%s) \n", status, wvl_filename, function_name, file_name);
      return status;
    }

    if (min_columns < 2) {
      fprintf (stderr, "Error, wavelength_file %s must have at least two columns! \n", wvl_filename);
      return -1;
    }

    if ((status = ASCII_calloc_string (&string, *nlambda, max_columns, FILENAME_MAX + 1)) != 0) {
      fprintf (stderr, "Error %d allocating string in %s (%s)\n", status, function_name, file_name);
      return status;
    }

    if ((status = ASCII_readfile (wvl_filename, string)) != 0) {
      fprintf (stderr, "Error %d reading wavelength-file: %s in %s (%s) \n", status, wvl_filename, function_name, file_name);
      return status;
    }

    if ((*lambda = calloc (*nlambda, sizeof (double))) == NULL) {
      fprintf (stderr, "Error %d allocating output->wl.lambda in %s (%s) \n", status, function_name, file_name);
      return status;
    }

    /* second column is center wavelength */
    for (iv = 0; iv < *nlambda; iv++) {
      (*lambda)[iv] = (double) strtod (string[iv][1], &dummy) * unit_factor;
    }

    if (verbose)
      for (iv = 0; iv < *nlambda; iv++)
        fprintf (stderr, "   lambda[%3d] = %10.3f nm\n", iv + 1, (*lambda)[iv]);

    /* Check validity of specified wavelength indices */
    if (lambda_i_start != -999 && lambda_i_start < 1) {
      fprintf (stderr, "  wavelength_index_start = %d has to be 1 or larger \n", lambda_i_start);
      status = -1;
    }
    if (*nlambda < lambda_i_start - 1) {
      fprintf (stderr, "  wavelength_index_start = %d, but there are only %d wavelength bands \n", lambda_i_start, *nlambda);
      status = -2;
    }

    if (lambda_i_end != -999 && lambda_i_end < lambda_i_start) {
      fprintf (stderr,
               "  wavelength_index_end = %d has to be larger than wavelength_index_start = %d\n",
               lambda_i_end,
               lambda_i_start);
      status = -3;
    }
    if (*nlambda < lambda_i_end) {
      fprintf (stderr, "  wavelength_index_end = %d, but there are only %d wavelength bands \n", lambda_i_end, *nlambda);
      status = -4;
    }
  }
  return status;
}

/**
 * calc_distribution_radius_grid
 * 
 * Calculate a radius grid for the size distribution over which the
 * optical properties are intgrated. This grid must be sufficiently
 * fine, otherwise "ripples" in the phase function are not averaged
 * out. Here the sampling is done in the same way as in cloudprp by
 * F. Evans.
 *
 * @author ???
 * @date   ???
 */
int calc_distribution_radius_grid (/* Input */
                                   int     distribution,
                                   char*   sd_filename,
                                   double  r_eff_min,
                                   double  r_eff_max,
                                   double* reff,
                                   int     nreff,
                                   double  alpha,
                                   int     nalpha,
                                   double  step,
                                   double  lambda, /*wavelength in micrometer*/
                                   double  dx_max,
                                   double  n_r_max,
                                   int     verbose,
                                   /* Output */
                                   dist_struct** dist) {
  double  radmin = 0.0, radmax = 0.0, rad = 0.0, delta_rad = 0.0;
  double  x = 0.0, delta_x = 0.0, dx = 0.0;
  int     is     = 0;
  int     status = 0, counter = 0, nint = 0;
  double* tmprad = NULL;
  int     nsteps = 0;

  double *radius_tmp = NULL, *dens_tmp = NULL;
  int     n_tmp = 0;

  switch (distribution) {

  case DIST_NONE: /* single radius calculation, nothing to do */

    // probably not needed but since we free everything at the end we should also allocate it */
    for (int ia = 0; ia < nalpha; ia++)
      for (int ir = 0; ir < nreff; ir++) {
        dist[ia][ir].n_dens      = 1;
        dist[ia][ir].radius      = calloc (1, sizeof (double));
        dist[ia][ir].number_dens = calloc (1, sizeof (double));

        dist[ia][ir].radius[0]      = 1;
        dist[ia][ir].number_dens[0] = 1;
      }
    break;

  case DIST_FILE:
    /* read size distribution from file */
    status = read_2c_file (sd_filename, &radius_tmp, &dens_tmp, &n_tmp);

    if (status != 0) {
      fprintf (stderr, "error %d reading %s\n", status, sd_filename);
      return status;
    }

    if (verbose)
      fprintf (stderr, " ... read %d data points from %s\n", n_tmp, sd_filename);

    for (int ia = 0; ia < nalpha; ia++)
      for (int ir = 0; ir < nreff; ir++) {
        dist[ia][ir].n_dens      = n_tmp;
        dist[ia][ir].radius      = calloc (n_tmp, sizeof (double));
        dist[ia][ir].number_dens = calloc (n_tmp, sizeof (double));

        for (int is = 0; is < n_tmp; is++) {
          dist[ia][ir].radius[is]      = radius_tmp[is];
          dist[ia][ir].number_dens[is] = dens_tmp[is];
        }
      }

    break;
  case DIST_GAMMA:
  case DIST_LOGNORMAL:
    /* ????? are these reasonable values? Evans used 0.02 reff_min - 5 reff_max; ????? */
    /* ????? BM found that 5 is not enough for small particles - but somebody    ????? */
    /* ????? increased the upper boundary to 8 reff_max which is VERY EXPENSIVE  ????? */

    radmin = 0.02 * r_eff_min;    /* cover the whole droplet spectrum   */
    radmax = n_r_max * r_eff_max; /* in order to save mie-calculations  */

    /* CE used n_r_max = 5 to calculate new mie data */
    /* radmax = 5.00 * input.r_eff_max;*/

    /* CE ????? */
    /* dx is the maximum bin size (size parameter) in size distribution grid */
    /* Original value is 0.03, this seems to be unsufficient, especially */
    /* for backward scattering region, also in rainbow region lots of wiggles. */
    /* May be we should think of another way to sample the size distribution */
    dx = dx_max;

    /* count number of radii */
    rad = radmin;
    is  = 0;
    while (rad <= radmax) {
      x         = 2 * M_PI * rad / lambda;
      delta_x   = double_max (dx, dx * sqrt (x));
      delta_rad = delta_x * lambda / (2 * M_PI);
      rad       = rad + delta_rad;
      is        = is + 1;
    }
    n_tmp = is;

    if (verbose)
      fprintf (stderr, "size distributions: r_min %g um, r_max %g um, number of grid points %d \n", radmin, radmax, n_tmp);

    for (int ia = 0; ia < nalpha; ia++)
      for (int ir = 0; ir < nreff; ir++) {
        dist[ia][ir].n_dens      = n_tmp;
        dist[ia][ir].radius      = calloc (n_tmp, sizeof (double));
        dist[ia][ir].number_dens = calloc (n_tmp, sizeof (double));
      }

    /* calculate radius grid */
    rad = radmin;
    is  = 0;
    while (rad <= radmax) {

      for (int ia = 0; ia < nalpha; ia++)
        for (int ir = 0; ir < nreff; ir++)
          dist[ia][ir].radius[is] = rad;

      x         = 2 * M_PI * rad / lambda;
      delta_x   = double_max (dx, dx * sqrt (x));
      delta_rad = delta_x * lambda / (2 * M_PI);
      rad       = rad + delta_rad;
      is        = is + 1;
    }
    break;

  case DIST_TRIANGULAR:
    nint   = (int) (alpha / step);
    nsteps = 2 * nint + 1;
    tmprad = calloc (nreff * nsteps, sizeof (double));
    for (int ireff = 0; ireff < nreff; ireff++)
      for (int jreff = -nint; jreff <= nint; jreff++)
        tmprad[counter++] = reff[ireff] + (double) jreff * step;

    // remove duplicates by sorting the array and keeping only those that differ by more than EPSILON
    qsort (tmprad, nreff * nsteps, sizeof (double), sortasc);

    counter = 1;
    for (int i = 0; i < nreff * nsteps - 1; i++)
      if (fabs (tmprad[i + 1] - tmprad[i]) > 10.0 * DBL_EPSILON * tmprad[i + 1])
        counter++;

    n_tmp = counter;
    fprintf (stderr, " ... found %d grid points for triangular calculation\n", n_tmp);
    fflush (stderr);

    for (int ia = 0; ia < nalpha; ia++)
      for (int ir = 0; ir < nreff; ir++) {
        dist[ia][ir].n_dens      = n_tmp;
        dist[ia][ir].radius      = calloc (n_tmp, sizeof (double));
        dist[ia][ir].number_dens = calloc (n_tmp, sizeof (double));
      }

    for (int ia = 0; ia < nalpha; ia++)
      for (int ir = 0; ir < nreff; ir++)
        dist[ia][ir].radius[0] = tmprad[0];

    counter = 0;
    for (int i = 0; i < nreff * nsteps - 1; i++)
      if (fabs (tmprad[i + 1] - tmprad[i]) > 10.0 * DBL_EPSILON * tmprad[i + 1]) {
        for (int ia = 0; ia < nalpha; ia++)
          for (int ir = 0; ir < nreff; ir++)
            dist[ia][ir].radius[++counter] = tmprad[i + 1];
      }
    free (tmprad);

    for (int i = 0; i < n_tmp; i++)
      fprintf (stderr, "%f\n", dist[0][0].radius[i]);
    fflush (stderr);

    break;
  case DIST_AER:
    /* Do nothing here. Radius grid and number densities are calculated in
       aerosol_size_distribution. */
    break;
  default:
    fprintf (stderr, "Error in calc_distribution_radius_grid: unrecognized distribution %d\n", distribution);
    return -1;
    break;
  }

  /* check if data points are sorted in ascending order */
  switch (distribution) {
  case DIST_NONE:
    break;
  case DIST_FILE:
  case DIST_GAMMA:
  case DIST_LOGNORMAL:
  case DIST_TRIANGULAR:
    /* for gamma and lognormal this is not nessesary, but it costs almost no time */
    for (int ia = 0; ia < nalpha; ia++)
      for (int ir = 0; ir < nreff; ir++) {
        for (is = 0; is < dist[ia][ir].n_dens - 1; is++)
          if (dist[ia][ir].radius[is + 1] < dist[ia][ir].radius[is]) {
            if (distribution == DIST_FILE)
              fprintf (stderr, "Error: %s not sorted in ascending order\n", sd_filename);
            else
              fprintf (stderr, "Error: radius in distribution not sorted in ascending order, internal program bug!!!\n");
            return -5;
          }
      }

    for (int ia = 0; ia < nalpha; ia++)
      for (int ir = 0; ir < nreff; ir++)
        if (dist[ia][ir].radius[0] == 0) {
          dist[ia][ir].radius[0] = 0.1 * dist[ia][ir].radius[1];
          fprintf (stderr, " ... cannot handle r=0; setting first data point to r=%g\n", dist[ia][ir].radius[0]);
        }

    break;
  default:
    fprintf (stderr, "Error in calc_distribution_radius_grid: unrecognized distribution %d\n", distribution);
    break;
  }

  return status;
}

/**
 * calc_distribution_number_dens
 *
 * Calculate particle number densities for each radius in the size
 * distribution.
 *
 * @author ???
 * @date   ???
 */
int calc_distribution_number_dens (/* Input */
                                   int     distribution,
                                   double  alpha,
                                   double  r_eff,
                                   double  rho_medium,
                                   int     n_dens,
                                   double* radius,
                                   int     verbose,
                                   /* Output */
                                   double* number_dens) {
  double A = 0.0, B = 0.0;
  float  LWC    = 1.0;
  int    is     = 0;
  int    status = 0;
  double lnA    = 0;

  /* fprintf(stderr,"calculate droplet radius distribution\n"); */

  if (rho_medium <= 0.0) {
    fprintf (stderr, "Error, rho_medium = %f is 0 or negative!\n", rho_medium);
    return -1;
  }

  switch (distribution) {
  case DIST_NONE:
  case DIST_FILE:
    /* single radius calculation, nothing to do */
    break;

  case DIST_GAMMA:

    /* if (verbose) */
    /*       fprintf(stderr, "Gamma size distribution, alpha=%f \n", alpha); */

    B = (alpha + 3.0) / r_eff;
    /* original version by Claudia Emde */
    /* A = (0.75/M_PI)*pow(B,(alpha+4))*1000./exp(gamma_ln(alpha+4.)); */
    /* A = A*LWC/rho_medium;                                         */

    /* new version by Bernhard Mayer, 31.8.2016 */
    lnA = log (0.75 / M_PI * 1000.0 * LWC / rho_medium) + (alpha + 4) * log (B) - gamma_ln (alpha + 4.);

    /* calculate number density */
    for (is = 0; is < n_dens; is++) {
      /* replaced gamma function by exp (ln(gamma)) since otherwise we often get a numerical overflow or NAN.         */
      /* instead of multiplying a very large number by a very small number, we calculate the logarithm and then exp() */
      /* original version by Claudia Emde */
      // number_dens[is] = A * pow(radius[is], alpha) * exp(-B*radius[is]);

      /* new version by Bernhard Mayer, 31.8.2016 */
      number_dens[is] = exp (lnA + alpha * log (radius[is]) - B * radius[is]);

      /* if(verbose) */
      /* 	fprintf(stderr," %3d %g %g \n", is, radius[is], number_dens[is]); */
    }
    break;

  case DIST_LOGNORMAL:

    if (verbose)
      fprintf (stderr, "log-normal distribution, alpha=%f, r_mod=%f \n", alpha, r_eff);

    /* specified r_eff is taken as r_mod for log normal distribution */
    for (is = 0; is < n_dens; is++) {
      A               = (log ((radius)[is]) - log (r_eff)) / log (alpha);
      number_dens[is] = 1.0 / (sqrt (2.0 * M_PI) * log (alpha) * radius[is]) * exp (-0.5 * A * A);
    }

    /* old implementation:  */
    /*   CE: this was not consistent with the commonly used log-normal size distribution as for  */
    /*       instance in OPAC (Hess et al. 1998). Now made implementation consistent with OPAC. */
    /* B = r_eff*exp(-2.5*alpha*alpha); */
    /*     A = 1000./((4./3.*M_PI)* sqrt(2*M_PI)*alpha * B*B*B *exp(4.5*alpha*alpha)); */

    /*     A = A*LWC/rho_medium; */

    /* calculate number density */
    /* for (is=0;is<n_dens;is++){ */
    /*       number_dens[is] = A/radius[is] *  */
    /*         exp(-0.5* pow(log(radius[is]/B),2) / pow(alpha,2)); */
    /*       /\* if(verbose ) *\/ */
    /*         fprintf(stderr," %3d %g %g \n", is, radius[is], number_dens[is]);    */
    /*     } */
    break;
  case DIST_TRIANGULAR:
    for (is = 0; is < n_dens; is++) {
      if (fabs (radius[is] - r_eff) < alpha)
        number_dens[is] = 1.0 - fabs (radius[is] - r_eff) / alpha;
      else
        number_dens[is] = 0;
    }

    for (is = 0; is < n_dens; is++)
      fprintf (stderr, "nd %f %f\n", radius[is], number_dens[is]);

    fflush (stderr);

    break;
  default:
    fprintf (stderr, "Error in make_distribution : unrecognized distribution %d\n", distribution);
    status = -1;
    break;
  } /* switch(distribution) */

  return status;
}

double gamma_ln (double xx) {
  int     j = 0;
  double* cof;
  double  stp  = 2.50662827465;
  double  half = 0.5, one = 1.0, fpf = 5.5, x, tmp, ser;
  double  gammaln = 0.0;

  cof    = calloc (6, sizeof (double));
  cof[0] = 76.18009173;
  cof[1] = -86.50532033;
  cof[2] = 24.01409822;
  cof[3] = -1.231739516;
  cof[4] = 0.120858003;
  cof[5] = -0.53638200;

  x   = xx - one;
  tmp = x + fpf;
  tmp = (x + half) * log (tmp) - tmp;
  ser = one;
  for (j = 0; j <= 5; j++) {
    x   = x + one;
    ser = ser + cof[j] / x;
  }
  gammaln = tmp + log (stp * ser);

  free (cof);

  return gammaln;
}

double double_min (double x, double y) {
  if (x <= y)
    return x;
  else
    return y;
}

double double_max (double x, double y) {
  if (x >= y)
    return x;
  else
    return y;
}

/**
 * mie_calc_all_radii
 * 
 * Perform Mie calculations for all radii of all size
 * distributions. The resulting optical properties are afterwards
 * integrated over the size distributions. This function needs to be
 * called only once for each wavelength and optical properties for all
 * specified effective radii can be obtained by integrating the
 * resulting optical properties arrays over the various size
 * distributions. 
 *
 * @author ...
 * @date   2010-02-11 restructured by C. Emde
 */
int mie_calc_all_radii (/* Input */
                        int            program,
                        int            medium,
                        int            water_properties,
                        mie_complex    crefin,
                        double         wavelength,
                        float          temperature,
                        int            nstokes,
                        int            numang,
                        dist_struct**  dist,
                        int            ia,
                        int            ir,
                        int            nalpha,
                        int            n_r_eff,
                        mie_inp_struct input,
                        int            ncid_phase,
                        int            varid,
                        int            dimids[3],
                        int            write_all_phasefunctions,
                        int            verbose,
                        /* Output */
                        mie_out_struct* output,
                        mie_complex*    ref,
                        float**         qext_array,
                        float**         qsca_array,
                        float**         gsca_array,
                        float**         qback_array,
                        float****       pmom_array,
                        double****      phase_amp,
                        double**        phase_normalization) {
  int    status  = 0;
  double p11_tmp = 0, p12_tmp = 0, p33_tmp = 0, p43_tmp = 0;
  float  phase_tmp[nstokes][numang];  // for experimental output to netcdf file

  // for aerosol because the size distribution changes with ia and ir
  for (int is = 0; is < dist[ia][ir].n_dens; is++) {

    if (is % 100 == 0 && verbose) {
      fprintf (stderr, "%d ", is);
      fflush (stderr);
    }

    status = mie_calc (input,
                       output,
                       program,
                       medium,
                       water_properties,
                       crefin,
                       wavelength,
                       dist[ia][ir].radius[is],
                       temperature,
                       nstokes,
                       numang,
                       ref);
    if (status != 0) {
      fprintf (stderr, "Error %d returned by mie_calc()\n", status);
      return status;
    }

    /* save results for further calculations */
    (*qext_array)[is]  = output->qext;
    (*qsca_array)[is]  = output->qsca;
    (*gsca_array)[is]  = output->gsca;
    (*qback_array)[is] = output->qback;

    if (!(output->qext > 0) && verbose) {
      fprintf (stderr, "Warning, something strange happened in Mie calculation: \n");
      fprintf (stderr, "is %d wavelength %g x_size %g qext %g \n", is, wavelength, dist[ia][ir].radius[is], output->qext);
    }

    if (input.nmom > 0)
      for (int ip = 0; ip < nstokes; ip++)
        for (int im = 0; im <= input.nmom; im++)
          (*pmom_array)[ip][im][is] = output->pmom[ip][im];

    if (numang > 1) {
      for (int iang = 0; iang < numang; iang++) {
        // we need to integrate on the fly to avoid arrays of dimension numang * n_size
        // which may have up to 1E9 elements and more for large particles

        phase_tmp[0][iang] = 0.5 * ((output->s2[iang].re * output->s2[iang].re + output->s2[iang].im * output->s2[iang].im) +
                                    (output->s1[iang].re * output->s1[iang].re + output->s1[iang].im * output->s1[iang].im));

        phase_tmp[1][iang] = 0.5 * ((output->s2[iang].re * output->s2[iang].re + output->s2[iang].im * output->s2[iang].im) -
                                    (output->s1[iang].re * output->s1[iang].re + output->s1[iang].im * output->s1[iang].im));

        phase_tmp[2][iang] = (output->s1[iang].re * output->s2[iang].re + output->s1[iang].im * output->s2[iang].im);

        phase_tmp[3][iang] = -(output->s1[iang].re * output->s2[iang].im - output->s2[iang].re * output->s1[iang].im);

        for (int ia = 0; ia < nalpha; ia++) {
          for (int ir = 0; ir < n_r_eff; ir++) {

            p11_tmp = phase_tmp[0][iang] * dist[ia][ir].number_dens[is];
            p12_tmp = phase_tmp[1][iang] * dist[ia][ir].number_dens[is];
            p33_tmp = phase_tmp[2][iang] * dist[ia][ir].number_dens[is];
            p43_tmp = phase_tmp[3][iang] * dist[ia][ir].number_dens[is];

            int     n   = dist[ia][ir].n_dens;
            double* x   = dist[ia][ir].radius;
            double  d2x = 0;
            double  X   = (2 * M_PI * dist[ia][ir].radius[is] / wavelength);  // size parameter

            if (is > 0 && is < n - 1)
              d2x = dist[ia][ir].radius[is + 1] - dist[ia][ir].radius[is - 1];

            if (is == 0) {
              phase_amp[0][ia][ir][iang] = 0.5 * p11_tmp * (x[1] - x[0]);
              if (iang == 0)
                phase_normalization[ia][ir] = X * X * output->qsca / 4.0 * dist[ia][ir].number_dens[is] * (x[1] - x[0]);
              if (nstokes == 4) {
                phase_amp[1][ia][ir][iang] = 0.5 * p12_tmp * (x[1] - x[0]);
                phase_amp[2][ia][ir][iang] = 0.5 * p33_tmp * (x[1] - x[0]);
                phase_amp[3][ia][ir][iang] = 0.5 * p43_tmp * (x[1] - x[0]);
              }
            }

            if (is > 0 && is < n - 1) {
              phase_amp[0][ia][ir][iang] += 0.5 * p11_tmp * d2x;

              if (iang == 0)
                phase_normalization[ia][ir] += X * X * output->qsca / 4.0 * dist[ia][ir].number_dens[is] * d2x;

              if (nstokes == 4) {
                phase_amp[1][ia][ir][iang] += 0.5 * p12_tmp * d2x;
                phase_amp[2][ia][ir][iang] += 0.5 * p33_tmp * d2x;
                phase_amp[3][ia][ir][iang] += 0.5 * p43_tmp * d2x;
              }
            }

            if (is == n - 1) {
              phase_amp[0][ia][ir][iang] += 0.5 * p11_tmp * (x[n - 1] - x[n - 2]);

              if (iang == 0)
                phase_normalization[ia][ir] += X * X * output->qsca / 4.0 * dist[ia][ir].number_dens[is] * (x[n - 1] - x[n - 2]);

              if (nstokes == 4) {
                phase_amp[1][ia][ir][iang] += 0.5 * p12_tmp * (x[n - 1] - x[n - 2]);
                phase_amp[2][ia][ir][iang] += 0.5 * p33_tmp * (x[n - 1] - x[n - 2]);
                phase_amp[3][ia][ir][iang] += 0.5 * p43_tmp * (x[n - 1] - x[n - 2]);
              }
            }
          }
        }  // for (int is = 0; is < dist[ia][ir].n_dens; is++) { ...
      }

      // Write the data
      if (write_all_phasefunctions) {
        size_t start[3] = {0, 0, 0};             // Start from the beginning of the array
        size_t count[3] = {1, nstokes, numang};  // Write one "slice" at a time
        start[0]        = is;
        if ((status = nc_put_vara_float (ncid_phase, varid, start, count, &phase_tmp[0][0])) != NC_NOERR) {
          fprintf (stderr, "Error writing data: %s\n", nc_strerror (status));
          return status;
        }
      }
    }  //  for (int iang = 0; iang < numang; iang++)
  }

  // normalize to integral over size distribution
  if (numang > 1)
    for (ia = 0; ia < nalpha; ia++) 
      for (ir = 0; ir < n_r_eff; ir++) 
        for (int ip = 0; ip < nstokes; ip++)
          for (int iang = 0; iang < numang; iang++)
            phase_amp[ip][ia][ir][iang] /= (0.5 * phase_normalization[ia][ir]);

  if (verbose)
    fprintf (stderr, "\n");

  return status;
}

/**
 * integrate_sizedist
 * 
 * Integrate optical properties arrays (output of mie_calc_all_radii)
 * over a size distribution. 
 * 
 * @author ...
 * @date   2010-02-11 restructured by C.Emde
 */
int integrate_sizedist (/* Input */
                        float*   qext_array,
                        float*   qsca_array,
                        float*   gsca_array,
                        float*   qback_array,
                        float*** pmom_array,
                        int      nstokes,
                        int      nmom,
                        int      numang,
                        double*  x_size,
                        double*  y_size,
                        int      n_size,
                        int      verbose,
                        /* Output */
                        double*  beta,
                        double*  omega,
                        double*  g,
                        double*  back,
                        float*** pmom) {

  double xsquared = 0, xcubed = 0, norm = 0;
  double ext = 0, sca = 0, asy = 0, vol = 0, backsca = 0;

  /* allocate memory for fields */
  double* extinction     = calloc (n_size, sizeof (double));
  double* scattering     = calloc (n_size, sizeof (double));
  double* asymmetry      = calloc (n_size, sizeof (double));
  double* volume         = calloc (n_size, sizeof (double));
  double* backscattering = calloc (n_size, sizeof (double));

  double*** pmom_tmp = NULL;

  if (nmom > 0) {

    pmom_tmp = calloc (nstokes, sizeof (double**));
    for (int ip = 0; ip < nstokes; ip++) {
      pmom_tmp[ip] = calloc (nmom + 1, sizeof (double*));
      for (int im = 0; im <= nmom; im++)
        pmom_tmp[ip][im] = calloc (n_size, sizeof (double));
    }
  }

  for (int is = 0; is < n_size; is++) {
    xsquared = x_size[is] * x_size[is] * 1e-12;
    xcubed   = xsquared * x_size[is] * 1e-6;

    extinction[is]     = qext_array[is] * xsquared * y_size[is];
    scattering[is]     = qsca_array[is] * xsquared * y_size[is];
    backscattering[is] = qback_array[is] * xsquared * y_size[is];
    asymmetry[is]      = gsca_array[is] * scattering[is];
    volume[is]         = xcubed * y_size[is];

    if (nmom > 0)
      for (int ip = 0; ip < nstokes; ip++)
        for (int im = 0; im <= nmom; im++)
          pmom_tmp[ip][im][is] = pmom_array[ip][im][is] * y_size[is];
  }

  norm = integrate (x_size, y_size, n_size);

  ext     = integrate (x_size, extinction, n_size) / norm;
  sca     = integrate (x_size, scattering, n_size) / norm;
  backsca = integrate (x_size, backscattering, n_size) / norm;
  asy     = integrate (x_size, asymmetry, n_size) / norm;
  vol     = integrate (x_size, volume, n_size) / norm;

  if (nmom > 0) {
    for (int ip = 0; ip < nstokes; ip++)
      for (int im = 0; im <= nmom; im++)
        (*pmom)[ip][im] = integrate (x_size, pmom_tmp[ip][im], n_size) / norm;

    if (nstokes == 4) {
      /* Copy moments */
      for (int ip = 0; ip < nstokes; ip++)
        for (int im = 0; im <= nmom; im++)
          pmom_tmp[ip][im][0] = (*pmom)[ip][im];

      /*  Convert from Wiscombe definition to Evans definition of Mueller matrix */
      norm = pmom_tmp[0][0][0] + pmom_tmp[1][0][0];
      for (int im = 0; im <= nmom; im++) {
        (*pmom)[0][im] = (pmom_tmp[0][im][0] + pmom_tmp[1][im][0]) / norm;
        (*pmom)[1][im] = (pmom_tmp[1][im][0] - pmom_tmp[0][im][0]) / norm;
        (*pmom)[2][im] = 2. * pmom_tmp[2][im][0] / norm;  // why factor 2? Also needed for phase_from_amplitudes!
        (*pmom)[3][im] = 2. * pmom_tmp[3][im][0] / norm;
      }
    } else if (nstokes == 1) {
      norm = (*pmom)[0][0];
      for (int im = 0; im <= nmom; im++)
        (*pmom)[0][im] = (*pmom)[0][im] / norm;
    } else {
      fprintf (stderr, "Error, nstokes must be either 1 or 4\n");
      return -1;
    }
  }

  /* extinction coefficient */
  *beta = 3.0 * ext / (4.0 * vol) / 1000.0;
  /* Single scattering albedo */
  *omega = sca / ext;
  /* Asymmetry parameter */
  *g = asy / sca;
  /* Backscattering coefficient */
  *back = 3.0 * backsca / (4.0 * vol) / 1000.0;

  /* free memory */
  free (extinction);
  free (scattering);
  free (asymmetry);
  free (volume);

  if (nmom > 0) {
    for (int ip = 0; ip < nstokes; ip++) {
      for (int im = 0; im <= nmom; im++)
        free (pmom_tmp[ip][im]);
      free (pmom_tmp[ip]);
    }
    free (pmom_tmp);
  }

  return 0;
}

/** 
 * aerosol_input_data
 *
 * Reads aeosol input data from OPAC datafile 
 * data/aerosol/OPAC/size_distr.cfg. 
 * The data includes size distribution parameters 
 * and the density of the medium. 
 * 
 * @author Claudia Emde
 */
int aerosol_input_data (/* Input */
                        char* data_files_path,
                        int   aerosol_type,
                        int   verbose,
                        /* Output */
                        aerdist_struct* aerdistr) {

  int     i = 0, ih = 0;
  int     status;
  int     n_dist, n_hum;
  double* aerosol           = NULL;
  double* rel_hum           = NULL;
  double* r_min             = NULL;
  double* r_max             = NULL;
  double* rmod              = NULL;
  double* rho               = NULL;
  double* sigma             = NULL;
  char*   aerdistr_filename = NULL;

  aerdistr_filename = calloc (FILENAME_MAX, sizeof (char));

  strcpy (aerdistr_filename, data_files_path);
  strcat (aerdistr_filename, "aerosol/OPAC/size_distr.cfg");

  /* allocate memory for size distribution data */
  n_dist = 39; /* number of distributions in datafile */

  n_hum  = 0;
  status = read_7c_file (aerdistr_filename, &aerosol, &rel_hum, &r_min, &r_max, &rmod, &rho, &sigma, &n_dist);
  if (status != 0) {
    fprintf (stderr, "error %d reading %s\n", status, aerdistr_filename);
    return status;
  }

  /* Check whether aerosol type is water soluble or dry. In the first case size distributions
     for 8 relative humidities are needed; in the second case only 1. */
  if (aerosol_type == TYPE_INSO || aerosol_type == TYPE_MIAM || aerosol_type == TYPE_MICM || aerosol_type == TYPE_MINM ||
      aerosol_type == TYPE_MITR || aerosol_type == TYPE_SOOT) {
    n_hum = 1;
  } else {
    n_hum = 8;
  }

  /* allocate aerdistr_struct */
  (*aerdistr).n_hum   = n_hum;
  (*aerdistr).rel_hum = calloc (n_hum, sizeof (double));
  (*aerdistr).r_min   = calloc (n_hum, sizeof (double));
  (*aerdistr).r_max   = calloc (n_hum, sizeof (double));
  (*aerdistr).rmod    = calloc (n_hum, sizeof (double));
  (*aerdistr).rho     = calloc (n_hum, sizeof (double));

  (*aerdistr).type = aerosol_type;
  while (i < n_dist) {
    if ((int) aerosol[i] == aerosol_type) {
      (*aerdistr).rel_hum[ih] = rel_hum[i];
      (*aerdistr).r_min[ih]   = r_min[i];
      (*aerdistr).r_max[ih]   = r_max[i];
      (*aerdistr).rmod[ih]    = rmod[i];
      (*aerdistr).rho[ih]     = rho[i];
      if (ih == 0)
        (*aerdistr).sigma = sigma[i];
      ih++;
    }
    i++;
  }

  /* Check, whether all relative humidities have been  read. */
  if (ih != (*aerdistr).n_hum) {
    fprintf (stderr, "Something went wrong while reading aerosol input.ih %d nhum %d \n", ih, (*aerdistr).n_hum);
    return -1;
  }

  free (aerdistr_filename);
  free (aerosol);
  free (rel_hum);
  free (r_min);
  free (r_max);
  free (rmod);
  free (rho);
  free (sigma);
  return 0;
}

/** 
 * aerosol_size_distribution
 *
 * Aerosol size distributions are generated using the OPAC parameters. The 
 * distributions are log-normal distributions. The radius grid is generated 
 * using the method proposed by Evans (cloudprp.f), which seems to be an efficient 
 * and accurate method.
 *
 * @author Claudia Emde
 */
int aerosol_size_distribution (/* Input */
                               double rmin,
                               double rmax,
                               double rmod,
                               double sigma,
                               double lambda,
                               int    verbose,
                               /* Output */
                               int*     ndens,
                               double** radius,
                               double** number_dens) {

  int    i     = 1;
  double rstep = 0.0, a = 0.0, xstep = 0.0, x = 0.0;
  double rad = 0.0, nd = 0.0;

  /* Calculate number of radius grid points */
  rad = rmin;
  /* The next points are determined iteratively */
  while (rad <= rmax &&
         /* Cutoff for small particle number densities. */
         /* Assure that calculation is not stopped too early.*/
         (nd > 1e-12 || rad < rmod)) {

    /* Size parameter */
    x     = 2 * M_PI * rad / (lambda * 1.e-03); /* 1E-3 = nm -> micro m */
    xstep = double_max (0.03, 0.03 * sqrt (x));
    /* This is the same rule as proposed by Evans in cloudprp.f */
    rstep = xstep * lambda * 1.e-03 / (2.0 * M_PI);
    rad += rstep;
    a  = (log (rad) - log (rmod)) / log (sigma);
    nd = 1.0 / (sqrt (2.0 * M_PI) * log (sigma) * log (10.0) * rad) * exp (-0.5 * a * a);
    i++;
  }
  *ndens = i;

  /* Allocate memory */
  *radius      = calloc (*ndens, sizeof (double));
  *number_dens = calloc (*ndens, sizeof (double));

  /* Calculate size distribution, same as above, but now we may fill
     the arrays */
  (*radius)[0]      = rmin;
  a                 = (log ((*radius)[0]) - log (rmod)) / log (sigma);
  (*number_dens)[0] = 1.0 / (sqrt (2.0 * M_PI) * log (sigma) / log (10)) * exp (-0.5 * a * a);

  if (verbose) /*     fprintf(stderr, "aerosol size distribution: \n");  */
    fprintf (stderr, " ... %d radius grid points in aerosol size distribution \n", *ndens);

  i = 1;
  /* The next points are determined iteratively */
  while (i < *ndens) {
    /* Size parameter */
    x     = 2 * M_PI * (*radius)[i - 1] / (lambda * 1.e-03); /* 1E-3 = nm -> micro m */
    xstep = double_max (0.03, 0.03 * sqrt (x));
    /* This is the same rule as proposed by Evans in cloudprp.f */
    rstep             = xstep * lambda * 1.e-03 / (2.0 * M_PI);
    (*radius)[i]      = (*radius)[i - 1] + rstep;
    a                 = (log ((*radius)[i]) - log (rmod)) / log (sigma);
    (*number_dens)[i] = 1.0 / (sqrt (2.0 * M_PI) * log (sigma) * log (10.0) * (*radius)[i]) * exp (-0.5 * a * a);

    /* if(verbose) */
    /*       fprintf(stderr, "i %d %g   %g \n", i, (*radius)[i], (*number_dens)[i]); */

    i++;
  }

  /* An alternative would be to use a logarithmic radius. This is done in OPAC.f. */
  /* The solution above should be more accurate */
  /*  int n;  */
  /*   double rlog; */
  /*   n=1000;  */

  /*   rstep=(log(rmax)-log(rmin))/(double)n;  */

  /*   (*radius)[0]=rmin;   */
  /*   rlog=log(rmin);   */
  /*   i=1;    */
  /*   while (i<n+1){   */
  /*     rlog+=rstep;   */
  /*     /\* radius grid *\/  */
  /*     (*radius)[i]=exp(rlog);   */
  /*     /\* particle number densities, total number density is set to 1.*\/   */
  /*     a=(log((*radius[i])-log(rmod)))/  */
  /*       log(sigma);   */
  /*     (*number_dens)[i]=1/(sqrt(2*M_PI)*log(sigma)*log(10))*  */
  /*       exp(-0.5*a*a);    */
  /*     i++;   */
  /*   }   */
  /*   n_dens=n+1;   */

  return 0;
}

/** 
 * aerosol_refractive_index
 *
 * This function reads the refractive index from the OPAC data 
 * in data/aerosol/OPAC for the requested aerosol type and 
 * interpolates it on the required wavelength.
 * 
 * 
 * @author Claudia Emde
 */
int aerosol_refractive_index (/* Input */
                              double* lambda,
                              int     nlambda,
                              char*   data_files_path,
                              int     aerosol_type,
                              int     n_hum,
                              int     verbose,
                              /* Output */
                              double*** ref_real,
                              double*** ref_imag)

{
  int     status        = 0;
  double* lambda_data   = NULL;
  double* ref_real_data = NULL;
  double* ref_imag_data = NULL;
  int     nlambda_data = 0, ih = 0, iv = 0;
  char*   type;
  char*   rel_hum;
  char*   filename;

  /* Coefficents for interpolation */
  double *a0 = NULL, *a1 = NULL, *a2 = NULL, *a3 = NULL;

  if (strlen (data_files_path) == 0)
    data_files_path = "../data/";

  filename = calloc (FILENAME_MAX, sizeof (char));

  if (aerosol_type == TYPE_INSO)
    type = "inso";
  else if (aerosol_type == TYPE_WASO)
    type = "waso";
  else if (aerosol_type == TYPE_SOOT)
    type = "soot";
  else if (aerosol_type == TYPE_SSAM)
    type = "ssam";
  else if (aerosol_type == TYPE_SSCM)
    type = "sscm";
  else if (aerosol_type == TYPE_MINM)
    type = "minm";
  else if (aerosol_type == TYPE_MIAM)
    type = "miam";
  else if (aerosol_type == TYPE_MICM)
    type = "micm";
  else if (aerosol_type == TYPE_MITR)
    type = "mitr";
  else if (aerosol_type == TYPE_SUSO)
    type = "suso";
  else {
    fprintf (stderr, "Error: Unknown aerosol type \n");
    return -1;
  }

  /* Allocate refractive index array*/
  (*ref_real) = calloc (n_hum, sizeof (double*));
  (*ref_imag) = calloc (n_hum, sizeof (double*));
  for (ih = 0; ih < n_hum; ih++) {
    (*ref_real)[ih] = calloc (nlambda, sizeof (double));
    (*ref_imag)[ih] = calloc (nlambda, sizeof (double));
  }

  for (ih = 0; ih < n_hum; ih++) {

    if (ih == 0)
      rel_hum = "00";
    else if (ih == 1)
      rel_hum = "50";
    else if (ih == 2)
      rel_hum = "70";
    else if (ih == 3)
      rel_hum = "80";
    else if (ih == 4)
      rel_hum = "90";
    else if (ih == 5)
      rel_hum = "95";
    else if (ih == 6)
      rel_hum = "98";
    else if (ih == 7)
      rel_hum = "99";
    else {
      fprintf (stderr, "Error: Humidity index ih=%d is wrong. \n", ih);
      return -1;
    }

    /* Create datafile name */
    strcpy (filename, data_files_path);
    strcat (filename, "aerosol/OPAC/refractive_indices/");
    strcat (filename, type);
    strcat (filename, rel_hum);
    strcat (filename, "_refr.dat");

    status = read_3c_file (filename, &lambda_data, &(ref_real_data), &(ref_imag_data), &nlambda_data);
    if (status != 0) {
      fprintf (stderr, "error %d reading %s\n", status, filename);
      return status;
    }

    /* for (iv=0; iv<nlambda_data; iv++) */
    /*       fprintf(stderr, "%g %g %g \n", lambda_data[iv], ref_real_data[iv], ref_imag_data[iv]); */

    /* Interpolate refractive index on the required wavelength. */
    /* real part*/
    status = linear_coeffc (lambda_data, ref_real_data, nlambda_data, &a0, &a1, &a2, &a3);
    for (iv = 0; iv < nlambda; iv++) {
      status = calc_splined_value (lambda[iv] * 1e-3, &(*ref_real)[ih][iv], lambda_data, nlambda_data, a0, a1, a2, a3);
      if (status != 0) {
        fprintf (stderr, "Sorry cannot interpolate refractive index (real)\n");
        return status;
      }
    }
    free (a0);
    free (a1);
    free (a2);
    free (a3);

    /* imaginary part */
    status = linear_coeffc (lambda_data, ref_imag_data, nlambda_data, &a0, &a1, &a2, &a3);
    for (iv = 0; iv < nlambda; iv++) {
      status = calc_splined_value (lambda[iv] * 1e-3, &(*ref_imag)[ih][iv], lambda_data, nlambda_data, a0, a1, a2, a3);
      if (status != 0) {
        fprintf (stderr, "Sorry cannot interpolate refractive index (imaginary)\n");
        return status;
      }
      (*ref_imag)[ih][iv] *= -1; /* Mie tool requires positive values*/
    }
    free (a0);
    free (a1);
    free (a2);
    free (a3);
  }

  /* Free variables */

  free (filename);
  free (lambda_data);

  free (ref_real_data);
  free (ref_imag_data);

  return 0;
}

/**
 * write_output
 * 
 * Write output to the screen or into output files. The output file
 * format corresponds to the output of the cloudprp program by
 * F. Evans. The screen output can be specified by the user. 
 * 
 * @author ???, C. Emde (phase matrix output) 
 * @date ???
 */

int write_output (input_struct input, output_struct* output, int ir, int iv) {

  int  nmom = -999, nstokes = 1;
  char outputfile[255];
  char ck_band_str[14];
  int  status         = 0;
  int  nphamat        = 0;
  int  default_output = 0;

  if (input.output_format == OUTPUT_FORMAT_NORMAL) {

    if (input.n_output_user == 0) {
      /* default output */
      input.n_output_user  = 8;
      input.output_user    = calloc (input.n_output_user, sizeof (int));
      input.output_user[0] = OUTPUT_USER_WAVE;
      input.output_user[1] = OUTPUT_USER_REFRAC_REAL;
      input.output_user[2] = OUTPUT_USER_REFRAC_IMAG;
      input.output_user[3] = OUTPUT_USER_QEXT;
      input.output_user[4] = OUTPUT_USER_OMEGA;
      input.output_user[5] = OUTPUT_USER_GG;
      input.output_user[6] = OUTPUT_USER_SPIKE;
      input.output_user[7] = OUTPUT_USER_PMOM;
      default_output       = 1;
    }

    for (int io = 0; io < input.n_output_user; io++) {
      switch (input.output_user[io]) {
      case OUTPUT_USER_NONE:
        fprintf (stdout, "%9s ", " ");
        break;
      case OUTPUT_USER_WAVE:
        fprintf (stdout, "%9.3f ", output->wl.lambda[iv]);
        break;
      case OUTPUT_USER_WAVENUMBER:
        fprintf (stdout, "%9.3f ", 1E7 / output->wl.lambda[iv]);
        break;
      case OUTPUT_USER_R_EFF:
        fprintf (stdout, "%12.6f ", output->r_eff[ir]);
        break;
      case OUTPUT_USER_REFRAC_REAL:
        fprintf (stdout, "%12.6e ", output->ref.re);
        break;
      case OUTPUT_USER_REFRAC_IMAG:
        fprintf (stdout, "%12.6e ", output->ref.im);
        break;
      case OUTPUT_USER_QEXT:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e ", output->mie.qext);
        else
          fprintf (stdout, "%12.6e ", output->beta);
        break;
      case OUTPUT_USER_QSCA:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e ", output->mie.qsca);
        else
          fprintf (stdout, "%12.6e ", output->omega * output->beta);
        break;
      case OUTPUT_USER_OMEGA:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e ", output->mie.qsca / output->mie.qext);
        else
          fprintf (stdout, "%12.6e ", output->omega);
        break;
      case OUTPUT_USER_GG:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e ", output->mie.gsca);
        else
          fprintf (stdout, "%12.6e ", output->g);
        break;
      case OUTPUT_USER_SPIKE:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e ", output->mie.spike);
        else
          fprintf (stdout, "%12.6e ", -1.0);
        break;
      case OUTPUT_USER_SFORW:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e %12.6e ", output->mie.sforw.re, output->mie.sforw.im);
        else
          fprintf (stdout, "%12.6e %12.6e ", NAN, NAN); /* not yet implemented*/
        break;
      case OUTPUT_USER_SBACK:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e %12.6e ", output->mie.sback.re, output->mie.sback.im);
        else
          fprintf (stdout, "%12.6e %12.6e ", NAN, NAN); /* not yet implemented*/
        break;
      case OUTPUT_USER_QBACK:
        if (input.program == BHMIE) {
          if (input.distribution == DIST_NONE)
            fprintf (stdout, "%12.6e ", output->mie.qback);
          else
            fprintf (stdout, "%12.6e ", output->back);
        } else
          fprintf (stdout, "%12.6e ", NAN); /* not yet implemented */
        break;
      case OUTPUT_USER_PMOM:
        if (input.mie.nmom > 0) {
          if (input.nstokes > 1 || default_output == 0) {
            fprintf (stdout, "\n");
            for (int im = 0; im <= input.mie.nmom; im++) {
              for (int ip = 0; ip < input.nstokes; ip++)
                fprintf (stdout, " %14.6e", output->pmom[ip][im]);
              fprintf (stdout, "\n");
            }
          } else {
            for (int im = 0; im <= input.mie.nmom; im++)
              fprintf (stdout, " %14.6e", output->pmom[0][im]);
            fprintf (stdout, "\n");
          }
        }
        break;
      default:
        fprintf (stderr, "Error, unknown user output %d\n", input.output_user[io]);
        break;
      } /* switch (input.output_user[io]) { */
    }   /* for (io=0; io<input.n_output_user; io++) { */
    fprintf (stdout, "\n");
    fflush (stdout);
  } /* if (input.output_user[0] != OUTPUT_USER_CLOUDPROP) { */

  else if (input.output_format == OUTPUT_FORMAT_CLOUDPROP || input.output_format == OUTPUT_FORMAT_AERPROP) {

    if (input.distribution == DIST_NONE) {
      fprintf (stderr, "The option output_user cloudprop works only for size distributions, not for mono-disperse particles. \n");
      return -1;
    }

    /* For polarized calculations, all four elements are calculated for mie files */
    if (input.nstokes > 1)
      nstokes = 4;

    if (ir == 0) {
      /* generate file name */
      switch (input.medium) {
      case WATER:
        strcpy (outputfile, "wc.");
        break;
      case ICE:
        strcpy (outputfile, "ic.");
        break;
      case USER:
        strcpy (outputfile, "uc.");
        break;
      case AEROSOL:
        if (input.aerosol_type == TYPE_INSO)
          strcpy (outputfile, "inso.");
        else if (input.aerosol_type == TYPE_WASO)
          strcpy (outputfile, "waso.");
        else if (input.aerosol_type == TYPE_SOOT)
          strcpy (outputfile, "soot.");
        else if (input.aerosol_type == TYPE_SSAM)
          strcpy (outputfile, "ssam.");
        else if (input.aerosol_type == TYPE_SSCM)
          strcpy (outputfile, "sscm.");
        else if (input.aerosol_type == TYPE_MINM)
          strcpy (outputfile, "minm.");
        else if (input.aerosol_type == TYPE_MIAM)
          strcpy (outputfile, "miam.");
        else if (input.aerosol_type == TYPE_MICM)
          strcpy (outputfile, "micm.");
        else if (input.aerosol_type == TYPE_MITR)
          strcpy (outputfile, "mitr.");
        else if (input.aerosol_type == TYPE_SUSO)
          strcpy (outputfile, "suso.");
        else {
          fprintf (stderr, "Error: Unknown aerosol type \n");
          return -1;
        }
        break;
      }

      if (iv < -1 || iv > 998)
        return err_out ("band has to be positive and smaller than 999, otherwise doesnt fit the string here", iv);
      snprintf (ck_band_str, strlen (ck_band_str), "%03d", iv + 1);
      strcat (outputfile, ck_band_str);
      strcat (outputfile, ".mie");
      fprintf (stderr, " ... generating file: %s\n", outputfile);

      if ((output->ptr_output = fopen (outputfile, "w")) == NULL) {
        fprintf (stderr, "Error, cannot open file: %s in write_output! \n", outputfile);
        return -1;
      }
    }

    if (input.medium != AEROSOL) {
      fprintf (stderr,
               "   lambda = %10.3f nm (%3d of %3d), r_eff = %8.3f micro meter (%3d of %3d) \n",
               output->wl.lambda[iv],
               iv + 1,
               output->wl.nlambda + 1,
               output->r_eff[ir],
               ir + 1,
               output->n_r_eff + 1);
    } else {
      fprintf (stderr,
               "   lambda = %10.3f nm (%3d of %3d), rel_hum = %4.0f percent (%3d of %3d) \n",
               output->wl.lambda[iv],
               iv + 1,
               output->wl.nlambda + 1,
               output->aerdistr.rel_hum[ir],
               ir + 1,
               output->aerdistr.n_hum);
    }

    if (ir == 0) {
      switch (input.medium) {
      case WATER:
        fprintf (output->ptr_output, "! Mie table vs. effective radius (LWC=1 g/m^3)\n");
        break;
      case ICE:
        fprintf (output->ptr_output, "! Mie table vs. effective radius (IWC=1 g/m^3)\n");
        break;
      case USER:
        fprintf (output->ptr_output,
                 "! Mie table vs. effective radius (ATTENTION!!! assume rho_medium = %7.4f g/cm3)\n",
                 output->rho_medium);
        break;
      case AEROSOL:
        fprintf (output->ptr_output, "! Mie table vs. relative humidity (N=1, rho_medium = %7.4f g/cm3)\n", output->rho_medium);
        break;
      }

      fprintf (output->ptr_output, " %12.6e  wavelength (micron) \n", output->wl.lambda[iv] / 1000.);
      fprintf (output->ptr_output, " %12.6e  %12.6e  index of refraction \n", output->ref.re, output->ref.im);
      if (input.distribution == DIST_NONE) {
        fprintf (output->ptr_output, "S  distribution type ! one Single radius  \n");
        fprintf (output->ptr_output, "-999.99  distribution shape parameter !(dummy to fit the format)\n");
      } else if (input.distribution == DIST_GAMMA) {
        fprintf (output->ptr_output, "G  distribution type\n");
        fprintf (output->ptr_output, "%7.4f  distribution shape parameter\n", input.alpha);
      } else if (input.distribution == DIST_LOGNORMAL) {
        fprintf (output->ptr_output, "L  distribution type\n");
        fprintf (output->ptr_output, "%7.4f  distribution shape parameter\n", input.alpha);
      } else if (input.distribution == DIST_FILE) {
        fprintf (output->ptr_output, "U  distribution type ! User defined, from file \n");
        fprintf (output->ptr_output, "-999.99  distribution shape parameter !(dummy to fit the format)\n");
      } else if (input.distribution == DIST_AER) {
        fprintf (output->ptr_output, "A  distribution type \n");
        fprintf (output->ptr_output, "Size distribution parameters from OPAC database \n");
      }

      if (input.medium != AEROSOL) {
        fprintf (output->ptr_output,
                 "%4d %8.4f %8.4f  number, starting, ending effective radius \n",
                 output->n_r_eff + 1,
                 output->r_eff[0],
                 output->r_eff[output->n_r_eff]);
      } else {
        fprintf (output->ptr_output,
                 "%4d %8.4f %8.4f  number, starting, ending relative humidity \n",
                 output->aerdistr.n_hum,
                 output->aerdistr.rel_hum[0],
                 output->aerdistr.rel_hum[output->aerdistr.n_hum - 1]);
      }

    } /* if (ir == 0) { */

    nmom = 0;
    for (int ip = 0; ip < nstokes; ip++) {
      for (int im = 0; im <= input.mie.nmom; im++) {
        if (output->pmom[ip][im] == FLT_MAX) {
          if (im > nmom)
            nmom = im;
          break;
        }
      }
    }

    /* remember to write a warning, if there could be more Legendre coefficients */
    if (nmom == input.mie.nmom)
      output->warnings_nmom[iv][ir] = 1;

    if (nstokes == 1)
      nphamat = 1;
    else
      nphamat = 4;

    if (input.medium != AEROSOL) {
      fprintf (output->ptr_output,
               "%9.4f  %12.5e  %12.5e %6d  %4d Reff  Ext  Alb  Nleg Nphamat\n",
               output->r_eff[ir],
               output->beta / output->rho_medium,
               output->omega,
               nmom,
               nphamat);
    } else {
      fprintf (output->ptr_output,
               "%9.4f  %12.5e  %12.5e %6d  %4d RelHum  Ext  Alb  Nleg Nphamat\n",
               output->aerdistr.rel_hum[ir],
               output->beta / output->rho_medium,
               output->omega,
               nmom,
               nphamat);
    }

    if (input.mie.nmom > 0) {
      fprintf (output->ptr_output, "  ");
      for (int ip = 0; ip < input.nstokes; ip++) {
        for (int im = 0; im <= nmom; im++)
          fprintf (output->ptr_output, " %9.5f", (2 * im + 1) * output->pmom[ip][im]);
        fprintf (output->ptr_output, "\n");
      }
    }
    fprintf (output->ptr_output, "\n");
    fflush (output->ptr_output);

    if ((ir == output->n_r_eff && input.output_format == OUTPUT_FORMAT_CLOUDPROP) ||
        (ir == output->aerdistr.n_hum && input.output_format == OUTPUT_FORMAT_AERPROP))
      fclose (output->ptr_output);

  } /* else if (input.output_format == OUTPUT_FORMAT_CLOUDPROP || OUTPUT_FORMAT_AERPROP) ... */

  return status;
}

#if HAVE_NETCDF4

/** 
 * create_netcdf_file 
 *
 * Create netcdf file.
 * This can directly be used as input for uvspec. 
 * 
 *
 * @author  Claudia Emde
 * @date    2010-01-09 Created
 */
int create_netcdf_file (/* Input */
                        char*   basename,
                        int     medium,
                        int     aerosol_type,
                        int     distribution,
                        double  alpha,
                        double  step,
                        int     nlam,
                        double* lambda,
                        int     nreff,
                        double* reff,
                        int     nmommax,
                        int     nstokes,
                        int     nthetamax,
                        int     nrho,
                        int     complete_moments,
                        double  n_r_max,
                        double  dx_max,
                        char*   sd_filename,
                        /* Output */
                        netcdf_ids* ncid,
                        char**      outputfile) {
  /* Indices*/
  int retval = 0, status = 0, i = 0;

  /* Netcdf dimension IDs */
  int dimids2D[2], dimids3D[3], dimids4D[4];

  float fill_value[1];
  char  name[MAXLINE], unit[23], varname[10];

  char   medium_str[MAXLINE];
  char   alpha_str[MAXLINE];
  char   distr[MAXLINE], param[MAXLINE];
  double distr_param[1] = {0};
  double distr_step[1]  = {0};
  long   version[1];
  int    comp_moments[1];

  int    nphamat = nstokes;
  double lambda_tmp[nlam];

  /*-------------------------------------------------------------------*/
  /* Generate filename                                                 */
  /*-------------------------------------------------------------------*/

  if (strlen (basename) != 0) {
    strcpy (*outputfile, basename);
  } else {
    switch (medium) {
    case WATER:
      strcpy (*outputfile, "wc");
      strcpy (medium_str, "water");
      break;
    case ICE:
      strcpy (*outputfile, "ic");
      strcpy (medium_str, "ice");
      break;
    case USER:
      strcpy (*outputfile, "uc");
      strcpy (medium_str, "user defined refractive index");
      break;
    case AEROSOL:

      if (aerosol_type == TYPE_INSO) {
        strcpy (*outputfile, "inso");
        strcpy (medium_str, "Aerosol: insoluble");
      } else if (aerosol_type == TYPE_WASO) {
        strcpy (*outputfile, "waso");
        strcpy (medium_str, "Aerosol: water soluble");
      } else if (aerosol_type == TYPE_SOOT) {
        strcpy (*outputfile, "soot");
        strcpy (medium_str, "Aerosol: soot");
      } else if (aerosol_type == TYPE_SSAM) {
        strcpy (*outputfile, "ssam");
        strcpy (medium_str, "Aerosol: sea salt accumulated mode");
      } else if (aerosol_type == TYPE_SSCM) {
        strcpy (*outputfile, "sscm");
        strcpy (medium_str, "Aerosol: sea salt coarse mode");
      } else if (aerosol_type == TYPE_MINM) {
        strcpy (*outputfile, "minm");
        strcpy (medium_str, "Aerosol: mineral nucleation mode");
      } else if (aerosol_type == TYPE_MIAM) {
        strcpy (*outputfile, "miam");
        strcpy (medium_str, "Aerosol: mineral accumulated mode");
      } else if (aerosol_type == TYPE_MICM) {
        strcpy (*outputfile, "micm");
        strcpy (medium_str, "Aerosol: mineral coarse mode");
      } else if (aerosol_type == TYPE_MITR) {
        strcpy (*outputfile, "mitr");
        strcpy (medium_str, "Aerosol: mineral transported");
      } else if (aerosol_type == TYPE_SUSO) {
        strcpy (*outputfile, "suso");
        strcpy (medium_str, "Aerosol: sulfate");
      } else {
        fprintf (stderr, "Error: Unknown aerosol type \n");
        return -1;
      }
      break;
    default:
      fprintf (stderr, "Error in write_output_netcdf : unrecognized medium %d\n", medium);
      return -1;
    }
  }

  switch (distribution) {
  case DIST_NONE:
    strcpy (distr, "Calculation for single radius.");
    distr_param[0] = -999.;
    distr_step[0]  = -999.;
    break;
  case DIST_GAMMA:
    strcpy (distr, "Gamma distribution.");
    distr_param[0] = alpha;
    snprintf (alpha_str, sizeof (alpha_str), ".gamma_%05.1f", alpha);
    strcat (*outputfile, alpha_str);
    break;
  case DIST_LOGNORMAL:
    strcpy (distr, "Lognormal distribution.");
    distr_param[0] = alpha;
    snprintf (alpha_str, sizeof (alpha_str), ".lognormal_%05.1f", alpha);
    strcat (*outputfile, alpha_str);
    break;
  case DIST_TRIANGULAR:
    strcpy (distr, "Triangular distribution.");
    distr_param[0] = alpha;
    distr_step[0]  = step;
    snprintf (alpha_str, sizeof (alpha_str), ".triangular_%05.1f", alpha);
    strcat (*outputfile, alpha_str);
    break;
  case DIST_FILE:
    strcpy (distr, "User defined size distribution.");
    distr_param[0] = -999.;
    break;
  case DIST_AER:
    strcpy (distr, "Aerosol size distribution,\nlog-normal with parameters from OPAC");
    distr_param[0] = -999.;
    break;
  default:
    fprintf (stderr, "Error in write_output_netcdf : unrecognized distribution %d\n", distribution);
    return -1;
    break;
  }

  strcat (*outputfile, ".mie.cdf");
  fprintf (stderr, " ... generating file: %s\n", *outputfile);

  if ((retval = nc_create (*outputfile, NC_CLOBBER | NC_NETCDF4, &(*ncid).ncid))) {
    fprintf (stderr, "Error: creating netcdf file");
    return -1;
  }

  /* ----------------------------------------------------------------*/
  /* Set global attributes                                           */
  /* ----------------------------------------------------------------*/

  version[0] = 20110916;
  strcpy (param, "mie");
  comp_moments[0] = complete_moments;

  strcpy (name, "Netcdf file created using libRadtran-Mie-tool");
  nc_put_att_text ((*ncid).ncid, NC_GLOBAL, "file_info", strlen (name), name);

  nc_put_att_long ((*ncid).ncid, NC_GLOBAL, "version", NC_LONG, 1, version);

  nc_put_att_int ((*ncid).ncid, NC_GLOBAL, "complete_moments", NC_INT, 1, comp_moments);

  nc_put_att_text ((*ncid).ncid, NC_GLOBAL, "parameterization", strlen (param), param);

  nc_put_att_text ((*ncid).ncid, NC_GLOBAL, "size_distr", strlen (distr), distr);

  nc_put_att_double ((*ncid).ncid, NC_GLOBAL, "param_alpha", NC_DOUBLE, 1, distr_param);

  nc_put_att_double ((*ncid).ncid, NC_GLOBAL, "param_step", NC_DOUBLE, 1, distr_step);

  nc_put_att_double ((*ncid).ncid, NC_GLOBAL, "n_r_max", NC_DOUBLE, 1, &n_r_max);

  nc_put_att_double ((*ncid).ncid, NC_GLOBAL, "dx_max", NC_DOUBLE, 1, &dx_max);

  if (strcmp (sd_filename, "") == 1)
    nc_put_att_text ((*ncid).ncid, NC_GLOBAL, "size distribution file", strlen (sd_filename), sd_filename);

  /* ----------------------------------------------------------------*/
  /* Define dimensions.                                              */
  /* ----------------------------------------------------------------*/

  if ((retval = nc_def_dim ((*ncid).ncid, "nlam", nlam, &(*ncid).id_nlam))) {
    fprintf (stderr, "Error creating dimension nlam.\n");
    return -1;
  }

  if ((retval = nc_def_dim ((*ncid).ncid, "nmommax", nmommax, &(*ncid).id_nmommax))) {
    fprintf (stderr, "Error creating dimension nmommax.\n");
    return -1;
  }

  if ((retval = nc_def_dim ((*ncid).ncid, "nphamat", nphamat, &(*ncid).id_nphamat))) {
    fprintf (stderr, "Error creating dimension nphase.\n");
  }

  if (medium == AEROSOL) {
    if ((retval = nc_def_dim ((*ncid).ncid, "nhum", nreff, &(*ncid).id_nreff))) {
      fprintf (stderr, "Error creating dimension nreff.\n");
      return -1;
    }
  } else {
    if ((retval = nc_def_dim ((*ncid).ncid, "nreff", nreff, &(*ncid).id_nreff))) {
      fprintf (stderr, "Error creating dimension nreff.\n");
      return -1;
    }
  }

  if (nthetamax > 0) {
    if ((retval = nc_def_dim ((*ncid).ncid, "nthetamax", nthetamax, &(*ncid).id_nthetamax))) {
      fprintf (stderr, "Error creating dimension nthetamax.\n");
      return -1;
    }
  }

  if (medium != AEROSOL)
    nrho = 1;

  if ((retval = nc_def_dim ((*ncid).ncid, "nrho", nrho, &(*ncid).id_nrho))) {
    fprintf (stderr, "Error creating dimension nrho.\n");
    return -1;
  }

  /* ---------------------------------------------------------------*/
  /* Define variables                                               */
  /* ---------------------------------------------------------------*/

  /* Wavelength*/

  if ((retval = nc_def_var ((*ncid).ncid, "wavelen", NC_DOUBLE, 1, &(*ncid).id_nlam, &(*ncid).id_lam))) {
    fprintf (stderr, "Error: defining variable wavelen. \n");
    return -1;
  }
  strcpy (name, "wavelength");
  strcpy (unit, "micrometer");
  nc_put_att_text ((*ncid).ncid, (*ncid).id_lam, "long_name", strlen (name), name);
  nc_put_att_text ((*ncid).ncid, (*ncid).id_lam, "units", strlen (unit), unit);

  /* effective radius or relative humidity */

  if (medium == AEROSOL) {
    strcpy (varname, "hum");
    strcpy (name, "relative humidity");
    strcpy (unit, "per cent");
  } else {
    strcpy (varname, "reff");
    strcpy (name, "effective radius");
    strcpy (unit, "micrometer");
  }
  if ((retval = nc_def_var ((*ncid).ncid, varname, NC_DOUBLE, 1, &(*ncid).id_nreff, &(*ncid).id_reff))) {
    fprintf (stderr, "Error defining variable reff. \n");
    return -1;
  }
  nc_put_att_text ((*ncid).ncid, (*ncid).id_reff, "long_name", strlen (name), name);
  nc_put_att_text ((*ncid).ncid, (*ncid).id_reff, "units", strlen (unit), unit);

  /* Number of scattering angles */
  if (nthetamax > 0) {

    dimids3D[0] = (*ncid).id_nlam;
    dimids3D[1] = (*ncid).id_nreff;
    dimids3D[2] = (*ncid).id_nphamat;

    if ((retval = nc_def_var ((*ncid).ncid, "ntheta", NC_INT, 3, dimids3D, &(*ncid).id_ntheta))) {
      fprintf (stderr, "Error defining variable ntheta. \n");
      return -1;
    }
    strcpy (name, "number of scattering angles");
    strcpy (unit, "-");
    nc_put_att_text ((*ncid).ncid, (*ncid).id_ntheta, "long_name", strlen (name), name);
    nc_put_att_text ((*ncid).ncid, (*ncid).id_ntheta, "units", strlen (unit), unit);

    /* Scattering angle grid */

    fill_value[0] = -999.0;

    /* The dimids array is used to pass the IDs of the dimensions of
     * the variable. */
    dimids4D[0] = (*ncid).id_nlam;
    dimids4D[1] = (*ncid).id_nreff;
    dimids4D[2] = (*ncid).id_nphamat;
    dimids4D[3] = (*ncid).id_nthetamax;

    fprintf (stderr, "nthetamax = %d %d %d\n", (*ncid).id_nthetamax, dimids4D[3], nthetamax);

    if ((retval = nc_def_var ((*ncid).ncid, "theta", NC_FLOAT, 4, dimids4D, &(*ncid).id_theta))) {
      fprintf (stderr, "Error defining variable theta. \n");
      return -1;
    }
    strcpy (name, "theta");
    strcpy (unit, "degrees");
    nc_put_att_text ((*ncid).ncid, (*ncid).id_theta, "long_name", strlen (name), name);
    nc_put_att_text ((*ncid).ncid, (*ncid).id_theta, "units", strlen (unit), unit);
    nc_put_att_float ((*ncid).ncid, (*ncid).id_theta, "_FillValue", NC_FLOAT, 1, fill_value);

    /* Phase matrix */

    if ((retval = nc_def_var ((*ncid).ncid, "phase", NC_FLOAT, 4, dimids4D, &(*ncid).id_phase))) {
      fprintf (stderr, "Error defining variable phase. \n");
      return -1;
    }
    strcpy (name, "phase");
    strcpy (unit, "-");
    nc_put_att_text ((*ncid).ncid, (*ncid).id_phase, "long_name", strlen (name), name);
    nc_put_att_text ((*ncid).ncid, (*ncid).id_phase, "units", strlen (unit), unit);
    nc_put_att_float ((*ncid).ncid, (*ncid).id_phase, "_FillValue", NC_FLOAT, 1, fill_value);
  }

  /* Number of Legendre polynomials */
  dimids3D[0] = (*ncid).id_nlam;
  dimids3D[1] = (*ncid).id_nreff;
  dimids3D[2] = (*ncid).id_nphamat;

  if ((retval = nc_def_var ((*ncid).ncid, "nmom", NC_INT, 3, dimids3D, &(*ncid).id_nmom))) {
    fprintf (stderr, "Error defining variable nmom. \n");
    return -1;
  }
  strcpy (name, "number of Legendre polynomials");
  strcpy (unit, "-");
  nc_put_att_text ((*ncid).ncid, (*ncid).id_nmom, "long_name", strlen (name), name);
  nc_put_att_text ((*ncid).ncid, (*ncid).id_nmom, "units", strlen (unit), unit);

  /* Legendre moments */
  fill_value[0] = FLT_MAX;

  dimids4D[0] = (*ncid).id_nlam;
  dimids4D[1] = (*ncid).id_nreff;
  dimids4D[2] = (*ncid).id_nphamat;
  dimids4D[3] = (*ncid).id_nmommax;

  if ((retval = nc_def_var ((*ncid).ncid, "pmom", NC_FLOAT, 4, dimids4D, &(*ncid).id_pmom))) {
    fprintf (stderr, "Error defining variable pmom. \n");
    return -1;
  }
  strcpy (name, "Legendre polynomials");
  strcpy (unit, "including factor 2*l+1");
  nc_put_att_text ((*ncid).ncid, (*ncid).id_pmom, "long_name", strlen (name), name);
  nc_put_att_text ((*ncid).ncid, (*ncid).id_pmom, "units", strlen (unit), unit);
  nc_put_att_float ((*ncid).ncid, (*ncid).id_pmom, "_FillValue", NC_FLOAT, 1, fill_value);

  /* Extinction coefficient */

  dimids2D[0] = (*ncid).id_nlam;
  dimids2D[1] = (*ncid).id_nreff;

  if ((retval = nc_def_var ((*ncid).ncid, "ext", NC_DOUBLE, 2, dimids2D, &(*ncid).id_ext))) {
    fprintf (stderr, "Error defining variable ext. \n");
    return -1;
  }
  strcpy (name, "extinction coefficient");
  strcpy (unit, "km^-1/(g/m^3)");
  nc_put_att_text ((*ncid).ncid, (*ncid).id_ext, "long_name", strlen (name), name);
  nc_put_att_text ((*ncid).ncid, (*ncid).id_ext, "units", strlen (unit), unit);

  /* Single scattering albedo */

  if ((retval = nc_def_var ((*ncid).ncid, "ssa", NC_DOUBLE, 2, dimids2D, &(*ncid).id_ssa))) {
    fprintf (stderr, "Error defining variable ssa. \n");
    return -1;
  }
  strcpy (name, "single scattering albedo");
  strcpy (unit, "-");
  nc_put_att_text ((*ncid).ncid, (*ncid).id_ssa, "long_name", strlen (name), name);
  nc_put_att_text ((*ncid).ncid, (*ncid).id_ssa, "units", strlen (unit), unit);

  /* Asymmetry factor */

  if ((retval = nc_def_var ((*ncid).ncid, "gg", NC_DOUBLE, 2, dimids2D, &(*ncid).id_gg))) {
    fprintf (stderr, "Error defining variable gg. \n");
    return -1;
  }
  strcpy (name, "Asymmetry factor");
  strcpy (unit, "-");
  nc_put_att_text ((*ncid).ncid, (*ncid).id_gg, "long_name", strlen (name), name);
  nc_put_att_text ((*ncid).ncid, (*ncid).id_gg, "units", strlen (unit), unit);

  /* Refractive index (real) */

  if (medium != AEROSOL) {
    if ((retval = nc_def_var ((*ncid).ncid, "refre", NC_DOUBLE, 1, &(*ncid).id_nlam, &(*ncid).id_refre))) {
      fprintf (stderr, "Error defining variable refre. \n");
      return -1;
    }
  } else {
    if ((retval = nc_def_var ((*ncid).ncid, "refre", NC_DOUBLE, 2, dimids2D, &(*ncid).id_refre))) {
      fprintf (stderr, "Error defining variable refre. \n");
      return -1;
    }
  }

  strcpy (name, "refractive index (real)");
  strcpy (unit, " ");
  nc_put_att_text ((*ncid).ncid, (*ncid).id_refre, "long_name", strlen (name), name);
  nc_put_att_text ((*ncid).ncid, (*ncid).id_refre, "units", strlen (unit), unit);

  /* Refractive index (imaginary) */
  if (medium != AEROSOL) {
    if ((retval = nc_def_var ((*ncid).ncid, "refim", NC_DOUBLE, 1, &(*ncid).id_nlam, &(*ncid).id_refim))) {
      fprintf (stderr, "Error defining variable refre. \n");
      return -1;
    }
  } else {
    if ((retval = nc_def_var ((*ncid).ncid, "refim", NC_DOUBLE, 2, dimids2D, &(*ncid).id_refim))) {
      fprintf (stderr, "Error defining variable refre. \n");
      return -1;
    }
  }

  strcpy (name, "refractive index (imaginary)");
  strcpy (unit, " ");
  nc_put_att_text ((*ncid).ncid, (*ncid).id_refim, "long_name", strlen (name), name);
  nc_put_att_text ((*ncid).ncid, (*ncid).id_refim, "units", strlen (unit), unit);

  if ((retval = nc_def_var ((*ncid).ncid, "rho", NC_DOUBLE, 2, dimids2D, &(*ncid).id_rho))) {
    fprintf (stderr, "Error defining variable refre. \n");
    return -1;
  }
  strcpy (name, "density of medium");
  strcpy (unit, "g/cm^3");
  nc_put_att_text ((*ncid).ncid, (*ncid).id_rho, "long_name", strlen (name), name);
  nc_put_att_text ((*ncid).ncid, (*ncid).id_rho, "units", strlen (unit), unit);

  /* End definitions */
  if ((retval = nc_enddef ((*ncid).ncid))) {
    fprintf (stderr, "Error netcdf definition end:\n");
    fprintf (stderr, "%s\n", nc_strerror (retval));
    fprintf (stderr, "nlam * nreff * nphamat * nthetamax = %d \n", nlam * nreff * nphamat * nthetamax);
    return -1;
  }

  /* Variables that are not depedent on ir, iv  */
  /* wavelengths */
  for (i = 0; i < nlam; i++)
    lambda_tmp[i] = lambda[i] * 1e-3;
  if ((retval = nc_put_var_double ((*ncid).ncid, (*ncid).id_lam, &lambda_tmp[0]))) {
    fprintf (stderr, "Error writing variable lambda.\n");
    return -1;
  }

  /* effective radius */
  if ((retval = nc_put_var_double ((*ncid).ncid, (*ncid).id_reff, &reff[0]))) {
    fprintf (stderr, "Error writing variable reff.\n");
    return -1;
  }

  /* Close netcdf file. This frees up any internal netCDF resources
   * associated with the file, and flushes any buffers. */
  if ((status = nc_close ((*ncid).ncid))) {
    fprintf (stderr, "Error closing netcdf file.\n");
    return -1;
  }

  return 0;
}

/** 
 * write_output_netcdf 
 *
 * Write the output of the Mie calculation in netcdf format. 
 * This can directly be used as input for uvspec. 
 * 
 *
 * @author Claudia Emde
 * @date   ????
 *         2009-06-29 (modifications to adapt to new format of uvspec files)
 *         2010-02-10 (rewrite)
 */
int write_output_netcdf (netcdf_ids ncid, input_struct input, output_struct* output, int iv, int ir, int ia) {
  /* Indices*/
  int i      = 0;
  int retval = 0, status = 0;

  int nphamat = input.nstokes;

  /* Temporary 1D arrays for netcdf output */
  int    int_array_tmp[1];
  int    id;
  double double_array_tmp[1];
  float* pmom_tmp;
  int    nmoms = 0, nmommax = 0;
  int    comp_moments[1];

  size_t    start1D[1], start2D[2], start3D[3], start4D[4];
  ptrdiff_t stride1D[1], stride2D[2], stride3D[3], stride4D[4];
  size_t    count1D[1], count2D[2], count3D[3], count4D[4];

  /* Initialization*/
  start1D[0]  = iv;
  count1D[0]  = 1;
  stride1D[0] = 1;

  start2D[0] = iv;
  start2D[1] = ir;
  for (i = 0; i < 2; i++) {
    count2D[i]  = 1; /* number of indices to be written*/
    stride2D[i] = 1; /* index interval */
  }

  start3D[0] = iv;
  start3D[1] = ir;
  start3D[2] = 0;
  for (i = 0; i < 3; i++) {
    count3D[i]  = 1;
    stride3D[i] = 1;
  }

  start4D[0] = iv;
  start4D[1] = ir;
  start4D[2] = 0;
  start4D[3] = 0;
  for (i = 0; i < 4; i++) {
    count4D[i]  = 1;
    stride4D[i] = 1;
  }
  if (input.nmom_netcdf != 0) {
    nmommax = input.nmom_netcdf;
    nmoms   = nmommax;
  } else {
    nmommax = input.mie.nmom;
  }

  pmom_tmp = calloc (nmommax, sizeof (float));

  /* Open netcdf file. */
  if ((status = nc_open (output->outputfile, NC_WRITE, &id))) {
    fprintf (stderr, "Error open netcdf file.\n");
    return -1;
  }

  /* ---------------------------------------------------------------*/
  /* Write variables                                                */
  /* ---------------------------------------------------------------*/

  if (input.verbose)
    fprintf (stderr, "... write data to netcdf file for iv %d, ir %d, ia %d \n", iv, ir, ia);

  for (int ip = 0; ip < nphamat; ip++) {

    start3D[2] = ip;
    start4D[2] = ip;

    if (input.nthetamax > 0) {
      /* ntheta */
      int_array_tmp[0] = output->ntheta[ip];
      if ((retval = nc_put_vars_int (id, ncid.id_ntheta, start3D, count3D, stride3D, int_array_tmp))) {
        fprintf (stderr, "Error writing variable ntheta.\n");
        return -1;
      }

      /* theta */
      count4D[3] = output->ntheta[ip];
      if ((retval = nc_put_vars_float (id, ncid.id_theta, start4D, count4D, stride4D, output->theta[ip]))) {
        fprintf (stderr, "Error writing variable theta.\n");
        return -1;
      }

      /* phase functions */
      if ((retval = nc_put_vars_float (id, ncid.id_phase, start4D, count4D, stride4D, output->phase[ip]))) {
        fprintf (stderr, "Error writing variable phase.\n");
        return -1;
      }
    }

    /* Number of Legendre polynomials in netcdf file*/
    nmoms = output->nmommax;

    if (nmoms == 0) {
      comp_moments[0] = 0;
      nc_put_att_int (id, NC_GLOBAL, "complete_moments", NC_INT, 1, comp_moments);
    }

    /* Number of Legendre polynomials  */
    int_array_tmp[0] = nmoms;
    if ((retval = nc_put_vars_int (id, ncid.id_nmom, start3D, count3D, stride3D, int_array_tmp))) {
      fprintf (stderr, "Error writing variable nmom.\n");
      return -1;
    }

    /* Legendre polynomials */
    count4D[3] = nmommax;
    for (i = 0; i < nmommax; i++)
      pmom_tmp[i] = output->pmom[ip][i] * (2.0 * i + 1.0);
    if ((retval = nc_put_vars_float (id, ncid.id_pmom, start4D, count4D, stride4D, pmom_tmp))) {
      fprintf (stderr, "Error writing variable pmom.\n");
      return -1;
    }
  }

  /* extinction coeffidient */
  double_array_tmp[0] = output->beta / output->rho_medium;
  if ((retval = nc_put_vars_double (id, ncid.id_ext, start2D, count2D, stride2D, double_array_tmp))) {
    fprintf (stderr, "Error writing variable ext.\n");
    return -1;
  }

  /* single scattering albedo */
  double_array_tmp[0] = output->omega;
  if ((retval = nc_put_vars_double (id, ncid.id_ssa, start2D, count2D, stride2D, double_array_tmp))) {
    fprintf (stderr, "Error writing variable ssa.\n");
    return -1;
  }

  /* asymmetry factor */
  if (input.distribution == DIST_NONE)
    double_array_tmp[0] = output->mie.gsca;
  else
    double_array_tmp[0] = output->g;
  if ((retval = nc_put_vars_double (id, ncid.id_gg, start2D, count2D, stride2D, double_array_tmp))) {
    fprintf (stderr, "Error writing variable gg.\n");
    return -1;
  }

  /* refractive index (real) */
  double_array_tmp[0] = output->ref.re;
  if (input.medium != AEROSOL) {
    if ((retval = nc_put_vars_double (id, ncid.id_refre, start1D, count1D, stride1D, double_array_tmp))) {
      fprintf (stderr, "Error writing variable refre.\n");
      return -1;
    }
  } else {
    if ((retval = nc_put_vars_double (id, ncid.id_refre, start2D, count2D, stride2D, double_array_tmp))) {
      fprintf (stderr, "Error writing variable refre.\n");
      return -1;
    }
  }

  /* refractive index (imag) */
  double_array_tmp[0] = output->ref.im;
  if (input.medium != AEROSOL) {
    if ((retval = nc_put_vars_double (id, ncid.id_refim, start1D, count1D, stride1D, double_array_tmp))) {
      fprintf (stderr, "Error writing variable refim.\n");
      return -1;
    }
  } else {
    if ((retval = nc_put_vars_double (id, ncid.id_refim, start2D, count2D, stride2D, double_array_tmp))) {
      fprintf (stderr, "Error writing variable refim.\n");
      return -1;
    }
  }

  /* density */
  double_array_tmp[0] = output->rho_medium;
  if ((retval = nc_put_vars_double (id, ncid.id_rho, start2D, count2D, stride2D, double_array_tmp))) {
    fprintf (stderr, "Error writing variable rho.\n");
    return -1;
  }

  if ((status = nc_close (id))) {
    fprintf (stderr, "Error closing netcdf file.\n");
    return -1;
  }

  return 0;
}

#endif

/**
 * calc_phase_functiom
 *
 * Calculate the phase function from Legendre moments. The phase
 * function is stored on an optimized scattering angle grid. If the
 * accuracy of the phase function is less than 1% (when a specified
 * maximum number of grid points - nthetamax - is used) a warning is
 * printed to the screen. 
 *
 * @author Claudia Emde
 * @date   2010-02-11 restructured 
 */
int calc_phase_function (/* Input */
                         int      nmom,
                         double   lambda,
                         double   reff,
                         int      nstokes,
                         int      nthetamax,
                         int      phase_from_amplitudes,
                         double** phase_in,
                         double*  xmu_in,
                         int      ntheta_in,
                         int      ntheta_user,
                         double   accuracy,
                         /* Output */
                         float*** pmom,
                         float*** phase_out,
                         float*** theta,
                         int**    ntheta,
                         int*     warnings_ntheta) {
  int     ip = 0, it = 0, im = 0;
  double  theta_i = 0.0, dtheta = 0.0; /* accuracy=0.0; */
  double *pmom_tmp = NULL, *mu_tmp = NULL, **phase_tmp = NULL;
  double *mu_opt = NULL, **phase_opt = NULL;
  int     Nopt = 0, Ntheta = 0;

  int status = 0;

  /* check size parameter, if larger than 250, 2 digits in angular 
     grid are not sufficient */

  if (!phase_from_amplitudes) {
    if (2 * M_PI * reff / lambda > 250.0) {
      Ntheta = 180001;
      dtheta = 0.001;
      fprintf (stderr, "3 digits optimization \n");
    } else {
      Ntheta = 18001;
      dtheta = 0.01;
      fprintf (stderr, "2 digits optimization \n");
    }

    /* Calculate phase functions */
    pmom_tmp  = calloc (nmom, sizeof (double));
    mu_tmp    = calloc (Ntheta, sizeof (double));
    phase_tmp = calloc (nstokes, sizeof (double*));
    phase_opt = calloc (nstokes, sizeof (double*));
    /* Allocate maximum allowed number of grid points */
    mu_opt = calloc (nthetamax, sizeof (double));
    for (ip = 0; ip < nstokes; ip++) {
      phase_tmp[ip] = calloc (Ntheta, sizeof (double));
      phase_opt[ip] = calloc (nthetamax, sizeof (double));
    }

    for (ip = 0; ip < nstokes; ip++) {

      for (im = 0; im < nmom; im++)
        pmom_tmp[im] = (2 * im + 1) * (*pmom)[ip][im];

      for (it = 0; it < Ntheta; it++) {
        theta_i           = 180.0 - it * dtheta;
        mu_tmp[it]        = cos (theta_i / 180.0 * M_PI);
        phase_tmp[ip][it] = mom2phase (cos (theta_i / 180.0 * M_PI), pmom_tmp, nmom);
      }
    }
    free (pmom_tmp);
  } else {
    mu_tmp = calloc (ntheta_in, sizeof (double));
    mu_opt = calloc (nthetamax, sizeof (double));

    phase_tmp = calloc (nstokes, sizeof (double*));
    phase_opt = calloc (nstokes, sizeof (double*));

    for (ip = 0; ip < nstokes; ip++) {
      phase_tmp[ip] = calloc (ntheta_in, sizeof (double));
      phase_opt[ip] = calloc (nthetamax, sizeof (double));
    }

    for (int iang = 0; iang < ntheta_in; iang++) {
      mu_tmp[iang] = xmu_in[iang];
      for (ip = 0; ip < nstokes; ip++)
        phase_tmp[ip][iang] = phase_in[ip][iang];
    }

    Ntheta = ntheta_in;
  }

  /* Accuracy of phase function: 1% */
  /* accuracy=0.01; */
  /* Optimize angular grid for phase function */
  if (ntheta_user == 0) {
    status = optimize_theta_grid_phamat (&mu_opt,
                                         &phase_opt,
                                         &Nopt,
                                         mu_tmp,
                                         phase_tmp,
                                         Ntheta,
                                         accuracy,
                                         nthetamax,
                                         nstokes,
                                         &(*warnings_ntheta));

    if (status != 0) {
      fprintf (stderr, "Error %d returned by optimize_theta_grid()\n", status);
      return status;
    }
  } else {
    Nopt = Ntheta;
    for (int i = 0; i < Nopt; i++) {
      for (ip = 0; ip < nstokes; ip++)
        phase_opt[ip][i] = phase_tmp[ip][i];
      mu_opt[i] = mu_tmp[i];
    }
  }

  for (ip = 0; ip < nstokes; ip++) {
    (*ntheta)[ip] = Nopt;
    for (it = 0; it < Nopt; it++) {
      /* At the moment the maximum accuracy is 3 digits. */
      /* Print a warning if this accuracy is not sufficient */
      (*theta)[ip][it]     = acos (mu_opt[it]) * 180.0 / M_PI;
      (*phase_out)[ip][it] = phase_opt[ip][it];
      //      if ((*theta)[ip][it] == 0.001 && ip == 0)    // BM, 24.3.2024; replaced because compiler warned "floating-point comparison is always false"
      if (fabs ((*theta)[ip][it] - 0.001) < 1E-6 && ip == 0)
        fprintf (stderr, "Warning: The accuracy of the angular grid is not sufficient !!!!\n)");
    }
  }

  for (ip = 0; ip < nstokes; ip++) {
    free (phase_tmp[ip]);
    free (phase_opt[ip]);
  }

  free (mu_tmp);
  free (phase_tmp);
  free (mu_opt);
  free (phase_opt);

  return 0;
}

int setup_alpha (double** alpha, int* nalpha, char* alpha_filename, double input_alpha, int verbose) {
  int status = 0, ia = 0;

  if (strlen (alpha_filename) != 0) {

    status = read_1c_file (alpha_filename, &(*alpha), &(*nalpha));
    if (verbose)
      fprintf (stderr, "... read %d data points from %s\n", *nalpha, alpha_filename);

    if (status != 0) {
      fprintf (stderr, "error %d reading %s\n", status, alpha_filename);
      return status;
    }
  } else {
    *alpha      = calloc (1, sizeof (double));
    (*alpha)[0] = input_alpha;
    *nalpha     = 1;
  }

  if (verbose) {
    fprintf (stderr, "Size distibution parameter alpha: ");

    for (ia = 0; ia < (*nalpha); ia++)
      fprintf (stderr, "%.1f ", (*alpha)[ia]);
    fprintf (stderr, "\n");
  }

  return 0;
}

static int sortasc (const void* x1, const void* x2) {
  double* fx1 = (double*) x1;
  double* fx2 = (double*) x2;

  if (*fx1 > *fx2)
    return 1;

  if (*fx1 < *fx2)
    return -1;

  return 0;
}
