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
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "common_math.h"
#include "equation.h"
#include "twostrebe.h"
#include "dynamic_twostream.h"
#include "solver.h"

/* Function that calculates the up- and downward thermal emission for
   all atmospheric layers similar to the way in Bernhard Mayer's twostrebe
   algorithm
*/
static void calculate_thermal_emission (double  B0,
                                        double  B1,
                                        double  dtau_ext,
                                        double  g_Parameter,
                                        double  omega_0,
                                        double* B_up,
                                        double* B_down) {

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

  *B_up   = -a11 * (alpha7 + alpha8 * dtau_ext) - a12 * (alpha9) + alpha7;
  *B_down = -a12 * (alpha7 + alpha8 * dtau_ext) - a11 * (alpha9) + alpha9 + alpha10 * dtau_ext;
}

/* Function that applies one Gauss-Seidel iteration to update the radiative quantities
*/
static void one_iteration_forward (int     it,
                                   int     nlyr,
                                   float*  dtau_org,
                                   float*  omega0_org,
                                   float*  g_org,
                                   double  S0,
                                   double  mu0,
                                   double  Ag,
                                   int     planck,
                                   int     delta,
                                   double* B,
                                   double  Bg,
                                   double* E_dir,
                                   double* E_down,
                                   double* E_up,
                                   double* E_down_hr,
                                   double* divE) {

  /*
  Definition of fields for layer reflectivities (r) and transmissions (t), upward (r_dir) and downward (s_dir) 
  scattered as well as transmitted (t_dir) fractions of direct solar radiation
  */
  double r[nlyr], t[nlyr], r_dir[nlyr], t_dir[nlyr], s_dir[nlyr];
  double B_up[nlyr], B_down[nlyr];
  {
    const double taumax = 100.0;  // maximum layer optical thickness allowed to avoid problems with the Eddington code

    // Boundary Conditions I: Top of the Atmosphere
    E_dir[0]  = S0 * mu0;
    E_down[0] = 0.0;

    /* Swap of the iteration direction in every Gauß-Seidel iteration: If we are in the thermal spectral range, 
       we start to iterate from bottom to top and alternate that iteration direction from there on; in the 
       solar spectral range, we start iterating from top to bottom.
    */
    unsigned int i_start, nit, i_increment;
    if ((planck && it % 2 == 1) || (!planck && it % 2 == 0)) {
      i_start     = 0;
      i_increment = 1;
    } else {
      i_start     = nlyr - 1;
      i_increment = -1;
    }
    nit = nlyr;

    for (unsigned int i = i_start; i < nit; i += i_increment) {

      double dtau = MIN (taumax, dtau_org[i]);
      double w0   = omega0_org[i];
      double g    = g_org[i];

      /*
      Get values for the layer reflectivities (r) and transmissions (t), upward (r_dir) and downward (s_dir) scattered 
      as well as transmitted (t_dir) fractions of direct solar radiation using the eddington_v2 code
      */
      if (delta) {
        delta_scale_hg (dtau, w0, g, &dtau, &w0, &g);
      }
      eddington_v2 (dtau, g, w0, mu0, &(t[i]), &(r[i]), &(r_dir[i]), &(s_dir[i]), &(t_dir[i]));

      // Calculate up- and downward facing thermal emission for all atmopsheric layers
      if (planck) {
        double B0;
        double B1;
        if (dtau > 0.01) {
          B1 = (B[i + 1] - B[i]) / dtau;
          B0 = B[i];
        } else {
          B1 = 0.0;
          B0 = (B[i + 1] + B[i]) * 0.5;
        }
        calculate_thermal_emission (B0, B1, dtau, g, w0, &(B_up[i]), &(B_down[i]));
      } else {
        B_up[i]   = 0.0;
        B_down[i] = 0.0;
      }

      // Dynamic Radiative Transfer: For every grid box, new output irradiances fluxes are calculated using the input irradiances of the previous timestep
      E_dir[i + 1]  = E_dir[i] * t_dir[i];
      E_down[i + 1] = t[i] * E_down[i] + r[i] * E_up[i + 1] + s_dir[i] * E_dir[i] + B_down[i];
      E_up[i]       = t[i] * E_up[i + 1] + r[i] * E_down[i] + r_dir[i] * E_dir[i] + B_up[i];
    }

    // Boundary Conditions II: Surface
    E_up[nlyr] = Ag * (E_down[nlyr] + E_dir[nlyr]) + (1.0 - Ag) * M_PI * Bg;
  }

  // Calculation of physically incorrect downward irradiances that can be used to physically correct calculate heating rates using the classical
  // forward differences method implemented in libRadtran
  E_down_hr[0] = 0.0;
  for (int i = 0; i < nlyr; i++) {
    E_down_hr[i + 1] = E_down_hr[i] + (E_dir[i] - E_dir[i + 1]) + (E_up[i + 1] - E_up[i]) -
                       (1.0 - t[i] - r[i]) * (E_down[i] + E_up[i + 1]) - (1.0 - t_dir[i] - r_dir[i] - s_dir[i]) * E_dir[i] +
                       (B_down[i] + B_up[i]);

    divE[i] =
      (1.0 - t[i] - r[i]) * (E_down[i] + E_up[i + 1]) + (1.0 - t_dir[i] - r_dir[i] - s_dir[i]) * E_dir[i] - (B_down[i] + B_up[i]);
  }
}

int dynamic_twostream (int    n_iterations,
                       float* dtau_org,
                       float* omega0_org,
                       float* g_org,
                       float* dtau_clr_org,
                       float* omega0_clr_org,
                       float* g_clr_org,
                       float* cf,
                       int    nlev,
                       double S0,
                       double mu0,
                       double Ag,
                       int    planck,
                       int    delta,
                       int    nzout,
                       float* zd,
                       float* temper,
                       float  btemp,
                       float  wvnmlo,
                       float  wvnmhi,
                       float* zout,
                       float* fldn,
                       float* flup,
                       float* fldir,
                       float* uavg) {

  fprintf (stderr,
           "ATTENTION: In the moment, the DYNAMIC TWOSTREAM solver IGNORES ANY CLOUD FRACTIONS specified. Hence, all clouds extend "
           "throughout the whole respective model layer. \n");

  const int nlyr = nlev - 1;  // number of atmospheric model layers

  // Setup of Radiative and Optical Quantities
  double E_down[nlev], E_up[nlev], E_dir[nlev], E_down_hr[nlev], divE[nlyr];

  for (int i = 0; i < nlev; i++) {
    E_down[i] = 0.0;
    E_up[i]   = 0.0;
    E_dir[i]  = 0.0;
  }

  // Planck Setup
  double B[nlev];   // emitted black body radiance at all atmospheric levels within the given wavelength spectrum
  double Bg = 0.0;  // black body radiance by the ground within the given wavelength spectrum

  for (int i = 0; i < nlev; i++) {
    float plkavg = 0.0;
    if (planck) {
      F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &temper[i], &plkavg);
    }
    B[i] = plkavg;
  }
  if (planck) {
    float plkavg;
    F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &btemp, &plkavg);
    Bg = plkavg;
  }

  for (int k = 0; k < n_iterations; k++) {
    one_iteration_forward (k,
                           nlyr,
                           dtau_org,
                           omega0_org,
                           g_org,
                           S0,
                           mu0,
                           Ag,
                           planck,
                           delta,
                           B,
                           Bg,
                           E_dir,
                           E_down,
                           E_up,
                           E_down_hr,
                           divE);
  }

  /* copy results to final fields */
  /* When heating_rate layer_fd is specified, heating rates are calculated via flux divergence and thus have to use
     the adjusted downward flux E_down_hr. heating_rate local works as well, as libRadtran then expects divE in uavg
     for rte_solver dynamic_twostream
  */
  int i_last_zout = 0;
  for (int i = 0; i < nlev; i++) {
    for (int j = 0; j < nzout; j++) {
      if (zout[j] == zd[i]) {
        fldn[j]  = E_down_hr[i];
        flup[j]  = E_up[i];
        fldir[j] = E_dir[i];
        uavg[j]  = 0.0;
        for (int l = i_last_zout; l < i; l++) {
          uavg[j] += divE[l];
        }
        i_last_zout = i;
      }
    }
  }

  fprintf (stderr, " ... *** Richard Maier's dynamic twostream solver with %d iterations\n", n_iterations);

  return 0;
}
