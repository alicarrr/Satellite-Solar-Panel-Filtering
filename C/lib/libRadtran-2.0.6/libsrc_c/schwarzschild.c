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

/* full solution of the Schwarzschild equation */

#include <math.h>

#include "cdisort.h"
#include "schwarzschild.h"

// solve schwarzschild equation along a path tau where L is the initial radiance, B1 is planck emission at the start and B2 is the emission at the end of the path
double schwarzschild_radiance (double L, double tau, double B1, double B2) {
  if (tau > 1e-30) {
    const double tm1 = expm1 (-tau);
    return L * (tm1 + 1) + (B2 - B1) - (B1 - (B2 - B1) / tau) * tm1;
  } else {
    return (B1 + B2) * .5 * tau + L * (1 - tau);
  }
}

int schwarzschild (const int   Nstreams,
                   const float Ag,
                   float*      dtau,
                   float*      ssa,
                   int         nlev,
                   int         planck,
                   int         nzout,
                   float*      zd,
                   float*      temper,
                   float       btemp,
                   float       wvnmlo,
                   float       wvnmhi,
                   float*      zout,
                   float*      fldn,
                   float*      flup,
                   float*      fldir,
                   float*      uavg) {

  if (planck) {
    const int nlyr = nlev - 1;
    double    tauabs[nlev - 1];

    double B[nlev];
    double E_up[nlev];
    double E_dn[nlev];

    const double dmu = 1.0 / (double) Nstreams;

    for (int k = 0; k < nlyr; k++) {  // only consider absorption optical depth
      tauabs[k] = (dtau[k] * (1 - ssa[k]));
    }

    for (int k = 0; k < nlev; k++) { /* level temperatures and Planck functions */
      E_up[k] = 0;
      E_dn[k] = 0;
      B[k]    = c_planck_func1 (wvnmlo, wvnmhi, temper[k]);
    }
    /* surface temperature and Planck function */
    const double Bg = c_planck_func1 (wvnmlo, wvnmhi, btemp);

    for (int imu = 0; imu < Nstreams; imu++) {
      const double mu = (imu + 0.5) * dmu;

      /* ## Boundary conditions at TOA ## */
      double L_dn = 0;
      // E_dn[0] = 0;

      /* ## Calculate downwelling radiation ## */
      for (int k = 0; k < nlyr; k++) {
        L_dn = schwarzschild_radiance (L_dn, tauabs[k] / mu, B[k], B[k + 1]);
        E_dn[k + 1] += L_dn * mu;
      }
    }

    for (int imu = 0; imu < Nstreams; imu++) {
      const double mu = (imu + 0.5) * dmu;
      /* ## Boundary conditions at Ground, assuming lambertian surface ## */
      double L_up = (1 - Ag) * Bg + Ag * E_dn[nlev - 1] * 2 * dmu;
      E_up[nlev - 1] += L_up * mu;

      /* ## Calculate upwelling radiation ## */
      for (int k = (nlyr - 1); k >= 0; k--) {
        L_up = schwarzschild_radiance (L_up, tauabs[k] / mu, B[k + 1], B[k]);
        E_up[k] += L_up * mu;
      }
    }

    /*copy results to final fields*/
    for (int k = 0; k < nlev; k++) {
      for (int lu = 0; lu < nzout; lu++) {
        if (fabs (zout[lu] - zd[k]) <= 1E-6 * zout[lu]) {
          fldn[lu]  = 2 * M_PI * dmu * E_dn[k];
          flup[lu]  = 2 * M_PI * dmu * E_up[k];
          fldir[lu] = 0;
          uavg[lu]  = 0;
        }
      }
    }

  } else {
    for (int lu = 0; lu < nzout; lu++) {
      fldn[lu]  = 0.0;
      flup[lu]  = 0.0;
      fldir[lu] = 0;
      uavg[lu]  = 0;
    }
  }

  return 0;
}
