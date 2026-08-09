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

#include "uvspec.h"
#include "fortran_and_c.h"
#include "f77-uscore.h"
#include "errors.h"
#include "ckdfu.h"

/* molecular constants */
#define M_AIR 28.977
#define M_O3 48.000
#define M_H2O 18.015

#define N_MOL_CKDFU 6.0235e23

int ckdfu (float       umco2,
           float       umch4,
           float       umn2o,
           float       umf11,
           float       umf12,
           float       umf22,
           float***    temper,
           float***    press,
           float*      z,
           float****   dens,
           float****   dens_avg,
           int         nlev,
           int         Nx,
           int         Ny,
           int         h2o_cont,
           int         lower_wl_id,
           int         upper_wl_id,
           int         quiet,
           ck_profile* ck) {

#include "fl_radparams.cinc"

  float* tmp_ck   = calloc (nv1x * mbx * mg, sizeof (float));
  float* tmp_delg = calloc (mbx * mg, sizeof (float));

  int* kg = calloc (mbx, sizeof (int));

  float *rhoair = NULL, *rhoo3 = NULL, *rhoh2o = NULL;
  float *pres = NULL, *temp = NULL, *davg = NULL;

  int first = 1;

  /* allocate memory for density profiles */
  rhoair = calloc (nlev, sizeof (float));
  rhoh2o = calloc (nlev, sizeof (float));
  rhoo3  = calloc (nlev, sizeof (float));
  temp   = calloc (nlev, sizeof (float));
  pres   = calloc (nlev, sizeof (float));
  davg   = calloc (nlev, sizeof (float));

  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      for (int lc = 0; lc < nlev; lc++) {
        rhoair[lc] = dens[MOL_AIR][ix][iy][lc] * M_AIR / N_MOL_CKDFU * 1.0e6;
        rhoh2o[lc] = dens[MOL_H2O][ix][iy][lc] * M_H2O / N_MOL_CKDFU * 1.0e6;
        rhoo3[lc]  = dens[MOL_O3][ix][iy][lc] * M_O3 / N_MOL_CKDFU * 1.0e6;
        temp[lc]   = temper[ix][iy][lc];
        pres[lc]   = press[ix][iy][lc];
      }

      for (int lc = 0; lc < nlev - 1; lc++)
        davg[lc] = dens_avg[MOL_AIR][ix][iy][lc] * M_AIR / N_MOL_CKDFU * 1.0e6 * (z[lc] - z[lc + 1]);

      void F77_FUNC (ckdfuf, CKDFUF) (float* umco2,
                                      float* umch4,
                                      float* umn2o,
                                      float* umf11,
                                      float* umf12,
                                      float* umf22,
                                      float* rhoair,
                                      float* rhoh2o,
                                      float* rhoo3,
                                      float* pres,
                                      float* temp,
                                      float* davg,
                                      float* z,
                                      int*   nlev,
                                      int*   quiet,
                                      int*   h2o_cont,
                                      int*   kg,
                                      float* tmp_delg,
                                      float* tmp_ck);

      F77_FUNC (ckdfuf, CKDFUF)
      (&umco2,
       &umch4,
       &umn2o,
       &umf11,
       &umf12,
       &umf22,
       rhoair,
       rhoh2o,
       rhoo3,
       pres,
       temp,
       davg,
       z,
       &nlev,
       &h2o_cont,
       &quiet,
       kg,
       tmp_delg,
       tmp_ck);

      /* we don't want first to be static because we want to allocate memory each time ckdfu() is called  */
      /* BM 1.4.2022: Do we really want that in the loop? Why not move it to front and remove if (first)? */
      if (first) {
        first = 0;

        for (int ib = lower_wl_id; ib <= upper_wl_id; ib++) {
          ck[ib].weight = calloc (kg[ib], sizeof (double));
          ck[ib].crs    = calloc (Nx, sizeof (float***));
          for (int ixx = 0; ixx < Nx; ixx++) {
            ck[ib].crs[ixx] = calloc (Ny, sizeof (float**));
            for (int iyy = 0; iyy < Ny; iyy++) {
              ck[ib].crs[ixx][iyy] = calloc (nlev, sizeof (float*));
              for (int lc = 0; lc < nlev; lc++)
                ck[ib].crs[ixx][iyy][lc] = calloc (kg[ib], sizeof (float));
            }
          }
        }
      }

      /* Loop over bands */
      for (int ib = lower_wl_id; ib <= upper_wl_id; ib++) {

        /* copy number of subbands to final destination */
        ck[ib].ngauss = kg[ib];

        for (int ig = 0; ig < kg[ib]; ig++) {

          /* copy weight of subbands to final destination */
          ck[ib].weight[ig] = tmp_delg[ig * mbx + ib];

          /* copy absorption coefficient to final destination */
          for (int lu = 0; lu < nlev - 1; lu++) {
            ck[ib].crs[ix][iy][lu][ig] = tmp_ck[ig * nvx * mbx + lu * mbx + ib];
            if (isinf (ck[ib].crs[ix][iy][lu][ig])) {
              ck[ib].crs[ix][iy][lu][ig] = 0.0;

              if (!quiet)
                fprintf (stderr, "Warning, ckdfu() returned INF in band %d, setting to 0\n", ib + 1);
            }
          }
        }
      }
    }
  }

  free (tmp_ck);
  free (tmp_delg);
  free (kg);

  free (rhoair);
  free (rhoh2o);
  free (rhoo3);
  free (temp);
  free (pres);
  free (davg);

  return 0;
}

int determine_ckdfu_band_loop_indices (int ib, int* ib_start, int* ib_end) {
#include "fl_radparams.cinc"

  /* the following statement is meaningless but it prevents compiler warnings */
  int null = nvx + nv1x + mbx + mg - nvx - nv1x - mbx - mg;
  null++;

  if (ib == NOT_DEFINED_INTEGER) {
    *ib_start = 1;
    *ib_end   = mbx;
  } else {
    if (ib < 0 || ib >= mbx) {
      int status = ib;
      fct_err_out (status, "wrong band ib to compute.. has to be between 0 and 17", ERROR_POSITION);
      return ib;
    }
    *ib_start = ib + 1;  // note here that this is a fortran loop index which begins at 1 and includes the upper bound
    *ib_end   = ib + 1;
  }

  return 0;
}

/***********************************************************************************/
/* Function: ckdfucld                                                              */
/* Description:                                                                    */
/*  compute cloud liquid water optical properties with code from Fu&Liou           */
/* Parameters: ib specifies a particular correlated-k band.                        */
/*             If you provide ib == NOT_DEFINED_INTEGER,                           */
/*             then we compute optical properties for all bands                    */
/***********************************************************************************/
int ckdfucld (float* z, float* reff, float* lwc, int nlev, int ib, float** tau, float** gg, float** ssa) {

#include "fl_radparams.cinc"

  /* the following statement is meaningless but it prevents compiler warnings */
  /* "warning: variable ‘mg’ set but not used" etc                            */
  int null = nvx + nv1x + mbx + mg - nvx - nv1x - mbx - mg;
  null++;

  float tmp_tau[nv1x * mbx];
  float tmp_ssa[nv1x * mbx];
  float tmp_ww[nv1x * mbx * 4];

  int ib_start = 0, ib_end = 0;

  int ierr = determine_ckdfu_band_loop_indices (ib, &ib_start, &ib_end);
  if (ierr)
    fct_err_out (ierr, "ckdfucld -- wrong band ib to compute..", ERROR_POSITION);

  void F77_FUNC (ckdcldf, CKDCLDF) (float* z,
                                    float* reff,
                                    float* lwc,
                                    int*   nlev,
                                    int*   ib_start,
                                    int*   ib_end,
                                    float* tmp_tau,
                                    float* tmp_ww,
                                    float* tmp_ssa);

  for (int lu = 0; lu < nlev; lu++) {
    if (lwc[lu] > 0) { /* need only test if lwc>0; otherwise optical properties are not required */

      /* "experimentally" determined valid ranges by comparison with Mie and geometrical optics, BM 16.1.2018 */
      /* different subroutines are used for solar and thermal: fl_misc_subs.f:water() for thermal,            */
      /* fl_water_hu2.f:water_hu() for solar. Limits:                                                         */
      /* solar: 2.5 - 40 micron, thermal 4.18 - 31.23 micron (for thermal see comment in fl_misc_subs.f)      */

      if (reff[lu] < 4.18 || reff[lu] > 31.23) {
        fprintf (stderr, "Error, effective radius %f um not covered by Fu and Liou model.\n", reff[lu]);
        fprintf (stderr, "Allowed range is %7.3f - %7.3f um\n", 4.18, 31.23);
        return -1;
      }
    }
  }

  F77_FUNC (ckdcldf, CKDCLDF) (z, reff, lwc, &nlev, &ib_start, &ib_end, tmp_tau, tmp_ww, tmp_ssa); /* in libsrc_c/ckdcldf.f */

  /* Loop over bands */
  for (int ibb = ib_start - 1; ibb < ib_end; ibb++) {

    /* copy data to final destination */
    for (int lu = 0; lu < nlev - 1; lu++) {
      tau[ibb][lu] = tmp_tau[lu * mbx + ibb];
      ssa[ibb][lu] = tmp_ssa[lu * mbx + ibb];
      gg[ibb][lu]  = tmp_ww[lu * mbx * 4 + 4 * ibb + 0] / 3.0; /* first moment of the phase function */
    }
  }

  return 0;
}

/***********************************************************************************/
/* Function: ckdfuice                                                              */
/* Description:                                                                    */
/*  compute cloud ice water optical properties with code from Fu&Liou              */
/* Parameters: ib specifies a particular correlated-k band.                        */
/*             If you provide ib == NOT_DEFINED_INTEGER,                           */
/*             then we compute optical properties for all bands                    */
/***********************************************************************************/
int ckdfuice (float* z, float* reff, float* lwc, int nlev, int ib, float** tau, float** gg, float** ssa, int unscaled) {

#include "fl_radparams.cinc"

  /* the following statement is meaningless but it prevents compiler warnings */
  /* "warning: variable ‘mg’ set but not used" etc                            */
  int null = nvx + nv1x + mbx + mg - nvx - nv1x - mbx - mg;
  null++;

  float tmp_tau[nv1x * mbx];
  float tmp_ssa[nv1x * mbx];
  float tmp_ww[nv1x * mbx * 4];
  float deff[nlev];

  int ib_start = 0, ib_end = 0;

  int ierr = determine_ckdfu_band_loop_indices (ib, &ib_start, &ib_end);
  if (ierr)
    fct_err_out (ierr, "ckdfuice -- wrong band ib to compute..", ERROR_POSITION);

  void F77_FUNC (ckdicef, CKDICEF) (float* z,
                                    float* reff,
                                    float* lwc,
                                    int*   nlev,
                                    int*   ib_start,
                                    int*   ib_end,
                                    float* tmp_tau,
                                    float* tmp_ww,
                                    float* tmp_ssa,
                                    int*   unscaled); /* in src/ckdicef.f */

  /* Fu et al. requires effective diameter */
  for (int lu = 0; lu < nlev; lu++) {
    deff[lu] = 2.0 * reff[lu];

    /* test if effective diameter within allowed range; unfortunately we need */
    /* to exclude more than necessary: shortwave and longwave properties are  */
    /* valid for slightly different size ranges; however, both are calculated */
    /* at the same time, hence we can only use the region which is available  */
    /* in both parameterizations.                                             */

    if (lwc[lu] > 0) { /* need only test if lwc>0; otherwise optical properties are not required */
      if (deff[lu] < MIN_DEFF_FU96 || deff[lu] > MAX_DEFF_FU96) {
        fprintf (stderr, "Error, effective radius %f um not covered by Fu [1996]\n", deff[lu] / 2.0);
        fprintf (stderr, "Allowed range is %7.3f - %7.3f um\n", MIN_DEFF_FU96 / 2.0, MAX_DEFF_FU96 / 2.0);
        return -1;
      }

      if (deff[lu] < MIN_DEFF_FU98 || deff[lu] > MAX_DEFF_FU98) {
        fprintf (stderr, "Error, effective radius %f um not covered by Fu et al. [1998]\n", deff[lu] / 2.0);
        fprintf (stderr, "Allowed range is %7.3f - %7.3f um\n", MIN_DEFF_FU98 / 2.0, MAX_DEFF_FU98 / 2.0);
        return -1;
      }
    }
  }

  F77_FUNC (ckdicef, CKDICEF) (z, deff, lwc, &nlev, &ib_start, &ib_end, tmp_tau, tmp_ww, tmp_ssa, &unscaled);

  /* Loop over bands */
  for (int ibb = ib_start - 1; ibb < ib_end; ibb++) {

    /* copy data to final destination */
    for (int lu = 0; lu < nlev - 1; lu++) {
      tau[ibb][lu] = tmp_tau[lu * mbx + ibb];
      ssa[ibb][lu] = tmp_ssa[lu * mbx + ibb];
      gg[ibb][lu]  = tmp_ww[lu * mbx * 4 + 4 * ibb + 0] / 3.0; /* first moment of the phase function */
    }
  }

  return 0;
}

int ckdfuray (int ib, int ig, float u0, int nlev, float* tau) {

#include "fl_radparams.cinc"

  /* the following statement is meaningless but it prevents compiler warnings */
  /* "warning: variable ‘mg’ set but not used" etc                            */
  int null = nvx + nv1x + mbx + mg - nvx - nv1x - mbx - mg;

  void F77_FUNC (ckdrayf, CKDRAYF) (int* ib, int* ig, float* u0, int* nlev, float* tau);

  F77_FUNC (ckdrayf, CKDRAYF) (&ib, &ig, &u0, &nlev, tau);

  return null;
}
