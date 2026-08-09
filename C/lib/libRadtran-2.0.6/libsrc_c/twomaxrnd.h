/*--------------------------------------------------------------------
 * $Id: twomaxrnd.c 2623 2011-12-23 10:52:38Z nina.crnivec $
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

/* full solution of the multi-layer twostream equation solar + thermal */
/* with maximum-random overlap assumption for partial cloudiness */
/* written by Nina Crnivec */

#ifndef __twomaxrnd_h
#define __twomaxrnd_h

#if defined(__cplusplus)
extern "C" {
#endif

int twomaxrnd (float* dtau_org,
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
               float* uavg);

void eddington_coeffc (double  dtau,
                       double  g,
                       double  omega0,
                       double  mu0,
                       double* a11,
                       double* a12,
                       double* a13,
                       double* a23,
                       double* a33);

void calcThermalComponents (int ilyr, double* B, double dtau, double omega0, double g, double* theComp1, double* theComp2);

int makeMatrixA2 (int nlev, int nFluxPerLevel, double** matrixA);

int makeVectorMinusB (int nlev, int nFluxPerLevel, double* vectB);

#if defined(__cplusplus)
extern "C" {
#endif

#endif /* _twomaxrnd_h */
