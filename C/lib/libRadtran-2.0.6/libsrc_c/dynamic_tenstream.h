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

#ifndef __dynamic_tenstream_h
#define __dynamic_tenstream_h

#if defined(__cplusplus)
extern "C" {
#endif

#include "mystic.h"

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
                       float***            abs3d);

/* Function returning ecRad Eddington coefficients as implemented in eddington.F90 in the
   TenStream solver
 */
int eddington_ec (double dtau, double g, double omega0, double mu0, double* t, double* r, double* rdir, double* sdir, double* tdir);

#if defined(__cplusplus)
extern "C" {
#endif

#endif /* __dynamic_tenstream_h */
