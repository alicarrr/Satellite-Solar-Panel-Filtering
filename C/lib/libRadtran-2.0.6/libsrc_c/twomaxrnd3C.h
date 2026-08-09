/*--------------------------------------------------------------------
 * $Id: twostrebe.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __twomaxrnd3c_h
#define __twomaxrnd3c_h

#if defined(__cplusplus)
extern "C" {
#endif

int twomaxrnd3C (float* dtau_cldk,
                 float* omega0_cldk,
                 float* g_cldk,
                 float* dtau_cldn,
                 float* omega0_cldn,
                 float* g_cldn,
                 float* dtau_clr,
                 float* omega0_clr,
                 float* g_clr,
                 float* cf,
                 float  scale_cf,
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
                 float  wvnmup,
                 float* zout,
                 float* fldn,
                 float* flup,
                 float* fldir,
                 float* uavg);

#if defined(__cplusplus)
extern "C" {
#endif

#endif /* _twomaxrnd3c_h */
