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

#ifndef __specrend_uvspec_h
#define __specrend_uvspec_h

#if defined(__cplusplus)
extern "C" {
#endif

/* prototypes */
int spectrum_to_rgb_overwrite (float* spec, int nspec);
int spectrum_to_rgb (float* spec, int nspec, double* r, double* g, double* b);

int image_wide_rgb_correction1d (float** image, int nzout);
int image_wide_rgb_correction2d (float*** image, int nzout, int iupper);
int image_wide_rgb_correction3d (float**** image, int nzout, int ilower, int iupper, int istep, int jlower, int jupper, int jstep);
int image_wide_rgb_correction3drad (float****** image,
                                    int         nzout,
                                    int         ilower,
                                    int         iupper,
                                    int         istep,
                                    int         jlower,
                                    int         jupper,
                                    int         jstep);
#if defined(__cplusplus)
}
#endif

#endif
