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

#ifndef __uvspecrandom_h
#define __uvspecrandom_h

#if defined(__cplusplus)
extern "C" {
#endif

#if HAVE_LIBGSL
#include <gsl/gsl_rng.h>
#endif

/* for debugging. This mode writes the randomstatus file, with */
/* help of it one can start shortly before a buggy photon      */
/* instead of having to simulate many good photons before      */
/* reaching the problematic one                                */
#undef WRITERANDOMSTATUS

double uvspec_random();
double uvspec_random_gauss (double sigma);
int    init_uvspec_random (int* randomseed, int need_mt19937, int quiet);

int read_and_write_random_status (char* randomstatusfilename,
                                  int*  ili,
                                  int*  backward_is,
                                  int*  backward_js,
                                  int   photoncounter,
                                  int*  readrandomstatus);

#if defined(__cplusplus)
}
#endif

#endif
