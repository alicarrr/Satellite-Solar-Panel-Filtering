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

#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "errors.h"
#include "uvspecrandom.h"

#if HAVE_LIBGSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#endif

/* global variable: random number generator */
#if HAVE_LIBGSL
gsl_rng* uvspecrng;
#endif

/*****************************************************************/
/* uvspec random number generator; either the highly-recommended */
/* MT19937 if GSL is available, or random() otherwise.           */
/*****************************************************************/

inline double uvspec_random() {
#if HAVE_LIBGSL
  return gsl_rng_uniform (uvspecrng);
#else
  return (double) random();
#endif
}

double uvspec_random_gauss (double sigma) {
#if HAVE_LIBGSL
  return gsl_ran_gaussian (uvspecrng, sigma);
#else
  return (double) random();
#endif
}

int init_uvspec_random (int* randomseed, int need_mt19937, int quiet) {
  int rseed = 0;

/* global variable, defined in uvspec_lex.l */
#if HAVE_LIBGSL
  /*  extern gsl_rng *uvspecrng; */ /* global random number generator */
  const gsl_rng_type* T;
#endif

  /* always initialize standard rng */
  srandom (rseed);

  if (*randomseed > 0)
    rseed = *randomseed;
  else
    rseed = (int) time (NULL) + (int) getpid();

#if HAVE_LIBGSL
  T                    = gsl_rng_mt19937;
  gsl_rng_default_seed = rseed;
  uvspecrng            = gsl_rng_alloc (T);
  if (!quiet && need_mt19937)
    fprintf (stderr, " ... using %s random number generator\n", gsl_rng_name (uvspecrng));
#else
  if (need_mt19937) {
    fprintf (stderr, "Error, GSL not found. The MT19937 random number generator is\n");
    fprintf (stderr, "required for MYSTIC simulations. Please install GSL and recompile!\n");
    return -1;
  }
#endif

  return rseed;
}

/***********************************************************************************/
/* Function: read_and_write_random_status                                 @62_30i@ */
/* Description:                                                                    */
/*  check whether to read/write random status.                                     */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

#if HAVE_LIBGSL
int read_and_write_random_status (char* randomstatusfilename,
                                  int*  ili,
                                  int*  backward_is,
                                  int*  backward_js,
                                  int   photoncounter,
                                  int*  readrandomstatus) {
  FILE* randomstatusfile = NULL;
  int   status           = 0;

  /************************************+**************************/
  /* starting with the randomseed is nice, but starting with the */
  /* actual photon that produced the problem is even nicer...    */
  /* here we read and write the status of the random number      */
  /* generator                                                   */
  /* this only works with the gsl !!!                            */
  /************************************+**************************/

  /* read randomstatusfile if asked for */
  if (*readrandomstatus && photoncounter == 1) {
    if ((randomstatusfile = fopen (randomstatusfilename, "r")) != NULL) {
      status = gsl_rng_fread (randomstatusfile, uvspecrng);
      if (status)
        return err_out ("Error %d returned by gsl_rng_fread()\n", status);
      if (!fscanf (randomstatusfile, "%d %d %d", ili, backward_is, backward_js)) {
        fprintf (stderr, "Error reading file %s\n", randomstatusfilename);
        return -1;
      }
      (void) fclose (randomstatusfile);
      fprintf (stderr, "reading random status for photon...\n");
    } else
      return err_out ("Error %d, could not open randomstatusfile\n", -1);
  }

#ifdef WRITERANDOMSTATUS
  /* write status every 10 photons, starting with photon 1 */
  if (photoncounter % 10 == 1 && !(*readrandomstatus)) {
    if ((randomstatusfile = fopen (randomstatusfilename, "w")) != NULL) {
      status = gsl_rng_fwrite (randomstatusfile, uvspecrng);
      if (status)
        return err_out ("Error %d returned by gsl_rng_fwrite()\n", status);
      fprintf (randomstatusfile, "%d %d %d\n", *ili, *backward_is, *backward_js);
      (void) fclose (randomstatusfile);
    } else
      return err_out ("Error %d, could not open randomstatusfile\n", -1);
  }
#endif
  *readrandomstatus = 0;

  return 0;
}
#endif
