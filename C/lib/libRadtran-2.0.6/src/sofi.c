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
#include <string.h>
#include <math.h>
#include <getopt.h>

#include "ascii.h"
#include "numeric.h"
#include "sofi.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define PROGRAM "sofi"
#define VERS "0.3"

/*************************************************************/
/* print usage information                                   */
/*************************************************************/

static void print_usage (char* filename) {
  fprintf (stderr, "%s %s - Calculate irradiance or diffuse flux during a total eclipse. \n", PROGRAM, VERS);
  fprintf (stderr, "Per default it is assumed that the centre of the eclipse is at the \n");
  fprintf (stderr, "observer position. \n\n");
  fprintf (stderr, "written by Claudia Emde and Paul Ockenfuss,\n");
  fprintf (stderr, "           LMU, e-Mail claudia.emde@.de\n");
  fprintf (stderr, "Version %s finished May 30, 2006\n\n", "0.1");
  fprintf (stderr, "Version %s finished December 3, 2017\n\n", "0.2");
  fprintf (stderr, "Version %s finished October 8, 2018\n\n", "0.3");
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "-h            display help message\n");
  fprintf (stderr, "-n            ID number\n");
  fprintf (stderr, "-z <degrees>  solar zenith angle\n");
  fprintf (stderr, "-a <degrees>  solar azimuth angle \n");
  fprintf (stderr, "-v <degrees>  lunar zenith angle \n");
  fprintf (stderr, "-w <degrees>  lunar azimuth angle \n");
  fprintf (stderr, "-d <km>       distance Sun-Earth \n");
  fprintf (stderr, "-e <km>       distance Moon-Earth \n");
  fprintf (stderr, "-l <nm>       wavelength \n");
  /* fprintf (stderr, "-x <km>       horizontal shift of eclipse shadow \n"); */
  /* fprintf (stderr, "-y <km>       vertical shift of eclipse shadow \n"); */
  fprintf (stderr, "-b <km>       height of atmosphere \n");
  fprintf (stderr, "-s            DeltaX ??? \n");
  fprintf (stderr, "-r            apparent ratio moon/sun \n");
  fprintf (stderr, "-q            Sun radius \n");
  fprintf (stderr, "-i            limb darkening function \n                1-Pierce  2-Neckel  3-Koepke \n");
  fprintf (stderr, "-c            calculate clearsky \n");
  fprintf (stderr, "-p            enable polarization \n");
  fprintf (stderr, "USAGE: %s [options] <filename>\n\n", filename);
  fprintf (stderr, "%s calculates the radiance/irradiance  \n", PROGRAM);
  fprintf (stderr, "during a total solar eclipse.\n");
  fprintf (stderr, "The input file is the photon distribution\n");
  fprintf (stderr, "at the top of the atmosphere calculated using  \n");
  fprintf (stderr, "mystic (mc_backward_writeallpixels), e.g. mc35.rad. \n");
  fprintf (stderr, "The tool can also be used to calculte clearsky \n");
  fprintf (stderr, "radiance/irradiance provided the mystic file above. \n");
  fprintf (stderr, "\n");
}

static int get_options (int     argc,
                        char**  argv,
                        char*   programname,
                        char*   infilename,
                        int*    idNumber,
                        double* sza,
                        double* saa,
                        double* mza,
                        double* maa,
                        double* sdist,
                        double* mdist,
                        double* lambda,
                        double* z,
                        double* DeltaX,
                        double* ratio,
                        double* sunrad,
                        int*    limb,
                        int*    clearsky,
                        int*    pol) {
  int   c     = 0;
  char* dummy = NULL;

  /* save name of program */
  strncpy (programname, argv[0], FILENAME_MAX - 1);

  /* set defaults */
  strcpy (infilename, "");
  *idNumber = 0;
  *sza      = 0.0;
  *saa      = 0.0;
  *maa      = 0.0;
  *mza      = 0.0;
  *lambda   = 400.0;
  *z        = 50.0;
  *ratio    = 1.05;
  *DeltaX   = 6.e-6;
  *limb     = 2;

  while ((c = getopt (argc, argv, "n:z:a:v:w:d:e:l:b:s:r:q:i:cp")) != EOF) {
    switch (c) {
    case 'n':
      *idNumber = strtol (optarg, &dummy, 0);
      break;
    case 'z':
      *sza = strtod (optarg, &dummy);
      break;
    case 'a':
      *saa = strtod (optarg, &dummy);
      break;
    case 'v':
      *mza = strtod (optarg, &dummy);
      break;
    case 'w':
      *maa = strtod (optarg, &dummy);
      break;
    case 'd':
      *sdist = strtod (optarg, &dummy);
      break;
    case 'e':
      *mdist = strtod (optarg, &dummy);
      break;
    case 'l':
      *lambda = strtod (optarg, &dummy);
      break;
    case 'b':
      *z = strtod (optarg, &dummy);
      break;
    case 's':
      *DeltaX = strtod (optarg, &dummy);
      break;
    case 'r':
      *ratio = strtod (optarg, &dummy);
      break;
    case 'q':
      *sunrad = strtod (optarg, &dummy);
      break;
    case 'i':
      *limb = strtol (optarg, &dummy, 0);
      break;
    case 'c':
      *clearsky = 1;
      break;
    case 'p':
      *pol = 1;
      break;
    default:
      print_usage (programname);
      return (-1);
    }
  }

  /* check number of remaining command line arguments */
  if (argc - optind != 1) {
    print_usage (programname);
    return -1;
  }

  strncpy (infilename, argv[optind + 0], FILENAME_MAX - 1);

  return 0; /* if o.k. */
}

int main (int argc, char** argv) {
  int    i = 0, j = 0, il = 0, ip = 0, status = 0;
  int    idNumber = 0;
  int    limb     = 0;
  char   programname[FILENAME_MAX];
  char   infilename[FILENAME_MAX];
  double sza    = 0.0;
  double saa    = 0.0;
  double maa    = 0.0;
  double mza    = 0.0;
  double sdist  = 0.0;
  double mdist  = 0.0;
  double sunrad = 0.0;
  double lambda = 0.0;
  double z_toa  = 50.0;
  /* double dx_moon=0.0; */
  /* double dy_moon=0.0; */
  double DeltaX   = 0.0;
  double ratioObs = 0.0;
  int    clearsky = 0;
  int    pol      = 0;
  double r_earth  = 6371211.0;
  /* double LengthX=2000000.0; FIXCE: not needed?? Side length of the sampling grid */
  /* double LengthY=2000000.0; */
  double* irr;
  double* irr_sum;
  double  r = 0.0, t = 0.0, alpha = 0.0, x = 0.0, y = 0.0;
  int     i_re = 0;
  /* double x2=0.0, y2=0.0;   */
  /* int dx=0.0, dy=0.0; */
  double cossaa = 0.0, sinsaa = 0.0, cossza = 0.0, tansza = 0.0, sinsza = 0.0;
  int    Nx   = 4000; /*Number of points, for which distr. func. at TOA is
		 calculated. DeltaX*Nx must be bigger than ratio+1*/
  FILE * file = NULL, *file2 = NULL, *file3 = NULL;
  int    nstokes = 1;
  double sum_Q2  = 0.;
  /* Variables to read in MYSTIC file (mc_backward_writeallpixels) */
  double *fi = NULL, *se = NULL, *th = NULL, *fo = NULL, *fv = NULL, *si = NULL, *sv = NULL, *ei = NULL;
  int     L = 0;

  double S[3]; /*direction vector pointing towards sun*/
  double M[3]; /*direction vector pointing towards moon*/
  double P[3];
  double Ns[3];
  double Nm[3];
  double dM = 0;

  /* lunar radius will be approximated at these distances from
     observer in SP in direction of the moon's azimuth */
  int    distances[] = {0, -100000, 100000, -200000, 200000, -500000, 500000, -1000000, 1000000};
  int    Ndistances  = 9;
  double ratios[Ndistances];

  /* pd - weighting function */
  double** pd     = calloc (Ndistances, sizeof (double*));
  double   pd_dir = 0.0, gamma = 0.0, cosgamma = 0.0, X = 0.0, ratio = 0.0;
  int      indexRatio = 0, f = 0;
  double   nearestratio = 0.0;

  /* Specify additional files to be written */
  int write_weighted_contribution_function = 0;
  int write_weighting_function             = 0;

  /* Solar eclipse weighting function , see Eq. 4, Emde and Mayer 2007 */
  if ((file = fopen ("sofi_weighting_function.dat", "w")) == NULL)
    return ASCIIFILE_NOT_FOUND;

  /* Result file, includes (polarized) radiance/irradiance during
     eclipse */
  if ((file2 = fopen ("sofi_results.dat", "w")) == NULL)
    return ASCIIFILE_NOT_FOUND;

  status = get_options (argc,
                        argv,
                        programname,
                        infilename,
                        &idNumber,
                        &sza,
                        &saa,
                        &mza,
                        &maa,
                        &sdist,
                        &mdist,
                        &lambda,
                        &z_toa,
                        &DeltaX,
                        &ratioObs,
                        &sunrad,
                        &limb,
                        &clearsky,
                        &pol);

  if (status != 0)
    return status;

  /* enable polarization */
  if (pol == 1)
    nstokes = 4;

  /* Define sun and moon direction vectors */
  cossaa = cos (saa * M_PI / 180.0);
  sinsaa = sin (saa * M_PI / 180.0);
  tansza = tan (sza * M_PI / 180.0);
  cossza = cos (sza * M_PI / 180.0);
  sinsza = sin (sza * M_PI / 180.0);

  maa    = maa * M_PI / 180.0;
  mza    = mza * M_PI / 180.0;
  sunrad = sunrad * M_PI / 180.0;
  z_toa *= 1000;
  sdist *= 1000;
  mdist *= 1000;

  S[0] = sdist * sinsaa * sinsza;
  S[1] = sdist * cossaa * sinsza;
  S[2] = sdist * cossza;
  M[0] = mdist * sin (maa) * sin (mza);
  M[1] = mdist * cos (maa) * sin (mza);
  M[2] = mdist * cos (mza);

  irr_sum = calloc (nstokes, sizeof (double));

  status = read_8c_file (infilename, &fi, &se, &th, &fo, &fv, &si, &sv, &ei, &L);
  fprintf (stderr, " ... read %d data points from %s\n", L, infilename);

  if (status != 0) {
    fprintf (stderr, "error %d reading %s\n", status, infilename);
    return status;
  }

  if (sza >= 90.0) {
    fprintf (stderr, "Solar eclipse calculations not implemented for \n");
    fprintf (stderr, "solar zenith angles greater or equal 90 degrees. \n");
    return -1;
  }

  /* Calculate probability density function. */

  for (int i = 0; i < Ndistances; i++) {
    P[0]      = z_toa * tansza * sinsaa + distances[i] * sinsaa;
    P[1]      = z_toa * tansza * cossaa + distances[i] * cossaa;
    P[2]      = z_toa;
    Ns[0]     = S[0] - P[0];
    Ns[1]     = S[1] - P[1];
    Ns[2]     = S[2] - P[2];
    Nm[0]     = M[0] - P[0];
    Nm[1]     = M[1] - P[1];
    Nm[2]     = M[2] - P[2];
    dM        = sqrt (Nm[0] * Nm[0] + Nm[1] * Nm[1] + Nm[2] * Nm[2]);
    ratios[i] = asin (1737400 / dM) / sunrad;

    pd[i]  = calloc (Nx, sizeof (double));
    status = sample_photons_sofi (lambda, DeltaX, ratios[i], sdist, Nx, limb, pd[i]);
    if (status != 0)
      return status;
    fprintf (stdout, " %.6f", ratios[i]);
  }
  fprintf (stdout, "\n");

  if (write_weighting_function) {
    fprintf (file, "#X and pd(X,ratio) for different ratios\n #X ");
    for (j = 0; j < Ndistances; j++) {
      fprintf (file, " %.6f", ratios[j]);
    }
    fprintf (file, "\n");
    for (i = 0; i < Nx; i++) {
      fprintf (file, "%g", i * DeltaX);
      for (j = 0; j < Ndistances; j++) {
        fprintf (file, " %g", pd[j][i]);
      }
      fprintf (file, "\n");
    }
  }

  /* old method, shift shadow over contribution function given dx_moon and dy_moom,
     replaced by more accurate method */
  /*height of the model atmosphere, convert to meter*/
  /* z_toa*=1000; */
  /* dx_moon*=1000; */
  /* dy_moon*=1000;  */
  /* /\* calculate centre of umbral shadow at TOA *\/ */
  /* r=z_toa*tansza;  */
  /* dx=r*sinsaa;  */
  /* dy=r*cossaa; */

  /* Sampled irradiance */
  irr     = calloc (L, sizeof (double));
  irr_sum = calloc (nstokes, sizeof (double));

  /* Darkning of direct radiance */
  /* calculate sunspot at TOA */
  P[0]  = z_toa * tansza * sinsaa;
  P[1]  = z_toa * tansza * cossaa;
  P[2]  = z_toa;
  Ns[0] = S[0] - P[0];
  Ns[1] = S[1] - P[1];
  Ns[2] = S[2] - P[2];
  Nm[0] = M[0] - P[0];
  Nm[1] = M[1] - P[1];
  Nm[2] = M[2] - P[2];
  dM    = sqrt (Nm[0] * Nm[0] + Nm[1] * Nm[1] + Nm[2] * Nm[2]);

  /* Detailed documentation of following equations in Paul
     Ockenfuss' Bachelor thesis */
  cosgamma = (Ns[0] * Nm[0] + Ns[1] * Nm[1] + Ns[2] * Nm[2]) / (sqrt (Ns[0] * Ns[0] + Ns[1] * Ns[1] + Ns[2] * Ns[2]) * dM);
  if (cosgamma < 1)
    gamma = acos (cosgamma);
  X            = gamma / sunrad;
  ratio        = asin (1737400 / dM) / sunrad;
  i_re         = (int) floor (X / DeltaX);
  nearestratio = 100000.;

  for (f = 0; f < Ndistances; f++) {
    if (fabs (ratio - ratios[f]) < nearestratio) {
      indexRatio   = f;
      nearestratio = fabs (ratio - ratios[f]);
    }
  }
  if (i_re > Nx - 1)
    i_re = Nx - 1;

  /* direct radiation during solar eclipse event */
  pd_dir = pd[indexRatio][i_re];
  fprintf (stdout, "Gamma at observer: %.5e\n", gamma);
  fprintf (stdout, "Calculated ratio at sunspot in sampling plane: %.8e \n", ratio);

  if (!clearsky) {
    for (il = 0; il < L / nstokes; il++)
      for (ip = 0; ip < nstokes; ip++) {
        i = nstokes * il + ip;

        /* Coordinates on spherical sampling area */
        /* CE ??? LengthX and LengthY not needed? */
        x = fi[i] - fi[0] - (fi[L - 1] + fi[0]) / 2;
        y = se[i] - se[0] - (se[L - 1] + se[0]) / 2;

        /* Calculate shift from spherical sampling domain to plane*/
        r     = z_toa + r_earth;
        alpha = asin (sqrt (x * x + y * y) / r);
        t     = r * (1.0 - cos (alpha)) / cossza;

        /* coordinates on plane */
        x += t * sinsza * sinsaa;
        y += t * sinsza * cossaa;

        /* following lines are replaced by more accurate method by Paul */
        /* Shift shadow to intersection point sun - TOA, 
	     needs to be done in plane, since dx, dy are calculated in plane */
        /* x+=(dx+dx_moon); */
        /* y+=(dy+dy_moon); */

        /* /\* Rotation and dilation in y direction *\/ */
        /* x2=cossaa*x-sinsaa*y; */
        /* y2=(sinsaa*x+cossaa*y)*cossza; */

        /* i_re= (int)fabs((sqrt(x2*x2+y2*y2))/1000.0); */

        P[0]     = -x;
        P[1]     = -y;
        P[2]     = z_toa;
        Ns[0]    = S[0] - P[0];
        Ns[1]    = S[1] - P[1];
        Ns[2]    = S[2] - P[2];
        Nm[0]    = M[0] - P[0];
        Nm[1]    = M[1] - P[1];
        Nm[2]    = M[2] - P[2];
        dM       = sqrt (Nm[0] * Nm[0] + Nm[1] * Nm[1] + Nm[2] * Nm[2]);
        cosgamma = (Ns[0] * Nm[0] + Ns[1] * Nm[1] + Ns[2] * Nm[2]) / (sqrt (Ns[0] * Ns[0] + Ns[1] * Ns[1] + Ns[2] * Ns[2]) * dM);
        if (cosgamma < 1)
          gamma = acos (cosgamma);
        else {
          gamma = 0.;
        }
        X            = gamma / sunrad;
        ratio        = asin (1737400 / dM) / sunrad;
        indexRatio   = 0;
        nearestratio = 100000.;
        for (f = 0; f < Ndistances; f++) {
          if (fabs (ratio - ratios[f]) < nearestratio) {
            indexRatio   = f;
            nearestratio = fabs (ratio - ratios[f]);
          }
        }

        i_re = (int) floor (X / DeltaX);

        /* Out of range of weighting function, set i_re to last value (1) */
        if (i_re > Nx - 1)
          i_re = Nx - 1;

        irr[i] = ei[i] * pd[indexRatio][i_re];

        irr_sum[ip] += irr[i];

        /* write weighted contribution function, see Fig. 8, Emde
	     and Mayer 2007. */
        if (write_weighted_contribution_function) {

          if ((file3 = fopen ("sofi_weighted_contribution.dat", "w")) == NULL)
            return ASCIIFILE_NOT_FOUND;

          fprintf (file, "%f  %f %g \n ", fi[i], se[i], irr[i]);

          /*the following line would print the shadow*/
          /*fprintf(stdout,  "%f  %f %g \n " ,fi[i], se[i], pd[i_re]);*/
        }
      }
  } else

    for (il = 0; il < L / nstokes; il++)
      for (ip = 0; ip < nstokes; ip++) {
        i = nstokes * il + ip;
        irr_sum[ip] += ei[i];
      }

  if (pol == 0) {
    fprintf (stdout, "Total radiance/flux I: %.5e \n", irr_sum[0]);
    fprintf (stdout, "Direct irradiance:     %.5e \n", pd_dir);
  }

  else {

    for (ip = 1; ip < nstokes; ip++)
      sum_Q2 += irr_sum[ip] * irr_sum[ip];

    fprintf (stdout, "Total radiance/flux (I,Q,U,V): %.4e %.4e %.4e %.4e\n", irr_sum[0], irr_sum[1], irr_sum[2], irr_sum[3]);
    fprintf (stdout, "Degree of polarization [per cent]: %.4e  \n ", sqrt (sum_Q2) / irr_sum[0] * 100);
    fprintf (stdout, "Direct irradiance:     %.5e \n", pd_dir);

    fprintf (file2,
             "%12.4e %12.4e %12.4e %12.4e %10.2f \n",
             irr_sum[0],
             irr_sum[1],
             irr_sum[2],
             irr_sum[3],
             sqrt (sum_Q2) / irr_sum[0] * 100);
  }

  fclose (file);
  fclose (file2);
  return 0;
}
