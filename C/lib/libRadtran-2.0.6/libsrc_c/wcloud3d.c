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
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <float.h>
#include "ascii.h"
#include "errors.h"
#include "wcloud3d.h"
#include "allocnd.h"

#if HAVE_NETCDF4
#include <netcdf.h>
#define RUN_NC(command, msg)                                                                                                       \
  do {                                                                                                                             \
    int status = (command);                                                                                                        \
    if (status != NC_NOERR) {                                                                                                      \
      fprintf (stderr, "netCDF error %d: \"%s\" %s from %s\n", status, nc_strerror (status), (msg), filename);                     \
      return status;                                                                                                               \
    }                                                                                                                              \
  } while (0)
#endif

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

/* Maximum number of characters per line */
#define MAX_LENGTH_OF_LINE 65536

/** \brief compute readout slice from atmos_region
 *
 * \note all pixel numbers are inclusive
 *
 * \param ixmin [in] smallest pixel index to take in x direction
 * \param ixmax [in] largest pixel index to take in x direction
 * \param iymin [in] smallest pixel index to take in y direction
 * \param iymax [in] largest pixel index to take in y direction
 * \param filename [in] name of input file (for error message only)
 * \param Nx [in/out] number of pixels in x direction (input: total number, output: used number)
 * \param Ny [in/out] number of pixels in y direction (input: total number, output: used number)
 * \param ixstart [out] smallest pixel of selected region in x direction
 * \param iystart [out] smallest pixel of selected region in y direction
 *
 * \return 0 on sucess, -1 otherwise
 */
int get_atmos_slice (int         ixmin,
                     int         ixmax,
                     int         iymin,
                     int         iymax,
                     const char* filename,
                     int*        Nx,
                     int*        Ny,
                     int*        ixstart,
                     int*        iystart) {
  if (ixmin >= 0 && ixmax >= 0 && iymin >= 0 && iymax >= 0) {
    if (ixmax >= *Nx || iymax >= *Ny) {
      fprintf (stderr,
               "Error! You have specified a region with `atmos_region` (%d,%d,%d,%d) which is larger than the file %s\n",
               ixmin,
               ixmax,
               iymin,
               iymax,
               filename);
      return -1;
    }
    *Nx      = ixmax - ixmin + 1;
    *Ny      = iymax - iymin + 1;
    *ixstart = ixmin;
    *iystart = iymin;
  } else {
    ixstart = 0;
    iystart = 0;
  }
  return 0;
}

double kilometer_to_meter (double value) {
  return value * 1000.;
}

void check_threed (float*** data, size_t Nz, size_t Nxy, float invalid_value, int** threed) {
  for (size_t i = 0; i < Nz; i++) {
    if ((*threed)[i])
      continue;
    float* layer_data = *(data[i]);
    for (size_t j = 0; j < Nxy; ++j) {
      if (layer_data[j] != invalid_value) {
        (*threed)[i] = 1;
        break;
      }
    }
  }
}

void transpose_nc_data (float* src, size_t Nz, size_t Nx, size_t Ny, float*** dst) {
  for (size_t l = 0, j = 0; j < Ny; ++j)
    for (size_t i = 0; i < Nx; ++i)
      for (size_t k = 0; k < Nz; ++k, ++l)
        dst[k][i][j] = src[l];
}

/*******************************************************/
/* Read data from a 3D cloud description file.         */
/*******************************************************/
int read_3D_caoth_header (char*   filename,
                          int*    Nx,
                          int*    Ny,
                          int*    Nz,
                          int*    cldproperties,
                          double* delX,
                          double* delY,
                          float** z,
                          int*    rows) {
  int   i = 0, number = 0, status = 0, min_columns = 0, max_columns = 0, max_length = 0;
  char  line[MAX_LENGTH_OF_LINE + 1] = "";
  char* string                       = NULL;

  FILE*  file  = NULL;
  char** array = NULL;
  char*  dummy = NULL;
  /* check input file */
  status = ASCII_checkfile (filename, rows, &min_columns, &max_columns, &max_length);
  if (status != 0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  if (*rows < 2) {
    fprintf (stderr, "Error: found less than two rows in %s\n", filename);
    return -1;
  }

  if (min_columns < 3) {
    fprintf (stderr, "Error: found less than three columns in %s\n", filename);
    return -1;
  }

  /* reset string to the beginning of the memory block */
  string = line;

  /* open file */
  if ((file = fopen (filename, "r")) == NULL)
    return ASCIIFILE_NOT_FOUND;

  /* first line */
  number = 0;
  while (number == 0) {
    dummy  = fgets (string, MAX_LENGTH_OF_LINE, file);
    status = ASCII_parsestring (string, &array, &number);

    if (status != 0) {
      fprintf (stderr, "Error %d reading 1st line of %s\n", status, filename);
      return status;
    }
  }
  if (number < 3) {
    fprintf (stderr, "Error: expected at least 3 values in 1st line, found only %d\n", number);
    return -1;
  }

  /* first 3 values are number of cells */
  *Nx = strtol (array[0], &dummy, 0);
  *Ny = strtol (array[1], NULL, 0);
  *Nz = strtol (array[2], NULL, 0);

  if (number > 3)
    *cldproperties = strtol (array[3], NULL, 0);

  if (number > 4) {
    fprintf (stderr, "*** WARNING: Found a 5th number in the first line of %s.\n", filename);
    fprintf (stderr, "*** Possibly you use the old format. Please be aware that wspec\n");
    fprintf (stderr, "*** doesn't exist anymore and the 5th number is ignored!\n");
  }

  free (array);

  /* 2nd line */
  number = 0;
  while (number == 0) {
    dummy  = fgets (string, MAX_LENGTH_OF_LINE, file);
    status = ASCII_parsestring (string, &array, &number);

    if (status != 0) {
      fprintf (stderr, "Error %d reading 2n line of %s\n", status, filename);
      return status;
    }
  }
  if (number != *Nz + 1 + 2) {
    fprintf (stderr, "Error: found %d z levels, expected %d\n", number - 2, *Nz + 1);
    return -1;
  }
  /* horizontal grid sizes, in km */
  *delX = strtod (array[0], NULL);
  *delY = strtod (array[1], NULL);
  /* rest of line 2 */
  /* set altitude levels */
  if (!(*z = calloc_float_1D (*Nz + 1, "z")))
    return -1;
  for (i = 2; i <= *Nz + 2; i++)
    (*z)[i - 2] = (float) strtod (array[i], NULL);

  free (array);

  fclose (file);
  return 0;
}

int read_3D_caoth_data (const char* filename,
                        const int   Nx,
                        const int   Ny,
                        const int   Nz,
                        const int   rows,
                        float****   data1,
                        float****   data2,
                        float****   data3,
                        int*        indx,
                        int*        indy,
                        int*        indz) {
  FILE* file = NULL;
  int   i = 0, ix = 0, iy = 0, iz = 0, status = 0, number = 0;
  char  line[MAX_LENGTH_OF_LINE + 1] = "";
  char* string                       = NULL;
  string                             = line;

  /* count how many data columns are specified */
  int ndata = data1 != NULL;
  ndata += data2 != NULL;
  ndata += data3 != NULL;

  /* open file */
  if ((file = fopen (filename, "r")) == NULL)
    CHKERROUT (ASCIIFILE_NOT_FOUND, "ASCIIFILE_NOT_FOUND");

  for (i = 0; i < rows; i++) {
    number       = 0;
    char** array = NULL;
    while (number == 0) {
      char* dummy = fgets (string, MAX_LENGTH_OF_LINE, file);
      if (dummy == NULL) {
        fprintf (stderr, "Error reading file %s\n", filename);
        CHKERROUT (-1, "could not determine lenght of file");
      }
      status = ASCII_parsestring (string, &array, &number);
      if (status != 0) {
        fprintf (stderr, "Error %d reading line %d of %s\n", status, i, filename);
        CHKERR (status);
      }
    }
    if (i < 2) { /* jump over the first two non-comment lines (header)*/
      free (array);
      continue;
    }
    if (number < ndata + 3) {
      fprintf (stderr, "Error: found %d columns in row %d, expected %d\n", number, i, ndata + 3);
      CHKERR (-1);
    }
    /* cell indices */
    ix = strtol (array[0], NULL, 0) - 1;
    iy = strtol (array[1], NULL, 0) - 1;
    iz = strtol (array[2], NULL, 0) - 1;
    if (ix < 0 || ix >= Nx) {
      fprintf (stderr, "Error, ix = %d out of bounds in line %d\n", ix + 1, i + 1);
      CHKERROUT (-1, "ix out of bounds");
    }

    if (iy < 0 || iy >= Ny) {
      fprintf (stderr, "Error, iy = %d out of bounds in line %d\n", iy + 1, i + 1);
      CHKERROUT (-1, "iy out of bounds");
    }

    if (iz < 0 || iz >= Nz) {
      fprintf (stderr, "Error, iz = %d out of bounds in line %d\n", iz + 1, i + 1);
      CHKERROUT (-1, "iz out of bounds");
    }
    if (data1) {
      (*data1)[iz][ix][iy] = strtod (array[3], NULL);
      if (data2) {
        (*data2)[iz][ix][iy] = strtod (array[4], NULL);
        if (data3) {
          (*data3)[iz][ix][iy] = strtod (array[5], NULL);
        }
      }
    }
    if (indx) {
      indx[i - 2] = ix;
    }
    if (indy) {
      indy[i - 2] = iy;
    }
    if (indz) {
      indz[i - 2] = iz;
    }
    free (array);
  }
  fclose (file);

  return 0;
}

int read_3D_caoth (char*     filename,
                   int*      Nx,
                   int*      Ny,
                   int*      Nz,
                   double*   delX,
                   double*   delY,
                   float**   z,
                   float**** lwc,
                   float**** reff,
                   float**** cf,
                   float**** ext,
                   float**** g1,
                   float**** g2,
                   float**** ff,
                   float**** f,
                   float**** ssa,
                   float**** dscale,
                   float*    rmin,
                   float*    rmax,
                   int*      cldproperties,
                   int**     threed,
                   int       ixmin,
                   int       ixmax,
                   int       iymin,
                   int       iymax,
                   int       quiet) {
  int status = 0;

  int rows = 0;  // min_columns=0, max_columns=0, max_length=0;

#if HAVE_NETCDF4
  int ixstart = 0, iystart = 0;

  size_t n           = 0;
  int    read_netcdf = 0;

  int get_ext = 0, get_g1 = 0, get_ssa = 0, get_lwc = 0, get_cf = 0, get_reff = 0;

  int ncid = 0, idd_nx = 0, idd_ny = 0, idd_nz = 0;
  int id_z = 0, id_ext = 0, id_lwc = 0, id_cf = 0, id_reff = 0, id_g1 = 0, id_ssa = 0;
#endif

#if HAVE_NETCDF4
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status == NC_NOERR) {
    read_netcdf = 1;

    if (!quiet)
      fprintf (stderr, " ... reading Cloud data from netCDF file %s\n", filename);

    RUN_NC (nc_inq_dimid (ncid, "nx", &idd_nx), "reading nx");
    RUN_NC (nc_inq_dimlen (ncid, idd_nx, &n), "reading nx");
    *Nx = n;

    RUN_NC (nc_inq_dimid (ncid, "ny", &idd_ny), "reading ny");
    RUN_NC (nc_inq_dimlen (ncid, idd_ny, &n), "reading ny");
    *Ny = n;

    RUN_NC (nc_inq_dimid (ncid, "nz", &idd_nz), "reading nz");
    RUN_NC (nc_inq_dimlen (ncid, idd_nz, &n), "reading nz");
    *Nz = n;

    *z = calloc_float_1D (*Nz + 1, "z");
    CHKPOINTER (*z);

    RUN_NC (nc_get_att_int (ncid, NC_GLOBAL, "cldproperties", cldproperties), "reading cldproperties");
    RUN_NC (nc_get_att_double (ncid, NC_GLOBAL, "dx", delX), "reading dx");
    RUN_NC (nc_get_att_double (ncid, NC_GLOBAL, "dy", delY), "reading dy");

    if (get_atmos_slice (ixmin, ixmax, iymin, iymax, filename, Nx, Ny, &ixstart, &iystart) != 0)
      CHKERR (-1);

  } else {
    if (!quiet)
      fprintf (stderr, " ... %s not in netCDF format, trying to open as ASCII\n", filename);
#endif

    status = read_3D_caoth_header (filename, Nx, Ny, Nz, cldproperties, delX, delY, z, &rows);
    CHKERR (status);

#if HAVE_NETCDF4
  }
#endif

  if (!quiet) {
    fprintf (stderr, " ... reading %d x %d x %d data points from %s\n", *Nx, *Ny, *Nz, filename);
    fprintf (stderr, " ... cldproperties = %d\n", *cldproperties);
    fprintf (stderr, " ... cloud grid size in x direction: %g km\n", *delX);
    fprintf (stderr, " ... cloud grid size in y direction: %g km\n", *delY);
  }

  /* allocate memory for altitude levels */
  *threed = calloc_int_1D (*Nz, "threed");
  CHKPOINTER (*threed);

  /* allocate memory for optical depth, asymmetry factor, and single scattering albedo */
  *ext = calloc_float_3D (*Nz, *Nx, *Ny, "ext");
  CHKPOINTER (*ext);
  *g1 = calloc_float_3D (*Nz, *Nx, *Ny, "g1");
  CHKPOINTER (*g1);
  *g2 = calloc_float_3D (*Nz, *Nx, *Ny, "g2");
  CHKPOINTER (*g2);
  *ff = calloc_float_3D (*Nz, *Nx, *Ny, "ff");
  CHKPOINTER (*ff);
  *ssa = calloc_float_3D (*Nz, *Nx, *Ny, "ssa");
  CHKPOINTER (*ssa);
  *f = calloc_float_3D (*Nz, *Nx, *Ny, "f");
  CHKPOINTER (*f);
  *dscale = calloc_float_3D (*Nz, *Nx, *Ny, "dscale");
  CHKPOINTER (*dscale);
  *lwc = calloc_float_3D (*Nz, *Nx, *Ny, "lwc");
  CHKPOINTER (*lwc);
  *cf = calloc_float_3D (*Nz, *Nx, *Ny, "cf");
  CHKPOINTER (*cf);
  *reff = calloc_float_3D (*Nz, *Nx, *Ny, "reff");
  CHKPOINTER (*reff);

  /* initialize forward HG fraction with 1 */
  {
    size_t N = (*Nz) * (*Nx) * (*Ny);
    for (size_t i = 0; i < N; ++i) {
      (***ff)[i] = 1.0;
    }
  }

  /* initialize rmin, rmax */
  switch (*cldproperties) {
  case CLD_OPTPROP:
    *rmin = NAN;
    *rmax = NAN; /* rmin, rmax not used */
    break;

  case CLD_EXTREFF:
  case CLD_LWCREFF:
  case CLD_LWCREFFCF:
    *rmin = +FLT_MAX;
    *rmax = -FLT_MAX;
    break;

  default:
    fprintf (stderr, "Error, unknown cloud properties %d in %s\n", *cldproperties, filename);
    CHKERROUT (-1, "Error, unknown cloud properties");
  }

#if HAVE_NETCDF4
  if (read_netcdf == 1) {

    /* read data */

    /* set which parameters to read */
    switch (*cldproperties) {
    case CLD_OPTPROP:
      get_ext = 1;
      get_g1  = 1;
      get_ssa = 1;
      break;
    case CLD_EXTREFF:
      get_ext  = 1;
      get_reff = 1;
      break;
    case CLD_LWCREFF:
      get_lwc  = 1;
      get_reff = 1;
      break;
    case CLD_LWCREFFCF:
      get_lwc  = 1;
      get_reff = 1;
      get_cf   = 1;
      break;
    default:
      fprintf (stderr, "Error, unknown cloud properties %d in %s\n", *cldproperties, filename);
      CHKERROUT (-1, "Error, unknown cloud properties");
    }

    size_t zstart = 0;
    size_t nz_lev = *Nz + 1;
    RUN_NC (nc_inq_varid (ncid, "z", &id_z), "reading z");
    RUN_NC (nc_get_vara_float (ncid, id_z, &zstart, &nz_lev, *z), "reading z");

    int    use_direct_nc_read = ixstart != 0 || iystart != 0 || ixmax + 1 != *Nx || iymax + 1 != *Ny;
    float* tmp_data           = NULL;
    if (use_direct_nc_read) {
      if (!quiet) {
        fprintf (stderr, "Using direct, consecutive netcdf read mode, not the mapped version because it is way faster\n");
      }
      tmp_data = malloc (sizeof (float) * (*Nz) * (*Nx) * (*Ny));
    }

    size_t    N          = (*Nz) * (*Nx) * (*Ny);
    size_t    tstart[3]  = {iystart, ixstart, 0};
    size_t    tcount[3]  = {*Ny, *Nx, *Nz};
    ptrdiff_t tstride[3] = {1, 1, 1};
    ptrdiff_t timap[3]   = {1, *Ny, (*Ny) * (*Nx)};
    if (get_ext) {
      RUN_NC (nc_inq_varid (ncid, "ext", &id_ext), "reading ext");
      if (use_direct_nc_read) {
        RUN_NC (nc_get_var_float (ncid, id_ext, tmp_data), "reading ext");
        transpose_nc_data (tmp_data, *Nz, *Nx, *Ny, *ext);
      } else {
        RUN_NC (nc_get_varm_float (ncid, id_ext, tstart, tcount, tstride, timap, ***ext), "reading ext");
      }
      for (size_t i = 0; i < N; ++i) {
        (***ext)[i] /= 1000.0; /* convert from km-1 to m-1 */
      }
      check_threed (*ext, *Nz, (*Nx) * (*Ny), 0.0, threed);
    }

    if (get_lwc) {
      RUN_NC (nc_inq_varid (ncid, "lwc", &id_lwc), "reading lwc");
      if (use_direct_nc_read) {
        RUN_NC (nc_get_var_float (ncid, id_lwc, tmp_data), "reading lwc");
        transpose_nc_data (tmp_data, *Nz, *Nx, *Ny, *lwc);
      } else {
        RUN_NC (nc_get_varm_float (ncid, id_lwc, tstart, tcount, tstride, timap, ***lwc), "reading lwc");
      }
      check_threed (*lwc, *Nz, (*Nx) * (*Ny), 0.0, threed);
    }

    if (get_cf) {
      RUN_NC (nc_inq_varid (ncid, "cf", &id_cf), "reading cloud fraction");
      if (use_direct_nc_read) {
        RUN_NC (nc_get_var_float (ncid, id_cf, tmp_data), "reading cloud fraction");
        transpose_nc_data (tmp_data, *Nz, *Nx, *Ny, *cf);
      } else {
        RUN_NC (nc_get_varm_float (ncid, id_cf, tstart, tcount, tstride, timap, ***cf), "reading cloud fraction");
      }
      check_threed (*cf, *Nz, (*Nx) * (*Ny), 0.0, threed);
    }

    if (get_reff) {
      RUN_NC (nc_inq_varid (ncid, "reff", &id_reff), "reading reff");
      if (use_direct_nc_read) {
        RUN_NC (nc_get_var_float (ncid, id_reff, tmp_data), "reading reff");
        transpose_nc_data (tmp_data, *Nz, *Nx, *Ny, *reff);
      } else {
        RUN_NC (nc_get_varm_float (ncid, id_reff, tstart, tcount, tstride, timap, ***reff), "reading reff");
      }
      for (size_t i = 0; i < N; ++i) {
        if (get_lwc) {
          if ((***lwc)[i] > 0) {
            *rmin = MIN (*rmin, (***reff)[i]);
            *rmax = MAX (*rmax, (***reff)[i]);
          }
        }
      }
    }

    if (get_g1) {
      RUN_NC (nc_inq_varid (ncid, "g1", &id_g1), "reading g1");
      if (use_direct_nc_read) {
        RUN_NC (nc_get_var_float (ncid, id_g1, tmp_data), "reading g1");
        transpose_nc_data (tmp_data, *Nz, *Nx, *Ny, *g1);
      } else {
        RUN_NC (nc_get_varm_float (ncid, id_g1, tstart, tcount, tstride, timap, ***g1), "reading g1");
      }
    }

    if (get_ssa) {
      RUN_NC (nc_inq_varid (ncid, "ssa", &id_ssa), "reading ssa");
      if (use_direct_nc_read) {
        RUN_NC (nc_get_var_float (ncid, id_ssa, tmp_data), "reading ssa");
        transpose_nc_data (tmp_data, *Nz, *Nx, *Ny, *ssa);
      } else {
        RUN_NC (nc_get_varm_float (ncid, id_ssa, tstart, tcount, tstride, timap, ***ssa), "reading ssa");
      }
    }

    nc_close (ncid);

    if (use_direct_nc_read)
      free (tmp_data);

    if (!quiet)
      fprintf (stderr, "\n");
  } else {
#endif

    {
      int* indx = calloc (rows - 2, sizeof (int));
      int* indy = calloc (rows - 2, sizeof (int));
      int* indz = calloc (rows - 2, sizeof (int));
      /* copy data to result arrays */
      switch (*cldproperties) {
      case CLD_OPTPROP:
        status = read_3D_caoth_data (filename, *Nx, *Ny, *Nz, rows, ext, g1, ssa, indx, indy, indz);
        CHKERROUT (status, "bad CLD_OPTPROP file");
        for (size_t i = 0; i < rows - 2; i++) {
          (*ext)[indz[i]][indx[i]][indy[i]] /= 1000.0; /* convert from km-1 to m-1 */
          (*ff)[indz[i]][indx[i]][indy[i]] = 1.0;
        }

        break;

      case CLD_EXTREFF:
        status = read_3D_caoth_data (filename, *Nx, *Ny, *Nz, rows, ext, reff, NULL, indx, indy, indz);
        CHKERROUT (status, "bad CLD_EXTREFF file");
        for (size_t j = 0; j < rows - 2; j++) {
          (*ext)[indz[j]][indx[j]][indy[j]] /= 1000.0; /* convert from km-1 to m-1 */
          if ((*ext)[indz[j]][indx[j]][indy[j]] > 0) {
            *rmin = MIN (*rmin, (*reff)[indz[j]][indx[j]][indy[j]]);
            *rmax = MAX (*rmax, (*reff)[indz[j]][indx[j]][indy[j]]);
          }
        }
        break;

      case CLD_LWCREFF:
        status = read_3D_caoth_data (filename, *Nx, *Ny, *Nz, rows, lwc, reff, NULL, indx, indy, indz);
        CHKERROUT (status, "bad CLD_LWCREFF file");
        for (size_t j = 0; j < rows - 2; j++) {
          if ((*lwc)[indz[j]][indx[j]][indy[j]] > 0) {
            *rmin = MIN (*rmin, (*reff)[indz[j]][indx[j]][indy[j]]);
            *rmax = MAX (*rmax, (*reff)[indz[j]][indx[j]][indy[j]]);
          }
        }
        break;

      case CLD_LWCREFFCF:
        status = read_3D_caoth_data (filename, *Nx, *Ny, *Nz, rows, lwc, reff, cf, indx, indy, indz);
        CHKERROUT (status, "bad CLD_LWCREFF_CF file");
        for (size_t j = 0; j < rows - 2; j++) {
          if ((*lwc)[indz[j]][indx[j]][indy[j]] > 0) {
            *rmin = MIN (*rmin, (*reff)[indz[j]][indx[j]][indy[j]]);
            *rmax = MAX (*rmax, (*reff)[indz[j]][indx[j]][indy[j]]);
          }
        }
        break;

      default:
        fprintf (stderr, "Error, unknown cloud properties %d in %s\n", *cldproperties, filename);
        status = -1;
        CHKERROUT (status, "Error, unknown cloud properties");
      }

      /* set flag indicating that this layer is 3D */
      for (size_t j = 0; j < rows - 2; j++) {
        (*threed)[indz[j]] = 1;
      }
      UVSPEC_FREE (indx);
      UVSPEC_FREE (indy);
      UVSPEC_FREE (indz);
    }

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", rows - 2, filename);

#if HAVE_NETCDF4
  }
#endif

  // convert units
  *delX = kilometer_to_meter (*delX);
  *delY = kilometer_to_meter (*delY);

  // check for valid effective radius read ins
  const float reff_limits[2] = {0, FLT_MAX};
  if (*rmin <= reff_limits[0] || *rmax >= reff_limits[1]) {
    fprintf (stderr,
             "\nCaution, found invalid effective radius!\n"
             " Expect reff to be in range ] %g, %g [\n but found reff_min/max [ %g, %g ]\n",
             reff_limits[0],
             reff_limits[1],
             *rmin,
             *rmax);
    status = -1;
    if (*rmin <= reff_limits[0])
      CHKERROUT (status, "Error, Found too small effective radius");
    if (*rmin >= reff_limits[0])
      CHKERROUT (status, "Error, Found too large effective radius");
  }

  return 0;
}

/*******************************************************/
/* Read data from a 2D albedo description file.        */
/*******************************************************/

int read_2D_albedo (char*     filename,
                    int*      Nx,
                    int*      Ny,
                    double*   delX,
                    double*   delY,
                    double*** albedo,
                    int       ixmin,
                    int       ixmax,
                    int       iymin,
                    int       iymax,
                    int       quiet) {
  int i = 0, status = 0;
  int ix = 0, iy = 0;

  double** value = NULL;
  int      rows = 0, min_columns = 0, max_columns = 0;

#if HAVE_NETCDF4
  int ixstart = 0, iystart = 0;

  struct stat buf;
  size_t      n = 0;

  int ncid = 0, idd_nx = 0, idd_ny = 0;
  int id_type = 0;
#endif

  /* try to open as netCDF file */

#if HAVE_NETCDF4
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status == NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading 2D albedo data from netCDF file %s\n", filename);

    /* determine file date and stop if the file was created before */
    /* September 7, 2007 where the format was changed              */
    stat (filename, &buf);

    /* CP: Commented out since it's not needed */
    /*tmstruct=gmtime(&buf.st_mtime);*/

    if (buf.st_mtime < 1189187082) {
      fprintf (stderr, "\n");
      fprintf (stderr, "*** Error %s was last changed %s", filename, ctime (&buf.st_mtime));
      fprintf (stderr, "*** and the MYSTIC albedo2D convention has been changed on\n");
      fprintf (stderr, "*** September 7, 2007. It is likely that you use the \n");
      fprintf (stderr, "*** old convention! Please check the documentation!\n");
      fprintf (stderr, "\n");

      return -1;
    }

    fprintf (stderr, "\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!! Careful! The netcdf file format for 2D albedo was !!!\n");
    fprintf (stderr, "!!! changed on September 7, 2007. In particular,      !!!\n");
    fprintf (stderr, "!!! x and y were interchanged so that an albedo file  !!!\n");
    fprintf (stderr, "!!! viewed with ncview should look like a map now,    !!!\n");
    fprintf (stderr, "!!! with (x,y)=(1,1) in the lower left corner.        !!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "\n");

    RUN_NC (nc_inq_dimid (ncid, "nx", &idd_nx), "reading nx");
    RUN_NC (nc_inq_dimlen (ncid, idd_nx, &n), "reading nx");
    *Nx = n;

    RUN_NC (nc_inq_dimid (ncid, "ny", &idd_ny), "reading ny");
    RUN_NC (nc_inq_dimlen (ncid, idd_ny, &n), "reading ny");
    *Ny = n;

    RUN_NC (nc_get_att_double (ncid, NC_GLOBAL, "dx", delX), "reading dx");
    RUN_NC (nc_get_att_double (ncid, NC_GLOBAL, "dy", delY), "reading dy");

    if (get_atmos_slice (ixmin, ixmax, iymin, iymax, filename, Nx, Ny, &ixstart, &iystart) != 0)
      return -1;

    /* allocate memory for 2D albedo field */
    if (!(*albedo = calloc_double_2D (*Nx, *Ny, "albedo")))
      return -1;

    /* read data */
    size_t    tstart[2]  = {iystart, ixstart};
    size_t    tcount[2]  = {*Ny, *Nx};
    ptrdiff_t tstride[2] = {1, 1};
    ptrdiff_t timap[2]   = {1, *Ny};

    RUN_NC (nc_inq_varid (ncid, "albedo", &id_type), "reading albedo");
    RUN_NC (nc_get_varm_double (ncid, id_type, tstart, tcount, tstride, timap, **albedo), "reading albedo");
    nc_close (ncid);

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", (*Nx) * (*Ny), filename);
  } else {
    if (!quiet)
      fprintf (stderr, " ... %s not in netCDF format, trying to open as ASCII\n", filename);
#endif

    /* read file; ASCII_file2double is a waste of memory in this */
    /* case, but it is very convenient because comments and      */
    /* everything are handled correctly.                         */

    status = ASCII_file2double (filename, &rows, &max_columns, &min_columns, &value);

    if (status != 0) {
      fprintf (stderr, "Error %d reading file %s\n", status, filename);
      return status;
    }

    if (max_columns < 4) {
      fprintf (stderr, "Error, found less than four columns in %s\n", filename);
      return -1;
    }

    if (rows < 1) {
      if (!quiet)
        fprintf (stderr, "Error, no header in %s\n", filename);
      return -1;
    }

    /* 1st line, number of cells */
    *Nx = (int) (value[0][0] + 0.5);
    *Ny = (int) (value[0][1] + 0.5);

    *delX = value[0][2];
    *delY = value[0][3];

    if (!quiet) {
      fprintf (stderr, " ... x-distance: %g km\n", *delX);
      fprintf (stderr, " ... y-distance: %g km\n", *delY);
    }

    /* allocate memory for 2D albedo field */
    if (!(*albedo = calloc_double_2D (*Nx, *Ny, "albedo")))
      return -1;

    for (i = 1; i < rows; i++) {
      ix = (int) (value[i][0] + 0.5) - 1;
      iy = (int) (value[i][1] + 0.5) - 1;

      if (ix < 0 || ix >= *Nx) {
        fprintf (stderr, "Error, ix = %d out of bounds in line %d\n", ix + 1, i + 1);
        free_double_2D (*albedo);
        *albedo = NULL;
        return -1;
      }

      if (iy < 0 || iy >= *Ny) {
        fprintf (stderr, "Error, iy = %d out of bounds in line %d\n", iy + 1, i + 1);
        free_double_2D (*albedo);
        *albedo = NULL;
        return -1;
      }

      /* copy data to result arrays */
      (*albedo)[ix][iy] = value[i][2];
    }

    /* free memory */
    (void) ASCII_free_double (value, rows);

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", rows - 1, filename);

#if HAVE_NETCDF4
  }
#endif

  // convert units
  *delX = kilometer_to_meter (*delX);
  *delY = kilometer_to_meter (*delY);

  return 0;
}

/*******************************************************/
/* Read data from a 2D umu/phi description file.        */
/*******************************************************/

int read_2D_umu (char*     filename,
                 int*      Nx,
                 int*      Ny,
                 double*   delX,
                 double*   delY,
                 double*** umu,
                 double*** phi,
                 int       ixmin,
                 int       ixmax,
                 int       iymin,
                 int       iymax,
                 int       quiet) {
  int i = 0, status = 0;
  int ix = 0, iy = 0;

  double** value = NULL;
  int      rows = 0, min_columns = 0, max_columns = 0;

#if HAVE_NETCDF4
  int ixstart = 0, iystart = 0;

  size_t n = 0;

  int ncid = 0, idd_nx = 0, idd_ny = 0;
  int id_type = 0;
#endif

  /* try to open as netCDF file */

#if HAVE_NETCDF4
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status == NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading 2D umu data from netCDF file %s\n", filename);

    RUN_NC (nc_inq_dimid (ncid, "nx", &idd_nx), "reading nx");
    RUN_NC (nc_inq_dimlen (ncid, idd_nx, &n), "reading nx");
    *Nx = n;

    RUN_NC (nc_inq_dimid (ncid, "ny", &idd_ny), "reading ny");
    RUN_NC (nc_inq_dimlen (ncid, idd_ny, &n), "reading ny");
    *Ny = n;

    RUN_NC (nc_get_att_double (ncid, NC_GLOBAL, "dx", delX), "reading dx");
    RUN_NC (nc_get_att_double (ncid, NC_GLOBAL, "dy", delY), "reading dy");

    if (get_atmos_slice (ixmin, ixmax, iymin, iymax, filename, Nx, Ny, &ixstart, &iystart) != 0)
      return -1;

    /* allocate memory for 2D umu/phi fields */
    if (!(*umu = calloc_double_2D (*Nx, *Ny, "umu")))
      return -1;
    if (!(*phi = calloc_double_2D (*Nx, *Ny, "phi"))) {
      free_double_2D (*umu);
      return -1;
    }

    /* read data */
    size_t    tstart[2]  = {iystart, ixstart};
    size_t    tcount[2]  = {*Ny, *Nx};
    ptrdiff_t tstride[2] = {1, 1};
    ptrdiff_t timap[2]   = {1, *Ny};

    RUN_NC (nc_inq_varid (ncid, "umu", &id_type), "reading umu");
    RUN_NC (nc_get_varm_double (ncid, id_type, tstart, tcount, tstride, timap, **umu), "reading umu");

    RUN_NC (nc_inq_varid (ncid, "phi", &id_type), "reading phi");
    RUN_NC (nc_get_varm_double (ncid, id_type, tstart, tcount, tstride, timap, **phi), "reading phi");
    for (size_t i = 0; i < (*Nx) * (*Ny); ++i) {
      /* Switch from DISORT to MYSTIC convention */
      (**phi)[i] += 180.0;
    }

    nc_close (ncid);

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", (*Nx) * (*Ny), filename);
  } else {
    if (!quiet)
      fprintf (stderr, " ... %s not in netCDF format, trying to open as ASCII\n", filename);
#endif

    /* read file; ASCII_file2double is a waste of memory in this */
    /* case, but it is very convenient because comments and      */
    /* everything are handled correctly.                         */

    status = ASCII_file2double (filename, &rows, &max_columns, &min_columns, &value);

    if (status != 0) {
      fprintf (stderr, "Error %d reading file %s\n", status, filename);
      return status;
    }

    if (max_columns < 4) {
      fprintf (stderr, "Error, found less than four columns in %s\n", filename);
      return -1;
    }

    if (rows < 1) {
      if (!quiet)
        fprintf (stderr, "Error, no header in %s\n", filename);
      return -1;
    }

    /* 1st line, number of cells */
    *Nx = (int) (value[0][0] + 0.5);
    *Ny = (int) (value[0][1] + 0.5);

    /* 1st line, horizontal distances */
    *delX = value[0][2];
    *delY = value[0][3];

    if (!quiet) {
      fprintf (stderr, " ... x-distance: %g km\n", *delX);
      fprintf (stderr, " ... y-distance: %g km\n", *delY);
    }

    /* allocate memory for 2D umu/phi fields */
    if (!(*umu = calloc_double_2D (*Nx, *Ny, "umu")))
      return -1;
    if (!(*phi = calloc_double_2D (*Nx, *Ny, "phi"))) {
      free_double_2D (*umu);
      *umu = NULL;
      return -1;
    }

    for (i = 1; i < rows; i++) {
      ix = (int) (value[i][0] + 0.5) - 1;
      iy = (int) (value[i][1] + 0.5) - 1;

      if (ix < 0 || ix >= *Nx) {
        fprintf (stderr, "Error, ix = %d out of bounds in line %d\n", ix + 1, i + 1);
        free_double_2D (*umu);
        *umu = NULL;
        free_double_2D (*phi);
        *phi = NULL;
        return -1;
      }

      if (iy < 0 || iy >= *Ny) {
        fprintf (stderr, "Error, iy = %d out of bounds in line %d\n", iy + 1, i + 1);
        free_double_2D (*umu);
        *umu = NULL;
        free_double_2D (*phi);
        *phi = NULL;
        return -1;
      }

      /* copy data to result arrays */
      (*umu)[ix][iy] = value[i][2];
      /* Switch from DISORT to MYSTIC convention */
      (*phi)[ix][iy] = value[i][3] + 180.0;
    }

    /* free memory */
    (void) ASCII_free_double (value, rows);

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", rows - 1, filename);

#if HAVE_NETCDF4
  }
#endif

  // convert units
  *delX = kilometer_to_meter (*delX);
  *delY = kilometer_to_meter (*delY);

  return 0;
}

/******************************************************************/
/* Read data from a 2D surface (albedo or RPV) description file.  */
/******************************************************************/

int read_2D_surface_labels (char*            filename,
                            char**           label_lib,
                            int              nlabel_lib,
                            int*             Nx,
                            int*             Ny,
                            double*          delX,
                            double*          delY,
                            unsigned char*** label,
                            int**            nlabel,
                            int*             therewereCaMs,
                            int              quiet) {
  int i = 0, il = 0, status = 0;
  int ix = 0, iy = 0;

  int rows = 0, min_columns = 0, max_columns = 0, max_length = 0;

  char*** string = NULL;

#if HAVE_NETCDF4
  size_t tstart[2] = {0, 0};
  size_t tcount[2] = {0, 0};

  size_t n = 0;

  int ncid   = 0;
  int idd_nx = 0, idd_ny = 0;
  int id_type = 0;

  unsigned char* temp       = NULL;
  char           tempstr[4] = "";
#endif

  *therewereCaMs = 0;

  /* try open as netCDF file */

#if HAVE_NETCDF4
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status == NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading surface data from netCDF file %s\n", filename);

    RUN_NC (nc_inq_dimid (ncid, "nx", &idd_nx), "reading nx");
    RUN_NC (nc_inq_dimlen (ncid, idd_nx, &n), "reading nx");
    *Nx = n;

    RUN_NC (nc_inq_dimid (ncid, "ny", &idd_ny), "reading ny");
    RUN_NC (nc_inq_dimlen (ncid, idd_ny, &n), "reading ny");
    *Ny = n;

    RUN_NC (nc_get_att_double (ncid, NC_GLOBAL, "dx", delX), "reading dx");
    RUN_NC (nc_get_att_double (ncid, NC_GLOBAL, "dy", delY), "reading dy");

    /* allocate memory for 2D BRDF field */
    if (!(*label = calloc_uchar_2D (*Nx, *Ny, "label")))
      return -1;
    if (!(*nlabel = calloc_int_1D (nlabel_lib, "nlabel"))) {
      free_uchar_2D (*label);
      *label = NULL;
      return -1;
    }

    /* read data */

    RUN_NC (nc_inq_varid (ncid, "type", &id_type), "reading type");

    temp = calloc (*Nx, sizeof (unsigned char));

    /* read type */
    tstart[1] = 0;   /* start with first x element    */
    tcount[1] = *Nx; /* read *Nx elements for each iy */

    tcount[0] = 1; /* read one y element at a time  */

    for (iy = 0; iy < *Ny; iy++) {
      tstart[0] = iy;

      if (!quiet)
        fprintf (stderr, ".");

      RUN_NC (nc_get_vara_uchar (ncid, id_type, tstart, tcount, temp), "reading type");

      for (ix = 0; ix < *Nx; ix++) {

        sprintf (tempstr, "%d", temp[ix]);

        for (il = 0; il < nlabel_lib; il++) {
          if (!strcmp (tempstr, label_lib[il]))
            break;
        }

        if (il == nlabel_lib) {
          fprintf (stderr, "Error, did not find an entry for %s\n", tempstr);
          return -1;
        }

        if (il > UCHAR_MAX) {
          fprintf (stderr, "Error, index %d larger than %d\n", il, UCHAR_MAX);
          return -1;
        }

        (*label)[ix][iy] = il;
        (*nlabel)[il]++;
      }
    }
    nc_close (ncid);

    free (temp);

    if (!quiet)
      fprintf (stderr, "\n");
  } else {
    if (!quiet)
      fprintf (stderr, " ... %s not in netCDF format, trying to open as ASCII\n", filename);
#endif

    /* read file; ASCII_file2double is a waste of memory in this */
    /* case, but it is very convenient because comments and      */
    /* everything are handled correctly.                         */

    status = ASCII_checkfile (filename, &rows, &min_columns, &max_columns, &max_length);

    if (status != 0) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }

    if (max_columns < 4) {
      fprintf (stderr, "Error, found less than four columns in %s\n", filename);
      return -1;
    }

    if (rows < 1) {
      if (!quiet)
        fprintf (stderr, "Error, no header in %s\n", filename);
      return -1;
    }

    /* allocate memory */
    status = ASCII_calloc_string (&string, rows, max_columns, max_length);

    if (status != 0) {
      fprintf (stderr, "Error %d allocating memory\n", status);
      return status;
    }

    /* read file to string array */
    status = ASCII_readfile (filename, string);

    if (status != 0) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }

    /* 1st line, number of cells */
    *Nx = strtol (string[0][0], NULL, 0);
    *Ny = strtol (string[0][1], NULL, 0);

    /* 1st line, horizontal distances */
    *delX = strtod (string[0][2], NULL);
    *delY = strtod (string[0][3], NULL);

    if (!quiet) {
      fprintf (stderr, " ... x-distance: %g km\n", *delX);
      fprintf (stderr, " ... y-distance: %g km\n", *delY);
    }

    /* allocate memory for 2D BRDF field */
    if (!(*label = calloc_uchar_2D (*Nx, *Ny, "label")))
      return -1;
    if (!(*nlabel = calloc_int_1D (nlabel_lib, "nlabel"))) {
      free_uchar_2D (*label);
      *label = NULL;
      return -1;
    }

    for (i = 1; i < rows; i++) {
      ix = strtol (string[i][0], NULL, 0) - 1;
      iy = strtol (string[i][1], NULL, 0) - 1;

      if (ix < 0 || ix >= *Nx) {
        fprintf (stderr, "Error, ix = %d out of bounds in line %d\n", ix + 1, i + 1);
        return -1;
      }

      if (iy < 0 || iy >= *Ny) {
        fprintf (stderr, "Error, iy = %d out of bounds in line %d\n", iy + 1, i + 1);
        return -1;
      }

      /* copy data to result arrays */
      for (il = 0; il < nlabel_lib; il++) {
        if (!strcmp (string[i][2], label_lib[il]))
          break;
      }

      if (il == nlabel_lib) {
        fprintf (stderr, "Error, did not find an entry for %s\n", string[i][2]);
        return -1;
      }

      /* if (il>UCHAR_MAX) { */
      /* 	fprintf (stderr, "Error, index %d larger than %d\n",  */
      /* 		 il, UCHAR_MAX); */
      /* 	return -1; */
      /* } */

      (*label)[ix][iy] = il;
      (*nlabel)[il]++;
    }

    /* free memory */
    (void) ASCII_free_string (string, rows, max_columns);

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", rows - 1, filename);

    if ((*nlabel)[0] > 0)
      *therewereCaMs = 1;

#if HAVE_NETCDF4
  }
#endif

  // convert units
  *delX = kilometer_to_meter (*delX);
  *delY = kilometer_to_meter (*delY);

  return 0;
}

/*******************************************************/
/* Read data from a 2D Ross Li BRDF description file.  */
/*******************************************************/

int read_2D_rossli (const char*         filename,
                    int*                Nx,
                    int*                Ny,
                    double*             delX,
                    double*             delY,
                    rossli_brdf_spec*** rossli,
                    int*                therewereCaMs,
                    int                 hotspot,
                    int                 isAmbralsFile,
                    int                 ixmin,
                    int                 ixmax,
                    int                 iymin,
                    int                 iymax,
                    int                 quiet) {
  int i = 0, status = 0;
  int ix = 0, iy = 0;
  int contains_CaM = 0, isCaM = 0;

#if HAVE_NETCDF4
  size_t tstart[2] = {0, 0};
  size_t tcount[2] = {0, 0};
  int    ixstart = 0, iystart = 0;

  size_t n = 0;

  int ncid = 0, idd_nx = 0, idd_ny = 0;
  int id_iscam = 0;
  int id_iso = 0, id_geo = 0, id_vol = 0;

  int*   isCaMv = NULL;
  float* temp   = NULL;
#endif

  double** value = NULL;
  int      rows = 0, min_columns = 0, max_columns = 0;

  double vol_fac = 1.0, iso_fac = 1.0, geo_fac = 1.0;

  *therewereCaMs = 0;

  if (isAmbralsFile) {
    vol_fac = 3. / 4.;
    iso_fac = 1. / PI;
    geo_fac = 1. / PI;
  }

  /* try open as netCDF file */

#if HAVE_NETCDF4
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status == NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading Ross-Li data from netCDF file %s\n", filename);

    RUN_NC (nc_inq_dimid (ncid, "nx", &idd_nx), "reading nx");
    RUN_NC (nc_inq_dimlen (ncid, idd_nx, &n), "reading nx");
    *Nx = n;

    RUN_NC (nc_inq_dimid (ncid, "ny", &idd_ny), "reading ny");
    RUN_NC (nc_inq_dimlen (ncid, idd_ny, &n), "reading ny");
    *Ny = n;

    RUN_NC (nc_get_att_double (ncid, NC_GLOBAL, "dx", delX), "reading dx");
    RUN_NC (nc_get_att_double (ncid, NC_GLOBAL, "dy", delY), "reading dy");

    if (get_atmos_slice (ixmin, ixmax, iymin, iymax, filename, Nx, Ny, &ixstart, &iystart) != 0)
      return -1;

    /* allocate memory for albedo */
    *rossli = calloc ((size_t) *Nx, sizeof (rossli_brdf_spec*));

    for (ix = 0; ix < *Nx; ix++)
      (*rossli)[ix] = calloc ((size_t) *Ny, sizeof (rossli_brdf_spec));

    /* read data */

    RUN_NC (nc_inq_varid (ncid, "iso", &id_iso), "reading iso");
    RUN_NC (nc_inq_varid (ncid, "geo", &id_geo), "reading geo");
    RUN_NC (nc_inq_varid (ncid, "vol", &id_vol), "reading vol");
    RUN_NC (nc_inq_varid (ncid, "isCaM", &id_iscam), "reading isCaM");

    temp   = calloc (*Nx, sizeof (float));
    isCaMv = calloc (*Nx, sizeof (int));

    /* read type */
    tstart[1] = ixstart; /* start with first x element    */
    tcount[1] = *Nx;     /* read *Nx elements for each iy */

    tcount[0] = 1; /* read one y element at a time  */

    for (iy = 0; iy < *Ny; iy++) {
      tstart[0] = iy + iystart;

      if (!quiet)
        fprintf (stderr, ".");

      RUN_NC (nc_get_vara_int (ncid, id_iscam, tstart, tcount, isCaMv), "reading isCaM");

      for (ix = 0; ix < *Nx; ix++)
        if (isCaMv[ix] == 1) {
          (*rossli)[ix][iy].isCaM = 1;
          *therewereCaMs          = 1;
        } else
          (*rossli)[ix][iy].isCaM = 0;

      RUN_NC (nc_get_vara_float (ncid, id_iso, tstart, tcount, temp), "reading iso");

      for (ix = 0; ix < *Nx; ix++)
        if (isCaMv[ix] == 1)
          (*rossli)[ix][iy].iso = 0.0;
        else
          (*rossli)[ix][iy].iso = iso_fac * temp[ix];

      RUN_NC (nc_get_vara_float (ncid, id_geo, tstart, tcount, temp), "reading geo");

      for (ix = 0; ix < *Nx; ix++)
        if (isCaMv[ix] == 1)
          (*rossli)[ix][iy].geo = 0.0;
        else
          (*rossli)[ix][iy].geo = geo_fac * temp[ix];

      RUN_NC (nc_get_vara_float (ncid, id_vol, tstart, tcount, temp), "reading vol");

      for (ix = 0; ix < *Nx; ix++)
        if (isCaMv[ix] == 1)
          (*rossli)[ix][iy].vol = 0.0;
        else
          (*rossli)[ix][iy].vol = vol_fac * temp[ix];

      for (ix = 0; ix < *Nx; ix++)
        if (isCaMv[ix] == 1)
          (*rossli)[ix][iy].hotspot = 0;
        else
          (*rossli)[ix][iy].hotspot = hotspot;
    }
    nc_close (ncid);

    free (temp);
    free (isCaMv);

    if (!quiet)
      fprintf (stderr, "\n");
  } else {
    if (!quiet)
      fprintf (stderr, " ... %s not in netCDF format, trying to open as ASCII\n", filename);
#endif

    /* read file; ASCII_file2double is a waste of memory in this */
    /* case, but it is very convenient because comments and      */
    /* everything are handled correctly.                         */

    status = ASCII_file2double (filename, &rows, &max_columns, &min_columns, &value);

    if (status != 0) {
      fprintf (stderr, "Error %d reading file %s\n", status, filename);
      return status;
    }

    if (max_columns < 5) {
      fprintf (stderr, "Error, found less than five columns in %s\n", filename);
      return -1;
    }

    /* if more than 5 columns, the 6th column defines whether CaM or not */
    if (max_columns > 5)
      contains_CaM = 1;

    if (rows < 1) {
      if (!quiet)
        fprintf (stderr, "Error, no header in %s\n", filename);
      return -1;
    }

    /* 1st line, number of cells */
    *Nx = (int) (value[0][0] + 0.5);
    *Ny = (int) (value[0][1] + 0.5);

    /* 1st line, horizontal distances */
    *delX = value[0][2];
    *delY = value[0][3];

    if (!quiet) {
      fprintf (stderr, " ... x-distance: %g km\n", *delX);
      fprintf (stderr, " ... y-distance: %g km\n", *delY);
    }

    /* allocate memory for albedo */
    *rossli = calloc ((size_t) *Nx, sizeof (rossli_brdf_spec*));

    for (ix = 0; ix < *Nx; ix++)
      (*rossli)[ix] = calloc ((size_t) *Ny, sizeof (rossli_brdf_spec));

    for (i = 1; i < rows; i++) {
      ix = (int) (value[i][0] + 0.5) - 1;
      iy = (int) (value[i][1] + 0.5) - 1;

      if (ix < 0 || ix >= *Nx) {
        fprintf (stderr, "Error, ix = %d out of bounds in line %d\n", ix + 1, i + 1);
        return -1;
      }

      if (iy < 0 || iy >= *Ny) {
        fprintf (stderr, "Error, iy = %d out of bounds in line %d\n", iy + 1, i + 1);
        return -1;
      }

      if (contains_CaM)
        isCaM = (value[i][5] == 1);

      if (isCaM) {
        /* CaM */
        (*rossli)[ix][iy].iso     = 0.0;
        (*rossli)[ix][iy].vol     = 0.0;
        (*rossli)[ix][iy].geo     = 0.0;
        (*rossli)[ix][iy].hotspot = 0;
        (*rossli)[ix][iy].isCaM   = 1;
        *therewereCaMs            = 1;
      } else {
        /* copy data to result arrays */
        (*rossli)[ix][iy].iso     = iso_fac * value[i][2];
        (*rossli)[ix][iy].vol     = geo_fac * value[i][3];
        (*rossli)[ix][iy].geo     = vol_fac * value[i][4];
        (*rossli)[ix][iy].hotspot = hotspot;
        (*rossli)[ix][iy].isCaM   = 0;
      }
    }

    /* free memory */
    (void) ASCII_free_double (value, rows);

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", rows - 1, filename);

#if HAVE_NETCDF4
  }
#endif

  // convert units
  *delX = kilometer_to_meter (*delX);
  *delY = kilometer_to_meter (*delY);

  return 0;
}
