#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ascii.h"
#include "cdisort.h"
#include "errors.h"
#include "mystic.h"
#include "nca.h"
#include "solver.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

/* ############################################################################### */
/* ############################## Function - NCA 3D ############################## */
/* ##### This function calculates the radiative transfer via a modified        ### */
/* ######################### schwarzschild solver. ############################### */
/* ###################### using forward mode of mystic ########################### */
/* ############################################################################### */

int nca_like_mystic (/* input */
                     char*              nca_options,
                     char*              datapath,
                     double             dx,
                     double             dy,
                     float*             zprof,
                     atmosphere_struct* atmos,
                     surftemp_struct*   surftemp,
                     float***           temperature,
                     float              wvnmlo,
                     float              wvnmhi,
                     int                source,
                     int                temper3d,
                     /* output */
                     float*** abs3d) {

  //local variables
  int       ix = 0, iy = 0, iz = 0, ilyr = 0;
  double*** Eup_3d   = NULL;
  double*** Edn_3d   = NULL;
  double*** abs_nca  = NULL;
  double*** planck   = NULL;
  double*** kabs_nca = NULL;

  // schwarzschild variables
  double mu      = 0.0;
  double dmu     = 0.01;
  double ds      = 0.0;
  double z_start = 0.0;
  double z1      = 0.0;
  double z2      = 0.0;
  double B1 = 0.0, B2 = 0.0;
  double alpha = 0.0;
  double beta  = 0.0;

  double*** L_up_schw = NULL;
  double*** L_dn_schw = NULL;

  double* L_up = calloc (atmos->Nz + 1, sizeof (double));
  double* L_dn = calloc (atmos->Nz + 1, sizeof (double));

  Eup_3d = calloc (atmos->Nz + 1, sizeof (double**));
  for (iz = 0; iz < atmos->Nz + 1; iz++) {
    Eup_3d[iz] = calloc (atmos->Nx, sizeof (double*));
    for (ix = 0; ix < atmos->Nx; ix++)
      Eup_3d[iz][ix] = calloc (atmos->Ny, sizeof (double));
  }
  Edn_3d = calloc (atmos->Nz + 1, sizeof (double**));
  for (iz = 0; iz < atmos->Nz + 1; iz++) {
    Edn_3d[iz] = calloc (atmos->Nx, sizeof (double*));
    for (ix = 0; ix < atmos->Nx; ix++)
      Edn_3d[iz][ix] = calloc (atmos->Ny, sizeof (double));
  }

  L_dn_schw = calloc (atmos->Nx, sizeof (double**));
  for (ix = 0; ix < atmos->Nx; ix++) {
    L_dn_schw[ix] = calloc (atmos->Ny, sizeof (double*));
    for (iy = 0; iy < atmos->Ny; iy++)
      L_dn_schw[ix][iy] = calloc (atmos->Nz + 1, sizeof (double));
  }

  L_up_schw = calloc (atmos->Nx, sizeof (double**));
  for (ix = 0; ix < atmos->Nx; ix++) {
    L_up_schw[ix] = calloc (atmos->Ny, sizeof (double*));
    for (iy = 0; iy < atmos->Ny; iy++)
      L_up_schw[ix][iy] = calloc (atmos->Nz + 1, sizeof (double));
  }

  abs_nca = calloc (atmos->Nz, sizeof (double**));
  for (iz = 0; iz < atmos->Nz; iz++) {
    abs_nca[iz] = calloc (atmos->Nx, sizeof (double*));
    for (ix = 0; ix < atmos->Nx; ix++)
      abs_nca[iz][ix] = calloc (atmos->Ny, sizeof (double));
  }

  planck = calloc (atmos->Nz + 1, sizeof (double**));
  for (iz = 0; iz < atmos->Nz + 1; iz++) {
    planck[iz] = calloc (atmos->Nx, sizeof (double*));
    for (ix = 0; ix < atmos->Nx; ix++)
      planck[iz][ix] = calloc (atmos->Ny, sizeof (double));
  }

  kabs_nca = calloc (atmos->Nz, sizeof (double**));
  for (iz = 0; iz < atmos->Nz; iz++) {
    kabs_nca[iz] = calloc (atmos->Nx, sizeof (double*));
    for (ix = 0; ix < atmos->Nx; ix++)
      kabs_nca[iz][ix] = calloc (atmos->Ny, sizeof (double));
  }

  for (ix = 0; ix < atmos->Nx; ix++) {
    for (iy = 0; iy < atmos->Ny; iy++) {

      // Set Planck in atmosphere
      for (ilyr = 0; ilyr < atmos->Nz + 1; ilyr++) {

        if (source == MCSRC_THERMAL_FORWARD) {
          if (temper3d) {
            planck[ilyr][ix][iy] = c_planck_func1 (wvnmlo, wvnmhi, temperature[ix][iy][ilyr]);
          } else {
            planck[ilyr][ix][iy] = c_planck_func1 (wvnmlo, wvnmhi, temperature[0][0][ilyr]);
          }
        } else {
          planck[ilyr][ix][iy] = 0.;
        }
      }
      if (surftemp->surf2D && source == MCSRC_THERMAL_FORWARD) {
        planck[ilyr][ix][iy] = c_planck_func1 (wvnmlo, wvnmhi, surftemp->temp2D[ix][iy]);
      }

      //get optical properties
      for (ilyr = 0; ilyr < atmos->Nz; ilyr++) {
        if (atmos->threed[MCCAOTH_TOT][ilyr] >= 1) {  // 3d atmosphere
          kabs_nca[ilyr][ix][iy] = ((atmos->kabs3D->prof[MCCAOTH_TOT])[ilyr][ix][iy]) * 1000.;
        } else {
          kabs_nca[ilyr][ix][iy] = ((atmos->kabs->prof[MCCAOTH_TOT])[ilyr]) * 1000.;
        }
      }
    }
  }

  for (ix = 0; ix < atmos->Nx; ix++) {   /* loop over all x gridboxes */
    for (iy = 0; iy < atmos->Ny; iy++) { /* loop over all y gridboxes */

      /* ## Integrate from mu=0 to mu=1 in steps of dmu ## */
      for (mu = dmu / 2; mu <= 1; mu += dmu) {

        /* ## Boundary conditions at TOA ## */
        L_dn[atmos->Nz] = 0;
        L_dn_schw[ix][iy][atmos->Nz] += L_dn[atmos->Nz] * mu;

        /* ## Calculate downwelling radiation ## */
        for (ilyr = atmos->Nz - 1; ilyr >= 0; ilyr--) {

          z1      = atmos->Z[ilyr] / 1000.;
          z2      = atmos->Z[ilyr + 1] / 1000.;
          ds      = fabs (z2 - z1) / fabs (mu);
          z_start = z1;

          B1 = planck[ilyr + 1][ix][iy];
          B2 = planck[ilyr][ix][iy];

          alpha = B2 + ((B1 - B2) / fabs (z1 - z2) *
                        fabs (z_start - z1)); /* ## in this case z_start and z2 are equal and therefore alpha = B2 ## */
          beta  = (B1 - B2) / fabs (z1 - z2) * fabs (mu);

          if (
            kabs_nca[ilyr][ix][iy] * ds <
            1e-8) /* ## possible numerical instability for very small kabs! This is solved by Taylor approximating of the last two terms ## */
            L_dn[ilyr] = L_dn[ilyr + 1] * exp (-kabs_nca[ilyr][ix][iy] * ds) + kabs_nca[ilyr][ix][iy] * ds * (alpha + beta * ds);
          else
            L_dn[ilyr] = L_dn[ilyr + 1] * exp (-kabs_nca[ilyr][ix][iy] * ds) +
                         (beta / kabs_nca[ilyr][ix][iy] + alpha) * (1 - exp (-kabs_nca[ilyr][ix][iy] * ds)) -
                         beta * ds * exp (-kabs_nca[ilyr][ix][iy] * ds);
          L_dn_schw[ix][iy][ilyr] += L_dn[ilyr] * mu * dmu * 2;
        }

        /* ## Boundary conditions at Ground ## */
        L_up[0] = planck[0][ix][iy];
        L_up_schw[ix][iy][0] += L_up[0] * mu * dmu * 2;

        /* ## Calculate upwelling radiation ## */
        for (ilyr = 1; ilyr < atmos->Nz + 1; ilyr++) {
          z1      = atmos->Z[ilyr] / 1000.;
          z2      = atmos->Z[ilyr - 1] / 1000.;
          ds      = fabs (z2 - z1) / fabs (mu);
          z_start = z2;

          B1 = planck[ilyr - 1][ix][iy];
          B2 = planck[ilyr][ix][iy];

          alpha = B2 + ((B1 - B2) / fabs (z1 - z2) *
                        fabs (z_start - z2)); /* ## in this case z_start and z2 are equal and therefore alpha = B2 ## */
          beta  = (B1 - B2) / fabs (z1 - z2) * fabs (mu);

          if (
            kabs_nca[ilyr - 1][ix][iy] * ds <
            1e-8) /* ## possible numerical instability for very small kabs! This is solved by Taylor approximating of the last two terms ## */
            L_up[ilyr] =
              L_up[ilyr - 1] * exp (-kabs_nca[ilyr - 1][ix][iy] * ds) + kabs_nca[ilyr - 1][ix][iy] * ds * (alpha + beta * ds);
          else
            L_up[ilyr] = L_up[ilyr - 1] * exp (-kabs_nca[ilyr - 1][ix][iy] * ds) +
                         (beta / kabs_nca[ilyr - 1][ix][iy] + alpha) * (1 - exp (-kabs_nca[ilyr - 1][ix][iy] * ds)) -
                         beta * ds * exp (-kabs_nca[ilyr - 1][ix][iy] * ds);

          L_up_schw[ix][iy][ilyr] += L_up[ilyr] * mu * dmu * 2;

        } /* end ilyr */
      }   /* end imu */

      for (ilyr = 0; ilyr < atmos->Nz + 1; ilyr++) {
        Eup_3d[ilyr][ix][iy] = L_up_schw[ix][iy][ilyr];
        Edn_3d[ilyr][ix][iy] = L_dn_schw[ix][iy][ilyr];
      }

    } /* ix */
  }   /* iy */

  // bring dx, dy and kabs to km grid!
  dx = dx / 1000.;
  dy = dy / 1000.;

  if (strcmp (nca_options, "nca2016") == 0) {
    int status = nca3d_1 (datapath, dx, dy, atmos, Edn_3d, Eup_3d, planck, kabs_nca, abs_nca);
    if (status != 0) {
      fprintf (stderr, "Error %d returned by nca_like_mystic with nca_option version1 ()\n", status);
      return status;
    }
  } else if (strcmp (nca_options, "nca2019_cuboid") == 0) {
    int status = nca3d_2 (datapath, dx, dy, atmos, Edn_3d, Eup_3d, planck, kabs_nca, abs_nca);
    if (status != 0) {
      fprintf (stderr, "Error %d returned by nca_like_mystic with nca_option version2_cuboid ()\n", status);
      return status;
    }
  } else if (strcmp (nca_options, "nca2019_triangle") == 0) {
    int status = nca3d_3 (datapath, dx, dy, atmos, Edn_3d, Eup_3d, planck, kabs_nca, abs_nca);
    if (status != 0) {
      fprintf (stderr, "Error %d returned by nca_like_mystic with nca_option version2_triangle ()\n", status);
      return status;
    }
  } else {
    printf ("nca_options <%s> is not valid, choose e.g. mc_nca {nca2016, nca2019_cuboid, nca2019_triangle}\n", nca_options);
    return -2;
  }

  //return nca heating
  for (ix = 0; ix < atmos->Nx; ix++) {   /* loop over all x gridboxes */
    for (iy = 0; iy < atmos->Ny; iy++) { /* loop over all y gridboxes */
      for (ilyr = 0; ilyr < atmos->Nz; ilyr++) {
        if (atmos->threed[MCCAOTH_TOT][ilyr] >= 1) {
          abs3d[ilyr][ix][iy] = abs_nca[ilyr][ix][iy];
        }
      }
    }
  }

  for (ix = 0; ix < atmos->Nx; ix++) {
    for (iy = 0; iy < atmos->Ny; iy++)
      free (L_up_schw[ix][iy]);
    free (L_up_schw[ix]);
  }
  free (L_up_schw);

  for (ix = 0; ix < atmos->Nx; ix++) {
    for (iy = 0; iy < atmos->Ny; iy++)
      free (L_dn_schw[ix][iy]);
    free (L_dn_schw[ix]);
  }
  free (L_dn_schw);

  for (iz = 0; iz < atmos->Nz; iz++) {
    for (ix = 0; ix < atmos->Nx; ix++)
      free (Edn_3d[iz][ix]);
    free (Edn_3d[iz]);
  }
  free (Edn_3d);

  for (iz = 0; iz < atmos->Nz; iz++) {
    for (ix = 0; ix < atmos->Nx; ix++)
      free (Eup_3d[iz][ix]);
    free (Eup_3d[iz]);
  }
  free (Eup_3d);

  for (iz = 0; iz < atmos->Nz; iz++) {
    for (ix = 0; ix < atmos->Nx; ix++)
      free (abs_nca[iz][ix]);
    free (abs_nca[iz]);
  }
  free (abs_nca);

  return 0;
}

/* ############################################################################ */
/* ########### 2016 Version of NCA. Please cite Klinger and Mayer 2016   ###### */
/* ########### when using this method .                                  ###### */
/* ########### NCA calculates thermal heating rates by only using input  ###### */
/* ########### directly neighboring columns, neglect of scattering and   ###### */
/* ########### by using an angle of 45 degree for radiation propagation. ###### */
/* ############################################################################ */

int nca3d_1 (char*              datapath,
             double             dx,
             double             dy,
             atmosphere_struct* atmos,
             double***          Edn_3d,
             double***          Eup_3d,
             double***          planck,
             double***          kabs_nca,
             double***          abs_nca) {

  //local variables

  double* HR_up   = calloc (atmos->Nz + 2, sizeof (double));
  double* HR_dn   = calloc (atmos->Nz + 2, sizeof (double));
  double* HR_up_s = calloc (atmos->Nz + 2, sizeof (double));
  double* HR_dn_s = calloc (atmos->Nz + 2, sizeof (double));

  int ilyr = 0;  // counter z levels
  int ix   = 0;  // counter x grid boxes
  int iy   = 0;  // counter y grid boxes

  double Absup  = 0.0;  // Upwelling, Absorption, lower/upper face
  double Absdn  = 0.0;  // Downwelling Absorption, lower/upper face
  double Absups = 0.0;  // Upwelling, Absorption, face 1
  double Absdns = 0.0;  // Downwelling Absorption, face 1

  double Emup  = 0.0;  // Upwelling, Emission, lower/upper face
  double Emdn  = 0.0;  // Downwelling Emission, lower/upper face
  double Emups = 0.0;  // Upwelling, Emission, face 1
  double Emdns = 0.0;  // Downwelling Emission, face 1

  double z1 = 0.0;  // height level 1
  double z2 = 0.0;  // height level 2

  double mu = 0.0;  // zenith angle

  double B1 = 0.0;  // Planck emission at layer 1
  double B2 = 0.0;  // Planck emission at layer 2

  double ax = 0.0;  // integration boarder for side contributions
  double bx = 0.0;  // integration boarder for side contributions
  double cx = 0.0;  // integration boarder for side contributions
  double ay = 0.0;  // integration boarder for side contributions
  double by = 0.0;  // integration boarder for side contributions
  double cy = 0.0;  // integration boarder for side contributions
  double az = 0.0;  // integration boarder for side contributions
  double bz = 0.0;  // integration boarder for side contributions
  double cz = 0.0;  // integration boarder for side contributions

  double az1 = 0.0;  // integration boarder for side contributions
  double bz1 = 0.0;  // integration boarder for side contributions
  double cz1 = 0.0;  // integration boarder for side contributions
  double z11 = 0.0;  // height level 1
  double z22 = 0.0;  // height level 2

  double l     = 0.0;
  double Trans = 0.0;
  double B1_1  = 0.0;  // Planck emission at layer 1
  double B2_1  = 0.0;  // Planck emission at layer 2

  double afit1   = 0.0;
  double bfit1   = 0.0;
  double bfit2   = 0.0;
  double afit2   = 0.0;
  double cfit2   = 0.0;
  double factor2 = 0.0;
  double factor1 = 0.0;
  double eps     = 0.000001;

  for (ix = 0; ix < atmos->Nx; ix++) {                /* loop over all x gridboxes */
    for (iy = 0; iy < atmos->Ny; iy++) {              /* loop over all y gridboxes */
      for (ilyr = 0; ilyr <= atmos->Nz - 1; ilyr++) { /* loop over all height levels */

        if (atmos->threed[MCCAOTH_TOT][ilyr] >= 1) { /* only use NCA in 3D layers */

          mu = cos (45. * PI / 180.);

          if (mu < 0)
            mu = mu * (-1.);

          z1 = atmos->Z[ilyr] / 1000.;
          z2 = atmos->Z[ilyr + 1] / 1000.;
          ax = 0.0;
          ay = 0.0;
          az = 0.0;
          cx = dx;
          cy = dy;
          cz = fabs (z2 - z1);

          // ## set parameters for integration
          bx = cx - (fabs (z2 - z1) / mu) * sqrt (1. - mu * mu);
          by = cy - (fabs (z2 - z1) / mu) * sqrt (1. - mu * mu);
          bz = cz - (dx / sqrt (1. - mu * mu)) * mu;

          if (fabs (bx - cx) < eps || fabs (by - cy) < eps) {
            bx = ax;
            by = ay;
          }
          if (fabs (bz - cz) < eps)
            bz = az;

          if (bx < 0 || by < 0) {
            bx = ax;
            by = ay;
          }
          if (bz < 0)
            bz = az;

          if (bx < 0.0001 || by < 0.0001) {
            bx = 0.0;
            by = 0.0;
          }
          if (bz < 0.0001)
            bz = 0.0;

          B1 = planck[ilyr][ix][iy];
          B2 = planck[ilyr + 1][ix][iy];

          // #### top/bottom ####
          if (ilyr < atmos->Nz) {

            B1_1 = planck[ilyr + 1][ix][iy];
            B2_1 = planck[ilyr + 2][ix][iy];

            z11 = atmos->Z[ilyr] / 1000.;
            z22 = atmos->Z[ilyr + 1] / 1000.;

            az1 = 0.0;
            cz1 = fabs (z22 - z11);

            // ## set parameters for integration

            bz1 = cz1 - (dx / sqrt (1. - mu * mu)) * mu;

            if (fabs (bz1 - cz1) < eps)
              bz1 = az1;

            if (bz1 < 0)
              bz1 = az1;

            if (bz1 < 0.0001)
              bz1 = 0.0;
            if (atmos->threed[MCCAOTH_TOT][ilyr] >= 1) {  // 3d atmosphere
              Trans = integrate_flux ((B1_1 + B2_1) / 2.,
                                      (B1_1 + B2_1) / 2.,
                                      Edn_3d[ilyr + 2][get_index (ix, atmos->Nx)][iy],
                                      kabs_nca[ilyr + 1][get_index (ix, atmos->Nx)][iy],
                                      kabs_nca[ilyr + 1][ix][iy],
                                      bz1,
                                      cz1,
                                      fabs (z22 - z11),
                                      dy,
                                      mu) +
                      integrate_flux ((B1_1 + B2_1) / 2.,
                                      (B1_1 + B2_1) / 2.,
                                      Edn_3d[ilyr + 2][get_index (ix - 1, atmos->Nx)][iy],
                                      kabs_nca[ilyr + 1][get_index (ix - 1, atmos->Nx)][iy],
                                      kabs_nca[ilyr + 1][ix][iy],
                                      bz1,
                                      cz1,
                                      fabs (z22 - z11),
                                      dy,
                                      mu) +
                      integrate_flux ((B1_1 + B2_1) / 2.,
                                      (B1_1 + B2_1) / 2.,
                                      Edn_3d[ilyr + 2][ix][get_index (iy - 1, atmos->Ny)],
                                      kabs_nca[ilyr + 1][ix][get_index (iy - 1, atmos->Ny)],
                                      kabs_nca[ilyr + 1][ix][iy],
                                      bz1,
                                      cz1,
                                      fabs (z22 - z11),
                                      dx,
                                      mu) +
                      integrate_flux ((B1_1 + B2_1) / 2.,
                                      (B1_1 + B2_1) / 2.,
                                      Edn_3d[ilyr + 2][ix][get_index (iy, atmos->Ny)],
                                      kabs_nca[ilyr + 1][ix][get_index (iy, atmos->Ny)],
                                      kabs_nca[ilyr + 1][ix][iy],
                                      bz1,
                                      cz1,
                                      fabs (z22 - z11),
                                      dx,
                                      mu);

              l = (Trans + Edn_3d[ilyr + 1][ix][iy]) / 5.;
            } else
              l = Edn_3d[ilyr + 1][ix][iy];
          } else
            l = 0;

          Absdn = l * integrate_emis (kabs_nca[ilyr][ix][iy], ax, bx, cx, dx, fabs (z2 - z1), mu);
          Emdn  = -((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], ax, bx, cx, dx, fabs (z2 - z1), mu);
          Absdn = Absdn + l * integrate_emis (kabs_nca[ilyr][ix][iy], ay, by, cy, dy, fabs (z2 - z1), mu);
          Emdn  = Emdn - ((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], ay, by, cy, dy, fabs (z2 - z1), mu);

          // #### go to left ####
          Absdns = integrate_abs ((B1 + B2) / 2.,
                                  Edn_3d[ilyr + 1][get_index (ix - 1, atmos->Nx)][iy],
                                  kabs_nca[ilyr][get_index (ix - 1, atmos->Nx)][iy],
                                  kabs_nca[ilyr][ix][iy],
                                  az,
                                  bz,
                                  cz,
                                  fabs (z2 - z1),
                                  dy,
                                  mu);
          Emdns  = -((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], az, bz, cz, fabs (z2 - z1), dy, mu);

          // #### go to right ####
          Absdns = Absdns + integrate_abs ((B1 + B2) / 2.,
                                           Edn_3d[ilyr + 1][get_index (ix + 1, atmos->Nx)][iy],
                                           kabs_nca[ilyr][get_index (ix + 1, atmos->Nx)][iy],
                                           kabs_nca[ilyr][ix][iy],
                                           az,
                                           bz,
                                           cz,
                                           fabs (z2 - z1),
                                           dy,
                                           mu);
          Emdns  = Emdns - ((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], az, bz, cz, fabs (z2 - z1), dy, mu);

          // #### go to front ####
          Absdns = Absdns + integrate_abs ((B1 + B2) / 2.,
                                           Edn_3d[ilyr + 1][ix][get_index (iy - 1, atmos->Ny)],
                                           kabs_nca[ilyr][ix][get_index (iy - 1, atmos->Ny)],
                                           kabs_nca[ilyr][ix][iy],
                                           az,
                                           bz,
                                           cz,
                                           fabs (z2 - z1),
                                           dx,
                                           mu);
          Emdns  = Emdns - ((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], az, bz, cz, fabs (z2 - z1), dx, mu);

          //##### go to back ####
          Absdns = Absdns + integrate_abs ((B1 + B2) / 2.,
                                           Edn_3d[ilyr + 1][ix][get_index (iy + 1, atmos->Ny)],
                                           kabs_nca[ilyr][ix][get_index (iy + 1, atmos->Ny)],
                                           kabs_nca[ilyr][ix][iy],
                                           az,
                                           bz,
                                           cz,
                                           fabs (z2 - z1),
                                           dx,
                                           mu);
          Emdns  = Emdns - ((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], az, bz, cz, fabs (z2 - z1), dx, mu);

          // ## upwelling
          // bottom face of gridbox
          if (ilyr > 0) {
            B2_1 = planck[ilyr - 1][ix][iy];  //atmos->Bplanck->prof  [MCCAOTH_TOT][ilyr+1][ix][iy];
            B1_1 = planck[ilyr][ix][iy];      //atmos->Bplanck->prof  [MCCAOTH_TOT][ilyr+2][ix][iy];

            z11 = atmos->Z[ilyr - 1] / 1000.;
            z22 = atmos->Z[ilyr] / 1000.;

            az1 = 0.0;
            cz1 = z22 - z11;

            // ## set parameters for integration
            bz1 = cz1 - (dx / sqrt (1. - mu * mu)) * mu;

            if (fabs (bz1 - cz1) < eps)
              bz1 = az1;

            if (bz1 < 0)
              bz1 = az1;

            if (bz1 < 1e-7)
              bz1 = 0.0;

            if (atmos->threed[MCCAOTH_TOT][ilyr] >= 1) {  // 3d atmosphere
              Trans = integrate_flux ((B1_1 + B2_1) / 2.,
                                      (B1_1 + B2_1) / 2.,
                                      Eup_3d[ilyr - 1][get_index (ix + 1, atmos->Nx)][iy],
                                      kabs_nca[ilyr - 1][get_index (ix + 1, atmos->Nx)][iy],
                                      kabs_nca[ilyr - 1][ix][iy],
                                      bz1,
                                      cz1,
                                      fabs (z22 - z11),
                                      dy,
                                      mu) +
                      integrate_flux ((B1_1 + B2_1) / 2.,
                                      (B1_1 + B2_1) / 2.,
                                      Eup_3d[ilyr - 1][get_index (ix - 1, atmos->Nx)][iy],
                                      kabs_nca[ilyr - 1][get_index (ix - 1, atmos->Nx)][iy],
                                      kabs_nca[ilyr - 1][ix][iy],
                                      bz1,
                                      cz1,
                                      fabs (z22 - z11),
                                      dy,
                                      mu) +
                      integrate_flux ((B1_1 + B2_1) / 2.,
                                      (B1_1 + B2_1) / 2.,
                                      Eup_3d[ilyr - 1][ix][get_index (iy + 1, atmos->Ny)],
                                      kabs_nca[ilyr - 1][ix][get_index (iy + 1, atmos->Ny)],
                                      kabs_nca[ilyr - 1][ix][iy],
                                      bz1,
                                      cz1,
                                      fabs (z22 - z11),
                                      dx,
                                      mu) +
                      integrate_flux ((B1_1 + B2_1) / 2.,
                                      (B1_1 + B2_1) / 2.,
                                      Eup_3d[ilyr - 1][ix][get_index (iy - 1, atmos->Ny)],
                                      kabs_nca[ilyr - 1][ix][get_index (iy - 1, atmos->Ny)],
                                      kabs_nca[ilyr - 1][ix][iy],
                                      bz1,
                                      cz1,
                                      fabs (z22 - z11),
                                      dx,
                                      mu);
              l = (Trans + Eup_3d[ilyr][ix][iy]) / 5.;
            } else
              l = Eup_3d[ilyr][ix][iy];
          } else
            l = Eup_3d[0][ix][iy];

          Absup = l * integrate_emis (kabs_nca[ilyr][ix][iy], ax, bx, cx, dx, fabs (z2 - z1), mu);
          Emup  = -((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], ax, bx, cx, dx, fabs (z2 - z1), mu);
          Absup = Absup + l * integrate_emis (kabs_nca[ilyr][ix][iy], ay, by, cy, dy, fabs (z2 - z1), mu);
          Emup  = Emup - ((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], ay, by, cy, dy, fabs (z2 - z1), mu);

          // #### go to left ####
          Absups = integrate_abs ((B1 + B2) / 2.,
                                  Eup_3d[ilyr][get_index (ix - 1, atmos->Nx)][iy],
                                  kabs_nca[ilyr][get_index (ix - 1, atmos->Nx)][iy],
                                  kabs_nca[ilyr][ix][iy],
                                  az,
                                  bz,
                                  cz,
                                  fabs (z2 - z1),
                                  dy,
                                  mu);
          Emups  = -((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], az, bz, cz, fabs (z2 - z1), dy, mu);

          // #### go to right ####
          Absups = Absups + integrate_abs ((B1 + B2) / 2.,
                                           Eup_3d[ilyr][get_index (ix + 1, atmos->Nx)][iy],
                                           kabs_nca[ilyr][get_index (ix + 1, atmos->Nx)][iy],
                                           kabs_nca[ilyr][ix][iy],
                                           az,
                                           bz,
                                           cz,
                                           fabs (z2 - z1),
                                           dy,
                                           mu);
          Emups  = Emups - ((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], az, bz, cz, fabs (z2 - z1), dy, mu);

          // #### go to front ####
          Absups = Absups + integrate_abs ((B1 + B2) / 2.,
                                           Eup_3d[ilyr][ix][get_index (iy - 1, atmos->Ny)],
                                           kabs_nca[ilyr][ix][get_index (iy - 1, atmos->Ny)],
                                           kabs_nca[ilyr][ix][iy],
                                           az,
                                           bz,
                                           cz,
                                           fabs (z2 - z1),
                                           dx,
                                           mu);
          Emups  = Emups - ((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], az, bz, cz, fabs (z2 - z1), dx, mu);

          // ##### go to back ####
          Absups = Absups + integrate_abs ((B1 + B2) / 2.,
                                           Eup_3d[ilyr][ix][get_index (iy + 1, atmos->Ny)],
                                           kabs_nca[ilyr][ix][get_index (iy + 1, atmos->Ny)],
                                           kabs_nca[ilyr][ix][iy],
                                           az,
                                           bz,
                                           cz,
                                           fabs (z2 - z1),
                                           dx,
                                           mu);
          Emups  = Emups - ((B1 + B2) / 2.) * integrate_emis (kabs_nca[ilyr][ix][iy], az, bz, cz, fabs (z2 - z1), dx, mu);

          afit1 = atan (log (dx / fabs (z2 - z1)) - 0.24) * (0.39);
          bfit1 = atan (log (dx / fabs (z2 - z1)) - 0.3) * (-0.67) + 1.28;
          afit2 = atan (dx / fabs (z2 - z1) * 7.35) * 0.6 - 0.35;
          cfit2 = atan (dx / fabs (z2 - z1) - 0.36) * 1.14 + (1 * 2.02) / exp (dx / fabs (z2 - z1) * (0.48)) - 0.84;
          bfit2 = (1 - cfit2) / (PI / 2.);

          factor2 = atan (kabs_nca[ilyr][ix][iy] * dx * afit2) * bfit2 + cfit2;

          if (dx < fabs (z2 - z1))
            factor1 = atan (kabs_nca[ilyr][ix][iy] * dx) * afit1 + bfit1;
          else if (dx > fabs (z2 - z1))
            factor1 = atan (kabs_nca[ilyr][ix][iy] * fabs (z2 - z1)) * afit1 + bfit1;
          else
            factor1 = atan (kabs_nca[ilyr][ix][iy] * dx) * afit1 + bfit1;

          HR_up[ilyr]           = (Absup + Emup) / 2 / factor2 * dy / (dy) / (dx);    ///(fabs(z2-z1));
          HR_dn[ilyr]           = (Absdn + Emdn) / 2 / factor2 * dy / (dy) / (dx);    ///(fabs(z2-z1));
          HR_up_s[ilyr]         = (Absups + Emups) / 2 / factor1 * dy / (dy) / (dx);  ///(fabs(z2-z1));
          HR_dn_s[ilyr]         = (Absdns + Emdns) / 2 / factor1 * dy / (dy) / (dx);  ///(fabs(z2-z1));
          abs_nca[ilyr][ix][iy] = (HR_up[ilyr] + HR_dn[ilyr] + HR_up_s[ilyr] + HR_dn_s[ilyr]) * PI;

        } else if (ilyr < atmos->Nz - 1)
          abs_nca[ilyr][ix][iy] =
            (Eup_3d[ilyr][ix][iy] - Eup_3d[ilyr + 1][ix][iy] + Edn_3d[ilyr + 1][ix][iy] - Edn_3d[ilyr][ix][iy]) * PI;

      }  //ilyr
    }    // iy
  }      // ix

  free (HR_up);
  free (HR_dn);
  free (HR_up_s);
  free (HR_dn_s);

  return 0;

}  // end nca function

/* ############################################################################ */
/* ########### 2019 Version of NCA. Please cite Klinger and Mayer 2019   ###### */
/* ########### when using this method .                                  ###### */
/* ########### NCA calculates thermal heating rates by only using input  ###### */
/* ########### directly neighboring columns. Emissivities and correction ###### */
/* ########### factors are tabulated and stored in lookup tables.        ###### */
/* ############################################################################ */
/* ############      This is the version for a cuboid/cubic grid    ########### */
/* ############################################################################ */

int nca3d_2 (char*              datapath,
             double             dx,
             double             dy,
             atmosphere_struct* atmos,
             double***          Edn_3d,
             double***          Eup_3d,
             double***          planck,
             double***          kabs_nca,
             double***          abs_nca) {

  /* local variables */
  double* HR_up   = calloc (atmos->Nz + 2, sizeof (double));
  double* HR_dn   = calloc (atmos->Nz + 2, sizeof (double));
  double* HR_up_s = calloc (atmos->Nz + 2, sizeof (double));
  double* HR_dn_s = calloc (atmos->Nz + 2, sizeof (double));

  int ilyr = 0;  // counter z levels
  int ix   = 0;  // counter x grid boxes
  int iy   = 0;  // counter y grid boxes

  double Absup  = 0.0;  // Upwelling, Absorption, lower/upper face
  double Absdn  = 0.0;  // Downwelling Absorption, lower/upper face
  double Absups = 0.0;  // Upwelling, Absorption, face 1
  double Absdns = 0.0;  // Downwelling Absorption, face 1

  double Emup  = 0.0;  // Upwelling, Emission, lower/upper face
  double Emdn  = 0.0;  // Downwelling Emission, lower/upper face
  double Emups = 0.0;  // Upwelling, Emission, face 1
  double Emdns = 0.0;  // Downwelling Emission, face 1

  double z1 = 0.0;  // height level 1
  double z2 = 0.0;  // height level 2

  double B1 = 0.0;  // Planck emission at layer 1
  double B2 = 0.0;  // Planck emission at layer 2

  double cz = 0.0;  // integration boarder for side contributions

  double l     = 0.0;
  double Trans = 0.0;

  double A_top  = 0.0;
  double A_side = 0.0;
  double V      = 0.0;
  double taux   = 0.0;
  double tauy   = 0.0;
  double tauz   = 0.0;
  double asp    = 0.0;

  double eps_side = 0.0;
  double eps_top  = 0.0;

  double f_final_t = 0.0;
  double f_final_s = 0.0;
  double f1 = 0.0, f2 = 0.0;

  double wa = 0.0, wb = 0.0, wc = 0.0, w1 = 0.0, w2 = 0.0;

  /* ################################################################################## */
  /* ######################## lookup table for emissivity ############################# */
  /* ################################################################################## */

  double** tmp    = NULL;
  int*     nn     = NULL;
  int      status = 0, n = 0, counter = 0, col = 0, i = 0, j = 0, k = 0;

  tmp = calloc (4, sizeof (double*));
  nn  = calloc (4, sizeof (int));

  double*   tau1    = NULL;
  double*   tau2    = NULL;
  double*   tau3    = NULL;
  double*** eps_tab = NULL;

  int ntaux = 0, ntauy = 0, ntauz = 0;

  char* filename = malloc (strlen (datapath) + strlen ("nca_lookup/lookup_tau_cuboid_sort.dat") + 1);
  strcpy (filename, datapath);
  strcat (filename, "nca_lookup/lookup_tau_cuboid_sort.dat");

  status = read_4c_file (filename, &(tmp[0]), &(tmp[1]), &(tmp[2]), &(tmp[3]), &n);
  if (status != 0) {
    fprintf (stderr, "Error %d reading Lookup Table of NCA 2019 Cubiod  %s\n", status, "lookup_tau_cuboid.dat");
    return status;
  }

  for (col = 0; col < 3; col++) {
    counter = 1;
    for (i = 1; i < n; i++) {
      if (tmp[col][i] > tmp[col][i - 1])
        counter++;

      if (tmp[col][i] < tmp[col][i - 1])
        break;
    }
    nn[col] = counter;
  }

  tau1    = calloc (nn[0], sizeof (double));
  tau2    = calloc (nn[1], sizeof (double));
  tau3    = calloc (nn[2], sizeof (double));
  eps_tab = calloc (nn[2], sizeof (double**));

  for (ix = 0; ix < nn[0]; ix++) {
    eps_tab[ix] = calloc (nn[1], sizeof (double*));
    for (iy = 0; iy < nn[1]; iy++)
      eps_tab[ix][iy] = calloc (nn[2], sizeof (double));
  }
  /* write values to arrays */

  for (i = 0; i < nn[0]; i++) {
    tau1[i] = tmp[0][i];
    ntaux   = nn[0];
  }
  for (i = 0; i < nn[1]; i++) {
    tau2[i] = tmp[1][i * nn[0]];
    ntauy   = nn[1];
  }
  for (i = 0; i < nn[2]; i++) {
    tau3[i] = tmp[2][i * nn[1] * nn[0]];
    ntauz   = nn[2];
  }

  for (i = 0; i < nn[0]; i++) {
    for (j = 0; j < nn[1]; j++) {
      for (k = 0; k < nn[2]; k++) {
        eps_tab[i][j][k] = tmp[3][i + j * nn[0] + k * nn[1] * nn[0]];
      }
    }
  }

  free (tmp);
  free (nn);
  free (filename);

  /* ################################################################################## */
  /* ######################## lookup table for correction ############################# */
  /* ################################################################################## */

  int* nn1 = NULL;
  tmp      = calloc (3, sizeof (double*));
  nn1      = calloc (3, sizeof (int));

  double*  asp_arr1     = NULL;
  double*  tau_top      = NULL;
  double** corr_tab_top = NULL;

  int ntau_t = 0;
  int nasp_t = 0;

  filename = malloc (strlen (datapath) + strlen ("/nca_lookup/lookup_correct_cuboid_top_sort.dat") + 1);
  strcpy (filename, datapath);
  strcat (filename, "/nca_lookup/lookup_correct_cuboid_top_sort.dat");

  status = read_3c_file (filename, &(tmp[0]), &(tmp[1]), &(tmp[2]), &n);
  if (status != 0) {
    fprintf (stderr, "Error %d reading Lookup Table of NCA 2019 Cubiod  %s\n", status, "lookup_correct_cuboid_top.dat");
    return status;
  }

  for (col = 0; col < 2; col++) {
    counter = 1;
    for (i = 1; i < n; i++) {
      if (tmp[col][i] > tmp[col][i - 1])
        counter++;

      if (tmp[col][i] < tmp[col][i - 1])
        break;
    }
    nn1[col] = counter;
  }

  asp_arr1     = calloc (nn1[0], sizeof (double));
  tau_top      = calloc (nn1[1], sizeof (double));
  corr_tab_top = calloc (nn1[1], sizeof (double*));

  for (ix = 0; ix < nn1[0]; ix++) {
    corr_tab_top[ix] = calloc (nn1[1], sizeof (double));
  }

  /* write values to arrays */
  for (i = 0; i < nn1[0]; i++) {
    asp_arr1[i] = tmp[0][i];
    nasp_t      = nn1[0];
  }
  for (i = 0; i < nn1[1]; i++) {
    tau_top[i] = tmp[1][i * nn1[0]];
    ntau_t     = nn1[1];
  }
  for (i = 0; i < nn1[0]; i++) {
    for (j = 0; j < nn1[1]; j++) {
      corr_tab_top[i][j] = tmp[2][i + nn1[0] * j];
    }
  }

  free (nn1);
  free (tmp);
  free (filename);

  /* ################################################################## */

  int* nn2 = NULL;
  tmp      = calloc (3, sizeof (double*));
  nn2      = calloc (3, sizeof (int));

  double*  asp_arr2      = NULL;
  double*  tau_side      = NULL;
  double** corr_tab_side = NULL;
  int      nasp_s        = 0;
  int      ntau_s        = 0;

  filename = malloc (strlen (datapath) + strlen ("/nca_lookup/lookup_correct_cuboid_side_sort.dat") + 1);
  strcpy (filename, datapath);
  strcat (filename, "/nca_lookup/lookup_correct_cuboid_side_sort.dat");

  status = read_3c_file (filename, &(tmp[0]), &(tmp[1]), &(tmp[2]), &n);
  if (status != 0) {
    fprintf (stderr, "Error %d reading Lookup Table of NCA 2019 Cubiod  %s\n", status, "lookup_correct_cuboid_side.dat");
    return status;
  }

  for (col = 0; col < 2; col++) {
    counter = 1;
    for (i = 1; i < n; i++) {
      if (tmp[col][i] > tmp[col][i - 1])
        counter++;

      if (tmp[col][i] < tmp[col][i - 1])
        break;
    }

    nn2[col] = counter;
  }

  asp_arr2      = calloc (nn2[0], sizeof (double));
  tau_side      = calloc (nn2[1], sizeof (double));
  corr_tab_side = calloc (nn2[1], sizeof (double*));

  for (ix = 0; ix < nn2[0]; ix++) {
    corr_tab_side[ix] = calloc (nn2[1], sizeof (double));
  }

  /* write values to arrays */

  for (i = 0; i < nn2[0]; i++) {
    asp_arr2[i] = tmp[0][i];
    nasp_s      = nn2[0];
  }
  for (i = 0; i < nn2[1]; i++) {
    tau_side[i] = tmp[1][i * nn2[0]];
    ntau_s      = nn2[1];
  }
  for (i = 0; i < nn2[0]; i++) {
    for (j = 0; j < nn2[1]; j++) {
      corr_tab_side[i][j] = tmp[2][i + nn2[0] * j];
    }
  }

  free (nn2);
  free (tmp);

  /* ####################################################################### */
  /* ########## Start NCA    ############################################### */
  /* ####################################################################### */

  for (ix = 0; ix < atmos->Nx; ix++) {                /* loop over all x gridboxes */
    for (iy = 0; iy < atmos->Ny; iy++) {              /* loop over all y gridboxes */
      for (ilyr = 0; ilyr <= atmos->Nz - 1; ilyr++) { /* loop over all height levels */

        if (atmos->threed[MCCAOTH_TOT][ilyr] >= 1) { /* only use NCA in 3D layers */

          z1 = atmos->Z[ilyr] / 1000.;
          z2 = atmos->Z[ilyr + 1] / 1000.;
          cz = fabs (z2 - z1);

          B1 = planck[ilyr][ix][iy];
          B2 = planck[ilyr + 1][ix][iy];

          A_top  = dx * dy;
          A_side = dx * cz;
          V      = dx * dy * cz;

          /* ### get emission from lookup table */

          taux    = kabs_nca[ilyr][ix][iy] * dx;
          tauy    = kabs_nca[ilyr][ix][iy] * dy;
          tauz    = kabs_nca[ilyr][ix][iy] * cz;
          eps_top = pol_3d (taux, tauy, tauz, eps_tab, tau1, tau2, tau3, ntaux, ntauy, ntauz);

          tauz     = kabs_nca[ilyr][ix][iy] * dx;
          tauy     = kabs_nca[ilyr][ix][iy] * dy;
          taux     = kabs_nca[ilyr][ix][iy] * cz;
          eps_side = pol_3d (taux, tauy, tauz, eps_tab, tau1, tau2, tau3, ntaux, ntauy, ntauz);

          asp = cz / dx;
          if (asp > 10)
            asp = 10;
          else if (asp < 0.1)
            asp = 0.1;

          taux = kabs_nca[ilyr][ix][iy] * dx;
          tauz = kabs_nca[ilyr][ix][iy] * cz;

          f_final_t = pol_2d (asp, tauz, corr_tab_top, asp_arr1, tau_top, nasp_t, ntau_t, 1);
          f_final_s = pol_2d (asp, taux, corr_tab_side, asp_arr2, tau_side, nasp_s, ntau_s, 1);

          /* #### top/bottom #### */
          if (ilyr < atmos->Nz) {
            wa = atan (asp * 1.21) * (-0.76) + 1.31;
            wb = pow (asp, 0.028) * (-8.08) + asp * 0.01 + atan (asp * 0.13) + 7.49;
            wc = pow (asp, 0.44) * (1.55) + asp * (-0.25) + atan (asp * (-0.3)) - 0.31;
            w1 = atan (kabs_nca[ilyr + 1][ix][iy] * cz * wa) * wb + wc;
            w2 = 1 - w1;

            Trans = Edn_3d[ilyr + 1][get_index (ix + 1, atmos->Nx)][iy] + Edn_3d[ilyr + 1][get_index (ix - 1, atmos->Nx)][iy] +
                    Edn_3d[ilyr + 1][ix][get_index (iy - 1, atmos->Ny)] + Edn_3d[ilyr + 1][ix][get_index (iy + 1, atmos->Ny)];

            l = Trans / 4. * w1 + Edn_3d[ilyr + 1][ix][iy] * w2;
          } else {
            l = Edn_3d[ilyr + 1][ix][iy];
          }

          Absdn = l;
          Emdn  = -(B1 + B2) / 2.;

          /* #### go to right #### */
          f1 = atan (kabs_nca[ilyr][get_index (ix + 1, atmos->Nx)][iy] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absdns =
            (f1 * Edn_3d[ilyr + 1][get_index (ix + 1, atmos->Nx)][iy] + f2 * Edn_3d[ilyr][get_index (ix + 1, atmos->Nx)][iy]);
          Emdns = -(B1 + B2) / 2.;

          /* #### go to left #### */
          f1 = atan (kabs_nca[ilyr][get_index (ix - 1, atmos->Nx)][iy] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absdns = Absdns + (f1 * Edn_3d[ilyr + 1][get_index (ix - 1, atmos->Nx)][iy] +
                             f2 * Edn_3d[ilyr][get_index (ix - 1, atmos->Nx)][iy]);
          Emdns  = Emdns - (B1 + B2) / 2.;

          /* #### go to back #### */
          f1 = atan (kabs_nca[ilyr][ix][get_index (iy - 1, atmos->Ny)] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absdns = Absdns + (f1 * Edn_3d[ilyr + 1][ix][get_index (iy - 1, atmos->Ny)] +
                             f2 * Edn_3d[ilyr][ix][get_index (iy - 1, atmos->Ny)]);
          Emdns  = Emdns - (B1 + B2) / 2.;

          /* #### go to front #### */
          f1 = atan (kabs_nca[ilyr][ix][get_index (iy + 1, atmos->Ny)] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absdns = Absdns + (f1 * Edn_3d[ilyr + 1][ix][get_index (iy + 1, atmos->Ny)] +
                             f2 * Edn_3d[ilyr][ix][get_index (iy + 1, atmos->Ny)]);
          Emdns  = Emdns - (B1 + B2) / 2.;

          /* ## upwelling */
          /* bottom face of gridbox */
          if (ilyr > 0) {

            wa = atan (asp * 1.21) * (-0.76) + 1.31;
            wb = pow (asp, 0.028) * (-8.08) + asp * 0.01 + atan (asp * 0.13) + 7.49;
            wc = pow (asp, 0.44) * (1.55) + asp * (-0.25) + atan (asp * (-0.3)) - 0.31;
            w1 = atan (kabs_nca[ilyr - 1][ix][iy] * cz * wa) * wb + wc;
            w2 = 1 - w1;

            Trans = Eup_3d[ilyr][get_index (ix + 1, atmos->Nx)][iy] + Eup_3d[ilyr][get_index (ix - 1, atmos->Nx)][iy] +
                    Eup_3d[ilyr][ix][get_index (iy - 1, atmos->Ny)] + Eup_3d[ilyr][ix][get_index (iy + 1, atmos->Ny)];
            l = Trans / 4. * w1 + w2 * Eup_3d[ilyr][ix][iy];
          } else
            l = Eup_3d[0][ix][iy];

          Absup = l;
          Emup  = -(B1 + B2) / 2.;

          // #### go to right ####
          f1 = atan (kabs_nca[ilyr][get_index (ix + 1, atmos->Nx)][iy] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absups =
            (f1 * Eup_3d[ilyr][get_index (ix + 1, atmos->Nx)][iy] + f2 * Eup_3d[ilyr + 1][get_index (ix + 1, atmos->Nx)][iy]);
          Emups = -(B1 + B2) / 2.;

          // #### go to left ####
          f1 = atan (kabs_nca[ilyr][get_index (ix - 1, atmos->Nx)][iy] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absups = Absups + (f1 * Eup_3d[ilyr][get_index (ix - 1, atmos->Nx)][iy] +
                             f2 * Eup_3d[ilyr + 1][get_index (ix - 1, atmos->Nx)][iy]);
          Emups  = Emups - (B1 + B2) / 2.;

          // #### go to back ####
          f1 = atan (kabs_nca[ilyr][ix][get_index (iy - 1, atmos->Ny)] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absups = Absups + (f1 * Eup_3d[ilyr][ix][get_index (iy - 1, atmos->Ny)] +
                             f2 * Eup_3d[ilyr + 1][ix][get_index (iy - 1, atmos->Ny)]);
          Emups  = Emups - (B1 + B2) / 2.;

          // #### go to front ####
          f1 = atan (kabs_nca[ilyr][ix][get_index (iy + 1, atmos->Ny)] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absups = Absups + (f1 * Eup_3d[ilyr][ix][get_index (iy + 1, atmos->Ny)] +
                             f2 * Eup_3d[ilyr + 1][ix][get_index (iy + 1, atmos->Ny)]);
          Emups  = Emups - (B1 + B2) / 2.;

          HR_up[ilyr]   = (Absup + Emup) * A_top / V * eps_top * f_final_t;
          HR_dn[ilyr]   = (Absdn + Emdn) * A_top / V * eps_top * f_final_t;
          HR_up_s[ilyr] = (Emups + Absups) / 2 * A_side / V * eps_side * f_final_s;
          HR_dn_s[ilyr] = (Emdns + Absdns) / 2 * A_side / V * eps_side * f_final_s;

          abs_nca[ilyr][ix][iy] = (HR_up[ilyr] + HR_dn[ilyr] + HR_up_s[ilyr] + HR_dn_s[ilyr]) * PI * (fabs (z2 - z1));
          ;
        } else if (ilyr < atmos->Nz - 1)
          abs_nca[ilyr][ix][iy] =
            (Eup_3d[ilyr][ix][iy] - Eup_3d[ilyr + 1][ix][iy] + Edn_3d[ilyr + 1][ix][iy] - Edn_3d[ilyr][ix][iy]) * PI *
            (fabs (z2 - z1));
        ;
      }  //ilyr
    }    // iy
  }      // ix

  free (HR_up);
  free (HR_dn);
  free (HR_up_s);
  free (HR_dn_s);

  return 0;

}  // end nca 2019 cubiod function

/* ############################################################################ */
/* ########### 2019 Version of NCA. Please cite Klinger and Mayer 2019   ###### */
/* ########### when using this method .                                  ###### */
/* ########### NCA calculates thermal heating rates by only using input  ###### */
/* ########### directly neighboring columns. Emissivities and correction ###### */
/* ########### factors are tabulated and stored in lookup tables.        ###### */
/* ############################################################################ */
/* ####### This is the version for a triangular grid based on cubic grid ###### */
/* ############################################################################ */

int nca3d_3 (char*              datapath,
             double             dx,
             double             dy,
             atmosphere_struct* atmos,
             double***          Edn_3d,
             double***          Eup_3d,
             double***          planck,
             double***          kabs_nca,
             double***          abs_nca) {

  /* local variables */
  double* HR_up   = calloc (atmos->Nz + 2, sizeof (double));
  double* HR_dn   = calloc (atmos->Nz + 2, sizeof (double));
  double* HR_up_s = calloc (atmos->Nz + 2, sizeof (double));
  double* HR_dn_s = calloc (atmos->Nz + 2, sizeof (double));

  int ilyr = 0;  // counter z levels
  int ix   = 0;  // counter x grid boxes
  int iy   = 0;  // counter y grid boxes

  double Absup  = 0.0;  // Upwelling, Absorption, lower/upper face
  double Absdn  = 0.0;  // Downwelling Absorption, lower/upper face
  double Absups = 0.0;  // Upwelling, Absorption, face 1
  double Absdns = 0.0;  // Downwelling Absorption, face 1

  double Emup  = 0.0;  // Upwelling, Emission, lower/upper face
  double Emdn  = 0.0;  // Downwelling Emission, lower/upper face
  double Emups = 0.0;  // Upwelling, Emission, face 1
  double Emdns = 0.0;  // Downwelling Emission, face 1

  double z1 = 0.0;  // height level 1
  double z2 = 0.0;  // height level 2

  double B1 = 0.0;  // Planck emission at layer 1
  double B2 = 0.0;  // Planck emission at layer 2

  double cz = 0.0;  // integration boarder for side contributions

  double l     = 0.0;
  double Trans = 0.0;

  double A_top   = 0.0;
  double A_side1 = 0.0;
  double A_side2 = 0.0;
  double A_side3 = 0.0;
  double V       = 0.0;
  double tauz    = 0.0;
  double tauhx   = 0.0;
  double asp     = 0.0;

  double eps_top   = 0.0;
  double eps_side1 = 0.0;
  double eps_side2 = 0.0;
  double eps_side3 = 0.0;

  double f_final  = 0.0;
  double f_final1 = 0.0;
  double f_final2 = 0.0;
  double f_final3 = 0.0;
  double f1 = 0.0, f2 = 0.0;

  double s = 0.0, dx2 = 0.0;

  double wa = 0.0, wb = 0.0, wc = 0.0, w1 = 0.0, w2 = 0.0;

  /* ################################################################################## */
  /* ######################## lookup table for emissivity ############################# */
  /* ################################################################################## */

  double** tmp    = NULL;
  int*     nn1    = NULL;
  tmp             = calloc (3, sizeof (double*));
  nn1             = calloc (3, sizeof (int));
  int      status = 0, n = 0, counter = 0, col = 0, i = 0, j = 0;
  double*  tau_top1    = NULL;
  double*  tau_top2    = NULL;
  double** eps_top_tab = NULL;

  int ntau1_t = 0;
  int ntau2_t = 0;

  char* filename = malloc (strlen (datapath) + strlen ("/nca_lookup/lookup_top_triangle_sort.dat") + 1);
  strcpy (filename, datapath);
  strcat (filename, "/nca_lookup/lookup_top_triangle_sort.dat");

  status = read_3c_file (filename, &(tmp[0]), &(tmp[1]), &(tmp[2]), &n);
  if (status != 0) {
    fprintf (stderr, "Error %d reading Lookup Table of NCA 2019 Triangle  %s\n", status, "lookup_top_triangle_sort.dat");
    return status;
  }

  for (col = 0; col < 2; col++) {
    counter = 1;
    for (i = 1; i < n; i++) {
      if (tmp[col][i] > tmp[col][i - 1])
        counter++;

      if (tmp[col][i] < tmp[col][i - 1])
        break;
    }
    nn1[col] = counter;
  }

  tau_top1    = calloc (nn1[0], sizeof (double));
  tau_top2    = calloc (nn1[1], sizeof (double));
  eps_top_tab = calloc (nn1[1], sizeof (double*));

  for (ix = 0; ix < nn1[0]; ix++) {
    eps_top_tab[ix] = calloc (nn1[1], sizeof (double));
  }

  /* write values to arrays */
  for (i = 0; i < nn1[0]; i++) {
    tau_top1[i] = tmp[0][i];
    ntau1_t     = nn1[0];
  }
  for (i = 0; i < nn1[1]; i++) {
    tau_top2[i] = tmp[1][i * nn1[0]];
    ntau2_t     = nn1[1];
  }
  for (i = 0; i < nn1[0]; i++) {
    for (j = 0; j < nn1[1]; j++) {
      eps_top_tab[i][j] = tmp[2][i + nn1[0] * j];
    }
  }

  free (nn1);
  free (tmp);
  free (filename);

  /* ################################################################################## */
  /* ######################## lookup table for correction ############################# */
  /* ################################################################################## */

  /* int *nn1=NULL; */
  tmp = calloc (3, sizeof (double*));
  nn1 = calloc (3, sizeof (int));

  double*  tau_side1    = NULL;
  double*  tau_side2    = NULL;
  double** eps_side_tab = NULL;

  int ntau1_s = 0;
  int ntau2_s = 0;

  filename = malloc (strlen (datapath) + strlen ("/nca_lookup/lookup_side_triangle_sort.dat") + 1);
  strcpy (filename, datapath);
  strcat (filename, "/nca_lookup/lookup_side_triangle_sort.dat");

  status = read_3c_file (filename, &(tmp[0]), &(tmp[1]), &(tmp[2]), &n);
  if (status != 0) {
    fprintf (stderr, "Error %d reading Lookup Table of NCA 2019 Triangle  %s\n", status, "lookup_side_triangle_sort.dat");
    return status;
  }

  for (col = 0; col < 2; col++) {
    counter = 1;
    for (i = 1; i < n; i++) {
      if (tmp[col][i] > tmp[col][i - 1])
        counter++;

      if (tmp[col][i] < tmp[col][i - 1])
        break;
    }
    nn1[col] = counter;
  }

  tau_side1    = calloc (nn1[0], sizeof (double));
  tau_side2    = calloc (nn1[1], sizeof (double));
  eps_side_tab = calloc (nn1[1], sizeof (double*));

  for (ix = 0; ix < nn1[0]; ix++) {
    eps_side_tab[ix] = calloc (nn1[1], sizeof (double));
  }

  /* write values to arrays */
  for (i = 0; i < nn1[0]; i++) {
    tau_side1[i] = tmp[0][i];
    ntau1_s      = nn1[0];
  }
  for (i = 0; i < nn1[1]; i++) {
    tau_side2[i] = tmp[1][i * nn1[0]];
    ntau2_s      = nn1[1];
  }
  for (i = 0; i < nn1[0]; i++) {
    for (j = 0; j < nn1[1]; j++) {
      eps_side_tab[i][j] = tmp[2][i + nn1[0] * j];
    }
  }

  free (nn1);
  free (tmp);
  free (filename);

  /* ################################################################## */
  /* int *nn1=NULL; */
  tmp = calloc (3, sizeof (double*));
  nn1 = calloc (3, sizeof (int));

  double*  asp_top11 = NULL;
  double*  tau_top22 = NULL;
  double** corr_top  = NULL;

  int nasp11_t = 0;
  int ntau22_t = 0;

  filename = malloc (strlen (datapath) + strlen ("/nca_lookup/lookup_correct_triangle_top_sort.dat") + 1);
  strcpy (filename, datapath);
  strcat (filename, "/nca_lookup/lookup_correct_triangle_top_sort.dat");

  status = read_3c_file (filename, &(tmp[0]), &(tmp[1]), &(tmp[2]), &n);
  if (status != 0) {
    fprintf (stderr, "Error %d reading Lookup Table of NCA 2019 Triangle  %s\n", status, "lookup_correct_triangle_top_sort.dat");
    return status;
  }

  for (col = 0; col < 2; col++) {
    counter = 1;
    for (i = 1; i < n; i++) {
      if (tmp[col][i] > tmp[col][i - 1])
        counter++;

      if (tmp[col][i] < tmp[col][i - 1])
        break;
    }
    nn1[col] = counter;
  }

  asp_top11 = calloc (nn1[0], sizeof (double));
  tau_top22 = calloc (nn1[1], sizeof (double));
  corr_top  = calloc (nn1[1], sizeof (double*));

  for (ix = 0; ix < nn1[0]; ix++) {
    corr_top[ix] = calloc (nn1[1], sizeof (double));
  }

  /* write values to arrays */
  for (i = 0; i < nn1[0]; i++) {
    asp_top11[i] = tmp[0][i];
    nasp11_t     = nn1[0];
  }
  for (i = 0; i < nn1[1]; i++) {
    tau_top22[i] = tmp[1][i * nn1[0]];
    ntau22_t     = nn1[1];
  }
  for (i = 0; i < nn1[0]; i++) {
    for (j = 0; j < nn1[1]; j++) {
      corr_top[i][j] = tmp[2][i + nn1[0] * j];
    }
  }

  free (nn1);
  free (tmp);
  free (filename);
  /* ################################################################## */
  /* int *nn1=NULL; */
  tmp = calloc (3, sizeof (double*));
  nn1 = calloc (3, sizeof (int));

  double*  asp_side11 = NULL;
  double*  tau_side22 = NULL;
  double** corr_side  = NULL;

  int nasp11_s = 0;
  int ntau22_s = 0;

  filename = malloc (strlen (datapath) + strlen ("/nca_lookup/lookup_correct_triangle_side_sort.dat") + 1);
  strcpy (filename, datapath);
  strcat (filename, "/nca_lookup/lookup_correct_triangle_side_sort.dat");

  status = read_3c_file (filename, &(tmp[0]), &(tmp[1]), &(tmp[2]), &n);
  if (status != 0) {
    fprintf (stderr, "Error %d reading Lookup Table of NCA 2019 Triangle  %s\n", status, "lookup_correct_triangle_side_sort.dat");
    return status;
  }

  for (col = 0; col < 2; col++) {
    counter = 1;
    for (i = 1; i < n; i++) {
      if (tmp[col][i] > tmp[col][i - 1])
        counter++;

      if (tmp[col][i] < tmp[col][i - 1])
        break;
    }
    nn1[col] = counter;
  }

  asp_side11 = calloc (nn1[0], sizeof (double));
  tau_side22 = calloc (nn1[1], sizeof (double));
  corr_side  = calloc (nn1[1], sizeof (double*));

  for (ix = 0; ix < nn1[0]; ix++) {
    corr_side[ix] = calloc (nn1[1], sizeof (double));
  }

  /* write values to arrays */
  for (i = 0; i < nn1[0]; i++) {
    asp_side11[i] = tmp[0][i];
    nasp11_s      = nn1[0];
  }
  for (i = 0; i < nn1[1]; i++) {
    tau_side22[i] = tmp[1][i * nn1[0]];
    ntau22_s      = nn1[1];
  }
  for (i = 0; i < nn1[0]; i++) {
    for (j = 0; j < nn1[1]; j++) {
      corr_side[i][j] = tmp[2][i + nn1[0] * j];
    }
  }

  free (nn1);
  free (tmp);
  free (filename);
  /* ####################################################################### */
  /* ########## Start NCA    ############################################### */
  /* ####################################################################### */

  for (ix = 0; ix < atmos->Nx; ix++) {                /* loop over all x gridboxes */
    for (iy = 0; iy < atmos->Ny; iy++) {              /* loop over all y gridboxes */
      for (ilyr = 0; ilyr <= atmos->Nz - 1; ilyr++) { /* loop over all height levels */

        if (atmos->threed[MCCAOTH_TOT][ilyr] >= 1) { /* only use NCA in 3D layers */
          z1 = atmos->Z[ilyr] / 1000.;
          z2 = atmos->Z[ilyr + 1] / 1000.;
          cz = fabs (z2 - z1);

          B1 = planck[ilyr][ix][iy];
          B2 = planck[ilyr + 1][ix][iy];

          dx2 = sqrt (dx * dx + dy * dy);
          //dx2=dx; single cloud solution

          //Check! This is valid as long as dx, dx2 and dy are similar! Otherwise new volumes and areas have to calcuated
          //Herons Formula

          s       = (dx + dx2 + dy) / 2;
          A_top   = sqrt (s * (s - dx) * (s - dx2) * (s - dy));
          A_side1 = dx * cz;
          A_side2 = dx2 * cz;
          A_side3 = dy * cz;

          V = A_top * cz;

          asp = cz / ((2 / ((dy + dx2 + dx) / 3)) * A_top);

          if (asp > 11)
            asp = 11;
          else if (asp < 0.118)
            asp = 0.118;

          /* ### get emission from lookup table */
          tauhx = kabs_nca[ilyr][ix][iy] * dx * 0.86603;
          tauz  = kabs_nca[ilyr][ix][iy] * cz;

          // find emissivity for top face
          eps_top = pol_2d (tauhx, tauz, eps_top_tab, tau_top1, tau_top2, ntau1_t, ntau2_t, 2);

          // find emissivity for side face
          eps_side1 = pol_2d (tauhx, tauz, eps_side_tab, tau_side1, tau_side2, ntau1_s, ntau2_s, 2);
          f_final1  = pol_2d (asp, tauhx, corr_side, asp_side11, tau_side22, nasp11_s, ntau22_s, 1);

          // find emissivity for side face
          tauhx     = kabs_nca[ilyr][ix][iy] * dx2 * 0.86603;
          eps_side2 = pol_2d (tauhx, tauz, eps_side_tab, tau_side1, tau_side2, ntau1_s, ntau2_s, 2);
          f_final2  = pol_2d (asp, tauhx, corr_side, asp_side11, tau_side22, nasp11_s, ntau22_s, 1);

          //find emissivity for side face
          tauhx     = kabs_nca[ilyr][ix][iy] * dy * 0.86603;
          eps_side3 = pol_2d (tauhx, tauz, eps_side_tab, tau_side1, tau_side2, ntau1_s, ntau2_s, 2);
          f_final3  = pol_2d (asp, tauhx, corr_side, asp_side11, tau_side22, nasp11_s, ntau22_s, 1);

          /* #### top/bottom #### */
          if (ilyr < atmos->Nz) {
            asp = cz / ((2 / ((dy + dx2 + dx) / 3)) * A_top);
            if (asp > 11)
              asp = 11;
            else if (asp < 0.118)
              asp = 0.118;

            wa = atan (asp * 1.29) * (-0.75) + 1.21;
            wb = pow (asp, 0.027) * (-7.98) + asp * (-0.01) + atan (asp * 0.11) + 7.36;
            wc = pow (asp, 0.49) * (1.46) + asp * (-0.25) + atan (asp * (-0.29)) - 0.12;
            w1 = atan (kabs_nca[ilyr + 1][ix][iy] * cz * wa) * wb + wc;
            w2 = 1 - w1;

            Trans = Edn_3d[ilyr + 1][get_index (ix + 1, atmos->Nx)][iy] + Edn_3d[ilyr + 1][get_index (ix - 1, atmos->Nx)][iy] +
                    Edn_3d[ilyr + 1][ix][get_index (iy - 1, atmos->Ny)] + Edn_3d[ilyr + 1][ix][get_index (iy + 1, atmos->Ny)] +
                    Edn_3d[ilyr + 1][ix][iy] * 2;

            l = (w1 * Trans / 6. + Edn_3d[ilyr + 1][ix][iy] * w2) * 2.;
            /* Singe Cloud */
            //l = Edn_3d[ilyr+1][ix][iy];
          } else {
            l = Edn_3d[ilyr + 1][ix][iy];
          }

          Absdn = l;
          Emdn  = -(B1 + B2);

          /* Singe Cloud */
          //Emdn = -(B1+B2)/2.;

          /*  #### CENTER ####  */
          f1 = atan (kabs_nca[ilyr][ix][iy] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absdns = 2 * (f1 * Edn_3d[ilyr + 1][ix][iy] + f2 * Edn_3d[ilyr][ix][iy]) * A_side2 * eps_side2 * f_final2;
          Emdns  = -(B1 + B2) * A_side2 * eps_side2 * f_final2;

          /*  #### go to right #### */
          f1 = atan (kabs_nca[ilyr][get_index (ix + 1, atmos->Nx)][iy] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;
          /* /\* for isolated single comparision *\/ */
          /* /\* for isolated single comparision  - comment out two of the side contributions. only 2 are needed in total*\/ */
          Absdns = Absdns + (f1 * Edn_3d[ilyr + 1][get_index (ix + 1, atmos->Nx)][iy] +
                             f2 * Edn_3d[ilyr][get_index (ix + 1, atmos->Nx)][iy]) *
                              A_side3 * eps_side3 * f_final3;
          Emdns = Emdns - (B1 + B2) / 2. * A_side3 * eps_side3 * f_final3;

          /* #### go to left #### */
          f1 = atan (kabs_nca[ilyr][get_index (ix - 1, atmos->Nx)][iy] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absdns = Absdns + (f1 * Edn_3d[ilyr + 1][get_index (ix - 1, atmos->Nx)][iy] +
                             f2 * Edn_3d[ilyr][get_index (ix - 1, atmos->Nx)][iy]) *
                              A_side3 * eps_side3 * f_final3;
          Emdns = Emdns - (B1 + B2) / 2. * A_side3 * eps_side3 * f_final3;

          /* #### go to back #### */
          f1 = atan (kabs_nca[ilyr][ix][get_index (iy - 1, atmos->Ny)] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absdns = Absdns + (f1 * Edn_3d[ilyr + 1][ix][get_index (iy - 1, atmos->Ny)] +
                             f2 * Edn_3d[ilyr][ix][get_index (iy - 1, atmos->Ny)]) *
                              A_side1 * eps_side1 * f_final1;
          Emdns = Emdns - (B1 + B2) / 2. * A_side1 * eps_side1 * f_final1;

          /* #### go to front #### */
          f1 = atan (kabs_nca[ilyr][ix][get_index (iy + 1, atmos->Ny)] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absdns = Absdns + (f1 * Edn_3d[ilyr + 1][ix][get_index (iy + 1, atmos->Ny)] +
                             f2 * Edn_3d[ilyr][ix][get_index (iy + 1, atmos->Ny)]) *
                              A_side1 * eps_side1 * f_final1;
          Emdns = Emdns - (B1 + B2) / 2. * A_side1 * eps_side1 * f_final1;

          /* ## upwelling */
          /* bottom face of gridbox */
          if (ilyr > 0) {
            asp = cz / ((2 / ((dx + dx2 + dy) / 3)) * A_top);
            if (asp > 11)
              asp = 11;
            else if (asp < 0.118)
              asp = 0.118;

            wa    = atan (asp * 1.29) * (-0.75) + 1.21;
            wb    = pow (asp, 0.027) * (-7.98) + asp * (-0.01) + atan (asp * 0.11) + 7.36;
            wc    = pow (asp, 0.49) * (1.46) + asp * (-0.25) + atan (asp * (-0.29)) - 0.12;
            w1    = atan (kabs_nca[ilyr - 1][ix][iy] * cz * wa) * wb + wc;
            w2    = 1 - w1;
            Trans = Eup_3d[ilyr][get_index (ix + 1, atmos->Nx)][iy] + Eup_3d[ilyr][get_index (ix - 1, atmos->Nx)][iy] +
                    Eup_3d[ilyr][ix][get_index (iy - 1, atmos->Ny)] + Eup_3d[ilyr][ix][get_index (iy + 1, atmos->Ny)] +
                    Eup_3d[ilyr][ix][iy] * 2;
            l = (w1 * Trans / 6. + Eup_3d[ilyr][ix][iy] * w2) * 2;
            /* Singe Cloud */
            //l= Eup_3d[ilyr][ix][iy];
          } else
            l = Eup_3d[0][ix][iy];

          Absup = l;

          Emup = -(B1 + B2);
          /* Singe Cloud */
          //Emup = -(B1+B2)/2;

          // #### center ####
          f1 = atan (kabs_nca[ilyr][ix][iy] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absups = 2 * (f1 * Eup_3d[ilyr][ix][iy] + f2 * Eup_3d[ilyr + 1][ix][iy]) * A_side2 * eps_side2 * f_final2;
          Emups  = -(B1 + B2) * A_side2 * eps_side2 * f_final2;

          // #### go to right ####
          f1 = atan (kabs_nca[ilyr][get_index (ix + 1, atmos->Nx)][iy] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;
          /* for isolated single comparision */
          /* for isolated single comparision  - comment out two of the side contributions. only 2 are needed in total*/
          Absups = Absups + (f1 * Eup_3d[ilyr][get_index (ix + 1, atmos->Nx)][iy] +
                             f2 * Eup_3d[ilyr + 1][get_index (ix + 1, atmos->Nx)][iy]) *
                              A_side3 * eps_side3 * f_final3;
          Emups = Emups - (B1 + B2) / 2. * A_side3 * eps_side3 * f_final3;

          // #### go to left ####
          f1 = atan (kabs_nca[ilyr][get_index (ix - 1, atmos->Nx)][iy] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absups = Absups + (f1 * Eup_3d[ilyr][get_index (ix - 1, atmos->Nx)][iy] +
                             f2 * Eup_3d[ilyr + 1][get_index (ix - 1, atmos->Nx)][iy]) *
                              A_side3 * eps_side3 * f_final3;
          Emups = Emups - (B1 + B2) / 2. * A_side3 * eps_side3 * f_final3;

          // #### go to back ####
          f1 = atan (kabs_nca[ilyr][ix][get_index (iy - 1, atmos->Ny)] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absups = Absups + (f1 * Eup_3d[ilyr][ix][get_index (iy - 1, atmos->Ny)] +
                             f2 * Eup_3d[ilyr + 1][ix][get_index (iy - 1, atmos->Ny)]) *
                              A_side3 * eps_side3 * f_final3;
          Emups = Emups - (B1 + B2) / 2. * A_side3 * eps_side3 * f_final3;

          // #### go to front ####
          f1 = atan (kabs_nca[ilyr][ix][get_index (iy + 1, atmos->Ny)] * cz * (-2.08 / asp)) * 0.31192 + 0.49;
          if (f1 < 0)
            f1 = 0.;
          f2 = 1 - f1;

          Absups = Absups + (f1 * Eup_3d[ilyr][ix][get_index (iy + 1, atmos->Ny)] +
                             f2 * Eup_3d[ilyr + 1][ix][get_index (iy + 1, atmos->Ny)]) *
                              A_side3 * eps_side3 * f_final3;
          Emups = Emups - (B1 + B2) / 2. * A_side3 * eps_side3 * f_final3;

          asp = cz / ((2 / ((dy + dx2 + dx) / 3)) * A_top);
          if (asp < 0.118)
            asp = 0.118;
          else if (asp > 11)
            asp = 11;

          f_final = pol_2d (asp, tauz, corr_top, asp_top11, tau_top22, nasp11_t, ntau22_t, 1);

          HR_up[ilyr]   = (Absup + Emup) * A_top / V * eps_top * f_final;
          HR_dn[ilyr]   = (Absdn + Emdn) * A_top / V * eps_top * f_final;
          HR_up_s[ilyr] = (Emups + Absups) / 2. / V;
          HR_dn_s[ilyr] = (Emdns + Absdns) / 2. / V;
          //fprintf(stderr, "%f %f %f %f %f %f %f %f \n", tauhx, tauz, HR_up[ilyr],HR_dn[ilyr],HR_up_s[ilyr],HR_dn_s[ilyr], eps_top, eps_side3 );
          abs_nca[ilyr][ix][iy] = (HR_up[ilyr] + HR_dn[ilyr] + HR_up_s[ilyr] + HR_dn_s[ilyr]) * PI * (fabs (z2 - z1)) / 2.;
          //abs_nca[ilyr][ix][iy] = (HR_up[ilyr]+HR_dn[ilyr]+HR_up_s[ilyr]+HR_dn_s[ilyr])*PI*(fabs(z2-z1));
        } else if (ilyr < atmos->Nz - 1)
          abs_nca[ilyr][ix][iy] =
            (Eup_3d[ilyr][ix][iy] - Eup_3d[ilyr + 1][ix][iy] + Edn_3d[ilyr + 1][ix][iy] - Edn_3d[ilyr][ix][iy]) * PI *
            (fabs (z2 - z1));
      }  //ilyr
    }    // iy
  }      // ix

  free (HR_up);
  free (HR_dn);
  free (HR_up_s);
  free (HR_dn_s);

  return 0;

}  // end nca 2019 triangle function

// ################################################################################
// ####################### Function - integrate_abs ###############################
// # This function integrates the contributiong emission/absorption of a grid box #
//################################################################################

double integrate_abs (double B_planck,
                      double L,
                      double kabs1,
                      double kabs2,
                      double a,
                      double b,
                      double c,
                      double delta_1,
                      double delta_2,
                      double mu) {

  const double sin_mu = sqrt (1. - mu * mu);
  const double LB     = L - B_planck;

  const double exp1 = exp (-kabs2 * delta_2 / sin_mu);
  const double exp2 = exp (-kabs2 * delta_1 / mu);

  if (kabs1 * c < 1e-4 || kabs1 * c < 1e-4) {    // Taylor solution,
    if (kabs2 * c < 1e-4 || kabs2 * c < 1e-4) {  //Taylor solution,
      const double integrate_abs = 0.0;
      return integrate_abs;
    } else {
      const double integrate_abs = L * ((1. - exp1) * (b - a) - (c - b) - exp2 * (c - b));
      // FJ In the original version, this was never happening...
      // the return statement was missing and 0 was always returned.
      // lets throw an error and see when that happens
      fprintf (stderr,
               "We ended up in a branch that was not expected. Returning 0 while it probably was meant to be %g. Please ask a "
               "developer to look into this (line %d, function '%s', file '%s')\n",
               integrate_abs,
               ERROR_POSITION);
      return 0;
    }
  } else {
    double factor_ab;
    double factor1_bc;
    if (kabs1 * c < 1e-4) {
      factor_ab  = (1. - exp1) * L * (b - a);
      factor1_bc = -(LB) * (c - b);
    } else {
      const double exp3 = exp (-kabs1 * b / mu);
      const double exp4 = exp (-kabs1 * a / mu);
      const double exp5 = exp (-kabs1 * c / mu);
      factor_ab         = (1. - exp1) * (B_planck * (b - a) - (LB) *mu / kabs1 * (exp3 - exp4));
      factor1_bc        = -(LB) *mu / kabs1 * (exp5 - exp3);
    }
    const double factor2_bc = B_planck * (c - b);

    double factor3_bc;
    if (kabs2 * c < 1e-4)
      factor3_bc = -B_planck * exp2 * (c - b);
    else
      factor3_bc = -B_planck * mu / kabs2 * (exp (kabs2 * (c - delta_1) / mu) - exp (kabs2 * (b - delta_1) / mu));

    double factor4_bc;
    if ((kabs1 * c < 1e-4 || fabs (kabs2 - kabs1) < 1e-4))
      factor4_bc = -(LB) *exp2 * (c - b);
    else
      factor4_bc = -(LB) *mu / (kabs2 - kabs1) *
                   (exp (((kabs2 - kabs1) * c - kabs2 * delta_1) / mu) - exp (((kabs2 - kabs1) * b - kabs2 * delta_1) / mu));

    const double integrate_abs = (factor_ab + factor1_bc + factor2_bc + factor3_bc + factor4_bc);
    return integrate_abs;
  }
}

// ################################################################################
// ####################### Function - integrate_flux ###############################
// # This function integrates the contributiong emission/absorption of a grid box #
// ################################################################################
double integrate_flux (double B_planck1,
                       double B_planck2,
                       double L,
                       double kabs1,
                       double kabs2,
                       double b,
                       double c,
                       double delta_1,
                       double delta_2,
                       double mu) {

  double factor2_bc;
  double factor3_bc;

  const double sin_mu = sqrt (1. - mu * mu);
  const double LB     = L - B_planck1;

  const double exp1 = exp (-kabs2 * delta_2 / sin_mu);

  const double factor1_bc = B_planck2 * (c - b);

  if (kabs2 * c < 1e-4)
    factor2_bc = (B_planck1 - B_planck2) * exp1 * (c - b);
  else if (kabs2 * c > 300)
    factor2_bc = 0.;
  else
    factor2_bc = (B_planck1 - B_planck2) * exp1 * mu / kabs2 * (exp (kabs2 * c / mu) - exp (kabs2 * b / mu));

  if (fabs (kabs2 - kabs1) * c < 1e-4)
    factor3_bc = (LB) * (c - b) * exp (-kabs2 * delta_1 / mu);
  else if (kabs2 * c > 300)
    factor3_bc = (LB);
  else
    factor3_bc = (LB) *mu / (kabs2 - kabs1) *
                 (exp (((kabs2 - kabs1) * c - (kabs2 * delta_1)) / mu) - exp (((kabs2 - kabs1) * b - (kabs2 * delta_1)) / mu));

  const double integrate_flux = (factor1_bc + factor2_bc + factor3_bc) / (c - b);

  return integrate_flux;
}

// ################################################################################
// ####################### Function - integrate_emis ##############################
// # This function integrates the contributiong emission/absorption of a grid box #
// ################################################################################

double integrate_emis (double kabs, double a, double b, double c, double delta_1, double delta_2, double mu) {

  double integrate_emis;

  const double sin_mu = sqrt (1. - mu * mu);
  const double exp1   = exp (-kabs * delta_2 / sin_mu);

  if (kabs < 1e-2)
    integrate_emis = (b - a) * (1. - exp1) + (c - b) * (1. - exp (-kabs * delta_1 / mu));
  else
    integrate_emis =
      ((b - a) * (1. - exp1) + (c - b) - mu / kabs * (exp (kabs * (c - delta_1) / mu) - exp (kabs * (b - delta_1) / mu)));

  return integrate_emis;
}

// ################################################################################
// ####################### Function - integrate_emis ##############################
// # This function integrates the contributiong emission/absorption of a grid box #
// ################################################################################

int get_index (double a, double b) {

  int back = 0;
  int frac = 0;

  frac = a / b;

  back = a - frac * b;

  if (back == -1)
    back = b - 1;

  return back;
}

// ################################################################################
// ################## Function - interpolate emissivity ###########################
// # This function interpolates in 3D space between different emissivities ########
// ################################################################################

double pol_3d (double    tauxx,
               double    tauyy,
               double    tauzz,
               double*** eps_tab,
               double*   tau_arrx,
               double*   tau_arry,
               double*   tau_arrz,
               int       ntaux,
               int       ntauy,
               int       ntauz) {

  int ix = 0, iy = 0, iz = 0, i = 0;
  int ixx = 0, iyy = 0, izz = 0;
  //  double taux=0.0, tauy=0.0, tauz=0.0;
  double delx = 0.0, dely = 0.0, delz = 0.0;
  double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0, temp4 = 0.0, tmp1 = 0.0, tmp2 = 0.0;
  double back_3d = 0.0;

  // find indeces of asp
  // at the table limits: set to lower or upper boundary
  for (i = 0; i < ntaux; i++) {
    if (tauxx > tau_arrx[i])
      ix = i;
    else if (tauxx < tau_arrx[0]) {
      ix    = 0;
      tauxx = tau_arrx[0];
    } else if (tauxx > tau_arrx[ntaux - 1]) {
      ix    = ntaux - 1;
      tauxx = tau_arrx[ntaux - 1];
    }
  }
  for (i = 0; i < ntauy; i++) {
    if (tauyy > tau_arry[i])
      iy = i;
    else if (tauyy < tau_arry[0]) {
      iy    = 0;
      tauyy = tau_arry[0];
    } else if (tauyy > tau_arry[ntauy - 1]) {
      iy    = ntauy - 1;
      tauyy = tau_arry[ntauy - 1];
    }
  }
  for (i = 0; i < ntauz; i++) {
    if (tauzz > tau_arrz[i])
      iz = i;
    else if (tauzz < tau_arrz[0]) {
      iz    = 0;
      tauzz = tau_arrz[0];
    } else if (tauyy > tau_arrz[ntauz - 1]) {
      iz    = ntauz - 1;
      tauzz = tau_arrz[ntauz - 1];
    }
  }

  if (tauxx < tau_arrx[0] || tauzz < tau_arrz[0] || tauyy < tau_arry[0]) {
    if (tauzz < tauxx || tauzz < tauyy)
      back_3d = 1 - (exp (-tauzz));
    else if (tauxx < tauzz || tauxx < tauyy)
      back_3d = 1 - (exp (-tauxx));
    else
      back_3d = 1 - (exp (-tauyy));
  } else {
    if (ix >= ntaux - 1)
      delx = 0;
    else {
      delx = (tauxx - tau_arrx[ix]) / (tau_arrx[ix + 1] - tau_arrx[ix]);
    }
    if (iy >= ntauy - 1)
      dely = 0;
    else
      dely = (tauyy - tau_arry[iy]) / (tau_arry[iy + 1] - tau_arry[iy]);
    if (iz >= ntauz - 1)
      delz = 0;
    else
      delz = (tauzz - tau_arrz[iz]) / (tau_arrz[iz + 1] - tau_arrz[iz]);

    if (ix == ntaux - 1)
      ixx = ntaux - 2;
    else
      ixx = ix;
    if (iy == ntauy - 1)
      iyy = ntauy - 2;
    else
      iyy = iy;
    if (iz == ntauz - 1)
      izz = ntauz - 2;
    else
      izz = iz;
    temp1 = eps_tab[ix][iy][iz] * (1 - delx) + eps_tab[ixx + 1][iy][iz] * delx;
    temp2 = eps_tab[ix][iy][izz + 1] * (1 - delx) + eps_tab[ixx + 1][iy][izz + 1] * delx;
    temp3 = eps_tab[ix][iyy + 1][iz] * (1 - delx) + eps_tab[ixx + 1][iyy + 1][iz] * delx;
    temp4 = eps_tab[ix][iyy + 1][izz + 1] * (1 - delx) + eps_tab[ixx + 1][iyy + 1][izz + 1] * delx;

    tmp1 = temp1 * (1 - dely) + temp3 * dely;
    tmp2 = temp2 * (1 - dely) + temp4 * dely;

    back_3d = tmp1 * (1 - delz) + tmp2 * delz;
  }
  //correct monte carlo noise in lookup table for high optical thickness (must be 1 in the limit)
  if (back_3d > 1)
    back_3d = 1.;

  return back_3d;
}

// ################################################################################
// ################## Function - interpolate emissivity ###########################
// # This function interpolates in 3D space between different emissivities ########
// ################################################################################

double pol_2d (double var1, double var2, double** tab1, double* asp, double* tau_var, int nasp, int ntau, int flag) {

  int    ix = 0, iy = 0, i = 0;
  double tmp1    = 0.0;
  double back_2d = 0.0;

  // find indeces of var1 and var2
  // at the table limits: set to lower or upper boundary
  for (i = 0; i < nasp; i++) {
    if (var1 >= asp[i] && var1 < asp[nasp - 1] && var1 > asp[0])
      ix = i;
    else if (var1 <= asp[0])
      ix = 0;
    else if (var1 >= asp[nasp - 1])
      ix = nasp - 1;
  }

  for (i = 0; i < ntau; i++) {
    if (var2 >= tau_var[i] && var2 < tau_var[ntau - 1] && var2 > tau_var[0])
      iy = i;
    else if (var2 <= tau_var[0])
      iy = 0;
    else if (var2 >= tau_var[ntau - 1])
      iy = ntau - 1;
  }

  //if optical depth is lower than lookup table, set to 1-exp(-tau) and exit program
  if (flag == 2) {
    if (var1 < asp[0] || var2 < tau_var[0]) {
      if (var2 < var1)
        back_2d = 1 - (exp (-var1));
      else
        back_2d = 1 - (exp (-var1));
    }
  }
  if (var1 < asp[nasp - 1] && var1 >= asp[0]) {
    if (var2 < tau_var[ntau - 1] && var2 >= tau_var[0]) {
      const double f1 = (asp[ix + 1] - var1) / (asp[ix + 1] - asp[ix]) * tab1[ix][iy] +
                        (var1 - asp[ix]) / (asp[ix + 1] - asp[ix]) * tab1[ix + 1][iy];
      const double f2 = (asp[ix + 1] - var1) / (asp[ix + 1] - asp[ix]) * tab1[ix][iy + 1] +
                        (var1 - asp[ix]) / (asp[ix + 1] - asp[ix]) * tab1[ix + 1][iy + 1];
      back_2d = (tau_var[iy + 1] - var2) / (tau_var[iy + 1] - tau_var[iy]) * f1 +
                (var2 - tau_var[iy]) / (tau_var[iy + 1] - tau_var[iy]) * f2;
    } else if (var2 >= tau_var[ntau - 1]) {
      back_2d = (asp[ix + 1] - var1) / (asp[ix + 1] - asp[ix]) * tab1[ix][iy] +
                (var1 - asp[ix]) / (asp[ix + 1] - asp[ix]) * tab1[ix + 1][iy];
    } else if (var2 < tau_var[0]) {
      tmp1    = tab1[ix][iy];
      back_2d = tmp1;
    }
  } else if (var1 >= asp[nasp - 1]) {
    if (var2 < tau_var[ntau - 1] && var2 > tau_var[0]) {
      tmp1 = (tau_var[iy + 1] - var2) / (tau_var[iy + 1] - tau_var[iy]) * tab1[ix][iy] +
             (var2 - tau_var[iy]) / (tau_var[iy + 1] - tau_var[iy]) * tab1[ix][iy + 1];
      back_2d = tmp1;
    } else if (var2 >= tau_var[ntau - 1]) {
      tmp1    = tab1[nasp - 1][ntau - 1];
      back_2d = tmp1;
    } else {
      fprintf (stderr, " %f\n", tab1[0][0]);
      tmp1    = tab1[0][0];
      back_2d = tmp1;
    }
  }

  return back_2d;
}
