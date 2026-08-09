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

#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>

#include "allocnd.h"
#include "ancillary.h"
#include "ascii.h"
#include "cdisort.h"
#include "cloud.h"
#include "common_math.h"
#include "Corefinder.h"
#include "c_tzs.h"
#include "dynamic_tenstream.h"
#include "dynamic_twostream.h"
#include "errors.h"
#include "fortran_and_c.h"
#include "LinArray.h"
#include "molecular.h"
#include "numeric.h"
#include "redistribute.h"
#include "rodents.h"
#include "solve_rte.h"
#include "sslidar.h"
#include "twomaxrnd.h"
#include "twostrebe.h"
#include "schwarzschild.h"
#include "uvspec.h"
#include "uvspecrandom.h"
#include "wcloud3d.h"

#if HAVE_TWOMAXRND3C
#include "twomaxrnd3C.h"
#endif
#if HAVE_MYSTIC
#include "mystic.h"
#endif
#if HAVE_TIPA
#include "tipa.h"
#endif
#if HAVE_SOS
#include "sos.h"
#endif

#include "f77-uscore.h"
#include "solver.h"

/* Definitions for numerical recipes functions */
#define SIGN(a, b) ((b) >= 0.0 ? fabs (a) : -fabs (a))

typedef struct {
  char header[127]; /* A 127- (or less) character header for prints         */

  float accur; /* Convergence criterion for azimuthal series.          */

  float fbeam; /* Intensity of incident parallel beam at top boundary. */
               /*  [same units as PLKAVG (default W/sq m) if thermal   */
               /*  sources active, otherwise arbitrary units].         */

  float fisot; /* Intensity of top-boundary isotropic illumination.    */
               /*  [same units as PLKAVG (default W/sq m) if thermal   */
               /*  sources active, otherwise arbitrary units].         */
               /*  Corresponding incident flux is pi (3.14159...)      */
               /*  times FISOT.                                        */

  float* hl; /* K = 0 to NSTR.  Coefficients in Legendre-            */
             /* polynomial expansion of bottom bidirectional         */
             /* reflectivity                                         */

  int planck; /* TRUE,  include thermal emission                      */
              /* FALSE, ignore thermal emission (saves computer time) */

  float btemp; /* Temperature of bottom boundary (K)                   */
               /* (bottom emissivity is calculated from ALBEDO or HL,  */
               /* so it need not be specified).                        */
               /* Needed only if PLANK is TRUE.                        */

  float ttemp; /* Temperature of top boundary (K)                      */
               /* Needed only if PLANK is TRUE.                        */

  float temis; /* Emissivity of top boundary.                          */
               /* Needed only if PLANK is TRUE.                        */

  float* utau;
  float  umu0;

  int ierror_d[47];
  int ierror_s[33];
  int ierror_t[22];

  int prndis[7];
  int prndis2[5];
  int prntwo[2];

  int ibcnd; /* 0 : General case.                                    */
             /* 1 : Return only albedo and transmissivity of the     */
             /*     entire medium vs. incident beam angle.           */

  int lamber; /* TRUE, isotropically reflecting bottom boundary.      */
              /* FALSE, bidirectionally reflecting bottom boundary.   */

  int onlyfl; /* TRUE, return fluxes, flux divergences, and mean      */
              /*       intensities.                                   */
              /* FALSE, return fluxes, flux divergences, mean         */
              /*        intensities, azimuthally averaged intensities */
              /*        (at the user angles) AND intensities.         */

  int quiet;
  int usrang;
  int usrtau;

  /* sdisort-specific variables */
  int nil;
  int newgeo;
  int spher;

  /* qdisort-specific variables */
  int       gsrc;  /* Flag for general source for qdisort */
  double*** qsrc;  /* The general source, in FORTRAN: REAL*4 qsrc( MXCMU, 0:MXULV, MXCMU )*/
                   /* At computational angles                                             */
  double*** qsrcu; /* The general source, in FORTRAN: REAL*4 qsrc( MXUMU, 0:MXULV, MXCMU )*/
                   /* At user angles                                                      */

} rte_input;

typedef struct {
  float*   dfdt;
  float*   flup;
  float*   rfldir;
  float*   rfldn;
  float*   uavgso;
  float*   uavgdn;
  float*   uavgup;
  float*   uavg;
  float*   heat;
  float*   emis;
  float*   w_zout;
  float**  u0u;
  float*** uu;
  float*** uum; /* Fourier components of intensities, returned by qdisort */
  float*   sslidar_nphot;
  float*   sslidar_nphot_q;
  float*   sslidar_ratio;

  int         alloc3D; /* flag if 3D was allocated */
  float***    rfldir3d;
  float***    rfldn3d;
  float***    flup3d;
  float*****  fl3d_is; /* importance sampling */
  float***    uavgso3d;
  float***    uavgdn3d;
  float***    uavgup3d;
  float***    abs3d;
  float***    absback3d;
  float****   radiance3d;
  float****** radiance3d_is; /*importance sampling*/
  float****** jacobian;

  /* corresponding variances */
  float***  rfldir3d_var;
  float***  rfldn3d_var;
  float***  flup3d_var;
  float***  uavgso3d_var;
  float***  uavgdn3d_var;
  float***  uavgup3d_var;
  float***  abs3d_var;
  float***  absback3d_var;
  float**** radiance3d_var;

  float* albmed;
  float* trnmed;

  // triangular surface output
  struct t_triangle_radiation_field* triangle_results;
} rte_output;

typedef struct {
  float* tauw;
  float* taui;

  float* g1d;
  float* g2d;
  float* fd;
  float* g1i;
  float* g2i;
  float* fi;

  float* ssaw;
  float* ssai;
} save_optprop;

typedef struct {
  double**   dtauc;
  double*    fbeam;
  double**** uum;
  double **  uavgso, **uavgdn, **uavgup;
  double **  rfldir, **rfldn, **flup;
  double***  u0u;
  double**** uu;
} raman_qsrc_components;

/* prototypes of internal functions */
static int reverse_profiles (input_struct input, output_struct* output);

static int calloc_rte_output (input_struct input,
                              int          nzout,
                              int          output3D,
                              int          Nxcld,
                              int          Nycld,
                              int          Nzcld,
                              int          Nxsample,
                              int          Nysample,
                              int          Ncsample,
                              int          Nlambda,
                              int          Nlyr,
                              int*         threed,
                              const size_t N_triangles,
                              rte_output** result);

static int reset_rte_output (rte_output** rte,
                             input_struct input,
                             int          nzout,
                             int          Nxcld,
                             int          Nycld,
                             int          Nzcld,
                             int          Nxsample,
                             int          Nysample,
                             int          Ncsample,
                             int          Nlambda,
                             int          Nlyr,
                             int*         threed,
                             const size_t N_triangles);

static int free_rte_output (rte_output*  result,
                            input_struct input,
                            int          nzout,
                            int          Nxcld,
                            int          Nycld,
                            int          Nzcld,
                            int          Nxsample,
                            int          Nysample,
                            int          Ncsample,
                            int          Nlyr,
                            int          Nlambda,
                            int*         threed);

static int add_rte_output (rte_output*        rte,
                           const rte_output*  add,
                           const double       factor,
                           const double*      factor_spectral,
                           const input_struct input,
                           const int          nzout,
                           const int          Nxcld,
                           const int          Nycld,
                           const int          Nzcld,
                           const int          Nc,
                           const int          Nlyr,
                           const int          Nlambda,
                           const int*         threed,
                           const int          islower,
                           const int          isupper,
                           const int          jslower,
                           const int          jsupper,
                           const int          isstep,
                           const int          jsstep);

static int init_rte_input (rte_input* rte, input_struct input, output_struct* output);

static int setup_and_call_solver (input_struct           input,
                                  output_struct*         output,
                                  rte_input*             rte_in,
                                  rte_output*            rte_out,
                                  raman_qsrc_components* raman_qsrc_components,
                                  int                    iv,
                                  int                    ib,
                                  int                    ir,
                                  int*                   threed,
                                  int                    mc_loaddata);

static int call_solver (input_struct           input,
                        output_struct*         output,
                        int                    rte_solver,
                        rte_input*             rte_in,
                        raman_qsrc_components* raman_qsrc_components,
                        int                    iv,
                        int                    ib,
                        int                    ir,
                        rte_output*            rte_out,
                        int*                   threed,
                        int                    mc_loaddata,
                        int                    verbose);

static int calc_spectral_heating (input_struct   input,
                                  output_struct* output,
                                  float*         dz,
                                  double*        rho_mass_zout,
                                  float*         k_abs,
                                  float*         k_abs_layer,
                                  rte_output*    rte_out,
                                  float*         heat,
                                  float*         emis,
                                  float*         w_zout,
                                  int            iv);

static float*** calloc_abs3d (const int Nx, const int Ny, const int Nz, const int* threed);

static void free_abs3d (float*** abs3d, const int Nx, const int Ny, const int Nz, const int* threed);

static double get_unit_factor (input_struct input, output_struct* output, int iv);

static int
generate_effective_cloud (input_struct input, output_struct* output, save_optprop* save_cloud, int iv, int iq, int verbose);

static int set_raman_source (double***              qsrc,
                             double***              qsrcu,
                             int                    maxphi,
                             int                    nlev,
                             int                    nzout,
                             int                    nstr,
                             int                    n_shifts,
                             float                  wanted_wl,
                             double*                wl_shifts,
                             float                  umu0,
                             float*                 zd,
                             float*                 zout,
                             float                  zenang,
                             float                  fbeam,
                             float                  radius,
                             float*                 dens,
                             double**               crs_RL,
                             double**               crs_RG,
                             float*                 ssalb,
                             int                    numu,
                             float*                 umu,
                             int                    usrang,
                             int*                   cmuind,
                             float***               pmom,
                             raman_qsrc_components* raman_qsrc_components,
                             int*                   zout_comp_index,
                             float                  altitude,
                             int                    last,
                             int                    verbose);

static save_optprop* calloc_save_optprop (int Nlev);

raman_qsrc_components* calloc_raman_qsrc_components (int raman_fast, int nzout, int maxumu, int nphi, int nstr, int nlambda_shift);

static void free_raman_qsrc_components (raman_qsrc_components* result, int raman_fast, int nzout, int maxumu, int nphi, int nstr);

void F77_FUNC (swde, SWDE) (float* g_scaled,
                            float* pref,
                            float* prmuz,
                            float* tau,
                            float* ssa_gas,
                            float* pre1,
                            float* pre2,
                            float* ptr1,
                            float* ptr2);

float zbrent_taueff (float mu_eff,
                     float g_scaled,
                     float ssa_gas,
                     float transmission_cloud,
                     float transmission_layer,
                     float x1,
                     float x2,
                     float tol);

static int find_veiled_core (const input_struct* input, const int iv_count, output_struct* output);

/*********************************************************************/
/* Main function. Loop over wavelengths or wavelength bands,         */
/* correlated-k quadrature points, and independent pixels.           */
/*********************************************************************/

int solve_rte (input_struct input, output_struct* output) {
  static int first  = TRUE;
  int        status = 0, add = 0;

  int       nr = 0, irs = 0;
  int       lower_wl_id = 0, upper_wl_id = 0, lower_iq_id = 0, upper_iq_id = 0, ivr = 0;
  int       nlambda       = 0;
  float*    dz            = NULL;
  double    weight        = 1;
  double*   rho_mass_zout = NULL;
  float*    k_abs         = NULL;
  float*    k_abs_layer   = NULL;
  float*    k_abs_outband = NULL;
  int       mc_loaddata   = 1;
  double    weight2       = 0;
  double    unit_factor   = 0;
  double    ffactor = 0.0, rfactor = 0.0, hfactor = 1.0;
  double**  u0u_raman    = NULL;
  double*   uavgso_raman = NULL;
  double*   uavgdn_raman = NULL;
  double*   uavgup_raman = NULL;
  double*   rfldir_raman = NULL;
  double*   rfldn_raman  = NULL;
  double*   flup_raman   = NULL;
  double*** uu_raman     = NULL; /* The radiance at user angles for the general source for qdisort */
                                 /* Only used if raman scattering is on                            */

  /* float heating_rate_emission = 0.0; */

  rte_input*  rte_in      = NULL;
  rte_output* rte_out     = NULL;
  rte_output* rte_outband = NULL;

  save_optprop*          save_cloud            = NULL;
  raman_qsrc_components* raman_qsrc_components = NULL;

  double* weight_spectral = NULL;

  //FIX 3DAbs not initialized when aerosol setup not done, should be redundant now
  // output->mc.alis.Nc = 1;

  if (first) {
    first = FALSE;

    /* this is not the right place to do this! This should be done somewhere else! Please clean up! BCA */
    if (input.ipa3d) {
      /* ipa3d and tipa were configured and tested only with 3D-cloudfiles containing lwc and reff, thus check: */
      for (int isp = 0; isp < input.n_caoth; isp++)
        if (!(output->caoth3d[isp].cldproperties == 3 || output->caoth3d[isp].cldproperties == 4))
          CHKERROUT (-1, "ipa3d/tipa does not work with cldproperties flag different from 3");

      /* finally, we add normalization to the number of pixels to ipaweight, BM 3.7.2020 */
      /* ipaweight[] either considers cloudcover or normalization for ipa3d              */
      for (int ipa = 0; ipa < output->nipa; ipa++)
        output->ipaweight[ipa] /= ((double) output->niipa * (double) output->njipa);
    }

    rte_in = calloc (1, sizeof (rte_input));

    status = calloc_rte_output (input,
                                output->atm.nzout,
                                input.output3D,
                                output->atm.Nxcld,
                                output->atm.Nycld,
                                output->atm.Nzcld,
                                output->mc.sample.Nx,
                                output->mc.sample.Ny,
                                output->mc.alis.Nc,
                                output->mc.alis.nlambda_abs,
                                output->atm.nlev - 1,
                                output->atm.threed,
                                output->mc.triangular_surface.N_triangles,
                                &rte_out);
    CHKERR (status);

    // in case of IPA3D we don't need 3D information for rte_outband
    int output3D = (input.ipa3d ? 0 : input.output3D);

    status = calloc_rte_output (input,
                                output->atm.nzout,
                                output3D,
                                output->atm.Nxcld,
                                output->atm.Nycld,
                                output->atm.Nzcld,
                                output->mc.sample.Nx,
                                output->mc.sample.Ny,
                                output->mc.alis.Nc,
                                output->mc.alis.nlambda_abs,
                                output->atm.nlev - 1,
                                output->atm.threed,
                                output->mc.triangular_surface.N_triangles,
                                &rte_outband);
    CHKERR (status);

    save_cloud = calloc_save_optprop (output->atm.nlev);

    /* initialize RTE input */
    status = init_rte_input (rte_in, input, output);
    CHKERR (status);

    /* allocate memory for the output result structures */
    status = setup_result (input, output, &dz, &rho_mass_zout, &k_abs_layer, &k_abs, &k_abs_outband);
    CHKERR (status);

    if (input.raman) {
      nr = 2;

      status = ASCII_calloc_double_3D (&uu_raman, output->atm.nzout, input.rte.nphi, input.rte.numu);
      CHKERR (status);
      status = ASCII_calloc_double (&u0u_raman, output->atm.nzout, input.rte.numu);
      CHKERR (status);
      uavgso_raman = calloc (output->atm.nzout, sizeof (double));
      uavgdn_raman = calloc (output->atm.nzout, sizeof (double));
      uavgup_raman = calloc (output->atm.nzout, sizeof (double));
      rfldir_raman = calloc (output->atm.nzout, sizeof (double));
      rfldn_raman  = calloc (output->atm.nzout, sizeof (double));
      flup_raman   = calloc (output->atm.nzout, sizeof (double));

      /* Number of wavelength shifts considered is independent of primary wavelength, */
      /* hence use output->atm.nq_t[0] below                                          */
      if (input.raman_fast) {
        nlambda = output->wl.nlambda_r;
      } else {
        nlambda = output->crs.number_of_ramanwavelengths;
      }
      raman_qsrc_components = calloc_raman_qsrc_components (input.raman_fast,
                                                            output->atm.nzout,
                                                            input.rte.maxumu,
                                                            input.rte.nphi,
                                                            input.rte.nstr,
                                                            nlambda);
    } else
      nr = 1;
  }

  if (input.raman) {
    /* For Raman scattering only include wavelengths that the user asked. */
    /* Internally we have to include more wavelengths to account for      */
    /* Raman scattered radiation, see loop over quadrature points below.  */

    lower_wl_id = output->wl.raman_start_id;
    upper_wl_id = output->wl.raman_end_id;
  } else {
    lower_wl_id = output->wl.iv_lower;
    upper_wl_id = output->wl.iv_upper;
  }

  /* concentration importance sampling */
  if (input.rte.mc.concentration_is) {
    weight_spectral    = calloc (1, sizeof (double));
    weight_spectral[0] = 1.0;
  }

  if (input.rte.mc.spectral_is) {

    weight_spectral = calloc (output->mc.alis.nlambda_abs, sizeof (double));
    /* Take wavelength in center of spectrum if not specified explicitly */
    /* FIXCE should also check whether absorption is not too high here   */
    if (input.rte.mc.spectral_is_wvl[0] == 0.) {
      lower_wl_id                    = (int) (0.5 * ((float) output->wl.iv_lower + (float) output->wl.iv_upper));
      output->mc.alis.nlambda_ref    = 1;
      output->mc.alis.ilambda_ref    = calloc (1, sizeof (int));
      output->mc.alis.ilambda_ref[0] = lower_wl_id;
    }

    else if (input.rte.mc.spectral_is_wvl[0] > 0.) {
      /* Find wavelength index for specified wavelength */
      lower_wl_id = 0;
      for (int iv_alis = 0; iv_alis < output->mc.alis.nlambda_abs; iv_alis++) {
        if (output->mc.alis.lambda[iv_alis] > input.rte.mc.spectral_is_wvl[0]) {
          lower_wl_id                    = iv_alis - 1;
          output->mc.alis.nlambda_ref    = 1;
          output->mc.alis.ilambda_ref    = calloc (1, sizeof (int));
          output->mc.alis.ilambda_ref[0] = lower_wl_id;
          break;
        }
      }
    } else {
      /* several calc wvls do not work so far */
      lower_wl_id                 = (int) (0.5 * ((float) output->wl.iv_lower + (float) output->wl.iv_upper));
      output->mc.alis.nlambda_ref = input.rte.mc.spectral_is_nwvl;
      output->mc.alis.ilambda_ref = calloc (output->mc.alis.nlambda_ref, sizeof (int));
      for (int iv_alis_ref = 0; iv_alis_ref < output->mc.alis.nlambda_ref; iv_alis_ref++) {
        for (int iv_alis = 0; iv_alis < output->mc.alis.nlambda_abs; iv_alis++) {
          if (output->mc.alis.lambda[iv_alis] > input.rte.mc.spectral_is_wvl[iv_alis_ref]) {
            output->mc.alis.ilambda_ref[iv_alis_ref] = iv_alis;
          }
        }
      }
    }

    upper_wl_id = lower_wl_id;
    if (!input.quiet) {
      fprintf (stderr, "... ALIS calculation wavelength: %g nm \n", output->mc.alis.lambda[lower_wl_id]);
      if (output->mc.alis.nlambda_ref > 1)
        for (int iv_alis_ref = 1; iv_alis_ref < output->mc.alis.nlambda_ref; iv_alis_ref++)
          fprintf (stderr,
                   "... helper ALIS wavelength: %g nm \n",
                   output->mc.alis.lambda[output->mc.alis.ilambda_ref[iv_alis_ref]]);
    }
  }

#if HAVE_TIPA
  /* ulrike: TIPA DIR. The "tilted cloud matrix" is used only for the
                       calculation of the DIRECT radiation Tilting for
                       every z-level is done here (outside the loop
                       over the wavelength)  */
  if (input.tipa == TIPA_DIR || input.rte.mc.tipa == TIPA_DIR) {
    for (int isp = 0; isp < input.n_caoth; isp++)
      if (input.caoth[isp].source == CAOTH_FROM_3D) {
        if (!input.quiet)
          fprintf (stderr, " ... performing the tilting for tipa dir (water clouds)\n");
        status = tipa_dirtilt (&(output->caoth3d[isp]),
                               output->atm,
                               output->alt,
                               &(output->caoth[isp].tipa),
                               input.tipa,
                               input.atm.sza,
                               input.atm.phi0,
                               lower_wl_id,
                               upper_wl_id,
                               input.rte.mc.tipa);
        if (status)
          return fct_err_out (status, "tipa_dirtilt", ERROR_POSITION);
      }
  }
#endif

  /************************/
  /* loop over wavelength */
  /************************/
  int iv_count = 0;
  for (int iv = lower_wl_id; iv <= upper_wl_id; iv++) {

    // ipa3d hack, BM 28.1.2023: We don't have a 1D polarisation solver, therefore we set nphamat to 1;
    // although polarization information might be provided in the optical property file, we don't want them
    // in the CAOTH to save memory; also, for the required freeing in the (iipa, jipa) loop it would cause segfaults
    for (int isp = 0; isp < input.n_caoth; isp++) {
      if (input.caoth[isp].ipa || (input.ipa3d && input.caoth[isp].source == CAOTH_FROM_3D)) {
        output->caoth[isp].ssprop.nphamat  = 1;
        output->caoth[isp].optprop.nphamat = 1;
      }
    }

    /***********************************************************/
    /* iterate over wavelengths, required for raman scattering */
    /***********************************************************/

    irs = 0;
    if (input.raman_fast && iv > lower_wl_id)
      irs = 1; /* AK20110407: All the elastic wavelengths are done for  */
               /* iv=lower_wl_id and stored. Thus only need to do       */
               /* the inelastic part for remaining wavelengths.         */

    for (int ir = irs; ir < nr; ir++) {

      /* solar zenith angle at this wavelength */
      rte_in->umu0 = cos (output->atm.sza_r[iv] * M_PI / 180.0);

      /* calculate 3D caoth properties for this wavelength */

      if (input.rte.solver == SOLVER_MONTECARLO || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) {
        for (int isp = 0; isp < input.n_caoth; isp++) {
          if (strcmp (input.caoth[isp].name, "molecular_3d") != 0) {
            status = convert_caoth_mystic (input, input.caoth[isp], output, &(output->caoth[isp]), &(output->caoth3d[isp]), iv);
            CHKERR (status);
          }
        }
      }

      if (input.rte.solver == SOLVER_MONTECARLO || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) {
        status = find_veiled_core (&input, iv_count, output);
        CHKERR (status);
      }

      /********************************/
      /* loop over independent pixels */
      /********************************/
      for (int ipa = 0; ipa < output->nipa; ipa++) {
        for (int iipa = 0; iipa < output->niipa; iipa++) {
          for (int jipa = 0; jipa < output->njipa; jipa++) {

            if (input.verbose && (output->niipa > 1 || output->njipa > 1))
              fprintf (stderr, " ... ipa loop over iipa=%d, jipa=%d\n", iipa, jipa);

            if (input.verbose && output->nipa == 1)
              fprintf (stderr,
                       "\n\n*** wavelength: iv = %d, %f nm, albedo = %f \n",
                       iv,
                       output->wl.lambda_r[iv],
                       output->alb.albedo_r[iv]);

            if (input.verbose && output->nipa > 1)
              fprintf (stderr,
                       "\n\n*** wavelength: iv = %d, %f nm, looking at column %d, albedo = %f\n",
                       iv,
                       output->wl.lambda_r[iv],
                       ipa,
                       output->alb.albedo_r[iv]);

            /* need to load data during first call of mystic() for each band/wavelength */
            mc_loaddata = 1;

            /* reset band integral */
            reset_rte_output (&rte_outband,
                              input,
                              output->atm.nzout,
                              output->atm.Nxcld,
                              output->atm.Nycld,
                              output->atm.Nzcld,
                              output->mc.sample.Nx,
                              output->mc.sample.Ny,
                              output->mc.alis.Nc,
                              output->mc.alis.nlambda_abs,
                              output->atm.nlev - 1,
                              output->atm.threed,
                              output->mc.triangular_surface.N_triangles);

            /* copy pixel number ipa to 1D caoth data, ulrike: allow ipa3d (and tipa dir) for caoth */
            for (int isp = 0; isp < input.n_caoth; isp++) {
              if (input.caoth[isp].ipa || (input.ipa3d && input.caoth[isp].source == CAOTH_FROM_3D)) {

                // BM, 8.9.2023, old code (manual free by UW) replaced by free_optprop() plus calloc_optprop()
                // Why do we need that? It seems to run without, but the computational cost is very small anyway.
                // BM 24.3.2024  commented out; obviously not needed; profiles of tau etc are zeroed anyway in caoth_prop_switch()
                /*
                if (input.ipa3d) {
                  free_optprop (&(output->caoth[isp].optprop));
                  calloc_optprop (&(output->caoth[isp].optprop), output->wl.nlambda_r, output->atm.nlev, 1);
                }
		*/

                /* copy everything except the single scattering properties */
                if (input.caoth[isp].ipa) {
                  status = cp_caoth_out (&(output->caoth[isp]), output->caoth_ipa[isp][ipa], 0, 0, input.quiet);
                  if (status)
                    return fct_err_out (status, "cp_cld_out", ERROR_POSITION);
                }

                if (input.ipa3d) {
                  if (input.verbose)
                    fprintf (stderr, " ... copying 3d to 1d for %s, pixel %d %d\n", output->caoth[isp].fullname, iipa, jipa);

                  if (output->caoth3d[isp].nthreed > 0) {  // only if we actually have 3D clouds, BM 30.7.2023
                    status = cp_caoth3d_out (&(output->caoth[isp]), output->caoth3d[isp], input.quiet, iipa, jipa);
                    CHKERR (status);
                  }

                  /* copy cloud fraction, if available */
                  if (output->caoth3d[isp].cldproperties == CLD_LWCREFFCF) {
                    if (input.verbose)
                      fprintf (stderr, " ... copying cloud fraction from caoth %d\n", isp);

                    /* need to allocate memory for cloud fraction profile? */
                    if (output->cf.nlev == 0) {
                      if (!input.quiet)
                        fprintf (stderr, " ... allocating memory for cloud fraction profile\n");

                      output->cf.nlev = output->atm.nlev - 1;
                      output->cf.zd   = calloc (output->cf.nlev, sizeof (float));
                      output->cf.cf   = calloc (output->cf.nlev, sizeof (float));

                      for (int lu = 0; lu < output->cf.nlev; lu++)
                        output->cf.zd[lu] = output->atm.zd[lu];
                    }

                    status = cp_caoth3d_cf_out (output->cf.cf, output->caoth3d[isp], input.quiet, iipa, jipa);
                    CHKERR (status);
                  }
                }  // end if (input.ipa3d)

#if HAVE_TIPA
                if (input.tipa == TIPA_DIR) { /* ulrike: calculate tilted dtau for caoth */
                  status = tipa_calcdtau (input.caoth[isp],
                                          output->caoth3d[isp],
                                          iipa,
                                          jipa,
                                          iv,
                                          input,
                                          output->wl,
                                          &(output->caoth[isp]),
                                          &(output->caoth[isp].tipa));
                  CHKERR (status);
                }
#endif

                /* calculate optical properties for the caoth properties */
                /* specified in the input-file (ulrike)                  */
                status = caoth_prop_switch (input, input.caoth[isp], output->wl, iv, &(output->caoth[isp]));
                CHKERR (status);

                /* overwrite these properties with user-defined optical thickness, ssa, etc */
                status = apply_user_defined_properties_to_caoth (input.caoth[isp],
                                                                 output->wl.nlambda_r,
                                                                 output->wl.lambda_r,
                                                                 output->alt.altitude,
                                                                 &(output->caoth[isp]));
                CHKERR (status);
              } /* end of "if (input.caoth[isp].ipa || input.ipa3d)" */
            }   /* end loop isp */

            if (input.rte.solver == SOLVER_TWOMAXRND && output->cf.nlev == 0) {
              CHKERROUT (
                -1,
                "rte_solver twomaxrnd makes only sense with cloud_fraction_file or cloud fraction defined in 3D cloud file.");
            }

            if (input.ipa) {
              /* copy cloud fraction structure, if needed */
              switch (input.cloud_overlap) {
              case CLOUD_OVERLAP_MAX:
              case CLOUD_OVERLAP_MAXRAND: /* target */ /* source */ /* alloc */
                if (input.rte.solver != SOLVER_TWOMAXRND && input.rte.solver != SOLVER_TWOMAXRND3C &&
                    input.rte.solver != SOLVER_DYNAMIC_TWOSTREAM) {
                  status = copy_cloud_fraction (&(output->cf), output->cfipa[ipa], FALSE); /* in cloud.c */
                  if (status != 0) {
                    fprintf (stderr, "Error %d copying output->cfipa[ipa] to output->cf\n", status);
                    return status;
                  }
                }
                break;
              case CLOUD_OVERLAP_RAND:
              case CLOUD_OVERLAP_OFF:
                /* nothing to do here */
                break;
              default:
                fprintf (stderr,
                         "Error, unknown cloud_overlap assumption %d. (line %d, function %s in %s)\n",
                         input.cloud_overlap,
                         __LINE__,
                         __func__,
                         __FILE__);
                return -1;
              }
            }

            /* IPA molecular absorption and aerosols */

            /* these lines also optimise the iq-loop for corr-k schemes, also if there is no ipa */
            /* CE: with spectral importance sampling number of calculations always corresponds   */
            /* to maximum number of bands                                                        */
            if (input.ck_scheme == CK_LOWTRAN && !input.rte.mc.spectral_is)
              output->atm.nq_r[iv] = output->crs_ck.profile[0][iv].ngauss;

            /* If only one subband is used in the LOWTRAN parameterization, */
            /* the photon weights of the three subbands are added; this is  */
            /* necessary because the number of subbands changes with        */
            /* concentration and is therefore not known beforehand.         */
            /* The correct use of mc_photons_file for LOWTRAN is then       */
            /* to always distribute the photons over three subbands;        */
            /* uvspec decides automatically if only one is needed           */

            if (input.ck_scheme == CK_LOWTRAN)
              if (output->atm.nq_r[iv] == 1)
                for (int iq = 1; iq < LOWTRAN_MAXINT; iq++)
                  output->mc_photons_r[iv][0] += output->mc_photons_r[iv][iq];

            if (input.raman) {
              if (ir == 0) {
                if (input.raman_fast) {
                  lower_iq_id          = output->wl.iv_lower;
                  output->atm.nq_r[iv] = output->wl.iv_upper;
                } else {
                  output->atm.nq_r[iv] = output->crs.number_of_ramanwavelengths;
                }
                upper_iq_id = output->atm.nq_r[iv];
              } else if (ir == 1) {
                output->atm.nq_r[iv] = 1;
                lower_iq_id          = 0;
                upper_iq_id          = output->atm.nq_r[iv];
              }
            } else {
              lower_iq_id = 0;
              upper_iq_id = output->atm.nq_r[iv];
            }

            /******************************************/
            /* loop over quadrature points (subbands) */
            /******************************************/
            for (int iq = lower_iq_id; iq < upper_iq_id; iq++) {

              if (input.verbose && (input.ck_scheme != CK_CRS))
                fprintf (stderr,
                         "\n*** wavelength: iv = %d, %f nm, looking at column %d, quadrature point nr %d, albedo = %f\n",
                         iv,
                         output->wl.lambda_r[iv],
                         ipa,
                         iq,
                         output->alb.albedo_r[iv]);

              if (input.verbose && (input.ck_scheme == CK_RAMAN))
                fprintf (stderr,
                         "\n*** wavelength: iv = %d, %f nm, looking at column %d, quadrature point nr %d, wvl = %f\n",
                         iv,
                         output->wl.lambda_r[iv],
                         ipa,
                         iq,
                         output->wl.lambda_r[iq]);

              /* set number of photons for this band */
              if (input.rte.solver == SOLVER_MONTECARLO && !input.rte.mc.spectral_is)
                output->mc_photons = (long int) (output->mc_photons_r[iv][iq] * (double) input.rte.mc.photons + 0.5);

              /* run at least MIN_MCPHOTONS for each band */
              /* no need to increase for backward direct because backward direct is (nearly) exact */
              if (output->mc.sample.backward != MCBACKWARD_EDIR && output->mc.sample.backward != MCBACKWARD_FDIR) {
                if (input.rte.mc.minphotons) { /* set by user */
                  if (output->mc_photons < (long int) input.rte.mc.minphotons)
                    output->mc_photons = input.rte.mc.minphotons;
                } else { /* default */
                  if (output->mc_photons < MIN_MCPHOTONS)
                    output->mc_photons = MIN_MCPHOTONS;
                }
              } else { /* however, we need at least one photon for direct */
                if (output->mc_photons < 1)
                  output->mc_photons = 1;
              }
              /* For spectral importance sampling, only one wavelength is */
              /* calculated, therefore no distribution of photons required.  */
              if (input.rte.mc.spectral_is)
                output->mc_photons = input.rte.mc.photons;

              if (input.verbose)
                fprintf (stderr, " ... ck weight %9.7f\n", output->atm.wght_r[iv][iq]);

              /* if the level number of cloud fraction data is more than 0, than ... */
              if (input.cloud_overlap != CLOUD_OVERLAP_OFF && input.rte.solver != SOLVER_TWOMAXRND &&
                  input.rte.solver != SOLVER_TWOMAXRND3C && input.rte.solver != SOLVER_DYNAMIC_TWOSTREAM) {

                /* Save optical properties for wavelength iv. This is necessary because averaged optical properties 
	       are calulated for each subband and put into output->wc.optprop.... */
                if (iq == 0) {

                  for (int lc = 0; lc < output->atm.nlev - 1; lc++) {
                    if (input.i_wc != -1) {
                      /* optical depth */
                      save_cloud->tauw[lc] = output->caoth[input.i_wc].optprop.dtau[iv][lc];
                      /* asymmetry parameter */
                      save_cloud->g1d[lc] = output->caoth[input.i_wc].optprop.g1[iv][lc];
                      save_cloud->g2d[lc] = output->caoth[input.i_wc].optprop.g2[iv][lc];
                      save_cloud->fd[lc]  = output->caoth[input.i_wc].optprop.ff[iv][lc];
                      /* single scattering albedo */
                      save_cloud->ssaw[lc] = output->caoth[input.i_wc].optprop.ssa[iv][lc];
                    } else {
                      save_cloud->tauw[lc] = 0.0;
                      save_cloud->g1d[lc]  = 0.0;
                      save_cloud->g2d[lc]  = 0.0;
                      save_cloud->fd[lc]   = 0.0;
                      save_cloud->ssaw[lc] = 0.0;
                    }

                    if (input.i_ic != -1) {
                      /* optical depth */
                      save_cloud->taui[lc] = output->caoth[input.i_ic].optprop.dtau[iv][lc];
                      /* asymmetry parameter */
                      save_cloud->g1i[lc] = output->caoth[input.i_ic].optprop.g1[iv][lc];
                      save_cloud->g2i[lc] = output->caoth[input.i_ic].optprop.g2[iv][lc];
                      save_cloud->fi[lc]  = output->caoth[input.i_ic].optprop.ff[iv][lc];
                      /* single scattering albedo */
                      save_cloud->ssai[lc] = output->caoth[input.i_ic].optprop.ssa[iv][lc];
                    } else {
                      save_cloud->taui[lc] = 0.0;
                      save_cloud->g1i[lc]  = 0.0;
                      save_cloud->g2i[lc]  = 0.0;
                      save_cloud->fi[lc]   = 0.0;
                      save_cloud->ssai[lc] = 0.0;
                    }
                  }
                }

                /* calculate effective cloud optical properties for fractional cloud cover */
                status = generate_effective_cloud (input, output, save_cloud, iv, iq, input.verbose);
                CHKERR (status);
              }

              /* 3DAbs include caoth3d for 3D molecular atmosphere, right place here ??? */
              if (input.rte.solver == SOLVER_MONTECARLO && output->molecular3d)
                optical_properties_molecular3d (input, output, &(output->caoth3d[CAOTH_FIR]), iv, iq);

              /* setup optical properties and call the RTE solver */
              status = setup_and_call_solver (input,
                                              output,
                                              rte_in,
                                              rte_out,
                                              raman_qsrc_components,
                                              iv,
                                              iq,
                                              ir,
                                              output->atm.threed,
                                              mc_loaddata);

              CHKERR (status);

              /* verbose output */
              if (input.verbose) {
                fprintf (stderr,
                         "  iv = %d, %f nm, iq = %d, flux_dir[lu=0] = %13.7e, flux_dn[lu=0] = %13.7e, flux_up[lu=0] = %13.7e, "
                         "weight_r = %13.7e \n",
                         iv,
                         output->wl.lambda_r[iv],
                         iq,
                         rte_out->rfldir[0],
                         rte_out->rfldn[0],
                         rte_out->flup[0],
                         output->atm.wght_r[iv][iq]);
              }

              /* data need to be loaded only once per iv */
              mc_loaddata = 0;

              if (input.heating != HEAT_NONE)
                calc_spectral_heating (input,
                                       output,
                                       dz,
                                       rho_mass_zout,
                                       k_abs,
                                       k_abs_layer,
                                       rte_out,
                                       rte_out->heat,
                                       rte_out->emis,
                                       rte_out->w_zout,
                                       iv);

              /**********************************************************************/
              /* Store intensities for later use in second round of raman iteration */
              /**********************************************************************/
              if (input.raman) {
                if (ir == 0) {
                  if (input.raman_fast) {
                    for (int lu = 0; lu < output->atm.nlyr; lu++)
                      raman_qsrc_components->dtauc[lu][iq] = (double) output->dtauc[lu];
                    raman_qsrc_components->fbeam[iq] = rte_in->fbeam;
                    for (int lu = 0; lu < output->atm.nzout; lu++) {
                      raman_qsrc_components->uavgso[lu][iq] = (double) rte_out->uavgso[lu];
                      raman_qsrc_components->uavgdn[lu][iq] = (double) rte_out->uavgdn[lu];
                      raman_qsrc_components->uavgup[lu][iq] = (double) rte_out->uavgup[lu];
                      raman_qsrc_components->rfldir[lu][iq] = (double) rte_out->rfldir[lu];
                      raman_qsrc_components->rfldn[lu][iq]  = (double) rte_out->rfldn[lu];
                      raman_qsrc_components->flup[lu][iq]   = (double) rte_out->flup[lu];
                      for (int iu = 0; iu < input.rte.nstr; iu++) {
                        raman_qsrc_components->u0u[lu][input.rte.cmuind[iu]][iq] = (double) rte_out->u0u[lu][input.rte.cmuind[iu]];
                        for (int j = 0; j < input.rte.nphi; j++) {
                          raman_qsrc_components->uu[lu][j][input.rte.cmuind[iu]][iq] =
                            (double) rte_out->uu[j][lu][input.rte.cmuind[iu]];
                        }
                      }
                      for (int iu = 0; iu < input.rte.numu - input.rte.nstr; iu++) {
                        raman_qsrc_components->u0u[lu][input.rte.umuind[iu]][iq] = (double) rte_out->u0u[lu][input.rte.umuind[iu]];
                        for (int j = 0; j < input.rte.nphi; j++) {
                          raman_qsrc_components->uu[lu][j][input.rte.umuind[iu]][iq] =
                            (double) rte_out->uu[j][lu][input.rte.umuind[iu]];
                        }
                      }
                    }

                    for (int lu = 0; lu < output->atm.nzout; lu++)
                      for (int j = 0; j < input.rte.nstr; j++)
                        for (int iu = 0; iu < input.rte.numu; iu++)
                          raman_qsrc_components->uum[lu][j][iu][iq] = (double) rte_out->uum[j][lu][iu];
                  } else {
                    raman_qsrc_components->fbeam[iq] = (double) rte_in->fbeam;
                    if (input.verbose)
                      fprintf (stderr,
                               "Storing Raman quantities for iq: %3d out of %3d.\n",
                               iq,
                               output->crs.number_of_ramanwavelengths - 1);
                    if (iq == output->crs.number_of_ramanwavelengths - 1) {
                      /* Only store radiation for the wanted wavelength which should be at the last index */
                      for (int lu = 0; lu < output->atm.nlyr; lu++)
                        raman_qsrc_components->dtauc[lu][iq] = (double) output->dtauc[lu];
                      for (int lu = 0; lu < output->atm.nzout; lu++) {
                        uavgso_raman[lu]    = rte_out->uavgso[lu];
                        uavgdn_raman[lu]    = rte_out->uavgdn[lu];
                        uavgup_raman[lu]    = rte_out->uavgup[lu];
                        rfldir_raman[lu]    = rte_out->rfldir[lu];
                        rfldn_raman[lu]     = rte_out->rfldn[lu];
                        flup_raman[lu]      = rte_out->flup[lu];
                        rte_out->rfldir[lu] = 0;
                        rte_out->rfldn[lu]  = 0;
                        rte_out->flup[lu]   = 0;
                        for (int iu = 0; iu < input.rte.nstr; iu++) {
                          u0u_raman[lu][input.rte.cmuind[iu]] = rte_out->u0u[lu][input.rte.cmuind[iu]];
                          for (int j = 0; j < input.rte.nphi; j++)
                            uu_raman[lu][j][input.rte.cmuind[iu]] = rte_out->uu[j][lu][input.rte.cmuind[iu]];
                        }
                        for (int iu = 0; iu < input.rte.numu - input.rte.nstr; iu++) {
                          u0u_raman[lu][input.rte.umuind[iu]] = rte_out->u0u[lu][input.rte.umuind[iu]];
                          for (int j = 0; j < input.rte.nphi; j++)
                            uu_raman[lu][j][input.rte.umuind[iu]] = rte_out->uu[j][lu][input.rte.umuind[iu]];
                        }
                      }
                    }

                    /* Store source components at all shifted wavelengths. Store uum for all wavelengths, */
                    /* including the last index which contains the wanted wavelength                      */

                    for (int lu = 0; lu < output->atm.nlyr; lu++)
                      raman_qsrc_components->dtauc[lu][iq] = (double) output->dtauc[lu];
                    for (int lu = 0; lu < output->atm.nzout; lu++)
                      for (int j = 0; j < input.rte.nstr; j++)
                        for (int iu = 0; iu < input.rte.numu; iu++)
                          raman_qsrc_components->uum[lu][j][iu][iq] = (double) rte_out->uum[j][lu][iu];
                  }
                } else if (ir == 1) {

                  add = 1;
                  if (input.raman_fast) {
                    ivr = iv + output->wl.iv_lower;
                    for (int lu = 0; lu < output->atm.nzout; lu++) {
                      rte_out->uavgso[lu] += (double) raman_qsrc_components->uavgso[lu][ivr];
                      rte_out->uavgup[lu] += (double) raman_qsrc_components->uavgup[lu][ivr];
                      rte_out->uavgdn[lu] += (double) raman_qsrc_components->uavgdn[lu][ivr];
                      rte_out->rfldir[lu] += (double) raman_qsrc_components->rfldir[lu][ivr];
                      rte_out->rfldn[lu] += (double) raman_qsrc_components->rfldn[lu][ivr];
                      rte_out->flup[lu] += (double) raman_qsrc_components->flup[lu][ivr];
                      for (int iu = 0; iu < input.rte.nstr; iu++) {
                        rte_out->u0u[lu][input.rte.cmuind[iu]] +=
                          (double) raman_qsrc_components->u0u[lu][input.rte.cmuind[iu]][ivr];
                        for (int j = 0; j < input.rte.nphi; j++) {
                          rte_out->uu[j][lu][input.rte.cmuind[iu]] +=
                            (double) raman_qsrc_components->uu[lu][j][input.rte.cmuind[iu]][ivr];
                        }
                      }
                      for (int iu = 0; iu < input.rte.numu - input.rte.nstr; iu++) {
                        rte_out->u0u[lu][input.rte.umuind[iu]] +=
                          (double) raman_qsrc_components->u0u[lu][input.rte.umuind[iu]][ivr];
                        for (int j = 0; j < input.rte.nphi; j++) {
                          rte_out->uu[j][lu][input.rte.umuind[iu]] +=
                            (double) raman_qsrc_components->uu[lu][j][input.rte.umuind[iu]][ivr];
                        }
                      }
                    }
                  } else {
                    for (int lu = 0; lu < output->atm.nzout; lu++) {
                      if (add) {
                        rte_out->uavgso[lu] += uavgso_raman[lu];
                        rte_out->uavgdn[lu] += uavgdn_raman[lu];
                        rte_out->uavgup[lu] += uavgup_raman[lu];
                        rte_out->rfldir[lu] += rfldir_raman[lu];
                        rte_out->rfldn[lu] += rfldn_raman[lu];
                        rte_out->flup[lu] += flup_raman[lu];
                        for (int iu = 0; iu < input.rte.nstr; iu++) {
                          rte_out->u0u[lu][input.rte.cmuind[iu]] += u0u_raman[lu][input.rte.cmuind[iu]];
                          for (int j = 0; j < input.rte.nphi; j++) {
                            rte_out->uu[j][lu][input.rte.cmuind[iu]] += uu_raman[lu][j][input.rte.cmuind[iu]];
                          }
                        }
                        for (int iu = 0; iu < input.rte.numu - input.rte.nstr; iu++) {
                          rte_out->u0u[lu][input.rte.umuind[iu]] += u0u_raman[lu][input.rte.umuind[iu]];
                          for (int j = 0; j < input.rte.nphi; j++) {
                            rte_out->uu[j][lu][input.rte.umuind[iu]] += uu_raman[lu][j][input.rte.umuind[iu]];
                          }
                        }
                      }
                    }
                  }
                }
              }

              /* add result for the current quadrature point considering quadrature weight */
              if (!input.raman || (input.raman && ir == 0 && iq == output->crs.number_of_ramanwavelengths - 1) ||
                  (input.raman && ir == 1)) {

                if (input.raman) {
                  if (ir == 0)
                    weight = 0;
                  else if (ir == 1)
                    weight = 1;
                } else
                  weight = output->atm.wght_r[iv][iq];

                if (input.rte.mc.spectral_is)
                  for (int iv_alis = 0; iv_alis < output->mc.alis.nlambda_abs; iv_alis++)
                    weight_spectral[iv_alis] = output->atm.wght_r[iv_alis][iq];

                status = add_rte_output (rte_outband,
                                         rte_out,
                                         weight,
                                         weight_spectral,
                                         input,
                                         output->atm.nzout,
                                         output->atm.Nxcld,
                                         output->atm.Nycld,
                                         output->atm.Nzcld,
                                         output->mc.alis.Nc,
                                         output->atm.nlev - 1,
                                         output->mc.alis.nlambda_abs,
                                         output->atm.threed,
                                         output->islower,
                                         output->isupper,
                                         output->jslower,
                                         output->jsupper,
                                         output->isstep,
                                         output->jsstep);
                CHKERR (status);
              }

            } /* for (iq=0; iq<output->atm.nq_r[iv]; iq++) { ...   == 'loop over quadrature points' */

            /************************************************/
            /* add result for the current independent pixel */
            /************************************************/

            if (rte_in->ibcnd) {
              for (int iu = 0; iu < input.rte.numu; iu++) {
                output->albmed_r[iu][iv] += rte_outband->albmed[iu];
                output->trnmed_r[iu][iv] += rte_outband->trnmed[iu];
              }
            } else {
              for (int lu = 0; lu < output->atm.nzout; lu++) {
                output->rfldir_r[lu][iv] += output->ipaweight[ipa] * rte_outband->rfldir[lu];
                output->rfldn_r[lu][iv] += output->ipaweight[ipa] * rte_outband->rfldn[lu];
                output->flup_r[lu][iv] += output->ipaweight[ipa] * rte_outband->flup[lu];
                output->uavg_r[lu][iv] += output->ipaweight[ipa] * rte_outband->uavg[lu];
                output->uavgdn_r[lu][iv] += output->ipaweight[ipa] * rte_outband->uavgdn[lu];
                output->uavgso_r[lu][iv] += output->ipaweight[ipa] * rte_outband->uavgso[lu];
                output->uavgup_r[lu][iv] += output->ipaweight[ipa] * rte_outband->uavgup[lu];
                output->heat_r[lu][iv] += output->ipaweight[ipa] * rte_outband->heat[lu];
                output->emis_r[lu][iv] += output->ipaweight[ipa] * rte_outband->emis[lu];
                output->w_zout_r[lu][iv] += output->ipaweight[ipa] * rte_outband->w_zout[lu];
                output->sslidar_nphot_r[lu][iv] += output->ipaweight[ipa] * rte_outband->sslidar_nphot[lu];
                output->sslidar_nphot_q_r[lu][iv] += output->ipaweight[ipa] * rte_outband->sslidar_nphot_q[lu];
                output->sslidar_ratio_r[lu][iv] += output->ipaweight[ipa] * rte_outband->sslidar_ratio[lu];

                /* intensities */
                for (int iu = 0; iu < input.rte.numu; iu++) {
                  output->u0u_r[lu][iu][iv] += output->ipaweight[ipa] * rte_outband->u0u[lu][iu];

                  for (int j = 0; j < input.rte.nphi; j++)
                    output->uu_r[lu][j][iu][iv] += output->ipaweight[ipa] * rte_outband->uu[j][lu][iu];
                }

                /* 3D fields */
                if (input.output3D && (input.rte.solver == SOLVER_MONTECARLO || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM)) {

                  for (int is = output->islower; is <= output->isupper; is += output->isstep) {
                    for (int js = output->jslower; js <= output->jsupper; js += output->jsstep) {

                      output->rfldir3d_r[lu][is][js][iv] += output->ipaweight[ipa] * rte_outband->rfldir3d[lu][is][js];
                      output->rfldn3d_r[lu][is][js][iv] += output->ipaweight[ipa] * rte_outband->rfldn3d[lu][is][js];
                      output->flup3d_r[lu][is][js][iv] += output->ipaweight[ipa] * rte_outband->flup3d[lu][is][js];
                      output->uavgso3d_r[lu][is][js][iv] += output->ipaweight[ipa] * rte_outband->uavgso3d[lu][is][js];
                      output->uavgdn3d_r[lu][is][js][iv] += output->ipaweight[ipa] * rte_outband->uavgdn3d[lu][is][js];
                      output->uavgup3d_r[lu][is][js][iv] += output->ipaweight[ipa] * rte_outband->uavgup3d[lu][is][js];

                      if (output->mc.sample.spectral_is)
                        for (int ivs = 0; ivs < output->mc.alis.nlambda_abs; ivs++)
                          output->fl3d_is_r[lu][is][js][0][ivs] +=
                            output->ipaweight[ipa] * rte_outband->fl3d_is[lu][0][is][js][ivs];

                      if (output->mc.sample.concentration_is)
                        for (int ic = 0; ic < output->mc.alis.Nc; ic++)
                          output->fl3d_is_r[lu][is][js][ic][iv] += output->ipaweight[ipa] * rte_outband->fl3d_is[lu][ic][is][js][0];

                      if (input.output3Dradiances) {
                        if (output->mc.sample.spectral_is || output->mc.sample.concentration_is) {
                          for (int ip = 0; ip < input.rte.mc.nstokes; ip++)
                            for (int ivs = 0; ivs < output->mc.alis.nlambda_abs; ivs++)
                              for (int ic = 0; ic < output->mc.alis.Nc; ic++)
                                output->radiance3d_r[lu][is][js][ip][ic][ivs] +=
                                  output->ipaweight[ipa] * rte_outband->radiance3d_is[lu][ic][is][js][ip][ivs];
                        } else {
                          for (int ip = 0; ip < input.rte.mc.nstokes; ip++)
                            output->radiance3d_r[lu][is][js][ip][0][iv] +=
                              output->ipaweight[ipa] * rte_outband->radiance3d[lu][is][js][ip];
                        }
                      }

                      if (input.rte.mc.jacobian[DIM_1D]) {
                        for (int isp = 0; isp < input.n_caoth + 2; isp++)
                          for (int ijac = 0; ijac < 2; ijac++)
                            for (int lc = 0; lc < output->atm.nlyr; lc++) /* scattering and absorption */
                              output->jacobian_r[lu][is][js][isp][ijac][lc][iv] +=
                                output->ipaweight[ipa] * rte_outband->jacobian[lu][is][js][isp][ijac][lc];
                      }

                      if (input.rte.mc.backward.absorption)
                        output->absback3d_r[lu][is][js][iv] += output->ipaweight[ipa] * rte_outband->absback3d[lu][is][js];

                      /* variances */
                      if (input.rte.mc.std) {

                        /* variance is weighted with square of weight */
                        weight2 = output->ipaweight[ipa] * output->ipaweight[ipa];

                        output->rfldir3d_var_r[lu][is][js][iv] += weight2 * rte_outband->rfldir3d_var[lu][is][js];
                        output->rfldn3d_var_r[lu][is][js][iv] += weight2 * rte_outband->rfldn3d_var[lu][is][js];
                        output->flup3d_var_r[lu][is][js][iv] += weight2 * rte_outband->flup3d_var[lu][is][js];
                        output->uavgso3d_var_r[lu][is][js][iv] += weight2 * rte_outband->uavgso3d_var[lu][is][js];
                        output->uavgdn3d_var_r[lu][is][js][iv] += weight2 * rte_outband->uavgdn3d_var[lu][is][js];
                        output->uavgup3d_var_r[lu][is][js][iv] += weight2 * rte_outband->uavgup3d_var[lu][is][js];

                        if (input.output3Dradiances) {
                          for (int ip = 0; ip < input.rte.mc.nstokes; ip++)
                            output->radiance3d_var_r[lu][is][js][ip][iv] += weight2 * rte_outband->radiance3d_var[lu][is][js][ip];
                        }

                        if (input.rte.mc.backward.absorption)
                          output->absback3d_var_r[lu][is][js][iv] += weight2 * rte_outband->absback3d_var[lu][is][js];
                      }
                    }
                  }
                } else if (input.ipa3d) {
                  /*ulrike: 3d-fields (without _r) are not needed*/

                  /* BM 3.7.2020: removed factor ipaweight[]; see comment above */
                  output->rfldir3d_r[lu][iipa][jipa][iv] += rte_outband->rfldir[lu];
                  output->rfldn3d_r[lu][iipa][jipa][iv] += rte_outband->rfldn[lu];
                  output->flup3d_r[lu][iipa][jipa][iv] += rte_outband->flup[lu];
                  output->uavgso3d_r[lu][iipa][jipa][iv] += rte_outband->uavgso[lu];
                  output->uavgdn3d_r[lu][iipa][jipa][iv] += rte_outband->uavgdn[lu];
                  output->uavgup3d_r[lu][iipa][jipa][iv] += rte_outband->uavgup[lu];

                } /*ulrike: end of: else if (input.ipa3d)*/
              }   /*ulrike: end for-loop over lev lu*/

              /* 3D absorption */
              if (input.output3D && ((input.rte.mc.absorption != MCFORWARD_ABS_NONE && input.rte.solver == SOLVER_MONTECARLO) ||
                                     input.rte.solver == SOLVER_DYNAMIC_TENSTREAM)) {
                size_t k3d = 0;  // count 3d layers, because output->abs3d_r has only entries for 3d layers
                for (int ks = 0; ks < output->atm.Nzcld; ks++) {
                  if (output->atm.threed[ks]) { /* only for 3D layers */
                    for (int is = 0; is < output->atm.Nxcld; is++) {
                      for (int js = 0; js < output->atm.Nycld; js++) {
                        output->abs3d_r[k3d][is][js][iv] += output->ipaweight[ipa] * rte_outband->abs3d[ks][is][js];
                        if (input.rte.mc.std)
                          output->abs3d_var_r[k3d][is][js][iv] += output->ipaweight[ipa] * rte_outband->abs3d_var[ks][is][js];
                      }
                    }
                    k3d++;
                  }
                }
              } else if (input.ipa3d && input.heating != HEAT_NONE) {

                /* calculation of net flux differences is a bit tricky since we need   */
                /* to figure out the correct levels above and below each 3D layer;     */
                /* for OUTPUT_3D_LEVELS we selected exactly all levels above and below */
                /* all 3D layers.                                                      */

                size_t k3d  = 0;  // count 3d layers, because output->abs3d_r has only entries for 3d layers
                int    klev = 0;  /* level counter */

                int notfirst = 0;
                for (int ks = 0; ks < output->atm.Nzcld; ks++) {
                  if (output->atm.threed[ks]) {

                    /* for all 3D layers except the first one we need to add 1 to */
                    /* the level index if there is a non-3D layer below           */
                    if (notfirst && output->atm.threed[ks - 1] == 0)
                      klev++;

                    notfirst = 1;

                    /* heating rate = Enet[k3d+1] - Enet[k3d] */
                    double dEnet = (rte_outband->rfldir[klev + 1] + rte_outband->rfldn[klev + 1] - rte_outband->flup[klev + 1]) -
                                   (rte_outband->rfldir[klev] + rte_outband->rfldn[klev] - rte_outband->flup[klev]);

                    output->abs3d_r[k3d][iipa][jipa][iv] += dEnet;

                    k3d++;
                    klev++;
                  }
                }
              } /* else if (input.ipa3d) */
            }   /* endof of if {rte_in->ibcnd} else {} */

            if (rte_outband->triangle_results) {  // Add triangle result -> output->result
              const double factor = output->ipaweight[ipa];
              status = add_triangular_surface_result (factor, rte_outband->triangle_results, output->triangle_results_r[iv]);
              CHKERR (status);
            }

            /* verbose output */
            if (input.verbose) {
              fprintf (stderr,
                       "  iv = %d, %f nm, sum iq, flux_dir[lu=0] = %13.7e, flux_dn[lu=0] = %13.7e, flux_up[lu=0] = %13.7e \n",
                       iv,
                       output->wl.lambda_r[iv],
                       output->rfldir_r[0][iv],
                       output->rfldn_r[0][iv],
                       output->flup_r[0][iv]);
            }
          }  // for (jipa=0; ipa<output->njipa; jipa++) independent pixel
        }    // for (iipa=0; ipa<output->niipa; iipa++) independent pixel
      }      // for (ipa=0; ipa<output->nipa; ipa++) independent pixel

      /* change unit of the solar spectrum [e.g. W/(m2 nm)] or terrestral spectrum [e.g. W/(m2 cm-1)]  */
      /* to output units wanted by the user: 'output per_nm', 'output per_cm', or 'output per_band' */
      /* but only, when dealing with unit (not transmission or reflectivity).                          */
      /* Unit conversion must happen before interpolate transmittance, as some unit conversions        */
      /* use the internal thermal bandwidths or correlated-k bandwidth.                                */
      /* UH 2006-03                                                                                    */

      if (output->wl.use_reptran)
        unit_factor = 1; /* conversion is done in internal_to_transmittance_grid() */
      else
        unit_factor = get_unit_factor (input, output, iv);

      if (unit_factor <= 0) {
        fprintf (stderr, "Error, calculating unit_factor = %f\n", unit_factor);
        CHKERR (-1);
      }

      switch (input.source) {
      case SRC_THERMAL:
        ffactor = unit_factor;
        rfactor = unit_factor;
        break;

      case SRC_SOLAR:
      case SRC_BLITZ: /* BCA */
      case SRC_LIDAR: /* BCA */
        switch (input.processing) {
        case PROCESS_INT:
        case PROCESS_SUM:
        case PROCESS_RGB:
        case PROCESS_RGBRAW:
          ffactor = unit_factor;
          rfactor = unit_factor;
          break;

        case PROCESS_NONE:
        case PROCESS_RAMAN:
          switch (input.calibration) {
          case OUTCAL_ABSOLUTE:
            ffactor = unit_factor;
            rfactor = unit_factor;
            break;

          case OUTCAL_TRANSMITTANCE:
            ffactor = 1.0;
            rfactor = 1.0;
            break;

          case OUTCAL_REFLECTIVITY:
            ffactor = 1.0;
            rfactor = 1.0;
            break;

          default:
            fprintf (stderr, "Error, unknown output calibration %d\n", input.calibration);
            return -1;
          }

          break;

        default:
          fprintf (stderr, "Error, unknown output processing %d\n", input.processing);
          return -1;
        }

        break;

      default:
        fprintf (stderr, "Error, unknown source %d\n", input.source);
        return -1;
      }

      hfactor = unit_factor;

      /*****************************************************************************************/
      /* now scale irradiances with ffactor, radiances with rfactor, heating rate with hfactor */
      /*****************************************************************************************/
      status = scale_output (input,
                             &(output->rfldir_r),
                             &(output->rfldn_r),
                             &(output->flup_r),
                             &(output->albmed_r),
                             &(output->trnmed_r),
                             &(output->uavgso_r),
                             &(output->uavgdn_r),
                             &(output->uavgup_r),
                             &(output->uavg_r),
                             &(output->u0u_r),
                             &(output->uu_r),
                             &(output->heat_r),
                             &(output->emis_r),
                             &(output->w_zout_r),
                             &(output->rfldir3d_r),
                             &(output->rfldn3d_r),
                             &(output->flup3d_r),
                             &(output->fl3d_is_r),
                             &(output->uavgso3d_r),
                             &(output->uavgdn3d_r),
                             &(output->uavgup3d_r),
                             &(output->radiance3d_r),
                             &(output->jacobian_r),
                             &(output->absback3d_r),
                             &(output->rfldir3d_var_r),
                             &(output->rfldn3d_var_r),
                             &(output->flup3d_var_r),
                             &(output->uavgso3d_var_r),
                             &(output->uavgdn3d_var_r),
                             &(output->uavgup3d_var_r),
                             &(output->radiance3d_var_r),
                             &(output->abs3d_var_r),
                             &(output->absback3d_var_r),
                             output->atm.nzout,
                             output->atm.Nxcld,
                             output->atm.Nycld,
                             output->atm.Nzcld,
                             output->mc.alis.Nc,
                             output->atm.nlev - 1,
                             output->atm.threed,
                             output->islower,
                             output->isupper,
                             output->jslower,
                             output->jsupper,
                             output->isstep,
                             output->jsstep,
                             &(output->abs3d_r),
                             output->triangle_results_r,
                             ffactor,
                             rfactor,
                             hfactor,
                             iv);
      /* in ancillary.c */
      CHKERR (status);

      /* free 3D cloud properties */
      if (input.rte.solver == SOLVER_MONTECARLO)
        for (int isp = 0; isp < input.n_caoth; isp++) {
          status = free_caoth_mystic (input.caoth[isp].properties, &(output->caoth3d[isp]));
          CHKERR (status);
        }
    } /*   for (ir=0; ir<nr; ir++) */
    iv_count++;

  } /* for (iv=output->wl.iv_lower; iv<=output->wl.iv_upper; iv++) */

  /* ulrike: free msorted and zsorted (for tipa dir!!! for tipa
             dirdiff msorted is freed already),
             output->(w/i)c.tipa.taudircld, ... */
  if (input.tipa == TIPA_DIR || input.rte.mc.tipa == TIPA_DIR) {
    for (int isp = 0; isp < input.n_caoth; isp++)
      if (input.caoth[isp].source == CAOTH_FROM_3D) {
        /* free m-and z-sorted for caoth */
        for (int iz = 0; iz < (output->caoth[isp].tipa.nztilt); iz++) {
          for (int ks = 0; ks < (output->caoth[isp].tipa.totlev[iz]); ks++)
            free ((output->caoth[isp].tipa.msorted)[iz][ks]);
          free ((output->caoth[isp].tipa.msorted)[iz]);
          free ((output->caoth[isp].tipa.zsorted)[iz]);
        }
        free (output->caoth[isp].tipa.msorted);
        free (output->caoth[isp].tipa.zsorted);

        for (int js = 0; js < (upper_wl_id - lower_wl_id + 1); js++) /* free taudircld for wc */
          free ((output->caoth[isp].tipa.taudircld)[js]);
        free (output->caoth[isp].tipa.taudircld);
      }
  }

  /* free temporary memory */

  status = free_rte_output (rte_out,
                            input,
                            output->atm.nzout,
                            output->atm.Nxcld,
                            output->atm.Nycld,
                            output->atm.Nzcld,
                            output->mc.sample.Nx,
                            output->mc.sample.Ny,
                            output->mc.alis.Nc,
                            output->atm.nlyr - 1,
                            output->mc.alis.nlambda_abs,
                            output->atm.threed);
  CHKERR (status);

  status = free_rte_output (rte_outband,
                            input,
                            output->atm.nzout,
                            output->atm.Nxcld,
                            output->atm.Nycld,
                            output->atm.Nzcld,
                            output->mc.sample.Nx,
                            output->mc.sample.Ny,
                            output->mc.alis.Nc,
                            output->atm.nlyr - 1,
                            output->mc.alis.nlambda_abs,
                            output->atm.threed);
  CHKERR (status);

  free (rte_in->hl);
  free (rte_in->utau);
  free (rte_in);
  if (input.rte.solver == SOLVER_DISORT && input.raman) {
    if (uu_raman != NULL)
      ASCII_free_double_3D (uu_raman, output->atm.nzout, input.rte.nphi);
    ASCII_free_double (u0u_raman, output->atm.nzout);
    free (uavgso_raman);
    free (uavgdn_raman);
    free (uavgup_raman);
    free (rfldir_raman);
    free (rfldn_raman);
    free (flup_raman);
  }

  if (input.raman) {
    free_raman_qsrc_components (raman_qsrc_components,
                                input.raman_fast,
                                output->atm.nzout,
                                input.rte.maxumu,
                                input.rte.nphi,
                                input.rte.nstr);
  }

  if (input.heating != HEAT_NONE) {
    free (dz);
    free (k_abs_layer);
    free (k_abs);
  }

  if (save_cloud != NULL) {
    free (save_cloud->tauw);
    free (save_cloud->taui);
    free (save_cloud->g1d);
    free (save_cloud->g2d);
    free (save_cloud->fd);
    free (save_cloud->g1i);
    free (save_cloud->g2i);
    free (save_cloud->fi);
    free (save_cloud->ssaw);
    free (save_cloud->ssai);
    free (save_cloud);
  }

#if HAVE_LIBGSL
#ifdef WRITERANDOMSTATUS
  if (remove (input.filename[FN_RANDOMSTATUS]) != 0)
    fprintf (stderr, "Error deleting randomstatusfile");
#endif
#endif

  return 0;
}  // end of solve_rte()

/* get factor for unit conversion */
static double get_unit_factor (input_struct input, output_struct* output, int iv) {
  double unit_factor = 0.0;

  switch (output->spectrum_unit) {
  case UNIT_PER_NM:
    switch (input.output_unit) {
    case UNIT_PER_NM:
      unit_factor = 1.0;
      break;
    case UNIT_PER_CM_1:
      /* unit_factor = (lambda/k) */ /* (lambda/k) = (lambda**2) / 1.0e+7 */ /* 1.0e+7 == cm -> nm; */
      unit_factor = (output->wl.lambda_r[iv] * output->wl.lambda_r[iv]) / 1.0e+7;
      break;
    case UNIT_PER_BAND:
      /* unit_factor = delta_lambda */ /* lambda_max = 1.0e+7 / k_lower;  lambda_min = 1.0e+7 / k_upper */
      unit_factor = 1.0e+7 / output->wl.wvnmlo_r[iv] - 1.0e+7 / output->wl.wvnmhi_r[iv];
      break;
    case UNIT_NOT_DEFINED:
      unit_factor = 1.0;
      break;
    default:
      fprintf (stderr, "Error, unknown output unit %d\n", input.output_unit);
      CHKERR (-1);
    }
    break;
  case UNIT_PER_CM_1:
    switch (input.output_unit) {
    case UNIT_PER_NM:
      /* unit_factor = (k/lambda) */ /* (k/lambda)= 1.0e+7 / (lambda**2) */ /* 1.0e+7 == cm -> nm; */
      /* k wavenumber in 1/cm**-1, lambda in nm */
      unit_factor = 1.0e+7 / (output->wl.lambda_r[iv] * output->wl.lambda_r[iv]);
      break;
    case UNIT_PER_CM_1:
      unit_factor = 1.0;
      break;
    case UNIT_PER_BAND:
      /* unit_factor = delta_k */
      unit_factor = output->wl.wvnmhi_r[iv] - output->wl.wvnmlo_r[iv];
      break;
    case UNIT_NOT_DEFINED:
      unit_factor = 1.0;
      break;
    default:
      fprintf (stderr, "Error: unknown output unit %d\n", input.output_unit);
      CHKERR (-1);
    }
    break;
  case UNIT_PER_BAND:
    switch (input.output_unit) {
    case UNIT_PER_NM:
      unit_factor = 1.0 / (1.0e7 / output->wl.wvnmlo_r[iv] - 1.0e7 / output->wl.wvnmhi_r[iv]);
      break;
    case UNIT_PER_CM_1:
      unit_factor = 1.0 / (output->wl.wvnmhi_r[iv] - output->wl.wvnmlo_r[iv]);
      break;
    case UNIT_PER_BAND:
      unit_factor = 1.0;
      break;
    case UNIT_NOT_DEFINED: /* not defined */
      unit_factor = 1.0;
      break;
    default:
      fprintf (stderr, "Error: unknown output unit %d\n", input.output_unit);
      CHKERR (-1);
    }
    break;
  case UNIT_NOT_DEFINED:
    switch (input.output_unit) {
    case UNIT_PER_NM:
    case UNIT_PER_CM_1:
    case UNIT_PER_BAND:
      CHKERROUT (-1,
                 "Error, cannot convert undefined solar spectrum to output with units. "
                 "please use 'solar_file filename unit' in order to "
                 "specify the unit of the spectrum");
      break;
    case UNIT_NOT_DEFINED: /* not defined */
      unit_factor = 1.0;
      break;
    default:
      fprintf (stderr, "Error, unknown output unit %d\n", input.output_unit);
      CHKERR (-1);
    }
    break;
  default:
    fprintf (stderr, "Error, unknown output unit %d\n", output->spectrum_unit);
    CHKERR (-1);
  }

  return unit_factor;
}  // end of get_unit_factor()

int setup_result (input_struct   input,
                  output_struct* output,
                  float**        p_dz,
                  double**       p_rho_mass_zout,
                  float**        p_k_abs_layer,
                  float**        p_k_abs,
                  float**        p_k_abs_outband) {
  int status;

  /* FIX 3DAbs need to be initialized when aerosol is not set up, now redundant ??? */
  /* output->mc.alis.Nc=1; */

  status = ASCII_calloc_float (&output->flup_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->rfldir_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->rfldn_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->uavg_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->uavgdn_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->uavgso_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->uavgup_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->heat_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->emis_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->w_zout_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->sslidar_nphot_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->sslidar_nphot_q_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->sslidar_ratio_r, output->atm.nzout, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->albmed_r, input.rte.numu, output->wl.nlambda_r);
  CHKERR (status);
  status = ASCII_calloc_float (&output->trnmed_r, input.rte.numu, output->wl.nlambda_r);
  CHKERR (status);

  /* variables in order to calculate heating rates (by actinic flux) */
  if (input.heating != HEAT_NONE) {

    *p_dz = calloc (output->atm.nlyr, sizeof (float)); /* dz in m for all (nlyr) layers */
    CHKPOINTER (*p_dz);

    /* Initialisation */
    for (int lc = 0; lc < output->atm.nlyr; lc++)
      (*p_dz)[lc] = (output->atm.zd[lc] - output->atm.zd[lc + 1]) * 1000.0; /* km -> m */

    *p_rho_mass_zout = calloc (output->atm.nzout, sizeof (double));
    CHKPOINTER (*p_rho_mass_zout);

    *p_k_abs_layer = calloc (output->atm.nlyr, sizeof (float));
    CHKPOINTER (*p_k_abs_layer);

    *p_k_abs = calloc (output->atm.nzout, sizeof (float));
    CHKPOINTER (*p_k_abs);

    *p_k_abs_outband = calloc (output->atm.nzout, sizeof (float));
    CHKPOINTER (*p_k_abs_outband);
  }

  if (input.rte.numu > 0) {
    status = ASCII_calloc_float_3D (&output->u0u_r, output->atm.nzout, input.rte.numu, output->wl.nlambda_r);
    CHKERR (status);
  }

  if (input.rte.numu > 0 && input.rte.nphi > 0) {
    status = ASCII_calloc_float_4D (&output->uu_r, output->atm.nzout, input.rte.nphi, input.rte.numu, output->wl.nlambda_r);
    CHKERR (status);
  }

  if (input.output3D) {
    const size_t Nz        = output->atm.nzout;
    const size_t Nx        = output->isupper - output->islower + 1;
    const size_t Ny        = output->jsupper - output->jslower + 1;
    const size_t Nlambda_r = output->wl.iv_upper - output->wl.iv_lower + 1;
    const size_t Nalis     = output->mc.alis.Nc;
    const size_t Nstokes   = input.rte.mc.nstokes;

    { /** Allocate memory for 3D results output->*_r */

      {  // alloc 4d vars like flxs and absorption
        const size_t dims[4]    = {Nz, Nx, Ny, Nlambda_r};
        const int    offsets[4] = {0, -output->islower, -output->jslower, -output->wl.iv_lower};

        if (!input.quiet)
          fprintf (stderr,
                   " ... allocating Nz %lu x Nx %lu x Ny %lu x Nwvl %lu = %lu entries (%.3f Gb) foreach 3D flux var\n",
                   dims[0],
                   dims[1],
                   dims[2],
                   dims[3],
                   dims[0] * dims[1] * dims[2] * dims[3],
                   dims[0] * dims[1] * dims[2] * dims[3] * sizeof (float) / 1024. / 1024. / 1024.);

        status = alloc_array_with_offsets_4d (dims, offsets, "rfldir3d_r", &output->rfldir3d_r);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "rfldn3d_r", &output->rfldn3d_r);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "flup3d_r", &output->flup3d_r);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "uavgso3d_r", &output->uavgso3d_r);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "uavgdn3d_r", &output->uavgdn3d_r);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "uavgup3d_r", &output->uavgup3d_r);
        CHKERR (status);

        if (input.rte.mc.std) { /* variances */
          status = alloc_array_with_offsets_4d (dims, offsets, "rfldir3d_var_r", &output->rfldir3d_var_r);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "rfldn3d_var_r", &output->rfldn3d_var_r);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "flup3d_var_r", &output->flup3d_var_r);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "uavgso3d_var_r", &output->uavgso3d_var_r);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "uavgdn3d_var_r", &output->uavgdn3d_var_r);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "uavgup3d_var_r", &output->uavgup3d_var_r);
          CHKERR (status);
        }

        /* backward absorption */
        if (input.rte.mc.backward.absorption) {
          status = alloc_array_with_offsets_4d (dims, offsets, "absback3d_r", &output->absback3d_r);
          CHKERR (status);

          if (input.rte.mc.std) {
            status = alloc_array_with_offsets_4d (dims, offsets, "absback3d_var_r", &output->absback3d_var_r);
            CHKERR (status);
          }
        }

        /* 3d absorption */
        if (input.rte.mc.absorption != MCFORWARD_ABS_NONE || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d) {
          const size_t N3d_layers = vec_count_true (output->atm.Nzcld, output->atm.threed);
          const size_t dims[4]    = {N3d_layers, output->atm.Nxcld, output->atm.Nycld, Nlambda_r};
          const int    offsets[4] = {0, 0, 0, -output->wl.iv_lower};
          if (!input.quiet)
            fprintf (stderr,
                     " ... allocating Nz3D %lu x Nx %lu x Ny %lu x Nwvl %lu = %lu entries (%.3f Gb) foreach 3D absorption var\n",
                     dims[0],
                     dims[1],
                     dims[2],
                     dims[3],
                     dims[0] * dims[1] * dims[2] * dims[3],
                     dims[0] * dims[1] * dims[2] * dims[3] * sizeof (float) / 1024. / 1024. / 1024.);

          status = alloc_array_with_offsets_4d (dims, offsets, "abs3d_r", &output->abs3d_r);
          CHKERR (status);

          if (input.rte.mc.std) {
            status = alloc_array_with_offsets_4d (dims, offsets, "abs3d_var_r", &output->abs3d_var_r);
            CHKERR (status);
          }
        }
      }

      if (input.output3Dradiances) {
        const size_t dims[6]    = {Nz, Nx, Ny, Nstokes, Nalis, Nlambda_r};
        const int    offsets[6] = {0, -output->islower, -output->jslower, 0, 0, -output->wl.iv_lower};

        if (!input.quiet)
          fprintf (stderr,
                   " ... allocating Nz3D %lu x Nx %lu x Ny %lu x Nstokes %lu x Nalis %lu x Nwvl %lu = %lu entries (%.3f Gb) "
                   "foreach 3D radiance var\n",
                   dims[0],
                   dims[1],
                   dims[2],
                   dims[3],
                   dims[4],
                   dims[5],
                   dims[0] * dims[1] * dims[2] * dims[3] * dims[4] * dims[5],
                   dims[0] * dims[1] * dims[2] * dims[3] * dims[4] * dims[5] * sizeof (float) / 1024. / 1024. / 1024.);

        status = alloc_array_with_offsets_6d (dims, offsets, "radiance3d_r", &output->radiance3d_r);
        CHKERR (status);
      }

      if (input.output3Dradiances) {
        if (input.rte.mc.std) {
          const size_t dims[5]    = {Nz, Nx, Ny, Nstokes, Nlambda_r};
          const int    offsets[5] = {0, -output->islower, -output->jslower, 0, -output->wl.iv_lower};

          status = alloc_array_with_offsets_5d (dims, offsets, "radiance3d_var_r", &output->radiance3d_var_r);
          CHKERR (status);
        }
      }

      if (input.rte.mc.spectral_is || input.rte.mc.concentration_is) {
        const size_t dims[5]    = {Nz, Nx, Ny, Nalis, Nlambda_r};
        const int    offsets[5] = {0, -output->islower, -output->jslower, 0, -output->wl.iv_lower};

        status = alloc_array_with_offsets_5d (dims, offsets, "fl3d_is_r", &output->fl3d_is_r);
        CHKERR (status);
      }

      /* So far we allow only 1D output for postprocessing */
      if (input.rte.mc.jacobian[DIM_1D]) {
        const size_t dims[7]    = {Nz, 1, 1, input.n_caoth + 2, 2, output->atm.nlev - 1, Nlambda_r};
        const int    offsets[7] = {0, 0, 0, 0, -output->wl.iv_lower};

        status = alloc_array_with_offsets_7d (dims, offsets, "jacobian_r", &output->jacobian_r);
        CHKERR (status);
      }
    }

  }  // endif input.output3D

  /* need to allocate enough memory for both output->wl.nlambda_s */
  /* and output->wl.nlambda_h, hence using whichever is larger    */
  const int nlambda = MAX (output->wl.nlambda_h, output->wl.nlambda_s);
  {
    status = ASCII_calloc_float (&output->flup, output->atm.nzout, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->rfldn, output->atm.nzout, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->rfldir, output->atm.nzout, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->uavg, output->atm.nzout, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->uavgdn, output->atm.nzout, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->uavgso, output->atm.nzout, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->uavgup, output->atm.nzout, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->heat, output->atm.nzout, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->emis, output->atm.nzout, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->albmed, input.rte.numu, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->trnmed, input.rte.numu, nlambda);
    CHKERR (status);

    status = ASCII_calloc_float (&output->w_zout, output->atm.nzout, nlambda);
    CHKERR (status);

    if (input.rte.numu > 0) {
      status = ASCII_calloc_float_3D (&output->u0u, output->atm.nzout, input.rte.numu, nlambda);
      CHKERR (status);
    }

    if (input.rte.numu > 0 && input.rte.nphi > 0) {
      status = ASCII_calloc_float_4D (&output->uu, output->atm.nzout, input.rte.nphi, input.rte.numu, nlambda);
      CHKERR (status);
    }

    status = ASCII_calloc_float (&output->sslidar_nphot, output->atm.nzout, nlambda);
    CHKERR (status);
    status = ASCII_calloc_float (&output->sslidar_nphot_q, output->atm.nzout, nlambda);
    CHKERR (status);
    status = ASCII_calloc_float (&output->sslidar_ratio, output->atm.nzout, nlambda);
    CHKERR (status);

    output->sza_h = calloc (output->wl.nlambda_h, sizeof (float));
    CHKPOINTER (output->sza_h);
  }

  if (input.output3D) {
    /* allocate only output pixels which are actually required  */
    /* (defined by mc_backward islower jslower isupper jsupper) */
    {
      const size_t Nz      = output->atm.nzout;
      const size_t Nx      = output->isupper - output->islower + 1;
      const size_t Ny      = output->jsupper - output->jslower + 1;
      const size_t Nstokes = input.rte.mc.nstokes;
      const size_t Nalis   = output->mc.alis.Nc;

      {  // alloc 4d vars like flxs and absorption
        const size_t dims[4]    = {Nz, Nx, Ny, nlambda};
        const int    offsets[4] = {0, -output->islower, -output->jslower, 0};

        status = alloc_array_with_offsets_4d (dims, offsets, "rfldir3d", &output->rfldir3d);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "rfldn3d", &output->rfldn3d);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "flup3d", &output->flup3d);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "uavgso3d", &output->uavgso3d);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "uavgdn3d", &output->uavgdn3d);
        CHKERR (status);
        status = alloc_array_with_offsets_4d (dims, offsets, "uavgup3d", &output->uavgup3d);
        CHKERR (status);

        if (input.rte.mc.std) { /* variances */
          status = alloc_array_with_offsets_4d (dims, offsets, "rfldir3d_var", &output->rfldir3d_var);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "rfldn3d_var", &output->rfldn3d_var);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "flup3d_var", &output->flup3d_var);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "uavgso3d_var", &output->uavgso3d_var);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "uavgdn3d_var", &output->uavgdn3d_var);
          CHKERR (status);
          status = alloc_array_with_offsets_4d (dims, offsets, "uavgup3d_var", &output->uavgup3d_var);
          CHKERR (status);
        }

        /* backward absorption */
        if (input.rte.mc.backward.absorption) {
          status = alloc_array_with_offsets_4d (dims, offsets, "absback3d", &output->absback3d);
          CHKERR (status);

          if (input.rte.mc.std) {
            status = alloc_array_with_offsets_4d (dims, offsets, "absback3d_var", &output->absback3d_var);
            CHKERR (status);
          }
        }
      }

      if (input.output3Dradiances) {
        const size_t dims[6]    = {Nz, Nx, Ny, Nstokes, Nalis, nlambda};
        const int    offsets[6] = {0, -output->islower, -output->jslower, 0, 0, 0};

        status = alloc_array_with_offsets_6d (dims, offsets, "radiance3d", &output->radiance3d);
        CHKERR (status);
      }

      if (input.output3Dradiances) {
        if (input.rte.mc.std) {
          const size_t dims[5]    = {Nz, Nx, Ny, Nstokes, nlambda};
          const int    offsets[5] = {0, -output->islower, -output->jslower, 0, 0};

          status = alloc_array_with_offsets_5d (dims, offsets, "radiance3d_var", &output->radiance3d_var);
          CHKERR (status);
        }
      }

      /* 3D absorption */
      if (input.rte.mc.absorption != MCFORWARD_ABS_NONE || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d) {
        const size_t N3d_layers = vec_count_true (output->atm.Nzcld, output->atm.threed);
        const size_t Nx         = output->atm.Nxcld;
        const size_t Ny         = output->atm.Nycld;
        output->abs3d           = calloc_float_4D (N3d_layers, Nx, Ny, nlambda, "abs3d");
        CHKPOINTER (output->abs3d);

        if (input.rte.mc.std) {
          output->abs3d_var = calloc_float_4D (N3d_layers, Nx, Ny, nlambda, "abs3d_var");
          CHKPOINTER (output->abs3d_var);
        }
      }

      if (input.rte.mc.spectral_is || input.rte.mc.concentration_is) {
        output->fl3d_is = calloc_float_5D (Nz, output->mc.sample.Nx, output->mc.sample.Ny, Nalis, nlambda, "fl3d_is");
        CHKPOINTER (output->fl3d_is);
      }

      if (input.rte.mc.jacobian[DIM_1D]) {
        output->jacobian = calloc_float_7D (Nz, 1, 1, input.n_caoth + 2, 2, output->atm.nlev - 1, nlambda, "jacobian");
        CHKPOINTER (output->jacobian);
      }
    }
  }  // output3D

  {  // allocate triangle output
    output->triangle_results_r = NULL;
    output->triangle_results_t = NULL;
    output->triangle_results_o = NULL;

    const int ierr = init_spectral_triangular_surface_result_struct (output->wl.nlambda_r,
                                                                     output->mc.triangular_surface.N_triangles,
                                                                     &(output->triangle_results_r));
    CHKERR (ierr);
  }

  return 0;
}  // end of setup_result()

static int reverse_profiles (input_struct input, output_struct* output) {
  double tmp = 0;

  for (int lu = 0; lu < output->atm.nlyr / 2; lu++) {

    /* reverse profile of optical depth dtauc */
    tmp                                      = output->dtauc[output->atm.nlyr - lu - 1];
    output->dtauc[output->atm.nlyr - lu - 1] = output->dtauc[lu];
    output->dtauc[lu]                        = tmp;

    /* reverse profile of single scattering albedo ssalb */
    tmp                                      = output->ssalb[output->atm.nlyr - lu - 1];
    output->ssalb[output->atm.nlyr - lu - 1] = output->ssalb[lu];
    output->ssalb[lu]                        = tmp;

    /* reverse profile of phase function pmom */
    for (int iu = 0; iu <= input.rte.nstr; iu++) {
      tmp                                            = output->pmom[output->atm.nlyr - lu - 1][0][iu];
      output->pmom[output->atm.nlyr - lu - 1][0][iu] = output->pmom[lu][0][iu];
      output->pmom[lu][0][iu]                        = tmp;
    }
  }
  return 0;
}  // end of reverse_profiles ()

/* allocate memory for Raman scattering components needed for calculation of */
/* source function for second iteration                                      */
raman_qsrc_components* calloc_raman_qsrc_components (int raman_fast, int nzout, int maxumu, int nphi, int nstr, int nlambda) {
  raman_qsrc_components* result = calloc (1, sizeof (raman_qsrc_components));

  if (raman_fast) {
    ASCII_calloc_double (&result->uavgso, nzout, nlambda);
    ASCII_calloc_double (&result->uavgdn, nzout, nlambda);
    ASCII_calloc_double (&result->uavgup, nzout, nlambda);
    ASCII_calloc_double (&result->rfldir, nzout, nlambda);
    ASCII_calloc_double (&result->rfldn, nzout, nlambda);
    ASCII_calloc_double (&result->flup, nzout, nlambda);
    ASCII_calloc_double_3D (&result->u0u, nzout, maxumu, nlambda);
    ASCII_calloc_double_4D (&result->uu, nzout, maxumu, maxumu, nlambda);
  }
  result->fbeam = calloc (nlambda, sizeof (double));
  ASCII_calloc_double (&result->dtauc, nzout, nlambda);
  ASCII_calloc_double_4D (&result->uum, nzout, maxumu, maxumu, nlambda);

  return result;
}  // end of calloc_raman_qsrc_components()

/* free memory for Raman scattering components needed for calculation of */
/* source function for second iteration                                      */
static void free_raman_qsrc_components (raman_qsrc_components* result, int raman_fast, int nzout, int maxumu, int nphi, int nstr) {
  if (raman_fast) {
    if (result->uavgso != NULL)
      ASCII_free_double (result->uavgso, nzout);
    if (result->uavgup != NULL)
      ASCII_free_double (result->uavgup, nzout);
    if (result->uavgdn != NULL)
      ASCII_free_double (result->uavgdn, nzout);
    if (result->rfldir != NULL)
      ASCII_free_double (result->rfldir, nzout);
    if (result->rfldn != NULL)
      ASCII_free_double (result->rfldn, nzout);
    if (result->flup != NULL)
      ASCII_free_double (result->flup, nzout);
    if (result->u0u != NULL)
      ASCII_free_double_3D (result->u0u, nzout, maxumu);
    if (result->uu != NULL)
      ASCII_free_double_4D (result->uu, nzout, maxumu, maxumu);
  }
  if (result->fbeam != NULL)
    free (result->fbeam);
  if (result->dtauc != NULL)
    ASCII_free_double (result->dtauc, nzout);
  if (result->uum != NULL)
    ASCII_free_double_4D (result->uum, nzout, maxumu, maxumu);

  free (result);
}  // end of free_raman_qsrc_components ()

/* allocate temporary memory optical properties */
static save_optprop* calloc_save_optprop (int Nlev) {
  save_optprop* result = calloc (1, sizeof (save_optprop));

  result->tauw = calloc (Nlev, sizeof (float));
  result->taui = calloc (Nlev, sizeof (float));
  result->g1d  = calloc (Nlev, sizeof (float));
  result->g2d  = calloc (Nlev, sizeof (float));
  result->fd   = calloc (Nlev, sizeof (float));
  result->g1i  = calloc (Nlev, sizeof (float));
  result->g2i  = calloc (Nlev, sizeof (float));
  result->fi   = calloc (Nlev, sizeof (float));
  result->ssaw = calloc (Nlev, sizeof (float));
  result->ssai = calloc (Nlev, sizeof (float));

  return result;
}

/* allocate temporary memory for the RTE solvers */
static int calloc_rte_output (input_struct input,
                              int          nzout,
                              int          output3D,
                              int          Nxcld,
                              int          Nycld,
                              int          Nzcld,
                              int          Nxsample,
                              int          Nysample,
                              int          Ncsample,
                              int          Nlambda,
                              int          Nlyr,
                              int*         threed,
                              const size_t N_triangles,
                              rte_output** result) {
  int status = 0;

  *result = calloc (1, sizeof (rte_output));

  (*result)->albmed          = calloc (input.rte.maxumu, sizeof (float));
  (*result)->trnmed          = calloc (input.rte.maxumu, sizeof (float));
  (*result)->dfdt            = calloc (nzout, sizeof (float));
  (*result)->flup            = calloc (nzout, sizeof (float));
  (*result)->rfldir          = calloc (nzout, sizeof (float));
  (*result)->rfldn           = calloc (nzout, sizeof (float));
  (*result)->uavg            = calloc (nzout, sizeof (float));
  (*result)->uavgdn          = calloc (nzout, sizeof (float));
  (*result)->uavgso          = calloc (nzout, sizeof (float));
  (*result)->uavgup          = calloc (nzout, sizeof (float));
  (*result)->heat            = calloc (nzout, sizeof (float));
  (*result)->emis            = calloc (nzout, sizeof (float));
  (*result)->w_zout          = calloc (nzout, sizeof (float));
  (*result)->sslidar_nphot   = calloc (nzout, sizeof (float));
  (*result)->sslidar_nphot_q = calloc (nzout, sizeof (float));
  (*result)->sslidar_ratio   = calloc (nzout, sizeof (float));

  if (input.rte.maxumu > 0) {
    status = ASCII_calloc_float (&((*result)->u0u), nzout, input.rte.maxumu);
    CHKERR (status);
  }

  if (input.rte.maxphi > 0 && input.rte.maxumu > 0) {
    status = ASCII_calloc_float_3D (&((*result)->uu), input.rte.maxphi, nzout, input.rte.maxumu);
    CHKERR (status);
  }

  if (input.rte.solver == SOLVER_DISORT && input.raman) {
    status = ASCII_calloc_float_3D (&((*result)->uum), input.rte.nstr, nzout, input.rte.maxumu);
    CHKERR (status);
  }

  /* 3d fields */
  if (output3D) {
    (*result)->alloc3D = 1;

    status = ASCII_calloc_float_3D (&((*result)->rfldir3d), nzout, Nxsample, Nysample);
    CHKERR (status);
    status = ASCII_calloc_float_3D (&((*result)->rfldn3d), nzout, Nxsample, Nysample);
    CHKERR (status);
    status = ASCII_calloc_float_3D (&((*result)->flup3d), nzout, Nxsample, Nysample);
    CHKERR (status);
    status = ASCII_calloc_float_3D (&((*result)->uavgso3d), nzout, Nxsample, Nysample);
    CHKERR (status);
    status = ASCII_calloc_float_3D (&((*result)->uavgdn3d), nzout, Nxsample, Nysample);
    CHKERR (status);
    status = ASCII_calloc_float_3D (&((*result)->uavgup3d), nzout, Nxsample, Nysample);
    CHKERR (status);

    /* variances */
    if (input.rte.mc.std) {
      status = ASCII_calloc_float_3D (&((*result)->rfldir3d_var), nzout, Nxsample, Nysample);
      CHKERR (status);
      status = ASCII_calloc_float_3D (&((*result)->rfldn3d_var), nzout, Nxsample, Nysample);
      CHKERR (status);
      status = ASCII_calloc_float_3D (&((*result)->flup3d_var), nzout, Nxsample, Nysample);
      CHKERR (status);
      status = ASCII_calloc_float_3D (&((*result)->uavgso3d_var), nzout, Nxsample, Nysample);
      CHKERR (status);
      status = ASCII_calloc_float_3D (&((*result)->uavgdn3d_var), nzout, Nxsample, Nysample);
      CHKERR (status);
      status = ASCII_calloc_float_3D (&((*result)->uavgup3d_var), nzout, Nxsample, Nysample);
      CHKERR (status);
    }

    if (input.output3Dradiances) {
      status = ASCII_calloc_float_4D (&((*result)->radiance3d), nzout, Nxsample, Nysample, input.rte.mc.nstokes);
      CHKERR (status);
      status =
        ASCII_calloc_float_6D (&((*result)->radiance3d_is), nzout, Ncsample, Nxsample, Nysample, input.rte.mc.nstokes, Nlambda);
      CHKERR (status);
      if (input.rte.mc.std) {
        status = ASCII_calloc_float_4D (&((*result)->radiance3d_var), nzout, Nxsample, Nysample, input.rte.mc.nstokes);
        CHKERR (status);
      }
    }

    if (input.rte.mc.spectral_is || input.rte.mc.concentration_is) {
      status = ASCII_calloc_float_5D (&((*result)->fl3d_is), nzout, Ncsample, Nxsample, Nysample, Nlambda);
      CHKERR (status);
    }

    if (input.rte.mc.jacobian[DIM_1D]) {
      /* fprintf(stderr, "calloc jacobian nzout %d  Nxsample %d Nysample %d input.n_caoth %d abs/sca %d Nzcld %d \n", nzout, Nxsample, Nysample, input.n_caoth, 2, Nlyr);  */
      status = ASCII_calloc_float_6D (&((*result)->jacobian), nzout, 1, 1, input.n_caoth + 2, 2, Nlyr);
      CHKERR (status);
    }

    if (input.rte.mc.backward.absorption) {
      status = ASCII_calloc_float_3D (&((*result)->absback3d), nzout, Nxsample, Nysample);
      CHKERR (status);
      if (input.rte.mc.std) {
        status = ASCII_calloc_float_3D (&((*result)->absback3d_var), nzout, Nxsample, Nysample);
        CHKERR (status);
      }
    }

    if (input.rte.mc.absorption != MCFORWARD_ABS_NONE || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d) {
      (*result)->abs3d = calloc_abs3d (Nxcld, Nycld, Nzcld, threed);
      CHKPOINTER ((*result)->abs3d);
      if (input.rte.mc.std) {
        (*result)->abs3d_var = calloc_abs3d (Nxcld, Nycld, Nzcld, threed);
        CHKPOINTER ((*result)->abs3d_var);
      }
    }
  }

  (*result)->triangle_results = NULL;
  status                      = init_triangular_surface_result_struct (N_triangles, &((*result)->triangle_results));
  CHKERR (status);

  return 0;
}  // end of calloc_rte_output()

/* reset rte_output structure */
static int reset_rte_output (rte_output** rte,
                             input_struct input,
                             int          nzout,
                             int          Nxcld,
                             int          Nycld,
                             int          Nzcld,
                             int          Nxsample,
                             int          Nysample,
                             int          Ncsample,
                             int          Nlambda,
                             int          Nlyr,
                             int*         threed,
                             const size_t N_triangles) {
  int status;
  int output3D = (*rte)->alloc3D;

  status = free_rte_output (*rte, input, nzout, Nxcld, Nycld, Nzcld, Nxsample, Nysample, Ncsample, Nlyr, Nlambda, threed);
  CHKERR (status);

  status = calloc_rte_output (input,
                              nzout,
                              output3D,
                              Nxcld,
                              Nycld,
                              Nzcld,
                              Nxsample,
                              Nysample,
                              Ncsample,
                              Nlambda,
                              Nlyr,
                              threed,
                              N_triangles,
                              rte);
  CHKERR (status);

  return 0;
}  // reset_rte_output ()

static int add_rte_output (rte_output*        rte,
                           const rte_output*  add,
                           const double       factor,
                           const double*      factor_spectral,
                           const input_struct input,
                           const int          nzout,
                           const int          Nxcld,
                           const int          Nycld,
                           const int          Nzcld,
                           const int          Nc,
                           const int          Nlyr,
                           const int          Nlambda,
                           const int*         threed,
                           const int          islower,
                           const int          isupper,
                           const int          jslower,
                           const int          jsupper,
                           const int          isstep,
                           const int          jsstep) {
  const double factor2 = factor * factor;

  for (int lu = 0; lu < nzout; lu++) {
    rte->rfldir[lu] += factor * add->rfldir[lu];
    rte->rfldn[lu] += factor * add->rfldn[lu];
    rte->flup[lu] += factor * add->flup[lu];
    rte->uavg[lu] += factor * add->uavg[lu];
    rte->uavgdn[lu] += factor * add->uavgdn[lu];
    rte->uavgso[lu] += factor * add->uavgso[lu];
    rte->uavgup[lu] += factor * add->uavgup[lu];
    rte->dfdt[lu] += factor * add->dfdt[lu];
    rte->heat[lu] += factor * add->heat[lu];
    rte->emis[lu] += factor * add->emis[lu];
    rte->w_zout[lu] += factor * add->w_zout[lu];
    rte->sslidar_nphot[lu] += factor * add->sslidar_nphot[lu];
    rte->sslidar_nphot_q[lu] += factor * add->sslidar_nphot_q[lu];
    rte->sslidar_ratio[lu] += factor * add->sslidar_ratio[lu];

    for (int iu = 0; iu < input.rte.numu; iu++) {
      rte->u0u[lu][iu] += factor * add->u0u[lu][iu];

      for (int j = 0; j < input.rte.nphi; j++)
        rte->uu[j][lu][iu] += factor * add->uu[j][lu][iu];
    }
  }

  for (int iu = 0; iu < input.rte.numu; iu++) {
    rte->albmed[iu] += factor * add->albmed[iu];
    rte->trnmed[iu] += factor * add->trnmed[iu];
  }

  if (rte->alloc3D && add->alloc3D) {
    for (int ks = 0; ks < nzout; ks++) {
      for (int is = islower; is <= isupper; is += isstep) {
        for (int js = jslower; js <= jsupper; js += jsstep) {

          rte->rfldir3d[ks][is][js] += factor * add->rfldir3d[ks][is][js];
          rte->rfldn3d[ks][is][js] += factor * add->rfldn3d[ks][is][js];
          rte->flup3d[ks][is][js] += factor * add->flup3d[ks][is][js];
          rte->uavgso3d[ks][is][js] += factor * add->uavgso3d[ks][is][js];
          rte->uavgdn3d[ks][is][js] += factor * add->uavgdn3d[ks][is][js];
          rte->uavgup3d[ks][is][js] += factor * add->uavgup3d[ks][is][js];

          /* variances */
          if (input.rte.mc.std) {
            rte->rfldir3d_var[ks][is][js] += factor2 * add->rfldir3d_var[ks][is][js];
            rte->rfldn3d_var[ks][is][js] += factor2 * add->rfldn3d_var[ks][is][js];
            rte->flup3d_var[ks][is][js] += factor2 * add->flup3d_var[ks][is][js];
            rte->uavgso3d_var[ks][is][js] += factor2 * add->uavgso3d_var[ks][is][js];
            rte->uavgdn3d_var[ks][is][js] += factor2 * add->uavgdn3d_var[ks][is][js];
            rte->uavgup3d_var[ks][is][js] += factor2 * add->uavgup3d_var[ks][is][js];
          }

          if (input.rte.mc.concentration_is) {
            for (int ic = 0; ic < Nc; ic++)
              rte->fl3d_is[ks][ic][is][js][0] += factor * add->fl3d_is[ks][ic][is][js][0];
          } else if (input.rte.mc.spectral_is) {
            for (int iv_alis = 0; iv_alis < Nlambda; iv_alis++)
              rte->fl3d_is[ks][0][is][js][iv_alis] += factor_spectral[iv_alis] * add->fl3d_is[ks][0][is][js][iv_alis];
          }

          if (input.output3Dradiances) {
            for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
              rte->radiance3d[ks][is][js][ip] += factor * add->radiance3d[ks][is][js][ip];

              /* FIXCE spectral and concentration importance sampling together not */
              /* yet working correctly */
              if (input.rte.mc.concentration_is) {
                for (int ic = 0; ic < Nc; ic++) {
                  rte->radiance3d_is[ks][ic][is][js][ip][0] += factor * add->radiance3d_is[ks][ic][is][js][ip][0];
                }
              }

              if (input.rte.mc.spectral_is) {
                for (int ic = 0; ic < Nc; ic++)
                  for (int iv_alis = 0; iv_alis < Nlambda; iv_alis++)
                    rte->radiance3d_is[ks][ic][is][js][ip][iv_alis] +=
                      factor_spectral[iv_alis] * add->radiance3d_is[ks][ic][is][js][ip][iv_alis];
              }
            }

            if (input.rte.mc.std) {
              for (int ip = 0; ip < input.rte.mc.nstokes; ip++) {
                rte->radiance3d_var[ks][is][js][ip] += factor2 * add->radiance3d_var[ks][is][js][ip];
              }
            }
          }

          if (input.rte.mc.backward.absorption) {
            rte->absback3d[ks][is][js] += factor * add->absback3d[ks][is][js];
            if (input.rte.mc.std) {
              rte->absback3d_var[ks][is][js] += factor2 * add->absback3d_var[ks][is][js];
            }
          }
        }
      }
    }

    if ((input.rte.solver == SOLVER_MONTECARLO && input.rte.mc.absorption != MCFORWARD_ABS_NONE) ||
        input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) {
      for (int ks = 0; ks < Nzcld; ks++) {
        if (threed[ks]) { /* only for 3D layers */
          for (int is = 0; is < Nxcld; is++) {
            for (int js = 0; js < Nycld; js++) {
              rte->abs3d[ks][is][js] += factor * add->abs3d[ks][is][js];
              if (input.rte.mc.std)
                rte->abs3d_var[ks][is][js] += factor2 * add->abs3d_var[ks][is][js];
            }
          }
        }
      }
    }
  }

  if (rte->triangle_results) {
    const int ierr = add_triangular_surface_result (factor, add->triangle_results, rte->triangle_results);
    CHKERR (ierr);
  }

  return 0;
}  // end if add_rte_output ()

/* free temporary memory for the RTE solvers */
static int free_rte_output (rte_output*  result,
                            input_struct input,
                            int          nzout,
                            int          Nxcld,
                            int          Nycld,
                            int          Nzcld,
                            int          Nxsample,
                            int          Nysample,
                            int          Ncsample,
                            int          Nlyr,
                            int          Nlambda,
                            int*         threed) {
  CHKPOINTEROUT (result, "rte_output cannot be freed, it is not allocated!");

  free (result->albmed);
  result->albmed = NULL;
  free (result->trnmed);
  result->trnmed = NULL;
  free (result->dfdt);
  result->dfdt = NULL;
  free (result->flup);
  result->flup = NULL;
  free (result->rfldir);
  result->rfldir = NULL;
  free (result->rfldn);
  result->rfldn = NULL;
  free (result->uavg);
  result->uavg = NULL;
  free (result->uavgdn);
  result->uavgdn = NULL;
  free (result->uavgso);
  result->uavgso = NULL;
  free (result->uavgup);
  result->uavgup = NULL;
  free (result->heat);
  result->heat = NULL;
  free (result->emis);
  result->emis = NULL;
  free (result->w_zout);
  result->w_zout = NULL;
  free (result->sslidar_nphot);
  result->sslidar_nphot = NULL;
  free (result->sslidar_nphot_q);
  result->sslidar_nphot_q = NULL;
  free (result->sslidar_ratio);
  result->sslidar_ratio = NULL;

  if (result->u0u != NULL) {
    ASCII_free_float (result->u0u, nzout);
    result->u0u = NULL;
  }

  if (result->uu != NULL) {
    ASCII_free_float_3D (result->uu, input.rte.nphi, nzout);
    result->uu = NULL;
  }

  if (result->uum != NULL) {
    ASCII_free_float_3D (result->uum, input.rte.nstr, nzout);
    result->uum = NULL;
  }

  /* free 3d fields */
  if (result->alloc3D) {

    if (result->rfldir3d != NULL) {
      ASCII_free_float_3D (result->rfldir3d, nzout, Nxsample);
      result->rfldir3d = NULL;
    }

    if (result->rfldn3d != NULL) {
      ASCII_free_float_3D (result->rfldn3d, nzout, Nxsample);
      result->rfldn3d = NULL;
    }

    if (result->flup3d != NULL) {
      ASCII_free_float_3D (result->flup3d, nzout, Nxsample);
      result->flup3d = NULL;
    }

    if (result->uavgso3d != NULL) {
      ASCII_free_float_3D (result->uavgso3d, nzout, Nxsample);
      result->uavgso3d = NULL;
    }

    if (result->uavgdn3d != NULL) {
      ASCII_free_float_3D (result->uavgdn3d, nzout, Nxsample);
      result->uavgdn3d = NULL;
    }

    if (result->uavgup3d != NULL) {
      ASCII_free_float_3D (result->uavgup3d, nzout, Nxsample);
      result->uavgup3d = NULL;
    }

    if (result->radiance3d != NULL) {
      ASCII_free_float_4D (result->radiance3d, nzout, Nxsample, Nysample);
      result->radiance3d = NULL;
    }

    if (input.rte.mc.spectral_is || input.rte.mc.concentration_is)
      if (result->fl3d_is != NULL) {
        ASCII_free_float_5D (result->fl3d_is, nzout, Ncsample, Nxsample, Nysample);
        result->fl3d_is = NULL;
      }

    if (result->radiance3d_is != NULL) {
      ASCII_free_float_6D (result->radiance3d_is, nzout, Ncsample, Nxsample, Nysample, input.rte.mc.nstokes);
      result->radiance3d_is = NULL;
    }

    if (input.rte.mc.jacobian[DIM_1D])
      if (result->jacobian != NULL) {
        ASCII_free_float_6D (result->jacobian, nzout, Nxsample, Nysample, input.n_caoth + 2, 2);
        result->jacobian = NULL;
      }

    if (input.rte.mc.backward.absorption)
      if (result->absback3d != NULL) {
        ASCII_free_float_3D (result->absback3d, nzout, Nxsample);
        result->absback3d = NULL;
      }

    if (input.rte.mc.absorption != MCFORWARD_ABS_NONE || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d)
      if (result->abs3d != NULL) {
        free_abs3d (result->abs3d, Nxcld, Nycld, Nzcld, threed);
        result->abs3d = NULL;
      }

    /* variances */
    if (input.rte.mc.std) {

      if (result->rfldir3d_var != NULL) {
        ASCII_free_float_3D (result->rfldir3d_var, nzout, Nxsample);
        result->rfldir3d_var = NULL;
      }
      if (result->rfldn3d_var != NULL) {
        ASCII_free_float_3D (result->rfldn3d_var, nzout, Nxsample);
        result->rfldn3d_var = NULL;
      }
      if (result->flup3d_var != NULL) {
        ASCII_free_float_3D (result->flup3d_var, nzout, Nxsample);
        result->flup3d_var = NULL;
      }
      if (result->uavgso3d_var != NULL) {
        ASCII_free_float_3D (result->uavgso3d_var, nzout, Nxsample);
        result->uavgso3d_var = NULL;
      }
      if (result->uavgdn3d_var != NULL) {
        ASCII_free_float_3D (result->uavgdn3d_var, nzout, Nxsample);
        result->uavgdn3d_var = NULL;
      }
      if (result->uavgup3d_var != NULL) {
        ASCII_free_float_3D (result->uavgup3d_var, nzout, Nxsample);
        result->uavgup3d_var = NULL;
      }
      if (result->radiance3d_var != NULL) {
        ASCII_free_float_4D (result->radiance3d_var, nzout, Nxsample, Nysample);
        result->radiance3d_var = NULL;
      }

      if (input.rte.mc.backward.absorption)
        if (result->absback3d_var != NULL) {
          ASCII_free_float_3D (result->absback3d_var, nzout, Nxsample);
          result->absback3d_var = NULL;
        }

      if (input.rte.mc.absorption != MCFORWARD_ABS_NONE || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM || input.ipa3d)
        if (result->abs3d_var != NULL) {
          free_abs3d (result->abs3d_var, Nxcld, Nycld, Nzcld, threed);
          result->abs3d_var = NULL;
        }
    }
  }

  const int ierr = free_triangular_surface_result_struct (result->triangle_results);
  CHKERR (ierr);

  free (result);
  return 0;
}  // end of free_rte_output ()

/******************************************************/
/* setup optical properties, do some initializations, */
/* and call the RTE solver.                           */
/******************************************************/

static int setup_and_call_solver (input_struct           input,
                                  output_struct*         output,
                                  rte_input*             rte_in,
                                  rte_output*            rte_out,
                                  raman_qsrc_components* raman_qsrc_components,
                                  int                    iv,
                                  int                    ib,
                                  int                    ir,
                                  int*                   threed,
                                  int                    mc_loaddata) {
  int        status = 0, iv1 = 0, iv2 = 0;
  static int first                   = 1;
  static int verbose                 = 0;
  int        rte_solver              = NOT_DEFINED_INTEGER;
  int        skip_optical_properties = FALSE;

  /* change rte solver to NULL solver for                                   */
  /* solar simulations with sza > 90 degrees done by plane paralell solvers */

  rte_solver = input.rte.solver;

  switch (input.source) {
  case SRC_THERMAL:
    /* do not change solver */
    break;
  case SRC_SOLAR:
  case SRC_BLITZ: /* BCA */
  case SRC_LIDAR: /* BCA */
    if (output->atm.sza_r[iv] >= 90.0) {
      switch (input.rte.solver) {

      /* plane parallel solvers */
      case SOLVER_FDISORT1:
      case SOLVER_FDISORT2:
      case SOLVER_RODENTS:
      case SOLVER_SCHWARZSCHILD:
      case SOLVER_TWOSTREBE:
      case SOLVER_TWOMAXRND:
      case SOLVER_TWOMAXRND3C:
      case SOLVER_DYNAMIC_TWOSTREAM:
      case SOLVER_DYNAMIC_TENSTREAM:
      case SOLVER_SSS:
      case SOLVER_SSSI:
      case SOLVER_DISORT:  // aky 24022011, cdisort may produce results for sza>90 if intensity correction
                           // are turned off. But results are a bit dubious if compared with mystic.
                           /* change solver to SOLVER_NULL, as output is 0 anyway */
        /* It is only 0 if pseudospherical is off aky 05022014 */
        if (!input.rte.pseudospherical) {
          rte_solver              = SOLVER_NULL;
          skip_optical_properties = TRUE;
          if (input.verbose)
            fprintf (stderr, " ... solar calculation with SZA>90 and pp solver, switch solver to NULL-solver \n");
        }
        break;

      /* spherical solvers */
      case SOLVER_SDISORT:
      case SOLVER_SPSDISORT:
      case SOLVER_FTWOSTR:
      case SOLVER_TWOSTR:
      case SOLVER_SOS:
        /* do not change solver as these solvers should be   */
        /* able to do solar calculations for sza > 90 degree */
        break;

      /* special solvers */
      case SOLVER_MONTECARLO: /* user must know it */
      case SOLVER_TZS:        /* only thermal anyway */
      case SOLVER_SSLIDAR:
        /* do nothing */
        break;

      case SOLVER_NULL:
        break;
      default:
        fprintf (stderr, "Error, unknown RTE solver %d\n", input.rte.solver);
        break;
      }
    }
    break;
  default:
    fprintf (stderr, "Error, unknown source %d\n", input.source);
    return -1;
  }  // switch (input.source)

  if (input.raman) {
    /* Optical properties are calculated just before calling qdisort in solve_rte.c */
    skip_optical_properties = TRUE;
    /* For raman_fast optical properties are calculated as usual */
    if (input.raman_fast && ir == 0)
      skip_optical_properties = FALSE;
  }

  verbose = input.verbose;

  if (input.verbose && input.ck_scheme == CK_LOWTRAN)
    verbose = 0; /* no verbose output for the following call of optical properties */

  /* calculate optical properties from model input */
  if (input.raman_fast) {
    iv1 = ib;
    iv2 = ib;
  } else {
    iv1 = iv;
    iv2 = iv;
  }

  if (input.rte.mc.spectral_is) {
    if (ib == 0) {
      /* fprintf(stderr, "alloc ALIS struct iv %d\n", iv); */
      output->mc.alis.dt = calloc (output->mc.alis.nlambda_abs, sizeof (double**));
      output->mc.alis.om = calloc (output->mc.alis.nlambda_abs, sizeof (double**));
      for (int iv_alis = 0; iv_alis < output->mc.alis.nlambda_abs; iv_alis++) {
        output->mc.alis.dt[iv_alis] = calloc (input.n_caoth + 2, sizeof (double*));
        output->mc.alis.om[iv_alis] = calloc (input.n_caoth + 2, sizeof (double*));
      }
      for (int iv_alis = 0; iv_alis < output->mc.alis.nlambda_abs; iv_alis++) {
        for (int isp = 0; isp < input.n_caoth + 2; isp++) {
          output->mc.alis.dt[iv_alis][isp] = calloc (output->atm.nlev_common - 1, sizeof (double));
          output->mc.alis.om[iv_alis][isp] = calloc (output->atm.nlev_common - 1, sizeof (double));
        }
      }
    }

    for (int iv_alis = 0; iv_alis < output->mc.alis.nlambda_abs; iv_alis++) {
      status = optical_properties (input, output, 0.0, ir, iv_alis, iv_alis, ib, verbose, skip_optical_properties);
      CHKERR (status);
    }

    /* spatially constant spectral Lambertian albedo */
    output->mc.alis.albedo = calloc (output->mc.alis.nlambda_abs, sizeof (double));
    for (int iv_alis = 0; iv_alis < output->mc.alis.nlambda_abs; iv_alis++)
      output->mc.alis.albedo[iv_alis] = output->alb.albedo_r[iv_alis];

    /* 2D spectral Lambertian albedo */
    if (output->surfaces.n > 0 && output->surfaces.albedo_r != NULL) {
      output->mc.alis.alb_type = calloc (output->surfaces.n, sizeof (double*));

      for (int il = 0; il < output->surfaces.n; il++) {
        output->mc.alis.alb_type[il] = calloc (output->mc.alis.nlambda_abs, sizeof (double));

        for (int iv_alis = 0; iv_alis < output->mc.alis.nlambda_abs; iv_alis++)
          output->mc.alis.alb_type[il][iv_alis] = output->surfaces.albedo_r[il][iv_alis];
      }
    }

    if (output->mc.alis.dt[iv][0][output->atm.nlev_common - 2] * (1.0 - output->mc.alis.om[iv][0][output->atm.nlev_common - 2]) >
        0.5) {
      fprintf (stderr, " ... Warning: Absorption is very high at the calculation wavelength for\n");
      fprintf (stderr, " ...          spectral importance sampling. In order to improve the result \n");
      fprintf (stderr, " ...          please select another calculation wavelength using the option\n");
      fprintf (stderr, " ...          *mc_spectral_is_wvl*. \n");
      fprintf (stderr,
               " ...          (Absorption optical depth in lowest layer is %g) \n",
               output->mc.alis.dt[iv][0][output->atm.nlev_common - 2] *
                 (1.0 - output->mc.alis.om[iv][0][output->atm.nlev_common - 2]));
    }
  }

  //TODO: else statement? call optical properties twice -> double allocating etc
  status = optical_properties (input, output, 0.0, ir, iv1, iv2, ib, verbose, skip_optical_properties); /* in ancillary.c */
  CHKERR (status);

  if (input.ck_scheme == CK_LOWTRAN) {

    /* this is a little inefficient as we need first the scattering optical */
    /* properties to calculate the absorption coefficients and then         */
    /* recalculate the optical properties; the reason is that SBDART        */
    /* requires the total scattering cross section as input, and this       */
    /* quantity is only available after the call to optical_properties      */

    /* ??? where to get the mixing ratio of N2 from ??? */
    /* ??? setting this value to -1                 ??? */
    status = sbdart_profile (output->atm.microphys.temper,
                             output->atm.microphys.press,
                             output->atm.zd,
                             output->atm.microphys.dens[MOL_H2O],
                             output->atm.microphys.dens[MOL_O3],
                             -1.0,
                             output->mixing_ratio[MX_O2],
                             output->mixing_ratio[MX_CO2],
                             output->mixing_ratio[MX_CH4],
                             output->mixing_ratio[MX_N2O],
                             output->mixing_ratio[MX_NO2],
                             input.ck_abs[CK_ABS_O4],
                             input.ck_abs[CK_ABS_N2],
                             input.ck_abs[CK_ABS_CO],
                             input.ck_abs[CK_ABS_SO2],
                             input.ck_abs[CK_ABS_NH3],
                             input.ck_abs[CK_ABS_NO],
                             input.ck_abs[CK_ABS_HNO3],
                             output->dtauc,
                             output->ssalb,
                             output->atm.nlev,
                             output->atm.sza_r[iv],
                             output->wl.lambda_r[iv],
                             iv,
                             output->crs_ck.profile[0],
                             input.rte.mc.spectral_is);
    CHKERR (status);

    /* allocate memory for absorption coefficient profile */
    if (first) {
      status = ASCII_calloc_float (&output->kabs.kabs, output->atm.nlev, output->wl.nlambda_r);
      CHKERR (status);

      first = 0;
    }

    /* set absorption coefficient */
    for (int lu = 0; lu < output->atm.nlyr; lu++)
      output->kabs.kabs[lu][iv] = output->crs_ck.profile[0][iv].crs[0][0][lu][ib];

    if (input.rte.mc.spectral_is) {
      for (int iv_alis = 0; iv_alis < output->mc.alis.nlambda_abs; iv_alis++) {
        status = optical_properties (input, output, 0.0, ir, iv_alis, iv_alis, ib, verbose, skip_optical_properties);
        CHKERR (status);
      }
    }

    /* calculate optical properties from model input */
    status = optical_properties (input, output, 0.0, ir, iv, iv, ib, input.verbose, skip_optical_properties); /* in ancillaries.c */
    CHKERR (status);
  }

  /* translate output altitudes (zout) to optical depths (utau) */

  switch (input.rte.solver) {
  case SOLVER_FDISORT1:
  case SOLVER_SDISORT:
  case SOLVER_FTWOSTR:
  case SOLVER_SOS:
  case SOLVER_FDISORT2:
  case SOLVER_SPSDISORT:
  case SOLVER_TZS:
  case SOLVER_SSS:
  case SOLVER_SSSI:
    /* "very old" version, recycled */
    F77_FUNC (setout, SETOUT)
    (output->dtauc, &(output->atm.nlyr), &(output->atm.nzout), rte_in->utau, output->atm.zd, output->atm.zout_sur);

    break;
  case SOLVER_TWOSTR:
  case SOLVER_SCHWARZSCHILD:
  case SOLVER_TWOSTREBE:
  case SOLVER_TWOMAXRND:
  case SOLVER_TWOMAXRND3C:
  case SOLVER_DYNAMIC_TWOSTREAM:
  case SOLVER_DYNAMIC_TENSTREAM:
  case SOLVER_RODENTS:
  case SOLVER_DISORT:
    /* "old" version */
    /*
    status = set_out (output->atm.zd, output->dtauc, output->atm.nlyr, 
		      input.atm.zout, output->atm.nzout, rte_in->utau);

    if (status!=0) {
      fprintf (stderr, "Error %d returned by set_out()\n", status);
      return status;
    }
    */

    status = c_setout (output->dtauc, output->atm.nlyr, output->atm.nzout, rte_in->utau, output->atm.zd, output->atm.zout_sur);
    CHKERR (status);
    break;
  case SOLVER_MONTECARLO:
  case SOLVER_SSLIDAR:
  case SOLVER_NULL:
    break;
  default:
    fprintf (stderr, "Error, unknown rte_solver %d\n", input.rte.solver);
    CHKERR (-1);
  }

  /* reverse profiles, if required */
  if (input.rte.reverse) {
    status = reverse_profiles (input, output);
    CHKERR (status);
  }

  /**** Switch to SOLVER_NULL to run optical_properties tests without solving RTE, Bettina Richter ****/
  if (input.write_optical_properties) {
    fprintf (stderr, " ... switch rte_solver to null_solver\n");
    input.rte.solver = SOLVER_NULL;
    rte_solver       = SOLVER_NULL;
  }

  /* call RTE solver */
  status =
    call_solver (input, output, rte_solver, rte_in, raman_qsrc_components, iv, ib, ir, rte_out, threed, mc_loaddata, input.verbose);
  CHKERR (status);

  return 0;
}  // end of setup_and_call_solver ()

static int call_solver (input_struct           input,
                        output_struct*         output,
                        int                    rte_solver,
                        rte_input*             rte_in,
                        raman_qsrc_components* raman_qsrc_comp,
                        int                    iv,
                        int                    ib,
                        int                    ir,
                        rte_output*            rte_out,
                        int*                   threed,
                        int                    mc_loaddata,
                        int                    verbose) {
  int                    status = 0;
  int                    ivi    = 0;
  double                 start = 0, end = 0, last = 0;
  raman_qsrc_components* tmp_raman_qsrc_comp = NULL;

  double*** tmp_crs = NULL;

  /* CE: commented since I have introduced the option earth_radius, default value is 6370.0 */
  /* float radius= 6370.0; */ /* Earth's radius in km */

  /* c_disort and c_twostr double declarations as they expect all input in double */
  /* and f77 disort and twostr is all float. AK 23.09.2010                        */

  disort_state  ds_in, twostr_ds;
  disort_output ds_out, twostr_out;
  double*       c_twostr_gg = NULL;
  double*       c_zd        = NULL;
  double        c_r_earth   = (double) input.r_earth;
  int           maz         = 0;

  float*  twostr_gg         = NULL;
  float*  twostr_gg_clr     = NULL;
  float*  twostr_gg_cldk    = NULL;
  float*  twostr_gg_cldn    = NULL;
  float*  twostr_cf         = NULL;
  float*  twostr_ff         = NULL;
  float*  sdisort_beta      = NULL;
  float*  sdisort_sig       = NULL;
  float** sdisort_denstab   = NULL;
  float*  tosdisort_denstab = NULL;
  float * disort_pmom = NULL, *disort2_pmom = NULL, *sss_pmom = NULL, *disort2_phaso = NULL;
  int     disort2_ntheta_default   = 0;
  int*    disort2_ntheta           = &disort2_ntheta_default;
  double* disort2_mup              = NULL;
  int     intensity_correction     = TRUE;
  int     old_intensity_correction = FALSE;
  int     rodents_delta_method     = 0;

  float qwanted_wvl[1];
  float qfbeam[1];
  float qalbedo[1];
  int   iv1 = 0, iv2 = 0;
  int   skip_optical_properties = FALSE;
  int   planck_tempoff          = 0;

  double** phase_back = NULL;

#if HAVE_SOS
  float** pmom_sos = NULL;
#endif

  /* temporary Fortran arrays */
  float *disort_u0u = NULL, *disort_uu = NULL;
  float *sss_u0u = NULL, *sss_uu = NULL;

#if HAVE_MYSTIC
  int      source          = 0;
  long int thermal_photons = 0;

  /* rpv arrays */
  float* alb_type  = NULL;
  float* rpv_rho0  = NULL;
  float* rpv_k     = NULL;
  float* rpv_theta = NULL;
  float* rpv_scale = NULL;
  float* rpv_sigma = NULL;
  float* rpv_t1    = NULL;
  float* rpv_t2    = NULL;

  float* rossli_iso = NULL;
  float* rossli_vol = NULL;
  float* rossli_geo = NULL;

  float* hapke_h  = NULL;
  float* hapke_b0 = NULL;
  float* hapke_w  = NULL;

  int write_files = 0;

  /* write MYSTIC monochromatic output only if monochromatic, non-ck uvspec calculation */
  write_files = (output->wl.nlambda_r * output->atm.nq_r[output->wl.iv_lower] > 1 ? 0 : 1);

  /* RPB very dirty trick, please dont kill me for this!!! */
  if (output->mc.sample.LidarLocEst)
    write_files = 1;
  /* write files for spectral/concentration importance sampling, but only if not spectrally post-processed */
  if ((input.rte.mc.spectral_is || input.rte.mc.concentration_is) && input.processing == PROCESS_NONE)
    write_files = 1;
#endif

  if (input.rte.mc.backward.writeback == 1)
    write_files = 1;

  if (verbose)
    start = clock();

  switch (rte_solver) { /* this is the ONLY use of rte_solver, everywhere else still input.rte.solver is used */

  case SOLVER_MONTECARLO:
#if HAVE_MYSTIC

    if (ib == 0)
      if (input.verbose)
        fprintf (stderr, " ... start Monte Carlo simulation for lambda = %10.2f nm \n", output->wl.lambda_r[iv]);

    /* create rpv arrays */
    /* BCA: this is not nice, especially case surf.n=0 and il=0, affects mystic.c and albedo.c */
    if (output->surfaces.n > 0 && output->surfaces.rpv != NULL) {

      rpv_rho0  = calloc (output->surfaces.n, sizeof (float));
      rpv_k     = calloc (output->surfaces.n, sizeof (float));
      rpv_theta = calloc (output->surfaces.n, sizeof (float));
      rpv_scale = calloc (output->surfaces.n, sizeof (float));
      rpv_sigma = calloc (output->surfaces.n, sizeof (float));
      rpv_t1    = calloc (output->surfaces.n, sizeof (float));
      rpv_t2    = calloc (output->surfaces.n, sizeof (float));

      for (int il = 0; il < output->surfaces.n; il++) {
        rpv_rho0[il]  = output->surfaces.rpv[il].rho0_r[iv];
        rpv_k[il]     = output->surfaces.rpv[il].k_r[iv];
        rpv_theta[il] = output->surfaces.rpv[il].theta_r[iv];
        rpv_scale[il] = output->surfaces.rpv[il].scale_r[iv];
        rpv_sigma[il] = output->surfaces.rpv[il].sigma_r[iv];
        rpv_t1[il]    = output->surfaces.rpv[il].t1_r[iv];
        rpv_t2[il]    = output->surfaces.rpv[il].t2_r[iv];
      }
    } else {
      rpv_rho0  = calloc (1, sizeof (float));
      rpv_k     = calloc (1, sizeof (float));
      rpv_theta = calloc (1, sizeof (float));
      rpv_scale = calloc (1, sizeof (float));
      rpv_sigma = calloc (1, sizeof (float));
      rpv_t1    = calloc (1, sizeof (float));
      rpv_t2    = calloc (1, sizeof (float));

      rpv_rho0[0]  = output->rpv.rho0_r[iv];
      rpv_k[0]     = output->rpv.k_r[iv];
      rpv_theta[0] = output->rpv.theta_r[iv];
      rpv_scale[0] = output->rpv.scale_r[iv];
      rpv_sigma[0] = output->rpv.sigma_r[iv];
      rpv_t1[0]    = output->rpv.t1_r[iv];
      rpv_t2[0]    = output->rpv.t2_r[iv];
    }

    if (output->surfaces.n > 0 && output->surfaces.albedo_r != NULL) {
      alb_type = calloc (output->surfaces.n, sizeof (float));

      for (int il = 0; il < output->surfaces.n; il++)
        alb_type[il] = output->surfaces.albedo_r[il][iv];
    } else {
      alb_type    = calloc (1, sizeof (float));
      alb_type[0] = output->alb.albedo_r[iv];
    }

    /* for rossli we don't have a surface type yet, but that */
    /* would be straightforward to implement                  */
    if (output->surfaces.n > 0 && output->surfaces.rossli != NULL) {

      rossli_iso = calloc (output->surfaces.n, sizeof (float));
      rossli_vol = calloc (output->surfaces.n, sizeof (float));
      rossli_geo = calloc (output->surfaces.n, sizeof (float));

      for (int il = 0; il < output->surfaces.n; il++) {
        rossli_iso[il] = output->surfaces.rossli[il].iso_r[iv];
        rossli_vol[il] = output->surfaces.rossli[il].vol_r[iv];
        rossli_geo[il] = output->surfaces.rossli[il].geo_r[iv];
      }
    } else {
      rossli_iso = calloc (1, sizeof (float));
      rossli_vol = calloc (1, sizeof (float));
      rossli_geo = calloc (1, sizeof (float));

      rossli_iso[0] = output->rossli.iso_r[iv];
      rossli_vol[0] = output->rossli.vol_r[iv];
      rossli_geo[0] = output->rossli.geo_r[iv];
    }

    /* for hapke we don't have a surface type yet, but that */
    /* would be straightforward to implement                  */
    hapke_w  = calloc (1, sizeof (float));
    hapke_b0 = calloc (1, sizeof (float));
    hapke_h  = calloc (1, sizeof (float));

    /* also, wavelength dependence is missing, but that */
    /* would also be straightforward to implement       */
    hapke_w[0]  = output->hapke.w_r[iv];
    hapke_b0[0] = output->hapke.b0_r[iv];
    hapke_h[0]  = output->hapke.h_r[iv];

    switch (input.source) {
    case SRC_SOLAR: /* solar source */
    case SRC_BLITZ: /* blitz source */
    case SRC_LIDAR: /* lidar source */

      switch (input.source) {
      case SRC_SOLAR: /* solar source */
        source = MCSRC_SOLAR;
        break;
      case SRC_BLITZ: /* blitz source */
        source = MCSRC_BLITZ;
        break;
      case SRC_LIDAR: /* lidar source */
        source = MCSRC_LIDAR;
        break;
      default:
        fprintf (stderr, "Error, source %d should not turn up here!\n", input.source);
        return -1;
      }

      status = mystic (&output->atm.nlyr,
                       threed,
                       input.n_caoth + 2,
                       output->mc.dt,
                       output->mc.om,
                       output->mc.g1,
                       output->mc.g2,
                       output->mc.ff,
                       output->mc.ds,
                       &(output->mc.alis),
                       output->mc.refind,
                       input.r_earth * 1000.0,
                       input.rte.mc.refractive_index_pv,
                       &(output->rayleigh_depol[iv]),
                       output->caoth3d,
                       output->mc.re,
                       output->mc.temper,
                       input.atmosphere3d,
                       output->mc.z,
                       output->mc.momaer,
                       output->mc.nmomaer,
                       output->mc.nthetaaer,
                       output->mc.thetaaer,
                       output->mc.muaer,
                       output->mc.phaseaer,
                       output->mc.nphamataer,
                       &output->alb.albedo_r[iv],
                       alb_type,
                       output->atm.sza_r[iv],
                       output->atm.phi0_r[iv],
                       input.atm.sza_spher,
                       input.atm.phi0_spher,
                       output->mc_photons,
                       &source,
                       &(output->wl.wvnmlo_r[iv]),
                       &(output->wl.wvnmhi_r[iv]),
                       &(output->wl.lambda_r[iv]),
                       iv,
                       output->atm.zout_sur,
                       &(output->atm.nzout),
                       rpv_rho0,
                       rpv_k,
                       rpv_theta,
                       rpv_scale,
                       rpv_sigma,
                       rpv_t1,
                       rpv_t2,
                       hapke_h,
                       hapke_b0,
                       hapke_w,
                       rossli_iso,
                       rossli_vol,
                       rossli_geo,
                       input.rossli.hotspot,
                       output->surfaces.mixed_surface,
                       &(input.cm.param[BRDF_CAM_U10]),
                       &(input.cm.param[BRDF_CAM_PCL]),
                       &(input.cm.param[BRDF_CAM_SAL]),
                       &(input.cm.param[BRDF_CAM_UPHI]),
                       &(input.cm.solar_wind),
                       &(input.bpdf),
                       &(input.rte.mc.tenstream),
                       input.rte.mc.tenstream_options,
                       &(input.rte.mc.nca), /* Carolin Klinger 2019 */
                       input.rte.mc.nca_options,
                       &(input.rte.mc.ipa),
                       &(input.rte.mc.absorption),
                       &(input.rte.mc.backward.thermal_heating_method),
                       &mc_loaddata,
                       &(output->mc.sample),
                       &(output->mc.elev),
                       &(output->mc.triangular_surface),
                       &(output->mc.surftemp),
                       input.rte.mc.filename[FN_MC_BASENAME],
                       input.rte.mc.filename[FN_MC_UMU],
                       input.rte.mc.filename[FN_MC_SUNSHAPE_FILE],
                       input.rte.mc.filename[FN_MC_ALBEDO],
                       input.rte.mc.filename[FN_MC_AMBRALS],
                       input.rte.mc.filename[FN_MC_ROSSLI],
                       input.rte.mc.filename[FN_MC_ALBEDO_TYPE],
                       input.rte.mc.filename[FN_MC_RPV2D_TYPE],
                       input.rte.mc.filename[FN_MC_AMBRALS_TYPE],
                       input.rte.mc.netcdf_surface_files,
                       output->surfaces.label,
                       &(output->surfaces.n),
                       &(input.rte.mc.delta_scaling_mucut), /*TZ ds*/
                       &(input.rte.mc.truncate),
                       &(input.rte.mc.reflectalways),
                       input.quiet,
                       rte_out->rfldir,
                       rte_out->rfldn,
                       rte_out->flup,
                       rte_out->uavgso,
                       rte_out->uavgdn,
                       rte_out->uavgup,
                       rte_out->rfldir3d,
                       rte_out->rfldn3d,
                       rte_out->flup3d,
                       rte_out->fl3d_is,
                       rte_out->uavgso3d,
                       rte_out->uavgdn3d,
                       rte_out->uavgup3d,
                       rte_out->radiance3d,
                       rte_out->absback3d,
                       rte_out->abs3d,
                       rte_out->radiance3d_is,
                       rte_out->jacobian,
                       rte_out->rfldir3d_var,
                       rte_out->rfldn3d_var,
                       rte_out->flup3d_var,
                       rte_out->uavgso3d_var,
                       rte_out->uavgdn3d_var,
                       rte_out->uavgup3d_var,
                       rte_out->radiance3d_var,
                       rte_out->absback3d_var,
                       rte_out->abs3d_var,
                       rte_out->triangle_results,
                       output->atm.Nxcld,
                       output->atm.Nycld,
                       output->atm.Nzcld,
                       output->atm.dxcld,
                       output->atm.dycld,
                       input.filename[FN_PATH],
                       input.filename[FN_RANDOMSTATUS],
                       input.rte.mc.readrandomstatus,
                       input.write_output_as_netcdf,
                       write_files,
                       input.rte.mc.visualize,
                       input.rte.mc.progressbar);

      CHKERR (status);

      break;

    case SRC_THERMAL: /* thermal source */
      /* normal atmospheric + surface emission */
      switch (input.rte.mc.absorption) {
      case MCFORWARD_ABS_NONE:
      case MCFORWARD_ABS_HEATING:
        if (!input.rte.mc.backward.yes) { /*TZ bt*/

          thermal_photons = output->mc_photons;

          /* thermal emission of the atmosphere */
          if (ib == 0 && !input.quiet)
            fprintf (stderr, " ... thermal forward\n");

        } else {                                /*TZ bt*/
          thermal_photons = output->mc_photons; /*TZ bt*/

          /* CE for spectral importance sampling, we run a MC */
          /* calculation at only one wavelengths, even if input */
          /* wavelength is not exactly included in */
          /* molecular_tau_file   */
          if (input.rte.mc.spectral_is)
            thermal_photons = input.rte.mc.photons;

          /* thermal emission into the atmosphere TZ bt*/
          if (ib == 0 && !input.quiet)
            fprintf (stderr, " ... thermal backward"); /*TZ*/
        }                                              /*TZ bt*/

        if (!input.rte.mc.backward.yes) /*TZ bt*/
          source = MCSRC_THERMAL_FORWARD;
        else                               /*TZ bt*/
          source = MCSRC_THERMAL_BACKWARD; /*TZ bt*/

        status = mystic (&output->atm.nlyr,
                         threed,
                         input.n_caoth + 2,
                         output->mc.dt,
                         output->mc.om,
                         output->mc.g1,
                         output->mc.g2,
                         output->mc.ff,
                         output->mc.ds,
                         &(output->mc.alis),
                         output->mc.refind,
                         input.r_earth * 1000.0,
                         input.rte.mc.refractive_index_pv,
                         &(output->rayleigh_depol[iv]),
                         output->caoth3d,
                         output->mc.re,
                         output->mc.temper,
                         input.atmosphere3d,
                         output->mc.z,
                         output->mc.momaer,
                         output->mc.nmomaer,
                         output->mc.nthetaaer,
                         output->mc.thetaaer,
                         output->mc.muaer,
                         output->mc.phaseaer,
                         output->mc.nphamataer,
                         &output->alb.albedo_r[iv],
                         alb_type,
                         output->atm.sza_r[iv],
                         output->atm.phi0_r[iv],
                         input.atm.sza_spher,
                         input.atm.phi0_spher,
                         thermal_photons,
                         &source,
                         &(output->wl.wvnmlo_r[iv]),
                         &(output->wl.wvnmhi_r[iv]),
                         &(output->wl.lambda_r[iv]),
                         iv,
                         output->atm.zout_sur,
                         &(output->atm.nzout),
                         rpv_rho0,
                         rpv_k,
                         rpv_theta,
                         rpv_scale,
                         rpv_sigma,
                         rpv_t1,
                         rpv_t2,
                         hapke_h,
                         hapke_b0,
                         hapke_w,
                         rossli_iso,
                         rossli_vol,
                         rossli_geo,
                         input.rossli.hotspot,
                         output->surfaces.mixed_surface,
                         &(input.cm.param[BRDF_CAM_U10]),
                         &(input.cm.param[BRDF_CAM_PCL]),
                         &(input.cm.param[BRDF_CAM_SAL]),
                         &(input.cm.param[BRDF_CAM_UPHI]),
                         &(input.cm.solar_wind),
                         &(input.bpdf),
                         &(input.rte.mc.tenstream),
                         input.rte.mc.tenstream_options,
                         &(input.rte.mc.nca), /* Carolin Klinger 2019 */
                         input.rte.mc.nca_options,
                         &(input.rte.mc.ipa),
                         &(input.rte.mc.absorption),
                         &(input.rte.mc.backward.thermal_heating_method),
                         &mc_loaddata,
                         &(output->mc.sample),
                         &(output->mc.elev),
                         &(output->mc.triangular_surface),
                         &(output->mc.surftemp),
                         input.rte.mc.filename[FN_MC_BASENAME],
                         input.rte.mc.filename[FN_MC_UMU],
                         input.rte.mc.filename[FN_MC_SUNSHAPE_FILE],
                         input.rte.mc.filename[FN_MC_ALBEDO],
                         input.rte.mc.filename[FN_MC_AMBRALS],
                         input.rte.mc.filename[FN_MC_ROSSLI],
                         input.rte.mc.filename[FN_MC_ALBEDO_TYPE],
                         input.rte.mc.filename[FN_MC_RPV2D_TYPE],
                         input.rte.mc.filename[FN_MC_AMBRALS_TYPE],
                         input.rte.mc.netcdf_surface_files,
                         output->surfaces.label,
                         &(output->surfaces.n),
                         &(input.rte.mc.delta_scaling_mucut), /*TZ ds*/
                         &(input.rte.mc.truncate),
                         &(input.rte.mc.reflectalways),
                         input.quiet,
                         rte_out->rfldir,
                         rte_out->rfldn,
                         rte_out->flup,
                         rte_out->uavgso,
                         rte_out->uavgdn,
                         rte_out->uavgup,
                         rte_out->rfldir3d,
                         rte_out->rfldn3d,
                         rte_out->flup3d,
                         rte_out->fl3d_is,
                         rte_out->uavgso3d,
                         rte_out->uavgdn3d,
                         rte_out->uavgup3d,
                         rte_out->radiance3d,
                         rte_out->absback3d,
                         rte_out->abs3d,
                         rte_out->radiance3d_is,
                         rte_out->jacobian,
                         rte_out->rfldir3d_var,
                         rte_out->rfldn3d_var,
                         rte_out->flup3d_var,
                         rte_out->uavgso3d_var,
                         rte_out->uavgdn3d_var,
                         rte_out->uavgup3d_var,
                         rte_out->radiance3d_var,
                         rte_out->absback3d_var,
                         rte_out->abs3d_var,
                         rte_out->triangle_results,
                         output->atm.Nxcld,
                         output->atm.Nycld,
                         output->atm.Nzcld,
                         output->atm.dxcld,
                         output->atm.dycld,
                         input.filename[FN_PATH],
                         input.filename[FN_RANDOMSTATUS],
                         input.rte.mc.readrandomstatus,
                         input.write_output_as_netcdf,
                         write_files,
                         input.rte.mc.visualize,
                         input.rte.mc.progressbar);
        CHKERR (status);
        break;

      default:
        fprintf (stderr, "Error, unknown absorption type %d\n", input.rte.mc.absorption);
        CHKERR (-1);
      }

      break;

    default:
      CHKERROUT (-1, "unknown source");
    }

    /* free RPV arrays */
    free (alb_type);
    free (rpv_rho0);
    free (rpv_k);
    free (rpv_theta);
    free (rpv_scale);
    free (rpv_sigma);
    free (rpv_t1);
    free (rpv_t2);

    /* free ROSSLI arrays */
    free (rossli_iso);
    free (rossli_vol);
    free (rossli_geo);

    /* free HAPKE arrays */
    free (hapke_h);
    free (hapke_b0);
    free (hapke_w);

    break;

#else
    fprintf (stderr, "Error: MYSTIC solver not included in uvspec build.\n");
    fprintf (stderr, "Error: Please contact bernhard.mayer@lmu.de\n");
    CHKERR (-1);
#endif

  case SOLVER_FDISORT1:
  case SOLVER_FDISORT2:
  case SOLVER_DISORT:

    if (rte_solver != SOLVER_DISORT) {
      disort_u0u = calloc (output->atm.nzout * input.rte.maxumu, sizeof (float));
      disort_uu  = calloc (output->atm.nzout * input.rte.maxumu * input.rte.maxphi, sizeof (float));
    }

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    if (input.rte.solver == SOLVER_FDISORT1) {

      if (iv == output->wl.iv_lower && ib == 0) {
        status = F77_FUNC (dcheck, DCHECK) (&output->atm.nlyr,
                                            &output->atm.nzout,
                                            &input.rte.nstr,
                                            &input.rte.numu,
                                            &input.rte.nphi,
                                            &input.optimize_fortran,
                                            &input.optimize_delta);
        CHKERR (status);
      }

      disort_pmom = c2fortran_3D_float_ary (output->atm.nlyr, 1, input.rte.nstr + 1, output->pmom);

      if (((input.source == SRC_SOLAR) && (rte_in->umu0 > 0)) || (input.source == SRC_THERMAL)) {
        /* no need to call disort otherwise as sun below horizon */
        F77_FUNC (disort, DISORT)
        (&output->atm.nlyr,
         output->dtauc,
         output->ssalb,
         disort_pmom,
         output->atm.microphys.temper[0][0],
         &(output->wl.wvnmlo_r[iv]),
         &(output->wl.wvnmhi_r[iv]),
         &(rte_in->usrtau),
         &output->atm.nzout,
         rte_in->utau,
         &input.rte.nstr,
         &(rte_in->usrang),
         &input.rte.numu,
         input.rte.umu,
         &input.rte.nphi,
         input.rte.phi,
         &(rte_in->ibcnd),
         &(rte_in->fbeam),
         &(rte_in->umu0),
         &output->atm.phi0_r[iv],
         &(rte_in->fisot),
         &(rte_in->lamber),
         &output->alb.albedo_r[iv],
         rte_in->hl,
         &(rte_in->btemp),
         &(rte_in->ttemp),
         &(rte_in->temis),
         &input.rte.deltam,
         &(rte_in->planck),
         &(rte_in->onlyfl),
         &(rte_in->accur),
         rte_in->prndis,
         rte_in->header,
         &output->atm.nlyr,
         &output->atm.nzout,
         &input.rte.maxumu,
         &input.rte.nstr,
         &input.rte.maxphi,
         rte_out->rfldir,
         rte_out->rfldn,
         rte_out->flup,
         rte_out->dfdt,
         rte_out->uavg,
         disort_uu,
         disort_u0u,
         rte_out->albmed,
         rte_out->trnmed,
         rte_out->uavgdn,
         rte_out->uavgso,
         rte_out->uavgup,
         &(input.quiet));

        for (int lev = 0; lev < output->atm.nzout; lev++)
          rte_out->uavg[lev] = rte_out->uavgso[lev] + rte_out->uavgdn[lev] + rte_out->uavgup[lev];
      }
    } else if (input.rte.solver == SOLVER_FDISORT2) {

      if (iv == output->wl.iv_lower && ib == 0) {
        status = F77_FUNC (dcheck, DCHECK) (&output->atm.nlyr,
                                            &output->atm.nzout,
                                            &input.rte.nstr,
                                            &input.rte.numu,
                                            &input.rte.nphi,
                                            &input.optimize_fortran,
                                            &input.optimize_delta);
        CHKERR (status);
      }

      /* BRDF or Lambertian albedo */
      if (input.disort2_brdf != BRDF_NONE)
        rte_in->lamber = 0;
      else
        rte_in->lamber = 1;

      disort2_pmom = c2fortran_3D_float_ary (output->atm.nlyr, 1, output->atm.nmom + 1, output->pmom);

      switch (input.rte.disort_icm) {
      case DISORT_ICM_OFF:
        intensity_correction = FALSE;
        break;
      case DISORT_ICM_MOMENTS:
        intensity_correction     = TRUE;
        old_intensity_correction = TRUE;
        break;
      case DISORT_ICM_PHASE:
        intensity_correction     = TRUE;
        old_intensity_correction = FALSE;
        disort2_ntheta           = &(output->ntheta[0][0]);
        disort2_phaso            = c2fortran_3D_float_ary (output->atm.nlyr, 1, output->ntheta[0][0], output->phase);
        disort2_mup              = c2fortran_3D_double_ary (1, 1, output->ntheta[0][0], output->mu);
        break;
      default:
        fprintf (stderr, "Error: unknown disort_icm %d\n", input.rte.disort_icm);
        fprintf (stderr, "       (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
        return -1;
      }

      if (((input.source == SRC_SOLAR) && (rte_in->umu0 > 0)) || (input.source == SRC_THERMAL)) {
        /* no need to call disort2 otherwise, as sun below horizon */
        F77_FUNC (disort2, DISORT2)
        (&output->atm.nlyr,
         output->dtauc,
         output->ssalb,
         &output->atm.nmom,
         disort2_pmom,
         disort2_ntheta,
         disort2_phaso,
         disort2_mup,
         output->atm.microphys.temper[0][0],
         &output->wl.wvnmlo_r[iv],
         &output->wl.wvnmhi_r[iv],
         &rte_in->usrtau,
         &output->atm.nzout,
         rte_in->utau,
         &input.rte.nstr,
         &rte_in->usrang,
         &input.rte.numu,
         input.rte.umu,
         &input.rte.nphi,
         input.rte.phi,
         &rte_in->ibcnd,
         &rte_in->fbeam,
         &rte_in->umu0,
         &output->atm.phi0_r[iv],
         &rte_in->fisot,
         &rte_in->lamber,
         &output->alb.albedo_r[iv],
         &rte_in->btemp,
         &rte_in->ttemp,
         &rte_in->temis,
         &rte_in->planck,
         &rte_in->onlyfl,
         &rte_in->accur,
         rte_in->prndis2,
         rte_in->header,
         &output->atm.nlyr,
         &output->atm.nzout,
         &input.rte.maxumu,
         &input.rte.maxphi,
         &output->atm.nmom,
         rte_out->rfldir,
         rte_out->rfldn,
         rte_out->flup,
         rte_out->dfdt,
         rte_out->uavg,
         disort_uu,
         disort_u0u,
         rte_out->albmed,
         rte_out->trnmed,
         rte_out->uavgdn,
         rte_out->uavgso,
         rte_out->uavgup,
         &(input.disort2_brdf),
         &(output->rpv.rho0_r[iv]),
         &(output->rpv.k_r[iv]),
         &(output->rpv.theta_r[iv]),
         &(output->rpv.sigma_r[iv]),
         &(output->rpv.t1_r[iv]),
         &(output->rpv.t2_r[iv]),
         &(output->rpv.scale_r[iv]),
         &(output->rossli.iso_r[iv]),
         &(output->rossli.vol_r[iv]),
         &(output->rossli.geo_r[iv]),
         &(input.cm.param[BRDF_CAM_U10]),
         &(input.cm.param[BRDF_CAM_PCL]),
         &(input.cm.param[BRDF_CAM_SAL]),
         &intensity_correction,
         &old_intensity_correction,
         &(input.quiet));
        for (int lev = 0; lev < output->atm.nzout; lev++)
          rte_out->uavg[lev] = rte_out->uavgso[lev] + rte_out->uavgdn[lev] + rte_out->uavgup[lev];
      }
    } else if (input.rte.solver == SOLVER_DISORT) {

      if (input.flu.source != NOT_DEFINED_INTEGER) {  // We have included fluorescence as radiation source
        rte_in->fbeam = output->wl.fbeam[iv];         // Need to do calibrated simulation.
      }

      if (input.raman) {

        rte_in->gsrc              = 0;
        ds_in.flag.general_source = FALSE; /* No extra source term for zero Raman scattering. */

        if (ir == 0 && input.raman_fast) {
          rte_in->fbeam = output->wl.fbeam[ib];
        } else if (ir == 0 && !input.raman_fast) {

          if (ib == 0) {
            status = ASCII_calloc_double_3D (&tmp_crs, output->atm.nlev, output->crs.number_of_ramanwavelengths, 3);
            CHKERR (status);

            /* First time for each user wavelength calculate Raman shifted wavelengths and Raman cross section */
            for (int lu = 0; lu < output->atm.nlev; lu++) {
              ivi    = 0;
              status = crs_raman_N2 (output->wl.lambda_r[iv],
                                     input.n_raman_transitions_N2,
                                     output->atm.microphys.temper[lu][0][0],
                                     &tmp_crs[lu],
                                     &ivi,
                                     verbose);
              CHKERR (status);
              status = crs_raman_O2 (output->wl.lambda_r[iv],
                                     input.n_raman_transitions_O2,
                                     output->atm.microphys.temper[lu][0][0],
                                     &tmp_crs[lu],
                                     &ivi,
                                     verbose);
              CHKERR (status);

              /* Put the wanted wavelength last in the tmp_crs by doing the following */
              /* and set correctly after sort                                         */
              tmp_crs[lu][output->crs.number_of_ramanwavelengths - 1][0] = 999e+9;

              /* sort cross section data in ascending order */
              status = ASCII_sortarray (tmp_crs[lu], output->crs.number_of_ramanwavelengths, 3, 0, 0);
              CHKERR (status);

              for (ivi = 0; ivi < output->crs.number_of_ramanshifts; ivi++) {
                output->crs.wvl_of_ramanshifts[ivi] = tmp_crs[0][ivi][0];
                output->crs.crs_raman_RL[lu][ivi]   = tmp_crs[lu][ivi][1];
                output->crs.crs_raman_RG[lu][ivi]   = tmp_crs[lu][ivi][2];
              }
              ivi                                 = output->crs.number_of_ramanwavelengths - 1;
              output->crs.wvl_of_ramanshifts[ivi] = output->wl.lambda_r[iv];
              output->crs.crs_raman_RL[lu][ivi]   = 0.0;
              output->crs.crs_raman_RG[lu][ivi]   = 0.0;
            }
            if (tmp_crs != NULL)
              ASCII_free_double_3D (tmp_crs, output->atm.nlev, output->crs.number_of_ramanwavelengths);
          }

          /* Interpolate all optical quantities to the Raman shifted wavelength */

          qwanted_wvl[0] = (float) output->crs.wvl_of_ramanshifts[ib];
          status =
            arb_wvn (output->wl.nlambda_r, output->wl.lambda_r, output->wl.fbeam, 1, qwanted_wvl, qfbeam, INTERP_METHOD_LINEAR, 0);
          CHKERR (status);

          /* FIXME, redistribute photons instead of interpolating spectrum */
          rte_in->fbeam = qfbeam[0];
          status        = arb_wvn (output->wl.nlambda_r,
                            output->wl.lambda_r,
                            output->alb.albedo_r,
                            1,
                            qwanted_wvl,
                            qalbedo,
                            INTERP_METHOD_LINEAR,
                            0);
          CHKERR (status);

          /* Find iv indices closest to ib wavelength above and below */
          iv1 = closest_above (qwanted_wvl[0], output->wl.lambda_r, output->wl.nlambda_r);
          iv2 = closest_below (qwanted_wvl[0], output->wl.lambda_r, output->wl.nlambda_r);
          /* Calculate optical properties for closest wavelength above and below wanted wavelength,*/
          /* and interpolate optical properties to wanted wavelength.       */
          status = optical_properties (input,
                                       output,
                                       qwanted_wvl[0],
                                       ir,
                                       iv1,
                                       iv2,
                                       0,
                                       verbose,
                                       skip_optical_properties); /* in ancillary.c */
          CHKERR (status);

          /* We also need utau at the "new" interpolated optical depth. */
          /* "very old" version, recycled */
          /* zd is altitude which defines  dtauc */
          /* zout is the altitudes at which tau is wanted */
          F77_FUNC (setout, SETOUT)
          (output->dtauc, &(output->atm.nlyr), &(output->atm.nzout), rte_in->utau, output->atm.zd, output->atm.zout_sur);
        }

        else if (ir == 1) {

          rte_in->fbeam             = 0.0;
          rte_in->gsrc              = 1;
          ds_in.flag.general_source = TRUE; /* Include extra source term for first order Raman scattering. */

          status = ASCII_calloc_double_3D (&(rte_in->qsrc), input.rte.nstr, output->atm.nzout, input.rte.nstr);
          CHKERR (status);

          if (rte_in->usrang) {
            status = ASCII_calloc_double_3D (&(rte_in->qsrcu), input.rte.nstr, output->atm.nzout, input.rte.numu);
            CHKERR (status);
          }

          if (input.raman_fast) {

            /* Calculate Raman crs at wanted wavelength. */
            status = ASCII_calloc_double_3D (&tmp_crs, output->atm.nlev, output->crs.number_of_ramanwavelengths, 3);
            CHKERR (status);

            /* For each user wavelength calculate Raman shifted wavelengths and Raman cross section */
            for (int lu = 0; lu < output->atm.nlev; lu++) {
              ivi    = 0;
              status = crs_raman_N2 (output->wl.lambda_r[iv + output->wl.iv_lower],
                                     input.n_raman_transitions_N2,
                                     output->atm.microphys.temper[0][0][lu],
                                     &tmp_crs[lu],
                                     &ivi,
                                     verbose);
              CHKERR (status);
              status = crs_raman_O2 (output->wl.lambda_r[iv + output->wl.iv_lower],
                                     input.n_raman_transitions_O2,
                                     output->atm.microphys.temper[0][0][lu],
                                     &tmp_crs[lu],
                                     &ivi,
                                     verbose);
              CHKERR (status);

              /* Put the wanted wavelength last in the tmp_crs by doing the following */
              /* and set correctly after sort                                         */
              tmp_crs[lu][output->crs.number_of_ramanwavelengths - 1][0] = 999e+9;

              /* sort cross section data in ascending order */
              status = ASCII_sortarray (tmp_crs[lu], output->crs.number_of_ramanwavelengths, 3, 0, 0);
              CHKERR (status);

              for (ivi = 0; ivi < output->crs.number_of_ramanshifts; ivi++) {
                output->crs.wvl_of_ramanshifts[ivi] = tmp_crs[0][ivi][0];
                output->crs.crs_raman_RL[lu][ivi]   = tmp_crs[lu][ivi][1];
                // Cross sections are reversed in wavelength. This because photons are scattered into the
                // wavelength of interest, lambda_1, from these wavelengths (lambda_j). If lambda_j >
                // lambda_1, then the photon gains energy. Thus, the original cross section indices
                // must be reversed to account for this, and vice versa. AK, 20130513.
                output->crs.crs_raman_RG[lu][output->crs.number_of_ramanshifts - ivi] = tmp_crs[lu][ivi][2];
              }
              ivi                                 = output->crs.number_of_ramanwavelengths - 1;
              output->crs.wvl_of_ramanshifts[ivi] = output->wl.lambda_r[iv];
              output->crs.crs_raman_RL[lu][ivi]   = 0.0;
              output->crs.crs_raman_RG[lu][ivi]   = 0.0;
            }
            if (tmp_crs != NULL)
              ASCII_free_double_3D (tmp_crs, output->atm.nlev, output->crs.number_of_ramanwavelengths);

            /* Interpolate raman_qsrc_comp at wl.lambda_r resolution to raman_wavelength grid */

            tmp_raman_qsrc_comp = calloc_raman_qsrc_components (input.raman_fast,
                                                                output->atm.nzout,
                                                                input.rte.maxumu,
                                                                input.rte.nphi,
                                                                input.rte.nstr,
                                                                output->crs.number_of_ramanwavelengths);

            double tmp_out_int[output->crs.number_of_ramanshifts + 1];
            double tmp_in_int[output->wl.nlambda_r];
            double tmp_wl[output->wl.nlambda_r];

            for (ivi = 0; ivi < output->wl.nlambda_r; ivi++)
              tmp_wl[ivi] = output->wl.lambda_r[ivi];

            /* Interpolate dtauc */
            for (int lu = 0; lu < output->atm.nzout - 1; lu++) {

              for (ivi = 0; ivi < output->wl.nlambda_r; ivi++)
                tmp_in_int[ivi] = raman_qsrc_comp->dtauc[lu][ivi];

              status = arb_wvn_double (output->wl.nlambda_r,
                                       tmp_wl,
                                       tmp_in_int,
                                       output->crs.number_of_ramanshifts,
                                       output->crs.wvl_of_ramanshifts,
                                       tmp_out_int,
                                       INTERP_METHOD_LINEAR,
                                       0);
              CHKERR (status);

              for (ivi = 0; ivi < output->crs.number_of_ramanshifts; ivi++)
                tmp_raman_qsrc_comp->dtauc[lu][ivi] = tmp_out_int[ivi];

              /* The wavelength we are calculating is stored last in the ary */
              tmp_raman_qsrc_comp->dtauc[lu][output->crs.number_of_ramanshifts] = raman_qsrc_comp->dtauc[lu][iv];
            }

            /* Interpolate fbeam */
            for (ivi = 0; ivi < output->wl.nlambda_r; ivi++)
              tmp_in_int[ivi] = raman_qsrc_comp->fbeam[ivi];

            status = arb_wvn_double (output->wl.nlambda_r,
                                     tmp_wl,
                                     tmp_in_int,
                                     output->crs.number_of_ramanshifts,
                                     output->crs.wvl_of_ramanshifts,
                                     tmp_out_int,
                                     INTERP_METHOD_LINEAR,
                                     0);
            CHKERR (status);

            for (ivi = 0; ivi < output->crs.number_of_ramanshifts; ivi++)
              tmp_raman_qsrc_comp->fbeam[ivi] = tmp_out_int[ivi];

            /* The wavelength we are calculating is stored last in the ary */
            tmp_raman_qsrc_comp->fbeam[output->crs.number_of_ramanshifts] = raman_qsrc_comp->fbeam[iv];

            /* Interpolate uum */
            for (int lu = 0; lu < output->atm.nzout; lu++) {
              for (maz = 0; maz < input.rte.nstr; maz++) {
                for (int iq = 0; iq < input.rte.numu; iq++) {

                  for (ivi = 0; ivi < output->wl.nlambda_r; ivi++)
                    tmp_in_int[ivi] = raman_qsrc_comp->uum[lu][maz][iq][ivi];

                  status = arb_wvn_double (output->wl.nlambda_r,
                                           tmp_wl,
                                           tmp_in_int,
                                           output->crs.number_of_ramanshifts,
                                           output->crs.wvl_of_ramanshifts,
                                           tmp_out_int,
                                           INTERP_METHOD_LINEAR,
                                           0);
                  CHKERR (status);

                  for (ivi = 0; ivi < output->crs.number_of_ramanshifts; ivi++)
                    tmp_raman_qsrc_comp->uum[lu][maz][iq][ivi] = tmp_out_int[ivi];

                  /* The wavelength we are calculating is stored last in the ary */
                  tmp_raman_qsrc_comp->uum[lu][maz][iq][output->crs.number_of_ramanshifts] = raman_qsrc_comp->uum[lu][maz][iq][iv];
                }
              }
            }
            qwanted_wvl[0] = tmp_wl[iv]; /* Needed later in set_raman_source */

            status = optical_properties (input,
                                         output,
                                         output->wl.lambda_r[iv],
                                         ir,
                                         iv,
                                         iv,
                                         0,
                                         verbose,
                                         skip_optical_properties); /* in ancillary.c */
            CHKERR (status);

            /* We also need utau. */
            /* "very old" version, recycled */
            F77_FUNC (setout, SETOUT)
            (output->dtauc, &(output->atm.nlyr), &(output->atm.nzout), rte_in->utau, output->atm.zd, output->atm.zout_sur);

            /* Set the general inhomogeneous source applicable for Raman scattering   */
            if (iv == output->wl.raman_end_id)
              last = 1;
            status = set_raman_source (rte_in->qsrc,
                                       rte_in->qsrcu,
                                       input.rte.maxphi,
                                       output->atm.nlyr + 1,
                                       output->atm.nzout,
                                       input.rte.nstr,
                                       output->crs.number_of_ramanshifts,
                                       qwanted_wvl[0],
                                       output->crs.wvl_of_ramanshifts,
                                       rte_in->umu0,
                                       output->atm.zd,
                                       output->atm.zout_sur,
                                       output->atm.sza_r[iv],
                                       output->wl.fbeam[iv],
                                       input.r_earth,
                                       output->atm.microphys.dens[MOL_AIR][0][0],
                                       output->crs.crs_raman_RL,
                                       output->crs.crs_raman_RG,
                                       output->ssalb,
                                       input.rte.numu,
                                       input.rte.umu,
                                       rte_in->usrang,
                                       input.rte.cmuind,
                                       output->pmom,
                                       tmp_raman_qsrc_comp,
                                       output->atm.zout_comp_index,
                                       output->alt.altitude,
                                       last,
                                       verbose);
            CHKERR (status);

            free_raman_qsrc_components (tmp_raman_qsrc_comp,
                                        input.raman_fast,
                                        output->atm.nzout,
                                        input.rte.maxumu,
                                        input.rte.nphi,
                                        input.rte.nstr);

          } /* END if ( input.raman_fast )  */
          else {

            /* Find iv indices closest to ib wavelength above and below */
            qwanted_wvl[0] = (float) output->crs.wvl_of_ramanshifts[output->crs.number_of_ramanwavelengths - 1];
            iv1            = closest_above (qwanted_wvl[0], output->wl.lambda_r, output->wl.nlambda_r);
            iv2            = closest_below (qwanted_wvl[0], output->wl.lambda_r, output->wl.nlambda_r);
            /* Calculate optical properties for closest wavelength above and below wanted wavelength,*/
            /* and interpolate optical properties to wanted wavelength.       */
            status = optical_properties (input,
                                         output,
                                         qwanted_wvl[0],
                                         ir,
                                         iv1,
                                         iv2,
                                         0,
                                         verbose,
                                         skip_optical_properties); /* in ancillary.c */
            CHKERR (status);

            /* We also need utau at the "new" interpolated optical depth. */
            /* "very old" version, recycled */

            F77_FUNC (setout, SETOUT)
            (output->dtauc, &(output->atm.nlyr), &(output->atm.nzout), rte_in->utau, output->atm.zd, output->atm.zout_sur);

            /* Set the general inhomogeneous source applicable for Raman scattering   */

            status = set_raman_source (rte_in->qsrc,
                                       rte_in->qsrcu,
                                       input.rte.maxphi,
                                       output->atm.nlyr + 1,
                                       output->atm.nzout,
                                       input.rte.nstr,
                                       output->crs.number_of_ramanshifts,
                                       qwanted_wvl[0],
                                       output->crs.wvl_of_ramanshifts,
                                       rte_in->umu0,
                                       output->atm.zd,
                                       output->atm.zout_sur,
                                       output->atm.sza_r[iv],
                                       output->wl.fbeam[iv],
                                       input.r_earth,
                                       output->atm.microphys.dens[MOL_AIR][0][0],
                                       output->crs.crs_raman_RL,
                                       output->crs.crs_raman_RG,
                                       output->ssalb,
                                       input.rte.numu,
                                       input.rte.umu,
                                       rte_in->usrang,
                                       input.rte.cmuind,
                                       output->pmom,
                                       raman_qsrc_comp,
                                       output->atm.zout_comp_index,
                                       output->alt.altitude,
                                       last,
                                       verbose);
            CHKERR (status);
          }
        } /*  END    if ( input.raman_fast ) {} else    */
      }   /*  END    if ( input.raman )    */

      /* BRDF or Lambertian albedo */
      if (input.disort2_brdf != BRDF_NONE)
        rte_in->lamber = 0;
      else
        rte_in->lamber = 1;

      ds_in.nlyr  = output->atm.nlyr;
      ds_in.ntau  = output->atm.nzout;
      ds_in.nstr  = input.rte.nstr;
      ds_in.numu  = input.rte.numu;
      ds_in.nmom  = output->atm.nmom;
      ds_in.nphi  = input.rte.nphi;
      ds_in.accur = rte_in->accur;
      if (input.rte.disort_icm == DISORT_ICM_PHASE)
        ds_in.nphase = output->ntheta[0][0];

      /* choose how to do intensity correction */
      switch (input.rte.disort_icm) {
      case DISORT_ICM_OFF:
        ds_in.flag.intensity_correction = FALSE;
        break;
      case DISORT_ICM_MOMENTS:
        ds_in.flag.intensity_correction     = TRUE;
        ds_in.flag.old_intensity_correction = TRUE;
        break;
      case DISORT_ICM_PHASE:
        ds_in.flag.intensity_correction     = TRUE;
        ds_in.flag.old_intensity_correction = FALSE;
        break;
      default:
        fprintf (stderr, "Error: unknown disort_icm %d\n", input.rte.disort_icm);
        CHKERR (-1);
      }

      ds_in.flag.quiet  = input.quiet;
      ds_in.flag.ibcnd  = rte_in->ibcnd;
      ds_in.flag.planck = rte_in->planck;
      ds_in.flag.lamber = rte_in->lamber;
      ds_in.flag.usrtau = rte_in->usrtau;
      ds_in.flag.usrang = rte_in->usrang;
      ds_in.flag.onlyfl = rte_in->onlyfl;
      if (input.raman) {
        ds_in.flag.output_uum = TRUE;
        if (ir == 1) {
          ds_in.flag.general_source = TRUE;
          ds_in.bc.fluor            = 0.0;  // Only fluorescence source for the zeroth iteration
                                            // Otherwise source is included twice.
        } else {
          ds_in.flag.general_source = FALSE;
          ds_in.bc.fluor            = (double) output->flu.fluorescence_r[ib];
        }
        ds_in.bc.albedo = (double) output->alb.albedo_r[ib];
      } else {
        ds_in.flag.output_uum     = FALSE;
        ds_in.flag.general_source = FALSE;
        ds_in.bc.albedo           = (double) output->alb.albedo_r[iv];
        ds_in.bc.fluor            = (double) output->flu.fluorescence_r[iv];
      }
      ds_in.flag.spher = rte_in->spher;
      ds_in.radius     = c_r_earth;

      /* ds_in.flag.spher   = input.rte.pseudospherical; */

      ds_in.bc.btemp       = (double) rte_in->btemp;
      ds_in.bc.fbeam       = (double) rte_in->fbeam;
      ds_in.bc.fisot       = (double) rte_in->fisot;
      ds_in.bc.temis       = (double) rte_in->temis;
      ds_in.bc.ttemp       = (double) rte_in->ttemp;
      ds_in.bc.umu0        = (double) rte_in->umu0;
      ds_in.bc.phi0        = (double) output->atm.phi0_r[iv];
      ds_in.wvnmlo         = output->wl.wvnmlo_r[iv];
      ds_in.wvnmhi         = output->wl.wvnmhi_r[iv];
      ds_in.flag.prnt[0]   = rte_in->prndis2[0];
      ds_in.flag.prnt[1]   = rte_in->prndis2[1];
      ds_in.flag.prnt[2]   = rte_in->prndis2[2];
      ds_in.flag.prnt[3]   = rte_in->prndis2[3];
      ds_in.flag.prnt[4]   = rte_in->prndis2[4];
      ds_in.flag.brdf_type = input.disort2_brdf;

      c_disort_state_alloc (&ds_in);
      c_disort_out_alloc (&ds_in, &ds_out);

      for (int lc = 0; lc < output->atm.nlyr; lc++) {
        ds_in.dtauc[lc] = (double) output->dtauc[lc];
        ds_in.ssalb[lc] = (double) output->ssalb[lc];

        for (int k = 0; k <= output->atm.nmom; k++) {
          ds_in.pmom[k + lc * (ds_in.nmom_nstr + 1)] = (double) output->pmom[lc][0][k];
        }
      }
      for (int lc = 0; lc < output->atm.nzout; lc++)
        ds_in.utau[lc] = (double) rte_in->utau[lc];

      for (size_t iu = 0; iu < input.rte.numu; iu++)
        ds_in.umu[iu] = (double) input.rte.umu[iu];

      for (size_t iu = 0; iu < input.rte.nphi; iu++)
        ds_in.phi[iu] = (double) input.rte.phi[iu];

      for (int lu = 0; lu <= output->atm.nlyr; lu++) {
        if (rte_in->planck)
          ds_in.temper[lu] = (double) output->atm.microphys.temper[0][0][lu];
      }
      if (ds_in.flag.spher) {
        for (int lu = 0; lu <= output->atm.nlyr; lu++) {
          ds_in.zd[lu] = (double) output->atm.zd[lu];
        }
      }
      if (input.rte.disort_icm == DISORT_ICM_PHASE) {
        for (int imu = 0; imu < ds_in.nphase; imu++)
          ds_in.mu_phase[imu] = output->mu[0][0][imu];
        for (int imu = 0; imu < ds_in.nphase; imu++)
          for (int lc = 0; lc < ds_in.nlyr; lc++)
            ds_in.phase[imu + lc * (ds_in.nphase)] = (double) output->phase[lc][0][imu];
      }

      /* albedo stuff */
      switch (ds_in.flag.brdf_type) {
      case BRDF_RPV:
        ds_in.brdf.rpv->rho0  = output->rpv.rho0_r[iv];
        ds_in.brdf.rpv->k     = output->rpv.k_r[iv];
        ds_in.brdf.rpv->theta = output->rpv.theta_r[iv];
        ds_in.brdf.rpv->sigma = output->rpv.sigma_r[iv];
        ds_in.brdf.rpv->t1    = output->rpv.t1_r[iv];
        ds_in.brdf.rpv->t2    = output->rpv.t2_r[iv];
        ds_in.brdf.rpv->scale = output->rpv.scale_r[iv];
        break;
      case BRDF_CAM:
        ds_in.brdf.cam->u10  = input.cm.param[BRDF_CAM_U10];
        ds_in.brdf.cam->pcl  = input.cm.param[BRDF_CAM_PCL];
        ds_in.brdf.cam->xsal = input.cm.param[BRDF_CAM_SAL];
        break;
      case BRDF_HAPKE:
        ds_in.brdf.hapke->b0 = output->hapke.b0_r[iv];
        ds_in.brdf.hapke->h  = output->hapke.h_r[iv];
        ds_in.brdf.hapke->w  = output->hapke.w_r[iv];
        break;
      case BRDF_ROSSLI:
        ds_in.brdf.rossli->iso     = output->rossli.iso_r[iv];
        ds_in.brdf.rossli->vol     = output->rossli.vol_r[iv];
        ds_in.brdf.rossli->geo     = output->rossli.geo_r[iv];
        ds_in.brdf.rossli->hotspot = input.rossli.hotspot;
        break;
      default:
        break;
      }

      if (input.raman && ir == 1) {
        for (maz = 0; maz < ds_in.nstr; maz++)
          for (int lc = 0; lc < ds_in.nlyr; lc++)
            for (int iq = 0; iq < ds_in.nstr; iq++)
              ds_in.gensrc[iq + (lc + maz * ds_in.nlyr) * ds_in.nstr] = (double) rte_in->qsrc[maz][lc][iq];
        for (maz = 0; maz < ds_in.nstr; maz++)
          for (int lc = 0; lc < ds_in.nlyr; lc++)
            for (size_t iu = 0; iu < ds_in.numu; iu++)
              ds_in.gensrcu[iu + (lc + maz * ds_in.nlyr) * ds_in.numu] = (double) rte_in->qsrcu[maz][lc][iu];
      }

      if (ir == 1) {
        if (rte_in->qsrc != NULL)
          ASCII_free_double_3D (rte_in->qsrc, input.rte.nstr, output->atm.nzout);
        if (rte_in->qsrcu != NULL)
          ASCII_free_double_3D (rte_in->qsrcu, input.rte.nstr, output->atm.nzout);
      }

      if (((input.source == SRC_SOLAR) && (rte_in->umu0 > 0)) /* no need to call plane-parralel cdisort, as sun below horizon */
          || ((input.source == SRC_SOLAR) && (input.rte.pseudospherical)) /* pseudo-spherical cdisort  */
          || (input.source == SRC_THERMAL)) {
        c_disort (&ds_in, &ds_out);
      }

      if (rte_in->ibcnd) {
        for (size_t iu = 0; iu < input.rte.numu; iu++) {
          rte_out->albmed[iu] = (float) ds_out.albmed[iu];
          rte_out->trnmed[iu] = (float) ds_out.trnmed[iu];
        }
      }
      for (int lu = 0; lu < output->atm.nzout; lu++) {
        rte_out->dfdt[lu]   = (float) ds_out.rad[lu].dfdt;
        rte_out->rfldir[lu] = (float) ds_out.rad[lu].rfldir;
        rte_out->rfldn[lu]  = (float) ds_out.rad[lu].rfldn;
        rte_out->flup[lu]   = (float) ds_out.rad[lu].flup;
        rte_out->uavg[lu]   = (float) ds_out.rad[lu].uavg;
        rte_out->uavgdn[lu] = (float) ds_out.rad[lu].uavgdn;
        rte_out->uavgup[lu] = (float) ds_out.rad[lu].uavgup;
        rte_out->uavgso[lu] = (float) ds_out.rad[lu].uavgso;

        for (int j = 0; j < input.rte.nphi; j++)
          for (int iu = 0; iu < input.rte.numu; iu++)
            rte_out->uu[j][lu][iu] = ds_out.uu[iu + (lu + j * ds_in.ntau) * ds_in.numu];

        if (input.raman)
          for (int j = 0; j < input.rte.nstr; j++)  // j is the same as mazim in c_disort
            for (int iu = 0; iu < input.rte.numu; iu++)
              rte_out->uum[j][lu][iu] = ds_out.uum[iu + (lu + j * ds_in.ntau) * ds_in.numu];

        for (size_t iu = 0; iu < input.rte.numu; iu++)
          rte_out->u0u[lu][iu] = ds_out.u0u[iu + lu * ds_in.numu];
      }
      c_disort_out_free (&ds_in, &ds_out);
      c_disort_state_free (&ds_in);
    }

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    if (rte_solver != SOLVER_DISORT) {
      /* convert temporary Fortran arrays to permanent result for the fortran solver, not cdisort*/
      fortran2c_2D_float_ary_noalloc (output->atm.nzout, input.rte.maxumu, disort_u0u, rte_out->u0u);

      fortran2c_3D_float_ary_noalloc (input.rte.maxphi, output->atm.nzout, input.rte.maxumu, disort_uu, rte_out->uu);

      free (disort_pmom);
      free (disort2_pmom);
      free (disort2_phaso);
      free (disort2_mup);
      free (disort_u0u);
      free (disort_uu);
    }

    break;

  case SOLVER_TZS:

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    // BM, 11.2.2022: Is this check needed for the c_tzs solver?
    // It was needed for the Fortran version
    if (iv == output->wl.iv_lower && ib == 0) {
      status = F77_FUNC (dcheck, DCHECK) (&output->atm.nlyr,
                                          &output->atm.nzout,
                                          &input.rte.nstr,
                                          &input.rte.numu,
                                          &input.rte.nphi,
                                          &input.optimize_fortran,
                                          &input.optimize_delta);
      CHKERR (status);
    }

    /* call RTE solver */
    status = c_tzs (output->atm.nlyr,
                    output->dtauc,
                    output->atm.nlev_common,
                    output->atm.zd_common,
                    output->atm.nzout,
                    output->atm.zout_sur,
                    output->ssalb,
                    output->atm.microphys.temper[0][0],
                    output->wl.wvnmlo_r[iv],
                    output->wl.wvnmhi_r[iv],
                    rte_in->usrtau,
                    output->atm.nzout,
                    rte_in->utau,
                    rte_in->usrang,
                    input.rte.numu,
                    input.rte.umu,
                    input.rte.nphi,
                    input.rte.phi,
                    output->alb.albedo_r[iv],
                    rte_in->btemp,
                    rte_in->ttemp,
                    rte_in->temis,
                    rte_in->planck,
                    rte_in->prndis,
                    rte_in->header,
                    rte_out->rfldir,
                    rte_out->rfldn,
                    rte_out->flup,
                    rte_out->dfdt,
                    rte_out->uavg,
                    rte_out->uu,
                    input.quiet);
    CHKERR (status);

    if (planck_tempoff)
      rte_in->planck = 1;

    break;

  case SOLVER_SSS:

    sss_u0u = calloc (output->atm.nzout * input.rte.maxumu, sizeof (float));
    sss_uu  = calloc (output->atm.nzout * input.rte.maxumu * input.rte.maxphi, sizeof (float));

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    if (iv == output->wl.iv_lower && ib == 0) {
      status = F77_FUNC (dcheck, DCHECK) (&output->atm.nlyr,
                                          &output->atm.nzout,
                                          &input.rte.nstr,
                                          &input.rte.numu,
                                          &input.rte.nphi,
                                          &input.optimize_fortran,
                                          &input.optimize_delta);
      CHKERR (status);
    }

    /* call RTE solver */
    sss_pmom = c2fortran_3D_float_ary (output->atm.nlyr, 1, output->atm.nmom + 1, output->pmom);
    F77_FUNC (sss, SSS)
    (&output->atm.nlyr,
     output->dtauc,
     output->ssalb,
     &output->atm.nmom,
     sss_pmom,
     output->atm.microphys.temper[0][0],
     &output->wl.wvnmlo_r[iv],
     &output->wl.wvnmhi_r[iv],
     &rte_in->usrtau,
     &output->atm.nzout,
     rte_in->utau,
     &input.rte.nstr,
     &rte_in->usrang,
     &input.rte.numu,
     input.rte.umu,
     &input.rte.nphi,
     input.rte.phi,
     &rte_in->ibcnd,
     &rte_in->fbeam,
     &rte_in->umu0,
     &output->atm.phi0_r[iv],
     &rte_in->fisot,
     &rte_in->lamber,
     &output->alb.albedo_r[iv],
     &rte_in->btemp,
     &rte_in->ttemp,
     &rte_in->temis,
     &rte_in->planck,
     &rte_in->onlyfl,
     &rte_in->accur,
     rte_in->prndis2,
     rte_in->header,
     &output->atm.nlyr,
     &output->atm.nzout,
     &input.rte.maxumu,
     &input.rte.maxphi,
     &output->atm.nmom,
     rte_out->rfldir,
     rte_out->rfldn,
     rte_out->flup,
     rte_out->dfdt,
     rte_out->uavg,
     sss_uu,
     rte_out->albmed,
     rte_out->trnmed,
     rte_out->uavgdn,
     rte_out->uavgso,
     rte_out->uavgup,
     &(input.disort2_brdf),
     &(output->rpv.rho0_r[iv]),
     &(output->rpv.k_r[iv]),
     &(output->rpv.theta_r[iv]),
     &(input.cm.param[BRDF_CAM_U10]),
     &(input.cm.param[BRDF_CAM_PCL]),
     &(input.cm.param[BRDF_CAM_SAL]));

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    /* convert temporary Fortran arrays to permanent result */
    fortran2c_2D_float_ary_noalloc (output->atm.nzout, input.rte.maxumu, sss_u0u, rte_out->u0u);

    fortran2c_3D_float_ary_noalloc (input.rte.maxphi, output->atm.nzout, input.rte.maxumu, sss_uu, rte_out->uu);

    free (sss_pmom);
    free (sss_u0u);
    free (sss_uu);

    break;

  case SOLVER_SSSI:

#if HAVE_SSSI

    /* get cloud reflectivity from lookup table */
    status = read_isccp_reflectivity (output->sssi.type,
                                      output->atm.sza_r[iv],
                                      output->sssi.tautot,
                                      input.filename[FN_PATH],
                                      input.quiet,
                                      &(output->sssi.ref));
    CHKERR (status);

    if (input.verbose) {
      fprintf (stderr, "*** SSSI cloud properties:\n");
      fprintf (stderr, "    top level:          zd[%d] = %f\n", output->sssi.lctop, output->atm.zd[output->sssi.lctop]);
      switch (output->sssi.type) {
      case ISCCP_WATER:
        fprintf (stderr, "    type:               water\n");
        break;
      case ISCCP_ICE:
        fprintf (stderr, "    type:               ice\n");
        break;
      default:
        fprintf (stderr, "Error, unknown ISCCP cloud type %d\n", output->sssi.type);
        return -1;
      }
      fprintf (stderr, "    optical thickness:  %f\n", output->sssi.tautot);
      fprintf (stderr, "    reflectivity:       %f\n", output->sssi.ref);
    }
    fflush (stderr);

    sss_u0u = calloc (output->atm.nzout * input.rte.maxumu, sizeof (float));
    sss_uu  = calloc (output->atm.nzout * input.rte.maxumu * input.rte.maxphi, sizeof (float));

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    if (iv == output->wl.iv_lower && ib == 0) {
      status = F77_FUNC (dcheck, DCHECK) (&output->atm.nlyr,
                                          &output->atm.nzout,
                                          &input.rte.nstr,
                                          &input.rte.numu,
                                          &input.rte.nphi,
                                          &input.optimize_fortran,
                                          &input.optimize_delta);
      CHKERR (status);
    }

    /* call RTE solver */
    sss_pmom = c2fortran_3D_float_ary (output->atm.nlyr, 1, output->atm.nmom + 1, output->pmom);
    F77_FUNC (sssi, SSSI)
    (&output->atm.nlyr,
     output->dtauc,
     output->ssalb,
     &output->atm.nmom,
     sss_pmom,
     output->atm.microphys.temper[0][0],
     &output->wl.wvnmlo_r[iv],
     &output->wl.wvnmhi_r[iv],
     &rte_in->usrtau,
     &output->atm.nzout,
     rte_in->utau,
     &input.rte.nstr,
     &rte_in->usrang,
     &input.rte.numu,
     input.rte.umu,
     &input.rte.nphi,
     input.rte.phi,
     &rte_in->ibcnd,
     &rte_in->fbeam,
     &rte_in->umu0,
     &output->atm.phi0_r[iv],
     &rte_in->fisot,
     &rte_in->lamber,
     &output->alb.albedo_r[iv],
     &rte_in->btemp,
     &rte_in->ttemp,
     &rte_in->temis,
     &rte_in->planck,
     &rte_in->onlyfl,
     &rte_in->accur,
     rte_in->prndis2,
     rte_in->header,
     &output->atm.nlyr,
     &output->atm.nzout,
     &input.rte.maxumu,
     &input.rte.maxphi,
     &output->atm.nmom,
     rte_out->rfldir,
     rte_out->rfldn,
     rte_out->flup,
     rte_out->dfdt,
     rte_out->uavg,
     sss_uu,
     rte_out->albmed,
     rte_out->trnmed,
     rte_out->uavgdn,
     rte_out->uavgso,
     rte_out->uavgup,
     &(input.disort2_brdf),
     &(output->rpv.rho0_r[iv]),
     &(output->rpv.k_r[iv]),
     &(output->rpv.theta_r[iv]),
     &(input.cm.param[BRDF_CAM_U10]),
     &(input.cm.param[BRDF_CAM_PCL]),
     &(input.cm.param[BRDF_CAM_SAL]),
     &(output->sssi.ref),
     &(output->sssi.lctop));

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    /* convert temporary Fortran arrays to permanent result */
    fortran2c_2D_float_ary_noalloc (output->atm.nzout, input.rte.maxumu, sss_u0u, rte_out->u0u);

    fortran2c_3D_float_ary_noalloc (input.rte.maxphi, output->atm.nzout, input.rte.maxumu, sss_uu, rte_out->uu);

    free (sss_pmom);
    free (sss_u0u);
    free (sss_uu);
#else
    CHKERROUT (-1, "SSSI solver not available!");
#endif

    break;

  case SOLVER_SDISORT:

    if (iv == output->wl.iv_lower && ib == 0) {
      status = F77_FUNC (dcheck, DCHECK) (&output->atm.nlyr,
                                          &output->atm.nzout,
                                          &input.rte.nstr,
                                          &input.rte.numu,
                                          &input.rte.nphi,
                                          &input.optimize_fortran,
                                          &input.optimize_delta);
      CHKERR (status);
    }

    sdisort_beta = calloc (output->atm.nlyr + 1, sizeof (float));
    disort_pmom  = c2fortran_3D_float_ary (output->atm.nlyr, 1, input.rte.nstr + 1, output->pmom);
    disort_u0u   = calloc (output->atm.nzout * input.rte.maxumu, sizeof (float));
    disort_uu    = calloc (output->atm.nzout * input.rte.maxumu * input.rte.maxphi, sizeof (float));
    sdisort_sig  = calloc (output->atm.nlyr + 1, sizeof (float));
    if ((status = ASCII_calloc_float (&sdisort_denstab, output->atm.microphys.nsza_denstab, output->atm.nlyr + 1)) != 0)
      return status;
    if (output->atm.microphys.denstab_id > 0) {
      for (int lu = 0; lu <= output->atm.nlyr; lu++) {
        sdisort_sig[lu] = output->crs.crs_amf[iv][lu];
        for (int is = 0; is < output->atm.microphys.nsza_denstab; is++)
          sdisort_denstab[is][lu] = output->atm.microphys.denstab_amf[is][lu];
      }

      tosdisort_denstab = c2fortran_2D_float_ary (output->atm.microphys.nsza_denstab, output->atm.nlyr + 1, sdisort_denstab);
      if (sdisort_denstab != NULL)
        ASCII_free_float (sdisort_denstab, output->atm.microphys.nsza_denstab);
    }

    /* Definition of refind in (refractive index - 1), function SOLVER_SDISORT takes refractive index */
    for (int lu = 0; lu <= output->atm.nlyr; lu++)
      output->atm.microphys.refind[iv][lu] += 1.;

    F77_FUNC (sdisort, SDISORT)
    (&output->atm.nlyr,
     output->dtauc,
     output->ssalb,
     disort_pmom,
     output->atm.microphys.temper[0][0],
     &(output->wl.wvnmlo_r[iv]),
     &(output->wl.wvnmhi_r[iv]),
     &(rte_in->usrtau),
     &output->atm.nzout,
     rte_in->utau,
     &input.rte.nstr,
     &(rte_in->usrang),
     &input.rte.numu,
     input.rte.umu,
     &input.rte.nphi,
     input.rte.phi,
     &(rte_in->fbeam),
     sdisort_beta,
     &(rte_in->nil),
     &(rte_in->umu0),
     &output->atm.phi0_r[iv],
     &(rte_in->newgeo),
     output->atm.zd,
     &(rte_in->spher),
     &input.r_earth,
     &(rte_in->fisot),
     &output->alb.albedo_r[iv],
     &(rte_in->btemp),
     &(rte_in->ttemp),
     &(rte_in->temis),
     &input.rte.deltam,
     &(rte_in->planck),
     &(rte_in->onlyfl),
     &(rte_in->accur),
     &(rte_in->quiet),
     rte_in->ierror_s,
     rte_in->prndis,
     rte_in->header,
     &output->atm.nlyr,
     &output->atm.nzout,
     &input.rte.maxumu,
     &input.rte.nstr,
     &input.rte.maxphi,
     rte_out->rfldir,
     rte_out->rfldn,
     rte_out->flup,
     rte_out->dfdt,
     rte_out->uavg,
     disort_uu,
     disort_u0u,
     &input.rte.sdisort[SDISORT_NSCAT],
     rte_out->uavgdn,
     rte_out->uavgso,
     rte_out->uavgup,
     &input.rte.sdisort[SDISORT_NREFRAC],
     &input.rte.sdisort[SDISORT_ICHAPMAN],
     output->atm.microphys.refind[iv],
     &output->atm.microphys.nsza_denstab,
     output->atm.microphys.sza_denstab,
     sdisort_sig,
     tosdisort_denstab,
     output->dtauc_md);

    /* convert temporary Fortran arrays to permanent result */

    fortran2c_2D_float_ary_noalloc (output->atm.nzout, input.rte.maxumu, disort_u0u, rte_out->u0u);

    fortran2c_3D_float_ary_noalloc (input.rte.maxphi, output->atm.nzout, input.rte.maxumu, disort_uu, rte_out->uu);

    free (sdisort_beta);
    free (disort_pmom);
    free (disort_u0u);
    free (disort_uu);

    break;
  case SOLVER_SPSDISORT:

    if (iv == output->wl.iv_lower && ib == 0) {
      status = F77_FUNC (dcheck, DCHECK) (&output->atm.nlyr,
                                          &output->atm.nzout,
                                          &input.rte.nstr,
                                          &input.rte.numu,
                                          &input.rte.nphi,
                                          &input.optimize_fortran,
                                          &input.optimize_delta);
      CHKERR (status);
    }

    sdisort_beta = calloc (output->atm.nlyr + 1, sizeof (float));
    disort_pmom  = c2fortran_3D_float_ary (output->atm.nlyr, 1, input.rte.nstr + 1, output->pmom);
    disort_u0u   = calloc (output->atm.nzout * input.rte.maxumu, sizeof (float));
    disort_uu    = calloc (output->atm.nzout * input.rte.maxumu * input.rte.maxphi, sizeof (float));

    F77_FUNC (spsdisort, SPSDISORT)
    (&output->atm.nlyr,
     output->dtauc,
     output->ssalb,
     disort_pmom,
     output->atm.microphys.temper[0][0],
     &(output->wl.wvnmlo_r[iv]),
     &(output->wl.wvnmhi_r[iv]),
     &(rte_in->usrtau),
     &output->atm.nzout,
     rte_in->utau,
     &input.rte.nstr,
     &(rte_in->usrang),
     &input.rte.numu,
     input.rte.umu,
     &input.rte.nphi,
     input.rte.phi,
     &(rte_in->fbeam),
     sdisort_beta,
     &(rte_in->nil),
     &(rte_in->umu0),
     &output->atm.phi0_r[iv],
     &(rte_in->newgeo),
     output->atm.zd,
     &(rte_in->spher),
     &input.r_earth,
     &(rte_in->fisot),
     &output->alb.albedo_r[iv],
     &(rte_in->btemp),
     &(rte_in->ttemp),
     &(rte_in->temis),
     &input.rte.deltam,
     &(rte_in->planck),
     &(rte_in->onlyfl),
     &(rte_in->accur),
     &(rte_in->quiet),
     rte_in->ierror_s,
     rte_in->prndis,
     rte_in->header,
     &output->atm.nlyr,
     &output->atm.nzout,
     &input.rte.maxumu,
     &input.rte.nstr,
     &input.rte.maxphi,
     rte_out->rfldir,
     rte_out->rfldn,
     rte_out->flup,
     rte_out->dfdt,
     rte_out->uavg,
     disort_uu,
     disort_u0u,
     rte_out->uavgdn,
     rte_out->uavgso,
     rte_out->uavgup);

    /* convert temporary Fortran arrays to permanent result */

    fortran2c_2D_float_ary_noalloc (output->atm.nzout, input.rte.maxumu, disort_u0u, rte_out->u0u);

    fortran2c_3D_float_ary_noalloc (input.rte.maxphi, output->atm.nzout, input.rte.maxumu, disort_uu, rte_out->uu);

    free (sdisort_beta);
    free (disort_pmom);
    free (disort_u0u);
    free (disort_uu);

    break;
  case SOLVER_FTWOSTR:

    if (iv == output->wl.iv_lower && ib == 0) {
      status = F77_FUNC (tcheck, TCHECK) (&output->atm.nlyr,
                                          &output->atm.nzout,
                                          &input.rte.nstr,
                                          &input.rte.numu,
                                          &input.rte.nphi,
                                          &input.optimize_fortran,
                                          &input.optimize_delta);
      if (status != 0)
        return status;
    }

    twostr_gg = calloc (output->atm.nlyr, sizeof (float));
    for (int lu = 0; lu < output->atm.nlyr; lu++)
      twostr_gg[lu] = output->pmom[lu][0][1];

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    F77_FUNC (twostr, TWOSTR)
    (&output->alb.albedo_r[iv],
     &(rte_in->btemp),
     &input.rte.deltam,
     output->dtauc,
     &(rte_in->fbeam),
     &(rte_in->fisot),
     twostr_gg,
     rte_in->header,
     rte_in->ierror_t,
     &output->atm.nlyr,
     &output->atm.nzout,
     &(rte_in->newgeo),
     &output->atm.nlyr,
     &(rte_in->planck),
     &output->atm.nzout,
     rte_in->prntwo,
     &(rte_in->quiet),
     &input.r_earth,
     &(rte_in->spher),
     output->ssalb,
     &(rte_in->temis),
     output->atm.microphys.temper[0][0],
     &(rte_in->ttemp),
     &(rte_in->umu0),
     &(rte_in->usrtau),
     rte_in->utau,
     &(output->wl.wvnmlo_r[iv]),
     &(output->wl.wvnmhi_r[iv]),
     output->atm.zd,
     rte_out->dfdt,
     rte_out->flup,
     rte_out->rfldir,
     rte_out->rfldn,
     rte_out->uavg);
    free (twostr_gg);

    for (int lev = 0; lev < output->atm.nzout; lev++) {
      rte_out->uavgso[lev] = NAN;
      rte_out->uavgdn[lev] = NAN;
      rte_out->uavgup[lev] = NAN;
    }

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    break;

  case SOLVER_TWOSTR:

    twostr_ds.nlyr        = output->atm.nlyr;
    twostr_ds.ntau        = output->atm.nzout;
    twostr_ds.flag.planck = rte_in->planck;
    twostr_ds.flag.quiet  = rte_in->quiet;
    twostr_ds.flag.spher  = rte_in->spher;
    /* twostr_ds.flag.spher = input.rte.pseudospherical; */
    twostr_ds.flag.usrtau = rte_in->usrtau;
    c_twostr_state_alloc (&twostr_ds);
    c_twostr_out_alloc (&twostr_ds, &twostr_out);

    c_twostr_gg = calloc (output->atm.nlyr, sizeof (double));
    for (int lu = 0; lu < output->atm.nlyr; lu++) {
      c_twostr_gg[lu]     = (double) output->pmom[lu][0][1];
      twostr_ds.dtauc[lu] = (double) output->dtauc[lu];
      twostr_ds.ssalb[lu] = (double) output->ssalb[lu];
    }
    for (int lu = 0; lu < output->atm.nzout; lu++)
      twostr_ds.utau[lu] = (double) rte_in->utau[lu];

    for (int lu = 0; lu <= output->atm.nlyr; lu++) {
      twostr_ds.zd[lu] = (double) output->atm.zd[lu];
      if (rte_in->planck)
        twostr_ds.temper[lu] = (double) output->atm.microphys.temper[0][0][lu];
    }
    twostr_ds.bc.albedo    = (double) output->alb.albedo_r[iv];
    twostr_ds.bc.btemp     = (double) rte_in->btemp;
    twostr_ds.bc.fbeam     = (double) rte_in->fbeam;
    twostr_ds.bc.fisot     = (double) rte_in->fisot;
    twostr_ds.bc.temis     = (double) rte_in->temis;
    twostr_ds.bc.ttemp     = (double) rte_in->ttemp;
    twostr_ds.bc.umu0      = (double) rte_in->umu0;
    twostr_ds.flag.prnt[0] = rte_in->prntwo[0];
    twostr_ds.flag.prnt[1] = rte_in->prntwo[1];
    twostr_ds.wvnmlo       = output->wl.wvnmlo_r[iv];
    twostr_ds.wvnmhi       = output->wl.wvnmhi_r[iv];

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    c_twostr (&twostr_ds, &twostr_out, input.rte.deltam, c_twostr_gg, rte_in->ierror_t, c_r_earth);

    for (int lu = 0; lu < output->atm.nzout; lu++) {
      rte_out->dfdt[lu]   = (float) twostr_out.rad[lu].dfdt;
      rte_out->rfldir[lu] = (float) twostr_out.rad[lu].rfldir;
      rte_out->rfldn[lu]  = (float) twostr_out.rad[lu].rfldn;
      rte_out->flup[lu]   = (float) twostr_out.rad[lu].flup;
      rte_out->uavg[lu]   = (float) twostr_out.rad[lu].uavg;
    }

    free (c_twostr_gg);
    free (c_zd);
    c_twostr_state_free (&twostr_ds);
    c_twostr_out_free (&twostr_ds, &twostr_out);

    for (int lev = 0; lev < output->atm.nzout; lev++) {
      rte_out->uavgso[lev] = NAN;
      rte_out->uavgdn[lev] = NAN;
      rte_out->uavgup[lev] = NAN;
    }

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    break;

  case SOLVER_RODENTS: /* ulrike, Robert Buras' two-stream model */

    if (input.tipa == TIPA_DIR) /* BCA this should be somewhere else */
      /* for tipa_dir, delta-scaling is not yet implemented!!! */
      rodents_delta_method = RODENTS_DELTA_METHOD_OFF;
    else
      /* test showed that f=g*g is better than f=p2; master thesis to improve this! */
      rodents_delta_method = RODENTS_DELTA_METHOD_HG;

    twostr_gg = calloc (output->atm.nlyr, sizeof (float));
    for (int lu = 0; lu < output->atm.nlyr; lu++)
      twostr_gg[lu] = output->pmom[lu][0][1];
    twostr_ff = calloc (output->atm.nlyr, sizeof (float));
    if (rodents_delta_method ==
        RODENTS_DELTA_METHOD_ON) /* use second moment for delta-scaling */ /* BCA this should be somewhere else */
      for (int lu = 0; lu < output->atm.nlyr; lu++)
        twostr_ff[lu] = output->pmom[lu][0][2];
    else /* f is set to zero, and evtl set to g*g internally */
      for (int lu = 0; lu < output->atm.nlyr; lu++)
        twostr_ff[lu] = 0.0;

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */

    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    status = rodents (/* INPUT */
                      output->atm.nlyr,
                      output->dtauc,
                      output->ssalb,
                      twostr_gg,
                      twostr_ff,
                      rodents_delta_method,
                      output->atm.microphys.temper[0][0],
                      output->wl.wvnmlo_r[iv],
                      output->wl.wvnmhi_r[iv],
                      rte_in->usrtau,
                      output->atm.nzout,
                      rte_in->utau,
                      rte_in->fbeam,
                      rte_in->umu0,
                      output->alb.albedo_r[iv],
                      rte_in->btemp,
                      rte_in->planck,
                      /* NECESSARY FOR TIPA DIR */
                      input.tipa, /* if ==2, then tipa dir */
                      output->tausol,
                      /* OUTPUT */
                      rte_out->rfldn,  /* e_minus */
                      rte_out->flup,   /* e_plus */
                      rte_out->rfldir, /* s_direct */
                      rte_out->uavg);  /* KST ???? */

    if (status != 0) {
      fprintf (stderr, "Error %d returned by rodents()\n", status);
      return status;
    }

    free (twostr_gg);
    free (twostr_ff);

    for (int lev = 0; lev < output->atm.nzout; lev++) {
      rte_out->uavgso[lev] = NAN;
      rte_out->uavgdn[lev] = NAN;
      rte_out->uavgup[lev] = NAN;
    }

    if (planck_tempoff) {
      rte_in->planck = 1;
    }
    break; /* END of rodents */

  case SOLVER_TWOSTREBE:

    twostr_gg = calloc (output->atm.nlyr, sizeof (float));
    for (int lu = 0; lu < output->atm.nlyr; lu++)
      twostr_gg[lu] = output->pmom[lu][0][1];

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    status = twostrebe (output->dtauc,            /* dtau (rodents) = dtau_org (twostrebe) */
                        output->ssalb,            /* omega_0 */
                        twostr_gg,                /* g (rodents) = g_org (twostrebe) */
                        output->atm.nlyr + 1,     /* nlev */
                        rte_in->fbeam,            /* S_0 */
                        rte_in->umu0,             /* mu_0 */
                        output->alb.albedo_r[iv], /* surface albedo */
                        rte_in->planck,           /* whether to use planck */
                        input.rte.deltam,         /* delta scaling */
                        output->atm.nzout,        /* nzout */
                        output->atm.zd,           /* z-levels */
                        output->atm.microphys.temper[0][0],
                        rte_in->btemp, /* surface temperature */
                        output->wl.wvnmlo_r[iv],
                        output->wl.wvnmhi_r[iv],
                        output->atm.zout_sur, /* zout's (in km) */
                        /* output */
                        rte_out->rfldn,  /* e_minus */
                        rte_out->flup,   /* e_plus */
                        rte_out->rfldir, /* s_direct */
                        rte_out->uavg);  /* KST ???? */

    if (status != 0) {
      fprintf (stderr, "Error %d returned by twostrebe()\n", status);
      return status;
    }

    free (twostr_gg);

    for (int lev = 0; lev < output->atm.nzout; lev++) {
      rte_out->uavgso[lev] = NAN;
      rte_out->uavgdn[lev] = NAN;
      rte_out->uavgup[lev] = NAN;
    }

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    break; /* END twostrebe */

  case SOLVER_SCHWARZSCHILD: /*start schwarzschild 01.03.22*/

    /*twostr_gg = calloc (output->atm.nlyr, sizeof (float));
    for (lu = 0; lu < output->atm.nlyr; lu++)
      twostr_gg[lu] = output->pmom[lu][0][1];*/

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    status = schwarzschild (input.rte.schwarzschild_streams,
                            output->alb.albedo_r[iv], /* surface albedo */
                            output->dtauc,            /* dtau (rodents) = dtau_org (twostrebe) */
                            output->ssalb,            /* omega_0 */
                            output->atm.nlyr + 1,     /* nlev */
                            rte_in->planck,           /* whether to use planck */
                            output->atm.nzout,        /* nzout */
                            output->atm.zd,           /* z-levels */
                            output->atm.microphys.temper[0][0],
                            rte_in->btemp, /* surface temperature */
                            output->wl.wvnmlo_r[iv],
                            output->wl.wvnmhi_r[iv],
                            output->atm.zout_sur, /* zout's (in km) */
                            /* output */
                            rte_out->rfldn,  /* e_minus */
                            rte_out->flup,   /* e_plus */
                            rte_out->rfldir, /* s_direct */
                            rte_out->uavg);  /* KST ???? */

    if (status != 0) {
      fprintf (stderr, "Error %d returned by schwarzschild()\n", status);
      return status;
    }

    /*free (twostr_gg);*/

    for (int lev = 0; lev < output->atm.nzout; lev++) {
      rte_out->uavgso[lev] = NAN;
      rte_out->uavgdn[lev] = NAN;
      rte_out->uavgup[lev] = NAN;
    }

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    break; /* END schwarzschild */

  case SOLVER_TWOMAXRND: /* Bernhard Mayer, 7.7.2016, Nina Crnivec twostream with Maximum Random Overlap */

    twostr_gg     = calloc (output->atm.nlyr, sizeof (float));
    twostr_gg_clr = calloc (output->atm.nlyr, sizeof (float));
    for (int lu = 0; lu < output->atm.nlyr; lu++) {
      twostr_gg[lu]     = output->pmom[lu][0][1];
      twostr_gg_clr[lu] = output->pmom01_clr[lu];
    }

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    twostr_cf = calloc (output->atm.nlyr, sizeof (float));

    if (output->cf.nlev != 0) {

      if (output->atm.nlyr != output->cf.nlev) {
        fprintf (stderr,
                 "Fatal error! Cloud fraction grid different from atmospheric grid. %d levels vs. %d levels\n",
                 output->cf.nlev + 1,
                 output->atm.nlyr + 1);
        return -1;
      } else {
        for (int lu = 0; lu < output->atm.nlyr; lu++)
          twostr_cf[lu] = output->cf.cf[lu];
      }
    }

    status = twomaxrnd (output->dtauc,            /* dtau (rodents) = dtau_org (twostrebe) */
                        output->ssalb,            /* omega_0 */
                        twostr_gg,                /* g (rodents) = g_org (twostrebe) */
                        output->dtauc_clr,        /* dtau (rodents) = dtau_org (twostrebe) */
                        output->ssalb_clr,        /* omega_0 */
                        twostr_gg_clr,            /* g (rodents) = g_org (twostrebe) */
                        twostr_cf,                /* cloud fraction */
                        output->atm.nlyr + 1,     /* nlev */
                        rte_in->fbeam,            /* S_0 */
                        rte_in->umu0,             /* mu_0 */
                        output->alb.albedo_r[iv], /* surface albedo */
                        rte_in->planck,           /* whether to use planck */
                        input.rte.deltam,         /* delta scaling */
                        output->atm.nzout,        /* nzout */
                        output->atm.zd,           /* z-levels */
                        output->atm.microphys.temper[0][0],
                        rte_in->btemp, /* surface temperature */
                        output->wl.wvnmlo_r[iv],
                        output->wl.wvnmhi_r[iv],
                        output->atm.zout_sur, /* zout's (in km) */
                        /* output */
                        rte_out->rfldn,  /* e_minus */
                        rte_out->flup,   /* e_plus */
                        rte_out->rfldir, /* s_direct */
                        rte_out->uavg);  /* KST ???? */

    CHKERR (status);

    free (twostr_gg);
    free (twostr_gg_clr);
    free (twostr_cf);

    for (int lev = 0; lev < output->atm.nzout; lev++) {
      rte_out->uavgso[lev] = NAN;
      rte_out->uavgdn[lev] = NAN;
      rte_out->uavgup[lev] = NAN;
    }

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    break; /* END twomaxrnd */

  case SOLVER_DYNAMIC_TWOSTREAM: /* Bernhard Mayer, 23.7.2020, Richard Maier, dynamic twostream */

    twostr_gg     = calloc (output->atm.nlyr, sizeof (float));
    twostr_gg_clr = calloc (output->atm.nlyr, sizeof (float));
    for (int lu = 0; lu < output->atm.nlyr; lu++) {
      twostr_gg[lu]     = output->pmom[lu][0][1];
      twostr_gg_clr[lu] = output->pmom01_clr[lu];
    }

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    twostr_cf = calloc (output->atm.nlyr, sizeof (float));

    if (output->cf.nlev != 0) {

      if (output->atm.nlyr != output->cf.nlev) {
        fprintf (stderr,
                 "Fatal error! Cloud fraction grid different from atmospheric grid. %d levels vs. %d levels\n",
                 output->cf.nlev,
                 output->atm.nlyr + 1);
        return -1;
      } else {
        for (int lu = 0; lu < output->atm.nlyr; lu++)
          twostr_cf[lu] = output->cf.cf[lu];
      }
    }

    status = dynamic_twostream (input.rte.dynamic_iterations, /* number of iterations */
                                output->dtauc,                /* dtau (rodents) = dtau_org (twostrebe) */
                                output->ssalb,                /* omega_0 */
                                twostr_gg,                    /* g (rodents) = g_org (twostrebe) */
                                output->dtauc_clr,            /* dtau (rodents) = dtau_org (twostrebe) */
                                output->ssalb_clr,            /* omega_0 */
                                twostr_gg_clr,                /* g (rodents) = g_org (twostrebe) */
                                twostr_cf,                    /* cloud fraction */
                                output->atm.nlyr + 1,         /* nlev */
                                rte_in->fbeam,                /* S_0 */
                                rte_in->umu0,                 /* mu_0 */
                                output->alb.albedo_r[iv],     /* surface albedo */
                                rte_in->planck,               /* whether to use planck */
                                input.rte.deltam,             /* delta scaling */
                                output->atm.nzout,            /* nzout */
                                output->atm.zd,               /* z-levels */
                                output->atm.microphys.temper[0][0],
                                rte_in->btemp, /* surface temperature */
                                output->wl.wvnmlo_r[iv],
                                output->wl.wvnmhi_r[iv],
                                output->atm.zout_sur, /* zout's (in km) */
                                /* output */
                                rte_out->rfldn,  /* e_minus */
                                rte_out->flup,   /* e_plus */
                                rte_out->rfldir, /* s_direct */
                                rte_out->uavg);  /* KST ???? */

    if (status != 0) {
      fprintf (stderr, "Error %d returned by dynamic_twostream()\n", status);
      return status;
    }

    free (twostr_gg);
    free (twostr_gg_clr);
    free (twostr_cf);

    for (int lev = 0; lev < output->atm.nzout; lev++) {
      rte_out->uavgso[lev] = NAN;
      rte_out->uavgdn[lev] = NAN;
      rte_out->uavgup[lev] = NAN;
    }

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    break; /* END dynamic_twostream */

  case SOLVER_DYNAMIC_TENSTREAM: /* Bernhard Mayer, 6.12.2020, Richard Maier, dynamic tenstream */

    twostr_gg     = calloc (output->atm.nlyr, sizeof (float));
    twostr_gg_clr = calloc (output->atm.nlyr, sizeof (float));
    for (int lu = 0; lu < output->atm.nlyr; lu++) {
      twostr_gg[lu]     = output->pmom[lu][0][1];
      twostr_gg_clr[lu] = output->pmom01_clr[lu];
    }

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    status = dynamic_tenstream (input.rte.dynamic_iterations, /* number of iterations */
                                input.rte.dynamic_history,    /* read history files   */
                                output->caoth3d,
                                input.n_caoth + 2,
                                output->dtauc,            /* dtau (rodents) = dtau_org (twostrebe) */
                                output->ssalb,            /* omega_0 */
                                twostr_gg,                /* g (rodents) = g_org (twostrebe) */
                                output->atm.nlyr + 1,     /* nlev */
                                rte_in->fbeam,            /* S_0 */
                                rte_in->umu0,             /* mu_0 */
                                output->atm.phi0_r[iv],   /* phi_0 */
                                output->alb.albedo_r[iv], /* surface albedo */
                                rte_in->planck,           /* whether to use planck */
                                input.rte.deltam,         /* delta scaling */
                                output->atm.nzout,        /* nzout */
                                output->atm.zd,           /* z-levels */
                                output->atm.microphys.temper[0][0],
                                rte_in->btemp, /* surface temperature */
                                output->wl.wvnmlo_r[iv],
                                output->wl.wvnmhi_r[iv],
                                output->atm.zout_sur, /* zout's (in km) */
                                iv,
                                ib,
                                input.filename[FN_PATH], /* datapath */
                                input.quiet,
                                /* output */
                                rte_out->rfldn,  /* e_minus */
                                rte_out->flup,   /* e_plus */
                                rte_out->rfldir, /* s_direct */
                                rte_out->uavg,
                                rte_out->rfldir3d,
                                rte_out->rfldn3d,
                                rte_out->flup3d,
                                rte_out->abs3d);
    CHKERR (status);

    free (twostr_gg);
    free (twostr_gg_clr);

    for (int lev = 0; lev < output->atm.nzout; lev++) {
      rte_out->uavgso[lev] = NAN;
      rte_out->uavgdn[lev] = NAN;
      rte_out->uavgup[lev] = NAN;
    }

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    break; /* END dynamic_tenstream */

  case SOLVER_TWOMAXRND3C: /* Bernhard Mayer, 11.4.2018, Nina Crnivec twostream with Maximum Random Overlap and tripleclouds */
    twostr_gg_cldk = calloc (output->atm.nlyr, sizeof (float));
    twostr_gg_cldn = calloc (output->atm.nlyr, sizeof (float));
    twostr_gg_clr  = calloc (output->atm.nlyr, sizeof (float));
    for (int lu = 0; lu < output->atm.nlyr; lu++) {
      twostr_gg_cldk[lu] = output->pmom01_cldk[lu];
      twostr_gg_cldn[lu] = output->pmom01_cldn[lu];
      twostr_gg_clr[lu]  = output->pmom01_clr[lu];
    }

    /* ??? no need for thermal below 2um;                 ??? */
    /* ??? need to do that to avoid numerical underflow;  ??? */
    /* ??? however, this could be done without actually   ??? */
    /* ??? calling the solver                             ??? */
    if (rte_in->planck && output->wl.lambda_r[iv] < 2000.0) {
      rte_in->planck = 0;
      planck_tempoff = 1;
    }

    twostr_cf = calloc (output->atm.nlyr, sizeof (float));

    if (output->cf.nlev != 0) {

      if (output->atm.nlyr != output->cf.nlev) {
        fprintf (stderr,
                 "Fatal error! Cloud fraction grid different from atmospheric grid. %d levels vs. %d levels\n",
                 output->cf.nlev,
                 output->atm.nlyr + 1);
        return -1;
      } else {
        for (int lu = 0; lu < output->atm.nlyr; lu++)
          twostr_cf[lu] = output->cf.cf[lu];
      }
    }

#if HAVE_TWOMAXRND3C
    status = twomaxrnd3C (output->dtauc_cldk, /* dtau (rodents) = dtau_org (twostrebe) */
                          output->ssalb_cldk, /* omega_0 */
                          twostr_gg_cldk,     /* g (rodents) = g_org (twostrebe) */
                          output->dtauc_cldn, /* dtau (rodents) = dtau_org (twostrebe) */
                          output->ssalb_cldn, /* omega_0 */
                          twostr_gg_cldn,     /* g (rodents) = g_org (twostrebe) */
                          output->dtauc_clr,  /* dtau (rodents) = dtau_org (twostrebe) */
                          output->ssalb_clr,  /* omega_0 */
                          twostr_gg_clr,      /* g (rodents) = g_org (twostrebe) */
                          twostr_cf,          /* cloud fraction */
                          input.rte.twomaxrnd3C_scale_cf,
                          output->atm.nlyr + 1,     /* nlev */
                          rte_in->fbeam,            /* S_0 */
                          rte_in->umu0,             /* mu_0 */
                          output->alb.albedo_r[iv], /* surface albedo */
                          rte_in->planck,           /* whether to use planck */
                          input.rte.deltam,         /* delta scaling */
                          output->atm.nzout,        /* nzout */
                          output->atm.zd,           /* z-levels */
                          output->atm.microphys.temper[0][0],
                          rte_in->btemp, /* surface temperature */
                          output->wl.wvnmlo_r[iv],
                          output->wl.wvnmhi_r[iv],
                          output->atm.zout_sur, /* zout's (in km) */
                          /* output */
                          rte_out->rfldn,  /* e_minus */
                          rte_out->flup,   /* e_plus */
                          rte_out->rfldir, /* s_direct */
                          rte_out->uavg);  /* KST ???? */
    CHKERR (status);

#else
    CHKERROUT (-1, "twomaxrnd3C solver not included in uvspec build.");
#endif

    free (twostr_gg_cldk);
    free (twostr_gg_cldn);
    free (twostr_gg_clr);
    free (twostr_cf);

    for (int lev = 0; lev < output->atm.nzout; lev++) {
      rte_out->uavgso[lev] = NAN;
      rte_out->uavgdn[lev] = NAN;
      rte_out->uavgup[lev] = NAN;
    }

    if (planck_tempoff) {
      rte_in->planck = 1;
    }

    break; /* END twomaxrnd3C */

  case SOLVER_SOS:
#if HAVE_SOS
    ASCII_calloc_float (&pmom_sos, output->atm.nlyr, output->atm.nmom + 1);

    for (int lu = 0; lu < output->atm.nlyr; lu++)
      for (int k = 0; k < output->atm.nmom + 1; k++)
        pmom_sos[lu][k] = output->pmom[lu][0][k];

    status = sos (output->atm.nlyr,
                  rte_in->newgeo,
                  input.rte.nstr,
                  input.rte.sos_nscat,
                  output->alb.albedo_r[iv],
                  input.r_earth,
                  output->atm.zd,
                  output->ssalb,
                  pmom_sos,
                  output->dtauc,
                  output->atm.sza_r[iv],
                  output->atm.nzout,
                  input.rte.numu,
                  input.rte.umu,
                  rte_in->utau,
                  rte_out->rfldir,
                  rte_out->rfldn,
                  rte_out->flup,
                  rte_out->uavgso,
                  rte_out->uavgdn,
                  rte_out->uavgup,
                  rte_out->u0u);
    CHKERR (status);

    break;
#else
    CHKERROUT (-1, "sos solver not included in uvspec build. Please contact arve.kylling@gmail.com");
#endif

  case SOLVER_SSLIDAR:

    /* NOTE! Lidar only uses one umu!!! */

    phase_back = calloc ((size_t) output->atm.nlyr, sizeof (double*));

    for (int lu = 0; lu < output->atm.nlyr; lu++)
      phase_back[lu] = calloc ((size_t) output->nphamat, sizeof (double));

    /* this only works because mu[0] = -1.0 */
    for (int lu = 0; lu < output->atm.nlyr; lu++)
      for (int is = 0; is < output->nphamat; is++)
        phase_back[lu][is] = output->phase[lu][is][0];

    status = ss_lidar (/* input: atmosphere */
                       output->atm.nlyr,
                       output->atm.zd, /* z-levels */
                       output->alt.altitude,
                       output->dtauc,            /* optical depth */
                       output->ssalb,            /* omega_0 */
                       phase_back,               /* phase function in backward direction */
                       output->alb.albedo_r[iv], /* albedo (rodents) = Ag (twostrebe) */
                       /* input: lidar */
                       output->wl.lambda_r[iv], /* lidar wavelength */
                       input.sslidar[SSLIDAR_E0],
                       input.sslidar[SSLIDAR_POSITION],
                       input.rte.umu[0], /* only first umu is lidar direction */
                       input.sslidar_nranges,
                       input.sslidar[SSLIDAR_RANGE],
                       output->atm.zout_sur, /* ranges (in km) */
                       input.sslidar[SSLIDAR_EFF],
                       input.sslidar[SSLIDAR_AREA],
                       input.sslidar_polarisation,
                       /* OUTPUT / RESULT */
                       rte_out->sslidar_nphot,
                       rte_out->sslidar_nphot_q,
                       rte_out->sslidar_ratio);
    CHKERR (status);

    for (int lu = 0; lu < output->atm.nlyr; lu++)
      free (phase_back[lu]);
    free (phase_back);
    break;
  case SOLVER_NULL:
    /* do nothing */
    break;

  default:
    fprintf (stderr, "Error: RTE solver %d not yet implemented, call_solver (solve_rte.c)\n", input.rte.solver);
    CHKERR (-1);
  }

  if (verbose) {
    end = clock();
    fprintf (stderr, "*** last solver call: %f seconds\n", ((double) (end - start)) / CLOCKS_PER_SEC);
  }

  return 0;
}  // end of call_solver ()

static int init_rte_input (rte_input* rte, input_struct input, output_struct* output) {
  int status = 0;

  strcpy (rte->header, "");
  rte->accur = 1.0e-5;

  switch (input.source) {
  case SRC_NONE:
  case SRC_BLITZ:
  case SRC_LIDAR:
    rte->planck = 0;
    rte->fbeam  = 0.0;
    rte->fisot  = 0.0;
    break;

  case SRC_SOLAR:
    rte->planck = 0;

    if (input.rte.fisot > 0) {
      rte->fbeam = 0.0;
      rte->fisot = 1.0;
    } else {
      rte->fbeam = 1.0;
      rte->fisot = 0.0;
    }

    break;

  case SRC_THERMAL:
    rte->planck = 1;
    rte->fbeam  = 0.0;
    rte->fisot  = 0.0;

    rte->btemp = output->surface_temperature;
    rte->ttemp = output->atm.microphys.temper[0][0][0];

    rte->temis = 0.0;

    break;
  default:
    fprintf (stderr, "Error, unknown source %d\n", input.source);
    CHKERR (-1);
  }

  rte->umu0 = 0.0;

  rte->ierror_d[0] = 0;
  rte->ierror_s[0] = 0;
  rte->ierror_t[0] = 0;

  for (size_t i = 0; i < 7; i++)
    rte->prndis[i] = 0;

  for (size_t i = 0; i < 5; i++)
    rte->prndis2[i] = 0;

  for (size_t i = 0; i < 2; i++)
    rte->prntwo[i] = 0;

  rte->ibcnd  = input.rte.ibcnd;
  rte->lamber = 1;
  rte->newgeo = 1;
  rte->nil    = 0;
  rte->onlyfl = 1;
  rte->quiet  = input.quiet;
  rte->usrang = 0;
  rte->usrtau = 1;

  if (input.rte.pseudospherical || input.rte.solver == SOLVER_SDISORT ||
      input.rte.solver == SOLVER_SPSDISORT) /* these solvers are spherical by default */
    rte->spher = 1;
  else
    rte->spher = 0;

  rte->hl   = calloc (input.rte.maxumu + 1, sizeof (float));
  rte->utau = calloc (output->atm.nzout, sizeof (float));

  if (input.rte.numu > 0) {
    rte->onlyfl = 0;
    rte->usrang = 1;
  }

  switch (input.rte.solver) {
  case SOLVER_MONTECARLO:
    /* set number of photons */
    output->mc_photons = input.rte.mc.photons;
    break;
  case SOLVER_SDISORT:
  case SOLVER_SPSDISORT:
  case SOLVER_FDISORT1:
  case SOLVER_SSS:
  case SOLVER_SSSI:
  case SOLVER_TZS:
    for (int ip = 0; ip < input.rte.nprndis; ip++)
      (rte->prndis)[input.rte.prndis[ip] - 1] = 1;
    break;
  case SOLVER_DISORT:
  case SOLVER_FDISORT2:
    for (int ip = 0; ip < input.rte.nprndis; ip++)
      (rte->prndis2)[input.rte.prndis[ip] - 1] = 1;
    break;
  case SOLVER_TWOSTR:
  case SOLVER_FTWOSTR:
    for (int ip = 0; ip < input.rte.nprndis; ip++)
      if (input.rte.prndis[ip] == 1 || input.rte.prndis[ip] == 2)
        (rte->prntwo)[input.rte.prndis[ip] - 1] = 1;
    break;
  case SOLVER_SOS:
  case SOLVER_NULL:
  case SOLVER_RODENTS:
  case SOLVER_SCHWARZSCHILD:
  case SOLVER_TWOSTREBE:
  case SOLVER_TWOMAXRND:
  case SOLVER_TWOMAXRND3C:
  case SOLVER_DYNAMIC_TWOSTREAM:
  case SOLVER_DYNAMIC_TENSTREAM:
  case SOLVER_SSLIDAR:
    break;
  default:
    fprintf (stderr, "Error: RTE solver %d not yet implemented, init_rte_input (solve_rte.c)\n", input.rte.solver);
    return -1;
    break;
  }

  return status;
}  // end of init_rte_input ()

/***************************************************************/
/* calc_spectral_heating calculates the divergence of the flux */
/* either by differences of the flux or                        */
/* with the help of the actinic flux                           */
/***************************************************************/

static int calc_spectral_heating (input_struct   input,
                                  output_struct* output,
                                  float*         dz,
                                  double*        rho_mass_zout,
                                  float*         k_abs,
                                  float*         k_abs_layer,
                                  rte_output*    rte_out,
                                  float*         heat,
                                  float*         emis,
                                  float*         w_zout,
                                  int            iv) {
  int status = 0;
  int nzout = 0, nlev = 0, nlyr = 0;

  float* F_net = NULL;
  float  dFdz  = 0;
  float* c_p   = 0;

  float* dtheta_dz_layer = NULL;
  float* dtheta_dz       = NULL;

  float* dFdz_array = NULL;
  int    lz1 = 0, lz2 = 0, n_lz = 0;
  float  M_AIR           = MOL_MASS_AIR / 1000.0; /* molecular weight of air (kg mol-1) */
  float  planck_radiance = 0.0;

  float* z_center  = NULL;
  float* zout_in_m = NULL;
  int    outside   = 0;
  int    start     = 0;

  nlev  = output->atm.nlev;
  nlyr  = nlev - 1;
  nzout = output->atm.nzout;

  /* center heights of the atmosphere layers in m */
  z_center = calloc (nlyr, sizeof (float));
  CHKPOINTER ((z_center));
  for (int lc = 0; lc < nlyr; lc++)
    z_center[lc] = output->atm.zd[lc] * 1000.0 - dz[lc] / 2; /* 1000 == km -> m */

  /* zout_in_m == levels (layer boundaries) of zout levels in m above surface */
  zout_in_m = calloc (nzout, sizeof (float));
  CHKPOINTER ((zout_in_m));
  for (int lz = 0; lz < nzout; lz++)
    zout_in_m[lz] = output->atm.zout_sur[lz] * 1000.0; /* 1000 == km -> m */

  dtheta_dz = calloc (nzout, sizeof (float));
  CHKPOINTER ((dtheta_dz));

  /* calculate normalized (1/incident_flux) spectral heating rate */
  switch (input.heating) {
  case HEAT_LAYER_CD:
  case HEAT_LAYER_FD:

    F_net = calloc (nzout, sizeof (float));
    CHKPOINTER (F_net);

    dFdz_array = calloc (nzout, sizeof (float));
    CHKPOINTER ((dFdz_array));

    c_p = calloc (nzout, sizeof (float));
    CHKPOINTER ((c_p));

    switch (input.rte.solver) {
    case SOLVER_FDISORT1:
    case SOLVER_SDISORT:
    case SOLVER_FTWOSTR:
    case SOLVER_TWOSTR:
    case SOLVER_RODENTS:
    case SOLVER_SCHWARZSCHILD:
    case SOLVER_TWOSTREBE:
    case SOLVER_TWOMAXRND:
    case SOLVER_TWOMAXRND3C:
    case SOLVER_DYNAMIC_TWOSTREAM:
    case SOLVER_DYNAMIC_TENSTREAM:
    case SOLVER_SOS:
    case SOLVER_MONTECARLO:
    case SOLVER_FDISORT2:
    case SOLVER_DISORT:
    case SOLVER_SPSDISORT:
    case SOLVER_TZS:
    case SOLVER_SSS:
    case SOLVER_SSSI:
    case SOLVER_NULL:
      for (int lz = 0; lz < nzout; lz++)
        F_net[lz] = (rte_out->rfldir[lz] + rte_out->rfldn[lz]) - rte_out->flup[lz];
      break;
    default:
      fprintf (stderr, "Error: unknown solver id number %d\n", input.rte.solver);
      CHKERR (-1);
      break;
    }

    if (input.heating == HEAT_LAYER_CD)
      n_lz = nzout;
    if (input.heating == HEAT_LAYER_FD)
      n_lz = nzout - 1;

    for (int lz = 0; lz < n_lz; lz++) {

      if (input.heating == HEAT_LAYER_CD) {
        if (lz == 0) {
          lz1 = lz;
          lz2 = lz + 1;
        } else if (lz == output->atm.nzout - 1) {
          lz1 = lz - 1;
          lz2 = lz;
        } else {
          lz1 = lz - 1;
          lz2 = lz + 1;
        }

        if (lz != 0 && lz != output->atm.nzout - 1) {
          /* centered difference */ /* 1.0e+6: convert from cm-3 to m-3 */
          rho_mass_zout[lz] = output->atm.microphys.dens_zout[MOL_AIR][lz] * 1.0e+6 * M_AIR / AVOGADRO;

          /* mass weighted mean of specific heating rates */
          c_p[lz] = output->atm.microphys.c_p[lz];
        } else {
          /* boundary, no centered difference possible, (log) average density for forward difference */
          rho_mass_zout[lz] =
            log_average (output->atm.microphys.dens_zout[MOL_AIR][lz1], output->atm.microphys.dens_zout[MOL_AIR][lz2]) * 1.0e+6 *
            M_AIR / AVOGADRO;

          /* mass weighted mean of specific heating rates, assuming exponential change of density and linear change of c_p */
          c_p[lz] = mass_weighted_average (output->atm.microphys.c_p[lz1],
                                           output->atm.microphys.c_p[lz2],
                                           output->atm.microphys.dens_zout[MOL_AIR][lz1],
                                           output->atm.microphys.dens_zout[MOL_AIR][lz2]);
        }
      } else if (input.heating == HEAT_LAYER_FD) {
        /* forward difference */
        lz1 = lz;
        lz2 = lz + 1;

        /* effective density for one layer (logarithmic) */ /* 1.0e+6: convert from cm-3 to m-3 */
        rho_mass_zout[lz] =
          log_average (output->atm.microphys.dens_zout[MOL_AIR][lz1], output->atm.microphys.dens_zout[MOL_AIR][lz2]) * 1.0e+6 *
          M_AIR / AVOGADRO;

        /* mass weighted mean of specific heating rates, assuming exponential change of density and linear change of c_p */
        c_p[lz] = mass_weighted_average (output->atm.microphys.c_p[lz1],
                                         output->atm.microphys.c_p[lz2],
                                         output->atm.microphys.dens_zout[MOL_AIR][lz1],
                                         output->atm.microphys.dens_zout[MOL_AIR][lz2]);
      } else {
        fprintf (stderr, "Error, unknown processing scheme %d\n", input.processing);
        CHKERR (-1);
      }

      /* Compute the derivative dF/dz using a two-point formula */

      if (fabs ((F_net[lz1] - F_net[lz2]) / F_net[lz1]) > 1.0E-6 ||
          fabs ((output->atm.zout_sur[lz1] - output->atm.zout_sur[lz2]) * rho_mass_zout[lz]) > 1.0E-6) {
        dFdz = (F_net[lz1] - F_net[lz2]) / (output->atm.zout_sur[lz1] - output->atm.zout_sur[lz2]);
      } else
        dFdz = NAN;

      dFdz_array[lz] = 0.001 * dFdz; /* 1/1000 = km -> m, for dz */

      heat[lz] = dFdz_array[lz] / (rho_mass_zout[lz] * c_p[lz]);

      /* calculate dtheta_dz from zout-levels */
      dtheta_dz[lz] = (output->atm.microphys.theta_zout[lz1] - output->atm.microphys.theta_zout[lz2]) /
                      (1000.0 * (output->atm.zout_sur[lz1] - output->atm.zout_sur[lz2]));

      w_zout[lz] = (output->atm.microphys.theta_zout[lz] / output->atm.microphys.temper_zout[lz]) * 1.0 / dtheta_dz[lz] * heat[lz];
    }

    free (F_net);
    free (dFdz_array);
    free (c_p);

    break;
  case HEAT_LOCAL:

    /* spline interpolation of all level EXEPT LOWEST AND HIGHEST LEVEL */
    /* they are outside of the range of layer midpoints and must therefor be extrapolated !!!! */

    dtheta_dz_layer = calloc (nlyr, sizeof (float));
    CHKPOINTER ((dtheta_dz_layer));

    outside = 0;
    /* k_abs_layer == absorption coefficient representative for one layer (layer midpoint is z_center) */
    for (int lc = 0; lc < nlyr; lc++) {
      k_abs_layer[lc]     = (1.0 - output->ssalb[lc]) * output->dtauc[lc] / dz[lc];
      dtheta_dz_layer[lc] = (output->atm.microphys.theta[lc + 1] - output->atm.microphys.theta[lc]) /
                            ((output->atm.zd[lc + 1] - output->atm.zd[lc]) * 1000.0);
    }

    /* if uppermost zout level is above the uppermost layer midpoint (e.g. zout TOA),  */
    /* than extrapolate k_abs from the uppermost 2 layers                              */
    if (zout_in_m[nzout - 1] > z_center[0]) {
      outside = 1;
      /* exponentiell extrapolation, as linear might cause negative values */
      if (k_abs_layer[1] != 0.0)
        k_abs[nzout - 1] =
          k_abs_layer[0] * pow (k_abs_layer[0] / k_abs_layer[1], dz[0] / ((output->atm.zd[0] - output->atm.zd[2]) * 1000.0));
      else
        /* if not possible, take value from the last layer, in most cases also 0.0  */
        k_abs[nzout - 1] = k_abs_layer[0];
      /* first difference */
      dtheta_dz[nzout - 1] = dtheta_dz_layer[0]; /* zout is sorted ascending, z_atm descending */
    }

    /* if lowermost zout level is below the lowest layer midpoint (e.g. zout surface),  */
    /* than extrapolate k_abs from the lowermost 2 layers                              */
    if (zout_in_m[0] < z_center[nlyr - 1]) {
      outside = outside + 1;
      start   = 1;
      /* linear extrapolation */
      k_abs[0] = k_abs_layer[nlyr - 1] - dz[nlyr - 1] * (k_abs_layer[nlyr - 1] - k_abs_layer[nlyr - 2]) /
                                           ((output->atm.zd[nlyr] - output->atm.zd[nlyr - 2]) * 1000.0); /* 1000 = km -> m */
      /* canceled 2/2 in (dz/2)/((z[2]-z[0])/2) */
      /* last difference */
      dtheta_dz[0] = dtheta_dz_layer[nlyr - 1]; /* zout is sorted ascending, z_atm descending */
    }

    /* interpolate the rest, attention pointer arithmetic, last argument 1 means descending order of z_center */
    /* in versions before Jan 2008 also INTERP_METHOD_SPLINE was tested                 */
    /* but this caused overshootings, if clouds are present.                            */
    /* thereforwe use linear PLUS additional levels around the zout level, UH Feb 2008 */
    status = arb_wvn (nlyr, z_center, k_abs_layer, nzout - outside, zout_in_m + start, k_abs + start, INTERP_METHOD_LINEAR, 1);
    if (status != 0) {
      fprintf (stderr, " Error, interpolation of 'k_abs'\n");
      CHKERR (-1);
    }

    /* interpolate the rest, attention pointer arithmetic, last argument 1 means descending order of z_center */
    status =
      arb_wvn (nlyr, z_center, dtheta_dz_layer, nzout - outside, zout_in_m + start, dtheta_dz + start, INTERP_METHOD_LINEAR, 1);
    if (status != 0) {
      fprintf (stderr, " Error, interpolation of 'dtheta_dz'\n");
      CHKERR (-1)
    }

    for (int lz = 0; lz < nzout; lz++) {

      /* level (local) property */ /* 1.0e+6: convert from cm-3 to m-3 */
      rho_mass_zout[lz] = output->atm.microphys.dens_zout[MOL_AIR][lz] * 1.0e+6 * M_AIR / AVOGADRO;

      /* correction of the emission term of the heating rate, when calculated with actinic flux */
      if (input.source == SRC_THERMAL)
        planck_radiance = c_planck_func1 (output->wl.wvnmlo_r[iv], output->wl.wvnmhi_r[iv], output->atm.microphys.temper_zout[lz]);

      /* (else (in the solar case) planck_radiance == 0) */

      heat[lz] =
        k_abs[lz] * 4.0 * M_PI * (rte_out->uavg[lz] - planck_radiance) / (rho_mass_zout[lz] * output->atm.microphys.c_p[lz]);
      emis[lz] = k_abs[lz] * 4.0 * M_PI * (-planck_radiance) / (rho_mass_zout[lz] * output->atm.microphys.c_p[lz]);

      w_zout[lz] = (output->atm.microphys.theta_zout[lz] / output->atm.microphys.temper_zout[lz]) * 1.0 / dtheta_dz[lz] * heat[lz];
    }

    /* overwrite heating rate in case of dynamic twostream; in that case uavg contains dEnet */
    if (input.rte.solver == SOLVER_DYNAMIC_TWOSTREAM || input.rte.solver == SOLVER_DYNAMIC_TENSTREAM) {

      if (!input.quiet)
        fprintf (stderr, " ... overwriting heating rate by uavg (dynamic_twostream special)\n");

      for (int lz = 0; lz < nzout - 1; lz++) {
        /* BM, 18.9.2020: copied from layer_fd above */
        double rho_mass_tmp =
          log_average (output->atm.microphys.dens_zout[MOL_AIR][lz], output->atm.microphys.dens_zout[MOL_AIR][lz + 1]) * 1.0e+6 *
          M_AIR / AVOGADRO;

        /* mass weighted mean of specific heating rates,                   */
        /* assuming exponential change of density and linear change of c_p */
        double c_p_tmp = mass_weighted_average (output->atm.microphys.c_p[lz],
                                                output->atm.microphys.c_p[lz + 1],
                                                output->atm.microphys.dens_zout[MOL_AIR][lz],
                                                output->atm.microphys.dens_zout[MOL_AIR][lz + 1]);

        // factor 1000.0 converts z from km to m
        heat[lz] =
          rte_out->uavg[lz] / (rho_mass_tmp * c_p_tmp) / (output->atm.zout_sur[lz + 1] - output->atm.zout_sur[lz]) / 1000.0;
        emis[lz] = NAN;
      }
    }

    break;

  case HEAT_IPA3D: /* do nothing, BM, 20.1.2022 */
    break;

  default:
    fprintf (stderr, "Error, unknown processing scheme %d\n", input.processing);
    CHKERR (-1);
  }

  free (z_center);
  free (zout_in_m);
  free (dtheta_dz);

  return status;
}  // end of calc_spectral_heating ()

static float*** calloc_abs3d (const int Nx, const int Ny, const int Nz, const int* threed) {
  float*** tmp;

  tmp = calloc (Nz, sizeof (float**));
  VALIDPTR (tmp);

  for (int lc = 0; lc < Nz; lc++)
    if (threed[lc]) {
      tmp[lc] = calloc (Nx, sizeof (float*));
      VALIDPTR (tmp[lc]);

      for (int lx = 0; lx < Nx; lx++) {
        tmp[lc][lx] = calloc (Ny, sizeof (float));
        VALIDPTR (tmp[lc][lx]);
      }
    }

  return tmp;
}

static void free_abs3d (float*** abs3d, const int Nx, const int Ny, const int Nz, const int* threed) {

  for (int lc = 0; lc < Nz; lc++)
    if (threed[lc]) {
      for (int lx = 0; lx < Nx; lx++)
        free (abs3d[lc][lx]);
      free (abs3d[lc]);
    }
  free (abs3d);
}

/***********************************************************************************/
/* Function: generate_effective_cloud                                              */
/*                                                                                 */
/* Description:                                                                    */
/*  * Generates an effective cloud assuming random overlap (as in ECHAM)           */
/*    given cloud profiles and cloud fraction.                                     */
/*  * save the effective cloud optical property on the wc-structure                */
/*  * ic-structure is set to 0.0, as wc-structure contrains both now               */
/*                                                                                 */
/* Parameters (input/output):                                                      */
/*   input          uvspec input structure                                         */
/*   output         uvspec output structure                                        */
/*   save_cloud     unmodified cloud properties for given wavelength               */
/*   iv             wavelength index                                               */
/*   iq             subband index                                                  */
/*   verbose        flag for verbose output                                        */
/*                                                                                 */
/*                                                                                 */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    solver_rte.c                                                          */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    xxx 200x   C. Emde       Created                                             */
/*                                                                                 */
/***********************************************************************************/

static int
generate_effective_cloud (input_struct input, output_struct* output, save_optprop* save_cloud, int iv, int iq, int verbose) {

  float  tauw = 0.0, taui = 0.0, taua = 0.0, taum = 0.0, taur = 0.0, tau_clear = 0.0, tau_cloud = 0.0;
  float  g1d = 0.0, g2d = 0.0, fd = 0.0, g1i = 0.0, g2i = 0.0, fi = 0.0, gi = 0.0, gw = 0.0;
  float  ssai = 0.0, ssaw = 0.0;
  float *ssa = NULL, *tau = NULL, *g = NULL;

  /* Reflectivity of underlying layer, should be zero, because here
     we calculate just the effective optical thickness of one
     layer. This optical thickness is used in the RTE solver, where
     then of course the reflectivity of the neighbour layers is
     considered. */
  float pref = 0.0;
  /* Solar zenith angle */
  float mu = 0.0;
  /* Effective solar zenith angle, accounts for the decrease of the
     direct solar beam and the corresponding increase of the diffuse
     part of the radiation (in ECHAM). Here always the solar zenith
     angle is used (see below). */
  float mu_eff = 0.0;
  /* Diffusivity factor, 1.66 in ECHAM */
  float r = 1.66;
  /* Effective cloudyness */
  float C_eff = 0.0, product = 1.0;
  /* Output variables of swde */
  float pre1 = 0.0, pre2 = 0.0, ptr2 = 0.0;
  /* Transmission of clear and cloudy parts and of the layer */
  float transmission_clear = 0.0, transmission_cloud = 0.0, transmission_layer = 0.0;
  /* Variables for iteratiom */
  float taueff = 0.0;
  /* cloud fraction to scale the optical thickness */
  /* Attention: ECHAM input is already scaled to cloudy part */
  float cf            = 1.0;
  int   first_verbose = TRUE;

  /*   if ( input.cloud_overlap == CLOUD_OVERLAP_OFF ) { */
  /*     fprintf (stderr, "Error: call of %s, but cloud overlap schema is switched off\n", __func__ ); */
  /*     return -1; */
  /*   } */

  if (verbose && output->cf.nlev > 0) {
    fprintf (stderr, " ... generate effective cloud\n");
  }

  ssa = calloc (output->atm.nlev - 1, sizeof (float));
  tau = calloc (output->atm.nlev - 1, sizeof (float));
  g   = calloc (output->atm.nlev - 1, sizeof (float));

  mu = cos (output->atm.sza_r[iv] * M_PI / 180.0);

  /* scattering properties */
  for (int lc = 0; lc < output->atm.nlev - 1; lc++) {

    /* Ice and water clouds */
    tauw = save_cloud->tauw[lc];
    taui = save_cloud->taui[lc];

    g1d = save_cloud->g1d[lc];
    g2d = save_cloud->g2d[lc];
    fd  = save_cloud->fd[lc];
    gw  = g1d * fd + (1.0 - fd) * g2d;

    g1i = save_cloud->g1i[lc];
    g2i = save_cloud->g2i[lc];
    fi  = save_cloud->fi[lc];
    gi  = g1i * fi + (1.0 - fi) * g2i;

    ssaw = save_cloud->ssaw[lc];
    ssai = save_cloud->ssai[lc];

    /* molecular absorption */
    taum = output->atm.optprop.tau_molabs_r[0][0][lc][iv][iq];

    /* aerosol */
    taua = output->aer.optprop.dtau[iv][lc];

    /* Rayleigh */
    switch (input.ck_scheme) {
    case CK_FU:
      taur = output->atm.optprop.tau_rayleigh_r[0][0][lc][iv][iq];
      break;

    case CK_KATO:
    case CK_KATO2:
    case CK_KATO2_96:
    case CK_KATO2ANDWANDJI:
    case CK_AVHRR_KRATZ:
    case CK_FILE:
    case CK_LOWTRAN:
    case CK_CRS:
    case CK_REPTRAN:
    case CK_REPTRAN_CHANNEL:
    case CK_RAMAN:
      taur = output->atm.optprop.tau_rayleigh_r[0][0][lc][iv][0];
      break;

    default:
      fprintf (stderr, "Error: unsupported correlated-k scheme %d\n", input.ck_scheme);
      return -1;

      break;
    }

    /* Calculate mean optical properties of the layer */

    cf = 0.0;
    for (int isp = 0; isp < input.n_caoth; isp++)
      if (input.caoth[isp].source == CAOTH_FROM_ECHAM || input.caoth[isp].source == CAOTH_FROM_1D)
        cf = 1.0;

    if (cf == 0.0) { /* else */
      if (output->cf.cf[lc] != 0.0)
        cf = output->cf.cf[lc];
      else
        cf = 1.0;
    }

    /* total optical thickness */
    tau_cloud = (tauw + taui) / cf; /* Water and Ice cloud */
    tau_clear = taum + taua + taur; /* Molecular, Aerosol, and Rayleigh */

    /* Calculate effective Cloudiness */

    /* Total optical thickness */
    tau[lc] = tau_cloud + tau_clear;

    if ((tauw || taui) != 0.0) {
      g[lc] = (tauw * ssaw * gw + taui * ssai * gi) / (tauw * ssaw + taui * ssai);
      /* effective single scattering albedo */
      ssa[lc] = (tauw * ssaw + taui * ssai) / (tauw + taui);
    } else {
      g[lc]   = 1.;
      ssa[lc] = 1.;
    }

    if (input.source == SRC_THERMAL)
      mu_eff = 1. / r;
    else {
      /*Calculate effective zenith angle */
      for (int i = lc; i >= 0; i--) {
        if (mu != 0.0)
          product *= 1.0 - output->cf.cf[i] * (1.0 - exp (-((1.0 - ssa[i] * g[i] * g[i]) * (tau[lc])) / mu));
        else
          product *= 1.0 - output->cf.cf[i];
      }

      C_eff  = 1.0 - product;
      mu_eff = mu / (1.0 - C_eff + mu * r * C_eff);
    }

    /* Call ECHAM radiation routine SWDE. */
    if ((tauw || taui) != 0.0) {

      /*Initialize inputs for swde.*/
      pre1               = 0.0;
      pre2               = 0.0;
      transmission_cloud = 0.0;
      transmission_clear = 0.0;
      ptr2               = 0.0;

      if (verbose) {
        if (first_verbose) {
          fprintf (
            stderr,
            "   lc      g        pref    theta_eff   tau        tau_clear      ptr1      ptr2     ssa       cf     taueff\n");
          first_verbose = FALSE;
        }
        fprintf (stderr,
                 " %4d %9.6f %9.6f %9.3f %12.6e %12.6e %9.6f %9.6f %9.6f %8.5f",
                 lc,
                 g[lc],
                 pref,
                 acos (mu_eff) * 180.0 / M_PI,
                 tau[lc],
                 tau_clear,
                 transmission_cloud,
                 ptr2,
                 ssa[lc],
                 output->cf.cf[lc]);
      }

      if (ssa[lc] == 1.0)
        ssa[lc] = 0.99999;

      /* Perform radiative transfer (twostream) with scaled optical properties
         to calculate effective optical thickness. */
      /* Scattering + sbsorption + rayleigh + aerosol*/
      F77_FUNC (swde, SWDE) (&g[lc], &pref, &mu_eff, &tau[lc], &ssa[lc], &pre1, &pre2, &transmission_cloud, &ptr2);

      ptr2 = 0.0;
      /* Only absorption, rayleigh, aerosol */
      F77_FUNC (swde, SWDE) (&g[lc], &pref, &mu_eff, &tau_clear, &ssa[lc], &pre1, &pre2, &transmission_clear, &ptr2);

      /* Transmission of the layer*/
      transmission_layer = output->cf.cf[lc] * transmission_cloud + (1.0 - output->cf.cf[lc]) * transmission_clear;

      /* Calculate effective optical thickness */
      taueff = zbrent_taueff (mu_eff,
                              g[lc],
                              ssa[lc],
                              transmission_cloud,
                              transmission_layer,
                              tau_clear,
                              tau[lc],
                              0.00001); /* in solve_rte.c, next function */

      if (verbose)
        fprintf (stderr, "  %g\n", taueff);

      /* Set the optical properties to be used in rte calculation.*/

      /* wc is now representing both, water and ice clouds */
      output->caoth[input.i_wc].optprop.dtau[iv][lc] = taueff;
      output->caoth[input.i_wc].optprop.g1[iv][lc]   = g[lc];
      output->caoth[input.i_wc].optprop.g2[iv][lc]   = 0.0;
      output->caoth[input.i_wc].optprop.ff[iv][lc]   = 1.0;
      output->caoth[input.i_wc].optprop.ssa[iv][lc]  = ssa[lc];

      /* ic is now not needed any more */
      output->caoth[input.i_ic].optprop.dtau[iv][lc] = 0.0;
      output->caoth[input.i_ic].optprop.g1[iv][lc]   = 0.0;
      output->caoth[input.i_ic].optprop.g2[iv][lc]   = 0.0;
      output->caoth[input.i_ic].optprop.ff[iv][lc]   = 0.0;
      output->caoth[input.i_ic].optprop.ssa[iv][lc]  = 0.0;

    } else {
      /* wc is now representing both, water and ice clouds */
      output->caoth[input.i_wc].optprop.dtau[iv][lc] = 0.0;
      output->caoth[input.i_wc].optprop.g1[iv][lc]   = 0.0;
      output->caoth[input.i_wc].optprop.g2[iv][lc]   = 0.0;
      output->caoth[input.i_wc].optprop.ff[iv][lc]   = 0.0;
      output->caoth[input.i_wc].optprop.ssa[iv][lc]  = 0.0;

      /* ic is now not needed any more */
      output->caoth[input.i_ic].optprop.dtau[iv][lc] = 0.0;
      output->caoth[input.i_ic].optprop.g1[iv][lc]   = 0.0;
      output->caoth[input.i_ic].optprop.g2[iv][lc]   = 0.0;
      output->caoth[input.i_ic].optprop.ff[iv][lc]   = 0.0;
      output->caoth[input.i_ic].optprop.ssa[iv][lc]  = 0.0;
    }
  }

  free (ssa);
  free (tau);
  free (g);
  return 0;
}  // end of generate_effective_cloud ()

/** 
 * *zbrent_taueff* returns the effective optical depth of a layer.
 *
 * The function is a slightly modified version of the "zbrent" function in the 
 * "Numerical Recipes in C" (p. 352 ff.) for finding roots of an arbitrary function. 
 * The function is here the ECHAM radiative tarnsfer routine swde.f and the root of it 
 * is the effective optical thickness.
 * 
 * @param mu_eff effective zenith angle
 * @param g asymmetry parameter
 * @param ssa single scattering albedo
 * @param transmission_cloud transmission cloudy part
 * @param transmission_layer total transmission
 * @param x1 clear optical depth
 * @param x2 cloudy optical depth
 * @param tol accuracy
 * 
 * @return effective tau
 */
float zbrent_taueff (float mu_eff,
                     float g,
                     float ssa,
                     float transmission_cloud,
                     float transmission_layer,
                     float x1,
                     float x2,
                     float tol) {
  int   iter = 0;
  float a = x1, b = x2, c = x2, d = 0, e = 0, min1 = 0, min2 = 0;
  float fa = 0, fb = 0;
  float fc = 0, p = 0, q = 0, r = 0, s = 0, tol1 = 0, xm = 0;
  int   itmax = 100;
  float eps   = 3.0e-8;
  float dummy = 0.0;
  float pref  = 0.0;

  F77_FUNC (swde, SWDE) (&g, &pref, &mu_eff, &x1, &ssa, &dummy, &dummy, &transmission_cloud, &dummy);
  fa = transmission_cloud - transmission_layer;

  dummy = 0.0;
  pref  = 0.0;
  F77_FUNC (swde, SWDE) (&g, &pref, &mu_eff, &x2, &ssa, &dummy, &dummy, &transmission_cloud, &dummy);
  fb = transmission_cloud - transmission_layer;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    fprintf (stderr, "Root must be bracketed in zbrent \n");
    fprintf (stderr, "Please check whether the cloud effective optical thickness has been \n");
    fprintf (stderr, "calculated correctly. \n");
  }
  fc = fb;
  for (iter = 1; iter <= itmax; iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c  = a;
      fc = fa;
      e = d = b - a;
    }
    if (fabs (fc) < fabs (fb)) {
      a  = b;
      b  = c;
      c  = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0 * eps * fabs (b) + 0.5 * tol;
    xm   = 0.5 * (c - b);
    if (fabs (xm) <= tol1 || fb == 0.0)
      return b;
    if (fabs (e) >= tol1 && fabs (fa) > fabs (fb)) {
      s = fb / fa;
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0)
        q = -q;
      p    = fabs (p);
      min1 = 3.0 * xm * q - fabs (tol1 * q);
      min2 = fabs (e * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2)) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a  = b;
    fa = fb;
    if (fabs (d) > tol1)
      b += d;
    else
      b += SIGN (tol1, xm);

    dummy = 0.0;
    pref  = 0.0;
    F77_FUNC (swde, SWDE) (&g, &dummy, &mu_eff, &b, &ssa, &dummy, &dummy, &transmission_cloud, &dummy);
    fb = transmission_cloud - transmission_layer;
  }
  fprintf (stderr, "Maximum number of iterations exceeded in zbrent \n");
  return 0.0;
}  // end of zbrent_taueff ()

static int set_raman_source (double***              qsrc,
                             double***              qsrcu,
                             int                    maxphi,
                             int                    nlev,
                             int                    nzout,
                             int                    nstr,
                             int                    n_shifts,
                             float                  wanted_wl,
                             double*                wl_shifts,
                             float                  umu0,
                             float*                 zd,
                             float*                 zout,
                             float                  zenang,
                             float                  fbeam,
                             float                  radius,
                             float*                 dens,
                             double**               crs_RL,
                             double**               crs_RG,
                             float*                 ssalb,
                             int                    numu,
                             float*                 umu,
                             int                    usrang,
                             int*                   cmuind,
                             float***               pmom,
                             raman_qsrc_components* raman_qsrc_comp,
                             int*                   zout_comp_index,
                             float                  altitude,
                             int                    last,
                             int                    verbose) {

  /* Equation numbers in this function refer to equations in ESAS-LIGHT */
  /* report for WP2200.                                                 */

  int status = 0, twonm1, test = 0;
#if HAVE_SOS
  int* nfac = NULL;
#endif

  int static first = 1, nlyr = 0;  // ????????
  int            lua = 0, lub = 0, lv = 0, iu = 0, jq = 0, k = 0, l = 0, nn = 0, maz = 0;
  double         sum = 0, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
  double static *cmu = NULL, *cwt = NULL, ***ylmc = NULL, ***ylmu = NULL, ***ylm0 = NULL, sgn = 0, *tmpumu = NULL;
  double *       tmp = NULL, *tmpylmc = NULL, *tmpylmu = NULL, *tmpylm0 = NULL;
  double static *trs = NULL, **trs_shifted = NULL, *tmp_shifted = NULL, *tmp_dtauc = NULL, *chtau = NULL;
  double static* tmp_dtauc_shifted = NULL;
  double static *tmp_dens = NULL, *tmp_dens_org = NULL, **tmp_crs_RG = NULL, **tmp_crs_RL = NULL;
  double static *tmp_in = NULL, *tmp_out = NULL;

  double static *tmp_cumtau = NULL, *tmp_cumtauint = NULL, **tmp_user_dtauc = NULL, *tmp_zd = NULL, *tmp_zd_org = NULL;
  int static ntmp_zd = 0;

#if HAVE_SOS
  double static** fac = NULL;
#endif

  double g2_R   = 1 / 100.;      /* Legendre expansion coefficients for Raman scattering */
  double PI_R_b = 0, PI_R_d = 0; /* PIs in Eq. (28)-(29) */
  double ssalbRL = 0;            /* ssalb for Raman loss, Eq. (40) in Spurr et al 2008 */
  double ssalbRG = 0;            /* ssalb for Raman gain, Eq. (39) in Spurr et al 2008 */
  double deltaz = 0, km2cm = 1E+5;

  double delm0       = 1;
  double crs_RL_tot  = 0;
  double crs_RG_tot  = 0;
  double lambda_fact = 0; /* Conversion factor lambdap**2/lambda**2, Eq. A3 Edgington et al. 1999 */

  if (first) {
    nlyr    = nlev - 1;
    ntmp_zd = nlev;

    /* Need quadrature angles and weights and Legendre polynomials */

    cmu = calloc (nstr, sizeof (double));
    cwt = calloc (nstr, sizeof (double));
    nn  = nstr / 2;

    c_gaussian_quadrature (nn, cmu, cwt);

    /* Rearrange cmu and cwt such that they are ascending order. qsrc  should */
    /* have cmu in ascending order, however, internally qdisort does not treat*/
    /* cmu in ascending order. This feature is inherited from disort.         */
    for (int iq = 0; iq < nn; iq++) {
      cmu[nn + iq] = cmu[iq];
      cwt[nn + iq] = cwt[iq];
    }
    for (int iq = 0; iq < nn; iq++) {
      cmu[iq] = -cmu[nstr - 1 - iq];
      cwt[iq] = cwt[nstr - 1 - iq];
    }

    /* Calculate Legendre polynomials for each m */
    if ((tmpylmc = calloc ((size_t) ((nstr + 1) * nstr), sizeof (double))) == NULL)
      return ASCII_NO_MEMORY;
    if ((tmpylm0 = calloc ((size_t) ((nstr + 1) * nstr), sizeof (double))) == NULL)
      return ASCII_NO_MEMORY;
    if ((status = ASCII_calloc_double_3D (&ylmc, nstr, nstr, nstr + 1)) != 0)
      return ASCII_NO_MEMORY;
    if ((status = ASCII_calloc_double_3D (&ylm0, nstr, nstr, nstr + 1)) != 0)
      return ASCII_NO_MEMORY;
    if (usrang) {
      if ((tmpylmu = calloc ((size_t) ((nstr + 1) * numu), sizeof (double))) == NULL)
        return ASCII_NO_MEMORY;
      if ((status = ASCII_calloc_double_3D (&ylmu, nstr, numu, nstr + 1)) != 0)
        return ASCII_NO_MEMORY;
    }
    for (maz = 0; maz < nstr; maz++) {
      twonm1 = nstr - 1;
      nn     = 1;
      if ((tmp = calloc ((size_t) (1), sizeof (double))) == NULL)
        return ASCII_NO_MEMORY;
      tmp[0] = -umu0;
      c_legendre_poly (nn, maz, nstr, twonm1, tmp, tmpylm0);
      free (tmp);
      fortran2c_2D_double_ary_noalloc (nstr, nstr + 1, tmpylm0, ylm0[maz]);

      nn = nstr / 2;
      c_legendre_poly (nn, maz, nstr, twonm1, cmu, tmpylmc);
      fortran2c_2D_double_ary_noalloc (nstr, nstr + 1, tmpylmc, ylmc[maz]);

      /* Evaluate Legendre polynomials with negative -cmu- from those with*/
      /* positive -cmu-;  Dave Armstrong Eq. (15) */
      sgn = -1.0;
      for (l = 0; l < nstr; l++) {
        sgn = -sgn;
        for (int iq = nn; iq < nstr; iq++) {
          ylmc[maz][iq][l] = sgn * ylmc[maz][nstr - 1 - iq][l];
        }
      }

      if (usrang) {
        if ((tmpumu = calloc ((size_t) (numu), sizeof (double))) == NULL) {
          status = ASCII_NO_MEMORY;
          fprintf (stderr,
                   "Unable to allocate memory for tmpumu, status: %d returned at (line %d, function %s in %s)\n",
                   status,
                   __LINE__,
                   __func__,
                   __FILE__);
          return status;
        }
        for (iu = 0; iu < numu; iu++)
          tmpumu[iu] = umu[iu];
        c_legendre_poly (numu, maz, nstr, twonm1, tmpumu, tmpylmu);
        fortran2c_2D_double_ary_noalloc (numu, nstr + 1, tmpylmu, ylmu[maz]);
        free (tmpumu);
      }
    }
    free (tmpylmc);
    free (tmpylm0);
    if (usrang) {
      free (tmpylmu);
    }

    if ((trs = calloc ((size_t) (ntmp_zd), sizeof (double))) == NULL) {
      status = ASCII_NO_MEMORY;
      fprintf (stderr,
               "Unable to allocate memory for trs, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    if ((status = ASCII_calloc_double (&trs_shifted, ntmp_zd, n_shifts)) != 0) {
      fprintf (stderr,
               "Unable to allocate memory for trs_shifted, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    if ((chtau = calloc ((size_t) (ntmp_zd), sizeof (double))) == NULL) {
      status = ASCII_NO_MEMORY;
      fprintf (stderr,
               "Unable to allocate memory for chtau, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    if ((tmp_dtauc = calloc ((size_t) (ntmp_zd), sizeof (double))) == NULL) {
      status = ASCII_NO_MEMORY;
      fprintf (stderr,
               "Unable to allocate memory for tmp_dtauc, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    if ((tmp_dens = calloc ((size_t) (ntmp_zd), sizeof (double))) == NULL) {
      status = ASCII_NO_MEMORY;
      fprintf (stderr,
               "Unable to allocate memory for tmp_dens, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    if ((tmp_dens_org = calloc ((size_t) (nlev), sizeof (double))) == NULL) {
      status = ASCII_NO_MEMORY;
      fprintf (stderr,
               "Unable to allocate memory for tmp_dens_org, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    if ((tmp_dtauc_shifted = calloc ((size_t) (ntmp_zd), sizeof (double))) == NULL) {
      status = ASCII_NO_MEMORY;
      fprintf (stderr,
               "Unable to allocate memory for tmp_dtauc_shifted, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    if ((tmp_in = calloc ((size_t) (nzout), sizeof (double))) == NULL) {
      status = ASCII_NO_MEMORY;
      fprintf (stderr,
               "Unable to allocate memory for tmp_in, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    if ((tmp_out = calloc ((size_t) (ntmp_zd), sizeof (double))) == NULL) {
      status = ASCII_NO_MEMORY;
      fprintf (stderr,
               "Unable to allocate memory for tmp_out, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    if ((tmp_shifted = calloc ((size_t) (ntmp_zd), sizeof (double))) == NULL) {
      status = ASCII_NO_MEMORY;
      fprintf (stderr,
               "Unable to allocate memory for tmp_shifted, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    if ((status = ASCII_calloc_double (&tmp_crs_RG, ntmp_zd, n_shifts)) != 0) {
      fprintf (stderr, "Error %d allocating memory for crs_RG \n", status);
      fprintf (stderr, "       (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
      return status;
    }

    if ((status = ASCII_calloc_double (&tmp_crs_RL, ntmp_zd, n_shifts)) != 0) {
      fprintf (stderr, "Error %d allocating memory for crs_RG \n", status);
      fprintf (stderr, "       (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
      return status;
    }

    if ((status = ASCII_calloc_double (&tmp_user_dtauc, ntmp_zd, n_shifts + 1)) != 0) {
      fprintf (stderr,
               "Unable to allocate memory for user_dtauc, status: %d returned at (line %d, function %s in %s)\n",
               status,
               __LINE__,
               __func__,
               __FILE__);
      return status;
    }

    first = 0;

  } /* if ( first ) */

  /* Find all unique altitudes */
  ntmp_zd = nlev;
  if ((tmp_zd = calloc (ntmp_zd, sizeof (double))) == NULL) {
    fprintf (stderr, "Error, allocating memory for tmp_zd\n");
    fprintf (stderr, "      (line %d, function %s in %s)  \n", __LINE__, __func__, __FILE__);
    return -1;
  }
  if ((tmp_zd_org = calloc (nlyr + 1, sizeof (double))) == NULL) {
    fprintf (stderr, "Error, allocating memory for tmp_zd_org\n");
    fprintf (stderr, "      (line %d, function %s in %s)  \n", __LINE__, __func__, __FILE__);
    return -1;
  }
  for (int lc = 0; lc <= nlyr; lc++) {
    tmp_zd_org[lc]   = (double) zd[lc];
    tmp_dens_org[lc] = (double) dens[lc];
  }

  lv = 0;
  for (int lu = 0; lu < nlev - 1; lu++)
    if (zd[lu] >= 0.0)
      tmp_zd[lv++] = (double) zd[lu];

  if (zd[nlev] < 0.0)
    tmp_zd[lv] = 0.0;
  else
    tmp_zd[lv] = (double) zd[nlev];

    // aky20042012 removed this, made source zero at bottom level if altitude was set different from level.
    // Probably a leftover from when source varied within layers and forgotten to clean up.....
    //  if ( altitude != 0) {
    //aky    ntmp_zd = lv;  /* Number of output levels may be reduced due to altitude option */
    // }

#if HAVE_SOS
  /* Calculate geometric correction factor needed for chapman function */
  if ((nfac = calloc ((size_t) (ntmp_zd - 1), sizeof (int))) == NULL)
    return ASCII_NO_MEMORY;
  if ((status = ASCII_calloc_double (&fac, ntmp_zd - 1, 2 * (ntmp_zd - 1))) != 0)
    return status;

  /* Share dtauc from zd layering to zout */
  for (int lu = 0; lu < ntmp_zd - 1; lu++)
    tmp_dtauc[lu] = raman_qsrc_comp->dtauc[lu][n_shifts];

  chtau[0] = 0;
  for (int lc = 1; lc <= ntmp_zd - 1; lc++)
    chtau[lc] = c_chapman_simpler (lc, 0.5, ntmp_zd, tmp_zd, tmp_dtauc, zenang, radius);

  /* Transmittance of the atmosphere at wanted wavelength */
  trans_double (ntmp_zd - 1, chtau, trs);

  for (int is = 0; is < n_shifts; is++) {

    for (int lu = 0; lu < ntmp_zd - 1; lu++)
      tmp_dtauc_shifted[lu] = raman_qsrc_comp->dtauc[lu][is];

    chtau[0] = 0;
    for (int lc = 1; lc <= ntmp_zd - 1; lc++)
      chtau[lc] = c_chapman_simpler (lc, 0.5, ntmp_zd, tmp_zd, tmp_dtauc_shifted, zenang, radius);
    trans_double (ntmp_zd - 1, chtau, tmp_shifted);
    for (int lu = 0; lu < ntmp_zd; lu++)
      trs_shifted[lu][is] = tmp_shifted[lu];
  }

#else
  fprintf (stderr, "Error, need SOS source code for Raman scattering!\n");
  return -1;
#endif

  /* Interpolate density and Raman cross sections from zd to zout grid */
  status = arb_wvn_double (nlev, tmp_zd_org, tmp_dens_org, ntmp_zd, tmp_zd, tmp_dens, INTERP_METHOD_LOG, 1);

  for (int is = 0; is < n_shifts; is++) {
    for (int lu = 0; lu < nlev; lu++)
      tmp_in[lu] = crs_RG[lu][is];
    status = arb_wvn_double (nlev, tmp_zd_org, tmp_in, ntmp_zd, tmp_zd, tmp_out, INTERP_METHOD_LINEAR, 1);
    for (int lu = 0; lu < ntmp_zd; lu++)
      tmp_crs_RG[lu][is] = tmp_out[lu];
  }
  for (int is = 0; is < n_shifts; is++) {
    for (int lu = 0; lu < nlev; lu++)
      tmp_in[lu] = crs_RL[lu][is];
    status = arb_wvn_double (nlev, tmp_zd_org, tmp_in, ntmp_zd, tmp_zd, tmp_out, INTERP_METHOD_LINEAR, 1);
    for (int lu = 0; lu < ntmp_zd; lu++)
      tmp_crs_RL[lu][is] = tmp_out[lu];
  }

  test = 0;
  if (test) {
    /*************************************************************/
    /* To check that all angles etc. are correctly treated test  */
    /* with the direct beam source. This should give identical   */
    /* results for the diffuse radiation as the first raman      */
    /* wavelength loop                                           */
    /*************************************************************/

    delm0 = 1;
    fbeam = 1;
    for (maz = 0; maz < nstr; maz++) {
      if (maz > 0)
        delm0 = 0;
      for (int lu = 0; lu < nzout - 1; lu++) {
        for (int iq = 0; iq < nstr; iq++) {
          sum = 0;
          for (k = maz; k < nstr; k++) {
            sum += (2 * k + 1) * ssalb[lu] * pmom[lu][0][k] * ylmc[maz][iq][k] * ylm0[maz][0][k];
          }
          qsrc[maz][lu][iq] = sum * (2 - delm0) * fbeam / (4.0 * M_PI);
          if (lu == nzout) { /* Set source at bottom level */
            qsrc[maz][lu + 1][iq] = qsrc[maz][lu][iq] * trs[lu + 1];
          }
          qsrc[maz][lu][iq] = qsrc[maz][lu][iq] * trs[lu];
        }
        if (usrang) {
          for (iu = 0; iu < numu; iu++) {
            sum = 0;
            for (k = maz; k < nstr; k++)
              sum += (2 * k + 1) * ssalb[lu] * pmom[lu][0][k] * ylmu[maz][iu][k] * ylm0[maz][0][k];

            qsrcu[maz][lu][iu] = sum * (2 - delm0) * fbeam / (4.0 * M_PI);
            if (lu == nzout - 2) /* Set source at bottom level */
              qsrcu[maz][lu + 1][iu] = qsrcu[maz][lu][iu] * trs[lu + 1];
            qsrcu[maz][lu][iu] = qsrcu[maz][lu][iu] * trs[lu];
          }
        }
      }
    }
  } else { /* Raman scattering source */

    if (verbose)
      fprintf (stderr,
               "Raman_src %1s  %2s %2s    %2s     %3s        %4s          %4s          %4s          %4s          %4s  %9s %9s "
               "%9s %9s )\n",
               "m",
               "lu",
               "iq",
               "zd",
               "cmu",
               "sum1",
               "sum2",
               "sum3",
               "sum4",
               "qsrc",
               "ssalbRL",
               "ssalbRG",
               "PI_R_b",
               "PI_R_d");

    delm0 = 1;

    /* dtauc is not necessarily at all user altitudes. So first interpolate to all user altitudes.... */
    tmp_cumtau = calloc ((size_t) (ntmp_zd), sizeof (double));
    CHKPOINTER (tmp_cumtau);
    tmp_cumtauint = calloc ((size_t) (ntmp_zd), sizeof (double));
    CHKPOINTER (tmp_cumtauint);

    for (int is = 0; is <= n_shifts; is++) {

      /* Start by calculating the cumulative optical depth. */
      tmp_cumtau[0] = 0.0;
      for (int lu = 1; lu < ntmp_zd; lu++)
        tmp_cumtau[lu] = tmp_cumtau[lu - 1] + raman_qsrc_comp->dtauc[lu - 1][is];

      /* Interpolate the cumulative optical depth to all user altitudes. */
      status = arb_wvn_double (ntmp_zd, tmp_zd, tmp_cumtau, ntmp_zd, tmp_zd, tmp_cumtauint, INTERP_METHOD_LINEAR, 1);
      /* Finally calculate the optical depth of each layer. */
      tmp_user_dtauc[0][is] = 0.0;
      for (int lu = 1; lu < ntmp_zd; lu++)
        tmp_user_dtauc[lu][is] = tmp_cumtauint[lu] - tmp_cumtauint[lu - 1];
    }

    free (tmp_cumtau);
    free (tmp_cumtauint);

    for (maz = 0; maz < nstr; maz++) {
      if (maz > 0)
        delm0 = 0;

      for (int lu = 1; lu < ntmp_zd; lu++) { /* No need to include first level as source is calculated at middle of layer */

        lua = lu - 1;
        lub = lu;

        deltaz  = (tmp_zd[lua] - tmp_zd[lub]) * km2cm;
        ssalbRL = 0;
        //aky	tmpsum=0.0;
        for (int is = 0; is < n_shifts; is++) {
          //aky	  ssalbRL += 0.5 * deltaz * (tmp_dens[lua]*tmp_crs_RL[lua][is] + tmp_dens[lub]*tmp_crs_RL[lub][is]) /
          ssalbRL += deltaz * dlog_average (tmp_dens[lua] * tmp_crs_RL[lua][is], tmp_dens[lub] * tmp_crs_RL[lub][is]) /
                     tmp_user_dtauc[lub][n_shifts];
          //aky	  tmpsum += dlog_average(tmp_crs_RL[lua][is], tmp_crs_RL[lub][is]);
        }
        verbose = 0;
        if (verbose && maz == 0) {
          crs_RL_tot = 0;
          crs_RG_tot = 0;
          for (int is = 0; is < n_shifts; is++) {
            crs_RL_tot += tmp_crs_RL[lu][is];
            crs_RG_tot += tmp_crs_RG[lu][is];
          }
          fprintf (stderr, "%3d, zout: %7.2f, crs_RL_tot: %13.6e, crs_RG_tot: %13.6e\n", lu, zout[lu], crs_RL_tot, crs_RG_tot);
        }

        for (int iq = 0; iq < nstr; iq++) {
          sum1 = 0;
          sum2 = 0;

          /* Calculate PI_R_b (Eq. 28) for Raman scattering */
          if (maz == 0)
            PI_R_b = ylmc[maz][iq][0] * ylm0[maz][0][0] + 5 * g2_R * ylmc[maz][iq][2] * ylm0[maz][0][2];
          else if (maz == 1)
            PI_R_b = -5 * g2_R * ylmc[maz][iq][2] * ylm0[maz][0][2];
          else if (maz == 2)
            PI_R_b = +5 * g2_R * ylmc[maz][iq][2] * ylm0[maz][0][2];
          else
            PI_R_b = 0;

          /* First part of Raman scattering source term, Eq. (31). */

          ssalbRG = 0;
          for (int is = 0; is < n_shifts; is++) {
            lambda_fact = (wl_shifts[is] * wl_shifts[is]) / (wanted_wl * wanted_wl);
            sum         = 0;
            for (jq = 0; jq < nstr; jq++) {
              /* Calculate  PI_R_d (Eq. 29) for Raman scattering */
              if (maz == 0)
                PI_R_d = ylmc[maz][iq][0] * ylmc[maz][jq][0] + 5 * g2_R * ylmc[maz][iq][2] * ylmc[maz][jq][2];
              else if (maz == 1)
                PI_R_d = +5 * g2_R * ylmc[maz][iq][2] * ylmc[maz][jq][2];
              else if (maz == 2)
                PI_R_d = +5 * g2_R * ylmc[maz][iq][2] * ylmc[maz][jq][2];
              else
                PI_R_d = 0;

              sum += cwt[jq] * PI_R_d * raman_qsrc_comp->uum[zout_comp_index[lub]][maz][cmuind[jq]][is];
            }
            //aky	    ssalbRG = 0.5 * deltaz * lambda_fact * (tmp_dens[lua] * tmp_crs_RG[lua][is] +
            //aky				      tmp_dens[lub] * tmp_crs_RG[lub][is]) /
            ssalbRG = deltaz * lambda_fact *
                      dlog_average (tmp_dens[lua] * tmp_crs_RG[lua][is], tmp_dens[lub] * tmp_crs_RG[lub][is]) /
                      tmp_user_dtauc[lub][n_shifts];

            sum1 += ssalbRG * sum;
          }
          sum1 *= +1. / 2.;

          /* Second part of Raman scattering source term, Eq. (31). */

          for (int is = 0; is < n_shifts; is++) {
            lambda_fact = (wl_shifts[is] * wl_shifts[is]) / (wanted_wl * wanted_wl);
            //aky	    ssalbRG = 0.5 * deltaz * (tmp_dens[lua] * tmp_crs_RG[lua][is] + tmp_dens[lub] *
            //aky				      tmp_crs_RG[lub][is]) / tmp_user_dtauc[lub][n_shifts];
            ssalbRG = deltaz * dlog_average (tmp_dens[lua] * tmp_crs_RG[lua][is], tmp_dens[lub] * tmp_crs_RG[lub][is]) /
                      tmp_user_dtauc[lub][n_shifts];

            sum2 += raman_qsrc_comp->fbeam[is] * ssalbRG * trs_shifted[lub][is] * lambda_fact;
          }
          sum2 *= +((2 - delm0) / (4.0 * M_PI)) * PI_R_b;

          /* Third part of Raman scattering source term, Eq. (31). */

          sum = 0;
          for (jq = 0; jq < nstr; jq++) {
            /* Calculate  PI_R_d (Eq. 29) for Raman scattering */
            if (maz == 0)
              PI_R_d = ylmc[maz][iq][0] * ylmc[maz][jq][0] + 5 * g2_R * ylmc[maz][iq][2] * ylmc[maz][jq][2];
            else if (maz == 1)
              PI_R_d = +5 * g2_R * ylmc[maz][iq][2] * ylmc[maz][jq][2];
            else if (maz == 2)
              PI_R_d = +5 * g2_R * ylmc[maz][iq][2] * ylmc[maz][jq][2];
            else
              PI_R_d = 0;

            sum += cwt[jq] * PI_R_d * raman_qsrc_comp->uum[zout_comp_index[lub]][maz][cmuind[jq]][n_shifts];
          }

          sum3 = -(ssalbRL / 2) * sum;

          /* Fourth part of Raman scattering source term, Eq. (31). */

          sum4 = -((ssalbRL * fbeam) / (4.0 * M_PI)) * (2.0 - delm0) * PI_R_b * trs[lub];

          qsrc[maz][lu - 1][iq] = sum1 + sum2 + sum3 + sum4; /* lu starts at 1, source index is one less */

          //aky	  verbose = 0;
          if (verbose) {
            if (maz == 0 && iq == 0)
              fprintf (stderr,
                       "Raman_src %2d %2d %2d %7.3f %7.3f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e  %7.3f %7.3f %13.6e\n",
                       maz,
                       lu,
                       iq,
                       zout[lu],
                       cmu[iq],
                       sum1,
                       sum2,
                       sum3,
                       sum4,
                       qsrc[maz][lu - 1][iq],
                       ssalbRL,
                       ssalbRG,
                       PI_R_b,
                       PI_R_d,
                       trs[lu]);
          }

        } /* for (iq=0;iq<nstr;iq++) */

        if (usrang) {
          /* Also have to calculate the source for all user angles */
          for (iu = 0; iu < numu; iu++) {
            sum1 = 0;
            sum2 = 0;

            /* Calculate PI_R_b (Eq. 28) for Raman scattering */
            if (maz == 0)
              PI_R_b = ylmu[maz][iu][0] * ylm0[maz][0][0] + 5 * g2_R * ylmu[maz][iu][2] * ylm0[maz][0][2];
            else if (maz == 1)
              PI_R_b = -5 * g2_R * ylmu[maz][iu][2] * ylm0[maz][0][2];
            else if (maz == 2)
              PI_R_b = +5 * g2_R * ylmu[maz][iu][2] * ylm0[maz][0][2];
            else
              PI_R_b = 0;

            /* First part of Raman scattering source term, Eq. (31). */

            ssalbRG = 0;
            for (int is = 0; is < n_shifts; is++) {
              sum = 0;
              for (jq = 0; jq < nstr; jq++) {
                /* Calculate  PI_R_d (Eq. 29) for Raman scattering */
                if (maz == 0)
                  PI_R_d = ylmu[maz][iu][0] * ylmc[maz][jq][0] + 5 * g2_R * ylmu[maz][iu][2] * ylmc[maz][jq][2];
                else if (maz == 1)
                  PI_R_d = +5 * g2_R * ylmu[maz][iu][2] * ylmc[maz][jq][2];
                else if (maz == 2)
                  PI_R_d = +5 * g2_R * ylmu[maz][iu][2] * ylmc[maz][jq][2];
                else
                  PI_R_d = 0;

                sum += cwt[jq] * PI_R_d * raman_qsrc_comp->uum[zout_comp_index[lub]][maz][cmuind[jq]][is];
              }
              //aky	      ssalbRG = 0.5 * deltaz * (tmp_dens[lua] * tmp_crs_RG[lua][is] + tmp_dens[lub] *
              //aky					tmp_crs_RG[lub][is])/
              ssalbRG = deltaz * dlog_average (tmp_dens[lua] * tmp_crs_RG[lua][is], tmp_dens[lub] * tmp_crs_RG[lub][is]) /
                        tmp_user_dtauc[lub][n_shifts];

              sum1 += ssalbRG * sum;
            }
            sum1 *= +1. / 2.;

            /* Second part of Raman scattering source term, Eq. (31). */

            for (int is = 0; is < n_shifts; is++) {
              //aky	      ssalbRG = 0.5 * deltaz * (tmp_dens[lua] * tmp_crs_RG[lua][is] + tmp_dens[lub] *
              //aky					tmp_crs_RG[lub][is])/tmp_user_dtauc[lub][n_shifts];
              ssalbRG = deltaz * dlog_average (tmp_dens[lua] * tmp_crs_RG[lua][is], tmp_dens[lub] * tmp_crs_RG[lub][is]) /
                        tmp_user_dtauc[lub][n_shifts];

              sum2 += raman_qsrc_comp->fbeam[is] * ssalbRG * trs_shifted[lub][is];
            }
            sum2 *= +((2 - delm0) / (4.0 * M_PI)) * PI_R_b;

            /* Third part of Raman scattering source term, Eq. (31). */

            sum = 0;
            for (jq = 0; jq < nstr; jq++) {
              /* Calculate  PI_R_d (Eq. 29) for Raman scattering */
              if (maz == 0)
                PI_R_d = ylmu[maz][iu][0] * ylmc[maz][jq][0] + 5 * g2_R * ylmu[maz][iu][2] * ylmc[maz][jq][2];
              else if (maz == 1)
                PI_R_d = +5 * g2_R * ylmu[maz][iu][2] * ylmc[maz][jq][2];
              else if (maz == 2)
                PI_R_d = +5 * g2_R * ylmu[maz][iu][2] * ylmc[maz][jq][2];
              else
                PI_R_d = 0;
              sum += cwt[jq] * PI_R_d * raman_qsrc_comp->uum[zout_comp_index[lub]][maz][cmuind[jq]][n_shifts];
            }
            sum3 = -(ssalbRL / 2) * sum;

            /* Fourth part of Raman scattering source term, Eq. (31). */

            sum4 = -((ssalbRL * fbeam) / (4.0 * M_PI)) * (2 - delm0) * PI_R_b * trs[lub];

            qsrcu[maz][lu - 1][iu] = sum1 + sum2 + sum3 + sum4; /* lu starts at 1, source index is one less */

            verbose = 0;
            if (verbose) {
              if (maz == 0 && iu == 0)
                fprintf (stderr,
                         "Raman_srcu %2d %2d %2d %7.3f %7.3f %13.6e %13.6e %13.6e %13.6e %13.6e  %13.6e %13.6e  %7.3f %7.3f\n",
                         maz,
                         lu,
                         iu,
                         zout[lu],
                         umu[iu],
                         sum1,
                         sum2,
                         sum3,
                         sum4,
                         qsrcu[maz][lu][iu],
                         ssalbRL,
                         ssalbRG,
                         PI_R_b,
                         PI_R_d);
            }
          } /* for (iq=0;iq<nstr;iq++) */
        }   /* if ( usrang ) {  */
      }     /* for (lu=0;lu<nzout;lu++) */
    }
  }

  free (tmp_zd);
  free (tmp_zd_org);
  free (nfac);
  if (last) {
    free (cmu);
    free (cwt);
  }

  return status;
}  // end of set_raman_source()

/** Section to find the veiled core. Author: Paul Ockenfuß, Bernhard Mayer
 *    Modifies the parameters ext, g1, g2, ssa and ff in the profiles "wc" and "hiddencore_dummy"
 *    Important developer information: Units in caoth3d: ext[...] 1/m; atm.dxcld & atm.dycld in m, atm.zd_common in km!
 *    -z-Index="0" means surface in caoth3d.
 *    -z-Index="0" means TOA in "atm"!
 **/
static int find_veiled_core (const input_struct* input, const int iv_count, output_struct* output) {
  int isp_hiddencore = -1, isp_wc3D = -1;
  for (int isp = 0; isp < input->n_caoth; isp++) {
    if (strcmp (input->caoth[isp].name, "hiddencore_dummy") == 0)
      isp_hiddencore = isp;

    if (strcmp (input->caoth[isp].name, "wc") == 0)
      isp_wc3D = isp;
  }

  if (input->rte.mc.core_isactive && (iv_count == 0 || output->caoth3d[isp_wc3D].cldproperties == CLD_LWCREFF ||
                                      output->caoth3d[isp_wc3D].cldproperties == CLD_LWCREFFCF)) {
    assert (isp_hiddencore >= 0);
    assert (isp_wc3D >= 0);
    assert ((output->caoth3d[isp_wc3D]).nthreed > 0);
    assert ((output->caoth3d[isp_hiddencore]).nthreed > 0);
    if (!input->quiet)
      fprintf (stderr, "Found wc3D; starting to modify core\n");
    int      nx     = (output->caoth3d[isp_wc3D]).Nx;
    int      ny     = (output->caoth3d[isp_wc3D]).Ny;
    int      nz     = (output->caoth3d[isp_wc3D]).nlyr;
    int*     threed = (output->caoth3d[isp_wc3D]).threed;
    float    dx = output->atm.dxcld, dy = output->atm.dycld;
    float*** core;
    float*** distances;
    distances = calloc_float_3D (nz, nx, ny, "distances");
    CHKPOINTER (distances);

    if (input->rte.mc.core_inputfile) { /*If the user specified the cloud core*/
      int    core_nx, core_ny, core_nz, core_rows, core_flag = 0, status = 0;
      double core_dx, core_dy;
      float* core_dz;
      status = read_3D_caoth_header (input->rte.mc.core_inputfile,
                                     &core_nx,
                                     &core_ny,
                                     &core_nz,
                                     &core_flag,
                                     &core_dx,
                                     &core_dy,
                                     &core_dz,
                                     &core_rows);
      CHKERR (status);
      if (core_flag != 4) {
        fprintf (stderr,
                 "!!Warning: when reading %s, another flag than 4 was found! Is this really a file describing a cloud core?\n",
                 input->rte.mc.core_inputfile);
      }
      core_dx *= 1000;
      core_dy *= 1000;
      if (nx != core_nx || ny != core_ny || dx != core_dx || dy != core_dy) {
        fprintf (stderr,
                 "Error when reading core profile %s: The horizontal grid of the core must be the same as in wc_file 3D!\n",
                 input->rte.mc.core_inputfile);
        CHKERR (-1);
      }
      float*** column4;
      int*     core_threed = calloc (core_nz, sizeof (int));
      column4              = calloc_float_3D (core_nz, core_nx, core_ny, "column4");
      CHKPOINTER (column4);

      int* indz = calloc (core_rows - 2, sizeof (int));
      status    = read_3D_caoth_data (input->rte.mc.core_inputfile,
                                   core_nx,
                                   core_ny,
                                   core_nz,
                                   core_rows,
                                   &column4,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   indz);
      CHKERR (status);

      for (size_t i = 0; i < core_rows - 2; i++) {
        core_threed[indz[i]] = 1;
      }

      status = redistribute_3D (&column4, nx, ny, core_dz, core_nz, output->atm.zd_common, nz, core_threed, threed, 0);
      CHKERR (status);
      core = column4;

      free (indz);
      free (core_threed);
    } else { /*If the user specified no cloud core, start the corefinder*/
      core = calloc_float_3D (nz, nx, ny, "core");
      CHKPOINTER (core);
      // Initialize array with scattering coefficients from wc3D profile
      // This profile can contain layers which are not 3D. They are coming from interpolation
      // to other grids (e.g. ic_file). Therefore, wc scattering is always zero in these layers.
      arr3d_f* k_scat = calloc3d_float (nx, ny, nz);
      for (size_t i = 0; i < k_scat->nx; i++) {
        for (size_t j = 0; j < k_scat->ny; j++) {
          for (size_t k = 0; k < k_scat->nz; k++) {
            if ((output->caoth3d[isp_wc3D]).threed[k])
              set3d_f (k_scat, i, j, k, ((output->caoth3d[isp_wc3D]).ext[k][i][j]) * ((output->caoth3d[isp_wc3D]).ssa[k][i][j]));
            else
              set3d_f (k_scat, i, j, k, 0.0);
          }
        }
      }
      //calculate layer distances in atmosphere
      float* deltaz = malloc (k_scat->nz * sizeof (float));
      for (size_t i = 0; i < k_scat->nz; i++)
        deltaz[i] = 1000 * (output->atm.zd_common[k_scat->nz - i - 1] - output->atm.zd_common[k_scat->nz - i]);  //convert to meters

      // Calculate starting points for the corefinder: every specified zout-level will create one xy-layer of starting points
      // Notebook entry 21
      size_t  total_starts = output->atm.nzout_user * k_scat->nx * k_scat->ny;
      size_t* xstart       = calloc (total_starts, sizeof (size_t));
      size_t* ystart       = calloc (total_starts, sizeof (size_t));
      size_t* zstart       = calloc (total_starts, sizeof (size_t));
      size_t* z_index      = calloc (output->atm.nzout_user, sizeof (size_t));

      int kc       = 0;
      int iz_index = 0;
      for (kc = 0; kc < (output->caoth3d[isp_wc3D]).nlyr; kc++)
        if (output->mc.sample.sample[kc])
          z_index[iz_index++] = kc;

      // At TOA, the first layer below TOA is used for starting points
      if (output->mc.sample.sample[kc])
        z_index[iz_index++] = kc - 1;

      if (!input->quiet) {
        fprintf (stderr, "Corefinder starting points:\n");
        for (size_t i = 0; i < output->atm.nzout_user; i++)
          fprintf (stderr, "%d\n", (int) z_index[i]);
      }

      // consistency check
      if (iz_index != output->atm.nzout_user)
        CHKERROUT (-1, "number of output levels not matching in corefinder");

      for (size_t k = 0; k < output->atm.nzout_user; k++)
        for (size_t j = 0; j < k_scat->ny; j++)
          for (size_t i = 0; i < k_scat->nx; i++) {
            xstart[k * k_scat->ny * k_scat->nx + j * k_scat->nx + i] = i;
            ystart[k * k_scat->ny * k_scat->nx + j * k_scat->nx + i] = j;
            zstart[k * k_scat->ny * k_scat->nx + j * k_scat->nx + i] = z_index[k];
          }

      //Start to find the core
      arr3d_i* core_linarray = calloc3d_int (k_scat->nx, k_scat->ny, k_scat->nz);
      if (!input->quiet)
        fprintf (stderr, "...start finding core. Threshold: %.2f\n", input->rte.mc.core_threshold);
      arr3d_f* distances_linarray =
        get_distances2 (*k_scat, dx, dy, deltaz, xstart, ystart, zstart, total_starts, input->rte.mc.core_threshold);
      black_white_filter (distances_linarray, input->rte.mc.core_threshold, core_linarray);

      for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
          for (size_t k = 0; k < nz; k++) {
            core[k][i][j]      = (float) get3d_i (core_linarray, i, j, k);
            distances[k][i][j] = get3d_f (distances_linarray, i, j, k);
          }
        }
      }

      free (xstart);
      free (ystart);
      free (zstart);
      free (z_index);
      free (deltaz);
      free3d_float (k_scat);
      free3d_int (core_linarray);
      free3d_float (distances_linarray);
    } /*END Corefinder*/

    //Optionally: Save the core to the file given by input->rte.mc.core_savefile
    if (input->rte.mc.core_savefile) {
      if (!input->quiet)
        fprintf (stderr, "...saving core to %s\n", input->rte.mc.core_savefile);

      FILE* fp;
      if ((fp = fopen (input->rte.mc.core_savefile, "w")) == NULL) {
        fprintf (stderr, "Could not open %s!\n", input->rte.mc.core_savefile);
        CHKERR (-1);
      }
      fprintf (fp, "%d %d %d 4\n", nx, ny, nz);
      fprintf (fp, "%g %g", output->atm.dxcld / 1000.0, output->atm.dycld / 1000.0);

      for (size_t i = 0; i < nz + 1; i++)
        fprintf (fp, " %g", output->atm.zd_common[nz - i]);
      fprintf (fp, "\n#IndexX IndexY IndexZ Core Distance\n");
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++) {
            if (threed[k]) {
              if (core[k][i][j] == 0.0 && input->rte.mc.core_inputfile == NULL) {
                fprintf (fp, "%d\t%d\t%d\t%.1f\t%g\n", i + 1, j + 1, k + 1, core[k][i][j], distances[k][i][j]);
              } else {
                fprintf (fp, "%d\t%d\t%d\t%.1f\tnan\n", i + 1, j + 1, k + 1, core[k][i][j]);
              }
            }
            // else //do also print the non-threed (non-cloud) zeros
            // {
            //   fprintf(fp, "%d\t%d\t%d\t%.1f\n", i + 1, j + 1, k + 1, 0.0);
            // }
          }
        }
      }
      fclose (fp);
    }
    //Start to create a delta-scaled profile for the area inside the core
    if (!input->quiet)
      fprintf (stderr, "...start delta scaling core\n");

    int counter   = 0;
    int counter3d = 0;
    for (size_t k = 0; k < nz; k++) {
      if (threed[k]) {
        counter3d++;
        for (size_t i = 0; i < nx; i++) {
          for (size_t j = 0; j < ny; j++) {
            if (core[k][i][j] > 0.5) {
              counter++;
              float g         = (output->caoth3d[isp_wc3D]).g1[k][i][j];
              float w0        = (output->caoth3d[isp_wc3D]).ssa[k][i][j];
              float ext       = (output->caoth3d[isp_wc3D]).ext[k][i][j];
              float f_scaling = input->rte.mc.core_scale;
              if (f_scaling < 0.0)
                f_scaling = g;

              (output->caoth3d[isp_wc3D]).ext[k][i][j]       = 0.0;
              (output->caoth3d[isp_hiddencore]).ext[k][i][j] = (1 - w0 * f_scaling) * ext;
              (output->caoth3d[isp_hiddencore]).ssa[k][i][j] = (1 - f_scaling) * w0 / (1 - w0 * f_scaling);
              (output->caoth3d[isp_hiddencore]).g1[k][i][j]  = (g - f_scaling) / (1 - f_scaling);
              (output->caoth3d[isp_hiddencore]).g2[k][i][j]  = 0.0;
              (output->caoth3d[isp_hiddencore]).ff[k][i][j]  = 0.0;
            } else {
              (output->caoth3d[isp_hiddencore]).ext[k][i][j] = 0.0;
            }
          }
        }
      }
    }

    if (!input->quiet)
      fprintf (stderr,
               "scaled %d pixels of %d 3d pixels! (%.2f%%)\n",
               counter,
               counter3d * nx * ny,
               100 * counter / ((float) counter3d * nx * ny));

    free_float_3D (core);
    free_float_3D (distances);
  } /* End Section to find hidden core */

  return 0;
}  // end of find_veiled_core()
