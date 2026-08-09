#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "bpdf.h"
#include "errors.h"
#include "mystic.h"

#define debugd(...)
#define debug(...)

/* Some of this code was based on FORTRAN code by Pavel Litvinov
 * Most formulas are explained in much greater detail in the following Bachelor thesis:
 * Relevance of Polarized Reflection at Land Surfaces
 * for the Interpretation of Earthshine Observations (Cedric Wehrum, 2023)
 */

/* Multiply two matrices (lhs and rhs) and store it in result
 * if result is NULL it will be automatically allocated
 * TODO simplify to only support 4x4 matrices and remove auto allocation
 */
static inline double* matMul (enum MatrixType matrixType, double* lhs, double* rhs, double* result) {
  if (result == NULL)
    result = malloc (sizeof (double) * matrixType * matrixType);
  MatIndexFunc ind = matIndexByType[matrixType];
  for (unsigned i = 0; i < matrixType; i++) {
    for (unsigned j = 0; j < matrixType; j++) {
      result[ind (i, j)] = 0.;
      for (unsigned k = 0; k < matrixType; k++)
        result[ind (i, j)] += lhs[ind (i, k)] * rhs[ind (k, j)];
    }
  }
  return result;
}

/* An older model of Litvinov
 * TODO either test or remove model
 */
int reflection_matrix_litvinov_sron (enum MatrixType matrixType,
                                     double          solar_za,
                                     double          view_za,
                                     double          solar_aa,
                                     double          view_aa,
                                     const double    params[3],
                                     double _Complex m_ind,
                                     double* rMat) {
  double alpha;
  double mu_0 = solar_za, mu_v = view_za;
  double theta_0 = acos (solar_za), phi_0 = solar_aa, theta_v = acos (view_za), phi_v = view_aa;
  double n_0[3]    = {sin (theta_0) * cos (phi_0), sin (theta_0) * sin (phi_0), cos (theta_0)};
  double n_v[3]    = {sin (fabs (theta_v)) * cos (phi_v), sin (fabs (theta_v)) * sin (phi_v), cos (theta_v)};
  double n_abs_sum = 0.;
  for (unsigned i = 0; i < 3; i++)
    n_abs_sum += (n_0[i] + n_v[i]) * (n_0[i] + n_v[i]);
  double mu_n = n_v[2] * n_v[2] + n_0[2] * n_0[2] / sqrt (n_abs_sum);

  double mu_scat = mu_0 * mu_v - sqrt (1. - mu_0 * mu_0) * sqrt (1. - mu_v * mu_v) * cos (phi_0 - phi_v);
  double theta_scat;
  if (mu_scat > 1.) {
    mu_scat    = 1.;
    theta_scat = 0.;
  } else if (mu_scat < -1.) {
    mu_scat    = -1.;
    theta_scat = M_PI;
  } else
    theta_scat = acos (mu_scat);
  double shadow_param = params[2] * theta_scat;
  double f_shadow     = ((1. + cos (shadow_param)) / 2.);
  f_shadow *= f_shadow * f_shadow;
  if (fabs (shadow_param) > M_PI)
    f_shadow = 0.;

  double mu_n_2 = mu_n * mu_n;
  double f_n    = 1. / (2. * mu_n_2 * mu_n * params[1]) * exp (-(1. - mu_n_2) / (2. * mu_n_2 * params[1]));

  reflFresnelMy (matrixType, solar_za, view_za, fabs (phi_0 - phi_v), m_ind, rMat, &alpha);
  double coefficient = params[0] / (4. * mu_n * (mu_0 + mu_v)) * f_shadow * f_n;

  for (unsigned i = 0; i < matrixType * matrixType; i++)
    rMat[i] = rMat[i] * coefficient;
  return 0;
}

/* N_BETA and N_ALPHA are the number of points used in the Gauss-Legendre quadrature
 * alpha refers to the azimuth angle and beta refers to the zenith angle of the integral
 */
#define N_BETA 9
#define N_ALPHA 22

// See header file
int reflection_matrix_litvinov (enum MatrixType      matrixType,
                                double               solar_za,
                                double               view_za,
                                double               solar_aa,
                                double               view_aa,
                                const unsigned short numParams,
                                const double         params[5],
                                double _Complex m_ind,
                                double* rMat) {
  // define some numerical constants
  double eps = 0.00000000001, sin_eps = 0.00000001, sig_min = 0.02;
  // define some constant model parameters
  // Note that sig2 can be overwritten by params[4]
  double par_free = 4. / 3., alb_m = 0.2;  //alb_m
  //double depol = 0.005, sig2 = 0.005; //for ocean
  double depol = 0.005, sig2 = 0.9f;  //for land

  // define all intermediate variables
  double alb_s, sig1, fren_fr, k_shad, mu0, muv, phii, del_phi, sin_ti, sin_tv, cos_tet, tet, shadow, r11_elem, pRel, sin_tet,
    depol_22, mu_n1, mubv, mubi, sinbv, sinbi, mu1v, mu2i, randF0;
  double randF1[N_ALPHA][N_BETA], randF1_fr[N_ALPHA][N_BETA];
  double intAlphaBeta1, intAlphaBeta1_fr, intBeta1, intBeta1_fr;
  double cos_etv, cos_eti, sin_etv, sin_eti, cos_2etv, sin_2etv, cos_2eti, sin_2eti;
  double f_pdf0 = 0., mu_n0 = 0.;

  // the previous calculation will be cached
  // in case this function is called multiple times with the same parameters and angles
  // the matrix will only be calculated once
  // this was adapted from Litvinov's FORTRAN code and might not be necessary for MYSTIC
  static bool   gaussInit = false;
  static double sza_old = 0., vza_old = 0., saa_old = 0., vaa_old = 0.;
  static double params_old[5] = {0.}, rMat_old[4 * 4] = {0.}, f_pdf1[N_BETA];
  static double alpha[N_ALPHA], mu_beta[N_BETA], beta[N_BETA], Gw_p[N_ALPHA], Gw_t[N_BETA];

  // azimuth bounds for the gaussian quadrature
  double low_lp = 0.;
  double up_lp  = 2. * M_PI;

  // We need to calculate the weights and roots of the Gaussian quadratures just one time
  // which is done at the first invocation of the function
  if (!gaussInit) {
    gaussLegendre (N_ALPHA, N_ALPHA, low_lp, up_lp, &alpha[0], &Gw_p[0]);
    gaussLegendre (N_BETA, N_BETA, 0., 1., &mu_beta[0], &Gw_t[0]);
    for (unsigned i = 0; i < N_BETA; i++)
      beta[i] = acos (mu_beta[i]);
    gaussInit = true;
  }

  // check if there are no significant changes and return previously calculated matrix if necessary
  bool noChanges = ((fabs (sza_old - solar_za) <= eps) && (fabs (vza_old - view_za) <= eps) && (fabs (saa_old - solar_aa) <= eps) &&
                    (fabs (vaa_old - view_aa) <= eps));

  for (unsigned short i = 0; i < numParams; i++)
    noChanges = noChanges && (fabs (params_old[i] - params[i]) <= eps);

  if (noChanges) {
    for (unsigned short i = 0; i < matrixType * matrixType; i++)
      rMat[i] = rMat_old[i];
  } else {
    // check which parameters were provided and use default values otherwise
    alb_s   = params[0];
    sig1    = params[1];
    fren_fr = params[2];
    switch (numParams) {
    case 3:
      k_shad = sig1;
      break;
    case 4:
      k_shad = params[3];
      break;
    case 5:
      k_shad = params[3];
      sig2   = params[4];
      break;
    default:
      printf ("invalid amount of parameters numParams=%d", numParams);
      exit (-1);
    }

    // some adjustments to how the angles are interpreted to make it work with MYSTIC's angular definitions
    // Note that Litvinov's code might have had a bug here that was corrected
    mu0 = solar_za;
    muv = view_za;

    view_aa = (solar_aa - view_aa) + 2. * M_PI;
    while (view_aa >= 2. * M_PI)
      view_aa = view_aa - 2. * M_PI;
    solar_aa       = 0.;
    double phi_raw = (view_aa - solar_aa);

    // phii: incident azimuth angle according to MYSTIC
    // del_phi: incident relative azimuth
    // sin_tv/sin_ti: sin(theta_view)/sin(theta_incident)
    //    theta is the zenith angle
    // tet: angle between incident and viewing angle
    phii    = solar_aa + M_PI;
    del_phi = view_aa - phii;
    sin_ti  = sqrt (1. - mu0 * mu0);
    sin_tv  = sqrt (1. - muv * muv);
    cos_tet = -(muv * mu0 + sin_tv * sin_ti * cos (view_aa - solar_aa));
    tet     = acos (cos_tet);

    // calculate the Fresnel matrix, shadow, p_rel and delta_p_22 (see Bachelor thesis)
    reflFresnelMy (matrixType, mu0, muv, phi_raw / (M_PI / 180.), m_ind, rMat, NULL);
    shadow   = shadowFunc (mu0, muv, k_shad * sig1);
    r11_elem = (1. - alb_m * sin (tet / par_free) * sin (tet / par_free));
    for (unsigned i = 0; i < matrixType * matrixType; i++)
      rMat[i] = fren_fr * rMat[i] * shadow / mu0 / muv;

    pRel    = -(1. - cos_tet * cos_tet);
    sin_tet = sqrt (1. - cos_tet * cos_tet);
    //depol_22 = 0.8*cos(tet/4.)*cos(tet/4.);
    depol_22 = 0.8 * (1. - sin (tet / 6.));

    // if the surface is not flat in scale 1
    if (sig1 >= sig_min) {
      // If sig1 is the same then f_pdf1 is not calculated and taken from previous
      // values since f_pdf1 is saved

      // perform Gaussian quadrature to approximate the integrals
      double sum_mu_beta = 0.;
      for (unsigned i = 0; i < N_BETA; i++)
        sum_mu_beta += mu_beta[i];

      if (fabs (params_old[1] - params[1]) > eps) {
        for (unsigned i = 0; i < N_BETA; i++) {
          mu_n1 = mu_beta[i];
          f_pdf1[i] =
            exp (-(1. - mu_n1 * mu_n1) / (2. * sig1) / mu_n1 / mu_n1) / M_PI / (2. * sig1) / mu_n1 / mu_n1 / mu_n1 / mu_n1;
        }
      }

      double sum_f_pdf1 = 0.;
      for (unsigned i = 0; i < N_BETA; i++)
        sum_f_pdf1 += f_pdf1[i];

      for (unsigned i = 0; i < N_BETA; i++) {
        mubv  = mu_beta[i] * muv;
        mubi  = mu_beta[i] * mu0;
        sinbv = sqrt (1. - mu_beta[i] * mu_beta[i]) * sin_tv;
        sinbi = sqrt (1. - mu_beta[i] * mu_beta[i]) * sin_ti;
        for (unsigned j = 0; j < N_ALPHA; j++) {
          mu1v = mubv + sinbv * cos (alpha[j] - view_aa);
          mu2i = mubi + sinbi * cos (alpha[j] - solar_aa);
          if ((mu1v >= 0.f) && (mu2i >= 0.f)) {
            randF1_fr[j][i] = pdfFresnel (alpha[j], beta[i], sig2, &mu0, phii, &muv, view_aa, mu1v, mu2i);
            randF0          = mu1v * mu2i * f_pdf1[i];
            randF1[j][i]    = randF0;
            randF1_fr[j][i] *= randF0;
          } else {
            randF1[j][i]    = 0.;
            randF1_fr[j][i] = 0.;
          }
        }
      }

      double sum_randF1 = 0.;
      for (unsigned i = 0; i < N_BETA; i++)
        for (unsigned j = 0; j < N_ALPHA; j++)
          sum_randF1 += randF1[j][i];
      double sum_randF1_fr = 0.;
      for (unsigned i = 0; i < N_BETA; i++)
        for (unsigned j = 0; j < N_ALPHA; j++)
          sum_randF1_fr += randF1_fr[j][i];

      intAlphaBeta1    = 0.;
      intAlphaBeta1_fr = 0.;
      for (unsigned i = 0; i < N_BETA; i++) {
        intBeta1    = 0;
        intBeta1_fr = 0;
        for (unsigned j = 0; j < N_ALPHA; j++) {
          intBeta1 += Gw_p[j] * randF1[j][i];
          intBeta1_fr += Gw_p[j] * randF1_fr[j][i];
        }
        intAlphaBeta1 += intBeta1 * Gw_t[i];
        intAlphaBeta1_fr += intBeta1_fr * Gw_t[i];
      }

      // calculate actual matrix elements
      for (unsigned i = 0; i < matrixType * matrixType; i++)
        rMat[i] *= intAlphaBeta1_fr;
      rMat[0] = alb_s * r11_elem * CB_peak_new (tet) * intAlphaBeta1 / mu0 / muv * shadow + rMat[0];
      if (matrixType >= MATRIX3) {
        if (fabs (sin_tet) < sin_eps)
          sin_tet = sin_eps;
        if (fabs (sin_tv) < sin_eps)
          sin_tv = sin_eps;
        if (fabs (sin_ti) < sin_eps)
          sin_ti = sin_eps;

        // These sines and cosines come from the rotation of the Stokes vector
        cos_etv = (-mu0 - muv * cos_tet) / sin_tet / sin_tv;
        cos_eti = (muv + mu0 * cos_tet) / sin_tet / sin_ti;  // checked "-"

        sin_etv = sin (del_phi) * sin_ti / sin_tet;
        sin_eti = sin (del_phi) * sin_tv / sin_tet;

        cos_2etv = cos_etv * cos_etv - sin_etv * sin_etv;
        sin_2etv = 2. * cos_etv * sin_etv;

        cos_2eti = cos_eti * cos_eti - sin_eti * sin_eti;
        sin_2eti = 2. * cos_eti * sin_eti;

        // put it all together
        MatIndexFunc ind = matIndexByType[matrixType];
        rMat[ind (1, 1)] =
          depol_22 * alb_s * r11_elem * CB_peak_new (tet) * intAlphaBeta1 / mu0 / muv * shadow * cos_2etv * cos_2eti +
          rMat[ind (1, 1)];
        rMat[ind (1, 0)] = depol * pRel * cos_2etv * intAlphaBeta1 / mu0 / muv * shadow + rMat[ind (1, 0)];
        rMat[ind (0, 1)] = depol * pRel * cos_2eti * intAlphaBeta1 / mu0 / muv * shadow + rMat[ind (0, 1)];
        rMat[ind (2, 0)] = -depol * pRel * sin_2etv * intAlphaBeta1 / mu0 / muv * shadow + rMat[ind (2, 0)];
        rMat[ind (0, 2)] = depol * pRel * sin_2eti * intAlphaBeta1 / mu0 / muv * shadow + rMat[ind (0, 2)];  // Checked "-"
      }
    } else {
      // Same with a completely flat surface in scale 1
      //TODO f_pdf1?
      for (unsigned i = 0; i < matrixType * matrixType; i++)
        rMat[i] *= f_pdf0 * M_PI / 4. / mu_n0;
      rMat[0] = alb_s * r11_elem * CB_peak_new (tet) * shadow + rMat[0];
      if (matrixType >= MATRIX3) {
        if (fabs (sin_tet) < sin_eps)
          sin_tet = sin_eps;
        if (fabs (sin_tv) < sin_eps)
          sin_tv = sin_eps;
        if (fabs (sin_ti) < sin_eps)
          sin_ti = sin_eps;

        cos_etv = (-mu0 - muv * cos_tet) / sin_tet / sin_tv;
        cos_eti = (muv + mu0 * cos_tet) / sin_tet / sin_ti;  // checked "-"

        sin_etv = sin (del_phi) * sin_ti / sin_tet;
        sin_eti = sin (del_phi) * sin_tv / sin_tet;

        cos_2etv = cos_etv * cos_etv - sin_etv * sin_etv;
        sin_2etv = 2. * cos_etv * sin_etv;

        cos_2eti = cos_eti * cos_eti - sin_eti * sin_eti;
        sin_2eti = 2. * cos_eti * sin_eti;

        MatIndexFunc ind = matIndexByType[matrixType];
        rMat[ind (1, 1)] = depol_22 * alb_s * r11_elem * CB_peak_new (tet) * shadow * cos_2etv * cos_2eti + rMat[ind (1, 1)];
        rMat[ind (1, 0)] = depol * pRel * cos_2etv * shadow + rMat[ind (1, 0)];
        rMat[ind (0, 1)] = depol * pRel * cos_2eti * shadow + rMat[ind (0, 1)];
        rMat[ind (2, 0)] = -depol * pRel * sin_2etv * shadow + rMat[ind (2, 0)];
        rMat[ind (0, 2)] = depol * pRel * sin_2eti * shadow + rMat[ind (0, 2)];  // Checked "-"
      }
    }
    // return matrix
    for (unsigned short i = 0; i < matrixType * matrixType; i++)
      rMat_old[i] = rMat[i];
    for (unsigned short i = 0; i < numParams; i++)
      params_old[i] = params[i];
    sza_old = solar_za;
    vza_old = view_za;
    saa_old = solar_aa;
    vaa_old = view_aa;
  }
  return 0;  //rMat;
}

// utilizes reflFresnelMy to calculate the Fresnel matrix and then multiplies with some factors (see Bachelor thesis)
int maign_breon_xh (enum MatrixType matrixType,
                    double          solar_za,
                    double          view_za,
                    double          phi,
                    double _Complex m_ind,
                    const double c_mgn,
                    const double arvi,
                    double*      rMat) {
  double alpha;
  reflFresnelMy (matrixType, solar_za, view_za, phi, m_ind, rMat, &alpha);
  const double coefficient = c_mgn * exp (-tan (alpha * 0.5)) * 0.25 / (solar_za + view_za) * exp (-arvi);

  for (size_t i = 0; i < matrixType * matrixType; i++)
    rMat[i] = rMat[i] * coefficient;
  return 0;
}

// This does exactly the same as maign_breon_xh, but incorporates reflFresnelMy
// Might be faster?
int maign_breon (double solar_za, double view_za, double phi, double _Complex m_ind, double c_mgn, double arvi, double* rMat) {
  double st_min = 0.001, sin_eps = 0.00000001;
  double zz, alpha, uu, coefficient, mu_alpha;
  double _Complex r1, r2, r_par1, r_par2, r_perp1, r_perp2, r_par, r_perp, r_par_per1, r_par_per2;
  double r_par_sq, r_perp_sq;
  double del_phi, sin_ti, sin_tv, sin_tet, cos_tet, cos_etv, cos_eti, sin_eti, sin_etv, cos_2etv, sin_2etv, cos_2eti, sin_2eti;
  double L_v[4 * 4], L_inc[4 * 4], F_fr[4 * 4], mu0, muv, tet0v, sin_phi;

  const double dpi_180 = M_PI / 180.0;
  mu0                  = solar_za;
  muv                  = view_za;

  sin_tv  = sqrt (1.0 - muv * muv);
  sin_ti  = sqrt (1.0 - mu0 * mu0);
  zz      = muv * mu0 + sin_ti * sin_tv * cos (phi * dpi_180);
  cos_tet = -zz;
  sin_tet = sqrt (1.0 - cos_tet * cos_tet);

  if (fabs (sin_tv) < sin_eps) {
    tet0v   = acos (muv);
    tet0v   = fabs (tet0v - st_min * dpi_180);
    muv     = cos (tet0v);
    sin_tv  = sqrt (1.0 - muv * muv);
    zz      = muv * mu0 + sin_ti * sin_tv * cos (phi * dpi_180);
    cos_tet = -zz;
    sin_tet = sqrt (1.0 - cos_tet * cos_tet);
  }

  if (fabs (sin_ti) < sin_eps) {
    tet0v   = acos (mu0);
    tet0v   = fabs (tet0v - st_min * dpi_180);
    mu0     = cos (tet0v);
    sin_ti  = sqrt (1.0 - mu0 * mu0);
    zz      = muv * mu0 + sin_ti * sin_tv * cos (phi * dpi_180);
    cos_tet = -zz;
    sin_tet = sqrt (1.0 - cos_tet * cos_tet);
  }

  if (fabs (sin_tet) < sin_eps) {
    tet0v   = acos (mu0);
    tet0v   = fabs (tet0v - st_min * dpi_180);
    muv     = cos (tet0v);
    sin_tv  = sqrt (1.0 - muv * muv);
    zz      = muv * mu0 + sin_ti * sin_tv * cos (phi * dpi_180);
    cos_tet = -zz;
    sin_tet = sqrt (1.0 - cos_tet * cos_tet);
  }

#if VERBOSE_LEVEL > 1

  if (fabs (sin_tet) < sin_eps)
    fprintf (stderr, "sin_tet: %f\n", sin_tet);
  if (fabs (sin_tv) < sin_eps)
    fprintf (stderr, "sin_tv: %f\n", sin_tv);
  if (fabs (sin_ti) < sin_eps)
    fprintf (stderr, "sin_ti: %f\n", sin_ti);

#endif

  alpha       = acos (zz);
  uu          = tan (alpha * 0.50);
  mu_alpha    = cos (alpha * 0.50);
  coefficient = c_mgn * exp (-uu) * 0.250 / (mu0 + muv) * exp (-arvi);

  r1     = m_ind * m_ind * mu_alpha;
  r2     = csqrt (m_ind * m_ind - 1.0 + mu_alpha * mu_alpha);
  r_par1 = r1 - r2;
  r_par2 = r1 + r2;

  r_perp1 = mu_alpha - r2;
  r_perp2 = mu_alpha + r2;

  r_par  = r_par1 / r_par2;
  r_perp = r_perp1 / r_perp2;

  r_par_sq   = creal (r_par * conj (r_par));
  r_perp_sq  = creal (r_perp * conj (r_perp));
  r_par_per1 = r_par * conj (r_perp);
  r_par_per2 = conj (r_par_per1);

  for (unsigned i = 0; i < 4 * 4; i++)
    F_fr[i] = 0.0;

  MatIndexFunc ind = matIndexByType[MATRIX4];
  F_fr[ind (0, 0)] = (r_par_sq + r_perp_sq) * coefficient * 0.5;
  F_fr[ind (1, 1)] = F_fr[ind (0, 0)];
  F_fr[ind (1, 0)] = (r_par_sq - r_perp_sq) * coefficient * 0.5;
  F_fr[ind (0, 1)] = F_fr[ind (1, 0)];
  F_fr[ind (2, 2)] = creal (r_par_per1 + r_par_per2) * coefficient * 0.5;
  F_fr[ind (3, 3)] = F_fr[ind (2, 2)];
  //TODO
  F_fr[ind (2, 3)] = creal (-I * (r_par_per1 - r_par_per2)) * coefficient * 0.5;
  F_fr[ind (3, 2)] = -F_fr[ind (2, 3)];

  del_phi = phi * dpi_180 - M_PI;
  cos_etv = (-mu0 - muv * cos_tet) / sin_tet / sin_tv;
  cos_eti = (muv + mu0 * cos_tet) / sin_tet / sin_ti;

  sin_etv = sin (del_phi) * sin_ti / sin_tet;
  sin_eti = sin (del_phi) * sin_tv / sin_tet;

  cos_2etv = cos_etv * cos_etv - sin_etv * sin_etv;
  sin_2etv = 2.0 * cos_etv * sin_etv;

  cos_2eti = cos_eti * cos_eti - sin_eti * sin_eti;
  sin_2eti = 2.0 * cos_eti * sin_eti;

  for (unsigned i = 0; i < 4 * 4; i++) {
    L_v[i]   = 0.0;
    L_inc[i] = 0.0;
  }

  sin_phi = sin (phi * dpi_180);
  if (fabs (sin_phi) < sin_eps / 10.0) {
    //TODO
    L_v[ind (0, 0)] = 1.0;
    L_v[ind (1, 1)] = 1.0;
    L_v[ind (1, 2)] = 0;
    L_v[ind (2, 2)] = 1.0;
    L_v[ind (2, 1)] = 0;
    L_v[ind (3, 3)] = 1.0;

    L_inc[ind (0, 0)] = 1.0;
    L_inc[ind (1, 1)] = 1.0;
    L_inc[ind (1, 2)] = 0;
    L_inc[ind (2, 2)] = 1.0;
    L_inc[ind (2, 1)] = 0;
    L_inc[ind (3, 3)] = 1.0;
  } else {
    L_v[ind (0, 0)] = 1.0;
    L_v[ind (1, 1)] = cos_2etv;
    L_v[ind (1, 2)] = sin_2etv;
    L_v[ind (2, 2)] = cos_2etv;
    L_v[ind (2, 1)] = -sin_2etv;
    L_v[ind (3, 3)] = 1.0;

    L_inc[ind (0, 0)] = 1.0;
    L_inc[ind (1, 1)] = cos_2eti;
    L_inc[ind (1, 2)] = sin_2eti;
    L_inc[ind (2, 2)] = cos_2eti;
    L_inc[ind (2, 1)] = -sin_2eti;
    L_inc[ind (3, 3)] = 1.0;
  }
  double* out = matMul (MATRIX4, F_fr, L_inc, NULL);
  matMul (MATRIX4, L_v, out, rMat);

  free (out);
  return 0;
}

// coherent backscattering for infinite medium (see Bachelor thesis)
double CB_peak_new (double zenithAngle) {
  double eps = 0.00001, kl_ext = 30., Coh_peak = 0.5;
  double q = 2. * kl_ext * cos (zenithAngle * 0.5);

  if (fabs (zenithAngle - M_PI) > eps)
    return 1. + Coh_peak * atan (q) / q;
  else
    return 1. + Coh_peak;
}

// Shadowing function (see Bachelor thesis)
double shadowFunc (double mu0, double muv, double sigma2) {
  double s1, s2, s3, xi, xxi, dcot, t1, t2, shadowi, shadowr;

  // another shadowing
  s1      = sqrt (2. * sigma2 / M_PI);
  s3      = 1. / (sqrt (2 * sigma2));
  s2      = s3 * s3;
  xi      = mu0;
  xxi     = xi * xi;
  dcot    = xi / sqrt (1. - xxi);
  t1      = exp (-dcot * dcot * s2);
  t2      = erfc (dcot * s3);
  shadowi = 0.5 * (s1 * t1 / dcot - t2);
  xi      = muv;
  xxi     = xi * xi;
  dcot    = xi / sqrt (1. - xxi);
  t1      = exp (-dcot * dcot * s2);
  t2      = erfc (dcot * s3);
  shadowr = 0.5 * (s1 * t1 / dcot - t2);
  return 1. / (1. + shadowi + shadowr);
}

// Calculate Fresnel matrix
double* reflFresnelMy (enum MatrixType matrixType,
                       double          mu01,
                       double          muv1,
                       double          phi,
                       double _Complex m,
                       double* rMat_fr,
                       double* alpha1) {
  static const double st_min = 0.001, sin_eps = 0.00000001;
  //sin_eps = sin(st_min*M_PI/180.);
  if (rMat_fr == NULL)
    rMat_fr = malloc (sizeof (double) * matrixType * matrixType);
  double mu0, muv, sin_tv, sin_ti, zz, cos_tet, sin_tet, tet0v, sin_phi, cos_etv, cos_eti, sin_eti, sin_etv, cos_2etv, sin_2etv,
    cos_2eti, sin_2eti;
  double alpha, mu_alpha;
  double _Complex r1, r2, r_par1, r_par2, r_perp1, r_perp2, r_par, r_perp, r_par_per1, r_par_per2;
  double r_par_sq, r_perp_sq;
  double l_v[4 * 4], l_inc[4 * 4], f_fr[4 * 4];
  double del_phi;

  mu0 = mu01;
  muv = muv1;

  sin_tv  = sqrt (1. - muv * muv);
  sin_ti  = sqrt (1. - mu0 * mu0);
  zz      = muv * mu0 + sin_ti * sin_tv * cos (phi * M_PI / 180.);
  cos_tet = -zz;
  sin_tet = sqrt (1. - cos_tet * cos_tet);

  if (fabs (sin_tv) < sin_eps) {
    tet0v   = acos (muv);
    tet0v   = fabs (tet0v - st_min * M_PI / 180.);
    muv     = cos (tet0v);
    sin_tv  = sqrt (1. - muv * muv);
    zz      = muv * mu0 + sin_ti * sin_tv * cos (phi * M_PI / 180.);
    cos_tet = -zz;
    sin_tet = sqrt (1. - cos_tet * cos_tet);
  }

  if (fabs (sin_ti) < sin_eps) {
    tet0v   = acos (mu0);
    tet0v   = fabs (tet0v - st_min * M_PI / 180.);
    mu0     = cos (tet0v);
    sin_ti  = sqrt (1. - mu0 * mu0);
    zz      = muv * mu0 + sin_ti * sin_tv * cos (phi * M_PI / 180.);
    cos_tet = -zz;
    sin_tet = sqrt (1. - cos_tet * cos_tet);
  }

  if (fabs (sin_tet) < sin_eps) {
    tet0v   = acos (mu0);
    tet0v   = fabs (tet0v - st_min * M_PI / 180.);
    muv     = cos (tet0v);
    sin_tv  = sqrt (1. - muv * muv);
    zz      = muv * mu0 + sin_ti * sin_tv * cos (phi * M_PI / 180.);
    cos_tet = -zz;
    sin_tet = sqrt (1. - cos_tet * cos_tet);
  }

#if VERBOSE_LEVEL > 1
  if (fabs (sin_tet) < sin_eps)
    printf ("\nsin_tet = %4.2f\n", sin_tet);
  if (fabs (sin_tv) < sin_eps)
    printf ("\nsin_tv = %4.2f\n", sin_tv);
  if (fabs (sin_ti) < sin_eps)
    printf ("\nsin_ti = %4.2f\n", sin_ti);
#endif

  alpha = acos (zz);
  if (alpha1 != NULL)
    *alpha1 = alpha;
  mu_alpha = cos (alpha * 0.5);

  r1     = m * m * mu_alpha;
  r2     = csqrt (m * m - 1. + mu_alpha * mu_alpha);
  r_par1 = r1 - r2;
  r_par2 = r1 + r2;

  r_perp1 = mu_alpha - r2;
  r_perp2 = mu_alpha + r2;

  r_par  = r_par1 / r_par2;
  r_perp = r_perp1 / r_perp2;

  r_par_sq   = (double) (r_par * conj (r_par));
  r_perp_sq  = (double) (r_perp * conj (r_perp));
  r_par_per1 = r_par * conj (r_perp);
  r_par_per2 = conj (r_par_per1);

  for (unsigned i = 0; i < 4 * 4; i++)
    f_fr[i] = 0.;
  f_fr[0] = (r_par_sq + r_perp_sq) * 0.5;

  if (matrixType >= MATRIX3) {
    MatIndexFunc ind = matIndexByType[matrixType];
    f_fr[ind (1, 1)] = f_fr[ind (0, 0)];
    f_fr[ind (1, 0)] = (r_par_sq - r_perp_sq) * 0.5;
    f_fr[ind (0, 1)] = f_fr[ind (1, 0)];
    f_fr[ind (2, 2)] = creal (r_par_per1 + r_par_per2) * 0.5;

    if (matrixType == MATRIX4) {
      f_fr[ind (3, 3)] = f_fr[ind (2, 2)];
      f_fr[ind (2, 3)] = creal (-I * (r_par_per1 - r_par_per2)) * 0.5;
      f_fr[ind (3, 2)] = -f_fr[ind (2, 3)];
    }

    del_phi = phi * M_PI / 180. - M_PI;
    cos_etv = (-mu0 - muv * cos_tet) / sin_tet / sin_tv;
    cos_eti = (muv + mu0 * cos_tet) / sin_tet / sin_ti;  // checked "-"

    sin_etv = sin (del_phi) * sin_ti / sin_tet;
    sin_eti = sin (del_phi) * sin_tv / sin_tet;

    cos_2etv = cos_etv * cos_etv - sin_etv * sin_etv;
    sin_2etv = 2. * cos_etv * sin_etv;

    cos_2eti = cos_eti * cos_eti - sin_eti * sin_eti;
    sin_2eti = 2. * cos_eti * sin_eti;

    for (unsigned i = 0; i < 4 * 4; i++) {
      l_v[i]   = 0.;
      l_inc[i] = 0.;
    }

    sin_phi = sin (phi * M_PI / 180.);
    if (fabs (sin_phi) < sin_eps / 10.) {
      MatIndexFunc ind = matIndexByType[matrixType];
      l_v[ind (0, 0)]  = 1.;
      l_v[ind (1, 1)]  = 1.;
      l_v[ind (1, 2)]  = 0.;
      l_v[ind (2, 2)]  = 1.;
      l_v[ind (2, 1)]  = 0.;

      l_inc[ind (0, 0)] = 1.;
      l_inc[ind (1, 1)] = 1.;
      l_inc[ind (1, 2)] = 0.;
      l_inc[ind (2, 2)] = 1.;
      l_inc[ind (2, 1)] = 0.;

      if (matrixType == MATRIX4) {
        l_v[ind (3, 3)]   = 1.;
        l_inc[ind (3, 3)] = 1.;
      }

    } else {
      l_v[ind (0, 0)] = 1.;
      l_v[ind (1, 1)] = cos_2etv;
      l_v[ind (1, 2)] = sin_2etv;
      l_v[ind (2, 2)] = cos_2etv;
      l_v[ind (2, 1)] = -sin_2etv;

      l_inc[ind (0, 0)] = 1.;
      l_inc[ind (1, 1)] = cos_2eti;
      l_inc[ind (1, 2)] = sin_2eti;
      l_inc[ind (2, 2)] = cos_2eti;
      l_inc[ind (2, 1)] = -sin_2eti;

      if (matrixType == MATRIX4) {
        l_v[ind (3, 3)]   = 1.;
        l_inc[ind (3, 3)] = 1.;
      }
    }
    double temp[4 * 4];
    rMat_fr = matMul (matrixType, l_v, matMul (matrixType, f_fr, l_inc, temp), rMat_fr);
  } else {
    for (unsigned i = 0; i < 4 * 4; i++) {
      rMat_fr[i] = f_fr[i];
    }
  }
  return rMat_fr;
}

double* gaussLegendre (unsigned numDimensions, unsigned numGaussLegPoints, double a, double b, double* abscissas, double* weights) {
  /***********************************************************************
* Given the lower and upper limits of integration a and b,
* and given the number of Gauss-Legendre points numGaussLegPoints,
* this routine returns through array abscissas the abscissas and through
* array weights the weights of the Gauss-Legendre quadrature formula.
* eps is the desired accuracy of the abscissas.
* This routine is documented further in:
*   W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)
*   page 125 ISBN 0-521-30811-9
**********************************************************************/
  //IMPLICIT Real(8) (a-h,o-z)
  double   eps = 1.e-14;
  unsigned m;
  if (abscissas == NULL)
    abscissas = malloc (sizeof (float) * numDimensions);
  if (weights == NULL)
    weights = malloc (sizeof (float) * numDimensions);
  double xm, xl, z, zdd, p1, p2, p3, pp, z1;

  m  = (numGaussLegPoints + 1) / 2;
  xm = 0.5 * (a + b);
  xl = 0.5 * (b - a);

  for (unsigned i = 1; i <= m; i++) {
    z   = cos (M_PI * ((double) i - 0.25) / ((double) numGaussLegPoints + 0.5));
    zdd = 1.;
    while (zdd > eps) {
      p1 = 1.;
      p2 = 0.;
      for (unsigned j = 1; j <= numGaussLegPoints; j++) {
        p3 = p2;
        p2 = p1;
        p1 = (((double) (2 * j) - 1.) * z * p2 - ((double) j - 1.) * p3) / (double) j;
      }
      pp  = numGaussLegPoints * (z * p1 - p2) / (z * z - 1.);
      z1  = z;
      z   = z1 - p1 / pp;
      zdd = fabs (z - z1);
    }

    abscissas[i - 1]                 = xm - xl * z;
    abscissas[numGaussLegPoints - i] = xm + xl * z;
    weights[i - 1]                   = 2. * xl / ((1. - z * z) * pp * pp);
    weights[numGaussLegPoints - i]   = weights[i - 1];
  }
  return abscissas;
}

double pdfFresnel (double  alpha,
                   double  beta,
                   double  sigma2,
                   double* dmui1,
                   double  phii,
                   double* dmur1,
                   double  phir,
                   double  mu1r,
                   double  mu2i) {
  double r_rot[3][3];

  if (fabs (*dmui1 - 1.) < 1e-10)
    *dmui1 = 0.999999999999;
  if (fabs (*dmur1 - 1.) < 1e-10)
    *dmur1 = 0.999999999999;

  double dcosi = cos (phii);
  double dsini = sin (phii);
  double dcosr = cos (phir);
  double dsinr = sin (phir);

  double dsi = sqrt (1. - *dmui1 * *dmui1);
  double dsr = sqrt (1. - *dmur1 * *dmur1);

  //coordinates of incident normal
  double vi1o = dsi * dcosi;
  double vi2o = dsi * dsini;
  double vi3o = -*dmui1;

  //coordinates of viewing normal
  double vr1o = dsr * dcosr;
  double vr2o = dsr * dsinr;
  double vr3o = *dmur1;

  //rotation matrix
  r_rot[0][0] = cos (alpha) * cos (beta);
  r_rot[0][1] = sin (alpha) * cos (beta);
  r_rot[0][2] = -sin (beta);

  r_rot[1][0] = -sin (alpha);
  r_rot[1][1] = cos (alpha);
  r_rot[1][2] = 0.;

  r_rot[2][0] = cos (alpha) * sin (beta);
  r_rot[2][1] = sin (alpha) * sin (beta);
  r_rot[2][2] = cos (beta);

  //transformation on rotation
  double vi1 = r_rot[0][0] * vi1o + r_rot[0][1] * vi2o + r_rot[0][2] * vi3o;
  double vi2 = r_rot[1][0] * vi1o + r_rot[1][1] * vi2o + r_rot[1][2] * vi3o;
  double vi3 = r_rot[2][0] * vi1o + r_rot[2][1] * vi2o + r_rot[2][2] * vi3o;

  // printf("%10.5f, %10.5f, %20.16f, %20.16f, %20.16f", acos(*dmui1)*180.f/M_PI, phii*180.f/M_PI, vi1o, vi2o, vi3o);
  //local solar zenith angle
  double dmui = -vi3;

  //*****************
  //  dmui=mu2i
  if (fabs (dmui - mu2i) > 0.00001)
    printf ("Error dmui 1 %4.2f %4.2f", dmui, mu2i);

  if (dmui < 0)
    printf ("Error dmui 2");
  if (fabs (dmui - 1.) < 1. - 10.)
    dmui = 0.999999999999;
  dsi = sqrt (1. - dmui * dmui);

  //local solar azimuth angle
  dcosi = vi1 / dsi;
  dsini = vi2 / dsi;
  //*****************

  if (fabs (fabs (dsini) - sqrt (1. - dcosi * dcosi)) > 0.00001)
    printf ("Error dsini %4.2f %4.2f", dsini, sqrt (1 - dcosi * dcosi));

  double vr1 = r_rot[0][0] * vr1o + r_rot[0][1] * vr2o + r_rot[0][2] * vr3o;
  double vr2 = r_rot[1][0] * vr1o + r_rot[1][1] * vr2o + r_rot[1][2] * vr3o;
  double vr3 = r_rot[2][0] * vr1o + r_rot[2][1] * vr2o + r_rot[2][2] * vr3o;

  //local viewing zenith angle
  double dmur = vr3;

  //*****************
  if (fabs (dmur - mu1r) > 0.00001)
    printf ("Error dmur %4.2f %4.2f", dmur, mu1r);

  if (dmur < 0)
    printf ("Error dmur");
  if (fabs (dmur - 1.) < 1. - 10.)
    dmur = 0.999999999999;
  dsr = sqrt (1. - dmur * dmur);
  //local viewing azimuth angle
  dcosr = vr1 / dsr;
  dsinr = vr2 / dsr;

  if (fabs (fabs (dsinr) - sqrt (1. - dcosr * dcosr)) > 0.0001) {
    printf ("Error dsinr %4.2f %4.2f %4.2f", acos (*dmur1) * 180.f / M_PI, dsinr, sqrt (1 - dcosr * dcosr));
    printf ("%16.5f %16.5f %16.5f %16.5f", vr2, vr1o, vr2o, vr3o);
  }

  //LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION
  double unit1 = vi1 - vr1;
  double unit2 = vi2 - vr2;
  double unit3 = vi3 - vr3;
  double fact1 = unit1 * unit1 + unit2 * unit2 + unit3 * unit3;
  //double factor = sqrt(1./fact1);

  double rdz2 = unit3 * unit3;
  double rdz4 = rdz2 * rdz2;

  double dcoeff = 1. / (4. * dmui * dmur * rdz4 * 2. * sigma2);  //original
  double dex    = -(unit1 * unit1 + unit2 * unit2) / (2. * sigma2 * rdz2);
  dex           = exp (dex);
  dcoeff        = dcoeff * fact1 * fact1 * dex;  // mu_n^4  //original

  //  SHADOWING
  double p    = acos (-1.);
  double s1   = sqrt (2. * sigma2 / p);
  double s3   = 1. / (sqrt (2. * sigma2));
  double s2   = s3 * s3;
  double xi   = dmui;
  double xxi  = xi * xi;
  double dcot = xi / sqrt (1. - xxi);
  double t1   = exp (-dcot * dcot * s2);

  //t2 = erfc_func(dcot*s3)

  double t2      = erfc (dcot * s3);
  double shadowi = 0.5 * (s1 * t1 / dcot - t2);
  xi             = dmur;
  xxi            = xi * xi;
  dcot           = xi / sqrt (1. - xxi);
  t1             = exp (-dcot * dcot * s2);

  //t2 = erfc_func(dcot*s3)
  t2 = erfc (dcot * s3);
  //t2 = 1.-erf(dcot*s3);
  double shadowr = 0.5 * (s1 * t1 / dcot - t2);
  double shadow  = 1. / (1. + shadowi + shadowr);

  return dcoeff * shadow;
}

int reflection_matrix_litvinov_mystic (double mu_esc,
                                       double mu_inc,
                                       double phi_inc,
                                       double phi_esc,
                                       double _Complex refr_index,
                                       double    albedo,
                                       double    rms_slope,
                                       double    fresnel,
                                       double    shadowing,
                                       double    rms_slope2,
                                       int       backward,
                                       double*** refl_mat) {
  // TODO remove index function and use C-style multidimensional arrays
  static MatIndexFunc ind = NULL;
  if (ind == NULL)
    ind = matIndexByType[MATRIX4];
  double litvinov_mat[4 * 4];
  double params[] = {albedo, rms_slope, fresnel, shadowing, rms_slope2};

  if (backward) {
    const int ierr =
      reflection_matrix_litvinov (MATRIX4, mu_esc, mu_inc, phi_esc, phi_inc + M_PI, 5, params, refr_index, litvinov_mat);
    CHKERR (ierr);
  } else {
    const int ierr =
      reflection_matrix_litvinov (MATRIX4, mu_inc, mu_esc, phi_inc + M_PI, phi_esc, 5, params, refr_index, litvinov_mat);
    CHKERR (ierr);
  }
  // TODO remove index function
  for (unsigned i = 0; i < 4; i++) {
    for (unsigned j = 0; j < 4; j++) {
      (*refl_mat)[i][j] = litvinov_mat[ind (i, j)];
    }
  }
  // TODO remove index function
  return 0;
}

int reflection_matrix_maignan_breon_mystic (double mu_esc,
                                            double mu_inc,
                                            double phi_inc,
                                            double phi_esc,
                                            double _Complex refr_index,
                                            double    c_maign,
                                            double    arvi,
                                            double    rpv_rho0,
                                            double    rpv_k,
                                            double    rpv_theta,
                                            double    rpv_scale,
                                            double    rpv_sigma,
                                            double    rpv_t1,
                                            double    rpv_t2,
                                            int       backward,
                                            double*** refl_mat) {
  // TODO remove index function and use C-style multidimensional arrays
  static MatIndexFunc ind = NULL;
  if (ind == NULL)
    ind = matIndexByType[MATRIX4];

  int     status;
  double* maign_mat;
  maign_mat = malloc (sizeof (double) * 4 * 4);
  if (backward) {
    status = maign_breon_xh (MATRIX4,
                             mu_esc,
                             mu_inc,
                             360. - ((phi_inc + M_PI) - phi_esc) / M_PI * 180.,
                             refr_index,
                             c_maign,
                             arvi,
                             maign_mat);
  } else {
    status = maign_breon_xh (MATRIX4,
                             mu_inc,
                             mu_esc,
                             360. - (phi_esc - (phi_inc + M_PI)) / M_PI * 180.,
                             refr_index,
                             c_maign,
                             arvi,
                             maign_mat);
  }
  maign_mat[ind (0, 0)] += rpv_brdf (rpv_rho0,
                                     rpv_k,
                                     rpv_theta,
                                     rpv_scale,
                                     rpv_sigma,
                                     rpv_t1,
                                     rpv_t2,
                                     mu_inc,
                                     mu_esc,
                                     fabs ((phi_esc - phi_inc) * 180.0 / M_PI));

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned j = 0; j < 4; j++) {
      (*refl_mat)[i][j] = maign_mat[ind (i, j)];
    }
  }
  free (maign_mat);
  return status;
}
