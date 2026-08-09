#include <mystic.h>

int nca_like_mystic (char*              nca_options,
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
                     float***           abs3d);

int nca3d_1 (char*              datapath,
             double             dx,
             double             dy,
             atmosphere_struct* atmos,
             double***          Edn_3d,
             double***          Eup_3d,
             double***          planck,
             double***          kabs_nca,
             double***          abs_nca);

int nca3d_2 (char*              datapath,
             double             dx,
             double             dy,
             atmosphere_struct* atmos,
             double***          Edn_3d,
             double***          Eup_3d,
             double***          planck,
             double***          kabs_nca,
             double***          abs_nca);

int nca3d_3 (char*              datapath,
             double             dx,
             double             dy,
             atmosphere_struct* atmos,
             double***          Edn_3d,
             double***          Eup_3d,
             double***          planck,
             double***          kabs_nca,
             double***          abs_nca);

double calc_emis (double kabs, double ds, double alpha, double beta);

double integrate_emis (double kabs, double a, double b, double c, double delta_1, double delta_2, double mu);

double integrate_abs (double B,
                      double L,
                      double kabs1,
                      double kabs2,
                      double a,
                      double b,
                      double c,
                      double delta_1,
                      double delta_2,
                      double mu);

double integrate_flux (double B_planck1,
                       double B_planck2,
                       double L,
                       double kabs1,
                       double kabs2,
                       double b,
                       double c,
                       double delta_1,
                       double delta_2,
                       double mu);

int get_index (double a, double b);

double pol_3d (double    tauxx,
               double    tauyy,
               double    tauzz,
               double*** esp_tab,
               double*   tau_arrx,
               double*   tau_arry,
               double*   tau_arrz,
               int       ntaux,
               int       ntauy,
               int       ntauz);

double pol_2d (double var1, double var2, double** tab1, double* asp, double* tau_var, int nasp, int ntau, int flag);
