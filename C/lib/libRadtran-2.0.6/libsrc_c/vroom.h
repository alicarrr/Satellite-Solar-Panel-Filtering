#ifndef q0
#define q0
#if defined(q1)
extern"\103"{
#endif
#define q2 10.0
typedef struct scadis_struct{double q3;double q4;double q5;double dirold_dx[3]
;double d_phi;double epsfac;double mu_max;double mu_min;double q6[99],F_min[99
];}scadis_struct;int set_vroom_settings(int vroom,sample_struct*sample,int q7)
;int mc_vroom_check_and_verbose(sample_struct*sample,int q7,int q8);int 
mc_vroom_prepare(sample_struct*sample,atmosphere_struct*q9,float q10,int q11,
int q12,int q7);int mc_vroom_cloning(photon_struct*p,atmosphere_struct*q9,
sample_struct*sample,result_struct*q13,elevation_struct*q14,albedo_struct*q15,
t_triangular_surface*q16,surftemp_struct*q17,int*q18,int*q19,int q20,int q21,
float*q22,int q7);int mc_vroom_splitting_and_rr(photon_struct*p,
atmosphere_struct*q9,sample_struct*sample,result_struct*q13,elevation_struct*
q14,albedo_struct*q15,t_triangular_surface*q16,surftemp_struct*q17,int*q18,int
*q19,int q20,int q21,float*q22,int q7);double get_phase_max(pft**phase_max,int
 n_phase_max,double q23,int SC_mode);void cp_locest(locest_struct*q24,
locest_struct*q19,sample_struct*sample,int n_caoth);int mc_vroom_prep_DDIS(
sample_struct*sample,photon_struct*p,atmosphere_struct*q9,int*q25,int*q26,int*
q27,locest_struct*lest,scadis_struct*q28);int mu_scatter_special(
atmosphere_struct*q9,photon_struct*p,pft**phase_max,int n_phase_max,int q29,
double*mu,scadis_struct q28,int q25,int q27);int random_reflection_special(
sample_struct*sample,atmosphere_struct*q9,elevation_struct*q14,photon_struct*p
,pft**phase_max,int n_phase_max,double*mu,double*phi,double*q30,scadis_struct*
q28,locest_struct lest,int q31,int q25,int q26,int q27);int 
mc_vroom_scattering_calc_phases_and_jacobians(sample_struct*sample,
photon_struct*p,atmosphere_struct*q9,double q32,int q33,double*q34,double*q35,
double*q36);int mc_vroom_DDIS_weight_and_prep_stuff(sample_struct*sample,
atmosphere_struct*q9,double q23,double q37,double q38,int q25,int q27,double 
q34,double q35,double q36,locest_struct lest,scadis_struct q28,photon_struct*p
);int mc_vroom_set_mus_and_phis(sample_struct*sample,photon_struct*p,int q25,
locest_struct lest,scadis_struct q28,double*mu,double*q23,double*q37,double 
phi,double*q38);int calloc_hybrid3D_field(int****q39,atmosphere_struct*q9);
void free_hybrid3D_field(int****q39,atmosphere_struct*q9);
#if defined(q1)
}
#endif
#endif

