#ifndef __GPE_ENGINE_API__
#define __GPE_ENGINE_API__

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

#include <complex.h>

EXTERN void gpe_get_lattice_api(int *_nx, int *_ny, int *_nz);
EXTERN int gpe_create_engine_api(double alpha, double beta, double dt, double npart);
EXTERN int gpe_set_user_params_api(int size, double *params);
EXTERN int gpe_set_psi_api(double t, double complex * psi);
EXTERN int gpe_normalize_psi_api();
EXTERN int gpe_energy_api(double *t, double *ekin, double *eint, double *eext);
EXTERN int gpe_evolve_api(int nt);
EXTERN int gpe_get_psi_api(double *t, double complex * psi);
EXTERN int gpe_destroy_engine_api();

#undef EXTERN

#endif