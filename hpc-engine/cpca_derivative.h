// File list functions included into pca_derivative.cu

#ifndef __PCA_DERIVATIVE__
#define __PCA_DERIVATIVE__

#ifdef TDWSLDA_MAIN
#define cufftDoubleComplex double complex
#else
#include <cufft.h>
#endif

extern void *pca_cufft_work_area;
int create_cufftPlans(int batch_size,  int nwfip, size_t *workSize);
int set_workspace_for_cufftPlan(void *workArea);
int compute_derivatives(int n, cufftDoubleComplex *wf, cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *wf_d_dz, cufftDoubleComplex *wf_laplace, int nthreads);
int compute_gradient_real_f(double *f, double *df_dx, double *df_dy, double *df_dz, int nthreads);
int compute_derivative_real_vector_f(double *fx, double *fy, double *fz, double *dfx_dx, double *dfy_dy, double *dfz_dz,int nthreads);
int compute_laplace_real_f(double *f, double *laplace_f, int nthreads);
int compute_laplace(int n, cufftDoubleComplex *wf, cufftDoubleComplex *wf_laplace, int nthreads);
int compute_divergence_real_vector_f(double *fx, double *fy, double *fz, double *divf, int nthreads);

#ifdef TDWSLDA_MAIN
#undef cufftDoubleComplex
#endif

#endif
