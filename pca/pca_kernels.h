// List of functions from pca_kernels.cu

#ifndef __PCA_KERNELS__
#define __PCA_KERNELS__

// EDF functions
double polarization_h(double n_a, double n_b);
double der_polarization__der_na_h(double n_a, double n_b);
double der_polarization__der_nb_h(double n_a, double n_b);
double alpha_h(double p);
double alpha_a_h(double p);
double alpha_b_h(double p);
double alpha_plus_h(double p);
double alpha_minus_h(double p);
double der_alpha_plus__der_na_h(double n_a, double n_b);
double der_alpha_plus__der_nb_h(double n_a, double n_b);
double der_alpha_minus__der_na_h(double n_a, double n_b);
double der_alpha_minus__der_nb_h(double n_a, double n_b);
double funG_h(double p);
double funD_h(double n_a, double n_b);
double der_funD__der_na_h(double n_a, double n_b);
double der_funD__der_nb_h(double n_a, double n_b);
double tildeC_h(double n_a, double n_b);
double der_tildeC__der_na_h(double n_a, double n_b);
double der_tildeC__der_nb_h(double n_a, double n_b);

// params handling
void process_params(double *params, double kF);
int memcopy_const_params(double *params);

#ifdef WORK_IN_ROTATING_FRAME
int memcopy_const_Omega(double Omega_a, double Omega_b);
int set_Omega(double *params, double kF, double time, double *Omega_a, double *Omega_b);
#endif
#ifdef BDG_MODE
int memcopy_const_BdG(double aBdG);
#endif

// Utility functions
int set_gpu(int device); 
/**
 * Allocates on gpu device memory of given size 
 * and sets pointer
 * @return 0 - OK, otherwise - ERROR
 * */
int gpu_malloc(size_t size, void ** pointer);
int gpu_free(void * pointer);

/**
 * Page locked variants for host allocation
 * */
int host_malloc_pl(size_t size, void ** pointer);
int host_free_pl(void * pointer);

/**
 * functions to copy data cpu<->gpu
 * */
int memcopy_host2gpu(void * host, void * gpu, size_t size);
int memcopy_gpu2host(void * gpu, void * host, size_t size);
int memcopy_gpu2gpu(void * gpusrc, void * gpudst, size_t size);
int memcopy_const(double mu_a, double mu_b, double ec, double t0, double dt, double kF);

// reduction
int local_reductionR(double *array, int size, double *partial_sums, int threads, int mode);

// physics
int compute_potentials(int it, double *d_densities, double *d_potentials, double cccoeff, int nthreads);
int compute_energy(int it, double *d_densities, double *d_potentials, double *d_workarea, int nthreads);
int apply_hamiltonian(int n, cufftDoubleComplex *wf_in, cufftDoubleComplex *wf_out, 
                            cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *wf_d_dz, cufftDoubleComplex *wf_laplace, cufftDoubleComplex *alphawf_laplace,
                            double *d_densities, double *d_potentials, double qfalpha, double *useqpe, double cccoeff, 
                            int nthreads, cudaStream_t* streams);
int compute_ovelap(int n, cufftDoubleComplex *wf1, cufftDoubleComplex *wf2, double *overlap_re, double *overlap_im, 
                              double *workarea, int nthreads);
int amb_step1(int n, cufftDoubleComplex *ykm1, 
                         cufftDoubleComplex *fkm1, cufftDoubleComplex *fkm2, cufftDoubleComplex *fkm3, 
                         int nthreads);
int amb_step4(int n, cufftDoubleComplex *ykm1_in, cufftDoubleComplex *ykm1_out, 
                         cufftDoubleComplex *fkm3, 
                         int nthreads);
int amb45_step1(int n, cufftDoubleComplex *ykm1, 
                         cufftDoubleComplex *fkm1, cufftDoubleComplex *fkm2, cufftDoubleComplex *fkm3, cufftDoubleComplex *fkm4, 
                         int nthreads);
int amb45_step4(int n, cufftDoubleComplex *ykm1_in, cufftDoubleComplex *ykm1_out, 
                         cufftDoubleComplex *fkm4, 
                         int nthreads);
int normalize_wf(int n, cufftDoubleComplex *wf, int nthreads, cudaStream_t* streams);
int multiply_wf_by_alpha(int n, cufftDoubleComplex *wf_in, cufftDoubleComplex *wf_out, double *d_densities, int nthreads);
int taylor_expansion_contribution(int it, double dt, int n, cufftDoubleComplex *wf_hpsi, 
                                             cufftDoubleComplex *wf_update, cufftDoubleComplex *wf_contr, int nthreads);

int creat_streams (cudaStream_t* streams);
#endif

