// List of functions from pca_kernels.cu

#ifndef __CPCA_KERNELS__
#define __CPCA_KERNELS__

#ifdef TDWSLDA_MAIN
#define cufftDoubleComplex double complex
#endif

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
// params handling
#ifdef TDWSLDA
void process_params(double *params, double *kF, double *mu, size_t extra_data_size, void *extra_data);
size_t get_extra_data_size(double *params);
int load_extra_data(size_t size, void *extra_data, double *params);
int memcopy_extra_data(size_t extra_data_size, void *extra_data);
#else
void process_params(double *params, double kF, double *mu);
#endif
int memcopy_const_params(double *params);
int memcopy_const_sclgth(double a);

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
int apply_hamiltonian(int it, int n, cufftDoubleComplex *wf_in, cufftDoubleComplex *wf_out, 
                            cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, double *d_kkz, cufftDoubleComplex *wf_laplace, cufftDoubleComplex *alphawf_laplace,
                            double *d_densities, double *d_potentials, double qfswitch, double *useqpe, double pccoeff,
                            int nthreads);
int compute_ovelap(int n, cufftDoubleComplex *wf1, cufftDoubleComplex *wf2, double *overlap_re, double *overlap_im, 
                              double *workarea, int nthreads);
int abm_step1(int n, cufftDoubleComplex *ykm1,
                         cufftDoubleComplex *fkm1,
                         cufftDoubleComplex *fkm2,
                         cufftDoubleComplex *fkm3,
                         cufftDoubleComplex *fkm4,
                         cufftDoubleComplex *fkm5,
                         int abm_scheme,
                         int nthreads);
int abm_step4(int n, cufftDoubleComplex *ykm1_in, cufftDoubleComplex *ykm1_out,
                         cufftDoubleComplex *fkm_last,
                         int abm_scheme,
                         int nthreads);
int normalize_wf(int n, cufftDoubleComplex *wf, int nthreads);
int multiply_wf_by_alpha(int n, cufftDoubleComplex *wf_in, cufftDoubleComplex *wf_out, double *d_densities, int nthreads);
int taylor_expansion_contribution(int it, double dt, int n, cufftDoubleComplex *wf_hpsi, 
                                             cufftDoubleComplex *wf_update, cufftDoubleComplex *wf_contr, int nthreads);

#ifdef TDWSLDA_MAIN
#undef cufftDoubleComplex
#endif

#ifdef ENABLE_MODIFY_DENSITIES
void modify_densities(int it, wslda_density h_densities, double *params, size_t extra_data_size, void *h_extra_data,
                              wslda_density d_densities,                                         void *d_extra_data);
#endif

#endif

