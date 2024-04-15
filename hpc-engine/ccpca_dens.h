#ifdef TDWSLDA_MAIN
#define cufftDoubleComplex double complex
#endif

int calculate_densities(int n, cufftDoubleComplex *wf,
                            cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *d_wf_laplace, double *kkyz,
                            double *d_fbetaEn, 
                            double *weights,
                            int *cnt,
                            double *d_densities,  
                            int gradients_computed, int nthreads);
int density_caculate_tau(double *d_densities, int nthreads);
int symmetrize_densities_device(double *d_densities);
int calculate_quantum_friction_densities(int n, Complex *wf,
                            Complex *d_wf_laplace,
                            double *kkyz,
                            double *d_fbetaEn,
                            double *weights,
                            int *cnt,
                            double *d_qf_densities,
                            int nthreads);

#ifdef TDWSLDA_MAIN
#undef cufftDoubleComplex
#endif
