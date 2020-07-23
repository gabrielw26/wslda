int calculate_densities(int n, cufftDoubleComplex *wf,
                            cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *d_wf_laplace, double *kkz, 
                            double *d_fbetaEn, 
                            double *d_densities,  
                            int gradients_computed, int nthreads);

int calculate_densities_weighted(int n, cufftDoubleComplex *wf,
                            double *kkz, 
                            double *d_fbetaEn, 
                            double *d_weights, 
                            double *d_densities,
                            int nthreads);
int density_caculate_tau(double *d_densities, int nthreads);
int symmetrize_densities_device(double *d_densities);

