int calculate_densities(int n, cufftDoubleComplex *wf,
                            cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *wf_d_dz, 
                            cufftDoubleComplex *d_wf_laplace, 
                            double *d_fbetaEn, 
                            double *d_densities,  
                            int gradients_computed, int nthreads);

int calculate_densities_weighted(int n, cufftDoubleComplex *wf,
                            double *d_fbetaEn, 
                            double *d_weights, 
                            double *d_densities,
                            int nthreads);
int density_caculate_tau(double *d_densities, int nthreads);
