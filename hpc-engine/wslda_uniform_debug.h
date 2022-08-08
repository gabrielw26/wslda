/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 28-07-2021
 * 
 * This file provides simple functions for testings corctness of results in case of uniform-testcase.
 * */  


int debug_densities_uniform(int blocksize, double *h_densities )
{
    // Set pointers for to simplify notation
    // densities 
    double complex *nu   =(double complex *)(h_densities +  0*blocksize);
    double *rho_a = (double *)(h_densities +  2*blocksize);
    double *tau_a = (double *)(h_densities +  3*blocksize);
    double *j_a_x = (double *)(h_densities +  4*blocksize);
    double *j_a_y = (double *)(h_densities +  5*blocksize);
    double *j_a_z = (double *)(h_densities +  6*blocksize);
    double *rho_b = (double *)(h_densities +  7*blocksize);
    double *tau_b = (double *)(h_densities +  8*blocksize);
    double *j_b_x = (double *)(h_densities +  9*blocksize);
    double *j_b_y = (double *)(h_densities + 10*blocksize);
    double *j_b_z = (double *)(h_densities + 11*blocksize);
    
    int i;
    double epsilon=1.0e-6;
    int ierr=0; // OK
    
    for(i=0; i<blocksize; i++) 
    {
        if(fabs(rho_a[i]-__md_pca_uniform.n0_a)>epsilon) 
        {
            ierr++; 
            printf("DEBUG: ERROR: rho_a: %f != %f\n", rho_a[i],__md_pca_uniform.n0_a);
            break;
        }
    }
    
    for(i=0; i<blocksize; i++) 
    {
        if(fabs(fabs(tau_a[i])-__md_pca_uniform.tau_a)>epsilon) 
        {
            ierr++; 
            printf("DEBUG: ERROR: tau_a: %f != %f\n", tau_a[i],__md_pca_uniform.tau_a);
            break;
        }
    }
    
    for(i=0; i<blocksize; i++) 
    {
        if(fabs(j_a_x[i]-0.0)>epsilon) 
        {
            ierr++; 
            printf("DEBUG: ERROR: j_a_x: %f != %f\n", j_a_x[i],0.0);
            break;
        }
    }
    
    for(i=0; i<blocksize; i++) 
    {
        if(fabs(j_a_y[i]-0.0)>epsilon) 
        {
            ierr++; 
            printf("DEBUG: ERROR: j_a_y: %f != %f\n", j_a_y[i],0.0);
            break;
        }
    }
    
    for(i=0; i<blocksize; i++) 
    {
        if(fabs(j_a_z[i]-0.0)>epsilon) 
        {
            ierr++; 
            printf("DEBUG: ERROR: j_a_z: %f != %f\n", j_a_z[i],0.0);
            break;
        }
    }
    
    for(i=0; i<blocksize; i++) 
    {
        if(fabs(cabs(nu[i])-fabs(__md_pca_uniform.nu))>epsilon) 
        {
            ierr++; 
            printf("DEBUG: ERROR: nu: %f != %f\n", cabs(nu[i]),fabs(__md_pca_uniform.nu));
            break;
        }
    }
    return ierr; 
}
