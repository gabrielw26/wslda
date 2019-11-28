#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/complex.h>
typedef thrust::complex<double> Complex;
#include "pca_settings.h"



// ================================================================================================
// ========================================= calculate_densities ==================================
// ================================================================================================
__global__ void kernel_calculate_densities(size_t n, Complex *wf, 
                                         Complex *wf_d_dx, Complex *wf_d_dy, double *kkz, 
                                         double *fbetaEn,
                                         double *rho_a, double *rho_b,
                                         double *tau_a, double *tau_b,
                                         Complex *nu,
                                         double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    double na=0.0, nb=0.0;
    double taua=0.0, taub=0.0;
    Complex _nu=Complex(0.0, 0.0);
    double jax=0.0, jay=0.0, jaz=0.0;
    double jbx=0.0, jby=0.0, jbz=0.0;

    Complex u, v, wfdx, wfdy, wfdz;
    double fbEn, fbmEn, kz, wcnt;
    #define DENS_FACTOR_M 10000.
    
    size_t iwf;
    
    if(ixyz<NXY)
    {
        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {    
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            kz = kkz[iwf];
            wcnt = 2.0; // take into account +kz and -kz
            if(fabs(kz)<1.0e-12) wcnt = 1.0; // except for kz=0.0
            if(fabs(kz+M_PI)<1.0e-12) kz = 0.0; // momentum for which I should kill contribution for gradients
            
            // read u and v
            u=wf[       iwf*NXY+ixyz];
            v=wf[n*NXY +iwf*NXY+ixyz];
            
            // form na, nb, nu
            na+=thrust::norm(u)*fbEn *wcnt;
            nb+=thrust::norm(v)*fbmEn*wcnt;
            _nu+=u*thrust::conj(v)*(fbmEn-fbEn)/2.0 *wcnt;
            
            // read derivatives from u 
            wfdx=wf_d_dx[       iwf*NXY+ixyz];
            wfdy=wf_d_dy[       iwf*NXY+ixyz];
            wfdz=Complex(0.0,kz)*u; // i*kz*u(x,y)
            
            // form taua and j_a
            taua+=(thrust::norm(wfdx)+thrust::norm(wfdy)+thrust::norm(wfdz))*fbEn*wcnt;
            jax+=(thrust::conj(u)*wfdx).imag()*fbEn*wcnt;
            jay+=(thrust::conj(u)*wfdy).imag()*fbEn*wcnt;
            // jaz+=(thrust::conj(u)*wfdz).imag()*fbEn*wcnt; // no currents along z direction
            
            // read derivatives from v 
            wfdx=wf_d_dx[n*NXY+iwf*NXY+ixyz];
            wfdy=wf_d_dy[n*NXY+iwf*NXY+ixyz];
            wfdz=Complex(0.0,kz)*v; // i*kz*v(x,y)
            
            // form taua and j_a
            taub+=(thrust::norm(wfdx)+thrust::norm(wfdy)+thrust::norm(wfdz))*fbmEn*wcnt;
            jbx-=(thrust::conj(v)*wfdx).imag()*fbmEn*wcnt;
            jby-=(thrust::conj(v)*wfdy).imag()*fbmEn*wcnt;
            // jbz-=(thrust::conj(v)*wfdz).imag()*fbmEn*wcnt; // no currents along z direction
        }
        
        // save result to global memory and add missing NZ factor from 1/sqrt(NZ) * exp(i*kz*z)
        rho_a[ixyz]=na/DENS_FACTOR_M/(double)NZ;
        rho_b[ixyz]=nb/DENS_FACTOR_M/(double)NZ;
        tau_a[ixyz]=taua/DENS_FACTOR_M/(double)NZ;
        tau_b[ixyz]=taub/DENS_FACTOR_M/(double)NZ;
        nu[ixyz]=_nu/DENS_FACTOR_M/(double)NZ;
        j_a_x[ixyz]=jax/DENS_FACTOR_M/(double)NZ;
        j_a_y[ixyz]=jay/DENS_FACTOR_M/(double)NZ;
        j_a_z[ixyz]=jaz/DENS_FACTOR_M/(double)NZ;
        j_b_x[ixyz]=jbx/DENS_FACTOR_M/(double)NZ;
        j_b_y[ixyz]=jby/DENS_FACTOR_M/(double)NZ;
        j_b_z[ixyz]=jbz/DENS_FACTOR_M/(double)NZ;
    }
}

__global__ void kernel_calculate_densities_limited(size_t n, Complex *wf, double *kkz,
                                         double *fbetaEn,
                                         double *rho_a, double *rho_b,
                                         double *tau_a, double *tau_b,
                                         Complex *nu,
                                         double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    double na=0.0, nb=0.0;
    Complex _nu=Complex(0.0, 0.0);

    Complex u, v;
    double fbEn, fbmEn,kz, wcnt;
    #define DENS_FACTOR_M 10000.
    
    size_t iwf;
    
    if(ixyz<NXY)
    {
        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {    
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            kz = kkz[iwf];
            wcnt = 2.0; // take into account +kz and -kz
            if(fabs(kz)<1.0e-12) wcnt = 1.0; // except for kz=0.0
            
            // read u and v
            u=wf[       iwf*NXY+ixyz];
            v=wf[n*NXY+iwf*NXY+ixyz];
            
            // form na, nb, nu
            na+=thrust::norm(u)*fbEn *wcnt;
            nb+=thrust::norm(v)*fbmEn*wcnt;
            _nu+=u*thrust::conj(v)*(fbmEn-fbEn)/2.0 *wcnt;
        }
        
        // save result to global memory
        rho_a[ixyz]=na/DENS_FACTOR_M/(double)NZ;
        rho_b[ixyz]=nb/DENS_FACTOR_M/(double)NZ;
        tau_a[ixyz]=0.0;
        tau_b[ixyz]=0.0;
        nu[ixyz]=_nu/DENS_FACTOR_M/(double)NZ;
        j_a_x[ixyz]=0.0;
        j_a_y[ixyz]=0.0;
        j_a_z[ixyz]=0.0;
        j_b_x[ixyz]=0.0;
        j_b_y[ixyz]=0.0;
        j_b_z[ixyz]=0.0;
    }
}
    
/**
 * Function computes density.
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf array with wave-functions (INPUT)
 * @param wf_d_dx derivative with respect to dx (INPUT)
 * @param wf_d_dy derivative with respect to dy (INPUT)
 * @param kkz value of kz (INPUT)
 * @param fbetaEn weight of wave-function (INPUT)
 * @param d_densites (OUTPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]  
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NXY,
 *                          nu - double complex array of size NXY,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NXY
 *                   In total size of d_densites is 12*NXY
 * @param gradients_computed if 1 then gradients are computed, otherwise only normal and anomalus density will be computed
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int calculate_densities(int n, cufftDoubleComplex *wf,
                            cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, double *kkz, 
                            double *d_fbetaEn, 
                            double *d_densities,
                            int gradients_computed, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXY/nthreads);
    
    // Set pointers for to simplify notation
    // densities 
    double *rho_a = (double *)(d_densities +  0*NXY);
    double *rho_b = (double *)(d_densities +  1*NXY);
    double *tau_a = (double *)(d_densities +  2*NXY);
    double *tau_b = (double *)(d_densities +  3*NXY);
    Complex *nu   =(Complex *)(d_densities +  4*NXY);
    double *j_a_x = (double *)(d_densities +  6*NXY);
    double *j_a_y = (double *)(d_densities +  7*NXY);
    double *j_a_z = (double *)(d_densities +  8*NXY);
    double *j_b_x = (double *)(d_densities +  9*NXY);
    double *j_b_y = (double *)(d_densities + 10*NXY);
    double *j_b_z = (double *)(d_densities + 11*NXY);
    
    
    if(gradients_computed) // computation of all densities
    {
        kernel_calculate_densities<<<nblocks, nthreads>>>(n, (Complex *)wf, 
                                            (Complex *)wf_d_dx, (Complex *)wf_d_dy, kkz, 
                                            d_fbetaEn,
                                            rho_a, rho_b,
                                            tau_a, tau_b,
                                            nu,
                                            j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z);
    }
    else // only normal and anomalus density will be computed, other are set to zero
    {
        kernel_calculate_densities_limited<<<nblocks, nthreads>>>(n, (Complex *)wf, kkz,
                                            d_fbetaEn,
                                            rho_a, rho_b,
                                            tau_a, tau_b,
                                            nu,
                                            j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z);        
    }

    return 0;
}

// ----------------------------------------------------------------------------------------
// ----------------------------------calculate_densities_weighted -------------------------
// ----------------------------------------------------------------------------------------
__global__ void kernel_calculate_densities_weighted(size_t n, Complex *wf, double *kkz,
                                         double *fbetaEn, double *weights,
                                         double *rho_a, double *rho_b,
                                         double *tau_a, double *tau_b,
                                         Complex *nu,
                                         double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    double na=0.0, nb=0.0;

    Complex u, v;
    double fbEn, fbmEn,kz, wcnt;
    #define DENS_FACTOR_M 10000.
    
    size_t iwf;
    
    if(ixyz<NXY)
    {
        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {    
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            kz = kkz[iwf];
            wcnt = 2.0; // take into account +kz and -kz
            if(fabs(kz)<1.0e-12) wcnt = 1.0; // except for kz=0.0
            
            wcnt*=weights[iwf]; // add external weight 
            
            // read u and v
            u=wf[       iwf*NXY+ixyz];
            v=wf[n*NXY+iwf*NXY+ixyz];
            
            // form na, nb, nu
            na+=thrust::norm(u)*fbEn *wcnt;
            nb+=thrust::norm(v)*fbmEn*wcnt;
        }
        
        // save result to global memory
        rho_a[ixyz]=na/DENS_FACTOR_M/(double)NZ;
        rho_b[ixyz]=nb/DENS_FACTOR_M/(double)NZ;
        tau_a[ixyz]=0.0;
        tau_b[ixyz]=0.0;
        nu[ixyz]=0.0;
        j_a_x[ixyz]=0.0;
        j_a_y[ixyz]=0.0;
        j_a_z[ixyz]=0.0;
        j_b_x[ixyz]=0.0;
        j_b_y[ixyz]=0.0;
        j_b_z[ixyz]=0.0;
    }
}


/**
 * Function computes density.
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf array with wave-functions (INPUT)
 * @param kkz value of kz (INPUT)
 * @param fbetaEn weight of wave-function (INPUT)
 * @param weights for density computation (INPUT)
 * @param d_densites (OUTPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]  
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NXY,
 *                          nu - double complex array of size NXY,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NXY
 *                   In total size of d_densites is 12*NXY
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int calculate_densities_weighted(int n, cufftDoubleComplex *wf,
                            double *kkz, 
                            double *d_fbetaEn, 
                            double *d_weights, 
                            double *d_densities,
                            int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXY/nthreads);
    
    // Set pointers for to simplify notation
    // densities 
    double *rho_a = (double *)(d_densities +  0*NXY);
    double *rho_b = (double *)(d_densities +  1*NXY);
    double *tau_a = (double *)(d_densities +  2*NXY);
    double *tau_b = (double *)(d_densities +  3*NXY);
    Complex *nu   =(Complex *)(d_densities +  4*NXY);
    double *j_a_x = (double *)(d_densities +  6*NXY);
    double *j_a_y = (double *)(d_densities +  7*NXY);
    double *j_a_z = (double *)(d_densities +  8*NXY);
    double *j_b_x = (double *)(d_densities +  9*NXY);
    double *j_b_y = (double *)(d_densities + 10*NXY);
    double *j_b_z = (double *)(d_densities + 11*NXY);
    

    kernel_calculate_densities_weighted<<<nblocks, nthreads>>>(n, (Complex *)wf, kkz,
                                        d_fbetaEn, d_weights, 
                                        rho_a, rho_b,
                                        tau_a, tau_b,
                                        nu,
                                        j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z);        


    return 0;
}
