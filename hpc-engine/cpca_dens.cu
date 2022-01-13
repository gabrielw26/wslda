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
                                         Complex *wf_d_dx, Complex *wf_d_dy, Complex *d_wf_laplace, double *kkz, 
                                         double *fbetaEn, double *d_weights,
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

#ifndef TAU_COMPUTATION_VIA_GRADIENTS
    double kz2;
#endif
    
    size_t iwf;
    double wght=1.0;
    
    if(ixyz<NXY)
    {
        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {    
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            if(d_weights!=NULL) wght=d_weights[iwf];
            fbEn*=wght; fbmEn*=wght; 
            kz = kkz[iwf];
#ifndef TAU_COMPUTATION_VIA_GRADIENTS
            kz2=kz*kz; // do not kill N/2 component in case of laplace
#endif
            wcnt = 2.0; // take into account +kz and -kz
            if(fabs(kz)<1.0e-12) wcnt = 1.0; // except for kz=0.0
#ifdef USE_CUBIC_CUTOFF
            if(fabs(kz+M_PI/DZ)<1.0e-12) {kz = 0.0; wcnt=1;}// momentum for which I should kill contribution for gradients
#else
            if(fabs(kz+M_PI/DZ)<1.0e-12) kz = 0.0; // momentum for which I should kill contribution for gradients
#endif            
            // read u and v
            u=wf[       iwf*NXY+ixyz];
            v=wf[n*NXY +iwf*NXY+ixyz];
            
            // form na, nb, nu
#ifdef SPINSYMMETRY_MODE
            // na+=... will be taken later as nb
            nb+=(thrust::norm(v)*fbmEn + thrust::norm(u)*fbEn)*wcnt;
            _nu+=u*thrust::conj(v)*(fbmEn-fbEn)*wcnt;
#else
            na+=thrust::norm(u)*fbEn *wcnt;
            nb+=thrust::norm(v)*fbmEn*wcnt;
            _nu+=u*thrust::conj(v)*(fbmEn-fbEn)/2.0 *wcnt;
#endif
            

            // read derivatives from u 
            wfdx=wf_d_dx[       iwf*NXY+ixyz];
            wfdy=wf_d_dy[       iwf*NXY+ixyz];
            wfdz=Complex(0.0,kz)*u; // i*kz*u(x,y)
            
            // form taua and j_a
#ifdef SPINSYMMETRY_MODE
            // do not compute taua and j_a - will taken from taub and j_b
            // compute only contributions to taub and j_b comming from derivatives of u
    #ifdef TAU_COMPUTATION_VIA_GRADIENTS
            taub+=(thrust::norm(wfdx)+thrust::norm(wfdy)+thrust::norm(wfdz))*fbEn*wcnt;
    #else
            taub+=(thrust::conj(u)*( d_wf_laplace[       iwf*NXY+ixyz] - u*kz2 )).real()*fbEn*wcnt; 
    #endif
            jbx+=(thrust::conj(u)*wfdx).imag()*fbEn*wcnt;
            jby+=(thrust::conj(u)*wfdy).imag()*fbEn*wcnt;
            // jbz+=(thrust::conj(u)*wfdz).imag()*fbEn*wcnt; // no currents along z direction
#else
          
    #ifdef TAU_COMPUTATION_VIA_GRADIENTS
            taua+=(thrust::norm(wfdx)+thrust::norm(wfdy)+thrust::norm(wfdz))*fbEn*wcnt;
    #else
            taua+=(thrust::conj(u)*( d_wf_laplace[       iwf*NXY+ixyz] - u*kz2 )).real()*fbEn*wcnt; 
    #endif
            jax+=(thrust::conj(u)*wfdx).imag()*fbEn*wcnt;
            jay+=(thrust::conj(u)*wfdy).imag()*fbEn*wcnt;
            // jaz+=(thrust::conj(u)*wfdz).imag()*fbEn*wcnt; // no currents along z direction
#endif
            
            // read derivatives from v 
            wfdx=wf_d_dx[n*NXY+iwf*NXY+ixyz];
            wfdy=wf_d_dy[n*NXY+iwf*NXY+ixyz];
            wfdz=Complex(0.0,kz)*v; // i*kz*v(x,y)
            
            // form taua and j_a
#ifdef TAU_COMPUTATION_VIA_GRADIENTS
            taub+=(thrust::norm(wfdx)+thrust::norm(wfdy)+thrust::norm(wfdz))*fbmEn*wcnt;
#else
            taub+=(thrust::conj(v)*( d_wf_laplace[n*NXY+iwf*NXY+ixyz] - v*kz2 )).real()*fbmEn*wcnt; 
#endif
            jbx-=(thrust::conj(v)*wfdx).imag()*fbmEn*wcnt;
            jby-=(thrust::conj(v)*wfdy).imag()*fbmEn*wcnt;
            // jbz-=(thrust::conj(v)*wfdz).imag()*fbmEn*wcnt; // no currents along z direction
        }
        
#ifdef SPINSYMMETRY_MODE
        na=nb;
        taua=taub;
        jax=jbx;
        jay=jby;
        jaz=jbz; 
#endif
        
        // save result to global memory and add missing LZ factor from 1/sqrt(LZ) * exp(i*kz*z)
        rho_a[ixyz]=na/DENS_FACTOR_M/(double)LZ;
        rho_b[ixyz]=nb/DENS_FACTOR_M/(double)LZ;
        tau_a[ixyz]=taua/DENS_FACTOR_M/(double)LZ;
        tau_b[ixyz]=taub/DENS_FACTOR_M/(double)LZ;
        nu[ixyz]=_nu/DENS_FACTOR_M/(double)LZ;
        j_a_x[ixyz]=jax/DENS_FACTOR_M/(double)LZ;
        j_a_y[ixyz]=jay/DENS_FACTOR_M/(double)LZ;
        j_a_z[ixyz]=jaz/DENS_FACTOR_M/(double)LZ;
        j_b_x[ixyz]=jbx/DENS_FACTOR_M/(double)LZ;
        j_b_y[ixyz]=jby/DENS_FACTOR_M/(double)LZ;
        j_b_z[ixyz]=jbz/DENS_FACTOR_M/(double)LZ;
    }
}

__global__ void kernel_calculate_densities_limited(size_t n, Complex *wf, double *kkz,
                                         double *fbetaEn, double *d_weights,
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
    
    double wght=1.0;
    
    if(ixyz<NXY)
    {
        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {    
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            if(d_weights!=NULL) wght=d_weights[iwf];
            fbEn*=wght; fbmEn*=wght; 
            kz = kkz[iwf];
            wcnt = 2.0; // take into account +kz and -kz
            if(fabs(kz)<1.0e-12) wcnt = 1.0; // except for kz=0.0
#ifdef USE_CUBIC_CUTOFF
            if(fabs(kz+M_PI/DZ)<1.0e-12) wcnt = 1.0;// momentum for which I should kill contribution for gradients
#endif            
            // read u and v
            u=wf[       iwf*NXY+ixyz];
            v=wf[n*NXY+iwf*NXY+ixyz];
                        
            // form na, nb, nu
#ifdef SPINSYMMETRY_MODE
            // na+=... will be taken later as nb
            nb+=(thrust::norm(v)*fbmEn + thrust::norm(u)*fbEn)*wcnt;
            _nu+=u*thrust::conj(v)*(fbmEn-fbEn)*wcnt;
#else
            na+=thrust::norm(u)*fbEn *wcnt;
            nb+=thrust::norm(v)*fbmEn*wcnt;
            _nu+=u*thrust::conj(v)*(fbmEn-fbEn)/2.0 *wcnt;
#endif
        }
        
#ifdef SPINSYMMETRY_MODE
        na=nb;
#endif
        
        // save result to global memory
        rho_a[ixyz]=na/DENS_FACTOR_M/(double)LZ;
        rho_b[ixyz]=nb/DENS_FACTOR_M/(double)LZ;
        tau_a[ixyz]=0.0;
        tau_b[ixyz]=0.0;
        nu[ixyz]=_nu/DENS_FACTOR_M/(double)LZ;
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
 * @param weights  extra weights, for computing subset densities, NULL - no weights (INPUT)
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
                            cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, 
                            cufftDoubleComplex *d_wf_laplace, double *kkz, 
                            double *d_fbetaEn, 
                            double *weights,
                            double *d_densities,
                            int gradients_computed, 
                            int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXY/nthreads);
    
    // Set pointers for to simplify notation
    // densities 
    Complex *nu   =(Complex *)(d_densities +  0*NXY);
    double *rho_a = (double *)(d_densities +  2*NXY);
    double *tau_a = (double *)(d_densities +  3*NXY);
    double *j_a_x = (double *)(d_densities +  4*NXY);
    double *j_a_y = (double *)(d_densities +  5*NXY);
    double *j_a_z = (double *)(d_densities +  6*NXY);
    double *rho_b = (double *)(d_densities +  7*NXY);
    double *tau_b = (double *)(d_densities +  8*NXY);
    double *j_b_x = (double *)(d_densities +  9*NXY);
    double *j_b_y = (double *)(d_densities + 10*NXY);
    double *j_b_z = (double *)(d_densities + 11*NXY);
    
    
    if(gradients_computed) // computation of all densities
    {
        kernel_calculate_densities<<<nblocks, nthreads>>>(n, (Complex *)wf, 
                                            (Complex *)wf_d_dx, (Complex *)wf_d_dy, (Complex *)d_wf_laplace, kkz, 
                                            d_fbetaEn, weights,
                                            rho_a, rho_b,
                                            tau_a, tau_b,
                                            nu,
                                            j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z);
    }
    else // only normal and anomalus density will be computed, other are set to zero
    {
        kernel_calculate_densities_limited<<<nblocks, nthreads>>>(n, (Complex *)wf, kkz,
                                            d_fbetaEn, weights,
                                            rho_a, rho_b,
                                            tau_a, tau_b,
                                            nu,
                                            j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z);        
    }

    return 0;
}


// ----------------------------------------------------------------------------------------
// ----------------------------------density_caculate_tau ---------------------------------
// ----------------------------------------------------------------------------------------
extern "C" void *pca_cufft_work_area;
extern "C" int compute_laplace_real_f(double *f, double *laplace_f, int nthreads);

__global__ void kernel_density_caculate_tau(double *laplace_rho, double *tau)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double _tau;
    
    if(ixyz<NXY)
    {
        _tau = 0.5*laplace_rho[ixyz] - tau[ixyz];
        
        // tau cannot be negative
        if(_tau<0.0) _tau=0.0;
        
        // save value
        tau[ixyz]=_tau;
    }
}

/**
 * Function computes tau according formula:
 * tau = (1/2)laplace(rho) - Re( sum_n Psi_n^* laplapce (Psi_n) )
 * See Eq.(28) in PHYSICAL REVIEW C 95, 044302 (2017)
 * @param d_densities array with densities (INPUT: `tau` array keeps only  Re( sum_n Psi_n^* laplapce (Psi_n) ), OUTPUT: kinetic energy density)
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR 
 */
extern "C" int density_caculate_tau(double *d_densities, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXY/nthreads);
    int ierr;
    
    // Set pointers for to simplify notation
    // densities 
//     Complex *nu   =(Complex *)(d_densities +  0*NXY);
    double *rho_a = (double *)(d_densities +  2*NXY);
    double *tau_a = (double *)(d_densities +  3*NXY);
//     double *j_a_x = (double *)(d_densities +  4*NXY);
//     double *j_a_y = (double *)(d_densities +  5*NXY);
//     double *j_a_z = (double *)(d_densities +  6*NXY);
    double *rho_b = (double *)(d_densities +  7*NXY);
    double *tau_b = (double *)(d_densities +  8*NXY);
//     double *j_b_x = (double *)(d_densities +  9*NXY);
//     double *j_b_y = (double *)(d_densities + 10*NXY);
//     double *j_b_z = (double *)(d_densities + 11*NXY);
    
    // I can use pca_cufft_work_area as working buffer for computation 
    double *laplace_rho=(double *)pca_cufft_work_area;
    
    // ------------- for b component -------------
    ierr = compute_laplace_real_f(rho_b, laplace_rho, nthreads);
    if(ierr!=0) return ierr - 101;
    kernel_density_caculate_tau<<<nblocks, nthreads>>>(laplace_rho,tau_b);
    
    // ------------- for a component ------------- 
#ifdef SPINSYMMETRY_MODE
    if( cudaMemcpy( tau_a , tau_b, sizeof(double)*NXY, cudaMemcpyDeviceToDevice )!= cudaSuccess ) return 100;
#else
    ierr = compute_laplace_real_f(rho_a, laplace_rho, nthreads);
    if(ierr!=0) return ierr - 102;
    kernel_density_caculate_tau<<<nblocks, nthreads>>>(laplace_rho,tau_a);
#endif
    
    return 0;
}

/**
 * Function sets the same values for densities "b" and for densities "a"
 * */
extern "C" int symmetrize_densities_device(double *d_densities)
{    
    // Set pointers for to simplify notation
    // densities 
//     Complex *nu   =(Complex *)(d_densities +  0*NXY);
    double *rho_a = (double *)(d_densities +  2*NXY);
//     double *tau_a = (double *)(d_densities +  3*NXY);
//     double *j_a_x = (double *)(d_densities +  4*NXY);
//     double *j_a_y = (double *)(d_densities +  5*NXY);
//     double *j_a_z = (double *)(d_densities +  6*NXY);
    double *rho_b = (double *)(d_densities +  7*NXY);
//     double *tau_b = (double *)(d_densities +  8*NXY);
//     double *j_b_x = (double *)(d_densities +  9*NXY);
//     double *j_b_y = (double *)(d_densities + 10*NXY);
//     double *j_b_z = (double *)(d_densities + 11*NXY);
    if( cudaMemcpy( rho_b , rho_a, sizeof(double)*NXY*5, cudaMemcpyDeviceToDevice )!= cudaSuccess ) return 100;
    return 0;
}
