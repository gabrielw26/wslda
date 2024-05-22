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
                                         Complex *wf_d_dx, Complex *d_wf_laplace, double *kky, double *kkz, 
                                         double *fbetaEn, double *d_weights,
                                         double *rho_a, double *rho_b,
                                         double *tau_a, double *tau_b,
                                         Complex *nu,
                                         double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                         int *cnt
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    double na=0.0, nb=0.0;
    double taua=0.0, taub=0.0;
    Complex _nu=Complex(0.0, 0.0);
    double jax=0.0, jay=0.0, jaz=0.0;
    double jbx=0.0, jby=0.0, jbz=0.0;

    Complex u, v, wfdx, wfdy, wfdz;
    double fbEn, fbmEn, ky, kz, wcnt;
    #define DENS_FACTOR_M 10000.

#ifndef TAU_COMPUTATION_VIA_GRADIENTS
    double kyz2;
#endif
    
    size_t iwf;
    
    double wght=1.0;
    
    if(ixyz<NX)
    {
        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {    
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            if(d_weights!=NULL) wght=d_weights[iwf];
            fbEn*=wght; fbmEn*=wght; 
            ky = kky[iwf];
            kz = kkz[iwf];
#ifndef TAU_COMPUTATION_VIA_GRADIENTS
            kyz2=kz*kz + ky*ky; // do not kill N/2 component in case of laplace
#endif

            wcnt = cnt[iwf];  
            if(fabs(ky+M_PI/DY)<1.0e-12) ky = 0.0; // momentum for which I should kill contribution for gradients
            if(fabs(kz+M_PI/DZ)<1.0e-12) kz = 0.0; // momentum for which I should kill contribution for gradients
            
            // read u and v
            u=wf[      iwf*NX+ixyz];
            v=wf[n*NX +iwf*NX+ixyz];
            
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
            wfdx=wf_d_dx[       iwf*NX+ixyz];
            wfdy=Complex(0.0,ky)*u; // i*ky*u(x)
            wfdz=Complex(0.0,kz)*u; // i*kz*u(x)
            
            // form taua and j_a
#ifdef SPINSYMMETRY_MODE
            // do not compute taua and j_a - will taken from taub and j_b
            // compute only contributions to taub and j_b comming from derivatives of u
    #ifdef TAU_COMPUTATION_VIA_GRADIENTS
            taub+=(thrust::norm(wfdx)+thrust::norm(wfdy)+thrust::norm(wfdz))*fbEn*wcnt;
    #else
            taub+=(thrust::conj(u)*( d_wf_laplace[       iwf*NX+ixyz] - u*kyz2 )).real()*fbEn*wcnt; 
    #endif
            jbx+=(thrust::conj(u)*wfdx).imag()*fbEn*wcnt;
            // jby+=(thrust::conj(u)*wfdy).imag()*fbEn*wcnt; // no currents along y direction
            // jbz+=(thrust::conj(u)*wfdz).imag()*fbEn*wcnt; // no currents along z direction
#else
          
    #ifdef TAU_COMPUTATION_VIA_GRADIENTS
            taua+=(thrust::norm(wfdx)+thrust::norm(wfdy)+thrust::norm(wfdz))*fbEn*wcnt;
    #else
            taua+=(thrust::conj(u)*( d_wf_laplace[       iwf*NX+ixyz] - u*kyz2 )).real()*fbEn*wcnt; 
    #endif
            jax+=(thrust::conj(u)*wfdx).imag()*fbEn*wcnt;
            // jay+=(thrust::conj(u)*wfdy).imag()*fbEn*wcnt; // no currents along y direction
            // jaz+=(thrust::conj(u)*wfdz).imag()*fbEn*wcnt; // no currents along z direction
#endif
            
            // read derivatives from v 
            wfdx=wf_d_dx[n*NX+iwf*NX+ixyz];
            wfdy=Complex(0.0,ky)*v; // i*ky*v(x)
            wfdz=Complex(0.0,kz)*v; // i*kz*v(x)
            
            // form taua and j_a
#ifdef TAU_COMPUTATION_VIA_GRADIENTS
            taub+=(thrust::norm(wfdx)+thrust::norm(wfdy)+thrust::norm(wfdz))*fbmEn*wcnt;
#else
            taub+=(thrust::conj(v)*( d_wf_laplace[n*NX+iwf*NX+ixyz] - v*kyz2 )).real()*fbmEn*wcnt; 
#endif
            jbx-=(thrust::conj(v)*wfdx).imag()*fbmEn*wcnt;
            // jby-=(thrust::conj(v)*wfdy).imag()*fbmEn*wcnt; // no currents along y direction
            // jbz-=(thrust::conj(v)*wfdz).imag()*fbmEn*wcnt; // no currents along z direction
        }
        
#ifdef SPINSYMMETRY_MODE
        na=nb;
        taua=taub;
        jax=jbx;
        jay=jby;
        jaz=jbz; 
#endif
        
        // save result to global memory and add missing (LY*LZ) factor from 1/sqrt((LY*LZ)) * exp(i*kz*z)
        rho_a[ixyz]=na/DENS_FACTOR_M/(double)(LY*LZ);
        rho_b[ixyz]=nb/DENS_FACTOR_M/(double)(LY*LZ);
        tau_a[ixyz]=taua/DENS_FACTOR_M/(double)(LY*LZ);
        tau_b[ixyz]=taub/DENS_FACTOR_M/(double)(LY*LZ);
        nu[ixyz]=_nu/DENS_FACTOR_M/(double)(LY*LZ);
        j_a_x[ixyz]=jax/DENS_FACTOR_M/(double)(LY*LZ);
        j_a_y[ixyz]=jay/DENS_FACTOR_M/(double)(LY*LZ);
        j_a_z[ixyz]=jaz/DENS_FACTOR_M/(double)(LY*LZ);
        j_b_x[ixyz]=jbx/DENS_FACTOR_M/(double)(LY*LZ);
        j_b_y[ixyz]=jby/DENS_FACTOR_M/(double)(LY*LZ);
        j_b_z[ixyz]=jbz/DENS_FACTOR_M/(double)(LY*LZ);
    }
}

__global__ void kernel_calculate_densities_limited(size_t n, Complex *wf, double *kky, double *kkz,
                                         double *fbetaEn, double *d_weights,
                                         double *rho_a, double *rho_b,
                                         double *tau_a, double *tau_b,
                                         Complex *nu,
                                         double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                         int *cnt
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    double na=0.0, nb=0.0;
    Complex _nu=Complex(0.0, 0.0);

    Complex u, v;
    double fbEn, fbmEn, /*ky, kz,*/ wcnt;
    #define DENS_FACTOR_M 10000.
    
    size_t iwf;
    
    double wght=1.0;
    
    if(ixyz<NX)
    {
        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {    
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            if(d_weights!=NULL) wght=d_weights[iwf];
            fbEn*=wght; fbmEn*=wght; 
//             ky = kky[iwf];
//             kz = kkz[iwf];
            wcnt = cnt[iwf];  
            
            // read u and v
            u=wf[     iwf*NX+ixyz];
            v=wf[n*NX+iwf*NX+ixyz];
                        
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
        rho_a[ixyz]=na/DENS_FACTOR_M/(double)(LY*LZ);
        rho_b[ixyz]=nb/DENS_FACTOR_M/(double)(LY*LZ);
        tau_a[ixyz]=0.0;
        tau_b[ixyz]=0.0;
        nu[ixyz]=_nu/DENS_FACTOR_M/(double)(LY*LZ);
        j_a_x[ixyz]=0.0;
        j_a_y[ixyz]=0.0;
        j_a_z[ixyz]=0.0;
        j_b_x[ixyz]=0.0;
        j_b_y[ixyz]=0.0;
        j_b_z[ixyz]=0.0;
    }
}
    
/**
 * Function computes densities.
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf array with wave-functions (INPUT)
 * @param wf_d_dx derivative with respect to dx (INPUT)
 * @param d_wf_laplace laplacian of wave-functions (INPUT)
 * @param kkyz value of ky and kz (INPUT)
 * @param fbetaEn weight of wave-function (INPUT)
 * @param weights extra weights, for computing subset densities, NULL - no weights (INPUT)
 * @param cnt degeneracy of states
 * @param d_densites (OUTPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]  
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NX,
 *                          nu - double complex array of size NX,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NX
 *                   In total size of d_densites is 12*NX
 * @param gradients_computed if 1 then gradients are computed, otherwise only normal and anomalus density will be computed
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int calculate_densities(int n, cufftDoubleComplex *wf,
                            cufftDoubleComplex *wf_d_dx,
                            cufftDoubleComplex *d_wf_laplace, double *kkyz,
                            double *d_fbetaEn, 
                            double *weights,
                            int *cnt,
                            double *d_densities,
                            int gradients_computed, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);
    
    // Set pointers for to simplify notation
    // densities 
    Complex *nu   =(Complex *)(d_densities +  0*NX);
    double *rho_a = (double *)(d_densities +  2*NX);
    double *tau_a = (double *)(d_densities +  3*NX);
    double *j_a_x = (double *)(d_densities +  4*NX);
    double *j_a_y = (double *)(d_densities +  5*NX);
    double *j_a_z = (double *)(d_densities +  6*NX);
    double *rho_b = (double *)(d_densities +  7*NX);
    double *tau_b = (double *)(d_densities +  8*NX);
    double *j_b_x = (double *)(d_densities +  9*NX);
    double *j_b_y = (double *)(d_densities + 10*NX);
    double *j_b_z = (double *)(d_densities + 11*NX);
    
    double *kky = kkyz;
    double *kkz = kkyz+n;
    
    if(gradients_computed) // computation of all densities
    {
        kernel_calculate_densities<<<nblocks, nthreads>>>(n, (Complex *)wf, 
                                            (Complex *)wf_d_dx, (Complex *)d_wf_laplace, kky, kkz, 
                                            d_fbetaEn, weights,
                                            rho_a, rho_b,
                                            tau_a, tau_b,
                                            nu,
                                            j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z, cnt);
    }
    else // only normal and anomalus density will be computed, other are set to zero
    {
        kernel_calculate_densities_limited<<<nblocks, nthreads>>>(n, (Complex *)wf, kky, kkz,
                                            d_fbetaEn, weights,
                                            rho_a, rho_b,
                                            tau_a, tau_b,
                                            nu,
                                            j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z, cnt);        
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
    
    if(ixyz<NX)
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
    int nblocks = (int)ceil((float)NX/nthreads);
    int ierr;
    
    // Set pointers for to simplify notation
    // densities 
//     Complex *nu   =(Complex *)(d_densities +  0*NX);
    double *rho_a = (double *)(d_densities +  2*NX);
    double *tau_a = (double *)(d_densities +  3*NX);
//     double *j_a_x = (double *)(d_densities +  4*NX);
//     double *j_a_y = (double *)(d_densities +  5*NX);
//     double *j_a_z = (double *)(d_densities +  6*NX);
    double *rho_b = (double *)(d_densities +  7*NX);
    double *tau_b = (double *)(d_densities +  8*NX);
//     double *j_b_x = (double *)(d_densities +  9*NX);
//     double *j_b_y = (double *)(d_densities + 10*NX);
//     double *j_b_z = (double *)(d_densities + 11*NX);
    
    // I can use pca_cufft_work_area as working buffer for computation 
    double *laplace_rho=(double *)pca_cufft_work_area;
    
    // ------------- for b component -------------
    ierr = compute_laplace_real_f(rho_b, laplace_rho, nthreads);
    if(ierr!=0) return ierr - 101;
    kernel_density_caculate_tau<<<nblocks, nthreads>>>(laplace_rho,tau_b);
    
    // ------------- for a component ------------- 
#ifdef SPINSYMMETRY_MODE
    if( cudaMemcpy( tau_a , tau_b, sizeof(double)*NX, cudaMemcpyDeviceToDevice )!= cudaSuccess ) return 100;
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
//     Complex *nu   =(Complex *)(d_densities +  0*NX);
    double *rho_a = (double *)(d_densities +  2*NX);
//     double *tau_a = (double *)(d_densities +  3*NX);
//     double *j_a_x = (double *)(d_densities +  4*NX);
//     double *j_a_y = (double *)(d_densities +  5*NX);
//     double *j_a_z = (double *)(d_densities +  6*NX);
    double *rho_b = (double *)(d_densities +  7*NX);
//     double *tau_b = (double *)(d_densities +  8*NX);
//     double *j_b_x = (double *)(d_densities +  9*NX);
//     double *j_b_y = (double *)(d_densities + 10*NX);
//     double *j_b_z = (double *)(d_densities + 11*NX);
    if( cudaMemcpy( rho_b , rho_a, sizeof(double)*NX*5, cudaMemcpyDeviceToDevice )!= cudaSuccess ) return 100;
    return 0;
}

// ================================================================================================
// ============================ calculate_quantum_friction_densities ==============================
// ================================================================================================
//---------EA VERSION
__global__ void kernel_calculate_quantum_friction_densities(size_t n, Complex *wf, Complex *d_wf_laplace, double *kky, double *kkz,
                                         double *fbetaEn, double *d_weights,
                                         double *d_qf_density_for_Ua,double *d_qf_density_for_Ub, Complex *d_qf_density_for_D,
                                         int *cnt
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex u, v, lap_v, lap_u;
    double fbEn, fbmEn;
    #define DENS_FACTOR_M 10000.

    size_t iwf;
    double /*kx, kz,*/ wcnt;

    double wght=1.0;

    if(ixyz<NX)
    {
        // Initialize variables for accumulating densities
        double  U_loc_a = 0.0;
        double  U_loc_b = 0.0;
        Complex D_loc = Complex(0.0, 0.0);

        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            if(d_weights!=NULL) wght=d_weights[iwf];
            fbEn*=wght; fbmEn*=wght; 

            // ky = kky[iwf];
            // kz = kkz[iwf];
            wcnt = cnt[iwf]; // degeneracy of the state

            // read u and v (from global memory)
            u=wf[     iwf*NX+ixyz];
            v=wf[n*NX+iwf*NX+ixyz];

            // read laplaces of u and v (from global memory)
            lap_u=d_wf_laplace[     iwf*NX+ixyz];
            lap_v=d_wf_laplace[n*NX+iwf*NX+ixyz];

            // TODO: EA: use kernel_calculate_densities for reference
            //  which implements computatoin of densities as presented
            //  https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Physical%20quantities#densities
            // ...
                U_loc_a -= (thrust::conj(lap_u)*u).imag()*fbEn;
                U_loc_b  = (thrust::conj(lap_v)*v).imag()*fbmEn;
                D_loc   -= lap_u*thrust::conj(v)+u*thrust::conj(lap_v);    //ojooooooooooooo D_loc doesn't include thermal factors
        }

    // send to global memory
        d_qf_density_for_Ua[ixyz] =  U_loc_a;   // TODO
        d_qf_density_for_Ub[ixyz] =  U_loc_b;   // TODO
        d_qf_density_for_D[ixyz] =   0.5*D_loc; // TODO
    }
}


// TODO: EA: Update of descriptions accordingly
/**
 * Function computes densities (generalzied) densities for quantum friction force.
 * U_a: -sum_n Im [u_n^* Laplace v_n] the cooling potencial for the current terms in species (a)
 * U_b:  sum_n Im [v_n^* Laplace v_n] the cooling potencial for the current terms in species (b)
 * D: sum_n [v_n^* Laplace u_n + u_n Laplace v_n^*] the cooling potential in pairing sector
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf array with wave-functions (INPUT)
 * @param d_wf_laplace laplacian of wave-functions (INPUT)
 * @param kkyz value of ky and kz (INPUT)
 * @param fbetaEn weight of wave-function (INPUT)
 * @param weights  extra weights, for computing subset densities, NULL - no weights (INPUT)
 * @param cnt degeneracy of states
 * @param d_qf_densities storage buffer for output densities (OUTPUT) +3 NXYZ
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 * */
extern "C" int calculate_quantum_friction_densities(int n, cufftDoubleComplex *wf, //---------EA VERSION
                            cufftDoubleComplex *d_wf_laplace,
                            double *kkyz,
                            double *d_fbetaEn,
                            double *weights,
                            int *cnt,
                            double *d_qf_densities,
                            int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);

    double *kky = kkyz;
    double *kkz = kkyz+n;

    // pointers algebra
    double * d_qf_density_for_Ua = (double *) d_qf_densities;        // density for diagonal part (U) of quantum friction force
    double * d_qf_density_for_Ub = (double *) d_qf_densities + NXYZ;        // density for diagonal part (U) of quantum friction force
    Complex *d_qf_density_for_D = (Complex *) d_qf_densities + 2*NXYZ; // density for off-diagonal part (Delta) of quantum friction force
    
    //kernel_calculate_quantum_friction_densities<<<nblocks, nthreads>>>(n, (Complex *)wf, (Complex *)d_wf_laplace, kky, kkz, d_fbetaEn, weights, d_qf_density_for_Ua, d_qf_density_for_Ub, d_qf_density_for_D, cnt);

    return 0;
}
