// Author: Gabriel Wlazlowski
// Date: 09-09-2016

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/complex.h>
#include <stdio.h>
typedef thrust::complex<double> Complex;

#include "pca_settings.h"
#define CODEDIM 3
#include "pca_macro.h"
#include "pca_edf.h"
#ifdef HIPMODE
#include "hip/wslda_cuda_utils.hpp"
#else
#include "wslda_cuda_utils.h"
#endif
#define double_complex Complex
#define __externc
#include "wslda_potdens.h"
#include "reduce_many.h"
#include "pca_utils.h" // access to input (md) structure
#include "tdwslda_memory_management.h"


#ifdef TDWSLDA

#include "problem-definition.h"

#ifdef ENABLE_V_EXT
#define u_ext(ix, iy, iz, it, spin) v_ext(ix, iy, iz, it, spin, dc_params, dc_extra_data_size, dc_extra_data)
#else
#define u_ext(ix, iy, iz, it, spin) 0.0
#endif

#ifdef ENABLE_DELTA_EXT
#define macro_delta_ext(ix, iy, iz, it, delta) delta_ext(ix, iy, iz, it, delta, dc_params, dc_extra_data_size, dc_extra_data)
#else
#define macro_delta_ext(ix, iy, iz, it, delta) Complex(0.0,0.0)
#endif

#endif

// Functionals
#include "tdwslda_functionals.h"

// ode integrator
#ifdef HIPMODE
#include "hip/tdwslda_ode_integrator.hpp"
#else
#include "tdwslda_ode_integrator.h"
#endif


// =======================================================================================
// ================================ compute_potentials ===================================
// =======================================================================================
/**
 * Function computes potentials V_a, V_b and delta
 * using formulas from section "9.3.2.2 Summary"
 * @param it index of time step, it is ised for proper evaluation of external potential
 * @param d_densites (INPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NXYZ,
 *                          nu - double complex array of size NXYZ,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NXYZ
 *                   In total size of d_densites is 12*NXYZ
 * @param d_potentials (INPUT/OUTPUT)
 *                     collective array with potentials [V_a, V_b, delta]
 *                     where: V_a, V_b - double arrays of size NXYZ
 *                            delta - double complex array of size NXYZ
 *                     In total size of d_potentials is 4*NXYZ
 *                     NOTE: I assume that d_potentials contains potentials from previous iteration,
 *                           i.e. they are good starting point for self-consistent process.
 * @param cccoeff the current corrections coefficient
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 **/
extern "C" int compute_potentials(int it, double *d_densities, double *d_potentials, double cccoeff, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXYZ/nthreads);
    wslda_density densall=convert_into_wslda_density(d_densities, NXYZ);
    wslda_potential potsall=convert_into_wslda_potential(d_potentials, NXYZ, NULL);

    tdwslda_compute_potentials<<<nblocks, nthreads>>>(it, densall, potsall, cccoeff);

#ifdef ENABLE_MODIFY_POTENTIALS
    modify_potentials<<<nblocks, nthreads>>>(it, densall, potsall);
#endif
    return 0;
}

// =======================================================================================
// ====================================== zero_array =====================================
// =======================================================================================
__global__ void kernel_zero_array(int asize, double *array)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    if(ixyz<asize) array[ixyz]=0.0;
}

/**
 * Fills array with zeros
 * @param asize size of array
 * @param array pointer to array
 * */
int zero_array(int asize, double *array)
{
    // number of blocks
    int nblocks = (int)ceil((float)asize/512);
    kernel_zero_array<<<nblocks, 512>>>(asize,array);
    return 0;
}

// =======================================================================================
// ================================== compute_energy =====================================
// =======================================================================================
__global__ void kernel_compute_energy_v_ext(int it,
                                      double *rho_a, double *rho_b, Complex *nu,
                                      double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                      double *tau_a, double *tau_b,
                                      Complex *delta,
                                      double *E_ext
                                      )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;

    if(ixyz<NXYZ)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); // decode cartesian coordinates

        // External potential energy
        E_ext[ixyz]=(rho_a[ixyz]*u_ext(ix,iy,iz,it,SPINA) + rho_b[ixyz]*u_ext(ix,iy,iz,it,SPINB))*DXYZ;
    }
}

__global__ void kernel_compute_energy_delta_ext(int it,
                                      double *rho_a, double *rho_b, Complex *nu,
                                      double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                      double *tau_a, double *tau_b,
                                      Complex *delta,
                                      double *E_ext
                                      )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;

    double na, nb;

    if(ixyz<NXYZ)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); // decode cartesian coordinates

        // External potential energy
        E_ext[ixyz]=(thrust::conj(nu[ixyz])*macro_delta_ext(ix, iy, iz, it, delta[ixyz])).real()*(-2.0)*DXYZ;
    }
}

__global__ void kernel_compute_energy_velocity_ext(int it,
                                      double *rho_a, double *rho_b, Complex *nu,
                                      double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                      double *tau_a, double *tau_b,
                                      Complex *delta,
                                      double *E_ext
                                      )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;

    if(ixyz<NXYZ)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); // decode cartesian coordinates

        // External potential energy
        E_ext[ixyz]=  (
                         j_a_x[ixyz]*velocity_ext(ix, iy, iz, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data)
                        +j_a_y[ixyz]*velocity_ext(ix, iy, iz, it, SPINA, YAXIS, dc_params, dc_extra_data_size, dc_extra_data)
                        +j_a_z[ixyz]*velocity_ext(ix, iy, iz, it, SPINA, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data)
                        +j_b_x[ixyz]*velocity_ext(ix, iy, iz, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data)
                        +j_b_y[ixyz]*velocity_ext(ix, iy, iz, it, SPINB, YAXIS, dc_params, dc_extra_data_size, dc_extra_data)
                        +j_b_z[ixyz]*velocity_ext(ix, iy, iz, it, SPINB, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data)
                       )*(-1.0)*DXYZ;

    }
}

__global__ void kernel_compute_angular_momentum_z(double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                      double *Laz, double *Lbz)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;

    double _x, _y;

    if(ixyz<NXYZ)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); // decode cartesian coordinates

        _x = (double)(ix-NX/2)*DX;
        _y = (double)(iy-NY/2)*DY;

        Laz[ixyz] = (_x*j_a_y[ixyz] - _y*j_a_x[ixyz])*DXYZ;
        Lbz[ixyz] = (_x*j_b_y[ixyz] - _y*j_b_x[ixyz])*DXYZ;
    }
}

// The main body of function computing the energy is moved to separate file
#ifdef HIPMODE
#include "hip/tdwslda_energy_generic.hpp"
#else
#include "tdwslda_energy_generic.h"
#endif
extern "C" int compute_energy(int it, double *d_densities, double *d_potentials, double *d_workarea, int nthreads)
{
    return compute_energy_generic(it, d_densities, d_potentials, d_workarea, nthreads);
}

// =======================================================================================
// ================================== apply_hamiltonian ==================================
// =======================================================================================
extern "C" int compute_gradient_real_f(double *f, double *df_dx, double *df_dy, double *df_dz, int nthreads);
extern "C" int compute_derivative_real_vector_f(double *fx, double *fy, double *fz, double *dfx_dx, double *dfy_dy, double *dfz_dz,int nthreads);
extern "C" int compute_laplace_real_f(double *f, double *laplace_f, int nthreads);
extern "C" int compute_divergence_real_vector_f(double *fx, double *fy, double *fz, double *divf, int nthreads);
extern "C" int high_frequency_filter_d(double *in, double *out, double fd_mu, double fd_T, int nthreads);
extern "C" int high_frequency_filter_c(Complex *in, Complex *out, double fd_mu, double fd_T, int nthreads);
extern "C" int high_frequency_filter_massive_d(int n, double *in, double *out, double fd_mu, double fd_T, int nthreads);

__global__ void kernel_apply_hamiltonian(int it, wslda_potential h_potentials,
                                         double * laplace_alpha_a, double *laplace_alpha_b,
                                         double *j_corr_a_x, double *j_corr_a_y, double *j_corr_a_z, double *j_corr_b_x, double *j_corr_b_y, double *j_corr_b_z,
                                         size_t n, Complex *wf_in, Complex *wf_out,
                                         Complex *wf_d_dx, Complex *wf_d_dy, Complex *wf_d_dz, Complex *wf_laplace, Complex *alphawf_laplace,
                                         double *vx_a, double *vy_a, double *vz_a, double *divv_a, double *vx_b, double *vy_b, double *vz_b, double *divv_b
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double aa, ab;
    double na, nb;
    double Va, Vb;
    Complex D;
    double cja=0.0, cjb=0.0;
#ifdef CURRENT_CORRECTIONS
    Complex gax, gay, gaz, gbx, gby, gbz;
    double ja, jb/*, jp*/;
    double /*fr,*/ fra, frb; // regularization functions
#endif
#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
    int ix, iy, iz, i; // need to decode coordinate
#endif

    size_t iwf;
    Complex u, v, tu, tv;

    if(ixyz<NXYZ)
    {
#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); // decode cartesian coordinates
#endif

        // effective mass
        aa=h_potentials.alpha_a[ixyz]; // will be multiplied by -1/2 later
        ab=h_potentials.alpha_b[ixyz]; // will be multiplied by -1/2 later

        // read potentials
        Va=h_potentials.V_a[ixyz];
        Vb=h_potentials.V_b[ixyz];
        D =h_potentials.delta[ixyz];

#ifdef ENABLE_DELTA_EXT
        D = D + macro_delta_ext(ix, iy, iz, it, D);
#endif

#ifdef CURRENT_CORRECTIONS
        // reset
        gax=Complex(0.0, 0.0);
        gay=Complex(0.0, 0.0);
        gaz=Complex(0.0, 0.0);
        gbx=Complex(0.0, 0.0);
        gby=Complex(0.0, 0.0);
        gbz=Complex(0.0, 0.0);
#endif

        // read gradient corrections
#ifdef CURRENT_CORRECTIONS
        cja+=-0.5*(j_corr_a_x[ixyz]+j_corr_a_y[ixyz]+j_corr_a_z[ixyz]);
        cjb+=-0.5*(j_corr_b_x[ixyz]+j_corr_b_y[ixyz]+j_corr_b_z[ixyz]);

        gax+=Complex(0.0,  -1.*h_potentials.A_a_x[ixyz]);
        gay+=Complex(0.0,  -1.*h_potentials.A_a_y[ixyz]);
        gaz+=Complex(0.0,  -1.*h_potentials.A_a_z[ixyz]);
        gbx+=Complex(0.0,   1.*h_potentials.A_b_x[ixyz]); // note conjugate of complex number (beacuse of "-h*" operator)
        gby+=Complex(0.0,   1.*h_potentials.A_b_y[ixyz]); // note conjugate of complex number (beacuse of "-h*" operator)
        gbz+=Complex(0.0,   1.*h_potentials.A_b_z[ixyz]); // note conjugate of complex number (beacuse of "-h*" operator)
#endif

#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
        // multiply mass by -1/4
        aa*=-0.5;
        ab*=-0.5;
#else
        // multiply mass by -1/4
        aa*=-0.25;
        ab*=-0.25;

        // Read laplace of effective mass - I use na and nb as working buffers
        na=0.25*laplace_alpha_a[ixyz];
        nb=0.25*laplace_alpha_b[ixyz];
#endif

#ifdef ENABLE_VELOCITY_EXT
        gax+=Complex(0.0, 1.0*vx_a[ixyz]);
        gbx+=Complex(0.0,-1.0*vx_b[ixyz]); // NOTE: complex conjugate included

        gay+=Complex(0.0, 1.0*vy_a[ixyz]);
        gby+=Complex(0.0,-1.0*vy_b[ixyz]); // NOTE: complex conjugate included

        gaz+=Complex(0.0, 1.0*vz_a[ixyz]);
        gbz+=Complex(0.0,-1.0*vz_b[ixyz]); // NOTE: complex conjugate included

        cja+=  0.5*divv_a[ixyz];
        cjb+=  0.5*divv_b[ixyz];
#endif

        // apply to each wave-function
        for(iwf=0; iwf<n; iwf++)
        {
            // reset
            u=Complex(0.0, 0.0);
            v=Complex(0.0, 0.0);

            // read wf
            tu=wf_in[       iwf*NXYZ+ixyz];
            tv=wf_in[n*NXYZ+iwf*NXYZ+ixyz];

            u+=tu*Complex(Va-dc_mu_a,     cja); // V_a*u
            v-=tv*Complex(Vb-dc_mu_b,-1.0*cjb); // V_b*v, note conjugate of complex number (beacuse of "-h*" operator)

            u+=tv*D;                // delta   * v
            v+=tu*thrust::conj(D);  // delta^* * u

#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // kinetic part - contribution: (1/4)*laplace(alpha)*u, see further
            u+=tu*na; // note: na keeps (1/4)*laplace(alpha_a)
            v-=tv*nb; // note: nb keeps (1/4)*laplace(alpha_b)
#endif

#ifdef CURRENT_CORRECTIONS
            // read gradients of wf: x-coordinate
            tu=wf_d_dx[       iwf*NXYZ+ixyz];
            tv=wf_d_dx[n*NXYZ+iwf*NXYZ+ixyz];
            u+=tu*gax;
            v-=tv*gbx;

            // read gradients of wf: y-coordinate
            tu=wf_d_dy[       iwf*NXYZ+ixyz];
            tv=wf_d_dy[n*NXYZ+iwf*NXYZ+ixyz];
            u+=tu*gay;
            v-=tv*gby;

            // read gradients of wf: z-coordinate
            tu=wf_d_dz[       iwf*NXYZ+ixyz];
            tv=wf_d_dz[n*NXYZ+iwf*NXYZ+ixyz];
            u+=tu*gaz;
            v-=tv*gbz;
#endif


            // kinetic part
            // note: (-1/2)\nabla(alpha * nabla u) = -(1/4)*laplace(alpha*u) - (1/4)*alpha*laplace(u) + (1/4)*laplace(alpha)*u
            // (1/4)*laplace(alpha)*u already done
            // read laplace of wf
            tu=wf_laplace[       iwf*NXYZ+ixyz];
            tv=wf_laplace[n*NXYZ+iwf*NXYZ+ixyz];
            u+=tu*aa; // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_a is given by aa
            v-=tv*ab; // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_b is given by ab
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // read laplace(alpha*u)
            u+=alphawf_laplace[       iwf*NXYZ+ixyz]*(-0.25);
            v-=alphawf_laplace[n*NXYZ+iwf*NXYZ+ixyz]*(-0.25);
#endif
            // save to global memory wf
            wf_out[       iwf*NXYZ+ixyz]=u;
            wf_out[n*NXYZ+iwf*NXYZ+ixyz]=v;
        }
    }
}

__global__ void kernel_apply_hamiltonian_bdg(int it,
                                         double *V_a, double *V_b, Complex *delta,
                                         size_t n, Complex *wf_in, Complex *wf_out,
                                         Complex *wf_d_dx, Complex *wf_d_dy, Complex *wf_d_dz, Complex *wf_laplace,
                                         double *vx_a, double *vy_a, double *vz_a, double *divv_a, double *vx_b, double *vy_b, double *vz_b, double *divv_b
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x *  blockDim.x; // compute for this point

    double Va, Vb;
    Complex D;

#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
    int ix, iy, iz, i; // need to decode coordinate
#endif

    double cja=0.0, cjb=0.0;
#ifdef CURRENT_CORRECTIONS
    Complex gax, gay, gaz;
    Complex gbx, gby, gbz;
#endif

    size_t iwf;
    Complex u, v, tu, tv;

    if(ixyz<NXYZ)
    {
#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); // decode cartesian coordinates
#endif
        // read potentials
        Va=V_a[ixyz];
        Vb=V_b[ixyz];
        D=delta[ixyz];

#ifdef ENABLE_DELTA_EXT
        D = D + macro_delta_ext(ix, iy, iz, it, D);
#endif

#ifdef ENABLE_VELOCITY_EXT
        gax=Complex(0.0, 1.0*vx_a[ixyz]);
        gbx=Complex(0.0,-1.0*vx_b[ixyz]); // NOTE: complex conjugate included

        gay=Complex(0.0, 1.0*vy_a[ixyz]);
        gby=Complex(0.0,-1.0*vy_b[ixyz]); // NOTE: complex conjugate included

        gaz=Complex(0.0, 1.0*vz_a[ixyz]);
        gbz=Complex(0.0,-1.0*vz_b[ixyz]); // NOTE: complex conjugate included

        cja=  0.5*divv_a[ixyz];
        cjb=  0.5*divv_b[ixyz];
#endif

        // apply to each wave-function
        for(iwf=0; iwf<n; iwf++)
        {
            // reset
            u=Complex(0.0, 0.0);
            v=Complex(0.0, 0.0);

            // read wf
            tu=wf_in[        iwf*NXYZ+ixyz];
            tv=wf_in[n*NXYZ +iwf*NXYZ+ixyz];

            u+=tu*Complex(Va-dc_mu_a,      cja); // V_a*u
            v-=tv*Complex(Vb-dc_mu_b, -1.0*cjb); // V_b*v, note conjugate of complex number (beacuse of "-h*" operator)

            u+=tv*D;                // delta   * v
            v+=tu*thrust::conj(D);  // delta^* * u

#ifdef CURRENT_CORRECTIONS
            // read gradients of wf: x-coordinate
            tu=wf_d_dx[       iwf*NXYZ+ixyz];
            tv=wf_d_dx[n*NXYZ+iwf*NXYZ+ixyz];
            u+=tu*gax;
            v-=tv*gbx;

            // read gradients of wf: y-coordinate
            tu=wf_d_dy[       iwf*NXYZ+ixyz];
            tv=wf_d_dy[n*NXYZ+iwf*NXYZ+ixyz];
            u+=tu*gay;
            v-=tv*gby;

            // read gradients of wf: z-coordinate
            tu=wf_d_dz[       iwf*NXYZ+ixyz];
            tv=wf_d_dz[n*NXYZ+iwf*NXYZ+ixyz];
            u+=tu*gaz;
            v-=tv*gbz;
#endif
            // kinetic part

            // read laplace of wf
            tu=wf_laplace[        iwf*NXYZ+ixyz];
            tv=wf_laplace[n*NXYZ +iwf*NXYZ+ixyz];
            u+=tu*(-0.5); // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_a is given by aa
            v-=tv*(-0.5); // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_b is given by ab

            // save to global memory wf
            wf_out[        iwf*NXYZ+ixyz]=u;
            wf_out[n*NXYZ +iwf*NXYZ+ixyz]=v;
        }
    }
}

 __global__ void kernel_add_quantum_friction(double *rho_a, double *rho_b, //add difference of particle number and desired particle number//add difference of particle number and desired particle number
                                             double* d_qf_density_for_Ua, double* d_qf_density_for_Ub, Complex*d_qf_density_for_D, 
                                             double *V_a, double *V_b, Complex *delta, double qfalpha, double qfbeta,double qfgamma)
{

    double PhDel,PhDen;
    double _V_a, _V_b;
    Complex _delta;
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    if(ixyz<NXYZ)
    {
        // see Eq.(3) in paper https://arxiv.org/abs/1305.6891
        PhDel         = thrust::arg(delta[ixyz]);
        PhDen         = thrust::arg(d_qf_density_for_D[ixyz]);

        _V_a = V_a[ixyz]; // keep the original value
        V_a[ixyz]   -=d_qf_density_for_Ua[ixyz]*qfalpha/dc_nF; // here I divide be reference density, to avoid problems of division by zero
        _V_b = V_b[ixyz]; // keep the original value
        V_b[ixyz]   -=d_qf_density_for_Ub[ixyz]*qfalpha/dc_nF; // here I divide be reference density, to avoid problems of division by zero
        _delta = delta[ixyz]; // keep the original value
        delta[ixyz] += qfbeta*thrust::abs(delta[ixyz])*Complex(cos(PhDel),sin(PhDel))*sin(PhDel - PhDen); //particle conserving part
        // particle control part
        delta[ixyz] += Complex(0.0,qfgamma)*_delta;

        // Save original values - it will be used by kernel_remove_quantum_friction(...)
        d_qf_density_for_Ua[ixyz]=_V_a;
        d_qf_density_for_Ub[ixyz]=_V_b;
        d_qf_density_for_D[ixyz]=_delta;
    }
}

__global__ void kernel_remove_quantum_friction(double* d_qf_density_for_Ua, double* d_qf_density_for_Ub, Complex*d_qf_density_for_D,
                                             double *V_a, double *V_b, Complex *delta)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    if(ixyz<NXYZ)
    {
        // restore values without quantum friction
        V_a[ixyz]=d_qf_density_for_Ua[ixyz];
        V_b[ixyz]=d_qf_density_for_Ub[ixyz];
        delta[ixyz]=d_qf_density_for_D[ixyz];
    }
}


__global__ void kernel_compute_qpe(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *re)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    if(ixyz<NXYZ)
    {
        p=thrust::conj(wf1_u[ixyz])*wf2_u[ixyz] + thrust::conj(wf1_v[ixyz])*wf2_v[ixyz];
        re[ixyz]=p.real()*DXYZ;
    }
}
__global__ void kernel_compute_qpe(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *re, int noAllElements)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    if(ixyz<noAllElements)
    {
        p=thrust::conj(wf1_u[ixyz])*wf2_u[ixyz] + thrust::conj(wf1_v[ixyz])*wf2_v[ixyz];
        re[ixyz]=p.real()*DXYZ;
    }
}

__global__ void kernel_compute_qpe_norm(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *qpe_re, double *norm_re)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    Complex _wf1_u, _wf1_v;
    if(ixyz<NXYZ)
    {
        // qpe
        _wf1_u=wf1_u[ixyz];
        _wf1_v=wf1_v[ixyz];
        p=thrust::conj(_wf1_u)*wf2_u[ixyz] + thrust::conj(_wf1_v)*wf2_v[ixyz];
        qpe_re[ixyz]=p.real()*DXYZ;
        // norm
        norm_re[ixyz]=(thrust::norm(_wf1_u)+thrust::norm(_wf1_v))*DXYZ;
    }
}

__global__ void kernel_subtruct_qpe(int n, Complex *wf, Complex *Hwf, double *qpe)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double e;
    int iwf;
    if(ixyz<NXYZ)
    {
        for(iwf=0; iwf<n; iwf++)
        {
            e=qpe[iwf];
            Hwf[ixyz+iwf*NXYZ       ]-=wf[ixyz+iwf*NXYZ       ]*e;
            Hwf[ixyz+iwf*NXYZ+n*NXYZ]-=wf[ixyz+iwf*NXYZ+n*NXYZ]*e;
        }
    }
}

__global__ void kernel_subtruct_qpe_norm(int n, Complex *wf, Complex *Hwf, double *qpe, double *norm)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double e, c;
    int iwf;
    Complex _u, _v;
    if(ixyz<NXYZ)
    {
        for(iwf=0; iwf<n; iwf++)
        {
            e=qpe[iwf];
            c=1.0/sqrt(norm[iwf]);

            // normalize wf
            _u=wf[ixyz+iwf*NXYZ       ]; _u*=c;
            _v=wf[ixyz+iwf*NXYZ+n*NXYZ]; _v*=c;

            //subtruct <H>, note: <c*psi|H|c*psi>=c*c*<H>
            Hwf[ixyz+iwf*NXYZ       ]=(Hwf[ixyz+iwf*NXYZ       ]*c - _u*e*c*c);
            Hwf[ixyz+iwf*NXYZ+n*NXYZ]=(Hwf[ixyz+iwf*NXYZ+n*NXYZ]*c - _v*e*c*c);

            // save normalized wave-functions
            wf[ixyz+iwf*NXYZ       ]=_u;
            wf[ixyz+iwf*NXYZ+n*NXYZ]=_v;
        }
    }
}

#ifdef ENABLE_VELOCITY_EXT
__global__ void kernel_get_vector_vext(int it, int spin, double *vx, double *vy, double *vz)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;

    if(ixyz<NXYZ)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); // decode cartesian coordinates

        vx[ixyz]=velocity_ext(ix, iy, iz, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        vy[ixyz]=velocity_ext(ix, iy, iz, it, spin, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        vz[ixyz]=velocity_ext(ix, iy, iz, it, spin, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data);

//         if(ixyz==0) printf("TEST1: %12.6f %12.6f %12.6f %12.6f\n", 0.1*it, vx[ixyz], vy[ixyz], vz[ixyz]);
    }
}
#endif

/**
 * Function applies hamiltonian (H-<H>)*Psi.
 * NOTE: this function executes cuFFT
 * @param it  iteration number
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf_in array with wave-functions (INPUT)
 * @param wf_out array with wave-functions (OUTPUT),
 *               It can be the same as wf_d_dx, ..., wf_d_laplace, but it CANNOT be wf_in
 * @param wf_d_dx derivative with respect to dx (INPUT)
 * @param wf_d_dy derivative with respect to dy (INPUT)
 * @param wf_d_dz derivative with respect to dz (INPUT)
 * @param wf_d_laplace laplace of wave-functions (INPUT)
 * @param d_densites (INPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NXYZ,
 *                          nu - double complex array of size NXYZ,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NXYZ
 *                   In total size of d_densites is 12*NXYZ
 * @param d_potentials (INPUT)
 *                     collective array with potentials [V_a, V_b, delta]
 *                     where: V_a, V_b - double arrays of size NXYZ
 *                            delta - double complex array of size NXYZ
 *                     In total size of d_potentials is 4*NXYZ
 * @param qfswitch switch coefficient for quantum friction term, if qfswitch=0.0 then quantum friction is NOT active
 * @param useqpe array of size [n]
 *               if NULL then quasiparticle energies will be computed from wf_in,
 *               otherwise given array will be used,
 * @param pccoeff the particle control coefficient: N(t) - N_req
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 * */
extern "C" int apply_hamiltonian(int it, int n, cufftDoubleComplex *wf_in, cufftDoubleComplex *wf_out,
                            cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *wf_d_dz, 
                            cufftDoubleComplex *wf_laplace, cufftDoubleComplex *alphawf_laplace,
                            double *d_densities, double *d_potentials, double qfswitch, double *useqpe, double pccoeff,
                            int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXYZ/nthreads);
    int ierr;

    wslda_density densall=convert_into_wslda_density(d_densities, NXYZ);
    wslda_potential potsall=convert_into_wslda_potential(d_potentials, NXYZ, NULL);

    // Set pointers for to simplify notation
    // densities
//     Complex *nu   =(Complex *)(d_densities +  0*NXYZ);
    double *rho_a = (double *)(d_densities +  2*NXYZ);
//     double *tau_a = (double *)(d_densities +  3*NXYZ);
    double *j_a_x = (double *)(d_densities +  4*NXYZ);
    double *j_a_y = (double *)(d_densities +  5*NXYZ);
    double *j_a_z = (double *)(d_densities +  6*NXYZ);
    double *rho_b = (double *)(d_densities +  7*NXYZ);
//     double *tau_b = (double *)(d_densities +  8*NXYZ);
    double *j_b_x = (double *)(d_densities +  9*NXYZ);
    double *j_b_y = (double *)(d_densities + 10*NXYZ);
    double *j_b_z = (double *)(d_densities + 11*NXYZ);
    // pontentials
    double *V_a = (double *)(d_potentials +  0*NXYZ);
    double *V_b = (double *)(d_potentials +  1*NXYZ);
    Complex *delta = (Complex *)(d_potentials +  2*NXYZ);

    // DERIVATIVES
    double * grad_alpha_a  = (double *)mm_get_pointer_to_kernels_workspace(NXYZ); // Use here work area of cuFFT
    double * grad_alpha_b  = grad_alpha_a + NXYZ*3;
    double * grad_j_corr_a = grad_alpha_a + NXYZ*6; // (j+/n+ - alpha_a*ja/na), see GW notes
    double * grad_j_corr_b = grad_alpha_a + NXYZ*9; // (j+/n+ - alpha_b*jb/nb), see GW notes
    double * laplace_alpha_a  = grad_alpha_a + NXYZ*12;
    double * laplace_alpha_b  = grad_alpha_a + NXYZ*13;

    // Step 1: if quantum friction is active, update mean-field potentials
    if(qfswitch>0.0)
    {
        double qfalpha = md.qfalpha*qfswitch;
        double qfbeta = md.qfbeta*qfswitch;
        double qfgamma = md.qfgamma*qfswitch*pccoeff;

        double * d_qf_density_for_Ua = (double *) (d_densities+12*NXYZ);  // density for diagonal part (U) of quantum friction force
        double * d_qf_density_for_Ub = (double *) (d_qf_density_for_Ua+1*NXYZ);
        Complex *d_qf_density_for_D = (Complex *) (d_qf_density_for_Ua+2*NXYZ); // density for off-diagonal part (Delta) of quantum friction force
    
        
        // update mean field potential by friction terms
        kernel_add_quantum_friction<<<nblocks, nthreads>>>(rho_a, rho_b,
                                                    d_qf_density_for_Ua, d_qf_density_for_Ub, d_qf_density_for_D,
                                                    V_a, V_b, delta, qfalpha, qfbeta, qfgamma);
        
    }

    // filtering of mean-fields
    if(md.hkf_mode>=1)
    {
        ierr = high_frequency_filter_massive_d(2, V_a, V_a, md.hkf_mu, md.hkf_T, nthreads);
        if(ierr!=0) return ierr;
    }

    double *vecvext_a    = NULL;
    double *divvext_a    = NULL;
    double *vecvext_b    = NULL;
    double *divvext_b    = NULL;
#ifdef ENABLE_VELOCITY_EXT
    // set pointers
    vecvext_a    = (double *)(grad_alpha_a + 14*NXYZ); // storage for keeping vext=[vx(r),vy(r),vz(r)]
    divvext_a    = (double *)(grad_alpha_a + 17*NXYZ); // storage for keeping div(vext)
    vecvext_b    = (double *)(grad_alpha_a + 18*NXYZ); // storage for keeping vext=[vx(r),vy(r),vz(r)]
    divvext_b    = (double *)(grad_alpha_a + 21*NXYZ); // storage for keeping div(vext)

    // fill arrays with data
    kernel_get_vector_vext<<<nblocks, nthreads>>>(it, SPINA, vecvext_a, vecvext_a+NXYZ, vecvext_a+2*NXYZ); // NOTE - only SPINA
    ierr=compute_divergence_real_vector_f(vecvext_a, vecvext_a+NXYZ, vecvext_a+2*NXYZ, divvext_a, nthreads);
    if(ierr!=0) return ierr+300;

    kernel_get_vector_vext<<<nblocks, nthreads>>>(it, SPINB, vecvext_b, vecvext_b+NXYZ, vecvext_b+2*NXYZ); // NOTE - only SPINB
    ierr=compute_divergence_real_vector_f(vecvext_b, vecvext_b+NXYZ, vecvext_b+2*NXYZ, divvext_b, nthreads);
    if(ierr!=0) return ierr+400;
#endif

    if(md.hkf_mode>=2)
    {
        ierr = high_frequency_filter_c(potsall.delta, potsall.delta, md.hkf_mu, md.hkf_T, nthreads);
        if(ierr!=0) return ierr;
    }

#ifdef BDG_MODE
    // Step 3: apply hamiltonian
    kernel_apply_hamiltonian_bdg<<<nblocks, nthreads>>>(it,
                                            V_a, V_b, delta,
                                            n, (Complex *)wf_in, (Complex *)wf_out,
                                            (Complex *)wf_d_dx, (Complex *)wf_d_dy, (Complex *)wf_d_dz, (Complex *)wf_laplace,
                                            vecvext_a, vecvext_a+NXYZ, vecvext_a+2*NXYZ, divvext_a, vecvext_b, vecvext_b+NXYZ, vecvext_b+2*NXYZ, divvext_b
                                                   );
#else

#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
    // and compute gradient and laplace of effective mass
    ierr=compute_laplace_real_f(potsall.alpha_a, laplace_alpha_a, nthreads);
    if(ierr!=0) return ierr;
    ierr=compute_laplace_real_f(potsall.alpha_b, laplace_alpha_b, nthreads);
    if(ierr!=0) return ierr;
#endif

#ifdef CURRENT_CORRECTIONS
    // filtering of vector potentials
    if(md.hkf_mode>=1)
    {
        ierr = high_frequency_filter_massive_d(6, potsall.A_a_x, potsall.A_a_x, md.hkf_mu, md.hkf_T, nthreads);
        if(ierr!=0) return ierr;
    }
    ierr=compute_derivative_real_vector_f(potsall.A_a_x, potsall.A_a_y, potsall.A_a_z, grad_j_corr_a, grad_j_corr_a+NXYZ, grad_j_corr_a+NXYZ*2, nthreads);
    if(ierr!=0) return ierr;
    ierr=compute_derivative_real_vector_f(potsall.A_b_x, potsall.A_b_y, potsall.A_b_z, grad_j_corr_b, grad_j_corr_b+NXYZ, grad_j_corr_b+NXYZ*2, nthreads);
    if(ierr!=0) return ierr;
#endif

    // Step 3: apply hamiltonian
    kernel_apply_hamiltonian<<<nblocks, nthreads>>>(it, potsall,
                                            laplace_alpha_a, laplace_alpha_b,
                                            grad_j_corr_a, grad_j_corr_a+NXYZ, grad_j_corr_a+NXYZ*2, grad_j_corr_b, grad_j_corr_b+NXYZ, grad_j_corr_b+NXYZ*2,
                                            n, (Complex *)wf_in, (Complex *)wf_out,
                                            (Complex *)wf_d_dx, (Complex *)wf_d_dy, (Complex *)wf_d_dz, (Complex *)wf_laplace, (Complex *)alphawf_laplace,
                                            vecvext_a, vecvext_a+NXYZ, vecvext_a+2*NXYZ, divvext_a, vecvext_b, vecvext_b+NXYZ, vecvext_b+2*NXYZ, divvext_b
                                                   );
#endif
    // Step 4: to increas stability - subtruct <H>*wf, where <H> is quasi particle energy
    // and normalize
    double *gpe;

    if(useqpe==NULL) // compute quasiparticle energies
    {
        gpe=(double *)mm_get_pointer_to_total_workspace(NXYZ); // for easier notation

        // compute quasi-particle energy for each wave-function
        int noAllElements = n*NXYZ;
        nblocks = (int)ceil((float)noAllElements/nthreads);
        kernel_compute_qpe<<<nblocks, nthreads>>>(  (Complex *)wf_in,               (Complex *)wf_out, 
                                                    (Complex *)wf_in+noAllElements, (Complex *)wf_out+noAllElements, 
                                                    gpe, noAllElements);

        ierr = local_reductions_many(n, NXYZ, gpe, gpe);
        if(ierr!=0) return ierr;

    }
    else // use given values
    {
        gpe=useqpe;
    }

    // subtruct <H>*wf
    nblocks = (int)ceil((float)NXYZ/nthreads);
    kernel_subtruct_qpe<<<nblocks, nthreads>>>(n, (Complex *)wf_in, (Complex *)wf_out, gpe);  
 

    // if quantum friction was active, remove contributions to mean-fields
    if(qfswitch>0.0)
    {
        double * d_qf_density_for_Ua = (double *) (d_densities+12*NXYZ);  // density for diagonal part (U) of quantum friction force
        double * d_qf_density_for_Ub = (double *) (d_qf_density_for_Ua+1*NXYZ);
        Complex *d_qf_density_for_D = (Complex *) (d_qf_density_for_Ua+2*NXYZ); // density for off-diagonal part (Delta) of quantum

        nblocks = (int)ceil((float)NXYZ/nthreads);
        kernel_remove_quantum_friction<<<nblocks, nthreads>>>(d_qf_density_for_Ua, d_qf_density_for_Ub, d_qf_density_for_D,
                                                              V_a, V_b, delta);
    }

    return 0;
}

// ================================================================================================
// ========================================= compute_ovelap =======================================
// ================================================================================================
__global__ void kernel_compute_ovelap(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *re, double *im)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    if(ixyz<NXYZ)
    {
        p=thrust::conj(wf1_u[ixyz])*wf2_u[ixyz] + thrust::conj(wf1_v[ixyz])*wf2_v[ixyz];
        re[ixyz]=p.real()*DXYZ;
        im[ixyz]=p.imag()*DXYZ;
    }
}
/**
 * Function computes overlaps between two wave-functions (wf1,wf2)
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf1 array with wave-functions (INPUT)
 * @param wf2 array with wave-functions (INPUT)
 * @param overlap_re computed values of overlaps, real parts, array of size n (INPUT/OUTPUT)
 *                   If pointer is set as NULL on input then computation of real part is skipped.
 * @param overlap_im computed values of overlaps, real parts, array of size n (INPUT/OUTPUT)
 *                   If pointer is set as NULL on input then computation of imaginary part is skipped.
 * @param workarea work space of size 2*NXYZ, it can be the same workspace as used by cuFFT
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 * */
extern "C" int compute_ovelap(int n, cufftDoubleComplex *wf1, cufftDoubleComplex *wf2, double *overlap_re, double *overlap_im,
                              double *workarea, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXYZ/nthreads);

    double *re = workarea;
    double *im = workarea + NXYZ;
    int iwf;
    size_t shift=0;
    int ierr;

    for(iwf=0; iwf<n; iwf++) // for each wave-function
    {
        kernel_compute_ovelap<<<nblocks, nthreads>>>((Complex *)wf1+shift, (Complex *)wf2+shift, (Complex *)wf1+shift+n*NXYZ, (Complex *)wf2+shift+n*NXYZ, re, im);
        shift+=NXYZ; // move pointer to next wf

        ierr = local_reductions_many(2, NXYZ, workarea, workarea);
        if(ierr!=0) return ierr;

        if(overlap_re!=NULL)
        {
            if( cudaMemcpy( overlap_re+iwf , workarea+0 , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -333;
        }

        if(overlap_im!=NULL)
        {
            if( cudaMemcpy( overlap_im+iwf , workarea+0 , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -334;
        }
    }
    return 0;
}

// =======================================================================================
// ==================================== normalize_wf =====================================
// =======================================================================================
__global__ void kernel_compute_norm(Complex *wf_u, Complex *wf_v, double *norm)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    if(ixyz<NXYZ)
    {
        norm[ixyz]=(thrust::norm(wf_u[ixyz])+thrust::norm(wf_v[ixyz]))*DXYZ;
    }
}
__global__ void kernel_compute_norm(Complex *wf_u, Complex *wf_v, double *norm, int noAllElements)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    if(ixyz<noAllElements)
    {
        norm[ixyz]=(thrust::norm(wf_u[ixyz])+thrust::norm(wf_v[ixyz]))*DXYZ;
    }
}

__global__ void kernel_normalize_wf(int n, Complex *wf, double *norm)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double c;
    int iwf;
    if(ixyz<NXYZ)
    {
        for(iwf=0; iwf<n; iwf++)
        {
            c=1.0/sqrt(norm[iwf]);

            // normalize wf
            wf[ixyz+iwf*NXYZ       ]*=c;
            wf[ixyz+iwf*NXYZ+n*NXYZ]*=c;
        }
    }
}


/**
 * Function normalizes wave-functions
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf pointer to wave-functions
 * @return 0 - OK, otherwise ERROR
 * */
extern "C" int normalize_wf(int n, cufftDoubleComplex *wf, int nthreads)
{
    // number of blocks
    int ierr;
    int noAllElements = NXYZ*n; 
    double * norm  = (double *)mm_get_pointer_to_total_workspace(NXYZ); // Use here work area of cuFFT as working buffer
    
    // compute norm for each wave-function
    int nblocks = (int)ceil((float)noAllElements/nthreads);
    kernel_compute_norm<<<nblocks, nthreads>>>((Complex *)wf, (Complex *)wf+noAllElements, norm, noAllElements);

    // reduce wave-functions pararell
    ierr = local_reductions_many(n, NXYZ, norm, norm);
    if(ierr != 0) return ierr;
    
    // normalize wf
    nblocks = (int)ceil((float)NXYZ/nthreads);
    kernel_normalize_wf<<<nblocks, nthreads>>>(n, (Complex *)wf, norm); 

    return 0;

}

// =======================================================================================
// ============================== multiply_wf_by_alpha ===================================
// =======================================================================================
__global__ void kernel_multiply_wf_by_alpha(int n, double * alpha_a, double * alpha_b, Complex *wf_in, Complex *wf_out)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    // registers
    int iwf;
    double na, nb;
    Complex u, v;

    if(ixyz<NXYZ)
    {
        na=alpha_a[ixyz];
        nb=alpha_b[ixyz];

        for(iwf=0; iwf<n; iwf++)
        {
            // read u and v
            u=wf_in[       iwf*NXYZ+ixyz];
            v=wf_in[n*NXYZ+iwf*NXYZ+ixyz];

            wf_out[       iwf*NXYZ+ixyz]=u*na; // multiply by effective mass and save
            wf_out[n*NXYZ+iwf*NXYZ+ixyz]=v*nb; // multiply by effective mass and save
        }
    }
}

/**
 * Function applies hamiltonian.
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf_in array with wave-functions (INPUT)
 * @param wf_out array with wave-functions (OUTPUT),
 * @param d_potentials collective array with potentials (INPUT)
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 **/
extern "C" int multiply_wf_by_alpha(int n, cufftDoubleComplex *wf_in, cufftDoubleComplex *wf_out, double *d_potentials, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXYZ/nthreads);

    wslda_potential potsall=convert_into_wslda_potential(d_potentials, NXYZ, NULL);

    kernel_multiply_wf_by_alpha<<<nblocks, nthreads>>>(n, potsall.alpha_a, potsall.alpha_b, (Complex *)wf_in, (Complex *)wf_out);

    return 0;
}


// =======================================================================================
// ======================== taylor_expansion_contribution ================================
// =======================================================================================
__global__ void kernel_taylor_expansion_contribution(int n, Complex *wf_hpsi, Complex *wf_update, Complex *wf_contr, double t_coeff)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    // registers
    int iwf;

    Complex u, v;

    if(ixyz<NXYZ)
    {

        for(iwf=0; iwf<n; iwf++)
        {
            // read u and v
            u=wf_hpsi[       iwf*NXYZ+ixyz];
            v=wf_hpsi[n*NXYZ+iwf*NXYZ+ixyz];

            // add Taylor coefficients
            u*=Complex(0.0, t_coeff);
            v*=Complex(0.0, t_coeff);


            wf_contr[       iwf*NXYZ+ixyz]=u; // save taylor contribution
            wf_contr[n*NXYZ+iwf*NXYZ+ixyz]=v; // save taylor contribution

            wf_update[       iwf*NXYZ+ixyz]+=u; // add to rest of Taylor expansion
            wf_update[n*NXYZ+iwf*NXYZ+ixyz]+=v; // add to rest of Taylor expansion
        }
    }
}

/**
 * Function add it-order conytribution from Taylor expansion: (1/it!)(-i*(H-<H>)dt)*psi
 * @param it Taylor order contribution
 * @param dt time step
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf_hpsi result of apply_hamiltonian (INPUT)
 * @param wf_update array with wave-functions constructed by Taylor expansion (OUTPUT),
 * @param wf_contr array with wave-functions constructed by Taylor expansion (OUTPUT),
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 **/
extern "C" int taylor_expansion_contribution(int it, double dt, int n, cufftDoubleComplex *wf_hpsi,
                                             cufftDoubleComplex *wf_update, cufftDoubleComplex *wf_contr, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXYZ/nthreads);

    kernel_taylor_expansion_contribution<<<nblocks, nthreads>>>(n, (Complex *)wf_hpsi, (Complex *)wf_update, (Complex *)wf_contr,
                                                                -1.0*dt/(double)(it));

    return 0;
}

// =======================================================================================
// ======================== get_ext_potentials - handlers ================================
// =======================================================================================

__global__ void kernel_get_v_ext(int it, int spin, double *data)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;

    if(ixyz<NXYZ)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); // decode cartesian coordinates
        data[ixyz]=u_ext(ix,iy,iz,it,spin);
    }
}

/**
 * @return array with external potential
 * */
extern "C" int get_v_ext(int datadim, int spin, int it, double *data)
{
    int nthreads = 256;
    int nblocks = (int)ceil((float)NXYZ/nthreads);

    // allocate cuda memory
    double *wrkspace = (double *)mm_get_pointer_to_kernels_workspace(NXYZ); // reuse workspace
    kernel_get_v_ext<<<nblocks, nthreads>>>(it, spin, wrkspace);

    if( cudaMemcpy( data , wrkspace, sizeof(double)*NXYZ, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;

    return 0;
}

__global__ void kernel_get_delta_ext(int it, Complex *deltain, Complex *data)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;

    if(ixyz<NXYZ)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); // decode cartesian coordinates
        data[ixyz]=macro_delta_ext(ix, iy, iz, it, deltain[ixyz]);
    }
}

/**
 * @return array with external potential
 * */
extern "C" int get_delta_ext(int datadim, int it, void *deltain, void *data)
{
    int nthreads = 256;
    int nblocks = (int)ceil((float)NXYZ/nthreads);

    // allocate cuda memory
    Complex *wrkspace = (Complex *)mm_get_pointer_to_kernels_workspace(NXYZ); // reuse workspace
    kernel_get_delta_ext<<<nblocks, nthreads>>>(it, (Complex *)deltain, (Complex *)wrkspace);

    if( cudaMemcpy( data , wrkspace, sizeof(Complex)*NXYZ, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;

    return 0;
}

__global__ void kernel_get_velocity_ext(int it, int spin, double *datax, double *datay, double *dataz)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;

    if(ixyz<NXYZ)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); // decode cartesian coordinates
        datax[ixyz]=velocity_ext(ix, iy, iz, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        datay[ixyz]=velocity_ext(ix, iy, iz, it, spin, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        dataz[ixyz]=velocity_ext(ix, iy, iz, it, spin, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data);
    }
}

/**
 * @return array with external potential
 * */
extern "C" int get_velocity_ext(int datadim, int spin, int it, double *data)
{
    int nthreads = 256;
    int nblocks = (int)ceil((float)NXYZ/nthreads);

    double *wrkspace = (double *)mm_get_pointer_to_kernels_workspace(NXYZ); // reuse workspace
    double *vx = wrkspace + 0*NXYZ;
    double *vy = wrkspace + 1*NXYZ;
    double *vz = wrkspace + 2*NXYZ;

    kernel_get_velocity_ext<<<nblocks, nthreads>>>(it, spin, vx, vy, vz);

    if( cudaMemcpy( data , wrkspace, sizeof(double)*NXYZ*3, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;

    return 0;
}
