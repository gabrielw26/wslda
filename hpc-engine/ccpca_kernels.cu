// Author: Gabriel Wlazlowski
// Date: 23-04-2016

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/complex.h>
#include <stdio.h>
typedef thrust::complex<double> Complex;

#include "pca_settings.h"
#define CODEDIM 1
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


#ifdef TDWSLDA

#include "problem-definition.h"

#ifdef ENABLE_V_EXT
#define u_ext(ix, iy, iz, it, spin) v_ext(ix, 0, 0, it, spin, dc_params, dc_extra_data_size, dc_extra_data)
#else
#define u_ext(ix, iy, iz, it, spin) 0.0
#endif

#ifdef ENABLE_DELTA_EXT
#define macro_delta_ext(ix, iy, iz, it, delta) delta_ext(ix, 0, 0, it, delta, dc_params, dc_extra_data_size, dc_extra_data)
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
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NX,
 *                          nu - double complex array of size NX,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NX
 *                   In total size of d_densites is 12*NX
 * @param d_potentials (INPUT/OUTPUT)
 *                     collective array with potentials [V_a, V_b, delta]
 *                     where: V_a, V_b - double arrays of size NX
 *                            delta - double complex array of size NX
 *                     In total size of d_potentials is 4*NX
 *                     NOTE: I assume that d_potentials contains potentials from previous iteration,
 *                           i.e. they are good starting point for self-consistent process.
 * @param cccoeff the current corrections coefficient
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 **/
extern "C" int compute_potentials(int it, double *d_densities, double *d_potentials, double cccoeff, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);
    wslda_density densall=convert_into_wslda_density(d_densities, NX);
    wslda_potential potsall=convert_into_wslda_potential(d_potentials, NX, NULL);

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
    int ix;

    if(ixyz<NXY)
    {
        ix = ixyz; // decode cartesian coordinates

        // External potential energy
        E_ext[ixyz]=(rho_a[ixyz]*u_ext(ix,0,0,it,SPINA) + rho_b[ixyz]*u_ext(ix,0,0,it,SPINB))*DXYZ*NY*NZ;
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
    int ix;

    double na, nb;

    if(ixyz<NXY)
    {
        ix = ixyz; // decode cartesian coordinates

        // External potential energy
        E_ext[ixyz]=(thrust::conj(nu[ixyz])*macro_delta_ext(ix, 0, 0, it, delta[ixyz])).real()*(-2.0)*DXYZ*NY*NZ;
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
    int ix;

    if(ixyz<NXY)
    {
        ix = ixyz; // decode cartesian coordinates

        // External potential energy
        E_ext[ixyz]=  (
                         j_a_x[ixyz]*velocity_ext(ix, 0, 0, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data)
                        +j_b_x[ixyz]*velocity_ext(ix, 0, 0, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data)
                       )*(-1.0)*DXYZ*NY*NZ;

    }
}


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
extern "C" void *pca_cufft_work_area;
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
                                         Complex *wf_d_dx, double *d_kky, double *d_kkz, Complex *wf_laplace, Complex *alphawf_laplace,
                                         double cccoeff,
                                         double *vx_a, double *divv_a, double *vx_b, double *divv_b
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double aa, ab;
    double na, nb;
    double Va, Vb;
    Complex D;
    double cja=0.0, cjb=0.0;
    double ja, jb/*, jp*/;
#ifdef CURRENT_CORRECTIONS
    Complex gax, gbx;
    double /*fr,*/ fra, frb; // regularization functions
#endif

#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
    int ix;
#endif
    size_t iwf;
    Complex u, v, tu, tv;

    if(ixyz<NX)
    {
#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
        ix=ixyz;
#endif
        // effective mass
        aa=h_potentials.alpha_a[ixyz]; // will be multiplied by -1/2 later
        ab=h_potentials.alpha_b[ixyz]; // will be multiplied by -1/2 later

        // read potentials
        Va=h_potentials.V_a[ixyz];
        Vb=h_potentials.V_b[ixyz];
        D =h_potentials.delta[ixyz];

#ifdef ENABLE_DELTA_EXT
        D = D + macro_delta_ext(ix, 0, 0, it, D);
#endif

#ifdef CURRENT_CORRECTIONS
        // reset
        gax=Complex(0.0, 0.0);
        gbx=Complex(0.0, 0.0);
#endif

        // read gradient corrections
#ifdef CURRENT_CORRECTIONS
        cja+=-0.5*(j_corr_a_x[ixyz]);
        cjb+=-0.5*(j_corr_b_x[ixyz]);

        gax+=Complex(0.0,  -1.*h_potentials.A_a_x[ixyz]);
        gbx+=Complex(0.0,   1.*h_potentials.A_b_x[ixyz]); // note conjugate of complex number (beacuse of "-h*" operator)
#endif

#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
        // multiply mass by -1/2
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
            tu=wf_in[      iwf*NX+ixyz];
            tv=wf_in[n*NX +iwf*NX+ixyz];

#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
            // ja = (1/2)*alpha_a*kz^2, NOTE aa has been multiplied by -0.5 already
            ja = d_kky[iwf]; // kz, jb as working buffer
            jb = d_kkz[iwf]; // kz, jb as working buffer
            jb = jb*jb + ja*ja;      // ky^2 + kz^2, jb as working buffer
            ja=-1.0*aa*jb;
            jb=-1.0*ab*jb;
#else
            // ja = (1/2)*alpha_a*kz^2, NOTE aa has been multiplied by -0.25 already
            ja = d_kky[iwf]; // kz, jb as working buffer
            jb = d_kkz[iwf]; // kz, jb as working buffer
            jb = jb*jb + ja*ja;      // ky^2 + kz^2, jb as working buffer
            ja=-2.0*aa*jb;
            jb=-2.0*ab*jb;
#endif
            u+=tu*Complex(Va-dc_mu_a+ja,     cja); // V_a*u
            v-=tv*Complex(Vb-dc_mu_b+jb,-1.0*cjb); // V_b*v, note conjugate of complex number (beacuse of "-h*" operator)

            u+=tv*D;                // delta   * v
            v+=tu*thrust::conj(D);  // delta^* * u

#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // kinetic part - contribution: (1/4)*laplace(alpha)*u, see further
            u+=tu*na; // note: na keeps (1/4)*laplace(alpha_a)
            v-=tv*nb; // note: nb keeps (1/4)*laplace(alpha_b)
#endif

#ifdef CURRENT_CORRECTIONS
            // read gradients of wf: x-coordinate
            tu=wf_d_dx[      iwf*NX+ixyz];
            tv=wf_d_dx[n*NX +iwf*NX+ixyz];
            u+=tu*gax;
            v-=tv*gbx;
#endif


            // kinetic part
            // note: (-1/2)\nabla(alpha * nabla u) = -(1/4)*laplace(alpha*u) - (1/4)*alpha*laplace(u) + (1/4)*laplace(alpha)*u
            // (1/4)*laplace(alpha)*u already done
            // read laplace of wf
            tu=wf_laplace[      iwf*NX+ixyz];
            tv=wf_laplace[n*NX +iwf*NX+ixyz];
            u+=tu*aa; // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_a is given by aa
            v-=tv*ab; // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_b is given by ab
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // read laplace(alpha*u)
            u+=alphawf_laplace[      iwf*NX+ixyz]*(-0.25);
            v-=alphawf_laplace[n*NX +iwf*NX+ixyz]*(-0.25);
#endif
            // save to global memory wf
            wf_out[      iwf*NX+ixyz]=u;
            wf_out[n*NX +iwf*NX+ixyz]=v;
        }
    }
}

__global__ void kernel_apply_hamiltonian_bdg(int it,
                                         double *V_a, double *V_b, Complex *delta,
                                         size_t n, Complex *wf_in, Complex *wf_out,
                                         Complex *wf_d_dx, double *d_kky, double *d_kkz, Complex *wf_laplace,
                                         double *vx_a, double *divv_a, double *vx_b, double *divv_b
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    double p,k;
    double Va, Vb;
    Complex D;

    size_t iwf;
    Complex u, v, tu, tv;

#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
    int ix;
#endif

    double cja=0.0, cjb=0.0;
#ifdef CURRENT_CORRECTIONS
    Complex gax;
    Complex gbx;
#endif

    if(ixyz<NX)
    {
#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
        ix=ixyz;
#endif
        // read potentials
        Va=V_a[ixyz];
        Vb=V_b[ixyz];
        D=delta[ixyz];

#ifdef ENABLE_DELTA_EXT
        D = D + macro_delta_ext(ix, 0, 0, it, D);
#endif

#ifdef ENABLE_VELOCITY_EXT
        gax=Complex(0.0, 1.0*vx_a[ixyz]);
        gbx=Complex(0.0,-1.0*vx_b[ixyz]); // NOTE: complex conjugate included

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
            tu=wf_in[      iwf*NX+ixyz];
            tv=wf_in[n*NX +iwf*NX+ixyz];

            // p = (1/2)*kz^2
            k = d_kky[iwf]; // ky
            p = d_kkz[iwf]; // kz
            p = 0.5*(p*p+k*k);    // kz^2/2 + ky^2/2

            u+=tu*Complex(Va-dc_mu_a+p,       cja); // V_a*u
            v-=tv*Complex(Vb-dc_mu_b+p,  -1.0*cjb); // V_b*v, note conjugate of complex number (beacuse of "-h*" operator)

            u+=tv*D;                // delta   * v
            v+=tu*thrust::conj(D);  // delta^* * u

#ifdef CURRENT_CORRECTIONS
            // read gradients of wf: x-coordinate
            tu=wf_d_dx[      iwf*NX+ixyz];
            tv=wf_d_dx[n*NX +iwf*NX+ixyz];
            u+=tu*gax;
            v-=tv*gbx;
#endif

            // kinetic part

            // read laplace of wf
            tu=wf_laplace[      iwf*NX+ixyz];
            tv=wf_laplace[n*NX +iwf*NX+ixyz];
            u+=tu*(-0.5); // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_a is given by aa
            v-=tv*(-0.5); // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_b is given by ab

            // save to global memory wf
            wf_out[      iwf*NX+ixyz]=u;
            wf_out[n*NX +iwf*NX+ixyz]=v;
        }
    }
}

__global__ void kernel_add_quantum_friction(double *rho_a, double *rho_b,
                                            double *djax_dx, double *djay_dy, double *djaz_dz, double *djbx_dx, double *djby_dy, double *djbz_dz,
                                            double *V_a, double *V_b, double qfalpha)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    int ix;
    double coeff=0.0, r;
    if(ixyz<NX)
    {
        // see Eq.(3) in paper https://arxiv.org/abs/1305.6891
        V_a[ixyz]-=qfalpha*(djax_dx[ixyz])/dc_nF; // here I divide be reference density, to avoid problems of division by zero
        V_b[ixyz]-=qfalpha*(djbx_dx[ixyz])/dc_nF; // here I divide be reference density, to avoid problems of division by zero
    }
}

__global__ void kernel_compute_qpe(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *re)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    if(ixyz<NX)
    {
        p=thrust::conj(wf1_u[ixyz])*wf2_u[ixyz] + thrust::conj(wf1_v[ixyz])*wf2_v[ixyz];
        re[ixyz]=p.real()*DX;
    }
}
__global__ void kernel_compute_qpe(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *re, double noAllElements)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    if(ixyz<noAllElements)
    {
        p=thrust::conj(wf1_u[ixyz])*wf2_u[ixyz] + thrust::conj(wf1_v[ixyz])*wf2_v[ixyz];
        re[ixyz]=p.real()*DX;
    }
}


__global__ void kernel_compute_qpe_norm(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *qpe_re, double *norm_re)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    Complex _wf1_u, _wf1_v;
    if(ixyz<NX)
    {
        // qpe
        _wf1_u=wf1_u[ixyz];
        _wf1_v=wf1_v[ixyz];
        p=thrust::conj(_wf1_u)*wf2_u[ixyz] + thrust::conj(_wf1_v)*wf2_v[ixyz];
        qpe_re[ixyz]=p.real()*DX;
        // norm
        norm_re[ixyz]=(thrust::norm(_wf1_u)+thrust::norm(_wf1_v))*DX;
    }
}

__global__ void kernel_subtruct_qpe(int n, Complex *wf, Complex *Hwf, double *qpe)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double e;
    int iwf;
    if(ixyz<NX)
    {
        for(iwf=0; iwf<n; iwf++)
        {
            e=qpe[iwf];
            Hwf[ixyz+iwf*NX      ]-=wf[ixyz+iwf*NX      ]*e;
            Hwf[ixyz+iwf*NX +n*NX]-=wf[ixyz+iwf*NX +n*NX]*e;
        }
    }
}

__global__ void kernel_subtruct_qpe_norm(int n, Complex *wf, Complex *Hwf, double *qpe, double *norm)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double e, c;
    int iwf;
    Complex _u, _v;
    if(ixyz<NX)
    {
        for(iwf=0; iwf<n; iwf++)
        {
            e=qpe[iwf];
            c=1.0/sqrt(norm[iwf]);

            // normalize wf
            _u=wf[ixyz+iwf*NX      ]; _u*=c;
            _v=wf[ixyz+iwf*NX +n*NX]; _v*=c;

            //subtruct <H>, note: <c*psi|H|c*psi>=c*c*<H>
            Hwf[ixyz+iwf*NX      ]=(Hwf[ixyz+iwf*NX      ]*c - _u*e*c*c);
            Hwf[ixyz+iwf*NX +n*NX]=(Hwf[ixyz+iwf*NX +n*NX]*c - _v*e*c*c);

            // save normalized wave-functions
            wf[ixyz+iwf*NX      ]=_u;
            wf[ixyz+iwf*NX +n*NX]=_v;
        }
    }
}

#ifdef ENABLE_VELOCITY_EXT
__global__ void kernel_get_vector_vext(int it, int spin, double *vx)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;

    if(ixyz<NX)
    {
        ix=ixyz;

        vx[ixyz]=velocity_ext(ix, 0, 0, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);

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
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NX,
 *                          nu - double complex array of size NX,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NX
 *                   In total size of d_densites is 12*NX
 * @param d_potentials (INPUT)
 *                     collective array with potentials [V_a, V_b, delta]
 *                     where: V_a, V_b - double arrays of size NX
 *                            delta - double complex array of size NX
 *                     In total size of d_potentials is 4*NX
 * @param qfswitch switch coefficient for quantum friction term, if qfswitch=0.0 then quantum friction is NOT active
 * @param useqpe array of size [n]
 *               if NULL then quasiparticle energies will be computed from wf_in,
 *               otherwise given array will be used,
 * @param cccoeff the current corrections coefficient
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 * */
extern "C" int apply_hamiltonian(int it, int n, cufftDoubleComplex *wf_in, cufftDoubleComplex *wf_out,
                            cufftDoubleComplex *wf_d_dx, double *d_kkyz, cufftDoubleComplex *wf_laplace, cufftDoubleComplex *alphawf_laplace,
                            double *d_densities, double *d_potentials, double qfswitch, double *useqpe, double cccoeff,
                            int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);
    int ierr;

    wslda_density densall=convert_into_wslda_density(d_densities, NX);
    wslda_potential potsall=convert_into_wslda_potential(d_potentials, NX, NULL);

    // Set pointers for to simplify notation
    // densities
//     Complex *nu   =(Complex *)(d_densities +  0*NX);
    double *rho_a = (double *)(d_densities +  2*NX);
//     double *tau_a = (double *)(d_densities +  3*NX);
    double *j_a_x = (double *)(d_densities +  4*NX);
    double *j_a_y = (double *)(d_densities +  5*NX);
    double *j_a_z = (double *)(d_densities +  6*NX);
    double *rho_b = (double *)(d_densities +  7*NX);
//     double *tau_b = (double *)(d_densities +  8*NX);
    double *j_b_x = (double *)(d_densities +  9*NX);
    double *j_b_y = (double *)(d_densities + 10*NX);
    double *j_b_z = (double *)(d_densities + 11*NX);

    // pontentials
    double *V_a = (double *)(d_potentials +  0*NX);
    double *V_b = (double *)(d_potentials +  1*NX);
    Complex *delta = (Complex *)(d_potentials +  2*NX);

    double *d_kky = d_kkyz;
    double *d_kkz = d_kkyz+n;

    double * grad_alpha_a  = (double *)pca_cufft_work_area; // Use here work area of cuFFT
    double * grad_alpha_b  = grad_alpha_a + NX*3;
    double * grad_j_corr_a = grad_alpha_a + NX*6; // (j+/n+ - alpha_a*ja/na), see GW notes
    double * grad_j_corr_b = grad_alpha_a + NX*9; // (j+/n+ - alpha_b*jb/nb), see GW notes
    double * laplace_alpha_a  = grad_alpha_a + NX*12;
    double * laplace_alpha_b  = grad_alpha_a + NX*13;

    // Step 1: if quantum friction is active, update mean-field potentials
    if(qfswitch>0.0)
    {
        double qfalpha = md.qfalpha*qfswitch;
        double qfbeta = md.qfbeta*qfswitch;
        double qfgamma = md.qfgamma*qfswitch;

        double * d_qf_density_for_U = (double *)  d_densities+12*NXYZ;;  // density for diagonal part (U) of quantum friction force
        Complex *d_qf_density_for_D = (Complex *) d_qf_density_for_U + NXYZ; // density for off-diagonal part (Delta) of quantum friction force

        // TODO: EA: Update this section
        // TODO: For now I leave the old method, but you should replace it with computation via second derivatives
        // TODO: Here you need to update mean-field potentials (V_a,V_b) and pairing potential (Delta)

        // compute nabla*j, use grad_j_corr_a and grad_j_corr_b as temporary buffers
        ierr=compute_derivative_real_vector_f(j_a_x, NULL, NULL, grad_j_corr_a, NULL, NULL, nthreads);
        if(ierr!=0) return ierr;
        ierr=compute_derivative_real_vector_f(j_b_x, NULL, NULL, grad_j_corr_b, NULL, NULL, nthreads);
        if(ierr!=0) return ierr;

        // update mean field potential by friction term
        kernel_add_quantum_friction<<<nblocks, nthreads>>>(rho_a, rho_b,
                                                    grad_j_corr_a, NULL, NULL, grad_j_corr_b, NULL, NULL,
                                                    V_a, V_b, qfalpha);
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
    vecvext_a    = (double *)(grad_alpha_a + 14*NX); // storage for keeping vext=[vx(r),vy(r),vz(r)]
    divvext_a    = (double *)(grad_alpha_a + 17*NX); // storage for keeping div(vext)
    vecvext_b    = (double *)(grad_alpha_a + 18*NX); // storage for keeping vext=[vx(r),vy(r),vz(r)]
    divvext_b    = (double *)(grad_alpha_a + 21*NX); // storage for keeping div(vext)

    // fill arrays with data
    kernel_get_vector_vext<<<nblocks, nthreads>>>(it, SPINA, vecvext_a); // NOTE - only SPINA
    ierr=compute_divergence_real_vector_f(vecvext_a, NULL, NULL, divvext_a, nthreads);
    if(ierr!=0) return ierr+300;

    kernel_get_vector_vext<<<nblocks, nthreads>>>(it, SPINB, vecvext_b); // NOTE - only SPINB
    ierr=compute_divergence_real_vector_f(vecvext_b, NULL, NULL, divvext_b, nthreads);
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
                                            (Complex *)wf_d_dx, d_kky, d_kkz, (Complex *)wf_laplace,
                                            vecvext_a, divvext_a, vecvext_b, divvext_b
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
    ierr=compute_derivative_real_vector_f(potsall.A_a_x, NULL, NULL, grad_j_corr_a, NULL, NULL, nthreads);
    if(ierr!=0) return ierr;
    ierr=compute_derivative_real_vector_f(potsall.A_b_x, NULL, NULL, grad_j_corr_b, NULL, NULL, nthreads);
    if(ierr!=0) return ierr;
#endif

    // Step 3: apply hamiltonian
    kernel_apply_hamiltonian<<<nblocks, nthreads>>>(it, potsall,
                                            laplace_alpha_a, laplace_alpha_b,
                                            grad_j_corr_a, NULL, NULL, grad_j_corr_b, NULL, NULL,
                                            n, (Complex *)wf_in, (Complex *)wf_out,
                                            (Complex *)wf_d_dx, d_kky, d_kkz, (Complex *)wf_laplace, (Complex *)alphawf_laplace,
                                            cccoeff,
                                            vecvext_a, divvext_a, vecvext_b, divvext_b
                                                   );
#endif

    // Step 4: to increas stability - subtruct <H>*wf, where <H> is quasi particle energy
    // and normalize
    // use grad_alpha_a as working buffer
    double *gpe;

    if(useqpe==NULL) // compute quasiparticle energies
    {
        gpe=grad_alpha_a; // for easier notation
        
        // compute quasi-particle energy for each wave-function
        int noAllElements = n*NX;
        nblocks = (int)ceil((float)noAllElements/nthreads);
        kernel_compute_qpe<<<nblocks, nthreads>>>(  (Complex *)wf_in,               (Complex *)wf_out, 
                                                    (Complex *)wf_in+noAllElements, (Complex *)wf_out+noAllElements, 
                                                    gpe, noAllElements);

        ierr = local_reductions_many(n, NX, gpe, gpe);
        if(ierr!=0) return ierr;

    }
    else // use given values
    {
        gpe=useqpe;
    }

    // subtruct <H>*wf
    nblocks = (int)ceil((float)NX/nthreads);
    kernel_subtruct_qpe<<<nblocks, nthreads>>>(n, (Complex *)wf_in, (Complex *)wf_out, gpe);  
 
    return 0;
}

// ================================================================================================
// ========================================= compute_ovelap =======================================
// ================================================================================================
__global__ void kernel_compute_ovelap(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *re, double *im)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    if(ixyz<NX)
    {
        p=thrust::conj(wf1_u[ixyz])*wf2_u[ixyz] + thrust::conj(wf1_v[ixyz])*wf2_v[ixyz];
        re[ixyz]=p.real()*DX;
        im[ixyz]=p.imag()*DX;
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
 * @param workarea work space of size 2*NX, it can be the same workspace as used by cuFFT
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 * */
extern "C" int compute_ovelap(int n, cufftDoubleComplex *wf1, cufftDoubleComplex *wf2, double *overlap_re, double *overlap_im,
                              double *workarea, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);

    double *re = workarea;
    double *im = workarea + NX;
    int iwf;
    size_t shift=0;
    int ierr;

    for(iwf=0; iwf<n; iwf++) // for each wave-function
    {
        kernel_compute_ovelap<<<nblocks, nthreads>>>((Complex *)wf1+shift, (Complex *)wf2+shift, (Complex *)wf1+shift+n*NX, (Complex *)wf2+shift+n*NX, re, im);
        shift+=NX; // move pointer to next wf
        
        ierr = local_reductions_many(2, NX, workarea, workarea);
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
    if(ixyz<NX)
    {
        norm[ixyz]=(thrust::norm(wf_u[ixyz])+thrust::norm(wf_v[ixyz]))*DX;
    }
}
__global__ void kernel_compute_norm(Complex *wf_u, Complex *wf_v, double *norm, int noAllElements)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    if(ixyz<noAllElements)
    {
        norm[ixyz]=(thrust::norm(wf_u[ixyz])+thrust::norm(wf_v[ixyz]))*DX;
    }
}


__global__ void kernel_normalize_wf(int n, Complex *wf, double *norm)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double c;
    int iwf;
    if(ixyz<NX)
    {
        for(iwf=0; iwf<n; iwf++)
        {
            c=1.0/sqrt(norm[iwf]);

            // normalize wf
            wf[ixyz+iwf*NX     ]*=c;
            wf[ixyz+iwf*NX+n*NX]*=c;
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
    int ierr;
    int noAllElements = NX*n;
    double * norm  = (double *)pca_cufft_work_area; // Use here work area of cuFFT as working buffer

    // compute norm
    int nblocks = (int)ceil((float)noAllElements/nthreads);
    kernel_compute_norm<<<nblocks, nthreads>>>((Complex *)wf, (Complex *)wf+noAllElements, norm, noAllElements);

    // reduce wave-functions pararell
    ierr = local_reductions_many(n, NX, norm, norm);
    if(ierr != 0) return ierr;
    
    // normalize wf
    nblocks = (int)ceil((float)NX/nthreads);
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

    if(ixyz<NX)
    {
        na=alpha_a[ixyz];
        nb=alpha_b[ixyz];

        for(iwf=0; iwf<n; iwf++)
        {
            // read u and v
            u=wf_in[       iwf*NX+ixyz];
            v=wf_in[n*NX +iwf*NX+ixyz];

            wf_out[       iwf*NX+ixyz]=u*na; // multiply by effective mass and save
            wf_out[n*NX +iwf*NX+ixyz]=v*nb; // multiply by effective mass and save
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
    int nblocks = (int)ceil((float)NX/nthreads);

    wslda_potential potsall=convert_into_wslda_potential(d_potentials, NX, NULL);

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

    if(ixyz<NX)
    {

        for(iwf=0; iwf<n; iwf++)
        {
            // read u and v
            u=wf_hpsi[       iwf*NX+ixyz];
            v=wf_hpsi[n*NX +iwf*NX+ixyz];

            // add Taylor coefficients
            u*=Complex(0.0, t_coeff);
            v*=Complex(0.0, t_coeff);


            wf_contr[       iwf*NX+ixyz]=u; // save taylor contribution
            wf_contr[n*NX +iwf*NX+ixyz]=v; // save taylor contribution

            wf_update[       iwf*NX+ixyz]+=u; // add to rest of Taylor expansion
            wf_update[n*NX +iwf*NX+ixyz]+=v; // add to rest of Taylor expansion
        }
    }
}

/**
 * Function add it-order conytribution from Taylor expansion: (1/it!)(-i*(H-<H>)dt)*psi
 * @param it Taylor order contribution
 * @param dt time step
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf_hpsi result of apply_hamiltonian (INPUT)
 * @param wf_update array with wave-functions contructed by Taylor expansion (OUTPUT),
 * @param wf_contr array with wave-functions contructed by Taylor expansion (OUTPUT),
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 **/
extern "C" int taylor_expansion_contribution(int it, double dt, int n, cufftDoubleComplex *wf_hpsi,
                                             cufftDoubleComplex *wf_update, cufftDoubleComplex *wf_contr, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);

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
    int ix;

    if(ixyz<NX)
    {
        ix = ixyz; // decode cartesian coordinates
        data[ixyz]=u_ext(ix,0,0,it,spin);
    }
}

/**
 * @return array with external potential
 * */
extern "C" int get_v_ext(int datadim, int spin, int it, double *data)
{
    int nthreads = 256;
    int nblocks = (int)ceil((float)NX/nthreads);

    // allocate cuda memory
    double *wrkspace = (double *)pca_cufft_work_area; // reuse workspace
    kernel_get_v_ext<<<nblocks, nthreads>>>(it, spin, wrkspace);

    if( cudaMemcpy( data , wrkspace, sizeof(double)*NX, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;

    return 0;
}

__global__ void kernel_get_delta_ext(int it, Complex *deltain, Complex *data)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;

    if(ixyz<NX)
    {
        ix = ixyz; // decode cartesian coordinates
        data[ixyz]=macro_delta_ext(ix, 0, 0, it, deltain[ixyz]);
    }
}

/**
 * @return array with external potential
 * */
extern "C" int get_delta_ext(int datadim, int it, void *deltain, void *data)
{
    int nthreads = 256;
    int nblocks = (int)ceil((float)NX/nthreads);

    // allocate cuda memory
    Complex *wrkspace = (Complex *)pca_cufft_work_area; // reuse workspace
    kernel_get_delta_ext<<<nblocks, nthreads>>>(it, (Complex *)deltain, (Complex *)wrkspace);

    if( cudaMemcpy( data , wrkspace, sizeof(Complex)*NX, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;

    return 0;
}

__global__ void kernel_get_velocity_ext(int it, int spin, double *datax, double *datay, double *dataz)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;

    if(ixyz<NX)
    {
        ix = ixyz; // decode cartesian coordinates
        datax[ixyz]=velocity_ext(ix, 0, 0, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        datay[ixyz]=velocity_ext(ix, 0, 0, it, spin, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        dataz[ixyz]=velocity_ext(ix, 0, 0, it, spin, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data);
    }
}

/**
 * @return array with external potential
 * */
extern "C" int get_velocity_ext(int datadim, int spin, int it, double *data)
{
    int nthreads = 256;
    int nblocks = (int)ceil((float)NX/nthreads);

    double *wrkspace = (double *)pca_cufft_work_area; // reuse workspace
    double *vx = wrkspace + 0*NX;
    double *vy = wrkspace + 1*NX;
    double *vz = wrkspace + 2*NX;

    kernel_get_velocity_ext<<<nblocks, nthreads>>>(it, spin, vx, vy, vz);

    if( cudaMemcpy( data , wrkspace, sizeof(double)*NX*3, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;

    return 0;
}
