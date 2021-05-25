// Author: Gabriel Wlazlowski

// Code for computing densities

#include <stdlib.h>
#include <stddef.h>
#include "pca_settings.h"
#include "pca_macro.h"
#include "s1dpca_fft.h"
#include "s1dpca_densities.h"


double fbeta(double E, double beta);

/**
 * Function computes contribution to the densities
 * @param nwf number of states for kz that contribute to the density (INPUT)
 * @param En eigen energies (INPUT)
 * @param psi eigenvectors (INPUT)
 * @param ecut eenergy cut-off (INPUT)
 * @param beta inverse of temperature
 * @param h_densities array with densities to be updated (INPUT/OUTPUT)
 * @param mdfft metadata for ffts plans execution (INPUT)
 * @param ky value of ky (INPUT) 
 * @param kz value of kz (INPUT) 
 * @param weight degeneracy of state (INPUT)
 * @param spinsymmetry (INPUT)
 * @param S buffer for entopy (INPUT/OUTPUT)
 * */
int compute_contribution_to_densities(int nwf, double *En, double complex *psi, double ecut, double beta, wslda_density h_densities, 
                                      metadata_s1dpca_fft *mdfft, double ky, double kz, int weight, int spinsymmetry, double *S)
{
    int ien; 
    int ix, iy, iz, ixyz;
    double complex *u, *v; // psi=(u,v)
    double kx, k2;
    double fbEn, fbmEn;
    double complex wfdx, wfdy, wfdz;
    double complex wfdxv, wfdyv, wfdzv;
    double complex wfdxu, wfdyu, wfdzu;
    double norm;
    
    // densities - decode 
    double *rho_a = h_densities.rho_a;
    double *rho_b = h_densities.rho_b;
    double *tau_a = h_densities.tau_a;
    double *tau_b = h_densities.tau_b;
    double complex *nu = h_densities.nu;
    double *j_a_x = h_densities.j_a_x;
    double *j_a_y = h_densities.j_a_y;
    double *j_a_z = h_densities.j_a_z;
    double *j_b_x = h_densities.j_b_x;
    double *j_b_y = h_densities.j_b_y;
    double *j_b_z = h_densities.j_b_z;
    
    double _kz=2.*M_PI/( double )LZ * ( double )(NZ/2-NZ); // momentum for which I should kill contribution for gradients
    if(fabs(_kz-kz)<1.0e-12) _kz = 0.0;
    else                     _kz = kz;
    
    double _ky=2.*M_PI/( double )LY * ( double )(NY/2-NY); // momentum for which I should kill contribution for gradients
    if(fabs(_ky-ky)<1.0e-12) _ky = 0.0;
    else                     _ky = ky;
        
    for(ien=0; ien<nwf; ien++) // for each eigen-energy 
    {
        
        if(fabs(En[ien])>ecut) continue; // skip states above the cut-off energy
        if(spinsymmetry>0 && En[ien]<0.0) continue; // take only positive states
        
        // docompose state
        u = psi + ien*2*NX;
        v = u + NX;
        
        // normalize eigen-vectors
        #define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a))
        norm=0.0;
        for(ixyz=0; ixyz<NX; ixyz++) norm+=cnorm(u[ixyz])+cnorm(v[ixyz]); // \int (|u|^2 + |v|^2)d^3r
        norm*=DX; // volume element
        norm=1./sqrt(norm); // normalization factor
        for(ixyz=0; ixyz<NX; ixyz++) u[ixyz]*=norm;
        for(ixyz=0; ixyz<NX; ixyz++) v[ixyz]*=norm;
        
        // compute gradients
        for(ixyz=0; ixyz<2*NX; ixyz++) mdfft->fft1uv[ixyz]=u[ixyz];
        fftw_execute(mdfft->plan_f_uv);
        
        ixyz=0;
        for(ix=0; ix<NX; ix++) 
        {
            // extract momentum
            if(ix<NX/2) kx=2.*M_PI/(( double )NX*DX) * ( double )(ix   ) / (NX);
            else        kx=2.*M_PI/(( double )NX*DX) * ( double )(ix-NX) / (NX);
                        
            k2 = -1.0*(kx*kx)*NX;
            
            // corrections for gradient computation
            if(ix==NX/2) kx=0.0;
            
            mdfft->fft1grad[ixyz     ] = mdfft->fft1uv[ixyz   ] * I * kx; // du/dx
            mdfft->fft1grad[ixyz+1*NX] = mdfft->fft1uv[ixyz+NX] * I * kx; // dv/dx
            
            // laplace
            mdfft->fft1uv[ixyz   ]*=k2;
            mdfft->fft1uv[ixyz+NX]*=k2;
            
            ixyz++;
        }
        
        fftw_execute(mdfft->plan_b_grad);
        fftw_execute(mdfft->plan_b_uv);
        
        // form densities 
        
        // weight
        #define DENS_FACTOR_M 1.
        #define Complex(a,b) (a + I*b)
        #define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a))
        
        fbEn=fbeta(En[ien], beta)*DENS_FACTOR_M;
        fbmEn = DENS_FACTOR_M - fbEn;
        
        // entropy 
        if(spinsymmetry>0)
        {
            if(fbEn >NUMERICAL_ZERO) S[0] -=  fbEn * log( fbEn) * weight * 2.0;
            if(fbmEn>NUMERICAL_ZERO) S[0] -= fbmEn * log(fbmEn) * weight * 2.0;
        }
        else
        {
            if(fbEn >NUMERICAL_ZERO) S[0] -=  fbEn * log( fbEn) * weight;
            if(fbmEn>NUMERICAL_ZERO) S[0] -= fbmEn * log(fbmEn) * weight;
        }

        if(spinsymmetry>0) for(ixyz=0; ixyz<NX; ixyz++)
        {
            // form na, nb, nu
            // rho_a[ixyz]+=... will be taken as rho_b later
            rho_b[ixyz]+=(cnorm(v[ixyz])*fbmEn + cnorm(u[ixyz])*fbEn) / LY / LZ * weight;
            nu[ixyz]+=u[ixyz]*conj(v[ixyz])*(fbmEn-fbEn) / LY / LZ * weight;     
            
            // form taua and j_a
            // skip, will be taken from b component later
            
            // form taub and j_b
            wfdxu = mdfft->fft1grad[ixyz     ]; // du/dx
            wfdyu = u[ixyz] * I * _ky; // du/dy
            wfdzu = u[ixyz] * I * _kz; // du/dz
            wfdxv = mdfft->fft1grad[ixyz+1*NX]; // dv/dx
            wfdyv = v[ixyz] * I * _ky; // dv/dy
            wfdzv = v[ixyz] * I * _kz; // dv/dz
#ifdef TAU_COMPUTATION_VIA_GRADIENTS
            tau_b[ixyz]+=(cnorm(wfdxv)+cnorm(wfdyv)+cnorm(wfdzv))*fbmEn / LY / LZ * weight 
                       + (cnorm(wfdxu)+cnorm(wfdyu)+cnorm(wfdzu))*fbEn  / LY / LZ * weight;
#else
            tau_b[ixyz]+=creal( conj(v[ixyz])*(mdfft->fft1uv[ixyz+NX] - v[ixyz]*ky*ky - v[ixyz]*kz*kz)*fbmEn / LY / LZ * weight 
                       + conj(u[ixyz])*(mdfft->fft1uv[ixyz   ] - u[ixyz]*ky*ky - u[ixyz]*kz*kz)*fbEn  / LY / LZ * weight );
#endif
            j_b_x[ixyz]-=cimag(conj(v[ixyz])*wfdxv)*fbmEn / LY / LZ * weight - cimag(conj(u[ixyz])*wfdxu)*fbEn / LY / LZ * weight;
            // j_b_y[ixyz]+=0.0;
            // j_b_z[ixyz]+=0.0;
        }
        else for(ixyz=0; ixyz<NX; ixyz++)
        {
            // form na, nb, nu
            rho_a[ixyz]+=cnorm(u[ixyz])*fbEn / LY / LZ * weight;
            rho_b[ixyz]+=cnorm(v[ixyz])*fbmEn / LY / LZ * weight;
            nu[ixyz]+=u[ixyz]*conj(v[ixyz])*(fbmEn-fbEn)/2.0 / LY / LZ * weight;     
            
            // form taua and j_a
            wfdx = mdfft->fft1grad[ixyz     ]; // du/dx
            wfdy = u[ixyz] * I * _ky; // du/dy
            wfdz = u[ixyz] * I * _kz; // du/dz
#ifdef TAU_COMPUTATION_VIA_GRADIENTS
            tau_a[ixyz]+=(cnorm(wfdx)+cnorm(wfdy)+cnorm(wfdz))*fbEn / LY / LZ * weight;
#else
            tau_a[ixyz]+=creal( conj(u[ixyz])*(mdfft->fft1uv[ixyz     ] - u[ixyz]*ky*ky - u[ixyz]*kz*kz)*fbEn / LY / LZ * weight );
#endif
            j_a_x[ixyz]+=cimag(conj(u[ixyz])*wfdx)*fbEn / LY / LZ * weight;
            // j_a_y[ixyz]+=0.0;
            // j_a_z[ixyz]+=0.0;
            
            // form taub and j_b
            wfdx = mdfft->fft1grad[ixyz+1*NX]; // dv/dx
            wfdy = v[ixyz] * I * _ky; // dv/dy
            wfdz = v[ixyz] * I * _kz; // dv/dz
#ifdef TAU_COMPUTATION_VIA_GRADIENTS
            tau_b[ixyz]+=(cnorm(wfdx)+cnorm(wfdy)+cnorm(wfdz))*fbmEn / LY / LZ * weight;
#else
            tau_b[ixyz]+=creal( conj(v[ixyz])*(mdfft->fft1uv[ixyz+NX] - v[ixyz]*ky*ky - v[ixyz]*kz*kz)*fbmEn / LY / LZ * weight );
#endif
            j_b_x[ixyz]-=cimag(conj(v[ixyz])*wfdx)*fbmEn / LY / LZ * weight;
            // j_b_y[ixyz]+=0.0;
            // j_b_z[ixyz]+=0.0;
        }
        
    } // for(ien=0; ien<2*NX; ien++)
    
//     double ttt=0.0;
//     for(ixyz=0; ixyz<NX; ixyz++) ttt+=rho_a[ixyz] * LZ;
//     wprintf("ttt=%f\n", ttt);
    
    return 0;
}

/**
 * Function computes tau according formula:
 * tau = (1/2)laplace(rho) - Re( sum_n Psi_n^* laplapce (Psi_n) )
 * See Eq.(28) in PHYSICAL REVIEW C 95, 044302 (2017)
 * @param h_densities array with densities (INPUT: `tau` array keeps only  Re( sum_n Psi_n^* laplapce (Psi_n) ), OUTPUT: kinetic energy density)
 * @param mdfft handler for fftw plans
 * @return 0 - OK, otherwise ERROR 
 */
int density_caculate_tau(wslda_density h_densities, metadata_s1dpca_fft *mdfft)
{
    // densities - decode 
    double *rho_a = h_densities.rho_a;
    double *rho_b = h_densities.rho_b;
    double *tau_a = h_densities.tau_a;
    double *tau_b = h_densities.tau_b;
    double complex *nu = h_densities.nu;
    double *j_a_x = h_densities.j_a_x;
    double *j_a_y = h_densities.j_a_y;
    double *j_a_z = h_densities.j_a_z;
    double *j_b_x = h_densities.j_b_x;
    double *j_b_y = h_densities.j_b_y;
    double *j_b_z = h_densities.j_b_z;
    
    int ierr, ixyz;
    double *laplace_rho;
    cppmallocl(laplace_rho,NX,double);
    
    // ------------- for a component -------------
    ierr = compute_laplace_real_f(rho_a, laplace_rho, mdfft);
    if(ierr!=0) return ierr - 101;
    for(ixyz=0; ixyz<NX; ixyz++)
    {
        tau_a[ixyz] = 0.5*laplace_rho[ixyz] - tau_a[ixyz];
        // tau cannot be negative
        if(tau_a[ixyz]<0.0) tau_a[ixyz]=0.0;    
    }
    
    // ------------- for b component -------------
    ierr = compute_laplace_real_f(rho_b, laplace_rho, mdfft);
    if(ierr!=0) return ierr - 102;
    for(ixyz=0; ixyz<NX; ixyz++)
    {
        tau_b[ixyz] = 0.5*laplace_rho[ixyz] - tau_b[ixyz];
        // tau cannot be negative
        if(tau_b[ixyz]<0.0) tau_b[ixyz]=0.0;    
    }
    
    // clear memory
    free(laplace_rho);
    
    return 0;
}
