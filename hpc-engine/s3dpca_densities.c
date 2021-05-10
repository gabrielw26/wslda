// Author: Gabriel Wlazlowski

// Code for computing densities

#include <stdlib.h>
#include <stddef.h>
#include "pca_settings.h"
#include "pca_macro.h"
#include "s3dpca_fft.h"
#include "s3dpca_densities.h"

double fbeta(double E, double beta);

/**
 * Function computes contribution to the densities
 * @param En eigen energies (INPUT)
 * @param psi eigenvectors (INPUT)
 * @param nwfip number of local wavefunctions (INPUT)
 * @param beta inverse of temperature
 * @param h_densities array with densities to be updated (INPUT/OUTPUT)
 * @param mdfft metadata for ffts plans execution (INPUT) 
 * @param spinsymmetry (INPUT)
 * @param S buffer for entopy (INPUT/OUTPUT)
 * */
int compute_contribution_to_densities(double *En, double complex *psi, int nwfip, double beta, wslda_density h_densities, 
                                      metadata_s3dpca_fft *mdfft, int spinsymmetry, double *S)
{
    int ien; 
    int ix, iy, iz, ixyz;
    double complex *u, *v; // psi=(u,v)
    double kx, ky, kz, k2;
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
    
    
    for(ien=0; ien<nwfip; ien++) // for each eigen-energy 
    {
        // docompose state
        u = psi + ien*2*NXYZ;
        v = u + NXYZ;
        
        // normalize eigen-vectors
        #define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a))
        norm=0.0;
        for(ixyz=0; ixyz<NXYZ; ixyz++) norm+=cnorm(u[ixyz])+cnorm(v[ixyz]); // \int (|u|^2 + |v|^2)d^3r
        norm*=DXYZ; // volume element
        norm=1./sqrt(norm); // normalization factor
        for(ixyz=0; ixyz<NXYZ; ixyz++) u[ixyz]*=norm;
        for(ixyz=0; ixyz<NXYZ; ixyz++) v[ixyz]*=norm;
        
            
            
        // compute gradients
        for(ixyz=0; ixyz<2*NXYZ; ixyz++) mdfft->fft3uv[ixyz]=u[ixyz];
        fftw_execute(mdfft->plan_f_uv);
        
        ixyz=0;
        for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++) 
        {
            // extract momentum
            if(ix<NX/2) kx=2.*M_PI/(( double )NX*DX) * ( double )(ix   ) / (NXYZ);
            else        kx=2.*M_PI/(( double )NX*DX) * ( double )(ix-NX) / (NXYZ);
            
            if(iy<NY/2) ky=2.*M_PI/(( double )NY*DY) * ( double )(iy   ) / (NXYZ);
            else        ky=2.*M_PI/(( double )NY*DY) * ( double )(iy-NY) / (NXYZ);
            
            if(iz<NZ/2) kz=2.*M_PI/(( double )NZ*DZ) * ( double )(iz   ) / (NXYZ);
            else        kz=2.*M_PI/(( double )NZ*DZ) * ( double )(iz-NZ) / (NXYZ);
            
            k2 = -1.0*(kx*kx + ky*ky + kz*kz)*NXYZ;
            
            // corrections for gradient computation
            if(ix==NX/2) kx=0.0;
            if(iy==NY/2) ky=0.0;    
            if(iz==NZ/2) kz=0.0;
            
            // gradients
            mdfft->fft3grad[ixyz       ] = mdfft->fft3uv[ixyz     ] * I * kx; // du/dx
            mdfft->fft3grad[ixyz+1*NXYZ] = mdfft->fft3uv[ixyz+NXYZ] * I * kx; // dv/dx
            mdfft->fft3grad[ixyz+2*NXYZ] = mdfft->fft3uv[ixyz     ] * I * ky; // du/dy
            mdfft->fft3grad[ixyz+3*NXYZ] = mdfft->fft3uv[ixyz+NXYZ] * I * ky; // dv/dy
            mdfft->fft3grad[ixyz+4*NXYZ] = mdfft->fft3uv[ixyz     ] * I * kz; // du/dz
            mdfft->fft3grad[ixyz+5*NXYZ] = mdfft->fft3uv[ixyz+NXYZ] * I * kz; // dv/dz          
            
            // laplace
            mdfft->fft3uv[ixyz     ]*=k2;
            mdfft->fft3uv[ixyz+NXYZ]*=k2;
            
            ixyz++;
        }
        
        fftw_execute(mdfft->plan_b_grad);
        fftw_execute(mdfft->plan_b_uv);
        
        // form densities 
        
        // weight
        #define DENS_FACTOR_M 1.
        #define Complex(a,b) (a + I*b)
        
        fbEn=fbeta(En[ien], beta)*DENS_FACTOR_M;
        fbmEn = DENS_FACTOR_M - fbEn;
        
        // entropy 
        if(spinsymmetry>0)
        {
            if(fbEn >NUMERICAL_ZERO) S[0] -=  fbEn * log( fbEn)*2.0;
            if(fbmEn>NUMERICAL_ZERO) S[0] -= fbmEn * log(fbmEn)*2.0;
        }
        else
        {
            if(fbEn >NUMERICAL_ZERO) S[0] -=  fbEn * log( fbEn);
            if(fbmEn>NUMERICAL_ZERO) S[0] -= fbmEn * log(fbmEn);
        }
            
        if(spinsymmetry>0) for(ixyz=0; ixyz<NXYZ; ixyz++) // spin symmetric case
        {
            // form na, nb, nu
            // rho_a[ixyz]+=... will be taken as rho_b later
            rho_b[ixyz]+=(cnorm(v[ixyz])*fbmEn + cnorm(u[ixyz])*fbEn);
            nu[ixyz]+=u[ixyz]*conj(v[ixyz])*(fbmEn-fbEn);     
            
            // form taua and j_a
            // skip, will be taken from b component later
            
            // form taub and j_b
            wfdxu = mdfft->fft3grad[ixyz       ]; // du/dx
            wfdyu = mdfft->fft3grad[ixyz+2*NXYZ]; // du/dy
            wfdzu = mdfft->fft3grad[ixyz+4*NXYZ]; // du/dz
            wfdxv = mdfft->fft3grad[ixyz+1*NXYZ]; // dv/dx
            wfdyv = mdfft->fft3grad[ixyz+3*NXYZ]; // dv/dy
            wfdzv = mdfft->fft3grad[ixyz+5*NXYZ]; // dv/dz
#ifdef TAU_COMPUTATION_VIA_GRADIENTS
            tau_b[ixyz]+=(cnorm(wfdxv)+cnorm(wfdyv)+cnorm(wfdzv))*fbmEn + (cnorm(wfdxu)+cnorm(wfdyu)+cnorm(wfdzu))*fbEn;
#else
            tau_b[ixyz]+=creal( conj(v[ixyz])*mdfft->fft3uv[ixyz+NXYZ]*fbmEn + conj(u[ixyz])*mdfft->fft3uv[ixyz     ]*fbEn );
#endif
            j_b_x[ixyz]-=cimag(conj(v[ixyz])*wfdxv)*fbmEn - cimag(conj(u[ixyz])*wfdxu)*fbEn;
            j_b_y[ixyz]-=cimag(conj(v[ixyz])*wfdyv)*fbmEn - cimag(conj(u[ixyz])*wfdyu)*fbEn;
            j_b_z[ixyz]-=cimag(conj(v[ixyz])*wfdzv)*fbmEn - cimag(conj(u[ixyz])*wfdzu)*fbEn;
            
        }
        else for(ixyz=0; ixyz<NXYZ; ixyz++) // spin-imbalanced case
        {
            // form na, nb, nu
            rho_a[ixyz]+=cnorm(u[ixyz])*fbEn;
            rho_b[ixyz]+=cnorm(v[ixyz])*fbmEn;
            nu[ixyz]+=u[ixyz]*conj(v[ixyz])*(fbmEn-fbEn)/2.0;     
            
            // form taua and j_a
            wfdx = mdfft->fft3grad[ixyz       ]; // du/dx
            wfdy = mdfft->fft3grad[ixyz+2*NXYZ]; // du/dy
            wfdz = mdfft->fft3grad[ixyz+4*NXYZ]; // du/dz
#ifdef TAU_COMPUTATION_VIA_GRADIENTS
            tau_a[ixyz]+=(cnorm(wfdx)+cnorm(wfdy)+cnorm(wfdz))*fbEn;
#else
            tau_a[ixyz]+=creal( conj(u[ixyz])*mdfft->fft3uv[ixyz     ]*fbEn );
#endif
            j_a_x[ixyz]+=cimag(conj(u[ixyz])*wfdx)*fbEn;
            j_a_y[ixyz]+=cimag(conj(u[ixyz])*wfdy)*fbEn;
            j_a_z[ixyz]+=cimag(conj(u[ixyz])*wfdz)*fbEn;
            
            // form taub and j_b
            wfdx = mdfft->fft3grad[ixyz+1*NXYZ]; // dv/dx
            wfdy = mdfft->fft3grad[ixyz+3*NXYZ]; // dv/dy
            wfdz = mdfft->fft3grad[ixyz+5*NXYZ]; // dv/dz
#ifdef TAU_COMPUTATION_VIA_GRADIENTS
            tau_b[ixyz]+=(cnorm(wfdx)+cnorm(wfdy)+cnorm(wfdz))*fbmEn;
#else
            tau_b[ixyz]+=creal( conj(v[ixyz])*mdfft->fft3uv[ixyz+NXYZ]*fbmEn );
#endif
            j_b_x[ixyz]-=cimag(conj(v[ixyz])*wfdx)*fbmEn;
            j_b_y[ixyz]-=cimag(conj(v[ixyz])*wfdy)*fbmEn;
            j_b_z[ixyz]-=cimag(conj(v[ixyz])*wfdz)*fbmEn;
            
        }
        
    } // for(ien=0; ien<nwfip; ien++)
    
    
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
int density_caculate_tau(wslda_density h_densities, metadata_s3dpca_fft *mdfft)
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
    cppmallocl(laplace_rho,NXYZ,double);
    
    // ------------- for a component -------------
    ierr = compute_laplace_real_f(rho_a, laplace_rho, mdfft);
    if(ierr!=0) return ierr - 101;
    for(ixyz=0; ixyz<NXYZ; ixyz++)
    {
        tau_a[ixyz] = 0.5*laplace_rho[ixyz] - tau_a[ixyz];
        // tau cannot be negative
        if(tau_a[ixyz]<0.0) tau_a[ixyz]=0.0;    
    }
    
    // ------------- for b component -------------
    ierr = compute_laplace_real_f(rho_b, laplace_rho, mdfft);
    if(ierr!=0) return ierr - 102;
    for(ixyz=0; ixyz<NXYZ; ixyz++)
    {
        tau_b[ixyz] = 0.5*laplace_rho[ixyz] - tau_b[ixyz];
        // tau cannot be negative
        if(tau_b[ixyz]<0.0) tau_b[ixyz]=0.0;    
    }
    
    // clear memory
    free(laplace_rho);
    
    return 0;
}
