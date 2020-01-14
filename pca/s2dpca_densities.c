// Author: Gabriel Wlazlowski

// Code for computing densities

#include <stdlib.h>
#include <stddef.h>
#include "pca_settings.h"
#include "pca_macro.h"
#include "s2dpca_fft.h"
#include "s2dpca_densities.h"

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
 * @param kz value of kz (INPUT)  
 * */
int compute_contribution_to_densities(int nwf, double *En, double complex *psi, double ecut, double beta, double *h_densities, 
                                      metadata_s2dpca_fft *mdfft, double kz)
{
    int ien; 
    int ix, iy, iz, ixyz;
    double complex *u, *v; // psi=(u,v)
    double kx, ky;
    double fbEn, fbmEn;
    double complex wfdx, wfdy, wfdz;
    double norm;
    
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    double _kz=2.*M_PI/( double )LZ * ( double )(NZ/2-NZ); // momentum for which I should kill contribution for gradients
    if(fabs(_kz-kz)<1.0e-12) _kz = 0.0;
    else                     _kz = kz;
    
    double weight=2.0; // all kz vectors are in +/- pairs, I compute for only single component of pair and double contribution
    if(kz==0.0) weight=1.0; // except for kz=0 compoment, which has no pair
    
    for(ien=0; ien<nwf; ien++) // for each eigen-energy 
    {
        
        if(fabs(En[ien])>ecut) continue; // skip states above the cut-off energy
        
        // docompose state
        u = psi + ien*2*NX*NY;
        v = u + NX*NY;
        
        // normalize eigen-vectors
        #define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a))
        norm=0.0;
        for(ixyz=0; ixyz<NX*NY; ixyz++) norm+=cnorm(u[ixyz])+cnorm(v[ixyz]); // \int (|u|^2 + |v|^2)d^3r
        norm*=DX*DY; // volume element
        norm=1./sqrt(norm); // normalization factor
        for(ixyz=0; ixyz<NX*NY; ixyz++) u[ixyz]*=norm;
        for(ixyz=0; ixyz<NX*NY; ixyz++) v[ixyz]*=norm;
        
        // compute gradients
        for(ixyz=0; ixyz<2*NX*NY; ixyz++) mdfft->fft2uv[ixyz]=u[ixyz];
        fftw_execute(mdfft->plan_f_uv);
        
        ixyz=0;
        for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) 
        {
            // extract momentum
            if(ix<NX/2) kx=2.*M_PI/(( double )NX*DX) * ( double )(ix   ) / (NX*NY);
            else        kx=2.*M_PI/(( double )NX*DX) * ( double )(ix-NX) / (NX*NY);
            
            if(iy<NY/2) ky=2.*M_PI/(( double )NY*DY) * ( double )(iy   ) / (NX*NY);
            else        ky=2.*M_PI/(( double )NY*DY) * ( double )(iy-NY) / (NX*NY);
            
            // corrections for gradient computation
            if(ix==NX/2) kx=0.0;
            if(iy==NY/2) ky=0.0;        
            
            mdfft->fft2grad[ixyz        ] = mdfft->fft2uv[ixyz      ] * I * kx; // du/dx
            mdfft->fft2grad[ixyz+1*NX*NY] = mdfft->fft2uv[ixyz+NX*NY] * I * kx; // dv/dx
            mdfft->fft2grad[ixyz+2*NX*NY] = mdfft->fft2uv[ixyz      ] * I * ky; // du/dy
            mdfft->fft2grad[ixyz+3*NX*NY] = mdfft->fft2uv[ixyz+NX*NY] * I * ky; // dv/dy
            
            ixyz++;
        }
        
        fftw_execute(mdfft->plan_b_grad);
        
        // form densities 
        
        // weight
        #define DENS_FACTOR_M 1.
        #define Complex(a,b) (a + I*b)
        #define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a))
        
        fbEn=fbeta(En[ien], beta)*DENS_FACTOR_M;
        fbmEn = DENS_FACTOR_M - fbEn;
            
        for(ixyz=0; ixyz<NX*NY; ixyz++)
        {
            // form na, nb, nu
            rho_a[ixyz]+=cnorm(u[ixyz])*fbEn / LZ * weight;
            rho_b[ixyz]+=cnorm(v[ixyz])*fbmEn / LZ * weight;
            nu[ixyz]+=u[ixyz]*conj(v[ixyz])*(fbmEn-fbEn)/2.0 / LZ * weight;     
            
            // form taua and j_a
            wfdx = mdfft->fft2grad[ixyz        ]; // du/dx
            wfdy = mdfft->fft2grad[ixyz+2*NX*NY]; // du/dy
            wfdz = u[ixyz] * I * _kz; // du/dz
            tau_a[ixyz]+=(cnorm(wfdx)+cnorm(wfdy)+cnorm(wfdz))*fbEn / LZ * weight;
            j_a_x[ixyz]+=cimag(conj(u[ixyz])*wfdx)*fbEn / LZ * weight;
            j_a_y[ixyz]+=cimag(conj(u[ixyz])*wfdy)*fbEn / LZ * weight;
            // j_a_z[ixyz]+=0.0;
            
            // form taub and j_b
            wfdx = mdfft->fft2grad[ixyz+1*NX*NY]; // dv/dx
            wfdy = mdfft->fft2grad[ixyz+3*NX*NY]; // dv/dy
            wfdz = v[ixyz] * I * _kz; // dv/dz
            tau_b[ixyz]+=(cnorm(wfdx)+cnorm(wfdy)+cnorm(wfdz))*fbmEn / LZ * weight;
            j_b_x[ixyz]-=cimag(conj(v[ixyz])*wfdx)*fbmEn / LZ * weight;
            j_b_y[ixyz]-=cimag(conj(v[ixyz])*wfdy)*fbmEn / LZ * weight;
            // j_b_z[ixyz]+=0.0;
        }
        
    } // for(ien=0; ien<2*NX*NY; ien++)
    
//     double ttt=0.0;
//     for(ixyz=0; ixyz<NX*NY; ixyz++) ttt+=rho_a[ixyz] * LZ;
//     printf("ttt=%f\n", ttt);
    
    return 0;
}
