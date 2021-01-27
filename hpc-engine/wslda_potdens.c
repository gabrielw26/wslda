/** 
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 * 
 * @author Gabriel Wlazlowski
 * @date 29.08.2020
 * */

#include "pca_settings.h"
#include "wslda_potdens.h"


wslda_density convert_into_wslda_density(double *h_densities, int blocklength)
{
    wslda_density d;
    
    if(blocklength==NX      ) {d.nx=NX; d.ny=1 ; d.nz=1 ; d.datadim=1;} // 1D code
    if(blocklength==NX*NY   ) {d.nx=NX; d.ny=NY; d.nz=1 ; d.datadim=2;} // 2D code
    if(blocklength==NX*NY*NZ) {d.nx=NX; d.ny=NY; d.nz=NZ; d.datadim=3;} // 3D code
    d.blocklength = d.nx*d.ny*d.nz;
    
    d.nu = (double complex *)(h_densities +  0*blocklength);
    
    d.rho_a = (double *)(h_densities +  2*blocklength);
    d.tau_a = (double *)(h_densities +  3*blocklength);
    d.j_a_x = (double *)(h_densities +  4*blocklength);
    d.j_a_y = (double *)(h_densities +  5*blocklength);
    d.j_a_z = (double *)(h_densities +  6*blocklength);
    
    d.rho_b = (double *)(h_densities +  7*blocklength);
    d.tau_b = (double *)(h_densities +  8*blocklength);
    d.j_b_x = (double *)(h_densities +  9*blocklength);
    d.j_b_y = (double *)(h_densities + 10*blocklength);
    d.j_b_z = (double *)(h_densities + 11*blocklength);
    
    return d;
}

/**
 * Converts array into wslda_potential structure
 * */
wslda_potential convert_into_wslda_potential(double *h_potentials, int blocklength, double *mu)
{
    wslda_potential d;
    
    if(blocklength==NX      ) {d.nx=NX; d.ny=1 ; d.nz=1 ; d.datadim=1;} // 1D code
    if(blocklength==NX*NY   ) {d.nx=NX; d.ny=NY; d.nz=1 ; d.datadim=2;} // 2D code
    if(blocklength==NX*NY*NZ) {d.nx=NX; d.ny=NY; d.nz=NZ; d.datadim=3;} // 3D code
    d.blocklength = d.nx*d.ny*d.nz;
        
    // DO NOT CHANGE ORDER - IT WILL BREAK TD CODES!
    d.V_a = (double *)(h_potentials +  0*blocklength);
    d.V_b = (double *)(h_potentials +  1*blocklength);
    d.delta = (double complex *)(h_potentials +  2*blocklength);
    
    d.alpha_a = (double *)(h_potentials +  4*blocklength);
    d.alpha_b = (double *)(h_potentials +  5*blocklength);
    
    d.A_a_x = (double *)(h_potentials +  6*blocklength);
    d.A_a_y = (double *)(h_potentials +  7*blocklength);
    d.A_a_z = (double *)(h_potentials +  8*blocklength);
    d.A_b_x = (double *)(h_potentials +  9*blocklength);
    d.A_b_y = (double *)(h_potentials + 10*blocklength);
    d.A_b_z = (double *)(h_potentials + 11*blocklength);
    
    d.mu=mu;
    
    return d;
}
