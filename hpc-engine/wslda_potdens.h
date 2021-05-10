/** 
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 * 
 * @author Gabriel Wlazlowski
 * @date 29.08.2020
 * */

#include <complex.h>

#ifndef _WSLDA_POTDENS_
#define _WSLDA_POTDENS_

// ordering arrays in h_densities, 
// c - complex, r- real
#define DENSCNT 12
#define DENSTYPE "crrrrrrrrrr"

// ordering arrays in h_potentials, 
// c - complex, r- real
#define POTCNT 12
#define POTTYPE "rrcrrrrrrrr"

#if CODEDIM==1

#define DENSDIM (DENSCNT*NX)
#define POTDIM  (POTCNT*NX)
#define BLOCKLENGTH (NX)

#elif CODEDIM==2

#define DENSDIM (DENSCNT*NX*NY)
#define POTDIM  (POTCNT*NX*NY)
#define BLOCKLENGTH (NX*NY)

#else

#define DENSDIM (DENSCNT*NXYZ)
#define POTDIM  (POTCNT*NXYZ)
#define BLOCKLENGTH (NXYZ) 

#endif


typedef struct
{
    int nx;             
    int ny;             /// for 1d code it is set to 1
    int nz;             /// for 1d and 2d code it is set to 1
    int datadim;        /// dimensonality of data
    int blocklength;    /// number of elements in potential array = nx*ny*nz
    
    // densities
    double complex *nu;
    double *rho_a;
    double *rho_b;
    double *tau_a;
    double *tau_b;
    double *j_a_x;
    double *j_a_y;
    double *j_a_z;
    double *j_b_x;
    double *j_b_y;
    double *j_b_z;  
} wslda_density; 

typedef struct
{
    int nx;             
    int ny;             /// for 1d code it is set to 1
    int nz;             /// for 1d and 2d code it is set to 1
    int datadim;        /// dimensonality of data
    int blocklength;    /// number of elements in potential array = nx*ny*nz
    
    // potentials
    double complex *delta;
    double *alpha_a; /// effective mass, spin-a
    double *alpha_b; /// effective mass, spin-b
    double *V_a;
    double *V_b;
    double *A_a_x;
    double *A_a_y;
    double *A_a_z;
    double *A_b_x;
    double *A_b_y;
    double *A_b_z;
    
    // chemical potentials
    double *mu;
} wslda_potential; 


/**
 * Converts array into wslda_density structure
 * */
wslda_density convert_into_wslda_density(double *h_densities, int blocklength);

/**
 * Converts array into wslda_potential structure
 * */
wslda_potential convert_into_wslda_potential(double *h_potentials, int blocklength, double *mu);

#endif
