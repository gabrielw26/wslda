// Author: Gabriel Wlazlowski

// Code for computing densities
#include "wslda_potdens.h"

#ifndef __KZPCA_DENSITIES__
#define __KZPCA_DENSITIES__
int compute_contribution_to_densities(double *En, double complex *psi, int nwfip, double beta, wslda_density h_densities, 
                                      metadata_s3dpca_fft *mdfft, int spinsymmetry);
int density_caculate_tau(wslda_density h_densities, metadata_s3dpca_fft *mdfft);
#endif
