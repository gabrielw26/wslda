// Author: Gabriel Wlazlowski

// Code for computing densities
#include "wslda_potdens.h"

#ifndef __KZPCA_DENSITIES__
#define __KZPCA_DENSITIES__
int compute_contribution_to_densities(int nwf, double *En, double complex *psi, double ecut, double beta, wslda_density h_densities, 
                                      metadata_s1dpca_fft *mdfft, double ky, double kz, int weight, int spinsymmetry);
int density_caculate_tau(wslda_density h_densities, metadata_s1dpca_fft *mdfft);
#endif
