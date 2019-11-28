// Author: Gabriel Wlazlowski

// Code for computing densities

#ifndef __KZPCA_DENSITIES__
#define __KZPCA_DENSITIES__
int compute_contribution_to_densities(int nwf, double *En, double complex *psi, double ecut, double beta, double *h_densities, 
                                      metadata_kzpca_fft *mdfft, double kz);
#endif
