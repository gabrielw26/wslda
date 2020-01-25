// Author: Gabriel Wlazlowski

// Code for computing densities

#ifndef __KZPCA_DENSITIES__
#define __KZPCA_DENSITIES__
int compute_contribution_to_densities(double *En, double complex *psi, int nwfip, double beta, double *h_densities, 
                                      metadata_s3dpca_fft *mdfft, int spinsymmetry);
int density_caculate_tau(double *h_densities, metadata_s3dpca_fft *mdfft);
#endif
