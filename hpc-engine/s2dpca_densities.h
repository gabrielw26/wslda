// Author: Gabriel Wlazlowski

// Code for computing densities

#ifndef __KZPCA_DENSITIES__
#define __KZPCA_DENSITIES__
int compute_contribution_to_densities(int nwf, double *En, double complex *psi, double ecut, double beta, double *h_densities, 
                                      metadata_s2dpca_fft *mdfft, double kz, int spinsymmetry);
int density_caculate_tau(double *h_densities, metadata_s2dpca_fft *mdfft);
#endif
