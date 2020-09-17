// Author: Gabriel Wlazlowski

// Code for computing matrix elements of hamiltonian
#include "s3dpca_grid.h"
#include "wslda_potdens.h"

#ifndef __KZPCA_ME__
#define __KZPCA_ME__

int compute_matrix_elements_2d(metadata_s3dpca_grid *bgrid, int it, wslda_density h_densities, wslda_potential h_potentials, metadata_s2dpca_fft *mdfft, double complex *h, double kz, double complex * me_d_dx, double complex * me_d_dy);
int compute_matrix_elements_of_momentum_operator(int nx, double dx, double complex *me);
int compute_angular_momentum_Lz(double *jx, double *jy, double *Lz);
int compute_vext_dot_j(int it, int spin, double *jx, double *jy, double *vext_dot_j);
// int test_me(int nx);

#endif

