// Author: Gabriel Wlazlowski

// Code for computing matrix elements of hamiltonian
#include "s3dpca_grid.h"
#include "wslda_potdens.h"

#ifndef __KZPCA_ME__
#define __KZPCA_ME__

int compute_matrix_elements_1d(metadata_s3dpca_grid *bgrid, int it, wslda_density h_densities, wslda_potential h_potentials, metadata_s1dpca_fft *mdfft, double complex *h, double ky, double kz, double complex * me_d_dx);
int compute_matrix_elements_of_momentum_operator(int nx, double dx, double complex *me);
int compute_vext_dot_j(int it, int spin, double *jx, double *vext_dot_j);
// int test_me(int nx);

#endif

