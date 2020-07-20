// Author: Gabriel Wlazlowski

// Code for computing matrix elements of hamiltonian

#ifndef __S3DPCA_ME__
#define __S3DPCA_ME__

#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "s3dpca_fft.h"
#include "s3dpca_grid.h"

int recompute_potentials(int it, double *h_densities, double *h_potentials, double *h_potentials_new);
int recompute_potentials_meanfield_only(int it, double *h_densities, double *h_potentials, double *h_potentials_new);
int compute_matrix_elements(metadata_s3dpca_grid *bgrid, int it, double *h_densities, double *h_potentials, metadata_s3dpca_fft *mdfft, double complex *h, double complex * me_d_dx, double complex * me_d_dy, double complex * me_d_dz);
int compute_energy(int it, double *h_densities, double *h_potentials, double *energy, double *npart);
int compute_matrix_elements_of_momentum_operator(int nx, double dx, double complex *me);
int compute_angular_momentum_Lz(double *jx, double *jy, double *Lz);
int compute_vext_dot_j(int it, int spin, double *jx, double *jy, double *jz, double *vext_dot_j);
// int test_me(int nx);
#endif

