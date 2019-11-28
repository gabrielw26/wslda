// Author: Gabriel Wlazlowski

// Code for computing matrix elements of hamiltonian

#ifndef __KZPCA_ME__
#define __KZPCA_ME__

void process_params(double *params, double kF);
int recompute_potentials(int it, double *h_densities, double *h_potentials, double *h_potentials_new);
int compute_matrix_elements(int it, double *h_densities, double *h_potentials, metadata_kzpca_fft *mdfft, double complex *h, double kz, double complex * me_d_dx, double complex * me_d_dy);
int compute_energy(int it, double *h_densities, double *h_potentials, double *energy, double *npart);
int compute_matrix_elements_of_momentum_operator(int nx, double complex *me);
int compute_angular_momentum_Lz(double *jx, double *jy, double *Lz);
// int test_me(int nx);
#endif

