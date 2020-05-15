/*
 * broyden.h
 *
 *  Created on: Oct 9, 2019
 *      Author: kaskadermike
 */

#ifndef BROYDEN_H_
#define BROYDEN_H_

void update_mu(double **dens_in, double **dens_out, double *Vin, double *Vout, int M, int dim, double mu_a, double mu_b, double mu_a_old, double mu_b_old);

int Broyden(double *h_dens, double **dens_in, double **dens_out, int M, int dim, double omega_0, double omega_n, double omega_k, double alpha);

int Broyden_mu(double *h_dens, double **dens_in, double **dens_out, int M, int dim, double omega_0, double omega_n, double omega_k, double alpha,
		double mu_a, double mu_b);


#endif /* BROYDEN_H_ */
