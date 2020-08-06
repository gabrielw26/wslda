/*
 * broyden.c
 *
 *  Created on: Oct 9, 2019
 *      Author: kaskadermike
 * 
 *  Modified by: Gabriel Wlazlowski, May 2020
 */


// #include "sxdpca_broyden.h"

#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#include "netlib-lapack.h"

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    }
    
// functions declaration
int ran(int min, int max);
void delta(double *dV, double *dF, double *Vn1, double *Vn, double *Fn1, double *Fn, int dim);
void calculate_F(double *F, double *Vin, double *Vout, int dim);
void update(double **dens_in, double **dens_out, double *Vin, double *Vout, int M, int dim);
void update_mu(double **dens_in, double **dens_out, double *Vin, double *Vout, int M, int dim, double mu_a, double mu_b, double mu_a_old, double mu_b_old);

// LAPACK ROUTINES FOR INVERSION
// call dgetrf( m, n, a, lda, ipiv, info )
extern void dgetrf_(int *m, int *n,  double *a, int *lda, int *ipiv, int *info );

// call dgetri( n, a, lda, ipiv, work, lwork, info )
extern void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info );

/**
 * Functions inverts square matrix a (NxN) and stores result into matrix b
 * */
int lapack_inversion(int M, double **a, double **b)
{
    int i,j; 
	double *A;
	cppmallocl(A, M*M, double);
	for(i=0; i<M; i++) for(j=0; j<M; j++) A[j*M +i]=a[i][j]; // copy to column-major format (Fortran style)
	
    // invert matrix
    // needed storage for ipiv matrix
    int *ipiv;
    cppmallocl(ipiv,M,int);
    
    int info;
    int lwork=-1;
    double testwork[1];
    
    // call for optimal workspace for _dgetri
    dgetri_(&M, A, &M, ipiv, testwork, &lwork, &info );
    if(info!=0) { printf("Error: memory query: _dgetri=%d\n", info); return 1;}
    lwork=(int)testwork[0];
//     printf("dgetri: request for memory: %d\n", lwork);
    double *work;
    cppmallocl(work,lwork,double);

    dgetrf_(&M, &M, A, &M, ipiv, &info );
    if(info!=0) { printf("Error: dgetrf_=%d\n", info); return info;}
    
    dgetri_(&M, A, &M, ipiv, work, &lwork, &info );
    if(info!=0) { printf("Error: dgetri_=%d\n", info); return info;}
    
	for(i=0; i<M; i++) for(j=0; j<M; j++) b[i][j] = A[j*M +i]; // copy to result table
	
	free(A);
    free(ipiv);
    free(work);
    
    return 0;
}
    
    
int ran(int min, int max) {
    int tmp;
    if (max >= min)
        max -= min;
    else {
        tmp = min - max;
        min = max;
        max = tmp;
    }
    return (rand() % (max - min + 1)) + min;
}

void calculate_F(double *F, double *Vin, double *Vout, int dim){
	int ixyz;
	for (ixyz = 0; ixyz < dim; ixyz++){
		F[ixyz] = Vout[ixyz] - Vin[ixyz];
	}
}

void delta(double *dV, double *dF, double *Vn1, double *Vn, double *Fn1, double *Fn, int dim){
	double denominator = 0.;
	double tmp;
	int ixyz;
	for (ixyz = 0; ixyz < dim; ixyz++){
		tmp = Fn1[ixyz] - Fn[ixyz];
		denominator += tmp * tmp;
	}
	denominator = 1. / sqrt(denominator);

	for(ixyz = 0; ixyz < dim; ixyz++){
		dV[ixyz] = (Vn1[ixyz] - Vn[ixyz]) * denominator;
		dF[ixyz] = (Fn1[ixyz] - Fn[ixyz]) * denominator;
	}
}

void update(double **dens_in, double **dens_out, double *Vin, double *Vout, int M, int dim){
	int i, j;
	for (i = 0; i < M; i++){
		for (j = 0; j < dim; j++){
			dens_in[i][j] = dens_in[i + 1][j];
			dens_out[i][j] = dens_out[i + 1][j];
		}
	}
	for (i = 0; i < dim; i++){
		dens_in[M][i] = Vin[i];
		dens_out[M][i] = Vout[i];
	}
}

void update_mu(double **dens_in, double **dens_out, double *Vin, double *Vout, int M, int dim, double mu_a, double mu_b, double mu_a_old, double mu_b_old){
	int i, j;
	for (i = 0; i < M; i++){
		for (j = 0; j < dim+2; j++){
			dens_in[i][j] = dens_in[i + 1][j];
			dens_out[i][j] = dens_out[i + 1][j];
		}
	}
	for (i = 0; i < dim; i++){
		dens_in[M][i] = Vin[i];
		dens_out[M][i] = Vout[i];
	}
	dens_in[M][i+0] = mu_a_old;
	dens_out[M][i+0] = mu_a;
	dens_in[M][i+1] = mu_b_old;
	dens_out[M][i+1] = mu_b;
}

int Broyden(double *h_dens, double **dens_in, double **dens_out, int M, int dim, double omega_0, double omega_n, double omega_k, double alpha){

	int ixyz, i, j, k, n;
	double ckm;
	double gamma_mn;

	double *u;
	double *delta_Vn;
	double *delta_Vk;
	double *delta_Fn;
	double *delta_Fk;
	double *Fn;
	double *Fn1;
	double *Fk;
	double *Fk1;
	double *sum;
	double **a;
	double **tmp;
	double **beta;

	cppmallocl(u, dim, double);
	cppmallocl(delta_Vn, dim, double);
	cppmallocl(delta_Vk, dim, double);
	cppmallocl(delta_Fn, dim, double);
	cppmallocl(delta_Fk, dim, double);
	cppmallocl(Fn, dim, double);
	cppmallocl(Fn1, dim, double);
	cppmallocl(Fk, dim, double);
	cppmallocl(Fk1, dim, double);
	cppmallocl(sum, dim, double);
	cppmallocl(beta, M, double*);
	cppmallocl(a, M, double*);
	cppmallocl(tmp, dim, double*);
	for (i = 0; i < M; i++){
		cppmallocl(a[i], M, double);
		cppmallocl(beta[i], M, double);
	}
	for(i = 0; i < dim; i++){
		cppmallocl(tmp[i], M, double);
	}
	for (i = 0; i < M; i++){
		for (j = 0; j < M; j++){
			a[i][j] = 0.;
			beta[i][j] = 0.;
		}
	}

	for (n = 0; n < M; n++){
		calculate_F(Fn, dens_in[n], dens_out[n], dim); 							// F(n) for delta F(n) fraction
		calculate_F(Fn1, dens_in[n+1], dens_out[n+1], dim);						// F(n+1) for delta F(n) fraction
		delta(delta_Vn, delta_Fn, dens_in[n+1], dens_in[n], Fn1, Fn, dim);	 	// both delta V(n) and delta F(n) fractions
		for (k = 0; k < M; k++){
			calculate_F(Fk, dens_in[k], dens_out[k], dim);						// F(k) for delta F(n) fraction
			calculate_F(Fk1, dens_in[k+1], dens_out[k+1], dim); 				// F(k+1) for delta F(n) fraction
			delta(delta_Vk, delta_Fk, dens_in[k+1], dens_in[k], Fk1, Fk, dim); 	// both delta V(k)(not used) and delta F(k) fractions
			for (ixyz = 0; ixyz < dim; ixyz++){ 								// a matrix
				a[k][n] += omega_n * omega_k * delta_Fn[ixyz] * delta_Fk[ixyz];
				sum[ixyz] = 0.;
			}
		}
	}

	// overwriting "a" matrix as sum of "a" and omega_0^2 times I
	for (n = 0; n < M; n++){
		for (k = 0; k < M; k++){
			if (k == n){
				a[k][n] = omega_0 * omega_0 + a[k][n];
			}
			else {
				continue;
			}
		}
	}

	// inversion to beta
	lapack_inversion(M, a, beta);
    

	for (n = 0; n < M; n++){
		ckm = 0.;
		gamma_mn = 0.;
		calculate_F(Fn, dens_in[n], dens_out[n], dim); 							// F(n) for delta F(n) fraction
		calculate_F(Fn1, dens_in[n+1], dens_out[n+1], dim);						// F(n+1) for delta F(n) fraction
		delta(delta_Vn, delta_Fn, dens_in[n+1], dens_in[n], Fn1, Fn, dim);	 	// both delta V(n) and delta F(n) fractions
		for (ixyz = 0; ixyz < dim; ixyz++){										// u vector
			u[ixyz] = alpha * delta_Fn[ixyz] + delta_Vn[ixyz];
		}
		for (k = 0; k < M; k++){
			ckm = 0.;
			calculate_F(Fk, dens_in[k], dens_out[k], dim);						// F(k) for delta F(n) fraction
			calculate_F(Fk1, dens_in[k+1], dens_out[k+1], dim); 				// F(k+1) for delta F(n) fraction
			delta(delta_Vk, delta_Fk, dens_in[k+1], dens_in[k], Fk1, Fk, dim); 	// both delta V(k)(not used) and delta F(k) fractions
			for (ixyz = 0; ixyz < dim; ixyz++){ 								// ckm scalar and a matrix
				ckm += delta_Fk[ixyz] * (dens_out[M][ixyz] - dens_in[M][ixyz]);
			}
			ckm = ckm * omega_k;
			gamma_mn += ckm * beta[k][n];
		}

		for (ixyz = 0; ixyz < dim; ixyz++){
			tmp[ixyz][n] = omega_n * gamma_mn * u[ixyz];
		}
	}

	// sum up
	for (ixyz = 0; ixyz < dim; ixyz++){
		for (n = 0; n < M; n++){
			sum[ixyz] += tmp[ixyz][n];
		}
	}

	// calculate final density vector
	for (ixyz = 0; ixyz < dim; ixyz++){
		h_dens[ixyz] = alpha * dens_out[M][ixyz] + (1. -  alpha) * dens_in[M][ixyz] - sum[ixyz];
	}
	
    // clear memory
	for (i = 0; i < M; i++){
		free(a[i]);
		free(beta[i]);
	}
	for(i = 0; i < dim; i++){
		free(tmp[i]);
	}
	free(u);
	free(delta_Vn);
	free(delta_Vk);
	free(delta_Fn);
	free(delta_Fk);
	free(Fn);
	free(Fn1);
	free(Fk);
	free(Fk1);
	free(sum);
	free(beta);
	free(a);
	free(tmp);

	return 1;
}

int Broyden_mu(double *h_dens, double **dens_in, double **dens_out, int M, int dim, double omega_0, double omega_n, double omega_k, double alpha,
		double *mu_a, double *mu_b){

	int ixyz, i, j, k, n;
	double ckm;
	double gamma_mn;

	double *u;
	double *delta_Vn;
	double *delta_Vk;
	double *delta_Fn;
	double *delta_Fk;
	double *Fn;
	double *Fn1;
	double *Fk;
	double *Fk1;
	double *sum;
	double **a;
	double **tmp;
	double **beta;

	cppmallocl(u, dim, double);
	cppmallocl(delta_Vn, dim, double);
	cppmallocl(delta_Vk, dim, double);
	cppmallocl(delta_Fn, dim, double);
	cppmallocl(delta_Fk, dim, double);
	cppmallocl(Fn, dim, double);
	cppmallocl(Fn1, dim, double);
	cppmallocl(Fk, dim, double);
	cppmallocl(Fk1, dim, double);
	cppmallocl(sum, dim, double);
	cppmallocl(beta, M, double*);
	cppmallocl(a, M, double*);
	cppmallocl(tmp, dim, double*);
	for (i = 0; i < M; i++){
		cppmallocl(a[i], M, double);
		cppmallocl(beta[i], M, double);
	}
	for(i = 0; i < dim; i++){
		cppmallocl(tmp[i], M, double);
	}
	for (i = 0; i < M; i++){
		for (j = 0; j < M; j++){
			a[i][j] = 0.;
			beta[i][j] = 0.;
		}
	}
	
	for (n = 0; n < M; n++){
		calculate_F(Fn, dens_in[n], dens_out[n], dim); 							// F(n) for delta F(n) fraction
		calculate_F(Fn1, dens_in[n+1], dens_out[n+1], dim);						// F(n+1) for delta F(n) fraction
		delta(delta_Vn, delta_Fn, dens_in[n+1], dens_in[n], Fn1, Fn, dim);	 	// both delta V(n) and delta F(n) fractions
		for (k = 0; k < M; k++){
			calculate_F(Fk, dens_in[k], dens_out[k], dim);						// F(k) for delta F(n) fraction
			calculate_F(Fk1, dens_in[k+1], dens_out[k+1], dim); 				// F(k+1) for delta F(n) fraction
			delta(delta_Vk, delta_Fk, dens_in[k+1], dens_in[k], Fk1, Fk, dim); 	// both delta V(k)(not used) and delta F(k) fractions
			for (ixyz = 0; ixyz < dim; ixyz++){ 								// a matrix
				a[k][n] += omega_n * omega_k * delta_Fn[ixyz] * delta_Fk[ixyz];
				sum[ixyz] = 0.;
			}
		}
	}

	// overwriting "a" matrix as sum of "a" and omega_0^2 times I
	for (n = 0; n < M; n++){
		for (k = 0; k < M; k++){
			if (k == n){
				a[k][n] = omega_0 * omega_0 + a[k][n];
			}
			else {
				continue;
			}
		}
	}
	
	
	// inversion to beta
	lapack_inversion(M, a, beta);

	for (n = 0; n < M; n++){
		ckm = 0.;
		gamma_mn = 0.;
		calculate_F(Fn, dens_in[n], dens_out[n], dim); 							// F(n) for delta F(n) fraction
		calculate_F(Fn1, dens_in[n+1], dens_out[n+1], dim);						// F(n+1) for delta F(n) fraction
		delta(delta_Vn, delta_Fn, dens_in[n+1], dens_in[n], Fn1, Fn, dim);	 	// both delta V(n) and delta F(n) fractions
		for (ixyz = 0; ixyz < dim; ixyz++){										// u vector
			u[ixyz] = alpha * delta_Fn[ixyz] + delta_Vn[ixyz];
		}
		for (k = 0; k < M; k++){
			ckm = 0.;
			calculate_F(Fk, dens_in[k], dens_out[k], dim);						// F(k) for delta F(n) fraction
			calculate_F(Fk1, dens_in[k+1], dens_out[k+1], dim); 				// F(k+1) for delta F(n) fraction
			delta(delta_Vk, delta_Fk, dens_in[k+1], dens_in[k], Fk1, Fk, dim); 	// both delta V(k)(not used) and delta F(k) fractions
			for (ixyz = 0; ixyz < dim; ixyz++){ 								// ckm scalar and a matrix
				ckm += delta_Fk[ixyz] * (dens_out[M][ixyz] - dens_in[M][ixyz]);
			}
			ckm = ckm * omega_k;
			gamma_mn += ckm * beta[k][n];
		}

		for (ixyz = 0; ixyz < dim; ixyz++){
			tmp[ixyz][n] = omega_n * gamma_mn * u[ixyz];
		}
	}

	// sum up
	for (ixyz = 0; ixyz < dim; ixyz++){
		for (n = 0; n < M; n++){
			sum[ixyz] += tmp[ixyz][n];
		}
	}

	// calculate final density vector
	for (ixyz = 0; ixyz < (dim-2); ixyz++){
		h_dens[ixyz] = alpha * dens_out[M][ixyz] + (1. -  alpha) * dens_in[M][ixyz] - sum[ixyz];
	}
	*mu_a = alpha * dens_out[M][dim-2] + (1. -  alpha) * dens_in[M][dim-2] - sum[dim-2];
	*mu_b = alpha * dens_out[M][dim-1] + (1. -  alpha) * dens_in[M][dim-1] - sum[dim-1];

    // clear memory
	for (i = 0; i < M; i++){
		free(a[i]);
		free(beta[i]);
	}
	for(i = 0; i < dim; i++){
		free(tmp[i]);
	}
	free(u);
	free(delta_Vn);
	free(delta_Vk);
	free(delta_Fn);
	free(delta_Fk);
	free(Fn);
	free(Fn1);
	free(Fk);
	free(Fk1);
	free(sum);
	free(beta);
	free(a);
	free(tmp);
    
	return 1;
}
