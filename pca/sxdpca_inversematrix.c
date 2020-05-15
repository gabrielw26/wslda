/*
 * inverse.c
 *
 *  Created on: Oct 15, 2019
 *      Author: kaskadermike
 */

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

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    }

//	calculate minor of matrix OR build new matrix : k-had = minor
void calculate_minor(double **b, double **a, int i, int M){
	int h = 0;
	int k = 0;
	int j, l;
	for(l = 1; l < M; l++)
		for(j = 0; j < M; j++){
			if(j == i)
				continue;
				b[h][k] = a[l][j];
				k++;
			if(k == (M -1)){
				h++;
				k = 0;
			}
		}
} // end function


//	calculate determinte of matrix
double determinant(double **a, int M){
	double **b;
	int i;
	cppmallocl(b, M, double*);
	for (i = 0; i < M; i++){
		cppmallocl(b[i], M, double);
	}
	double sum = 0.0;
	if (M == 1) {
		return a[0][0];
	}
	else if(M == 2) {
		return (a[0][0] * a[1][1] - a[0][1] * a[1][0]);
	}
	else {
		for(i = 0; i < M; i++){
			calculate_minor(b, a, i, M);	// read function
			if (i % 2 == 0) sum += a[0][i] * determinant(b, (M - 1));	// read function	// sum = determinte matrix
			else sum += a[0][i] * -1. * determinant(b, (M - 1));	// read function	// sum = determinte matrix
		}
	}
	return sum;
} // end function


//	calculate transpose of matrix
void transpose(double **c, double **d, int M, double det){
	double b[M][M];
	int i, j;
	for (i = 0; i < M; i++){
		for (j = 0; j < M; j++){
			b[i][j] = c[j][i];
		}
	}
	det = 1. / det;
	for (i = 0; i < M; i++){
		for (j = 0; j < M; j++) {
			d[i][j] = b[i][j] * det;
		}
	} // array d[][] = inverse matrix
} // end function


//	calculate cofactor of matrix
int cofactor(double **a, double **d, int M, double det){
	double **b;
	int i, j, h, l;
	cppmallocl(b, M, double*);
	for (i = 0; i < M; i++){
		cppmallocl(b[i], M, double);
	}
	double **c;
	cppmallocl(c, M, double*);
	for (i = 0; i < M; i++){
		cppmallocl(c[i], M, double);
	}
	int m, k;
	for (h = 0; h < M; h++){
		for (l = 0; l < M; l++){
			m = 0;
			k = 0;
			for (i = 0; i < M; i++)
			for (j = 0; j < M; j++)
			if (i != h && j != l){
				b[m][k] = a[i][j];
				if (k < (M - 2))
					k++;
				else {
					k = 0;
					m++;
				}
			}
			if ((h + l) % 2 == 0) c[h][l] = determinant(b, (M - 1));	// c = cofactor Matrix
			else c[h][l] = -1. * determinant(b, (M - 1));
		}
	}
	transpose(c, d, M, det);	// read function
	return 1;
} // end function


// calculate inverse of matrix
void inverse(double **a, double **d, int M, double det){
	if(det == 0) printf("\nInverse of Entered Matrix is not possible\n");
	else if(M == 1){
		d[0][0] = 1;
	}
	else
	cofactor(a, d, M, det);	// read function
} // end function

