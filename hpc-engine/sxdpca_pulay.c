/**
 * W-SLDA Toolkit
 *
 * Add here description of content
 * */
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <lapacke.h>
#include "sxdpca_pulay.h"
#include "wslda_errors.h"

// function calculates dot product of two vectors, it takes size of vectors and their addresses in the list, it also takes in the list
double dot2(int n, int v1, int v2, double *list) {
    double sum=0;
    for(int i=0; i<n; i++){
        sum+=(list[v1+i]*list[v2+i]);
    }
    return sum;
}

// function adds given pair of vectors to the list
void addtolist(int listlenght, double *list, int n, double *vect,double *R_vect){
    for(int i=0; i<listlenght-2; i+=2){
        for(int j=0; j<n; j++){
            list[(i)*n+j]=list[(i+2)*n+j];
            list[(i+1)*n+j]=list[(i+3)*n+j];
        }
    }
    for(int j=0; j<n; j++){
        list[(listlenght-2)*n+j]=vect[j];
        list[(listlenght-1)*n+j]=R_vect[j];
    }
}


size_t pulay_compute_workspace(int order, int vector_length)
{
    int m = order;
    int n = vector_length;
    size_t workspace_size = 0;

    // Calculate the maximum size of Rmatrix and b arrays
    size_t Rmatrix_size = (m + 2) * (m + 2) * sizeof(double);
    size_t b_size = (m + 2) * sizeof(double);

    // Calculate space for the other arrays
    size_t listofPairs_size = (2 * m) * n * sizeof(double);
    size_t R_vect_size = n * sizeof(double);
    size_t vect_x_size = n * sizeof(double);

    // Calculate size for variables like R_count
    size_t variables_size = sizeof(int);

    // Add some extra space for other variables and overhead
    size_t extra_space = 1024;

    // Calculate the total workspace size
    workspace_size = Rmatrix_size + b_size + listofPairs_size + R_vect_size + vect_x_size + variables_size + extra_space;

    return workspace_size;
}

int pulay_algorithm(int order, int vector_length, double *in, double *out, double alpha, double *workspace)
{
    int m = order;
    int n = vector_length;

    // Calculate the offset for the variables in the workspace
    int Rmatrix_offset = 0;
    int b_offset = (m + 2) * (m + 2);
    int listofPairs_offset = b_offset + (m + 2);
    int R_vect_offset = listofPairs_offset + (2 * m) * n;
    int vect_x_offset = R_vect_offset + n;
    int R_count_offset = vect_x_offset + n;

    // Get the pointers to the variables in the workspace
    double *Rmatrix = &workspace[Rmatrix_offset];
    double *b = &workspace[b_offset];
    double *listofPairs = &workspace[listofPairs_offset];
    double *R_vect = &workspace[R_vect_offset];
    double *vect_x = &workspace[vect_x_offset];
    int *R_count = (int *)&workspace[R_count_offset];
    
    double *vect_x = in; // in density vector
    double alpha = alpha; // alpha is a constant get from function argument

    addtolist(2*m, listofPairs, n, vect_x, R_vect);
    
    if(*R_count<m){ // if R_count is less than m we construct smaller matrix size of R_count+1 by R_count+1
        for(int i=0; i<*R_count+1; i++){
            b[i]=0;
        }
        b[*R_count+1]=1;
        for(int j=0; j<(*R_count+2)*(*R_count+2); j++){ // loop going through the matrix and assigning values
            if(j>((*R_count+2)*(*R_count+1)-1)){
                if(j==((*R_count+2)*(*R_count+2)-1)){
                    Rmatrix[j]=0;
                } else {
                    Rmatrix[j]=1;
                }
            } else if((j+1)%(*R_count+2)==0){
                Rmatrix[j]=1;
            } else {
                if(j<=*R_count){
                    Rmatrix[j]=dot2(n, 2*m-1-2*(*R_count), 2*m-1-2*(*R_count-(j%(*R_count+1))), listofPairs);
                } else {
                    Rmatrix[j]=dot2(n, 2*m-1-2*(*R_count-(j-j%(*R_count+2))/(*R_count+2)), 2*m-1-2*(*R_count-(j%(*R_count+2))), listofPairs);
                }
            }
        }

        if( LAPACKE_dposv( LAPACK_ROW_MAJOR, 'U', *R_count+2, 1, Rmatrix, *R_count+2, b, *R_count+2) > 0 ) { // calculating Rmatrix*b=c where c is the vector of weights for mixing and if there is an error we set R_count=0
                *R_count=0;
        } else {
            for(int i=0; i<n; i++){
                vect_x[i]=0;
            }

            for(int i=0; i<n; i++){ // calculating new density vector
                for(int j=0; j<*R_count+1; j++){
                    vect_x[i]+=b[j]*(listofPairs[(2*m-2*(*R_count+1)+2*j)*n+i]+alpha*listofPairs[(2*m-2*(*R_count+1)+1+2*j)*n+i]);
                }
            }
        }
    } else {
        for(int i=0; i<m; i++){
            b[i]=0;
        }
        b[m]=1;
        for(int j=0; j<(m+1)*(m+1); j++){
            if(j>((m+1)*m-1)){
                if(j==((m+1)*(m+1)-1)){
                    Rmatrix[j]=0;
                } else {
                    Rmatrix[j]=1;
                }
            } else if((j+1)%(m+1)==0){
                Rmatrix[j]=1;
            } else {
                Rmatrix[j]=dot2(n, 2*m-1-2*(m-1-(j-j%(m+1))/(m+1)), 2*m-1-2*(m-1-(j%(m+1))), listofPairs);
            }
        }
        if( LAPACKE_dposv( LAPACK_ROW_MAJOR, 'U', m+1, 1, Rmatrix, m+1, b, *R_count+2) > 0 ) { // calculating Rmatrix*b=c where c is the vector of weights for mixing and if there is an error we set R_count=0
                *R_count=0;
        } else {

            for(int i=0; i<n; i++){
                vect_x[i]=0;
            }

            for(int i=0; i<n; i++){ // calculating new density vector
                for(int j=0; j<m; j++){
                    vect_x[i]+=b[j]*(listofPairs[(2*j)*n+i]+alpha*listofPairs[(2*j+1)*n+i]);
                }
            }
        }
    }
    *R_count++;

    return WSLDA_OK; // wslda_errors.h for errors
}

// TODO: remove later this info
/**
 * To include into compilation process:
 * Edit mk.st file
 * Add to DEPS file sxdpca_pulay.h
 *    DEPS = ....
 *           wslda_api_version.h sxdpca_pulay.h ../VERSION.h
 *
 * Add to SRCS file sxdpca_pulay.o (not make will change .o into .c
 *     SRCS = .... sxdpca_pulay.o
 *
 * Add pulay binary to executable solve-uniform:
 * (BINDIR)solve-uniform: $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS) $(DEPSLOCAL) $(OBJDIR)sxdpca_pulay.o
	$(CXX) -o $(BINDIR)solve-uniform $(CFLAGS) $(WSLDADIR)extensions/solve-uniform.c $(OBJDIR)sxdpca_pulay.o -lm -lfftw3
 * */

