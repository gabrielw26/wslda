#ifndef _SXDPCA_PULAY_H_
#define _SXDPCA_PULAY_H_

/**
 * This function calculates the size of workspace needed for Pulay algorithm
 * @param order order of Pulay algorithm
 * @param vector_length length of vectors
 * @return size of workspace in bytes
 * */
size_t pulay_compute_workspace(int order, int vector_length);

/**
 * This function is the main function of Pulay algorithm, it computes the new density vector given one
 * @param order order of Pulay algorithm
 * @param vector_length length of vectors
 * @param in input data
 * @param output output data
 * @param alpha mixing parameter
 * @param workspace working buffer, of size as given by pulay_compute_workspace(...)
 * @return error code
 * */
int pulay_algorithm(int order, int vector_length, double *in, double *out, double alpha, double *workspace);

double dot2(int n, int v1, int v2, double *list);
void addtolist(int listlenght, double *list, int n, double *vect,double *R_vect);

#endif
