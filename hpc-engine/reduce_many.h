#ifndef REDUCE_MANY_H
#define REDUCE_MANY_H

template <unsigned int blockSize>
void warpReduce(volatile double *shm, int tid);

template <unsigned int blockSize>
void reduce_many(double *g_idata, double *g_odata, int distBetweenSums, int mmm, int nblocksPerArray, int endVal);
void reduce_many_32(double *g_idata, double *g_odata, int distBetweenSums, int mmm, int nblocksPerArray, int endVal);
void call_template_reduce_many_kernels(int dimGrid, int dimBlock, double *d_idata, double *d_odata, 
    int distBetweenSums, int nElementsInVector, int nBlocksPerArray, int endVal);

void write_at_begining(double *g_odata, int nElementsInVector, int nVectors, int idxFirstArrayToRewrite);

int getNumberOfBlockDim(int nElementsToSum, int optNumOperation);
int local_reductions_many(int n, int m, double *d_in, double *d_out);

#endif