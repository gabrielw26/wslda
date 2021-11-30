#ifndef REDUCE_MANY_H
#define REDUCE_MANY_H

template <unsigned int blockSize>
void warpReduce(volatile double *shm, int tid, int nEndValue);

template <unsigned int blockSize>
void reduce_many(double *g_idata, double *g_odata, unsigned int nElementsInVector, 
    unsigned int nblocksPerArray, unsigned int nMemoryBlocksPerArray, unsigned int gridDimPerArray);

int getEndValue(int currArrayIdx, int nElementsInVector, int nMemoryBlocksPerArray, int currMemoryBlockIdx);
void call_template_kernels(int dimGrid, int dimBlock, double *d_idata, double *d_odata, 
    int nElementsInVector, int nblocksPerArray, int nMemoryBlocksPerArray);
int getNumberOfGridDim(int m, int optNumOperation, int blockDim);
int getNumberOfBlockDim(int m, int optNumOperation);
int local_reductions_many(int n, int m, double *d_in, double *d_out);

#endif