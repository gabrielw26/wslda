#ifndef __PCA_REDSUF__
#define __PCA_REDSUF__

double warpReduceSum(double val); 

double blockReduceSum(double val); 

void deviceReduceKernel(double* in, double* out, int N); 

void deviceReduce_suffle(double *in, double* out, int N, int threads); 

void deviceReduce_suffle_st(double *in, double* out, int N, int threads, cudaStream_t stream); 

void deviceReduceKernel_MANY(double* in, double* out, int N, int M); 

void deviceReduce_suffle_MANY_st(double *in, double* out, int N, int threads, cudaStream_t stream, int M, double* work); 

void deviceReduce_suffle_MANY(double *in, double* out, int N, int threads, int M, double* work);

#endif
