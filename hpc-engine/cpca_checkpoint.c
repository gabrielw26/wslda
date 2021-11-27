#include <cuda.h>
#include <cufft.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <complex.h>

#include "pca_settings.h"
#include "pca_macro.h"
#include "pca_kernels.h"

#ifndef CODEDIM
    CODEDIM NOT DEFINED!
#endif

#if CODEDIM==1
#define CHUNKN NX
#elif CODEDIM==2
#define CHUNKN NXY
#else
    INCORRECT CODEDIM
#endif

void getnwfip( int ip , int np , int nwf , int * nwfip );


int load_nwf (MPI_Comm comm, char* inprefix,
                int* nwf, int* nwfip_out,
                int HowMany)
{
 char file_name[256];

 int comm_size, comm_rank;
 int nr, ierr;
 *nwf=-1;

 size_t shift_0;

 MPI_Comm comm_io;

 MPI_Status status;
 MPI_File in;

 MPI_Comm_size(comm, &comm_size);
 MPI_Comm_rank(comm, &comm_rank);
 nr = comm_rank / HowMany;

 sprintf(file_name, "%s/checkpoint_%d", inprefix, nr);

 MPI_Comm_split(comm, nr, comm_rank, &comm_io);

 MPI_File_open (comm_io, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &in);

 if ( in == NULL ) {
    return WSLDA_ERR_TD_CANNOT_LOAD_CHECKPOINT_NWF;
 }

 shift_0 = sizeof(double);
 shift_0 += 2 * sizeof(double);
 shift_0 += sizeof(double);
 shift_0 += sizeof(double);
 MPI_File_read_at(in, shift_0, nwf, 1, MPI_INT, &status);
//compute nwfip
 getnwfip( comm_rank , comm_size , *nwf , nwfip_out );

 MPI_File_close (&in);
 
 if(*nwf<=0) return WSLDA_ERR_TD_CANNOT_LOAD_CHECKPOINT_NWF;

 return 0;
}         
                        
int load_all (double complex * h_wavefun, MPI_Comm comm, char* inprefix, 
		cufftDoubleComplex *d_wf,
    		cufftDoubleComplex *d_fkm1, cufftDoubleComplex *d_fkm2, cufftDoubleComplex *d_fkm3,
		double* d_potentials, double* t0,
		int* nwf, int* nwfip_out,
		double* h_fbetaEn, double* h_kkz, double* mu,double* ec, double* kF, double* eF, double* Effg, double *beta,
		int HowMany)
{
 char file_name[256];
   
 int comm_size, comm_rank;
 int nr, ierr, ip;
 int nwfip, nwfp;

 size_t shift, shift_0;

 MPI_Comm comm_io;

 MPI_Status status;
 MPI_File in;

 MPI_Comm_size(comm, &comm_size);
 MPI_Comm_rank(comm, &comm_rank); 

 nr = comm_rank / HowMany;

 sprintf(file_name, "%s/checkpoint_%d", inprefix, nr);

 MPI_Comm_split(comm, nr, comm_rank, &comm_io);
 MPI_Comm_rank(comm_io, &ip);
 MPI_Comm_size(comm_io, &nr);

 MPI_File_open (comm_io, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &in);

 if ( in == NULL ) {
    return WSLDA_ERR_TD_CANNOT_LOAD_CHECKPOINT_DATA;
 }

 MPI_File_read_at(in, 0, t0, 1, MPI_DOUBLE, &status);
 shift_0 = sizeof(double);
 MPI_File_read_at(in, shift_0, mu, 2, MPI_DOUBLE, &status);
 shift_0 += 2 * sizeof(double);
 MPI_File_read_at(in, shift_0, ec, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
 MPI_File_read_at(in, shift_0, kF, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
 MPI_File_read_at(in, shift_0, nwf, 1, MPI_INT, &status);
 shift_0 += sizeof(int); 
 MPI_File_read_at(in, shift_0, eF, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
 MPI_File_read_at(in, shift_0, Effg, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
 MPI_File_read_at(in, shift_0, beta, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);

//compute nwfip
 getnwfip (comm_rank, comm_size, *nwf, nwfip_out);
 nwfip = *nwfip_out;

 if (comm_rank < *nwf%comm_size) nwfp = 0;
 else nwfp = 1;

 shift = shift_0 + sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_read_at(in, shift, h_fbetaEn, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;
 
 shift = shift_0;
 shift += sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_read_at(in, shift, h_kkz, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;
 
#if CODEDIM==1
 // h_kkz keeps both ky and kz
 shift = shift_0;
 shift += sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_read_at(in, shift, h_kkz+nwfip, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;
#endif
 
//d_wf 
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;//(comm_rank % HowMany);
 MPI_File_read_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 gpu_exec( memcopy_host2gpu(h_wavefun, d_wf,  (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr; 
//d_fkm1
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 MPI_File_read_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm1,  (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//_fkm2
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 MPI_File_read_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm2,  (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fkm3
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 MPI_File_read_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm3,  (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//V_a, V_b, DELTA
 shift = shift_0;
 // everyone reads
 {
     size_t _i;
     for(_i=0; _i<12; _i++) // because h_wavefun may not have enough space to store all elements
     {
        MPI_File_read_at(in, shift+_i*CHUNKN*sizeof(double), h_wavefun, (size_t)1*CHUNKN, MPI_DOUBLE, &status);
        gpu_exec( memcopy_host2gpu(h_wavefun, d_potentials+_i*CHUNKN, (size_t)1*CHUNKN*sizeof(double)) );
        
     }
 }

 MPI_File_close (&in);
 return 0;
}

int save_all(double complex * h_wavefun, MPI_Comm comm, char* outprefix,
                cufftDoubleComplex *d_wf,
                cufftDoubleComplex *d_fkm1, cufftDoubleComplex *d_fkm2, cufftDoubleComplex *d_fkm3,
                double* d_potentials, double *t0,
		int nwf, int nwfip,
                double* h_fbetaEn, double* h_kkz, double* mu, double *ec, double *kF, double* eF, double* Effg, double *beta,
                int HowMany)
{
 char file_name[256];

 int comm_size, comm_rank;
 int nr, ierr, nwfp, ip;

 size_t shift, shift_0;

 MPI_Comm comm_io;

 MPI_Status status; 
 MPI_File in; 

 MPI_Comm_size(comm, &comm_size);
 MPI_Comm_rank(comm, &comm_rank);
 ip = comm_rank;

 nr = comm_rank / HowMany;

 sprintf(file_name, "%s/checkpoint_%d", outprefix, nr);
    
 MPI_Comm_split(comm, nr, comm_rank, &comm_io);
 MPI_Comm_rank(comm_io, &ip);
 MPI_Comm_size(comm_io, &nr);

 MPI_File_open (comm_io, file_name, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &in);

 if ( in == NULL ) {
    perror ( "Unable to open the file" );
    MPI_Abort( MPI_COMM_WORLD , -1 );
    exit ( EXIT_FAILURE );
 }

if (ip==0) MPI_File_write_at(in, 0, t0, 1, MPI_DOUBLE, &status);
 shift_0 = sizeof(double);
if (ip==0) MPI_File_write_at(in, shift_0, mu, 2, MPI_DOUBLE, &status);
 shift_0 += 2 * sizeof(double);
if (ip==0) MPI_File_write_at(in, shift_0, ec, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
if (ip==0) MPI_File_write_at(in, shift_0, kF, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
if (ip==0) MPI_File_write_at(in, shift_0, &nwf, 1, MPI_INT, &status);
 shift_0 += sizeof(int);
if(ip==0) MPI_File_write_at(in, shift_0, eF, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
if(ip==0) MPI_File_write_at(in, shift_0, Effg, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
if(ip==0) MPI_File_write_at(in, shift_0, beta, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);


 if (comm_rank < nwf % comm_size) nwfp = 0;
 else nwfp = 1;

 shift = shift_0 + sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_write_at(in, shift, h_fbetaEn, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;

 shift = shift_0;
 shift += sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_write_at(in, shift, h_kkz, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;
 
#if CODEDIM==1
 // h_kkz keeps both ky and kz
 shift = shift_0;
 shift += sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_write_at(in, shift, h_kkz+nwfip, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;
#endif

//d_wf 
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 gpu_exec( memcopy_gpu2host(d_wf, h_wavefun, (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 MPI_File_write_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fkm1
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 gpu_exec( memcopy_gpu2host(d_fkm1, h_wavefun, (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 MPI_File_write_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fk2
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 gpu_exec( memcopy_gpu2host(d_fkm2, h_wavefun, (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 MPI_File_write_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fkm3
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 gpu_exec( memcopy_gpu2host(d_fkm3, h_wavefun, (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 MPI_File_write_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//V_a, V_b, DELTA
 shift = shift_0;
 if(ip==0) // writes only the ip==0
 {
     size_t _i;
     for(_i=0; _i<12; _i++) // because h_wavefun may not have enough space to store all elements
     {
        gpu_exec( memcopy_gpu2host(d_potentials+_i*CHUNKN, h_wavefun, (size_t)1*CHUNKN*sizeof(double)) );
        MPI_File_write_at(in, shift+_i*CHUNKN*sizeof(double), h_wavefun, (size_t)1*CHUNKN, MPI_DOUBLE, &status);
     }
 }
 
 MPI_File_close (&in);
 return 0;
}


int load_all_45 (double complex * h_wavefun, MPI_Comm comm, char* inprefix,
                cufftDoubleComplex *d_wf,
                cufftDoubleComplex *d_fkm1, cufftDoubleComplex *d_fkm2, 
		cufftDoubleComplex *d_fkm3, cufftDoubleComplex *d_fkm4,		
                double* d_potentials, double* t0,
                int* nwf, int* nwfip_out,
                double* h_fbetaEn, double* h_kkz, double* mu,double* ec, double* kF, double* eF, double* Effg, double *beta,
                int HowMany)
{
char file_name[256];

 int comm_size, comm_rank;
 int nr, ierr, ip;
 int nwfip, nwfp;

 size_t shift, shift_0;

 MPI_Comm comm_io;

 MPI_Status status;
 MPI_File in;

 MPI_Comm_size(comm, &comm_size);
 MPI_Comm_rank(comm, &comm_rank);

 nr = comm_rank / HowMany;

 sprintf(file_name, "%s/checkpoint_%d", inprefix, nr);

 MPI_Comm_split(comm, nr, comm_rank, &comm_io);
 MPI_Comm_rank(comm_io, &ip);
 MPI_Comm_size(comm_io, &nr);

 MPI_File_open (comm_io, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &in);

 if ( in == NULL ) {
    return WSLDA_ERR_TD_CANNOT_LOAD_CHECKPOINT_DATA;
 }

 MPI_File_read_at(in, 0, t0, 1, MPI_DOUBLE, &status);
 shift_0 = sizeof(double);
 MPI_File_read_at(in, shift_0, mu, 2, MPI_DOUBLE, &status);
 shift_0 += 2 * sizeof(double);
 MPI_File_read_at(in, shift_0, ec, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
 MPI_File_read_at(in, shift_0, kF, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
 MPI_File_read_at(in, shift_0, nwf, 1, MPI_INT, &status);
 shift_0 += sizeof(int);
 MPI_File_read_at(in, shift_0, eF, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
 MPI_File_read_at(in, shift_0, Effg, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
 MPI_File_read_at(in, shift_0, beta, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);

//compute nwfip
 getnwfip (comm_rank, comm_size, *nwf, nwfip_out);
 nwfip = *nwfip_out;

 if (comm_rank < *nwf%comm_size) nwfp = 0;
 else nwfp = 1;

 shift = shift_0 + sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_read_at(in, shift, h_fbetaEn, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;
 
 shift = shift_0;
 shift += sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_read_at(in, shift, h_kkz, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;
 
#if CODEDIM==1
 // h_kkz keeps both ky and kz
 shift = shift_0;
 shift += sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_read_at(in, shift, h_kkz+nwfip, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;
#endif 

//d_wf 
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;//(comm_rank % HowMany);
 MPI_File_read_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 gpu_exec( memcopy_host2gpu(h_wavefun, d_wf,  (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fkm1
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 MPI_File_read_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm1,  (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//_fkm2
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 MPI_File_read_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm2,  (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fkm3
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 MPI_File_read_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm3,  (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fkm4
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 MPI_File_read_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm4,  (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//V_a, V_b, DELTA
 shift = shift_0;
 // everyone reads
 {
     size_t _i;
     for(_i=0; _i<12; _i++) // because h_wavefun may not have enough space to store all elements
     {
        MPI_File_read_at(in, shift+_i*CHUNKN*sizeof(double), h_wavefun, (size_t)1*CHUNKN, MPI_DOUBLE, &status);
        gpu_exec( memcopy_host2gpu(h_wavefun, d_potentials+_i*CHUNKN, (size_t)1*CHUNKN*sizeof(double)) );
        
     }
 }

 MPI_File_close (&in);
 return 0;
}


int save_all_45(double complex * h_wavefun, MPI_Comm comm, char* outprefix,
                cufftDoubleComplex *d_wf,
                cufftDoubleComplex *d_fkm1, cufftDoubleComplex *d_fkm2, 
		cufftDoubleComplex *d_fkm3, cufftDoubleComplex *d_fkm4,
                double* d_potentials, double *t0,
                int nwf, int nwfip,
                double* h_fbetaEn, double* h_kkz, double* mu, double *ec, double *kF, double* eF, double* Effg, double *beta,
                int HowMany)
{
 char file_name[256];

 int comm_size, comm_rank;
 int nr, ierr, nwfp, ip;

 size_t shift, shift_0;

 MPI_Comm comm_io;

 MPI_Status status;
 MPI_File in;

 MPI_Comm_size(comm, &comm_size);
 MPI_Comm_rank(comm, &comm_rank);
 ip = comm_rank;

 nr = comm_rank / HowMany;

 sprintf(file_name, "%s/checkpoint_%d", outprefix, nr);

 MPI_Comm_split(comm, nr, comm_rank, &comm_io);
 MPI_Comm_rank(comm_io, &ip);
 MPI_Comm_size(comm_io, &nr);

 MPI_File_open (comm_io, file_name, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &in);

 if ( in == NULL ) {
    perror ( "Unable to open the file" );
    MPI_Abort( MPI_COMM_WORLD , -1 );
    exit ( EXIT_FAILURE );
 }

if (ip==0) MPI_File_write_at(in, 0, t0, 1, MPI_DOUBLE, &status);
 shift_0 = sizeof(double);
if (ip==0) MPI_File_write_at(in, shift_0, mu, 2, MPI_DOUBLE, &status);
 shift_0 += 2 * sizeof(double);
if (ip==0) MPI_File_write_at(in, shift_0, ec, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
if (ip==0) MPI_File_write_at(in, shift_0, kF, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
if (ip==0) MPI_File_write_at(in, shift_0, &nwf, 1, MPI_INT, &status);
 shift_0 += sizeof(int);
if(ip==0) MPI_File_write_at(in, shift_0, eF, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
if(ip==0) MPI_File_write_at(in, shift_0, Effg, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);
if(ip==0) MPI_File_write_at(in, shift_0, beta, 1, MPI_DOUBLE, &status);
 shift_0 += sizeof(double);


 if (comm_rank < nwf % comm_size) nwfp = 0;
 else nwfp = 1;

 shift = shift_0 + sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_write_at(in, shift, h_fbetaEn, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;

 shift = shift_0;
 shift += sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_write_at(in, shift, h_kkz, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;
 
#if CODEDIM==1
 // h_kkz keeps both ky and kz
 shift = shift_0;
 shift += sizeof(double) * (nwfip +nwfp) * ip;
 MPI_File_write_at(in, shift, h_kkz+nwfip, nwfip, MPI_DOUBLE, &status);
 shift_0 += sizeof(double) * (nwfip+nwfp) * nr ;
#endif 

//d_wf
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 gpu_exec( memcopy_gpu2host(d_wf, h_wavefun, (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 MPI_File_write_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fkm1
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 gpu_exec( memcopy_gpu2host(d_fkm1, h_wavefun, (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 MPI_File_write_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fk2
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 gpu_exec( memcopy_gpu2host(d_fkm2, h_wavefun, (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 MPI_File_write_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fkm3
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 gpu_exec( memcopy_gpu2host(d_fkm3, h_wavefun, (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 MPI_File_write_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;
//d_fkm4
 shift = shift_0;
 shift += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * ip;
 gpu_exec( memcopy_gpu2host(d_fkm4, h_wavefun, (size_t)2*nwfip*CHUNKN*sizeof(cufftDoubleComplex)) );
 MPI_File_write_at(in, shift, h_wavefun, (size_t)2*nwfip*CHUNKN*2, MPI_DOUBLE, &status);
 shift_0 += (size_t)2*(nwfip+nwfp)*CHUNKN*sizeof(cufftDoubleComplex) * nr;   
//V_a, V_b, DELTA
 shift = shift_0;
 if(ip==0) // writes only the ip==0
 {
     size_t _i;
     for(_i=0; _i<12; _i++) // because h_wavefun may not have enough space to store all elements
     {
        gpu_exec( memcopy_gpu2host(d_potentials+_i*CHUNKN, h_wavefun, (size_t)1*CHUNKN*sizeof(double)) );
        MPI_File_write_at(in, shift+_i*CHUNKN*sizeof(double), h_wavefun, (size_t)1*CHUNKN, MPI_DOUBLE, &status);
     }
 }
 
 MPI_File_close (&in);
 return 0;
}
