#ifndef __TDWSLDA_CHECKPOINT__
#define __TDWSLDA_CHECKPOINT__

int load_nwf (MPI_Comm comm, char* inprefix,
                int* nwf, int* nwfip_out,
                int HowMany);

int load_checkpoint(double complex * h_wavefun, MPI_Comm comm, char* inprefix,
                double complex *d_wf,
                double complex *d_fkm1,
                double complex *d_fkm2,
                double complex *d_fkm3,
                double complex *d_fkm4,
                double complex *d_fkm5,
                double* d_potentials, double* t0,
                int* nwf, int* nwfip_out,
                double* h_fbetaEn, double* h_kkz, double* mu,double* ec, double* kF, double* eF, double* Effg, double *beta,
                int HowMany, int abm_scheme);

int save_checkpoint(double complex * h_wavefun, MPI_Comm comm, char* outprefix,
                double complex *d_wf,
                double complex *d_fkm1,
                double complex *d_fkm2,
                double complex *d_fkm3,
                double complex *d_fkm4,
                double complex *d_fkm5,
                double* d_potentials, double *t0,
                int nwf, int nwfip,
                double* h_fbetaEn, double* h_kkz, double* mu, double *ec, double *kF, double* eF, double* Effg, double *beta,
                int HowMany, int abm_scheme);

#endif
