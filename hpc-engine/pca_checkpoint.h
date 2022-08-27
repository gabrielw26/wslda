int load_nwf (MPI_Comm comm, char* inprefix,
                int* nwf, int* nwfip_out,
                int HowMany);

int load_all(double complex * h_wavefun, MPI_Comm comm, char* inprefix, 
		double complex *d_wf,
    		double complex *d_fkm1, double complex *d_fkm2, double complex *d_fkm3,
		double* d_potentials, double *t0,
                int* nwf, int* nwfip,
                double* h_fbetaEn, double* mu,double* ec, double* kF, double* eF, double* Effg, double *beta,
		int HowMany);

int save_all(double complex * h_wavefun, MPI_Comm comm, char* outprefix,
                double complex *d_wf,
                double complex *d_fkm1, double complex *d_fkm2, double complex *d_fkm3,
                double* d_potentials, double *t0,
                int nwf, int nwfip,
                double* h_fbetaEn, double* mu,double *ec, double *kF, double* eF, double* Effg, double *beta,
                int HowMany);

int load_all_45(double complex * h_wavefun, MPI_Comm comm, char* inprefix,
                double complex *d_wf,
                double complex *d_fkm1, double complex *d_fkm2, 
		double complex *d_fkm3, double complex *d_fkm4,
                double* d_potentials, double *t0,
                int* nwf, int* nwfip,
                double* h_fbetaEn, double* mu,double* ec, double* kF, double* eF, double* Effg, double *beta,
                int HowMany);

int save_all_45(double complex * h_wavefun, MPI_Comm comm, char* outprefix,
                double complex *d_wf,
                double complex *d_fkm1, double complex *d_fkm2,
		double complex *d_fkm3, double complex *d_fkm4,
                double* d_potentials, double *t0,
                int nwf, int nwfip,
                double* h_fbetaEn, double* mu,double *ec, double *kF, double* eF, double* Effg, double *beta,
                int HowMany);

