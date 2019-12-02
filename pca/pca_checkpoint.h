int load_nwf (MPI_Comm comm, char* inprefix,
                int* nwf, int* nwfip_out,
                int HowMany);

int load_all(double complex * h_wavefun, MPI_Comm comm, char* inprefix, 
		cufftDoubleComplex *d_wf,
    		cufftDoubleComplex *d_fkm1, cufftDoubleComplex *d_fkm2, cufftDoubleComplex *d_fkm3,
		double* d_potentials, double *t0,
                int* nwf, int* nwfip,
                double* h_fbetaEn, double* mu,double* ec, double* kF, double* eF, double* Effg,
		int HowMany);

int save_all(double complex * h_wavefun, MPI_Comm comm, char* outprefix,
                cufftDoubleComplex *d_wf,
                cufftDoubleComplex *d_fkm1, cufftDoubleComplex *d_fkm2, cufftDoubleComplex *d_fkm3,
                double* d_potentials, double *t0,
                int nwf, int nwfip,
                double* h_fbetaEn, double* mu,double *ec, double *kF, double* eF, double* Effg,
                int HowMany);

int load_all_45(double complex * h_wavefun, MPI_Comm comm, char* inprefix,
                cufftDoubleComplex *d_wf,
                cufftDoubleComplex *d_fkm1, cufftDoubleComplex *d_fkm2, 
		cufftDoubleComplex *d_fkm3, cufftDoubleComplex *d_fkm4,
                double* d_potentials, double *t0,
                int* nwf, int* nwfip,
                double* h_fbetaEn, double* mu,double* ec, double* kF, double* eF, double* Effg,
                int HowMany);

int save_all_45(double complex * h_wavefun, MPI_Comm comm, char* outprefix,
                cufftDoubleComplex *d_wf,
                cufftDoubleComplex *d_fkm1, cufftDoubleComplex *d_fkm2,
		cufftDoubleComplex *d_fkm3, cufftDoubleComplex *d_fkm4,
                double* d_potentials, double *t0,
                int nwf, int nwfip,
                double* h_fbetaEn, double* mu,double *ec, double *kF, double* eF, double* Effg,
                int HowMany);

