// This is code for solving stationary DFT equations for polarized cold atoms 
// in 3D, without any symmetry restrictions

// Authors:
// Gabriel Wlazlowski <gabriel.wlazlowski@pw.edu.pl>

// NOTE: Inspect for NOTE flags the code before running it!!!

// #define VERBOSE
// #define S3DDEBUG


#if DIAGONALIZATION_ROUTINE==PZHEEVR
#define USE_SCALAPACK_PZHEEVR
#elif DIAGONALIZATION_ROUTINE==PZHEEVD
#define USE_SCALAPACK_PZHEEVD
#else
    select DIAGONALIZATION ROUTINE in pca_settings
    // #define USE_SCALAPACK_PZHEEV
#endif

// #define USE_ELPA
// #define USE_ELPA_NEV_FRACTION 0.775

#ifdef S3DDEBUG
#define ECHOLINE                                                                                                        \
    {                                                                                                                   \
        if(iam==0) {printf("# ECHOLINE REACHED LINE %d IN FILE %s\n", __LINE__ , __FILE__); fflush(stdout); }           \
        MPI_Barrier(MPI_COMM_WORLD);                                                                                    \
    }
#else
#define ECHOLINE
#endif

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>

#include "pca_settings.h"
#include "pca_macro.h"
#include "pca_utils.h"
#include "pca_io.h"
#include "s2dpca_edf.h"
#include "pca_uniform.h"
#include "pca_logger.h"
#include "s3dpca_fft.h"
#include "s3dpca_me.h"
#include "s3dpca_grid.h"
#include "s3dpca_densities.h"
#include "s3dpca_inittype4.h"
// #include "s3dpca_testfun.h"

#ifdef USE_SCALAPACK_PZHEEVR
/* PZHEEVR prototype */
extern void pzheevr_(char* jobz, char* range, char* uplo, int* n, double complex* a, int* ia, int* ja, int* desca, 
                    double* vl, double* vu, int* il, int* iu, int* m, int* nz, double* w, double complex* z, int* iz, int* jz, int* descz, 
                    double complex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info);
extern void pdsyevr_(char* jobz, char* range, char* uplo, int* n, double*         a, int* ia, int* ja, int* desca, 
                    double* vl, double* vu, int* il, int* iu, int* m, int* nz, double* w, double*         z, int* iz, int* jz, int* descz, 
                    double*         work, int* lwork,                             int* iwork, int* liwork, int* info);
#endif
#ifdef USE_SCALAPACK_PZHEEVD
/* PZHEEVD prototype */
extern void pzheevd_ (char *jobz , char *uplo , int *n , double complex *a , int *ia , int *ja , int *desca , 
             double *w , double complex *z , int *iz , int *jz , int *descz , 
             double complex *work , int *lwork , double *rwork , int *lrwork , 
             int *iwork, int *liwork, int *info );
extern void pdsyevd_ (char *jobz , char *uplo , int *n , double         *a , int *ia , int *ja , int *desca , 
             double *w , double         *z , int *iz , int *jz , int *descz , 
             double         *work , int *lwork , 
             int *iwork , int *liwork , int *info );
#endif
#ifdef USE_SCALAPACK_PZHEEV
/* PZHEEVD prototype */
extern void pzheev_ (char *jobz , char *uplo , int *n , double complex *a , int *ia , int *ja , int *desca , 
             double *w , double complex *z , int *iz , int *jz , int *descz , 
             double complex *work , int *lwork , double *rwork , int *lrwork , int *info );
#endif
#ifdef USE_ELPA
#include <elpa/elpa.h>

#define assert_elpa_ok(x) assert(x == ELPA_OK)
#endif

/* Matrix Redistribution function prototype */
void Cpzgemr2d(int m, int n, double complex *A, int IA, int JA, int *descA, 
               double complex *B, int IB, int JB, int *descB, int gcontext);

/* other routines used in the main driver */
void Cblacs_pinfo( int * , int * ) ;
int  Csys2blacs_handle( MPI_Comm comm ) ;
void Cblacs_setup( int * , int * ) ;
void Cblacs_get( int , int , int * ) ;
void Cblacs_gridinit( int * , char * , int , int ) ;
void Cblacs_gridinfo( int , int * , int * , int * , int * ) ;
void Cblacs_exit( int ) ;
void Cblacs_gridexit(int ConTxt);
int numroc_(int*, int*, int*, int*, int*);
void descinit_( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* ); 
int indxl2g_(int*, int*, int*, int*, int*);

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, double complex* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", creal(a[i+j*lda]) );
                printf( "\n" );
        }
}
void process_params(double *params, double kF, double *mu, size_t extra_data_size, void *extra_data);
void modify_densities(int it, double *h_densities, double *params, size_t extra_data_size, void *extra_data);
void modify_potentials(int it, double *h_densities, double *h_potentials, double *params, size_t extra_data_size, void *extra_data);
size_t get_extra_data_size(double *params);
int load_extra_data(size_t size, void *extra_data, double *params);

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


// ====================================================================================
// =================================== GLOBAL VARIABLES ===============================
// ====================================================================================
// This section reporduces variables keeped in GPU constant memory 
double *dc_params; /* Declaration of the variable */
size_t dc_extra_data_size;
void *dc_extra_data;
double dc_mu_a;
double dc_mu_b;
double dc_mu_a_old;		// for Broyden
double dc_mu_b_old;		// for Broyden
double dc_ec;

// rotating frame
double dc_Omega_a;
double dc_Omega_b;

// BdG mode
double aBdG;

typedef char * string;

// make -f Makefile.kzsolver

int main( int argc , char ** argv ) 
{    
    int i, j, k; // basic iterators
    int ix, iy, iz, ixyz; // lattice iterators
    int ierr; // error flag
    int iam, np; // basic MPI indicators
    int nwf; // number of wave-functions
    int nwfip=0; // number of wave-functions per process
    int iwf; // wave-function iterator
    int ikz; // index of kz vector
    // other technical variables
    int *wf_tbl, *wf_idx_tbl; // table of size np, keeps number of managed wf by each process
    int kziter=0; // local iteration number
    int it=0; // global iteration number
    double beta; // inverse of temperature
    double Lz_a=0.0, Lz_b=0.0, Lz=0.0;
    double Lz_a_old=0.0, Lz_b_old=0.0, Lz_old=0.0;
    double vextja=0.0, vextjb=0.0,  vextja_old=0.0, vextjb_old=0.0;
    
    double eF_a, eF_b, eF, Effg, kF;
    double mu[2]; // chemical potential
    double ec; // energy cut-off
    double rt_zheev, rt_dens, rt_pot, rt_other, rt_me, rt_tot=0.0, rt_redistrib; // run time
    
    // arrays
    double complex *h_wavefun; // pointer to wave-functions on host (cpu) side 
    double *h_fbetaEn; // pointer to weights of wave-functions on host (cpu) side
    double *h_En; // pointer to eigen-values of wave-functions on host (cpu) side
    double *h_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *h_densities_old; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *h_densities_partial; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *h_potentials; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    double *h_potentials_old; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    double *h_energy; // buffer for energies (CPU)
    double energy[5], energy_old[5];
    string energy_labels[5];
    energy_labels[0] = "E_kin";
    energy_labels[1] = "E_pot";
    energy_labels[2] = "E_pair";
    energy_labels[3] = "E_CM";
    energy_labels[4] = "E_ext";
    double E_tot, E_tot_old;
    double npart[2], npart_old[2], nparttest;
    char convstatus[2][6]; sprintf(convstatus[0], "FAIL"); sprintf(convstatus[1], "PASS");
    int is_converged, is_converged_local;
    int saving_iteration=0;
    dc_Omega_a=0.0;
    dc_Omega_b=0.0;
    
    char file_name[256];
    
    double *extra_data = NULL;
    size_t extra_data_size;
  
    /* start main */
    MPI_Init( &argc , &argv ) ; /* set up the parallel WORLD */
    MPI_Comm_size( MPI_COMM_WORLD , &np ) ; /* total number of processes */
    MPI_Comm_rank( MPI_COMM_WORLD , &iam ) ; /* id of process st 0 <= iam < np */
    
    // initial memory allocation
    cppmallocl( wf_tbl,np,int);
    cppmallocl( wf_idx_tbl,np,int);
    
    // Read of input parameters
    char execcmd[ 256 ] ;
    strcpy( execcmd , argv[ 0 ] ) ;
    for( i = 1 ; i < argc ; i++ ) 
    {
        strcat( execcmd , " " ) ; 
        strcat( execcmd , argv[ i ] ) ;
    }
    
    if( iam == 0 )
    {
        i = readcmd( argc , argv ) ;
        if( i == -1 )
        {
            printf( "TERMINATING! NO INPUT FILE.\n" ) ;
            ierr = -1 ;
            MPI_Abort( MPI_COMM_WORLD , ierr ) ;
            return( EXIT_FAILURE ) ;
        }
        
        // Read input file
        // Info from file is loaded into metadata structure
        j = parse_input_file(argv[i]);
        if ( j == 0 )
        {
            ierr = -1 ;
            printf("PROBLEM WITH INPUT FILE: `%s`.\n" , argv[ i ] ) ;
            MPI_Abort( MPI_COMM_WORLD , ierr ) ;
            return( EXIT_FAILURE ) ;      
        }
        
    }
    
    // Broadcast input parameter
    MPI_Bcast( &md , sizeof(md) , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
    if(iam==0) printf("# LATTICE: %d x %d x %d\n", NX, NY, NZ);
#ifdef USE_SCALAPACK_PZHEEVR
    if(iam==0) printf("# USING SCALAPACK WITH PZHEEVR.\n");
#endif
#ifdef USE_SCALAPACK_PZHEEVD
    if(iam==0) printf("# USING SCALAPACK WITH PZHEEVD.\n");
#endif
#ifdef USE_SCALAPACK_PZHEEV
    if(iam==0) printf("# USING SCALAPACK WITH PZHEEV.\n");
#endif
#ifdef USE_ELPA
    if(iam==0) printf("# USING ELPA.\n");
#endif
    
#if FUNCTIONAL==BDG
    aBdG = md.aBdG; // copy to global momeory
    if ( fabs(aBdG)<1.0e-12 )
    {
        ierr = -1 ;
        if(iam==0) printf("SET aBdG IN INPUT FILE!") ;
        MPI_Abort( MPI_COMM_WORLD , ierr ) ;
        return( EXIT_FAILURE ) ;      
    }
#else
    aBdG = 0.0; // deactivate BdG functional
#endif
    
    if(iam==0)
    {
        if(fabs(aBdG)<1.0e-12) printf("# ENERGY DENSITY FUNCTIONAL: ASLDA [UNITARITY]\n");
        else                   printf("# ENERGY DENSITY FUNCTIONAL: BdG [a=%16.8f]\n", aBdG);
    }
    
#ifdef SPINSYMMETRY_MODE
    md.spinsymmetry=1; // force spin symmetry mode
#endif
    if(md.spinsymmetry>0 && iam==0)  printf("# SPINSYMMETRY MODE IS ACTIVE.\n");
    if(md.nocurrents>0 && iam==0)  printf("# CODE IMPOSES NO CURRENTS FOR THE SOLUTION.\n");
    
#ifdef UNIFORM_TEST_MODE
    md.Na = ceil(1.0/(6.*M_PI*M_PI)*LXYZ);
    md.Nb = md.Na+1;
    if(md.spinsymmetry==1) md.Nb = md.Na;
    md.init0Na = md.Na;
    md.init0Nb = md.Nb;
    if(iam==0) printf("# UNIFORM_TEST_MODE: Setting number of particles to be: (%f,%f)\n", md.Na,md.Nb);
#endif
    
    // ====================================================================================
    // ==================================== BLACS GRID ====================================
    // ====================================================================================
    
    if(md.p==0 || md.q==0)
    { 
        if(iam==0) printf("# AUTOMATIC DIVISION OF WORK - MAY NOT BE OPTIMAL!\n"); fflush(stdout);
        int dims[2] = {0,0};
        MPI_Dims_create(np, 2, dims); // however, you can also set nprow and npcol by hand, keeping constraing nprow*npcol=np
        md.p = dims[0]; // cartesian direction 0
        md.q = dims[1]; // cartesian direction 1        
    }
    
    // for hamiltonian diagonalization
    int iam_blacs, nprocs_blacs, ictxt;
    char * b_order ;
    int MONE = -1 , ZERO = 0 , ONE = 1;
    int DESCA[ 9 ];
    int p, q, ip, iq, nip, niq;
    int info;
    int Hsize = NXYZ*2; // size of hamiltonian matrix
    
    if(iam==0) printf("# HAMILTONIAN SIZE: %d x %d\n", Hsize, Hsize);
#ifdef MATRIX_IS_REAL
    if(iam==0) printf("# HAMILTONIAN TOTAL STORAGE: %.3fGB\n", 1.0*sizeof(double)*Hsize*Hsize/1024/1024/1024);
#else
    if(iam==0) printf("# HAMILTONIAN TOTAL STORAGE: %.3fGB\n", 1.0*sizeof(double complex)*Hsize*Hsize/1024/1024/1024);
#endif 
    
    /* initialize the BLACS grid for hamiltonian diagonalizaion- a virtual rectangular grid */
    if(iam==0) printf("# CREATING CBLACS GRID OF SIZE (pzheev): [%d x %d]\n", md.p, md.q);
    b_order = "R" ;
    Cblacs_pinfo( &iam_blacs , &nprocs_blacs ) ;
    if ( nprocs_blacs < 1 ) 
        Cblacs_setup( &iam_blacs , &nprocs_blacs ) ;
    Cblacs_get( MONE , ZERO , &ictxt ) ; 
    Cblacs_gridinit( &ictxt , b_order , md.p , md.q ) ;  /* 'Row-Major' */
    Cblacs_gridinfo( ictxt , &p , &q , &ip , &iq ) ; /* ip,iq: the process row,column id */ 
    if(p!=md.p || q!=md.q) error_msg_mpi_abort(iam, p!=md.p || q!=md.q);       
    nip = numroc_( &Hsize, &md.mb, &ip, &ZERO, &p );
    niq = numroc_( &Hsize, &md.nb, &iq, &ZERO, &q );
#ifdef VERBOSE
    printf("# CHECK: iam=%d, ip=%d, iq=%d, nip=%d, niq=%d mb=%d nb=%d\n", iam, ip, iq, nip, niq, md.mb, md.nb);    
#endif
    /* Descriptor for the matrix */
    descinit_( DESCA, &Hsize, &Hsize, &md.mb, &md.nb, &ZERO, &ZERO, &ictxt, &nip, &info ); 
    if(info!=0) error_msg_mpi_abort(iam, info!=0);      
    
    // write info to structure metadata_s3dpca_grid
    metadata_s3dpca_grid bgrid;
    bgrid.ip=ip; // grid identifier
    bgrid.iq=iq; // grid identifier
    bgrid.nip=nip; // size of local matrix
    bgrid.niq=niq; // size of local matrix
    bgrid.p=p; // grid size
    bgrid.q=q; // grid size
    bgrid.mb=md.mb; // block size
    bgrid.nb=md.nb; // block size
     
    
    // ====================================================================================
    // ================================ ALLOCATE CPU BUFFERS ==============================
    // ====================================================================================
    
    // ---------------- POTENTIALS & DENSITIES----------------
    cppmallocl(h_densities,12*NXYZ,double);
    cppmallocl(h_densities_old,12*NXYZ,double);
    cppmallocl(h_densities_partial,12*NXYZ,double);
    cppmallocl(h_potentials,4*NXYZ,double);
    cppmallocl(h_potentials_old,4*NXYZ,double);
    cppmallocl(dc_params, MAX_USER_PARAMS, double);
    
    // For easier access to data
    // densities 
    double *rho_a = (double *)(h_densities +  0*NXYZ);
    double *rho_b = (double *)(h_densities +  1*NXYZ);
    double *tau_a = (double *)(h_densities +  2*NXYZ);
    double *tau_b = (double *)(h_densities +  3*NXYZ);
    double complex *nu = (double complex *)(h_densities +  4*NXYZ);
    double *j_a_x = (double *)(h_densities +  6*NXYZ);
    double *j_a_y = (double *)(h_densities +  7*NXYZ);
    double *j_a_z = (double *)(h_densities +  8*NXYZ);
    double *j_b_x = (double *)(h_densities +  9*NXYZ);
    double *j_b_y = (double *)(h_densities + 10*NXYZ);
    double *j_b_z = (double *)(h_densities + 11*NXYZ);
    
    // pontentials
    double *V_a = (double *)(h_potentials +  0*NXYZ);
    double *V_b = (double *)(h_potentials +  1*NXYZ);
    double complex *delta = (double complex *)(h_potentials +  2*NXYZ);
    
    cppmallocl(h_fbetaEn, 2*NXYZ,double);
    cppmallocl(h_En, 2*NXYZ,double);

    // ---------------- HAMILTONIAN ----------------
    // hamitonian
    double complex *h;
    cppmallocl(h, nip*niq,double complex); // allocate only local fraction of the hamiltonian matrix
    // eigen-vectors
    double complex *U;
#ifdef MATRIX_IS_REAL
    U = h; // in place working mode
#else
    cppmallocl(U, nip*niq,double complex); // allocate only local fraction of the matrix
#endif
    // eigen-values
    double * En; 
    cppmallocl(En, 2*NXYZ,double); // allocate space for the whole vector
    
#ifdef MATRIX_IS_REAL
    // use only h matrix as working area
    double *hR = (double *)h;
    double *UR = hR + nip*niq;
#endif
    
    // ====================================================================================
    // ================================ CREATE WAVE VECTORS ===============================
    // ====================================================================================
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }
    
    // ===================================================================================
    // ======================================== UEXT =====================================
    // ===================================================================================    
    dc_mu_a=__md_pca_uniform.mu_a;
    dc_mu_b=__md_pca_uniform.mu_b;
    dc_ec = __md_pca_uniform.ec;
    beta = __md_pca_uniform.beta;
    
    // reset energy and particle number
    for(ixyz=0; ixyz< 5; ixyz++) energy[ixyz] = 0.0;
    npart[SPINA]=0.0; npart[SPINB]=0.0;
    
    // ===================================================================================
    // =============================== INITIAL STATE =====================================
    // ===================================================================================
    if(md.inittype==0 || md.inittype==1) // Start from uniform solution
    {
        if(md.inittype==0)
        {
            if(iam==0) printf("# CREATING UNIFORM SOLUTION...\n");

            // Generate uniform initial 
            if(fabs(aBdG)<1.0e-12) solve_uniform_problem    (md.init0Na/LXYZ, md.init0Nb/LXYZ, &nwf, iam==0);
            else                   solve_uniform_problem_bdg(md.init0Na/LXYZ, md.init0Nb/LXYZ, &nwf, iam==0);
           
            // Save solution
            if(iam==0 && md.init0save)
            {
                cpu_exec( save_uniform() );
//                 ABORT;
            }  
        }
        else 
        {
            if(iam==0) printf("# READING UNIFORM SOLUTION...\n");
            if(iam==0) { cpu_exec( read_uniform(&nwf, iam==0) ); }
            MPI_Bcast( &__md_pca_uniform , sizeof(metadata_pca_uniform_t), MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
            MPI_Bcast( &nwf , 1, MPI_INT , 0 , MPI_COMM_WORLD ) ;
        }
                       
        // initialize potentials and densities
        for(ixyz=0; ixyz<NXYZ; ixyz++)
        {
            // potentials
            V_a[ixyz]=__md_pca_uniform.V_a; 
            V_b[ixyz]=__md_pca_uniform.V_b; 
            delta[ixyz]=__md_pca_uniform.delta + I*0.0;
            
            // densities
            rho_a[ixyz]=__md_pca_uniform.n0_a; 
            rho_b[ixyz]=__md_pca_uniform.n0_b; 
            tau_a[ixyz]=__md_pca_uniform.tau_a; 
            tau_b[ixyz]=__md_pca_uniform.tau_b; 
            nu[ixyz]=__md_pca_uniform.nu; 
            j_a_x[ixyz]=0.0;
            j_a_y[ixyz]=0.0;
            j_a_z[ixyz]=0.0;
            j_b_x[ixyz]=0.0;
            j_b_y[ixyz]=0.0;
            j_b_z[ixyz]=0.0;           
        }
        
        // set eF, kF and Effg
        eF=pow(3.0*M_PI*M_PI*(__md_pca_uniform.n0_a+__md_pca_uniform.n0_b), 2.0/3.0) / 2.0;
        kF=pow(3.0*M_PI*M_PI*(__md_pca_uniform.n0_a+__md_pca_uniform.n0_b), 1.0/3.0);
        eF_a=pow(6.0*M_PI*M_PI*__md_pca_uniform.n0_a, 2.0/3.0) / 2.0;
        eF_b=pow(6.0*M_PI*M_PI*__md_pca_uniform.n0_b, 2.0/3.0) / 2.0;
        Effg = 0.6*__md_pca_uniform.n0_a*eF_a*LXYZ + 0.6*__md_pca_uniform.n0_b*eF_b*LXYZ; 
        
        // set global variables
        dc_mu_a=__md_pca_uniform.mu_a;
        dc_mu_b=__md_pca_uniform.mu_b;
        dc_ec = __md_pca_uniform.ec;
        beta = __md_pca_uniform.beta;
    }
    else if(md.inittype==2 || md.inittype==22) // start from checkpoint
    {
        if(iam==0)
        {
            sprintf(file_name, "%s_checkpoint.s3dpca", md.inprefix);
            printf("# READING CHECKPOINT FILE `%s`\n", file_name);
            FILE * pFile = fopen(file_name, "rb");
            
            // write all nescesary data to file
            fread(&it          , sizeof(int)         , 1 , pFile); // iteration number
            fread(&dc_mu_a     , sizeof(double)      , 1 , pFile); 
            fread(&dc_mu_b     , sizeof(double)      , 1 , pFile); 
            fread(&dc_ec       , sizeof(double)      , 1 , pFile); 
            fread(&beta        , sizeof(double)      , 1 , pFile); 
            fread(&eF          , sizeof(double)      , 1 , pFile); 
            fread(&kF          , sizeof(double)      , 1 , pFile);
            fread(&Effg        , sizeof(double)      , 1 , pFile);
            fread(h_potentials , sizeof(double)*NXYZ , 4 , pFile);
            fread(h_densities  , sizeof(double)*NXYZ , 12, pFile);
            fread(energy       , sizeof(double)      , 5 , pFile);
            fread(npart        , sizeof(double)      , 2 , pFile);
            fread(&dc_mu_a_old , sizeof(double)      , 1 , pFile); 
            fread(&dc_mu_b_old , sizeof(double)      , 1 , pFile);
                  
            fclose(pFile);
            
            printf("# CHECKPOINT READ: it=%d\n", it);
            printf("# CHECKPOINT READ: dc_mu_a=%16.8g  dc_mu_b=%16.8g\n", dc_mu_a, dc_mu_b);
            printf("# CHECKPOINT READ: dc_ec=%16.8g  beta=%16.8g\n", dc_ec, beta);
            printf("# CHECKPOINT READ: eF=%16.8g  kF=%16.8g  Effg=%16.8g\n", eF, kF, Effg);
            printf("# CHECKPOINT READ: ------- NPART -------\n");
            printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
                "SPINA", npart[SPINA], npart[SPINA], (npart[SPINA]-npart[SPINA]));
            printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
                "SPINB", npart[SPINB], npart[SPINB], (npart[SPINB]-npart[SPINB]));  
            printf("# CHECKPOINT READ: ------- ENERGY -------\n");
            E_tot=0.0; E_tot_old=0.0;
            for(i=0; i<5; i++)
            {
                E_tot+=energy[i]; 
                
                printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
                    energy_labels[i], energy[i]/Effg, energy[i]/Effg, (energy[i]-energy[i])/Effg);
            }
            printf("  ------------------------------------------------------------------------\n");
            printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
                "E_tot", E_tot/Effg, E_tot/Effg, (E_tot-E_tot)/Effg);
        }
        
        // Send data to all processes
        MPI_Bcast(&it          , 1 , MPI_INT    , 0 , MPI_COMM_WORLD );
        MPI_Bcast(&dc_mu_a     , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&dc_mu_b     , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&dc_ec       , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&beta        , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&eF          , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&kF          , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(&Effg        , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(h_potentials , 4*NXYZ , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(h_densities  , 12*NXYZ, MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(energy       , 5 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(npart        , 2 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(&dc_mu_a_old , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&dc_mu_b_old , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    }
    else if(md.inittype==4) // start from cuGPE initial solutions
    {
        if(iam==0) 
        {
            printf("# INITTYPE[4]: DENSITIES FROM dpca FILES.\n");
            file_operation( inittype4_readdpca(md.inprefix, h_densities, h_potentials, &eF) );
            
            // recompute potentials
            kF=sqrt(2*eF);
            if(iam==0) printf("# EXECUTING: process_params(md.params, %f)\n", kF);
            for(i=0; i<MAX_USER_PARAMS; i++) dc_params[i]=md.params[i];
            mu[SPINA]=dc_mu_a; mu[SPINB]=dc_mu_b;
            process_params(dc_params, kF, mu, extra_data_size, extra_data);
            cpu_exec( recompute_potentials_meanfield_only(it, h_densities, h_potentials, h_potentials) ); 
            
            // set chemical potentials using TF approximation
            zbrent_data_t npart_mu;
            for(j=0; j<2; j++)
            {
                reset_zbrent_data(&npart_mu);
                double sol;
                double init_mu_1=-100;
                double init_mu_2= 100.;
                double *__U;
                double __N;
                if(j==SPINA) {__U = V_a; __N=md.Na;}
                else         {__U = V_b; __N=md.Nb;}
                
//                 for(ix=0; ix<NX; ix++) printf("UUU %d %f %f\n", ix, __U[(NZ/2)+NZ*(NY/2)+NZ*NY*ix], cabs(delta[(NZ/2)+NZ*(NY/2)+NZ*NY*ix]));
                
                while(1)
                {
                    sol=zbrent2(&npart_mu, &init_mu_1, &init_mu_2, 1.0e-12, &ierr);
        //             if(ip==0) printf("sol=%f, ierr=%d\n", sol, ierr);
                    if(ierr==0) break; // we converged
                    if(ierr>0) // request for point
                    {
                        double x = get_point_from_zbrent_data(&npart_mu);
        //                 if(ip==0) printf("zbrent needs value of function for x=%f\n", x);
                        double fx=npart_TF(NXYZ, __U, x, DXYZ);
//                         if(iam==0) printf("Setting: f(%f)=%f\n", x, fx);
                        fx=fx-__N;
//                         if(iam==0) printf("Setting: f(%f)=%f\n", x, fx);
                        set_point_in_zbrent_data(&npart_mu, fx);
                    }
                    if(ierr<0) ABORT_NOBARRIER;
                }
                mu[j]=sol;
            }
            
            dc_mu_a = mu[SPINA];
            dc_mu_b = mu[SPINB];
            dc_ec=M_PI*M_PI/(2.*DX*DX);
            Effg = 0.6 * md.Na * eF;
            beta = 1.0 / (md.kztemp * eF); 
        }
        
        MPI_Bcast(&it          , 1 , MPI_INT    , 0 , MPI_COMM_WORLD );
        MPI_Bcast(&dc_mu_a     , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&dc_mu_b     , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&dc_ec       , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&beta        , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&eF          , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&kF          , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(&Effg        , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(h_potentials , 4*NXYZ , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(h_densities  , 12*NXYZ, MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        

//         ABORT; // TODO
    }
    else
    {
        if(iam==0) printf("NOT SUPPORTED INITTYPE=%d!\n", md.inittype);
        ABORT;
    }
    
    // ===================================================================================
    // ================================== EXTRA DATA =====================================
    // ===================================================================================
    if(iam==0) extra_data_size = get_extra_data_size(md.params);
    MPI_Bcast( &extra_data_size , sizeof(size_t) , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
    if(extra_data_size>0)
    {
        if(iam==0) printf("# EXTRA_DATA IS ACTIVE.\n");
        if(iam==0) printf("# ALLOCATING EXTRA_DATA OF SIZE %ld B.\n", extra_data_size); fflush(stdout);
        if ( ( extra_data = (void *) malloc( extra_data_size ) ) == NULL  )
        {                                                             
            fprintf( stderr , "error: cannot malloc()! Exiting!\n") ; 
            fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; 
            MPI_Finalize() ;
            /* Arrays will be cleared automatically */
            return( EXIT_FAILURE ) ; 
        }
        
        if(iam==0) cpu_exec( load_extra_data(extra_data_size, extra_data, md.params) );
        MPI_Bcast( extra_data , extra_data_size , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
    }
    dc_extra_data_size=extra_data_size;
    dc_extra_data=extra_data;

    // ===================================================================================
    // ================================== FFTW PLANS =====================================
    // ===================================================================================
    if(iam==0) { printf("# CREATING FFTW PLANS...\n"); fflush(stdout); }
    metadata_s3dpca_fft mdfft; // keeps plans and buffers for fftw
    create_fft_plans(&mdfft, md.batch);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // ===================================================================================
    // ========================= MATRIX ELEMENTS OF GRADIENTS ============================
    // ===================================================================================    
    double complex * me_d_dx;
    double complex * me_d_dy;
    double complex * me_d_dz;
    cppmallocl(me_d_dx, NX*NX,double complex);
    cppmallocl(me_d_dy, NY*NY,double complex);
    cppmallocl(me_d_dz, NZ*NZ,double complex);
    
    cpu_exec( compute_matrix_elements_of_momentum_operator(NX, DX, me_d_dx) );
    cpu_exec( compute_matrix_elements_of_momentum_operator(NY, DY, me_d_dy) );
    cpu_exec( compute_matrix_elements_of_momentum_operator(NZ, DZ, me_d_dz) );
        
    // ===================================================================================
    // ======================================= LOGGER ====================================
    // =================================================================================== 
    
    // NOTE settings some variables
#ifndef UNIFORM_TEST_MODE
    if(md.referencekF>0.0) kF = md.referencekF;
    eF = 0.5*kF*kF;
    Effg = 0.6 * (md.Na+md.Nb) * eF;
    beta = 1.0 / (md.kztemp * eF);    
    if(md.ec>0.0) dc_ec = md.ec; 
    else          dc_ec = M_PI*M_PI/(2.*DX*DX);
#endif
    
    double dc_ec_l=-1.0*dc_ec; // lower bound for states extraction
    double dc_ec_u= 1.0*dc_ec; // upper bound for states exteraction
    if(md.spinsymmetry>0) dc_ec_l=0.0; // take only positive states
    
    // special case for saving
    if(md.writewf==1 && md.kzmaxiters==1) saving_iteration=1;
    
    // Create binary files and add initial measurement
    if(iam==0)
    {
        // Create empty files with headers - do it once
        sprintf(file_name, "%s_density_a.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, NZ, DX, DY, DZ, eF, 1.0*it, 1.0) );
        sprintf(file_name, "%s_density_b.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, NZ, DX, DY, DZ, eF, 1.0*it, 1.0) );
        sprintf(file_name, "%s_delta.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, NZ, DX, DY, DZ, eF, 1.0*it, 1.0) );    
        sprintf(file_name, "%s_current_a.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, NZ, DX, DY, DZ, eF, 1.0*it, 1.0) );
        sprintf(file_name, "%s_current_b.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, NZ, DX, DY, DZ, eF, 1.0*it, 1.0) );    
        
        // for each measurement add data to file
        sprintf(file_name, "%s_density_a.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, rho_a, sizeof(double)*NXYZ) );
        sprintf(file_name, "%s_density_b.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, rho_b, sizeof(double)*NXYZ) );
        sprintf(file_name, "%s_delta.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, delta, sizeof(double complex)*NXYZ) );
        sprintf(file_name, "%s_current_a.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, j_a_x, sizeof(double)*NXYZ*3) );
        sprintf(file_name, "%s_current_b.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, j_b_x, sizeof(double)*NXYZ*3) );
        
        // Create run log and add entry
        mu[SPINA]=dc_mu_a; mu[SPINB]=dc_mu_b;
        cpu_exec( create_header_of_runlog(execcmd, kF, Effg, mu, dc_ec, nwf, np, nwfip) );
    }
    
    // ====================================================================================
    // ====================================== WORKSPACE ===================================
    // ====================================================================================    
#if defined(USE_SCALAPACK_PZHEEVD) || defined(USE_SCALAPACK_PZHEEV) || defined(USE_SCALAPACK_PZHEEVR)
    /* Query and allocate the optimal workspace */
#ifdef USE_SCALAPACK_PZHEEVR
    if(iam==0) printf("# PREPARING WORKING BUFFERS FOR: `pzheevr`\n");
#endif
#ifdef USE_SCALAPACK_PZHEEVD
    if(iam==0) printf("# PREPARING WORKING BUFFERS FOR: `pzheevd`\n");
#endif
#ifdef USE_SCALAPACK_PZHEEV
    if(iam==0) printf("# PREPARING WORKING BUFFERS FOR: `pzheev`\n");
#endif    
#ifdef MATRIX_IS_REAL
    if(iam==0) printf("# MATRIX IS ASSUMED TO BE REAL! SWITCHING TO pdsyev* VERSION!\n");
#endif
    double complex * work , tw[ 2 ] = {1.0+I*0.0,1.0+I*0.0};
    double * rwork , tw_[2] = {1.0, 1.0} ;
    int lwork = -1, lrwork = -1, liwork = 7 * Hsize + 8 * q + 2 ; 
    int * iwork;

#ifdef USE_SCALAPACK_PZHEEVR
    double pzheevr_vl=dc_ec_l;
    double pzheevr_vu=dc_ec_u;

    int pzheevr_il=1;     // not referenced in this context
    int pzheevr_iu=Hsize; // not referenced in this context
    int pzheevr_m=-1;
    int pzheevr_nz=-1;
    cppmallocl(iwork, liwork, int);
#ifdef MATRIX_IS_REAL
    pdsyevr_( "V", "V", "U", &Hsize,hR, &ONE , &ONE , DESCA, &pzheevr_vl, &pzheevr_vu, &pzheevr_il, &pzheevr_iu, &pzheevr_m, &pzheevr_nz,
              En,UR, &ONE , &ONE, DESCA            , tw_, &lrwork, iwork, &liwork, &info );
#else
    pzheevr_( "V", "V", "U", &Hsize, h, &ONE , &ONE , DESCA, &pzheevr_vl, &pzheevr_vu, &pzheevr_il, &pzheevr_iu, &pzheevr_m, &pzheevr_nz,
              En, U, &ONE , &ONE, DESCA, tw, &lwork, tw_, &lrwork, iwork, &liwork, &info );
#endif
    liwork = iwork[ 0 ] ;
    free( iwork ) ;
    cppmallocl(iwork, liwork, int);
#endif
#ifdef USE_SCALAPACK_PZHEEVD
    cppmallocl(iwork, liwork, int);
#ifdef MATRIX_IS_REAL
    pdsyevd_( "V", "U", &Hsize,hR, &ONE , &ONE , DESCA, En,UR, &ONE , &ONE, DESCA            , tw_, &lrwork, iwork, &liwork, &info );
#else
    pzheevd_( "V", "U", &Hsize, h, &ONE , &ONE , DESCA, En, U, &ONE , &ONE, DESCA, tw, &lwork, tw_, &lrwork, iwork, &liwork, &info );
#endif
    liwork = iwork[ 0 ] ;
    free( iwork ) ;
    cppmallocl(iwork, liwork, int);
#endif
#ifdef USE_SCALAPACK_PZHEEV
    pzheev_( "V", "U", &Hsize, h, &ONE , &ONE , DESCA, En, U, &ONE , &ONE, DESCA, tw, &lwork, tw_, &lrwork, &info );
#endif
    if(info!=0) error_msg_mpi_abort(iam, info!=0);   

    lwork = ( int ) creal( tw[ 0 ] ) ;
    lrwork = ( int ) creal( tw_[ 0 ] ) ;
    
    cppmallocl(work, lwork, double complex);
    cppmallocl(rwork, lrwork, double);
    
#ifdef VERBOSE 
    printf("# pzheevd-set[%d]: info=%d  lwork=%d  lrwork=%d  liwork=%d\n",iam, info, lwork, lrwork, liwork);
    size_t size_workspace = sizeof(double complex)*lwork+sizeof(double)*lrwork;
    printf("# SIZE OF WORKSPACE [%d]: size_workspace=%.3fGB\n",iam, 1.0*(double)(size_workspace)/1024./1024./1024.);  
#endif
 
#endif
    
#ifdef USE_ELPA
    if(iam==0) printf("# SETTING UP ELPA...\n");
    if (elpa_init(20181112) != ELPA_OK) error_msg_mpi_abort(iam, ELPA API version not supported);
    
    elpa_t handle;
    handle = elpa_allocate(&info); 
    if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    
    /* Set parameters */
    elpa_set(handle, "na", Hsize, &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    int elpa_nev = (int)(USE_ELPA_NEV_FRACTION*Hsize);
    elpa_set(handle, "nev", elpa_nev, &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    elpa_set(handle, "local_nrows", nip, &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    elpa_set(handle, "local_ncols", niq, &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    elpa_set(handle, "nblk", md.nb, &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    elpa_set(handle, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    elpa_set(handle, "process_row",ip, &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    elpa_set(handle, "process_col", iq, &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);

    /* Setup */
    info=elpa_setup(handle);    
    if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);    
    
    elpa_set(handle, "solver", ELPA_SOLVER_2STAGE, &info);
    if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    if(iam==0) printf("# SETTING UP OF ELPA  DONE.\n");
#endif
        
    fflush(stdout);
         
    // ===================================================================================
    // ========================== SELF-CONSITENT LOOP ====================================
    // ===================================================================================     
    if(md.resetit) it=0; // reset iterator counter
    
    // special case - only one iteration for diagonalization
    if(md.kzmaxiters==1 && md.writewf==1) saving_iteration=1;
    
    while(1) // do until reached self-consitency
    {
        b_t();
        if((kziter+1)==md.kzmaxiters && md.writewf==1) 
        {
            if(iam==0) printf("# EXECUTING LAST ITERATION WITH SAVING DATA [md.writewf==1]\n");
            saving_iteration=1;
        }
        rt_zheev=0.0; rt_dens=0.0; rt_pot=0.0; rt_other=0.0; rt_me=0.0; rt_redistrib=0.0;
        b_t();
        // Make copy of potentials and densities
        for(ixyz=0; ixyz< 4*NXYZ; ixyz++) h_potentials_old[ixyz] = h_potentials[ixyz];
        for(ixyz=0; ixyz<12*NXYZ; ixyz++) h_densities_old[ixyz]  = h_densities[ixyz];
        for(ixyz=0; ixyz<12*NXYZ; ixyz++) h_densities_partial[ixyz] = 0.0; // reset
        for(i=0; i< 5; i++) energy_old[i] = energy[i]; // make copy
        npart_old[SPINA]=npart[SPINA]; npart_old[SPINB]=npart[SPINB]; // make copy 
        Lz_a_old=Lz_a; Lz_b_old=Lz_b; Lz_old=Lz; 
        vextja_old=vextja; vextjb_old=vextjb; 
        dc_mu_a_old = dc_mu_a; dc_mu_b_old = dc_mu_b;
        
        if(md.referencekF>0.0) kF = md.referencekF;
        else                   kF = pow(6.*M_PI*M_PI*rho_a[NZ/2 + NZ*NY/2 + NZ*NY*NX/2],1./3.); 
                               // take density in the center and use it for definition of the kF (for SPINA)
        
        eF = 0.5 * kF * kF;
        beta = 1.0 / (md.kztemp * eF);
        if(iam==0) printf("# EXECUTING: process_params(md.params, %f)\n", kF);
        for(i=0; i<MAX_USER_PARAMS; i++) dc_params[i]=md.params[i];
        mu[SPINA]=dc_mu_a; mu[SPINB]=dc_mu_b;
        process_params(dc_params, kF, mu, extra_data_size, extra_data);
        ECHOLINE;
        rt_other+=e_t(0);
        
        // matrix elements
        b_t();
        cpu_exec( compute_matrix_elements(&bgrid, it, h_densities, h_potentials, &mdfft, h, me_d_dx, me_d_dy, me_d_dz) );
        rt_me+=e_t(0);
        ECHOLINE;
        
        // diagonalize
        b_t();
        if(iam==0) printf("# DIAGONALIZATION %d...\n", it); fflush(stdout);
#ifdef MATRIX_IS_REAL
        // convert matrix to real version
//         for(ixyz=0; ixyz<nip*niq; ixyz++) if(fabs(cimag(h[ixyz]))>1.0e-12) {printf("# ERROR: matrix has imginary components!\n"); ABORT_NOBARRIER;}
        for(ixyz=0; ixyz<nip*niq; ixyz++) hR[ixyz] = creal(h[ixyz]);
#endif
#ifdef USE_SCALAPACK_PZHEEVR
        pzheevr_m=-1;
        pzheevr_nz=-1;
#ifdef MATRIX_IS_REAL
        pdsyevr_( "V", "V", "U", &Hsize,hR, &ONE , &ONE , DESCA, &pzheevr_vl, &pzheevr_vu, &pzheevr_il, &pzheevr_iu, &pzheevr_m, &pzheevr_nz,
                En,UR, &ONE , &ONE, DESCA              , rwork, &lrwork, iwork, &liwork, &info );
#else
        pzheevr_( "V", "V", "U", &Hsize, h, &ONE , &ONE , DESCA, &pzheevr_vl, &pzheevr_vu, &pzheevr_il, &pzheevr_iu, &pzheevr_m, &pzheevr_nz,
                En, U, &ONE , &ONE, DESCA, work, &lwork, rwork, &lrwork, iwork, &liwork, &info );
#endif
#endif
#ifdef USE_SCALAPACK_PZHEEVD
#ifdef MATRIX_IS_REAL
        pdsyevd_( "V", "U", &Hsize,hR, &ONE , &ONE , DESCA, En,UR, &ONE , &ONE, DESCA              , rwork, &lrwork, iwork, &liwork, &info );
#else
        pzheevd_( "V", "U", &Hsize, h, &ONE , &ONE , DESCA, En, U, &ONE , &ONE, DESCA, work, &lwork, rwork, &lrwork, iwork, &liwork, &info );
#endif
#endif
#ifdef USE_SCALAPACK_PZHEEV
#ifdef MATRIX_IS_REAL
        NOT IMPLEMENTED!
#else
        pzheev_( "V", "U", &Hsize, h, &ONE , &ONE , DESCA, En, U, &ONE , &ONE, DESCA, work, &lwork, rwork, &lrwork, &info );
#endif
#endif
#if defined(USE_SCALAPACK_PZHEEVD) || defined(USE_SCALAPACK_PZHEEV) || defined(USE_SCALAPACK_PZHEEVR)
        /* Check for convergence */
        if( info > 0 ) 
        {
            printf( "The algorithm failed to compute eigenvalues!\n" );
            ABORT_NOBARRIER;
        }
#endif
#ifdef USE_ELPA
        elpa_eigenvectors(handle, h, En, U, &info);
        if( info !=ELPA_OK ) 
        {
            printf( "The algorithm failed to compute eigenvalues!\n" );
            ABORT_NOBARRIER;
        }     
        if(En[elpa_nev-1]<dc_ec)
        {
            printf( "ELPA: fraction of extracted values is too small!\n" );
            ABORT_NOBARRIER;
        }
#endif
#ifdef MATRIX_IS_REAL
        // convert result back to complex
        for(ixyz=0; ixyz<nip*niq; ixyz++) hR[ixyz] = UR[ixyz]; // back to hR matrix
        for(ixyz=nip*niq-1; ixyz>=0; ixyz--) U[ixyz] = hR[ixyz] + I*0.0;
#endif
        rt_zheev+=e_t(0);
        if(iam==0) printf("# DIAGONALIZATION DONE [%.0f sec]\n", rt_zheev); fflush(stdout);
        ECHOLINE;
        
        // preparation for density computation
        b_t();
#ifdef USE_SCALAPACK_PZHEEVR
        nwf=pzheevr_m;
#else
        // do matrix redistribution for computation of densities
        for(i=0; i<Hsize-1; i++) if(En[i]>En[i+1]) 
        {
            printf( "Eigenvalues are not sorted correctly!\n" );
            ABORT_NOBARRIER;
        }
        // find min and max index energies in the interval E in [-ecut,+ecut]
        ix=-1; iy=Hsize;
        for(i=0; i<Hsize; i++)
        {
            if(En[i]<dc_ec_l) ix=MAX(i,ix);
            if(En[i]>dc_ec_u) iy=MIN(i,iy);
        }
        ix=ix+1; // shift to next eigenvalue
        nwf=iy-ix;
#endif            
        if(iam==0) printf("# Number of nwf in [-ecut,+ecut] to be extracted is: %d (%.1f%% of total number of states)\n", nwf, 100.0*nwf/Hsize);
        
        // Temporary grid for density computation
        int ictxt_d; // context for density computation
        int p_d, q_d, ip_d, iq_d, nip_d, niq_d;
        int DESCUD[ 9 ];
        Cblacs_get( MONE , ZERO , &ictxt_d ) ; 
        Cblacs_gridinit( &ictxt_d , b_order , 1 ,  np) ;  /* 'Row-Major' */
        Cblacs_gridinfo( ictxt_d , &p_d , &q_d , &ip_d , &iq_d ) ; /* ip,iq: the process row,column id */
        nip_d = numroc_( &Hsize, &md.mb, &ip_d, &ZERO, &p_d );
        niq_d = numroc_( &nwf, &ONE, &iq_d, &ZERO, &q_d );
#ifdef VERBOSE
        printf("# CHECK-D: iam=%d, ip_d=%d, iq_d=%d, nip_d=%d, niq_d=%d\n", iam, ip_d, iq_d, nip_d, niq_d);    
#endif
        /* Descriptor for the matrix */
        descinit_( DESCUD, &Hsize, &nwf, &md.mb, &ONE, &ZERO, &ZERO, &ictxt_d, &nip_d, &info ); 
        if(info!=0) error_msg_mpi_abort(iam, info!=0);
        
        // allocate memory for temporary matrix
        double complex *U_d; // eigen-vectors for density computation 
        double *En_d; // eigen-energies restricted to [-ecut,+ecut]
        double *En_d_local;
        cppmallocl(U_d, nip_d*niq_d,double complex); // allocate only local fraction of the matrix
        cppmallocl(En_d, nwf,double); // allocate only for fraction of eigenvalues
        cppmallocl(En_d_local, niq_d,double); // allocate only for local fraction of eigenvalues
#ifdef USE_SCALAPACK_PZHEEVR
        for(i=0; i<nwf; i++) En_d[i]=En[i];
        
        for(i=0; i<niq_d; i++)
        {
            j=i+1;
            k=indxl2g_( &j, &ONE, &iq_d, &ZERO, &np )-1; // back to C standard
            En_d_local[i]=En_d[k];
        }
        
        Cpzgemr2d(Hsize, nwf, U, 1, 1, DESCA, U_d, 1, 1, DESCUD, ictxt);
#else
        for(i=0; i<nwf; i++) En_d[i]=En[ix+i];
        
        for(i=0; i<niq_d; i++)
        {
            j=i+1;
            k=indxl2g_( &j, &ONE, &iq_d, &ZERO, &np )-1; // back to C standard
            En_d_local[i]=En_d[k];
        }
        
        Cpzgemr2d(Hsize, nwf, U, 1, ix+1, DESCA, U_d, 1, 1, DESCUD, ictxt);
#endif
        rt_redistrib+=e_t(0);
        
        // --------- DATA SAVING ----------
        if(saving_iteration==1) // save extracted wave-functions
        {
            b_t();
            if(iam==0) printf("# WRITING WAVE-FUNCTIONS WHERE |E_n/eF|<%.6g\n", md.writeecut); fflush(stdout);
            if(iam==0) printf("# WRITING USING %d I/O GROUPS\n", md.iogroups); fflush(stdout);
            
            // Creating IO groups
            int idgroup; // identifier of the group
            MPI_Comm mpi_comm_group;
            int gr_iam, gr_np;
            
            idgroup=iam % md.iogroups;
            MPI_Comm_split(MPI_COMM_WORLD, idgroup, iam, &mpi_comm_group);
            MPI_Comm_rank(mpi_comm_group, &gr_iam);
            MPI_Comm_size(mpi_comm_group, &gr_np);
            if(gr_iam==0) printf("# I/O GROUP %d WITH %d PROCESSES HAS BEEN SUCCESSFULLY CREATED.\n", idgroup, gr_np); fflush(stdout);
            
            // !!!!! I/O can be written better with MPI I/O - maybe it will be improved in future !!!!!
            MPI_Status MPIStat;
            int lastwf=0;
        
            // ----- part 0 ------
            if(gr_iam==0)
            {
                // Create empty files
                sprintf(file_name, "%s_s3dpca.%04d.wfu", md.outprefix, idgroup);
                file_operation( touch_file(file_name) );
            
                sprintf(file_name, "%s_s3dpca.%04d.wfv", md.outprefix, idgroup);
                file_operation( touch_file(file_name) );

                sprintf(file_name, "%s_s3dpca.%04d.en", md.outprefix, idgroup);
                file_operation( touch_file(file_name) );
            }
            
            // ----- part 1 ------
            // wait for my turn of writing event
            if(gr_iam!=0) // Wait until gr_iam-1 process finish his job
                j = MPI_Recv(&lastwf, 1, MPI_INT, gr_iam-1, 99, mpi_comm_group, &MPIStat);    
                
            // save my wave-functions if generated
            j=0;
            file_operation( append_wf_from_s3dpca_part1(md.outprefix, En_d_local, U_d, md.writeecut*eF, beta, idgroup, niq_d, &j) );
            lastwf+=j;
                
            if(gr_iam!=(gr_np-1)) // File is free, send info to next process
                j = MPI_Send(&lastwf, 1, MPI_INT, gr_iam+1, 99, mpi_comm_group); 
                
            // ----- part 2 ------
            // wait for my turn of writing event
            if(gr_iam!=0) // Wait until gr_iam-1 process finish his job
                j = MPI_Recv(&j, 1, MPI_INT, gr_iam-1, 99, mpi_comm_group, &MPIStat);    
            
            // save my wave-functions if generated
            file_operation( append_wf_from_s3dpca_part2(md.outprefix, En_d_local, U_d, md.writeecut*eF, beta, idgroup, niq_d, &j) );
            
            if(gr_iam!=(gr_np-1)) // File is free, send info to next process
                j = MPI_Send(&j, 1, MPI_INT, gr_iam+1, 99, mpi_comm_group);    
            
            // ----- part 3 ------
            // wait for my turn of writing event
            if(gr_iam!=0) // Wait until gr_iam-1 process finish his job
                j = MPI_Recv(&j, 1, MPI_INT, gr_iam-1, 99, mpi_comm_group, &MPIStat);    
            
            // save my wave-functions if generated
            file_operation( append_wf_from_s3dpca_part3(md.outprefix, En_d_local, U_d, md.writeecut*eF, beta, idgroup, niq_d, &j) );
            
            if(gr_iam!=(gr_np-1)) // File is free, send info to next process
                j = MPI_Send(&j, 1, MPI_INT, gr_iam+1, 99, mpi_comm_group);  
                
            // finilize I/O
            MPI_Bcast( &lastwf , 1, MPI_INT , gr_np-1 , mpi_comm_group ) ;
            
            // create info file by each group
            sprintf(file_name, "%s_s3dpca.%04d.info", md.outprefix, idgroup);
            mu[SPINA] = dc_mu_a; mu[SPINB] = dc_mu_b;
            if(gr_iam==0) file_operation( create_checkpoint_info_pca(file_name, lastwf, NX, NY, NZ, DX, DY, DZ, kF, mu, md.writeecut*eF, beta) );
                
            double rt = e_t(0);
            if(gr_iam==0) printf("# DATA WRITING BY I/O GROUP %d TOOK %.1f SEC. WRITTEN %.2fMB. WRITTEN STATES=%d\n", idgroup, rt, (double)1.*lastwf*NX*NY*NZ*2*16/1024./1024., lastwf); fflush(stdout);
            
            // recompute total nwf
            if(gr_iam!=0) lastwf=0;
            MPI_Allreduce( &lastwf, &nwf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if(iam==0) printf("# DATA SAVING DONE. IN TOTAL WRITTEN %d WAVE-FUNCTIONS.\n", nwf);
            rt_other+=rt;
        }
        
        // --------- DENSITIES ----------
        b_t();
        // compute contribution to the densities
        cpu_exec( compute_contribution_to_densities(En_d_local, U_d, niq_d, beta, h_densities_partial, &mdfft, md.spinsymmetry) );
        // compute densities as global reduction
        MPI_Allreduce( h_densities_partial, h_densities, 12*NXYZ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        if(md.spinsymmetry==1)  for(ixyz=0; ixyz<NXYZ; ixyz++) // impose by hand symmetry on densities
        {
            rho_a[ixyz]=rho_b[ixyz];
            tau_a[ixyz]=tau_b[ixyz];
            j_a_x[ixyz]=j_b_x[ixyz];
            j_a_y[ixyz]=j_b_y[ixyz];
            j_a_z[ixyz]=j_b_z[ixyz];
        }
        if(md.nocurrents) for(ixyz=0; ixyz<NXYZ; ixyz++) // impose by hand no currents
        {
            j_a_x[ixyz]=0.0; j_b_x[ixyz]=0.0;
            j_a_y[ixyz]=0.0; j_b_y[ixyz]=0.0;
            j_a_z[ixyz]=0.0; j_b_z[ixyz]=0.0;
        }
        
#ifndef TAU_COMPUTATION_VIA_GRADIENTS  
        // finalize computation of tau
        cpu_exec( density_caculate_tau(h_densities, &mdfft) );
#endif
        rt_dens+=e_t(0);
                
        // free temporary resources
        b_t();
        free(U_d); 
        free(En_d);
        free(En_d_local);
        Cblacs_gridexit( ictxt_d );
        rt_other+=e_t(0);

        
        // ------------------ mix densities ------------------
        b_t();
        if(md.inittype==22 && kziter==0) //special case - started from interpolated checkpoint
        {
            // some densities like tau and nu are cut-off dependent, and may be very different 
            // when moving onle from to another lattice
            // what matters is only E_kin+E_pair which is well defined
            
            // pass - do not mix
            if(iam==0) printf("# SPECIAL CASE: START FROM INTERPOLATED SOLUTION [md.inittype==22]! MIXING SKIPPED!\n");
        }
        else if(saving_iteration==1) //special case - saving interation
        {
            // typically results from this iteration are loded into dynamical code
            // for clear comparision of read corretness skip mixing here
            
            // pass - do not mix
            if(iam==0) printf("# SPECIAL CASE: SAVING ITERATION! MIXING SKIPPED!\n");
        }
        else
        {
            for(ixyz=0; ixyz<12*NXYZ; ixyz++) h_densities[ixyz] = md.linearmixing * h_densities[ixyz] + (1.0-md.linearmixing) * h_densities_old[ixyz];
        }
        // ------------------ update chemical potentials ------------------
        if(iam==0) printf("# MUCHANGE FROM: dc_mu_a=%16.8g  dc_mu_b=%16.8g\n", dc_mu_a, dc_mu_b);
        if(it>0 && saving_iteration==0) // skip upfating the potential is it is saving iteration
        {       
            npart[SPINA]=0.0; npart[SPINB]=0.0;
            for(ixyz=0; ixyz<NXYZ; ixyz++) {npart[SPINA]+=rho_a[ixyz]; npart[SPINB]+=rho_b[ixyz];}
            npart[SPINA]*=DXYZ; npart[SPINB]*=DXYZ; 
            double kzmuchange_a = md.kzmuchange*(npart[SPINA] - md.Na)/md.Na;
            double kzmuchange_b = md.kzmuchange*(npart[SPINB] - md.Nb)/md.Nb;
            if(fabs(kzmuchange_a)>md.mumaxchange*eF)
            {
                if(kzmuchange_a>0.0) kzmuchange_a=     md.mumaxchange*eF;
                else                 kzmuchange_a=-1.0*md.mumaxchange*eF;
            }
            if(fabs(kzmuchange_b)>md.mumaxchange*eF)
            {
                if(kzmuchange_b>0.0) kzmuchange_b=     md.mumaxchange*eF;
                else                 kzmuchange_b=-1.0*md.mumaxchange*eF;
            }
            dc_mu_a -= kzmuchange_a;
            dc_mu_b -= kzmuchange_b;  
            if(md.spinsymmetry==1) dc_mu_b=dc_mu_a; // activate constraint
        }
        if(iam==0) printf("# MUCHANGE TO  : dc_mu_a=%16.8g  dc_mu_b=%16.8g\n", dc_mu_a, dc_mu_b);
        rt_other+=e_t(0);
        
        modify_densities(it, h_densities, dc_params, extra_data_size, extra_data) ;
        
        // ------------------ compute new potentials ------------------
        b_t();
        cpu_exec( recompute_potentials(it, h_densities, h_potentials, h_potentials) ); 
        modify_potentials(it, h_densities, h_potentials, dc_params, extra_data_size, extra_data) ;
        rt_pot+=e_t(0);
        
        // ------------------ angular momentum ------------------
        b_t();
        cpu_exec( compute_angular_momentum_Lz(j_a_x, j_a_y, &Lz_a) );
        cpu_exec( compute_angular_momentum_Lz(j_b_x, j_b_y, &Lz_b) );
        Lz = Lz_a + Lz_b; // total angular momentum        
        if(iam==0) printf("# ANGULAR MOMENTUM: it=%d\n", it);
        if(iam==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "LZ_A", Lz_a/npart[SPINA], Lz_a_old/npart_old[SPINA], (Lz_a/npart[SPINA]-Lz_a_old/npart_old[SPINA]));
        if(iam==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "LZ_B", Lz_b/npart[SPINB], Lz_b_old/npart_old[SPINB], (Lz_b/npart[SPINB]-Lz_b_old/npart_old[SPINB]));  
        if(iam==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "LZ", Lz/(npart[SPINA]+npart[SPINB]), Lz_old/(npart_old[SPINA]+npart_old[SPINB]), (Lz/(npart[SPINA]+npart[SPINB])-Lz_old/(npart_old[SPINA]+npart_old[SPINB])));
        if(iam==0) printf("dc_Omega_a=%16.8g  dc_Omega_b=%16.8g\n",dc_Omega_a, dc_Omega_b);
        
        // ------------------ external velocity ------------------
        cpu_exec( compute_vext_dot_j(it, SPINA, j_a_x, j_a_y, j_a_z, &vextja) );
        cpu_exec( compute_vext_dot_j(it, SPINB, j_a_x, j_a_y, j_a_z, &vextjb) );  
        if(iam==0) printf("# EXTERNAL VELOCITY: it=%d\n", it);
        if(iam==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "VEXT_A*J", vextja, vextja_old, (vextja-vextja_old));    
        if(iam==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "VEXT_B*J", vextjb, vextjb_old, (vextjb-vextjb_old));
        
        // ------------------ check convergence ------------------
        is_converged=1;
        cpu_exec( compute_energy(it, h_densities, h_potentials, energy, npart) );
        if(iam==0) printf("# CONVERGENCE REPORT PARTICLE NUMBER: it=%d\n", it);
        nparttest=fabs(npart[SPINA]-md.Na)/(md.Na+md.Nb);
        if(nparttest>md.npartconveps) {is_converged=0; is_converged_local=0;} else {is_converged_local=1;}
        if(iam==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g CONVSTATUS=%6s NPARTCONV=%16.8g\n", 
            "SPINA", npart[SPINA], npart_old[SPINA], (npart[SPINA]-npart_old[SPINA]), convstatus[is_converged_local], nparttest);
        nparttest=fabs(npart[SPINB]-md.Nb)/(md.Na+md.Nb);
        if(nparttest>md.npartconveps) {is_converged=0; is_converged_local=0;} else {is_converged_local=1;}
        if(iam==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g CONVSTATUS=%6s NPARTCONV=%16.8g\n", 
            "SPINB", npart[SPINB], npart_old[SPINB], (npart[SPINB]-npart_old[SPINB]), convstatus[is_converged_local], nparttest);
            
        if(iam==0) printf("# CONVERGENCE REPORT ENERGY: it=%d\n", it);
        E_tot=0.0; E_tot_old=0.0;
        for(i=0; i<5; i++)
        {
            E_tot+=energy[i]; 
            E_tot_old+=energy_old[i];
            
            if(fabs((energy[i]-energy_old[i])/Effg)>md.energyconveps) {is_converged=0; is_converged_local=0;} else {is_converged_local=1;}
            
            if(iam==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g CONVSTATUS=%6s\n", 
                energy_labels[i], energy[i]/Effg, energy_old[i]/Effg, (energy[i]-energy_old[i])/Effg, convstatus[is_converged_local]);
            
        }
        if(iam==0) printf("  ------------------------------------------------------------------------\n");
        if(fabs((E_tot-E_tot_old)/Effg)>md.energyconveps) {is_converged=0; is_converged_local=0;} else {is_converged_local=1;}
        if(iam==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g CONVSTATUS=%6s\n", 
                "E_tot", E_tot/Effg, E_tot_old/Effg, (E_tot-E_tot_old)/Effg, convstatus[is_converged_local]);
        
        double minF_new = E_tot - dc_mu_a*npart[SPINA] - dc_mu_b*npart[SPINB] - dc_Omega_a*Lz_a - dc_Omega_b*Lz_b - vextja - vextjb;
        double minF_old = E_tot_old - dc_mu_a_old*npart_old[SPINA] - dc_mu_b_old*npart_old[SPINB] - dc_Omega_a*Lz_a_old - dc_Omega_b*Lz_b_old - vextja_old - vextjb_old;
        if(iam==0) printf("# MINIMIZATION FUNCTION: %16.8f\n", minF_new);
        if(iam==0) printf("# FUNCTION CHANGED BY: %16.8f\n", minF_new-minF_old);
        if(iam==0)
        {
            #define OUTPUT_ENTRIES 18
            double line_items[OUTPUT_ENTRIES]={     
                npart[SPINA], // 2
                npart[SPINB], // 3
                npart[SPINA]+npart[SPINB], // 4
                E_tot/Effg, // 5
                energy[0]/Effg, // 6
                energy[1]/Effg, // 7
                energy[2]/Effg, // 8
                energy[3]/Effg, // 9
                energy[4]/Effg, //10
                minF_new, // 11
                dc_mu_a, //12
                dc_mu_b, //13
                beta, // 14
                Lz_a/npart[SPINA], // 15
                Lz_b/npart[SPINB], // 16
                Lz/(npart[SPINA]+npart[SPINB]) ,  // 17
                vextja, // 18
                vextjb, // 19
            };
            cpu_exec( add_line_to_file(it, rt_tot, OUTPUT_ENTRIES, line_items) );
            
            // for each measurement add data to file
            sprintf(file_name, "%s_density_a.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, rho_a, sizeof(double)*NXYZ) );
            sprintf(file_name, "%s_density_b.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, rho_b, sizeof(double)*NXYZ) );
            sprintf(file_name, "%s_delta.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, delta, sizeof(double complex)*NXYZ) );
            sprintf(file_name, "%s_current_a.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, j_a_x, sizeof(double)*NXYZ*3) );
            sprintf(file_name, "%s_current_b.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, j_b_x, sizeof(double)*NXYZ*3) );
        }
                
        // checkpoint - only by iam==0
        if(md.checkpoint && iam==0)
        {
            sprintf(file_name, "%s_checkpoint.s3dpca", md.outprefix);
            printf("# CREATING CHECKPOINT FILE `%s`\n", file_name);
            FILE * pFile = fopen(file_name, "wb");
            
            // write all nescesary data to file
            i=it+1;
            fwrite(&i           , sizeof(int)         , 1 , pFile); // iteration number
            fwrite(&dc_mu_a     , sizeof(double)      , 1 , pFile); 
            fwrite(&dc_mu_b     , sizeof(double)      , 1 , pFile); 
            fwrite(&dc_ec       , sizeof(double)      , 1 , pFile); 
            fwrite(&beta        , sizeof(double)      , 1 , pFile); 
            fwrite(&eF          , sizeof(double)      , 1 , pFile); 
            fwrite(&kF          , sizeof(double)      , 1 , pFile);
            fwrite(&Effg        , sizeof(double)      , 1 , pFile);
            fwrite(h_potentials , sizeof(double)*NXYZ , 4 , pFile);
            fwrite(h_densities  , sizeof(double)*NXYZ , 12, pFile);
            fwrite(energy       , sizeof(double)      , 5 , pFile);
            fwrite(npart        , sizeof(double)      , 2 , pFile);
            fwrite(&dc_mu_a_old , sizeof(double)      , 1 , pFile); 
            fwrite(&dc_mu_b_old , sizeof(double)      , 1 , pFile);
                  
            fclose(pFile);
        }
        
        rt_other+=e_t(0);
        
        // ------------------ timing------------------
        rt_tot=rt_zheev+rt_dens+rt_pot+rt_other+rt_me+rt_redistrib; 
        if(iam==0) printf("# TIMING rt_tot=%8.2f: rt_zheev=%8.2f[%5.2f%%] rt_redistrib=%8.2f[%5.2f%%] rt_dens=%8.2f[%5.2f%%] rt_pot=%8.2f[%5.2f%%] rt_me=%8.2f[%5.2f%%] rt_other=%8.2f[%5.2f%%]\n", 
            rt_tot, rt_zheev, rt_zheev/rt_tot*100., rt_redistrib, rt_redistrib/rt_tot*100., rt_dens, rt_dens/rt_tot*100., rt_pot, rt_pot/rt_tot*100., rt_me, rt_me/rt_tot*100., rt_other, rt_other/rt_tot*100.);
        fflush(stdout);
        
        
        if(saving_iteration==1)
        {
            // Write missing files: info for the total set of wf and file with potentials
            if(iam==0)
            {
                // write info file
                sprintf(file_name, "%s_s3dpca.info", md.outprefix);
                mu[SPINA] = dc_mu_a; mu[SPINB] = dc_mu_b;
                file_operation( create_checkpoint_info_pca(file_name, nwf, NX, NY, NZ, DX, DY, DZ, kF, mu, dc_ec, beta) );
                
                // write potentials
                sprintf(file_name, "%s_s3dpca.pud", md.outprefix);
                file_operation( checkpoint_save_u_and_delta_kzpca(file_name, NX*NY*NZ, V_a, delta) );
            }
            
            // Create check.stamp
            if(iam==0)
            {
                // write check.stamp file
                sprintf(file_name, "%s_check.stamp", md.outprefix);
                printf("# CREATING CHECK STAMP FILE: `%s`\n",file_name);
                file_operation( touch_file(file_name) );
                file_operation( check_stamp_entry(file_name, 12, NXYZ, h_densities, 5, energy) );
            }
            
            if(iam==0) printf("# SAVING ITERATION DONE.\n");
            break;
        }
        
        if(is_converged && kziter>0)
        {
            if(iam==0) printf("# ALGORITHM CONVERGED!\n");
            
            if(md.writewf==0) break;
            else saving_iteration=1;
        }
        
        it++; // go to next iteration
        kziter++;
        if(kziter==md.kzmaxiters)
        {
            if(iam==0) printf("# MAXIMUM NUMBER OF ITERATIONS REACHED!\n"); fflush(stdout);
            
            if(md.writewf==0) break;
            else saving_iteration=1;
        }
           
//         break; //NOTE - for tests only
    } // while(1)

    /* messy exit here */
    destroy_fft_plans(&mdfft);
    
    MPI_Barrier( MPI_COMM_WORLD ) ;
    MPI_Finalize() ;
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
