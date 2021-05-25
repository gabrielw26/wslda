// This is code for solving stationary DFT equations for polarized cold atoms 
// in 3D, without any symmetry restrictions

// Authors:
// Gabriel Wlazlowski <gabriel.wlazlowski@pw.edu.pl>

// NOTE: Inspect for NOTE flags the code before running it!!!

// #define VERBOSE
// #define S3DDEBUG

#ifdef S3DDEBUG
#define ECHOLINE                                                                                                        \
    {                                                                                                                   \
        if(iam==0) {wprintf("# ECHOLINE REACHED LINE %d IN FILE %s\n", __LINE__ , __FILE__); fflush(stdout); }           \
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

#include "wdata.h"

#include "pca_settings.h"
#include "pca_macro.h"
#include "pca_utils.h"
#include "wslda_potdens.h"
#include "wslda_wavevectors.h"
#include "pca_io.h"
#include "s2dpca_edf.h"
#include "pca_uniform.h"
#include "pca_logger.h"
#include "s3dpca_fft.h"
#include "s3dpca_me.h"
#include "s3dpca_grid.h"
#include "s3dpca_densities.h"
#include "sxdpca_broyden.h"
#include "wslda_writevars.h"
#include "wslda_functionals.h"
#include "wslda_reproducibility.h"
#include "wslda_interpolation.h"
#include "wslda_resize.h"
#include "wslda_st_checkpoint.h"

#if DIAGONALIZATION_ROUTINE==PZHEEVR
#define USE_SCALAPACK_PZHEEVR
#elif DIAGONALIZATION_ROUTINE==PZHEEVD
#define USE_SCALAPACK_PZHEEVD
#elif DIAGONALIZATION_ROUTINE==ELPA
#define USE_ELPA
#else
#error "select DIAGONALIZATION ROUTINE in predifines.h"
#endif


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
        wprintf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) wprintf( " %6.2f", creal(a[i+j*lda]) );
                wprintf( "\n" );
        }
}
void process_params(double *params, double *kF, double *mu, size_t extra_data_size, void *extra_data);
void modify_densities(int it, wslda_density h_densities, double *params, size_t extra_data_size, void *extra_data);
void modify_potentials(int it, wslda_density h_densities, wslda_potential h_potentials, double *params, size_t extra_data_size, void *extra_data);
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

// BdG mode
double aBdG;

int wsldapid; // process id - global variable

#include "logger.h"

typedef char * string;

int main( int argc , char ** argv ) 
{    
    int i, j, k; // basic iterators
    int ix, iy, iz, ixyz; // lattice iterators
    int ierr; // error flag
    int iam, np, ip, iq ; // basic MPI indicators
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
    double S=0.0, S_old=0.0; // entropy
    
    double eF_a, eF_b, eF, Effg, kF;
    double mu[2]; // chemical potential
    double ec; // energy cut-off
    double rt_zheev, rt_dens, rt_pot, rt_other, rt_me, rt_tot=0.0, rt_redistrib; // run time
    
    // arrays
    double complex *h_wavefun; // pointer to wave-functions on host (cpu) side 
    double *h_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *h_densities_old; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *h_densities_partial; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *h_potentials; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    double *h_potentials_old; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    double *h_energy; // buffer for energies (CPU)
    double observables[WSLDAITEMS];
    double energy[ENERGYITEMS], energy_old[ENERGYITEMS];
    string energy_labels[ENERGYITEMS];
    energy_labels[EKIN] = "E_kin";
    energy_labels[EPOT] = "E_pot";
    energy_labels[EPAIR] = "E_pair";
    energy_labels[ECURRENT] = "E_curr";
    energy_labels[EPOTEXT] = "E_potext";
    energy_labels[EPAIREXT] = "E_pairext";
    energy_labels[EVELEXT] = "E_velext";
    double E_tot, E_tot_old;
    double npart[2]={NUMERICAL_ZERO, NUMERICAL_ZERO}, npart_old[2]={NUMERICAL_ZERO, NUMERICAL_ZERO}, nparttest;
    char convstatus[2][6]; sprintf(convstatus[0], "FAIL"); sprintf(convstatus[1], "PASS");
    int is_converged, is_converged_local;
    int saving_iteration=0;
    
    // parameters for Broyden method
	double omega_0 = md.omega0broyden;	// weight assigned to the error in the inverse Jacobian
	double omega_n = md.omeganbroyden;	// weight associated with each previous iteration
	double omega_k = md.omegakbroyden;	// 		---//---
    double **dens_in;		// pointer to array of arrays of densities
    double **dens_out;		// 		---//---
    
    char file_name[256];
    
    void *extra_data = NULL;
    size_t extra_data_size;
  
    /* start main */
    MPI_Init( &argc , &argv ) ; /* set up the parallel WORLD */
    MPI_Comm_size( MPI_COMM_WORLD , &np ) ; /* total number of processes */
    MPI_Comm_rank( MPI_COMM_WORLD , &iam ) ; /* id of process st 0 <= iam < np */
    wsldapid=iam; // save to global variable
    
    if(iam==0) wprintf("# START OF THE MAIN FUNCTION\n");
    
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
            wprintf( "TERMINATING! NO INPUT FILE.\n" ) ; something_to_cheer_you_up(stdout);
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
            wprintf("PROBLEM WITH INPUT FILE: `%s`.\n" , argv[ i ] ) ; something_to_cheer_you_up(stdout);
            MPI_Abort( MPI_COMM_WORLD , ierr ) ;
            return( EXIT_FAILURE ) ;      
        }
        
        // Make copy of input file
        ip=iam; file_operation( check_if_can_overwrite_files() ); // terminate if file exists, and input->overwrite==0
        sprintf(file_name, "%s_input.txt", md.outprefix);
        copy_input_file(argv[i],file_name) ; 
        assure_reproducibility(md.outprefix);
    }
    
    if(iam==0) wprintf("# CODE: ST-WSLDA-3D\n");
    if(iam==0) wprintf("# VERSION: %s\n", VERSION);
    
    // Broadcast input parameter
    MPI_Bcast( &md , sizeof(md) , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
    if(iam==0) wprintf("# LATTICE: %d x %d x %d\n", NX, NY, NZ);
    if(iam==0) wprintf("# SPACING: %f x %f x %f\n", DX, DY, DZ);
#ifdef USE_SCALAPACK_PZHEEVR
    if(iam==0) wprintf("# USING SCALAPACK WITH PZHEEVR.\n");
#endif
#ifdef USE_SCALAPACK_PZHEEVD
    if(iam==0) wprintf("# USING SCALAPACK WITH PZHEEVD.\n");
#endif
#ifdef USE_SCALAPACK_PZHEEV
    if(iam==0) wprintf("# USING SCALAPACK WITH PZHEEV.\n");
#endif
#ifdef USE_ELPA
    if(iam==0) wprintf("# USING ELPA.\n");
#endif
    
#if FUNCTIONAL==BDG
    aBdG = md.aBdG; // copy to global momeory
    if ( fabs(aBdG)<1.0e-12 )
    {
        ierr = -1 ;
        if(iam==0) wprintf("ERROR: SET aBdG IN INPUT FILE!\n"); 
        something_to_cheer_you_up_pid0(stdout);
        fflush(stdout);
        MPI_Abort( MPI_COMM_WORLD , ierr ) ;
        return( EXIT_FAILURE ) ;      
    }
#else
    aBdG = 0.0; // deactivate BdG functional
#endif

#if FUNCTIONAL==BDG
    if(iam==0) wprintf("# ENERGY DENSITY FUNCTIONAL: BDG\n");
#elif FUNCTIONAL==SLDA    
    if(iam==0) wprintf("# ENERGY DENSITY FUNCTIONAL: SLDA\n");
#elif FUNCTIONAL==ASLDA    
    if(iam==0) wprintf("# ENERGY DENSITY FUNCTIONAL: ASLDA\n");    
#elif FUNCTIONAL==CUSTOMEDF    
    if(iam==0) wprintf("# ENERGY DENSITY FUNCTIONAL: CUSTOMEDF\n"); 
#endif
    
#ifdef SPINSYMMETRY_MODE
    md.spinsymmetry=1; // force spin symmetry mode
#endif
    if(md.spinsymmetry>0 && iam==0)  wprintf("# SPINSYMMETRY MODE IS ACTIVE.\n");
    if(md.nocurrents>0 && iam==0)  wprintf("# CODE IMPOSES NO CURRENTS FOR THE SOLUTION.\n");
    
#ifdef UNIFORM_TEST_MODE
    md.Na = ceil(1.0/(6.*M_PI*M_PI)*LXYZ);
    md.Nb = md.Na+1;
    if(md.spinsymmetry==1) md.Nb = md.Na;
    md.init0Na = md.Na;
    md.init0Nb = md.Nb;
    if(iam==0) wprintf("# UNIFORM_TEST_MODE: Setting number of particles to be: (%f,%f)\n", md.Na,md.Nb);
#endif
    
    // ====================================================================================
    // ==================================== BLACS GRID ====================================
    // ====================================================================================
    
    if(md.p==0 || md.q==0)
    { 
        if(iam==0) wprintf("# AUTOMATIC DIVISION OF WORK - MAY NOT BE OPTIMAL!\n"); fflush(stdout);
        int dims[2] = {0,0};
        MPI_Dims_create(np, 2, dims); // however, you can also set nprow and npcol by hand, keeping constraing nprow*npcol=np
        if(dims[0]<dims[1])
        {
            md.p = dims[0]; // cartesian direction 0
            md.q = dims[1]; // cartesian direction 1
        }
        else
        {
            md.p = dims[1]; // cartesian direction 0
            md.q = dims[0]; // cartesian direction 1
        }
        
    }
    
    if(md.p==0 || md.q==0) 
    {
        if(iam==0) wprintf("ERROR: CANNOT SET p AND q VALUES! CHECK INPUT FILE SETTINGS!\n"); fflush(stdout);
        ABORTip(iam);
    }
    
    // check if input parameters are ok
    if(md.p*md.q!=np) 
    {
        if(iam==0) wprintf("ERROR: md.p*md.q!=np: CHECK INPUT FILE SETTINGS!\n"); fflush(stdout);
        ABORTip(iam);
    }
    
    // for hamiltonian diagonalization
    int iam_blacs, nprocs_blacs, ictxt;
    char * b_order ;
    int MONE = -1 , ZERO = 0 , ONE = 1;
    int DESCA[ 9 ];
    int p, q, nip, niq;
    int info;
    int Hsize = NXYZ*2; // size of hamiltonian matrix
    
    if(iam==0) wprintf("# HAMILTONIAN SIZE: %d x %d\n", Hsize, Hsize);
#ifdef MATRIX_IS_REAL
    if(iam==0) wprintf("# HAMILTONIAN TOTAL STORAGE: %.3fGB\n", 1.0*sizeof(double)*Hsize*Hsize/1024/1024/1024);
#else
    if(iam==0) wprintf("# HAMILTONIAN TOTAL STORAGE: %.3fGB\n", 1.0*sizeof(double complex)*Hsize*Hsize/1024/1024/1024);
#endif 
    
    /* initialize the BLACS grid for hamiltonian diagonalizaion- a virtual rectangular grid */
    if(iam==0) wprintf("# CREATING CBLACS GRID OF SIZE [%d x %d] WITH BLOCK SIZE [%d x %d]\n", md.p, md.q, md.mb, md.nb);
    b_order = "R" ;
    Cblacs_pinfo( &iam_blacs , &nprocs_blacs ) ;
    if ( nprocs_blacs < 1 ) 
        Cblacs_setup( &iam_blacs , &nprocs_blacs ) ;
    Cblacs_get( MONE , ZERO , &ictxt ) ; 
    Cblacs_gridinit( &ictxt , b_order , md.p , md.q ) ;  /* 'Row-Major' */
    Cblacs_gridinfo( ictxt , &p , &q , &ip , &iq ) ; /* ip,iq: the process row,column id */ 
    if(p!=md.p || q!=md.q) 
    {
        report_error(WSLDA_ERR_INOCRRECT_PQ, stderr);
        error_msg_mpi_abort(iam, p!=md.p || q!=md.q);
    }
    nip = numroc_( &Hsize, &md.mb, &ip, &ZERO, &p );
    niq = numroc_( &Hsize, &md.nb, &iq, &ZERO, &q );
#ifdef VERBOSE
    wprintf("# CHECK: iam=%d, ip=%d, iq=%d, nip=%d, niq=%d mb=%d nb=%d\n", iam, ip, iq, nip, niq, md.mb, md.nb);    
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
    cppmallocl(h_densities,DENSDIM,double);
    cppmallocl(h_densities_old,DENSDIM,double);
    cppmallocl(h_densities_partial,DENSDIM,double);
    cppmallocl(h_potentials,POTDIM,double);
    cppmallocl(h_potentials_old,POTDIM,double);
    cppmallocl(dc_params, MAX_USER_PARAMS, double);
    
    // solution buffers
    int SOLDIM;
    double *h_solution; // pointer to array with solution
    double *h_solution_old; // pointer to array with solution
    if(md.mixingtype=='d') { SOLDIM = DENSDIM; h_solution = h_densities ; h_solution_old = h_densities_old ; }
    else                   { SOLDIM = POTDIM ; h_solution = h_potentials; h_solution_old = h_potentials_old; }

    if(iam==0)
    {
        if(md.mixingtype=='d') wprintf("# MIXING TYPE: (d)ensities.\n");
        else                   wprintf("# MIXING TYPE: (p)otentials.\n");      
    }
    
    // For easier access to data
    wslda_density densall = convert_into_wslda_density(h_densities, NXYZ);
    wslda_density densall_partial = convert_into_wslda_density(h_densities_partial, NXYZ);
    wslda_potential potsall = convert_into_wslda_potential(h_potentials, NXYZ, mu);
            
    // For Broyden method
    cppmallocl(dens_in, (md.Mbroyden + 1), double*);
    cppmallocl(dens_out, (md.Mbroyden + 1), double*);
    for (i = 0; i < (md.Mbroyden + 1); i++){
        cppmallocl(dens_in[i], SOLDIM + 2, double);
        cppmallocl(dens_out[i], SOLDIM + 2, double);
        
        // reset values
        for(j=0; j<SOLDIM + 2; j++) dens_in[i][j]=0.0;
        for(j=0; j<SOLDIM + 2; j++) dens_out[i][j]=0.0;
    }

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
    
    if(iam==0) cpu_exec( create_directory(md.outprefix) );
    
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
    for(ixyz=0; ixyz< ENERGYITEMS; ixyz++) energy[ixyz] = 0.0;
    npart[SPINA]=NUMERICAL_ZERO; npart[SPINB]=NUMERICAL_ZERO;
    
    // ===================================================================================
    // =============================== INITIAL STATE =====================================
    // ===================================================================================
    if(md.inittype==0 || md.inittype==10) // Start from uniform solution
    {
        if(md.inittype==0)
        {
            if(iam==0) wprintf("# CREATING UNIFORM SOLUTION...\n");

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
            if(iam==0) wprintf("# READING UNIFORM SOLUTION...\n");
            if(iam==0) { cpu_exec( read_uniform(&nwf, iam==0) ); }
            MPI_Bcast( &__md_pca_uniform , sizeof(metadata_pca_uniform_t), MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
            MPI_Bcast( &nwf , 1, MPI_INT , 0 , MPI_COMM_WORLD ) ;
        }
                       
        // initialize potentials and densities
        for(ixyz=0; ixyz<NXYZ; ixyz++)
        {
            // potentials
            potsall.V_a[ixyz]=__md_pca_uniform.V_a; 
            potsall.V_b[ixyz]=__md_pca_uniform.V_b; 
            potsall.delta[ixyz]=__md_pca_uniform.delta + I*0.0;
            potsall.alpha_a[ixyz]=__md_pca_uniform.alph_a; 
            potsall.alpha_b[ixyz]=__md_pca_uniform.alph_b;
            potsall.A_a_x[ixyz]=0.0; 
            potsall.A_a_y[ixyz]=0.0;  
            potsall.A_a_z[ixyz]=0.0;  
            potsall.A_b_x[ixyz]=0.0; 
            potsall.A_b_y[ixyz]=0.0; 
            potsall.A_b_z[ixyz]=0.0; 
            
            // densities
            densall.rho_a[ixyz]=__md_pca_uniform.n0_a; 
            densall.rho_b[ixyz]=__md_pca_uniform.n0_b; 
            densall.tau_a[ixyz]=__md_pca_uniform.tau_a; 
            densall.tau_b[ixyz]=__md_pca_uniform.tau_b; 
            densall.nu[ixyz]=__md_pca_uniform.nu; 
            densall.j_a_x[ixyz]=0.0;
            densall.j_a_y[ixyz]=0.0;
            densall.j_a_z[ixyz]=0.0;
            densall.j_b_x[ixyz]=0.0;
            densall.j_b_y[ixyz]=0.0;
            densall.j_b_z[ixyz]=0.0;           
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
        
        // set energy buffers
        energy[EKIN]=__md_pca_uniform.ekin;
        energy[EPOT]=__md_pca_uniform.epot;
        energy[EPAIR]=__md_pca_uniform.epair;
    }
    else if(md.inittype==5) // start from checkpoint
    {
        if(iam==0)
        {
            // Check format of checkpoint
            i = wslda_stcheckpoint_format(3);
            if(i==WSLDA_ST_CHECKPOINT_DAT)
            {
                // this format supports extansions and interpolatons
                int _interop, _resop;
                file_operation( wslda_st_required_operations(3, &_interop, &_resop) );
                
                // prepare for reading
                double trd_consts[11];
                j=0; // file idx

                // convert checkpoint
                if(_interop>0)
                {
                    file_operation(
                        wslda_st_checkpoint_convert(ST_CHECKPOINT_RESIZE, j, 3, &it, 11, trd_consts, POTDIM, h_potentials, DENSDIM, h_densities, ENERGYITEMS, energy, SOLDIM + 2, dens_in, dens_out)
                    ); 
                    j++;
                }
                
                if(_resop==23)
                {
                    file_operation(
                        wslda_st_checkpoint_convert(ST_CHECKPOINT_2D_TO_3D, j, 3, &it, 11, trd_consts, POTDIM, h_potentials, DENSDIM, h_densities, ENERGYITEMS, energy, SOLDIM + 2, dens_in, dens_out)
                    ); 
                    j++;
                }
                if(_resop==13)
                {
                    file_operation(
                        wslda_st_checkpoint_convert(ST_CHECKPOINT_1D_TO_3D, j, 3, &it, 11, trd_consts, POTDIM, h_potentials, DENSDIM, h_densities, ENERGYITEMS, energy, SOLDIM + 2, dens_in, dens_out)
                    ); 
                    j++;
                }
                
                // read checkpoint
                file_operation(
                    wslda_st_read_checkpoint(j, 3, &it, 11, trd_consts, POTDIM, h_potentials, DENSDIM, h_densities, ENERGYITEMS, energy, SOLDIM + 2, dens_in, dens_out)
                );
                
                // decode constants
                dc_mu_a=trd_consts[0];  dc_mu_b=trd_consts[1];  dc_mu_a_old=trd_consts[2];  dc_mu_b_old=trd_consts[3];  dc_ec=trd_consts[4];  beta=trd_consts[5];  eF=trd_consts[6];  kF=trd_consts[7];  Effg=trd_consts[8];  npart[SPINA]=trd_consts[9];  npart[SPINB]=trd_consts[10]; 
                
                // copy init checkpoint
                copy_initcheckpoint();
            }
            else if(i==WSLDA_ST_CHECKPOINT_OLD)
            {
                // NOTE: It wll be removed in future
                // here I keep it only to be compatible with our past caculations
                sprintf(file_name, "%s/checkpoint.s3dpca", md.inprefix);
                wprintf("# READING CHECKPOINT FILE `%s`\n", file_name);
                wprintf("# !!! !!! YOU ARE USING OLD CHECKPOINT FORMAT !!! !!! SUPPORT OF THIS FORMAT WILL BE REMOVED IN FUTURE!\n");
                FILE * pFile = fopen(file_name, "rb");
                if(pFile==NULL)
                {
                    wprintf("# CANNOT FIND CHECKPOINT FILE: `%s`\n", file_name); fflush(stdout);
                    ABORT_NOBARRIER;
                }
                
                // write all nescesary data to file
                fread(&it          , sizeof(int)         , 1 , pFile); // iteration number
                fread(&dc_mu_a     , sizeof(double)      , 1 , pFile); 
                fread(&dc_mu_b     , sizeof(double)      , 1 , pFile); 
                fread(&dc_ec       , sizeof(double)      , 1 , pFile); 
                fread(&beta        , sizeof(double)      , 1 , pFile); 
                fread(&eF          , sizeof(double)      , 1 , pFile); 
                fread(&kF          , sizeof(double)      , 1 , pFile);
                fread(&Effg        , sizeof(double)      , 1 , pFile);
                fread(h_potentials , sizeof(double)*POTDIM , 1, pFile);
                fread(h_densities  , sizeof(double)*DENSDIM, 1, pFile);
                fread(energy       , sizeof(double)      , ENERGYITEMS , pFile);
                fread(npart        , sizeof(double)      , 2 , pFile);
                fread(&dc_mu_a_old , sizeof(double)      , 1 , pFile); 
                fread(&dc_mu_b_old , sizeof(double)      , 1 , pFile);
                for (i = 0; i < (md.Mbroyden + 1); i++) fread(dens_in[i]   , sizeof(double) , SOLDIM + 2 , pFile);
                for (i = 0; i < (md.Mbroyden + 1); i++) fread(dens_out[i]  , sizeof(double) , SOLDIM + 2 , pFile);
                    
                fclose(pFile);

            }
            else
            {
                sprintf(file_name, "%s_checkpoint.dat", md.inprefix);
                wprintf("# CANNOT FIND CHECKPOINT FILE: `%s`\n", file_name); fflush(stdout);
                ABORT_NOBARRIER;
            }
                        
            wprintf("# CHECKPOINT READ: it=%d\n", it);
            wprintf("# CHECKPOINT READ: dc_mu_a=%16.8g  dc_mu_b=%16.8g\n", dc_mu_a, dc_mu_b);
            wprintf("# CHECKPOINT READ: dc_ec=%16.8g  beta=%16.8g\n", dc_ec, beta);
            wprintf("# CHECKPOINT READ: eF=%16.8g  kF=%16.8g  Effg=%16.8g\n", eF, kF, Effg);
            wprintf("# CHECKPOINT READ: ------- NPART -------\n");
            wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
                "SPINA", npart[SPINA], npart[SPINA], (npart[SPINA]-npart[SPINA]));
            wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
                "SPINB", npart[SPINB], npart[SPINB], (npart[SPINB]-npart[SPINB]));  
            wprintf("# CHECKPOINT READ: ------- ENERGY -------\n");
            E_tot=0.0; E_tot_old=0.0;
            for(i=0; i<ENERGYITEMS; i++)
            {
                E_tot+=energy[i]; 
                
                wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
                    energy_labels[i], energy[i]/Effg, energy[i]/Effg, (energy[i]-energy[i])/Effg);
            }
            wprintf("  ------------------------------------------------------------------------\n");
            wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
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
        MPI_Bcast(h_potentials , POTDIM , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(h_densities  , DENSDIM, MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(energy       , ENERGYITEMS , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(npart        , 2 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(&dc_mu_a_old , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&dc_mu_b_old , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        for (i = 0; i < (md.Mbroyden + 1); i++) MPI_Bcast(dens_in[i] , SOLDIM + 2 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        for (i = 0; i < (md.Mbroyden + 1); i++) MPI_Bcast(dens_out[i], SOLDIM + 2 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    }
    else if(md.inittype==-1) // start from checkpoint
    {
        if(iam==0) wprintf("# BUFFERS WILL BE INITIALIZED VIA modify_densities AND modify_potentials FUNCTIONS\n");
        if(iam==0) wprintf("# FORCING: resetit=1 AND nomixstart=1\n");
        if(iam==0) wprintf("# FORCING: ITERATION NUMBER: it=-1\n");
        md.resetit=1;
        md.nomixstart=1;
    }
    else
    {
        if(iam==0) wprintf("NOT SUPPORTED INITTYPE=%d!\n", md.inittype);
        ABORT;
    }
    
    if(md.inittype>=0) cpu_exec( wslda_check_array_against_naninf(ENERGYITEMS, energy) );
    
    // ===================================================================================
    // ================================== EXTRA DATA =====================================
    // ===================================================================================
    if(iam==0) extra_data_size = get_extra_data_size(md.params);
    MPI_Bcast( &extra_data_size , sizeof(size_t) , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
    if(extra_data_size>0)
    {
        if(iam==0) wprintf("# EXTRA_DATA IS ACTIVE.\n");
        if(iam==0) wprintf("# ALLOCATING EXTRA_DATA OF SIZE %ld B.\n", extra_data_size); fflush(stdout);
        if ( ( extra_data = (void *) malloc( extra_data_size ) ) == NULL  )
        {                                                             
            wfprintf( stderr , "error: cannot malloc()! Exiting!\n") ; 
            wfprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; 
            something_to_cheer_you_up_pid0(stdout);
            MPI_Finalize() ;
            /* Arrays will be cleared automatically */
            return( EXIT_FAILURE ) ; 
        }
        if(iam==0) wprintf("# EXECUTING: load_extra_data(%zu, extra_data, input->params)\n", extra_data_size);
        if(iam==0) cpu_exec( load_extra_data(extra_data_size, extra_data, md.params) );
        MPI_Bcast( extra_data , extra_data_size , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;

        // reproducibility pack
        if(iam==0) save_extradata_to_file(extra_data_size, extra_data);
    }
    dc_extra_data_size=extra_data_size;
    dc_extra_data=extra_data;

    // ===================================================================================
    // ================================== FFTW PLANS =====================================
    // ===================================================================================
    if(iam==0) { wprintf("# CREATING FFTW PLANS...\n"); fflush(stdout); }
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
    if(md.resetit) it=0; // reset iterator counter
    
    // NOTE settings some variables
#ifndef UNIFORM_TEST_MODE
    if(md.referencekF>0.0) kF = md.referencekF;
    eF = 0.5*kF*kF;
    Effg = 0.6 * (md.Na+md.Nb) * eF;
    beta = 1.0 / (md.temperature * eF);    
    if(md.ec>0.0) dc_ec = md.ec; 
    else          dc_ec = M_PI*M_PI/(2.*DX*DX);
#endif
    
    mu[SPINA]=dc_mu_a; mu[SPINB]=dc_mu_b;
    
    wdata_metadata wdmd; 
    file_operation( create_wdata_metadata(&md, 3, 1.0*(it-1), 1.0, md.spinsymmetry, &wdmd) );
    
    // set constants
    wdata_setconst(&wdmd, "kF", kF);
    wdata_setconst(&wdmd, "eF", eF);
    wdata_setconst(&wdmd, "mu_a", mu[SPINA]);
    wdata_setconst(&wdmd, "mu_b", mu[SPINB]);
    
    // prepare database
    MPI_Barrier(MPI_COMM_WORLD);
    if(iam==0) 
    {
        file_operation( clear_files(&md, &wdmd) );
        file_operation( write_wdata_metadata_file(&md, &wdmd, "st-wslda-3d") );
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(iam==0) wprintf("# EXECUTING: process_params(input->params, [%f], [%f,%f], %zu, extra_data)\n", kF, mu[SPINA], mu[SPINB], extra_data_size);
    for(i=0; i<MAX_USER_PARAMS; i++) dc_params[i]=md.params[i];
    mu[SPINA]=dc_mu_a; mu[SPINB]=dc_mu_b;
    process_params(dc_params, &kF, mu, extra_data_size, extra_data);
    dc_mu_a=mu[SPINA]; dc_mu_b=mu[SPINB];
    modify_densities(it-1, densall, dc_params, extra_data_size, extra_data) ;
    modify_potentials(it-1, densall, potsall, dc_params, extra_data_size, extra_data) ;
    dc_mu_a=mu[SPINA]; dc_mu_b=mu[SPINB];
    file_operation( write_measurments(&wdmd, MPI_COMM_WORLD, "st", it-1, densall, potsall) );
    if(iam==0) file_operation( write_wdata_metadata_file(&md, &wdmd, "st-wslda-3d") );
    ECHOLINE;
    
    double dc_ec_l=-1.0*dc_ec; // lower bound for states extraction
    double dc_ec_u= 1.0*dc_ec; // upper bound for states exteraction
    if(md.spinsymmetry>0) dc_ec_l=0.0; // take only positive states
    
    // special case for saving
    if(md.writewf==1 && md.maxiters==1) saving_iteration=1;
    
    // Create run log and add entry
    if(iam==0) cpu_exec( logger_create_header(execcmd) );
    
    // ====================================================================================
    // ====================================== WORKSPACE ===================================
    // ====================================================================================    
#if defined(USE_SCALAPACK_PZHEEVD) || defined(USE_SCALAPACK_PZHEEV) || defined(USE_SCALAPACK_PZHEEVR)
    /* Query and allocate the optimal workspace */
#ifdef USE_SCALAPACK_PZHEEVR
    if(iam==0) wprintf("# PREPARING WORKING BUFFERS FOR: `pzheevr`\n");
#endif
#ifdef USE_SCALAPACK_PZHEEVD
    if(iam==0) wprintf("# PREPARING WORKING BUFFERS FOR: `pzheevd`\n");
#endif
#ifdef USE_SCALAPACK_PZHEEV
    if(iam==0) wprintf("# PREPARING WORKING BUFFERS FOR: `pzheev`\n");
#endif    
#ifdef MATRIX_IS_REAL
    if(iam==0) wprintf("# MATRIX IS ASSUMED TO BE REAL! SWITCHING TO pdsyev* VERSION!\n");
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
    wprintf("# pzheevd-set[%d]: info=%d  lwork=%d  lrwork=%d  liwork=%d\n",iam, info, lwork, lrwork, liwork);
    size_t size_workspace = sizeof(double complex)*lwork+sizeof(double)*lrwork;
    wprintf("# SIZE OF WORKSPACE [%d]: size_workspace=%.3fGB\n",iam, 1.0*(double)(size_workspace)/1024./1024./1024.);  
#endif
 
#endif
    
#ifdef USE_ELPA
    if(iam==0) wprintf("# SETTING UP ELPA...\n");
    if (elpa_init(20181112) != ELPA_OK) error_msg_mpi_abort(iam, ELPA API version not supported);
    
    elpa_t handle;
    handle = elpa_allocate(&info); 
    if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    
    /* Set parameters */
    elpa_set(handle, "na", Hsize, &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    int elpa_nev = (int)(ELPA_NEV_FRACTION*Hsize); if(elpa_nev>Hsize) elpa_nev=Hsize;
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
    
//     ISSUE: autotuning
//     elpa_autotune_t autotune_handle = elpa_autotune_setup(handle, ELPA_AUTOTUNE_FAST, ELPA_AUTOTUNE_DOMAIN_COMPLEX, &info); 
//     if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
//     int elpa_atotue_unfinished=1;
    
    elpa_set(handle, "solver", ELPA_USE_SOLVER, &info);  
    if(iam==0) wprintf("# ELPA: SETTINGS SOLVER: `%s`\n", STRINGIZE(ELPA_USE_SOLVER));
    if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
    
#ifdef ELPA_USE_GPU
    if(iam==0) wprintf("# ELPA: ACTIVATING GPUs\n");
    elpa_set(handle, "gpu", 1, &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
#else
    elpa_set(handle, "gpu", 0, &info); if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
#endif
#ifdef MATRIX_IS_REAL
    elpa_set(handle, "real_kernel", ELPA_USE_REAL_KERNEL, &info); 
    if(iam==0) wprintf("# ELPA: SETTINGS REAL KERNEL: `%s`\n", STRINGIZE(ELPA_USE_REAL_KERNEL));
#else
    elpa_set(handle, "complex_kernel", ELPA_USE_COMPLEX_KERNEL, &info); 
    if(iam==0) wprintf("# ELPA: SETTINGS COMPLEX KERNEL: `%s`\n", STRINGIZE(ELPA_USE_COMPLEX_KERNEL));
#endif
    if(info!=ELPA_OK) error_msg_mpi_abort(iam, info!=ELPA_OK);
            
    if(iam==0) wprintf("# SETTING UP OF ELPA  DONE.\n");
#endif
        
    fflush(stdout);
         
    // ===================================================================================
    // ========================== SELF-CONSITENT LOOP ====================================
    // ===================================================================================     
    
    // special case - only one iteration for diagonalization
    if(md.maxiters==1 && md.writewf==1) saving_iteration=1;
    
    while(1) // do until reached self-consitency
    {
        b_t();
        if((kziter+1)==md.maxiters && md.writewf==1) 
        {
            if(iam==0) wprintf("# EXECUTING LAST ITERATION WITH SAVING DATA [md.writewf==1]\n");
            saving_iteration=1;
        }
        rt_zheev=0.0; rt_dens=0.0; rt_pot=0.0; rt_other=0.0; rt_me=0.0; rt_redistrib=0.0;
        b_t();
        // Make copy of densities and potentials
        for(ixyz=0; ixyz<DENSDIM; ixyz++) h_densities_old[ixyz]  = h_densities[ixyz];
        for(ixyz=0; ixyz<POTDIM ; ixyz++) h_potentials_old[ixyz] = h_potentials[ixyz];
        for(ixyz=0; ixyz<DENSDIM; ixyz++) h_densities_partial[ixyz] = 0.0; // reset
        for(i=0; i< ENERGYITEMS; i++) energy_old[i] = energy[i]; // make copy
        npart_old[SPINA]=npart[SPINA]; npart_old[SPINB]=npart[SPINB]; // make copy 
        Lz_a_old=Lz_a; Lz_b_old=Lz_b; Lz_old=Lz; 
        S_old=S; S=0.0; // save value of entropy and reset buffer
        dc_mu_a_old = dc_mu_a; dc_mu_b_old = dc_mu_b;
        
        if(md.referencekF>0.0) kF = md.referencekF;
        else 
        {
            double _max_dens=0.0;
            for(ixyz=0; ixyz<BLOCKLENGTH; ixyz++) if(densall.rho_a[ixyz]+densall.rho_b[ixyz]>_max_dens) _max_dens=densall.rho_a[ixyz]+densall.rho_b[ixyz];
            kF = pow(3.*M_PI*M_PI*_max_dens,1./3.);
        }
#ifndef UNIFORM_TEST_MODE        
        eF = 0.5 * kF * kF;
        beta = 1.0 / (md.temperature * eF);
#endif
        if(iam==0) wprintf("# EXECUTING: process_params(input->params, [%f], [%f,%f], %zu, extra_data)\n", kF, mu[SPINA], mu[SPINB], extra_data_size);
        for(i=0; i<MAX_USER_PARAMS; i++) dc_params[i]=md.params[i];
        mu[SPINA]=dc_mu_a; mu[SPINB]=dc_mu_b;
        process_params(dc_params, &kF, mu, extra_data_size, extra_data);
        dc_mu_a=mu[SPINA]; dc_mu_b=mu[SPINB];
        ECHOLINE;
        rt_other+=e_t(0);
        
        // matrix elements
        b_t();
        cpu_exec( compute_matrix_elements_3d(&bgrid, it, densall, potsall, &mdfft, h, me_d_dx, me_d_dy, me_d_dz) );
        rt_me+=e_t(0);
        ECHOLINE;
        
        // diagonalize
        b_t();
        if(iam==0) wprintf("# DIAGONALIZATION %d...\n", it); fflush(stdout);
#ifdef MATRIX_IS_REAL
        // convert matrix to real version
//         for(ixyz=0; ixyz<nip*niq; ixyz++) if(fabs(cimag(h[ixyz]))>1.0e-12) {wprintf("# ERROR: matrix has imginary components!\n"); ABORT_NOBARRIER;}
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
            wprintf( "The algorithm failed to compute eigenvalues!\n" );
            ABORT_NOBARRIER;
        }
#endif
#ifdef USE_ELPA
//  ISSUE: autotue
//         if (elpa_atotue_unfinished != 0 ) elpa_atotue_unfinished = elpa_autotune_step(handle, autotune_handle, &info);
//         sprintf(file_name, "%s_elpa.autotune", md.outprefix);
//         if(iam==0) elpa_autotune_save_state(handle, autotune_handle, file_name, &info);
//         if(iam==0) wprintf("----> elpa_atotue_unfinished=%d info=%d\n", elpa_atotue_unfinished, info);
//         if (elpa_atotue_unfinished == 0) 
//         {
//             elpa_autotune_set_best(handle, autotune_handle, &info);  // from now on use values used by autotuning
//             elpa_autotune_deallocate(autotune_handle, &info);        // cleanup autotuning
//         }
//         if (elpa_atotue_unfinished == 0 && iam==0) 
//         {
//             wprintf("# ELPA autotuning finished in the %d th scf step \n",kziter);
//         }

#ifdef MATRIX_IS_REAL
        elpa_eigenvectors(handle, hR, En, UR, &info);
#else
        elpa_eigenvectors(handle, h , En, U , &info);
#endif
        if( info !=ELPA_OK ) 
        {
            wprintf( "The algorithm failed to compute eigenvalues!\n" );
            ABORT_NOBARRIER;
        }     
        if(elpa_nev<Hsize && En[elpa_nev-1]<dc_ec)
        {
            wprintf( "# !!!!!!! WARNING !!!!!!!: ELPA_NEV_FRACTION IS TOO SMALL!!!!!!!\n" );
        }
#endif
#ifdef MATRIX_IS_REAL
        // convert result back to complex
        for(ixyz=0; ixyz<nip*niq; ixyz++) hR[ixyz] = UR[ixyz]; // back to hR matrix
        for(ixyz=nip*niq-1; ixyz>=0; ixyz--) U[ixyz] = hR[ixyz] + I*0.0;
#endif
        rt_zheev+=e_t(0);
        if(iam==0) wprintf("# DIAGONALIZATION DONE [%.0f sec]\n", rt_zheev); fflush(stdout);
        ECHOLINE;
        
        // preparation for density computation
        b_t();
#ifdef USE_SCALAPACK_PZHEEVR
        nwf=pzheevr_m;
#else
        // do matrix redistribution for computation of densities
        for(i=0; i<Hsize-1; i++) if(En[i]>En[i+1]) 
        {
            wprintf( "Eigenvalues are not sorted correctly!\n" );
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
        if(iam==0) wprintf("# NUMBER OF EXTRACTED nwf IN ENERGY RANGE [-ecut,+ecut] IS %d (%.1f%% OF TOTAL NUMBER OF STATES)\n", nwf, 100.0*nwf/Hsize);
        
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
        wprintf("# CHECK-D: iam=%d, ip_d=%d, iq_d=%d, nip_d=%d, niq_d=%d\n", iam, ip_d, iq_d, nip_d, niq_d);    
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
            if(iam==0) wprintf("# WRITING WAVE-FUNCTIONS WHERE |E_n/eF|<%.6g\n", md.writeecut); fflush(stdout);
            if(iam==0) wprintf("# WRITING USING %d I/O GROUPS\n", md.iogroups); fflush(stdout);
            
            // Creating IO groups
            int idgroup; // identifier of the group
            MPI_Comm mpi_comm_group;
            int gr_iam, gr_np;
            
            idgroup=iam % md.iogroups;
            MPI_Comm_split(MPI_COMM_WORLD, idgroup, iam, &mpi_comm_group);
            MPI_Comm_rank(mpi_comm_group, &gr_iam);
            MPI_Comm_size(mpi_comm_group, &gr_np);
            if(gr_iam==0) wprintf("# I/O GROUP %d WITH %d PROCESSES HAS BEEN SUCCESSFULLY CREATED.\n", idgroup, gr_np); fflush(stdout);
            
            // !!!!! I/O can be written better with MPI I/O - maybe it will be improved in future !!!!!
            MPI_Status MPIStat;
            int lastwf=0;
        
            // ----- part 0 ------
            if(gr_iam==0)
            {
                // Create empty files
                sprintf(file_name, "%s/s3dpca.%04d.wfu", md.outprefix, idgroup);
                file_operation( touch_file(file_name) );
            
                sprintf(file_name, "%s/s3dpca.%04d.wfv", md.outprefix, idgroup);
                file_operation( touch_file(file_name) );

                sprintf(file_name, "%s/s3dpca.%04d.en", md.outprefix, idgroup);
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
            sprintf(file_name, "%s/s3dpca.%04d.info", md.outprefix, idgroup);
            mu[SPINA] = dc_mu_a; mu[SPINB] = dc_mu_b;
            if(gr_iam==0) file_operation( create_checkpoint_info_pca(file_name, lastwf, NX, NY, NZ, DX, DY, DZ, kF, mu, md.writeecut*eF, beta) );
                
            double rt = e_t(0);
            if(gr_iam==0) wprintf("# DATA WRITING BY I/O GROUP %d TOOK %.1f SEC. WRITTEN %.2fMB. WRITTEN STATES=%d\n", idgroup, rt, (double)1.*lastwf*NX*NY*NZ*2*16/1024./1024., lastwf); fflush(stdout);
            
            // recompute total nwf
            if(gr_iam!=0) lastwf=0;
            MPI_Allreduce( &lastwf, &nwf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if(iam==0) wprintf("# DATA SAVING DONE. IN TOTAL WRITTEN %d WAVE-FUNCTIONS.\n", nwf);
            rt_other+=rt;
        }
        
        // For wf-reproducibility pack - save present checkpoint file
        if(saving_iteration==1 && iam==0) copy_checkpoint();
        
        // --------- DENSITIES ----------
        b_t();
        // compute contribution to the densities
        cpu_exec( compute_contribution_to_densities(En_d_local, U_d, niq_d, beta, densall_partial, &mdfft, md.spinsymmetry, &S) );
        // compute densities as global reduction
        MPI_Allreduce( h_densities_partial, h_densities, DENSDIM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce( MPI_IN_PLACE, &S, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        if(md.spinsymmetry==1)  for(ixyz=0; ixyz<NXYZ; ixyz++) // impose by hand symmetry on densities
        {
            densall.rho_a[ixyz]=densall.rho_b[ixyz];
            densall.tau_a[ixyz]=densall.tau_b[ixyz];
            densall.j_a_x[ixyz]=densall.j_b_x[ixyz];
            densall.j_a_y[ixyz]=densall.j_b_y[ixyz];
            densall.j_a_z[ixyz]=densall.j_b_z[ixyz];
        }
        if(md.nocurrents) for(ixyz=0; ixyz<NXYZ; ixyz++) // impose by hand no currents
        {
            densall.j_a_x[ixyz]=0.0; densall.j_b_x[ixyz]=0.0;
            densall.j_a_y[ixyz]=0.0; densall.j_b_y[ixyz]=0.0;
            densall.j_a_z[ixyz]=0.0; densall.j_b_z[ixyz]=0.0;
        }
        
#ifndef TAU_COMPUTATION_VIA_GRADIENTS  
        // finalize computation of tau
        cpu_exec( density_caculate_tau(densall, &mdfft) );
#endif
        
        modify_densities(it, densall, dc_params, extra_data_size, extra_data) ;
        rt_dens+=e_t(0);
                        
        // free temporary resources
        b_t();
        free(U_d); 
        free(En_d);
        free(En_d_local);
        Cblacs_gridexit( ictxt_d );
        rt_other+=e_t(0);
        
        // ------------------ compute new potentials ------------------
        b_t();
        cpu_exec( compute_potentials(it, densall, potsall) ); 
        mu[SPINA]=dc_mu_a; mu[SPINB]=dc_mu_b;
        modify_potentials(it, densall, potsall, dc_params, extra_data_size, extra_data) ;
        dc_mu_a=mu[SPINA]; dc_mu_b=mu[SPINB];
        rt_pot+=e_t(0);
        
        // ------------------ update chemical potentials ------------------
        b_t();
        if(iam==0) wprintf("# MUCHANGE FROM: mu_a/eF=%16.8g  mu_b/eF=%16.8g\n", dc_mu_a/eF, dc_mu_b/eF);
        i=0; // as flag for broyden
        if(it>0 && saving_iteration==0) // skip updating the potential if it is saving iteration
        {       
            npart[SPINA]=0.0; npart[SPINB]=0.0;
            for(ixyz=0; ixyz<NXYZ; ixyz++) {npart[SPINA]+=densall.rho_a[ixyz]; npart[SPINB]+=densall.rho_b[ixyz];}
            npart[SPINA]*=DXYZ; npart[SPINB]*=DXYZ; 
            double muchange_a = md.muchange*(npart[SPINA] - md.Na)/md.Na;
            double muchange_b = md.muchange*(npart[SPINB] - md.Nb)/md.Nb;
            if(fabs(muchange_a)>md.mumaxchange*eF)
            {
                if(muchange_a>0.0) muchange_a=     md.mumaxchange*eF;
                else               muchange_a=-1.0*md.mumaxchange*eF;
                i=1; // deactivate broyden
            }
            if(fabs(muchange_b)>md.mumaxchange*eF)
            {
                if(muchange_b>0.0) muchange_b=     md.mumaxchange*eF;
                else               muchange_b=-1.0*md.mumaxchange*eF;
                i=1; // deactivate broyden
            }
            dc_mu_a -= muchange_a;
            dc_mu_b -= muchange_b;  
            if(md.spinsymmetry==1) dc_mu_b=dc_mu_a; // activate constraint
        }
        if(iam==0) wprintf("# MUCHANGE TO  : mu_a/eF=%16.8g  mu_b/eF=%16.8g\n", dc_mu_a/eF, dc_mu_b/eF);
        
        if(i==1 && it>md.startbroyden && it<md.stopbroyden && md.broyden == 1 && md.broydenautores==1)
        {
            md.startbroyden=it;
            if(iam==0) wprintf("# CHEMICAL POTENTIAL HAS CHANGED BY `mumaxchange`! RESTARTING BROYDEN TO AVOID INSTABILITY.\n");
        }
        rt_other+=e_t(0);
        
        // ------------------ mix densities ------------------
        b_t();
        if(md.nomixstart==1 && kziter==0) //special case - no mixing for the first iteration
        {
            // pass - do not mix
            if(iam==0) wprintf("# DENSITIES MIX: SPECIAL CASE: NO MIXING FOR STARTING ITERATION (nomixstart==1)! MIXING SKIPPED!\n");
        }
        else if(saving_iteration==1) //special case - saving interation
        {
            // typically results from this iteration are loded into dynamical code
            // for clear comparision of read corretness skip mixing here
            
            // pass - do not mix
            if(iam==0) wprintf("# DENSITIES MIX: SPECIAL CASE: SAVING ITERATION! MIXING SKIPPED!\n");
        }
        else if (((it-md.startbroyden) >= 0) && ((it-md.startbroyden) < (md.Mbroyden + 1)) && (md.broyden == 1) )
        {
            int rkziter=it-md.startbroyden;
            for(ixyz = 0; ixyz < SOLDIM; ixyz++) {
				dens_in[rkziter][ixyz] = h_solution_old[ixyz];
				dens_out[rkziter][ixyz] = h_solution[ixyz];
            	h_solution[ixyz] = md.linearmixing * h_solution[ixyz] + (1.0 - md.linearmixing) * h_solution_old[ixyz];
            }
			dens_in[rkziter][ixyz+0] = dc_mu_a_old;
			dens_out[rkziter][ixyz+0] = dc_mu_a;
			dens_in[rkziter][ixyz+1] = dc_mu_b_old;
			dens_out[rkziter][ixyz+1] = dc_mu_b;
            if(iam==0) wprintf("# DENSITIES MIX: BROYDEN IS STORING DATA, MIXING=LINEAR\n");
        }
        else if (((it-md.startbroyden) >= (md.Mbroyden+1)) && (it-md.stopbroyden)<=0 && (md.broyden == 1))
        {
        	update_mu(dens_in, dens_out, h_solution_old, h_solution, md.Mbroyden, SOLDIM, dc_mu_a, dc_mu_b, dc_mu_a_old, dc_mu_b_old);
        	Broyden_mu(h_solution, dens_in, dens_out, md.Mbroyden, SOLDIM+2, omega_0, omega_n, omega_k, md.broydenmixing, &dc_mu_a, &dc_mu_b);
            
            if     (dc_mu_a-dc_mu_a_old>md.mumaxchange*eF) dc_mu_a = dc_mu_a_old+md.mumaxchange*eF;
            else if(dc_mu_a_old-dc_mu_a>md.mumaxchange*eF) dc_mu_a = dc_mu_a_old-md.mumaxchange*eF;
            
            if     (dc_mu_b-dc_mu_b_old>md.mumaxchange*eF) dc_mu_b = dc_mu_b_old+md.mumaxchange*eF;
            else if(dc_mu_b_old-dc_mu_b>md.mumaxchange*eF) dc_mu_b = dc_mu_b_old-md.mumaxchange*eF;
            
            if(md.spinsymmetry==1) dc_mu_b=dc_mu_a; // activate constraint
            
            if(iam==0) wprintf("# MUCHANGE BROY: mu_a/eF=%16.8g  mu_b/eF=%16.8g\n", dc_mu_a/eF, dc_mu_b/eF);
            if(iam==0) wprintf("# DENSITIES MIX: BROYDEN MIXING\n");
        }
        else
        {
        	for(ixyz = 0; ixyz < SOLDIM; ixyz++) h_solution[ixyz] = md.linearmixing * h_solution[ixyz] + (1.0-md.linearmixing) * h_solution_old[ixyz];
            if(iam==0) wprintf("# DENSITIES MIX: LINEAR MIXING\n");
        }
        
        // impose by hand nonegativity of densities
        double dens_min = 1.0e-14;
        for (ixyz = 0; ixyz < NXYZ; ixyz++)
        {
        	if (densall.rho_a[ixyz] < 0.) densall.rho_a[ixyz] = dens_min;
        	if (densall.rho_b[ixyz] < 0.) densall.rho_b[ixyz] = dens_min;
        	if (densall.tau_a[ixyz] < 0.) densall.tau_a[ixyz] = dens_min;
        	if (densall.tau_b[ixyz] < 0.) densall.tau_b[ixyz] = dens_min;
        }
        rt_other+=e_t(0);
        
        // ------------------------ energy ----------------------
        cpu_exec( compute_energy(it, densall, potsall, energy, npart) );
        
        // ---------------------- entropy -----------------------
        if(iam==0) wprintf("# ENTROPY [T/eF=%16.8g]: it=%d\n", 1.0/(beta*eF), it);
        if(iam==0) wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "S/NkB", S/(npart[SPINA]+npart[SPINB]), S_old/(npart_old[SPINA]+npart_old[SPINB]), 
                           (S/(npart[SPINA]+npart[SPINB])-S_old/(npart_old[SPINA]+npart_old[SPINB])));
        
        // ------------------ angular momentum ------------------
        b_t();
        cpu_exec( compute_angular_momentum_Lz(densall.j_a_x, densall.j_a_y, &Lz_a) );
        cpu_exec( compute_angular_momentum_Lz(densall.j_b_x, densall.j_b_y, &Lz_b) );
        Lz = Lz_a + Lz_b; // total angular momentum        
        if(iam==0) wprintf("# ANGULAR MOMENTUM: it=%d\n", it);
        if(iam==0) wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "Lza/Na", Lz_a/npart[SPINA], Lz_a_old/npart_old[SPINA], (Lz_a/npart[SPINA]-Lz_a_old/npart_old[SPINA]));
        if(iam==0) wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "Lzb/Nb", Lz_b/npart[SPINB], Lz_b_old/npart_old[SPINB], (Lz_b/npart[SPINB]-Lz_b_old/npart_old[SPINB]));  
        if(iam==0) wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "Lzt/Nt", Lz/(npart[SPINA]+npart[SPINB]), Lz_old/(npart_old[SPINA]+npart_old[SPINB]), (Lz/(npart[SPINA]+npart[SPINB])-Lz_old/(npart_old[SPINA]+npart_old[SPINB])));
        
        // ------------------ check convergence ------------------
        is_converged=1;
        if(iam==0) wprintf("# CONVERGENCE REPORT PARTICLE NUMBER: it=%d\n", it);
        nparttest=fabs(npart[SPINA]-md.Na)/(md.Na+md.Nb);
        if(nparttest>md.npartconveps) {is_converged=0; is_converged_local=0;} else {is_converged_local=1;}
        if(iam==0) wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g CONVSTATUS=%6s NPARTCONV=%16.8g\n", 
            "SPINA", npart[SPINA], npart_old[SPINA], (npart[SPINA]-npart_old[SPINA]), convstatus[is_converged_local], nparttest);
        nparttest=fabs(npart[SPINB]-md.Nb)/(md.Na+md.Nb);
        if(nparttest>md.npartconveps) {is_converged=0; is_converged_local=0;} else {is_converged_local=1;}
        if(iam==0) wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g CONVSTATUS=%6s NPARTCONV=%16.8g\n", 
            "SPINB", npart[SPINB], npart_old[SPINB], (npart[SPINB]-npart_old[SPINB]), convstatus[is_converged_local], nparttest);
            
        if(iam==0) wprintf("# CONVERGENCE REPORT ENERGY: it=%d\n", it);
        Effg = 0.6 * (npart[SPINA]+npart[SPINB]) * eF;
        E_tot=0.0; E_tot_old=0.0;
        for(i=0; i<ENERGYITEMS; i++)
        {
            E_tot+=energy[i]; 
            E_tot_old+=energy_old[i];
            
            if(fabs((energy[i]-energy_old[i])/Effg)>md.energyconveps) {is_converged=0; is_converged_local=0;} else {is_converged_local=1;}
            
            if(iam==0) wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g CONVSTATUS=%6s\n", 
                energy_labels[i], energy[i]/Effg, energy_old[i]/Effg, (energy[i]-energy_old[i])/Effg, convstatus[is_converged_local]);
            
        }
        if(iam==0) wprintf("  ------------------------------------------------------------------------------------------\n");
        if(fabs((E_tot-E_tot_old)/Effg)>md.energyconveps) {is_converged=0; is_converged_local=0;} else {is_converged_local=1;}
        if(iam==0) wprintf("%9s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g CONVSTATUS=%6s\n", 
                "E_tot", E_tot/Effg, E_tot_old/Effg, (E_tot-E_tot_old)/Effg, convstatus[is_converged_local]);
        
        double minF_new = E_tot - dc_mu_a*npart[SPINA] - dc_mu_b*npart[SPINB];
        double minF_old = E_tot_old - dc_mu_a_old*npart_old[SPINA] - dc_mu_b_old*npart_old[SPINB];
        if(iam==0) wprintf("# MINIMIZATION FUNCTION: %16.8f\n", minF_new);
        if(iam==0) wprintf("# FUNCTION CHANGED BY: %16.8f\n", minF_new-minF_old);
        for(i=0; i<ENERGYITEMS; i++) observables[i]=energy[i]; observables[ENTROPY]=S;
        if(iam==0) cpu_exec( logger_add_entry(it, densall, potsall, kF, mu, observables, npart, dc_params, dc_extra_data_size, dc_extra_data) );
        cpu_exec( wslda_check_array_against_naninf(ENERGYITEMS, energy) );
        
        // broyden restarting
        i=0;
        if(fabs((E_tot-E_tot_old)/Effg)>md.broydenEmaxchg) i=1;
        if(i==1 && it+1-md.Mbroyden-md.broydenEdelay>md.startbroyden && it+1<md.stopbroyden && md.broyden == 1 && md.broydenautores==1)
        {             //md.Mbroyden+x - typically after a few iteraton once the Broyden starts we observer energy fluctuation that should vanish 
            md.startbroyden=it+1;
            if(iam==0) wprintf("# ENERGY HAS CHANGED MORE THAN broydenEmaxchg=%f! RESTARTING BROYDEN TO AVOID INSTABILITY.\n", md.broydenEmaxchg);
        }
        
        // set constants after update
        wdata_setconst(&wdmd, "kF", kF);
        wdata_setconst(&wdmd, "eF", eF);
        wdata_setconst(&wdmd, "mu_a", mu[SPINA]);
        wdata_setconst(&wdmd, "mu_b", mu[SPINB]);
        file_operation( write_measurments(&wdmd, MPI_COMM_WORLD, "st", it, densall, potsall) );
        if(iam==0) file_operation( write_wdata_metadata_file(&md, &wdmd, "st-wslda-3d") );
                
        // checkpoint - only by iam==0
        if(md.checkpoint && iam==0) 
        {
            // prepare data info for writing
            double twrt_consts[11] = {dc_mu_a, dc_mu_b, dc_mu_a_old, dc_mu_b_old, dc_ec, beta, eF, kF, Effg, npart[SPINA], npart[SPINB]};
            // write checkpoint
            file_operation(
                wslda_st_write_checkpoint(3, it, 11, twrt_consts, POTDIM, h_potentials, DENSDIM, h_densities, ENERGYITEMS, energy, SOLDIM + 2, dens_in, dens_out)
            );
        }
        
        
        rt_other+=e_t(0);
        
        // ------------------ timing------------------
        rt_tot=rt_zheev+rt_dens+rt_pot+rt_other+rt_me+rt_redistrib; 
        if(iam==0) wprintf("# TIMING rt_tot=%8.2f: rt_zheev=%8.2f[%5.2f%%] rt_redistrib=%8.2f[%5.2f%%] rt_dens=%8.2f[%5.2f%%] rt_pot=%8.2f[%5.2f%%] rt_me=%8.2f[%5.2f%%] rt_other=%8.2f[%5.2f%%]\n", 
            rt_tot, rt_zheev, rt_zheev/rt_tot*100., rt_redistrib, rt_redistrib/rt_tot*100., rt_dens, rt_dens/rt_tot*100., rt_pot, rt_pot/rt_tot*100., rt_me, rt_me/rt_tot*100., rt_other, rt_other/rt_tot*100.);
        fflush(stdout);
        
        
        if(saving_iteration==1)
        {
            // Write missing files: info for the total set of wf and file with potentials
            if(iam==0)
            {
                // write info file
                sprintf(file_name, "%s/s3dpca.info", md.outprefix);
                mu[SPINA] = dc_mu_a; mu[SPINB] = dc_mu_b;
                file_operation( create_checkpoint_info_pca(file_name, nwf, NX, NY, NZ, DX, DY, DZ, kF, mu, dc_ec, beta) );
                
                // write potentials
                sprintf(file_name, "%s/s3dpca.pud", md.outprefix);
                file_operation( checkpoint_save_u_and_delta_kzpca(file_name, NX*NY*NZ, potsall.V_a, potsall.delta) );
            }
            
            // Create check.stamp
            if(iam==0)
            {
                // write check.stamp file
                sprintf(file_name, "%s_check.stamp", md.outprefix);
                wprintf("# CREATING CHECK STAMP FILE: `%s`\n",file_name);
                file_operation( touch_file(file_name) );
                file_operation( check_stamp_entry(file_name, 12, NXYZ, h_densities, ENERGYITEMS, energy) );
            }
            
            if(iam==0)
            {
                wprintf("# CREATING WAVE-FUNCTIONS REPRODUCIBILITY PACK: %s/reprowf.tar\n", md.outprefix);
                create_reprowf_tar(extra_data_size);
            }
            
            if(iam==0) wprintf("# SAVING ITERATION DONE.\n");
            break;
        }
        
        it++; // go to next iteration - global counter
        kziter++; // go to next iteration - this run counter
        
        if(is_converged && kziter>0)
        {
            if(iam==0) wprintf("# ALGORITHM CONVERGED!\n");
            
            if(md.writewf==0) break;
            else saving_iteration=1;
        }
        
        if(kziter==md.maxiters)
        {
            if(iam==0) wprintf("# MAXIMUM NUMBER OF ITERATIONS REACHED!\n"); fflush(stdout);
            
            if(md.writewf==0) break;
            else saving_iteration=1;
        }
           
//         break; //NOTE - for tests only
    } // while(1)

    /* messy exit here */
    destroy_fft_plans(&mdfft);
    
#ifdef TESTSUITE
    if(iam==0) testsuite_ok();
#endif
    
    MPI_Barrier( MPI_COMM_WORLD ) ;
    MPI_Finalize() ;
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
