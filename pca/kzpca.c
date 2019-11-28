// This is code for solving stationary DFT equations for polarized cold atoms 

// Authors:
// Gabriel Wlazlowski <gabrielw@if.pw.edu.pl>

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
#include "kzpca_edf.h"
#include "pca_uniform.h"
#include "pca_logger.h"
#include "kzpca_fft.h"
#include "kzpca_me.h"
#include "kzpca_densities.h"

/* ZHEEV prototype */
extern void zheev_( char* jobz, char* uplo, int* n, double complex* a, int* lda,
                double* w, double complex* work, int* lwork, double* rwork, int* info );

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, double complex* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", creal(a[i+j*lda]) );
                printf( "\n" );
        }
}
double u_ext(int ix, int iy, int it, int spin);
void process_params(double *params, double kF);
void modify_potentials(int it, double *h_densities, double *h_potentials);


// ====================================================================================
// =================================== GLOBAL VARIABLES ===============================
// ====================================================================================
// This section reporduces variables keeped in GPU constant memory 
double *dc_params; /* Declaration of the variable */
double dc_mu_a;
double dc_mu_b;
double dc_ec;

// rotating frame
double dc_Omega_a;
double dc_Omega_b;

typedef char * string;

// make -f Makefile.kzsolver

int main( int argc , char ** argv ) 
{    
    int i, j, k; // basic iterators
    int ix, iy, iz, ixyz; // lattice iterators
    int ierr; // error flag
    int ip, np; // basic MPI indicators
    int nwf; // number of wave-functions
    int nwfip; // number of wave-functions per process
    int iwf; // wave-function iterator
    int ikz; // index of kz vector
    // other technical variables
    int *wf_tbl, *wf_idx_tbl; // table of size np, keeps number of managed wf by each process
    int kziter=0;
    int it=0; // global iteration number
    double beta; // inverse of temperature
    double Lz_a, Lz_b, Lz;
    double Lz_a_old=0.0, Lz_b_old=0.0, Lz_old=0.0;
    
    double eF_a, eF_b, eF, Effg, kF;
    double mu[2]; // chemical potential
    double ec; // energy cut-off
    double rt_zheev, rt_dens, rt_pot, rt_other, rt_me, rt_tot=0.0; // run time
    
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
    double npart[2], npart_old[2];
    int is_converged;
    dc_Omega_a=0.0;
    dc_Omega_b=0.0;
    
    char file_name[256];
  
    /* start main */
    MPI_Init( &argc , &argv ) ; /* set up the parallel WORLD */
    MPI_Comm_size( MPI_COMM_WORLD , &np ) ; /* total number of processes */
    MPI_Comm_rank( MPI_COMM_WORLD , &ip ) ; /* id of process st 0 <= ip < np */
    
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
    
    if( ip == 0 )
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
    if(ip==0) printf("# LATTICE: %d x %d x %d\n", NX, NY, NZ);
    
    // ====================================================================================
    // ================================ ALLOCATE CPU BUFFERS ==============================
    // ====================================================================================
    cppmallocl(h_densities,12*NX*NY,double);
    cppmallocl(h_densities_old,12*NX*NY,double);
    cppmallocl(h_densities_partial,12*NX*NY,double);
    cppmallocl(h_potentials,4*NX*NY,double);
    cppmallocl(h_potentials_old,4*NX*NY,double);
    cppmallocl(dc_params, MAX_USER_PARAMS, double);
    
    // For easier access to data
    // densities 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    cppmallocl(h_fbetaEn, 2*NXYZ,double);
    cppmallocl(h_En, 2*NXYZ,double);
    
    // hamitonian
    double complex *h;
    cppmallocl(h, (2*NX*NY)*(2*NX*NY),double complex);
    // eigen-values
    double * En; 
    cppmallocl(En, 2*NX*NY,double);
    
    // working space for zheev
    int n = 2*NX*NY, lda = 2*NX*NY, info, lwork;
    double complex wkopt;
    double complex * work;
    double *rwork;
    cppmallocl(rwork, 3*n-2, double);
    
    if(ip==0) printf("# MATRIX SIZE: %d x %d\n", n, n);
    /* Query and allocate the optimal workspace */
    lwork = -1;
    zheev_( "V", "U", &n, h, &lda, En, &wkopt, &lwork, rwork, &info );
    lwork = (int)creal(wkopt);
    cppmallocl(work, lwork, double complex);
    if(ip==0) printf("# SETTING UP ZHEEV: lwork = %d\n", lwork); fflush(stdout);
    

    // ====================================================================================
    // ================================ CREATE WAVE VECTORS ===============================
    // ====================================================================================
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / NX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / NX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / NY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / NY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / NZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / NZ * ( double ) j ; 
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
            if(ip==0) printf("# CREATING UNIFORM SOLUTION...\n");
            
            // Generate initial state for testing
//             cpu_exec( solve_uniform_problem(md.Na/NXYZ, md.Nb/NXYZ, &nwf, ip==0) ); // do not check
            solve_uniform_problem(md.Na/NXYZ, md.Nb/NXYZ, &nwf, ip==0);
            
            // Save solution
            if(ip==0 && md.init0save)
            {
                cpu_exec( save_uniform() );
//                 ABORT;
            }  
        }
        else 
        {
            if(ip==0) printf("# READING UNIFORM SOLUTION...\n");
            if(ip==0) { cpu_exec( read_uniform(&nwf, ip==0) ); }
            MPI_Bcast( &__md_pca_uniform , sizeof(metadata_pca_uniform_t), MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
            MPI_Bcast( &nwf , 1, MPI_INT , 0 , MPI_COMM_WORLD ) ;
        }
                       
        // initialize potentials and densities
        for(ixyz=0; ixyz<NX*NY; ixyz++)
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
        Effg = 0.6*__md_pca_uniform.n0_a*eF_a*NXYZ + 0.6*__md_pca_uniform.n0_b*eF_b*NXYZ; 
        
        // set global variables
        dc_mu_a=__md_pca_uniform.mu_a;
        dc_mu_b=__md_pca_uniform.mu_b;
        dc_ec = __md_pca_uniform.ec;
        beta = __md_pca_uniform.beta;
    }
    else if(md.inittype==2) // start from checkpoint
    {
        if(ip==0)
        {
            sprintf(file_name, "%s_checkpoint.kzpca", md.inprefix);
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
            fread(h_potentials , sizeof(double)*NX*NY, 4 , pFile);
            fread(h_densities  , sizeof(double)*NX*NY, 12, pFile);
            fread(energy       , sizeof(double)      , 5 , pFile);
            fread(npart        , sizeof(double)      , 2 , pFile);
                  
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
        MPI_Bcast( &it , 1, MPI_INT , 0 , MPI_COMM_WORLD ) ;
        MPI_Bcast(&dc_mu_a     , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&dc_mu_b     , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&dc_ec       , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&beta        , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&eF          , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
        MPI_Bcast(&kF          , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(&Effg        , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(h_potentials , 4*NX*NY , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(h_densities  , 12*NX*NY, MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(energy       , 5 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
        MPI_Bcast(npart        , 2 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    }
    else
    {
        if(ip==0) printf("NOT SUPPORTED INITTYPE=%d!\n", md.inittype);
        ABORT;
    }
    // ===================================================================================
    // =============================== WORK DIVISION =====================================
    // ===================================================================================
    
    // assing number of kz values taken by each process
    if ( np > NZ )
    {
        if(ip==0) printf("ERROR: TOO MUCH RESOURCES! (%d>%d)\n", np, NZ);
        ABORT;
    }
    getnwfip( ip , np , NZ , &nwfip ) ;
    MPI_Gather( &nwfip , 1 , MPI_INT , wf_tbl , 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
    MPI_Bcast( wf_tbl , np , MPI_INT , 0 , MPI_COMM_WORLD ) ; 
    
    // my range of wf to manage
    int myuidx=0, mylidx=0;
    for(i=0; i<=ip; i++)
        myuidx+=wf_tbl[i];
    mylidx=myuidx-nwfip;

    j=0;
    for(i=0; i<ip; i++) j+=wf_tbl[i];    
    MPI_Gather( &j , 1 , MPI_INT , wf_idx_tbl , 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
    MPI_Bcast( wf_idx_tbl , np , MPI_INT , 0 , MPI_COMM_WORLD ) ;
    
    printf("# PROCESS %d COMPUTES FOR %d kz-values [%d,%d)\n", ip, nwfip, mylidx,myuidx); fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // ===================================================================================
    // ================================== FFTW PLANS =====================================
    // ===================================================================================
    omp_set_num_threads(md.nthreads);
    printf("# PROCESS %d ACTIVATES %d THREADS FOR FFTW, BTACH=%d\n", ip, omp_get_max_threads(), md.batch); fflush(stdout);
    metadata_kzpca_fft mdfft; // keeps plans and buffers for fftw
    create_fft_plans(&mdfft, md.batch);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // ===================================================================================
    // ========================= MATRIX ELEMENTS OF GRADIENTS ============================
    // ===================================================================================    
    double complex * me_d_dx;
    double complex * me_d_dy;
    cppmallocl(me_d_dx, NX*NX,double complex);
    cppmallocl(me_d_dy, NY*NY,double complex);
    
    cpu_exec( compute_matrix_elements_of_momentum_operator(NX, me_d_dx) );
    cpu_exec( compute_matrix_elements_of_momentum_operator(NY, me_d_dy) );
        
    // ===================================================================================
    // ======================================= LOGGER ====================================
    // =================================================================================== 
    // Create binary files and add initial measurement
    kF=1.0;
    eF = 0.5*kF*kF;
    Effg = 0.6 * md.Na * eF;
    beta = 1.0 / (md.kztemp * eF);    
    dc_ec = M_PI*M_PI/2.; // TODO - set by hand
    if(ip==0)
    {
        // Create empty files with headers - do it once
        sprintf(file_name, "%s_density_a.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, 1.0*it, 1.0) );
        sprintf(file_name, "%s_density_b.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, 1.0*it, 1.0) );
        sprintf(file_name, "%s_delta.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, 1.0*it, 1.0) );    
        sprintf(file_name, "%s_current_a.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, 1.0*it, 1.0) );
        sprintf(file_name, "%s_current_b.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, 1.0*it, 1.0) );    
        
        // for each measurement add data to file
        sprintf(file_name, "%s_density_a.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, rho_a, sizeof(double)*NX*NY) );
        sprintf(file_name, "%s_density_b.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, rho_b, sizeof(double)*NX*NY) );
        sprintf(file_name, "%s_delta.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, delta, sizeof(double complex)*NX*NY) );
        sprintf(file_name, "%s_current_a.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, j_a_x, sizeof(double)*NX*NY*3) );
        sprintf(file_name, "%s_current_b.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, j_b_x, sizeof(double)*NX*NY*3) );
        
        // Create run log and add entry
        mu[SPINA]=dc_mu_a; mu[SPINB]=dc_mu_b;
        cpu_exec( create_header_of_runlog(execcmd, kF, Effg, mu, dc_ec, nwf, np, nwfip) );
    }
         
    // ===================================================================================
    // ========================== SELF-CONSITENT LOOP ====================================
    // ===================================================================================     
    it=0; // TODO - reset
    while(1) // do until reached self-consitency
    {
        b_t();
        rt_zheev=0.0; rt_dens=0.0; rt_pot=0.0; rt_other=0.0; rt_me=0.0;
        // Make copy of potentials and densities
        for(ixyz=0; ixyz< 4*NX*NY; ixyz++) h_potentials_old[ixyz] = h_potentials[ixyz];
        for(ixyz=0; ixyz<12*NX*NY; ixyz++) h_densities_old[ixyz]  = h_densities[ixyz];
        for(ixyz=0; ixyz<12*NX*NY; ixyz++) h_densities_partial[ixyz] = 0.0; // reset
        for(i=0; i< 5; i++) energy_old[i] = energy[i]; // make copy
        npart_old[SPINA]=npart[SPINA]; npart_old[SPINB]=npart[SPINB]; // make copy 
        Lz_a_old=Lz_a; Lz_b_old=Lz_b; Lz_old=Lz; 
        
        // take density in the center and use it for definition of the kF (for SPINA)
        kF = pow(6.*M_PI*M_PI*rho_a[NY/2 + NX/2*NY],1./3.);
        // kF = 1.000; // TODO: overwrite
        eF = 0.5 * kF * kF;
        beta = 1.0 / (md.kztemp * eF);
        if(ip==0) printf("# EXECUTING: process_params(md.params, %f)\n", kF);
        for(i=0; i<MAX_USER_PARAMS; i++) dc_params[i]=md.params[i];
        process_params(dc_params, kF);
        
	/*
        // ajusting particle number
        double tkF = 1.0; // TODO
        double tol = 0.005*tkF;
        kF = pow(6.*M_PI*M_PI*rho_a[NY/2 + NX/2*NY],1./3.);
        double diff = kF-tkF;
        if(fabs(diff)>tol) 
        {
            if(diff>0.0) md.Na-=1.0; 
            else         md.Na+=1.0;
            if(ip==0) printf("# CHANGING PARTICLE NUMBER TO md.Na=%f\n", md.Na);
        }
	*/
        md.Nb=md.Na;
        if(ip==0) printf("# LOCAL FERMI MOMENTUM kF=%f\n", kF);
        kF=1.0; // set by hand to 
        
        // plot potential along x and y direction
        sprintf(file_name, "%s_potential.txt", md.outprefix);
        FILE * fpot = fopen(file_name, "w");
        fprintf(fpot, "# ix u_ext(ix,NY/2,it,SPINA) u_ext(ix,NY/2,it,SPINB)\n");
        for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) if(iy==NY/2) fprintf(fpot, "%8d %16.8g %16.8g\n", ix, u_ext(ix,iy,it,SPINA), u_ext(ix,iy,it,SPINB));
        fprintf(fpot, "\n# iy u_ext(NX/2,iy,it,SPINA) u_ext(NX/2,iy,it,SPINB)\n");
        for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) if(ix==NX/2) fprintf(fpot, "%8d %16.8g %16.8g\n", iy, u_ext(ix,iy,it,SPINA), u_ext(ix,iy,it,SPINB));

        fclose(fpot);
        rt_other+=e_t(0);
        
        // ------------------ diagonalize for each kz ------------------
        nwf=0;
        for(ikz=mylidx; ikz<myuidx; ikz++) // for each kz
        {
            // matrix elements
            b_t();
            cpu_exec( compute_matrix_elements(it, h_densities, h_potentials, &mdfft, h, kkz[ikz], me_d_dx, me_d_dy) );
            rt_me+=e_t(0);
            
            // diagonalize
            b_t();
            zheev_( "V", "U", &n, h, &lda, En, work, &lwork, rwork, &info );
            /* Check for convergence */
            if( info > 0 ) 
            {
                printf( "The algorithm failed to compute eigenvalues!\n" );
                ABORT_NOBARRIER;
            }
            rt_zheev+=e_t(0);
            
            // compute contribution to the densities
            b_t();
            cpu_exec( compute_contribution_to_densities(En, h, dc_ec, beta, h_densities_partial, &mdfft, kkz[ikz], &nwf) );
            rt_dens+=e_t(0);
            
        } // for(ikz=mylidx; ikz<myuidx; ikz++)

        
        // global reduction
        b_t();
        MPI_Allreduce( h_densities_partial, h_densities, 12*NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce( MPI_IN_PLACE, &nwf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if(ip==0) printf("# NWF=%d\n", nwf);
        rt_dens+=e_t(0);
        
        // ------------------ mix densities ------------------
        b_t();
        for(ixyz=0; ixyz<12*NX*NY; ixyz++) h_densities[ixyz] = md.kzmixparam * h_densities[ixyz] + (1.0-md.kzmixparam) * h_densities_old[ixyz];
        
        // ------------------ update chemical potentials ------------------
        if(ip==0) printf("# MUCHNAGE FROM: dc_mu_a=%16.8g  dc_mu_b=%16.8g\n", dc_mu_a, dc_mu_b);
        if(it>0)
        {       
            npart[SPINA]=0.0; npart[SPINB]=0.0;
            for(ixyz=0; ixyz<NX*NY; ixyz++) {npart[SPINA]+=rho_a[ixyz]*NZ; npart[SPINB]+=rho_b[ixyz]*NZ;}
	    double kzmuchange_a = md.kzmuchange*(npart[SPINA] - md.Na);
	    double kzmuchange_b = md.kzmuchange*(npart[SPINB] - md.Nb);
	    // double leF_a = pow(6.*M_PI*M_PI*rho_a[NY/2 + NX/2*NY],2./3.) / 2.0;
	    // double leF_b = pow(6.*M_PI*M_PI*rho_b[NY/2 + NX/2*NY],2./3.) / 2.0;
	    double max_change = 0.1;
	    if(fabs(kzmuchange_a)>max_change)
	    {
	      if(kzmuchange_a>0.0) kzmuchange_a=     max_change;
	      else                 kzmuchange_a=-1.0*max_change;
	    }
            dc_mu_a -= kzmuchange_a;
            dc_mu_b -= kzmuchange_b;  
            dc_mu_b=dc_mu_a; // TODO
        }
        if(ip==0) printf("# MUCHNAGE TO  : dc_mu_a=%16.8g  dc_mu_b=%16.8g\n", dc_mu_a, dc_mu_b);
        rt_other+=e_t(0);
        
        // ------------------ compute new potentials ------------------
        b_t();
        cpu_exec( recompute_potentials(it, h_densities, h_potentials, h_potentials) ); 
        modify_potentials(it, h_densities, h_potentials) ;
        rt_pot+=e_t(0);
        
        // ------------------ angular momentum ------------------
        b_t();
        cpu_exec( compute_angular_momentum_Lz(j_a_x, j_a_y, &Lz_a) );
        cpu_exec( compute_angular_momentum_Lz(j_b_x, j_b_y, &Lz_b) );
        Lz = Lz_a + Lz_b; // total angular momentum
        if(ip==0) printf("# ANGULAR MOMENTUM: it=%d\n", it);
        if(ip==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "LZ_A", Lz_a/npart[SPINA], Lz_a_old/npart_old[SPINA], (Lz_a/npart[SPINA]-Lz_a_old/npart_old[SPINA]));
        if(ip==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "LZ_B", Lz_b/npart[SPINB], Lz_b_old/npart_old[SPINB], (Lz_b/npart[SPINB]-Lz_b_old/npart_old[SPINB]));  
        if(ip==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "LZ", Lz/(npart[SPINA]+npart[SPINB]), Lz_old/(npart_old[SPINA]+npart_old[SPINB]), (Lz/(npart[SPINA]+npart[SPINB])-Lz_old/(npart_old[SPINA]+npart_old[SPINB])));
        if(ip==0) printf("dc_Omega_a=%16.8g  dc_Omega_b=%16.8g\n",dc_Omega_a, dc_Omega_b);
        
        // ------------------ check convergence ------------------
        cpu_exec( compute_energy(it, h_densities, h_potentials, energy, npart) );
        if(ip==0) printf("# PARTICLE NUMBER: it=%d\n", it);
        if(ip==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "SPINA", npart[SPINA], npart_old[SPINA], (npart[SPINA]-npart_old[SPINA]));
        if(ip==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "SPINB", npart[SPINB], npart_old[SPINB], (npart[SPINB]-npart_old[SPINB]));
            
        is_converged=1;
        if(ip==0) printf("# CONVERGENCE REPORT: it=%d\n", it);
        E_tot=0.0; E_tot_old=0.0;
        for(i=0; i<5; i++)
        {
            E_tot+=energy[i]; 
            E_tot_old+=energy_old[i];
            
            if(ip==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
                energy_labels[i], energy[i]/Effg, energy_old[i]/Effg, (energy[i]-energy_old[i])/Effg);
            
            if(fabs((energy[i]-energy_old[i])/Effg)>md.kzconveps) is_converged=0;
        }
        if(ip==0) printf("  ------------------------------------------------------------------------\n");
        if(ip==0) printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
                "E_tot", E_tot/Effg, E_tot_old/Effg, (E_tot-E_tot_old)/Effg);
        if(ip==0)
        {
            #define OUTPUT_ENTRIES 15
            double line_items[OUTPUT_ENTRIES]={     
                npart[SPINA], // 1
                npart[SPINB], // 2
                npart[SPINA]+npart[SPINB], // 3
                E_tot/Effg, // 4
                energy[0]/Effg, // 5
                energy[1]/Effg, // 6
                energy[2]/Effg, // 7
                energy[3]/Effg, // 8
                energy[4]/Effg, //9
                dc_mu_a, //10
                dc_mu_b, //11
                beta, // 12
                Lz_a/npart[SPINA], // 13
                Lz_b/npart[SPINB], // 14
                Lz/(npart[SPINA]+npart[SPINB]) ,  // 15
            };
            cpu_exec( add_line_to_file(it, rt_tot, OUTPUT_ENTRIES, line_items) );
            
            // for each measurement add data to file
            sprintf(file_name, "%s_density_a.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, rho_a, sizeof(double)*NX*NY) );
            sprintf(file_name, "%s_density_b.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, rho_b, sizeof(double)*NX*NY) );
            sprintf(file_name, "%s_delta.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, delta, sizeof(double complex)*NX*NY) );
            sprintf(file_name, "%s_current_a.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, j_a_x, sizeof(double)*NX*NY*3) );
            sprintf(file_name, "%s_current_b.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, j_b_x, sizeof(double)*NX*NY*3) );
        }
        
        // checkpoint - only by ip==0
        if(md.checkpoint && ip==0)
        {
            sprintf(file_name, "%s_checkpoint.kzpca", md.outprefix);
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
            fwrite(h_potentials , sizeof(double)*NX*NY, 4 , pFile);
            fwrite(h_densities  , sizeof(double)*NX*NY, 12, pFile);
            fwrite(energy       , sizeof(double)      , 5 , pFile);
            fwrite(npart        , sizeof(double)      , 2 , pFile);
                  
            fclose(pFile);
        }
        
        rt_other+=e_t(0);
        
        if(is_converged && kziter>0)
        {
            if(ip==0) printf("# ALGORITHM CONVERGED!\n");
            break;
        }
        
        it++; // go to next iteration
        kziter++;
        if(kziter==md.kzmaxiters)
        {
            if(ip==0) printf("# MAXIMUM NUMBER OF ITERATIONS REACHED!\n"); fflush(stdout);
            break;            
        }
        
        // ------------------ timing------------------
        rt_tot=rt_zheev+rt_dens+rt_pot+rt_other+rt_me; 
        if(ip==0) printf("# TIMING rt_tot=%8.2f: rt_zheev=%8.2f[%5.2f%%] rt_dens=%8.2f[%5.2f%%] rt_pot=%8.2f[%5.2f%%] rt_me=%8.2f[%5.2f%%] rt_other=%8.2f[%5.2f%%]\n", 
            rt_tot, rt_zheev, rt_zheev/rt_tot*100., rt_dens, rt_dens/rt_tot*100., rt_pot, rt_pot/rt_tot*100., rt_me, rt_me/rt_tot*100., rt_other, rt_other/rt_tot*100.);
        fflush(stdout);
        
    } // while(1)
    
    
    if(is_converged || kziter==md.kzmaxiters)
    {
        if(ip==0) printf("# WRITING WAVE-FUNCTIONS TO THE  FILE!\n");
        rt_zheev=0.0; rt_dens=0.0; rt_pot=0.0; rt_other=0.0; rt_me=0.0;
        MPI_Status MPIStat;
        int lastwf=0;
        
        b_t();
        // set new reference energy
        // take density in the center and use it for definition of the kF
        kF = pow(6.*M_PI*M_PI*rho_a[NY/2 + NX/2*NY],1./3.);    
        kF = 1.00; // TODO: overwrited
        eF = 0.5*kF*kF;
        Effg = 0.6 * npart[SPINA] * eF;
        beta = 1.0 / (md.kztemp * eF);
        
        for(ixyz=0; ixyz<12*NX*NY; ixyz++) h_densities_partial[ixyz] = 0.0; // reset
        for(i=0; i< 5; i++) energy_old[i] = energy[i]; // make copy
        
        // clear files
        if(ip==0)
        {
            sprintf(file_name, "%s_kzpca.wfu", md.outprefix);
            file_operation( touch_file(file_name) );
           
            sprintf(file_name, "%s_kzpca.wfv", md.outprefix);
            file_operation( touch_file(file_name) );
            
            sprintf(file_name, "%s_kzpca.kkz", md.outprefix);
            file_operation( touch_file(file_name) );

            sprintf(file_name, "%s_kzpca.en", md.outprefix);
            file_operation( touch_file(file_name) );
            
            sprintf(file_name, "%s_kzpca.info", md.outprefix);
            file_operation( touch_file(file_name) );
            
            sprintf(file_name, "%s_kzpca.pud", md.outprefix);
            file_operation( touch_file(file_name) );
        }
        rt_other+=e_t(0);
        
        // ------------------ diagonalize for each kz ------------------
        nwf=0;
        for(ikz=0; ikz<wf_tbl[0]; ikz++) // for each kz
        {
            i = ikz + mylidx;
            if(i<myuidx) // I have job to do
            {
                // matrix elements
                b_t();
                cpu_exec( compute_matrix_elements(it, h_densities, h_potentials, &mdfft, h, kkz[i], me_d_dx, me_d_dy) );
                rt_me+=e_t(0);
                
                // diagonalize
                b_t();
                zheev_( "V", "U", &n, h, &lda, En, work, &lwork, rwork, &info );
                /* Check for convergence */
                if( info > 0 ) 
                {
                    printf( "The algorithm failed to compute eigenvalues!\n" );
                    ABORT_NOBARRIER;
                }
                rt_zheev+=e_t(0);
                
                // compute contribution to the densities
                b_t();
                cpu_exec( compute_contribution_to_densities(En, h, dc_ec, beta, h_densities_partial, &mdfft, kkz[i], &nwf) );
                rt_dens+=e_t(0);  
            }
          
            b_t();
            
            // wait for my turn of writing event
            if(ip!=0) // Wait until ip-1 process finish his job
                j = MPI_Recv(&lastwf, 1, MPI_INT, ip-1, 99, MPI_COMM_WORLD, &MPIStat);    
            
            // save my wave-functions if generated
            nwf=0;
            if(i<myuidx) // I have job to do
            {
                file_operation( append_wf_from_kzpca(md.outprefix, En, h, dc_ec, beta, kkz[i], &nwf) );
            }
            lastwf+=nwf;
            
            if(ip!=(np-1)) // File is free, send info to next process
                j = MPI_Send(&lastwf, 1, MPI_INT, ip+1, 99, MPI_COMM_WORLD); 
            
            MPI_Bcast( &lastwf , 1, MPI_INT , np-1 , MPI_COMM_WORLD ) ;
            rt_other+=e_t(0);
        }

        // write missing stuff
        b_t();
        if(ip==0)
        {
            printf("# WRITTEN NWF=%d wave-functions\n", lastwf);
            
            // write info file
            sprintf(file_name, "%s_kzpca.info", md.outprefix);
            mu[SPINA] = dc_mu_a; mu[SPINB] = dc_mu_b;
            file_operation( create_checkpoint_info_pca(file_name, lastwf, NX, NY, NZ, 1.0, 1.0, 1.0, kF, mu, dc_ec, beta) );
            
            // write potentials
            sprintf(file_name, "%s_kzpca.pud", md.outprefix);
            file_operation( checkpoint_save_u_and_delta_kzpca(file_name, NX*NY, V_a, delta) );
        }
        rt_other+=e_t(0);
        
        // global reduction
        b_t();
        MPI_Allreduce( h_densities_partial, h_densities, 12*NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        rt_dens+=e_t(0);
        
        // ------------------ angular momentum ------------------
        b_t();
        cpu_exec( compute_angular_momentum_Lz(j_a_x, j_a_y, &Lz_a) );
        cpu_exec( compute_angular_momentum_Lz(j_b_x, j_b_y, &Lz_b) );
        Lz = Lz_a + Lz_b; // total angular momentum
        
        // ------------------ check convergence ------------------
        cpu_exec( compute_energy(it, h_densities, h_potentials, energy, npart) );
        if(ip==0) printf("# PARTICLE NUMBER: it=%d\n", it);
        if(ip==0) printf("%8s: NEW=%16.8g OLD=none DIFF=none\n", 
            "SPINA", npart[SPINA]);
        if(ip==0) printf("%8s: NEW=%16.8g OLD=none DIFF=none\n", 
            "SPINB", npart[SPINB]);
            

        if(ip==0) printf("# ENERGY REPORT: it=%d\n", it);
        E_tot=0.0; E_tot_old=0.0;
        for(i=0; i<5; i++)
        {
            E_tot+=energy[i]; 
            E_tot_old+=energy_old[i];
            
            if(ip==0) printf("%8s: NEW=%16.8g OLD=none DIFF=none\n", 
                energy_labels[i], energy[i]/Effg);
        }
        if(ip==0) printf("  ------------------------------------------------------------------------\n");
        if(ip==0) printf("%8s: NEW=%16.8g OLD=none DIFF=none\n", 
                "E_tot", E_tot/Effg);
        if(ip==0)
        {
            #define OUTPUT_ENTRIES 15
            double line_items[OUTPUT_ENTRIES]={     
                npart[SPINA], // 1
                npart[SPINB], // 2
                npart[SPINA]+npart[SPINB], // 3
                E_tot/Effg, // 4
                energy[0]/Effg, // 5
                energy[1]/Effg, // 6
                energy[2]/Effg, // 7
                energy[3]/Effg, // 8
                energy[4]/Effg, //9
                dc_mu_a, //10
                dc_mu_b, //11
                beta, // 12
                Lz_a/npart[SPINA], // 13
                Lz_b/npart[SPINB], // 14
                Lz/(npart[SPINA]+npart[SPINB]) ,  // 15
            };
            cpu_exec( add_line_to_file(it, rt_tot, OUTPUT_ENTRIES, line_items) );
            
            // for each measurement add data to file
            sprintf(file_name, "%s_density_a.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, rho_a, sizeof(double)*NX*NY) );
            sprintf(file_name, "%s_density_b.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, rho_b, sizeof(double)*NX*NY) );
            sprintf(file_name, "%s_delta.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, delta, sizeof(double complex)*NX*NY) );
            sprintf(file_name, "%s_current_a.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, j_a_x, sizeof(double)*NX*NY*3) );
            sprintf(file_name, "%s_current_b.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, j_b_x, sizeof(double)*NX*NY*3) );
            
            // write full density array
            sprintf(file_name, "%s_kzpca.dens", md.outprefix);
            if(ip==0) printf("# Writing density array to file: `%s`\n", file_name);
            FILE * fdensarray = fopen(file_name, "wb");
            fwrite(h_densities  , sizeof(double)*NX*NY, 12, fdensarray);
            fclose(fdensarray);
        }
        rt_other+=e_t(0);
        
        // ------------------ timing------------------
        rt_tot=rt_zheev+rt_dens+rt_pot+rt_other+rt_me; 
        if(ip==0) printf("# TIMING rt_tot=%8.2f: rt_zheev=%8.2f[%5.2f%%] rt_dens=%8.2f[%5.2f%%] rt_pot=%8.2f[%5.2f%%] rt_me=%8.2f[%5.2f%%] rt_other=%8.2f[%5.2f%%]\n", 
            rt_tot, rt_zheev, rt_zheev/rt_tot*100., rt_dens, rt_dens/rt_tot*100., rt_pot, rt_pot/rt_tot*100., rt_me, rt_me/rt_tot*100., rt_other, rt_other/rt_tot*100.);
        fflush(stdout);
    }
        
    /* messy exit here */
    destroy_fft_plans(&mdfft);
    
    MPI_Barrier( MPI_COMM_WORLD ) ;
    MPI_Finalize() ;
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
