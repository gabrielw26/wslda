// This is code for solving DFT equations for polarized cold atoms (pca)
// Version constrained polarized cold atoms (cpca)
// Along z and y direction translational symmetry is assumeed.
// 
// Authors:
// Gabriel Wlazlowski <gabrielw@if.pw.edu.pl>

#define TDWSLDA_MAIN

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
#include "cpca_derivative.h"
#include "ccpca_kernels.h"
#include "ccpca_dens.h"
#include "pca_utils.h"
#include "wslda_potdens.h"
#include "wslda_wavevectors.h"
#include "pca_io.h"
#include "pca_uniform.h"
#include "pca_logger.h"
#include "cpca_checkpoint.h"
#include "wslda_writevars.h"
#include "wslda_reproducibility.h"

static double dc_ec;
static double dc_t0;
static int dc_np;
static int dc_nwfip;
#include "logger.h"

int main( int argc , char ** argv ) 
{
    int i, j, k; // basic iterators
    int forceCP = 0; // force checkpointing flag
    int ix, iy, iz, ixyz; // lattice iterators
    int ierr; // error flag
    int ip, np; // basic MPI indicators
    int nwf; // number of wave-functions
    int nwfip; // number of wave-functions per process
    int iwf; // wave-function iterator
    
    // physical variables
    double mu[2]; // chemical potential
    double ec; // energy cut-off
    double t0=0.0, dt; 
    int it=0;
    double kF;
    double beta;

    int HowMany = 24;
#ifdef MPI_NP_PER_IO_GROUP
    HowMany=MPI_NP_PER_IO_GROUP;
#endif
    int gradients_computed = 1; // flag indicating if code uses gradients in computaation, by default equal 1
    
    // arrays
    double complex *h_wavefun; // pointer to wave-functions on host (cpu) side 
    double *h_fbetaEn; // pointer to weights of wave-functions on host (cpu) side
    double *d_fbetaEn; // pointer to weights of wave-functions on device (gpu) side
    double *h_kkyz; // pointer to ky and kz wave-vector values on host (cpu) side
    double *d_kkyz; // pointer to ky and kz wave-vector values on device (gpu) side
    double *h_En; // pointer to eigen-values of wave-functions on host (cpu) side
    double *h_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *d_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (GPU)
    double *h_potentials; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    double *d_potentials; // pointer to array with potentials [V_a, V_b, delta] (GPU)
    double *d_workarea; // pointer to working area, also used by cufft (GPU)
    double *h_energy; // buffer for energies (CPU)
    cufftDoubleComplex *d_wf; // pointer to wave-functions (GPU)
    cufftDoubleComplex *d_fkm1; // pointer to f_k-1 (GPU)
    cufftDoubleComplex *d_fkm2; // pointer to f_k-2 (GPU)
    cufftDoubleComplex *d_fkm3; // pointer to f_k-3 (GPU)
#if INTEGRATION_SCHEME==AB4AM5
    cufftDoubleComplex *d_fkm4; // pointer to f_k-3 (GPU)
#endif
    cufftDoubleComplex *d_wf_d_dx; // pointer to derivative of wave-function d/dx (GPU)
    cufftDoubleComplex *d_wf_laplace; // pointer to laplace of wave-function (d^2/dx^2 + d^2/dy^2 + d^2/dz^2) (GPU)
    cufftDoubleComplex *d_alphawf_laplace=NULL; // pointer to laplace of alpha*wave-function (d^2/dx^2 + d^2/dy^2 + d^2/dz^2) (GPU)
    cufftDoubleComplex *d_tmp_ptr;
    double *h_qpe_nwfip, *h_qpe_nwf; // buffers for quasiparticle energies

    // other technical variables
    int *wf_tbl, *wf_idx_tbl; // table of size np, keeps number of managed wf by each process
    size_t  workarea_size=(size_t)PCA_WORKSPACE_SHIFT*NX*sizeof(double)*2; // minimal size of workarea 
    int mpipackagesize;
    
    void *extra_data = NULL, *d_extra_data = NULL;
    size_t extra_data_size;
    
    // for reporting
    double eF_a, eF_b, eF, Effg;
    double time;
    double energy_kin, energy_pot, energy_pair, energy_current, energy_uext, energy_dext, energy_vext, energy_tot;   
    double Na, Nb;
    double Laz, Lbz; // angular momentum
    char file_name[256];
    
    // timing
    double rt;    
    
    // quantum friction
    double qfalpha=0.0;
    
    // current corrections
    double cccoeff=0.0;
    
    /* start main */
    wt_b_t(); // tag init time
    
    MPI_Init( &argc , &argv ) ; /* set up the parallel WORLD */
    MPI_Comm_size( MPI_COMM_WORLD , &np ) ; /* total number of processes */
    MPI_Comm_rank( MPI_COMM_WORLD , &ip ) ; /* id of process st 0 <= ip < np */
    
    if(ip==0) printf("# CODE: TD-WSLDA-1D\n");

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
       
        // Make copy of input file
        sprintf(file_name, "%s_input.txt", md.outprefix);
        file_operation( copy_input_file(argv[i],file_name) ); 
        file_operation( assure_reproducibility(md.outprefix) );
        file_operation( create_directory(md.outprefix) );
    }
    
    // Broadcast input parameter
    MPI_Bcast( &md , sizeof(md) , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
    
    // variables
    dt= md.dt ;
#ifndef TDWSLDA
    double Emax =  M_PI*M_PI/2.; // E_max = p_max^2 / 2m, where: p_max is maximum momentum on the lattice, p_max=M_PI (if lattice spacing is 1.0)
    dt/=Emax; // time step
#endif
    
#ifdef SPINSYMMETRY_MODE
    if(ip==0 && md.spinsymmetry==0) print_warning(WSLDA_WRN_SPINSYMMETRY0);
    md.spinsymmetry=1;
#else
    if(ip==0 && md.spinsymmetry==1) print_warning(WSLDA_WRN_SPINSYMMETRY1);
    md.spinsymmetry=0;
#endif
    
#ifdef UNIFORM_TEST_MODE
    md.Na = ceil(1.0/(6.0*M_PI*M_PI) * LXYZ);
    
    if(md.spinsymmetry==1) md.Nb = md.Na;
    else md.Nb = md.Na +1; 

    if(ip==0) printf("# UNIFORM_TEST_MODE: SETTING NUMBER OF PARTICLES Na=%f Nb=%f\n", md.Na, md.Nb);
    md.init0Na=md.Na; md.init0Nb=md.Nb;
#endif
    
    if(ip==0) printf("# MPI EXCHANGE PACKAGE SIZE=%.3f MB [%d]\n", 1.0*EXCHANGE_SIZE*NX*sizeof(double)/pow(2,20), EXCHANGE_SIZE);
    
    cpu_exec( wslda_check_settings() );
    
    // ====================================================================================
    // ============================= INITIALIZE GPU =======================================
    // ====================================================================================    
    int deviceId=0;
#ifdef CUSTOM_GPU_DISTRIBUTION
    if(ip==0) printf("# EXECUTING `assign_deviceid_to_mpi_process` TO GET GPUS DISTRIBUTION ACROSS THE SYSTEM\n");
    deviceId = assign_deviceid_to_mpi_process(MPI_COMM_WORLD);
#else
    if(ip==0) printf("# ASSUMING STANDARD DISTRIBUTION OF GPUS ACROSS THE SYSTEM [gpuspernode=%d]\n", md.gpuspernode);
    deviceId = ip % md.gpuspernode;
#endif
    
#ifdef PRINT_GPU_DISTRIBUTION
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("# PROCESS ip=%d RUNNING ON NODE %s USES device-id=%d\n", ip, processor_name, deviceId);
#endif
    gpu_exec( set_gpu(deviceId) );
        
    if(ip==0) printf("# LATTICE: %d x %d x %d\n", NX, NY, NZ);
    if(ip==0) printf("# SPACING: %f x %f x %f\n", DX, DY, DZ);
#if FUNCTIONAL==BDG
    if(ip==0) printf("# ENERGY DENSITY FUNCTIONAL: BDG\n");
#elif FUNCTIONAL==SLDA    
    if(ip==0) printf("# ENERGY DENSITY FUNCTIONAL: SLDA\n");
#elif FUNCTIONAL==ASLDA    
    if(ip==0) printf("# ENERGY DENSITY FUNCTIONAL: ASLDA\n");    
#elif FUNCTIONAL==CUSTOMEDF    
    if(ip==0) printf("# ENERGY DENSITY FUNCTIONAL: CUSTOMEDF\n"); 
#endif
    
    // ====================================================================================
    // ======================== ALLOCATE GPU AND CPU BUFFERS ==============================
    // ====================================================================================
    gpu_exec( host_malloc_pl((size_t)12*NX*sizeof(double), (void **)&h_densities) );
    gpu_exec(     gpu_malloc((size_t)12*NX*sizeof(double), (void **)&d_densities) );
    
    // potentials
    gpu_exec( host_malloc_pl((size_t)12*NX*sizeof(double), (void **)&h_potentials) ); // FIXME: only 4*NX is in use!
    for(i=0; i<12*NX; i++) h_potentials[i]=0.0; // reset
    gpu_exec(     gpu_malloc((size_t)4*NX*sizeof(double), (void **)&d_potentials) );   
    
    // energy
    gpu_exec( host_malloc_pl((size_t)TDWSLDAITEMS*sizeof(double), (void **)&h_energy) );
    
    // For easier access to data
    wslda_density densall = convert_into_wslda_density(h_densities, NX);
    wslda_potential potsall = convert_into_wslda_potential(h_potentials, NX, mu);
        
    // ====================================================================================
    // ================================ INITIAL STATE =====================================
    // ====================================================================================
    if(md.inittype==0 || md.inittype==10) // Start from uniform solution
    {
        if(md.inittype==0)
        {
            if(ip==0) printf("# CREATING UNIFORM SOLUTION...\n");
            
            // Generate initial state for testing
#ifdef BDG_MODE
            cpu_exec( solve_uniform_problem_bdg(md.init0Na/LXYZ, md.init0Nb/LXYZ, &nwf, ip==0) );
#else
            cpu_exec( solve_uniform_problem(md.init0Na/LXYZ, md.init0Nb/LXYZ, &nwf, ip==0) );
#endif
            cpu_exec( get_nwf_to_evolve_1d(&nwf) ); // correct number of states to evolve
            
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
            if(ip==0) { cpu_exec( read_uniform(&nwf, ip==0) ); cpu_exec( get_nwf_to_evolve_1d(&nwf) ); }
            MPI_Bcast( &__md_pca_uniform , sizeof(metadata_pca_uniform_t), MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
            MPI_Bcast( &nwf , 1, MPI_INT , 0 , MPI_COMM_WORLD ) ;
        }
        
        // divide wf over processes
        if(ip==0) printf("# INIT0-1: nwf=%d wave-functions to scatter\n", nwf);
//         ABORT;
        if ( np > nwf )
        {
            if(ip==0) printf("# INIT0-1: np[%d] > nwf[%d]!\n", np, nwf);
            ABORT;
        }
        getnwfip( ip , np , nwf , &nwfip ) ;
        MPI_Gather( &nwfip , 1 , MPI_INT , wf_tbl , 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
        MPI_Bcast( wf_tbl , np , MPI_INT , 0 , MPI_COMM_WORLD ) ; 
        
        // my range of wf to manage
        int myuidx=0, mylidx=0;
        for(i=0; i<=ip; i++)
            myuidx+=wf_tbl[i];
        mylidx=myuidx-nwfip;
        
        // allocate memory for my wf
        cppmallocl(h_wavefun, NX*nwfip*2,double complex);
        cppmallocl(h_fbetaEn, nwfip,double);
        cppmallocl(h_kkyz, nwfip*2,double); // 2 accounts that we have ky and kz
        cppmallocl(h_En, nwfip,double);
        
        // initialize wf
        cpu_exec( create_uniform_wf_1d(mylidx, myuidx, h_wavefun, &mu[SPINA], &mu[SPINB], &ec, h_fbetaEn, h_kkyz, h_En, ip==0) );
        
        // initialize potentials
        for(ixyz=0; ixyz<NX; ixyz++)
        {
            potsall.V_a[ixyz]=__md_pca_uniform.V_a; 
            potsall.V_b[ixyz]=__md_pca_uniform.V_b; 
            potsall.delta[ixyz]=__md_pca_uniform.delta + I*0.0;
        }
        
        // set eF, kF and Effg
        eF=pow(3.0*M_PI*M_PI*(__md_pca_uniform.n0_a+__md_pca_uniform.n0_b), 2.0/3.0) / 2.0;
        kF=pow(3.0*M_PI*M_PI*(__md_pca_uniform.n0_a+__md_pca_uniform.n0_b), 1.0/3.0);
        eF_a=pow(6.0*M_PI*M_PI*__md_pca_uniform.n0_a, 2.0/3.0) / 2.0;
        eF_b=pow(6.0*M_PI*M_PI*__md_pca_uniform.n0_b, 2.0/3.0) / 2.0;
        Effg = 0.6*__md_pca_uniform.n0_a*eF_a*LXYZ + 0.6*__md_pca_uniform.n0_b*eF_b*LXYZ; 
    }
    else if(md.inittype==5) 
    {
        // allocate memory for my wf
        cpu_exec( load_nwf (MPI_COMM_WORLD, md.inprefix, &nwf, &nwfip, HowMany) );
        cppmallocl(h_wavefun, NX*nwfip*2,double complex);
        cppmallocl(h_fbetaEn, nwfip,double);
        cppmallocl(h_kkyz, nwfip*2,double); // 2 accounts that we have ky and kz
//         printf("# WF SCATTER: ip=%d processes nwfip=%d wave-functions\n", ip, nwfip);

    }        
    else if(md.inittype==1) // Start from solution of st-wslda-1d solver 
    {
        // Load data from info file
        int _nx, _ny, _nz;
        double _dx, _dy, _dz;
        int *nwf_per_kyz;
        int nwf_s1dpca;
        int kvecs_to_consder = 0;    
        double * kkx , * kky , * kkz ; /* values */
        cppmallocl(kkx,NX,double); create_kkx(kkx);
        cppmallocl(kky,NY,double); create_kky(kky);
        cppmallocl(kkz,NZ,double); create_kkz(kkz);
        count_number_of_k_modes(kkx, kky, kkz, 1, &kvecs_to_consder);
        wslda_kmode *kvecs;
        cppmallocl(kvecs,kvecs_to_consder,wslda_kmode);
        create_k_modes(kkx, kky, kkz, 1, kvecs);
        cppmallocl(nwf_per_kyz, kvecs_to_consder,int);
        time=0.0;
        sprintf(file_name, "%s/s1dpca.info", md.inprefix);
        if(ip==0)
        {
            file_operation( read_checkpoint_info_pca(file_name, &nwf, &_nx, &_ny, &_nz, &_dx, &_dy, &_dz, &kF, &mu[0], &ec, &beta) );
            if(_nx!=NX || _ny!=NY || _nz!=NZ || _dx!=DX || _dy!=DY || _dz!=DZ)
            {
                printf("# ST-WSLDA-1D INFO FILE NOT CONSISTENT GIVEN SETTINGS\n");
                printf("# ST-WSLDA-1D: nx=%d, ny=%d, nz=%d, dx=%f, dy=%f, dz=%f\n", _nx, _ny, _nz, _dx, _dy, _dz);
                printf("# SETTINGS   : nx=%d, ny=%d, nz=%d, dx=%f, dy=%f, dz=%f\n", NX, NY, NZ, DX, DY, DZ);
                ABORT_NOBARRIER;
            }
            printf("# ST-WSLDA-1D: file_name=`%s`\n",file_name);
            printf("# ST-WSLDA-1D: nwf (all modes)=%d\n",nwf);
            printf("# ST-WSLDA-1D: nx=%d, ny=%d, nz=%d, dx=%f, dy=%f, dz=%f\n", _nx, _ny, _nz, _dx, _dy, _dz);
            printf("# ST-WSLDA-1D: kF=%f, mu_a=%f, mu_b=%f, ec=%f, beta=%f\n", kF, mu[SPINA], mu[SPINB], ec, beta);
            fflush(stdout);
            
            // scan files and determine nwf in each of them
            nwf_s1dpca=nwf;
            file_operation( scan_stwslda1d_info_files(md.inprefix, 1, kvecs_to_consder, kvecs, &nwf_s1dpca, nwf_per_kyz) );
            printf("# ST-WSLDA-1D: nwf in binary files=%d\n",nwf_s1dpca);
            nwf=nwf_s1dpca;
        }
        MPI_Bcast( &nwf , 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
        MPI_Bcast( &kF , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ; eF=0.5*kF*kF;
        MPI_Bcast( mu , 2 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ;
        MPI_Bcast( &ec , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ;
        MPI_Bcast( &beta , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ;  it=0; t0=0.0;
        MPI_Bcast( nwf_per_kyz , kvecs_to_consder , MPI_INT , 0 , MPI_COMM_WORLD ) ; 
        MPI_Bcast( kvecs , sizeof(wslda_kmode)*kvecs_to_consder , MPI_BYTE , 0 , MPI_COMM_WORLD ) ; 
                
        // divide wf over processes
        if(ip==0) printf("# INIT1: nwf=%d wave-functions to scatter\n", nwf);
        if ( np > nwf )
        {
            if(ip==0) printf("# INIT1: np[%d] > nwf[%d]!\n", np, nwf);
            ABORT;
        }
        getnwfip( ip , np , nwf , &nwfip ) ;
        MPI_Gather( &nwfip , 1 , MPI_INT , wf_tbl , 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
        MPI_Bcast( wf_tbl , np , MPI_INT , 0 , MPI_COMM_WORLD ) ; 
        
        // my range of wf to manage
        int myuidx=0, mylidx=0;
        for(i=0; i<=ip; i++)
            myuidx+=wf_tbl[i];
        mylidx=myuidx-nwfip;
        
        // allocate memory for my wf
        cppmallocl(h_wavefun, NX*nwfip*2,double complex);
        cppmallocl(h_fbetaEn, nwfip,double);
        cppmallocl(h_kkyz,    nwfip*2,double);
        
        // read data
        int _4_max_readers = md.iogroups;
        int _4_nblocks = (int)ceil((float)(np)/_4_max_readers);
        for(i=0; i<_4_nblocks; i++)
        {
            if(ip==0) { printf("# INIT1: BLOCK ID[%d] CONSITING WITH %d PROCESSES READS DATA...\n", i, _4_max_readers); fflush(stdout);}
            if(ip%_4_nblocks == i) file_operation( read_stwslda1d_wf(md.inprefix, 1, kvecs_to_consder, kvecs, nwf_per_kyz, mylidx, myuidx, h_wavefun, h_fbetaEn, h_kkyz) );
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        for(i=0; i<nwfip; i++) h_fbetaEn[i]=fbeta(h_fbetaEn[i],beta); // convert quasiparticle energies into weights
        
        // load u and delta
        if(ip==0)
        {
            sprintf(file_name, "%s/s1dpca.pud", md.inprefix);
            printf("# INIT1: LOADING POTENTIALS `%s`...\n", file_name);
            file_operation( read_binary_file(file_name, NX*4*sizeof(double), 0, h_potentials) );           
        }
        
        MPI_Bcast(h_potentials, 4*NX, MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ;  
        
        // compute particle number and set Effg;
        double Ntota=0.0, Nmya=0.0;
        double Ntotb=0.0, Nmyb=0.0; 
        for(iwf=0; iwf<nwfip; iwf++)
        {
            int wcnt=get_weight_1d(h_kkyz[iwf], h_kkyz[iwf+nwfip]);
            ixyz=0;
            for ( ix = 0 ; ix < NX ; ix++ )
            {
                Nmya+=(pow(creal(h_wavefun[         iwf*NX+ixyz]),2)+pow(cimag(h_wavefun[         iwf*NX+ixyz]),2))*h_fbetaEn[iwf]*wcnt;
                Nmyb+=(pow(creal(h_wavefun[nwfip*NX+iwf*NX+ixyz]),2)+pow(cimag(h_wavefun[nwfip*NX+iwf*NX+ixyz]),2))*(1.0-h_fbetaEn[iwf])*wcnt;
                ixyz++;
            }  
        }
        MPI_Allreduce( &Nmya, &Ntota, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce( &Nmyb, &Ntotb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if(ip==0) printf("# INIT1: TOTAL NUMBER OF PARTICLES: SPIN_A=%16.8g SPIN_B=%16.8g TOTAL=%16.8g\n", Ntota, Ntotb, Ntota+Ntotb); 
        if(md.spinsymmetry==1) Ntota = Ntota+Ntotb; 
        Effg = 0.6 * Ntota * eF;
        fflush(stdout);
        
        // free memory
        free(kkx); free(kky); free(kkz);
        free(nwf_per_kyz); free(kvecs);
    }
    else
    {
        if(ip==0) printf("NOT SUPPORTED INITTYPE=%d!\n", md.inittype);
        ABORT;
    }
    
    // wait till loading is done
    MPI_Barrier(MPI_COMM_WORLD);

    // ====================================================================================
    // ======================== ALLOCATE GPU AND CPU BUFFERS ==============================
    // ====================================================================================    
    // Allocate memory for wave-functions and derivatives
    gpu_exec( gpu_malloc(NX*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_wf) );
    gpu_exec( gpu_malloc(NX*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_fkm1) );
    gpu_exec( gpu_malloc(NX*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_fkm2) );
    gpu_exec( gpu_malloc(NX*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_fkm3) );
#if INTEGRATION_SCHEME==AB4AM5
    gpu_exec( gpu_malloc(NX*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_fkm4) );
#endif
    gpu_exec( gpu_malloc(NX*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_wf_d_dx) );
    gpu_exec( gpu_malloc(NX*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_wf_laplace) );
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
    gpu_exec( gpu_malloc(NX*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_alphawf_laplace) );
#endif
    gpu_exec( gpu_malloc(nwfip*sizeof(double), (void **)&d_fbetaEn) );
    gpu_exec( gpu_malloc(nwfip*2*sizeof(double), (void **)&d_kkyz) ); // 2 accounts for ky and kz
    gpu_exec( host_malloc_pl(nwfip*sizeof(double), (void **)&h_qpe_nwfip) );
    cppmallocl(h_qpe_nwf, nwf,double);
    // aditional data needed for saving qpe
    MPI_Gather( &nwfip , 1 , MPI_INT , wf_tbl , 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
    MPI_Bcast( wf_tbl , np , MPI_INT , 0 , MPI_COMM_WORLD ) ;
    j=0;
    for(i=0; i<ip; i++) j+=wf_tbl[i];    
    MPI_Gather( &j , 1 , MPI_INT , wf_idx_tbl , 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
    MPI_Bcast( wf_idx_tbl , np , MPI_INT , 0 , MPI_COMM_WORLD ) ;
    
    
    // ====================================================================================
    // ==================================== COPY DATA TO GPU ==============================
    // ====================================================================================    
    if(md.inittype==5){ 

        if(ip==0) printf("# LOADING CHECKPOINT\n");
        b_t();
        size_t memsize;
#if INTEGRATION_SCHEME==AB3AM4
        cpu_exec( load_all (h_wavefun, MPI_COMM_WORLD, md.inprefix,
                  d_wf, d_fkm1, d_fkm2, d_fkm3, 
		  d_potentials, &t0, 
                  &nwf, &nwfip,
                  h_fbetaEn, h_kkyz, mu, &ec, &kF, &eF, &Effg,
		  HowMany) );
        memsize = (size_t)(nwf)*(NX)*2*4*16;
#elif INTEGRATION_SCHEME==AB4AM5
        cpu_exec( load_all_45 (h_wavefun, MPI_COMM_WORLD, md.inprefix,
                     d_wf, d_fkm1, d_fkm2, d_fkm3, d_fkm4,
                     d_potentials, &t0,
                     &nwf, &nwfip,
                     h_fbetaEn, h_kkyz, mu, &ec, &kF, &eF, &Effg,
                     HowMany) );
        memsize = (size_t)(nwf)*(NX)*2*5*16;
#else
        CHECK PCA_SETTINGS.H
#endif

        MPI_Barrier( MPI_COMM_WORLD ) ;
        rt = e_t(0);
        if(ip==0)
        {
            double memsize_gb = (double)(memsize) / pow(2,30);
            printf("# CHECKPOINT INFO: MODE=READ: DATA SIZE=%12.2f GB\n",  memsize_gb);
            printf("# CHECKPOINT INFO: MPI_NP_PER_IO_GROUP=%d.\n", HowMany);
            printf("# CHECKPOINT INFO: READ TIME=%12.2f sec\n", rt);
            printf("# CHECKPOINT INFO: READ SPEED=%12.3f GB/sec\n", memsize_gb/rt);
        }
    }
    
#ifdef TDWSLDA
    dt/=eF; // time step
#endif
    
    if(ip==0) printf("# INITIALIZING GPU BUFFERS OF ABM ALGORITHM...\n");
    
    if(md.inittype!=5)
    { 
        // copy wave-functions
        gpu_exec( memcopy_host2gpu(h_wavefun, d_wf,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) );   

        // we start from eigenstates - then fkm1, fkm2, fkm3 are zero
        for(i=0; i<NX*nwfip*2; i++) h_wavefun[i]=0.0 + I*0.0;
        gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm1,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm2,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm3,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
#if INTEGRATION_SCHEME==AB4AM5
        gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm4,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
#endif

        // copy potentials
        gpu_exec( memcopy_host2gpu(h_potentials, d_potentials,  (size_t)4*NX*sizeof(double)) );  
    }
    // copy weights
    gpu_exec( memcopy_host2gpu(h_fbetaEn, d_fbetaEn,  (size_t)nwfip  *sizeof(double)) ); 
    gpu_exec( memcopy_host2gpu(h_kkyz   , d_kkyz   ,  (size_t)nwfip*2*sizeof(double)) );
    
    // Set constants
    gpu_exec( memcopy_const(mu[SPINA], mu[SPINB], ec, t0, dt, kF) );    
    md.ec=ec; dc_ec=ec; dc_t0=t0; dc_np=np; dc_nwfip=nwfip;
    
    // ===================================================================================
    // ================================== EXTRA DATA =====================================
    // ===================================================================================
    if(ip==0) extra_data_size = get_extra_data_size(md.params);
    MPI_Bcast( &extra_data_size , sizeof(size_t) , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
    if(extra_data_size>0)
    {
        if(ip==0) printf("# EXTRA_DATA IS ACTIVE.\n");
        if(ip==0) printf("# ALLOCATING EXTRA_DATA OF SIZE %ld B.\n", extra_data_size); fflush(stdout);
        if ( ( extra_data = (void *) malloc( extra_data_size ) ) == NULL  )
        {                                                             
            fprintf( stderr , "error: cannot malloc()! Exiting!\n") ; 
            fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; 
            MPI_Finalize() ;
            /* Arrays will be cleared automatically */
            return( EXIT_FAILURE ) ; 
        }
        
        if(ip==0) cpu_exec( load_extra_data(extra_data_size, extra_data, md.params) );
        MPI_Bcast( extra_data , extra_data_size , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
        
        // copy extra data to GPU
        gpu_exec( gpu_malloc(extra_data_size, (void **)&d_extra_data) );
        gpu_exec( memcopy_host2gpu(extra_data, d_extra_data,  extra_data_size) ); 
        gpu_exec( memcopy_extra_data(extra_data_size, d_extra_data) );
    }
    
    // Process params and copy them to gpu;
#ifdef TDWSLDA
    process_params(md.params, kF, mu, extra_data_size, extra_data);
#else
    process_params(md.params, kF, mu);
#endif
    gpu_exec( memcopy_const_params(md.params) );
    
#ifdef BDG_MODE   
    gpu_exec( memcopy_const_BdG(md.aBdG) );
#endif
    
    if(ip==0) printf("# DONE.\n");
    
    
    // CUFFT plans
    size_t cufft_workSize;
    size_t wf_size = (size_t)2*NX*nwfip*sizeof(double complex);
    if(md.batch>2*nwfip) md.batch=2*nwfip;
 
    // Create plans
    gpu_exec( create_cufftPlans(md.batch, nwfip, &cufft_workSize) );
    
    
    // Allocate memory for plans
    if(ip==0) printf("# CUFFT[ip=%d]: cufft_workSize=%.2f times space of wf (%.2fMB)\n", ip, (double)cufft_workSize/(double)wf_size, (double)wf_size/pow(2.,20));
    if(workarea_size<cufft_workSize) workarea_size=cufft_workSize;
    gpu_exec( gpu_malloc(workarea_size, (void **)&d_workarea) );
    
    // Assign work space with plans
    gpu_exec( set_workspace_for_cufftPlan(d_workarea) );

    // ====================================================================================
    // ================================= INITIAL MEASUREMENT ==============================
    // ====================================================================================
    if(ip==0) printf("# INITIAL MEASUREMENT\n");
    // normalize wf 
    gpu_exec( normalize_wf(nwfip, d_wf, md.nthreads) );      
    // derivatives
    gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, NULL, NULL, d_wf_laplace, md.nthreads) ); 
    // densities - local reduction
    gpu_exec( calculate_densities(nwfip, d_wf, d_wf_d_dx, d_wf_laplace, d_kkyz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
    // densities - global reduction
    gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NX*sizeof(double)) ); 
    MPI_Allreduce( MPI_IN_PLACE, h_densities, 12*NX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)12*NX*sizeof(double)) ); 
    if(md.spinsymmetry>0) symmetrize_densities_device(d_densities); // special calse: spin-symmetric system
#ifndef TAU_COMPUTATION_VIA_GRADIENTS
    if(gradients_computed) density_caculate_tau(d_densities, md.nthreads);
#endif
    // potentials
//     if(md.inittype!=2) gpu_exec( compute_potentials(it, d_densities, d_potentials, cccoeff, md.nthreads) );
    // energy
    gpu_exec( compute_energy(it, d_densities, d_potentials, d_workarea, md.nthreads) ); 
    
    // get data for reporting
    // energy
    gpu_exec( memcopy_gpu2host(d_workarea, h_energy,  (size_t)TDWSLDAITEMS*sizeof(double)) );   
    // potentials
    gpu_exec( memcopy_gpu2host(d_potentials, h_potentials,  (size_t)4*NX*sizeof(double)) );     
    // densities - they are in h_densities
    double N_tot_init = h_energy[NPARTA]+h_energy[NPARTB]; // save initial value of particle number
    if(md.inittype!=5) Effg = 0.6 * N_tot_init * eF; // set correct value of Effg
    
    // report result
    if(ip==0)
    {
        time=t0+it*dt;
        energy_kin = h_energy[EKIN];
        energy_pot = h_energy[EPOT];
        energy_pair = h_energy[EPAIR];
        energy_current = h_energy[ECURRENT];
        energy_uext = h_energy[EPOTEXT];
        energy_dext = h_energy[EPAIREXT];
        energy_vext = h_energy[EVELEXT];
        energy_tot = 0.0;
        for(i=0;i<=EVELEXT;i++) energy_tot+=h_energy[i];
        Na=h_energy[NPARTA];
        Nb=h_energy[NPARTB];
        Laz=h_energy[LZA];
        Lbz=h_energy[LZB];
        
        printf("# GPU ENERGY     : ETOT=%12.8f, EKIN=%12.8f, EPOT=%12.8f, EPAIR=%12.8f, ECURRENT=%12.8f, EPOTEXT=%12.8f, EPAIREXT=%12.8f, EVELEXT=%12.8f\n", energy_tot/Effg, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_current/Effg, energy_uext/Effg, energy_dext/Effg, energy_vext/Effg); 
        
        // Create check stamp file
        sprintf(file_name, "%s_check.stamp", md.outprefix);
        printf("# CREATING CHECK STAMP FILE: `%s`\n",file_name);
        file_operation( touch_file(file_name) );
        // Take densities from device
        gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NX*sizeof(double)) );
        file_operation( check_stamp_entry_coeff(file_name, 12, NX, h_densities, TDWSLDAITEMS, h_energy, LY*LZ) ); 
        
        printf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %8s\n", "time*eF", "Na", "Nb", "Na+Nb", "ETOT", "EKIN", "EPOT", "EPAIR", "ECURRENT", "EPOTEXT", "EPAIREXT", "EVELEXT", "rt"); 
        printf("%12.4f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", time*eF, Na, Nb, Na+Nb, energy_tot/Effg, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_current/Effg, energy_uext/Effg, energy_dext/Effg, energy_vext/Effg);     
        
        // Create run log and add entry
        cpu_exec( logger_create_header(execcmd) );
        double _npart[2]={h_energy[NPARTA],h_energy[NPARTB]};
        cpu_exec( logger_add_entry(0, densall, potsall, kF, mu, h_energy, _npart, md.params, extra_data_size, extra_data) );
    }    
    
    cpu_exec( wslda_check_array_against_naninf(TDWSLDAITEMS, h_energy) );

    // Create binary files and add initial measurement
    wdata_metadata wdmd; 
    file_operation( create_wdata_metadata(&md, 1, t0, md.timesteps*dt, md.spinsymmetry, &wdmd) );
    
    // set constants
    wdata_setconst(&wdmd, "kF", kF);
    wdata_setconst(&wdmd, "eF", eF);
    wdata_setconst(&wdmd, "mu_a", mu[SPINA]);
    wdata_setconst(&wdmd, "mu_b", mu[SPINB]);
    
    // prepare database
    MPI_Barrier(MPI_COMM_WORLD);
    if(ip==0) 
    {
        file_operation( clear_files(&md, &wdmd) );
        file_operation( write_wdata_metadata_file(&md, &wdmd, "td-wslda-1d") );
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NX*sizeof(double)) );
    gpu_exec( memcopy_gpu2host(d_potentials, h_potentials,  (size_t)4*NX*sizeof(double)) );
    set_ptr_d_delta(d_potentials+2*NXY);
    file_operation( write_measurments(&wdmd, MPI_COMM_WORLD, "td", it, densall, potsall) );
    if(ip==0) file_operation( write_wdata_metadata_file(&md, &wdmd, "td-wslda-1d") );
    if(ip==0)
    {  
#ifdef STORE_QPE
        sprintf(file_name, "%s_qpe.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, 1, 1, 1.0, 1.0, 1.0, eF, t0, md.timesteps*dt) ); 
        // save zeros for qpe for initial measurement - to avoid expensive computation of qpe
        for(i=0; i<nwf; i++) h_qpe_nwf[i]=0.0;
        sprintf(file_name, "%s_qpe.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, h_qpe_nwf, sizeof(double)*nwf) );  
#endif
    }
    
    // ====================================================================================
    // ================================ REAL TIME EVOLUTION  ==============================
    // ====================================================================================
    int i_meas, i_step;
    
    if(md.inittype!=5 && md.selfstart==1)
    { 
        // NOTE: I assume that potential is constant during first steps
        // NOTE: I assume there is no quantum friction during the first steps
        if(ip==0) printf("# SELFSTART: EXECUTING TAYLOR EXPANSION OF THE EVOLUTION OPERATOR.\n");
        
        double *d_qpe; // buffers for quasiparticle energies
        gpu_exec( gpu_malloc(nwfip*sizeof(double), (void **)&d_qpe) );
              
        int selfstart_steps, exp_iters;
#if INTEGRATION_SCHEME==AB3AM4
        selfstart_steps=3;
        exp_iters = 4;
#elif INTEGRATION_SCHEME==AB4AM5
        selfstart_steps=4;
        exp_iters = 5;
#endif  
        double complex *h_fkm;
        cppmallocl(h_fkm, NX*nwfip*2*selfstart_steps,double complex);

        for(i_step=0; i_step<selfstart_steps; i_step++)
        {
            // ----------------------------- predictor -----------------------------------
            // compute value of quantum friction coefficient
            qfalpha = 0.0; // NOTE: I assume there is no quantum friction during the first steps
            cccoeff = h_smooth_step(t0+it*dt, md.ccstart/eF,  md.ccstop/eF,  md.ccswitch/eF, 1.0);
            
            // derivatives
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
            if(qfalpha>0.0) gradients_computed=1; else gradients_computed=0;
#endif
            if(gradients_computed)
            {
                gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, NULL, NULL, d_wf_laplace, md.nthreads) ); 
            }
            else
            {
                gpu_exec( compute_laplace(2*nwfip, d_wf, d_wf_laplace, md.nthreads) );
            }
            // densities - local reduction
            gpu_exec( calculate_densities(nwfip, d_wf, d_wf_d_dx, d_wf_laplace, d_kkyz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
            // densities - global reduction
            gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NX*sizeof(double)) ); 
            MPI_Allreduce( MPI_IN_PLACE, h_densities, 12*NX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)12*NX*sizeof(double)) );
            if(md.spinsymmetry>0) symmetrize_densities_device(d_densities); // special calse: spin-symmetric system
#ifndef TAU_COMPUTATION_VIA_GRADIENTS
            if(gradients_computed) density_caculate_tau(d_densities, md.nthreads);
#endif
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // effective mass correction
            gpu_exec( multiply_wf_by_alpha(nwfip, d_wf, d_alphawf_laplace, d_densities, md.nthreads) );
            gpu_exec( compute_laplace(2*nwfip, d_alphawf_laplace, d_alphawf_laplace, md.nthreads) );
#endif
            // potentials - NOTE: it=0!
            gpu_exec( compute_potentials(0, d_densities, d_potentials, cccoeff, md.nthreads) );
            
            // executing exp[-i*H(t)*dt]*psi
            // H*psi - first execution, d_fkm3 as working buffer
            gpu_exec( memcopy_gpu2gpu(d_wf, d_fkm3, (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) );
            gpu_exec( apply_hamiltonian(0, nwfip, d_fkm3, d_fkm1, /* NOTE - d_fkm1 as output buffer  */
                                    d_wf_d_dx, d_kkyz, d_wf_laplace, d_alphawf_laplace,
                                    d_densities, d_potentials, qfalpha, NULL, cccoeff, 
                                    md.nthreads) );
            // Make copy of qpe
            gpu_exec( memcopy_gpu2gpu(d_workarea, d_qpe, (size_t)nwfip*sizeof(double)) );
            
            // Store H*Psi
            gpu_exec( memcopy_gpu2host(d_fkm1, h_fkm+i_step*2*nwfip*NX,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
                        
            // Add contribution from Taylor expansion
            gpu_exec( taylor_expansion_contribution(1, 0.5*dt, nwfip, d_fkm1, d_fkm3, d_fkm2, md.nthreads) );
                            
            // H*psi - remaining executions
            for(i_meas=1; i_meas<exp_iters; i_meas++)
            {
                // derivatives
                if(gradients_computed)
                {
                    gpu_exec( compute_derivatives(2*nwfip, d_fkm2, d_wf_d_dx, NULL, NULL, d_wf_laplace, md.nthreads) ); 
                }
                else
                {
                    gpu_exec( compute_laplace(2*nwfip, d_fkm2, d_wf_laplace, md.nthreads) );
                }
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
                // effective mass correction
                gpu_exec( multiply_wf_by_alpha(nwfip, d_fkm2, d_alphawf_laplace, d_densities, md.nthreads) );
                gpu_exec( compute_laplace(2*nwfip, d_alphawf_laplace, d_alphawf_laplace, md.nthreads) );
#endif
            
                // H*psi
                gpu_exec( apply_hamiltonian(0, nwfip, d_fkm2, d_fkm1,
                                        d_wf_d_dx, d_kkyz, d_wf_laplace, d_alphawf_laplace,
                                        d_densities, d_potentials, qfalpha, d_qpe, cccoeff, 
                                        md.nthreads) );
                // Add contribution from Taylor expansion
                gpu_exec( taylor_expansion_contribution(i_meas+1, 0.5*dt, nwfip, d_fkm1, d_fkm3, d_fkm2, md.nthreads) );
            }
            
            // normalize wf 
            gpu_exec( normalize_wf(nwfip, d_fkm3, md.nthreads) );
            
            // NOTE: d_fkm3 keeps prediction of wave-function for midpoint
            
            // ----------------------------- corrector -----------------------------------
            // compute value of quantum friction coefficient
            qfalpha = 0.0; // NOTE: I assume there is no quantum friction during the first steps
            cccoeff = h_smooth_step(t0+it*dt, md.ccstart/eF,  md.ccstop/eF,  md.ccswitch/eF, 1.0);

            // derivatives
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
            if(qfalpha>0.0) gradients_computed=1; else gradients_computed=0;
#endif
            if(gradients_computed)
            {
                gpu_exec( compute_derivatives(2*nwfip, d_fkm3, d_wf_d_dx, NULL, NULL, d_wf_laplace, md.nthreads) ); 
            }
            else
            {
                gpu_exec( compute_laplace(2*nwfip, d_fkm3, d_wf_laplace, md.nthreads) );
            }
            // densities - local reduction
            gpu_exec( calculate_densities(nwfip, d_fkm3, d_wf_d_dx, d_wf_laplace, d_kkyz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
            // densities - global reduction
            gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NX*sizeof(double)) ); 
            MPI_Allreduce( MPI_IN_PLACE, h_densities, 12*NX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)12*NX*sizeof(double)) );
            if(md.spinsymmetry>0) symmetrize_densities_device(d_densities); // special calse: spin-symmetric system
#ifndef TAU_COMPUTATION_VIA_GRADIENTS
            if(gradients_computed) density_caculate_tau(d_densities, md.nthreads);
#endif
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // effective mass correction
            gpu_exec( multiply_wf_by_alpha(nwfip, d_fkm3, d_alphawf_laplace, d_densities, md.nthreads) );
            gpu_exec( compute_laplace(2*nwfip, d_alphawf_laplace, d_alphawf_laplace, md.nthreads) );
#endif
            // potentials - NOTE: it=0!
            gpu_exec( compute_potentials(0, d_densities, d_potentials, cccoeff, md.nthreads) );
            // NOTE - densities and potentials are computed for midpoint 
            
            // recompute derivatives for d_wf
            if(gradients_computed)
            {
                gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, NULL, NULL, d_wf_laplace, md.nthreads) ); 
            }
            else
            {
                gpu_exec( compute_laplace(2*nwfip, d_wf, d_wf_laplace, md.nthreads) );
            }            
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // effective mass correction
            gpu_exec( multiply_wf_by_alpha(nwfip, d_wf, d_alphawf_laplace, d_densities, md.nthreads) );
            gpu_exec( compute_laplace(2*nwfip, d_alphawf_laplace, d_alphawf_laplace, md.nthreads) );
#endif
           
            // executing exp[-i*H(t+dt/2)*dt]*psi
            // H*psi - first execution
            gpu_exec( apply_hamiltonian(0, nwfip, d_wf, d_fkm1, /* NOTE - d_fkm1 as output buffer  */
                                    d_wf_d_dx, d_kkyz, d_wf_laplace, d_alphawf_laplace,
                                    d_densities, d_potentials, qfalpha, NULL, cccoeff, 
                                    md.nthreads) );
            // Make copy of qpe
            gpu_exec( memcopy_gpu2gpu(d_workarea, d_qpe, (size_t)nwfip*sizeof(double)) );
            
            // Add contribution from Taylor expansion
            gpu_exec( taylor_expansion_contribution(1, dt, nwfip, d_fkm1, d_wf, d_fkm2, md.nthreads) );
                
            // H*psi - remaining executions
            for(i_meas=1; i_meas<exp_iters; i_meas++)
            {
                // derivatives
                if(gradients_computed)
                {
                    gpu_exec( compute_derivatives(2*nwfip, d_fkm2, d_wf_d_dx, NULL, NULL, d_wf_laplace, md.nthreads) ); 
                }
                else
                {
                    gpu_exec( compute_laplace(2*nwfip, d_fkm2, d_wf_laplace, md.nthreads) );
                }
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
                // effective mass correction
                gpu_exec( multiply_wf_by_alpha(nwfip, d_fkm2, d_alphawf_laplace, d_densities, md.nthreads) );
                gpu_exec( compute_laplace(2*nwfip, d_alphawf_laplace, d_alphawf_laplace, md.nthreads) );
#endif                
           
                // H*psi
                gpu_exec( apply_hamiltonian(0, nwfip, d_fkm2, d_fkm1,
                                        d_wf_d_dx, d_kkyz, d_wf_laplace, d_alphawf_laplace,
                                        d_densities, d_potentials, qfalpha, d_qpe, cccoeff, 
                                        md.nthreads) );
                // Add contribution from Taylor expansion
                gpu_exec( taylor_expansion_contribution(i_meas+1, dt, nwfip, d_fkm1, d_wf, d_fkm2, md.nthreads) );
            }
            
            // normalize wf 
            gpu_exec( normalize_wf(nwfip, d_wf, md.nthreads) );
            
            // NOTE: d_wf keeps wave-function for t+dt
            if(ip==0) { printf("# SELFSTART: i_step=%d\n", i_step); fflush(stdout); }
        }
        
        // Copy fkm1, ..., fkm4 back to gpu
#if INTEGRATION_SCHEME==AB3AM4
        gpu_exec( memcopy_host2gpu(h_fkm+2*2*nwfip*NX, d_fkm1,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_fkm+1*2*nwfip*NX, d_fkm2,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_fkm+0*2*nwfip*NX, d_fkm3,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
#elif INTEGRATION_SCHEME==AB4AM5
        gpu_exec( memcopy_host2gpu(h_fkm+3*2*nwfip*NX, d_fkm1,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_fkm+2*2*nwfip*NX, d_fkm2,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_fkm+1*2*nwfip*NX, d_fkm3,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_fkm+0*2*nwfip*NX, d_fkm4,  (size_t)2*nwfip*NX*sizeof(cufftDoubleComplex)) ); 
#endif        
        
        // clear memory
        gpu_exec( gpu_free(d_qpe) );
        free(h_fkm);
        if(ip==0) printf("# SELFSTART: DONE.\n");
         
        // derivatives
        gradients_computed=1;
        gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, NULL, NULL, d_wf_laplace, md.nthreads) ); 
        // densities - local reduction
        gpu_exec( calculate_densities(nwfip, d_wf, d_wf_d_dx, d_wf_laplace, d_kkyz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
        // densities - global reduction
        gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NX*sizeof(double)) ); 
        MPI_Allreduce( MPI_IN_PLACE, h_densities, 12*NX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)12*NX*sizeof(double)) ); 
        if(md.spinsymmetry>0) symmetrize_densities_device(d_densities); // special calse: spin-symmetric system
#ifndef TAU_COMPUTATION_VIA_GRADIENTS
        if(gradients_computed) density_caculate_tau(d_densities, md.nthreads);
#endif
        // potentials
        gpu_exec( compute_potentials(it, d_densities, d_potentials, cccoeff, md.nthreads) );
        // energy
        gpu_exec( compute_energy(it, d_densities, d_potentials, d_workarea, md.nthreads) );
    
        // get data for reporting
        // energy
        gpu_exec( memcopy_gpu2host(d_workarea, h_energy,  (size_t)TDWSLDAITEMS*sizeof(double)) );   
        // potentials
        gpu_exec( memcopy_gpu2host(d_potentials, h_potentials,  (size_t)4*NX*sizeof(double)) );     
        // densities - they are in h_densities
    
        // report result
        if(ip==0)
        {
            time=t0+it*dt;
            energy_kin = h_energy[EKIN];
            energy_pot = h_energy[EPOT];
            energy_pair = h_energy[EPAIR];
            energy_current = h_energy[ECURRENT];
            energy_uext = h_energy[EPOTEXT];
            energy_dext = h_energy[EPAIREXT];
            energy_vext = h_energy[EVELEXT];
            energy_tot = 0.0;
            for(i=0;i<=EVELEXT;i++) energy_tot+=h_energy[i];
            Na=h_energy[NPARTA];
            Nb=h_energy[NPARTB];
            Laz=h_energy[LZA];
            Lbz=h_energy[LZB];
            
            printf("# AFTER SELFSTART: ETOT=%12.8f, EKIN=%12.8f, EPOT=%12.8f, EPAIR=%12.8f, ECURRENT=%12.8f, EPOTEXT=%12.8f, EPAIREXT=%12.8f, EVELEXT=%12.8f\n", energy_tot/Effg, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_current/Effg, energy_uext/Effg, energy_dext/Effg, energy_vext/Effg);  
             
            printf("%12.4f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", time*eF, Na, Nb, Na+Nb, energy_tot/Effg, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_current/Effg, energy_uext/Effg, energy_dext/Effg, energy_vext/Effg);     
        }    
    
    }
    
    // ====================================================================================
    // ========================================= ABM  =====================================
    // ====================================================================================
    for (i_meas=0; i_meas<md.measurements; i_meas++)
    {
        b_t(); // reset timer
        
        for(i_step=0; i_step<md.timesteps; i_step++)
        {
            // ----------------------------- predictor -----------------------------------
            // compute value of quantum friction coefficient and current corrections coeff
            qfalpha = md.qfalpha*h_smooth_step(t0+(it+1)*dt, md.qfstart/eF,  md.qfstop/eF,  md.qfswitch/eF, 1.0);
            cccoeff = h_smooth_step(t0+(it+1)*dt, md.ccstart/eF,  md.ccstop/eF,  md.ccswitch/eF, 1.0);

#if INTEGRATION_SCHEME==AB3AM4
            gpu_exec( amb_step1(nwfip, d_wf, d_fkm1, d_fkm2, d_fkm3, md.nthreads) );
#elif INTEGRATION_SCHEME==AB4AM5
            gpu_exec( amb45_step1(nwfip, d_wf, d_fkm1, d_fkm2, d_fkm3, d_fkm4, md.nthreads) );
#else
            CHECK PCA_SETTINGS.H
#endif
            // normalize wf 
            gpu_exec( normalize_wf(nwfip, d_wf, md.nthreads) );
            // derivatives
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
            if(qfalpha>0.0) gradients_computed=1; else gradients_computed=0;
#endif
            if(gradients_computed)
            {
                gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, NULL, NULL, d_wf_laplace, md.nthreads) ); 
            }
            else
            {
                gpu_exec( compute_laplace(2*nwfip, d_wf, d_wf_laplace, md.nthreads) );
            }
            // densities - local reduction
            gpu_exec( calculate_densities(nwfip, d_wf, d_wf_d_dx, d_wf_laplace, d_kkyz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
            // densities - global reduction
            mpipackagesize = EXCHANGE_SIZE;
            gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)mpipackagesize*NX*sizeof(double)) ); 
            MPI_Allreduce( MPI_IN_PLACE, h_densities, mpipackagesize*NX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)mpipackagesize*NX*sizeof(double)) );
            if(md.spinsymmetry>0) symmetrize_densities_device(d_densities); // special calse: spin-symmetric system
#ifndef TAU_COMPUTATION_VIA_GRADIENTS
            if(gradients_computed) density_caculate_tau(d_densities, md.nthreads);
#endif
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // effective mass correction
            gpu_exec( multiply_wf_by_alpha(nwfip, d_wf, d_alphawf_laplace, d_densities, md.nthreads) );
            gpu_exec( compute_laplace(2*nwfip, d_alphawf_laplace, d_alphawf_laplace, md.nthreads) );
#endif
            // potentials
            gpu_exec( compute_potentials(it+1, d_densities, d_potentials, cccoeff, md.nthreads) );
            // H*psi
            gpu_exec( apply_hamiltonian(it, nwfip, d_wf, d_wf_laplace, /* NOTE - d_wf_laplace as output buffer  */
                                    d_wf_d_dx, d_kkyz, d_wf_laplace, d_alphawf_laplace,
                                    d_densities, d_potentials, qfalpha, NULL, cccoeff, 
                                    md.nthreads) ); 
            
            // ----------------------------- corrector -----------------------------------
            // compute value of quantum friction coefficient  and current corrections coeff
            qfalpha = md.qfalpha*h_smooth_step(t0+(it+1)*dt, md.qfstart/eF,  md.qfstop/eF,  md.qfswitch/eF, 1.0);
            cccoeff = h_smooth_step(t0+(it+1)*dt, md.ccstart/eF,  md.ccstop/eF,  md.ccswitch/eF, 1.0);
#if INTEGRATION_SCHEME==AB3AM4
            gpu_exec( amb_step4(nwfip, d_wf_laplace, /* NOTE - d_wf_laplace as itermiediate buffer  */
                         d_wf, d_fkm3, md.nthreads) );  
#elif INTEGRATION_SCHEME==AB4AM5
            gpu_exec( amb45_step4(nwfip, d_wf_laplace, /* NOTE - d_wf_laplace as itermiediate buffer  */
                         d_wf, d_fkm4, md.nthreads) );
#else
            CHECK PCA_SETTINGS.H
#endif
            // normalize wf 
            gpu_exec( normalize_wf(nwfip, d_wf, md.nthreads) );
            // derivatives
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
            if(qfalpha>0.0 || i_step==md.timesteps-1) gradients_computed=1; else gradients_computed=0;
#endif
            if(gradients_computed)
            {
                gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, NULL, NULL, d_wf_laplace, md.nthreads) ); 
            }
            else
            {
                gpu_exec( compute_laplace(2*nwfip, d_wf, d_wf_laplace, md.nthreads) );
            }
            // densities - local reduction
            gpu_exec( calculate_densities(nwfip, d_wf, d_wf_d_dx, d_wf_laplace, d_kkyz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
            // densities - global reduction
            if(i_step==md.timesteps-1) mpipackagesize = 12; else mpipackagesize = EXCHANGE_SIZE;
            gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)mpipackagesize*NX*sizeof(double)) ); 
            MPI_Allreduce( MPI_IN_PLACE, h_densities, mpipackagesize*NX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)mpipackagesize*NX*sizeof(double)) ); 
            if(md.spinsymmetry>0) symmetrize_densities_device(d_densities); // special calse: spin-symmetric system
#ifndef TAU_COMPUTATION_VIA_GRADIENTS
            if(gradients_computed) density_caculate_tau(d_densities, md.nthreads);
#endif
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // effective mass correction
            gpu_exec( multiply_wf_by_alpha(nwfip, d_wf, d_alphawf_laplace, d_densities, md.nthreads) );
            gpu_exec( compute_laplace(2*nwfip, d_alphawf_laplace, d_alphawf_laplace, md.nthreads) );
#endif
            // potentials
            gpu_exec( compute_potentials(it+1, d_densities, d_potentials, cccoeff, md.nthreads) );  
            
            // Preparation of buffers for next step
#if INTEGRATION_SCHEME==AB3AM4
            d_tmp_ptr=d_fkm3;
            d_fkm3=d_fkm2;
            d_fkm2=d_fkm1;
            d_fkm1=d_tmp_ptr;
#elif INTEGRATION_SCHEME==AB4AM5            
            d_tmp_ptr=d_fkm4;
            d_fkm4=d_fkm3;
            d_fkm3=d_fkm2;
            d_fkm2=d_fkm1;
            d_fkm1=d_tmp_ptr;
#else
            CHECK PCA_SETTINGS.H
#endif
            // H*psi
            gpu_exec( apply_hamiltonian(it, nwfip, d_wf, d_fkm1, 
                                    d_wf_d_dx, d_kkyz, d_wf_laplace, d_alphawf_laplace, 
                                    d_densities, d_potentials, qfalpha, NULL, cccoeff, 
                                    md.nthreads) );            
            
            it++; // update global time counter
        }
        
        // ----------------------------- measurement -------------------------------------
#ifdef STORE_QPE
        // Save quasiparticle energies - computation of energy will destroy them
        gpu_exec( memcopy_gpu2host(d_workarea, h_qpe_nwfip,  (size_t)nwfip*sizeof(double)) );
        MPI_Gatherv(h_qpe_nwfip,nwfip,MPI_DOUBLE,h_qpe_nwf,wf_tbl,wf_idx_tbl,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if(ip==0)
        {
            sprintf(file_name, "%s_qpe.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, h_qpe_nwf, sizeof(double)*nwf) ); 
        }
#endif
        
        // energy
        gpu_exec( compute_energy(it, d_densities, d_potentials, d_workarea, md.nthreads) );
        
        // get data for reporting
        // energy
        gpu_exec( memcopy_gpu2host(d_workarea, h_energy,  (size_t)TDWSLDAITEMS*sizeof(double)) );   
        // potentials
        gpu_exec( memcopy_gpu2host(d_potentials, h_potentials,  (size_t)4*NX*sizeof(double)) );     
        // densities - they are in h_densities
        
        rt=e_t(0); // get timing
        
        // report result
        if(ip==0)
        {
            time=t0+it*dt;
            energy_kin = h_energy[EKIN];
            energy_pot = h_energy[EPOT];
            energy_pair = h_energy[EPAIR];
            energy_current = h_energy[ECURRENT];
            energy_uext = h_energy[EPOTEXT];
            energy_dext = h_energy[EPAIREXT];
            energy_vext = h_energy[EVELEXT];
            energy_tot = 0.0;
            for(i=0;i<=EVELEXT;i++) energy_tot+=h_energy[i];
            Na=h_energy[NPARTA];
            Nb=h_energy[NPARTB];
            Laz=h_energy[LZA];
            Lbz=h_energy[LZB];
            
            printf("%12.4f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %8.2f\n", time*eF, Na, Nb, Na+Nb, energy_tot/Effg, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_current/Effg, energy_uext/Effg, energy_dext/Effg, energy_vext/Effg, rt);        
            
            double _npart[2]={h_energy[NPARTA],h_energy[NPARTB]};
            cpu_exec( logger_add_entry(i_meas+1, densall, potsall, kF, mu, h_energy, _npart, md.params, extra_data_size, extra_data) );
        } 
        
        cpu_exec( wslda_check_array_against_naninf(TDWSLDAITEMS, h_energy) );
        
        // add binary data
        gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NX*sizeof(double)) );
        gpu_exec( memcopy_gpu2host(d_potentials, h_potentials,  (size_t)4*NX*sizeof(double)) );
        file_operation( write_measurments(&wdmd, MPI_COMM_WORLD, "td", it, densall, potsall) );
        if(ip==0) file_operation( write_wdata_metadata_file(&md, &wdmd, "td-wslda-1d") );
        
        if(ip==0) 
        {
            sprintf(file_name, "%s_checkpoint.make", md.outprefix);
            if(exists(file_name)) forceCP = 1;
            
            // another possibility - we approach walltime
            rt = wt_e_t();
//             printf("# TIME TO WALLTIME: %8.2f [h]\n", md.walltime-rt/3600.);
            if(rt>md.walltime*3600) 
            {
                forceCP = 1;
                printf("# WALLTIME REACHED!\n");
            }
        }
        MPI_Bcast(&forceCP, 1, MPI_INT , 0 , MPI_COMM_WORLD );
        if(forceCP==1) md.checkpoint = 1;
        if(forceCP==1 && ip==0) printf("# CONDUCTING EMERGENCY CHECKPOINT!\n");
        if(forceCP==1) break;
        
        // Check if siulation is stable
        if( fabs( (h_energy[NPARTA]+h_energy[NPARTB]-N_tot_init)/N_tot_init )>N_STABILITY_CRITERIA )
        {
            if(ip==0) printf("# SIMULATION INSTABILITY CRITERIA MET!!! BREAKING!!!\n");
            break;
        }

        fflush(stdout); // clear output
    }
    time=t0+it*dt;
    //check point
    if (md.checkpoint)
    {
        b_t(); // start measureing time of writing
        size_t memsize;
#if INTEGRATION_SCHEME==AB3AM4
        save_all(h_wavefun, MPI_COMM_WORLD, md.outprefix,
                    d_wf, d_fkm1, d_fkm2, d_fkm3,
                    d_potentials, &time, 
            nwf, nwfip,
            h_fbetaEn, h_kkyz, mu, &ec, &kF, & eF, &Effg, 
            HowMany);
        memsize = (size_t)(nwf)*(NX)*2*4*16;
#elif INTEGRATION_SCHEME==AB4AM5            
        save_all_45(h_wavefun, MPI_COMM_WORLD, md.outprefix,
                    d_wf, d_fkm1, d_fkm2, d_fkm3, d_fkm4,
                    d_potentials, &time,
                    nwf, nwfip,
                    h_fbetaEn, h_kkyz, mu, &ec, &kF, & eF, &Effg, 
                    HowMany);
        memsize = (size_t)(nwf)*(NX)*2*5*16;
#else
            CHECK PCA_SETTINGS.H
#endif
        MPI_Barrier( MPI_COMM_WORLD ) ;
        rt = e_t(0);
        if(ip==0)
        {
            double memsize_gb = (double)(memsize) / pow(2,30);
            printf("# CHECKPOINT INFO: MODE=WRITE: DATA SIZE=%12.2f GB\n",  memsize_gb);
            printf("# CHECKPOINT INFO: MPI_NP_PER_IO_GROUP=%d.\n", HowMany);
            printf("# CHECKPOINT INFO: WRITE TIME=%12.2f sec\n", rt);
            printf("# CHECKPOINT INFO: WRITE SPEED=%12.3f GB/sec\n", memsize_gb/rt);
            sprintf(file_name, "%s_check.stamp", md.outprefix);
            printf("# CREATING CHECK STAMP: `%s`\n",file_name);
            // Take densities from device
            gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NX*sizeof(double)) );
            file_operation( check_stamp_entry_coeff(file_name, 12, NX, h_densities, TDWSLDAITEMS, h_energy, LY*LZ) );   
        }
    }
    /* messy exit here */
    MPI_Barrier( MPI_COMM_WORLD ) ;
    MPI_Finalize() ;
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
