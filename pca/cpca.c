// This is code for solving DFT equations for polarized cold atoms (pca)
// Version constrained polarized cold atoms (cpca)
// Along z-direction translational symmetry is assumeed.

// Authors:
// Gabriel Wlazlowski <gabrielw@if.pw.edu.pl>

// compile:
// mpicc cpca.c -o cpca -lm

// TODO:
// 1. SPINSYMMETRY_MODE

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
#include "cpca_derivative.h"
#include "cpca_kernels.h"
#include "cpca_dens.h"
#include "pca_utils.h"
#include "pca_io.h"
#include "pca_uniform.h"
#include "pca_logger.h"
#include "cpca_checkpoint.h"

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
    int gradients_computed = 1; // flag indicating if code uses gradients in computaation, by default equal 1
    
    // arrays
    double complex *h_wavefun; // pointer to wave-functions on host (cpu) side 
    double *h_fbetaEn; // pointer to weights of wave-functions on host (cpu) side
    double *d_fbetaEn; // pointer to weights of wave-functions on device (gpu) side
    double *h_kkz; // pointer to kz wave-vector values on host (cpu) side
    double *d_kkz; // pointer to kz wave-vector values on device (gpu) side
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
    cufftDoubleComplex *d_wf_d_dy; // pointer to derivative of wave-function d/dy (GPU)
    cufftDoubleComplex *d_wf_laplace; // pointer to laplace of wave-function (d^2/dx^2 + d^2/dy^2 + d^2/dz^2) (GPU)
    cufftDoubleComplex *d_alphawf_laplace=NULL; // pointer to laplace of alpha*wave-function (d^2/dx^2 + d^2/dy^2 + d^2/dz^2) (GPU)
    cufftDoubleComplex *d_tmp_ptr;
    double *h_qpe_nwfip, *h_qpe_nwf; // buffers for quasiparticle energies

    // other technical variables
    int *wf_tbl, *wf_idx_tbl; // table of size np, keeps number of managed wf by each process
    size_t  workarea_size=(size_t)5*NXY*sizeof(double); // minimal size of workarea 
    
    // for reporting
    double eF_a, eF_b, eF, Effg;
    double time;
    double energy_kin, energy_pot, energy_pair, energy_CM, energy_uext, energy_tot;   
    double Na, Nb;
    double Laz, Lbz; // angular momentum
    char file_name[256];
    
    // timing
    double rt;    
    
    // quantum friction
    double qfalpha=0.0;
    
    // current corrections
    double cccoeff=0.0;
    
//     // TODO - ad hoc section
//     double *h_weight_andreev, *d_weight_andreev;
//     double *h_densities_andreev; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
//     double *d_densities_andreev; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (GPU)
//     gpu_exec( host_malloc_pl((size_t)12*NXY*sizeof(double), (void **)&h_densities_andreev) );
//     gpu_exec(     gpu_malloc((size_t)12*NXY*sizeof(double), (void **)&d_densities_andreev) );
//     // TODO - end of section
    
    /* start main */
    wt_b_t(); // tag init time
    
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
    
    // variables
    dt= md.dt ;
    double Emax =  M_PI*M_PI/2.; // E_max = p_max^2 / 2m, where: p_max is maximum momentum on the lattice, p_max=M_PI (if lattice spacing is 1.0)
    dt/=Emax; // time step
    
#ifdef UNIFORM_TEST_MODE
    md.Na = ceil(1.0/(6.0*M_PI*M_PI) * NXYZ);
#ifdef SPINSYMMETRY_MODE
    md.Nb = md.Na;
#else
    md.Nb = md.Na+1;
#endif
    if(ip==0) printf("# UNIFORM_TEST_MODE: SETTING NUMBER OF PARTICLES Na=%f\n", md.Na);
#endif
    
    // ====================================================================================
    // ============================= INITIALIZE GPU =======================================
    // ====================================================================================    
#if TARGET_MACHINE == TITAN
    if(ip==0) printf("# MACHINE: TITAN\n");
    
    int deviceId=0;
    printf("# Process ip=%d uses deviceId=%d\n", ip, deviceId);
    gpu_exec( set_gpu(deviceId) );
#endif
    
#if TARGET_MACHINE == TSUBAME
    if(ip==0) printf("# MACHINE: TSUBAME\n");
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    int deviceId=ip / md.tsubamenodes;
    printf("# Process ip=%d on node %s uses deviceId=%d\n", ip, processor_name, deviceId);
    gpu_exec( set_gpu(deviceId) );
#endif
    
#if TARGET_MACHINE == DWARF
    if(ip==0) printf("# MACHINE: DWARF\n");
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    int *ompi_local_rank;
    ompi_local_rank = (int *)malloc(sizeof(int)*np);
    int ompi_ppn=4;
    if(strcmp (processor_name,"node2061.grid4cern.if.pw.edu.pl")==0) ompi_ppn=8;
    if(strcmp (processor_name,"node2062.grid4cern.if.pw.edu.pl")==0) ompi_ppn=8;
    MPI_Allgather(&ompi_ppn,1,MPI_INT,ompi_local_rank,1,MPI_INT,MPI_COMM_WORLD);
    int ompi_i=0, ompi_j;
    while(ompi_i<np)
    {
        if(ompi_local_rank[ompi_i]==8)
        {
            for(ompi_j=0; ompi_j<8; ompi_j++) ompi_local_rank[ompi_i+ompi_j]=ompi_j;
            ompi_i+=8;
        }
        else
        {
            for(ompi_j=0; ompi_j<4; ompi_j++) ompi_local_rank[ompi_i+ompi_j]=ompi_j;
            ompi_i+=4;
        }
    }

    int deviceId=ompi_local_rank[ip];
    printf("ip=%d running on machine `%s` uses device number %d\n",ip, processor_name, deviceId);
    free(ompi_local_rank);

    // now you can use deviceId for set device
    gpu_exec( set_gpu(deviceId) );
#endif
    
    if(ip==0) printf("# LATTICE: %d x %d x %d\n", NX, NY, NZ);
    
    // ====================================================================================
    // ======================== ALLOCATE GPU AND CPU BUFFERS ==============================
    // ====================================================================================
    gpu_exec( host_malloc_pl((size_t)12*NXY*sizeof(double), (void **)&h_densities) );
    gpu_exec(     gpu_malloc((size_t)12*NXY*sizeof(double), (void **)&d_densities) );
    
    // potentials
    gpu_exec( host_malloc_pl((size_t)4*NXY*sizeof(double), (void **)&h_potentials) );
    gpu_exec(     gpu_malloc((size_t)4*NXY*sizeof(double), (void **)&d_potentials) );   
    
    // energy
    gpu_exec( host_malloc_pl((size_t)9*sizeof(double), (void **)&h_energy) );
    
    // For easier access to data
    // densities 
    double *rho_a = (double *)(h_densities +  0*NXY);
    double *rho_b = (double *)(h_densities +  1*NXY);
    double *tau_a = (double *)(h_densities +  2*NXY);
    double *tau_b = (double *)(h_densities +  3*NXY);
    double complex *nu = (double complex *)(h_densities +  4*NXY);
    double *j_a_x = (double *)(h_densities +  6*NXY);
    double *j_a_y = (double *)(h_densities +  7*NXY);
    double *j_a_z = (double *)(h_densities +  8*NXY);
    double *j_b_x = (double *)(h_densities +  9*NXY);
    double *j_b_y = (double *)(h_densities + 10*NXY);
    double *j_b_z = (double *)(h_densities + 11*NXY);
    
    // pontentials
    double *V_a = (double *)(h_potentials +  0*NXY);
    double *V_b = (double *)(h_potentials +  1*NXY);
    double complex *delta = (double complex *)(h_potentials +  2*NXY);
        
    // ====================================================================================
    // ================================ INITIAL STATE =====================================
    // ====================================================================================
    if(md.inittype==0 || md.inittype==1) // Start from uniform solution
    {
        if(md.inittype==0)
        {
            if(ip==0) printf("# CREATING UNIFORM SOLUTION...\n");
            
            // Generate initial state for testing
#ifdef BDG_MODE
            cpu_exec( solve_uniform_problem_bdg(md.Na/NXYZ, md.Nb/NXYZ, &nwf, ip==0) );
#else
            cpu_exec( solve_uniform_problem(md.Na/NXYZ, md.Nb/NXYZ, &nwf, ip==0) );
#endif
            cpu_exec( get_nwf_to_evolve_2d(&nwf) ); // correct number of states to evolve
            
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
            if(ip==0) { cpu_exec( read_uniform(&nwf, ip==0) ); cpu_exec( get_nwf_to_evolve_2d(&nwf) ); }
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
        cppmallocl(h_wavefun, NXY*nwfip*2,double complex);
        cppmallocl(h_fbetaEn, nwfip,double);
        cppmallocl(h_kkz, nwfip,double);
        cppmallocl(h_En, nwfip,double);
        
        // initialize wf
        cpu_exec( create_uniform_wf_2d(mylidx, myuidx, h_wavefun, &mu[SPINA], &mu[SPINB], &ec, h_fbetaEn, h_kkz, h_En, ip==0) );
        
        // initialize potentials
        for(ixyz=0; ixyz<NXY; ixyz++)
        {
            V_a[ixyz]=__md_pca_uniform.V_a; 
            V_b[ixyz]=__md_pca_uniform.V_b; 
            delta[ixyz]=__md_pca_uniform.delta + I*0.0;
        }
        
        // set eF, kF and Effg
        eF=pow(3.0*M_PI*M_PI*(__md_pca_uniform.n0_a+__md_pca_uniform.n0_b), 2.0/3.0) / 2.0;
        kF=pow(3.0*M_PI*M_PI*(__md_pca_uniform.n0_a+__md_pca_uniform.n0_b), 1.0/3.0);
        eF_a=pow(6.0*M_PI*M_PI*__md_pca_uniform.n0_a, 2.0/3.0) / 2.0;
        eF_b=pow(6.0*M_PI*M_PI*__md_pca_uniform.n0_b, 2.0/3.0) / 2.0;
        Effg = 0.6*__md_pca_uniform.n0_a*eF_a*NXYZ + 0.6*__md_pca_uniform.n0_b*eF_b*NXYZ; 
    }
    else if(md.inittype==3) 
    {
        // allocate memory for my wf
        load_nwf (MPI_COMM_WORLD, md.inprefix, &nwf, &nwfip, HowMany);
        cppmallocl(h_wavefun, NXYZ*nwfip*2,double complex);
        cppmallocl(h_fbetaEn, nwfip,double);
        cppmallocl(h_kkz, nwfip,double);
        printf("# WF SCATTER: ip=%d processes nwfip=%d wave-functions\n", ip, nwfip);

    }        
    else if(md.inittype==4) // Start from solution of kzSLpca solver 
    {
        // Load data from info file
        int _nx, _ny, _nz;
        double _dx, _dy, _dz;
        int *nwf_per_kz;
        cppmallocl(nwf_per_kz, NZ/2,int);
        
        time=0.0;
        sprintf(file_name, "%s_kzpca.info", md.inprefix);
        if(ip==0)
        {
            file_operation( read_checkpoint_info_pca(file_name, &nwf, &_nx, &_ny, &_nz, &_dx, &_dy, &_dz, &kF, &mu[0], &ec, &beta) );
            if(_nx!=NX || _ny!=NY || _nz!=NZ || _dx!=1.0 || _dy!=1.0 || _dz!=1.0)
            {
                printf("KZ-SOLVER INFO FILE NOT CONSISTENT GIVEN SETTINGS\n");
                printf("KZ-SOLVER: nx=%d, ny=%d, nz=%d, dx=%f, dy=%f, dz=%f\n", _nx, _ny, _nz, _dx, _dy, _dz);
                printf("SETTINGS  : nx=%d, ny=%d, nz=%d, dx=%f, dy=%f, dz=%f\n", NX, NY, NZ, 1.0, 1.0, 1.0);
                ABORT;
            }
            printf("KZ-SOLVER: file_name=`%s`\n",file_name);
            printf("KZ-SOLVER: nwf(including -kz's)=%d\n",nwf);
            printf("KZ-SOLVER: nx=%d, ny=%d, nz=%d, dx=%f, dy=%f, dz=%f\n", _nx, _ny, _nz, _dx, _dy, _dz);
            printf("KZ-SOLVER: kF=%f, mu_a=%f, mu_b=%f, ec=%f, beta=%f\n", kF, mu[SPINA], mu[SPINB], ec, beta);
            fflush(stdout);
            
            // scan files and determine nwf in each of them
            file_operation( scan_kzpca_info_files(md.inprefix, NZ, &nwf, nwf_per_kz) );
            printf("KZ-SOLVER: nwf(to evolve)=%d\n",nwf);
        }
        MPI_Bcast( &nwf , 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
        MPI_Bcast( &kF , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ; eF=0.5*kF*kF;
        MPI_Bcast( mu , 2 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ;
        MPI_Bcast( &ec , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ;
        MPI_Bcast( &beta , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ;  it=0; t0=0.0;
        MPI_Bcast( nwf_per_kz , NZ/2 , MPI_INT , 0 , MPI_COMM_WORLD ) ; 
                
        // divide wf over processes
        if(ip==0) printf("# INIT4: nwf=%d wave-functions to scatter\n", nwf);
        if ( np > nwf )
        {
            if(ip==0) printf("# INIT4: np[%d] > nwf[%d]!\n", np, nwf);
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
        cppmallocl(h_wavefun, NXY*nwfip*2,double complex);
        cppmallocl(h_fbetaEn, nwfip,double);
        cppmallocl(h_kkz,     nwfip,double);
        
        // read data
        int _4_max_readers = 32;
        int _4_nblocks = (int)ceil((float)(np)/_4_max_readers);
        for(i=0; i<_4_nblocks; i++)
        {
            if(ip==0) { printf("# INIT4: BLOCK ID[%d] CONSITING WITH %d PROCESSES READS DATA...\n", i, _4_max_readers); fflush(stdout);}
            if(ip%_4_nblocks == i) file_operation( read_kzSLpca_wf(md.inprefix, NZ, nwf_per_kz, mylidx, myuidx, h_wavefun, h_fbetaEn, h_kkz) );
            MPI_Barrier(MPI_COMM_WORLD);
        }
            
//         // TODO - ad hoc section - creates weights that select only Andreev states
//         cppmallocl(h_weight_andreev, nwfip, double);
//         #define ANDREEV_MAX_EN_CUTOFF 0.4
//         for(i=0; i<nwfip; i++) 
//             if(fabs(h_fbetaEn[i]/eF)<ANDREEV_MAX_EN_CUTOFF) h_weight_andreev[i]=1.0; 
//             else                                            h_weight_andreev[i]=0.0;
//             
//         gpu_exec( gpu_malloc(nwfip*sizeof(double), (void **)&d_weight_andreev) );
//         gpu_exec( memcopy_host2gpu(h_weight_andreev, d_weight_andreev,  (size_t)nwfip*sizeof(double)) ); 
//         
//         double andreev_states=0.0, andreev_states_t=0.0;
//         for(i=0; i<nwfip; i++) andreev_states+=h_weight_andreev[i];
//         MPI_Allreduce( &andreev_states, &andreev_states_t, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//         if(ip==0) printf("# ANDREEV # number of states=%.1f\n", andreev_states_t);
//         // TODO - end of ad hoc section
        
        for(i=0; i<nwfip; i++) h_fbetaEn[i]=fbeta(h_fbetaEn[i],beta); // convert quasiparticle energies into weights
        
        free(nwf_per_kz);
                
        // load u and delta
        if(ip==0)
        {
            sprintf(file_name, "%s_kzpca.pud", md.inprefix);
            printf("# INIT4: LOADING POTENTIALS `%s`...\n", file_name);
            file_operation( read_binary_file(file_name, NXY*4*sizeof(double), 0, h_potentials) );           
        }
        
        MPI_Bcast(h_potentials, 4*NXY, MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ;  
        
        // compute particle number and set Effg;
        double Ntota=0.0, Nmya=0.0;
        double Ntotb=0.0, Nmyb=0.0; 
        for(iwf=0; iwf<nwfip; iwf++)
        {
            double wcnt=2.0;
            if(fabs(h_kkz[iwf])<1.0e-12) wcnt = 1.0; // except for kz=0.0
            ixyz=0;
            for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ )
            {
                Nmya+=(pow(creal(h_wavefun[          iwf*NXY+ixyz]),2)+pow(cimag(h_wavefun[          iwf*NXY+ixyz]),2))*h_fbetaEn[iwf]*wcnt;
                Nmyb+=(pow(creal(h_wavefun[nwfip*NXY+iwf*NXY+ixyz]),2)+pow(cimag(h_wavefun[nwfip*NXY+iwf*NXY+ixyz]),2))*(1.0-h_fbetaEn[iwf])*wcnt;
                ixyz++;
            }  
        }
        MPI_Allreduce( &Nmya, &Ntota, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce( &Nmyb, &Ntotb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if(ip==0) printf("# INIT4: TOTAL NUMBER OF PARTICLES: SPIN_A=%16.8g SPIN_B=%16.8g\n", Ntota, Ntotb); 
        Effg = 0.6 * Ntota * eF;
        fflush(stdout);
    }
    else
    {
        if(ip==0) printf("NOT SUPPORTED INITTYPE=%d!\n", md.inittype);
        ABORT;
    }
    
    // Write info about distribution of wf
    MPI_Barrier(MPI_COMM_WORLD);
    printf("# WF SCATTER: ip=%d processes nwfip=%d wave-functions\n", ip, nwfip);
    MPI_Barrier(MPI_COMM_WORLD);

    // ====================================================================================
    // ======================== ALLOCATE GPU AND CPU BUFFERS ==============================
    // ====================================================================================    
    // Allocate memory for wave-functions and derivatives
    gpu_exec( gpu_malloc(NXY*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_wf) );
    gpu_exec( gpu_malloc(NXY*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_fkm1) );
    gpu_exec( gpu_malloc(NXY*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_fkm2) );
    gpu_exec( gpu_malloc(NXY*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_fkm3) );
#if INTEGRATION_SCHEME==AB4AM5
    gpu_exec( gpu_malloc(NXY*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_fkm4) );
#endif
    gpu_exec( gpu_malloc(NXY*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_wf_d_dx) );
    gpu_exec( gpu_malloc(NXY*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_wf_d_dy) );
    gpu_exec( gpu_malloc(NXY*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_wf_laplace) );
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
    gpu_exec( gpu_malloc(NXY*nwfip*2*sizeof(cufftDoubleComplex), (void **)&d_alphawf_laplace) );
#endif
    gpu_exec( gpu_malloc(nwfip*sizeof(double), (void **)&d_fbetaEn) );
    gpu_exec( gpu_malloc(nwfip*sizeof(double), (void **)&d_kkz) );
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
    if(md.inittype==3){ 

        if(ip==0) printf("# LOADING CHECKPOINT\n");
        b_t();
        size_t memsize;
#if INTEGRATION_SCHEME==AB3AM4
        load_all (h_wavefun, MPI_COMM_WORLD, md.inprefix,
                  d_wf, d_fkm1, d_fkm2, d_fkm3, 
		  d_potentials, &t0, 
                  &nwf, &nwfip,
                  h_fbetaEn, h_kkz, mu, &ec, &kF, &eF, &Effg,		  
		  HowMany);
        memsize = (size_t)(nwf)*(NX*NY*NZ)*2*4*16;
#elif INTEGRATION_SCHEME==AB4AM5
        load_all_45 (h_wavefun, MPI_COMM_WORLD, md.inprefix,
                     d_wf, d_fkm1, d_fkm2, d_fkm3, d_fkm4,
                     d_potentials, &t0,
                     &nwf, &nwfip,
                     h_fbetaEn, h_kkz, mu, &ec, &kF, &eF, &Effg,
                     HowMany);
        memsize = (size_t)(nwf)*(NX*NY*NZ)*2*5*16;
#else
        CHECK PCA_SETTINGS.H
#endif

        MPI_Barrier( MPI_COMM_WORLD ) ;
        rt = e_t(0);
        if(ip==0)
        {
            double memsize_gb = (double)(memsize) / pow(2,30);
            printf("# CHECKPOINT INFO: MODE=READ: DATA SIZE=%12.2f GB\n",  memsize_gb);
            printf("# CHECKPOINT INFO: HowMany=%d.\n", HowMany);
            printf("# CHECKPOINT INFO: READ TIME=%12.2f sec\n", rt);
            printf("# CHECKPOINT INFO: READ SPEED=%12.3f GB/sec\n", memsize_gb/rt);
        }
    }
    
    if(ip==0) printf("# INITIALIZING GPU BUFFERS OF ABM ALGORITHM...\n");
    
    if(md.inittype!=3)
    { 
        // copy wave-functions
        gpu_exec( memcopy_host2gpu(h_wavefun, d_wf,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) );   

        // we start from eigenstates - then fkm1, fkm2, fkm3 are zero
        for(i=0; i<NXY*nwfip*2; i++) h_wavefun[i]=0.0 + I*0.0;
        gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm1,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm2,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm3,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
#if INTEGRATION_SCHEME==AB4AM5
        gpu_exec( memcopy_host2gpu(h_wavefun, d_fkm4,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
#endif

        // copy potentials
        gpu_exec( memcopy_host2gpu(h_potentials, d_potentials,  (size_t)4*NXY*sizeof(double)) );  
    }
    // copy weights
    gpu_exec( memcopy_host2gpu(h_fbetaEn, d_fbetaEn,  (size_t)nwfip*sizeof(double)) ); 
    gpu_exec( memcopy_host2gpu(h_kkz    , d_kkz    ,  (size_t)nwfip*sizeof(double)) );
    
    // Set constants
    gpu_exec( memcopy_const(mu[SPINA], mu[SPINB], ec, t0, dt, kF) );    
    
    // Process params and copy them to gpu;
    process_params(md.params, kF);
    gpu_exec( memcopy_const_params(md.params) );
    
#ifdef WORK_IN_ROTATING_FRAME
    double Omega_a=0.0; // reset angular momentum
    double Omega_b=0.0; // reset angular momentum
    gpu_exec( memcopy_const_Omega(Omega_a, Omega_b) );
#endif
#ifdef BDG_MODE   
    gpu_exec( memcopy_const_BdG(md.aBdG) );
#endif
    
    if(ip==0) printf("# DONE.\n");
    
    
    // CUFFT plans
    size_t cufft_workSize;
    size_t wf_size = (size_t)2*NXY*nwfip*sizeof(double complex);
    if(md.batch>2*nwfip) md.batch=2*nwfip;
 
    // Create plans
    gpu_exec( create_cufftPlans(md.batch, nwfip, &cufft_workSize) );
    
    
    // Allocate memory for plans
    printf("# CUFFT[ip=%d]: cufft_workSize=%.2f times space of wf (%.2fMB)\n", ip, (double)cufft_workSize/(double)wf_size, (double)wf_size/pow(2.,20));
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
    gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, NULL, d_wf_laplace, md.nthreads) ); 
    // densities - local reduction
    gpu_exec( calculate_densities(nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, d_kkz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
    // densities - global reduction
    gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NXY*sizeof(double)) ); 
    MPI_Allreduce( MPI_IN_PLACE, h_densities, 12*NXY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)12*NXY*sizeof(double)) ); 
    // potentials
//     if(md.inittype!=3) gpu_exec( compute_potentials(it, d_densities, d_potentials, cccoeff, md.nthreads) );
    // energy
    gpu_exec( compute_energy(it, d_densities, d_potentials, d_workarea, md.nthreads) ); 
    
    // get data for reporting
    // energy
    gpu_exec( memcopy_gpu2host(d_workarea, h_energy,  (size_t)9*sizeof(double)) );   
    // potentials
    gpu_exec( memcopy_gpu2host(d_potentials, h_potentials,  (size_t)4*NXY*sizeof(double)) );     
    // densities - they are in h_densities
    double N_tot_init = h_energy[5]+h_energy[6]; // save initial value of particle number
    
    // report result
    if(ip==0)
    {
        time=t0+it*dt;
        energy_kin = h_energy[0];
        energy_pot = h_energy[1];
        energy_pair = h_energy[2];
        energy_CM = h_energy[3];
        energy_uext = h_energy[4];
        energy_tot = energy_kin+energy_pot+energy_pair+energy_CM+energy_uext;
        Na=h_energy[5]*NZ;
        Nb=h_energy[6]*NZ;
        Laz=h_energy[7];
        Lbz=h_energy[8];
        
        printf("# GPU ENERGY      : energy_kin=%16.12f, energy_pot=%16.12f, energy_pair=%16.12f, energy_tot=%16.12f, energy_CM=%16.12f, energy_uext=%16.12f\n", energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg, energy_CM/Effg, energy_uext/Effg);  
        printf("%12.4f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", time*eF, Na, Nb, Na+Nb, energy_tot/Effg, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_CM/Effg, energy_uext/Effg, Laz/Na, Lbz/Nb);     
        
        // Create run log and add entry
        cpu_exec( create_header_of_runlog(execcmd, kF, Effg, mu, ec, nwf, np, nwfip) );
#ifdef WORK_IN_ROTATING_FRAME
        #define OUTPUT_ENTRIES 16
#else
        #define OUTPUT_ENTRIES 14
#endif
        double line_items[OUTPUT_ENTRIES]={     
            time*eF, // 1
            Na, // 2
            Nb, // 3
            Na+Nb, // 4
            energy_tot/Effg, // 5
            energy_kin/Effg, // 6
            energy_pot/Effg, // 7
            energy_pair/Effg, // 8
            energy_CM/Effg, // 9
            energy_uext/Effg, //10
            qfalpha, //11
            Laz/Na, // 12
            Lbz/Nb, // 13
            (Laz+Lbz)/(Na+Nb), // 14
#ifdef WORK_IN_ROTATING_FRAME
            Omega_a, // 15
            Omega_b, // 16
#endif
        };
        cpu_exec( add_line_to_file(0, 0.0, OUTPUT_ENTRIES, line_items) );
    }    

    // Create binary files and add initial measurement
    if(ip==0)
    {
        // Create empty files with headers - do it once
        sprintf(file_name, "%s_density_a.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, t0, md.timesteps*dt) );
        sprintf(file_name, "%s_density_b.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, t0, md.timesteps*dt) );
        sprintf(file_name, "%s_delta.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, t0, md.timesteps*dt) );    
        sprintf(file_name, "%s_current_a.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, t0, md.timesteps*dt) );
        sprintf(file_name, "%s_current_b.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, t0, md.timesteps*dt) );   
        sprintf(file_name, "%s_qpe.dpca", md.outprefix);
        file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, t0, md.timesteps*dt) ); 
        
        // for each measurement add data to file
        sprintf(file_name, "%s_density_a.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, rho_a, sizeof(double)*NXY) );
        sprintf(file_name, "%s_density_b.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, rho_b, sizeof(double)*NXY) );
        sprintf(file_name, "%s_delta.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, delta, sizeof(double complex)*NXY) );
        sprintf(file_name, "%s_current_a.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, j_a_x, sizeof(double)*NXY*3) );
        sprintf(file_name, "%s_current_b.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, j_b_x, sizeof(double)*NXY*3) );
        // save zeros for qpe for initial measurement - to avoid expensive computation of qpe
        for(i=0; i<nwf; i++) h_qpe_nwf[i]=0.0;
        sprintf(file_name, "%s_qpe.dpca", md.outprefix);
        file_operation( add_measurement_entry(file_name, h_qpe_nwf, sizeof(double)*nwf) );   
    }

//     // TODO - ad hoc section
//     if(ip==0)
//     {
//         // Create empty files with headers - do it once
//         sprintf(file_name, "%s_density_andreev_a.dpca", md.outprefix);
//         file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, t0, md.timesteps*dt) );
//         sprintf(file_name, "%s_density_andreev_b.dpca", md.outprefix);
//         file_operation( create_measurement_file_with_header(file_name, NX, NY, 1, 1.0, 1.0, 1.0, eF, t0, md.timesteps*dt) );
//     }
//     
//     // densities - for andreev states only
//     gpu_exec( calculate_densities_weighted(nwfip, d_wf, d_kkz, d_fbetaEn, d_weight_andreev, d_densities_andreev, md.nthreads) );
//     // densities - global reduction
//     gpu_exec( memcopy_gpu2host(d_densities_andreev, h_densities_andreev,  (size_t)12*NXY*sizeof(double)) ); 
//     MPI_Allreduce( MPI_IN_PLACE, h_densities_andreev, 12*NXY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     
//     double andreev_particles[2], andreev_particles_total[2];
//     if(ip==0)
//     {
//         // for each measurement add data to file
//         sprintf(file_name, "%s_density_andreev_a.dpca", md.outprefix);
//         file_operation( add_measurement_entry(file_name, h_densities_andreev+0*NXY, sizeof(double)*NXY) );
//         sprintf(file_name, "%s_density_andreev_b.dpca", md.outprefix);
//         file_operation( add_measurement_entry(file_name, h_densities_andreev+1*NXY, sizeof(double)*NXY) );
//         
//         for(j=0; j<2; j++) // spin up and down
//         {
//             andreev_particles[j]=0.0;
//             for(ixyz=0; ixyz<NXY; ixyz++) andreev_particles[j]+=h_densities_andreev[j*NXY + ixyz];
//         }
//         
//         printf("# ANDREEV %12.4f %12.8f %12.8f\n", time*eF, andreev_particles[0], andreev_particles[1]);
//     }
//     // TODO - end of ad hoc section
    
    // ====================================================================================
    // ================================ REAL TIME EVOLUTION  ==============================
    // ====================================================================================
    int i_meas, i_step;
    
    if(md.inittype!=3 && md.selfstart==1)
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
        cppmallocl(h_fkm, NXY*nwfip*2*selfstart_steps,double complex);

        for(i_step=0; i_step<selfstart_steps; i_step++)
        {
            // ----------------------------- predictor -----------------------------------
            // compute value of quantum friction coefficient
            qfalpha = 0.0; // NOTE: I assume there is no quantum friction during the first steps
            cccoeff = h_smooth_step(t0+it*dt, md.ccstart/eF,  md.ccstop/eF,  md.ccswitch/eF, 1.0);
#ifdef WORK_IN_ROTATING_FRAME
            // compute rotation velocity
            gpu_exec( set_Omega(md.params, kF, t0, &Omega_a, &Omega_b) );
#endif
            // derivatives
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
            if(qfalpha>0.0) gradients_computed=1; else gradients_computed=0;
#endif
            if(gradients_computed)
            {
                gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, NULL, d_wf_laplace, md.nthreads) ); 
            }
            else
            {
                gpu_exec( compute_laplace(2*nwfip, d_wf, d_wf_laplace, md.nthreads) );
            }
            // densities - local reduction
            gpu_exec( calculate_densities(nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, d_kkz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
            // densities - global reduction
            gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NXY*sizeof(double)) ); 
            MPI_Allreduce( MPI_IN_PLACE, h_densities, 12*NXY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)12*NXY*sizeof(double)) );
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // effective mass correction
            gpu_exec( multiply_wf_by_alpha(nwfip, d_wf, d_alphawf_laplace, d_densities, md.nthreads) );
            gpu_exec( compute_laplace(2*nwfip, d_alphawf_laplace, d_alphawf_laplace, md.nthreads) );
#endif
            // potentials - NOTE: it=0!
            gpu_exec( compute_potentials(0, d_densities, d_potentials, cccoeff, md.nthreads) );
            
            // executing exp[-i*H(t)*dt]*psi
            // H*psi - first execution, d_fkm3 as working buffer
            gpu_exec( memcopy_gpu2gpu(d_wf, d_fkm3, (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) );
            gpu_exec( apply_hamiltonian(nwfip, d_fkm3, d_fkm1, /* NOTE - d_fkm1 as output buffer  */
                                    d_wf_d_dx, d_wf_d_dy, d_kkz, d_wf_laplace, d_alphawf_laplace,
                                    d_densities, d_potentials, qfalpha, NULL, cccoeff, 
                                    md.nthreads) );
            // Make copy of qpe
            gpu_exec( memcopy_gpu2gpu(d_workarea, d_qpe, (size_t)nwfip*sizeof(double)) );
            
            // Store H*Psi
            gpu_exec( memcopy_gpu2host(d_fkm1, h_fkm+i_step*2*nwfip*NXY,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
                        
            // Add contribution from Taylor expansion
            gpu_exec( taylor_expansion_contribution(1, 0.5*dt, nwfip, d_fkm1, d_fkm3, d_fkm2, md.nthreads) );
                            
            // H*psi - remaining executions
            for(i_meas=1; i_meas<exp_iters; i_meas++)
            {
                // derivatives
                if(gradients_computed)
                {
                    gpu_exec( compute_derivatives(2*nwfip, d_fkm2, d_wf_d_dx, d_wf_d_dy, NULL, d_wf_laplace, md.nthreads) ); 
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
                gpu_exec( apply_hamiltonian(nwfip, d_fkm2, d_fkm1,
                                        d_wf_d_dx, d_wf_d_dy, d_kkz, d_wf_laplace, d_alphawf_laplace,
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
#ifdef WORK_IN_ROTATING_FRAME
            // compute rotation velocity
            gpu_exec( set_Omega(md.params, kF, t0, &Omega_a, &Omega_b) );
#endif
            // derivatives
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
            if(qfalpha>0.0) gradients_computed=1; else gradients_computed=0;
#endif
            if(gradients_computed)
            {
                gpu_exec( compute_derivatives(2*nwfip, d_fkm3, d_wf_d_dx, d_wf_d_dy, NULL, d_wf_laplace, md.nthreads) ); 
            }
            else
            {
                gpu_exec( compute_laplace(2*nwfip, d_fkm3, d_wf_laplace, md.nthreads) );
            }
            // densities - local reduction
            gpu_exec( calculate_densities(nwfip, d_fkm3, d_wf_d_dx, d_wf_d_dy, d_kkz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
            // densities - global reduction
            gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NXY*sizeof(double)) ); 
            MPI_Allreduce( MPI_IN_PLACE, h_densities, 12*NXY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)12*NXY*sizeof(double)) );
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
                gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, NULL, d_wf_laplace, md.nthreads) ); 
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
            gpu_exec( apply_hamiltonian(nwfip, d_wf, d_fkm1, /* NOTE - d_fkm1 as output buffer  */
                                    d_wf_d_dx, d_wf_d_dy, d_kkz, d_wf_laplace, d_alphawf_laplace,
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
                    gpu_exec( compute_derivatives(2*nwfip, d_fkm2, d_wf_d_dx, d_wf_d_dy, NULL, d_wf_laplace, md.nthreads) ); 
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
                gpu_exec( apply_hamiltonian(nwfip, d_fkm2, d_fkm1,
                                        d_wf_d_dx, d_wf_d_dy, d_kkz, d_wf_laplace, d_alphawf_laplace,
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
        gpu_exec( memcopy_host2gpu(h_fkm+2*2*nwfip*NXY, d_fkm1,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_fkm+1*2*nwfip*NXY, d_fkm2,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_fkm+0*2*nwfip*NXY, d_fkm3,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
#elif INTEGRATION_SCHEME==AB4AM5
        gpu_exec( memcopy_host2gpu(h_fkm+3*2*nwfip*NXY, d_fkm1,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_fkm+2*2*nwfip*NXY, d_fkm2,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_fkm+1*2*nwfip*NXY, d_fkm3,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
        gpu_exec( memcopy_host2gpu(h_fkm+0*2*nwfip*NXY, d_fkm4,  (size_t)2*nwfip*NXY*sizeof(cufftDoubleComplex)) ); 
#endif        
        
        // clear memory
        gpu_exec( gpu_free(d_qpe) );
        free(h_fkm);
        if(ip==0) printf("# SELFSTART: DONE.\n");
         
        // derivatives
        gradients_computed=1;
        gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, NULL, d_wf_laplace, md.nthreads) ); 
        // densities - local reduction
        gpu_exec( calculate_densities(nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, d_kkz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
        // densities - global reduction
        gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NXY*sizeof(double)) ); 
        MPI_Allreduce( MPI_IN_PLACE, h_densities, 12*NXY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)12*NXY*sizeof(double)) ); 
        // potentials
        gpu_exec( compute_potentials(it, d_densities, d_potentials, cccoeff, md.nthreads) );
        // energy
        gpu_exec( compute_energy(it, d_densities, d_potentials, d_workarea, md.nthreads) );
    
        // get data for reporting
        // energy
        gpu_exec( memcopy_gpu2host(d_workarea, h_energy,  (size_t)9*sizeof(double)) );   
        // potentials
        gpu_exec( memcopy_gpu2host(d_potentials, h_potentials,  (size_t)4*NXY*sizeof(double)) );     
        // densities - they are in h_densities
    
        // report result
        if(ip==0)
        {
            time=t0+it*dt;
            energy_kin = h_energy[0];
            energy_pot = h_energy[1];
            energy_pair = h_energy[2];
            energy_CM = h_energy[3];
            energy_uext = h_energy[4];
            energy_tot = energy_kin+energy_pot+energy_pair+energy_CM+energy_uext;
            Na=h_energy[5]*NZ;
            Nb=h_energy[6]*NZ;
            Laz=h_energy[7];
            Lbz=h_energy[8];
            
            printf("# AFTER SELFSTART : energy_kin=%16.12f, energy_pot=%16.12f, energy_pair=%16.12f, energy_tot=%16.12f, energy_CM=%16.12f, energy_uext=%16.12f\n", energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg, energy_CM/Effg, energy_uext/Effg);  
            printf("%12.4f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", time*eF, Na, Nb, Na+Nb, energy_tot/Effg, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_CM/Effg, energy_uext/Effg, Laz/Na, Lbz/Nb);     
            
            // Create run log and add entry
            cpu_exec( create_header_of_runlog(execcmd, kF, Effg, mu, ec, nwf, np, nwfip) );

            double line_items[OUTPUT_ENTRIES]={     
                time*eF, // 1
                Na, // 2
                Nb, // 3
                Na+Nb, // 4
                energy_tot/Effg, // 5
                energy_kin/Effg, // 6
                energy_pot/Effg, // 7
                energy_pair/Effg, // 8
                energy_CM/Effg, // 9
                energy_uext/Effg, //10
                qfalpha, //11
                Laz/Na, // 12
                Lbz/Nb, // 13
                (Laz+Lbz)/(Na+Nb), // 14
#ifdef WORK_IN_ROTATING_FRAME
                Omega_a, // 15
                Omega_b, // 16
#endif
            };
            cpu_exec( add_line_to_file(0, 0.0, OUTPUT_ENTRIES, line_items) );
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
#ifdef WORK_IN_ROTATING_FRAME
            // compute rotation velocity
            gpu_exec( set_Omega(md.params, kF, t0+(it+1)*dt, &Omega_a, &Omega_b) );
#endif
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
                gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, NULL, d_wf_laplace, md.nthreads) ); 
            }
            else
            {
                gpu_exec( compute_laplace(2*nwfip, d_wf, d_wf_laplace, md.nthreads) );
            }
            // densities - local reduction
            gpu_exec( calculate_densities(nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, d_kkz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
            // densities - global reduction
            gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NXY*sizeof(double)) ); 
            MPI_Allreduce( MPI_IN_PLACE, h_densities, 12*NXY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)12*NXY*sizeof(double)) ); 
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // effective mass correction
            gpu_exec( multiply_wf_by_alpha(nwfip, d_wf, d_alphawf_laplace, d_densities, md.nthreads) );
            gpu_exec( compute_laplace(2*nwfip, d_alphawf_laplace, d_alphawf_laplace, md.nthreads) );
#endif
            // potentials
            gpu_exec( compute_potentials(it+1, d_densities, d_potentials, cccoeff, md.nthreads) );
            // H*psi
            gpu_exec( apply_hamiltonian(nwfip, d_wf, d_wf_laplace, /* NOTE - d_wf_laplace as output buffer  */
                                    d_wf_d_dx, d_wf_d_dy, d_kkz, d_wf_laplace, d_alphawf_laplace,
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
                gpu_exec( compute_derivatives(2*nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, NULL, d_wf_laplace, md.nthreads) ); 
            }
            else
            {
                gpu_exec( compute_laplace(2*nwfip, d_wf, d_wf_laplace, md.nthreads) );
            }
            // densities - local reduction
            gpu_exec( calculate_densities(nwfip, d_wf, d_wf_d_dx, d_wf_d_dy, d_kkz, d_fbetaEn, d_densities, gradients_computed, md.nthreads) );
            // densities - global reduction
            gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)12*NXY*sizeof(double)) ); 
            MPI_Allreduce( MPI_IN_PLACE, h_densities, 12*NXY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)12*NXY*sizeof(double)) ); 
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
            gpu_exec( apply_hamiltonian(nwfip, d_wf, d_fkm1, 
                                    d_wf_d_dx, d_wf_d_dy, d_kkz, d_wf_laplace, d_alphawf_laplace, 
                                    d_densities, d_potentials, qfalpha, NULL, cccoeff, 
                                    md.nthreads) );            
            
            it++; // update global time counter
        }
        
        // ----------------------------- measurement -------------------------------------
        // Save quasiparticle energies - computation of energy will destroy them
        gpu_exec( memcopy_gpu2host(d_workarea, h_qpe_nwfip,  (size_t)nwfip*sizeof(double)) );
        MPI_Gatherv(h_qpe_nwfip,nwfip,MPI_DOUBLE,h_qpe_nwf,wf_tbl,wf_idx_tbl,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if(ip==0)
        {
            sprintf(file_name, "%s_qpe.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, h_qpe_nwf, sizeof(double)*nwf) ); 
        }
        
        // energy
        gpu_exec( compute_energy(it, d_densities, d_potentials, d_workarea, md.nthreads) );
        
        // get data for reporting
        // energy
        gpu_exec( memcopy_gpu2host(d_workarea, h_energy,  (size_t)9*sizeof(double)) );   
        // potentials
        gpu_exec( memcopy_gpu2host(d_potentials, h_potentials,  (size_t)4*NXY*sizeof(double)) );     
        // densities - they are in h_densities
        
        rt=e_t(0); // get timing
        
        // report result
        if(ip==0)
        {
            time=t0+it*dt;
            energy_kin = h_energy[0];
            energy_pot = h_energy[1];
            energy_pair = h_energy[2];
            energy_CM = h_energy[3];
            energy_uext = h_energy[4];
            energy_tot = energy_kin+energy_pot+energy_pair+energy_CM+energy_uext;
            Na=h_energy[5]*NZ;
            Nb=h_energy[6]*NZ;
            Laz=h_energy[7];
            Lbz=h_energy[8];
            
            printf("%12.4f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %6.3f %8.2f\n", time*eF, Na, Nb, Na+Nb, energy_tot/Effg, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_CM/Effg, energy_uext/Effg, Laz/Na, Lbz/Nb, qfalpha, rt);        
            
            double line_items[OUTPUT_ENTRIES]={     
                time*eF, // 1
                Na, // 2
                Nb, // 3
                Na+Nb, // 4
                energy_tot/Effg, // 5
                energy_kin/Effg, // 6
                energy_pot/Effg, // 7
                energy_pair/Effg, // 8
                energy_CM/Effg, // 9
                energy_uext/Effg, //10
                qfalpha, //11
                Laz/Na, // 12
                Lbz/Nb, // 13
                (Laz+Lbz)/(Na+Nb), // 14
#ifdef WORK_IN_ROTATING_FRAME
                Omega_a, // 15
                Omega_b, // 16
#endif
            };
            cpu_exec( add_line_to_file(i_meas+1, rt, OUTPUT_ENTRIES, line_items) );
        } 
        
        // add binary data
        if(ip==0)
        {
            // for each measurement add data to file
            sprintf(file_name, "%s_density_a.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, rho_a, sizeof(double)*NXY) );
            sprintf(file_name, "%s_density_b.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, rho_b, sizeof(double)*NXY) );
            sprintf(file_name, "%s_delta.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, delta, sizeof(double complex)*NXY) );
            sprintf(file_name, "%s_current_a.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, j_a_x, sizeof(double)*NXY*3) );
            sprintf(file_name, "%s_current_b.dpca", md.outprefix);
            file_operation( add_measurement_entry(file_name, j_b_x, sizeof(double)*NXY*3) );
        }
        
//         // TODO ad hoc section
//         // densities - for andreev states only
//         gpu_exec( calculate_densities_weighted(nwfip, d_wf, d_kkz, d_fbetaEn, d_weight_andreev, d_densities_andreev, md.nthreads) );
//         // densities - global reduction
//         gpu_exec( memcopy_gpu2host(d_densities_andreev, h_densities_andreev,  (size_t)12*NXY*sizeof(double)) ); 
//         MPI_Allreduce( MPI_IN_PLACE, h_densities_andreev, 12*NXY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//         
//         if(ip==0)
//         {
//             // for each measurement add data to file
//             sprintf(file_name, "%s_density_andreev_a.dpca", md.outprefix);
//             file_operation( add_measurement_entry(file_name, h_densities_andreev+0*NXY, sizeof(double)*NXY) );
//             sprintf(file_name, "%s_density_andreev_b.dpca", md.outprefix);
//             file_operation( add_measurement_entry(file_name, h_densities_andreev+1*NXY, sizeof(double)*NXY) );
//             
//             for(j=0; j<2; j++) // spin up and down
//             {
//                 andreev_particles[j]=0.0;
//                 for(ixyz=0; ixyz<NXY; ixyz++) andreev_particles[j]+=h_densities_andreev[j*NXY + ixyz];
//             }
//             
//             printf("# ANDREEV %12.4f %12.8f %12.8f\n", time*eF, andreev_particles[0], andreev_particles[1]);
//         }
//         // TODO end of ad hoc section

        // TODO ad hoc section
        if(ip==0)
        {
            sprintf(file_name, "%s_checkpoint.kzpca", md.outprefix);
//             printf("# CREATING CHECKPOINT FILE `%s`\n", file_name);
            FILE * pFile = fopen(file_name, "wb");
            double npart[2];
            npart[SPINA]=Na;
            npart[SPINB]=Nb;
            beta=1./(0.01*eF); // set by hand
            
            // write all nescesary data to file
            i=0;
            fwrite(&i           , sizeof(int)         , 1 , pFile); // iteration number
            fwrite(&mu[SPINA]   , sizeof(double)      , 1 , pFile); 
            fwrite(&mu[SPINB]   , sizeof(double)      , 1 , pFile); 
            fwrite(&ec          , sizeof(double)      , 1 , pFile); 
            fwrite(&beta        , sizeof(double)      , 1 , pFile); 
            fwrite(&eF          , sizeof(double)      , 1 , pFile); 
            fwrite(&kF          , sizeof(double)      , 1 , pFile);
            fwrite(&Effg        , sizeof(double)      , 1 , pFile);
            fwrite(h_potentials , sizeof(double)*NX*NY, 4 , pFile);
            fwrite(h_densities  , sizeof(double)*NX*NY, 12, pFile);
            fwrite(h_energy     , sizeof(double)      , 5 , pFile);
            fwrite(npart        , sizeof(double)      , 2 , pFile);
                  
            fclose(pFile);
        }
        // TODO end of ad hoc section
        
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
        if( fabs( (h_energy[5]+h_energy[6]-N_tot_init)/N_tot_init )>N_STABILITY_CRITERIA )
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
		h_fbetaEn, h_kkz, mu, &ec, &kF, & eF, &Effg,
		HowMany);
       memsize = (size_t)(nwf)*(NX*NY*NZ)*2*4*16;
#elif INTEGRATION_SCHEME==AB4AM5            
       save_all_45(h_wavefun, MPI_COMM_WORLD, md.outprefix,
                   d_wf, d_fkm1, d_fkm2, d_fkm3, d_fkm4,
                   d_potentials, &time,
                   nwf, nwfip,
                   h_fbetaEn, h_kkz, mu, &ec, &kF, & eF, &Effg,
                   HowMany);
       memsize = (size_t)(nwf)*(NX*NY*NZ)*2*5*16;
#else
            CHECK PCA_SETTINGS.H
#endif
    MPI_Barrier( MPI_COMM_WORLD ) ;
    rt = e_t(0);
    if(ip==0)
    {
        double memsize_gb = (double)(memsize) / pow(2,30);
        printf("# CHECKPOINT INFO: MODE=WRITE: DATA SIZE=%12.2f GB\n",  memsize_gb);
        printf("# CHECKPOINT INFO: HowMany=%d.\n", HowMany);
        printf("# CHECKPOINT INFO: WRITE TIME=%12.2f sec\n", rt);
        printf("# CHECKPOINT INFO: WRITE SPEED=%12.3f GB/sec\n", memsize_gb/rt);
    }
}
    /* messy exit here */
    MPI_Barrier( MPI_COMM_WORLD ) ;
    MPI_Finalize() ;
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
