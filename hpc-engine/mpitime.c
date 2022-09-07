/*
 * W-SLDA Toolkit
 * This code is for testing performance of MPI communication
 * */

#define TDWSLDA_MAIN

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
#include "pca_derivative.h"
#include "pca_kernels.h"
#include "pca_utils.h"
#include "wslda_potdens.h"

int wsldapid; // process id - global variable
int wsldapnp; // total number of processes - global variable
#include "tdwslda_static_vars.h"
#define printf wprintf

int main( int argc , char ** argv )
{
    int ierr; // error flag
    int ip, np; // basic MPI indicators
    int i,j;

    MPI_Init( &argc , &argv ) ; /* set up the parallel WORLD */
    MPI_Comm_size( MPI_COMM_WORLD , &np ) ; /* total number of processes */
    MPI_Comm_rank( MPI_COMM_WORLD , &ip ) ; /* id of process st 0 <= ip < np */
    wsldapid=ip; // save to global variable
    wsldapnp=np; // save to global variable

    if(ip==0) wprintf("# START OF THE MPITIME\n");

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
    }

    // Broadcast input parameter
    MPI_Bcast( &md , sizeof(md) , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;

    int dim=atoi(argv[1]);
    int ixyz, nxyz;


    // header
    if(dim==1)
    {
        if(ip==0) printf("# LATTICE 1D WILL BE TESTED: %d\n", NX);
        nxyz=NX;
    }
    else if(dim==2)
    {
        if(ip==0) printf("# LATTICE 2D WILL BE TESTED: %d x %d = %d\n", NX, NY, NXY);
        nxyz=NXY;
    }
    else
    {
        if(ip==0) printf("# LATTICE 3D WILL BE TESTED: %d x %d x %d = %d\n", NX, NY, NZ, NXYZ);
        nxyz=NXYZ;
    }

    int mpipackagesize = EXCHANGE_SIZE;
    if(ip==0) printf("# EXCHANGE PACKAGE SIZE: %d x %d x %ld B= %.3f MB\n", mpipackagesize, nxyz, sizeof(double), 1.0*mpipackagesize*nxyz*sizeof(double)/1024/1024);
#ifdef USE_GPU_AWARE_MPI
    if(ip==0) printf("# USE_GPU_AWARE_MPI: YES\n");
#else
    if(ip==0) printf("# USE_GPU_AWARE_MPI: NO\n");
#endif
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // ====================================================================================
    // ============================= INITIALIZE GPU =======================================
    // ====================================================================================
    int deviceId=0;
#ifdef CUSTOM_GPU_DISTRIBUTION
    if(ip==0) wprintf("# EXECUTING `assign_deviceid_to_mpi_process` TO GET GPUS DISTRIBUTION ACROSS THE SYSTEM\n");
    deviceId = assign_deviceid_to_mpi_process(MPI_COMM_WORLD);
#else
    if(ip==0) wprintf("# ASSUMING STANDARD DISTRIBUTION OF GPUS ACROSS THE SYSTEM [gpuspernode=%d]\n", md.gpuspernode);
    deviceId = ip % md.gpuspernode;
#endif

#ifdef PRINT_GPU_DISTRIBUTION
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    wprintf("# PROCESS ip=%d RUNNING ON NODE %s USES device-id=%d\n", ip, processor_name, deviceId);
#endif
    gpu_exec( set_gpu(deviceId) );
    fflush(stdout);

    // ====================================================================================
    // ======================== ALLOCATE GPU AND CPU BUFFERS ==============================
    // ====================================================================================
    double *h_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *d_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (GPU)
    gpu_exec( host_malloc_pl((size_t)12*nxyz*sizeof(double), (void **)&h_densities) );
    gpu_exec(     gpu_malloc((size_t)12*nxyz*sizeof(double), (void **)&d_densities) );


    for(ixyz=0; ixyz<12*nxyz; ixyz++) h_densities[i]=1.0;
    gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)mpipackagesize*nxyz*sizeof(double)) ); // set array on gpu

    // ---------------------- TESTING LOOP---------------------------
    int i_meas, i_step;
    double rt_tot=0.0, rt_tot2=0.0;
    if(ip==0) printf("# NUMBER OF TESTS: measurements=%d WITH timesteps=%d\n", md.measurements, md.timesteps);
    if(ip==0) printf("# %4s %12s\n", "it", "t [sec]");
    if(ip==0) printf("# ------------------\n");
    for (i_meas=0; i_meas<md.measurements; i_meas++)
    {
        b_t(); // reset timer
        for(i_step=0; i_step<md.timesteps*2; i_step++) // factor 2-(predictor, corrector)
        {
            // densities - global reduction
#ifdef USE_GPU_AWARE_MPI
            MPI_Allreduce( MPI_IN_PLACE, d_densities, mpipackagesize*nxyz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
            gpu_exec( memcopy_gpu2host(d_densities, h_densities,  (size_t)mpipackagesize*nxyz*sizeof(double)) );
            MPI_Allreduce( MPI_IN_PLACE, h_densities, mpipackagesize*nxyz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            gpu_exec( memcopy_host2gpu(h_densities, d_densities,  (size_t)mpipackagesize*nxyz*sizeof(double)) );
#endif
        }

        double rt=e_t(0); // get timing
        rt_tot+=rt;
        rt_tot2+=rt*rt;

        if(ip==0) printf("%6d %12.3f\n", i_meas, rt);
    }

    if(ip==0) printf("# ------------------\n");
    double avg = rt_tot / md.measurements;
    double std = sqrt(rt_tot2/md.measurements - avg*avg);
    if(ip==0) printf("# TIMING PER MEASURMENT [sec]: %12.4f +/- %12.4f\n", avg,std);
    fflush(stdout);

    /* messy exit here */
    MPI_Barrier( MPI_COMM_WORLD ) ;
    MPI_Finalize() ;
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
