#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "predefines.h"
#include "gpe_utils.h"
#include "pca_utils.h"
#include "gpe_engine_api.h"
#include "wslda_reproducibility.h"

int wsldapid;

void read_of_input_parameters(int argc , char ** argv)
{
    int i;
    char execcmd[ 256 ];
    strcpy( execcmd , argv[ 0 ] ) ;
    for( i = 1 ; i < argc ; i++ ) 
    {
        strcat( execcmd , " " ) ; 
        strcat( execcmd , argv[ i ] ) ;
    }
}

int parse_command_line_and_get_idx_of_input_file(int argc , char ** argv)
{
    int i = readcmd( argc , argv ) ;
    if( i == -1 )
    {
        printf( "TERMINATING! NO INPUT FILE.\n" ) ;
        return(EXIT_FAILURE);
    }
    return i;
}

void read_input_file(int idx, char ** argv)
{
    // Read input file
    // Info from file is loaded into metadata structure
    int j = parse_input_file(argv[idx]);
    if ( j == 0 )
    {
        printf("PROBLEM WITH INPUT FILE: `%s`.\n" , argv[ idx ] ) ;
        exit(EXIT_FAILURE);      
    }
        
    // Input file tags are accessible through pointer `input`
}

int main( int argc , char ** argv ) 
{
    read_of_input_parameters(argc, argv);
    int input_idx = parse_command_line_and_get_idx_of_input_file(argc, argv);
    read_input_file(input_idx, argv);
    
//     file_operation( check_if_can_overwrite_files() ); // terminate if file exists, and input->overwrite==0
//     sprintf(file_name, "%s_input.txt", md.outprefix);
//     file_operation( copy_input_file(argv[i],file_name) ); 
    
    double ekin, eint, eext, etot, etot_prev, time, rt;
    double alpha=input->alpha_real;
    double beta=input->beta_real;
    double dt=input->dt;
    double npart=input->npart;
    const int device=input->gpuDevice;  
    
    set_gpu_device(device);

    int nx, ny, nz, ierr;
    gpe_get_lattice_api(&nx, &ny, &nz);
    printf("# GPE engine compiled for lattice: %d x %d x %d\n", nx, ny, nz);
    printf("# REAL TIME EVOLUTION\n");

    double complex *psi;
    uint nxyz=nx*ny*nz;
    alloc_host_memory(nxyz, &psi);
    // TODO
    if(input->iniitype==0)
    {
        set_initial_wave_function( nxyz, psi);
    }
    else if (input->iniitype==5)
    {
        // TODO: read from `inprefix`_psi.wdat
        read_initial_wave_function( nxyz, psi);
    }
    else
    {
        // error
    }
    gpe_create_engine_api(alpha, beta, dt, npart);
    gpe_set_user_params_api(MAX_USER_PARAMS, input->params);
    gpe_set_psi_api(0.0, psi);

    print_header();
    gpe_energy_api(&time, &ekin, &eint, &eext);

    etot = ekin + eint + eext;
    print_intial_results(time, npart, etot, ekin, eint, eext);
    
    // Create empty WDATA set
    // use example:
    // https://gitlab.fizyka.pw.edu.pl/wtools/wdata/-/blob/master/c-examples/example-write.c

    // TODO
    // if(mode=...)
    while(time <= input->timesteps)
    {
        b_t(); // reset timer
        // Psi(t)
        gpe_evolve_api(input->timesteps);
        // Psi(t+timesteps*dt)
        
        gpe_energy_api(&time, &ekin, &eint, &eext);
        rt = e_t(0); // get time

        etot = ekin + eint + eext;
        print_intial_results(time, npart, etot, ekin, eint, eext);
        
        gpe_get_psi_api(&time, psi);
        
        // Add new data to WDATA set
        // outprefix_psi.wdat        
//         gpe_get_density(&time, psi);
        // outprefix_density.wdat
//         gpe_get_current(&time, psi);
        // outprefix_current.wdat
        
        // add enetry to logger
    }

    // close files
    
    gpe_destroy_engine_api();
    free_host_memory(psi);
    return 0;
}
