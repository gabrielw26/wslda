#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "predefines.h"
#include "gpe_utils.h"
#include "pca_utils.h"
#include "gpe_engine_api.h"
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
    
    double ekin, eint, eext, etot, etot_prev, time, rt;
    double alpha=input->alpha_imag;
    double beta=input->beta_imag;
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
    set_initial_wave_function( nxyz, psi);
    gpe_create_engine_api(alpha, beta, dt, npart);
    gpe_set_user_params_api(3, input->params);
    gpe_set_psi_api(0.0, psi);
    gpe_normalize_psi_api();

    print_header();
    gpe_energy_api(&time, &ekin, &eint, &eext);

    etot = ekin + eint + eext;
    print_intial_results(time, npart, etot, ekin, eint, eext);

    while(1)
    {
        b_t(); // reset timer
        gpe_evolve_api(input->timesteps);
        gpe_energy_api(&time, &ekin, &eint, &eext);
        rt = e_t(0); // get time

        etot_prev=etot;
        etot = ekin + eint + eext;
        double diff=(etot_prev-etot)/npart; // diference in energy per particle
        print_results(time, npart, etot, ekin, eint, eext, diff, rt);
        
        if(fabs(diff) < input->energyconveps) break;
    }

    gpe_get_psi_api(&time, psi);
    write_to_binary_file(nxyz, psi);
    write_to_txt_file(nx, ny, nz, psi);
    gpe_destroy_engine_api();
    free_host_memory(psi);
    return 0;
}