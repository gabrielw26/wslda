#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "predefines.h"
#include "gpe_utils.h"
#include "pca_utils.h"
#include "gpe_engine_api.h"
#include "wslda_reproducibility.h"
int wsldapid;

int main( int argc , char ** argv ) 
{
    read_of_input_parameters(argc, argv);
    int input_idx = parse_command_line_and_get_idx_of_input_file(argc, argv);
    read_input_file(input_idx, argv);
    
//     file_operation( check_if_can_overwrite_files() ); // terminate if file exists, and input->overwrite==0
//     sprintf(file_name, "%s_input.txt", md.outprefix);
//     file_operation( copy_input_file(argv[i],file_name) ); 
    
    double ekin, eint, eext, etot, etot_prev, time, rt, diff=-1;
    double alpha=input->alpha;
    double beta=input->beta;
    double dt=input->dt;
    double npart=input->npart;
    const int device=input->gpuDevice;
    const int mode=input->gpe_mode;
    
    set_gpu_device(device);

    int nx, ny, nz, ierr;
    gpe_get_lattice_api(&nx, &ny, &nz);
    printf("# GPE engine compiled for lattice: %d x %d x %d\n", nx, ny, nz);
    if(mode==0)  printf("# IMAGINARY TIME PROJECTION\n");
    if(mode==1)  printf("# REAL TIME EVOLUTION\n");

    Complex *psi;
    uint nxyz=nx*ny*nz;
    alloc_host_memory(nxyz, &psi);
    // TODO
    if(mode==0) {
        set_initial_wave_function( nxyz, psi);
    } else if (mode==1) {
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
    if(mode==0) gpe_normalize_psi_api();

    print_header();
    gpe_energy_api(&time, &ekin, &eint, &eext);

    etot = ekin + eint + eext;
    print_intial_results(time, npart, etot, ekin, eint, eext);
    
    // Create empty WDATA set
    // use example:
    // https://gitlab.fizyka.pw.edu.pl/wtools/wdata/-/blob/master/c-examples/example-write.c

    // TODO
    // if(mode=...)
    while(1)
    {
        b_t(); // reset timer
        // Psi(t)
        gpe_evolve_api(input->timesteps);
        // Psi(t+timesteps*dt)
        
        gpe_energy_api(&time, &ekin, &eint, &eext);
        rt = e_t(0); // get time

        if(mode==0) { //TOREMOVE
            etot_prev=etot;
            etot = ekin + eint + eext;
            diff=(etot_prev-etot)/npart; // diference in energy per particle
            print_results(time, npart, etot, ekin, eint, eext, diff, rt);
        } else if(mode==1) { //TOREMOVE
            etot = ekin + eint + eext;
            print_intial_results(time, npart, etot, ekin, eint, eext);
        }

        
        // Add new data to WDATA set
        // outprefix_psi.wdat        
//         gpe_get_density(&time, psi);
        // outprefix_density.wdat
//         gpe_get_current(&time, psi);
        // outprefix_current.wdat
        
        // add enetry to logger

        if(mode==0 && fabs(diff) < input->energyconveps) break;
        if(mode==1 && time > input->timesteps) break;
    }

    // close files
    
    gpe_get_psi_api(&time, psi); //TOREMOVE
    write_to_binary_file(nxyz, psi); //TOREMOVE
    write_to_txt_file(nx, ny, nz, psi); //TOREMOVE
    gpe_destroy_engine_api();
    free_host_memory(psi);

    return 0;
}
