#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "wdata.h"


// #include "pca_settings.h"
// #include "wslda_potdens.h"


// #include "pca_macro.h"

#include "wslda_toolkit.h"
#include "pca_logger.h"
#include "logger.h"
#include "predefines.h"
#include "gpe_utils.h"
// #include "pca_utils.h"
#include "gpe_engine_api.h"
#include "wslda_reproducibility.h"

int main( int argc , char ** argv ) 
{
    char execcmd[ 256 ];
    read_of_input_parameters(execcmd, argc, argv);
    int input_idx = parse_command_line_and_get_idx_of_input_file(argc, argv);
    read_input_file(input_idx, argv);
    
//     file_operation( check_if_can_overwrite_files() ); // terminate if file exists, and input->overwrite==0
//     sprintf(file_name, "%s_input.txt", md.outprefix);
//     file_operation( copy_input_file(argv[i],file_name) ); 
    
    double ekin, eint, eext, etot, etot_prev, time, rt, diff=-1;
    const int device=input->gpuDevice;
    const int mode=input->gpe_mode;
    const double alpha=input->alpha;
    const double beta=input->beta;
    const double npart=input->npart;
    const double dt=input->dt;
    const double time0 = 0.0;
    int inittype = input->inittype;
    if(mode==1) inittype = 5;
    
    set_gpu_device(device);

    int nx, ny, nz, ierr, it = 0;
    gpe_get_lattice_api(&nx, &ny, &nz);
    printf("# GPE engine compiled for lattice: %d x %d x %d\n", nx, ny, nz);

    switch (mode)
    {
    case 0:
        printf("# IMAGINARY TIME PROJECTION\n");
        break;
    case 1:
        printf("# REAL TIME EVOLUTION\n");
        break;
    default:
        break;
    }

    Complex *psi;
    double *density, *currents, *energy;
    wslda_potential nullPotential;
    uint nxyz=nx*ny*nz;
    alloc_host_memory(nxyz, &psi, &density, &currents);

    switch (inittype)
    {
    case 0:
        set_initial_wave_function( nxyz, psi);
        break;
    case 5:
        read_initial_wave_function( nxyz, psi);
        break;
    default:
        break;
    }

    gpe_create_engine_api(alpha, beta, dt, npart);
    gpe_set_user_params_api(MAX_USER_PARAMS, input->params);
    gpe_set_psi_api(time0, psi);
    if(mode==0) gpe_normalize_psi_api();

    print_header();
    gpe_energy_api(&time, &ekin, &eint, &eext);

    etot = ekin + eint + eext;
    print_intial_results(time, npart, etot, ekin, eint, eext);
    cppmallocl(energy, ENERGYITEMS, double);

    
    // Create empty WDATA set
    // use example:
    // https://gitlab.fizyka.pw.edu.pl/wtools/wdata/-/blob/master/c-examples/example-write.c
    // create metadata handler
    wdata_metadata wmd;

    // Lattice
    wmd.datadim = 3; // 3D data
    wmd.nx = NX;
    wmd.ny = NY;
    wmd.nz = NZ;
    wmd.dx = DX;
    wmd.dy = DY;
    wmd.dz = DZ;
    strcpy(wmd.prefix, input->outprefix);
    wmd.t0 = 0.0;
    wmd.dt = input->dt*input->timesteps;
    
    // add variables to data set
    // for each variable binary file of name `prefix_`varname`.wdat will be created
    wdata_variable vdensity_a = {"density_a", "real", "none", "wdat"};
    wdata_add_variable(&wmd, &vdensity_a);
    wdata_variable vdelta = {"psi", "complex", "none", "wdat"};
    wdata_add_variable(&wmd, &vdelta);
    wdata_variable vcurrent_a = {"current_a", "vector", "none", "wdat"};
    wdata_add_variable(&wmd, &vcurrent_a);

    // add links to data sets
    // links are alternative names of the same variable
    wdata_link ldensity_b = {"density_b", "density_a"};
    wdata_add_link(&wmd, &ldensity_b);
    wdata_link lcurrent_b = {"current_b", "current_a"};
    wdata_add_link(&wmd, &lcurrent_b);
    wdata_link ldelta = {"delta", "psi"};
    wdata_add_link(&wmd, &ldelta);

    // add constants
//     wdata_const lconst_eF = {"eF", 0.5, "MeV"};
//     wdata_add_const(&wmd, &lconst_eF);
// 
//     wdata_const lconst_kF = {"kF", 1.0, "1/fm"};
//     wdata_add_const(&wmd, &lconst_kF);

    // just in case - clear data sets
    // it removes binary files if they alredy exists
    wdata_clear_database(&wmd);
    gpe_get_psi_api(&time, psi);
    wdata_write_cycle(&wmd, "psi", psi);
    
    gpe_get_density(&time, density);
    wdata_write_cycle(&wmd, "density_a", density);

    gpe_get_currents(&time, currents);
    wdata_write_cycle(&wmd, "current_a", currents);

    wdata_add_cycle(&wmd);
    wdata_write_metadata_to_file(&wmd, "");
    
    logger_create_header(execcmd);

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
        
        // TODO: other energies???
        energy[EKIN] = ekin;


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
        gpe_get_psi_api(&time, psi);
        wdata_write_cycle(&wmd, "psi", psi);
        
        gpe_get_density(&time, density);
        wdata_write_cycle(&wmd, "density_a", density);

        gpe_get_currents(&time, currents);
        wdata_write_cycle(&wmd, "current_a", currents);
        
        wdata_add_cycle(&wmd);
        wdata_write_metadata_to_file(&wmd, "");
        
        logger_add_entry(it++, convert_into_wslda_density(density, nxyz), nullPotential, input->referencekF, NULL, energy, &npart, NULL, 0, NULL);
        

        if(mode==0 && fabs(diff) < input->energyconveps) break; // algorithm converged
        if(time > dt*input->timesteps*input->measurements) 
        {
            if(mode==0) printf("WARNING: ...\n");
            break; // do not allow to iterate infinitly long
        }
    }

    // close files
    
//     gpe_get_psi_api(&time, psi); //TOREMOVE
//     write_to_binary_file(nxyz, psi); //TOREMOVE
//     write_to_txt_file(nx, ny, nz, psi); //TOREMOVE

    gpe_destroy_engine_api();
    free_host_memory(psi);

    return 0;
}
