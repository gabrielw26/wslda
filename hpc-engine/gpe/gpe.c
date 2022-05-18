#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "wdata.h"

#include "wslda_functionals.h"

#include "pca_settings.h"
#include "wslda_potdens.h"

#include "pca_utils.h"
#include "pca_logger.h"

double dc_ec;

#include "predefines.h"
#include "gpe_utils.h"

#include "gpe_engine_api.h"
#include "wslda_reproducibility.h"

int wsldapid; // process id - global variable
int wsldapnp; // total number of processes - global variable
#define printf wprintf
#include "logger.h"
#undef printf


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
    double npart=input->npart;
    const double dt=input->dt;
    const double time0 = 0.0;
    int inittype = input->inittype;
    if(mode==1) inittype = 5;
    
    set_gpu_device(device);

    int nx, ny, nz, ierr, it = 0, i;
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

    wslda_density densall;
    densall.nx=NX; densall.ny=NY; densall.nz=NZ; 
    densall.datadim=3;
    densall.blocklength=NX*NY*NZ;
    densall.nu=(double_complex *)psi;
    densall.rho_a=density;
    densall.rho_b=density;
    densall.j_a_x=currents; densall.j_a_y=currents+nxyz; densall.j_a_z=currents+nxyz*2;
    densall.j_b_x=currents; densall.j_b_y=currents+nxyz; densall.j_b_z=currents+nxyz*2;
    
    gpe_create_engine_api(alpha, beta, dt, npart);
    gpe_set_user_params_api(MAX_USER_PARAMS, input->params);
    gpe_set_sclgth_api(input->sclgth);
    gpe_set_psi_api(time0, psi);
    if(mode==0) gpe_normalize_psi_api();

    print_header();
    gpe_energy_api(&time, &ekin, &eint, &eext);

    etot = ekin + eint + eext;
    print_intial_results(time, npart, etot, ekin, eint, eext);
    cppmallocl(energy, ENERGYITEMS, double);

    
//     // Load data
//     if(ip==0) extra_data_size = get_extra_data_size(md.params);
//     MPI_Bcast( &extra_data_size , sizeof(size_t) , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
//     if(extra_data_size>0)
//     {
//         if(ip==0) wprintf("# EXTRA_DATA IS ACTIVE.\n");
//         if(ip==0) wprintf("# ALLOCATING EXTRA_DATA OF SIZE %ld B.\n", extra_data_size); fflush(stdout);
//         if ( ( extra_data = (void *) malloc( extra_data_size ) ) == NULL  )
//         {                                                             
//             wfprintf( stderr , "error: cannot malloc()! Exiting!\n") ; 
//             wfprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; 
//             MPI_Finalize() ;
//             /* Arrays will be cleared automatically */
//             return( EXIT_FAILURE ) ; 
//         }
//         if(ip==0) wprintf("# EXECUTING: load_extra_data(%zu, extra_data, input->params)\n", extra_data_size);
//         if(ip==0) cpu_exec( load_extra_data(extra_data_size, extra_data, md.params) );
//         MPI_Bcast( extra_data , extra_data_size , MPI_BYTE , 0 , MPI_COMM_WORLD ) ;
//         
//         // copy extra data to GPU
//         gpu_exec( gpu_malloc(extra_data_size, (void **)&d_extra_data) );
//         gpu_exec( memcopy_host2gpu(extra_data, d_extra_data,  extra_data_size) ); 
//         gpu_exec( memcopy_extra_data(extra_data_size, d_extra_data) );
//         
//         // reproducibility pack
//         if(ip==0) save_extradata_to_file(extra_data_size, extra_data);
//     }
//     wprintf("# EXECUTING: process_params(input->params, [%f], NULL, %zu, extra_data)\n", kF, extra_data_size);
//     process_params(md.params, &kF, NULL, extra_data_size, extra_data);
    
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
    
    gpe_get_density(&time, density); for(i=0; i<nxyz; i++)  density[i]*=0.5;
    wdata_write_cycle(&wmd, "density_a", density);

    gpe_get_currents(&time, currents); for(i=0; i<nxyz*3; i++)  currents[i]*=0.5;
    wdata_write_cycle(&wmd, "current_a", currents);

    wdata_add_cycle(&wmd);
    wdata_write_metadata_to_file(&wmd, "");
    
    logger_create_header(execcmd);

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
        
        gpe_get_density(&time, density); for(i=0; i<nxyz; i++)  density[i]*=0.5;
        wdata_write_cycle(&wmd, "density_a", density);

        gpe_get_currents(&time, currents); for(i=0; i<nxyz*3; i++)  currents[i]*=0.5;
        wdata_write_cycle(&wmd, "current_a", currents);
        
        wdata_add_cycle(&wmd);
        wdata_write_metadata_to_file(&wmd, "");
        
        logger_add_entry(it++, densall, nullPotential, input->referencekF, NULL, energy, &npart, NULL, 0, NULL);
        

        if(mode==0 && fabs(diff) < input->energyconveps) break; // algorithm converged
        if(time > dt*input->timesteps*input->measurements) 
        {
            if(mode==0) printf("WARNING: Program has executed %d steps and still doesn't converge.\n", input->measurements);
            break; // do not allow to iterate infinitly long
        }
    }

    // close files

    gpe_destroy_engine_api();
    free_host_memory(psi);

    return 0;
}
