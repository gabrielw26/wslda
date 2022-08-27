#ifndef __TDWSLDA_ENERGY_GENERIC__
#define __TDWSLDA_ENERGY_GENERIC__

// DO NOT REMOVE!
#include "tdwslda_functionals_framework_enable.h"

/**
 * Function computes energy, particle number and angular momentum
 * NOTE: this function executes cudaMemcpy, thus it is BLOCKING!
 * @param it index of time step, it is ised for proper evaluation of external potential
 * @param d_densites (INPUT)
 *                   collective array with densities
 * @param d_potentials (INPUT)
 *                     collective array with potentials
 * @param d_workarea (INPUT/OUTPUT)
 *                   working buffer of size PCA_WORKSPACE_SHIFT*NUMBER_ELEMENT*sizeof(double complex),
 *                   On OUTPUT first TDWSLDAITEMS elements contain energies and particle number, as specifie in pca_settings.h
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 * */
extern "C" int compute_energy_generic(int it, double *d_densities, double *d_potentials, double *d_workarea, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NUMBER_ELEMENT/nthreads);
    wslda_density densall=convert_into_wslda_density(d_densities, NUMBER_ELEMENT);
    wslda_potential potsall=convert_into_wslda_potential(d_potentials, NUMBER_ELEMENT, NULL);

    // Step 1: prepare buffers for local reductions
    // intrinsic energy
    double *E_kin = (double *)(d_workarea +  EKIN*NUMBER_ELEMENT);
    double *E_pot = (double *)(d_workarea +  EPOT*NUMBER_ELEMENT);
    double *E_pair= (double *)(d_workarea +  EPAIR*NUMBER_ELEMENT);
    double *E_Curr= (double *)(d_workarea +  ECURRENT*NUMBER_ELEMENT);
    tdwslda_compute_energy<<<nblocks, nthreads>>>(it, densall, potsall,
                                                 E_kin, E_pot, E_pair, E_Curr
                                                 );
    
    // external energy - potential
    double *E_ext = (double *)(d_workarea +  EPOTEXT*NUMBER_ELEMENT);
#ifdef ENABLE_V_EXT
    kernel_compute_energy_v_ext<<<nblocks, nthreads>>>(it,
                                                    densall.rho_a, densall.rho_b, densall.nu,
                                                    densall.j_a_x, densall.j_a_y, densall.j_a_z, densall.j_b_x, densall.j_b_y,densall. j_b_z,
                                                    densall.tau_a, densall.tau_b,
                                                    potsall.delta,
                                                    E_ext
                                                    );
#else
    cuda_set_array_elements(NUMBER_ELEMENT, E_ext,  0.0, nthreads); // set this contribution to 0.0
#endif
    
    // external energy - pairing
    double *E_ext_pair = (double *)(d_workarea +  EPAIREXT*NUMBER_ELEMENT);
#ifdef ENABLE_DELTA_EXT
    kernel_compute_energy_delta_ext<<<nblocks, nthreads>>>(it,
                                                    densall.rho_a, densall.rho_b, densall.nu,
                                                    densall.j_a_x, densall.j_a_y, densall.j_a_z, densall.j_b_x, densall.j_b_y, densall.j_b_z,
                                                    densall.tau_a, densall.tau_b,
                                                    potsall.delta,
                                                    E_ext_pair
                                                    );
#else
    cuda_set_array_elements(NUMBER_ELEMENT, E_ext_pair,  0.0, nthreads); // set this contribution to 0.0
#endif

    // external energy - velocity
    double *E_ext_vel = (double *)(d_workarea +  EVELEXT*NUMBER_ELEMENT);
#ifdef ENABLE_VELOCITY_EXT
    kernel_compute_energy_velocity_ext<<<nblocks, nthreads>>>(it,
                                                    densall.rho_a, densall.rho_b, densall.nu,
                                                    densall.j_a_x, densall.j_a_y, densall.j_a_z, densall.j_b_x, densall.j_b_y, densall.j_b_z,
                                                    densall.tau_a, densall.tau_b,
                                                    potsall.delta,
                                                    E_ext_vel
                                                    );
#else
    cuda_set_array_elements(NUMBER_ELEMENT, E_ext_vel,  0.0, nthreads); // set this contribution to 0.0
#endif

    // partile number
    double *Na = (double *)(d_workarea +  NPARTA*NUMBER_ELEMENT);
    double *Nb = (double *)(d_workarea +  NPARTB*NUMBER_ELEMENT);
    if( cudaMemcpy(Na, densall.rho_a, sizeof(double)*NUMBER_ELEMENT, cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -100;
    if( cudaMemcpy(Nb, densall.rho_b, sizeof(double)*NUMBER_ELEMENT, cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -100;
    
    // angular momentum
    double *La = (double *)(d_workarea +  LZA*NUMBER_ELEMENT);
    double *Lb = (double *)(d_workarea +  LZB*NUMBER_ELEMENT);
#if CODEDIM==1
    // no angular momentum for 1D case
    cuda_set_array_elements(NUMBER_ELEMENT, La,  0.0, nthreads); // set this contribution to 0.0
    cuda_set_array_elements(NUMBER_ELEMENT, Lb,  0.0, nthreads); // set this contribution to 0.0
#else
    kernel_compute_angular_momentum_z<<<nblocks, nthreads>>>(densall.j_a_x, densall.j_a_y, densall.j_a_z, 
                                                             densall.j_b_x, densall.j_b_y, densall.j_b_z, 
                                                             La, Lb);
#endif

    // Step 2: massive reductions
    int ierr = local_reductions_many(TDWSLDAITEMS, NUMBER_ELEMENT, d_workarea, d_workarea);
    if(ierr!=0) return ierr;
    
    // add missing volume element
    cuda_scale_array_elements(2, d_workarea +  NPARTA,  VOLUME_ELEMENT, 1); // add missing volume element
    
    return 0;
}


// DO NOT REMOVE!
#include "tdwslda_functionals_framework_disable.h"

#endif
