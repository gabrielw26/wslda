/**
 * Define lattice size and lattice spacing.
 * */
#define NX 8
#define NY 10
#define NZ 12

#define DX 1.0
#define DY 1.0
#define DZ 1.0

/**
 * Select functional:
 *  - SLDA: 
 *      for simulating unitary Fermi gas, 
 *      it sets effective mass of particles to 1.0 which assures better convergence properties,
 *      in case of time time-dependent calculations SLDA is about 2x faster than ASLDA.
 *  - ASLDA:
 *      for simulating unitary Fermi gas,
 *      at qualitative level it produces results compatible with SLDA, however it is more accurate,
 *      due to presence of current terms in the functional it has worse convergence properties.
 *  - BDG:
 *      for simulating systems in BCS regime,
 *      equations of motion are equivalent to Bogoliubov-de-Gennes equations,
 *      you MUST set aBdG value in input file when using this functional.
 * */
// #define FUNCTIONAL SLDA
#define FUNCTIONAL ASLDA
// #define FUNCTIONAL BDG

/**
 * Select which external potentials you want to use in simulations.
 *   - ENABLE_V_EXT:
 *       function v_ext(...) from problem definition will be called in each iteration.
 *   - ENABLE_DELTA_EXT:
 *       function delta_ext(...) from problem definition will be called in each iteration.
 *   - ENABLE_VELOCITY_EXT:
 *       function velocity_ext(...) from problem definition will be called in each iteration.
 * In order to achieve best performance disable call of empty functions.
 * */
#define ENABLE_V_EXT
// #define ENABLE_DELTA_EXT
// #define ENABLE_VELOCITY_EXT

/**
 * Maximal number of parameters in params array
 * */
#define MAX_USER_PARAMS 32 

/**
 * Minimal density to avoid numerical problems
 * below this treshold density is regarded as zero
 * */
#define DENSEPSILON 1.0e-8

/**
 * Active this flag in case of calculations for problems that preserve symmetry between spin up and down components.
 * In such case, the evolves only wave-functions for single spin component
 * and in consequence computing time decreases by factor of two
 * */
// #define SPINSYMMETRY_MODE

/**
 * Meaningful only in case of ASLDA.
 * Parameters defining stabilization procedure of ASLDA functional. 
 * For regions with density smaller than ASLDA_STABILIZATION_EXCLUDE_BELOW_DENISTY 
 * contribution from current term j^2/2n is assumed to be zero. 
 * For regions with density above ASLDA_STABILIZATION_RETAIN_ABOVE_DENSITY 
 * the contribution is assumed to be intact by stabilization procedure. 
 * For more info see: 
 * https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/Functionals#stabilization-of-aslda-functional
 * */
#define ASLDA_STABILIZATION_RETAIN_ABOVE_DENSITY  1.0e-5
#define ASLDA_STABILIZATION_EXCLUDE_BELOW_DENISTY 1.0e-7

/**
 * Number of mpi processes per IO group used for collective (parallel) writing of checkpoint files.
 * Performance of read/write checkpoint depends on the number of writes involved in IO process,
 * and optimal value depends on the computer. 
 * Use the default value (24) unless you are not satisfied with IO performance. 
 * */
#define MPI_NP_PER_IO_GROUP 24

/**
 * Active this flag in order to store quasi-particle energies for each measurement.
 * Note that in case of 1d or 2d codes it can require much more space than measurements itself.
 * */
// #define STORE_QPE

/**
 * Activate this flag in order to print to stdout
 * applied mapping mpi-process <==> device-id.
 * */
// #define PRINT_GPU_DISTRIBUTION

/**
 * Activate this flag if target machine has non-standard distribution of GPUs. 
 * In such case you need to provide body of function `assign_deviceid_to_mpi_process`.
 * If this flag is commented-out it is assumed that code is running on a machine 
 * with uniformly distributed GPU cards across the nodes, 
 * and each node has `gpuspernode` (input file parameter) cards.
 * */
// #define CUSTOM_GPU_DISTRIBUTION

/**
 * This function is used to assign unique device-id to mpi process.
 * @param comm MPI communicator
 * @return device-id assign to the process extracted by function MPI_Comm_rank(...)
 * DO NOT REMOVE STATEMENT `#if ... BELOW !!!
 * */
#if defined(CUSTOM_GPU_DISTRIBUTION) && defined(TDWSLDA_MAIN)
int assign_deviceid_to_mpi_process(MPI_Comm comm)
{
    int np, ip;
    MPI_Comm_size(comm, &np);
    MPI_Comm_rank(comm, &ip);
    
    // assign here deviceid to process with ip=iam
    int deviceid=0;
    
    return deviceid;
}
#endif
