#ifdef WSLDA
/**
 * Settings for static codes.
 * */

/**
 * Select diagonalization routine
 * ELPA demonstrates the best performance, use it if target system supports this lib.
 * Otherwise use standard ScaLapack lib (PZHEEV?) .
 * In case of ScaLapack it is recommended to use PZHEEVR, unless this routine does not work correctly (it may happen on some systems)
 * For more info see: https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Setting%20up%20diagonalization%20engine
 * */
// #define DIAGONALIZATION_ROUTINE PZHEEVR
// #define DIAGONALIZATION_ROUTINE PZHEEVD
#define DIAGONALIZATION_ROUTINE ELPA

/**
 * ---------------------- ELPA SETTINGS ---------------------------
 * Fill this part only if ELPA library is used for diagonalization
 * 
 * Default settings are: ELPA_SOLVER_1STAGE
 * but you can overwrite using options below
 * */

/**
 * uncomment it if you want to activate GPUs for diagonalizations 
 * */
#define ELPA_USE_GPU

/**
 * Select ELPA kernels,
 * for more info see documentation of ELPA lib
 * */
// #define ELPA_USE_SOLVER ELPA_SOLVER_2STAGE
// #define ELPA_USE_COMPLEX_KERNEL ELPA_2STAGE_COMPLEX_GPU
// #define ELPA_USE_REAL_KERNEL ELPA_2STAGE_REAL_GPU

/**
 * Fraction of eigenvectors to be extracted in each cycle.
 * 1.0 corresponds to extraction of all eigenvectors (default)
 * NOTE: value of this parameter should assure that all eigenstates below requested Ec are extracted.  
 * NOTE: For 3D case this value typically can be set to 0.78, for 1D and 2D casese 1.0 is recommended.
 * */
// #define ELPA_NEV_FRACTION 1.0

#endif

#ifdef TDWSLDA
/**
 * Settings for time-dependent codes.
 * Time-dependent codes require GPUs. 
 * */

/**
 * Number of mpi processes per IO group used for collective (parallel) writing of checkpoint files.
 * Performance of read/write checkpoint depends on the number of writes involved in IO process,
 * and optimal value depends on the computer. 
 * Use the default value (24) unless you are not satisfied with IO performance. 
 * */
#define MPI_NP_PER_IO_GROUP 24

/**
 * Default number of GPUs per node.
 * You can overwrite this value by using gpuspernode tag in the input file.
 * The flag is ingored if CUSTOM_GPU_DISTRIBUTION is selected.
 * */
#define GPUS_PER_NODE 8

/**
 * Activate this flag in order to print to stdout
 * applied mapping mpi-process <==> device-id.
 * */
#define PRINT_GPU_DISTRIBUTION

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

    return 0;
}
#endif

#endif

