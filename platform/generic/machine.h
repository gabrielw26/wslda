#ifdef WSLDA
/**
 * Settings for static codes.
 * */

/**
 * Select diagonalization routine
 * ELPA demonstrates the best performance; use it if the target system supports this library.
 * Otherwise use standard ScaLapack lib (PZHEEV?).
 * In case of ScaLapack, it is recommended to use PZHEEVR, unless this routine does not work correctly (it may happen on some systems)
 * For more info, see: Wiki -> Setting up diagonalization engine
 * */
#define DIAGONALIZATION_ROUTINE PZHEEVR
// #define DIAGONALIZATION_ROUTINE PZHEEVD
// #define DIAGONALIZATION_ROUTINE ELPA

/**
 * ---------------------- ELPA SETTINGS ---------------------------
 * Fill this part only if the ELPA library is used for diagonalization
 * 
 * Default settings are: ELPA_SOLVER_1STAGE
 * but you can overwrite using the options below
 * */

/**
 * Uncomment it if you want to activate GPUs for diagonalizations 
 * */
// #define ELPA_USE_GPU

/**
 * Select ELPA kernels,
 * for more info, see the  documentation of the ELPA lib
 * */
// #define ELPA_USE_SOLVER ELPA_SOLVER_2STAGE
// #define ELPA_USE_COMPLEX_KERNEL ELPA_2STAGE_COMPLEX_GPU
// #define ELPA_USE_REAL_KERNEL ELPA_2STAGE_REAL_GPU

/**
 * Fraction of eigenvectors to be extracted in each cycle.
 * 1.0 corresponds to the extraction of all eigenvectors (default)
 * NOTE: The value of this parameter should ensure that all eigenstates below the requested Ec are extracted.  
 * NOTE: For 3D cases, this value typically can be set to 0.78; for 1D and 2D cases, 1.0 is recommended.
 * */
// #define ELPA_NEV_FRACTION 1.0

#endif

#ifdef TDWSLDA
/**
 * Settings for time-dependent codes.
 * Time-dependent codes require GPUs. 
 * */

/**
 * Use this option if the machine has a GPU-aware MPI implementation,
 * i.e, the machine can handle a buffer regardless of whether it resides in host or device memory.
 * Usage of GPU-aware MPI can boost the performance of the computation.
 * */
// #define USE_GPU_AWARE_MPI

/**
 * Number of MPI processes per IO group used for collective (parallel) writing of checkpoint files.
 * Performance of read/write checkpoint depends on the number of writes involved in the IO process,
 * and the optimal value depends on the computer. 
 * Use the default value (24) unless you are not satisfied with the IO performance. 
 * */
#define MPI_NP_PER_IO_GROUP 24

/**
 * Default number of GPUs per node.
 * You can overwrite this value by using the gpuspernode tag in the input file.
 * The flag is ignored if CUSTOM_GPU_DISTRIBUTION is selected.
 * */
#define GPUS_PER_NODE 1

/**
 * Activate this flag in order to print to stdout
 * applied mapping mpi-process <==> device-id.
 * */
// #define PRINT_GPU_DISTRIBUTION

/**
 * Activate this flag if the target machine has a non-standard distribution of GPUs. 
 * In such a case, you need to provide the body of the function `assign_deviceid_to_mpi_process`.
 * If this flag is commented out, it is assumed that the code is running on a machine 
 * with uniformly distributed GPU cards across the nodes, 
 * controlled by GPUS_PER_NODE or `gpuspernode` input file tag.
 * */
// #define CUSTOM_GPU_DISTRIBUTION

/**
 * This function is used to assign a unique device ID to the MPI process.
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

#endif
