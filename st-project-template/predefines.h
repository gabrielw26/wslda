/**
 * Define lattice size and lattice spacing.
 * Keep them inside parentheses.
 * */
#define NX (8)
#define NY (10)
#define NZ (12)

#define DX (1.0)
#define DY (1.0)
#define DZ (1.0)

/**
 * Select functional:
 *  - SLDA: 
 *      For simulating a unitary Fermi gas. 
 *      It is equivalent to selecting ASLDA functional with SLDA_FORCE_A1 option.
 *  - ASLDA:
 *      For simulating a unitary Fermi gas.
 *      At the qualitative level, it produces results compatible with SLDA; however, it is more accurate.
 *      Due to the presence of current terms in the functional, it has worse convergence properties.
 *      For more info, see: https://arxiv.org/abs/1008.3933
 *  - SLDAE:
 *      For simulating a Fermi gas for an arbitrary value of akF.
 *      For small and negative akF, the functional is compatible with BDG, while for large akF, it is compatible with ASLDA.
 *      For more info, see: https://arxiv.org/abs/2201.07626
 *  - BDG:
 *      for simulating systems in the BCS regime,
 *      equations of motion are equivalent to Bogoliubov-de-Gennes equations,
 *  - CUSTOMEDF:
 *      Use this option to define your custom functional.
 *      Then you need to provide the body of functions: compute_energy_custom( ) and compute_potentials_custom( )
 *      in problem-definition.h file
 * */
#define FUNCTIONAL SLDA
// #define FUNCTIONAL ASLDA
// #define FUNCTIONAL SLDAE
// #define FUNCTIONAL BDG
// #define FUNCTIONAL CUSTOMEDF

/**
 * Sets the effective mass to be equal to 1.
 * Computation converges faster, but it also decreases the accuracy of the functional.
 * If functional is SLDA or BDG, this option is activated automatically.
 * */
// #define SLDA_FORCE_A1

/**
 * Activate this if you know that the Hamiltonian matrix is real. 
 * The code will utilize it in order to speed up the calculations by a factor of 4x (approximately)
 * */
// #define HAMILTONIAN_IS_REAL

/**
 * Scheme of pairing field renormalization procedure. 
 * For more info see: Wiki -> Regularization schemes of the pairing field
 * Select one:
 * SPHERICAL_CUTOFF: use spherical momentum space cutoff; in this case, you need to set the `ec` variable in the input file (default).
 * CUBIC_CUTOFF: use cubic momentum space cutoff; in this case, `ec` will be set to infinity automatically.
 * */
#define REGULARIZATION_SCHEME SPHERICAL_CUTOFF
// #define REGULARIZATION_SCHEME CUBIC_CUTOFF

/**
 * Maximal number of parameters in params array
 * */
#define MAX_USER_PARAMS 32 

/**
 * Minimal density to avoid numerical problems.
 * Below this threshold density is regarded as zero
 * */
#define DENSEPSILON 1.0e-8

/**
 * Meaningful only in the case of ASLDA and SLDAE.
 * Parameters defining the stabilization procedure of the ASLDA functional. 
 * For regions with density smaller than SLDA_STABILIZATION_EXCLUDE_BELOW_DENISTY, 
 * contribution from the current term j^2/2n is assumed to be zero. 
 * For regions with density above SLDA_STABILIZATION_RETAIN_ABOVE_DENSITY, 
 * the contribution is assumed to be intact by the stabilization procedure. 
 * For more info see: Wiki -> Stabilization of ASLDA functional
 * */
#define SLDA_STABILIZATION_RETAIN_ABOVE_DENSITY  1.0e-5
#define SLDA_STABILIZATION_EXCLUDE_BELOW_DENISTY 1.0e-7


/**
 * Machine file. 
 * This file contains information about the machine that will be used in the computation process.
 * You can specify the file in the following ways: 
 * - copy `machine.h` file to the current directory, see templates folder for various examples,
 * - specify the folder with the `machine.h` file via -I option in Makefile  
 * - use system variable WSLDA_MACHINE to specify the folder with the `machine.h` file, for example:
 *   export WSLDA_MACHINE=...
 **/  
#include "machine.h"

/**
 * Files: predefines.h, problem-definition.h, logger.h are assumed to be compatible with this API version
 * For a list of API versions, see Wiki -> API version
 * */
#define API_VERSION 20240131

