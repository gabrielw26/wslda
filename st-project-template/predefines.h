/**
 * Define lattice size and lattice spacing
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
 *  - SLDAE:
 *      for simulating Fermi gas for an arbitrary value of akF,
 *      for small and negative akF the functional is equivalent to BDG, while for large akF is equivalent to ASLDA.
 *  - BDG:
 *      for simulating systems in BCS regime,
 *      equations of motion are equivalent to Bogoliubov-de-Gennes equations,
 *      you MUST set aBdG value in input file when using this functional 
 *  - CUSTOMEDF:
 *      use this option to define your custom functional,
 *      then you need to provide body of functions: compute_energy_custom( ) and compute_potentials_custom( )
 *      in problem-definition.h file
 * */
// #define FUNCTIONAL SLDA
#define FUNCTIONAL ASLDA
// #define FUNCTIONAL SLDAE
// #define FUNCTIONAL BDG
// #define FUNCTIONAL CUSTOMEDF

/**
 * activate this if you know that Hamiltonian matrix is real, 
 * the code will utilize it in order to speed-up the calculations by factor 4x (approximately)
 * */
// #define MATRIX_IS_REAL

/**
 * Scheme of pairing field renormalization procedure. 
 * For more info see: https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/Regularization%20schemes%20of%20the%20pairing%20field
 * Select one:
 * SPHERICAL_CUTOFF: use spherical momentum space cutoff, in this case you need to set `ec` variable in input file (default).
 * CUBIC_CUTOFF: use cubic momentum space cutoff, in this case `ec` will be set to infinity automatically.
 * */
#define REGULARIZATION_SCHEME SPHERICAL_CUTOFF
// #define REGULARIZATION_SCHEME CUBIC_CUTOFF

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
 * Machine file. 
 * This file contains info about machine that will be used in the computation process.
 * You can specify the file in the following ways: 
 * - copy `machine.h` file to the current directory, see templates folder for various examples,
 * - specify the folder with `machine.h` file via -I option in Makefile  
 * - use system variable WSLDA_MACHINE to specify the folder with `machine.h` file, for example:
 *   export WSLDA_MACHINE=...
 **/ 
#include "machine.h"
