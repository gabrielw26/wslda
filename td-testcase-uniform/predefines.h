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
 *      for simulating unitary Fermi gas,
 *      it is equivalent to selecting ASLDA functional with SLDA_FORCE_A1 option.
 *  - ASLDA:
 *      for simulating unitary Fermi gas,
 *      at qualitative level it produces results compatible with SLDA, however it is more accurate,
 *      due to presence of current terms in the functional, but it has worse convergence properties.
 *      For more info see: https://arxiv.org/abs/1008.3933
 *  - SLDAE:
 *      for simulating Fermi gas for an arbitrary value of akF,
 *      for small and negative akF the functional is compatible with BDG, while for large akF is compatible with ASLDA.
 *      For more info see: https://arxiv.org/abs/2201.07626
 *  - BDG:
 *      for simulating systems in BCS regime,
 *      equations of motion are equivalent to Bogoliubov-de-Gennes equations,
 *  - CUSTOMEDF:
 *      use this option to define your custom functional,
 *      then you need to provide body of functions: compute_energy_custom( ) and compute_potentials_custom( )
 *      in problem-definition.h file
 * */
#define FUNCTIONAL SLDA
// #define FUNCTIONAL ASLDA
// #define FUNCTIONAL SLDAE
// #define FUNCTIONAL BDG
// #define FUNCTIONAL CUSTOMEDF

/**
 * Sets the effective mass to be equal 1.
 * Computation converges faster, but it also decreases the accuracy of the functional.
 * If functional is SLDA or BDG, this option is activated automatically.
 * */
// #define SLDA_FORCE_A1

/**
 * Select which external potentials you want to use in simulations.
 *   - ENABLE_V_EXT:
 *       function v_ext(...) from problem definition will be called in each iteration.
 *   - ENABLE_DELTA_EXT:
 *       function delta_ext(...) from problem definition will be called in each iteration.
 *   - ENABLE_VELOCITY_EXT:
 *       function velocity_ext(...) from problem definition will be called in each iteration.
 *   - ENABLE_MODIFY_PDENSITIES:
 *       function modify_densities(...) from problem definition will be called in each iteration.
 *   - ENABLE_MODIFY_POTENTIALS:
 *       function modify_potentials(...) from problem definition will be called in each iteration.
 *   - ENABLE_MODIFY_ENERGIES:
 *       function modify_energies(...) from problem definition will be called in each iteration.
 * In order to achieve best performance disable call of empty functions.
 * */
#define ENABLE_V_EXT
// #define ENABLE_DELTA_EXT
// #define ENABLE_VELOCITY_EXT
// #define ENABLE_MODIFY_DENSITIES
// #define ENABLE_MODIFY_POTENTIALS
// #define ENABLE_MODIFY_ENERGIES

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
 * In such case, the code evolves only wave functions for single spin component
 * and in consequence computing time decreases by a factor of two.
 * */
// #define SPINSYMMETRY_MODE

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
 * Meaningful only in case of ASLDA and SLDAE.
 * Parameters defining stabilization procedure of ASLDA functional. 
 * For regions with density smaller than SLDA_STABILIZATION_EXCLUDE_BELOW_DENISTY 
 * contribution from current term j^2/2n is assumed to be zero. 
 * For regions with density above SLDA_STABILIZATION_RETAIN_ABOVE_DENSITY 
 * the contribution is assumed to be intact by stabilization procedure. 
 * For more info see: 
 * https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/Functionals#stabilization-of-aslda-functional
 * */
#define SLDA_STABILIZATION_RETAIN_ABOVE_DENSITY  1.0e-5
#define SLDA_STABILIZATION_EXCLUDE_BELOW_DENISTY 1.0e-7

/**
 * Active this flag in order to store quasi-particle energies for each measurement.
 * Note that in case of 1d or 2d codes it can require much more space than measurements itself.
 * */
// #define STORE_QPE

/**
 * Integration scheme ABxAMy, where
 *   x - order of predictor,
 *   y - order of corrector.
 * Choose one: AB3AM4, AB4AM5, AB5AM5
 * */
// #define INTEGRATION_SCHEME AB3AM4
#define INTEGRATION_SCHEME AB4AM5
// #define INTEGRATION_SCHEME AB5AM5

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

/**
 * Files: predefines.h, problem-definition.h, logger.h are assumed to be compatible with this API version
 * For list of API versions see: https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/API-version
 * */
#define API_VERSION 20240131

/**
 * activate this flag for setting code in testing mode with uniform system
 * */
#define UNIFORM_TEST_MODE
