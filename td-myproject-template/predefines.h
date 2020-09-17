/**
 * Define lattice size and lattice spacing.
 * NOTE: presently only DX=DY=DZ=1 is implemented. 
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
 * Active this flag in order to store quasi-particle energies for each measurement.
 * Note that in case of 1d or 2d codes it can require much more space than measurements itself.
 * */
// #define STORE_QPE
