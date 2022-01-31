/**
 * Define lattice size and lattice spacing
 * */
#define NX 12
#define NY 14
#define NZ 16

#define DX 0.9
#define DY 0.9
#define DZ 0.9

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
 *      you MUST set aBdG value in input file when using this functional 
 *  - CUSTOMEDF:
 *      use this option to define your custom functional,
 *      then you need to provide body of functions: compute_energy_custom( ) and compute_potentials_custom( )
 *      in problem-definition.h file
 * */
#define FUNCTIONAL SLDAE
// #define FUNCTIONAL ASLDA
// #define FUNCTIONAL BDG
// #define FUNCTIONAL CUSTOMEDF

/**
 * Meaningful only in case SLDAE.
 * Sets effective mass to be equal, and speeds-up computation (approximately by a factor of two)
 * */
#define SLDAE_FORCE_A1

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
// #define REGULARIZATION_SCHEME SPHERICAL_CUTOFF
// #define REGULARIZATION_SCHEME CUBIC_CUTOFF

/**
 * Select diagonalization routine
 * ELPA demonstrates the best performance, use it if target system supports this lib.
 * Otherwise use standard ScaLapack lib (PZHEEV?) .
 * In case of ScaLapack it is recommended to use PZHEEVR, unless this routine does not work correctly (it may happen on some systems)
 * For more info see: http://git2.if.pw.edu.pl/gabrielw/cold-atoms/wikis/Parallelization-settings
 * */
// #define DIAGONALIZATION_ROUTINE PZHEEVR
// #define DIAGONALIZATION_ROUTINE PZHEEVD
#define DIAGONALIZATION_ROUTINE ELPA

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
 * activate this flag for setting code in testing mode with uniform system
 * */
#define UNIFORM_TEST_MODE

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
// #define ELPA_USE_GPU

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
