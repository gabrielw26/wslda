// Lattice
#define NX 8
#define NY 10
#define NZ 12

// Lattice spacing
// NOTE: presently only DX=DY=DZ=1 is implemented. 
#define DX 1.0
#define DY 1.0                                                                                                                                       
#define DZ 1.0


// #define FUNCTIONAL SLDA
#define FUNCTIONAL ASLDA
// #define FUNCTIONAL BDG

// Enable computation with extarnal potential field
// If this flag is active, you must provide body of delta_ext(...) function in pca_uext.h file
// #define ENABLE_V_EXT

// Enable computation with extarnal delta field
// If this flag is active, you must provide body of delta_ext(...) function in pca_uext.h file
// #define ENABLE_DELTA_EXT

// Enable computation with extarnal velocity field
// If this flag is active, you must provide body of vector_vext(...) function in pca_uext.h file
// #define ENABLE_VELOCITY_EXT

// Maximal number of parameters in params array
#define MAX_USER_PARAMS 32 

// Minimal density to avoid numerical problems
// below this treshold density is regarded as zero
#define DENSEPSILON 1.0e-8

// Active it when you do calculations for system that preserves symmetry between spin components
// spin-symmetric mode decreases computing time for factor about two
// #define SPINSYMMETRY_MODE

// // To switch to cubic cut-off mode
// #define USE_CUBIC_CUTOFF

// active this flag in order to store quasi-particle energies for each measurement
// note that in case of 1d or 2d codes this can require much more space than measurements itself
// meaningful only for dynamic codes
// #define STORE_QPE

// setting code in testing mode with uniform system
#define UNIFORM_TEST_MODE
