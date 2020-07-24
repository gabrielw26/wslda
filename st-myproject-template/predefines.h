// Lattice
#define NX 8
#define NY 10
#define NZ 12

#define DX 1.0
#define DY 1.0                                                                                                                                       
#define DZ 1.0

// #define FUNCTIONAL SLDA
#define FUNCTIONAL ASLDA
// #define FUNCTIONAL BDG

// activate this if you know that HFB matrix is real
// the code will utilize it in roder to speed-up the calculations
// meaningful only for static codes
// #define MATRIX_IS_REAL

// spin-symmetric mode decreases computing time for factor about two
// #define SPINSYMMETRY_MODE

// select diagonalization routine
// it is recommended to use PZHEEVR, unless this routine does not work correctly (it may happen on some systems)
#define DIAGONALIZATION_ROUTINE PZHEEVR
// #define DIAGONALIZATION_ROUTINE PZHEEVD

// Maximal number of parameters in params array
#define MAX_USER_PARAMS 32 

// Minimal density to avoid numerical problems
// below this threshold density is regarded as zero
#define DENSEPSILON 1.0e-8

