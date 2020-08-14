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
// #define MATRIX_IS_REAL

// Select diagonalization routine
// ELPA demonstrates the best performance, use it if target system supports this lib.
// Otherwise use standard ScaLapack lib (PZHEEV?) .
// In case of ScaLapack it is recommended to use PZHEEVR, unless this routine does not work correctly (it may happen on some systems)
// For more info see: http://git2.if.pw.edu.pl/gabrielw/cold-atoms/wikis/Parallelization-settings
#define DIAGONALIZATION_ROUTINE PZHEEVR
// #define DIAGONALIZATION_ROUTINE PZHEEVD
// #define DIAGONALIZATION_ROUTINE ELPA

// Maximal number of parameters in params array
#define MAX_USER_PARAMS 32 

// Minimal density to avoid numerical problems
// below this threshold density is regarded as zero
#define DENSEPSILON 1.0e-8


// ---------------------- ELPA SETTINGS ---------------------------
// Fill this part only if ELPA library is used for diagonnalization

// uncomment it if you want to activate GPU for diagonalizations 
#define ELPA_USE_GPU

// Select ELPA kernels
#define ELPS_USE_SOLVER ELPA_SOLVER_1STAGE
#define ELPA_USE_COMPLEX_KERNEL ELPA_2STAGE_COMPLEX_DEFAULT
#define ELPA_USE_REAL_KERNEL ELPA_2STAGE_REAL_DEFAULT

// Fraction of eigenvectors to be extracted in each cycle.
// 1.0 corresponds to extraction if all eigenvectors (USE IT IF YOU YOU ARE NOT SURE)
// NOTE: value of this parameter should assure that all eigenstates below requested Ec are extracted.  
// NOTE: For 3D case this value typically can be set to 0.78
#define ELPA_NEV_FRACTION 1.0
