// Author: Gabriel Wlazlowski

#ifndef __PCA_SETTINGS__
#define __PCA_SETTINGS__

// ===============================================================================
// ================================ GLOBAL SETTINGS ==============================
// ===============================================================================

#define CODE PCA_ASLDA
#define VERSION "1.10"

// Lattice

// Lattice
#define NX 40
#define NY 40
#define NZ 40

#define DX 1.0
#define DY 1.0                                                                                                                                       
#define DZ 1.0

// #define NX 96
// #define NY 64
// #define NZ 48
// 
// // Lattice spacing - NOTE only supported in s3dpca code!!!
// #define DX (96./NX)
// #define DY (64./NY)                                                                                                                                       
// #define DZ (48./NZ)


#define DXYZ (DX*DY*DZ)

// Volume settings
#define LX (DX*NX)
#define LY (DY*NY)
#define LZ (DZ*NZ)
#define LXYZ (LX*LY*LZ)

#define NXYZ (NX*NY*NZ)
 
// Settings for 2D calculations
#define NXY (NX*NY)
#define LXY (LX*LY)

#define SPINA 0
#define SPINB 1

// Integration scheme AB - predictor, AM - correctior, number specify order
// #define ITEGRATION_SCHEME AB3AM4
#define INTEGRATION_SCHEME AB4AM5

// Maximal number of parameters in params array
#define MAX_USER_PARAMS 32 

// Number of self-consistent iterations for U and delta computation
#define UD_SCITERS 15
// Mixing parameter for self-consistent algorithm - fraction of new solution used for mixing
#define UD_MIX_COEFF 0.75

// Minimal density to avoid numerical problems
#define DENSEPSILON 1.0e-8

// if particle number changed by this percentage then break the simulation
#define N_STABILITY_CRITERIA 0.25 

// spin-symmetric mode decreases computing time for factor about two
#define SPINSYMMETRY_MODE

// if spin-symmetric system is assumed (spinsymmetry==1) then only states where fbeta(ek,beta)>SPINSYMMETRY_CUTOFF are considered
#define SPINSYMMETRY_CUTOFF 1.0e-18

// // To switch to cubic cut-off mode
// #define USE_CUBIC_CUTOFF

// activate this flag for setting code in testing mode with uniform system
// #define UNIFORM_TEST_MODE

// active rotating frame framework
// NOTE: in rotating frame the code uses constant variables dc_Omega_a and dc_Omega_b
// NOTE: the systems rotates along z-axis
// #define WORK_IN_ROTATING_FRAME

// activate parallel method for local reductions
// #define PCA_REDUCE_MANY
// if yes fill correctly these values
#define blockSize_d 256
#define threads_red_d 1024

#define HOWMANY 16
#define streams_d 8

// ===================================================================================
// =========================== PARAMETERS OF EDF =====================================
// ===================================================================================

// if BDG_MODE then BdG functional is activated and aBdG parameter is active in dynamical codes
// #define BDG_MODE

// // effective mass - not 1.0 then current corrections are needed!
// #define CURRENT_CORRECTIONS
// #define A0 1.094
// #define A1 0.156
// #define A2 -0.532

// effective mass - equal 1.0 then no current corrections
// activate this flag to skip computation of gradients of wf - significant spped up
#define FAST_CONST_EFFECTIVE_MASS_MODE
#define A0 1.000
#define A1 0.0
#define A2 0.0

// normal part
#define G0 0.357
#define G1 0.642

// pairing
#define GAMMA0 -11.11
// #define GAMMA0 -0.000001

// regularization function parameters
#define P_NMIN 1.0e-7
#define P_NMAX 1.0e-5
#define P_ALPHA 1.0

// #define P_NMIN 1.0e-8
// #define P_NMAX 1.0e-6
// #define P_ALPHA 1.0

// ===================================================================================
// =================================== TECHNICAL =====================================
// ===================================================================================

// cufft plans
#define CUFFT_NUMBER_OF_PLANS 4
#define PLAN_Z2Z_BATCH 0
#define PLAN_Z2Z_ONE 1
#define PLAN_D2Z_ONE 2
#define PLAN_Z2D_ONE 3

// Target machine
// #define TARGET_MACHINE TITAN
// #define TARGET_MACHINE DWARF

#define TITAN 0
#define DWARF 1
#define TSUBAME 2
#define DWARF_ONE_NODE 3

#define AB3AM4 34
#define AB4AM5 45

#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
#undef CURRENT_CORRECTIONS
#endif

#ifdef BDG_MODE
#undef CURRENT_CORRECTIONS
#define FAST_CONST_EFFECTIVE_MASS_MODE
#endif

#ifdef WORK_IN_ROTATING_FRAME
#define CURRENT_CORRECTIONS
#undef FAST_CONST_EFFECTIVE_MASS_MODE
#endif

#endif
