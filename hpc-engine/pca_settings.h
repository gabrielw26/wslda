// Author: Gabriel Wlazlowski

#ifndef __PCA_SETTINGS__
#define __PCA_SETTINGS__

// ===============================================================================
// ================================ GLOBAL SETTINGS ==============================
// ===============================================================================

#ifdef WSLDA
// STATIC CODE
#define CODE "W-SLDA-TOOLKIT"
#define VERSION "0.1dev"
#include "netlib-lapack.h"
#include "predefines.h"

#elif TDWSLDA
// DYNAMIC  CODE
#define CODE "W-SLDA-TOOLKIT"
#define VERSION "0.1dev"
#include "predefines.h"

// HARD SET - TODO 
#define DX 1.0
#define DY 1.0                                                                                                                                       
#define DZ 1.0

#else
// DYNAMIC CODE - LEGACY MODE
#define CODE PCA_ASLDA
#define VERSION "1.10"

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

// Maximal number of parameters in params array
#define MAX_USER_PARAMS 32 

// Minimal density to avoid numerical problems
#define DENSEPSILON 1.0e-8

// spin-symmetric mode decreases computing time for factor about two
// #define SPINSYMMETRY_MODE

// active this flag in order to store quasi-particle energies for each measurment
// note that in case of 1d or 2d codes this can require much more space than measurments itself
// meaningful only for dynamic codes
// #define STORE_QPE

// Enable computation with extarnal delta field
// If this flag is active, you must provide body of delta_ext(...) function in pca_uext.h file
#define ENABLE_DELTA_EXT

// Enable computation with extarnal velocity field
// If this flag is active, you must provide body of vector_vext(...) function in pca_uext.h file
#define ENABLE_VELOCITY_EXT

// // To switch to cubic cut-off mode
// #define USE_CUBIC_CUTOFF

// activate this flag for setting code in testing mode with uniform system
#define UNIFORM_TEST_MODE

// compute kinetic energy density using formula tau ~ |nabla Psi|^2
// This is less acurate method than default, but we keep it for compability with older results
// #define TAU_COMPUTATION_VIA_GRADIENTS

// active rotating frame framework
// NOTE: in rotating frame the code uses constant variables dc_Omega_a and dc_Omega_b
// NOTE: the systems rotates along z-axis
// #define WORK_IN_ROTATING_FRAME

#endif

// Number of self-consistent iterations for U and delta computation
#define UD_SCITERS 15
// Mixing parameter for self-consistent algorithm - fraction of new solution used for mixing
#define UD_MIX_COEFF 0.75

// if particle number changed by this percentage then break the simulation
#define N_STABILITY_CRITERIA 0.25 

// Integration scheme AB - predictor, AM - correctior, number specify order
// #define ITEGRATION_SCHEME AB3AM4
#define INTEGRATION_SCHEME AB4AM5

// ===================================================================================
// =========================== PARAMETERS OF EDF =====================================
// ===================================================================================

#define SLDA 111
#define ASLDA 112
#define BDG 113
#define CUSTOMEDF 114

// if BDG_MODE then BdG functional is activated and aBdG parameter is active in dynamical codes
#if FUNCTIONAL==BDG
#define BDG_MODE
#define FAST_CONST_EFFECTIVE_MASS_MODE
#define A0 1.000
#define A1 0.0
#define A2 0.0
#endif

#if FUNCTIONAL==ASLDA
// effective mass - not 1.0 then current corrections are needed!
#define CURRENT_CORRECTIONS
#define A0 1.094
#define A1 0.156
#define A2 -0.532
#endif

#if FUNCTIONAL==SLDA
// effective mass - equal 1.0 then no current corrections
// activate this flag to skip computation of gradients of wf - significant spped up
#define FAST_CONST_EFFECTIVE_MASS_MODE
#define A0 1.000
#define A1 0.0
#define A2 0.0
#endif

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

// coordiates
#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2

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

// Package size
#ifdef SPINSYMMETRY_MODE

    #if FUNCTIONAL==ASLDA
    #define EXCHANGE_SIZE   7
    #endif
    #if FUNCTIONAL==SLDA
    #define EXCHANGE_SIZE   3
    #endif
    #if FUNCTIONAL==BDG
    #define EXCHANGE_SIZE   2
    #endif
    
#else

    #if FUNCTIONAL==ASLDA
    #define EXCHANGE_SIZE   12
    #endif
    #if FUNCTIONAL==SLDA
    #define EXCHANGE_SIZE   8
    #endif
    #if FUNCTIONAL==BDG
    #define EXCHANGE_SIZE   2
    #endif
    
#endif

#define PZHEEVR 1
#define PZHEEVD 2
#define ELPA 3

#endif
