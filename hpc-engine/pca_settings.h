// Author: Gabriel Wlazlowski

#ifndef __PCA_SETTINGS__
#define __PCA_SETTINGS__

// ===============================================================================
// ================================ GLOBAL SETTINGS ==============================
// ===============================================================================
#include "../VERSION.h"

#ifdef WSLDA
// STATIC CODE
#include "netlib-lapack.h"
#include "predefines.h"

#elif TDWSLDA
// DYNAMIC  CODE
#include "predefines.h"

// uncomment this if you want to assure that data is copied to the working buffer 
// before execution of derivative computation
// it is related to issue encountered for hipfft implementation
// #define DERIVATIVE_COPY_DATA_MODE
#else

#error "You need to select WSLDA or TDWSLDA!"

#endif

// check predefines
#include "wslda_predefines_test.h"

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
#define SLDAE 115

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

#if FUNCTIONAL==SLDAE
// effective mass - not 1.0 then current corrections are needed!
#define CURRENT_CORRECTIONS
#endif


#if FUNCTIONAL==SLDA
// effective mass - equal 1.0 then no current corrections
// activate this flag to skip computation of gradients of wf - significant spped up
#define FAST_CONST_EFFECTIVE_MASS_MODE
#define A0 1.000
#define A1 0.0
#define A2 0.0
#endif

#ifndef A0
#define A0 1.000
#define A1 0.0
#define A2 0.0
#endif

// normal part
#define G0 0.357
#define G1 0.642

// regularization function parameters
#ifdef ASLDA_STABILIZATION_EXCLUDE_BELOW_DENISTY
#define P_NMIN ASLDA_STABILIZATION_EXCLUDE_BELOW_DENISTY
#else
#define P_NMIN 1.0e-7
#endif

#ifdef ASLDA_STABILIZATION_RETAIN_ABOVE_DENSITY
#define P_NMAX ASLDA_STABILIZATION_RETAIN_ABOVE_DENSITY
#else
#define P_NMAX 1.0e-5
#endif

#define P_ALPHA 1.0

// #define P_NMIN 1.0e-8
// #define P_NMAX 1.0e-6
// #define P_ALPHA 1.0

// ===================================================================================
// =================================== TECHNICAL =====================================
// ===================================================================================

#define NUMERICAL_ZERO 1.0e-16

#define DXYZ (DX*DY*DZ)
#define DXY (DX*DY)

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

// ------------- observables -----------------
// energy contributions
#define ENERGYITEMS 7
#define EKIN        0
#define EPOT        1
#define EPAIR       2
#define ECURRENT    3
#define EPOTEXT     4
#define EPAIREXT    5
#define EVELEXT     6
// other contributions
// total number of items in TDWSLDA codes
#define TDWSLDAITEMS (ENERGYITEMS+4)
#define NPARTA      7
#define NPARTB      8
#define LZA         9
#define LZB         10

// total number of items in WSLDA codes
#define WSLDAITEMS (ENERGYITEMS+1)
#define ENTROPY     7

// cufft plans
#define CUFFT_NUMBER_OF_PLANS 4
#define PLAN_Z2Z_BATCH 0
#define PLAN_Z2Z_ONE 1
#define PLAN_D2Z_ONE 2
#define PLAN_Z2D_ONE 3

// Technical variable - amount of memory that is locked for axiliary array used in apply_hamiltonian
// i.e: sizeof(cufftDoubleComplex)*NXYZ*PCA_WORKSPACE_SHIFT
#define PCA_WORKSPACE_SHIFT 7
#ifdef ENABLE_VELOCITY_EXT
#undef PCA_WORKSPACE_SHIFT
#define PCA_WORKSPACE_SHIFT 11
#endif

// Target machine
// #define TARGET_MACHINE TITAN
// #define TARGET_MACHINE DWARF

#define TITAN 0
#define DWARF 1
#define TSUBAME 2
#define DWARF_ONE_NODE 3

#define AB3AM4 34
#define AB4AM5 45

#ifdef ENABLE_MODIFY_POTENTIALS
#undef BDG_MODE
#endif

#ifdef SLDAE_FORCE_A1
#define FAST_CONST_EFFECTIVE_MASS_MODE
#endif 

#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
#undef CURRENT_CORRECTIONS
#endif

#ifdef BDG_MODE
#undef CURRENT_CORRECTIONS
#define FAST_CONST_EFFECTIVE_MASS_MODE
#endif

#ifdef ENABLE_VELOCITY_EXT
#define CURRENT_CORRECTIONS
#undef FAST_CONST_EFFECTIVE_MASS_MODE
#endif

// Package size
#ifdef SPINSYMMETRY_MODE

    #if FUNCTIONAL==ASLDA
    #define EXCHANGE_SIZE   7
    #elif FUNCTIONAL==SLDA
    #define EXCHANGE_SIZE   3
    #elif FUNCTIONAL==BDG
    #define EXCHANGE_SIZE   2
    #else
    #define EXCHANGE_SIZE   7
    #endif
    
    #ifdef CURRENT_CORRECTIONS
    #undef EXCHANGE_SIZE
    #define EXCHANGE_SIZE 7
    #endif

#else

    #if FUNCTIONAL==ASLDA
    #define EXCHANGE_SIZE   12
    #elif FUNCTIONAL==SLDA
    #define EXCHANGE_SIZE   8
    #elif FUNCTIONAL==BDG
    #define EXCHANGE_SIZE   2
    #else
    #define EXCHANGE_SIZE   12
    #endif
    
    #ifdef CURRENT_CORRECTIONS
    #undef EXCHANGE_SIZE
    #define EXCHANGE_SIZE 12
    #endif

#endif

#define PZHEEVR 1
#define PZHEEVD 2
#define ELPA 3

#define SPHERICAL_CUTOFF 88
#define CUBIC_CUTOFF 89

#if REGULARIZATION_SCHEME==CUBIC_CUTOFF
#define USE_CUBIC_CUTOFF
#endif
// otherwise use speherical cutoff
#define REGULARIZATION_SCHEME_K_CONST 2.442749607806335

// pairing
#if REGULARIZATION_SCHEME==CUBIC_CUTOFF
#define GAMMA0 (-11.11*1.60)
#else
#define GAMMA0 -11.11
#endif

// defaults for ELPA
#ifndef ELPA_USE_SOLVER
#define ELPA_USE_SOLVER ELPA_SOLVER_1STAGE
#endif
#ifndef ELPA_USE_COMPLEX_KERNEL
#define ELPA_USE_COMPLEX_KERNEL ELPA_2STAGE_COMPLEX_DEFAULT
#endif
#ifndef ELPA_USE_REAL_KERNEL
#define ELPA_USE_REAL_KERNEL ELPA_2STAGE_REAL_DEFAULT
#endif
#ifndef ELPA_NEV_FRACTION
#define ELPA_NEV_FRACTION 1.0
#endif

#ifndef GPUS_PER_NODE
#define GPUS_PER_NODE 1
#endif

// ----- for TESTSUITE -----
// default energy error
#ifndef TS_EERR
#define TS_EERR 1.0e-4
#endif
// default particle error
#ifndef TS_NERR
#define TS_NERR 1.0e-3
#endif
// default chemical potential error
#ifndef TS_MUERR
#define TS_MUERR 1.0e-3
#endif
// default entropy error
#ifndef TS_SERR
#define TS_SERR 1.0e-3
#endif

// ------------ math -------------
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.570796326794896558
#endif

// to mantain legcy
#ifdef HAMILTONIAN_IS_REAL
#define MATRIX_IS_REAL
#endif

#ifndef API_VERSION
#define API_VERSION 20220218
#endif

#endif
