/**
 * W-SLDA Toolkit
 * Author: Antoine Boulet
 * Creation date: 2021.09.21
 *
 * td update (AB): 2021.11.25
 * Main change for td udate:
 * - simplification of functions (SLDA and functional parameters)
     by Padé[4/4] approximations
 * - avoid issue due to the use of recursive function in cuda
 * - the "old" functions are defined in sldae_dev.h
 * */

#ifndef _SLDAE_FUNCTIONAL_
#define _SLDAE_FUNCTIONAL_

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#ifdef TDWSLDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#define FDECORATOR __host__ __device__
#else
#define FDECORATOR
#endif // #ifdef TDWSLDA

// ---------------------------------------------------------------------------
#define LN2 log(2.)
#define LN2_2 pow(LN2, 2)
#define LN2_3 pow(LN2, 3)
#define LN2_4 pow(LN2, 4)
#define LN2_5 pow(LN2, 5)
#define M_PI_SQ pow(M_PI, 2)
// ---------------------------------------------------------------------------
#define AF_APS_UFG (0.8403540)
#define BF_APS_UFG (-0.281887)
#define CF_APS_UFG (-14.95850)
// ---------------------------------------------------------------------------
// Padé[4/4] functions (used by default in st and td code)
FDECORATOR double alpha_parameter_d0 (double _x);
FDECORATOR double beta_parameter_d0 (double _x);
FDECORATOR double inverse_gamma_parameter_d0 (double _x);
FDECORATOR double alpha_parameter_d1 (double _x);
FDECORATOR double beta_parameter_d1 (double _x);
FDECORATOR double inverse_gamma_parameter_d1 (double _x);
FDECORATOR double alpha_parameter_d2 (double _x);
FDECORATOR double beta_parameter_d2 (double _x);
FDECORATOR double inverse_gamma_parameter_d2 (double _x);
FDECORATOR double a_functional_d0 (double _x);
FDECORATOR double b_functional_d0 (double _x);
FDECORATOR double c_functional_d0 (double _x);
// ---------------------------------------------------------------------------

#ifdef SLDAE_DEV_MODE
// ---------------------------------------------------------------------------
/**
   ================================= TOOLS ===================================
 * */
// ---------------------------------------------------------------------------
#define HFB_ORDER 8
#define FUNCTIONAL_ORDER 2
// ---------------------------------------------------------------------------
#define DDCC_EPSILON 1.0e-32
// minimal _x to avoide divergences
// ---------------------------------------------------------------------------
FDECORATOR double rising_factorial (double x, int in);
FDECORATOR double falling_factorial (double x, int in);
FDECORATOR double binomial_coefficient (int in, int ik);
FDECORATOR double exponential_integral (double x);
FDECORATOR double pexp_bell_polynomial (int in, int ik, double * g);
FDECORATOR double xm (int dn, double x, int * ip);
// ---------------------------------------------------------------------------
FDECORATOR double power_lrule (int dn, double _x, int i, int dk, int * id, double (*f) (int, double, int *));
FDECORATOR double product_lrule (int dn, double _x, int dkf, int dkg, int * idf, int * idg, double (*f) (int, double, int *), double (*g) (int, double, int *));
FDECORATOR double inverse_lrule (int dn, double _x, int dk, int * id, double (*f) (int, double, int *));
FDECORATOR double composed_lrule (int dn, double _x, int dkf, int dkg, int * idf, int * idg, double (*f) (int, double, int *), double (*g) (int, double, int *));
// ---------------------------------------------------------------------------
FDECORATOR double b_expansion (int dp, double y_x, int * idx);
FDECORATOR double c_expansion (int dp, double y_x, int * idx);
// ---------------------------------------------------------------------------
/**
   =============================== FUNCTIONAL ================================
**/
// ---------------------------------------------------------------------------
#define GSE_UFG 0.3582341
#define PGF_UFG 0.4600000
#define IEM_UFG 0.8403361
// --------------------------------------------------------------------------
#define U_APS (5. / 24.)
#define V_APS (U_APS / tan(3. * M_PI / 16. * (1. - GSE_UFG)))
#define V_OVER_U_APS (V_APS / ((6. / (35. * M_PI)) * (11. - 2. * LN2)))
#define W_APS (V_OVER_U_APS * (24. / (35. * M_PI)) * (1. - 7. * LN2))
#define Y_APS (4. / 5.)
#define Z_APS (8. / exp(2.) / PGF_UFG)
#define X_APS ((9. * M_PI * pow(U_APS, 2)) / (7. * U_APS * V_APS) * ((5. * U_APS * W_APS) / (9. * M_PI * pow(U_APS, 2)) - (pow(U_APS, 2) + pow(V_APS, 2)) / U_APS * (IEM_UFG - 1.)))
#define A_APS ((5. * U_APS * W_APS - 7. * X_APS * U_APS * V_APS) / (9. * M_PI * pow(U_APS, 2)))
#define B_APS ((2. * W_APS + 7. * X_APS * V_APS) / (9. * M_PI * pow(U_APS, 2)))
// ---------------------------------------------------------------------------
FDECORATOR double s_aps (int dn, double _x, int * id);
FDECORATOR double s_aps_2d1 (int dn, double _x, int * id);
FDECORATOR double pairing_bcs (int dn, double _x, int * id);
FDECORATOR double pairing_fit (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
FDECORATOR double ground_state_energy (int dn, double _x, int * id);
FDECORATOR double chemical_potential (int dn, double _x, int * id);
FDECORATOR double inverse_effective_mass (int dn, double _x, int * id);
FDECORATOR double pairing_gap (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
FDECORATOR double one_over_a (int dn, double _x, int * id);
FDECORATOR double one_over_h (int dn, double _x, int * id);
FDECORATOR double h_sq (int dn, double _x, int * id);
FDECORATOR double h_over_a (int dn, double _x, int * id);
FDECORATOR double a_over_h (int dn, double _x, int * id);
FDECORATOR double h_over_a_pdk (int dn, double _x, int * idpk);
FDECORATOR double log_h_over_a (int dn, double _x, int * id);
FDECORATOR double b_series_coefficient (int dn, double _x, int * idp);
FDECORATOR double c_series_coefficient (int dn, double _x, int * idp);
FDECORATOR double b_series (int dn, double _x, int * id);
FDECORATOR double c_series (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
/**
   =============================== PARAMETERS ================================
**/
// ---------------------------------------------------------------------------
FDECORATOR double a_hfb (int dn, double _x, int * id);
FDECORATOR double b_hfb (int dn, double _x, int * id);
FDECORATOR double inverse_c_hfb (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
FDECORATOR double alpha_parameter (int dn, double _x, int * id);
FDECORATOR double beta_parameter (int dn, double _x, int * id);
FDECORATOR double inverse_gamma_parameter (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
FDECORATOR double a_functional (double _x, int * id);
FDECORATOR double b_functional (double _x, int * id);
FDECORATOR double c_functional (double _x, int * id);
// ---------------------------------------------------------------------------
#endif // #ifdef SLDAE_DEV_MODE

#endif // #ifndef _SLDAE_FUNCTIONAL_
