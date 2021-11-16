/**
 * W-SLDA Toolkit
 * Author: Antoine Boulet
 * Creation date: 2021.09.21
 * */

#ifndef _SLDAE_FUNCTIONAL_
#define _SLDAE_FUNCTIONAL_



#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <complex.h>


#ifdef TDWSLDA
#define FDECORATOR __host__ __device__
#else
#define FDECORATOR 
#endif


/**
  ================================== TOOLS ===================================
 * */
// ---------------------------------------------------------------------------
#define HFB_ORDER 8
#define FUNCTIONAL_ORDER 2
// ---------------------------------------------------------------------------
#define LN2 log(2.)
#define LN2_2 pow(LN2, 2)
#define LN2_3 pow(LN2, 3)
#define LN2_4 pow(LN2, 4)
#define LN2_5 pow(LN2, 5)
#define M_PI_SQ pow(M_PI, 2)
// ---------------------------------------------------------------------------
#define DDCC_EPSILON 1.0e-32
// minimal _x to avoide divergences
// ---------------------------------------------------------------------------
FDECORATOR double rising_factorial (double x, int in);
FDECORATOR double falling_factorial (double x, int in);
double binomial_coefficient (int in, int ik);
double exponential_integral (double x);
double pexp_bell_polynomial (int in, int ik, double * g);
double xm (int dn, double x, int * ip);
// ---------------------------------------------------------------------------
double power_lrule (int dn, double _x, int i, int dk, int * id, double (*f) (int, double, int *));
double product_lrule (int dn, double _x, int dkf, int dkg, int * idf, int * idg, double (*f) (int, double, int *), double (*g) (int, double, int *));
double inverse_lrule (int dn, double _x, int dk, int * id, double (*f) (int, double, int *));
double composed_lrule (int dn, double _x, int dkf, int dkg, int * idf, int * idg, double (*f) (int, double, int *), double (*g) (int, double, int *));
// ---------------------------------------------------------------------------
double b_expansion (int dp, double y_x, int * idx);
double c_expansion (int dp, double y_x, int * idx);
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------


/**
  ================================ FUNCTIONAL ================================
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
double s_aps (int dn, double _x, int * id);
double s_aps_2d1 (int dn, double _x, int * id);
double pairing_bcs (int dn, double _x, int * id);
double pairing_fit (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
double ground_state_energy (int dn, double _x, int * id);
double chemical_potential (int dn, double _x, int * id);
double inverse_effective_mass (int dn, double _x, int * id);
double pairing_gap (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
double one_over_a (int dn, double _x, int * id);
double one_over_h (int dn, double _x, int * id);
double h_sq (int dn, double _x, int * id);
double h_over_a (int dn, double _x, int * id);
double a_over_h (int dn, double _x, int * id);
double h_over_a_pdk (int dn, double _x, int * idpk);
double log_h_over_a (int dn, double _x, int * id);
double b_series_coefficient (int dn, double _x, int * idp);
double c_series_coefficient (int dn, double _x, int * idp);
double b_series (int dn, double _x, int * id);
double c_series (int dn, double _x, int * id);
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------


/**
  ================================ PARAMETERS ================================
**/
// ---------------------------------------------------------------------------
double a_hfb (int dn, double _x, int * id);
double b_hfb (int dn, double _x, int * id);
double inverse_c_hfb (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
double alpha_parameter (int dn, double _x, int * id);
double beta_parameter (int dn, double _x, int * id);
double inverse_gamma_parameter (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
double a_functional (double _x, int * id);
double b_functional (double _x, int * id);
double b_functional_aps (double _x, int * id);
double c_functional (double _x, int * id);
// ---------------------------------------------------------------------------
double pcc_renormalization(double _x, double lmu, int * id);
// ---------------------------------------------------------------------------

#endif
