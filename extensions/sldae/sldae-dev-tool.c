// ###########################################################################
// ---------------------------------------------------------------------------
//
//  sldae-dev-tool.c - version 1 (2022/01/27)
//
//  requirment - gsl library
//  usage > gcc/clang sldae-dev-tool.c -lm -lgsl -lgslcblas -o sldae-dev-tool
//
//  Antoine Boulet (antoine.boulet@protonmail.com)
//
//  attachment - arXiv:2201.07626
//
// ---------------------------------------------------------------------------
/*
 *
 The main function interpolate sldae parameters p = {beta, gamma, B, C}
 with a padé[4/4] function as used in the W_SLDA Toolkit (x = |akF|):
   p       a1 x + a2 x^2 + a3 x^3 + c x^4
 ----- = ----------------------------------
 p_UFG   1 + b1 x + b2 x^2 + b3 x^3 + c x^4
 [implementation in W-SLDA @ hpc-engine/sldae_functional.c]

 By default, the APS[x,y,z] functional is implemented, but user can
 redefine the ground state energy, the chemical potential,
 the inverse effective mass and the pairing gap (see bellow for more details).
 *
 */
// ---------------------------------------------------------------------------
// ###########################################################################

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

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


// ###########################################################################
// ---------------------------------------------------------------------------
/*
    //////////////////////////////////////////////////////////////////////
    ===  IMPLEMENTATION OF FUNCTIONS FOR THE SLDA EXTENDED FUNCTIONAL  ===
    //////////////////////////////////////////////////////////////////////

  # -- NOTATIONS AND CONVENTIONS
       =========================

  1/ variables
     _x = |akF|
          [density dependent coupling constant (DDCC)]
    y_x = \ln (\eta_x / \alpha_x)
          [appearing in b and c function expansions]
     dn = (\partial / \partial  _x)^n := \partial_n
          [functional derivative of order n according to the DDCC]
     dp = (\partial / \partial y_x)^n
          [derivative of coefficients appearing in b and c function expansions]
     id [2] = {FUNCTIONAL_ID, PAIRING_ID} := {fid, pid}
          [for general purpose]

  2/ functions
     f(int dk, double x, int id[]) = \partial_k f_id(x)

  # -- REQUIRMENTS (defined by user)
       ===========
  1/ ground state energy (and derivatives) as a function of _x
    [ground_state_energy(dn, _x, id)]
  2/ chemical potential (and derivatives) as a function of _x
    [chemical_potential(dn, _x, id)]
  3/ inverse effective mass (and derivatives) as a function of _x
    [inverse_effective_mass(dn, _x, id)]
  4/ pairing gap (and derivatives) as a function of _x
    [pairing_gap(dn, _x, id)]

  * remark 1 *
  The rest of the function provide automatic results for the SLDA
  extended methods as described in the attachment.
  Eventually, higher order correction, necesarry in case of
  (\eta_x / \alpha_x) > 0.5, can be implemented by adding
  higer order in functions:
        - b_expansion(int dp, double y_x, int idx[])
        - c_expansion(int dp, double y_x, int idx[])
          * notations *
              idx[0]: coefficient indices
              y_x = \ln (\eta_x / \alpha_x)
              dp correspond to the pth derivative according to y_x
  The trucation of such expansion is controled by HFB_ORDER corresponding
  to the tructation of the HFB series (set by user).
  This improvment require full version of Mathematica (free Wolfram cloud is
  unable to provide higher terms).

  * remark 2 *
  The derivatives (of product, compostition, etc.) are obtained using
  general Leibniz rule such that, only the derivatives of the required
  function are needed up to FUNCTIONAL_ORDER (see bellow) in order to
  express functional parameters as a series of the SLDA parameters, i.e.
  FUNCTIONAL_ORDER correspond to the tructation of the functional series.

  # -- functions defiend (follow the notes)
       =================
      [to be called in the main code to calculate energy and potentials]
  1/ HFB parameters: a_x, b_x, c_x
      - a_hfb(dn, _x, id)
      - b_hfb(dn, _x, id)
      - inverse_c_hfb(dn, _x, id)
  2/ SLDA parameters: \alpha_x, \beta_x, \gamma_x
      - alpha_parameter(dn, _x, id)
      - beta_parameter(dn, _x, id)
      - inverse_gamma_parameter(dn, _x, id)
  3/ functional parameters: A_x, B_x, C_x
      - a_functional(_x, id)
      - b_functional(_x, id)
      - c_functional(_x, id)
  These functions are general and do not need to be modified if the
  functional or the parameters FUNCTIONAL_ORDER and HFB_ORDER are changed.
*/
// ---------------------------------------------------------------------------
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ---------------------------------------------------------------------------
// Declaration of the main functions
// ---------------------------------------------------------------------------
// ###########################################################################
// phyical quantities
double ground_state_energy (int, double, int *);
double chemical_potential (int, double, int *);
double inverse_effective_mass (int, double, int *);
double pairing_gap (int, double, int *);
// HFB parameters
double a_hfb (int, double, int *);
double b_hfb (int, double, int *);
double inverse_c_hfb (int, double, int *);
// SLDA parameters
double alpha_parameter (int, double, int *);
double beta_parameter (int, double, int *);
double inverse_gamma_parameter (int, double, int *);
// functional parameters
double a_functional (double, int *);
double b_functional (double, int *);
double c_functional (double, int *);
// others
double xm (int, double, int *);
double product_lrule (int, double, int, int, int *, int *,
                      double (*) (int, double, int *),
                      double (*) (int, double, int *));
// APS[x,y,z] functional related
double s_aps (int, double, int *);
double s_aps_2d1 (int, double, int *);
double pairing_bcs (int, double, int *);
double pairing_fit (int, double, int *);
// ---------------------------------------------------------------------------
#define GSE_UFG 0.3582341
#define PGF_UFG 0.4600000
#define IEM_UFG 0.8403361
// --------------------------------------------------------------------------
#define U_APS 0.208333
#define V_APS 0.524595
#define W_APS -0.840781
#define Y_APS 0.800000
#define Z_APS 2.353657
#define X_APS -0.753128
#define A_APS -0.244173
#define B_APS -3.623881
// ---------------------------------------------------------------------------
// ###########################################################################



// ---------------------------------------------------------------------------
#define HFB_ORDER 8
#define FUNCTIONAL_ORDER 2
// ---------------------------------------------------------------------------



// ###########################################################################
// ---------------------------------------------------------------------------
// ground state energy, chemical potential, inverse effective mass,
// and pairing gap as a function of _x = |akF| (can be re-defiend by user)
// ---------------------------------------------------------------------------

// ground state energy
// ===================
// \parial_n \xi_x
double
ground_state_energy (int dn, double _x, int id [])
{
    int ip0 [1] = {0}; double x0 = xm (dn, _x, ip0);
    // int ip1 [1] = {1}; // double x1 = xm (dn, _x, ip1);
    // int ip2 [1] = {2}; // double x2 = xm (dn, _x, ip2);
    // ...
    return x0 - 16. / (3. * M_PI) * s_aps (dn, _x, id);
}

// chemical potential
// ==================
// \parial_n \zeta_x
double
chemical_potential (int dn, double _x, int id [])
{
    int ip0 [1] = {0}; double x0 = xm (dn, _x, ip0);
    int ip1 [1] = {1}; // double x1 = xm (dn, _x, ip1);
    // int ip2 [1] = {2}; // double x2 = xm (dn, _x, ip2);
    // ...
    return x0 - 16. / (3. * M_PI) * s_aps (dn, _x, id) -
    16. / (15. * M_PI) * product_lrule (dn, _x, 0, 1, ip1, id, xm, s_aps);
}

// inverse effective mass
// ======================
// \parial_n \alpha_x
double
inverse_effective_mass (int dn, double _x, int id [])
{
    int ip0 [1] = {0}; double x0 = xm (dn, _x, ip0);
    // int ip1 [1] = {1}; // double x1 = xm (dn, _x, ip1);
    int ip2 [1] = {2}; // double x2 = xm (dn, _x, ip2);
    int ip3 [1] = {3}; // double x3 = xm (dn, _x, ip3);
    // ...
    return x0 + A_APS * product_lrule (dn, _x, 0, 1, ip2, id, xm, s_aps) +
                B_APS * product_lrule (dn, _x, 0, 0, ip2, id, xm, s_aps_2d1) +
                B_APS * V_APS * product_lrule (dn, _x, 0, 0, ip3, id, xm, s_aps_2d1);
}

// pairing gap
// ===========
// \parial_n \eta_x
double
pairing_gap (int dn, double _x, int id [])
{
    if (pairing_bcs (0, _x, id) == 0) {
      return 0.;
      // avoid instability of derivative when vanishing pairing
    } else {
      return product_lrule (dn, _x, 0, 0, id, id, pairing_bcs, pairing_fit);
    }
}

// ---------------------------------------------------------------------------
// ###########################################################################




// ---------------------------------------------------------------------------



/* *
  ================================== TOOLS ===================================
* */

// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// usefull mathematical functions                        [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// rising factorials
// =================
// (x)^n = \prod_{k=0}^{n-1} (x + k)
// used to calulatate derivatives of the BCS pairing gap function
double
rising_factorial (double x, int in)
{
  if (in == 0) {
    return 1.;
  } else if (in > 0) {
    int ik; double r_ = 1.;
    for (ik = 0; ik < in; ik++) {
      r_ *= (x + ik);
    }
    return r_;
  } else {
    return 0.;
  }
}

// falling factorials
// ==================
// (x)_n = \prod_{k=0}^{n-1} (x - k)
double
falling_factorial (double x, int in)
{
  if (in == 0) {
    return 1.;
  } else if (in > 0) {
    int ik; double r_ = 1.;
    for (ik = 0; ik < in; ik++) {
      r_ *= (x - ik);
    }
    return r_;
  } else {
    return 0.;
  }
}

// binomial coefficients
// =====================
// C_{n,p} = n! / k! / (n-k)!
double
binomial_coefficient (int in, int ik)
{
  return tgamma (in + 1.) / tgamma (ik + 1.) / tgamma (in - ik + 1.);
}

// partial exponential Bell polynomials
// ====================================
// B_{n,k}(g_1, g_2, ..., g_{n-k+1})
// https://en.wikipedia.org/wiki/Bell_polynomials
// used for the Faà di Bruno's formula
double
pexp_bell_polynomial (int in, int ik, double g[])
{
  if (in == 0 && ik == 0) {
    return 1.;
  } else if (in == 0 && ik > 0) {
    return 0.;
  } else if (in > 0 && ik == 0) {
    return 0.;
  } else if (in > 0 && ik > 0) {
    int ij; double r_ = 0;
    for (ij = 1; ij <= in - ik + 1; ij++) {
      r_ += binomial_coefficient (in - 1, ij - 1) *
            g[ij - 1] *
            pexp_bell_polynomial (in - ij, ik - 1, g);
    }
    return r_;
  } else {
    return 0.;
  }
}

// unit polynomial functions
// =========================
// \partial_n x^p
// [implmented following the convention: ip[0] correspond to the pow of x]
double
xm (int dn, double x, int ip[])
{
  if (ip[0] == 0 && dn != 0) {
    return 0.;
  } else if (dn == 0) {
    return pow (x, ip[0]);
  } else if (dn > 0) {
    int ip_tmp[1] = {ip[0] - 1};
    return ip[0] * xm (dn - 1, x, ip_tmp);
    // polynomial Leibniz rule
  } else {
    return 0.;
  }
}
// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// general derivatives using Leibniz rules               [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// pow Leibniz rule
// ==================
// \parial_n [f^(k)(x)]^i =
// \sum_{p=0}^n C_{n,p} f^(k+p) \partial_{n-p} [f^(k)(x)]^{i-1}
double
pow_lrule
(
  int dn, double _x,
  int i, int dk,
  int id [],
  double (*f) (int, double, int *)
)
{
  if (i == 0 && dn != 0) {
    return 0.;
  } else if (dn == 0) {
    return pow ((*f) (dk, _x, id), i);
  } else if (dn > 0) {
    int dp; double r_ = 0;
    for (dp = 0; dp <= dn; dp++) {
      r_ += binomial_coefficient (dn, dp) *
            (*f) (dk + dp, _x, id) *
            pow_lrule (dn - dp, _x, i - 1, dk, id, (*f));
    }
    return r_;
  } else {
    return 0.;
  }
}

// product Leibniz rule
// ====================
// \parial_n [f^(kf)(x) * g^(kg)(x)] =
// \sum_{p=0}^n C_{n,p} f^(kf+p) * g^(kg+n-p)(x)
double
product_lrule
(
  int dn, double _x,
  int dkf, int dkg,
  int idf [], int idg [],
  double (*f) (int, double, int *),
  double (*g) (int, double, int *)
)
{
  if (dn == 0) {
    return (*f) (dkf, _x, idf) * (*g) (dkg, _x, idg);
  } else if (dn > 0) {
    int dp; double r_ = 0;
    for (dp = 0; dp <= dn; dp++) {
      r_ += binomial_coefficient (dn, dp) *
            (*f) (dkf + dp, _x, idf) *
            (*g) (dkg + dn - dp, _x, idg);
    }
    return r_;
  } else {
    return 0.;
  }
}

// inverse Leibniz rule
// ====================
// \parial_n [1 / f^(k)(x)]=
// \sum_{p=0}^n (-1)^p C_{n,p} (n+1) / (p+1)
//              (\partial_n [f^(k)(x)]^p) / [f^(k)(x)]^{p+1}
double
inverse_lrule
(
  int dn, double _x,
  int dk,
  int id [],
  double (*f) (int, double, int *)
)
{
  if (dn == 0) {
    return 1. / (*f) (dk, _x, id);
  } else if (dn > 0) {
    int dp; double r_ = 0;
    for (dp = 0; dp <= dn; dp++) {
      r_ += pow(-1., dp) * binomial_coefficient (dn, dp) *
            (dn + 1.) / (dp + 1.) *
            pow_lrule (dn, _x, dp, dk, id, (*f)) /
            pow((*f) (dk, _x, id), dp + 1);
      }
      return r_;
  } else {
      return 0.;
  }
}

// Faà di Bruno's formula
// ======================
// \partial_n [f^(kf)(g^(kg)(x))]
double
composed_lrule
(
  int dn, double _x,
  int dkf, int dkg,
  int idf [], int idg [],
  double (*f) (int, double, int *),
  double (*g) (int, double, int *)
)
{
  if (dn == 0) {
    return (*f) (dkf, (*g) (dkg, _x, idg), idf);
  } else if (dn > 0) {
    // vector of sucessive derivative of g
    int dj;
    // double dg [dn]; dn_max = FUNCTIONAL_ORDER
    double dg [FUNCTIONAL_ORDER];
    for (dj = 0; dj < dn; dj++) {
      dg [dj] = (*g) (dkg + dj + 1, _x, idg);
    }
    // Faà di Bruno's formula
    int dp; double r_ = 0.;
    for (dp = 1; dp <= dn; dp++) {
      r_ += (*f) (dkf + dp, (*g) (dkg, _x, idg), idf) *
            pexp_bell_polynomial (dn, dp, dg);
    }
    return r_;
  } else {
    return 0.;
  }
}

// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ---------------------------------------------------------------------------
// functions Bn and Cn for the expansion of the SLDA parameters b and c
// (higher order can be added systematically by user)
// ---------------------------------------------------------------------------

// (\partial / \partial y)^n \mathcal{B}_p(y)
double
b_expansion (int dp, double y_x, int idx []) // idx [0]: index of coefficient p
{
  int ip0 [1] = {0}; double y0 = xm (dp, y_x, ip0);
  int ip1 [1] = {1}; double y1 = xm (dp, y_x, ip1);
  int ip2 [1] = {2}; double y2 = xm (dp, y_x, ip2);
  int ip3 [1] = {3}; double y3 = xm (dp, y_x, ip3);
  int ip4 [1] = {4}; double y4 = xm (dp, y_x, ip4);
  // int ip5 [1] = {5}; double y5 = xm (dp, y_x, ip5);
  // ...

  if (idx [0] > HFB_ORDER) {
    return 0.;
  } else if (idx [0] == 0) {
    return -y0;
  } else if (idx [0] == 1) {
    return 0.;
  } else if (idx [0] == 2) {
    return -y1 / 4. +
            y0 / 8. +
            y0 * 3. * LN2 / 4.;
  } else if (idx [0] == 3) {
    return 0.;
  } else if (idx [0] == 4) {
    return +y2 * (3. / 64.) +
            y1 * (1. / 128.) * (5. - 36. * LN2) +
            y0 * (7. / 512.) +
            y0 * (3. / 128.) * LN2 * (18. * LN2 - 5.);
  } else if (idx [0] == 5) {
    return 0.;
  } else if (idx [0] == 6) {
    return -y3 * (7. / 384.) +
            y2 * ((168. * LN2 - 37.) / 1024.) -
            y1 * (3. * (23. + 8. * LN2 * (84. * LN2 - 37.)) / 4096.) +
            y0 * (209. / 24576.) +
            y0 * (9. * LN2 * (56. * LN2 - 23.) * (4. * LN2 - 1) / 4096.);
  } else if (idx [0] == 7) {
    return 0.;
  } else if (idx [0] == 8) {
    return +y4 * (55. / 6144.) +
            y3 * (43. / 1536.) -
            y3 * (55. / 512.) * LN2 +
            y2 * (1395. / 65536.) +
            y2 * (495. * LN2_2 / 1024.) -
            y2 * (129. / 512.) * LN2 -
            y1 * (2141. / 393216.) -
            y1 * (495. / 512.) * LN2_3 +
            y1 * (387. / 512.) * LN2_2 -
            y1 * (4185. * LN2 / 32768.) +
            y0 * (8657. / 3145728.) +
            y0 * (1485. * LN2_4 / 2048.) -
            y0 * (387. * LN2_3 / 512.) +
            y0 * (12555. * LN2_2 / 65536.) +
            y0 * (2141. * LN2 / 131072.);
  } else {
    return 0.;
    // expansion can be improved by adding higher order terms
  }
}

// (\partial / \partial y)^n \mathcal{C}_p(y)
double
c_expansion (int dp, double y_x, int idx []) // idx [0]: index of coefficient p
{
  int ip0 [1] = {0}; double y0 = xm (dp, y_x, ip0);
  int ip1 [1] = {1}; double y1 = xm (dp, y_x, ip1);
  int ip2 [1] = {2}; double y2 = xm (dp, y_x, ip2);
  int ip3 [1] = {3}; double y3 = xm (dp, y_x, ip3);
  int ip4 [1] = {4}; double y4 = xm (dp, y_x, ip4);
  int ip5 [1] = {5}; double y5 = xm (dp, y_x, ip5);
  // ...

  if (idx [0] > HFB_ORDER) {
    return 0.;
  } else if (idx [0] == 0) {
    return +y1 / 4. +
            y0 / 2. -
            y0 * 3. * LN2 / 4.;
  } else if (idx [0] == 1) {
    return 0.;
  } else if (idx [0] == 2) {
    return +y2 * (1. / 32.) -
            y1 * (3. / 16.) * LN2 +
            y0 * (1. / 64.) +
            y0 * (9. * LN2_2 / 32.);
  } else if (idx [0] == 3) {
    return 0.;
  } else if (idx [0] == 4) {
    return -y3 * (1. / 128.) +
            y2 * ((72. * LN2 - 13) / 1024.) -
            y1 * (3. * (3. + 2. * (36. * LN2 - 13.) * LN2) / 1024.) -
            y0 * (15. / 8192.) +
            y0 * (9. * LN2 * (3. + LN2 * (24. * LN2 - 13.)) / 1024.);
  } else if (idx [0] == 5) {
    return 0.;
  } else if (idx [0] == 6) {
    return +y4 * (5. / 1536.) +
            y3 * (29. / 3072.) -
            y3 * (5. / 128.) * LN2 +
            y2 * (315. / 32768.) +
            y2 * (45. / 256.) * LN2_2 -
            y2 * (87. * LN2 / 1024.) +
            y1 * (385. / 196608.) -
            y1 * (45. / 128.) * LN2_3 +
            y1 * (261. * LN2_2 / 1024.) -
            y1 * (945. * LN2 / 16384.) -
            y0 * (701. / 393216.) +
            y0 * (135. * LN2_4 / 512.) -
            y0 * (261. * LN2_3 / 1024.) +
            y0 * (2835. * LN2_2 / 32768.) -
            y0 * (385. * LN2 / 65536.);
  } else if (idx [0] == 7) {
    return 0.;
  } else if (idx [0] == 8) {
    return -y5 * (27. / 16384.) -
            y4 * (111. / 16384.) +
            y4 * (405. * LN2 / 16384.) -
            y3 * (1257. / 131072.) -
            y3 * (1215. * LN2_2 / 8192.) +
            y3 * (333. * LN2 / 4096.) -
            y2 * (3695. / 1048576.) +
            y2 * (3645. * LN2_3 / 8192.) -
            y2 * (2997. * LN2_2 / 8192.) +
            y2 * (11313. * LN2 / 131072.) +
            y1 * (5251. / 4194304.) -
            y1 * (10935. * LN2_4 / 16384.) +
            y1 * (2997. * LN2_3 / 4096.) -
            y1 * (33939. * LN2_2 / 131072.) +
            y1 * (11085. * LN2 / 524288.) -
            y0 * (6223. / 8388608.) +
            y0 * (6561. * LN2_5 / 16384.) -
            y0 * (8991. * LN2_4 / 16384.) +
            y0 * (33939. * LN2_3 / 131072.) -
            y0 * (33255. * LN2_2 / 1048576.) -
            y0 * (15753. * LN2 / 4194304.);
  } else {
    return 0.;
    // expansion can be improved by adding higher order terms
  }
}

// ---------------------------------------------------------------------------
// ###########################################################################


// ---------------------------------------------------------------------------


/* *
  ================================ FUNCTIONAL ================================
 * */


// ###########################################################################
// ---------------------------------------------------------------------------
// functions of _x used in APS[x,y,z] parametrization
// ---------------------------------------------------------------------------

// s_aps(n, x) =
// \partial_n arctan (u*x / (1 + v*x))
double
s_aps (int dn, double _x, int id [])
{
    if (dn == 0) {
      return atan (_x * U_APS / (1. + _x * V_APS));
    } else if (dn > 0) {
      return tgamma (dn + 0.) / 2. *
        creal(I * cpow (+(U_APS - I * V_APS) / (I * (1. + V_APS * _x) - U_APS * _x), dn) -
              I * cpow (-(U_APS + I * V_APS) / (I * (1. + V_APS * _x) + U_APS * _x), dn));
    } else {
        return 0.00;
    }
}

// s_aps_2d1(n, x) =
// \partial_n [s_aps(1, x)]^2
// [used for derivatives of the inverse effective mass]
double
s_aps_2d1 (int dn, double _x, int id [])
{
    return pow_lrule(dn, _x, 2, 1, id, s_aps);
}

// BCS pairing gap formula
// =======================
// \partial_n [(8/e^2)*\exp(-\pi/(2*_x))]
double
pairing_bcs (int dn, double _x, int id [])
{
    double a = 8. / exp (2.); double b = -M_PI / 2.;
    int j, k; double factor_ = 0.;
    for (k = 0; k <= dn; k++) for (j = 0; j <= k; j++) {
        factor_ += pow (-1., j) * pow (b / _x, k) *
                rising_factorial (1. + j - k - dn, dn) /
                tgamma (1. + j) / tgamma (1. + k - j);
    }
    return a * exp (b / _x) * pow (_x, -dn) * factor_;
}

// Padé[1/1] approximation of the pairing field
// ============================================
// \parial_n [(1 + y*x) / (1 + p*y*x)]
double
pairing_fit (int dn, double _x, int id [])
{
    if (dn == 0) {
      return (1. + Y_APS * _x) / (1. + Z_APS * Y_APS * _x);
    } else if (dn > 0) {
      return tgamma (dn + 1.) * (1. - Z_APS) / (pow (Z_APS, 2) * Y_APS) *
            pow(-Z_APS * Y_APS / (1. + Z_APS * Y_APS * _x), dn + 1);
    } else {
      return 0.;
    }

}

// ---------------------------------------------------------------------------
// ###########################################################################




// ---------------------------------------------------------------------------




// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// usefull mathematical tools for calculation            [! DO NOT MODIFY >>>]
// of HFB parameters a, b, and 1/c                       [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// \partial_n [1 / \alpha_x]
double
one_over_a (int dn, double _x, int id [])
{
    return inverse_lrule (dn, _x, 0, id, inverse_effective_mass);
}

// \partial_n [1 / \eta_x]
double
one_over_h (int dn, double _x, int id [])
{
    return inverse_lrule (dn, _x, 0, id, pairing_gap);
}

// \partial_n [\eta_x^2]
double
h_sq (int dn, double _x, int id [])
{
    return pow_lrule (dn, _x, 2, 0, id, pairing_gap);
}

// \partial_n [\eta_x / \alpha_x]
double
h_over_a (int dn, double _x, int id [])
{
    return product_lrule (dn, _x, 0, 0, id, id, pairing_gap, one_over_a);
}

// \partial_n [\alpha_x / \eta_x]
double
a_over_h (int dn, double _x, int id [])
{
    return product_lrule (dn, _x, 0, 0, id, id, one_over_h, inverse_effective_mass);
}

// \partial_n [(\partial_k \eta_x / \alpha_x)^p]
double
h_over_a_pdk (int dn, double _x, int idpk []) // idpk = {id} + {p} + {k}
{
    int id [2] = {idpk [0],idpk [1]};
    int ip [1] = {idpk [2]}; int ik [1] = {idpk [3]};
    return pow_lrule (dn, _x, ip [0], ik [0], id, h_over_a);
}

// \partial_n \ln [\eta_x / \alpha_x]
double
log_h_over_a (int dn, double _x, int id [])
{
    if (dn == 0) {
      return log (h_over_a (0, _x, id));
    } else if (dn == 1) {
      return h_over_a (1, _x, id) * a_over_h (0, _x, id);
    } else if (dn > 1) {
      return product_lrule (dn - 1, _x, 1, 0, id, id, h_over_a, a_over_h);
    } else {
      return 0.;
    }
}

// (\partial / \partial _x)^n \mathcal{B}_p(y(x))
double
b_series_coefficient (int dn, double _x, int idp []) // idp = {id} + {p}
{
    int id [2] = {idp [0],idp [1]}; int ip [1] = {idp [2]};
    return composed_lrule (dn, _x, 0, 0, ip, id, b_expansion, log_h_over_a);
}

// (\partial / \partial _x)^n \mathcal{C}_p(y(x))
double
c_series_coefficient (int dn, double _x, int idp []) // idp = {id} + {p}
{
    int id [2] = {idp [0],idp [1]}; int ip [1] = {idp [2]};
    return composed_lrule (dn, _x, 0, 0, ip, id, c_expansion, log_h_over_a);
}

// \partial_n \sum_p \mathcal{B}_p(y(x)) (\eta_x / \alpha_x)^p
double
b_series (int dn, double _x, int id [])
{
    int idx [3] = {id [0],id [1],0}; int ip [4] = {id [0],id [1],0,0};
    int ij; double r_ = 0.;
    for(ij = 0; ij <= HFB_ORDER; ij++) {
      idx [2] = ij; ip [2] = ij; ip [3] = 0;
      r_ += product_lrule (dn, _x, 0, 0, idx, ip, b_series_coefficient, h_over_a_pdk);
    }
    return r_;
}

// \partial_n \sum_p \mathcal{C}_p(y(x)) (\eta_x / \alpha_x)^p
double
c_series (int dn, double _x, int id [])
{
    int idx[3] = {id [0],id [1],0}; int ip [4] = {id [0],id [1],0,0};
    int ij; double r_ = 0;
    for(ij = 0; ij <= HFB_ORDER; ij++) {
      idx [2] = ij; ip [2] = ij; ip [3] = 0;
      r_ += product_lrule (dn, _x, 0, 0, idx, ip, c_series_coefficient, h_over_a_pdk);
    }
    return r_;
}

// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################


// ---------------------------------------------------------------------------


/* *
  ================================ PARAMETERS ================================
 * */

// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// HFB parameters a, b, and 1/c                          [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// \partial_n a_x
double
a_hfb (int dn, double _x, int id [])
{
    return inverse_effective_mass (dn, _x, id);
}

// \partial_n [a_x * b_x]
double
b_hfb (int dn, double _x, int id [])
{

    if (pairing_gap (0, _x, id) == 0) {
      return -a_hfb (dn, _x, id); // BCS limit
      // avoid instability of derivative when vanishing pairing
    } else {
      return product_lrule (dn, _x, 0, 0, id, id, a_hfb, b_series);
    }

}

// \partial_n [c_x / a_x]
double
inverse_c_hfb (int dn, double _x, int id [])
{
    if (pairing_gap (0, _x, id) == 0.) {
      if (_x < DDCC_EPSILON) {
        return inverse_c_hfb (dn, DDCC_EPSILON, id);
      } else {
        int ip[1] = {-1};
        return -(M_PI / 8.) * product_lrule (dn, _x, 0, 0, id, ip, one_over_a, xm); // BCS limit
      }
      // avoid instability of derivative when vanishing pairing
    } else {
      return product_lrule (dn, _x, 0, 0, id, id, one_over_a, c_series);
    }
}

// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################


// ---------------------------------------------------------------------------


// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// kinetic density                                       [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// \partial_n [1 / b_x]
double
one_over_b (int dn, double _x, int id [])
{
    return inverse_lrule (dn, _x, 0, id, b_hfb);
}

// \partial_n [\eta_x / b_x]
double
h_over_b (int dn, double _x, int id [])
{
    return product_lrule (dn, _x, 0, 0, id, id, pairing_gap, one_over_b);
}

// \partial_n [b_x / a_x]
double
b_over_a (int dn, double _x, int id [])
{
    return product_lrule (dn, _x, 0, 0, id, id, b_hfb, one_over_a);
}

// \partial_n [|b_x / a_x|]
double
abs_b_over_a (int dn, double _x, int id [])
{
    return -product_lrule (dn, _x, 0, 0, id, id, b_hfb, one_over_a);
}

// \partial_n [|b_x / a_x|^{5/2}]
double
b_over_a_pow5over2 (int dn, double _x, int id []) // take absolute value
{
    // return realpow_lrule (dn, _x, 5./2., 0, id, abs_b_over_a); // too time consuming
    if (dn == 0) {
        return pow(abs_b_over_a(0, _x, id), 5./2.);
    } else if (dn == 1) {
        return 5./2. * abs_b_over_a(1, _x, id) * pow(abs_b_over_a(0, _x, id), 3./2.);
    } else if (dn == 2) {
        return 5./2. * abs_b_over_a(2, _x, id) * pow(abs_b_over_a(0, _x, id), 3./2.)
        + 5./2. * 3./2. * pow(abs_b_over_a(1, _x, id), 2) * pow(abs_b_over_a(0, _x, id), 1./2.);
    } else if (dn == 3) {
        return 5./2. * abs_b_over_a(3, _x, id) * pow(abs_b_over_a(0, _x, id), 3./2.)
        + 3. * 5./2. * 3./2. * abs_b_over_a(2, _x, id) * abs_b_over_a(1, _x, id) * pow(abs_b_over_a(0, _x, id), 1./2.)
        + 5./2. * 3./2. * 1./2. * pow(abs_b_over_a(1, _x, id), 3) / pow(abs_b_over_a(0, _x, id), 1./2.);
    } else {
        return 0.;
    }
}

// \partial_n [|eta_x / b_x|^2]
double
h_over_b_pow2 (int dn, double _x, int id [])
{
    return product_lrule (dn, _x, 0, 0, id, id, h_over_b, h_over_b);
}

// \partial_n [|eta_x / b_x|^4]
double
h_over_b_pow4 (int dn, double _x, int id [])
{
    return product_lrule (dn, _x, 0, 0, id, id, h_over_b_pow2, h_over_b_pow2);
}

// \partial_n [b_x / eta_x]
double
b_over_h (int dn, double _x, int id [])
{
    return product_lrule (dn, _x, 0, 0, id, id, b_hfb, one_over_h);
}

// \partial_n \ln [\eta_x / b_x]
double
log_h_over_b (int dn, double _x, int id []) // take absolute value
{
    if (dn == 0) {
      return log (-h_over_b (0, _x, id));
    } else if (dn == 1) {
      return h_over_b (1, _x, id) * b_over_h (0, _x, id);
    } else if (dn > 1) {
      return product_lrule (dn - 1, _x, 1, 0, id, id, h_over_b, b_over_h);
    } else {
      return 0.;
    }
}

// auxiliary function used to define integral_bar bellow
double
integral_bar_0 (int dn, double _x, int id [])
{
  double l_ = id[0]/2.;
  int ip0 [1] = {0}; double x0 = xm (dn, _x, ip0);
  return 2. * x0 - l_ / 2. * h_over_b_pow2 (dn, _x, id) * (l_-1.)
    + l_ / 32. * h_over_b_pow4 (dn, _x, id) * (l_-1.) * (l_-2.) * (l_-3.);
}

// auxiliary function used to define integral_bar bellow
double
integral_bar_1 (int dn, double _x, int id [])
{
  double l_ = id[0]/2.;
  int ip0 [1] = {0}; double x0 = xm (dn, _x, ip0);
  double H_;
  if (fabs(l_ - 1./2.) < 0.01) {
    H_ = 2. - 2. * LN2;
  } else if (fabs(l_ - 3./2.) < 0.01) {
    H_ = 8./3. - 2. * LN2;
  } else if (fabs(l_ - 5./2.) < 0.01) {
    H_ = 46./15. - 2. * LN2;
  } else {
    H_ = 0.;
  }
  return (H_ - LN2)*x0 + log_h_over_b (dn, _x, id);
}

// auxiliary function used to define integral_bar bellow
double
integral_bar_2 (int dn, double _x, int id [])
{
  double l_ = id[0]/2.;
  return - 1. / 2. * h_over_b_pow2 (dn, _x, id) * (1. - l_ - pow(l_,2))
         + 1. / 32. * h_over_b_pow4 (dn, _x, id) * (6. - 13.*l_ + 3./2.*pow(l_,2) + 5.*pow(l_,3) - 3./2.*pow(l_,4));
}

// see attachment for full definition (appendix A)
double
integral_bar (int dn, double _x, int id [])
{
  return product_lrule (dn, _x, 0, 0, id, id, integral_bar_0, integral_bar_1) + integral_bar_2 (dn, _x, id);
}

// see attacment for full definition (appendix A), \propto \tau_x
double
integral_bar_kin (int dn, double _x, int id [])
{
  int id5[1] = {5}; int id3[1] = {3};
  return (integral_bar (dn, _x, id5) - integral_bar (dn, _x, id3));
}


// \partial_n [\tau_x]
double
kin_density (int dn, double _x, int id [])
{
    return 3. / 2. * product_lrule (dn, _x, 0, 0, id, id, b_over_a_pow5over2, integral_bar_kin);
}

// \partial_n [_x * \partial_1 a_x]
double
x_times_ad1  (int dn, double _x, int id [])
{
    int ip1 [1] = {1};
    return product_lrule (dn, _x, 0, 1, ip1, id, xm, a_hfb);
}
// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################


// ---------------------------------------------------------------------------


// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// SLDA parameters alpha, beta, and 1/gamma              [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// \partial_n \alpha_x = \partial_n a_x
double
alpha_parameter (int dn, double _x, int id [])
{
    return a_hfb (dn, _x, id);
}


// \partial_n \gamma_x^{-1} = \partial_n [6/(3\pi^2)^{2/3} / c_x]
double
inverse_gamma_parameter (int dn, double _x, int id [])
{
    int ip[1] = {1};
    double cst_ = 6. / pow (3. * pow (M_PI, 2), 2. / 3.);
    return cst_ * inverse_c_hfb (dn, _x, id)
         + cst_ * product_lrule (dn, _x, 1, 0, id, ip, inverse_c_hfb, xm);

}

// \partial_n \beta_x = \partial_n [b_x + \zeta_x + \eta_x^2 / c_x - partial_1 (\alpha_x) * tau / 2]
double
beta_parameter (int dn, double _x, int id [])
{
    double pairing_correction, kin_correction;
    double cst_ = 6. / pow (3. * pow (M_PI, 2), 2. / 3.);
      if (pairing_gap (0, _x, id) == 0) {
        pairing_correction = 0.;
        kin_correction = 0.;
      } else {
        pairing_correction = product_lrule (dn, _x, 0, 0, id, id, h_sq, inverse_gamma_parameter) / cst_;
        kin_correction = 1./2. * 1./3. * product_lrule (dn, _x, 0, 0, id, id, x_times_ad1, kin_density);
      }
      return b_hfb (dn, _x, id) + chemical_potential (dn, _x, id) + pairing_correction - kin_correction;
}



// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// functional parameter A, B, and C                      [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// A_x = \alpha_x = a_x
double
a_functional (double _x, int id [])
{
    return alpha_parameter (0, _x, id);
}

// B_x = \sum_n (-x)^n (5!/(5+n)!) \partial_n \beta_x
double
b_functional (double _x, int id [])
{
    double d_ = 5.;
    int dn; double r_ = 0.;
    for (dn = 0; dn <= FUNCTIONAL_ORDER; dn++) {
      r_ += pow (-_x, dn) * beta_parameter (dn, _x, id) *
            tgamma (d_ + 1.) / tgamma (d_ + 1. + dn);
    }
    return r_;
}


// C_x \propto c_x
double
c_functional (double _x, int id [])
{
    double cst_ = 6. / pow (3. * pow (M_PI, 2), 2. / 3.);
    return 1. / (cst_ * inverse_c_hfb (0, _x, id));
}

// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################


// ---------------------------------------------------------------------------


// ###########################################################################
// ---------------------------------------------------------------------------
// interpolation of parameters
// ---------------------------------------------------------------------------
struct rparams {
    double *X;
    double *F;
    double F_UFG;
};

int
beta_fit (const gsl_vector * x, void *params, gsl_vector * f)
{
    const double a1 = gsl_vector_get (x, 0);
    const double a2 = gsl_vector_get (x, 1);
    const double a3 = gsl_vector_get (x, 2);
    const double b1 = gsl_vector_get (x, 3);
    const double b2 = gsl_vector_get (x, 4);
    const double b3 = gsl_vector_get (x, 5);
    const double c0 = gsl_vector_get (x, 6);

    double _x, _f; int n; double _Y[7];
    double _f_ufg = ((struct rparams *) params)->F_UFG;

    for (n = 0; n < 7; n++) {
        _x = ((struct rparams *) params)->X[n];
        _f = ((struct rparams *) params)->F[n];
        _Y[n] = _f / _f_ufg * (1.0 + (b1) * _x + (b2) * pow(_x, 2) + (b3) * pow(_x, 3) + (c0) * pow(_x, 4))
        - (0.0 + (a1) * _x + (a2) * pow(_x, 2) + (a3) * pow(_x, 3) + (c0) * pow(_x, 4));
        gsl_vector_set (f, n, _Y[n]);
    }
    return GSL_SUCCESS;
}

int
b_functional_fit (const gsl_vector * x, void *params, gsl_vector * f)
{
    const double a1 = gsl_vector_get (x, 0);
    const double a2 = gsl_vector_get (x, 1);
    const double a3 = gsl_vector_get (x, 2);
    const double b1 = gsl_vector_get (x, 3);
    const double b2 = gsl_vector_get (x, 4);
    const double b3 = gsl_vector_get (x, 5);
    const double c0 = gsl_vector_get (x, 6);

    double _x, _f; int n; double _Y[7];
    double _f_ufg = ((struct rparams *) params)->F_UFG;

    for (n = 0; n < 7; n++) {
        _x = ((struct rparams *) params)->X[n];
        _f = ((struct rparams *) params)->F[n];
        _Y[n] = _f / _f_ufg * (1.0 + (b1) * _x + (b2) * pow(_x, 2) + (b3) * pow(_x, 3) + (c0) * pow(_x, 4))
        - (0.0 + (a1) * _x + (a2) * pow(_x, 2) + (a3) * pow(_x, 3) + (c0) * pow(_x, 4));
        gsl_vector_set (f, n, _Y[n]);
    }
    return GSL_SUCCESS;
}


int
inverse_gamma_fit (const gsl_vector * x, void *params, gsl_vector * f)
{
    const double a1 = gsl_vector_get (x, 0);
    const double a2 = gsl_vector_get (x, 1);
    const double a3 = gsl_vector_get (x, 2);
    const double b1 = gsl_vector_get (x, 3);
    const double b2 = gsl_vector_get (x, 4);
    const double b3 = gsl_vector_get (x, 5);
    const double c0 = gsl_vector_get (x, 6);

    double _x, _f; int n;
    double _f_ufg = ((struct rparams *) params)->F_UFG;
    double _Y[7];

    for (n = 0; n < 7; n++) {
        _x = ((struct rparams *) params)->X[n];
        _f = ((struct rparams *) params)->F[n];
        _Y[n] = _f / _f_ufg * (1.0 + (b1) * _x + (b2) * pow(_x, 2) + (b3) * pow(_x, 3) + (c0) * pow(_x, 4))
        - (0.0 + (a1) * _x + (a2) * pow(_x, 2) + (a3) * pow(_x, 3) + (c0) * pow(_x, 4));
        gsl_vector_set (f, n, _Y[n]);
    }
    return GSL_SUCCESS;
}

int
c_functional_fit (const gsl_vector * x, void *params, gsl_vector * f)
{
    const double a1 = gsl_vector_get (x, 0);
    const double a2 = gsl_vector_get (x, 1);
    const double a3 = gsl_vector_get (x, 2);
    const double b1 = gsl_vector_get (x, 3);
    const double b2 = gsl_vector_get (x, 4);
    const double b3 = gsl_vector_get (x, 5);
    const double c0 = gsl_vector_get (x, 6);

    double _x, _f; int n;
    double _f_ufg = ((struct rparams *) params)->F_UFG;
    double _Y[7];

    for (n = 0; n < 7; n++) {
        _x = ((struct rparams *) params)->X[n];
        _f = ((struct rparams *) params)->F[n];
        _Y[n] = _f / _f_ufg * (1.0 + (b1) * _x + (b2) * pow(_x, 2) + (b3) * pow(_x, 3) + (c0) * pow(_x, 4))
        - (0.0 + (a1) * _x + (a2) * pow(_x, 2) + (a3) * pow(_x, 3) + (c0) * pow(_x, 4));
        gsl_vector_set (f, n, _Y[n]);
    }
    return GSL_SUCCESS;
}


// ---------------------------------------------------------------------------


void print_sldae_parameters (double _x, int id [])
{
    printf("\n_x = %f\n\n", _x);
    printf("\t xi \t\t ms \t\t zeta \t\t eta\n");
    printf("\t %+lf \t %+lf \t %+lf \t %+lf\n",
           ground_state_energy(0, _x, id),
           inverse_effective_mass(0, _x, id),
           chemical_potential(0, _x, id),
           pairing_gap(0, _x, id));
    printf("\t %+lf \t %+lf \t %+lf \t %+lf\n",
           ground_state_energy(1, _x, id),
           inverse_effective_mass(1, _x, id),
           chemical_potential(1, _x, id),
           pairing_gap(1, _x, id));
    printf("\t %+lf \t %+lf \t %+lf \t %+lf\n",
           ground_state_energy(2, _x, id),
           inverse_effective_mass(2, _x, id),
           chemical_potential(2, _x, id),
           pairing_gap(2, _x, id));

    printf("\n");

    printf("\t a \t\t\t b \t\t\t inverse_c\n");
    printf("\t %+lf \t %+lf \t %+lf\n",
           a_hfb(0, _x, id),
           b_hfb(0, _x, id),
           inverse_c_hfb(0, _x, id));
    printf("\t %+lf \t %+lf \t %+lf\n",
           a_hfb(1, _x, id),
           b_hfb(1, _x, id),
           inverse_c_hfb(1, _x, id));
    printf("\t %+lf \t %+lf \t %+lf\n",
           a_hfb(2, _x, id),
           b_hfb(2, _x, id),
           inverse_c_hfb(2, _x, id));

    printf("\n");

    printf("\t alpha \t\t beta \t\t inverse_gamma\n");
    printf("\t %+lf \t %+lf \t %+lf\n",
           alpha_parameter(0, _x, id),
           beta_parameter(0, _x, id),
           inverse_gamma_parameter(0, _x, id));
    printf("\t %+lf \t %+lf \t %+lf\n",
           alpha_parameter(1, _x, id),
           beta_parameter(1, _x, id),
           inverse_gamma_parameter(1, _x, id));
    printf("\t %+lf \t %+lf \t %+lf\n",
           alpha_parameter(2, _x, id),
           beta_parameter(2, _x, id),
           inverse_gamma_parameter(2, _x, id));

    printf("\n");

    printf("\t A \t\t\t B \t\t\t C\n");
    printf("\t %+lf \t %+lf \t %+lf\n",
           a_functional(_x, id),
           b_functional(_x, id),
           c_functional(_x, id));

    printf("\n");

    printf("\t A_UFG \t\t B_UFG \t\t C_UFG\n");
    printf("\t %+lf \t %+lf \t %+lf\n", a_functional(1e9, id), b_functional(1e9, id), c_functional(1e9, id));

    printf("\n");

}

void print_interpolation_parameters (gsl_multiroot_fsolver *s, char *prefix, FILE *f)
{
    int iter;
    for (iter = 0; iter < 100; iter++) {
        if (gsl_multiroot_fsolver_iterate (s)) break;
    }
    printf ("\t iter = %3u \t f(x) = % .3e \n",
                iter,
                fabs(gsl_vector_get (s->f, 0)) +
                fabs(gsl_vector_get (s->f, 1)) +
                fabs(gsl_vector_get (s->f, 2)) +
                fabs(gsl_vector_get (s->f, 3)) +
                fabs(gsl_vector_get (s->f, 4)) +
                fabs(gsl_vector_get (s->f, 5)) +
                fabs(gsl_vector_get (s->f, 6)));
    printf("\t a1 = %+lf;\n", (gsl_vector_get (s->x, 0)));
    printf("\t a2 = %+lf;\n", (gsl_vector_get (s->x, 1)));
    printf("\t a3 = %+lf;\n", (gsl_vector_get (s->x, 2)));
    printf("\t b1 = %+lf;\n", (gsl_vector_get (s->x, 3)));
    printf("\t b2 = %+lf;\n", (gsl_vector_get (s->x, 4)));
    printf("\t b3 = %+lf;\n", (gsl_vector_get (s->x, 5)));
    printf("\t c  = %+lf;\n", (gsl_vector_get (s->x, 6)));
    printf("\n");
    
    // write to file
    fprintf(f,"#define SLDAE_%s_a1 %+lf\n", prefix, (gsl_vector_get (s->x, 0)));
    fprintf(f,"#define SLDAE_%s_a2 %+lf\n", prefix, (gsl_vector_get (s->x, 1)));
    fprintf(f,"#define SLDAE_%s_a3 %+lf\n", prefix, (gsl_vector_get (s->x, 2)));
    fprintf(f,"#define SLDAE_%s_b1 %+lf\n", prefix, (gsl_vector_get (s->x, 3)));
    fprintf(f,"#define SLDAE_%s_b2 %+lf\n", prefix, (gsl_vector_get (s->x, 4)));
    fprintf(f,"#define SLDAE_%s_b3 %+lf\n", prefix, (gsl_vector_get (s->x, 5)));
    fprintf(f,"#define SLDAE_%s_c  %+lf\n", prefix, (gsl_vector_get (s->x, 6)));
    
}

double
pade_function (gsl_multiroot_fsolver *s, double lambda, double ufg_val)
{
    const double a1 = (gsl_vector_get (s->x, 0));
    const double a2 = (gsl_vector_get (s->x, 1));
    const double a3 = (gsl_vector_get (s->x, 2));
    const double b1 = (gsl_vector_get (s->x, 3));
    const double b2 = (gsl_vector_get (s->x, 4));
    const double b3 = (gsl_vector_get (s->x, 5));
    const double c0 = (gsl_vector_get (s->x, 6));

    double _x=lambda; 
    return ufg_val * 
        (0.0 + (a1) * _x + (a2) * pow(_x, 2) + (a3) * pow(_x, 3) + (c0) * pow(_x, 4))
            /
        (1.0 + (b1) * _x + (b2) * pow(_x, 2) + (b3) * pow(_x, 3) + (c0) * pow(_x, 4));
}
// ---------------------------------------------------------------------------
// ###########################################################################


// ---------------------------------------------------------------------------


// ###########################################################################
// ---------------------------------------------------------------------------
// main interpolation (padé[4/4]) routine of sldae functional
// ---------------------------------------------------------------------------
int main() {

    int id[2] = {0,0}; // choose sldae functional
//     double _x = 1.000;
//     print_sldae_parameters (_x, id); // print sldae parameters at _x

    // text file for checking
    double *lamarr = (double *)malloc(sizeof(double)*1000);
    double lam, ufg_val;
    int lamidx=0, idx;
    for(lam=0.01; lam<0.10; lam+=0.01) {lamarr[lamidx]=lam; lamidx++;}
    for(lam=0.10; lam<1.0 ; lam+=0.1) {lamarr[lamidx]=lam; lamidx++;}
    for(lam= 1.0; lam<10.0; lam+=1.0) {lamarr[lamidx]=lam; lamidx++;}
    for(lam=10.0; lam<100; lam+=10) {lamarr[lamidx]=lam; lamidx++;}
    
    
    printf(" SLDAE DEV TOOLS\n");
    printf(" [padé[4/4] interpolation of sldae parameters p = {beta, gamma, B, C}]\n");
    printf("    p       a1 x + a2 x^2 + a3 x^3 + c x^4  \n");
    printf("  ----- = ----------------------------------\n");
    printf("  p_ufg   1 + b1 x + b2 x^2 + b3 x^3 + c x^4\n");
    printf("\n");


    // Create empty file
    char fname[128]="sldae_parameters.h";
    FILE * out_h = fopen(fname, "w");
    
    // initialization of interpolation routine
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (T, 7);

    // interpolation points (7 requiered)
    double X[7] = {0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0};
    double F_beta[7] = {
        beta_parameter(0, X[0], id),
        beta_parameter(0, X[1], id),
        beta_parameter(0, X[2], id),
        beta_parameter(0, X[3], id),
        beta_parameter(0, X[4], id),
        beta_parameter(0, X[5], id),
        beta_parameter(0, X[6], id)
    };
    double F_inverse_gamma[7] = {
        inverse_gamma_parameter(0, X[0], id),
        inverse_gamma_parameter(0, X[1], id),
        inverse_gamma_parameter(0, X[2], id),
        inverse_gamma_parameter(0, X[3], id),
        inverse_gamma_parameter(0, X[4], id),
        inverse_gamma_parameter(0, X[5], id),
        inverse_gamma_parameter(0, X[6], id)
    };
    double F_bf[7] = {
        b_functional(X[0], id),
        b_functional(X[1], id),
        b_functional(X[2], id),
        b_functional(X[3], id),
        b_functional(X[4], id),
        b_functional(X[5], id),
        b_functional(X[6], id)
    };
    double F_cf[7] = {
        c_functional(X[0], id),
        c_functional(X[1], id),
        c_functional(X[2], id),
        c_functional(X[3], id),
        c_functional(X[4], id),
        c_functional(X[5], id),
        c_functional(X[6], id)
    };
    struct rparams p_beta = {X, F_beta, beta_parameter(0, 1e9, id)};
    struct rparams p_inverse_gamma = {X, F_inverse_gamma, inverse_gamma_parameter(0, 1e9, id)};
    struct rparams p_bf = {X, F_bf, b_functional(1e9, id)};
    struct rparams p_cf = {X, F_cf, c_functional(1e9, id)};

    // initialization of the solution
    gsl_vector *x = gsl_vector_alloc (7);
    gsl_vector_set (x, 0, 1);
    gsl_vector_set (x, 1, 1);
    gsl_vector_set (x, 2, 1);
    gsl_vector_set (x, 3, 1);
    gsl_vector_set (x, 4, 1);
    gsl_vector_set (x, 5, 1);
    gsl_vector_set (x, 6, 1);

    // interpolation of beta_parameters
    gsl_multiroot_function beta = {&beta_fit, 7, &p_beta};
    gsl_multiroot_fsolver_set (s, &beta, x);
    ufg_val=beta_parameter(0, 1e9, id);
    printf ("\n beta interpolation (beta_ufg = %+lf)\n", ufg_val);
    print_interpolation_parameters (s, "beta", out_h);
    double *beta_fit_arr=(double *)malloc(sizeof(double)*lamidx);
    for(idx=0; idx<lamidx; idx++) beta_fit_arr[idx]=pade_function (s,  lamarr[idx], ufg_val);

    // interpolation of b_functional
    gsl_multiroot_function bf = {&b_functional_fit, 7, &p_bf};
    gsl_multiroot_fsolver_set (s, &bf, x);
    ufg_val=b_functional(1e9, id);
    printf ("\n b_functional interpolation (b_functional_ufg = %+lf)\n", ufg_val);
    print_interpolation_parameters (s, "b_functional", out_h);
    double *b_functional_fit_arr=(double *)malloc(sizeof(double)*lamidx);
    for(idx=0; idx<lamidx; idx++) b_functional_fit_arr[idx]=pade_function (s,  lamarr[idx], ufg_val);

    // interpolation of inverse_gamma_parameters
    gsl_multiroot_function inverse_gamma = {&inverse_gamma_fit, 7, &p_inverse_gamma};
    gsl_multiroot_fsolver_set (s, &inverse_gamma, x);
    ufg_val=inverse_gamma_parameter(0, 1e9, id);
    printf ("\n inverse_gamma interpolation (inverse_gamma_ufg = %+lf)\n", ufg_val);
    print_interpolation_parameters (s, "inverse_gamma", out_h);
    double *inverse_gamma_fit_arr=(double *)malloc(sizeof(double)*lamidx);
    for(idx=0; idx<lamidx; idx++) inverse_gamma_fit_arr[idx]=pade_function (s,  lamarr[idx], ufg_val);

    // interpolation of c_functional
    gsl_multiroot_function cf = {&c_functional_fit, 7, &p_cf};
    gsl_multiroot_fsolver_set (s, &cf, x);
    ufg_val=c_functional(1e9, id);
    printf ("\n c_functional interpolation (c_functional_ufg = %+lf)\n", ufg_val);
    print_interpolation_parameters (s, "c_functional", out_h);
    double *c_functional_fit_arr=(double *)malloc(sizeof(double)*lamidx);
    for(idx=0; idx<lamidx; idx++) c_functional_fit_arr[idx]=pade_function (s,  lamarr[idx], ufg_val);

    // reset memory of solver
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    
    printf("Parameters written to: %s\n", fname);
    fclose(out_h);
        
    printf("Creating file: sldae_parameters.txt\n");
    FILE * txt = fopen("sldae_parameters.txt", "w");
    fprintf(txt, "# 1: lambda\n");
    fprintf(txt, "# 2: ground_state_energy\n");
    fprintf(txt, "# 3: inverse_effective_mass\n");
    fprintf(txt, "# 4: chemical_potential\n");
    fprintf(txt, "# 5: pairing_gap\n");
    fprintf(txt, "# 6: beta_parameter\n");
    fprintf(txt, "# 7: beta_parameter - PADE\n");
    fprintf(txt, "# 8: b_functional\n");
    fprintf(txt, "# 9: b_functional - PADE\n");
    fprintf(txt, "#10: inverse_gamma_parameter\n");
    fprintf(txt, "#11: inverse_gamma_parameter - PADE\n");
    fprintf(txt, "#12: c_functional\n");
    fprintf(txt, "#13: c_functional - PADE\n");
    for(idx=0; idx<lamidx; idx++)
    {
        fprintf(txt, "%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
            lamarr[idx],            
            ground_state_energy(0, lamarr[idx], id),
            inverse_effective_mass(0, lamarr[idx], id),
            chemical_potential(0, lamarr[idx], id),
            pairing_gap(0, lamarr[idx], id),
            beta_parameter(0, lamarr[idx], id),
            beta_fit_arr[idx],
            b_functional(lamarr[idx], id),
            b_functional_fit_arr[idx],
            inverse_gamma_parameter(0, lamarr[idx], id),
            inverse_gamma_fit_arr[idx],
            c_functional(lamarr[idx], id),
            c_functional_fit_arr[idx]
        );
        fflush(txt);
    }
    fclose(txt);
  return 0;
}
// ---------------------------------------------------------------------------
// ###########################################################################
