/**
 * W-SLDA Toolkit
 * Author: Antoine Boulet
 * Creation date: 2021.09.21
 * */



#include "sldae_functional.h"




/**
    //////////////////////////////////////////////////////////////////////
    ===  IMPLEMENTATION OF FUNCTIONS FOR THE SLDA extended FUNCTIONAL  ===
    //////////////////////////////////////////////////////////////////////

  # -- NOTATIONS AND CONVENTIONS
       =========================
       FUNCTIONAL_ID / PAIRING_ID: select functional
       Note: only APS[x,y,z] is implemented

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
  User can define its own functional through the custom_* functions
  using FUNCTIONAL_ID = -1 and/or PAIRING_ID = -1 (to be done)

  * remark 1 *
  The rest of the function provide automatic results for the SLDA
  extended methods as described in the notes.
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
  unable to provide higher order terms).

  * remark 2 *
  The derivatives (of product, compostition, etc.) are obtained using
  general Leibniz rules such that, only the derivatives of the required
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
**/
// ---------------------------------------------------------------------------


/**
  ================================== TOOLS ===================================
**/

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
  return tgamma (in + 1) / tgamma (ik + 1) / tgamma (in - ik + 1);
}

// exponential integral
// ====================
// Ei(x) = - \int_{-x}^\infty dt e^{-t} / t
double
exponential_integral (double x)
{
  // return std::expint(x); // does not work...
  return 0.;
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
// [implmented following the convention: ip[0] correspond to the power of x]
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

// power Leibniz rule
// ==================
// \parial_n [f^(k)(x)]^i =
// \sum_{p=0}^n C_{n,p} f^(k+p) \partial_{n-p} [f^(k)(x)]^{i-1}
double
power_lrule
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
            power_lrule (dn - dp, _x, i - 1, dk, id, (*f));
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
      r_ += pow(-1, dp) * binomial_coefficient (dn, dp) *
            (dn + 1.) / (dp + 1.) *
            power_lrule (dn, _x, dp, dk, id, (*f)) /
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
    int dj; double dg [dn];
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
  int ip5 [1] = {5}; double y5 = xm (dp, y_x, ip5);
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


/**
  ================================ FUNCTIONAL ================================
**/


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
    double complex z_one = 1.00 + I * 0.00;
    return tgamma (dn) / 2. *
           creal(I * cpow (+(U_APS - I * V_APS) / (I * (1. + V_APS * _x) - U_APS * _x), dn) -
                 I * cpow (-(U_APS + I * V_APS) / (I * (1. + V_APS * _x) + U_APS * _x), dn));
  } else {
    return 0.;
  }
}

// s_aps_2d1(n, x) =
// \partial_n [s_aps(1, x)]^2
// [used for derivatives of the inverse effective mass]
double
s_aps_2d1 (int dn, double _x, int id [])
{
  return power_lrule(dn, _x, 2, 1, id, s_aps);
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
    factor_ += pow (-1, j) * pow (b / _x, k) *
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
    return tgamma (dn + 1) * (1. - Z_APS) / (pow (Z_APS, 2) * Y_APS) *
           pow(-Z_APS * Y_APS / (1. + Z_APS * Y_APS * _x), dn + 1);
  } else {
    return 0.;
  }
}

// ---------------------------------------------------------------------------
// ###########################################################################



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
  int ip1 [1] = {1}; double x1 = xm (dn, _x, ip1);
  int ip2 [1] = {2}; double x2 = xm (dn, _x, ip2);
  // ...
  int fid = id [0];       // FUNCTIONAL_ID
  int pid = id [1];       // PAIRING_ID

  double cst0_ = 16. / (3. * M_PI);
  return x0 - cst0_ * s_aps (dn, _x, id);
}

// chemical potential
// ==================
// \parial_n \zeta_x
double
chemical_potential (int dn, double _x, int id [])
{
  int ip0 [1] = {0}; double x0 = xm (dn, _x, ip0);
  int ip1 [1] = {1}; double x1 = xm (dn, _x, ip1);
  int ip2 [1] = {2}; double x2 = xm (dn, _x, ip2);
  // ...
  int fid = id [0];       // FUNCTIONAL_ID
  int pid = id [1];       // PAIRING_ID

  double cst0_ = 16. / (3. * M_PI);
  double cst1_ = 16. / (15. * M_PI);
  return x0 - cst0_ * s_aps (dn, _x, id) -
              cst1_ * product_lrule (dn, _x, 0, 1, ip1, id, xm, s_aps);
}

// inverse effective mass
// ======================
// \parial_n \alpha_x
double
inverse_effective_mass (int dn, double _x, int id [])
{
  int ip0 [1] = {0}; double x0 = xm (dn, _x, ip0);
  int ip1 [1] = {1}; double x1 = xm (dn, _x, ip1);
  int ip2 [1] = {2}; double x2 = xm (dn, _x, ip2);
  int ip3 [1] = {3}; double x3 = xm (dn, _x, ip3);
  // ...
  int fid = id [0];       // FUNCTIONAL_ID
  int pid = id [1];       // PAIRING_ID

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
  int ip0 [1] = {0}; double x0 = xm (dn, _x, ip0);
  int ip1 [1] = {1}; double x1 = xm (dn, _x, ip1);
  int ip2 [1] = {2}; double x2 = xm (dn, _x, ip2);
  // ...
  int fid = id [0];       // FUNCTIONAL_ID
  int pid = id [1];       // PAIRING_ID

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
  return power_lrule (dn, _x, 2, 0, id, pairing_gap);
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
  return power_lrule (dn, _x, ip [0], ik [0], id, h_over_a);
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

  int ij; double r_ = 0;
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


/**
  ================================ PARAMETERS ================================
**/

// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// HFB parameters a, b, and 1/c                          [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// \partial_n a_x
double
a_hfb (int dn, double _x, int id [])
{
  int fid = id [0];       // FUNCTIONAL_ID
  int pid = id [1];       // PAIRING_ID

  return inverse_effective_mass (dn, _x, id);
}

// \partial_n [a_x * b_x]
double
b_hfb (int dn, double _x, int id [])
{
  int fid = id [0];       // FUNCTIONAL_ID
  int pid = id [1];       // PAIRING_ID

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
  int fid = id [0];       // FUNCTIONAL_ID
  int pid = id [1];       // PAIRING_ID

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
// SLDA parameters alpha, beta, and 1/gamma              [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// \partial_n \alpha_x = \partial_n a_x
double
alpha_parameter (int dn, double _x, int id [])
{
  int fid = id [0]; // FUNCTIONAL_ID
  int pid = id [1]; // PAIRING_ID

  return a_hfb (dn, _x, id);
}

// \partial_n \beta_x = \partial_n [b_x + \zeta_x + \eta_x^2 / c_x]
double
beta_parameter (int dn, double _x, int id [])
{
  int fid = id [0]; // FUNCTIONAL_ID
  int pid = id [1]; // PAIRING_ID

  double pairing_correction;
  if (pairing_gap (0, _x, id) == 0) {
    pairing_correction = 0.;
  } else {
    pairing_correction = product_lrule (dn, _x, 0, 0, id, id, h_sq, inverse_c_hfb);
  }
  return b_hfb (dn, _x, id) +
    chemical_potential (dn, _x, id) +
    pairing_correction;
}

// \partial_n \gamma_x^{-1} = \partial_n [6/(3\pi^2)^{2/3} / c_x]
double
inverse_gamma_parameter (int dn, double _x, int id [])
{
  int fid = id [0]; // FUNCTIONAL_ID
  int pid = id [1]; // PAIRING_ID

  double cst_ = 6. / pow (3. * pow (M_PI, 2), 2. / 3.);
  return cst_ * inverse_c_hfb (dn, _x, id);
}

// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// functional parameter A, B, and C                      [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// A_x = \sum_n (-x)^n (3!/(3+n)!) \partial_n \alpha_x
double
a_functional (double _x, int id [])
{
  int fid = id [0]; // FUNCTIONAL_ID
  int pid = id [1]; // PAIRING_ID

  double d_ = 3.;
  int dn; double r_ = 0.;
  for (dn = 0; dn <= FUNCTIONAL_ORDER; dn++) {
    r_ += pow (-_x, dn) * alpha_parameter (dn, _x, id) *
          tgamma (d_ + 1) / tgamma (d_ + 1 + dn);
    }
  return r_;
}

// B_x = \sum_n (-x)^n (5!/(5+n)!) \partial_n \beta_x
double
b_functional (double _x, int id [])
{
  int fid = id [0]; // FUNCTIONAL_ID
  int pid = id [1]; // PAIRING_ID

  double d_ = 5.;
  int dn; double r_ = 0.;
  for (dn = 0; dn <= FUNCTIONAL_ORDER; dn++) {
    r_ += pow (-_x, dn) * beta_parameter (dn, _x, id) *
          tgamma (d_ + 1) / tgamma (d_ + 1 + dn);
  }
  return r_;
}

// C_x^{-1} = \sum_n (-x)^n (1!/(1+n)!) \partial_n [1 / \gamma_x]
double
c_functional (double _x, int id [])
{
  int fid = id [0]; // FUNCTIONAL_ID
  int pid = id [1]; // PAIRING_ID

  double d_ = 1.;
  if (_x == 0.) {
    return 0.;
  } else {
    int n; double r_ = 0.;
    for (n = 0; n <= FUNCTIONAL_ORDER; n++) {
      r_ += pow (-_x, n) * inverse_gamma_parameter (n, _x, id) *
            tgamma (d_ + 1) / tgamma (d_ + 1 + n);
    }
    return 1 / r_;
  }
}

// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################



// ---------------------------------------------------------------------------



/// ###########################################################################
// ---------------------------------------------------------------------------
// pairing coupling constant regularization routine
// ---------------------------------------------------------------------------
extern double dc_ec;

double
pcc_renormalization(double _x, double lmu, int id [])
{
  double alpha_ = inverse_effective_mass (0, _x, id);
  //    alpha k_0^2 / 2 - lmu = E_0
  //    alpha k_c^2 / 2 - lmu = E_c

  double k0, kc, lambda_;
  k0 = sqrt (fabs (2. * (0.000 + lmu) / alpha_));
  kc = sqrt (fabs (2. * (dc_ec + lmu) / alpha_));
  if (lmu >= 0) {
    lambda_ = (kc + k0) / (kc - k0);
    lambda_ = 1. - k0 / (2. * kc) * log(lambda_);
    lambda_ *= kc / (2. * M_PI_SQ);
  } else {
    lambda_ = k0 / kc;
    lambda_ = 1. + k0 / kc * atan(lambda_);
    lambda_ *= kc / (2. * M_PI_SQ);
  }
  return lambda_;
}

// ---------------------------------------------------------------------------
// ###########################################################################
