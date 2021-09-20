#include "sldae_parameters.h"
#include "sldae_functional.h"
#include "sldae_tools.h" // -> called in aps_functional.h (I don't know if it is required)
// ---------------------------------------------------------------------------



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
