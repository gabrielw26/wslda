#include "sldae_functional.h"
#include "sldae_tools.h"
// ---------------------------------------------------------------------------



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
    return tgamma (dn) / 2. *
           creal(I * pow (+(U_APS - I * V_APS) / (I * (1. + V_APS * _x) - U_APS * _x), dn) -
                 I * pow (-(U_APS + I * V_APS) / (I * (1. + V_APS * _x) + U_APS * _x), dn));
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

  return x0 + A_APS *
              product_lrule (dn, _x, 0, 1, ip2, id, xm, s_aps) +
              B_APS *
              product_lrule (dn, _x, 0, 0, ip2, id, xm, s_aps_2d1) +
              B_APS * V_APS *
              product_lrule (dn, _x, 0, 0, ip3, id, xm, s_aps_2d1);
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
