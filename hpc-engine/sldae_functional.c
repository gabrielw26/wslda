/**
 * W-SLDA Toolkit
 * Author: Antoine Boulet
 * Creation date: 2021.09.21
 *
 * td update (AB): 2021.11.25
 **/

#include "sldae_functional.h"

// ###########################################################################
// ---------------------------------------------------------------------------
// s_aps_dx function
// ==========================
FDECORATOR double
s_aps_d0 (double _x)
{
  return atan (_x * U_APS / (1. + _x * V_APS));
}
FDECORATOR double
s_aps_d1 (double _x)
{
  double U2plusV2 = (pow(U_APS, 2) + pow(V_APS, 2));
  double tmp_denominator = (1. + 2. * V_APS * _x + U2plusV2 * pow(_x, 2));
  return U_APS / tmp_denominator;
}
FDECORATOR double
s_aps_d2 (double _x)
{
  double U2plusV2 = (pow(U_APS, 2) + pow(V_APS, 2));
  double tmp_denominator = (1. + 2. * V_APS * _x + U2plusV2 * pow(_x, 2));
  return -2. * U_APS * (V_APS + U2plusV2 * _x) / pow(tmp_denominator, 2);
}
FDECORATOR double
s_aps_d3 (double _x)
{
  double U2plusV2 = (pow(U_APS, 2) + pow(V_APS, 2));
  double tmp_denominator = (1. + 2. * V_APS * _x + U2plusV2 * pow(_x, 2));
  return (6. * pow(U_APS, 5) * pow(_x, 2)
          + 6. * U_APS * pow(V_APS, 2) * pow(1. + V_APS * _x, 2)
          + 2. * pow(U_APS, 3) * (-1. + 6. * V_APS * _x * (1. + V_APS * _x)))
        / pow(tmp_denominator, 3);
}

// physical quantites
// ==========================
FDECORATOR double
ground_state_energy_d0 (double _x)
{
  return 1. - 16. / (3. * M_PI) * s_aps_d0 (_x);
}
FDECORATOR double
chemical_potential_d0 (double _x)
{
  return 1. - 16. / (3. * M_PI) * s_aps_d0 (_x)
            - 16. / (15. * M_PI) * _x * s_aps_d1 (_x);
}
FDECORATOR double
inverse_effective_mass_d0 (double _x)
{
  return 1. + A_APS * pow(_x, 2) * s_aps_d1 (_x)
            + B_APS * (1. + V_APS * _x) * pow(_x * s_aps_d1 (_x), 2);
}
FDECORATOR double
pairing_gap_d0 (double _x)
{
  return 8. / exp (2.) * exp (-M_PI / (2. * _x)) * (1. + Y_APS * _x) / (1. + Z_APS * Y_APS * _x);
}

// Padé[4/4] fitting function
// ==========================
FDECORATOR double
pade_d0 (double _x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4)
{
  return (a1*_x + a2*pow(_x,2) + a3*pow(_x,3) + a4*pow(_x,4))/ (1. + b1*_x + b2*pow(_x,2) + b3*pow(_x,3) + b4*pow(_x,4));
}

FDECORATOR double
pade_d1 (double _x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4)
{
  return -(((b1 + 2.*b2*_x + 3.*b3*pow(_x,2) + 4.*b4*pow(_x,3))* (a1*_x + a2*pow(_x,2) + a3*pow(_x,3) + a4*pow(_x,4)))/ pow(1. + b1*_x + b2*pow(_x,2) + b3*pow(_x,3) + b4*pow(_x,4),2)) + (a1 + 2.*a2*_x + 3.*a3*pow(_x,2) + 4.*a4*pow(_x,3))/ (1. + b1*_x + b2*pow(_x,2) + b3*pow(_x,3) + b4*pow(_x,4));
}

FDECORATOR double
pade_d2 (double _x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4)
{
  return (-2.*(a1 + 2.*a2*_x + 3.*a3*pow(_x,2) + 4.*a4*pow(_x,3))* (b1 + 2.*b2*_x + 3.*b3*pow(_x,2) + 4.*b4*pow(_x,3)))/ pow(1. + b1*_x + b2*pow(_x,2) + b3*pow(_x,3) + b4*pow(_x,4),2) + (2.*a2 + 6.*a3*_x + 12.*a4*pow(_x,2))/ (1. + b1*_x + b2*pow(_x,2) + b3*pow(_x,3) + b4*pow(_x,4)) + (a1*_x + a2*pow(_x,2) + a3*pow(_x,3) + a4*pow(_x,4))* ((2.*pow(b1 + 2.*b2*_x + 3.*b3*pow(_x,2) + 4.*b4*pow(_x,3),2))/ pow(1. + b1*_x + b2*pow(_x,2) + b3*pow(_x,3) + b4*pow(_x,4),3) - (2.*b2 + 6.*b3*_x + 12.*b4*pow(_x,2))/ pow(1. + b1*_x + b2*pow(_x,2) + b3*pow(_x,3) + b4*pow(_x,4),2));
}
// ---------------------------------------------------------------------------
// ###########################################################################


// ---------------------------------------------------------------------------


// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// SLDA parameters alpha, beta, and 1/gamma              [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// \partial_n \alpha_x = \partial_n a_x
FDECORATOR double
alpha_parameter_d0 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  a1 = 0.00000000; // = 0.;
  a2 = -0.2081560; // = U_APS * (A_APS + B_APS * U_APS);
  a3 = -0.1358830; // = U_APS * (2. * A_APS * V_APS + B_APS * U_APS * V_APS);
  a4 = -0.0162071; // = A_APS * U_APS * (pow(U_APS, 2) + pow(V_APS, 2));
  b1 = 2.09838000; // = 4. * V_APS;
  b2 = 1.73800000; // = (2. * pow(U_APS, 2) + 6. * pow(V_APS, 2));
  b3 = 0.66854900; // = 4. * V_APS * (pow(U_APS, 2) + pow(V_APS, 2));
  b4 = 0.10150800; // = pow(pow(U_APS, 2) + pow(V_APS, 2), 2);
  return 1. + pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
alpha_parameter_d1 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  a1 = 0.00000000; // = 0.;
  a2 = -0.2081560; // = U_APS * (A_APS + B_APS * U_APS);
  a3 = -0.1358830; // = U_APS * (2. * A_APS * V_APS + B_APS * U_APS * V_APS);
  a4 = -0.0162071; // = A_APS * U_APS * (pow(U_APS, 2) + pow(V_APS, 2));
  b1 = 2.09838000; // = 4. * V_APS;
  b2 = 1.73800000; // = (2. * pow(U_APS, 2) + 6. * pow(V_APS, 2));
  b3 = 0.66854900; // = 4. * V_APS * (pow(U_APS, 2) + pow(V_APS, 2));
  b4 = 0.10150800; // = pow(pow(U_APS, 2) + pow(V_APS, 2), 2);
  return pade_d1 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
alpha_parameter_d2 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  a1 = 0.00000000; // = 0.;
  a2 = -0.2081560; // = U_APS * (A_APS + B_APS * U_APS);
  a3 = -0.1358830; // = U_APS * (2. * A_APS * V_APS + B_APS * U_APS * V_APS);
  a4 = -0.0162071; // = A_APS * U_APS * (pow(U_APS, 2) + pow(V_APS, 2));
  b1 = 2.09838000; // = 4. * V_APS;
  b2 = 1.73800000; // = (2. * pow(U_APS, 2) + 6. * pow(V_APS, 2));
  b3 = 0.66854900; // = 4. * V_APS * (pow(U_APS, 2) + pow(V_APS, 2));
  b4 = 0.10150800; // = pow(pow(U_APS, 2) + pow(V_APS, 2), 2);
  return pade_d2 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

// \partial_n \beta_x = \partial_n [b_x + \zeta_x + \eta_x^2 / c_x]
FDECORATOR double
beta_parameter_d0 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 1.55368;                        // +/- 0.005789     (0.3726%)
  b1 = 7.96800;                        // +/- 0.6046       (7.588%)
  a2 = 9.47166;                        // +/- 0.8251       (8.712%)
  b2 = 6.06735;                        // +/- 0.4984       (8.214%)
  a3 = 1.91935;                        // +/- 0.1414       (7.369%)
  b3 = 5.60946;                        // +/- 0.4538       (8.09%)
  c  = 1.16711;                        // +/- 0.09957      (8.532%)
  a4 = c; b4 = c;
  return BF_APS_UFG * pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
beta_parameter_d1 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 1.55368;                        // +/- 0.005789     (0.3726%)
  b1 = 7.96800;                        // +/- 0.6046       (7.588%)
  a2 = 9.47166;                        // +/- 0.8251       (8.712%)
  b2 = 6.06735;                        // +/- 0.4984       (8.214%)
  a3 = 1.91935;                        // +/- 0.1414       (7.369%)
  b3 = 5.60946;                        // +/- 0.4538       (8.09%)
  c  = 1.16711;                        // +/- 0.09957      (8.532%)
  a4 = c; b4 = c;
  return BF_APS_UFG * pade_d1 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
beta_parameter_d2 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 1.55368;                        // +/- 0.005789     (0.3726%)
  b1 = 7.96800;                        // +/- 0.6046       (7.588%)
  a2 = 9.47166;                        // +/- 0.8251       (8.712%)
  b2 = 6.06735;                        // +/- 0.4984       (8.214%)
  a3 = 1.91935;                        // +/- 0.1414       (7.369%)
  b3 = 5.60946;                        // +/- 0.4538       (8.09%)
  c  = 1.16711;                        // +/- 0.09957      (8.532%)
  a4 = c; b4 = c;
  return BF_APS_UFG * pade_d2 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

// \partial_n \gamma_x^{-1} = \partial_n [6/(3\pi^2)^{2/3} / c_x]
FDECORATOR double
inverse_gamma_parameter_d0 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 0.2680370;                      // +/- 0.0001483    (0.05535%)
  b1 = 2.5914700;                      // +/- 0.02606      (1.005%)
  a2 = 0.7472010;                      // +/- 0.008348     (1.117%)
  b2 = 2.7034900;                      // +/- 0.02437      (0.9013%)
  a3 = 0.3300800;                      // +/- 0.00252      (0.7635%)
  b3 = 0.5475050;                      // +/- 0.004573     (0.8353%)
  c  = 0.0892438;                      // +/- 0.0008447    (0.9465%)
  a4 = c; b4 = c;
  return (1. / CF_APS_UFG) / pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
inverse_gamma_parameter_d1 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 0.2680370;                      // +/- 0.0001483    (0.05535%)
  b1 = 2.5914700;                      // +/- 0.02606      (1.005%)
  a2 = 0.7472010;                      // +/- 0.008348     (1.117%)
  b2 = 2.7034900;                      // +/- 0.02437      (0.9013%)
  a3 = 0.3300800;                      // +/- 0.00252      (0.7635%)
  b3 = 0.5475050;                      // +/- 0.004573     (0.8353%)
  c  = 0.0892438;                      // +/- 0.0008447    (0.9465%)
  a4 = c; b4 = c;
  return (1. / CF_APS_UFG)
        * pade_d1 (_x, a1, a2, a3, a4, b1, b2, b3, b4)
        / pow(pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4), 2);
}

FDECORATOR double
inverse_gamma_parameter_d2 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 0.2680370;                      // +/- 0.0001483    (0.05535%)
  b1 = 2.5914700;                      // +/- 0.02606      (1.005%)
  a2 = 0.7472010;                      // +/- 0.008348     (1.117%)
  b2 = 2.7034900;                      // +/- 0.02437      (0.9013%)
  a3 = 0.3300800;                      // +/- 0.00252      (0.7635%)
  b3 = 0.5475050;                      // +/- 0.004573     (0.8353%)
  c  = 0.0892438;                      // +/- 0.0008447    (0.9465%)
  a4 = c; b4 = c;
  return (1. / CF_APS_UFG)
        * (pade_d2 (_x, a1, a2, a3, a4, b1, b2, b3, b4)
          * pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4)
          - pow(pade_d1 (_x, a1, a2, a3, a4, b1, b2, b3, b4), 2))
        / pow(pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4), 4);
}

// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// functional parameter A, B, and C                      [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// A_x = \sum_n (-x)^n (3!/(3+n)!) \partial_n \alpha_x
FDECORATOR double
a_functional_d0 (double _x)
{
  double d_ = 3.;
  double r_ = 0.;
  r_ += alpha_parameter_d0 (_x);
  r_ += (-_x) * alpha_parameter_d1 (_x) * tgamma (d_ + 1.) / tgamma (d_ + 1. + 1.);
  r_ += pow(-_x, 2) * alpha_parameter_d2 (_x) * tgamma (d_ + 1.) / tgamma (d_ + 1. + 2.);
  return r_;
}

// B_x = \sum_n (-x)^n (5!/(5+n)!) \partial_n \beta_x
FDECORATOR double
b_functional_d0 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 1.283420;                       // +/- 0.0029       (0.2259%)
  b1 = 4.810160;                       // +/- 0.2614       (5.435%)
  a2 = 4.443550;                       // +/- 0.2959       (6.659%)
  b2 = 2.976660;                       // +/- 0.1833       (6.159%)
  a3 = 0.773896;                       // +/- 0.04053      (5.238%)
  b3 = 2.177010;                       // +/- 0.1298       (5.964%)
  c  = 0.365191;                       // +/- 0.0235       (6.434%)
  a4 = c; b4 = c;
  return BF_APS_UFG * pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

// C_x^{-1} = \sum_n (-x)^n (1!/(1+n)!) \partial_n [1 / \gamma_x]
FDECORATOR double
c_functional_d0 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 0.14594800;                     // +/- 7.88e-05     (0.05399%)
  b1 = 1.08153000;                     // +/- 0.01173      (1.085%)
  a2 = 0.17443100;                     // +/- 0.002099     (1.204%)
  b2 = 0.68441200;                     // +/- 0.005664     (0.8276%)
  a3 = 0.04190220;                     // +/- 0.0002529    (0.6036%)
  b3 = 0.08023010;                     // +/- 0.0005894    (0.7347%)
  c  = 0.00855043;                     // +/- 7.706e-05    (0.9012%)
  a4 = c; b4 = c;
  return CF_APS_UFG * pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################
