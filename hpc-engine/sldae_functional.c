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
  a1 = 1.58442 ;                     // +/- 0.01443      (0.9109%)
  b1 = 13.9448 ;                     // +/- 1.86         (13.34%)
  a2 = 17.2973 ;                     // +/- 2.525        (14.6%)
  b2 = 10.1914 ;                     // +/- 1.431        (14.04%)
  a3 = 0.820671;                     // +/- 0.07035      (8.573%)
  b3 = 11.9429 ;                     // +/- 1.636        (13.7%)
  c  = 3.48966 ;                     // +/- 0.4957       (14.21%)
  a4 = c; b4 = c;
  return BF_APS_UFG * pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
beta_parameter_d1 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 1.58442 ;                     // +/- 0.01443      (0.9109%)
  b1 = 13.9448 ;                     // +/- 1.86         (13.34%)
  a2 = 17.2973 ;                     // +/- 2.525        (14.6%)
  b2 = 10.1914 ;                     // +/- 1.431        (14.04%)
  a3 = 0.820671;                     // +/- 0.07035      (8.573%)
  b3 = 11.9429 ;                     // +/- 1.636        (13.7%)
  c  = 3.48966 ;                     // +/- 0.4957       (14.21%)
  a4 = c; b4 = c;
  return BF_APS_UFG * pade_d1 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
beta_parameter_d2 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 1.58442 ;                     // +/- 0.01443      (0.9109%)
  b1 = 13.9448 ;                     // +/- 1.86         (13.34%)
  a2 = 17.2973 ;                     // +/- 2.525        (14.6%)
  b2 = 10.1914 ;                     // +/- 1.431        (14.04%)
  a3 = 0.820671;                     // +/- 0.07035      (8.573%)
  b3 = 11.9429 ;                     // +/- 1.636        (13.7%)
  c  = 3.48966 ;                     // +/- 0.4957       (14.21%)
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
  a1 = 6.15853  ;                    // +/- 0.01528      (0.2482%)
  b1 = 1.48431  ;                    // +/- 0.01418      (0.9551%)
  a2 = 0.156722 ;                    // +/- 0.04688      (29.91%)
  b2 = 2.07468  ;                    // +/- 0.03981      (1.919%)
  a3 = 0.76668  ;                    // +/- 0.02634      (3.436%)
  b3 = 0.737836 ;                    // +/- 0.02402      (3.255%)
  c  = 0.0823379;                    // +/- 0.006872     (8.346%)
  a4 = c; b4 = c;
  return (1. / CF_APS_UFG) * pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
inverse_gamma_parameter_d1 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 6.15853  ;                    // +/- 0.01528      (0.2482%)
  b1 = 1.48431  ;                    // +/- 0.01418      (0.9551%)
  a2 = 0.156722 ;                    // +/- 0.04688      (29.91%)
  b2 = 2.07468  ;                    // +/- 0.03981      (1.919%)
  a3 = 0.76668  ;                    // +/- 0.02634      (3.436%)
  b3 = 0.737836 ;                    // +/- 0.02402      (3.255%)
  c  = 0.0823379;                    // +/- 0.006872     (8.346%)
  a4 = c; b4 = c;
  return (1. / CF_APS_UFG) * pade_d1 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
inverse_gamma_parameter_d2 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 6.15853  ;                    // +/- 0.01528      (0.2482%)
  b1 = 1.48431  ;                    // +/- 0.01418      (0.9551%)
  a2 = 0.156722 ;                    // +/- 0.04688      (29.91%)
  b2 = 2.07468  ;                    // +/- 0.03981      (1.919%)
  a3 = 0.76668  ;                    // +/- 0.02634      (3.436%)
  b3 = 0.737836 ;                    // +/- 0.02402      (3.255%)
  c  = 0.0823379;                    // +/- 0.006872     (8.346%)
  a4 = c; b4 = c;
  return (1. / CF_APS_UFG) * pade_d2 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// functional parameter A, B, and C                      [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]

// A_x = alpha_x
FDECORATOR double
a_functional_d0 (double _x)
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

// B_x = \sum_n (-x)^n (5!/(5+n)!) \partial_n \beta_x
FDECORATOR double
b_functional_d0 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 1.29871 ;                     // +/- 0.005204     (0.4007%)
  b1 = 6.47079 ;                     // +/- 0.462        (7.14%)
  a2 = 6.02212 ;                     // +/- 0.5124       (8.508%)
  b2 = 3.73389 ;                     // +/- 0.2985       (7.993%)
  a3 = 0.376541;                     // +/- 0.01686      (4.477%)
  b3 = 3.47301 ;                     // +/- 0.2627       (7.565%)
  c  = 0.810972;                     // +/- 0.06537      (8.06%)
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
  a1 = 0.267952 ;                    // +/- 0.0001541    (0.05752%)
  b1 = 2.61695  ;                    // +/- 0.02719      (1.039%)
  a2 = 0.754923 ;                    // +/- 0.008709     (1.154%)
  b2 = 2.71855  ;                    // +/- 0.0254       (0.9345%)
  a3 = 0.330699 ;                    // +/- 0.002623     (0.7931%)
  b3 = 0.550587 ;                    // +/- 0.004768     (0.866%)
  c  = 0.0900897;                    // +/- 0.0008812    (0.9781%)
  a4 = c; b4 = c;
  return CF_APS_UFG * pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################
