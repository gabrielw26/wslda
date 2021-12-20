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
  a1 = 1.61463;                      // +/- 0.03887      (2.407%)
  b1 = 31.2584;                      // +/- 9.376        (29.99%)
  a2 = 42.5454;                      // +/- 13.18        (30.97%)
  b2 = 23.7481;                      // +/- 7.214        (30.38%)
  a3 = 3.85825;                      // +/- 1.08         (27.99%)
  b3 = 25.1379;                      // +/- 7.611        (30.28%)
  c  = 6.2833 ;                      // +/- 1.934        (30.79%)
  a4 = c; b4 = c;
  return BF_APS_UFG * pade_d0 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
beta_parameter_d1 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 1.61463;                      // +/- 0.03887      (2.407%)
  b1 = 31.2584;                      // +/- 9.376        (29.99%)
  a2 = 42.5454;                      // +/- 13.18        (30.97%)
  b2 = 23.7481;                      // +/- 7.214        (30.38%)
  a3 = 3.85825;                      // +/- 1.08         (27.99%)
  b3 = 25.1379;                      // +/- 7.611        (30.28%)
  c  = 6.2833 ;                      // +/- 1.934        (30.79%)
  a4 = c; b4 = c;
  return BF_APS_UFG * pade_d1 (_x, a1, a2, a3, a4, b1, b2, b3, b4);
}

FDECORATOR double
beta_parameter_d2 (double _x)
{
  double a1, a2, a3, a4, b1, b2, b3, b4, c;
  //Final set of parameters            Asymptotic Standard Error
  //=======================            ==========================
  a1 = 1.61463;                      // +/- 0.03887      (2.407%)
  b1 = 31.2584;                      // +/- 9.376        (29.99%)
  a2 = 42.5454;                      // +/- 13.18        (30.97%)
  b2 = 23.7481;                      // +/- 7.214        (30.38%)
  a3 = 3.85825;                      // +/- 1.08         (27.99%)
  b3 = 25.1379;                      // +/- 7.611        (30.28%)
  c  = 6.2833 ;                      // +/- 1.934        (30.79%)
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
  a1 = 1.30703;                      // +/- 0.009127     (0.6983%)
  b1 = 10.1728;                      // +/- 1.253        (12.32%)
  a2 = 10.6772;                      // +/- 1.437        (13.46%)
  b2 = 6.15932;                      // +/- 0.7924       (12.87%)
  a3 = 1.02277;                      // +/- 0.1118       (10.93%)
  b3 = 5.37162;                      // +/- 0.682        (12.7%)
  c  = 1.07145;                      // +/- 0.1414       (13.2%)
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
