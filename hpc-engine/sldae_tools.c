#include "sldae_tools.h"
// ---------------------------------------------------------------------------



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
