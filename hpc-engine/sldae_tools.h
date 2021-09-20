// ###########################################################################
// ---------------------------------------------------------------------------
// Truncation parameters of series expansion
//      HFB_ORDER(MAX) = 8 (higher order can be added systematically)
//      FUNCTIONAL_ORDER(MAX) = maximal derivatives provided
//                              in the definition of physical quantites
// ---------------------------------------------------------------------------
#define HFB_ORDER 8
#define FUNCTIONAL_ORDER 2
// ---------------------------------------------------------------------------
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// defined mathematical and numerical constants          [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]
#define LN2 log(2.)
#define LN2_2 pow(LN2, 2)
#define LN2_3 pow(LN2, 3)
#define LN2_4 pow(LN2, 4)
#define LN2_5 pow(LN2, 5)
#define M_PI_SQ pow(M_PI, 2)

#define DDCC_EPSILON 1.0e-32
// minimal _x to avoide divergences
// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################



// ---------------------------------------------------------------------------


// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// usefull mathematical functions                        [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]
double rising_factorial (double x, int in);
double falling_factorial (double x, int in);
double binomial_coefficient (int in, int ik);
double exponential_integral (double x);
double pexp_bell_polynomial (int in, int ik, double * g);
double xm (int dn, double x, int * ip);
// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// general derivatives using Leibniz rules               [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]
double power_lrule (int dn, double _x, int i, int dk, int * id, double (*f) (int, double, int *));
double product_lrule (int dn, double _x, int dkf, int dkg, int * idf, int * idg, double (*f) (int, double, int *), double (*g) (int, double, int *));
double inverse_lrule (int dn, double _x, int dk, int * id, double (*f) (int, double, int *));
double composed_lrule (int dn, double _x, int dkf, int dkg, int * idf, int * idg, double (*f) (int, double, int *), double (*g) (int, double, int *));
// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ---------------------------------------------------------------------------
// functions Bn and Cn for the expansion of the SLDA parameters b and c
// (higher order can be added systematically by user)
// ---------------------------------------------------------------------------
double b_expansion (int dp, double y_x, int * idx);
double c_expansion (int dp, double y_x, int * idx);
// ---------------------------------------------------------------------------
// ###########################################################################
