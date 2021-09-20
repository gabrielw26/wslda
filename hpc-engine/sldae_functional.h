// ###########################################################################
// ---------------------------------------------------------------------------
// Value of dimensionless physical quantity of unitary gas at zero temperature
// [default: APS[x,y,z] parametrization of the functional]
//      GSE_UFG = ground state energy
//      PGE_UFG = pairing gap function
//      IEM_UFG = inverse effective mass
//
// GSE_UFG = 1-16/(3*PI)*atan((5/ 24) / ((6/(35*PI))*(11-2*ln(2))))
// from APS_X functional, in aggrement with MBPT(2)
//
// PGE_UFG and IEM_UFG are adjusted on Gorkov Green's function caluclation of
// [Phys. Rev. A80, 063612 (2009)] defining the APS[x,y,z] parametrization.
// ---------------------------------------------------------------------------
#define GSE_UFG 0.3582341
#define PGF_UFG 0.4600000
#define IEM_UFG 0.8403361
// ---------------------------------------------------------------------------
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ---------------------------------------------------------------------------
// constants used in APS[x,y,z] parametrization
// ---------------------------------------------------------------------------

// APS[x,y,z] constants
// ====================
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
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ---------------------------------------------------------------------------
// functions of _x used in APS[x,y,z] parametrization
// ---------------------------------------------------------------------------
double s_aps (int dn, double _x, int * id);
double s_aps_2d1 (int dn, double _x, int * id);
double pairing_bcs (int dn, double _x, int * id);
double pairing_fit (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
// ###########################################################################



// ---------------------------------------------------------------------------



// ###########################################################################
// ---------------------------------------------------------------------------
// ground state energy, chemical potential, inverse effective mass,
// and pairing gap as a function of _x = |akF| (can be re-defiend by user)
// ---------------------------------------------------------------------------
double ground_state_energy (int dn, double _x, int * id);
double chemical_potential (int dn, double _x, int * id);
double inverse_effective_mass (int dn, double _x, int * id);
double pairing_gap (int dn, double _x, int * id);
// ---------------------------------------------------------------------------
// ###########################################################################




// ---------------------------------------------------------------------------




// ###########################################################################
// ------------------------------------------------------[! DO NOT MODIFY >>>]
// usefull mathematical tools for calculation            [! DO NOT MODIFY >>>]
// of HFB parameters a, b, and 1/c                       [! DO NOT MODIFY >>>]
// ------------------------------------------------------[! DO NOT MODIFY >>>]
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
// ------------------------------------------------------[<<< DO NOT MODIFY !]
// ###########################################################################
