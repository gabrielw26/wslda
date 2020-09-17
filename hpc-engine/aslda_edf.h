/** 
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 * 
 * This header contains list of functions for ASLDA functional
 * 
 * The functions are implemented in:
 * td --> CUDA: pca_edf.h
 * st --> C   : s2dpca_edf.h
 * 
 * @author Gabriel Wlazlowski
 * @date 04.09.2020
 * */ 

#ifndef _ASLDA_EDF_FUNCTIONS_
#define _ASLDA_EDF_FUNCTIONS_

// EDF functions
double polarization_h(double n_a, double n_b);
double der_polarization__der_na_h(double n_a, double n_b);
double der_polarization__der_nb_h(double n_a, double n_b);
double alpha_h(double p);
double alpha_a_h(double p);
double alpha_b_h(double p);
double alpha_plus_h(double p);
double alpha_minus_h(double p);
double der_alpha_plus__der_na_h(double n_a, double n_b);
double der_alpha_plus__der_nb_h(double n_a, double n_b);
double der_alpha_minus__der_na_h(double n_a, double n_b);
double der_alpha_minus__der_nb_h(double n_a, double n_b);
double funG_h(double p);
double funD_h(double n_a, double n_b);
double der_funD__der_na_h(double n_a, double n_b);
double der_funD__der_nb_h(double n_a, double n_b);
double tildeC_h(double n_a, double n_b);
double der_tildeC__der_na_h(double n_a, double n_b);
double der_tildeC__der_nb_h(double n_a, double n_b);

double polarization(double n_a, double n_b);
double der_polarization__der_na(double n_a, double n_b);
double der_polarization__der_nb(double n_a, double n_b);
double alpha(double p);
double alpha_a(double p);
double alpha_b(double p);
double alpha_plus(double p);
double alpha_minus(double p);
double der_alpha_plus__der_na(double n_a, double n_b);
double der_alpha_plus__der_nb(double n_a, double n_b);
double der_alpha_minus__der_na(double n_a, double n_b);
double der_alpha_minus__der_nb(double n_a, double n_b);
double funG(double p);
double funD(double n_a, double n_b);
double der_funD__der_na(double n_a, double n_b);
double der_funD__der_nb(double n_a, double n_b);
double tildeC(double n_a, double n_b);
double der_tildeC__der_na(double n_a, double n_b);
double der_tildeC__der_nb(double n_a, double n_b);
double p_regularization(double n);
double der_p_regularization(double n);

#endif
