/** 
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 * 
 * The file implements routines for computation of potentials from densities
 * These routines are utilized by ST codes only!
 *  
 * @author Gabriel Wlazlowski
 * @date 04.09.2020
 * */ 

#ifndef _WSLDA_FUNCTIONALS_
#define _WSLDA_FUNCTIONALS_
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <complex.h>

#include "pca_settings.h"
#include "pca_macro.h"
#include "wslda_potdens.h"

#include "aslda_edf.h"

/**
 * Function computes potentials for selected functional
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials potentials from PREVIOUS iteration as input, 
 *                     updated values as output (INPUT/OUTPUT) 
 * */
int compute_potentials(int it, wslda_density h_densities, wslda_potential h_potentials);

#endif

