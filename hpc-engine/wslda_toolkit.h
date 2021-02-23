/**
 * W-SLDA Toolkit
 * 
 * This file includes all header files, needed for auxliary tools.
 * */  

#ifndef __WSLDA_TOOLKIT__
#define __WSLDA_TOOLKIT__

#ifndef TDWSLDA
#define WSLDA
#endif

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include "pca_settings.h"
#include "wslda_errors.h"
#include "pca_macro.h"
#ifdef WSLDA_NO_MAIN_FUNCTION
#include "pca_utils.h"
#else
#include "pca_utils.c"
#endif
#include "wslda_potdens.h"
#include "wslda_wavevectors.h"
#include "pca_io.h"
#include "s2dpca_edf.h"
#include "pca_uniform.h"
#include "wslda_interpolation.h"
#include "wslda_resize.h"


#endif
