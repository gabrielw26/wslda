/**
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 *
 * This file provides functions for the memory management.
 * In the long term, it should store all functions that manage memory for the td code.
 *
 * @author Gabriel Wlazlowski
 * @date 08.01.2025
 * */

#ifndef __TDWSLDA_MEMORY_MANAGEMENT__
#define __TDWSLDA_MEMORY_MANAGEMENT__

#ifdef __externc
#define __decorator extern "C"
#else
#define __decorator
#endif

#include <stddef.h>

// declarations
__decorator size_t mm_get_size_of_kernels_workspace(size_t n);
__decorator void * mm_get_pointer_to_kernels_workspace(size_t n);

__decorator size_t mm_get_size_of_derivatives_workspace(size_t n);
__decorator void * mm_get_pointer_to_derivatives_workspace(size_t n);

__decorator size_t mm_get_size_of_total_workspace(size_t n);
__decorator void * mm_get_pointer_to_total_workspace(size_t n);

__decorator void * mm_get_pointer_to_fftONE_workspace(size_t n);
#endif
