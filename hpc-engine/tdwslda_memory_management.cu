/**
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 *
 * This file provides functions for the memory management.
 * In the long term, it should store all functions that manage memory for the td code.
 *
 * NOTE: This file is closely related to *derivative.cu files!
 *
 * @author Gabriel Wlazlowski
 * @date 08.01.2025
 * */

#include "pca_settings.h"
#include "pca_macro.h"
// #include "tdwslda_memory_management.h"

// global pointer to workspace
extern "C" void *pca_cufft_work_area;


// --------------------------- Intrinsic functions ------------------------------

// --------------------------- End of Intrinsic functions ------------------------------

/**
 * The function returns the workspace size that can be used during the execution of kernels from *kernels.cu files.
 * @return size in Bytes.
 * */
extern "C" size_t mm_get_size_of_kernels_workspace(size_t n)
{
    return sizeof(double)*n*PCA_WORKSPACE_SHIFT;
}

/**
 * The function returns the workspace size that can be used during the execution of functions from *derivative.cu files.
 * @return size in Bytes.
 * */
extern "C" size_t mm_get_size_of_derivatives_workspace(size_t n)
{
    return sizeof(double)*n*PCA_DERIVATIVE_SHIFT;
}

/**
 * The function returns the total workspace size that can be used during the execution of functions from *derivative.cu and *kernels.cu files.
 * @return size in Bytes.
 * */
extern "C" size_t mm_get_size_of_total_workspace(size_t n)
{
    return mm_get_size_of_derivatives_workspace(n)+mm_get_size_of_kernels_workspace(n);
}

/**
 * The function returns the pointer to workspace that can be used during the execution of kernels from *kernels.cu files.
 * @return pointer.
 * */
extern "C" void * mm_get_pointer_to_kernels_workspace(size_t n)
{
    return pca_cufft_work_area;
}


/**
 * The function returns the pointer workspace that can be used during the execution of functions from *derivative.cu files.
 * @return pointer.
 * */
extern "C" void * mm_get_pointer_to_derivatives_workspace(size_t n)
{
    double *ptr = (double *)pca_cufft_work_area;
    ptr+=PCA_WORKSPACE_SHIFT*n;
    return (void *)ptr;
}

/**
 * The function returns the pointer to total workspace that can be used during the execution of kernels from *derivative.cu and *kernels.cu files.
 * @return pointer.
 * */
extern "C" void * mm_get_pointer_to_total_workspace(size_t n)
{
    return pca_cufft_work_area;
}

/**
 * The function returns the pointer to workspace that can be used by type ONE fft plans.
 * @return pointer.
 * */
extern "C" void * mm_get_pointer_to_fftONE_workspace(size_t n)
{
    double *ptr = (double *)pca_cufft_work_area;
    ptr+=PCA_WORKSPACE_SHIFT*n;
    ptr+=PCA_DERIVATIVE_SHIFT*n;
    return (void *)ptr;
}
