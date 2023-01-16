/**
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 *
 * @author Gabriel Wlazlowski
 *
 * The file contains basic routines for integration of ODE
 * */
#ifndef __TDWSLDA_ODE_INTEGRATOR__
#define __TDWSLDA_ODE_INTEGRATOR__

#if CODEDIM==1
#define ODE_NUMBER_ELEMENT NX
#elif CODEDIM==2
#define ODE_NUMBER_ELEMENT NXY
#elif CODEDIM==3
#define ODE_NUMBER_ELEMENT NXYZ
#else
#error CODEDIM NOT SPECIFIED
#endif

// =======================================================================================
// ====================================== abm_step1 ======================================
// =======================================================================================
__global__ void kernel_abm34_step1(size_t n, Complex *ykm1,
                         Complex *fkm1, Complex *fkm2, Complex *fkm3)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex _ykm1, _fkm1, _fkm2, _fkm3;
    Complex dti;

    if(ixyz<n*2*ODE_NUMBER_ELEMENT)
    {
        // read data
        _ykm1=ykm1[ixyz];
        _fkm1=fkm1[ixyz];
        _fkm2=fkm2[ixyz];
        _fkm3=fkm3[ixyz];

        dti = Complex(0.0, -1.0*dc_dt);

        // save data
        ykm1[ixyz] = _ykm1 + dti*(_fkm1*(23./12.) - _fkm2*(16./12.) + _fkm3*(5./12.));
        fkm3[ixyz] = _ykm1 + dti*(_fkm1*(19./24.) - _fkm2*( 5./24.) + _fkm3*(1./24.));

    }
}

__global__ void kernel_abm45_step1(size_t n, Complex *ykm1,
                         Complex *fkm1, Complex *fkm2, Complex *fkm3, Complex *fkm4)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex _ykm1, _fkm1, _fkm2, _fkm3, _fkm4;
    Complex dti;

    if(ixyz<n*2*ODE_NUMBER_ELEMENT)
    {
        // read data
        _ykm1=ykm1[ixyz];
        _fkm1=fkm1[ixyz];
        _fkm2=fkm2[ixyz];
        _fkm3=fkm3[ixyz];
        _fkm4=fkm4[ixyz];

        dti = Complex(0.0, -1.0*dc_dt);

        // save data
        ykm1[ixyz] = _ykm1 + dti*(_fkm1*(55./24.  ) - _fkm2*(59./24.  ) + _fkm3*(37./24.  ) - _fkm4*(9./24.  ));
        fkm4[ixyz] = _ykm1 + dti*(_fkm1*(646./720.) - _fkm2*(264./720.) + _fkm3*(106./720.) - _fkm4*(19./720.));
    }
}

__global__ void kernel_abm55_step1(size_t n, Complex *ykm1,
                         Complex *fkm1, Complex *fkm2, Complex *fkm3, Complex *fkm4, Complex *fkm5)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex _ykm1, _fkm1, _fkm2, _fkm3, _fkm4, _fkm5;
    Complex dti;

    if(ixyz<n*2*ODE_NUMBER_ELEMENT)
    {
        // read data
        _ykm1=ykm1[ixyz];
        _fkm1=fkm1[ixyz];
        _fkm2=fkm2[ixyz];
        _fkm3=fkm3[ixyz];
        _fkm4=fkm4[ixyz];
        _fkm5=fkm5[ixyz];

        dti = Complex(0.0, -1.0*dc_dt);

        // save data
        ykm1[ixyz] = _ykm1 + dti*(_fkm1*(1901./720.) - _fkm2*(2774./720.) + _fkm3*(2616./720.) - _fkm4*(1274./720.) + _fkm5*( 251./720.));
        fkm5[ixyz] = _ykm1 + dti*(_fkm1*( 646./720.) - _fkm2*( 264./720.) + _fkm3*( 106./720.) - _fkm4*(  19./720.));
    }
}

/**
 * Functions perform step 1 from intgration.pdf
 * @param n  number of wave-functions (u,v pairs) to process
 * @param ykm1-ykm5 historical data, array y_{k-1} of size 2*n*ODE_NUMBER_ELEMENT
 * @param abm_scheme integration scheme
 * @return 0 - OK, otherwise ERROR
 * */
extern "C" int abm_step1(int n, cufftDoubleComplex *ykm1,
                         cufftDoubleComplex *fkm1,
                         cufftDoubleComplex *fkm2,
                         cufftDoubleComplex *fkm3,
                         cufftDoubleComplex *fkm4,
                         cufftDoubleComplex *fkm5,
                         int abm_scheme,
                         int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)2*ODE_NUMBER_ELEMENT*n/nthreads);
    if(abm_scheme==AB3AM4)
    {
        kernel_abm34_step1<<<nblocks, nthreads>>>(n, (Complex *)ykm1, (Complex *)fkm1, (Complex *)fkm2, (Complex *)fkm3);
    }
    else if(abm_scheme==AB4AM5)
    {
        kernel_abm45_step1<<<nblocks, nthreads>>>(n, (Complex *)ykm1, (Complex *)fkm1, (Complex *)fkm2, (Complex *)fkm3, (Complex *)fkm4);
    }
    else if(abm_scheme==AB5AM5)
    {
        kernel_abm55_step1<<<nblocks, nthreads>>>(n, (Complex *)ykm1, (Complex *)fkm1, (Complex *)fkm2, (Complex *)fkm3, (Complex *)fkm4, (Complex *)fkm5);
    }
    else
    {
        return -999999; // should never happen
    }

    return 0;
}


__global__ void kernel_abm34_step4(size_t n, Complex *ykm1_in, Complex *ykm1_out, Complex *fkm3)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex _ykm1, _fkm3;
    Complex dti;

    if(ixyz<n*2*ODE_NUMBER_ELEMENT)
    {
        // read data
        _ykm1=ykm1_in[ixyz];
        _fkm3=fkm3[ixyz];

        dti = Complex(0.0, -1.0*dc_dt);

        // save data
        ykm1_out[ixyz] = dti*_ykm1*(9./24.) + _fkm3;
    }
}

__global__ void kernel_abm45_step4(size_t n, Complex *ykm1_in, Complex *ykm1_out, Complex *fkm4)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex _ykm1, _fkm4;
    Complex dti;

    if(ixyz<n*2*ODE_NUMBER_ELEMENT)
    {
        // read data
        _ykm1=ykm1_in[ixyz];
        _fkm4=fkm4[ixyz];

        dti = Complex(0.0, -1.0*dc_dt);

        // save data
        ykm1_out[ixyz] = dti*_ykm1*(251./720.) + _fkm4;
    }
}

__global__ void kernel_abm55_step4(size_t n, Complex *ykm1_in, Complex *ykm1_out, Complex *fkm5)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex _ykm1, _fkm5;
    Complex dti;

    if(ixyz<n*2*ODE_NUMBER_ELEMENT)
    {
        // read data
        _ykm1=ykm1_in[ixyz];
        _fkm5=fkm5[ixyz];

        dti = Complex(0.0, -1.0*dc_dt);

        // save data
        ykm1_out[ixyz] = dti*_ykm1*(251./720.) + _fkm5;
    }
}

extern "C" int abm_step4(int n, cufftDoubleComplex *ykm1_in, cufftDoubleComplex *ykm1_out,
                         cufftDoubleComplex *fkm_last,
                         int abm_scheme,
                         int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)2*ODE_NUMBER_ELEMENT*n/nthreads);
    if(abm_scheme==AB3AM4)
    {
        kernel_abm34_step4<<<nblocks, nthreads>>>(n, (Complex *)ykm1_in, (Complex *)ykm1_out, (Complex *)fkm_last);
    }
    else if(abm_scheme==AB4AM5)
    {
        kernel_abm45_step4<<<nblocks, nthreads>>>(n, (Complex *)ykm1_in, (Complex *)ykm1_out, (Complex *)fkm_last);
    }
    else if(abm_scheme==AB5AM5)
    {
        kernel_abm55_step4<<<nblocks, nthreads>>>(n, (Complex *)ykm1_in, (Complex *)ykm1_out, (Complex *)fkm_last);
    }
    else
    {
        return -999999; // should never happen
    }

    return 0;
}

#undef ODE_NUMBER_ELEMENT
#endif
