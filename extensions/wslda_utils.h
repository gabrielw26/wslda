/**
 * W-SLDA Toolkit
 *
 * This file contains collections of varius usfull functions
 * that faciliate specyfing the physical problem.
 *
 * */
#ifndef EXTENSION_WSLDA_UTILS
#define EXTENSION_WSLDA_UTILS

#ifdef TDWSLDA
#define CDECORATOR __device__ __host__
#else
#define CDECORATOR
#endif

/**
 * Switch function.
 * It rises smoothly from zero to one in range interval [0,1]
 * @param x
 * @param alpha controls the smoothness of the step function
 * */
CDECORATOR
double smooth_from0to1(double x, double alpha)
{
    if(x<=0.0) return 0.0;
    if(x>=1.0) return 1.0;
    return 0.5*( 1.0+tanh( alpha*tan( M_PI_2*( 2.0*x-1.0 ) ) ) );
}

/**
 * Function of form:
 * f(x) = lv for x < xl;
 *      = smooth switch from lv to rv for x in [xl,xr]
 *      = rv for x > xr;
 * @param x
 * @param xl x coordinate for start of switching
 * @param xr x coordinate for stop of switching
 * @param lv left value, for x<xl
 * @param rv right value, for x>xr
 * @param alpha controls the smoothness of the step function
 * */
CDECORATOR
double smooth_fromatob(double x, double xl, double xr, double lv, double rv)
{
    return lv + (rv-lv)*smooth_from0to1((x-xl)/(xr-xl),1.0);
}

/**
 * Smooth step function:
 *   [-infty,    x11) : 0
 *   [   x11,    x12) : increases smoothly from 0 to 1
 *   [   x12,    x21) : 1
 *   [   x21,    x22) : decreases smoothly from 1 to 0
 *   [   x22, +infty] : 0
 * @param x
 * @param alpha controls the smoothness of the step function
 * */
CDECORATOR
double smooth_step(double x, double x11, double x12, double x21, double x22, double alpha)
{
    if     (x<=x11) return 0.0;
    else if(x<=x12) return smooth_from0to1((x-x11)/(x12-x11),alpha);
    else if(x<=x21) return 1.0;
    else if(x<=x22) return 1.0-smooth_from0to1((x-x21)/(x22-x21),alpha);
    else            return 0.0;
}

/**
 * See Fig.4 of https://arxiv.org/pdf/1711.05803.pdf
 * @param x coordinate
 * @param omega frequency
 * @param x1 starting coordinate for smoothing, for x<=x1 value of the potential is omega^2*x^2/2
 * @param x2 for x>=x2 value of the potential is omega^2*x2^2/2
 * @param alpha controls smoothness of the step function
 * */
CDECORATOR
double harmonic_oscillator_smooth_edges(double x, double omega, double x1, double x2, double alpha)
{
    double xx=x*x, xx1=x1*x1, xx2=x2*x2, o2=0.5*omega*omega;

    double V=o2*xx;
    if(xx<=xx1) return V;

    double V2=o2*xx2;
    if(xx>=xx2) return V2;

    double s = smooth_from0to1((xx-xx1)/(xx2-xx1), alpha);
    return V*(1.0-s) + V2*s;
}

#undef CDECORATOR
#endif
