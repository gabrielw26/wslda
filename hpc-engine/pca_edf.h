// Author: Gabriel Wlazlowski
// 
// This header implements energy density functional of ultra-cold atoms
// and functional derivaties

#ifndef __PCA_EDF__
#define __PCA_EDF__

// ===================================================================================
// ================================= FUNCTIONS =======================================
// ===================================================================================


/**
 * Regularization funtion
 * */
__device__ __host__ inline double p_regularization(double n)
{
    if(n<=P_NMIN) return 0.0;
    else if(n>=P_NMAX) return 1.0;
    else return 0.5 + 0.5*tanh(P_ALPHA * tan( (M_PI*(n-P_NMIN))/(P_NMAX-P_NMIN) - M_PI/2.0 ) );
}

/**
 * Regularization funtion - derivative
 * */
__device__ __host__ inline double der_p_regularization(double n)
{
    if(n<=P_NMIN) return 0.0;
    else if(n>=P_NMAX) return 0.0;
    else return 0.5*P_ALPHA*M_PI / pow(cosh(P_ALPHA * tan( (M_PI*(n-P_NMIN))/(P_NMAX-P_NMIN) - M_PI/2.0 ) )*cos( (M_PI*(n-P_NMIN))/(P_NMAX-P_NMIN) - M_PI/2.0 ),2) / (P_NMAX-P_NMIN);
}

/**
 * Local polarization
 * @return (n_a-n_b)/(n_a+n_b)
 * */
__device__ __host__ inline double polarization(double n_a, double n_b)
{
//     return (n_a-n_b)/(n_a+n_b+DENSEPSILON);
    
    double pr = p_regularization(n_a+n_b);
    if(pr==0.0) return 0.0;
    else return pr*(n_a-n_b)/(n_a+n_b);
}
extern "C" double polarization_h(double n_a, double n_b)
{
    return polarization(n_a, n_b);
}

/**
 * derivative dp/dn_a
 * */
__device__ __host__ inline double der_polarization__der_na(double n_a, double n_b)
{
//     return (2.0*n_b+DENSEPSILON)/( (n_a+n_b+DENSEPSILON)*(n_a+n_b+DENSEPSILON) );
    
    double np=n_a+n_b;
    double pr = p_regularization(np);
    double dpr = der_p_regularization(np);
    
    double r=0.0;
    if(pr!=0.0) r+=pr*(2.0*n_b)/( (np)*(np) );
    if(dpr!=0.0) r+=dpr*(n_a-n_b)/np;
    return r;
}
extern "C" double der_polarization__der_na_h(double n_a, double n_b)
{
    return  der_polarization__der_na( n_a,  n_b);
}

/**
 * derivative dp/dn_b
 * */
__device__ __host__ inline double der_polarization__der_nb(double n_a, double n_b)
{
//     return (-2.0*n_a-DENSEPSILON)/( (n_a+n_b+DENSEPSILON)*(n_a+n_b+DENSEPSILON) );
    
    double np=n_a+n_b;
    double pr = p_regularization(np);
    double dpr = der_p_regularization(np);
    
    double r=0.0;
    if(pr!=0.0) r+=pr*(-2.0*n_a)/( (np)*(np) );
    if(dpr!=0.0) r+=dpr*(n_a-n_b)/np;
    return r;
}
extern "C" double der_polarization__der_nb_h(double n_a, double n_b)
{
    return  der_polarization__der_nb( n_a,  n_b);
}

/*
 * Equation: (9.98) from BFM review
 * */
__device__ __host__ inline double alpha(double p)
{
    double p2=p*p;
    double p4=p2*p2;
    return A0 + A1*p*(1.0 - (2.0/3.0)*p2 + (1.0/5.0)*p4) + A2*p2*(1.0 - p2 + (1.0/3.0)*p4);
}
extern "C" double alpha_h(double p)
{
    return  alpha( p);
}

/*
 * Equation: (9.98) from BFM review
 * */
__device__ __host__ inline double alpha_a(double p)
{
    return alpha(p);
}
extern "C" double alpha_a_h(double p)
{
    return  alpha_a( p);
}

/*
 * Equation: (9.98) from BFM review
 * */
__device__ __host__ inline double alpha_b(double p)
{
    return alpha(-1.0*p);
}
extern "C" double alpha_b_h(double p)
{
    return  alpha_b( p);
}

__device__ __host__ inline double alpha_plus(double p)
{
    return 0.5*(alpha_a(p)+alpha_b(p));
}
extern "C" double alpha_plus_h(double p)
{
    return  alpha_plus( p);
}

__device__ __host__ inline double alpha_minus(double p)
{
    return 0.5*(alpha_a(p)-alpha_b(p));
}
extern "C" double alpha_minus_h(double p)
{
    return  alpha_minus( p);
}

/**
 * derivative dalpha_plus/dn_a
 * */
__device__ __host__ inline double der_alpha_plus__der_na(double n_a, double n_b)
{
    double p=polarization(n_a,n_b);
    double p2=p*p;
    double p4=p2*p2;
    
    return ( 2.0*A2*p*(1.0 - p2 + (1.0/3.0)*p4) - 2.0*A2*p*p2*(1.0 - (2.0/3.0)*p2) ) * der_polarization__der_na(n_a, n_b);
}
extern "C" double der_alpha_plus__der_na_h(double n_a, double n_b)
{
    return  der_alpha_plus__der_na( n_a,  n_b);
}

/**
 * derivative dalpha_plus/dn_b
 * */
__device__ __host__ inline double der_alpha_plus__der_nb(double n_a, double n_b)
{
    double p=polarization(n_a,n_b);
    double p2=p*p;
    double p4=p2*p2;
    
    return ( 2.0*A2*p*(1.0 - p2 + (1.0/3.0)*p4) - 2.0*A2*p*p2*(1.0 - (2.0/3.0)*p2) ) * der_polarization__der_nb(n_a, n_b);
}
extern "C" double der_alpha_plus__der_nb_h(double n_a, double n_b)
{
    return  der_alpha_plus__der_nb( n_a,  n_b);
}

/**
 * derivative dalpha_minus/dn_a
 * */
__device__ __host__ inline double der_alpha_minus__der_na(double n_a, double n_b)
{
    double p=polarization(n_a,n_b);
    double p2=p*p;
    double p4=p2*p2;
    
    return ( A1*(1.0 - (2.0/3.0)*p2 + (1.0/5.0)*p4) - 4.0*A1*p2*((1.0/3.0) - (1.0/5.0)*p2) ) * der_polarization__der_na(n_a, n_b);
}
extern "C" double der_alpha_minus__der_na_h(double n_a, double n_b)
{
    return  der_alpha_minus__der_na( n_a,  n_b);
}

/**
 * derivative dalpha_minus/dn_b
 * */
__device__ __host__ inline double der_alpha_minus__der_nb(double n_a, double n_b)
{
    double p=polarization(n_a,n_b);
    double p2=p*p;
    double p4=p2*p2;
    
    return ( A1*(1.0 - (2.0/3.0)*p2 + (1.0/5.0)*p4) - 4.0*A1*p2*((1.0/3.0) - (1.0/5.0)*p2) ) * der_polarization__der_nb(n_a, n_b);
}
extern "C" double der_alpha_minus__der_nb_h(double n_a, double n_b)
{
    return  der_alpha_minus__der_nb( n_a,  n_b);
}

/*
 * Equation: (9.91) from BFM review
 * */
__device__ __host__ inline double funG(double p)
{
    return G0 + G1*p*p;
}
extern "C" double funG_h(double p)
{
    return  funG( p);
}

/**
 * normal phase part
 * */
__device__ __host__ inline double funD(double n_a, double n_b)
{
    double p=polarization(n_a,n_b);
    return (1.0 / (20.0*M_PI*M_PI))*pow(6.0*M_PI*M_PI*(n_a+n_b),5.0/3.0) *
            (
              funG(p) - alpha(p)*pow((1.0+p)/2.0,5.0/3.0) - alpha(-1.0*p)*pow((1.0-p)/2.0,5.0/3.0) 
            );
}
extern "C" double funD_h(double n_a, double n_b)
{
    return  funD( n_a,  n_b);
}

/**
 * normal phase part
 * derivative: dD/dna
 * */
__device__ __host__ inline double der_funD__der_na(double n_a, double n_b)
{
    double p=polarization(n_a,n_b);
    double p2=p*p;
    double p4=p2*p2;
    
    return (1.0 / (2.0))*pow(6.0*M_PI*M_PI*(n_a+n_b),2.0/3.0) *
            (
              funG(p) - alpha(p)*pow((1.0+p)/2.0,5.0/3.0) - alpha(-1.0*p)*pow((1.0-p)/2.0,5.0/3.0) 
            )  +
           (1.0 / (20.0*M_PI*M_PI))*pow(6.0*M_PI*M_PI*(n_a+n_b),5.0/3.0) *
            (
               2.0*G1*p 
               - ( // dalpha(p)/dna * ((1+p)/2)^(5/3)
                   A1*(1.0 - (2.0/3.0)*p2 + (1.0/5.0)*p4) - 4.0*A1*p2*((1.0/3.0) - (1.0/5.0)*p2) + 
                   2.0*A2*p*(1.0 - p2 + (1.0/3.0)*p4) - 2.0*A2*p*p2*(1.0 - (2.0/3.0)*p2)
                 )*pow((1.0+p)/2.0,5.0/3.0)
               - (5.0/6.0)*alpha(p)*pow((1.0+p)/2.0,2.0/3.0) // alpha(p) * d[((1+p)/2)^(5/3)]/dn_a
             
               - ( // dalpha(-p)/dna * ((1-p)/2)^(5/3)
                   -1.0*A1*(1.0 - (2.0/3.0)*p2 + (1.0/5.0)*p4) + 4.0*A1*p2*((1.0/3.0) - (1.0/5.0)*p2) + 
                   2.0*A2*p*(1.0 - p2 + (1.0/3.0)*p4) - 2.0*A2*p*p2*(1.0 - (2.0/3.0)*p2)
                 )*pow((1.0-p)/2.0,5.0/3.0)
               + (5.0/6.0)*alpha(-1.0*p)*pow((1.0-p)/2.0,2.0/3.0) // alpha(-p) * d[((1-p)/2)^(5/3)]/dn_a
            ) * der_polarization__der_na(n_a, n_b)
            ;
}
extern "C" double der_funD__der_na_h(double n_a, double n_b)
{
    return  der_funD__der_na( n_a,  n_b);
}

/**
 * normal phase part
 * derivative: dD/dnb
 * */
__device__ __host__ inline double der_funD__der_nb(double n_a, double n_b)
{
    double p=polarization(n_a,n_b);
    double p2=p*p;
    double p4=p2*p2;
    
    return (1.0 / (2.0))*pow(6.0*M_PI*M_PI*(n_a+n_b),2.0/3.0) *
            (
              funG(p) - alpha(p)*pow((1.0+p)/2.0,5.0/3.0) - alpha(-1.0*p)*pow((1.0-p)/2.0,5.0/3.0) 
            )  +
           (1.0 / (20.0*M_PI*M_PI))*pow(6.0*M_PI*M_PI*(n_a+n_b),5.0/3.0) *
            (
               2.0*G1*p 
               - ( // dalpha(p)/dna * ((1+p)/2)^(5/3)
                   A1*(1.0 - (2.0/3.0)*p2 + (1.0/5.0)*p4) - 4.0*A1*p2*((1.0/3.0) - (1.0/5.0)*p2) + 
                   2.0*A2*p*(1.0 - p2 + (1.0/3.0)*p4) - 2.0*A2*p*p2*(1.0 - (2.0/3.0)*p2)
                 )*pow((1.0+p)/2.0,5.0/3.0)
               - (5.0/6.0)*alpha(p)*pow((1.0+p)/2.0,2.0/3.0) // alpha(p) * d[((1+p)/2)^(5/3)]/dn_a
             
               - ( // dalpha(-p)/dna * ((1-p)/2)^(5/3)
                   -1.0*A1*(1.0 - (2.0/3.0)*p2 + (1.0/5.0)*p4) + 4.0*A1*p2*((1.0/3.0) - (1.0/5.0)*p2) + 
                   2.0*A2*p*(1.0 - p2 + (1.0/3.0)*p4) - 2.0*A2*p*p2*(1.0 - (2.0/3.0)*p2)
                 )*pow((1.0-p)/2.0,5.0/3.0)
               + (5.0/6.0)*alpha(-1.0*p)*pow((1.0-p)/2.0,2.0/3.0) // alpha(-p) * d[((1-p)/2)^(5/3)]/dn_a
            ) * der_polarization__der_nb(n_a, n_b)
            ;
}
extern "C" double der_funD__der_nb_h(double n_a, double n_b)
{
    return  der_funD__der_nb( n_a,  n_b);
}

/**
 * pairing, equation (9.97) from BFM review
 * */
__device__ __host__ inline double tildeC(double n_a, double n_b, double ccoeff)
{
    double p=polarization(n_a,n_b);
    return alpha_plus(p)*pow(n_a+n_b+DENSEPSILON,1.0/3.0) / (GAMMA0*ccoeff);
}
extern "C" double tildeC_h(double n_a, double n_b)
{
    return  tildeC( n_a,  n_b, 1.0);
}

/**
 * derivative dtildeC/dna
 * */
__device__ __host__ inline double der_tildeC__der_na(double n_a, double n_b, double ccoeff)
{
    double p=polarization(n_a,n_b);
    return ( der_alpha_plus__der_na(n_a, n_b)*pow(n_a+n_b+DENSEPSILON,1.0/3.0) + (1.0/3.0)*alpha_plus(p)*pow(n_a+n_b+DENSEPSILON,-2.0/3.0) )/ (GAMMA0*ccoeff);
}
extern "C" double der_tildeC__der_na_h(double n_a, double n_b)
{
    return  der_tildeC__der_na( n_a,  n_b, 1.0);
}

/**
 * derivative dtildeC/dnb
 * */
__device__ __host__ inline double der_tildeC__der_nb(double n_a, double n_b, double ccoeff)
{
    double p=polarization(n_a,n_b);
    return ( der_alpha_plus__der_nb(n_a, n_b)*pow(n_a+n_b+DENSEPSILON,1.0/3.0) + (1.0/3.0)*alpha_plus(p)*pow(n_a+n_b+DENSEPSILON,-2.0/3.0) )/ (GAMMA0*ccoeff);
}
extern "C" double der_tildeC__der_nb_h(double n_a, double n_b)
{
    return  der_tildeC__der_nb( n_a,  n_b, 1.0);
}

#endif
