/***************************************************************************
 *   Copyright (C) 2015 by                                                 *
 *   WARSAW UNIVERSITY OF TECHNOLOGY                                       *
 *   FACULTY OF PHYSICS                                                    *
 *   NUCLEAR THEORY GROUP                                                  *
 *   See also AUTHORS file                                                 *
 *                                                                         *
 *   This file is a part of GPE for GPU project.                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/ 
/**
 * @file
 * @brief cuGPE library - user defined functions
 * */

#ifndef __GPE_USER_DEFINED__
#define __GPE_USER_DEFINED__

/***************************************************************************/ 
/**************************** USER DEFINED *********************************/
/***************************************************************************/

/**
 * GPE_FOR can be either PARTICLES or DIMERS. If PARTICLES then \f$\kappa=1\f$, if DIMERS then \f$\kappa=2\f$
 * */
// #define GPE_FOR PARTICLES
#define GPE_FOR DIMERS

// VARIABLES ACCESSIBLE IN THIS FILE FOR USER
// !!!!!!!!!! DO NOT MODIFY IT !!!!!!!!!
#define MAX_USER_PARAMS 32
__constant__ double d_user_param[MAX_USER_PARAMS];
__constant__ uint d_nx; // lattice size in x direction
__constant__ uint d_ny; // lattice size in y direction
__constant__ uint d_nz; // lattice size in z direction
__constant__ double d_dt;
__constant__ double d_t0;
__constant__ double d_npart;
__constant__ Complex *d_psi_ref; // pointer to reference psi on device - use gpe_set_psi_ref() function to set it

/**
 * Function changes wave function.
 * This function is called before each integration step.
 * NOTE: This function assumes that norm is not changed after modification. 
 * @param ix - x coordinate, ix=0,1,...,d_nx-1, where d_nx is global variable 
 * @param iy - y coordinate, iy=0,1,...,d_ny-1, where d_ny is global variable 
 * @param iz - z coordinate, iy=0,1,...,d_nz-1, where d_ny is global variable 
 * @param it - time value, ie. time = d_t0 + it*d_dt, d_t0 and d_dt are global variables
 * @param psi - psi(ix, iy, iz, it) - psi is normalized, i.e. int n(r) d^3r = npart, where n(r) computed according gpe_density(psi)
 * @return value of wave function after modification
 * */
inline __device__  Complex gpe_modify_psi(uint ix, uint iy, uint iz, uint it, Complex psi)
{
    return psi; // no change
}

/**
 * Switch function - performs switch in time interval [0-T]
 * */
__device__ __host__ inline double switch_function(double t, double T, double alpha)
{
    return 0.5*( 1.0+tanh( alpha*tan( M_PI_2*( 2.0*t/T-1.0 ) ) ) );
}

/**
 * Function computes value of external potential V_ext(x,y,z,t)
 * @param ix - x coordinate, ix=0,1,...,d_nx-1, where d_nx is global variable 
 * @param iy - y coordinate, iy=0,1,...,d_ny-1, where d_ny is global variable 
 * @param iz - z coordinate, iy=0,1,...,d_nz-1, where d_ny is global variable 
 * @param it - time value, ie. time = d_t0 + it*d_dt, d_t0 and d_dt are global variables
 * @return value of external potential 
 * */
inline __device__  double gpe_external_potential(uint ix, uint iy, uint iz, uint it)
{
    // potential for Josephson effect
    // see: http://arxiv.org/pdf/1508.00733v1.pdf 
#define OMEGA_YX (148.0/15.0)
#define OMEGA_ZX (187.0/15.0) 
    
    double _ix = (double)(ix) - 1.0*(NX/2)+0.5;
    double _iy = (double)(iy) - 1.0*(NY/2)+0.5;
    double _iz = (double)(iz) - 1.0*(NZ/2)+0.5;
    
    double s;
    
#define LX (1.0*NX)
#define LY (1.0*NY)
#define xs 10.0
#define ys 10.0
#define HO(x, omega2p2) (omega2p2*x*x)
    double omega_x2 = d_user_param[0];            // passed from main part - this is 0.5*omega_x*omega_x
    double omega_y2 = omega_x2*OMEGA_YX*OMEGA_YX; // passed from main part - this is 0.5*omega_x*omega_x
    double hox;
    if(_ix>=(double)(-LX/2)+xs && _ix<=(double)(LX/2)-xs) 
    {
        hox = HO(_ix, omega_x2);
    }
    else if(_ix<(double)(-LX/2)+xs) 
    {   
        s = switch_function(_ix+(double)(LX/2), xs, 1.0);
        hox = s*HO(_ix, omega_x2) + (1.0-s)*d_user_param[12];
    }
    else
    {
        s = switch_function(_ix-((double)(LX/2)-xs), xs, 1.0);
        hox=(1.0-s)*HO(_ix, omega_x2) + s*d_user_param[12];
    }
    
    double hoy;
    if(_iy>=(double)(-LY/2)+ys && _iy<=(double)(LY/2)-ys) 
    {
        hoy = HO(_iy, omega_y2);
    }
    else if(_iy<(double)(-LY/2)+ys) 
    {   
        s = switch_function(_iy+(double)(LY/2), ys, 1.0);
        hoy = s*HO(_iy, omega_y2) + (1.0-s)*d_user_param[13];
    }
    else
    {
        s = switch_function(_iy-((double)(LY/2)-ys), ys, 1.0);
        hoy=(1.0-s)*HO(_iy, omega_y2) + s*d_user_param[13];
    }
    
   
    double trap = hox + hoy;
//     double trap = omega_x2*(_ix*_ix + OMEGA_YX*OMEGA_YX*_iy*_iy + d_user_param[4]*OMEGA_ZX*OMEGA_ZX*_iz*_iz);
    
    double barrier = 0.0;
    // 1 - target value of barrier
    // 5 - activate barrier - dynamic calculations=1
    // 6 - barrier heigh for static calculaions
    // 7 - change time from value 7 to value 1
    double swtch=0.0;
    double swtch2=0.0;
    if(d_user_param[5]>0.5)
    {
        double time = d_t0 + it*d_dt;
        
        // dynamical barrier
        if(time<d_user_param[7]) swtch=switch_function(time, d_user_param[7], 1.0);
        else                     swtch=1.0;
    
        double shift=0.0;
//         barrier = (d_user_param[6] - swtch*(d_user_param[6]-d_user_param[1])) * exp(-1.0*(_ix-shift)*(_ix-shift) * (d_user_param[2]));
        barrier = (d_user_param[6] - swtch*(d_user_param[6]-d_user_param[1])) * exp(-1.0*(_ix-shift)*(_ix-shift) * (d_user_param[2] - swtch*(d_user_param[2]-d_user_param[8])));
        
        if(time<d_user_param[7]) swtch2 = 1.0; // keep tilt
        else if(time<(d_user_param[7]+d_user_param[9])) swtch2 = 1.0 -switch_function(time-d_user_param[7], d_user_param[9], 1.0); // remove tilt
        else swtch2=0.0; // no tilt
    
//         swtch2 = 1.0;
    }
    else // static calculations
    {
        double shift=0.0;
        barrier = d_user_param[6] * exp(-1.0*(_ix-shift)*(_ix-shift) * d_user_param[2]);
        swtch2 = 1.0; // tilt is on
    }
    
    return trap + barrier -_ix*d_user_param[3]*swtch2;
           
}
#undef LX 
#undef LY 
#undef xs 
#undef ys 
#undef HO


inline __device__  double gpe_external_potential_old3(uint ix, uint iy, uint iz, uint it)
{
    // potential for Josephson effect
    // see: http://arxiv.org/pdf/1508.00733v1.pdf 
#define OMEGA_YX (148.0/15.0)
#define OMEGA_ZX (187.0/15.0) 
    
    double _ix = (double)(ix) - 1.0*(NX/2)+0.5;
    double _iy = (double)(iy) - 1.0*(NY/2)+0.5;
    double _iz = (double)(iz) - 1.0*(NZ/2)+0.5;
    double omega_x2 = d_user_param[0]; // passed from main part - this is 0.5*omega_x*omega_x

    double trap = omega_x2*(_ix*_ix + OMEGA_YX*OMEGA_YX*_iy*_iy + d_user_param[4]*OMEGA_ZX*OMEGA_ZX*_iz*_iz);
    
    double barrier = 0.0;
    // 5 - activate barrier
    // 6 - start to rise of barrier
    // 7 - stop to raise barrier
    // 8 - start to rise tilt
    // 9 - stop to rise tilt
    // 10 - start to remove tilt
    // 11 - stop to remove tilt
    double swtch=0.0;
    double swtch2=0.0;
    if(d_user_param[5]>0.5)
    {
        double time = d_t0 + it*d_dt;
        
        // dynamical barrier
        if(time<d_user_param[6])      swtch=0.0;
        else if(time<d_user_param[7]) swtch=switch_function(time-d_user_param[6], d_user_param[7]-d_user_param[6], 1.0);
        else                          swtch=1.0;
        
        // tilt
        if(time<d_user_param[8])       swtch2=0.0;
        else if(time<d_user_param[9])  swtch2=switch_function(time-d_user_param[8], d_user_param[9]-d_user_param[8], 1.0);
        else if(time<d_user_param[10]) swtch2=1.0;
        else if(time<d_user_param[11]) swtch2=1.0-switch_function(time-d_user_param[10], d_user_param[11]-d_user_param[10], 1.0);
        else                           swtch2=0.0;
            
        double shift=0.0;
        barrier = swtch*d_user_param[1] * exp(-1.0*(_ix-shift)*(_ix-shift) * d_user_param[2]);
    }
    
    return trap + barrier -_ix*d_user_param[3]*swtch2;
           
}

inline __device__  double gpe_external_potential_old2(uint ix, uint iy, uint iz, uint it)
{
    // potential for Josephson effect
    // see: http://arxiv.org/pdf/1508.00733v1.pdf 
#define OMEGA_YX (148.0/15.0)
#define OMEGA_ZX (187.0/15.0) 
    
    double _ix = (double)(ix) - 1.0*(NX/2)+0.5;
    double _iy = (double)(iy) - 1.0*(NY/2)+0.5;
    double _iz = (double)(iz) - 1.0*(NZ/2)+0.5;
    double omega_x2 = d_user_param[0]; // passed from main part - this is 0.5*omega_x*omega_x

    double trap = omega_x2*(_ix*_ix + OMEGA_YX*OMEGA_YX*_iy*_iy + d_user_param[4]*OMEGA_ZX*OMEGA_ZX*_iz*_iz);
    
    double barrier = 0.0;
    // 5 - activate barrier
    // 6 - rise of barrier
    // 7 - intial shift of the barrier
    // 8 - start to remove shift
    // 9 - stop to remove shift
    if(d_user_param[5]>0.5)
    {
        double time = d_t0 + it*d_dt;
        
        // dynamical barrier
        double swtch=1.0;
        if(time<d_user_param[6]) swtch=switch_function(time, d_user_param[6], 1.0);
        
        double shift = d_user_param[7];
        if(time<d_user_param[8])
        {
            shift *= 1.0; // keep intial value
        }
        else if(time<d_user_param[9])
        {
            shift *= (1.0-switch_function(time-d_user_param[8], d_user_param[9]-d_user_param[8], 1.0));
        }
        else
        {
            shift *= 0.0;
        }
            
        barrier = swtch*d_user_param[1] * exp(-1.0*(_ix-shift)*(_ix-shift) * d_user_param[2]);
    }
    
    return trap + barrier -_ix*d_user_param[3];
           
}

inline __device__  double gpe_external_potential_old(uint ix, uint iy, uint iz, uint it)
{
    // potential for Josephson effect
    // see: http://arxiv.org/pdf/1508.00733v1.pdf 
#define OMEGA_YX (148.0/15.0)
#define OMEGA_ZX (187.0/15.0) 
    
    double _ix = (double)(ix) - 1.0*(NX/2)+0.5;
    double _iy = (double)(iy) - 1.0*(NY/2)+0.5;
    double _iz = (double)(iz) - 1.0*(NZ/2)+0.5;
    double omega_x2 = d_user_param[0]; // passed from main part - this is 0.5*omega_x*omega_x

    double trap = omega_x2*(_ix*_ix + OMEGA_YX*OMEGA_YX*_iy*_iy + d_user_param[4]*OMEGA_ZX*OMEGA_ZX*_iz*_iz);
    
//     double _tilt=(_iy*TILT_COS + _iz*TILT_SIN)*TILT_ANGLE_TAN;
//     _ix+=_tilt;
    double barrier = d_user_param[1] * exp(-1.0*(_ix)*(_ix) * d_user_param[2]);
   
    return trap + barrier -_ix*d_user_param[3];
           
}

/**
 * Function returns value of energy density functional EDF
 * @param rho - density, computed according gpe_density(psi)
 * @param it - time value, ie. time = d_t0 + it*d_dt, d_t0 and d_dt are global variables
 * @return value of energy density functional
 * */
inline __device__  double gpe_EDF(double rho, uint it)
{
    // Density energy functional for unitary Fermi gas
    // see: Phys. Rev. A 90, 043638 (2014)
    return 0.46*0.6*rho*pow(3.0*M_PI*M_PI*rho, 2.0/3.0)/2.; // unitary limit
    
//     // Density energy functional for fermionic cold atoms
//     // see: Phys. Rev. Lett. 112, 025301 (2014)
//     double a = d_user_param[10]; // scattering length
//     #define XI 0.37
//     #define CONTACT 0.901
//     double kF=pow(3.0*M_PI*M_PI*rho , 1.0/3.0);
//     if(kF<1.0e-12) return 0.0;
//     double eF=0.5*kF*kF;
//     double x=1.0/(a*kF);
//     return 0.6*eF*rho*XI*(XI+x) / ( XI + x*(1.0+CONTACT) + 3.0*M_PI*XI*x*x ) - rho/(2.0*a*a);
}

/**
 * Function returns value of mean field, i.e U= d_EDF / dn - variational derivative of EDF with respect to density
 * @param rho - density, computed according gpe_density(psi)
 * @param it - time value, ie. time = d_t0 + it*d_dt, d_t0 and d_dt are global variables
 * @return value of mean field
 * */
inline __device__  double gpe_dEDFdn(double rho, uint it)
{
    // see: Phys. Rev. A 90, 043638 (2014)
    return 0.46*pow(3.0*M_PI*M_PI*rho, 2.0/3.0)/2.0; // unitary limit
    
//     // Density energy functional for fermionic cold atoms
//     // see: Phys. Rev. Lett. 112, 025301 (2014)
//     double a = d_user_param[10]; // scattering length
//     double kF=pow(3.0*M_PI*M_PI*rho , 1.0/3.0);
//     if(kF<1.0e-12) return 0.0;
//     double eF=0.5*kF*kF;
//     double x=1.0/(a*kF);
//     double D = ( XI + x*(1.0+CONTACT) + 3.0*M_PI*XI*x*x);
//     return XI*eF*(XI+0.8*x)/D + 0.2*XI*eF*(XI+x)*x*( (1.0+CONTACT)+6.0*M_PI*XI*x )/(D*D) - 1.0/(2.0*a*a);
}


#endif
