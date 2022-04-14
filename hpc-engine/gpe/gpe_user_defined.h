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
 * Function computes value of external potential V_ext(x,y,z,t)
 * @param ix - x coordinate, ix=0,1,...,d_nx-1, where d_nx is global variable 
 * @param iy - y coordinate, iy=0,1,...,d_ny-1, where d_ny is global variable 
 * @param iz - z coordinate, iy=0,1,...,d_nz-1, where d_ny is global variable 
 * @param it - time value, ie. time = d_t0 + it*d_dt, d_t0 and d_dt are global variables
 * @return value of external potential 
 * */
inline __device__  double gpe_external_potential(uint ix, uint iy, uint iz, uint it)
{
    
    // harmonic trap:
    // V(x,y,z) = 0.5*(omega_x*x)^2 + 0.5*(omega_y*y)^2 + 0.5*(omega_z*z)^2
    
    // frequencies are passed through d_user_param array
    double omega_x = d_user_param[0];
    double omega_y = d_user_param[1];
    double omega_z = d_user_param[2];
    
    // coordinate with respect to center of the box
    double _ix = (double)(ix) - 1.0*(NX/2);
    double _iy = (double)(iy) - 1.0*(NY/2);
    double _iz = (double)(iz) - 1.0*(NZ/2);


    double trap =   0.5*_ix*_ix*omega_x*omega_x
                  + 0.5*_iy*_iy*omega_y*omega_y 
                  + 0.5*_iz*_iz*omega_z*omega_z;
    
    return trap;
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
    return 0.37*0.6*rho*pow(3.0*M_PI*M_PI*rho, 2.0/3.0)/2.; // unitary limit
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
    return 0.37*pow(3.0*M_PI*M_PI*rho, 2.0/3.0)/2.0; // unitary limit
}


#endif
