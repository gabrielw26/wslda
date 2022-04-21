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
 * @brief cuGPE library
 * */

#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>

#ifndef __GPE_ENGINE__
#define __GPE_ENGINE__

typedef unsigned int uint;
typedef cufftDoubleComplex Complex;

// Macro allocates memory and check final status
#define gpemalloc(pointer,size,type)                                            \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        exit(1) ;                                                               \
    }
    
#define gpe_exec( cmd, ierr )                                                   \
    { ierr=cmd;                                                                 \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "error: cannot execute: %s\n", #cmd) ;                \
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        exit(1) ;                                                               \
    } }


/***************************************************************************/ 
/****************************** GPE ENGINE *********************************/
/***************************************************************************/
/**
 * Function provides sizes of mesh 
 * provided in compilation process.
 * @param *_nx size of mesh in x direction [OUTPUT]
 * @param *_ny size of mesh in y direction [OUTPUT]
 * @param *_nz size of mesh in z direction [OUTPUT]
 * */
void gpe_get_lattice(int *_nx, int *_ny, int *_nz);

/**
 * Function creates GPE engine with 1024 nthreads per block.
 * @param alpha \f$\alpha\f$ parameter of GPE equation [INPUT]
 * @param beta \f$\beta\f$ of GPE equation [INPUT]
 * @param dt integration step [INPUT]
 * @param npart number of particles [INPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_create_engine(double alpha, double beta, double dt, double npart , int nthreads=1024);

/**
 * Function changes values of GPE parameters
 * @param alpha \f$\alpha\f$ parameter of GPE equation [INPUT]
 * @param beta \f$\beta\f$ of GPE equation [INPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_change_alpha_beta(double alpha, double beta);

/**
 * Function sets new value of time for present wave function, i.e.  \f$\Psi(t)\rightarrow\Psi(t_0)\f$
 * @param t0 new value of time [INPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_set_time(double t0);

/**
 * Function sets user parameters. Maximal number of user parameters is 32.
 * User parameters are accessible in functions gpe_modify_psi(), gpe_external_potential(), gpe_EDF() and gpe_dEDFdn()
 * through array d_user_param.
 * @param size number of parameters, not bigger than 32 [INPUT]
 * @param params pointer to array with parameters [INPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_set_user_params(int size, double *params);

/**
 * Functions sets value of quantum friction coefficient.
 * Note - quantum friction produces extra cost in computation process - overhead is about 50%! 
 * @param qfcoeff - quantum friction coefficient \f$\gamma\f$. If qfcoeff==0.0 then quantum friction is deactivated. [INPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_set_quantum_friction_coeff(double qfcoeff);

/**
 * Function sets wave function for specified time i.e. \f$\Psi(t)\f$
 * @param t value of time [INPUT]
 * @param psi corresponding wave function [INPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_set_psi(double t, Complex * psi);

/**
 * Function returns wave function and corresponding time i.e. \f$\Psi(t)\f$
 * @param t value of time [OUTPUT]
 * @param psi corresponding wave function [OUTPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_get_psi(double *t, Complex * psi);

/**
 * Function copies wave function to d_psi_ref buffer.
 * This is read only buffer, and it is not modified during computation process.
 * User has access to psi_ref through d_psi_ref that is visible in each function listed in gpe_user_defined.
 * @param psi wave function [INPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_set_psi_ref(Complex * psi);

/**
 * Function copies present wave function (holds by evolver) to d_psi_ref buffer.
 * This is read only buffer, and it is not modified during computation process.
 * User has access to psi_ref through d_psi_ref that is visible in each function listed in gpe_user_defined.
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_set_psi_ref_from_present_state();

/**
 * Function removes all data from buffer d_psi_ref.
 * After this operation buffer d_psi_ref is set to NULL.
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_clear_psi_ref();

/**
 * Function returns density computed from wave function \f$n(\vec{r}, t)=\kappa |\Psi(\vec{r}, t)|^{2}\f$ and corresponding time.
 * \f$\kappa\f$ is equal 1 for particles and 2 for dimers.
 * @param t value of time [OUTPUT]
 * @param density corresponding density [OUTPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_get_density(double *t, double * density);

/**
 * Function returns currents computed from wave function \f$\vec{j}(\vec{r}, t)=\frac{1}{\kappa}\textrm{Im}[\Psi^*(\vec{r}, t)\vec{\nabla}\Psi(\vec{r}, t)]\f$ 
 * and corresponding time.
 * \f$\kappa\f$ is equal 1 for particles and 2 for dimers.
 * @param t value of time [OUTPUT]
 * @param jx x coordinate of currents i.e. j_x(r,t) [OUTPUT]
 * @param jy y coordinate of currents i.e. j_y(r,t) [OUTPUT]
 * @param jz z coordinate of currents i.e. j_z(r,t) [OUTPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_get_currents(double *t, double * jx, double * jy, double * jz);

/**
 * Function normalizes state
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_normalize_psi();

/**
 * Function evolves state nt steps in time i.e. \f$\Psi(t)\rightarrow\Psi(t+n_t dt)\f$
 * @param nt number of steps to evolve [INPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_evolve(int nt);

/**
 * Function evolves state nt steps in imaginary time using Taylor expansion up to given order.
 * In addition term -Omega*Lz is added to the hamiltonian (support for rotating frames)
 * @param nt number of steps to evolve [INPUT]
 * @param order order of Taylor expansion [INPUT]
 * @param Omega angular momentum for frame rotation
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_evolve_taylor(int nt, int order, double Omega);

/**
 * Function returns energy computed from wave function and corresponding time.
 * Total energy is equal etot=ekin+eint+eext.
 * @param t value of time [OUTPUT]
 * @param ekin expectation value of kinetic energy operator [OUTPUT]
 * @param eint expectation value of internal energy [OUTPUT]
 * @param eext expectation value of external energy [OUTPUT]
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_energy(double *t, double *ekin, double *eint, double *eext);

/**
 * Function destroys engine.
 * @return It returns 0 if success otherwise error code is returned.
 * */
int gpe_destroy_engine();
#endif
