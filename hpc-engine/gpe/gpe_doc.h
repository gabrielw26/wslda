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

// This is only page with main documentation code
// ---------------------------------------------------------------------------
/*! \mainpage cuGPE - accelerated engine for solving GPE-like equation

<h2>Description.</h2>

cuGPE library solves dimensionless (\f$\hbar=m=1\f$) GPE-like equation:
\f[
(i\alpha - \beta)\frac{\partial\Psi(\vec{r}, t)}{\partial t} =  \left(  -\frac{1}{2\kappa}\nabla^2 + \kappa\frac{\delta\mathcal{E}(n,t)}{\delta n} + \kappa V_{ext}(\vec{r}, t)\right) \Psi(\vec{r}, t)
\f]
where \f$\Psi\f$ is collective wave function of system, \f$\mathcal{E}\f$ is energy density functional and \f$V_{ext}\f$ is external potential.

Energy density functional defines internal energy of the system \f$E_{int}=\int d^{3}\vec{r}\;\mathcal{E}(n(\vec{r}),t)\f$

Parameter \f$\kappa\f$ can be:
- \f$\kappa=1\f$: \f$\Psi\f$  represents collective wave function of particles and density is defined as \f$n(\vec{r}, t)=|\Psi(\vec{r}, t)|^{2}\f$.
- \f$\kappa=2\f$: \f$\Psi\f$  represents collective wave function of dimers with effective mass equal 2 and density is defined as \f$n(\vec{r}, t)=2|\Psi(\vec{r}, t)|^{2}\f$.

In order to perform simulation user has to modify gpe_user_defined.h file and provide:
- energy density functional \f$\mathcal{E}(n,t)\f$: see function gpe_EDF()
- derivative of energy density functional \f$\frac{\delta\mathcal{E}(n,t)}{\delta n}\f$: see function gpe_dEDFdn()
- external potential \f$V_{ext}(\vec{r}, t)\f$: see function gpe_external_potential()

The GPE-like problem is solved on 3D spacial lattice of size \f$n_x \times n_y \times n_z\f$.
Presently the lattice size  is limited to \f$512^3\f$ (due to memory limitation).

<hr>

<h2>License</h2>
cuGPE library is released under the <a href="http://www.gnu.org/copyleft/gpl.html">GNU General Public License</a>.
It means that the library is totally <b>free</b>!

<hr>

<h2>Authors</h2>
- <b>Gabriel Wlazlowski</b>, 
<a href="http://eng.pw.edu.pl/">Warsaw University of Technology</a>, 
<a href="http://www.if.pw.edu.pl/www_en/">Faculty of Physics</a>,
<a href="mailto:gabrielw@if.pw.edu.pl">gabrielw@if.pw.edu.pl</a> 
(main developer) 
- <a href="http://www.if.pw.edu.pl/~magiersk/"><b>Piotr Magierski</b></a>,
<a href="http://eng.pw.edu.pl/">Warsaw University of Technology</a>, 
<a href="http://www.if.pw.edu.pl/www_en/">Faculty of Physics</a>, 
<a href="mailto:magiersk@if.pw.edu.pl">magiersk@if.pw.edu.pl</a> 
(meritorical support) 

<hr>

<h2>Requirements</h2>
cuGPE library uses CUDA technology. 
To download CUDA Toolkit visit <a href="https://developer.nvidia.com/cuda-zone">CUDA Zone</a> web page. 

<hr>

<h2>User's guide</h2>
In order to build a code following steps have to be accomplished:
1. Modify gpe_user_defined.h file
2. Create Main code
3. Compile cuGPE library and Main code

Examples below present how to prepare cuGPE engine, compile and run the code.
All examples are provided in Examples directory.

\subpage user_defined page presents content of gpe_user_defined.h file.

1. Imaginary time projection: \subpage gpe_imag 
2. Real time evolution: \subpage gpe_real
3. Ground state projection via real time evolution: \subpage gpe_real_quantum_friction
4. Real time evolution with removing energy of phonons: \subpage gpe_real_psi_ref

 */

//-----------------------------------------------------------
/*! \page user_defined User defined functions

For presented examples gpe_user_defined.h file has follwing content:

\code 
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
 * GPE_FOR can be either PARTICLES or DIMERS. If PARTICLES then kappa=1, if DIMERS then kappa=2
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

\endcode

 */

//-----------------------------------------------------------
/*! \page gpe_imag gpe_imag.cu

\section gpe_imag_main Content of the main code file gpe_imag.cu

\code
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
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <complex>      // std::complex

/***************************************************************************/ 
/**************************** GPE HEADERS **********************************/
/***************************************************************************/
#include "gpe_engine.h"
// Only to do timing
#include "gpe_timing.h" 
#include "gpe_user_defined.h"

/***************************************************************************/ 
/************************* MAIN FUNCTION  **********************************/
/***************************************************************************/
int main( int argc , char ** argv ) 
{
    // SETTINGS
    double alpha=0.0;
    double beta=1.0;
    double dt=0.025;
    double npart=1000.0;
    int device=0;
    
    int ierr;
    cudaError err;
    
    err=cudaSetDevice( device );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot cudaSetDevice(%d)!\n", device);
        return 1;
    }
    
    int nx, ny, nz;
    gpe_get_lattice(&nx, &ny, &nz);
    printf("# GPE engine compiled for lattice: %d x %d x %d\n", nx, ny, nz);
    
    uint nxyz=nx*ny*nz;
    uint ix, iy, iz, ixyz;
    
    double ekin, eint, eext, etot, etot_prev;
    double time;
    double rt;
  
    // ***************** Imaginary time projection **************************
    printf("# IMAGINARY TIME PROJECTION\n");
    
    // CPU memory for wave function - pinned for fast transfers
    Complex *psi; // Complex type defined in gpe_engine.h - structure with two doubles x and y for real and imaginary parts
    err=cudaHostAlloc( &psi , sizeof(Complex)*nxyz, cudaHostAllocDefault );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot allocate memory!\n");
        return 1;
    }
    
    // Set initial wave function - constant
    for(ixyz=0; ixyz<nxyz; ixyz++) { psi[ixyz].x = 1.0; psi[ixyz].y = 0.0; } // Fill with initial values
    
    // Create engine
    gpe_exec( gpe_create_engine(alpha, beta, dt, npart), ierr );
    
    // Prepare user defined parameters
    double params[3];
    params[0] = 0.01; // omega_x
    params[1] = 0.10; // omega_y
    params[2] = 0.11; // omega_z
    
    // Copy parameters to engine
    gpe_exec( gpe_set_user_params(3, params), ierr );
    
    // Copy wave function to GPU and normalize
    gpe_exec( gpe_set_psi(0.0, psi), ierr ) ;
    gpe_exec( gpe_normalize_psi(), ierr );

    // For nice printing
    printf("#%7s %12s %12s %12s %12s %12s %14s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "diff", "comp.time");
    
    // initial energy
    gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
    etot = ekin + eint + eext;
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart);  
    
    // evolve in imaginary time until convergence is achieved
    while(1)
    {
        b_t(); // reset timer
        
        // Evolve 100 steps forward
        gpe_exec( gpe_evolve(100), ierr ); 
        
        // Compute energy 
        gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
        
        rt = e_t(); // get time
        
        etot_prev=etot;
        etot = ekin + eint + eext;
        double diff=(etot_prev-etot)/npart; // diference in energy per particle
        printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.6g %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, diff, rt);  
        
        if(fabs(diff)<1.0e-9) break;
    }    
    
    // Get wave function and save it to file
    gpe_exec( gpe_get_psi(&time, psi), ierr ) ;    
            
    // write to binary file
    printf("# Writing psi to file\n");
    FILE * psiFile;
    psiFile = fopen ("psi.dat", "wb");
    fwrite (psi , sizeof(Complex)*nxyz, 1, psiFile);
    fclose (psiFile);   
    
    // Write to txt file |Psi(x,y,x)|^2 - sections along x, y and z axis
    FILE * fout = fopen("psi.txt", "w");
    iy=ny/2;
    iz=nz/2;
    for(ix=0; ix<nx; ix++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", ix-nx/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n");
    ix=nx/2;
    iz=nz/2;
    for(iy=0; iy<ny; iy++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iy-ny/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n"); 
    ix=nx/2;
    iy=ny/2;
    for(iz=0; iz<nz; iz++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iz-nz/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fclose(fout);
    
    // Destroy engine
    gpe_exec( gpe_destroy_engine(), ierr) ;
    
    // Clear memory
    cudaFreeHost(psi);
    
    return 0;
}
\endcode

\section gpe_imag_compile Compilation

Compilation commands are presented belowe.
\code
nvcc -c gpe_engine.cu -o gpe_engine.o -O3 -arch sm_35 -Xcudafe '--diag_suppress=declared_but_not_referenced --diag_suppress=set_but_not_used' -DNX=256 -DNY=32 -DNZ=32                                            
nvcc -c gpe_imag.cu -o gpe_imag.o -O3 -arch sm_35 -Xcudafe '--diag_suppress=declared_but_not_referenced --diag_suppress=set_but_not_used' -DNX=256 -DNY=32 -DNZ=32                                                
nvcc gpe_imag.o gpe_engine.o -o gpe_imag -lcudart -lcufft -lm
\endcode

Important remarks:
- Lattice size has to be provieded during compilation procees (performace related issue). 
In presented example lattice will have dimension 256 x 32 x 32 (note flags -DNX=256 -DNY=32 -DNZ=32)
- It is important to provide correct `compute capability` of target device - in presented example the code is compiled for GPU cards 
with compute capability 3.5 (note flag -arch sm_35)

\section gpe_imag_output Output of execution
Ouput of the code executed on K80 gpu card
\code
[gabrielw@node2062 gpe]$ ./gpe_imag 
# GPE engine compiled for lattice: 256 x 32 x 32
# IMAGINARY TIME PROJECTION
#   time         etot         ekin         eint         eext  (eint+eext)           diff    comp.time
    0.00   1.24378644   0.00000000   0.02593644   1.21785000   1.24378644
    2.50   0.34147110   0.00488094   0.11813150   0.21845867   0.33659017     0.902315       0.1814
    5.00   0.32571292   0.00650935   0.14323551   0.17596806   0.31920357    0.0157582       0.1639
    7.50   0.32416932   0.00715695   0.14974779   0.16726457   0.31701236   0.00154361       0.1553
   10.00   0.32388702   0.00743251   0.15212530   0.16432921   0.31645451  0.000282296       0.1495
   12.50   0.32381498   0.00754983   0.15315475   0.16311041   0.31626515  7.20384e-05       0.1459
   15.00   0.32379110   0.00759931   0.15364293   0.16254886   0.31619179  2.38784e-05       0.1434
   17.50   0.32378133   0.00761978   0.15389013   0.16227141   0.31616155  9.77498e-06       0.1435
   20.00   0.32377666   0.00762790   0.15402315   0.16212561   0.31614876  4.66764e-06       0.1435
   22.50   0.32377418   0.00763082   0.15409921   0.16204416   0.31614337  2.47579e-06       0.1435
   25.00   0.32377277   0.00763160   0.15414539   0.16199579   0.31614117  1.41083e-06       0.1435
   27.50   0.32377193   0.00763156   0.15417501   0.16196535   0.31614037  8.45833e-07       0.1437
   30.00   0.32377140   0.00763123   0.15419498   0.16194519   0.31614017  5.26448e-07       0.1435
   32.50   0.32377106   0.00763085   0.15420899   0.16193122   0.31614021  3.37123e-07       0.1435
   35.00   0.32377084   0.00763049   0.15421917   0.16192119   0.31614036  2.20707e-07       0.1435
   37.50   0.32377070   0.00763017   0.15422676   0.16191376   0.31614052  1.47026e-07       0.1435
   40.00   0.32377060   0.00762991   0.15423254   0.16190814   0.31614068   9.9304e-08       0.1435
   42.50   0.32377053   0.00762970   0.15423703   0.16190380   0.31614083  6.78172e-08       0.1439
   45.00   0.32377048   0.00762952   0.15424057   0.16190039   0.31614096   4.6728e-08       0.1436
   47.50   0.32377045   0.00762938   0.15424338   0.16189769   0.31614107  3.24298e-08       0.1435
   50.00   0.32377043   0.00762926   0.15424565   0.16189551   0.31614116  2.26389e-08       0.1435
   52.50   0.32377041   0.00762917   0.15424749   0.16189375   0.31614124  1.58797e-08       0.1435
   55.00   0.32377040   0.00762909   0.15424899   0.16189232   0.31614131  1.11825e-08       0.1436
   57.50   0.32377039   0.00762903   0.15425022   0.16189114   0.31614137  7.90012e-09       0.1435
   60.00   0.32377039   0.00762897   0.15425124   0.16189017   0.31614141  5.59613e-09       0.1437
   62.50   0.32377038   0.00762893   0.15425209   0.16188937   0.31614145  3.97285e-09       0.1435
   65.00   0.32377038   0.00762889   0.15425279   0.16188870   0.31614149  2.82561e-09       0.1435
   67.50   0.32377038   0.00762886   0.15425338   0.16188814   0.31614152  2.01273e-09       0.1435
   70.00   0.32377038   0.00762883   0.15425387   0.16188767   0.31614154  1.43554e-09       0.1435
   72.50   0.32377038   0.00762881   0.15425428   0.16188728   0.31614156  1.02496e-09       0.1435
   75.00   0.32377037   0.00762879   0.15425463   0.16188695   0.31614158  7.32459e-10       0.1438
# Writing psi to file

\endcode
 */

//-----------------------------------------------------------
/*! \page gpe_real gpe_real.cu

\section gpe_real_main Content of the main code file gpe_real.cu

\code
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
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <complex>      // std::complex

/***************************************************************************/ 
/**************************** GPE HEADERS **********************************/
/***************************************************************************/
#include "gpe_engine.h"
// Only to do timing
#include "gpe_timing.h" 
#include "gpe_user_defined.h"

/***************************************************************************/ 
/************************* MAIN FUNCTION  **********************************/
/***************************************************************************/
int main( int argc , char ** argv ) 
{
    // SETTINGS
    double alpha=1.0;
    double beta=0.0;
    double dt=0.025;
    double npart=1000.0;
    int device=0;
    
    int ierr;
    cudaError err;
    
    err=cudaSetDevice( device );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot cudaSetDevice(%d)!\n", device);
        return 1;
    }
    
    int nx, ny, nz;
    gpe_get_lattice(&nx, &ny, &nz);
    printf("# GPE engine compiled for lattice: %d x %d x %d\n", nx, ny, nz);
    
    uint nxyz=nx*ny*nz;
    uint ix, iy, iz, ixyz;
    
    double ekin, eint, eext, etot;
    double time;
    double rt;
  
    // ***************** Imaginary time projection **************************
    printf("# REAL TIME EVOLUTION\n");
    
    // CPU memory for wave function - pinned for fast transfers
    Complex *psi; // Complex type defined in gpe_engine.h - structure with two doubles x and y for real and imaginary parts
    err=cudaHostAlloc( &psi , sizeof(Complex)*nxyz, cudaHostAllocDefault );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot allocate memory!\n");
        return 1;
    }
    
    // Set initial wave function - read it from file psi.dat (see gpe_imag.cu)
    printf("# Reading psi from file\n");
    FILE * psiFile;
    psiFile = fopen ("psi.dat", "rb");
    size_t readok = fread (psi , sizeof(Complex)*nxyz, 1, psiFile);
    if (readok != 1)
    {
        printf("Reading error\n");
        return 1;
    }
    fclose (psiFile);  
    
    // Create engine
    gpe_exec( gpe_create_engine(alpha, beta, dt, npart), ierr );
    
    // Prepare user defined parameters
    double params[3];
    params[0] = 0.01; // omega_x
    params[1] = 0.10; // omega_y
    params[2] = 0.11; // omega_z
    
    // Copy parameters to engine
    gpe_exec( gpe_set_user_params(3, params), ierr );
    
    // Copy wave function to GPU
    gpe_exec( gpe_set_psi(0.0, psi), ierr ) ;

    // For nice printing
    printf("#%7s %12s %12s %12s %12s %12s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "comp.time");
    
    // initial energy
    gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
    etot = ekin + eint + eext;
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart);  
    
    // evolve in real time
    while(1)
    {
        b_t(); // reset timer
        
        // Evolve 100 steps forward
        gpe_exec( gpe_evolve(100), ierr ); 
        
        // Compute energy 
        gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
        
        rt = e_t(); // get time
        
        etot = ekin + eint + eext;

        printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, rt);  
        
        if(time>100.) break;
    }    
    
    // Get wave function and save it to file
    gpe_exec( gpe_get_psi(&time, psi), ierr ) ;    
            
    // write to binary file
    printf("# Writing psi to file\n");
    psiFile = fopen ("psi2.dat", "wb");
    fwrite (psi , sizeof(Complex)*nxyz, 1, psiFile);
    fclose (psiFile);   
    
    // Write to txt file |Psi(x,y,x)|^2 - sections along x, y and z axis
    FILE * fout = fopen("psi2.txt", "w");
    iy=ny/2;
    iz=nz/2;
    for(ix=0; ix<nx; ix++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", ix-nx/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n");
    ix=nx/2;
    iz=nz/2;
    for(iy=0; iy<ny; iy++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iy-ny/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n"); 
    ix=nx/2;
    iy=ny/2;
    for(iz=0; iz<nz; iz++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iz-nz/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fclose(fout);
    
    // Destroy engine
    gpe_exec( gpe_destroy_engine(), ierr) ;
    
    // Clear memory
    cudaFreeHost(psi);
    
    return 0;
}
\endcode

\section gpe_real_compile Compilation

Compilation commands are presented belowe.
\code
nvcc -c gpe_engine.cu -o gpe_engine.o -O3 -arch sm_35 -Xcudafe '--diag_suppress=declared_but_not_referenced --diag_suppress=set_but_not_used' -DNX=256 -DNY=32 -DNZ=32
nvcc -c gpe_real.cu -o gpe_real.o -O3 -arch sm_35 -Xcudafe '--diag_suppress=declared_but_not_referenced --diag_suppress=set_but_not_used' -DNX=256 -DNY=32 -DNZ=32
nvcc gpe_real.o gpe_engine.o -o gpe_real -lcudart -lcufft -lm
\endcode
See also remarks: \ref gpe_imag_compile

\section gpe_real_output Output of execution
Ouput of the code executed on K80 gpu card
\code
[gabrielw@node2062 gpe]$ ./gpe_real 
# GPE engine compiled for lattice: 256 x 32 x 32
# REAL TIME EVOLUTION
# Reading psi from file
#   time         etot         ekin         eint         eext  (eint+eext)    comp.time
    0.00   0.32377037   0.00762879   0.15425463   0.16188695   0.31614158
    2.50   0.32377037   0.00762879   0.15425462   0.16188696   0.31614158       0.1555
    5.00   0.32377037   0.00762879   0.15425461   0.16188698   0.31614159       0.1509
    7.50   0.32377037   0.00762879   0.15425459   0.16188700   0.31614159       0.1395
   10.00   0.32377037   0.00762879   0.15425456   0.16188702   0.31614158       0.1296
   12.50   0.32377037   0.00762880   0.15425454   0.16188704   0.31614158       0.1241
   15.00   0.32377037   0.00762880   0.15425454   0.16188703   0.31614157       0.1212
   17.50   0.32377037   0.00762880   0.15425456   0.16188701   0.31614157       0.1186
   20.00   0.32377037   0.00762880   0.15425460   0.16188698   0.31614158       0.1163
   22.50   0.32377037   0.00762879   0.15425466   0.16188692   0.31614158       0.1166
   25.00   0.32377037   0.00762879   0.15425472   0.16188686   0.31614158       0.1155
   27.50   0.32377037   0.00762879   0.15425479   0.16188680   0.31614159       0.1154
   30.00   0.32377037   0.00762879   0.15425485   0.16188674   0.31614159       0.1155
   32.50   0.32377037   0.00762879   0.15425489   0.16188670   0.31614159       0.1154
   35.00   0.32377037   0.00762879   0.15425492   0.16188666   0.31614158       0.1154
   37.50   0.32377037   0.00762879   0.15425494   0.16188664   0.31614158       0.1158
   40.00   0.32377037   0.00762879   0.15425496   0.16188663   0.31614159       0.1155
   42.50   0.32377037   0.00762878   0.15425498   0.16188662   0.31614160       0.1154
   45.00   0.32377037   0.00762877   0.15425501   0.16188660   0.31614161       0.1155
   47.50   0.32377037   0.00762876   0.15425505   0.16188657   0.31614161       0.1158
   50.00   0.32377037   0.00762876   0.15425510   0.16188652   0.31614162       0.1145
   52.50   0.32377037   0.00762875   0.15425517   0.16188645   0.31614162       0.1155
   55.00   0.32377037   0.00762874   0.15425526   0.16188637   0.31614163       0.1158
   57.50   0.32377037   0.00762873   0.15425536   0.16188628   0.31614164       0.1156
   60.00   0.32377037   0.00762872   0.15425545   0.16188620   0.31614165       0.1154
   62.50   0.32377037   0.00762872   0.15425554   0.16188612   0.31614166       0.1155
   65.00   0.32377037   0.00762871   0.15425561   0.16188605   0.31614166       0.1154
   67.50   0.32377037   0.00762871   0.15425567   0.16188600   0.31614167       0.1155
   70.00   0.32377037   0.00762870   0.15425571   0.16188596   0.31614167       0.1155
   72.50   0.32377037   0.00762870   0.15425575   0.16188593   0.31614168       0.1154
   75.00   0.32377037   0.00762869   0.15425579   0.16188590   0.31614169       0.1155
   77.50   0.32377037   0.00762868   0.15425584   0.16188586   0.31614170       0.1158
   80.00   0.32377037   0.00762867   0.15425590   0.16188581   0.31614171       0.1158
   82.50   0.32377037   0.00762866   0.15425598   0.16188574   0.31614172       0.1145
   85.00   0.32377037   0.00762865   0.15425606   0.16188566   0.31614172       0.1154
   87.50   0.32377037   0.00762865   0.15425616   0.16188556   0.31614172       0.1155
   90.00   0.32377037   0.00762865   0.15425627   0.16188546   0.31614172       0.1155
   92.50   0.32377037   0.00762865   0.15425637   0.16188536   0.31614173       0.1155
   95.00   0.32377037   0.00762865   0.15425646   0.16188527   0.31614173       0.1158
   97.50   0.32377037   0.00762865   0.15425654   0.16188519   0.31614173       0.1155
  100.00   0.32377037   0.00762864   0.15425659   0.16188514   0.31614173       0.1154
  102.50   0.32377037   0.00762864   0.15425664   0.16188510   0.31614173       0.1155
# Writing psi to file

\endcode
 */

//-----------------------------------------------------------
/*! \page gpe_real_quantum_friction gpe_real_quantum_friction.cu
\section gpe_real_quantum_friction_thoery Quantum friction - theory
Quantum friction potential is given by:
\f[
 U_{qf}(\vec{r},t) = -\gamma\frac{\vec{\nabla}\cdot\vec{j}(\vec{r},t)}{n(\vec{r},t)}
\f]
where \f$n(\vec{r},t)\f$ and \f$\vec{j}(\vec{r},t)\f$ stand for particle density and particle current respectively.

To activate quantm friction set value of \f$\gamma\f$ to non-zero value. To deactive quantum friction set \f$\gamma\f$ to zero - see: gpe_set_quantum_friction_coeff()

To learn more about quantum friction see:
- Aurel Bulgac, Michael McNeil Forbes, Kenneth J. Roche, Gabriel Wlazlowski, <i>Quantum Friction: Cooling Quantum Systems with Unitary Time Evolution</i>
<a href="http://arxiv.org/abs/1305.6891">arXiv:1305.6891</a>

\section gpe_real_quantum_friction_main Content of the main code file gpe_real_quantum_friction.cu

\code
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
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <complex>      // std::complex

/***************************************************************************/ 
/**************************** GPE HEADERS **********************************/
/***************************************************************************/
#include "gpe_engine.h"
// Only to do timing
#include "gpe_timing.h" 
#include "gpe_user_defined.h"

/***************************************************************************/ 
/************************* MAIN FUNCTION  **********************************/
/***************************************************************************/
int main( int argc , char ** argv ) 
{
    // SETTINGS
    double alpha=1.0;
    double beta=0.0;
    double dt=0.025;
    double npart=1000.0;
    int device=0;
    
    int ierr;
    cudaError err;
    
    err=cudaSetDevice( device );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot cudaSetDevice(%d)!\n", device);
        return 1;
    }
    
    int nx, ny, nz;
    gpe_get_lattice(&nx, &ny, &nz);
    printf("# GPE engine compiled for lattice: %d x %d x %d\n", nx, ny, nz);
    
    uint nxyz=nx*ny*nz;
    uint ix, iy, iz, ixyz;
    
    double ekin, eint, eext, etot, etot_prev;
    double time;
    double rt;
  
    // ***************** Imaginary time projection **************************
    printf("# GROUND STATE PROJECTION VIA REAL TIME EVOLUTION\n");
    
    // CPU memory for wave function - pinned for fast transfers
    Complex *psi; // Complex type defined in gpe_engine.h - structure with two doubles x and y for real and imaginary parts
    err=cudaHostAlloc( &psi , sizeof(Complex)*nxyz, cudaHostAllocDefault );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot allocate memory!\n");
        return 1;
    }
    
    // Set initial wave function - constant
    for(ixyz=0; ixyz<nxyz; ixyz++) { psi[ixyz].x = 1.0; psi[ixyz].y = 0.0; } // Fill with initial values
    
    // Create engine
    gpe_exec( gpe_create_engine(alpha, beta, dt, npart), ierr );
    
    // Prepare user defined parameters
    double params[3];
    params[0] = 0.01; // omega_x
    params[1] = 0.10; // omega_y
    params[2] = 0.11; // omega_z
    
    // Copy parameters to engine
    gpe_exec( gpe_set_user_params(3, params), ierr );
    
    // Copy wave function to GPU and normalize
    gpe_exec( gpe_set_psi(0.0, psi), ierr ) ;
    gpe_exec( gpe_normalize_psi(), ierr );
    
    // Activate quantum friction
    gpe_exec( gpe_set_quantum_friction_coeff (5.0), ierr );

    // For nice printing
    printf("#%7s %12s %12s %12s %12s %12s %14s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "diff", "comp.time");
    
    // initial energy
    gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
    etot = ekin + eint + eext;
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart);  
    
    // evolve in imaginary time until convergence is achieved
    while(1)
    {
        b_t(); // reset timer
        
        // Evolve 1000 steps forward
        gpe_exec( gpe_evolve(1000), ierr ); 
        
        // Compute energy 
        gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
        
        rt = e_t(); // get time
        
        etot_prev=etot;
        etot = ekin + eint + eext;
        double diff=(etot_prev-etot)/npart; // diference in energy per particle
        printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.6g %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, diff, rt);  
        
        if(fabs(diff)<1.0e-9) break;
    }    
    
    // Get wave function and save it to file
    gpe_exec( gpe_get_psi(&time, psi), ierr ) ;    
            
    // write to binary file
    printf("# Writing psi to file\n");
    FILE * psiFile;
    psiFile = fopen ("psiqf.dat", "wb");
    fwrite (psi , sizeof(Complex)*nxyz, 1, psiFile);
    fclose (psiFile);   
    
    // Write to txt file |Psi(x,y,x)|^2 - sections along x, y and z axis
    FILE * fout = fopen("psiqf.txt", "w");
    iy=ny/2;
    iz=nz/2;
    for(ix=0; ix<nx; ix++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", ix-nx/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n");
    ix=nx/2;
    iz=nz/2;
    for(iy=0; iy<ny; iy++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iy-ny/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n"); 
    ix=nx/2;
    iy=ny/2;
    for(iz=0; iz<nz; iz++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iz-nz/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fclose(fout);
    
    // Destroy engine
    gpe_exec( gpe_destroy_engine(), ierr) ;
    
    // Clear memory
    cudaFreeHost(psi);
    
    return 0;
}
\endcode

\section gpe_real_quantum_friction_compile Compilation

Compilation commands are presented belowe.
\code
nvcc -c gpe_engine.cu -o gpe_engine.o -O3 -arch sm_35 -Xcudafe '--diag_suppress=declared_but_not_referenced --diag_suppress=set_but_not_used' -DNX=256 -DNY=32 -DNZ=32
nvcc -c gpe_real_quantum_friction.cu -o gpe_real_qf.o -O3 -arch sm_35 -Xcudafe '--diag_suppress=declared_but_not_referenced --diag_suppress=set_but_not_used' -DNX=256 -DNY=32 -DNZ=32
nvcc gpe_real_qf.o gpe_engine.o -o gpe_real_qf -lcudart -lcufft -lm
\endcode
See also remarks: \ref gpe_imag_compile

\section gpe_real_quantum_friction_output Output of execution
Ouput of the code executed on K80 gpu card
\code
[gabrielw@node2062 gpe]$ ./gpe_real_qf 
# GPE engine compiled for lattice: 256 x 32 x 32
# GROUND STATE PROJECTION VIA REAL TIME EVOLUTION
#   time         etot         ekin         eint         eext  (eint+eext)           diff    comp.time
    0.00   1.24378644   0.00000000   0.02593644   1.21785000   1.24378644
   25.00   0.48578325   0.02201187   0.10627596   0.35749542   0.46377137     0.758003       2.2036
   50.00   0.47574053   0.05211403   0.11399788   0.30962862   0.42362650    0.0100427       2.1064
   75.00   0.46704667   0.08757238   0.12143789   0.25803640   0.37947429   0.00869386       2.1077
  100.00   0.45384939   0.11288317   0.13324496   0.20772125   0.34096622    0.0131973       2.1078
  125.00   0.43467789   0.11251891   0.15277489   0.16938409   0.32215898    0.0191715       2.1083
  150.00   0.41201298   0.07948604   0.18159876   0.15092818   0.33252694    0.0226649       2.1083
  175.00   0.39424576   0.03796596   0.20612113   0.15015867   0.35627980    0.0177672       2.1084
  200.00   0.37756614   0.01072003   0.21255218   0.15429393   0.36684611    0.0166796       2.1088
  225.00   0.36641927   0.01383815   0.19862689   0.15395423   0.35258112    0.0111469       2.1088
  250.00   0.35413241   0.01778710   0.18487965   0.15146566   0.33634531    0.0122869       2.1086
  275.00   0.34289881   0.01704471   0.17435156   0.15150254   0.32585410    0.0112336       2.1086
  300.00   0.33470700   0.01392432   0.16660198   0.15418070   0.32078269   0.00819181       2.1095
  325.00   0.33006086   0.01151410   0.16061579   0.15793097   0.31854676   0.00464614       2.1093
  350.00   0.32733507   0.00986607   0.15607393   0.16139508   0.31746901   0.00272579       2.1089
  375.00   0.32567387   0.00877891   0.15300790   0.16388705   0.31689496    0.0016612       2.1090
  400.00   0.32480603   0.00816223   0.15132975   0.16531405   0.31664380  0.000867835       2.1095
  425.00   0.32444986   0.00791717   0.15076761   0.16576508   0.31653269  0.000356172       2.1090
  450.00   0.32432479   0.00790968   0.15101790   0.16539720   0.31641511  0.000125069       2.1094
  475.00   0.32424462   0.00797026   0.15180286   0.16447151   0.31627437  8.01671e-05       2.1091
  500.00   0.32414652   0.00797122   0.15287301   0.16330228   0.31617530  9.81034e-05       2.1096
  525.00   0.32404084   0.00788772   0.15398890   0.16216422   0.31615312  0.000105681       2.1094
  550.00   0.32395272   0.00776101   0.15494957   0.16124214   0.31619171  8.81192e-05       2.1091
  575.00   0.32389732   0.00764797   0.15562092   0.16062843   0.31624935  5.53971e-05       2.1097
  600.00   0.32387272   0.00758594   0.15594684   0.16033994   0.31628678  2.46031e-05       2.1095
  625.00   0.32386467   0.00757836   0.15594478   0.16034153   0.31628632  8.04546e-06       2.1095
  650.00   0.32385813   0.00760372   0.15568609   0.16056832   0.31625441  6.54754e-06       2.1096
  675.00   0.32384566   0.00763527   0.15526860   0.16094179   0.31621038  1.24694e-05       2.1100
  700.00   0.32382818   0.00765594   0.15479175   0.16138048   0.31617224  1.74773e-05       2.1091
  725.00   0.32381061   0.00766210   0.15433987   0.16180864   0.31614851  1.75678e-05       2.1098
  750.00   0.32379728   0.00765904   0.15397366   0.16216458   0.31613824  1.33357e-05       2.1094
  775.00   0.32378963   0.00765401   0.15372792   0.16240770   0.31613562  7.64899e-06       2.1098
  800.00   0.32378644   0.00765129   0.15361316   0.16252198   0.31613514  3.18911e-06       2.1137
  825.00   0.32378530   0.00765105   0.15361944   0.16251482   0.31613426   1.1343e-06       2.1143
  850.00   0.32378418   0.00765087   0.15372166   0.16241164   0.31613331  1.12903e-06       2.1141
  875.00   0.32378215   0.00764834   0.15388562   0.16224819   0.31613381    2.023e-06       2.1096
  900.00   0.32377943   0.00764270   0.15407435   0.16206239   0.31613674  2.72103e-06       2.1099
  925.00   0.32377674   0.00763509   0.15425400   0.16188765   0.31614165  2.69361e-06       2.1142
  950.00   0.32377471   0.00762765   0.15439852   0.16174854   0.31614706  2.02575e-06       2.1172
  975.00   0.32377357   0.00762230   0.15449236   0.16165890   0.31615127  1.14395e-06       2.1102
 1000.00   0.32377310   0.00761992   0.15453104   0.16162214   0.31615318  4.66819e-07       2.1199
 1025.00   0.32377292   0.00762024   0.15451971   0.16163297   0.31615268  1.81643e-07       2.1214
 1050.00   0.32377270   0.00762231   0.15447036   0.16168003   0.31615039  2.20457e-07       2.1123
 1075.00   0.32377232   0.00762505   0.15439844   0.16174883   0.31614727  3.82107e-07       2.1189
 1100.00   0.32377183   0.00762765   0.15431967   0.16182452   0.31614418  4.84731e-07       2.1151
 1125.00   0.32377138   0.00762973   0.15424749   0.16189416   0.31614165   4.5426e-07       2.1234
 1150.00   0.32377106   0.00763119   0.15419153   0.16194834   0.31613987  3.22828e-07       2.1220
 1175.00   0.32377088   0.00763206   0.15415697   0.16198186   0.31613883  1.70836e-07       2.1165
 1200.00   0.32377082   0.00763239   0.15414466   0.16199377   0.31613843  6.50968e-08       2.1184
 1225.00   0.32377079   0.00763223   0.15415194   0.16198662   0.31613856  2.76243e-08       2.1160
 1250.00   0.32377075   0.00763163   0.15417373   0.16196539   0.31613912  3.98928e-08       2.1222
 1275.00   0.32377068   0.00763072   0.15420377   0.16193620   0.31613997  6.68608e-08       2.1203
 1300.00   0.32377060   0.00762966   0.15423580   0.16190515   0.31614095  8.08061e-08       2.1223
 1325.00   0.32377053   0.00762864   0.15426451   0.16187738   0.31614189  7.26678e-08       2.1225
 1350.00   0.32377048   0.00762783   0.15428618   0.16185647   0.31614265  4.96062e-08       2.1211
 1375.00   0.32377046   0.00762734   0.15429890   0.16184421   0.31614312   2.5065e-08       2.1199
 1400.00   0.32377045   0.00762718   0.15430255   0.16184071   0.31614327  9.21247e-09       2.1220
 1425.00   0.32377044   0.00762731   0.15429841   0.16184473   0.31614314  4.61039e-09       2.1208
 1450.00   0.32377044   0.00762763   0.15428870   0.16185411   0.31614281  7.48996e-09       2.1204
 1475.00   0.32377042   0.00762805   0.15427604   0.16186634   0.31614237  1.19336e-08       2.1222
 1500.00   0.32377041   0.00762848   0.15426297   0.16187896   0.31614193  1.36953e-08       2.1177
 1525.00   0.32377040   0.00762885   0.15425158   0.16188998   0.31614155  1.17639e-08       2.1218
 1550.00   0.32377039   0.00762912   0.15424325   0.16189802   0.31614127  7.65988e-09       2.1218
 1575.00   0.32377039   0.00762927   0.15423865   0.16190246   0.31614111  3.66353e-09       2.1226
 1600.00   0.32377039   0.00762931   0.15423771   0.16190336   0.31614108  1.30093e-09       2.1239
 1625.00   0.32377038   0.00762925   0.15423984   0.16190130   0.31614114  7.89323e-10       2.1214
# Writing psi to file

\endcode
 */


//-----------------------------------------------------------
/*! \page gpe_real_psi_ref gpe_real_psi_ref.cu
\section gpe_real_psi_ref_thoery Removing energy of phonons
In GPE the velocity field is fully determined by phase of wave-function:
\f$v(\vec{r},t) = \frac{\hbar}{m}\nabla\phi(\vec{r},t)\f$, where \f$\phi\f$ is phase i.e. \f$\psi=\sqrt{|\psi|^2}e^{i\phi}\f$.
Thus, it possible to generate from given state \f$\psi_{\textrm{ref}}\f$ a new state \f$\psi\f$ that it is characterized 
by EXACTLY the same velocity field as \f$\psi_{\textrm{ref}}\f$ but with lowest possible energy.
Physically, it means that from \f$\psi_{\textrm{ref}}\f$ we remove all phonon excitations.

From technical point of view procedure reduces to:
1. Fill buffer d_psi_ref with \f$\psi_{\textrm{ref}}\f$. To do this use gpe_set_psi_ref() or gpe_set_psi_ref_from_present_state().
2. Switch evolver to imaginary time projection.
3. In each step erase phase of present wave function \f$\psi\f$ and replace it by phase of \f$\psi_{\textrm{ref}}\f$. To do this use gpe_modify_psi().

Example below shows time evolution of quantized vortex line which is initially bend (it has shape of Gaussian function).
During time evolution vortex emits phonons - they are removed after emission thus vortex dynamics is not affected by interaction with phonons.

\section gpe_real_psi_ref_user_def Definition of user defined functions

The vortex in confined in external trap of Woods-Saxon shape:
\code
inline __device__  double gpe_external_potential(uint ix, uint iy, uint iz, uint it)
{
    
    // Woods-Saxon potential:
    // V(r) = -V_0 / ( 1+exp( (r-R)/a ) )
    
    // parameters V_0, R and a are passed from main code
    double V_0 = d_user_param[0];
    double R   = d_user_param[1];
    double a   = d_user_param[2];
    
    // coordinate with respect to center of the box
    double _ix = (double)(ix) - 1.0*(NX/2);
    double _iy = (double)(iy) - 1.0*(NY/2);
    double _iz = (double)(iz) - 1.0*(NZ/2);
    
    double _r = sqrt(_ix*_ix + _iy*_iy); // distance from center of external potential in plane XY


    double trap =  -1.0*V_0 / ( 1+exp( (_r-R)/a ) );
    
    return trap;
}
\endcode

gpe_modify_psi function allows for three modes of working. 
1. No change of wave function - it is used for standard real time evolution
2. Phase of wave-function is set to rotate around line with Gaussian shape. This mode is used to generate initial state (via imaginary time evolution) for real time dynamics.
3. Phase of wave-function is set to be identical with phase of psi_ref. 

To decide what mode should be used function uses d_user_param[3] element.

\code
inline __device__  Complex gpe_modify_psi(uint ix, uint iy, uint iz, uint it, Complex psi)
{
    // NOTE: I use d_user_param[3] element to switch between different options
    
    if(d_user_param[3]==0.0) // real time evolution, no modification of psi
    {
        return psi; // no change
    }
    else if(d_user_param[3]==1.0) // imaginary time evolution - I imprint vortex along z of Gaussian shape
    {
        double psi_abs = sqrt(psi.x*psi.x + psi.y*psi.y); // absolute value
        double _ix = (double)(ix) - 1.0*(NX/2);
        double _iz = (double)(iz) - 1.0*(NZ/2);
        
        // NOTE: d_user_param[4] and d_user_param[5] are passed from main code
        double _iy_shift=d_user_param[4] * exp(-1.0 * pow(_iz/d_user_param[5],2)); //Gaussian shape
        double _iy = (double)(iy) - 1.0*(NY/2) + _iy_shift;
        
        double psi_arg = atan2(_iy, _ix);
        
        // phase imprint
        psi.x = psi_abs*cos(psi_arg);
        psi.y = psi_abs*sin(psi_arg);
        
        return psi;
    }
    else // imaginary time evolution with frozen phase
    {
        uint ixyz= iz + NZ*iy + NZ*NY*ix;
        double psi_abs = sqrt(psi.x*psi.x + psi.y*psi.y); // absolute value
        double psi_arg = atan2(d_psi_ref[ixyz].y, d_psi_ref[ixyz].x); // take phase from psi_ref
        
        // phase imprint
        psi.x = psi_abs*cos(psi_arg);
        psi.y = psi_abs*sin(psi_arg);
        
        return psi;
    }
}
\endcode

\section gpe_real_psi_ref_main Content of the main code file gpe_real_psi_ref.cu

\code
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
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <complex>      // std::complex

/***************************************************************************/ 
/**************************** GPE HEADERS **********************************/
/***************************************************************************/
#include "gpe_engine.h"
// Only to do timing
#include "gpe_timing.h" 
#include "gpe_user_defined.h"


/***************************************************************************/ 
/************************* MAIN FUNCTION  **********************************/
/***************************************************************************/
int main( int argc , char ** argv ) 
{
    // SETTINGS
    double alpha=0.0;
    double beta=1.0;
    double dt=0.025;
    double npart=3000.0;
    int device=0;
    
    int ierr;
    cudaError err;
    
    err=cudaSetDevice( device );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot cudaSetDevice(%d)!\n", device);
        return 1;
    }
    
    int nx, ny, nz;
    gpe_get_lattice(&nx, &ny, &nz);
    printf("# GPE engine compiled for lattice: %d x %d x %d\n", nx, ny, nz);
    
    uint nxyz=nx*ny*nz;
    uint ix, iy, iz, ixyz;
    
    double ekin, eint, eext, etot, etot_prev, etot_phonons=0.0;
    double time;
    double rt;

    // ***************** Imaginary time projection **************************
    printf("# IMAGINARY TIME PROJECTION\n");
    
    // CPU memory for wave function - pinned for fast transfers
    Complex *psi; // Complex type defined in gpe_engine.h - structure with two doubles x and y for real and imaginary parts
    err=cudaHostAlloc( &psi , sizeof(Complex)*nxyz, cudaHostAllocDefault );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot allocate memory!\n");
        return 1;
    }
    
    // Set initial wave function - constant
    for(ixyz=0; ixyz<nxyz; ixyz++) { psi[ixyz].x = 1.0; psi[ixyz].y = 0.0; } // Fill with initial values
    
    // Create engine
    gpe_exec( gpe_create_engine(alpha, beta, dt, npart), ierr );
    
    // Prepare user defined parameters
    double params[6];
    params[0] = 10.0;    // Woods-Saxon potential: V_0
    params[1] = 0.90*nx/2; // Woods-Saxon potential: R
    params[2] = 0.02*nx/2; // Woods-Saxon potential: a
    params[3] = 1.0; // gpe_modify_psi working mode - imaginary time evolution with vortex along z of Gaussian shape
    params[4] = 10.0;// Gaussian shape: amplitude
    params[5] = 10.0; // Gaussian shape: width
    
    // Copy parameters to engine
    gpe_exec( gpe_set_user_params(6, params), ierr );
    
    // Copy wave function to GPU and normalize
    gpe_exec( gpe_set_psi(0.0, psi), ierr ) ;
    gpe_exec( gpe_normalize_psi(), ierr );

    // For nice printing
    printf("#%7s %12s %12s %12s %12s %12s %14s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "diff", "comp.time");
    
    // initial energy
    gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
    etot = ekin + eint + eext;
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart);  
    
    // evolve in imaginary time until convergence is achieved
    while(1)
    {
        b_t(); // reset timer
        
        // Evolve 100 steps forward
        gpe_exec( gpe_evolve(100), ierr ); 
        
        // Compute energy 
        gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
        
        rt = e_t(); // get time
        
        etot_prev=etot;
        etot = ekin + eint + eext;
        double diff=(etot_prev-etot)/npart; // diference in energy per particle
        printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.6g %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, diff, rt);  
        
        if(fabs(diff)<1.0e-9) break;
    }    
    
    // Get wave function and save it to file
    gpe_exec( gpe_get_psi(&time, psi), ierr ) ;    
            
    // ***************** Real time evolution **************************
    printf("# REAL TIME EVOLUTION\n");
        
    // reset time;
    gpe_exec( gpe_set_time(0.0), ierr );
    
    // For nice printing
    printf("#%7s %12s %12s %12s %12s\n", "time", "etot", "vor.lgth", "ephonons", "comp.time");
    
    // initial energy
    gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
    etot = ekin + eint + eext;
    
    // Compute vortex length - CAN BE IMPROVED!!! COMPUTATION WITH LOW ACCURACY!!
    int vorx[nz]; // cordinate of vortex - x coordinate
    int vory[nz]; // cordinate of vortex - x coordinate
    double vorlgth;
    double vortmp; // temporary variables for vortex length computation
    
    for(iz=0; iz<nz; iz++) // for each plane find point with minimal density - inside tube of given radius
    {
        vorx[iz]=-1; vory[iz]=-1;
        vortmp=999999.0; // minimal density
        for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++)
        {
            if( sqrt(pow((1.0*ix-nx/2),2) +  pow((1.0*iy-ny/2),2))>0.70*nx/2 ) continue; // out of tube
            
            ixyz= iz + nz*iy + nz*ny*ix;
            if( (pow(psi[ixyz].x,2) +  pow(psi[ixyz].y,2))<vortmp )
            {
                vortmp=pow(psi[ixyz].x,2) +  pow(psi[ixyz].y,2);
                vorx[iz]=ix; vory[iz]=iy;
            }
        }
       
        if(vorx[iz]==-1 || vory[iz]==-1) return -1; // error! 
    }
    vorlgth=0.0;
    for(iz=0; iz<nz; iz++) vorlgth+=sqrt( pow(1.0*(vorx[(iz+1)%nz]-vorx[iz]), 2 ) + pow(1.0*(vory[(iz+1)%nz]-vory[iz]), 2 ) + 1.0 );
        
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, vorlgth, 0.0 , 0.0);  
    
    // evolve in real time
    while(1)
    {
        b_t(); // reset timer
        
        // --------- REAL TIME EVOLUTION ---------
        // Change engine - switch to real time evolution
        alpha=1.0;
        beta=0.0;
        gpe_exec( gpe_change_alpha_beta(alpha, beta), ierr );
        
        // Copy parameters to engine
        params[3] = 0.0; // gpe_modify_psi working mode - real time evolution, no modification of psi
        gpe_exec( gpe_set_user_params(6, params), ierr );
            
        // Evolve 100 steps forward
        gpe_exec( gpe_evolve(100), ierr ); 
        
        // Compute energy 
        gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
        etot_prev = ekin + eint + eext;
        
        // --------- REMOVE PHONONS ---------
        // set psi_ref
        gpe_exec( gpe_set_psi_ref_from_present_state(), ierr ); 
        
        // Change engine - switch to imaginary time evolution
        alpha=0.0;
        beta=1.0;
        gpe_exec( gpe_change_alpha_beta(alpha, beta), ierr );
        
        // Copy parameters to engine
        params[3] = 2.0; // gpe_modify_psi working mode - // imaginary time evolution with frozen phase
        gpe_exec( gpe_set_user_params(6, params), ierr );
        
        gpe_exec( gpe_evolve(100), ierr ); 
        // reset time;
        gpe_exec( gpe_set_time(time), ierr );
        
        // Compute energy 
        gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
        
        etot = ekin + eint + eext;
        etot_phonons+=(etot_prev-etot); // amount of removed energy
        
        // --------- COMPUTE OBSERVABLES AND REPORT IT ---------

        // Get wave function
        gpe_exec( gpe_get_psi(&time, psi), ierr ) ;  
    
        // compute vortet length - CAN BE IMPROVED!!! COMPUTATION WITH LOW ACCURACY!!
        for(iz=0; iz<nz; iz++) // for each plane find point with minimal density - inside tube of given radius
        {
            vorx[iz]=-1; vory[iz]=-1;
            vortmp=999999.0; // minimal density
            for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++)
            {
                if( sqrt(pow((1.0*ix-nx/2),2) +  pow((1.0*iy-ny/2),2))>0.70*nx/2 ) continue; // out of tube
                
                ixyz= iz + nz*iy + nz*ny*ix;
                if( (pow(psi[ixyz].x,2) +  pow(psi[ixyz].y,2))<vortmp )
                {
                    vortmp=pow(psi[ixyz].x,2) +  pow(psi[ixyz].y,2);
                    vorx[iz]=ix; vory[iz]=iy;
                }
            }
        
            if(vorx[iz]==-1 || vory[iz]==-1) return -1; // error! 
        }
        vorlgth=0.0;
        for(iz=0; iz<nz; iz++) vorlgth+=sqrt( pow(1.0*(vorx[(iz+1)%nz]-vorx[iz]), 2 ) + pow(1.0*(vory[(iz+1)%nz]-vory[iz]), 2 ) + 1.0 );
    
        rt = e_t(); // get time

        // print
        printf("%8.2f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, vorlgth, etot_phonons/npart , rt);  
        
        if(time>100.) break;
    }    
    
    // Get wave function and save it to file
    gpe_exec( gpe_get_psi(&time, psi), ierr ) ;  
    
    // write to binary file
    FILE * psiFile;
    psiFile = fopen ("psi.dat", "wb");
    fwrite (psi , sizeof(Complex)*nxyz, 1, psiFile);
    fclose (psiFile); 
    
    // Destroy engine
    gpe_exec( gpe_destroy_engine(), ierr) ;
    
    // Clear memory
    cudaFreeHost(psi);
    
    return 0;
}
\endcode

\section gpe_real_psi_ref_compile Compilation

Compilation commands are presented belowe.
\code
nvcc -c gpe_engine.cu -o gpe_engine.o -O3 -arch sm_35 -Xcudafe '--diag_suppress=declared_but_not_referenced --diag_suppress=set_but_not_used' -DNX=64 -DNY=64 -DNZ=64
nvcc -c gpe_real_psi_ref.cu -o gpe_real_psi_ref.o -O3 -arch sm_35 -Xcudafe '--diag_suppress=declared_but_not_referenced --diag_suppress=set_but_not_used' -DNX=64 -DNY=64 -DNZ=64
nvcc gpe_real_psi_ref.o gpe_engine.o -o gpe_real_psi_ref -lcudart -lcufft -lm
\endcode
See also remarks: \ref gpe_imag_compile

\section gpe_real_psi_ref_output Output of execution
Ouput of the code executed on K80 gpu card
\code
gabrielw@w540:~/MyProjects/gpe/cugpe$ ./gpe_real_psi_ref 0
# GPE engine compiled for lattice: 64 x 64 x 64
# IMAGINARY TIME PROJECTION
#   time         etot         ekin         eint         eext  (eint+eext)           diff    comp.time
    0.00  -6.31753201   0.00000000   0.05394997  -6.37148198  -6.31753201
    2.50  -9.90604430   0.00442490   0.08335542  -9.99382461  -9.91046919      3.58851       4.6916
    5.00  -9.90701418   0.00417225   0.08478067  -9.99596711  -9.91118644  0.000969887       4.4904
    7.50  -9.90713887   0.00406551   0.08532607  -9.99653045  -9.91120438  0.000124685       4.5772
   10.00  -9.90716040   0.00401907   0.08555679  -9.99673626  -9.91117947  2.15304e-05       4.3940
   12.50  -9.90716449   0.00399861   0.08565744  -9.99682054  -9.91116310  4.08731e-06       4.4167
   15.00  -9.90716529   0.00398956   0.08570192  -9.99685678  -9.91115485  8.06186e-07       4.6579
   17.50  -9.90716546   0.00398555   0.08572171  -9.99687271  -9.91115101  1.61817e-07       4.6365
   20.00  -9.90716549   0.00398377   0.08573054  -9.99687980  -9.91114926  3.27737e-08       4.4153
   22.50  -9.90716550   0.00398299   0.08573450  -9.99688298  -9.91114848  6.67535e-09       4.4714
   25.00  -9.90716550   0.00398264   0.08573627  -9.99688440  -9.91114813  1.36645e-09       4.4406
   27.50  -9.90716550   0.00398248   0.08573707  -9.99688504  -9.91114798  2.81615e-10       4.4317
# REAL TIME EVOLUTION
#   time         etot     vor.lgth     ephonons    comp.time
    0.00  -9.90716550  72.28427125   0.00000000   0.00000000
    2.50  -9.90716810  72.28427125   0.00000261   7.36917210
    5.00  -9.90717191  72.28427125   0.00000642   7.38789415
    7.50  -9.90717578  72.28427125   0.00001029   7.40928316
   10.00  -9.90717937  72.28427125   0.00001388   7.43843508
   12.50  -9.90718268  72.28427125   0.00001718   7.40415907
   15.00  -9.90718573  72.28427125   0.00002023   7.41789412
   17.50  -9.90718857  72.28427125   0.00002307   7.43635702
   20.00  -9.90719122  72.28427125   0.00002572   7.51943588
   22.50  -9.90719372  72.09151861   0.00002823   7.60622287
   25.00  -9.90719609  72.09151861   0.00003059   8.05484700
   27.50  -9.90719834  72.28427125   0.00003285   7.41869903
   30.00  -9.90720050  72.09151861   0.00003500   7.35465884
   32.50  -9.90720256  72.09151861   0.00003706   7.33631301
   35.00  -9.90720452  72.09151861   0.00003903   7.40490294
   37.50  -9.90720641  72.09151861   0.00004091   7.44458199
   40.00  -9.90720822  72.09151861   0.00004272   7.62996197
   42.50  -9.90720995  71.26309149   0.00004446   7.54241300
   45.00  -9.90721163  71.45584412   0.00004613   7.52452588
   47.50  -9.90721323  71.45584412   0.00004774   7.62443304
   50.00  -9.90721478  71.45584412   0.00004929   7.77522206
   52.50  -9.90721628  71.45584412   0.00005078   7.68651390
   55.00  -9.90721773  71.26309149   0.00005223   7.68222690
   57.50  -9.90721913  71.26309149   0.00005363   7.88170815
   60.00  -9.90722047  71.26309149   0.00005497   7.43006992
   62.50  -9.90722176  72.72719310   0.00005627   7.42997503
   65.00  -9.90722302  72.53444047   0.00005752   7.42708707
   67.50  -9.90722423  72.53444047   0.00005874   7.40821600
   70.00  -9.90722541  72.72719310   0.00005991   7.44067621
   72.50  -9.90722655  72.91994574   0.00006105   7.47543311
   75.00  -9.90722766  72.91994574   0.00006216   7.43699312
   77.50  -9.90722873  72.09151861   0.00006323   7.46207309
   80.00  -9.90722977  72.09151861   0.00006427   7.66153908
   82.50  -9.90723077  72.09151861   0.00006528   7.41722894
   85.00  -9.90723175  72.09151861   0.00006625   7.70406914
   87.50  -9.90723270  72.09151861   0.00006720   7.46598911
   90.00  -9.90723361  72.09151861   0.00006812   7.77018404
   92.50  -9.90723451  72.09151861   0.00006901   7.35911393
   95.00  -9.90723537  72.09151861   0.00006987   7.35751104
   97.50  -9.90723621  72.09151861   0.00007071   7.35847187
  100.00  -9.90723702  72.09151861   0.00007152   7.35284185
  102.50  -9.90723781  71.89876598   0.00007231   7.35719395
\endcode

 */