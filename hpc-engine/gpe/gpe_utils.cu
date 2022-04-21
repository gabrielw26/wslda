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
#include <complex>

#include "predefines.h"
#include "pca_utils.h"
#include "gpe_engine.h"

/**
 * Function set gpu device.
 * @param device It is gpu machine number on which the program will be executed.
 * */
 extern "C" void set_gpu_device(int device)
{
    cudaError err=cudaSetDevice( device );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot cudaSetDevice(%d)!\n", device);
        exit(err);
    }
}

/**
 * Function allocate CPU memory for wave function. As a result it is pinned for fast transfers.
 * @param nxyz It is product of nx, ny and nz lattice.
 * @param psi It is structure with two doubles x and y for real and imaginary parts.
 * */
 extern "C" void alloc_host_memory(uint nxyz, Complex **psi)
{
    cudaError err=cudaHostAlloc( psi , sizeof(Complex)*nxyz, cudaHostAllocDefault );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot allocate memory!\n");
        exit(err);
    }
}

extern "C" void free_host_memory(void *psi)
{
    cudaFreeHost(psi);
}

extern "C" void set_initial_wave_function(uint nxyz, Complex *psi)
{
    uint ixyz;
    for(ixyz=0; ixyz<nxyz; ixyz++) 
    { 
        psi[ixyz].x = 1.0; psi[ixyz].y = 0.0; 
    }
}

extern "C" void read_initial_wave_function(uint nxyz, Complex *psi)
{
    FILE * psiFile;
    printf("# Reading psi from file\n");
    psiFile = fopen ("psi.dat", "rb");
    size_t readok = fread (psi , sizeof(Complex)*nxyz, 1, psiFile);
    if (readok != 1)
    {
        printf("Reading error\n");
        exit(1);
    }
    fclose (psiFile); 
}

extern "C" void write_to_binary_file(uint nxyz, Complex *psi)
{
    FILE * psiFile;
    printf("# Writing psi to file\n");

    //TODO 
    // imag -> "psi.dat"
    // real -> "psi2.dat"
    psiFile = fopen ("psi.dat", "wb");
    fwrite (psi , sizeof(Complex)*nxyz, 1, psiFile);
    fclose (psiFile);    
}

extern "C" void write_to_txt_file(uint nx, uint ny, uint nz, Complex *psi)
{
    FILE * fout;
    uint ix, iy, iz, ixyz;

    //TODO 
    // imag -> "psi.dat"
    // real -> "psi2.dat"
    fout = fopen("psi.txt", "w");

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
}

extern "C" void print_header()
{
    printf("#%7s %12s %12s %12s %12s %12s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "comp.time");
}

extern "C" void print_intial_results(double time, double npart, double etot, double ekin, double eint, double eext)
{
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart);  
}

extern "C" void print_results(double time, double npart, double etot, double ekin, double eint, double eext, double diff, double rt)
{
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.6g %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, diff, rt);
}