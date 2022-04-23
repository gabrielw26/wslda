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
#include <string.h>
#include <tgmath.h>
#include <complex.h>

#include "predefines.h"
#include "pca_utils.h"
#include "gpe_utils.h"


#define Complex double complex

void set_initial_wave_function(uint nxyz, Complex *psi)
{
    uint ixyz;
    for(ixyz=0; ixyz<nxyz; ixyz++) 
    { 
        psi[ixyz] = 0.0+0.0*I; 
    }
}

void read_initial_wave_function(uint nxyz, Complex *psi)
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

void write_to_binary_file(uint nxyz, Complex *psi)
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

void write_to_txt_file(uint nx, uint ny, uint nz, Complex *psi)
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
        fprintf(fout, "%6d %12.6g\n", ix-nx/2, (creal(psi[ixyz])*creal(psi[ixyz]) + cimag(psi[ixyz])*cimag(psi[ixyz])) );
    }
    
    fprintf(fout, "\n\n");
    ix=nx/2;
    iz=nz/2;
    for(iy=0; iy<ny; iy++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iy-ny/2, (creal(psi[ixyz])*creal(psi[ixyz]) + cimag(psi[ixyz])*cimag(psi[ixyz])) );
    }
    
    fprintf(fout, "\n\n"); 
    ix=nx/2;
    iy=ny/2;
    for(iz=0; iz<nz; iz++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iz-nz/2, (creal(psi[ixyz])*creal(psi[ixyz]) + cimag(psi[ixyz])*cimag(psi[ixyz])) );
    }
    
    fclose(fout);
}

void print_header()
{
    printf("#%7s %12s %12s %12s %12s %12s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "comp.time");
}

void print_intial_results(double time, double npart, double etot, double ekin, double eint, double eext)
{
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart);  
}

void print_results(double time, double npart, double etot, double ekin, double eint, double eext, double diff, double rt)
{
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.6g %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, diff, rt);
}

void read_of_input_parameters(int argc , char ** argv)
{
    int i;
    char execcmd[ 256 ];
    strcpy( execcmd , argv[ 0 ] ) ;
    for( i = 1 ; i < argc ; i++ ) 
    {
        strcat( execcmd , " " ) ; 
        strcat( execcmd , argv[ i ] ) ;
    }
}

int parse_command_line_and_get_idx_of_input_file(int argc , char ** argv)
{
    int i = readcmd( argc , argv ) ;
    if( i == -1 )
    {
        printf( "TERMINATING! NO INPUT FILE.\n" ) ;
        return(EXIT_FAILURE);
    }
    return i;
}

void read_input_file(int idx, char ** argv)
{
    // Read input file
    // Info from file is loaded into metadata structure
    int j = parse_input_file(argv[idx]);
    if ( j == 0 )
    {
        printf("PROBLEM WITH INPUT FILE: `%s`.\n" , argv[ idx ] ) ;
        exit(EXIT_FAILURE);      
    }
        
    // Input file tags are accessible through pointer `input`
}