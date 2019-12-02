/**
 * This is ToolKit for pca-type codes.
 * The tool interpoates checkpoint from 3D static code to a new lattice checkpoint
* 
 * COMPILATION:
 * GW laptop
 *      mpicc tkpca_3dcheckint.c -o tkpca_3dcheckint -O3 -lm -lfftw3
 * */

#define file_operationl( cmd )                                                  \
    { ierr=cmd;                                                                 \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "FILE ERROR:: cannot execute: %s\n" , #cmd);          \
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        return( EXIT_FAILURE ) ;                                                \
    } }

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <mpi.h>

#include "pca_settings.h"
#include "pca_macro.h"
#include "pca_utils.h"
#include "pca_io.h"

// -------------------------- IMPLEMENT POTENTIAL HERE!!! ----------------------------
#undef NX
#undef NY
#undef NZ
#undef DX
#undef DY
#undef DZ

double dc_params[32];

inline double switch_function(double t, double T, double alpha)
{
    return 0.5*( 1.0+tanh( alpha*tan( M_PI_2*( 2.0*t/T-1.0 ) ) ) );
}
double borderfun(double v0, double alpha, double width, double xmax, double x) {
    if(fabs(x) <= xmax - (width + 2.)) return 0.;
    else if(fabs(x) >= xmax - 2.) return v0;
    else {
        if(x >= 0) return ( v0 * ( .5 + .5 * tanh(alpha * tan(M_PI * (x - xmax + width + 2) / width - M_PI / 2.)) ) );
        else return ( v0 * ( .5 + .5 * tanh(alpha * tan(M_PI * (-x - xmax + width + 2) / width - M_PI / 2.)) ) );
    }
}
// double u_ext(int ix, int iy, int iz, int it, int spin, int NX, int NY, int NZ, double DX, double DY, double DZ)
// {
//     double _ix = ((double)(ix) - 1.0*(NX/2))*DX;
//     double _iy = ((double)(iy) - 1.0*(NY/2))*DY;
//     double _iz = ((double)(iz) - 1.0*(NZ/2))*DZ;
//     
//     double borderpot = 0.;
//     
//     borderpot += borderfun(2., 2., 7., DX*NX / 2., _ix);
//     borderpot += borderfun(2., 2., 7., DY*NY / 2., _iy);
//     borderpot += borderfun(2., 2., 7., DZ*NZ / 2., _iz);
//     
//     // add part chemical potantial dependent
//     double mu_avg = 0.5*(dc_params[10]+dc_params[11]); // average chemical potential
//     double mu_add;
//     if(spin==SPINA)
//         mu_add=dc_params[10]-mu_avg;
//     else
//         mu_add=dc_params[11]-mu_avg;
//     
//     #define SWTCH 5.0
//     if(_ix<-SWTCH) // no modify at left
//         return borderpot;
//     else if(_ix<SWTCH) // switch region
//         return borderpot + mu_add*switch_function(_ix+SWTCH, 2.*SWTCH, 1.0);
//     else
//         return borderpot + mu_add;
//     
//     return borderpot;
//     
// //     return 0.0;
// 
// }

#define xs (0.08*(double)(LX))
#define ys (0.08*(double)(LY))
#define zs (0.08*(double)(LZ))
#define HO(x, omega2p2) (omega2p2*x*x)
double u_ext(int ix, int iy, int iz, int it, int spin, int NX, int NY, int NZ, double DX, double DY, double DZ)
{    
    double LX = DX*NX;
    double LY = DY*NY;
    double LZ = DZ*NZ;
    
    double kF=  1.0;
    double eF = 0.5*kF*kF;
    
    double Rx = 0.87*LX/2;
    double omega_x = sqrt(2.*eF / (Rx*Rx) );
    double omega_y = omega_x * (151./91);
    double omega_z = omega_x * (235./91);

    params[0] = omega_x * omega_x / 2.0;
    params[1] = omega_y * omega_y / 2.0;
    params[2] = omega_z * omega_z / 2.0;
    
//     Umax=HO(xs) + 1.0*(HO(0)-HO(xs))
    params[3] = HO((double)(-LX/2)+xs, params[0]) + 1.0*(HO((double)(-LX/2), params[0])-HO((double)(-LX/2)+xs, params[0]));
    params[4] = HO((double)(-LY/2)+ys, params[1]) + 1.0*(HO((double)(-LY/2), params[1])-HO((double)(-LY/2)+ys, params[1]));
    params[5] = HO((double)(-LZ/2)+zs, params[2]) + 1.0*(HO((double)(-LZ/2), params[2])-HO((double)(-LZ/2)+zs, params[2]));
    
    
    double _x = (double)(ix)*DX - 1.0*(LX/2);
    double _y = (double)(iy)*DY - 1.0*(LY/2);
    double _z = (double)(iz)*DZ - 1.0*(LZ/2);
    
    // ---------------------------- HO ----------------------------
    double s;
    
    double hox;
    if(_x>=(double)(-LX/2)+xs && _x<=(double)(LX/2)-xs) 
    {
        hox = HO(_x, dc_params[0]);
    }
    else if(_x<(double)(-LX/2)+xs) 
    {   
        s = switch_function(_x+(double)(LX/2), xs, 1.0);
        hox = s*HO(_x, dc_params[0]) + (1.0-s)*dc_params[3];
    }
    else
    {
        s = switch_function(_x-((double)(LX/2)-xs), xs, 1.0);
        hox=(1.0-s)*HO(_x, dc_params[0]) + s*dc_params[3];
    }
    
    double hoy;
    if(_y>=(double)(-LY/2)+ys && _y<=(double)(LY/2)-ys) 
    {
        hoy = HO(_y, dc_params[1]);
    }
    else if(_y<(double)(-LY/2)+ys) 
    {   
        s = switch_function(_y+(double)(LY/2), ys, 1.0);
        hoy = s*HO(_y, dc_params[1]) + (1.0-s)*dc_params[4];
    }
    else
    {
        s = switch_function(_y-((double)(LY/2)-ys), ys, 1.0);
        hoy=(1.0-s)*HO(_y, dc_params[1]) + s*dc_params[4];
    }
    
    double hoz;
    if(_z>=(double)(-LZ/2)+zs && _z<=(double)(LZ/2)-zs) 
    {
        hoz = HO(_z, dc_params[2]);
    }
    else if(_z<(double)(-LZ/2)+zs) 
    {   
        s = switch_function(_z+(double)(LZ/2), zs, 1.0);
        hoz = s*HO(_z, dc_params[2]) + (1.0-s)*dc_params[5];
    }
    else
    {
        s = switch_function(_z-((double)(LZ/2)-zs), zs, 1.0);
        hoz=(1.0-s)*HO(_z, dc_params[2]) + s*dc_params[5];
    }
        
    return hox + hoy + hoz;

}
#undef xs
#undef ys
#undef zs
#undef HO

// -----------------------------------------------------------------------------------

#define ZERO_NOISE 1.0e-16
#define ZERO_NOISE_CORR(array) \
    for(ixyz=0; ixyz<(outNX*outNY*outNZ); ixyz++) if(array[ixyz]<ZERO_NOISE) array[ixyz]=ZERO_NOISE;

typedef char * string;
int interpolate_array_R(int inNX, int inNY, int inNZ, double inDX, double inDY, double inDZ, double *in_array, 
                        int outNX, int outNY, int outNZ, double outDX, double outDY, double outDZ, double *out_array,
                        MPI_Comm comm
                       );
int interpolate_array_C(int inNX, int inNY, int inNZ, double inDX, double inDY, double inDZ, double complex *in_array, 
                        int outNX, int outNY, int outNZ, double outDX, double outDY, double outDZ, double complex *out_array,
                        MPI_Comm comm
                       );

int main( int argc , char ** argv ) 
{
    
    // variables
    int i, j, k; // basic iterators
    int it;
    int ix, iy, iz, ixyz; // lattice iterators
    int ierr; // error flag
    char file_name[256];
    double eF, kF, Effg, beta;
    double energy[5];
    string energy_labels[5];
    energy_labels[0] = "E_kin";
    energy_labels[1] = "E_pot";
    energy_labels[2] = "E_pair";
    energy_labels[3] = "E_CM";
    energy_labels[4] = "E_ext";
    double E_tot;
    double npart[2];
    
    int iam, np; // basic MPI indicators
    double dc_mu_a;
    double dc_mu_b;
    double dc_ec;
    
    md.overwrite=1;
    
    /* start main */
    MPI_Init( &argc , &argv ) ; /* set up the parallel WORLD */
    MPI_Comm_size( MPI_COMM_WORLD , &np ) ; /* total number of processes */
    MPI_Comm_rank( MPI_COMM_WORLD , &iam ) ; /* id of process st 0 <= iam < np */
    
    
    // INPUT OUTPUT
    if(argc!=6) if(iam==0)
    {
        printf("Usage: %s infile_name outfile_name target_lattice_size\n", argv[0]);
        printf("    inprefix - [string] input prefix\n");
        printf("    outprefix - [string] output prefix\n");
        printf("    inNX - [int] input CUBIC lattice size`\n");
        printf("    outNX - [int] output CUBIC lattice size`\n");
        printf("    LZ - [double] total length of CUBIC lattice`\n");
    }
    if(argc!=6)
    {
        /* messy exit here */
        MPI_Barrier( MPI_COMM_WORLD ) ;
        MPI_Finalize() ;
    }
    
    
    // INPUT
    char *inprefix = argv[1];
    int inNX=atoi(argv[3]);
    int inNY=inNX;
    int inNZ=inNX;
    double inDX=atof(argv[5])/inNX;
    double inDY=inDX;
    double inDZ=inDX;
    
    // OUTPUT
    char *outprefix = argv[2];
    int outNX=atoi(argv[4]);
    int outNY=outNX;
    int outNZ=outNX;
    double outDX=atof(argv[5])/outNX;
    double outDY=outDX;
    double outDZ=outDX;  


    // arrays for input
    double *h_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *h_potentials; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    cppmallocl(h_densities,12*inNX*inNY*inNZ,double);
    cppmallocl(h_potentials,4*inNX*inNY*inNZ,double);
    
    // For easier access to data
    // densities 
    double *rho_a = (double *)(h_densities +  0*inNX*inNY*inNZ);
    double *rho_b = (double *)(h_densities +  1*inNX*inNY*inNZ);
    double *tau_a = (double *)(h_densities +  2*inNX*inNY*inNZ);
    double *tau_b = (double *)(h_densities +  3*inNX*inNY*inNZ);
    double complex *nu = (double complex *)(h_densities +  4*inNX*inNY*inNZ);
    double *j_a_x = (double *)(h_densities +  6*inNX*inNY*inNZ);
    double *j_a_y = (double *)(h_densities +  7*inNX*inNY*inNZ);
    double *j_a_z = (double *)(h_densities +  8*inNX*inNY*inNZ);
    double *j_b_x = (double *)(h_densities +  9*inNX*inNY*inNZ);
    double *j_b_y = (double *)(h_densities + 10*inNX*inNY*inNZ);
    double *j_b_z = (double *)(h_densities + 11*inNX*inNY*inNZ);
    
    // pontentials
    double *V_a = (double *)(h_potentials +  0*inNX*inNY*inNZ);
    double *V_b = (double *)(h_potentials +  1*inNX*inNY*inNZ);
    double complex *delta = (double complex *)(h_potentials +  2*inNX*inNY*inNZ);
    
    
    // arrays for output
    double *out_h_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *out_h_potentials; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    cppmallocl(out_h_densities,12*outNX*outNY*outNZ,double);
    cppmallocl(out_h_potentials,4*outNX*outNY*outNZ,double);
    
    // For easier access to data
    // densities 
    double *out_rho_a = (double *)(out_h_densities +  0*outNX*outNY*outNZ);
    double *out_rho_b = (double *)(out_h_densities +  1*outNX*outNY*outNZ);
    double *out_tau_a = (double *)(out_h_densities +  2*outNX*outNY*outNZ);
    double *out_tau_b = (double *)(out_h_densities +  3*outNX*outNY*outNZ);
    double complex *out_nu = (double complex *)(out_h_densities +  4*outNX*outNY*outNZ);
    double *out_j_a_x = (double *)(out_h_densities +  6*outNX*outNY*outNZ);
    double *out_j_a_y = (double *)(out_h_densities +  7*outNX*outNY*outNZ);
    double *out_j_a_z = (double *)(out_h_densities +  8*outNX*outNY*outNZ);
    double *out_j_b_x = (double *)(out_h_densities +  9*outNX*outNY*outNZ);
    double *out_j_b_y = (double *)(out_h_densities + 10*outNX*outNY*outNZ);
    double *out_j_b_z = (double *)(out_h_densities + 11*outNX*outNY*outNZ);
    
    // pontentials
    double *out_V_a = (double *)(out_h_potentials +  0*outNX*outNY*outNZ);
    double *out_V_b = (double *)(out_h_potentials +  1*outNX*outNY*outNZ);
    double complex *out_delta = (double complex *)(out_h_potentials +  2*outNX*outNY*outNZ);

    if(iam==0)
    {
        sprintf(file_name, "%s_checkpoint.s3dpca", inprefix);
        printf("# READING CHECKPOINT FILE `%s`\n", file_name);
        printf("# INPUT LATTICE: %d x %d x %d\n", inNX, inNY, inNZ);
        printf("# INPUT SPACING: %.2f x %.2f x %.2f\n", inDX, inDY, inDZ);
        printf("# INPUT VOLUME : %.2f x %.2f x %.2f\n", inDX*inNX, inDY*inNY, inDZ*inNZ);
        
        FILE * pFile = fopen(file_name, "rb");
        
        // write all nescesary data to file
        fread(&it          , sizeof(int)         , 1 , pFile); // iteration number
        fread(&dc_mu_a     , sizeof(double)      , 1 , pFile); 
        fread(&dc_mu_b     , sizeof(double)      , 1 , pFile); 
        fread(&dc_ec       , sizeof(double)      , 1 , pFile); 
        fread(&beta        , sizeof(double)      , 1 , pFile); 
        fread(&eF          , sizeof(double)      , 1 , pFile); 
        fread(&kF          , sizeof(double)      , 1 , pFile);
        fread(&Effg        , sizeof(double)      , 1 , pFile);
        fread(h_potentials , sizeof(double)*inNX*inNY*inNZ, 4 , pFile);
        fread(h_densities  , sizeof(double)*inNX*inNY*inNZ, 12, pFile);
        fread(energy       , sizeof(double)      , 5 , pFile);
        fread(npart        , sizeof(double)      , 2 , pFile);
                
        fclose(pFile);
        
        printf("# CHECKPOINT READ: it=%d\n", it);
        printf("# CHECKPOINT READ: dc_mu_a=%16.8g  dc_mu_b=%16.8g\n", dc_mu_a, dc_mu_b);
        printf("# CHECKPOINT READ: dc_ec=%16.8g  beta=%16.8g\n", dc_ec, beta);
        printf("# CHECKPOINT READ: eF=%16.8g  kF=%16.8g  Effg=%16.8g\n", eF, kF, Effg);
        printf("# CHECKPOINT READ: ------- NPART -------\n");
        printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "SPINA", npart[SPINA], npart[SPINA], (npart[SPINA]-npart[SPINA]));
        printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "SPINB", npart[SPINB], npart[SPINB], (npart[SPINB]-npart[SPINB]));  
        printf("# CHECKPOINT READ: ------- ENERGY -------\n");
        E_tot=0.0;
        for(i=0; i<5; i++)
        {
            E_tot+=energy[i]; 
            
            printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
                energy_labels[i], energy[i]/Effg, energy[i]/Effg, (energy[i]-energy[i])/Effg);
        }
        printf("  ------------------------------------------------------------------------\n");
        printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            "E_tot", E_tot/Effg, E_tot/Effg, (E_tot-E_tot)/Effg);
    }
    
    // Send data to all processes
    MPI_Bcast(&it          , 1 , MPI_INT    , 0 , MPI_COMM_WORLD ) ;
    MPI_Bcast(&dc_mu_a     , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
    MPI_Bcast(&dc_mu_b     , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
    MPI_Bcast(&dc_ec       , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
    MPI_Bcast(&beta        , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
    MPI_Bcast(&eF          , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD ); 
    MPI_Bcast(&kF          , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    MPI_Bcast(&Effg        , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    MPI_Bcast(h_potentials , 4*inNX*inNY*inNZ , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    MPI_Bcast(h_densities  , 12*inNX*inNY*inNZ, MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    MPI_Bcast(energy       , 5 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    MPI_Bcast(npart        , 2 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    
    if(iam==0)
    {
        printf("#\n");
        printf("# TARGET LATTICE: %d x %d x %d\n", outNX, outNY, outNZ);
        printf("# TARGET SPACING: %.2f x %.2f x %.2f\n", outDX, outDY, outDZ);
        printf("# TARGET VOLUME : %.2f x %.2f x %.2f\n", outDX*outNX, outDY*outNY, outDZ*outNZ);
        printf("#\n");
    }
    
    // test - particle number
    for(i=0; i<2; i++)
    {
        double _npart=0.0;
        double *rho = (h_densities +  i*inNX*inNY*inNZ);
        for(ixyz=0; ixyz<(inNX*inNY*inNZ); ixyz++) _npart+=rho[ixyz]*inDX*inDY*inDZ;
        if(iam==0) printf("# CHECK IN: NPART[%d]=%f\n", i, _npart);
    }
    
    // ---------------------- interpolating arrays -----------------------------
    if(iam==0) printf("#\n");
    if(iam==0) printf("# INTERPOLATING TO THE NEW TARGET LATTICE...\n");
    
    // V_a
    if(iam==0) printf("# INTERPOLATING: V_a...\n");
    ixyz=0;
    for(ix=0; ix<inNX; ix++) for(iy=0; iy<inNY; iy++) for(iz=0; iz<inNZ; iz++) // subtruct external potential
    {
        V_a[ixyz]-=u_ext(ix, iy, iz, it, SPINA, inNX, inNY, inNZ, inDX, inDY, inDZ);
        ixyz++;
    }
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , V_a, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_V_a,
                        MPI_COMM_WORLD);
    ixyz=0;
    for(ix=0; ix<outNX; ix++) for(iy=0; iy<outNY; iy++) for(iz=0; iz<outNZ; iz++) // add back external potential
    {
        out_V_a[ixyz]+=u_ext(ix, iy, iz, it, SPINA, outNX, outNY, outNZ, outDX, outDY, outDZ);
        ixyz++;
    }
    
    // V_b
    if(iam==0) printf("# INTERPOLATING: V_b...\n");
    ixyz=0;
    for(ix=0; ix<inNX; ix++) for(iy=0; iy<inNY; iy++) for(iz=0; iz<inNZ; iz++) // subtruct external potential
    {
        V_b[ixyz]-=u_ext(ix, iy, iz, it, SPINB, inNX, inNY, inNZ, inDX, inDY, inDZ);
        ixyz++;
    }
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , V_b, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_V_b,
                        MPI_COMM_WORLD);
    ixyz=0;
    for(ix=0; ix<outNX; ix++) for(iy=0; iy<outNY; iy++) for(iz=0; iz<outNZ; iz++) // add back external potential
    {
        out_V_b[ixyz]+=u_ext(ix, iy, iz, it, SPINB, outNX, outNY, outNZ, outDX, outDY, outDZ);
        ixyz++;
    }
    
    // delta
    if(iam==0) printf("# INTERPOLATING: delta...\n");
    interpolate_array_C(inNX , inNY , inNZ , inDX , inDY , inDZ , delta, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_delta,
                        MPI_COMM_WORLD);

    // rho_a
    if(iam==0) printf("# INTERPOLATING: rho_a...\n");
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , rho_a, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_rho_a,
                        MPI_COMM_WORLD);
    ZERO_NOISE_CORR(out_rho_a);
    
    // rho_b
    if(iam==0) printf("# INTERPOLATING: rho_b...\n");
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , rho_b, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_rho_b,
                        MPI_COMM_WORLD);
    ZERO_NOISE_CORR(out_rho_b);
    
    // tau_a
    if(iam==0) printf("# INTERPOLATING: tau_a...\n");
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , tau_a, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_tau_a,
                        MPI_COMM_WORLD);
    ZERO_NOISE_CORR(out_tau_a);
    
    // tau_b
    if(iam==0) printf("# INTERPOLATING: tau_b...\n");
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , tau_b, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_tau_b,
                        MPI_COMM_WORLD);
    ZERO_NOISE_CORR(out_tau_b);

    // nu
    if(iam==0) printf("# INTERPOLATING: nu...\n");
    interpolate_array_C(inNX , inNY , inNZ , inDX , inDY , inDZ , nu, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_nu,
                        MPI_COMM_WORLD);
    
    // j_a_x
    if(iam==0) printf("# INTERPOLATING: j_a_x...\n");
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , j_a_x, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_j_a_x,
                        MPI_COMM_WORLD);
    
    // j_a_y
    if(iam==0) printf("# INTERPOLATING: j_a_y...\n");
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , j_a_y, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_j_a_y,
                        MPI_COMM_WORLD);
    
    // j_a_z
    if(iam==0) printf("# INTERPOLATING: j_a_z...\n");
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , j_a_z, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_j_a_z,
                        MPI_COMM_WORLD);
    
    // j_b_x
    if(iam==0) printf("# INTERPOLATING: j_b_x...\n");
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , j_b_x, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_j_b_x,
                        MPI_COMM_WORLD);
    
    // j_b_y
    if(iam==0) printf("# INTERPOLATING: j_b_y...\n");
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , j_b_y, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_j_b_y,
                        MPI_COMM_WORLD);
    
    // j_b_z
    if(iam==0) printf("# INTERPOLATING: j_b_z...\n");
    interpolate_array_R(inNX , inNY , inNZ , inDX , inDY , inDZ , j_b_z, 
                        outNX, outNY, outNZ, outDX, outDY, outDZ, out_j_b_z,
                        MPI_COMM_WORLD);
    
    // test - particle number
    for(i=0; i<2; i++)
    {
        double _npart=0.0;
        double *rho = (out_h_densities +  i*outNX*outNY*outNZ);
        for(ixyz=0; ixyz<(outNX*outNY*outNZ); ixyz++) _npart+=rho[ixyz]*outDX*outDY*outDZ;
        if(iam==0) printf("# CHECK OUT: NPART[%d]=%f\n", i, _npart);
        npart[i]=_npart;
    }
    
    dc_ec = M_PI*M_PI/(2.*outDX*outDX); 
    if(iam==0) printf("# NEW VALUE FOR dc_ec=%f\n", dc_ec);
    
    // checkpoint - only by iam==0
    if(iam==0)
    {
        sprintf(file_name, "%s_checkpoint.s3dpca", outprefix);
        printf("# CREATING CHECKPOINT FILE `%s`\n", file_name);
        FILE * pFile = fopen(file_name, "wb");
        
        // write all nescesary data to file
        i=it;
        fwrite(&i           , sizeof(int)         , 1 , pFile); // iteration number
        fwrite(&dc_mu_a     , sizeof(double)      , 1 , pFile); 
        fwrite(&dc_mu_b     , sizeof(double)      , 1 , pFile); 
        fwrite(&dc_ec       , sizeof(double)      , 1 , pFile); 
        fwrite(&beta        , sizeof(double)      , 1 , pFile); 
        fwrite(&eF          , sizeof(double)      , 1 , pFile); 
        fwrite(&kF          , sizeof(double)      , 1 , pFile);
        fwrite(&Effg        , sizeof(double)      , 1 , pFile);
        fwrite(out_h_potentials , sizeof(double)*outNX*outNY*outNZ, 4 , pFile);
        fwrite(out_h_densities  , sizeof(double)*outNX*outNY*outNZ, 12, pFile);
        fwrite(energy       , sizeof(double)      , 5 , pFile);
        fwrite(npart        , sizeof(double)      , 2 , pFile);
                
        fclose(pFile);
    }
    
    // print sections allong directions

    if(iam==0)
    {
        sprintf(file_name, "%s_sectionX.txt", outprefix);
        FILE * pFile = fopen(file_name, "w");
        printf("# CREATING FILE: `%s`\n", file_name);
        
        int nx=outNX;
        int ny=outNY;
        int nz=outNZ;
    
        iy=ny/2;
        iz=nz/2;        
        for(ix=0; ix<nx; ix++)
        {
            ixyz=ix*ny*nz + iy*nz + iz;
            
            fprintf(pFile, "%4d %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g\n",
                ix, //1
                out_rho_a[ixyz], // 2
                out_rho_b[ixyz], // 3
                cabs(out_delta[ixyz]), // 4
                out_tau_a[ixyz], // 5
                out_tau_b[ixyz], // 6
                cabs(out_nu[ixyz]), // 7
                out_V_a[ixyz], // 8
                out_V_b[ixyz] // 9
            );
        }    
        
        fclose(pFile);
    }
        
    if(iam==0) printf("# TKPCA: DONE.\n");
    

    /* messy exit here */
    MPI_Barrier( MPI_COMM_WORLD ) ;
    MPI_Finalize() ;
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}


// -----------------------------------------------------------------------------------
// ----------------------- AUXLIARY FUNCTIONS ----------------------------------------
// -----------------------------------------------------------------------------------
/**
 * Funtion returns inerpolated value of function for coordinate (x,y,z).
 * Fourier transform coefficients are reqired
 * */
double complex interpolate(double x, double y, double z, double complex *fcoeffs, double *kkx, double *kky, double *kkz, int nx, int ny, int nz)
{
    int ix, iy, iz, ixyz;
    double complex r = 0.0 + I*0.0;
    ixyz=0;
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++)
    {
        r += fcoeffs[ixyz]*cexp(I*(kkx[ix]*x + kky[iy]*y + kkz[iz]*z));
        ixyz++;
    }
    
    return r;
}

int interpolate_array_R(int inNX, int inNY, int inNZ, double inDX, double inDY, double inDZ, double *in_array, 
                        int outNX, int outNY, int outNZ, double outDX, double outDY, double outDZ, double *out_array,
                        MPI_Comm comm
                       )
{
    int iam, np;
    MPI_Comm_size( comm , &np ) ; /* total number of processes */
    MPI_Comm_rank( comm , &iam ) ; /* id of process st 0 <= iam < np */    

    int i,j;
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,inNX,double);
    cppmallocl(kky,inNY,double);
    cppmallocl(kkz,inNZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= inNX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / (inDX*inNX) * ( double ) i ; }
    j = - i ;
    for ( i = inNX / 2 ; i < inNX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / (inDX*inNX) * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= inNY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / (inDY*inNY) * ( double ) i ; }
    j = - i ;
    for ( i = inNY / 2 ; i < inNY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / (inDY*inNY) * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= inNZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / (inDZ*inNZ) * ( double ) i ; }
    j = - i ;
    for ( i = inNZ / 2 ; i < inNZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / (inDZ*inNZ) * ( double ) j ; 
        j++ ;
    }

    // Creating FFTW plans needed for interpolations
    #define USE_FFTW_PLANNER FFTW_ESTIMATE
    double complex *fft3;
    cppmallocl(fft3,inNX*inNY*inNZ,double complex);
    
    fftw_plan plan_f = fftw_plan_dft_3d(inNX, inNY, inNZ, fft3, fft3, FFTW_FORWARD, USE_FFTW_PLANNER);

    int ix, iy, iz, ixyz;    
    for(ixyz=0; ixyz<(inNX*inNY*inNZ); ixyz++) fft3[ixyz] = in_array[ixyz] + I*0.0;
    fftw_execute(plan_f);
    for(ixyz=0; ixyz<(inNX*inNY*inNZ); ixyz++) fft3[ixyz]/=(inNX*inNY*inNZ); // introduce normalization factor here
    
    // reset output array;
    for(ixyz=0; ixyz<(outNX*outNY*outNZ); ixyz++) out_array[ixyz]=0.0;
    
    // interpolate for each point
    ixyz=-1;
    for(ix=0; ix<outNX; ix++) for(iy=0; iy<outNY; iy++) for(iz=0; iz<outNZ; iz++)
    {
        ixyz++;
        if(ixyz % np != iam) continue; // not my point
        
        out_array[ixyz] = creal( interpolate(outDX*ix, outDY*iy, outDZ*iz, fft3, kkx, kky, kkz, inNX, inNY, inNZ) );
    }
    
    // take results from all processes
    MPI_Allreduce( MPI_IN_PLACE, out_array, outNX*outNY*outNZ, MPI_DOUBLE, MPI_SUM, comm);

    free(fft3);
    free(kkx);
    free(kky);
    free(kkz);
    fftw_destroy_plan(plan_f);   
    return 0;
}

int interpolate_array_C(int inNX, int inNY, int inNZ, double inDX, double inDY, double inDZ, double complex *in_array, 
                        int outNX, int outNY, int outNZ, double outDX, double outDY, double outDZ, double complex *out_array,
                        MPI_Comm comm
                       )
{
    int iam, np;
    MPI_Comm_size( comm , &np ) ; /* total number of processes */
    MPI_Comm_rank( comm , &iam ) ; /* id of process st 0 <= iam < np */    

    int i,j;
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,inNX,double);
    cppmallocl(kky,inNY,double);
    cppmallocl(kkz,inNZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= inNX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / (inDX*inNX) * ( double ) i ; }
    j = - i ;
    for ( i = inNX / 2 ; i < inNX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / (inDX*inNX) * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= inNY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / (inDY*inNY) * ( double ) i ; }
    j = - i ;
    for ( i = inNY / 2 ; i < inNY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / (inDY*inNY) * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= inNZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / (inDZ*inNZ) * ( double ) i ; }
    j = - i ;
    for ( i = inNZ / 2 ; i < inNZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / (inDZ*inNZ) * ( double ) j ; 
        j++ ;
    }

    // Creating FFTW plans needed for interpolations
    #define USE_FFTW_PLANNER FFTW_ESTIMATE
    double complex *fft3;
    cppmallocl(fft3,inNX*inNY*inNZ,double complex);
    
    fftw_plan plan_f = fftw_plan_dft_3d(inNX, inNY, inNZ, fft3, fft3, FFTW_FORWARD, USE_FFTW_PLANNER);

    int ix, iy, iz, ixyz;    
    for(ixyz=0; ixyz<(inNX*inNY*inNZ); ixyz++) fft3[ixyz] = in_array[ixyz];
    fftw_execute(plan_f);
    for(ixyz=0; ixyz<(inNX*inNY*inNZ); ixyz++) fft3[ixyz]/=(inNX*inNY*inNZ); // introduce normalization factor here
    
    // reset output array;
    for(ixyz=0; ixyz<(outNX*outNY*outNZ); ixyz++) out_array[ixyz]=0.0 + I*0.0;
    
    // interpolate for each point
    ixyz=-1;
    for(ix=0; ix<outNX; ix++) for(iy=0; iy<outNY; iy++) for(iz=0; iz<outNZ; iz++)
    {
        ixyz++;
        if(ixyz % np != iam) continue; // not my point
        
        out_array[ixyz] = interpolate(outDX*ix, outDY*iy, outDZ*iz, fft3, kkx, kky, kkz, inNX, inNY, inNZ);
    }
    
    // take results from all processes
    MPI_Allreduce( MPI_IN_PLACE, out_array, outNX*outNY*outNZ, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);

    free(fft3);
    free(kkx);
    free(kky);
    free(kkz);
    fftw_destroy_plan(plan_f);   
    
    return 0;
}
