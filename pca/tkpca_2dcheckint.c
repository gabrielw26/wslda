/**
 * This is ToolKit for pca-type codes.
 * The tool interpoates checkpoint from 2D static code to a new lattice checkpoint
* 
 * COMPILATION:
 * GW laptop
 *      mpicc tkpca_2dcheckint.c -o tkpca_2dcheckint -O3 -lm -lfftw3
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

// Lattice spacing - NOTE only supported in s3dpca code!!!
#undef DX
#undef DY
#undef DZ

#undef LX
#undef LY
#undef LZ
inline double switch_function(double t, double T, double alpha)
{
    return 0.5*( 1.0+tanh( alpha*tan( M_PI_2*( 2.0*t/T-1.0 ) ) ) );
}

double u_ext(int ix, int iy, int it, int spin, int NX, int NY, double DX, double DY)
{
    return 0.0;
    /*
#define xs (0.08*(double)(LX))
#define ys (0.08*(double)(LY))
#define HO(x, omega2p2) (omega2p2*x*x)
    double kF = 0.938673;
    double eF = 0.5*kF*kF;
    double dc_params[32];
    double LX=DX*NX;
    double LY=DY*NY;
    
    dc_params[0]=                 0.75; //    # V(LX/2)/eF
    dc_params[1]=                 0.75; //    # V(LY/2)/eF
    dc_params[5]=                 -1;//  # only up to iter<dc_params5 do imprint, modify_potentials
    dc_params[6]=                 0.7;// # Omega=dc_params6 * omega_x for rotating frame system

    dc_params[0] = 4.0*dc_params[0]*eF/(DX*DX*NX*NX);
    dc_params[1] = 4.0*dc_params[1]*eF/(DY*DY*NY*NY);
    

 
//     Umax=HO(xs) + 1.0*(HO(0)-HO(xs))
    dc_params[2] = HO((double)(-LX/2)+xs, dc_params[0]) + 1.0*(HO((double)(-LX/2), dc_params[0])-HO((double)(-LX/2)+xs, dc_params[0]));
    dc_params[3] = HO((double)(-LY/2)+ys, dc_params[1]) + 1.0*(HO((double)(-LY/2), dc_params[1])-HO((double)(-LY/2)+ys, dc_params[1]));
    
    double _x = ((double)(ix) - 1.0*(NX/2))*DX;
    double _y = ((double)(iy) - 1.0*(NY/2))*DY;
        
    double swtch=1.0;
//     if(1.0*it<dc_params[4]) swtch=switch_function(1.0*it, dc_params[4], 1.0);
    
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
        hox = s*HO(_x, dc_params[0]) + (1.0-s)*dc_params[2];
    }
    else
    {
        s = switch_function(_x-((double)(LX/2)-xs), xs, 1.0);
        hox=(1.0-s)*HO(_x, dc_params[0]) + s*dc_params[2];
    }
    
    double hoy;
    if(_y>=(double)(-LY/2)+ys && _y<=(double)(LY/2)-ys) 
    {
        hoy = HO(_y, dc_params[1]);
    }
    else if(_y<(double)(-LY/2)+ys) 
    {   
        s = switch_function(_y+(double)(LY/2), ys, 1.0);
        hoy = s*HO(_y, dc_params[1]) + (1.0-s)*dc_params[3];
    }
    else
    {
        s = switch_function(_y-((double)(LY/2)-ys), ys, 1.0);
        hoy=(1.0-s)*HO(_y, dc_params[1]) + s*dc_params[3];
    }
        
    return swtch*(hox + hoy);  */  
    
}
// -----------------------------------------------------------------------------------

#define ZERO_NOISE 1.0e-16
#define ZERO_NOISE_CORR(array) \
    for(ixyz=0; ixyz<(outNX*outNY); ixyz++) if(array[ixyz]<ZERO_NOISE) array[ixyz]=ZERO_NOISE;

typedef char * string;
int interpolate_array_R(int inNX, int inNY, double inDX, double inDY, double *in_array, 
                        int outNX, int outNY, double outDX, double outDY, double *out_array,
                        MPI_Comm comm
                       );
int interpolate_array_C(int inNX, int inNY, double inDX, double inDY, double complex *in_array, 
                        int outNX, int outNY, double outDX, double outDY, double complex *out_array,
                        MPI_Comm comm
                       );



int main( int argc , char ** argv ) 
{
    // settings
    
/*    // INPUT
    char inprefix[256] =  "/home/gabrielw/MyProjects/dft/pca/tests/VRTEST/P0.48.i1";
    int inNX=48;
    int inNY=48;
    int inNZ=48;
    double inDX=(64./48.);
    double inDY=(64./48.);
    double inDZ=(64./48.);
    
    // OUTPUT
    char outprefix[256] = "/home/gabrielw/MyProjects/dft/pca/tests/VRTEST/P0.64.t0";
    int outNX=64;
    int outNY=64;
    int outNZ=64;
    double outDX=1.0;
    double outDY=1.0;
    double outDZ=1.0;  */  
    
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
    int inNZ=32;
    double inDX=atof(argv[5])/inNX;
    double inDY=inDX;
    double inDZ=1.0;
    
    // OUTPUT
    char *outprefix = argv[2];
    int outNX=atoi(argv[4]);
    int outNY=outNX;
    int outNZ=32;
    double outDX=atof(argv[5])/outNX;
    double outDY=outDX;
    double outDZ=1.0;  
    
    // arrays for input
    double *h_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *h_potentials; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    cppmallocl(h_densities,12*inNX*inNY,double);
    cppmallocl(h_potentials,4*inNX*inNY,double);
    
    // For easier access to data
    // densities 
    double *rho_a = (double *)(h_densities +  0*inNX*inNY);
    double *rho_b = (double *)(h_densities +  1*inNX*inNY);
    double *tau_a = (double *)(h_densities +  2*inNX*inNY);
    double *tau_b = (double *)(h_densities +  3*inNX*inNY);
    double complex *nu = (double complex *)(h_densities +  4*inNX*inNY);
    double *j_a_x = (double *)(h_densities +  6*inNX*inNY);
    double *j_a_y = (double *)(h_densities +  7*inNX*inNY);
    double *j_a_z = (double *)(h_densities +  8*inNX*inNY);
    double *j_b_x = (double *)(h_densities +  9*inNX*inNY);
    double *j_b_y = (double *)(h_densities + 10*inNX*inNY);
    double *j_b_z = (double *)(h_densities + 11*inNX*inNY);
    
    // pontentials
    double *V_a = (double *)(h_potentials +  0*inNX*inNY);
    double *V_b = (double *)(h_potentials +  1*inNX*inNY);
    double complex *delta = (double complex *)(h_potentials +  2*inNX*inNY);
    
    
    // arrays for output
    double *out_h_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *out_h_potentials; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    cppmallocl(out_h_densities,12*outNX*outNY,double);
    cppmallocl(out_h_potentials,4*outNX*outNY,double);
    
    // For easier access to data
    // densities 
    double *out_rho_a = (double *)(out_h_densities +  0*outNX*outNY);
    double *out_rho_b = (double *)(out_h_densities +  1*outNX*outNY);
    double *out_tau_a = (double *)(out_h_densities +  2*outNX*outNY);
    double *out_tau_b = (double *)(out_h_densities +  3*outNX*outNY);
    double complex *out_nu = (double complex *)(out_h_densities +  4*outNX*outNY);
    double *out_j_a_x = (double *)(out_h_densities +  6*outNX*outNY);
    double *out_j_a_y = (double *)(out_h_densities +  7*outNX*outNY);
    double *out_j_a_z = (double *)(out_h_densities +  8*outNX*outNY);
    double *out_j_b_x = (double *)(out_h_densities +  9*outNX*outNY);
    double *out_j_b_y = (double *)(out_h_densities + 10*outNX*outNY);
    double *out_j_b_z = (double *)(out_h_densities + 11*outNX*outNY);
    
    // pontentials
    double *out_V_a = (double *)(out_h_potentials +  0*outNX*outNY);
    double *out_V_b = (double *)(out_h_potentials +  1*outNX*outNY);
    double complex *out_delta = (double complex *)(out_h_potentials +  2*outNX*outNY);

    if(iam==0)
    {
        sprintf(file_name, "%s_checkpoint.kzpca", inprefix);
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
        fread(h_potentials , sizeof(double)*inNX*inNY, 4 , pFile);
        fread(h_densities  , sizeof(double)*inNX*inNY, 12, pFile);
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
    MPI_Bcast(h_potentials , 4*inNX*inNY , MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    MPI_Bcast(h_densities  , 12*inNX*inNY, MPI_DOUBLE , 0 , MPI_COMM_WORLD );
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
        double *rho = (h_densities +  i*inNX*inNY);
        for(ixyz=0; ixyz<(inNX*inNY); ixyz++) _npart+=rho[ixyz]*inDX*inDY;
        _npart*=inDZ*inNZ;
        if(iam==0) printf("# CHECK IN: NPART[%d]=%f\n", i, _npart);
    }
    
    // ---------------------- interpolating arrays -----------------------------
    if(iam==0) printf("#\n");
    if(iam==0) printf("# INTERPOLATING TO THE NEW TARGET LATTICE...\n");
    
    // V_a
    ixyz=0;
    for(ix=0; ix<inNX; ix++) for(iy=0; iy<inNY; iy++)
    {
        V_a[ixyz]-=u_ext(ix, iy, it, SPINA, inNX, inNY, inDX, inDY);
        ixyz++;
    }

    if(iam==0) printf("# INTERPOLATING: V_a...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , V_a, 
                        outNX, outNY, outDX, outDY, out_V_a,
                        MPI_COMM_WORLD);
    ixyz=0;
    for(ix=0; ix<outNX; ix++) for(iy=0; iy<outNY; iy++)
    {
        out_V_a[ixyz]+=u_ext(ix, iy, it, SPINA, outNX, outNY, outDX, outDY);
        ixyz++;
    }
    
    // V_b
    ixyz=0;
    for(ix=0; ix<inNX; ix++) for(iy=0; iy<inNY; iy++)
    {
        V_b[ixyz]-=u_ext(ix, iy, it, SPINB, inNX, inNY, inDX, inDY);
        ixyz++;
    }
    
    if(iam==0) printf("# INTERPOLATING: V_b...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , V_b, 
                        outNX, outNY, outDX, outDY, out_V_b,
                        MPI_COMM_WORLD);
    ixyz=0;
    for(ix=0; ix<outNX; ix++) for(iy=0; iy<outNY; iy++)
    {
        out_V_b[ixyz]+=u_ext(ix, iy, it, SPINB, outNX, outNY, outDX, outDY);
        ixyz++;
    }
    
    // delta
    if(iam==0) printf("# INTERPOLATING: delta...\n");
    interpolate_array_C(inNX , inNY , inDX , inDY , delta, 
                        outNX, outNY, outDX, outDY, out_delta,
                        MPI_COMM_WORLD);

    // rho_a
    if(iam==0) printf("# INTERPOLATING: rho_a...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , rho_a, 
                        outNX, outNY, outDX, outDY, out_rho_a,
                        MPI_COMM_WORLD);
    ZERO_NOISE_CORR(out_rho_a);
    
    // rho_b
    if(iam==0) printf("# INTERPOLATING: rho_b...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , rho_b, 
                        outNX, outNY, outDX, outDY, out_rho_b,
                        MPI_COMM_WORLD);
    ZERO_NOISE_CORR(out_rho_b);
    
    // tau_a
    if(iam==0) printf("# INTERPOLATING: tau_a...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , tau_a, 
                        outNX, outNY, outDX, outDY, out_tau_a,
                        MPI_COMM_WORLD);
    ZERO_NOISE_CORR(out_tau_a);
    
    // tau_b
    if(iam==0) printf("# INTERPOLATING: tau_b...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , tau_b, 
                        outNX, outNY, outDX, outDY, out_tau_b,
                        MPI_COMM_WORLD);
    ZERO_NOISE_CORR(out_tau_b);

    // nu
    if(iam==0) printf("# INTERPOLATING: nu...\n");
    interpolate_array_C(inNX , inNY , inDX , inDY , nu, 
                        outNX, outNY, outDX, outDY, out_nu,
                        MPI_COMM_WORLD);
    
    // j_a_x
    if(iam==0) printf("# INTERPOLATING: j_a_x...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , j_a_x, 
                        outNX, outNY, outDX, outDY, out_j_a_x,
                        MPI_COMM_WORLD);
    
    // j_a_y
    if(iam==0) printf("# INTERPOLATING: j_a_y...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , j_a_y, 
                        outNX, outNY, outDX, outDY, out_j_a_y,
                        MPI_COMM_WORLD);
    
    // j_a_z
    if(iam==0) printf("# INTERPOLATING: j_a_z...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , j_a_z, 
                        outNX, outNY, outDX, outDY, out_j_a_z,
                        MPI_COMM_WORLD);
    
    // j_b_x
    if(iam==0) printf("# INTERPOLATING: j_b_x...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , j_b_x, 
                        outNX, outNY, outDX, outDY, out_j_b_x,
                        MPI_COMM_WORLD);
    
    // j_b_y
    if(iam==0) printf("# INTERPOLATING: j_b_y...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , j_b_y, 
                        outNX, outNY, outDX, outDY, out_j_b_y,
                        MPI_COMM_WORLD);
    
    // j_b_z
    if(iam==0) printf("# INTERPOLATING: j_b_z...\n");
    interpolate_array_R(inNX , inNY , inDX , inDY , j_b_z, 
                        outNX, outNY, outDX, outDY, out_j_b_z,
                        MPI_COMM_WORLD);
    
    // test - particle number
    for(i=0; i<2; i++)
    {
        double _npart=0.0;
        double *rho = (out_h_densities +  i*outNX*outNY);
        for(ixyz=0; ixyz<(outNX*outNY); ixyz++) _npart+=rho[ixyz]*outDX*outDY;
        _npart*=outDZ*outNZ;
        if(iam==0) printf("# CHECK OUT: NPART[%d]=%f\n", i, _npart);
        npart[i]=_npart;
    }
    
    dc_ec = M_PI*M_PI/(2.*outDX*outDX); 
    if(iam==0) printf("# NEW VALUE FOR dc_ec=%f\n", dc_ec);
    
    // checkpoint - only by iam==0
    if(iam==0)
    {
        sprintf(file_name, "%s_checkpoint.kzpca", outprefix);
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
        fwrite(out_h_potentials , sizeof(double)*outNX*outNY, 4 , pFile);
        fwrite(out_h_densities  , sizeof(double)*outNX*outNY, 12, pFile);
        fwrite(energy       , sizeof(double)      , 5 , pFile);
        fwrite(npart        , sizeof(double)      , 2 , pFile);
                
        fclose(pFile);
    }

    if(iam==0)
    {
        sprintf(file_name, "%s_sectionX.out.txt", outprefix);
        FILE * pFile = fopen(file_name, "w");
        printf("# CREATING FILE: `%s`\n", file_name);
        
        int nx=outNX;
        int ny=outNY;
        int nz=1;
    
        iy=ny/2;
        iz=0;        
        for(ix=0; ix<nx; ix++)
        {
            ixyz=ix*ny*nz + iy*nz + iz;
            
            fprintf(pFile, "%14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g\n",
                outDX*ix, //1
                out_rho_a[ixyz], // 2
                out_rho_b[ixyz], // 3
                cabs(out_delta[ixyz]), // 4
                out_tau_a[ixyz], // 5
                out_tau_b[ixyz], // 6
                cabs(out_nu[ixyz]), // 7
                out_V_a[ixyz], // 8
                out_V_b[ixyz], // 9
                u_ext(ix, iy, it, SPINA, outNX, outNY, outDX, outDY), // 10
                u_ext(ix, iy, it, SPINB, outNX, outNY, outDX, outDY) // 11
            );
        }    
        
        // section for input
        sprintf(file_name, "%s_sectionX.in.txt", outprefix);
        pFile = fopen(file_name, "w");
        printf("# CREATING FILE: `%s`\n", file_name);
        
        nx=inNX;
        ny=inNY;
        nz=1;
    
        iy=ny/2;
        iz=0;        
        for(ix=0; ix<nx; ix++)
        {
            ixyz=ix*ny*nz + iy*nz + iz;
            
            fprintf(pFile, "%14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g\n",
                inDX*ix, //1
                rho_a[ixyz], // 2
                rho_b[ixyz], // 3
                cabs(delta[ixyz]), // 4
                tau_a[ixyz], // 5
                tau_b[ixyz], // 6
                cabs(nu[ixyz]), // 7
                V_a[ixyz], // 8
                V_b[ixyz], // 9
                u_ext(ix, iy, it, SPINA, inNX, inNY, inDX, inDY), // 10
                u_ext(ix, iy, it, SPINB, inNX, inNY, inDX, inDY) // 11
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
double complex interpolate(double x, double y, double complex *fcoeffs, double *kkx, double *kky, int nx, int ny)
{
    int ix, iy, ixyz;
    double complex r = 0.0 + I*0.0;
    ixyz=0;
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++)
    {
        r += fcoeffs[ixyz]*cexp(I*(kkx[ix]*x + kky[iy]*y));
        ixyz++;
    }
    
    return r;
}

int interpolate_array_R(int inNX, int inNY, double inDX, double inDY, double *in_array, 
                        int outNX, int outNY, double outDX, double outDY, double *out_array,
                        MPI_Comm comm
                       )
{
    int iam, np;
    MPI_Comm_size( comm , &np ) ; /* total number of processes */
    MPI_Comm_rank( comm , &iam ) ; /* id of process st 0 <= iam < np */    

    int i,j;
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky; /* values */
    cppmallocl(kkx,inNX,double);
    cppmallocl(kky,inNY,double);
    
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


    // Creating FFTW plans needed for interpolations
    #define USE_FFTW_PLANNER FFTW_ESTIMATE
    double complex *fft2;
    cppmallocl(fft2,inNX*inNY,double complex);
    
    fftw_plan plan_f = fftw_plan_dft_2d(inNX, inNY, fft2, fft2, FFTW_FORWARD, USE_FFTW_PLANNER);

    int ix, iy, ixyz;    
    for(ixyz=0; ixyz<(inNX*inNY); ixyz++) fft2[ixyz] = in_array[ixyz] + I*0.0;
    fftw_execute(plan_f);
    for(ixyz=0; ixyz<(inNX*inNY); ixyz++) fft2[ixyz]/=(inNX*inNY); // introduce normalization factor here
    
    // reset output array;
    for(ixyz=0; ixyz<(outNX*outNY); ixyz++) out_array[ixyz]=0.0;
    
    // interpolate for each point
    ixyz=-1;
    for(ix=0; ix<outNX; ix++) for(iy=0; iy<outNY; iy++)
    {
        ixyz++;
        if(ixyz % np != iam) continue; // not my point
        
        out_array[ixyz] = creal( interpolate(outDX*ix, outDY*iy, fft2, kkx, kky, inNX, inNY) );
    }
    
    // take results from all processes
    MPI_Allreduce( MPI_IN_PLACE, out_array, outNX*outNY, MPI_DOUBLE, MPI_SUM, comm);

    free(fft2);
    free(kkx);
    free(kky);
    fftw_destroy_plan(plan_f);   
    return 0;
}

int interpolate_array_C(int inNX, int inNY, double inDX, double inDY, double complex *in_array, 
                        int outNX, int outNY, double outDX, double outDY, double complex *out_array,
                        MPI_Comm comm
                       )
{
    int iam, np;
    MPI_Comm_size( comm , &np ) ; /* total number of processes */
    MPI_Comm_rank( comm , &iam ) ; /* id of process st 0 <= iam < np */    

    int i,j;
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky; /* values */
    cppmallocl(kkx,inNX,double);
    cppmallocl(kky,inNY,double);
    
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


    // Creating FFTW plans needed for interpolations
    #define USE_FFTW_PLANNER FFTW_ESTIMATE
    double complex *fft2;
    cppmallocl(fft2,inNX*inNY,double complex);
    
    fftw_plan plan_f = fftw_plan_dft_2d(inNX, inNY, fft2, fft2, FFTW_FORWARD, USE_FFTW_PLANNER);

    int ix, iy, ixyz;    
    for(ixyz=0; ixyz<(inNX*inNY); ixyz++) fft2[ixyz] = in_array[ixyz];
    fftw_execute(plan_f);
    for(ixyz=0; ixyz<(inNX*inNY); ixyz++) fft2[ixyz]/=(inNX*inNY); // introduce normalization factor here
    
    // reset output array;
    for(ixyz=0; ixyz<(outNX*outNY); ixyz++) out_array[ixyz]=0.0+I*0.0;
    
    // interpolate for each point
    ixyz=-1;
    for(ix=0; ix<outNX; ix++) for(iy=0; iy<outNY; iy++)
    {
        ixyz++;
        if(ixyz % np != iam) continue; // not my point
        
        out_array[ixyz] = interpolate(outDX*ix, outDY*iy, fft2, kkx, kky, inNX, inNY);
    }
    
    // take results from all processes
    MPI_Allreduce( MPI_IN_PLACE, out_array, outNX*outNY, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);

    free(fft2);
    free(kkx);
    free(kky);
    fftw_destroy_plan(plan_f);   
    return 0;
}
