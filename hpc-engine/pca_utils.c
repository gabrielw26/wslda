/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * */
#include "pca_settings.h"
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#define ALLOCATE_MD_STRUCTURE
#include "pca_utils.h"

// allocate global metadata structure
metadata_t md = 
{ // create with default values
-1, //inittype;             
1, //measurements;         
1, // timesteps;            
0.01, //dt;                
0.0, // kc;                
0.0, //ec;                
"none", // inprefix
"pca", // outprefix
512, // nthreads;             
1000000, // batch;                
0, // overwrite;            
0, // checkpoint;           
0, // selfstart;
0.0, //qfalpha;           
0.0, // qfbeta;            
0.0, //qfgamma;           
0.0, // qfstart;           
0.0, // qfstop;            
0.0, // qfswitch;          
100.0, // Na;                
100.0, // Nb;         
100.0, // init0Na;                
100.0, // init0Nb; 
1.0e-4, // muchange;   
0.2, // Tstart;
0.05, // Tstop;
0.01, // DeltaT;
1.0e-6, // init0eps;
0.25, // init0scmix
100000, // init0maxiter;
0, // init0debug
0, // init0save
0, // p;                    
0, // q;                    
32, // mb;                   
32, // nb;    
0, // tsubamenodes
1.0e-6, // energyconveps
1.0e-6, // npartconveps
0.5, // linearmixing
0.5, // kzmuchange
10000, // kzmaxiters
0.01, // kztemp
0.0, // referencekF
0, // spinsymmetry
0.1, // mumaxchange
1, // resetit
0, // writewf
1.0e12, // writeecut
0.0, // aBdG
0, // nocurrents
0, // broyden
5, // Mbroyden
0, // startbroyden
999999, // stopbroyden
0.75, // broydenmixing
0.01, //	omega0broyden
1.,	// omeganbroyden
1.,	// omegakbroyden
10000.0, // walltime
-10.0, // ccstart; 
99999.0, // ccstop;
10.0, // ccswitch;
1, // iogroups
};

/** 
 * Function reads input file 
 * and puts values into global struct `input`
 * @return 1 - success, 0 - fail
 * */
int parse_input_file(char * file_name)
{
    FILE *fp;
    fp=fopen(file_name, "r");
    if(fp==NULL)
        return 0;
    
    int i;
    for(i=0; i<MAX_USER_PARAMS; i++) md.params[i]=0.0; // reset parameters
    
    // reset list of variables
    md.nwritevars=0;
    // add default variables
    sprintf(md.writevars[md.nwritevars],"density"); md.nwritevars++;
    sprintf(md.writevars[md.nwritevars],"delta"); md.nwritevars++;
    sprintf(md.writevars[md.nwritevars],"current"); md.nwritevars++;
        
    char s[MAX_REC_LEN];
    char tag[MAX_REC_LEN];
    char ptag[MAX_REC_LEN];
    double tmpparam;
    while(fgets(s, MAX_REC_LEN, fp) != NULL)
    {
        // Read first element of line
        tag[0]='#'; tag[1]='\0';
        sscanf (s,"%s %*s",tag);
        
        // Loop over known tags;
        if(strcmp (tag,"#") == 0)
            continue;
        else if (strcmp (tag,"inittype") == 0)
            sscanf (s,"%s %d %*s",tag,&md.inittype);
        else if (strcmp (tag,"measurements") == 0)
            sscanf (s,"%s %d %*s",tag,&md.measurements);
        else if (strcmp (tag,"timesteps") == 0)
            sscanf (s,"%s %d %*s",tag,&md.timesteps);
        else if (strcmp (tag,"dt") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.dt);
        else if (strcmp (tag,"kc") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.kc);
        else if (strcmp (tag,"ec") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.ec);
        else if (strcmp (tag,"inprefix") == 0)
            sscanf (s,"%s %s %*s",tag,md.inprefix);
        else if (strcmp (tag,"outprefix") == 0)
            sscanf (s,"%s %s %*s",tag,md.outprefix);
        else if (strcmp (tag,"nthreads") == 0)
            sscanf (s,"%s %d %*s",tag,&md.nthreads);
        else if (strcmp (tag,"batch") == 0)
            sscanf (s,"%s %d %*s",tag,&md.batch); 
        else if (strcmp (tag,"overwrite") == 0)
            sscanf (s,"%s %d %*s",tag,&md.overwrite);
        else if (strcmp (tag,"checkpoint") == 0)
            sscanf (s,"%s %d %*s",tag,&md.checkpoint);
        else if (strcmp (tag,"selfstart") == 0)
            sscanf (s,"%s %d %*s",tag,&md.selfstart);        
        // QUANTUM FRICTION
        else if (strcmp (tag,"qfalpha") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.qfalpha);
        else if (strcmp (tag,"qfbeta") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.qfbeta);
        else if (strcmp (tag,"qfgamma") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.qfgamma);
        else if (strcmp (tag,"qfstart") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.qfstart);
        else if (strcmp (tag,"qfstop") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.qfstop);
        else if (strcmp (tag,"qfswitch") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.qfswitch);
        // PARTICLE NUMBER
        else if (strcmp (tag,"Na") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.Na);
        else if (strcmp (tag,"Nb") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.Nb);
        // INIT-0 parameters
        else if (strcmp (tag,"init0Na") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.init0Na);
        else if (strcmp (tag,"init0Nb") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.init0Nb);
        else if (strcmp (tag,"muchange") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.muchange);
        else if (strcmp (tag,"Tstart") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.Tstart);
        else if (strcmp (tag,"Tstop") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.Tstop);
        else if (strcmp (tag,"DeltaT") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.DeltaT);
        else if (strcmp (tag,"init0eps") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.init0eps);
        else if (strcmp (tag,"init0scmix") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.init0scmix);
        else if (strcmp (tag,"init0maxiter") == 0)
            sscanf (s,"%s %d %*s",tag,&md.init0maxiter);
        else if (strcmp (tag,"init0debug") == 0)
            sscanf (s,"%s %d %*s",tag,&md.init0debug);
        else if (strcmp (tag,"init0save") == 0)
            sscanf (s,"%s %d %*s",tag,&md.init0save);
        // SCALAPACK additional parameters
        else if (strcmp (tag,"p") == 0)
            sscanf (s,"%s %d %*s",tag,&md.p);
        else if (strcmp (tag,"q") == 0)
            sscanf (s,"%s %d %*s",tag,&md.q);
        else if (strcmp (tag,"mb") == 0)
            sscanf (s,"%s %d %*s",tag,&md.mb);
        else if (strcmp (tag,"nb") == 0)
            sscanf (s,"%s %d %*s",tag,&md.nb);
        else if (strcmp (tag,"tsubamenodes") == 0)
            sscanf (s,"%s %d %*s",tag,&md.tsubamenodes);
        // kz-solver
        else if (strcmp (tag,"energyconveps") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.energyconveps);
        else if (strcmp (tag,"npartconveps") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.npartconveps);
        else if (strcmp (tag,"linearmixing") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.linearmixing);
        else if (strcmp (tag,"kzmuchange") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.kzmuchange);
        else if (strcmp (tag,"kzmaxiters") == 0)
            sscanf (s,"%s %d %*s",tag,&md.kzmaxiters);
        else if (strcmp (tag,"kztemp") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.kztemp);
        else if (strcmp (tag,"referencekF") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.referencekF);
        else if (strcmp (tag,"spinsymmetry") == 0)
            sscanf (s,"%s %d %*s",tag,&md.spinsymmetry);   
        else if (strcmp (tag,"mumaxchange") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.mumaxchange);
        else if (strcmp (tag,"resetit") == 0)
            sscanf (s,"%s %d %*s",tag,&md.resetit);
        else if (strcmp (tag,"writewf") == 0)
            sscanf (s,"%s %d %*s",tag,&md.writewf);
        else if (strcmp (tag,"writeecut") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.writeecut);
        else if (strcmp (tag,"aBdG") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.aBdG);
        else if (strcmp (tag,"nocurrents") == 0)
            sscanf (s,"%s %d %*s",tag,&md.nocurrents);
        // broyden
        else if (strcmp (tag,"broyden") == 0)
            sscanf (s,"%s %d %*s",tag,&md.broyden);      
        else if (strcmp (tag,"Mbroyden") == 0)
            sscanf (s,"%s %d %*s",tag,&md.Mbroyden); 
        else if (strcmp (tag,"startbroyden") == 0)
            sscanf (s,"%s %d %*s",tag,&md.startbroyden); 
        else if (strcmp (tag,"stopbroyden") == 0)
            sscanf (s,"%s %d %*s",tag,&md.stopbroyden); 
        else if (strcmp (tag,"broydenmixing") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.broydenmixing);  
        else if (strcmp (tag,"omega0broyden") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.omega0broyden);     
        else if (strcmp (tag,"omegakbroyden") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.omegakbroyden);   
        else if (strcmp (tag,"omeganbroyden") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.omeganbroyden);   
        // technical
        else if (strcmp (tag,"walltime") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.walltime);
        // current corrections
        else if (strcmp (tag,"ccstart") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.ccstart);
        else if (strcmp (tag,"ccstop") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.ccstop);
        else if (strcmp (tag,"ccswitch") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.ccswitch);       
        // IO
        else if (strcmp (tag,"iogroups") == 0)
            sscanf (s,"%s %d %*s",tag,&md.iogroups);
        else
        {
            // POTENTIAL PARAMETERS
            for(i=0; i<MAX_USER_PARAMS; i++)
            {
                sprintf(ptag,"params%d",i);
                if (strcmp (tag,ptag) == 0)
                {
                    sscanf (s,"%s %lf %*s",tag,&tmpparam);
                    md.params[i]=(double)tmpparam;
                    break;
                }
            }
            
            // VARIABLES TO WITE 
            // TODO
        }
        
    }
        
    fclose(fp);
    return 1;
}
 

void b_t( void ) 
{ /* hack together a clock w/ microsecond resolution */
  gettimeofday( &s , NULL ) ;
  t_clock = clock() ;
  t_gettimeofday = s.tv_sec + 1e-6 * s.tv_usec ;
}

double e_t( int type ) 
{
  switch ( type ) 
    {
    case 0 :
      t_clock1 = clock() ;
      gettimeofday( &s , NULL ) ;
      t_clock = t_clock1 - t_clock ;
      t_gettimeofday = s.tv_sec + 1e-6 * s.tv_usec - t_gettimeofday ;
      return t_gettimeofday ;
    case 1 :
      return t_gettimeofday ;
    case 2 :
      return t_clock / ( double ) CLOCKS_PER_SEC ;
    }
  return t_gettimeofday ;
}


void wt_b_t( void ) 
{ /* hack together a clock w/ microsecond resolution */
  gettimeofday( &s , NULL ) ;
  wt_t_gettimeofday = s.tv_sec + 1e-6 * s.tv_usec ;
}

double wt_e_t( void ) 
{
    gettimeofday( &s , NULL ) ;
    return s.tv_sec + 1e-6 * s.tv_usec - wt_t_gettimeofday ;
}

void getnwfip( int ip , int np , int nwf , int * nwfip ) 
{
  *nwfip  = nwf / np ;
  if ( ip < (nwf % np) ) (*nwfip)++ ;
}

void print_help(char *progname)
{
  printf("USAGE:\n");
  printf("\t%s input_file_name\n\n",progname);
  printf("OR\n");
  printf("\t%s -v\n",progname);
  printf("\tfor printing information about the version.\n\n");
  
}

void print_version()
{
  printf("Code       : %s\n",STRINGIZE(CODE));
  printf("Version    : %s\n",VERSION);
  printf("Build time : %s, %s\n",__DATE__,__TIME__);
  printf("Lattice    : %d x %d x %d\n", NX, NY, NZ);
  printf("Defined macro-variables:\n");
  printf("\tNXYZ=%d\n", NXYZ);
#ifdef TARGET_MACHINE
  printf("\tTARGET_MACHINE=%s\n",STRINGIZE(TARGET_MACHINE));
#endif
#ifdef VERBOSE
  printf("\tVERBOSE\n");
#endif
#ifdef DEBUG
  printf("\tDEBUG\n");
#endif
#ifdef TIMING
  printf("\tTIMING\n");
#endif
#ifdef EPSILON
  printf("\tEPSILON=%g\n",EPSILON);  
#endif
#ifdef MAX_USER_PARAMS
    printf("\tMAX_USER_PARAMS=%d\n", MAX_USER_PARAMS);
#endif
    
    printf("\tUD_SCITERS=%d\n", UD_SCITERS);
    printf("\tUD_MIX_COEFF=%f\n", UD_MIX_COEFF);
    printf("\tDENSEPSILON=%g\n", DENSEPSILON);
    
#ifdef CURRENT_CORRECTIONS
    printf("\tCURRENT_CORRECTIONS\n");
#endif
    
}

/**
 * Function preprocess comand line
 * and returns index of file_name in array argv
 * If input file name not given then returns -1
 * */
int readcmd(int argc, char *argv[])
{
    static const char *optString = "vh";
    int opt = 0;
    do
    {
        opt = getopt( argc, argv, optString );
        switch( opt )
        {
            case 'v':
            print_version();
            break;
                        
            case 'h':
            print_help(argv[0]);
            break;
                        
            default:
            break;
        }
    }
    while(opt != -1);
        
    if (optind < argc)
        return optind;
    else
        return -1;
}

/**
 * Switch function - performs switch in time interval [0-T]
 * */
double h_switch_function(double t, double T, double alpha)
{
    return 0.5*( 1.0+tanh( alpha*tan( M_PI_2*( 2.0*t/T-1.0 ) ) ) );
}

double h_smooth_step(double t, double step_start, double step_stop, double T, double alpha)
{
    if(t<=step_start || t>=step_stop) return 0.0; 
    if(t>=step_start+T && t<=step_stop-T) return 1.0;
    if(t>step_start && t<step_start+T) return h_switch_function(t-step_start, T, alpha);
    else return 1.0-h_switch_function(t-step_stop+T, T, alpha);
}

void symmetrize_densities(double *h_densities)
{
    double complex *nu = (double complex *)(h_densities +  0*NXYZ);
    double *rho_a = (double *)(h_densities +  2*NXYZ);
    double *tau_a = (double *)(h_densities +  3*NXYZ);
    double *j_a_x = (double *)(h_densities +  4*NXYZ);
    double *j_a_y = (double *)(h_densities +  5*NXYZ);
    double *j_a_z = (double *)(h_densities +  6*NXYZ);    
    double *rho_b = (double *)(h_densities +  7*NXYZ);
    double *tau_b = (double *)(h_densities +  8*NXYZ);
    double *j_b_x = (double *)(h_densities +  9*NXYZ);
    double *j_b_y = (double *)(h_densities + 10*NXYZ);
    double *j_b_z = (double *)(h_densities + 11*NXYZ);
    
    int ixyz;
    for(ixyz=0; ixyz<  NXYZ; ixyz++) rho_a[ixyz]=rho_b[ixyz];
    for(ixyz=0; ixyz<  NXYZ; ixyz++) tau_a[ixyz]=tau_b[ixyz];
    for(ixyz=0; ixyz<3*NXYZ; ixyz++) j_a_x[ixyz]=j_b_x[ixyz];

}

int copy_input_file(char * input_file, char * file_name)
{
    FILE * log;
    
    // open file
    log = fopen (file_name,"w");
    if(log==NULL) // error - cannot create the file
        return 1; 
    
    FILE * inp;
    char s[MAX_REC_LEN];
    inp = fopen (input_file,"r");
    if(inp==NULL) return 2; // error - cannot open file
    while(fgets(s, MAX_REC_LEN, inp) != NULL) fprintf(log,"%s",s);
    fclose(inp);
    fclose(log);
    
    return 0;
}
