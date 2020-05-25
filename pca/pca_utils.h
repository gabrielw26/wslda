// Global buffer with all info from input file

#ifndef __PCA_UTILS__
#define __PCA_UTILS__


#define MD_CHAR_LGTH 256
typedef struct
{
    int inittype;                       // 0 - uniform (set INIT-0 parameters)
                                        // 1 - uniform solution but read it from file `inprefix`_uniform.solution
                                        // 2 - start from checkpoint
                                        // 22 - start from interpolated checkpoint, supported by s2dpca and s3dpca
                                        // 3 - start from s2dpca solver solution
                                        // 4 - start from dpca files (works in case of s3dpca) - deprecated!
                                        // 5 - start from s3dpca solver
    int measurements;                   // number of measurements
    int timesteps;                      // number of time steps between mesurements
    double dt;                          // time step
    double kc;                          // momentum cut-off 
    double ec;                          // energy cut-off 
    char inprefix[MD_CHAR_LGTH];        // prefix for input files
    char outprefix[MD_CHAR_LGTH];       // prefix for output files
    int nthreads;                       // number of gpu threads
    int batch;                          // batch size for cuFFT
    int overwrite;                      // eneable overwrite mode?
    int checkpoint;                     // do checkpoint 
    int selfstart;                      // 0 - assume that algorithm starts from eigenstate state (default)
                                        // 1 - for first steps use Taylor expansion of the evolution operator
    
    // QUANTUM FRICTION
    double qfalpha;                     // alpha parameter for quantum friction term
    double qfbeta;                      // beta parameter for quantum friction - parring channel
    double qfgamma;                     // gamma parameter for particle number control
    double qfstart;                     // start time for evolving with quantum friction, in units of eF
    double qfstop;                      // stop time for evolving with quantum friction, in units of eF
    double qfswitch;                    // time for switch function
    
    // PARTICLE NUMBER
    double Na;                 // Requested number of particles a-type
    double Nb;                 // Requested number of particles b-type
    
    // INIT-0 parameters
    double init0Na;            // Requested number of particles a-type - uniform solution
    double init0Nb;            // Requested number of particles b-type - uniform solution
    double muchange;           // Change rate of chemical potential
    double Tstart;             // Start temperature, in units of eF, default 0.2
    double Tstop;              // Stop temperature, in units of eF, default 0.05
    double DeltaT;             // change of temperature in units of eF, default 0.01
    double init0eps;           // epsilon for convergence, default 1.0e-6
    double init0scmix;         // mixing parameter in self-consitent process, default 0.25
    int init0maxiter;          // maximum number of iterations, default 100000
    int init0debug;            // debug level, default=0 (no debug info), 
    int init0save;             // save solution to file?, default=0, if 1 then solution is in file 'outprefix'_uniform.solution
    
    // SCLAPACK additional parameters
    int p;                              // CBLACS grid, defalt p=0 (atomatic assignment of the value)
    int q;                              // CBLACS grid, defalt q=0 (atomatic assignment of the value)
    int mb;                             // CBLACS grid, default mb=32
    int nb;                             // CBLACS grid, default nb=32
    
    // TSUBAME parameter
    int tsubamenodes; // number of nodes used in calculations on tsubame computer
    
    // kz-solver parameters
    double energyconveps; // convergence epsilon for energy- fraction of Effg needed to get convergence, default=1.0e-6
    double npartconveps; // convergence epsilon for particle number- fraction of N_tot=(Na+Nb) needed to get convergence, default=1.0e-6
    double linearmixing; // mixing parameter for linear algorithm, default=0.5
    double kzmuchange; // coefficient for changing chemical potential, default=0.5
    int kzmaxiters; // maximum number of iterations, default=10000
    double kztemp; // temperature in units of eF_a, default=0.01
    double referencekF; // value of reference kF used in calculations, if 0.0 then not set (default)
    int spinsymmetry; // impose spin symmetry, default 0 - no
    double mumaxchange; // maximal change of chemical potential per iteration, in units of Fermi energy, default 0.1
    int resetit; // if 1 set it=0, otherwise continue from value read from checkpoint file, default resetit=1
    int writewf; // if 1 the code will write wave-functions on exit, default writewf=0
    double writeecut; // only states with |E_n/eF|<writeecut will be written, default writeecut=INFINITY
    double aBdG; // scattering length for BdG mode, if aBdG=0.0 then ASLDA is activated, default aBdG=0.0
    int nocurrents; // if 1 then code imposes by hand no currents, default: nocurrents=0
    
    // broyden mixing parameters
    int broyden; // 0 - linear mixing, 1 - update densities with Broyden, default=0
    int Mbroyden; // number of previous iterations taken into account, default=5
    int startbroyden; // firts iteration for broyden activavtion, default=0
    int stopbroyden; // last iteration for broyden activavtion, default=999999
    double broydenmixing; // mixing parameter for Broyden algorithm, default=0.75
	double omega0broyden;	// weight assigned to the error in the inverse Jacobian, default=0.01
	double omeganbroyden;	  // weight associated with each previous iteration, default=1.0
	double omegakbroyden;	   // weight associated with each previous iteration, default=1.0
    
    // walltime
    double walltime; // after this time in hours the energency checkpoint will be executed, default=1000 
    
    // Current corrections
    double ccstart;                     // start time for turning on current corrections, in units of eF
    double ccstop;                      // stop time for the current corrections, in units of eF
    double ccswitch;                    // time for switch function
    
    // IO
    int iogroups;                       // number of IO groups used for wf writing, default=1

    // POTENTIAL PARAMETERS
    double params[MAX_USER_PARAMS];
} metadata_t;

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

#define MAX_REC_LEN 1024
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
        }
    }
        
    fclose(fp);
    return 1;
}
 
#include <sys/time.h>
#include <time.h> /* for ctime() */
static double t_gettimeofday ;
static clock_t t_clock , t_clock1 ;
static struct timeval s ;

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

/*
 * For walltime measuremnt 
 * */
static double wt_t_gettimeofday ;

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

#define STRINGIZE2(s) #s
#define STRINGIZE(s) STRINGIZE2(s)

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
#include <getopt.h>
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
#endif
