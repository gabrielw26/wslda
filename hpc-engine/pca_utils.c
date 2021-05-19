/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * */
#include "pca_settings.h"
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#define ALLOCATE_MD_STRUCTURE
#include "pca_utils.h"
#include "wslda_errors.h"

// allocate global metadata structure
metadata_t md = 
{ // create with default values
-1, //inittype;             
1, //measurements;         
1, // timesteps;            
0.01, //dt;                
0.999999*M_PI/DX, // kc;                
M_PI*M_PI/(2.*DX*DX), //ec;                
"none", // inprefix
"wslda", // outprefix
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
-1.0, // init0Na;                
-1.0, // init0Nb; 
-1.0, // init0muchange;   
-1.0, // init0Tstart;
-1.0, // init0Tstop;
0.01, // init0DeltaT;
1.0e-9, // init0eps;
-1.0, // init0scmix
-1, // init0maxiter;
0, // init0debug
0, // init0save
0, // p;                    
0, // q;                    
32, // mb;                   
32, // nb;    
1, // gpuspernode
1.0e-6, // energyconveps
1.0e-6, // npartconveps
0.5, // linearmixing
0.5, // muchange
10000, // maxiters
1.0e-9, // temperature
0.0, // referencekF
0, // spinsymmetry
0.1, // mumaxchange
1, // resetit
0, // writewf
1.0e12, // writeecut
0.0, // aBdG
0, // nocurrents
0, // nomixstart
'p', // mixingtype
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
"wdat", // dataformat
0, // initialized
"wslda.stdout", // stdoutfile
};

metadata_t *input = &md; // additional handler;

// Taken from:
// https://stackoverflow.com/questions/779875/what-function-is-to-replace-a-substring-from-a-string-in-c
void replace_str(char *str,char *org,char *rep)
{
    char *ToRep = strstr(str,org);
    char *Rest = (char*)malloc(strlen(ToRep));
    strcpy(Rest,((ToRep)+strlen(org)));

    strcpy(ToRep,rep);
    strcat(ToRep,Rest);

    free(Rest);
}

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
    md.nwritevar=0;
        
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
        else if (strcmp (tag,"init0muchange") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.init0muchange);
        else if (strcmp (tag,"init0Tstart") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.init0Tstart);
        else if (strcmp (tag,"init0Tstop") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.init0Tstop);
        else if (strcmp (tag,"init0DeltaT") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.init0DeltaT);
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
        else if (strcmp (tag,"gpuspernode") == 0)
            sscanf (s,"%s %d %*s",tag,&md.gpuspernode);
        // kz-solver
        else if (strcmp (tag,"energyconveps") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.energyconveps);
        else if (strcmp (tag,"npartconveps") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.npartconveps);
        else if (strcmp (tag,"linearmixing") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.linearmixing);
        else if (strcmp (tag,"muchange") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.muchange);
        else if (strcmp (tag,"maxiters") == 0)
            sscanf (s,"%s %d %*s",tag,&md.maxiters);
        else if (strcmp (tag,"temperature") == 0)
            sscanf (s,"%s %lf %*s",tag,&md.temperature);
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
        else if (strcmp (tag,"nomixstart") == 0)
            sscanf (s,"%s %d %*s",tag,&md.nomixstart);
        else if (strcmp (tag,"mixingtype") == 0)
            sscanf (s,"%s %c %*s",tag,&md.mixingtype);
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
        else if (strcmp (tag,"dataformat") == 0)
            sscanf (s,"%s %s %*s",tag,md.dataformat);
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
            if (strcmp (tag,"writevar") == 0)
            {
                int ivars, ierr;
                for(ivars=0; ivars<100; ivars++)
                {
//                     wprintf("[PARSER-i]: `%s`, `%s` `%s`\n", s, tag, ptag);
                    ierr = sscanf (s,"%s %s %*s",tag,ptag);
                    if(ierr<2) break;
                    
//                     if(ptag[0]=='#') wprintf("[PARSER-#]:\n");
                    if(ptag[0]=='#') break;
                    
                    if (strcmp (ptag,"all") == 0) 
                    {
#ifdef WSLDA
                        replace_str(s,ptag,"rho delta j nu tau V V_ext delta_ext velocity_ext alpha A");
#else
                        replace_str(s,ptag,"rho delta j nu tau V V_ext delta_ext velocity_ext");
#endif
//                         wprintf("[PARSER-R]: `%s`, `%s` `%s`\n", s, tag, ptag);
                        continue;
                    }
                    
                    if (strcmp (ptag,"default") == 0) 
                    {
                        replace_str(s,ptag,"rho delta j");
//                         wprintf("[PARSER-R]: `%s`, `%s` `%s`\n", s, tag, ptag);
                        continue;
                    }
                        
                    replace_str(s,ptag," ");
                    
                    // check if variable already added
                    ierr=0;
                    for(i=0; i<md.nwritevar; i++) if (strcmp (ptag,md.writevar[i]) == 0) {ierr=1; break;}
                    if(ierr==1) continue; // variable already added;
                    
                    // add variable
                    strcpy(md.writevar[md.nwritevar],ptag); md.nwritevar++;
//                     wprintf("[ADDING]: %d->%s\n", md.nwritevar-1, md.writevar[md.nwritevar-1]);
//                     wprintf("[PARSER-O]: `%s`, `%s` `%s`\n", s, tag, ptag);
                }
            }
        }
        
    }
    
    // prepare for wprintf()
    sprintf(md.stdoutfile, "%s.stdout", md.outprefix);
    FILE * f = fopen(md.stdoutfile, "w"); // clear file
    fclose(f);    
    md.initialized=1;
    
    // add default variables - if not added 
    if(md.nwritevar==0)
    {
        sprintf(md.writevar[md.nwritevar],"rho"); md.nwritevar++;
        sprintf(md.writevar[md.nwritevar],"delta"); md.nwritevar++;
        sprintf(md.writevar[md.nwritevar],"j"); md.nwritevar++;
    }
    
    // defult values
    if(md.temperature<1.0e-9) md.temperature=1.0e-9; // to avoid division by zero when computing beta=1/T
    if(md.init0Na<0.0) md.init0Na=md.Na;
    if(md.init0Nb<0.0) md.init0Nb=md.Nb;
    if(md.init0muchange<0.0) md.init0muchange=md.muchange;
    if(md.init0Tstop<0.0) md.init0Tstop = md.temperature; 
    if(md.init0Tstart<0.0) md.init0Tstart=md.init0Tstop;
    if(md.init0scmix<0.0) md.init0scmix=md.linearmixing;
    if(md.init0maxiter<0) md.init0maxiter=md.maxiters;
    
#ifdef WSLDA
    if(md.muchange>0 && md.Na!=md.Nb && md.spinsymmetry==1) 
    {
        warn_head(stdout);
        wfprintf(stdout, "#\tForcing spinsymmetry=0 since Na!=Nb and muchange>0 (fixed particle number mode). \n");
        warn_foot(stdout);
        md.spinsymmetry=0;
    }
#endif

#ifdef USE_CUBIC_CUTOFF
    wfprintf(stdout, "# CUBIC CUTOFF: RASING ec TO INFINITY!\n");
    md.ec=1.0e16;
    md.kc=1.0e16;
#endif
    
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
  wprintf("USAGE:\n");
  wprintf("\t%s input_file_name\n\n",progname);
  wprintf("OR\n");
  wprintf("\t%s -v\n",progname);
  wprintf("\tfor printing information about the version.\n\n");
  
}

void print_version()
{
  wprintf("Code       : %s\n",STRINGIZE(CODE));
  wprintf("Version    : %s\n",VERSION);
  wprintf("Build time : %s, %s\n",__DATE__,__TIME__);
  wprintf("Lattice    : %d x %d x %d\n", NX, NY, NZ);
  wprintf("Defined macro-variables:\n");
  wprintf("\tNXYZ=%d\n", NXYZ);
#ifdef TARGET_MACHINE
  wprintf("\tTARGET_MACHINE=%s\n",STRINGIZE(TARGET_MACHINE));
#endif
#ifdef VERBOSE
  wprintf("\tVERBOSE\n");
#endif
#ifdef DEBUG
  wprintf("\tDEBUG\n");
#endif
#ifdef TIMING
  wprintf("\tTIMING\n");
#endif
#ifdef EPSILON
  wprintf("\tEPSILON=%g\n",EPSILON);  
#endif
#ifdef MAX_USER_PARAMS
    wprintf("\tMAX_USER_PARAMS=%d\n", MAX_USER_PARAMS);
#endif
    
    wprintf("\tUD_SCITERS=%d\n", UD_SCITERS);
    wprintf("\tUD_MIX_COEFF=%f\n", UD_MIX_COEFF);
    wprintf("\tDENSEPSILON=%g\n", DENSEPSILON);
    
#ifdef CURRENT_CORRECTIONS
    wprintf("\tCURRENT_CORRECTIONS\n");
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

int wslda_check_settings()
{
#if FUNCTIONAL==BDG
    if(md.aBdG==0.0) return WSLDA_ERR_ABDG_NOT_SET;
#endif
    return 0;
}

/**
 * Function checks array againts NaN and Inf.
 * @return 0: WSLDA_OK, WSLDA_ERR_NAN_DETECTED, WSLDA_INF_NAN_DETECTED
 * */
int wslda_check_array_against_naninf(int n, double *array)
{
    int i;
    for(i=0; i<n; i++) 
    {
        if(isnan(array[i])) return WSLDA_ERR_NAN_DETECTED;
        if(isinf(array[i])) return WSLDA_ERR_INF_DETECTED;
    }
    return WSLDA_OK;
}

void wprintf( const char * format, ... )
{
  va_list args;
  va_start (args, format);
  vprintf (format, args); 
  va_end (args);
  if(md.initialized)
  {
      FILE * f = fopen(md.stdoutfile, "a");
      if(f==NULL) { printf("PROBLEM!\n"); fflush(stdout);}
      va_start (args, format);
      vfprintf (f, format, args);
      va_end (args);
      fclose(f);
  }
  fflush(stdout);
}

void wfprintf(FILE *stream,  const char * format, ... )
{
  va_list args;
  va_start (args, format);
  vfprintf (stream, format, args);
  va_end (args);
  if(md.initialized)
  {
      FILE * f = fopen(md.stdoutfile, "a");
      va_start (args, format);
      vfprintf (f, format, args);
      va_end (args);
      fclose(f);
  }
  fflush(stream); 
}


void testsuite_ok()
{
    char file_name[512];
    sprintf(file_name, "%s_testsuite.ok", md.outprefix);
    FILE * f = fopen(file_name, "w");
    fprintf(f,"%s\n",file_name);
    fclose(f);
}

void create_reprowf_tar()
{
    char cmd[2048];
    sprintf(cmd, 
        "tar -cf %s/reprowf.tar %s_predefines.h %s_problem-definition.h %s_logger.h %s/checkpoint.dat %s_input.txt %s.wlog %s.stdout",
        md.outprefix, md.outprefix, md.outprefix, md.outprefix, md.outprefix, md.outprefix, md.outprefix, md.outprefix);
    wprintf("# SYSTEM: %s\n", cmd);
    system(cmd);
}

void copy_checkpoint()
{
    char cmd[2048];
    sprintf(cmd, "cp -f %s_checkpoint.dat %s/checkpoint.dat", md.outprefix, md.outprefix);
    wprintf("# SYSTEM: %s\n", cmd);
    system(cmd);
}

void copy_initcheckpoint()
{
    char cmd[2048];
    sprintf(cmd, "cp -f %s_checkpoint.dat %s_checkpoint.dat.init", md.inprefix, md.outprefix);
    wprintf("# SYSTEM: %s\n", cmd);
    system(cmd);
}

void copy_reprowftar()
{
    char cmd[2048];
    sprintf(cmd, "cp -f %s/reprowf.tar %s_reprowf.tar", md.inprefix, md.outprefix);
    wprintf("# SYSTEM: %s\n", cmd);
    system(cmd);
}
    
