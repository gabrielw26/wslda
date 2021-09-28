/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 31-07-2020
 * */

// Global buffer with all info from input file

#ifndef __PCA_UTILS__
#define __PCA_UTILS__


#define MD_CHAR_LGTH 256
#define MAX_WRITEVARS 16
#define MAX_VARNAME_LGTH 32
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
    double init0muchange;           // Change rate of chemical potential
    double init0Tstart;             // Start temperature, in units of eF, default 0.2
    double init0Tstop;              // Stop temperature, in units of eF, default 0.05
    double init0DeltaT;             // change of temperature in units of eF, default 0.01
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

    // GPUS distribution
    int gpuspernode; // number of gpus per node, defualt=1

    // static solver parameters
    double energyconveps; // convergence epsilon for energy- fraction of Effg needed to get convergence, default=1.0e-6
    double npartconveps; // convergence epsilon for particle number- fraction of N_tot=(Na+Nb) needed to get convergence, default=1.0e-6
    double linearmixing; // mixing parameter for linear algorithm, default=0.5
    double muchange; // coefficient for changing chemical potential, default=0.5
    int maxiters; // maximum number of iterations, default=10000
    double temperature; // temperature in units of eF_a, default=0.01
    double referencekF; // value of reference kF used in calculations, if 0.0 then not set (default)
    int spinsymmetry; // impose spin symmetry, default 0 - no
    double mumaxchange; // maximal change of chemical potential per iteration, in units of Fermi energy, default 0.1
    int resetit; // if 1 set it=0, otherwise continue from value read from checkpoint file, default resetit=1
    int writewf; // if 1 the code will write wave-functions on exit, default writewf=0
    double writeecut; // only states with |E_n/eF|<writeecut will be written, default writeecut=INFINITY
    double aBdG; // scattering length for BdG mode, if aBdG=0.0 then ASLDA is activated, default aBdG=0.0
    int nocurrents; // if 1 then code imposes by hand no currents, default: nocurrents=0 - DEPRECATED
    int killcurrents; // if 1 then code imposes by hand no currents, default: nocurrents=0, the same as no currents;
    int nomixstart; // if 1 then in the first iteration do not do mixing, default nomixstart=0
    char mixingtype; // 'd' - mix densities, 'p' - mix potentials (default)

    // broyden mixing parameters
    int broyden; // 0 - linear mixing, 1 - update densities with Broyden, default=0
    int Mbroyden; // number of previous iterations taken into account, default=5
    int startbroyden; // firts iteration for broyden activavtion, default=0
    int stopbroyden; // last iteration for broyden activavtion, default=999999
    double broydenmixing; // mixing parameter for Broyden algorithm, default=0.75
	double omega0broyden;	// weight assigned to the error in the inverse Jacobian, default=0.01
	double omeganbroyden;	  // weight associated with each previous iteration, default=1.0
	double omegakbroyden;	   // weight associated with each previous iteration, default=1.0
	int broydenautores; // automatic restarts of Broyden algorithm: 0-no, 1-yes (default)
	double broydenEmaxchg; // if the total energy between iteration change be more than broydenEmaxchg then Broyden is restarted, (default=0.1)
	int broydenEdelay; // scan energy changes only after broydenEdelay with Broyden has been executed, typically, just after starting the Broyden energy fluctuations are observed which should decay within a few iterations (default=5)

    // walltime
    double walltime; // after this time in hours the energency checkpoint will be executed, default=1000

    // Current corrections
    double ccstart;                     // start time for turning on current corrections, in units of eF
    double ccstop;                      // stop time for the current corrections, in units of eF
    double ccswitch;                    // time for switch function

    // subset tracking
    double subsetMinEn;     // if subsetMinEn!=subsetMaxEn then td code track densities arising from states
    double subsetMaxEn;     // where En in [subsetMinEn,subsetMaxEn], default subsetMinEn=subsetMaxEn=0
    int subsetShiftDmu;     // if 1 then apply extra shift of quasiparticle energies by (mu_a-mu_b)/2, default=0

    // SLDAE
    double aSLDAe;  // s-wave scattering length, default = -1.0
    int pccrSLDAe;  // pairing coupling constant renormalization scheme
                    // 0: in-medium (default)
                    // 1: in-vacuum (Bulgac et al.)

    // IO
    int iogroups;                       // number of IO groups used for wf writing, default=1
    char dataformat[8];                 // format of produced files: wdat or npy, default=wdat
    char initialized; // technical variable, indicating that structure is initialized by the input file
    char stdoutfile[MD_CHAR_LGTH]; // technical variable,

    // POTENTIAL PARAMETERS
    double params[MAX_USER_PARAMS]; // double parameters
    char strings[MAX_USER_PARAMS][MD_CHAR_LGTH]; // strings

    // variables to write
    int nwritevar; // number of variables to write
    char writevar[MAX_WRITEVARS][MAX_VARNAME_LGTH]; // and their names
} metadata_t;

#ifndef ALLOCATE_MD_STRUCTURE
extern metadata_t md;
extern metadata_t *input;
#endif

#define MAX_REC_LEN 1024
/**
 * Function reads input file
 * and puts values into global struct `input`
 * @return 1 - success, 0 - fail
 * */
int parse_input_file(char * file_name);


#include <sys/time.h>
#include <time.h> /* for ctime() */
static double t_gettimeofday ;
static clock_t t_clock , t_clock1 ;
static struct timeval s ;

void b_t( void );

double e_t( int type );

/*
 * For walltime measuremnt
 * */
static double wt_t_gettimeofday ;

void wt_b_t( void );

double wt_e_t( void );

void getnwfip( int ip , int np , int nwf , int * nwfip );

void print_help(char *progname);

#define STRINGIZE2(s) #s
#define STRINGIZE(s) STRINGIZE2(s)

void print_version();

/**
 * Function preprocess comand line
 * and returns index of file_name in array argv
 * If input file name not given then returns -1
 * */
#include <getopt.h>
int readcmd(int argc, char *argv[]);

/**
 * Switch function - performs switch in time interval [0-T]
 * */
double h_switch_function(double t, double T, double alpha);

double h_smooth_step(double t, double step_start, double step_stop, double T, double alpha);

void symmetrize_densities(double *h_densities);

int copy_input_file(char * input_file, char * file_name);

int wslda_check_settings();

int wslda_check_array_against_naninf(int n, double *array);

void wprintf( const char * format, ... );
void wfprintf(FILE *stream,  const char * format, ... );

void testsuite_ok();

void create_reprowf_tar(size_t extra_data_size);
void copy_checkpoint();
void copy_initcheckpoint();
void copy_reprowftar();

void save_extradata_to_file(size_t size, void *extra_data);

#endif
