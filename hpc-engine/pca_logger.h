// Author: Gabriel Wlazlowski
// Date: 08-02-2017

int add_line_to_file(int imeasurement, double rt, int nr_of_items, double *items)
{
    FILE * log;
    char file_name[256];
    sprintf(file_name, "%s_run.log", md.outprefix);
        
    log = fopen (file_name,"a");
    if(log==NULL) // error - cannot create the file
        return 1;
    
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [20];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (buffer,20,"%x-%X",timeinfo);
    
    char line[1024] = "";
    char item[32];
    int i;
    sprintf(item, "%8d", imeasurement);
    strcat(line,item);
    for(i=0; i< nr_of_items; i++) 
    {
        sprintf(item, " %18.10g ", items[i]);
        strcat (line,item);
    }
    sprintf(item, " %8.3f ", rt);
    strcat(line,item);
    sprintf(item, " %20s", buffer);
    strcat(line,item);
    
    fprintf(log, "%s\n", line);
    
    fclose(log);
    
    return 0;   
}
// ------------------------------- CUSTOM LOGGER -------------------------------------
static double __logger_t_gettimeofday ;
static struct timeval __logger_s ;
int logger_create_header(const char *execcmd)
{
    FILE * log;
    char file_name[512];
    int i;
    sprintf(file_name, "%s.wlog", md.outprefix);
    
    if(md.overwrite==0) // break if file exists
    {
        log = fopen (file_name,"r"); // try to open
        if(log!=NULL)
        {
            fclose(log);
            return 1;
        }
    }
    
    // create file
    log = fopen (file_name,"w");
    if(log==NULL) // error - cannot create the file
        return 2;
    
    time_t rawtime;
    time ( &rawtime );
    fprintf(log,"# CREATION TIME OF THE LOG: %s", ctime (&rawtime));
    fprintf(log,"# EXECUTION COMMAND       : %s\n", execcmd);
    fprintf(log,"# CODE NAME               : %s\n",STRINGIZE(CODE));
    fprintf(log,"# VERSION OF THE CODE     : %s\n", VERSION);
    fprintf(log,"# COMPILATION DATE & TIME : %s, %s\n", __DATE__,__TIME__);
    
    fclose(log);
    
    // reset timers
    gettimeofday( &__logger_s , NULL ) ;
    __logger_t_gettimeofday = __logger_s.tv_sec + 1e-6 * __logger_s.tv_usec ;
    
    return 0;
}

double logger_get_time_from_last_entry()
{
    gettimeofday( &__logger_s , NULL ) ;
    double rt = __logger_s.tv_sec + 1e-6 * __logger_s.tv_usec - __logger_t_gettimeofday ;
    __logger_t_gettimeofday = __logger_s.tv_sec + 1e-6 * __logger_s.tv_usec ;
    return rt;
}

// implemented in logger.h
int logger(FILE *log, 
           int it, 
           wslda_density h_densities, wslda_potential h_potentials, 
           double kF, double *mu,
           double *energy, double *npart, 
           double *params, size_t extra_data_size, void *extra_data);


int logger_add_entry(int it, 
           wslda_density h_densities, wslda_potential h_potentials, 
           double kF, double *mu,
           double *energy, double *npart, 
           double *params, size_t extra_data_size, void *extra_data)
{
    FILE * log;
    char file_name[512];
    int i;
    sprintf(file_name, "%s.wlog", md.outprefix);
    
    // append
    log = fopen (file_name,"a");
    if(log==NULL) // error - cannot create the file
        return -1;
    
    int ierr = logger(log, it, h_densities, h_potentials, kF, mu, energy, npart, params, extra_data_size, extra_data);
    
    fclose(log);
    
    return ierr;
}
