/**
 * W-SLDA Toolkit
 * */

// #define WDATA_TESTING_MODE

// gcc -std=gnu99 wdata.c -o wdata.exe -lm -DWDATA_TESTING_MODE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include "wdata.h"

#define MAX_REC_LEN 1024

/** 
 * Function reads metadata
 * and puts values into struct wdata_metadata
 * @return 0: ok, 1: Cannot open metadata file,
 * */
int wdata_parse_metadata_file(const char * file_name, wdata_metadata *md)
{
    FILE *fp;
    fp=fopen(file_name, "r");
    if(fp==NULL)
        return 1; // Cannot open file
        
    // reset vars
    md->nvars=0;
    md->nlinks=0;
    md->nconsts=0;
    
    // buffers
    char s[MAX_REC_LEN];
    char tag[MAX_REC_LEN];
    char ptag[MAX_REC_LEN];
    
    while(fgets(s, MAX_REC_LEN, fp) != NULL)
    {
        // Read first element of line
        tag[0]='#'; tag[1]='\0';
        sscanf (s,"%s %*s",tag);
        
        // Loop over known tags;
        if(strcmp (tag,"#") == 0)
            continue;
        else if (strcmp (tag,"NX") == 0)
            sscanf (s,"%s %d %*s",tag,&md->NX);
        else if (strcmp (tag,"NY") == 0)
            sscanf (s,"%s %d %*s",tag,&md->NY);
        else if (strcmp (tag,"NZ") == 0)
            sscanf (s,"%s %d %*s",tag,&md->NZ);
        else if (strcmp (tag,"DX") == 0)
            sscanf (s,"%s %lf %*s",tag,&md->DX);
        else if (strcmp (tag,"DY") == 0)
            sscanf (s,"%s %lf %*s",tag,&md->DY);
        else if (strcmp (tag,"DZ") == 0)
            sscanf (s,"%s %lf %*s",tag,&md->DZ);
        else if (strcmp (tag,"datadim") == 0)
            sscanf (s,"%s %d %*s",tag,&md->datadim);
        else if (strcmp (tag,"prefix") == 0)
            sscanf (s,"%s %s %*s",tag,md->prefix);
        else if (strcmp (tag,"cycles") == 0)
            sscanf (s,"%s %d %*s",tag,&md->cycles);
        else if (strcmp (tag,"dt") == 0)
            sscanf (s,"%s %lf %*s",tag,&md->dt);
        else if (strcmp (tag,"t0") == 0)
            sscanf (s,"%s %lf %*s",tag,&md->t0);
        
        // variables
        else if (strcmp (tag,"var") == 0)
        {
            sscanf (s,"%s %s %s %s %lf %*s",tag, &md->vars[md->nvars].name, &md->vars[md->nvars].type, &md->vars[md->nvars].unit); 
            md->nvars++;       
        }
        
        // links
        else if (strcmp (tag,"link") == 0)
        {
            sscanf (s,"%s %s %s %*s",tag, &md->links[md->nlinks].name, &md->links[md->nlinks].linkto); 
            md->nlinks++;       
        }
        
        // consts
        else if (strcmp (tag,"const") == 0)
        {
            sscanf (s,"%s %s %lf %*s",tag, &md->consts[md->nconsts].name, &md->consts[md->nconsts].value); 
            md->nconsts++;       
        }
    }
    
    fclose(fp);
    return 0;
}

void wdata_print_metadata(wdata_metadata *md, FILE *out)
{
    int i;
    
    fprintf(out, "NX %24d   # lattice\n", md->NX);
    fprintf(out, "NY %24d   # lattice\n", md->NY);
    fprintf(out, "NZ %24d   # lattice\n", md->NZ);
    fprintf(out, "DX %24g   # spacing\n", md->DX);
    fprintf(out, "DY %24g   # spacing\n", md->DY);
    fprintf(out, "DZ %24g   # spacing\n", md->DZ);
    fprintf(out, "datadim %19d   # dimension of block size: 1=NX, 2=NX*NY, 3=NX*NY*NZ\n", md->datadim);
    fprintf(out, "prefix %20s   # prefix for files belonging to this data set, binary files have names prefix_variable.wdat\n", md->prefix);
    fprintf(out, "cycles %20d   # number of cycles (measurements)\n", md->cycles);
    fprintf(out, "t0 %24g   # time value for the first cycle\n", md->t0);
    fprintf(out, "dt %24g   # time interval between cycles\n", md->dt);

    // variables
    fprintf(out,"\n");
    fprintf(out,"# variables\n");
    fprintf(out,"# tag                  name                    type                    unit\n");
    for(i=0; i<md->nvars; i++) wdata_print_variable(&md->vars[i], out);
    
    // links
    fprintf(out,"\n");
    fprintf(out,"# links\n");
    fprintf(out,"# tag                  name                 link-to\n");
    for(i=0; i<md->nlinks; i++) wdata_print_link(&md->links[i], out);
    
    // links
    fprintf(out,"\n");
    fprintf(out,"# consts\n");
    fprintf(out,"# tag                  name                   value\n");
    for(i=0; i<md->nconsts; i++) wdata_print_const(&md->consts[i], out);
    
    fprintf(out,"\n");
}

void wdata_print_variable(wdata_variable *md, FILE *out)
{
    fprintf(out, "var%24s%24s%24s\n", md->name, md->type, md->unit);
}

void wdata_print_link(wdata_link *md, FILE *out)
{
    fprintf(out, "link%23s%24s\n", md->name, md->linkto);
}

void wdata_print_const(wdata_const *md, FILE *out)
{
    fprintf(out, "const%22s%24g\n", md->name, md->value);
}

void wdata_add_variable(wdata_metadata *md, wdata_variable *var)
{
    md->vars[md->nvars] = *var;
    md->nvars++;
}

void wdata_add_link(wdata_metadata *md, wdata_link *link)
{
    md->links[md->nlinks] = *link;
    md->nlinks++;
}


void wdata_add_const(wdata_metadata *md, wdata_const *_const)
{
    md->consts[md->nconsts] = *_const;
    md->nconsts++;
}

int wdata_get_blocksize(wdata_metadata *md)
{
    if(md->datadim==3) return md->NX*md->NY*md->NZ;
    if(md->datadim==2) return md->NX*md->NY       ;
    if(md->datadim==1) return md->NX              ;
    return -1; // error!!!
}

size_t wdata_get_blocksize_bytes(wdata_metadata *md, wdata_variable *var)
{
    if(strcmp(var->type, "real"   ) == 0) return sizeof(double)*wdata_get_blocksize(md)  ;
    if(strcmp(var->type, "complex") == 0) return sizeof(double)*wdata_get_blocksize(md)*2;    
    if(strcmp(var->type, "vector" ) == 0) return sizeof(double)*wdata_get_blocksize(md)*3;    
         
    return 0; // error
}
/**
 * Function adds new block to data file
 * @param md metadata for data set
 * @param var variable to be added to file with name `prefix`_`varname`.wdat
 * @param data pointer to binary data (INPUT)
 * @return 0: ok; 1: cannot open binary file; 2: cannot add datablock to file
 * */
int wdata_add_datablock(wdata_metadata *md, wdata_variable *var, void *data)
{
    char file_name[MD_CHAR_LGTH];
    wdata_get_filename(md, var, file_name);
//     printf("wdata_add_datablock: file_name=%s\n", file_name);
    
    
    FILE *pFile;
    
    pFile= fopen (file_name, "ab");
    if (pFile==NULL)  return 1; // cannot open    
        
    size_t test_ele = fwrite (data , wdata_get_blocksize_bytes(md,var), 1, pFile);
    if(test_ele!=1) return 2; // data not written 
    
    fclose(pFile);
    
    return 0;
}

/**
 * Function adds new block to data file
 * @param md metadata for data set
 * @param varname name of variable, can be from list of vars or links
 * @param data pointer to binary data (INPUT)
 * @return 0: ok; 1: cannot open binary file; 2: cannot add datablock to file, 11: variable is not defined
 * */
int wdata_write_cycle(wdata_metadata *md, const char *varname, void *data)
{
    int ierr;
    wdata_variable var;
    ierr = wdata_get_variable(md,varname, &var);
    if(ierr>0) return 10+ierr;
    
//     wdata_print_variable(&var, stdout); // for testing
    
    char file_name[MD_CHAR_LGTH];
    wdata_get_filename(md, &var, file_name);
//     printf("wdata_write_cycle: Writing to file %s\n", file_name);
    
    
    FILE *pFile;
    
    pFile= fopen (file_name, "ab");
    if (pFile==NULL)  return 1; // cannot open    
        
    size_t test_ele = fwrite (data , wdata_get_blocksize_bytes(md,&var), 1, pFile);
    if(test_ele!=1) return 2; // data not written 
    
    fclose(pFile);
    
    return 0;
}

/**
 * Function read block of data from file
 * @param md metadata for data set
 * @param varname name of variable, can be from list of vars or links
 * @param data pointer to binary data (OUTPUT)
 * @return 0: ok; 1: cannot open binary file; 2: cannot read data block from file; 3: cannot shift pointer;  11: variable is not defined
 * */
int wdata_read_cycle(wdata_metadata *md, const char *varname, int cycle, void *data)
{
    int ierr;
    wdata_variable var;
    ierr = wdata_get_variable(md,varname, &var);
    if(ierr>0) return 10+ierr;
    
//     wdata_print_variable(&var, stdout); // for testing
    
    char file_name[MD_CHAR_LGTH];
    wdata_get_filename(md, &var, file_name);
//     printf("wdata_read_cycle: Reading from file %s\n", file_name);
    
    FILE *pFile;
    
    pFile= fopen (file_name, "rb");
    if (pFile==NULL)  return 1; // cannot open    
        
    // set pointer to correct location
    if(fseek ( pFile, wdata_get_blocksize_bytes(md,&var)*cycle, SEEK_SET ) != 0 ) return 3; // cannot seek pointer
    
    size_t test_ele = fread (data , wdata_get_blocksize_bytes(md,&var), 1, pFile);
    if(test_ele!=1) return 2; // data not read
    
    fclose(pFile);
    
    return 0;
}

int wdata_add_cycle(wdata_metadata *md)
{
    md->cycles++;
    return 0;
}

void wdata_get_filename(wdata_metadata *md, wdata_variable *var, char *file_name)
{
    sprintf(file_name, "%s_%s.wdat", md->prefix, var->name);
}

/**
 * Functions extracts variable corresponding to given name
 * @param md metadata for data set
 * @param varname name of variable, can be from list of vars or links
 * @param var pointer to variable from md structure (OUTPUT)
 * @return 0: ok; 1: cannot find variable
 * */
int wdata_get_variable(wdata_metadata *md, const char *varname, wdata_variable *var)
{
    int i;
    char tvarname[MD_CHAR_LGTH]; // target variable name
    sprintf(tvarname,"%s", varname); // copy to tvarname
    
    // check is links redirects
    for(i=0; i<md->nlinks; i++) if(strcmp(md->links[i].name, varname) == 0) sprintf(tvarname,"%s", md->links[i].linkto);
    
    // find variable
    for(i=0; i<md->nvars; i++) if(strcmp(md->vars[i].name, tvarname) == 0)
    {
        *var = md->vars[i];
        return 0;
    }
    
    // cannot find variable
    return 1;
}

int wdata_get_const(wdata_metadata *md, const char *constname, wdata_const *_const)
{
    int i;
    for(i=0; i<md->nconsts; i++) if(strcmp(md->consts[i].name, constname) == 0) 
    {
        *_const = md->consts[i];
        return 0;
    }
    
    // cannot find const;
    return 1;
}

/**
 * Functions returns value of constant
 * If there is no constant with given name then 0.0 is returned 
 * */
double wdata_getconst(wdata_metadata *md, const char *constname)
{
    int i;
    for(i=0; i<md->nconsts; i++) if(strcmp(md->consts[i].name, constname) == 0) return md->consts[i].value;
    
    return 0.0;
}


/**
 * Functions sets value of constant. If constant was not added before it adds it and sets value.
 * */
void wdata_setconst(wdata_metadata *md, const char *constname, double constvalue)
{
    int i;
    for(i=0; i<md->nconsts; i++) if(strcmp(md->consts[i].name, constname) == 0) { md->consts[i].value=constvalue; return ;}
    
    wdata_const _const;
    strcpy(_const.name,constname);
    _const.value = constvalue;
    md->consts[md->nconsts] = _const;
    md->nconsts++;
    
    return ;
}

/**
 * Function checks if binary file exists for variable
 * @return 1 if binary file for this variable exists, otherwise 0
 */
int wdata_file_exists(wdata_metadata *md, const char *varname)
{
    int ierr;
    wdata_variable var;
    ierr = wdata_get_variable(md, varname, &var);
    if(ierr!=0) return 0;
    
    char file_name[MD_CHAR_LGTH];
    wdata_get_filename(md, &var, file_name);
    
    FILE *file;
    if ((file = fopen(file_name, "r")))
    {
        fclose(file);
        return 1;
    }
    return 0;
}

/**
 * Function removes binary file associated with variable varname
 * */
void wdata_clear_file(wdata_metadata *md, const char *varname)
{
    char file_name[MD_CHAR_LGTH];
    sprintf(file_name, "%s_%s.wdat", md->prefix, varname);
    
    remove(file_name);
}

/**
 * Removes all datafiles, except metadata file
 * */
void wdata_clear_database(wdata_metadata *md)
{
    char file_name[MD_CHAR_LGTH];
    int i;
    for(i=0; i<md->nvars; i++) 
    {
        wdata_get_filename(md, &md->vars[i], file_name);
        remove(file_name);
    }
    md->cycles=0;
}

#ifdef WDATA_TESTING_MODE

// ===========================================================================
// ============== FOR TESTING ONLY ===========================================
// ===========================================================================
#include <math.h>
#define SQRT_2PI 2.506628274631
double function_x(double x, double sigma)
{
    return 1./(sigma*SQRT_2PI) * exp(-0.5*x*x/(sigma*sigma));
}

double function_xyz(double x, double y, double z, double time)
{
    double val=0.0;
    double sigmax = 2.0 + 0.20*time;
    double sigmay = 3.0 + 0.15*time;
    double sigmaz = 4.0 + 0.10*time;
    
    val = function_x(x, sigmax)*function_x(y, sigmay)*function_x(z, sigmaz);
    
    return val;
}

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    }
    
int main()
{
    int err;
    int i, ix, iy, iz, ixyz;
    
    // create artificial data for visulisation in visit
    wdata_metadata md = {24, 28, 32, 1.0, 1.0, 1.0, 3, "testa", 0, 0.0, 1.0, 0, 0};
    
    wdata_variable vdensity_a = {"density_a", "real", "testunit"};
    wdata_add_variable(&md, &vdensity_a);
    
    wdata_variable vdelta = {"delta", "complex", "none"};
    wdata_add_variable(&md, &vdelta);
    
    wdata_variable vcurrent_a = {"current_a", "vector", "none"};
    wdata_add_variable(&md, &vcurrent_a);
    
    wdata_link ldensity_b = {"density_b", "density_a"};
    wdata_add_link(&md, &ldensity_b);
    
    wdata_link lcurrent_b = {"current_b", "current_a"};
    wdata_add_link(&md, &lcurrent_b);
    
    wdata_const lconst_eF = {"eF", 0.1};
    wdata_add_const(&md, &lconst_eF);
    
    wdata_const lconst_kF = {"kF", 1.1};
    wdata_add_const(&md, &lconst_kF);
    
    // just in case - clear data sets
    wdata_clear_database(&md);

    char file_name[256];

    
    // add artificial data to sets
    double *dataR;
    double complex *dataC;
    double *dataV;
    cppmallocl(dataR,md.NX*md.NY*md.NZ,double);
    cppmallocl(dataC,md.NX*md.NY*md.NZ,double complex);
    cppmallocl(dataV,md.NX*md.NY*md.NZ*3,double);
    
    int bdim = wdata_get_blocksize(&md);
    
    int ncycles=10;

    for(i=0; i<ncycles; i++)
    {
        ixyz=0;
        for(ix=0; ix<md.NX; ix++) for(iy=0; iy<md.NY; iy++) for(iz=0; iz<md.NZ; iz++)
        {
            double x = md.DX*(ix-md.NX/2);
            double y = md.DY*(iy-md.NY/2);
            double z = md.DZ*(iz-md.NZ/2);
            double time = md.t0 + md.dt*i;
            dataR[ixyz] = function_xyz(x,y,z,time);
            dataC[ixyz] = -1.0*dataR[ixyz] + I*0.0;
            
            // vector 
            dataV[ixyz+0*bdim] = -1.0*y;
            dataV[ixyz+1*bdim] =  1.0*x;
            dataV[ixyz+2*bdim] =  0.0;
            
            ixyz++;
        }
        
        // add cycle to binary sets
        err = wdata_add_datablock(&md, &vdensity_a, dataR);
        printf("A err=%d\n", err);
        err = wdata_write_cycle(&md, "delta", dataC);
        printf("B err=%d\n", err);
        err = wdata_write_cycle(&md, "current_a", dataV);
        printf("C err=%d\n", err);
        wdata_add_cycle(&md);
    }
    
    // write metadata file
    sprintf(file_name, "%s_info.wtxt", md.prefix);
    printf("creating %s\n", file_name);
    FILE * fout = fopen(file_name, "w");
    wdata_print_metadata(&md, fout);
    fclose(fout);


//     // test of reading
// //     for(i=0; i<ncycles+2; i++)
//     {
//         err = wdata_read_cycle(&md, "density_a", i, dataR);
//         printf("i=%d, err=%d\n", i, err);
//         
//         ixyz=0;
//         for(ix=0; ix<md.NX; ix++) for(iy=0; iy<md.NY; iy++) for(iz=0; iz<md.NZ; iz++)
//         {
//             printf("%d %d %d %f\n", ix, iy, iz, dataR[ixyz]);
//             ixyz++;
//         }
//     }
    
    return 0;
}

// int main()
// {
//     int err;
//     int i, ix, iy, iz, ixyz;
//     
//     // create artificial data for visulisation in visit
//     wdata_metadata md;
//     char file_name[256] = "testa_info.wtxt";
//     
//     // write metadata file
//     wdata_parse_metadata_file(file_name, &md);
//     
//     double *dataR;
//     double complex *dataC;
//     cppmallocl(dataR,md.NX*md.NY*md.NZ,double);
//     cppmallocl(dataC,md.NX*md.NY*md.NZ,double complex);
//     
//     wdata_print_metadata(&md,stdout);
//     
//     i=0;
//     
//     {
//         err = wdata_read_cycle(&md, "density_a", i, dataR);
//         printf("i=%d, err=%d\n", i, err);
//         
//         ixyz=0;
// //         for(ix=0; ix<md.NX; ix++) for(iy=0; iy<md.NY; iy++) for(iz=0; iz<md.NZ; iz++)
//         {
//             printf("%d %d %d %f\n", ix, iy, iz, dataR[ixyz]);
//             ixyz++;
//         }
//     }
// }

#endif

