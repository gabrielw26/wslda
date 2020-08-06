/**
 * W-SLDA Toolkit
 * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h> 
#include <libgen.h>
#include <complex.h>
#include "wdata.h"

#define MAX_REC_LEN 1024

#define WRKDIR_SET 123987

char __wdata__basedir[MAX_REC_LEN];
void wdata_goto_wrkdir(wdata_metadata *md)
{
    
    if(md->issetwrkdir==WRKDIR_SET)
    { 
        getcwd(__wdata__basedir, MAX_REC_LEN);
        chdir(md->wrkdir);
//         printf("changing to workdir=%s\n", md->wrkdir);
    }                                                                           
}

void wdata_goback_wrkdir(wdata_metadata *md)
{
    if(md->issetwrkdir==WRKDIR_SET) 
    {
        chdir(__wdata__basedir);
//         printf("changing to basedir=%s\n", __wdata__basedir);
    }
}

int wdata_parse_metadata_file(const char * file_name, wdata_metadata *md)
{
    FILE *fp;
    fp=fopen(file_name, "r");
    if(fp==NULL)
        return 1; // Cannot open file
        
    // reset vars
    md->nvar=0;
    md->nlink=0;
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
            sscanf (s,"%s %s %s %s %s %*s",tag, &md->var[md->nvar].name, &md->var[md->nvar].type, &md->var[md->nvar].unit, &md->var[md->nvar].format ); 
            md->nvar++;       
        }
        
        // links
        else if (strcmp (tag,"link") == 0)
        {
            sscanf (s,"%s %s %s %*s",tag, &md->link[md->nlink].name, &md->link[md->nlink].linkto); 
            md->nlink++;       
        }
        
        // consts
        else if (strcmp (tag,"const") == 0)
        {
            sscanf (s,"%s %s %lf %*s",tag, &md->consts[md->nconsts].name, &md->consts[md->nconsts].value); 
            md->nconsts++;       
        }
    }
    
    fclose(fp);
    
    // set working dir to be consistent with wtxt file location
    char ctmp[MD_CHAR_LGTH];
    strcpy(ctmp, file_name); 
    strcpy(md->wrkdir, dirname(ctmp)); 
    md->issetwrkdir = WRKDIR_SET;
    
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
    fprintf(out, "prefix %20s   # prefix for files belonging to this data set, binary files have names prefix_variable.format\n", md->prefix);
    fprintf(out, "cycles %20d   # number of cycles (measurements)\n", md->cycles);
    fprintf(out, "t0 %24g   # time value for the first cycle\n", md->t0);
    fprintf(out, "dt %24g   # time interval between cycles\n", md->dt);

    // variables
    fprintf(out,"\n");
    fprintf(out,"# variables\n");
    fprintf(out,"# tag                  name                    type                    unit                  format\n");
    for(i=0; i<md->nvar; i++) wdata_print_variable(&md->var[i], out);
    
    // links
    fprintf(out,"\n");
    fprintf(out,"# links\n");
    fprintf(out,"# tag                  name                 link-to\n");
    for(i=0; i<md->nlink; i++) wdata_print_link(&md->link[i], out);
    
    // consts
    fprintf(out,"\n");
    fprintf(out,"# consts\n");
    fprintf(out,"# tag                  name                   value\n");
    for(i=0; i<md->nconsts; i++) wdata_print_const(&md->consts[i], out);
    
    fprintf(out,"\n");
}

void wdata_print_variable(wdata_variable *md, FILE *out)
{
    fprintf(out, "var%24s%24s%24s%24s\n", md->name, md->type, md->unit, md->format);
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
    md->var[md->nvar] = *var;
    md->nvar++;
}

void wdata_add_link(wdata_metadata *md, wdata_link *link)
{
    md->link[md->nlink] = *link;
    md->nlink++;
}


void wdata_add_const(wdata_metadata *md, wdata_const *_const)
{
    md->consts[md->nconsts] = *_const;
    md->nconsts++;
}

int wdata_get_blocklength(wdata_metadata *md)
{
    if(md->datadim==3) return md->NX*md->NY*md->NZ;
    if(md->datadim==2) return md->NX*md->NY       ;
    if(md->datadim==1) return md->NX              ;
    return -1; // error!!!
}

size_t wdata_get_blocksize(wdata_metadata *md, wdata_variable *var)
{
    if(strcmp(var->type, "real"   ) == 0) return sizeof(double)*wdata_get_blocklength(md)  ;
    if(strcmp(var->type, "complex") == 0) return sizeof(double)*wdata_get_blocklength(md)*2;    
    if(strcmp(var->type, "vector" ) == 0) return sizeof(double)*wdata_get_blocklength(md)*3;    
         
    return 0; // error
}

int wdata_add_datablock(wdata_metadata *md, wdata_variable *var, void *data)
{
    wdata_goto_wrkdir(md);
    
    int ierr=-1;
    if     (strcmp(var->format, "wdat" ) == 0) ierr = wdata_add_datablock_wdat(md,var, data);
    else if(strcmp(var->format, "dpca" ) == 0) ierr = wdata_add_datablock_dpca(md,var, data);
    else if(strcmp(var->format, "npy"  ) == 0) ierr = wdata_add_datablock_npy(md,var, data);
    
    wdata_goback_wrkdir(md);
    
    return ierr; // error code
}

int wdata_add_datablock_wdat(wdata_metadata *md, wdata_variable *var, void *data)
{
    char file_name[MD_CHAR_LGTH];
    wdata_get_filename(md, var, file_name);
//     printf("wdata_add_datablock: file_name=%s\n", file_name);
        
    FILE *pFile;
    
    pFile= fopen (file_name, "ab");
    if (pFile==NULL)  return 1; // cannot open    
        
    size_t test_ele = fwrite (data , wdata_get_blocksize(md,var), 1, pFile);
    if(test_ele!=1) return 2; // data not written 
    
    fclose(pFile);
    
    return 0;
}

int wdata_add_datablock_dpca(wdata_metadata *md, wdata_variable *var, void *data)
{
    return -99; // pdca format is deprecated, thus writing files in this format is no longer supported!
}

int wdata_add_datablock_npy(wdata_metadata *md, wdata_variable *var, void *data)
{
    return -99; // TODO
}

int wdata_load_datablock(wdata_metadata *md, wdata_variable *var, int cycle, void *data)
{
    wdata_goto_wrkdir(md);
    
    int ierr=-1;
    if     (strcmp(var->format, "wdat" ) == 0) ierr = wdata_load_datablock_wdat(md,var, cycle, data);
    else if(strcmp(var->format, "dpca" ) == 0) ierr = wdata_load_datablock_dpca(md,var, cycle, data);
    else if(strcmp(var->format, "npy"  ) == 0) ierr = wdata_load_datablock_npy(md,var, cycle, data);
    
    wdata_goback_wrkdir(md);
    
    return ierr; // error code    
}

int wdata_load_datablock_wdat(wdata_metadata *md, wdata_variable *var, int cycle, void *data)
{
    char file_name[MD_CHAR_LGTH];
    wdata_get_filename(md, var, file_name);
//     printf("wdata_read_cycle: Reading from file %s\n", file_name);
    
    FILE *pFile;
    
    pFile= fopen (file_name, "rb");
    if (pFile==NULL)  return 1; // cannot open    
        
    // set pointer to correct location
    if(fseek ( pFile, wdata_get_blocksize(md,var)*cycle, SEEK_SET ) != 0 ) return 3; // cannot seek pointer
    
    size_t test_ele = fread (data , wdata_get_blocksize(md,var), 1, pFile);
    if(test_ele!=1) return 2; // data not read
    
    fclose(pFile);

    return 0; 
}

int wdata_load_datablock_dpca(wdata_metadata *md, wdata_variable *var, int cycle, void *data)
{
    char file_name[MD_CHAR_LGTH];
    wdata_get_filename(md, var, file_name);
//     printf("wdata_read_cycle: Reading from file %s\n", file_name);
    
    FILE *pFile;
    
    pFile= fopen (file_name, "rb");
    if (pFile==NULL)  return 1; // cannot open    
        
    // set pointer to correct location
    #define MIO_CNT_INTS 4
    #define MIO_CNT_double 6
    size_t header_size = (MIO_CNT_INTS+1)*sizeof(int) + MIO_CNT_double*sizeof(double);
    if(fseek ( pFile, header_size + wdata_get_blocksize(md,var)*cycle, SEEK_SET ) != 0 ) return 3; // cannot seek pointer
    
    size_t test_ele = fread (data , wdata_get_blocksize(md,var), 1, pFile);
    if(test_ele!=1) return 2; // data not read
    
    fclose(pFile);

    return 0;
}

int wdata_load_datablock_npy(wdata_metadata *md, wdata_variable *var, int cycle, void *data)
{
    return -99; // TODO
}

int wdata_write_cycle(wdata_metadata *md, const char *varname, void *data)
{
    
    int ierr;
    wdata_variable var;
    ierr = wdata_get_variable(md,varname, &var);
    if(ierr>0) return 10+ierr;
    
    return wdata_add_datablock(md, &var, data);
}

int wdata_read_cycle(wdata_metadata *md, const char *varname, int cycle, void *data)
{
    int ierr;
    wdata_variable var;
    ierr = wdata_get_variable(md,varname, &var);
    if(ierr>0) return 10+ierr;
    
    return wdata_load_datablock(md, &var, cycle, data);
}

int wdata_add_cycle(wdata_metadata *md)
{
    md->cycles++;
    return 0;
}

void wdata_get_filename(wdata_metadata *md, wdata_variable *var, char *file_name)
{
    if     (strcmp(var->format, "wdat" ) == 0) sprintf(file_name, "%s_%s.wdat", md->prefix, var->name);
    else if(strcmp(var->format, "dpca" ) == 0) sprintf(file_name, "%s_%s.dpca", md->prefix, var->name);
    else if(strcmp(var->format, "npy"  ) == 0) sprintf(file_name, "%s_%s.npy", md->prefix, var->name);
}

int wdata_get_variable(wdata_metadata *md, const char *varname, wdata_variable *var)
{
    int i;
    char tvarname[MD_CHAR_LGTH]; // target variable name
    sprintf(tvarname,"%s", varname); // copy to tvarname
    
    // check is links redirects
    for(i=0; i<md->nlink; i++) if(strcmp(md->link[i].name, varname) == 0) sprintf(tvarname,"%s", md->link[i].linkto);
    
    // find variable
    for(i=0; i<md->nvar; i++) if(strcmp(md->var[i].name, tvarname) == 0)
    {
        *var = md->var[i];
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

double wdata_getconst_value(wdata_metadata *md, const char *constname)
{
    int ierr=0;
    wdata_const _const;
    ierr = wdata_get_const(md, constname, &_const);
    if(ierr==0) return _const.value;
    else return 0.0;
}


double wdata_getconst(wdata_metadata *md, const char *constname)
{
    int i;
    for(i=0; i<md->nconsts; i++) if(strcmp(md->consts[i].name, constname) == 0) return md->consts[i].value;
    
    return 0.0;
}


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

int wdata_file_exists(wdata_metadata *md, const char *varname)
{
    wdata_goto_wrkdir(md);
    
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
    
    wdata_goback_wrkdir(md);
    
    return 0;
}


void wdata_clear_file(wdata_metadata *md, const char *varname)
{
    wdata_goto_wrkdir(md);
    
    char file_name[MD_CHAR_LGTH];
    sprintf(file_name, "%s_%s.wdat", md->prefix, varname);
    
    remove(file_name);
    
    wdata_goback_wrkdir(md);
}


void wdata_clear_database(wdata_metadata *md)
{
    wdata_goto_wrkdir(md);
    
    char file_name[MD_CHAR_LGTH];
    int i;
    for(i=0; i<md->nvar; i++) 
    {
        wdata_get_filename(md, &md->var[i], file_name);
        remove(file_name);
    }
    md->cycles=0;
    
    wdata_goback_wrkdir(md);
}

int wdata_write_metadata_to_file(wdata_metadata *md, const char * filename)
{
    char file_name[MD_CHAR_LGTH];
    if(strcmp(filename,"")==0) sprintf(file_name, "%s.wtxt", md->prefix);
    else                       sprintf(file_name, "%s", filename);
    
//     printf("creating %s\n", file_name);
    FILE * fout = fopen(file_name, "w");
    if(fout==NULL) return 1;
    wdata_print_metadata(md, fout);
    fclose(fout);
}

void wdata_set_working_dir(wdata_metadata *md, const char * wrkdir)
{
    strcpy(md->wrkdir, wrkdir);
    md->issetwrkdir=WRKDIR_SET;
}


int wdata_add_var_to_metadata_file(const char * file_name, wdata_variable *var)
{    
    FILE * fout = fopen(file_name, "a");
    if(fout==NULL) return 1;
    wdata_print_variable(var, fout);
    fclose(fout);
    return 0;
}


int wdata_add_link_to_metadata_file(const char * file_name, wdata_link *link)
{
    FILE * fout = fopen(file_name, "a");
    if(fout==NULL) return 1;
    wdata_print_link(link, fout);
    fclose(fout);
    return 0;
}


int wdata_add_const_to_metadata_file(const char * file_name, wdata_const *_const)
{
    FILE * fout = fopen(file_name, "a");
    if(fout==NULL) return 1;
    wdata_print_const(_const, fout);
    fclose(fout);
    return 0;
}
