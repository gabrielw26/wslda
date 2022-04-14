#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/***************************************************************************/ 
/**************************** GPE HEADERS **********************************/
/***************************************************************************/
#include "predefines.h"
#include "gpe_utils.h"
int wsldapid;

/***************************************************************************/ 
/***************************** FUNCTIONS ***********************************/
/***************************************************************************/
void help_msg(char *prog_name)
{
    printf("Not correct arguments. Try: \n%s real\n%s imag\n", prog_name, prog_name);
    exit(1);
}

/***************************************************************************/ 
/************************* MAIN FUNCTION  **********************************/
/***************************************************************************/
int main( int argc , char ** argv ) 
{
    int err = 0, type = -1;
    const char *imag = "imag";
    const char *real = "real";
    if(argc != 2) help_msg(argv[0]);

    if(0 == strcmp(imag, argv[1])) type = 0;
    else if(0 == strcmp(real, argv[1])) type = 1;
    else help_msg(argv[0]);

    err=gpe_compute(type);
    return err;
}