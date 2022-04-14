#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/***************************************************************************/ 
/**************************** GPE HEADERS **********************************/
/***************************************************************************/
//#include "gpe_engine.h"
#include "predefines.h"
#include "gpe_utils.h"
//#include "pca_utils.h"
//#include "wslda_toolkit.h"
// -> include "pca_utils.h", "gpe_user_defined.h"
//extern int wsldapid;
int wsldapid;

void help_msg()
{
    printf("Not correct arguments. Try: \n./gpe real\n./gpe imag\n");
    exit(1);
}

/***************************************************************************/ 
/************************* MAIN FUNCTION  **********************************/
/***************************************************************************/
int main( int argc , char ** argv ) 
{
    int err = 0;
    char *program_type;
    const char *imag = "imag";
    const char *real = "real";
    if(argc != 2) help_msg();

    if(0 == strcmp(imag, argv[1])) {
        err=gpe_imag();
    }
    else if(0 == strcmp(real, argv[1])) {
        err=gpe_real();
    }
    else {
        help_msg();
    } 
    return err;
}