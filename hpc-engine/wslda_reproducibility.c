/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 04-10-2020
 * */ 

#include <stdio.h>
#include "reproducibility.h"
#include "pca_settings.h"

int assure_reproducibility(const char *prefix)
{
    FILE * f;
    char file_name[512];
    
    sprintf(file_name, "%s_logger.h", prefix);
    f = fopen(file_name, "w");
    if(f==NULL) return 1;
    fprintf(f, "/**\n");
    fprintf(f, " * W-SLDA Toolkit\n");
    fprintf(f, " * Engine version: %s\n", VERSION);
    fprintf(f, " * */\n\n");
    logger_h(f);
    fclose(f);
    
    sprintf(file_name, "%s_problem-definition.h", prefix);
    f = fopen(file_name, "w");
    if(f==NULL) return 2;
    fprintf(f, "/**\n");
    fprintf(f, " * W-SLDA Toolkit\n");
    fprintf(f, " * Engine version: %s\n", VERSION);
    fprintf(f, " * */\n\n");
    problem_definition_h(f);
    fclose(f);
    
    sprintf(file_name, "%s_predefines.h", prefix);
    f = fopen(file_name, "w");
    if(f==NULL) return 3;
    fprintf(f, "/**\n");
    fprintf(f, " * W-SLDA Toolkit\n");
    fprintf(f, " * Engine version: %s\n", VERSION);
    fprintf(f, " * */\n\n");
    predefines_h(f);
    fclose(f); 
    
    return 0;
}
