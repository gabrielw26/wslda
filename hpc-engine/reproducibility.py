#!/usr/bin/python

# W-SLDA Toolkit
# @author Gabriel Wlazlowski

from string import *
import sys
import os

user_file=["predefines.h", "problem-definition.h", "logger.h"]
functions=["predefines_h", "problem_definition_h", "logger_h"]

targetdir = sys.argv[1]
targetfile = os.path.join(targetdir, "reproducibility.h")
#print targetfile

# create file
f = open(targetfile, "w")

for fin, funout in zip(user_file,functions):
    #print fin, funout
    f.write("void %s(FILE *fptr)\n" % funout)
    f.write("{\n")
    
    # insert file content
    src = open(fin, "r")
    for line in src:
        line = line.replace("\n", "")
        line = line.replace("\\", "\\\\")
        line = line.replace("\"", "\\\"")
        line = line.replace("\'", "\\\'")
        line = line.replace("%", "%%")
        
        f.write('fprintf(fptr, "%s\\n");\n' % line)
    src.close()
    
    f.write("}\n")

# close file
f.close()



