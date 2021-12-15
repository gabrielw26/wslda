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
    
    
# !!! special case: machine.h !!!
funout="machine_h"
f.write("void %s(FILE *fptr)\n" % funout)
f.write("{\n")

# test where is machine.h located
if os.path.isfile("machine.h"):
    fin="machine.h"
elif 'WSLDA_MACHINE' in os.environ:
    if os.path.isfile(os.environ['WSLDA_MACHINE']+"/machine.h"):
        fin=os.environ['WSLDA_MACHINE']+"/machine.h"
    else:
        fin="none"
else:
    fin="none"

if fin=="none":
    line="// CANNOT LOCATE machine.h. YOU NEED TO PROVIDE IT MANUALLY."
    f.write('fprintf(fptr, "%s\\n");\n' % line)
else:
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



