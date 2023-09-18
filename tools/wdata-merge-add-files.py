#!/usr/bin/python

# Script adds missing files after applying wdata-merge-tool

# Import packages
from math import *
from string import *
import os
import getopt, sys


# Remove 1st argument from the
# list of command line arguments
argumentList = sys.argv[1:]
options = "o:"
long_options = ["output="]

outprefix=None
try:
    # Parsing argument
    arguments, inprefixex = getopt.getopt(argumentList, options, long_options)
    # print(arguments)
    # print(inprefixex)

    # checking each argument
    for currentArgument, currentValue in arguments:

        if currentArgument in ("-o", "--output"):
            outprefix=currentValue

except getopt.error as err:
    # output error, and return with an error code
    print (str(err))
    sys.exit(1)

if outprefix==None or len(outprefix)==0:
    print("Usage: %s inprefix1 inprefix2 inprefix3 ... -o output_prefix" % sys.argv[0])
    sys.exit(1)


# Start merging of files
def runcmd(cmd):
    print(cmd)
    os.system(cmd)

def create_empty_file(fname):
    print("# CREATING EMPTY FILE: `%s`" % fname)
    f = open(fname, "w")
    f.close()

def add_to_file_another_file(fname, addfname, header='', start_line=''):
    print("# ADDING TO FILE: `%s` CONTENT OF `%s`" % (fname,addfname))
    f = open(fname, "a")

    # add header
    if header!='':
        f.write(header)

    # add content of another file
    fadd = open(addfname,"r")
    for l in fadd:
        f.write("%s%s" % (start_line,l))

    f.close()
    fadd.close()


# Merge wlog files
outfname=outprefix+'.wlog'
infnames=[s+'.wlog' for s in inprefixex]
create_empty_file(outfname)
for infname in infnames:
    add_to_file_another_file(outfname, infname, header='', start_line='')

# Merge stdout files
outfname=outprefix+'.stdout'
infnames=[s+'.stdout' for s in inprefixex]
create_empty_file(outfname)
for infname in infnames:
    add_to_file_another_file(outfname, infname, header='', start_line='')

# Merge input files
outfname=outprefix+'_input.txt'
infnames=[s+'_input.txt' for s in inprefixex]
create_empty_file(outfname)
for i,infname in enumerate(infnames):
    if i==0:
        add_to_file_another_file(outfname, infname, header='', start_line='')
    else:
        add_to_file_another_file(outfname, infname, header='\n\n### --- INPUT FOR CONTINUATION OF CALCULATIONS ---\n', start_line='#')


# Copy other files
sufixes=['predefines.h', 'problem-definition.h', 'logger.h', 'machine.h']
for s in sufixes:
    cmd = "cp %s_%s %s_%s" % (inprefixex[0],s,outprefix,s)
    runcmd(cmd)

print("Done.")

