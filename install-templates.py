#!/usr/bin/python

# Script for installing templates

import sys
import os

if __name__ == '__main__':
    
    WSLDA=os.environ['WSLDA']
    if WSLDA=="":
        print("WSLDA not set!")
        print("export WSLDA=...")
        print("before usage of this script!")
        sys.exit(1)
    
    argc=len(sys.argv)
    if argc!=2:
        print("Use:")
        print("%s system" % sys.argv[0])
        print("where system is:")
        dirs = os.listdir("%s/templates" % WSLDA)
        for d in dirs:
            if d in ["legacy", "README.md"]: continue
            print("\t%s" % d)
        print("\tclean <-- use this to clear templates")
        sys.exit(1)
        
    
    if sys.argv[1]=="clean":
        print("Cleaning...")
        ftorm=['env.sh', 'Makefile', 'job.sh', 'machine.h', 'README.md']
        folders=['st-project-template','st-testcase-uniform','td-project-template','td-testcase-uniform']
        for d in folders:
            for f in ftorm:
                cmd = "rm %s/%s/%s" % (WSLDA,d,f)
                print(cmd)
                os.system(cmd)
        
        print("Done.")
        sys.exit(1)
        
    
    # Copy templates
    for ss in ['st-project-template','st-testcase-uniform']:
        cmd = "cp -r %s/templates/%s/st/* %s/%s" % (WSLDA,sys.argv[1],WSLDA,ss)
        print(cmd)
        os.system(cmd)
        cmd = "cp %s/templates/%s/machine.h %s/%s" % (WSLDA,sys.argv[1],WSLDA,ss)
        print(cmd)
        os.system(cmd)
        
    for ss in ['td-project-template','td-testcase-uniform']:
        cmd = "cp -r %s/templates/%s/td/* %s/%s" % (WSLDA,sys.argv[1],WSLDA,ss)
        print(cmd)
        os.system(cmd)
        cmd = "cp %s/templates/%s/machine.h %s/%s" % (WSLDA,sys.argv[1],WSLDA,ss)
        print(cmd)
        os.system(cmd)
        
    print("Done.")

