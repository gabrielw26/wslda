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
            if d=="legacy": continue
            print("\t%s" % d)
        sys.exit(1)
        
    
    # Copy templates
    cmd = "cp -r %s/templates/%s/* %s" % (WSLDA,sys.argv[1],WSLDA)
    print cmd
    os.system(cmd)
    print("Done.")
    
