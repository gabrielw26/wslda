#!/usr/bin/python

# Import packages
from string import *
import os

for ss in [0,1]:
    for d in [1,2,3]:
        ifname="test.ss%d-%dd.cmp" % (ss,d)
        ofname=ifname+".ref"
        
        print("%s --> %s" % (ifname,ofname))
        
        fin = open(ifname, "r")
        fout= open(ofname, "w")
        
        for line in fin:
            line=line.replace("\n", "")
            line=line+"\t 1.0e-4"
            fout.write("%s\n" % line)
            
        fin.close()
        fout.close()
        
        cmd="rm %s" % (ifname)
        print(cmd)
        os.system(cmd)
        
        ifname="test.ss%d-%dd.wlog" % (ss,d)
        ofname=ifname+".ref"
        cmd="mv %s %s" % (ifname,ofname)
        print(cmd)
        os.system(cmd)
        
