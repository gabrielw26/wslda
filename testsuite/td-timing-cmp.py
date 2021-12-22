#!/usr/bin/python

# This script can be used to compare timings of td codes
# Instructions:
# grep 3.5000 testsuite.out > old.txt  # <-- do it for `old` version of the code
# grep 3.5000 testsuite.out > new.txt  # <-- do it for `new` version of the code
# ./td-timing-cmp.py

from string import *

new=[]
f=open("new.txt", "r")
for l in f:
    l=l.replace("\n", "")
    l=split(l," ")
    new.append(atof(l[-1]))
f.close()
print new

old=[]
f=open("old.txt", "r")
for l in f:
    l=l.replace("\n", "")
    l=split(l," ")
    old.append(atof(l[-1]))
f.close()
print old

print "| %8s | %8s | %8s | %8s |" % ("test", "old", "new", "old/new")
print "-" * 45
for i in range(len(old)):
    print "| %8d | %8.2f | %8.2f | %8.2f |" % (i+1, old[i], new[i], old[i]/new[i])
