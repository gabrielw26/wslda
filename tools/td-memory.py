#!/usr/bin/python

# Script computes memory utilization by bsk_dynamic3d code

# Import packages
from math import *
from string import *
import numpy as np
import matplotlib.pyplot as plt 

# SETTINGS
NX = 100
NY = 100
NZ = 100
codedim=3 # dimensonality of code
nwf=None # provide here number if you know it, otherwise the code will use simple estimate
mem_per_gpu = 16.0 # in GB
min_mem_utilization = 2.0 # in GB


# ----------------------------------------------------------------
# --------------------- DO NOT MODIFY HERE -----------------------
# ----------------------------------------------------------------
if nwf==None:
    if codedim==3: nwf = 0.5*NX*NY*NZ # or use simple estimate
    elif codedim==2: nwf = 0.5*NX*NY*NZ
    elif codedim==2: nwf = 0.5*NX*NY*NZ
    else: raise NameError('incorrect codedim')

if codedim==3:
    pass
elif codedim==2:
    NZ=1
elif codedim==2:
    NZ=1
    NY=1
else: 
    raise NameError('incorrect codedim')

def find_mem_per_gpu(ngpus):
    """
    returns in GB
    """
    nwfip = nwf / ngpus
    memwf = nwfip*2*NX*NY*NZ*16 # in bytes
    total = 11*memwf # contribution from main buffers
    total = total + 20*NX*NY*NZ*8 # densities storage
    total = total + 4*NX*NY*NZ*8  # potentials storage
    total = total + nwfip*8  # d_fbetaEn storage
    total = total + nwfip*8  # d_qpe storage
    
    return 1.0*total / 2**30


gpus=[]
mem=[]
ngpus=0
while 1:
    ngpus = ngpus+1
    m = find_mem_per_gpu(ngpus)
    
    if m<=mem_per_gpu:
        gpus.append(ngpus)
        mem.append(m)
        
    if m<min_mem_utilization:
        break
    
# print minimal request
print("MINIMAL NUMBER OF GPUs=%d" % gpus[0])

# Data for plotting

fig, ax = plt.subplots()
ax.plot(gpus, mem)

ax.set(xlabel='# GPUs', ylabel='memory per GPU [GB]',
       title="NX=%d, NY=%d, NZ=%d, nwf=%d" % (NX, NY, NZ, nwf))
ax.grid()

plt.show()

