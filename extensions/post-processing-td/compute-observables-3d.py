#!/bin/python
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from wdata.io import WData, Var

wtxt='./ufg/qtmrg.wtxt'
spindeg=2.0 # two-component system

data=WData.load(wtxt)

print("%12s %14s %14s %14s" % ("t*eF","N","Eflow","Econd"))
for it in range(data.Nt):
    t = data.t[it]*data.eF # time
    rho=data.rho_a[it,:,:,:]*spindeg # density
    delta=np.abs(data.delta[it,:,:,:]) # order parameter (absolute value)
    j2=(data.j_a[it,0,:,:,:]*spindeg)**2 + (data.j_a[it,1,:,:,:]*spindeg)**2 +  (data.j_a[it,2,:,:,:]*spindeg)**2 # current square
    N = np.sum(rho)
    Eflow = np.sum(j2/(2.*rho)) # flow energy
    kF = (3.*np.pi**2*rho)**(1./3.)
    eF = 0.5*kF**2
    Econd=(3./8.)*np.sum(delta**2 * rho/eF)
    print("%12.8f %14.8f %14.8f %14.8f" % (t,N,Eflow,Econd))


