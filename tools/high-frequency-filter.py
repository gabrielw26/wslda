#!/usr/bin/python

# Script computes memory utilization by bsk_dynamic3d code

# Import packages
from math import *
from string import *
import numpy as np
import matplotlib.pyplot as plt 
from scipy.fftpack import fft, fftfreq, fftshift, ifft, ifftshift


# Settings
NX=512
DX=1.0
hff_mu=0.9
hff_T=0.01
fname="high-frequency-filter"

# Domain
x = np.arange(-DX*NX/2, DX*NX/2, DX)

# test function
U = -1.0 / (np.exp((np.abs(x)-0.6*DX*NX/2)/5.0) + 1.0)

noise=(np.random.random(NX) - 0.5)*0.5
U = U + noise

# fourier space
Uk = fft(U) 
xk = fftfreq(NX, DX) * 2. * np.pi  # fftfreq does not multiply it by 2pi, do it by hand
Uk = fftshift(Uk)
xk = fftshift(xk)

# filter
kc=np.pi/DX
ec=0.5*kc*kc
mu=hff_mu*ec
T=hff_T*ec
ek=0.5*xk*xk
fil = 1.0/(np.exp((ek-mu)/T)+1.0)
Ukf=Uk*fil


# go back to coordinate space
Uf = ifft(ifftshift(Ukf))
Uf = np.real(Uf)


# make plot
plt.rcParams['font.size'] = 6
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(4.0,4.0))
ax1.plot(x,U, color="red", lw=2.0, label="input: U(x)")
ax1.plot(x,Uf, color="blue", lw=1.0, label="filtered: Uf(x)")
ax1.legend()
ax1.grid(True)

ax2.plot(xk,np.abs(Uk), color="green", lw=2.0, label="U(k)")
ax2.plot(xk,np.abs(Ukf), color="yellow", lw=1.0, label="U(k)*FD(k)")
ax2.legend()
ax2.grid(True)

ax3.plot(ek,fil, color="orange", lw=2.0, label="FD(%.2f,%.2f)" % (hff_mu,hff_T))
ax3.legend()
ax3.grid(True)

#plt.show()
fname=fname+".png"
print("Saving to file: "+fname)
#plt.tight_layout()
plt.savefig(fname, dpi=300)

