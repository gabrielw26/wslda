# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.1
#   kernelspec:
#     display_name: Python [conda env:work]
#     language: python
#     name: conda-env-work-py
# ---

# +
from importlib import reload
import numpy as np
import wdata.io;reload(wdata.io)

dim = 3
ws = (1.0,2.0,3.0)
Nxyz = (32, 33, 34)[:dim]
Lxyz = (2*np.pi,)*dim
dxyz = np.divide(Lxyz, Nxyz)

t0 = 0.1
dt = 1.1
Nt = 10

xyz = np.meshgrid(
    *((np.arange(_N)-_N/2)*_dx
      for _N, _dx in zip(Nxyz, dxyz)),
    indexing='ij', sparse=True)
ts = np.arange(Nt)*dt + t0

def f(ws, xyz, t, d=0):
    x, y, z = xyz
    wx, wy, wz = (w*t for w in ws)
    cos, sin = np.cos, np.sin
    if d == 0:
        return cos(wx*x)*cos(wy*y)*cos(wz*z)
    else:
        return [
            -wx*sin(wx*x)*cos(wy*y)*cos(wz*z),
            -wy*cos(wx*x)*sin(wy*y)*cos(wz*z),
            -wz*cos(wx*x)*cos(wy*y)*sin(wz*z)]
    
density = [f(ws, xyz, t) for t in ts]
delta = [1j*f(ws, xyz, t) for t in ts]
current = [f(ws, xyz, t, d=1) for t in ts]

# Automatic
res = wdata.io.WData(
    prefix="tmp",
    xyz=xyz,
    t=ts,
    variables=[
        wdata.io.Var(density=density), 
        wdata.io.Var(delta=delta),
        wdata.io.Var(current_a=current),
    ],
    aliases={'density_b': 'density_a',
             'current_b': 'current_a'},
    constants=dict(eF=0.5, kF=1),
)
print(res.get_metadata())
# -

data = wdata.io.WData.load('tmp.wtxt')
print(data.description)
data.density.shape
