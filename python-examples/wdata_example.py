import numpy as np
from wdata.io import WData, Var
from tqdm import tqdm

import sys
import os
winterp_dir = os.getcwd()+os.sep+os.pardir+os.sep+'python'
sys.path.append(winterp_dir)
import winterp as wi


def tfx(x, a, x0):
    return np.exp(a * np.power(x - x0, 2))


def main():
    # Initial lattice
    Nxyz = (20, 14, 32)
    Dxyz = (1, 1., 1.1)
    xx = np.linspace(0, (Nxyz[0]*Dxyz[0])-Dxyz[0], Nxyz[0])
    yy = np.linspace(0, (Nxyz[1]*Dxyz[1])-Dxyz[1], Nxyz[1])
    zz = np.linspace(0, (Nxyz[2]*Dxyz[2])-Dxyz[2], Nxyz[2])
    xyz = np.meshgrid(xx, yy, zz, indexing='ij', sparse=True)
    ts = np.arange(1)

    # Generate grids for ix, iy, iz
    ix_grid, iy_grid, iz_grid = np.meshgrid(
        np.arange(Nxyz[0]), np.arange(Nxyz[1]), np.arange(Nxyz[2]), indexing='ij')

    # Apply the tfx function to the distance grid to compute inputr and inputc
    ax, ay, az = -0.3, -0.18, -0.15
    # Fill input arrays with some function
    inputr = np.zeros((Nxyz[0], Nxyz[1], Nxyz[2]), dtype=np.float64)
    inputr = tfx(Dxyz[0]*ix_grid, ax, Dxyz[0]*Nxyz[0]/2)*\
             tfx(Dxyz[1]*iy_grid, ay, Dxyz[1]*Nxyz[1]/2)*\
             tfx(Dxyz[2]*iz_grid, az, Dxyz[2]*Nxyz[2]/2)
    
    wdata = WData(
        prefix=f"wdata_example",
        xyz=xyz,  # lattice
        dxyz=Dxyz,
        Nt=ts,  # number of cycles
        dt=1,  # time step
        variables=[
            Var(gauss=np.reshape(np.array(inputr), (len(ts), Nxyz[0], Nxyz[1], Nxyz[2])))
        ],
    )
    print(wdata.get_metadata())
    wdata.save(force=True)
    X_factor=4
    Y_factor=4
    Z_factor=2

    X=Nxyz[0]*X_factor
    Y=Nxyz[1]*Y_factor
    Z=Nxyz[2]*Z_factor

    xx = np.linspace(0, (Nxyz[0]*Dxyz[0])-Dxyz[0]/X_factor, X)
    yy = np.linspace(0, (Nxyz[1]*Dxyz[1])-Dxyz[1]/Y_factor, Y)
    zz = np.linspace(0, (Nxyz[2]*Dxyz[2])-Dxyz[2]/Z_factor, Z)
    xyz = np.meshgrid(xx, yy, zz, indexing='ij', sparse=True)

    interpr = wi.create_interpolator(inputr, Dxyz[0], Dxyz[1], Dxyz[2])
    outr = np.zeros((X,Y,Z), dtype=np.float64)    
    ix=0
    for x in tqdm(xx):
        iy=0
        for y in yy:
            iz=0
            for z in zz:
                outr[ix][iy][iz]=wi.interpolate(interpr, x, y, z)
                iz+=1
            iy+=1
        ix+=1
        
    wdata_interp = WData(
        prefix=f"wdata_example_interp",
        xyz=xyz,  # lattice
        dxyz=Dxyz,
        Nt=ts,  # number of cycles
        dt=1,  # time step
        variables=[
            Var(gauss=np.reshape(np.array(outr), (len(ts), X, Y, Z)))
        ],
    )
    print(wdata_interp.get_metadata())
    wdata_interp.save(force=True)
    
    
if __name__ == "__main__":
    main()
