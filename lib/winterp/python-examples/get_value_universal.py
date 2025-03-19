import numpy as np

import sys
import os
winterp_dir = os.getcwd()+os.sep+os.pardir+os.sep+'python'
sys.path.append(winterp_dir)
import winterp as wi


# Test function in Python

def tfx(x, a, x0):
    return np.exp(a * np.power(x - x0, 2))


def main():
    # Initial lattice
    inNX, inNY, inNZ = 32, 28, 40
    inDX, inDY, inDZ = 0.9, 1.1, 0.8

    # Fill input arrays with some function
    inputr = np.zeros((inNX, inNY, inNZ), dtype=np.float64)
    inputc = np.zeros((inNX, inNY, inNZ), dtype=np.complex128)

    ax, ay = -0.1, -0.08

    # Generate grids for ix, iy, iz
    ix_grid, iy_grid, iz_grid = np.meshgrid(
        np.arange(inNX), np.arange(inNY), np.arange(inNZ), indexing='ij')

    # Compute the distance grid _r using vectorized operations
    _r = np.sqrt((inDX * ix_grid - inDX * inNX / 2) ** 2 +
                 (inDY * iy_grid - inDY * inNY / 2) ** 2 +
                 (inDZ * iz_grid - inDZ * inNZ / 2) ** 2)

    # Apply the tfx function to the distance grid to compute inputr and inputc
    inputr = tfx(_r, ax, 0)
    inputc = inputr + 1j * tfx(_r, ay, 0)  # Gaussian complex

    x_array = np.linspace(0, inDX * inNX, 3200)
    y_array = np.linspace(0, inDY * inNY, 100)
    z_array = np.linspace(0, inDZ * inNZ, 100)

    interpr = wi.create_interpolator(inputr, inDX, inDY, inDZ)
    outr = np.array([wi.interpolate(interpr, 
                                         ix, inDY * inNY / 2, inDZ * inNZ / 2)
                     for ix in x_array])
    np.savetxt('get_value_3dr_uni.txt', outr)

    interpc = wi.create_interpolator(inputc, inDX, inDY, inDZ)
    outc = np.array([wi.interpolate(interpc, 
                                         ix, inDY * inNY / 2, inDZ * inNZ / 2)
                     for ix in x_array])
    np.savetxt('get_value_3dc_uni.txt', outc)

if __name__ == "__main__":
    main()
