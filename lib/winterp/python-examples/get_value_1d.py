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
    inNX = 32
    inDX = 1.0

    # Fill input arrays with some function
    inputr = np.zeros(inNX, dtype=np.float64)
    inputc = np.zeros(inNX, dtype=np.complex128)

    ax, ay = -0.02, -0.05

    # Compute real and complex values
    for ix in range(inNX):
        inputr[ix] = tfx(inDX * ix, ax, inDX * inNX / 2)
        inputc[ix] = tfx(inDX * ix, ax, inDX * inNX / 2) + 1j * tfx(inDX * ix, ay, inDX * inNX / 2)

    # Create and use interpolator for real data
    interp_r = wi.create_interpolator_1d_r(inNX, inDX, inputr)
    outr = np.array([wi.interpolate_1d_r(interp_r, inDX / 100 * ix) for ix in range(inNX * 100)])
    np.savetxt('get_value_1dr.txt', outr)
    wi.destroy_interpolator(interp_r)

    # Create and use interpolator for complex data
    interp_c = wi.create_interpolator_1d_c(inNX, inDX, inputc)
    outc = np.array([wi.interpolate_1d_c(interp_c, inDX * ix) for ix in range(inNX)])
    np.savetxt('get_value_1dc.txt', outc)
    wi.destroy_interpolator(interp_c)

if __name__ == "__main__":
    main()
