import numpy as np
import os

def read_checkpoint_info_pca(file_name):
    # This function reads checkpoint info from a binary file
    with open(file_name, 'rb') as file:
        prec = int(np.fromfile(file, dtype=np.int32, count=1)[0]) # skip
        nwf = int(np.fromfile(file, dtype=np.int32, count=1)[0])
        nx = int(np.fromfile(file, dtype=np.int32, count=1)[0])
        ny = int(np.fromfile(file, dtype=np.int32, count=1)[0])
        nz = int(np.fromfile(file, dtype=np.int32, count=1)[0])
        dx = np.fromfile(file, dtype=np.float64, count=1)[0]
        dy = np.fromfile(file, dtype=np.float64, count=1)[0]
        dz = np.fromfile(file, dtype=np.float64, count=1)[0]
        kF = np.fromfile(file, dtype=np.float64, count=1)[0]
        mu = np.fromfile(file, dtype=np.float64, count=2)
        ec = np.fromfile(file, dtype=np.float64, count=1)[0]
        beta = np.fromfile(file, dtype=np.float64, count=1)[0]
    return nwf, nx, ny, nz, dx, dy, dz, kF, mu, ec, beta

def process_wave_functions(prefix):
    # Read lattice settings from the info file
    file_name = os.path.join(prefix, 's2dpca.info')
    print(f"# READING INFO FILE: {file_name}")
    nwf, nx, ny, nz, dx, dy, dz, kF, mu, ec, beta = read_checkpoint_info_pca(file_name)

    print(f"# LATTICE: {nx} x {ny} x {nz}")
    print(f"# SPACING: {dx:.2f} x {dy:.2f} x {dz:.2f}")
    print(f"# VOLUME : {dx * nx:.2f} x {dy * ny:.2f} x {dz * nz:.2f}")

    # Allocate arrays for wave-functions
    u = np.zeros((nx, ny), dtype=np.complex128)
    v = np.zeros((nx, ny), dtype=np.complex128)

    N = 0.0

    for ikz in range(nz // 2):
        # Read file header
        file_name = os.path.join(prefix, f's2dpca.{ikz:04d}.info')
        print(f"# OPENING: {file_name}")
        nwf, nx, ny, nz, dx, dy, dz, kF, mu, ec, beta = read_checkpoint_info_pca(file_name)

        # Open files
        file_en = os.path.join(prefix, f's2dpca.{ikz:04d}.en')
        file_kkz = os.path.join(prefix, f's2dpca.{ikz:04d}.kkz')
        file_wfu = os.path.join(prefix, f's2dpca.{ikz:04d}.wfu')
        file_wfv = os.path.join(prefix, f's2dpca.{ikz:04d}.wfv')

        # degeneracy of the state
        wght=1
        if ikz>0: wght=2 #(+kz, -kz)

        try:
            with open(file_en, 'rb') as pFile_en, open(file_kkz, 'rb') as pFile_kkz, \
                 open(file_wfu, 'rb') as pFile_wfu, open(file_wfv, 'rb') as pFile_wfv:

                print(f"# PROCESSING WF[{nwf}] FOR ikz={ikz}")
                for iwf in range(nwf):
                    ei = np.fromfile(pFile_en, dtype=np.float64, count=1)[0]  # eigen energy
                    kz = np.fromfile(pFile_kkz, dtype=np.float64, count=1)[0]  # kz value
                    u = np.fromfile(pFile_wfu, dtype=np.complex128, count=nx*ny).reshape((nx, ny))  # u-component
                    v = np.fromfile(pFile_wfv, dtype=np.complex128, count=nx*ny).reshape((nx, ny))  # v-component

                    # Normalize
                    u /= np.sqrt(dx * dy)
                    v /= np.sqrt(dx * dy)

                    # ...
                    # ... WRITE HERE YOUR CODE ...
                    # ...
                    # ... like particle number (for spin-symmetric case)
                    # N += 2.0*np.sum(np.abs(v)**2)*dx*dy*wght

        except IOError as e:
            print(f"ERROR: Cannot open {e.filename}")
            return

    # Arrays will be cleared automatically
    print("Processing complete.")

# Example usage
prefix = "your_prefix"
process_wave_functions(prefix)
