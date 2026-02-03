import numpy as np
from wdata.io import WData, Var
from tqdm import tqdm
import matplotlib.pyplot as plt

import sys
import os

def extract_nodal_lines_tile(delta, tileN=10, refineGrid=5, minPoints=5, minPointsDense=0.75, phaseJump=np.pi, progress=1, debug=0):
    """
    Extract nodal lines from the delta data using a tile-based approach.
    For each tile in the x-y plane, sample points along the border of the tile,
    compute the phase difference of the delta values at these points, and determine
    if a nodal line passes through the tile based on the phase differences.
    Parameters:
    - delta: 2D numpy array of complex delta values.
    - tileN: Number of points to sample along each side of the tile.
    - refineGrid: Factor by which to refine the grid for interpolation before processing.
    - phaseJump: Phase jump to consider for nodal line detection (default is pi).
    - minPoints: Minimum number of points with phase difference close to pi to consider a nodal line.
    - minPointsDense: Minimum density of contiguous points with phase difference close to pi.
    - progress: If 1, print progress to the console.
    - debug: If 1, save debug plots for tiles that are close to the criteria.
    Returns:
    - nodal_data: 3D numpy array with shape (nx, ny, 6) containing nodal line information.
      If a nodal line is detected in a tile, then nodal_data[ix, iy, 0]>0 and the following entries contain:
        - nodal_data[ix, iy, 0]: Number of points on the nodal line.
        - nodal_data[ix, iy, 1]: x coordinate of the start of the nodal line segment.
        - nodal_data[ix, iy, 2]: y coordinate of the start of the nodal line segment.
        - nodal_data[ix, iy, 3]: x coordinate of the end of the nodal line segment.
        - nodal_data[ix, iy, 4]: y coordinate of the end of the nodal line segment.
        - nodal_data[ix, iy, 5]: Distance between the start and end points of the nodal line segment.
    """
    # check dimensions of delta
    dim=delta.shape
    if(len(dim)!=2):
        raise ValueError("delta must be a 2D array")
    
    # refine delta for better interpolation
    refine_dim = np.array(dim)*refineGrid

    # check consistency of refined dimensions
    # tileN must divide refine_dim
    if(refine_dim[0]%tileN!=0 or refine_dim[1]%tileN!=0):
        raise ValueError("refineGrid*dim must be divisible by tileN")
    dx=dim[0]/refine_dim[0]
    dy=dim[1]/refine_dim[1]

    if progress==1:
        print(f"\rInterpolating...", end='', flush=True)
    rdelta = wi.interpolation_2d_c(delta, refine_dim[0], refine_dim[1])
    #check new dimensions
    rdim = rdelta.shape
    if(rdim[0]!=refine_dim[0] or rdim[1]!=refine_dim[1]):
        raise ValueError("Refined delta has incorrect dimensions: Interpolation failed!")

    # buffer for results
    nodal_data = np.zeros((rdim[0]//tileN, rdim[1]//tileN, 6)) # [num_points_on_nodal_line, x1 y1, x2, y2, distance]

    # for each point in x-y plane
    for ix in range(0,rdim[0],tileN):
        if progress==1:
            print(f"\rProgress: {100*ix/rdim[0]:3.0f}%", end='', flush=True)

        for iy in range(0,rdim[1],tileN):
            x = ix*dx
            y = iy*dy

            # compute coodinates of points in tile with corners [x,y], [x+tileN*dx,y], [x,y+tileN*dy], [x+tileN*dx,y+tileN*dy] with tileN points per side
            tile_points = []
            tile_values = []
            tile_length = []
            tile_phasediff = []
            # border from [x,y] to [x+tileN*dx,y]
            lgth=0.0
            for i in range(tileN):
                tx = ix+i
                ty = iy
                tile_points.append((tx*dx, ty*dy))
                tile_values.append(rdelta[tx%rdim[0], ty%rdim[1]])
                tile_phasediff.append(np.angle(tile_values[-1]/tile_values[0]))
                lgth += dx
                tile_length.append(lgth)
            # border from [x+tileN*dx,y] to [x+tileN*dx,y+tileN*dy]
            for i in range(tileN):
                tx = ix + tileN
                ty = iy + i
                tile_points.append((tx*dx, ty*dy))
                tile_values.append(rdelta[tx%rdim[0], ty%rdim[1]])
                tile_phasediff.append(np.angle(tile_values[-1]/tile_values[0]))
                lgth += dy
                tile_length.append(lgth)
            # border from [x+tileN*dx,y+tileN*dy] to [x,y+tileN*dy]
            for i in range(tileN):
                tx = ix + tileN - i
                ty = iy + tileN
                tile_points.append((tx*dx, ty*dy))
                tile_values.append(rdelta[tx%rdim[0], ty%rdim[1]])
                tile_phasediff.append(np.angle(tile_values[-1]/tile_values[0]))
                lgth += dx
                tile_length.append(lgth)
            # border from [x,y+tileN*dy] to [x,y]
            for i in range(tileN):
                tx = ix 
                ty = iy + tileN - i 
                tile_points.append((tx*dx, ty*dy))
                tile_values.append(rdelta[tx%rdim[0], ty%rdim[1]])
                tile_phasediff.append(np.angle(tile_values[-1]/tile_values[0]))
                lgth += dy
                tile_length.append(lgth)

            tile_phasediff = np.abs(tile_phasediff) # take absolute value

            # compute max phase difference around the tile
            max_phasediff = np.max(tile_phasediff) - np.min(tile_phasediff)

            # compute optimal chi square value
            chi_square = 0.0
            tile_phasediff_01=[]
            for pd in tile_phasediff:
                t1,t2 = (pd-0.0)**2, (pd-phaseJump)**2
                if t1<t2:
                    tile_phasediff_01.append(0)
                    chi_square += t1
                else:
                    tile_phasediff_01.append(1)
                    chi_square += t2

            # convert to numpy array
            tile_phasediff_01 = np.array(tile_phasediff_01)   

            # decide if we have nodal line
            p1 = np.sum(tile_phasediff_01)
            if p1<minPoints:
                if debug==1 and p1>0:
                    # for testing
                    debug_file_png="debug_phassdiff.%03d_%03d.png"%(ix, iy)
                    plt.figure()
                    plt.plot(tile_length, tile_phasediff, 'o-')
                    plt.plot(tile_length, np.array(tile_phasediff_01)*phaseJump, 'o-')
                    plt.xlabel('Length')
                    plt.ylabel('Phase Difference')
                    plt.title('ix={}, iy={}: p1<minPoints'.format(ix, iy))
                    plt.grid()
                    plt.savefig(debug_file_png)
                    plt.close()
                continue # no nodal line
            # find min and max index where tile_phasediff_01 is one
            indices = np.where(tile_phasediff_01==1)[0]
            min_index = np.min(indices)
            max_index = np.max(indices)
            p2=max_index - min_index + 1
            if p1/p2 < minPointsDense:
                if debug==1:
                    # for testing
                    debug_file_png="debug_phassdiff.%03d_%03d.png"%(ix, iy)
                    plt.figure()
                    plt.plot(tile_length, tile_phasediff, 'o-')
                    plt.plot(tile_length, np.array(tile_phasediff_01)*phaseJump, 'o-')
                    plt.xlabel('Length')
                    plt.ylabel('Phase Difference')
                    plt.title('ix={}, iy={}: p1/p2 < minPointsDense'.format(ix, iy))
                    plt.grid()
                    plt.savefig(debug_file_png)
                    plt.close()
                continue # not enough contiguous points

            # store coordinates of nodal line point
            nix=ix//tileN
            niy=iy//tileN
            nodal_data[nix, niy, 0] = p1 # store number of points on nodal line
            nodal_data[nix, niy, 1] = tile_points[min_index][0] # x1 coord
            nodal_data[nix, niy, 2] = tile_points[min_index][1] # y1 coord
            nodal_data[nix, niy, 3] = tile_points[max_index][0] # x2 coord
            nodal_data[nix, niy, 4] = tile_points[max_index][1] # y2 coord
            # compute distance between points
            _dx = tile_points[max_index][0] - tile_points[min_index][0]
            _dy = tile_points[max_index][1] - tile_points[min_index][1]
            dist = np.sqrt(_dx*_dx + _dy*_dy)
            nodal_data[nix, niy, 5] = dist # distance between min and max point

    if progress==1:
        print("\rProgress: 100%")

    # return result
    return nodal_data

def compute_total_length_of_nodal_lines(nodal_data):
    """
    Compute the total length of nodal lines from the nodal data.
    Parameters:
    - nodal_data: 3D numpy array with shape (nx, ny, 7) containing nodal line information.
    Returns:
    - total_length: Total length of nodal lines.
    """
    total_length = 0.0
    dim = nodal_data.shape
    for ix in range(dim[0]):
        for iy in range(dim[1]):
            if nodal_data[ix, iy, 0] > 0:
                total_length += nodal_data[ix, iy, 5]
    return total_length

def make_plot_of_nodal_lines(delta,nodal_data, vmin=None, vmax=None, title=None, file_name=None):
    """
    Make a plot of the nodal lines from the nodal data.
    Parameters:
    - nodal_data: 3D numpy array with shape (nx, ny, 7) containing nodal line information.
    """
    ratio=delta.shape[0]/delta.shape[1]
    scale=5 # adjust scale as needed
    fig, ax = plt.subplots(figsize=(scale*ratio,1.2*scale))

    data_to_show=np.abs(delta)
    im = ax.imshow(
        np.transpose(data_to_show),
        interpolation='bilinear',
        extent=[0, delta.shape[0], 0, delta.shape[1]],
        aspect='equal',
        origin='lower',
        vmin=vmin if vmin is not None else np.nanmin(data_to_show),
        vmax=vmax if vmax is not None else np.nanmax(data_to_show),
    )
    # fig.colorbar(im, ax=ax, orientation="horizontal")

    for ix in range(nodal_data.shape[0]):
        for iy in range(nodal_data.shape[1]):
            npoints = nodal_data[ix, iy, 0]
            if npoints>0:
                x1 = nodal_data[ix, iy, 1]
                y1 = nodal_data[ix, iy, 2]
                x2 = nodal_data[ix, iy, 3]
                y2 = nodal_data[ix, iy, 4]
                ax.plot([x1, x2], [y1, y2], color='red', linewidth=1.0)

    # add title if provided
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title("Nodal Lines")

    plt.tight_layout()

    if file_name is not None:
        fig.savefig(file_name)
    else:
        plt.show()

    # close figure to free memory
    plt.close(fig)


# Load winterp python bindings
WSLDA_DIR = os.getenv('WSLDA') # get WSLDA directory from environment variable
WINTERP_DIR = WSLDA_DIR + "/lib/winterp" # default location of winterp
# WINTERP_DIR = "/home/gabrielw/MyProjects/winterp" # custom location for my system
WINTERP_DIR_PYTHON = WINTERP_DIR + "/python" # default location of winterp python bindings
sys.path.append(WINTERP_DIR_PYTHON)
import winterp as wi

# main function
if __name__ == "__main__":
    # set input and output files
    fwtxt="/home/gabrielw/MyProjects/winterp/nodal-lines/data/60_aslda_P35-v1x_v0y_TD_2.wtxt"
    outprefix="./result/v1.0"

    # load data
    data = WData.load(fwtxt)
    cycles = data.Nt
    print(f"Data has {cycles} time steps.")

    # arrays to store time and lengths
    time=[]
    lengths=[]

    # create empty txt file with header
    fout = open(outprefix + "_total_length_vs_time.txt", "w")
    fout.write("# Time\tTotal_Length_of_Nodal_Lines\n")
    fout.flush()

    # take from first time step min and max for plotting
    delta=data.delta[0] # sample data
    dmin=np.min(np.abs(delta))
    dmax=np.max(np.abs(delta))

    # loop over time steps
    for tstep in range(cycles):
        print(f"Processing time step {tstep+1}/{cycles}...")
        delta=data.delta[tstep] # sample data
        
        # if file f"./data/nodal_lines_tile_t{tstep:03d}.npy exists, skip
        if os.path.exists(outprefix + f"_nodal_lines_tile_t{tstep:03d}.npy"):
            print(f"Time step {tstep+1}/{cycles} already processed, skipping...")
            nodal_data = np.load(outprefix + f"_nodal_lines_tile_t{tstep:03d}.npy")
        else:
            nodal_data = extract_nodal_lines_tile(delta, tileN=15, refineGrid=5, minPoints=5, minPointsDense=0.75, phaseJump=np.pi*0.75)
            # save to npy file for later use
            np.save(outprefix + f"_nodal_lines_tile_t{tstep:03d}.npy", nodal_data)

        # continue work with nodal_data
        total_length = compute_total_length_of_nodal_lines(nodal_data)
        print(f"Total length of nodal lines: {total_length}")
        make_plot_of_nodal_lines(delta, nodal_data, vmin=dmin, vmax=dmax, title=f"Time={data.t0 + tstep*data.dt:8.2f}", file_name=outprefix + f"_nodal_lines_tile_t{tstep:03d}.png")

        time.append(data.t0 + tstep*data.dt)
        lengths.append(total_length)

        # append to txt file
        fout.write(f"{data.t0 + tstep*data.dt:8.2f}\t{total_length:12.6f}\n")
        fout.flush()

    # close txt file
    fout.close()

    # make plot of total length vs time
    plt.figure()
    plt.plot(time, lengths, marker='o', label='length from extractor', alpha=0.5)
    # compute moving average for smoothing
    window_size = 10
    lengths_smooth = np.convolve(lengths, np.ones(window_size)/window_size, mode='valid')
    plt.plot(time[window_size-1:], lengths_smooth, color='red', linewidth=2, label='smoothed')
    plt.xlabel("Time")
    plt.ylabel("Total Length of Nodal Lines")
    plt.title("Total Length of Nodal Lines vs Time")
    plt.grid()
    plt.legend()
    plt.savefig(outprefix + "_total_length_vs_time.png")
    