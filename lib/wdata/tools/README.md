
# wdata-copy
```
TOOL FOR COPYING WDATA SET
Usage: ../bin//wdata-copy -w wtxt -o outprefix ...
	 -w, --wtxt: wdata descriptor file if input dataset [REQUIRED]
	 -o, --outprefix: prefix for output dataset [REQUIRED]
	 -h, --help: print help
```
# wdata-subset
```
TOOL FOR SUBSET EXTRACTING FROM WDATA SET
Usage: ../bin//wdata-subset -w wtxt -o outprefix ...
	 -w, --wtxt: wdata descriptor file if input dataset [REQUIRED]
	 -o, --outprefix: prefix for output dataset [REQUIRED]
	 -b, --begin: first cycle of subset to extraction, default=0
	 -e, --end: last cycle of subset to extraction, default=-1 (the last cycle)
	 -s, --stride: every `stride` frame will be taken only, default=1
	 -h, --help: print help
```
# wdata-section
```
TOOL FOR GENERATING SECTION ALONG LINE (x1,y1,z1)--(x2,y2,z2) FOR SELECTED VARIABLE
Usage: ../bin//wdata-section -w wtxt -o output -v variable --x1 0 --y1 0 --z1 0 --x2 10 --y2 10 --z2 10 -p 100 -c 0
	 -w, --wtxt: wdata descriptor file [REQUIRED]
	 -v, --var: variable name [REQUIRED]
	 -o, --output: name of output file, default: prefix.var.cycle.txt
	 -x, --x1: x-coordinate for starting point, default=0
	 -y, --y1: y-coordinate for starting point, default=0, Ignore for 1D data
	 -z, --z1: z-coordinate for starting point, default=0, Ignore for 1D and 2D data
	 -i, --x2: x-coordinate for final point, default=nx*dx
	 -j, --y2: y-coordinate for final point, default=ny*dy, Ignore for 1D data
	 -k, --z2: z-coordinate for final point, default=nz*dz, Ignore for 1D and 2D data
	 -p, --points: number of sampling points along line (x1,y1,z1)-(x2,y2,z2), default=100
	 -c, --cycle: cycle id, default=0. Negative means take form the end, -1 is the last one.
	 -s, --screen: print on screen the cross-section values
	 -h, --help: print help
```
# wdata-append
```
TOOL FOR APPENDING WDATA SET TO EXISTING ONE
Usage: ../bin//wdata-append -w wtxt -a wtxt -t
	 -w, --wtxt: wdata descriptor file to which new data will be appended  [REQUIRED]
	 -a, --append: wdata descriptor file which will be appended [REQUIRED]
	 -t, --text: update text files, like wlog, stdout
	 -e, --eps: tolerance between least entry times in wdata sets, default=1.0e-6
	 -y, --yes: apply yes answer to prompt
	 -h, --help: print help
```
# wdata-datadim-up
```
WDATA TOOL: DATADIM MODIFIER
Usage: ../bin//wdata-datadim-up -w wtxt -o outprefix ...
	 -w, --wtxt: wdata descriptor file if input dataset [REQUIRED]
	 -o, --outprefix: prefix for output dataset and new metadata file outprefix.wtxt [REQUIRED]
	 -d, --dim: target datadim of output set, default=input datadim+1
	            uniformity along added dimension is assumed.
	 -h, --help: print help
```
# wdata-interpolate
```
WDATA TOOL FOR INTERPOLATING EXISTING DATASET TO NEW RESOLUTION
Usage: ../bin//wdata-interpolate -w wtxt -o outprefix ...
	 -w, --wtxt: wdata descriptor file if input dataset [REQUIRED]
	 -o, --outprefix: prefix for output dataset and new metadata file outprefix.wtxt [REQUIRED]
	 -x, --nx: number of points along x direction for new representation
	           dx will be ajusted such that to keep length along x fixed. Default=input.nx
	 -y, --ny: number of points along y direction for new representation
	           dy will be ajusted such that to keep length along y fixed. Default=input.ny
	 -z, --nz: number of points along z direction for new representation
	           dz will be ajusted such that to keep length along z fixed. Default=input.nz
	 -h, --help: print help
```
