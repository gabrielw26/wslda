# wdata-section
*Purpose*: extracts section of selected variable along line specified by two points in space.
```
TOOL FOR GENERATING SECTION ALONG LINE (x1,y1,z1)--(x2,y2,z2) FOR SELECTED VARIABLE
Usage: ./wdata-section -w wtxt -o output -v variable --x1 0 --y1 0 --z1 0 --x2 10 --y2 10 --z2 10 -p 100 -c 0
         -w, --wtxt: wdata descriptor file [REQUIRED]
         -v, --var: variable name [REQUIRED]
         -o, --output: name of output file, default: prefix.var.cycle.txt
         -x, --x1: x-coordinate for starting point, default=0
         -y, --y1: y-coordinate for starting point, default=0, Ignore for 1D data
         -z, --z1: z-coordinate for starting point, default=0, Ignore for 1D and 2D data
         -i, --x2: x-coordinate for final point, default=nx
         -j, --y2: y-coordinate for final point, default=ny, Ignore for 1D data
         -k, --z2: z-coordinate for final point, default=nz, Ignore for 1D and 2D data
         -p, --points: number of sampling points along line (x1,y1,z1)-(x2,y2,z2), default=100
         -c, --cycle: cycle id, default=0
         -h, --help: print help
```

# wdata-cut 
*Purpose*: extracts subset from existing wdata set.  

In many cases, for testing purposes, it is sufficient to have only a sample of data (for example to easy download to local computer and test it). You can use `wdata-cut` tool to extract a data sample. The syntax is following:  
```bash
WDATA SUBSET EXTRACTOR
Usage: ./wdata-cut file.wtxt outprefix start stop
        or
Usage: ./wdata-cut file.wtxt outprefix start stop stride
        file.wtxt    - metadata file
        outprefix    - for new files, metadata file will be written to outprefix.wtxt
        [start,stop) - range for subtructing cycles
        stride       - every `stride` frame will be taken only from given range, default stride=1
```
For example to extract the first cycle only:     
```
[wtools@dell c-examples]$ ../bin/wdata-cut ./test.wtxt sample 0 5 2
# WDATA SUBSET EXTRACTOR
# WORKING DIR: `.` --> `.`
# SUBTRUCTION RANGE: [0,5)
# STRIDE: 2
# READING FILE: `./test.wtxt`
# WRITING NEW METADAFILE: `sample.wtxt`
# DONE.
```

# wdata-interpolate
With this tool, you can interpolate the existing wdata set to a new resolution. Typically this tool is used to generate a figure of better accuracy.
```
[wtools@dell c-examples]$ ../bin/wdata-interpolate test.wtxt high-resolution 48 56 64 
# WDATA SET INTERPOLATOR
# READING INPUT DATA: `test.wtxt`
# **********************  INPUT LATTICE **********************
# LATTICE: 24 x 28 x 32
# SPACING: 1.000000 x 1.000000 x 1.000000
# VOLUME : 24.000000 x 28.000000 x 32.000000
# ********************** OUTPUT LATTICE **********************
# LATTICE: 48 x 56 x 64
# SPACING: 0.500000 x 0.500000 x 0.500000
# VOLUME : 24.000000 x 28.000000 x 32.000000
# WORKING DIR: `.` --> `.`
# INTERPOLATING `density_a`...
# INTERPOLATING `delta`...
# INTERPOLATING `current_a`...
# WRITING `high-resolution.wtxt`
# DONE.
```

# wdata-datadim-up
With this tool, you can increase data dimensionality of an existing data set. The system is assumed to be uniform along added dimension.
```
[wtools@dell c-examples]$ ../bin/wdata-datadim-up test.wtxt testup 3
# WDATA TOOL: DATADIM MODIFIER
# READING INPUT DATA: `test.wtxt`
# ************************ LATTICE ***************************
# LATTICE: 24 x 28 x 32
# SPACING: 1.000000 x 1.000000 x 1.000000
# VOLUME : 24.000000 x 28.000000 x 32.000000
# DIM-IN : 2
# DIM-OUT: 3
# RESIZING `density_a`...
# RESIZING `delta`...
# RESIZING `current_a`...
# WRITING `testup.wtxt`
# DONE.
```
# wdata-merge
With this tool, you can merge a few w-data sets into one.
The current version provides only basic functionality.
```
[gabrielw@wutdell st-testcase-uniform]$ wdata-merge test16-mrg test16 test16-2
# WDATA SETS MERGER
# WORKING WITH: test16.wtxt...
#       ADDING CYCLE 0->0 [time=-1.000000].
#       ADDING CYCLE 1->1 [time=0.000000].
#       ADDING CYCLE 2->2 [time=1.000000].
#       ADDING CYCLE 3->3 [time=2.000000].
#       ADDING CYCLE 4->4 [time=3.000000].
#       ADDING CYCLE 5->5 [time=4.000000].
#       ADDING CYCLE 6->6 [time=5.000000].
#       ADDING CYCLE 7->7 [time=6.000000].
#       ADDING CYCLE 8->8 [time=7.000000].
#       ADDING CYCLE 9->9 [time=8.000000].
#       ADDING CYCLE 10->10 [time=9.000000].
# WORKING WITH: test16-2.wtxt...
#       ADDING CYCLE 5->11 [time=10.000000].
#       ADDING CYCLE 6->12 [time=11.000000].
#       ADDING CYCLE 7->13 [time=12.000000].
#       ADDING CYCLE 8->14 [time=13.000000].
#       ADDING CYCLE 9->15 [time=14.000000].
# WRITING NEW METADAFILE: `test16-mrg.wtxt`
# DONE.
```
