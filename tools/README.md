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
