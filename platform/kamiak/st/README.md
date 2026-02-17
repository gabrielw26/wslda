# INFORMATION
Wiki pages:
https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/

The general Kamiak documentation is here:
https://hpc.wsu.edu/

CoCalc notes:
https://cocalc.com/projects/13fd6946-fdd9-474e-a824-5ba2b0d7efda/files/Computers/Kamiak-Notes.md

Discourse:
https://discourse.iscimath.org/t/using-the-kamiak-cluster/900

WSU CIRC:
https://forum.circ.wsu.edu/

# SETTING UP CALCULATIONS AND COMPILATION
The steps below assume that you have set system variable `WSLDA`.  
You can check if it is set by executing
```bash
echo $WSLDA
```

## STEP 1: Creating your project folder

```bash
cp -r $WSLDA/st-project-template/ name-of-your-project
cd name-of-your-project
```

## STEP 2: Set up compilation environment

```bash
source env.sh

make 1d 
# or 
make 2d 
# or 
make 3d
# or
make
# to make all codes
```

## STEP 3: Set up your problem

You need to modify:
```
    predefines.h
    problem-definition.h
    logger.h (optionally)
    machine.h (optionally)
    input.txt
```        
NOTE: You MUST recompile code after each modification of `*.h` files

NOTE: If you want to use GPUs for diagonalizations uncomment in `machine.h` line
```
#define ELPA_USE_GPU_NVIDIA
```
        
## STEP 4: Prepare job submission script & submit

This part strongly depends on the target system. Please, refer to your (super)computer documentation.
