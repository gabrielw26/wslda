# SETTING UP CALCULATIONS AND COMPILATION 

The steps below assume that you have set system variable `WSLDA`.  
You can check if it is set by executing
```bash
echo $WSLDA
```

## STEP 1: Creating your project folder

For calculations use $SCRATCH area

```bash
cd $SCRATCH
cp -r $WSLDA/st-project-template/ name-of-your-project
cd name-of-your-project
```

## STEP 2: Set up compilation environment

To compile you need first allocate the interactive session:

```bash
srun -p plgrid-gpu-a100 -A plginhsf-gpu-a100 -N 1 -n 1 --pty /bin/bash -l
```

Then, within the interactive session:
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
    input.txt
```        
NOTE: You MUST recompile code after each modification of `*.h` files

NOTE: You can use nano editor after loading module:
```bash
module load nano
```
        
## STEP 4: Prepare job submission script & submit

Use provided template `job.sh`
