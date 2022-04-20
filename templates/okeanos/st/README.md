# SETTING UP CALCULATIONS AND COMPILATION
(see also: https://gitlab.fizyka.pw.edu.pl/wtools/w-bsk/-/wikis/Setting%20up%20calculations)

The steps below assume that you have set system variable `BSK`.  
You can check if it is set by executing
```bash
echo $BSK
```

## STEP 1: Creating your project folder

```bash
cp -r $BSK/st-project-template/ name-of-your-project
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

## STEP 4: Prepare job submission script & submit

This part strongly depends on the target system. Please, refer to your (super)computer documentation.
