SETTING UP CALCULATIONS AND COMPILATION 
(see also: https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Setting%20up%20calculations)

It is recomended to set in your .bashrc

export WSLDA=/mnt/storage_2/project_data/grant_518/wslda
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/storage_2/project_data/grant_518/wslda/lib/wdata:/mnt/storage_2/project_data/grant_518/wslda/lib/wderiv:/mnt/storage_2/project_data/grant_518/wslda/lib/winterp
export PATH=$PATH:/mnt/storage_2/project_data/grant_518/wslda/lib/wdata/bin:/mnt/storage_2/project_data/grant_518/wslda/tools/bin


STEP 1: Creating your project folder
    
    cp -r $WSLDA/st-project-template name-of-your-project
    cd name-of-your-project
    
STEP 2: Set up compilation environment

    For compilation you need to call for interactive session:

    	srun --pty -n 1 /bin/bash

    and within such session you execute:

        source env.sh
        make 1d or make 2d or make 3d
        
    ALTERNATIVELY, you can compile your code using submission:
    
        sbatch job.compile.sh

STEP 3: Set up your problem

    You need to modify: 
        predefines.h (do not modify DIAGONALIZATION_ROUTINE)
        problem-definition.h
        logger.h (optionally)
        input.txt
	
    NOTE: after each modifcation of *.h file you MUST recompile the code
          to apply changes
        
STEP 4: Prepare job submission script & submit

    Edit script 
        nano job.sh

    Submit
        sbatch job.sh
        
    NOTE: sbatch can be executed only from login node
          
