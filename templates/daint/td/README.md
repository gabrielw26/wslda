
It is recomended to set in your .bashrc
export WSLDA=/project/pr125/share/wslda/


STEP 1: Creating your project folder
    
    cp -r /project/pr125/share/wslda/st-project-template name-of-your-project
    cd name-of-your-project
    
    For computation you must use $SCRATCH folder
        cd $SCRATCH
   NOTE: files older than 30days are automatically removed from  $SCRATCH
   For storing results use location:
        cd /project/pr125/
    
STEP 2: Set up compilation environment

    source env.sh
	make 1d or make 2d or make 3d

STEP 3: Set up your problem

    You need to modify: 
        predefines.h 
        problem-definition.h
        logger.h (optionally)
        input.txt
	
    NOTE: after each modifcation of *.h file you MUST recompile the code
          to apply changes
        
STEP 4: Prepare job submission script & submit

    Edit script 
        emacs job.sh

    Submit
        sbatch job.sh
        
    NOTE: sbatch can be executed only from login node
          
