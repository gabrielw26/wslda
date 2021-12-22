SETTING UP CALCULATIONS AND COMPILATION 
(see also: http://git2.if.pw.edu.pl/gabrielw/cold-atoms/wikis/setting-up-calculations-and-compilation)

STEP 1: Creating your project folder
    
    cp -r /gs/hs1/hp190063/share/wslda/st-project-template/ name-of-your-project
    cd name-of-your-project
    
    For compilation & computation it is recomended to use workspace
        /gs/hs1/hp190063/
        
STEP 2: Set up compilation environment

    cd name-of-your-project
    source env.sh
    make 1d or make 2d or make 3d

STEP 3: Set up your problem

    You need to modify: 
        predefines.h (do not modify DIAGONALIZATION_ROUTINE)
        problem-definition.h
        logger.h (optionally)
        input.txt
        
    NOTE: You MUST recompile code after each modification of *.h files
        
STEP 4: Prepare job submission script & submit

    Edit job script
        emacs job.sh

    Submit script
        qsub job.sh
                  
