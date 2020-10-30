SETTING UP CALCULATIONS AND COMPILATION 
(see also: http://git2.if.pw.edu.pl/gabrielw/cold-atoms/wikis/setting-up-calculations-and-compilation)

STEP 1: Creating your project folder
    
    cp -r /home2/archive/wslda/td-project-template/ name-of-your-project
    cd name-of-your-project
    
STEP 2: Set up compilation environment

    It is recomended to use one of computation nodes for compilation,
    for example node2066:

    ssh66
    cd name-of-your-project
    source env.sh

STEP 3: Set up your problem

    You need to modify: 
        predefines.h
        problem-definition.h
        logger.h (optionally)
        input.txt
        
STEP 4: Prepare job submission script & submit

    Edit
        job.sh

    Submit
        qsub job.sh
        
    NOTE: qsub can be executed only from login node (2072) 
          or from node2066 (ssh66) 
          
