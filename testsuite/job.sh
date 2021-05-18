#!/bin/bash
#PBS -N testsuite
#PBS -l nodes=1:ppn=40
## If you want to submit to specific node use this: 
##PBS -l nodes=node2067.grid4cern.if.pw.edu.pl:ppn=40
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -q long


# execute code
cd $PBS_O_WORKDIR

# execute testsutite 
./testsuiteRun.sh > testsuite.out 2>&1
