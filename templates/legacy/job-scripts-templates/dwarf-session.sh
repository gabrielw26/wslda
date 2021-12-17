#!/bin/bash

if [ $# -ne 2 ]
then
  echo "Usage: `basename $0` nodes ppn"
  exit $E_BADARGS
fi

cmd="qsub -q long -l walltime=12:00:00,nodes=$1:ppn=$2 -I"
echo $cmd
$cmd
