#!/bin/bash
nohup matlab -nojvm -r 'CrossVal;exit'



#This file can then be executed in the background with one of the three following commands, depending on whether you want to ignore all output, just stdout, or just stderr:

#$ ./run_script.sh > /dev/null 2>&1 &
#$ ./run_script.sh > /dev/null 2> err.txt &
#$ ./run_script.sh > log.txt 2> /dev/null &

