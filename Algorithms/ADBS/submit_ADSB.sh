#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J exometa_ADSB
### -- ask for number of cores (default: 1) -- 
#BSUB -n 1
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=16GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 16GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 24:00 
### -- set the email address -- 
### BSUB -u ericbaufa10@gmail.com
### -- send notification at completion -- 
### BSUB -N
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o output_exometa_ADSB_%J.out 
#BSUB -e output_exometa_ADSB_%J.err 

### Load modules
module load matlab/R2022a
module load gurobi/9.5.2

# here follow the commands you want to execute with input.in as the input file
matlab -nodisplay -nosplash -nodesktop -r "ADSB_parser('modelName','looplessOption','compactOption'); exit;"