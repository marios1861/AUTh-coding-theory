#!/bin/bash
#SBATCH -J LDPC-l-10-n-1000
#SBATCH -p batch
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1
#SBATCH -t 1:00:00

#### load module ####
module load matlab/R2022b

#### RUN MATLAB JOB ####
echo "Starting Matlab"
matlab -nodisplay << EOF &> matlab1.out
run("bonus_main_diff_graph.m");
run("bonus_main_same_graph.m");
exit
EOF
