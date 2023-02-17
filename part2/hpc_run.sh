#!/bin/bash
#SBATCH -J iLDPC-n-1000
#SBATCH -p batch
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1
#SBATCH -t 1:00:00

#### load module ####
module load matlab/R2022b

#### RUN MATLAB JOB ####
echo "Starting Matlab"
matlab -nodisplay << EOF &> matlab1.out
run("main_graph.m");
run("main.m");
run("main_naive.m");
run("main_b_graph.m");
exit
EOF
