#!/bin/bash
#SBATCH --job-name="one_rep2"
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p genacc_q
#SBATCH --mail-type="ALL"
#SBATCH -o one_rep2.out
#SBATCH -e one_rep2.err
module load matlab
matlab -nosplash -nojvm -nodesktop -r "one_rep2; exit"