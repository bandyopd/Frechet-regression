#!/bin/bash
#SBATCH --job-name="one_rep"
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p genacc_q
#SBATCH --mail-type="ALL"
#SBATCH -o one_rep.out
#SBATCH -e one_rep.err
module load matlab
matlab -nosplash -nojvm -nodesktop -r one_rep