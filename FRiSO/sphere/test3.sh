#!/bin/bash
#SBATCH --job-name="test3"
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --cores-per-socket=8
#SBATCH -p genacc_q
#SBATCH --mail-type="ALL"
#SBATCH -o test3.out
#SBATCH -e test3.err
module load matlab
matlab -nodisplay -nosplash -r test3