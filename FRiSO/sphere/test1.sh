#!/bin/bash
#SBATCH --job-name="test1"
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --cores-per-socket=8
#SBATCH -p genacc_q
#SBATCH --mail-type="ALL"
#SBATCH -o test1.out
#SBATCH -e test1.err
module load matlab
matlab -nodisplay -nosplash -r test1