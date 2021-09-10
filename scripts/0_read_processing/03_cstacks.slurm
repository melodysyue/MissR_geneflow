#!/bin/bash -l

#SBATCH --job-name=ersn_cstacks
#SBATCH -o ersn_cstacks.out
#SBATCH -e ersn_cstacks.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 31
#SBATCH --mail-user=yshi@uwsp.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish
#SBATCH --chdir=/scratch-lustre/wlarson/yue/missR/ersn

enable_lmod
module load stacks/2.3

cstacks -P ./stacks_clone_filtered -M ./popmaps/ERSN_popmap_cstacks.txt -n 3 -p 31 --disable_gapped
