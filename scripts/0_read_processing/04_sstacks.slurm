#!/bin/bash -l

#SBATCH --job-name=ersn_sstacks
#SBATCH -o ersn_sstacks.out
#SBATCH -e ersn_sstacks.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 31
#SBATCH --mail-user=melodysyue@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish
#SBATCH --chdir=/scratch-lustre/wlarson/yue/missR/ersn

enable_lmod
module load stacks/2.3

sstacks -P ./stacks_clone_filtered -M ./popmaps/ersn_popmap_all.txt -p 31 --disable_gapped

