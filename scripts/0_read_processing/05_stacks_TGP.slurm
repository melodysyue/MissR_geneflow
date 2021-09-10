#!/bin/bash -l

#SBATCH --job-name=ersn_stacks_TGP
#SBATCH -o ersn_stacks_TGP.out
#SBATCH -e ersn_stacks_TGP.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 31
#SBATCH --mail-user=melodysyue@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish
#SBATCH --chdir=/scratch-lustre/wlarson/yue/missR/ersn

enable_lmod
module load stacks/2.3

##remove the .1 from all the filenames in the directory of stacks. 
##rename .1. . *.tsv.gz

#run tsv2bam
tsv2bam -P ./stacks_clone_filtered -M ./popmaps/ersn_popmap_all_no1.txt -R ./clone_filtered -t 31
#run gstacks
gstacks -P ./stacks_clone_filtered -M ./popmaps/ersn_popmap_all_no1.txt -t 31
#run populations
mkdir ./stacks_clone_filtered/populations 
populations -P ./stacks_clone_filtered/ -M ./popmaps/ersn_popmap_all_no1.txt -R 0.3 -O ./stacks_clone_filtered/populations -t 31 --vcf --genepop


