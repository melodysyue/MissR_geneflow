#!/bin/bash -l

#SBATCH --job-name=ersn_ustacks
#SBATCH -o ersn_ustacks.out
#SBATCH -e ersn_ustacks.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 31
#SBATCH --mail-user=melodysyue@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish
#SBATCH --chdir=/scratch-lustre/wlarson/yue/missR/ersn

enable_lmod
module load stacks/2.3
##rename files in clone_filtered
##renmae .1.1. .1. ./clone_filtered/*.fq.gz
##renmae .2.2. .2. ./clone_filtered/*.fq.gz

i=1
samples=$(find ./clone_filtered -name '*.1.fq.gz')


for sample in $samples
do
   ustacks -t gzfastq -f $sample -o ./stacks_clone_filtered -i $i --disable-gapped -m 3 -M 3 -H -p 31 --max_locus_stacks 4 --model_type bounded --bound_high 0.05 
   let "i+=1";
done
