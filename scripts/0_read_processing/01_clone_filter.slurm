#!/bin/bash -l

#SBATCH --job-name=ersn_clone_filter
#SBATCH -o ersn_clone_filter.out
#SBATCH -e ersn_clone_filter.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 31
#SBATCH --mail-user=melodysyue@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish
#SBATCH --chdir=/scratch-lustre/wlarson/yue/missR/ersn/process_radtags


enable_lmod
module load stacks/2.3

forwards=$(find . -name '*.1.fq.gz' | grep -v '.rem.1.fq.gz')

for forward in $forwards
do
  reverse="$(basename $forward .1.fq.gz).2.fq.gz"
  clone_filter -1 $forward -2 $reverse -i gzfastq -o ../clone_filtered
done
