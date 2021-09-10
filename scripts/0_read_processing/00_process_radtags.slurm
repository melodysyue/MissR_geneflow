#!/bin/bash -l

#SBATCH --job-name=shiner_process_radtags
#SBATCH -o shiner_process_radtags.%A.%a.out
#SBATCH -e shiner_process_radtags.%A.%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 32
#SBATCH --array=1-3
#SBATCH --mail-user=yshi@uwsp.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish
#SBATCH --chdir=/scratch-lustre/wlarson/yue/missR/ersn


enable_lmod
module load stacks/2.3


forward=$(find . -name '*_1.fq.gz' | sort | awk "NR==$SLURM_ARRAY_TASK_ID")
reverse=$(find . -name '*_2.fq.gz' | sort | awk "NR==$SLURM_ARRAY_TASK_ID")
barcode=$(find . -name '*.barcode.txt' | sort | awk "NR==$SLURM_ARRAY_TASK_ID")

process_radtags -1 $forward -2 $reverse -b $barcode -e sbfI -i gzfastq -c -q -r --filter_illumina --bestrad -t 140 -o ./process_radtags
