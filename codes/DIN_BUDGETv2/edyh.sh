#!/bin/bash
#SBATCH --account=hbkoziel
#SBATCH -p standard96
#SBATCH -t 11:30:00
#SBATCH --mem=32G
#SBATCH --job-name=EDYH
#SBATCH --output=EDYH_%A_%a.out
#SBATCH --array=1-16


module load anaconda2

year=$(echo $SLURM_ARRAY_TASK_ID | awk '{print 1999+$1}')  # years from 1985 onwards
echo ${year}

srun python2 edyh.py ${year} > edyh_${year}.out