#!/bin/bash
#SBATCH --account=hbkoziel
#SBATCH -p standard96
#SBATCH -t 11:30:00
#SBATCH --mem=32G
#SBATCH --job-name=EDYV
#SBATCH --output=EDYV_%A_%a.out
#SBATCH --array=1-16

module load anaconda2

year=$(echo $SLURM_ARRAY_TASK_ID | awk '{print 1999+$1}')  # years from 1985 onwards
echo ${year}

srun python2 edyv.py ${year} > edyv_${year}.out
#srun python2 EDYH2.py > output.out