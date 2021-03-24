#!/bin/bash
#SBATCH -p standard96
#SBATCH -t 11:30:00
#SBATCH --mem=256G
#SBATCH --job-name=EDYHv2
#SBATCH --output=EDYHv2_%A_%a.out
#SBATCH --array=1-16
#SBATCH --cpus-per-task=16

module load anaconda2

year=$(echo $SLURM_ARRAY_TASK_ID | awk '{print 1999+$1}')  # years from 1985 onwards
echo ${year}

srun python2 edyhv2.py ${year} > edyhv2_${year}.out