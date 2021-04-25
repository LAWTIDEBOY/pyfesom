#!/bin/bash
#SBATCH -p standard96
#SBATCH --array=1-16

year=$(echo $SLURM_ARRAY_TASK_ID | awk '{print 1999+$1}')  # years from 1985 onwards
echo ${year}