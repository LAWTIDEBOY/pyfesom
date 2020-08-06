#!/bin/bash
#SBATCH -t 06:10:00
#SBATCH --cpus-per-task=2
#SBATCH -p standard96

module load anaconda2
python /home/hbkoziel/pyfesom/codes/export/export_map_mean_and_trend_SCM2.py >> testpy-output.txt && echo "Done with python_file.py"

