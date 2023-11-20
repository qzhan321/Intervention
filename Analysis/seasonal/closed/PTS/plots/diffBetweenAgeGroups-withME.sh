#!/bin/bash

#SBATCH --job-name=diffBetweenAgeGroups-withME.sh
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --array=1-2
#SBATCH --ntasks-per-node=1 
#SBATCH --mem-per-cpu=48000
#SBATCH --account=pi-jozik

seasonality=seasonal
openness=closed
state=withME
run=4
R CMD BATCH /home/qizhan/others/PhD/projects/intervention/analysis$run/scripts/actualRuns/$seasonality/$openness/PTS/plots/diffBetweenAgeGroups-$state$SLURM_ARRAY_TASK_ID.R /home/qizhan/others/PhD/projects/intervention/analysis$run/scripts/actualRuns/$seasonality/$openness/PTS/plots/diffBetweenAgeGroups-$state$SLURM_ARRAY_TASK_ID.Rout 