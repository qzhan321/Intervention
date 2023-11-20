#!/bin/bash

#SBATCH --job-name=epi-vs-time.sh
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --array=1
#SBATCH --ntasks-per-node=1 
#SBATCH --mem-per-cpu=38000
#SBATCH --account=pi-jozik

seasonality=seasonal
openness=closed
run=4
R CMD BATCH /home/qizhan/others/PhD/projects/intervention/analysis$run/scripts/actualRuns/$seasonality/$openness/epi/epi-vs-time.R /home/qizhan/others/PhD/projects/intervention/analysis$run/scripts/actualRuns/$seasonality/$openness/epi/epi-vs-time.Rout 