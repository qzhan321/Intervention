#!/bin/bash

#SBATCH --job-name=MOIInd.sh
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --array=1,4,5
#SBATCH --ntasks-per-node=1 
#SBATCH --mem-per-cpu=1500
#SBATCH --account=pi-pascualmm

prefix=survey
readDir0=/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-290923_NYU/MOI/MOIEstInputs
saveDir0=/project2/pascualmm/QZ/PhD/projects/intervention/natComRevision/files/figures/main/Fig5-290923_NYU/MOI/MOIEst
if [ -d "$saveDir0" ]; then
  echo "folder exists!"
else
  mkdir $saveDir0
fi

MOIEstDir=/home/qizhan/others/PhD/projects/intervention/natComRevision/analysis/files/Fig5-290923_NYU
repertoireSizeDistDir=/project2/pascualmm/QZ/PhD/projects/FOI/eLifeSubMarch2024/round2/files/exponential/FOI/inputs/MOIInd/repertoireSizeDist

Rscript $MOIEstDir/MOI_estimation.R -i $readDir0/${prefix}_${SLURM_ARRAY_TASK_ID}.csv -m 20 -r $repertoireSizeDistDir/repertoireSizeDistribution.csv -t "count" -p "uniform" -s "medium" -v TRUE -a "pool" -o $saveDir0/${prefix}_${SLURM_ARRAY_TASK_ID}.RData
