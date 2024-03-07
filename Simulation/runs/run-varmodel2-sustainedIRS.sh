#!/bin/bash
#SBATCH --job-name=Oct24th
#SBATCH --time=36:00:00
#SBATCH --output=/home/qizhan/others/PhD/projects/intervention/simulation4/scripts/actualRuns/seasonal/closed/outputAndErrors/sustainedIRS/Oct24th_%A_%a.out
#SBATCH --error=/home/qizhan/others/PhD/projects/intervention/simulation4/scripts/actualRuns/seasonal/closed/outputAndErrors/sustainedIRS/Oct24th_%A_%a.err
#SBATCH --array=101-114
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=68000
#SBATCH --partition=caslake
#SBATCH --account=pi-jozik
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qizhan@uchicago.edu

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
module load gcc/7.4.0 # gcc/6.1
module load python # python/cpython-3.7.0
module load R
# SLURM_ARRAY_TASK_ID=4
filePrefix=sim
seasonality=seasonal
openness=closed
code=varmodel2-master
runCategory=actualRuns
runPrepDir=runPrepFiles
modelDir=runModelsDir
sqliteDir=sqlitesDir
remoteDir=/scratch/midway2/qizhan/PhD/projects/intervention/simulation4/$runCategory/$seasonality/$openness/${filePrefix}_${SLURM_ARRAY_TASK_ID}
cp -r /home/qizhan/others/PhD/projects/intervention/simulation4/codes/$code $remoteDir
cp -r /home/qizhan/others/PhD/projects/intervention/simulation4/scripts/$runCategory/$seasonality/$openness/runInputFiles $remoteDir/$runPrepDir
cd $remoteDir
cd $runPrepDir
python writeParameters.py -p ${filePrefix}_param_sustainedIRS.csv -i parameters-template.py -n $SLURM_ARRAY_TASK_ID -r 1 -x $filePrefix -s
cd ..

mkdir $modelDir
mkdir $sqliteDir
# build the model, run preIRS first
./build.py -p $runPrepDir/${filePrefix}_${SLURM_ARRAY_TASK_ID}_r0_input.py -d $modelDir/s0
# execute the run
cd $sqliteDir
../$modelDir/s0/bin/varMig
cd ..
