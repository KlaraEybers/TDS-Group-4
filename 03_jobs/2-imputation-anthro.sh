#!/bin/bash
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=128:ompthreads=128:mem=200gb
#PBS -N imputation_anthro
#PBS -o 03_jobs/logs/
#PBS -e 03_jobs/logs/

set -e

export SCENARIO="sensitivity_anthro_confounder"

cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS-Group-4/

LOGFILE="03_jobs/logs/imputation_anthro_$(date +%Y%m%d_%H%M%S).log"
mkdir -p 03_jobs/logs

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate r413

module load NLopt/2.7.1-GCCcore-13.3.0

echo "Starting imputation: SCENARIO=${SCENARIO}" >> $LOGFILE 2>&1
echo "Started: $(date)" >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/12_imputation.R >> $LOGFILE 2>&1
echo "Finished: $(date)" >> $LOGFILE 2>&1