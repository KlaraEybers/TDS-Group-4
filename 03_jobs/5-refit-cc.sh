#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=32:mem=128gb
#PBS -N refit_cc
#PBS -o 03_jobs/logs/
#PBS -e 03_jobs/logs/
set -e

export SCENARIO="complete_case"

cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS-Group-4/

LOGFILE="03_jobs/logs/refit_cc_$(date +%Y%m%d_%H%M%S).log"
mkdir -p 03_jobs/logs

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate phd_r

echo "Starting refit: SCENARIO=${SCENARIO}" >> $LOGFILE 2>&1
echo "Started: $(date)" >> $LOGFILE 2>&1

echo "Starting mediation analysis..." >> $LOGFILE 2>&1
Rscript 01_scripts/02_analysis/04_multivariate/mediation_analysis.R >> $LOGFILE 2>&1

echo "Starting prediction analysis..." >> $LOGFILE 2>&1
Rscript 01_scripts/02_analysis/04_multivariate/prediction_analysis.R >> $LOGFILE 2>&1

echo "Starting sex stratified analysis..." >> $LOGFILE 2>&1
Rscript 01_scripts/02_analysis/04_multivariate/sex_stratified.R >> $LOGFILE 2>&1

echo "Finished: $(date)" >> $LOGFILE 2>&1