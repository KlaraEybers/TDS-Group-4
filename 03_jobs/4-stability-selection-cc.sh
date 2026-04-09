#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -N stabsel_cc
#PBS -o 03_jobs/logs/
#PBS -e 03_jobs/logs/

set -e

export SCENARIO="complete_case"

cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS-Group-4/

LOGFILE="03_jobs/logs/stab_sel_cc_$(date +%Y%m%d_%H%M%S).log"
mkdir -p 03_jobs/logs

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate r413

echo "Starting stability selection: SCENARIO=${SCENARIO}" >> $LOGFILE 2>&1
echo "Started: $(date)" >> $LOGFILE 2>&1

echo "Starting outcome stability selection..." >> $LOGFILE 2>&1
Rscript 01_scripts/02_analysis/03_stability_selection/01_stability_selection_outcome.R >> $LOGFILE 2>&1

echo "Starting mediator stability selection..." >> $LOGFILE 2>&1
Rscript 01_scripts/02_analysis/03_stability_selection/02_stability_selection_mediators1.R >> $LOGFILE 2>&1

echo "Starting incrementation analysis..." >> $LOGFILE 2>&1
Rscript 01_scripts/02_analysis/03_stability_selection/03_incrementation_analysis.R >> $LOGFILE 2>&1

echo "Finished: $(date)" >> $LOGFILE 2>&1