#!/bin/bash
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N exploratory
#PBS -o 03_jobs/logs/
#PBS -e 03_jobs/logs/

set -e

# Set project root path here
cd /rds/general/project/hda_25-26/live/TDS/TDS_Group4/

LOGFILE="03_jobs/logs/exploratory_$(date +%Y%m%d_%H%M%S).log"
mkdir -p 03_jobs/logs

# Set CVD outcome path here
cvd_events_path=/rds/general/project/hda_25-26/live/TDS/TDS_Group4/00_data/00_extracted/cvd_events.rds

#Set missingness thresholds
variable_threshold=0.3
participant_threshold=0.3


eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate r413

echo "Starting missingness visualisation script..." >> $LOGFILE 2>&1
Rscript 01_scripts/02_analysis/01_descriptive/00_missingness.R $variable_threshold $participant_threshold >> $LOGFILE 2>&1
echo "Missingness visualisation script done." >> $LOGFILE 2>&1

echo "Starting FAMD script..." >> $LOGFILE 2>&1
Rscript 01_scripts/02_analysis/01_descriptive/01_FAMD_and_Correlation.R >> $LOGFILE 2>&1 
echo "FAMD script done." >> $LOGFILE 2>&1

echo "Starting forest script..." >> $LOGFILE 2>&1
Rscript 01_scripts/02_analysis/02_univariate/00_forest_plot.R >> $LOGFILE 2>&1
echo "Forest script done." >> $LOGFILE 2>&1

echo "Starting manhattan script..." >> $LOGFILE 2>&1
Rscript 01_scripts/02_analysis/02_univariate/01_manhattan_plot.R >> $LOGFILE 2>&1
echo "Manhattan script done." >> $LOGFILE 2>&1

echo "Starting population flow diagram script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/14_population_flow_diagram.R $cvd_events_path >> $LOGFILE 2>&1
echo "Population flow diagram script done." >> $LOGFILE 2>&1

echo "Starting variable flow diagram script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/15_variable_flow_diagram.R >> $LOGFILE 2>&1
echo "Variable flow diagram script done." >> $LOGFILE 2>&1

echo "Starting variable flow diagram script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/16_full_table_1.R >> $LOGFILE 2>&1
echo "Variable flow diagram script done." >> $LOGFILE 2>&1

echo "Starting imputed table 1 script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/17_imputed_cohort_table_1.R >> $LOGFILE 2>&1
echo "Full imputed table 1 script done." >> $LOGFILE 2>&1



