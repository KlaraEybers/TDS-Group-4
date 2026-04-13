#!/bin/bash
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=4:mem=50gb
#PBS -N pre-processing
#PBS -o 03_jobs/logs/
#PBS -e 03_jobs/logs/

set -e

# Set project root path here
cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS-Group-4/

LOGFILE="03_jobs/logs/processing_$(date +%Y%m%d_%H%M%S).log"
mkdir -p 03_jobs/logs

# Set UK Biobank path here
ukb_path=/rds/general/project/hda_25-26/live/TDS/fg520/TDS-Group-4//extraction_and_recoding/outputs/ukb_recoded.rds
# Set CVD outcome path here
cvd_events_path=/rds/general/project/hda_25-26/live/TDS/fg520/TDS-Group-4/00_data/00_extracted/cvd_events_real.rds

#Set missingness thresholds
variable_threshold=0.3
participant_threshold=0.3


eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate phd_r

echo "Starting initial cleaning script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/00_initial_clean.R $ukb_path >> $LOGFILE 2>&1
echo "Initial cleaning script done." >> $LOGFILE 2>&1

echo "Starting sociodemographic script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/01_socio_demographics.R >> $LOGFILE 2>&1
echo "Sociodemographic script done." >> $LOGFILE 2>&1

echo "Starting physical measures script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/02_physical_measures.R >> $LOGFILE 2>&1
echo "Physical measures script done." >> $LOGFILE 2>&1

echo "Starting diet script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/03_diet.R >> $LOGFILE 2>&1
echo "Diet script done. Exit code:" >> $LOGFILE 2>&1

echo "Starting biological samples script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/04_biological_samples.R >> $LOGFILE 2>&1
echo "Biological samples script done." >> $LOGFILE 2>&1

echo "Starting medical history script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/05_health_and_medical_history.R >> $LOGFILE 2>&1
echo "Medical history script done." >> $LOGFILE 2>&1

echo "Starting early life script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/06_early_life_factors.R >> $LOGFILE 2>&1
echo "Early life script done." >> $LOGFILE 2>&1

echo "Starting environment script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/07_environment.R >> $LOGFILE 2>&1
echo "Environment script done." >> $LOGFILE 2>&1

echo "Starting lifestyle script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/08_lifestyle.R >> $LOGFILE 2>&1
echo "Lifestyle script done." >> $LOGFILE 2>&1

echo "Starting mental health script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/09_mental_health.R >> $LOGFILE 2>&1
echo "Mental health script done." >> $LOGFILE 2>&1

echo "Starting final cleaning script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/10_final_clean.R $cvd_events_path >> $LOGFILE 2>&1
echo "Final cleaning script done." >> $LOGFILE 2>&1

echo "Starting missingness script..." >> $LOGFILE 2>&1
Rscript 01_scripts/01_processing/11_remove_missing.R $variable_threshold $participant_threshold >> $LOGFILE 2>&1
echo "Missingness script done." >> $LOGFILE 2>&1

