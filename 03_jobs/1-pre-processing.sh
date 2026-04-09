#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=4:mem=50GB
#PBS -N pre-processing

set -e

# Set project root path here
cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS-Group-4/

# Set UK Biobank path here
ukb_path=/rds/general/project/hda_25-26/live/TDS/fg520/TDS-Group-4/extraction_and_recoding/outputs/ukb_recoded.rds

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate phd_r

echo "Starting initial cleaning script..."
Rscript 01_scripts/01_processing/00_initial_clean.R $ukb_path
echo "Initial cleaning script done. Exit code: $?"

echo "Starting sociodemographic script..."
Rscript 01_scripts/01_processing/01_socio_demographics.R
echo "Sociodemographic script done. Exit code: $?"

echo "Starting physical measures script..."
Rscript 01_scripts/01_processing/02_physical_measures.R
echo "Physical measures script done. Exit code: $?"

echo "Starting diet script..."
Rscript 01_scripts/01_processing/03_diet.R
echo "Diet script done. Exit code: $?"

echo "Starting biological samples script..."
Rscript 01_scripts/01_processing/04_biological_samples.R
echo "Biological samples script done. Exit code: $?"

echo "Starting medical history script..."
Rscript 01_scripts/01_processing/05_health_and_medical_history.R
echo "Medical history script done. Exit code: $?"

echo "Starting early life script..."
Rscript 01_scripts/01_processing/06_early_life_factors.R
echo "Early life script done. Exit code: $?"

echo "Starting environment script..."
Rscript 01_scripts/01_processing/07_environment.R
echo "Environment script done. Exit code: $?"

echo "Starting lifestyle script..."
Rscript 01_scripts/01_processing/08_lifestyle.R
echo "Lifestyle script done. Exit code: $?"

echo "Starting mental health script..."
Rscript 01_scripts/01_processing/09_mental_health.R
echo "Mental health script done. Exit code: $?"

echo "Starting final cleaning script..."
Rscript 01_scripts/01_processing/10_final_clean.R
echo "Final cleaning script done. Exit code: $?"

echo "Starting missingness script..."
Rscript 01_scripts/01_processing/11_remove_missing.R
echo "Missingness script done. Exit code: $?"

echo "Starting missingness visualisation script..."
# Rscript 01_scripts/03_report/00_missingness/missingness_final.R
echo "Missingness visualisation script done. Exit code: $?"
