#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N exploratory

set -e

# Set project root path here
cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS-Group-4/

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate phd_r

echo "Starting forest script..."
Rscript 01_scripts/02_analysis/02_univariate/00_forest_plot.R
echo "Forest script done. Exit code: $?"

echo "Starting manhattan script..."
Rscript 01_scripts/02_analysis/02_univariate/01_manhattan_plot.R
echo "Manhattan script done. Exit code: $?"

echo "Starting FAMD script..."
Rscript 01_scripts/02_analysis/03_exploratory/01_FAMD.R
echo "FAMD script done. Exit code: $?"
