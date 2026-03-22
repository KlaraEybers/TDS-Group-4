#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=128:ompthreads=128:mem=100gb
#PBS -N imputation

cd /rds/general/project/hda_25-26/live/TDS/TDS_Group4/

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate r413

module load NLopt/2.7.1-GCCcore-13.3.0

Rscript 01_scripts/01_processing/12_imputation.R 

