#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -N ml_models
#PBS -o 03_jobs/logs/
#PBS -e 03_jobs/logs/
set -e

cd /rds/general/project/hda_25-26/live/TDS/TDS_Group4/ # ADAPT THIS 

LOGFILE="03_jobs/logs/ml_models_$(date +%Y%m%d_%H%M%S).log"
mkdir -p 03_jobs/logs

source /rds/general/project/hda_25-26/live/TDS/TDS_Group4/group4_python_environment/bin/activate # ADAPT THIS

echo "Started: $(date)" >> $LOGFILE 2>&1

echo "Starting neural network..." >> $LOGFILE 2>&1
jupyter nbconvert --to notebook --execute 01_scripts/02_analysis/05_neural_network/cvd_nn.ipynb >> $LOGFILE 2>&1
echo "Neural network done." >> $LOGFILE 2>&1

echo "Starting decision trees (confounders)..." >> $LOGFILE 2>&1
jupyter nbconvert --to notebook --execute 01_scripts/02_analysis/06_decision_trees/decision_trees_confounders_w_shap.ipynb >> $LOGFILE 2>&1
echo "Decision trees (confounders) done." >> $LOGFILE 2>&1

echo "Starting decision trees (confounders selected)..." >> $LOGFILE 2>&1
jupyter nbconvert --to notebook --execute 01_scripts/02_analysis/06_decision_trees/decision_trees_confounders_selected_w_shap.ipynb >> $LOGFILE 2>&1
echo "Decision trees (confounders selected) done." >> $LOGFILE 2>&1

echo "Starting decision trees (all)..." >> $LOGFILE 2>&1
jupyter nbconvert --to notebook --execute 01_scripts/02_analysis/06_decision_trees/decision_trees_all_w_shap.ipynb >> $LOGFILE 2>&1
echo "Decision trees (all) done." >> $LOGFILE 2>&1

echo "Starting combined ROC curves..." >> $LOGFILE 2>&1
jupyter nbconvert --to notebook --execute 01_scripts/02_analysis/08_combined_roc_curves/combined_roc_curves.ipynb >> $LOGFILE 2>&1
echo "Combined ROC curves done." >> $LOGFILE 2>&1