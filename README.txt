# ============================================================
## PROJECT DESCRIPTION: TDS Group 4
# ============================================================
Research question:
Which external exposome factors are predictive of cardiovascular disease and through
which biological and anthropometric factors are their effects mediated?

Aims:
1) Identify external exposome factors associated with CVD using stability selection LASSO
2) Determine which biological markers and anthropometrics mediate the relationship
   between selected exposome factors and CVD risk
3) Evaluate the predictive performance of an interpretable CVD risk prediction model
   compared to flexible machine learning methods
4) Conduct sub-group analysis by sex to assess whether predictive risk factors differ
   between males and females

# ============================================================
## PROJECT STRUCTURE
# ============================================================

TDS_Group4/
|-- 00_data/
|-- 01_scripts/
|   |-- 00_variable_dictionary.R
|   |-- 01_processing/
|   |   |-- 00_initial_clean.R
|   |   |-- 01_socio_demographics.R
|   |   |-- 02_physical_measures.R
|   |   |-- 03_diet.R
|   |   |-- 04_biological_samples.R
|   |   |-- 05_health_and_medical_history.R
|   |   |-- 06_early_life_factors.R
|   |   |-- 07_environment.R
|   |   |-- 08_lifestyle.R
|   |   |-- 09_mental_health.R
|   |   |-- 10_final_clean.R
|   |   |-- 11_remove_missing.R
|   |   |-- 12_imputation_functions.R
|   |   `-- 12_imputation.R
|   |-- 02_analysis/
|   |   |-- 01_descriptive/
|   |   |   |-- 00_missingness.R
|   |   |   `-- 01_FAMD_and_Correlation.R
|   |   |-- 02_univariate/
|   |   |   |-- 00_forest_plot.R
|   |   |   `-- 01_manhattan_plot.R
|   |   |-- 03_stability_selection/
|   |   |   |-- 01_stability_selection_outcome.R
|   |   |   |-- 02_stability_selection_mediators1.R
|   |   |   `-- 03_incrementation_analysis.R
|   |   |-- 04_multivariate/
|   |   |   |-- mediation_analysis.R
|   |   |   `-- prediction_analysis.R
|   |   |-- 05_neural_network/
|   |   |   `-- cvd_nn.ipynb
|   |   |-- 06_decision_trees/
|   |   |   |-- decision_trees_all_w_shap.ipynb
|   |   |   |-- decision_trees_confounders_selected_w_shap.ipynb
|   |   |   `-- decision_trees_confounders_w_shap.ipynb
|   |   `-- 07_sensitivity_analysis/
|   |       |-- 01_sex_stratified.R
|   |       |-- 02_complete_case.R
|   |       |-- 03_subgroup_performance.R
|   |       `-- 04_anthropometrics_unpenalised.R
|   `-- install_dependencies.R
|-- 02_results/
|-- 03_jobs/
|-- config.py
|-- environment_r413.yml
|-- extraction_and_recoding/
|-- python_requirements.txt
|-- README.txt
|-- tds_group4/
`-- TDS_Group4.Rproj

# ============================================================
## SETUP (complete once before submitting jobs)
# ============================================================
 
1. Install R packages:
 
   eval "$(~/anaconda3/bin/conda shell.bash hook)"
   conda activate r413
   module load NLopt/2.7.1-GCCcore-13.3.0
   Rscript 01_scripts/install_dependencies.R
 
   Note: The NLopt module is required for miceRanger/nloptr to install correctly.
   
2. Transfer the outcomes file (cvd_events.rds) from the real UK Biobank
   to this folder: TDS_Group4/00_data/00_extracted keeping the same naming convention

# ============================================================
## HOW TO RUN JOB SCRIPTS
# ============================================================

All scripts are run via PBS job scripts from the project root. Submit them in 
order, adapting the working directory and environment in each job script to match yours. 

~/TDS_Group4/03_jobs/
  1-pre-processing.sh
  2-imputation.sh
  3-exploratory.sh
 
Default directory:   cd /rds/general/project/hda_25-26/live/TDS/TDS_Group4/
 
Step 1 — Data cleaning and preprocessing:
 
   qsub 03_jobs/01-preprocessing.sh
   
   Instruction: Adjust the ukb_path variable to reference the correct UK Biobank data
 
   Runs all scripts in 01_processing/ (00 through 11) plus missingness analysis.
   Input:  extraction_and_recoding/outputs/ukb_recoded.rds
   Output: 00_data/01_processed/ukb_final_merged.rds
 
Step 2 — Imputation:
 
   qsub 03_jobs/02-imputation.sh
 
   DAG-respecting random forest imputation via miceRanger.
   Stage 1: Exposures imputed using confounders + exposures only.
   Stage 2: Biomarkers imputed using confounders + imputed exposures + biomarkers.
   Output: 00_data/01_processed/03_imputation/
 
Step 3 — Descriptive analysis and univariate:
 
   qsub 03_jobs/03-exploratory.sh
 
   Missingness visualisation, FAMD/correlation, forest plots, Manhattan plots.
   Output: 02_results/02_missingness/
           02_figures/02_univariate/
           02_results/01_tables/00_flow_diagrams/

# ============================================================
## ENVIRONMENT DETAILS
# ============================================================

R version:      4.1.3
Python version: 3.11.13
 
R packages are managed via install_dependencies.R. The full list:
  tidyverse, dplyr, ggplot2, tidyr, purrr, stringr, tibble, readr,
  here, readxl, data.table, miceRanger, future, patchwork, corrplot,
  e1071, FactoMineR, sharp, boot, pROC, broom, scales
 
Python packages are managed via the Python conda environment. The full list:
  numpy, matplotlib, scikit-learn, torch, pyreadr, xgboost, shap, pandas 
 
HPC module required: NLopt/2.7.1-GCCcore-13.3.0 (for miceRanger/nloptr)

# ============================================================
## RESULTS STRUCTURE
# ============================================================

All outputs are written to 02_results/ and organised as follows:

02_results/
|-- 00_figures/
|   |-- 01_imputation/        Convergence diagnostics (numeric + categorical) for
|                             each data split (full, selection, train, test)
|   |-- 02_univariate/        Forest plots and Manhattan plots (raw, BH, Bonferroni)
|   |-- 03_multivariate/      Mediation forest plot, heatmap, ROC curves,
|                             calibration, DAG, odds ratio plots
|   |-- 04_exploratory/       FAMD plots, correlation heatmaps, density/boxplots
|                             by CVD status, incremental AUC analysis
|   |-- 05_sensitivity_sex/   Sex-stratified mediation, prediction, ROC curves
|   |-- 06_sensitivity_complete_cases/  Complete case incremental AUC
|   |-- 07_sensitivity_anthropometrics/ Anthropometric sensitivity incremental AUC
|   `-- overview_panel.png    Summary panel of key figures
|
|-- 01_tables/
|   |-- 00_flow_diagrams/     Population and variable flow numbers
|   |-- 01_multivariate/      AUC, confusion matrix, mediation tables,
|                             model objects (.rds), odds ratios
|   |-- 03_sensitivity_complete_case/  CC mediation effects and AUC
|   |-- 04_neural_network/    Best model weights, hyperparameter search, 
|                             three-model comparison
|   `-- 05_sensitivity_sex/   Sex-stratified mediation and odds ratios
|
|-- 02_missingness/           Variable- and participant-level missingness
|                             figures and summary CSVs
|
`-- 03_ml/                    Machine learning outputs per model configuration:
    |-- 00_neural_network/    NN evaluation plots
    |-- tree_models_python_confounders/           Confounders only
    |-- tree_models_python_all_predictors_.../    All predictors
    `-- tree_models_python_selected_.../          Selected + demographics
                              Each folder contains: metrics, ROC, confusion matrix,
                              feature importance (native, permutation, SHAP)

# ============================================================
## AUTHORS
# ============================================================

Group 4 — Klara Eybers, Hassan El-Maoula, Gong Vitayatanagorn, Shiqi Huang, Luane Branthonne
MSc Health Data Analytics and Machine Learning 
Imperial College London, 2025-26
