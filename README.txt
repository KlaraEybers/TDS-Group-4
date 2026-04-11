# ============================================================
## PROJECT DESCRIPTION: TDS Group 4
# ============================================================
Research question:
Which external exposome factors are predictive of cardiovascular disease and through
which biomarkers are their effects mediated?

Aims:
1) Identify external exposome factors associated with CVD using stability selection LASSO
2) Determine which biomarkers mediate the relationship between selected 
exposome factors and CVD risk
3) Evaluate the predictive performance of an interpretable CVD risk prediction model
   compared to flexible machine learning methods
4) Conduct sub-group analysis by sex to assess whether predictive risk factors differ
   between males and females

# ============================================================
## PROJECT STRUCTURE
# ============================================================

TDS_Group4/
|
|-- 00_data/                                  <- generated (empty on GitHub)
|
|-- 01_scripts/
|   |-- 00_variable_config.R                  # Central config (scenarios, paths, variable lists)
|   |-- 00_variable_dictionary.R              # Display name mappings
|   |-- 00_install_dependencies.R             # R package installation
|   |-- 01_processing/
|   |   |-- 00_initial_clean.R                # Initial data cleaning
|   |   |-- 01_socio_demographics.R           # Sociodemographics
|   |   |-- 02_physical_measures.R            # Anthropometrics
|   |   |-- 03_diet.R                         # Diet score and alcohol
|   |   |-- 04_biological_samples.R           # Blood biomarker cleaning
|   |   |-- 05_health_and_medical_history.R   # Medical history
|   |   |-- 06_early_life_factors.R           # Early life exposures
|   |   |-- 07_environment.R                  # Air pollution & urbanicity
|   |   |-- 08_lifestyle.R                    # Lifestyle variables
|   |   |-- 09_mental_health.R                # Psychosocial variables
|   |   |-- 10_final_clean.R                  # Merge & final variable set
|   |   |-- 11_remove_missing.R               # Missingness filtering
|   |   |-- 12_imputation.R                   # DAG-respecting RF imputation
|   |   |-- 12_imputation_functions.R         # Imputation helper functions
|   |   |-- 13_complete_case.R                # Complete case preprocessing
|   |   |-- 14_population_flow_diagram.R      # Sample flow diagram
|   |   |-- 15_variable_flow_diagram.R        # Variable inclusion flow
|   |   |-- 16_full_table_1.R                 # Table 1 (full cohort)
|   |   `-- 17_imputed_cohort_table_1.R       # Table 1 (imputed cohort)
|   `-- 02_analysis/
|   |   |-- 01_descriptive/
|   |   |   |-- 00_missingness.R              # Missingness visualisation
|   |   |   `-- 01_FAMD_and_Correlation.R     # FAMD & correlation analysis
|   |   |-- 02_univariate/
|   |   |   |-- 00_forest_plot.R              # Univariate OR forest plots
|   |   |   `-- 01_manhattan_plot.R           # Manhattan plot
|   |   |-- 03_stability_selection/
|   |   |   |-- 01_stability_selection_outcome.R    # LASSO 1: exposures -> CVD
|   |   |   |-- 02_stability_selection_mediators1.R # LASSO 2: exposures -> mediators
|   |   |   |-- 03_incrementation_analysis.R        # Incremental AUC analysis
|   |   |   `-- 04_duplication_script.R             # Stability selection diagnostics
|   |   |-- 04_multivariate/
|   |   |   |-- mediation_analysis.R          # Mediation analysis
|   |   |   |-- prediction_analysis.R         # Logistic regression & ROC
|   |   |   `-- sex_stratified.R              # Sex-stratified analysis
|   |   |-- 05_neural_network/
|   |   |   `-- cvd_nn.ipynb                  # Neural network comparator
|   |   |-- 06_decision_trees/
|   |   |   |-- decision_trees_all_w_shap.ipynb
|   |   |   |-- decision_trees_confounders_selected_w_shap.ipynb
|   |   |   `-- decision_trees_confounders_w_shap.ipynb
|   |   `-- 08_combined_roc_curves/
|   |       `-- combined_roc_curves.ipynb     # Cross-model ROC comparison
|   |`-- rds_file_conversion.R                # Converting/removing RDS files
|   `--- scenario_variable_names.R            # Script to extract and export variable names
|
|-- 02_results/                               <- generated (empty on GitHub)
|
|-- 03_jobs/
|   |-- 1-pre-processing.sh                   # Job 1: cleaning pipeline
|   |-- 2-imputation-main.sh                  # Job 2a: main imputation
|   |-- 2-imputation-anthro.sh                # Job 2b: anthropometrics sensitivity imputation
|   |-- 3-exploratory.sh                      # Job 3: EDA, FAMD, univariate
|   |-- 4-stability-selection-main.sh         # Job 4a: LASSO (main)
|   |-- 4-stability-selection-anthro.sh       # Job 4b: LASSO (anthropometrics)
|   |-- 4-stability-selection-cc.sh           # Job 4c: LASSO (complete case)
|   |-- 5-refit-main.sh                       # Job 5a: refit + mediation (main)
|   |-- 5-refit-anthro.sh                     # Job 5b: refit + mediation (anthropometrics)
|   |-- 5-refit-cc.sh                         # Job 5c: refit + mediation (complete case)
|   |-- 6-ml.sh                               # Job 6: ML models (Python)
|   `-- logs/                                 # PBS log files
|
|-- .gitignore
|-- config.py
|-- python_requirements.txt
|-- TDS_Group4.Rproj
`-- README.txt

# ============================================================
## SETUP (run once before submitting jobs)
# ============================================================

Package Installation:
 
1. Install R packages:
 
   eval "$(~/anaconda3/bin/conda shell.bash hook)"
   conda activate *your chosen environment*
   module load NLopt/2.7.1-GCCcore-13.3.0
   Rscript 01_scripts/install_dependencies.R
 
   Note: The NLopt module is required for miceRanger/nloptr to install correctly.

2. Install Python packages:

   eval "$(~/anaconda3/bin/conda shell.bash hook)"
   conda activate *your chosen environment*
   pip install -r python_requirements.txt  

# ============================================================
## HOW TO RUN
# ============================================================

All scripts are run via PBS job scripts from the project root. Submit them in 
order (scripts with the same numbering can be submitted together). 
The pipeline has three analysis scenarios controlled by the SCENARIO 
environment variable (set in relevant job scripts via 00_variable_config.R):

  - main:                          Anthropometrics excluded from analysis
  - sensitivity_anthro_confounder: Anthropometrics treated as confounders
  - complete_case:                 Complete case analysis (no imputation)

~/TDS_Group4/03_jobs/
  1-pre-processing.sh
  2-imputation-main.sh
  2-imputation-anthro.sh
  3-exploratory.sh
  4-stability-selection-main.sh
  4-stability-selection-anthro.sh
  4-stability-selection-cc.sh
  5-refit-main.sh
  5-refit-anthro.sh
  5-refit-cc.sh
  6-ml.sh

Default project root:   cd /rds/general/project/hda_25-26/live/TDS/TDS_Group4/

Step 0 - Before submitting, set the following variables in the job scripts:

  ukb_path                 Path to the raw UK Biobank recoded data file.
                           Set in: 1-pre-processing.sh

  cvd_events_path          Path to the CVD outcome events file.
                           Set in: 1-pre-processing.sh, 3-exploratory.sh

  variable_threshold       Maximum proportion of missingness allowed per variable
                           before it is dropped (default: 0.3).
                           Set in: 1-pre-processing.sh, 3-exploratory.sh

  participant_threshold    Maximum proportion of missingness allowed per 
                           participant before they are dropped (default: 0.3).
                           Set in: 1-pre-processing.sh, 3-exploratory.sh
                           
  cd path                  Project root path.
                           Set in: All scripts
                           
  Source                   R environment in job scripts needs to be replaced 
  (Environment)            with your chosen environment.
                           Set in: All job scripts except 6-ml.sh

Step 1 — Data cleaning and preprocessing:

   qsub 03_jobs/1-pre-processing.sh

   Runs all scripts in 01_processing/ 
   Input:  ukb_path, cvd_events_path
   Output: 00_data/01_processed/

Step 2 — Imputation (submit both simultaneously):

   qsub 03_jobs/2-imputation-main.sh
   qsub 03_jobs/2-imputation-anthro.sh

   DAG-respecting two-stage random forest imputation via miceRanger.
   Stage 1: Exposures imputed using confounders + exposures only.
   Stage 2: Biomarkers imputed using confounders + imputed exposures + biomarkers.
   The complete case scenario skips this step (uses 11b_complete_case.R output).
   Output: 00_data/01_processed/03_imputation/

Step 3 — Descriptive analysis and univariate:

   qsub 03_jobs/3-exploratory.sh

   Missingness visualisation, FAMD/correlation, forest plots, Manhattan plots.
   Output: 02_results/

Step 4 — LASSO variable selection (submit all three simultaneously):

   qsub 03_jobs/4-stability-selection-main.sh
   qsub 03_jobs/4-stability-selection-anthro.sh
   qsub 03_jobs/4-stability-selection-cc.sh

   LASSO stability selection using the sharp package on the selection set (50%).
   Separate runs for outcome model, mediator model, and incremental analysis.
   Output: 00_data/02_analysis/
           02_results/ 

Step 5 — Model refitting and mediation (submit all three simultaneously):

   qsub 03_jobs/5-refit-main.sh
   qsub 03_jobs/5-refit-anthro.sh
   qsub 03_jobs/5-refit-cc.sh

   Mediation analysis (total/direct/indirect effects with bootstrapped CIs),
   prediction model evaluation on the held-out test set, and sex-stratified
   subgroup analysis.
   Output: 02_results/

Step 6 — Machine learning models:

   qsub 03_jobs/6-ml.sh

   Uses the Python environment in group4_python_environment/.
   Decision trees (random forest, XGBoost) and feedforward neural network.
   Output: 02_results/
   
Step 7 - Converting/removing RDS files:

  rds_file_conversion.R
  
  Creates a copy of the 02_results folder (which has either converted RDS files
  into csv files or removed them entirely).
  Output: 02_resultsv2/

Step 8 - Extract Final Variable Names for Reporting
  
  scenario_variable_names.R
  
  Run the extraction script to output `scenario_variable_names.csv` containing the final variable names 
  used across the three different scenarios (Main, Anthropometric, and Complete Case).
  Output: 02_resultsv2/

# ============================================================
## ENVIRONMENT DETAILS
# ============================================================

R version:      4.1.3
Python version: 3.11.13
 
R packages are managed via install_dependencies.R. The full list:
  tidyverse, dplyr, ggplot2, tidyr, purrr, stringr, tibble, readr,
  here, readxl, data.table, miceRanger, future, patchwork, corrplot,
  e1071, FactoMineR, sharp, boot, pROC, broom, scales
 
Python packages are managed via the Python conda environment and python_requirements.txt. 
The full list:
  numpy, matplotlib, scikit-learn, torch, pyreadr, xgboost, shap, pandas 
 
HPC module required: NLopt/2.7.1-GCCcore-13.3.0 (for miceRanger/nloptr)

# ============================================================
## RESULTS STRUCTURE
# ============================================================

All outputs are written to 02_results/ and organised as follows:

02_results/                                   # <- generated (empty on GitHub)
|-- 00_figures/
|   |-- 01_imputation/                        # Convergence diagnostics for each data split
|   |-- 02_univariate/                        # Forest plots and Manhattan plots
|   |-- 03_multivariate/                      # Mediation, ROC, calibration, DAG, odds ratios
|   |-- 04_exploratory/                       # FAMD, correlation, stability selection, incremental plot
|   |-- 05_sensitivity_sex/                   # Sex-stratified results (main scenario)
|   |-- 06_sensitivity_complete_cases/        # Complete case figures
|   `-- 08_sensitivity_anthro_confounder/     # Anthro-as-confounder figures
|
|-- 01_tables/
|   |-- 00_flow_diagrams/                     # Population and variable flow numbers
|   |-- 01_multivariate/                      # AUC, confusion matrix, mediation, odds ratios
|   |-- 02_univariate/                        # Univariate regression results
|   |-- 03_sensitivity_complete_case/         # CC mediation, incremental, sex-stratified
|   |-- 04_neural_network/                    # NN hyperparameters and model comparison
|   |-- 05_sensitivity_sex/                   # Sex-stratified mediation and prediction
|   |-- 06_sensitivity_anthro_confounder/     # Anthro-as-confounder mediation, incremental
|   |-- 07_famd_correlation/                  # FAMD eigenvalues, correlation matrices
|   `-- 08_lasso/                             # LASSO selection tables, incremental AUC
|-- 02_missingness/                           # Missingness diagnostics
`-- 03_ml/                                    # ML model outputs (RF, XGBoost, NN, ROC)

# ============================================================
## AUTHORS
# ============================================================

Group 4 — Klara Eybers, Hassan El-Maoula, Gong Vitayatanagorn, Shiqi Huang, Luane Branthonne
MSc Health Data Analytics and Machine Learning 
Imperial College London, 2025-26
