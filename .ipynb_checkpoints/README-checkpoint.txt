PROJECT DESCRIPTION: TDS Group 4

Research question: 
Which external exposome factors are predictive of cardiovascular disease and through which biological factors do they exert their effects? 

Aims:
1) Identify external exposome factors associated with CVD using stability selection LASSO
2) Determine which biological markers and anthropometrics mediate the relationship between selected exposome factors and CVD risk
3) Evaluate the predictive performance of an interpretable CVD risk prediction model compared to flexible machine learning methods


STRUCTURE
rds_project/live/TDS/TDS_Group4
│
├── 00_data
│   ├── 00_extracted        # Raw input data (UK Biobank)
│   ├── 01_processed        # Cleaned and transformed datasets
│   └── 02_analysis         # Analysis output data
│
├── 01_scripts
│   ├── 00_extraction       # empty !!!!!
│   ├── 01_processing       # Data cleaning and feature engineering
│   ├── 02_analysis         # Statistical and ML analyses
│   │    └── 00_descriptive 
│   │    └── 01_imputation 
│   │    └── 02_univariate
│   │    └── 03_exploratory
│   │    └── 04_multivariate
│   │    └── 05_neural_network
│   │    └── 06_decision_trees
│   ├── 03_report           # Reporting scripts ???????????
│   ├── 04_visualisation    # Plotting and distributions
│   └── archive             # Deprecated / exploratory scripts
│
├── 02_results              # Outputs, figures, model results
└── 03_jobs                 # Job scripts (e.g. cluster execution)

110 directories, 999 files


## Requirements

**R** (version x.x.x) — main processing and analysis  
**Python** (version x.x.x) — machine learning models 

Key R packages:
```
tidyverse, mice, FactoMineR, survival, forestplot, ...
```

Key Python packages:
```
torch, scikit-learn, xgboost, pandas, numpy, ...
```

---

## Data Access

Raw data is stored in `00_data/00_extracted/` and sourced from the UK Biobank 
---

## How to Run

### 1. Processing
Run scripts in `01_scripts/01_processing/` sequentially:
```bash
Rscript 01_scripts/01_processing/00_initial_clean.R
Rscript 01_scripts/01_processing/01_socio_demographics.R
# ... through to ...
Rscript 01_scripts/01_processing/12_imputation.R
```

### 2. Analysis

**Descriptive EDA**
```bash
Rscript 01_scripts/02_analysis/00_descriptive/01_sociodemographic_eda.R
```

**Univariate analysis + plots**
```bash
Rscript 01_scripts/02_analysis/02_univariate/00_univariate_w_forest_plot.R
Rscript 01_scripts/02_analysis/02_univariate/01_manhattan_plot.R
```

**Multivariate analysis**
```bash
Rscript 01_scripts/02_analysis/04_multivariate/multivariate_regression_exposome_vs_biological.R
```

**Machine learning (Python)**
```bash
jupyter nbconvert --to notebook --execute 01_scripts/02_analysis/05_neural_network/cvd_nn.ipynb
jupyter nbconvert --to notebook --execute 01_scripts/02_analysis/06_decision_trees/random_forest_xgboost.ipynb
```

---

## Analysis Overview

| Step | Folder | Description |
|------|--------|-------------|
| 0 | `00_descriptive` | EDA across sociodemographics, physical measures, diet |
| 1 | `01_imputation` | Multiple imputation for missing data |
| 2 | `02_univariate` | Univariate associations + Manhattan/forest plots |
| 3 | `03_exploratory` | FAMD, stability selection, incremental analysis |
| 4 | `04_multivariate` | Regression, AUC, mediation analysis |
| 5 | `05_neural_network` | Neural network CVD prediction |
| 6 | `06_decision_trees` | Random forest and XGBoost |

---

## Notes

- `01_scripts/archive/` contains old exploratory scripts — **not part of the pipeline**

---

## Authors

Group 4 — Klara Eybers, Hassan El-Maoula, Gong Vitayatanagorn, Shiqi Huang, Luane Branthonne
MSc Health Data Analytics and Machine Learning 
Imperial College London
[Date]

