# UKB_PA_Fracture
UK Biobank: This repository contains the code for 3 linked studies of activity and fracture risk  

- Study 1: Cross sectional association of self-report PA and self-report fracture - Code starts wth CXA
- Study 2: Survival study of self-report PA and HES fracture - Code starts with TTE
- Study 3: Survival study of accelerometer PA and HES fracture

## How to reproduce

1. Place UKB raw data files in `/data_raw/` (not under version control).
2. For Study 1 Run scripts 00 through to 11
3. For Study 2 Run scripts 00 through to 07 (loads all helpers), then run 12 - 15 for survival analysis
4. Outputs (tables and figures) will appear in `/results/`.

## Dependencies
- R version and structure of folders is detailed in the config.env, acrtivate.R and Renv.lock files
- Key packages: tidyverse, survival, tableone, ggplot2, knitr, rmarkdown

## Data security
- No UK Biobank data had been committed to Git.
- Data storage complies with institutional and UKB data access rules.
- Files were loaded manually rather than by gitpush, a git.ignore file should be used if working straight from R


