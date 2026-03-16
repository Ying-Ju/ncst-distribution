# Non-Central Skew-t Distribution (NCST)

This repository contains the R code used to reproduce the results in

Hasan, A. & Chen, Y. J. (2026), *Flexible Modeling of Multivariate Skewed and Heavy-Tailed Data via a Non-Central Skew t Distribution with Application to Tumor Shape Data*.

## Repository structure

R/  
Implementation of NCST likelihood evaluation, simulation experiments,
and distribution fitting. This folder includes the following files:

- 00_ncst_core_functions.R — core functions
- 01_2D_ncst_example.R — 2D example used in Section 2
- 02_Quadratic_Form_Validation.R — reproduces Section 4.1
- 03_Skewness_and_Tail_Behavior.R — reproduces Section 4.2
- 04_Effect_of_Degrees_of_Freedom.R — reproduces Section 4.3
- 05_Comparison_Misspecified_Models.R — reproduces Section 4.4 
- 06_Read_Data.R — reproduces real-data results

data/  
Instructions for obtaining the dataset used in the empirical example.

figures/  
Figures generated from the simulations and model fitting.


