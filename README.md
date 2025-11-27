# PhyllosphereAMR

## Overview
PhyllosphereAMR is an code repository for the analysis and visualization of antimicrobial resistance (AMR) in phyllosphere. This repository provides scripts for data loading, processing, statistical evaluation, and graphical representation.

Project summary
---------------
PhyllosphereAMR implements end-to-end analyses for mapping microbiome- and MAG-derived signals related to antimicrobial resistance (AMR) across global / spatial datasets. The repository includes:
- geospatial extraction and raster handling,
- per-dataset preprocessing (16S and shotgun/metagenomic),
- variable selection (correlation filtering + mRMR),
- AutoML model training and spatial prediction,
- visualization and post-processing for MAG- and community-level results,
- spatial extrapolation / future-scenario evaluation and risk mapping.

Goals
-----
- Provide reproducible pipelines to extract spatial covariates, select robust predictors, train predictive models and produce maps.
- Make scripts modular, parameterizable, and portable across different environments.
- Replace user-specific absolute paths with configurable root variables (e.g., `/path/to/project`).
- Reduce duplication by consolidating repeated logic into shared utilities and parameterized scripts.

Repository layout
----------------
Primary subdirectories:
- 16s_rRNA_global_mapping/  — 16S workflows (geodata, feature selection, AutoML and mapping)
- Metagenomic_global_mapping/  — shotgun/metagenomic workflows (geodata, feature selection, AutoML and mapping)
- MAG/  — MAG-level analyses, clustering, plotting and statistical comparisons

Naming conventions 
------------------------------
- File names: <group>_<step>_<short_snake_case_description>.R  
  Examples:
  - 03_01_geodata_extraction.R
  - 02_05_auto_ml_train_optimal.R
- Keep variable names and environment variables in snake_case

High-level pipeline flow
------------------------
Common steps (applies to both 16S and metagenomic flows):
1. Geodata extraction & imputation: extract raster covariates at sample coordinates, compute NA diagnostics, remove poor predictors.
2. Feature selection: correlation-based pruning + mRMR (or equivalent) to generate compact predictor sets.
3. Model training (AutoML): seed-based exploration → optimal model training; maintain reproducibility via controlled seeds and env snapshots.
4. Global prediction: predict across the raster stack, combine results, evaluate and create maps/reports.
5. Extrapolation & scenario analysis: PCA convex hull, band-based extrapolation checks, future climate or landuse scenario mapping.

Files organized by group
------------------------
16s_rRNA_global_mapping/
- Geodata & imputation: 3_16s_rRNA_global_mapping_01_geodata.R
- Feature selection: 3_16s_rRNA_global_mapping_02_feature_selection_cor_mrmr.R
- AutoML training: 03–06 (seed exploration and optimal training for agri/non-agri)
- Prediction & mapping: 07–16 (prediction, combined maps, future scenario analyses)

Metagenomic_global_mapping/
- Mirror flow with 2_* naming: geodata extraction, 2_* feature selection, 2_* AutoML, 2_* prediction & visualization

MAG/
- 1_* scripts for MAG clustering, family-level visuals, pairwise testing and distribution analyses.
