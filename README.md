# PhyllosphereAMR

## Overview
PhyllosphereAMR is an code repository for the analysis and visualization of antimicrobial resistance (AMR) in phyllosphere. This repository provides scripts for data loading, processing, statistical evaluation, and graphical representation.


Version
-------
Refactored README — standardized naming, path placeholders and consolidated workflow descriptions.
All example file paths use the placeholder `/path/to/project`. Update `project_root` or environment variables before running scripts.

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
Top-level:
- 0_Raster_stack_initiate_and_check.R  — helper for building/validating raster stacks
- README.md, LICENSE

Primary subdirectories:
- 16s_rRNA_global_mapping/  — 16S workflows (geodata, feature selection, AutoML and mapping)
- Metagenomic_global_mapping/  — shotgun/metagenomic workflows (geodata, feature selection, AutoML and mapping)
- MAG/  — MAG-level analyses, clustering, plotting and statistical comparisons

Naming conventions (recommended)
------------------------------
- File names: <group>_<step>_<short_snake_case_description>.R  
  Examples:
  - 03_01_geodata_extraction.R
  - 02_05_auto_ml_train_optimal.R
- Use two-digit step numbers for consistent ordering (01, 02, 03...).
- Keep variable names and environment variables in snake_case: project_root, data_dir, output_dir, db_dir, conda_dir.
- Replace absolute, user-specific paths with a single configurable `project_root` or environment variable.

Configuration and environment
-----------------------------
- Use project-level configuration at each script top (set `project_root`):
  - project_root <- Sys.getenv("PROJECT_ROOT", "/path/to/project")
- Recommended R environment tooling: renv, packrat, or conda to pin package versions.
- Typical R packages used: dplyr, tidyr, ggplot2, sf, terra, raster, moments, caret/mRMRe/FSelector, H2O or AutoML wrappers.
- Store heavy artifacts (GIS stacks, trained model binaries) outside Git or using Git LFS.

Inputs / outputs (placeholders)
------------------------------
- Input data directory: /path/to/project/data/
- GIS raster stack placeholder: /path/to/project/tools/GIS_rasters/tif_stack_data_197all.rds
- Output / ML artifacts: /path/to/project/data/16S/ML/ or /path/to/project/data/metagenomic/ML/
- Always review the top of each script to set `project_root` and other environment variables prior to running.

High-level pipeline flow
------------------------
Common steps (applies to both 16S and metagenomic flows):
1. Geodata extraction & imputation: extract raster covariates at sample coordinates, compute NA diagnostics, impute (mean/mode), remove poor predictors.
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

Duplication & refactor notes
---------------------------
- Multiple scripts perform identical or near-identical tasks across groups (e.g., geodata extraction, feature selection, AutoML training).
- Recommended: consolidate into parameterized scripts + a shared utilities module (R/utils.R or R/pkg).
- Use YAML/JSON config files to differentiate dataset-specific settings (16S vs metagenomic, agri vs non-agri, seed/config choices).

Data hygiene & reproducibility
-----------------------------
- Keep `project_root` configurable; avoid hard-coded user paths.
- Add smoke-test datasets and CI checks that validate each pipeline step on a small subset.
- Save environment snapshots (renv.lock, environment.yml) next to the repository.
- Document model and dataset versions; store metadata about training runs.