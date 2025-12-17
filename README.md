# PhyllosphereAMR

## Project Overview
PhyllosphereAMR is an code repository for the analysis and visualization of antimicrobial resistance (AMR) in phyllosphere. This repository provides scripts for data loading, processing, statistical evaluation, and graphical representation.

## Repository Structure

### 1. `data/`
This directory contains all the input datasets used in the analysis. Key files include:
- `AMRcategory_abundance_tpm.csv`: Abundance of AMR categories in TPM.
- `ARG_abundance_tpm.csv`: Abundance of ARGs in TPM.
- `ARG_and_MGE_carrying_contig.csv`: Contigs carrying both ARGs and MGEs.
- `ARG_annotation_results.csv`: Annotation results for ARGs.
- `MAG_bigtable.csv`: Metadata for metagenome-assembled genomes (MAGs).
- `Metadata_metagenome.csv`: Metadata for metagenomic samples.

### 2. `figures/`
This directory stores all the generated figures from the analysis. Subdirectories include:
- `Fig2j/`: Contains pie charts for gene cluster proportions.

### 3. `scripts/`
This directory contains all the R scripts used for data analysis and visualization. Key scripts include:
- `Fig1a_resistome_patterns_circular.R`: Generates circular plots for resistome patterns.
- `Fig2f_heatmap_ARG_carrying_contig.R`: Creates heatmaps for ARG-carrying contigs.
- `Fig4ab_MAG_kmeans_umap.R`: Performs UMAP and K-means clustering on MAG data.
- `global_mapping/`: Contains scripts for global mapping and extrapolation analyses.

### 4. `global_mapping/`
This subdirectory includes scripts for predictive modeling and risk analysis:
- `1_prediction_agri.R`: Predictions for agricultural samples.
- `4_risk_minibacth_kmeans.R`: Mini-batch K-means clustering for risk analysis.
- `8_extrapolation_PCA_convex_hull.R`: PCA and convex hull analysis for extrapolation.

## Key Features
- **Data Processing**: Scripts for cleaning and merging datasets.
- **Visualization**: Heatmaps, chord diagrams, pie charts, and UMAP plots.
- **Statistical Analysis**: Kruskal-Wallis tests, clustering, and dimensionality reduction.

## How to Use
1. Clone the repository:
   ```bash
   git clone https://github.com/ShuleCai/PhyllosphereAMR.git
   ```
2. Navigate to the project directory:
   ```bash
   cd PhyllosphereAMR
   ```
3. Install required R packages (e.g., `dplyr`, `ggplot2`, `ComplexHeatmap`).
4. Run the scripts in the `scripts/` directory to reproduce the analyses.

## Dependencies
The project relies on the following R packages:
- `dplyr`
- `ggplot2`
- `ComplexHeatmap`
- `circlize`
- `tidyr`
- `umap`
- `cluster`
- `ggpubr`
- `agricolae`

## Outputs
The analysis generates:
- Heatmaps for ARG and MGE abundance.
- Chord diagrams for ARG-MGE co-occurrence.
- Pie charts for gene cluster proportions.
- UMAP and clustering visualizations.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## Acknowledgments
Special thanks to all contributors and collaborators who made this project possible.
