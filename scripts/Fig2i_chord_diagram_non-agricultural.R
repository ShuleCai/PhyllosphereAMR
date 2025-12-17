# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(reshape2)
library(dplyr)
library(tibble)
library(circlize)

# Read co-occurrence abundance table (ARG vs MGE categories)
both_df_taxo_abun <- read.csv("./data/Contig_co_occurrence_ARG_MGE.csv", check.names = F)

# Summarize non-agricultural abundance by ARG_group and MGE_major_category
# The summarise produces a table with columns: ARG_group, MGE_major_category, abun_nona
both_df_taxo_abun_nona <- both_df_taxo_abun %>% 
  # group rows by ARG_group and MGE_major_category and sum the 'Non.agricultural' values
  group_by(ARG_group, MGE_major_category) %>% 
  dplyr::summarise(abun_nona=sum(`Non.agricultural`)) %>% 
  # sort by descending abundance for inspection
  arrange(desc(abun_nona))

# Convert the long-format summary table to a wide matrix suitable for chordDiagram
# Rows = ARG_group, Columns = MGE_major_category, fill missing combinations with 0
mat_nona <- both_df_taxo_abun_nona %>% 
  dcast(ARG_group ~ MGE_major_category, fill = 0) %>% 
  tibble::column_to_rownames("ARG_group") %>% 
  as.matrix

# Define colors for each grid sector (ARG groups and MGE categories)
grid.col = c(
  multidrug = "#B2182B",
  `macrolide-lincosamide-streptogramin` = "#E69F00",
  aminoglycoside = "#56B4E9",
  bacitracin = "#009E73",
  beta_lactam = "#F0E442",
  chloramphenicol = "#CC6666",
  fosfomycin = "#D55E00",
  other_peptide_antibiotics = "#CC79A7",
  polymyxin = "#0072B2",
  quinolone = "#9999CC",
  tetracycline = "#66CC99",
  `integration/excision` = "grey",
  `phage` = "grey",
  `replication/recombination/repair` = "grey",
  `stability/transfer/defense` = "grey",
  `transfer` = "grey"
)

# Open PDF device and draw the chord diagram
# chordDiagram will draw links weighted by values in mat_nona
pdf("./figures/Fig2i.pdf", family = "ArialMT")
circlize::chordDiagram(
  mat_nona, 
  transparency = 0.4,
  grid.col = grid.col,
  annotationTrackHeight = c(0.05, 0.04),
  link.sort = TRUE
)
dev.off()
