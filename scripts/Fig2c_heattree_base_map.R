# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(dplyr)
library(stringr)
library(metacoder)
library(tidyr)

# Read ARG carrying contig (ACC) taxonomy data and exclude the 'Strain' column
acc_tax <- read.csv("data/ARG_carrying_contig_taxonomy.csv") %>% select(-Strain)

# Read project metadata
project_metadata <- read.csv("data/Metadata_metagenome.csv")

# Read ARG carrying contig quantification data (TPM)
acc_tpm <- read.csv("data/ARG_carrying_contig_TPM.csv")

# Format taxonomy columns by appending prefixes
acc_tax$Kingdom <- ifelse(is.na(acc_tax$Kingdom), NA, paste0("k__", acc_tax$Kingdom))
acc_tax$Phylum  <- ifelse(is.na(acc_tax$Phylum),  NA, paste0("p__", acc_tax$Phylum))
acc_tax$Class   <- ifelse(is.na(acc_tax$Class),   NA, paste0("c__", acc_tax$Class))
acc_tax$Order   <- ifelse(is.na(acc_tax$Order),   NA, paste0("o__", acc_tax$Order))
acc_tax$Family  <- ifelse(is.na(acc_tax$Family),  NA, paste0("f__", acc_tax$Family))
acc_tax$Genus   <- ifelse(is.na(acc_tax$Genus),   NA, paste0("g__", acc_tax$Genus))
acc_tax$Species <- ifelse(is.na(acc_tax$Species), NA, paste0("s__", acc_tax$Species))

# Concatenate taxonomy columns into a single string
acc_tax <- acc_tax %>% 
  unite("taxonomy", 2:ncol(acc_tax), sep = ";", na.rm = TRUE, remove = FALSE) %>% 
  dplyr::select(contig_name, taxonomy)

# Merge taxonomy data with TPM data
dat2 <- left_join(acc_tax, acc_tpm, by = c("contig_name" = "Contig"))

# Extract sample groups from metadata
groups <- project_metadata %>% mutate(Type = get("Agricultural"))

# Filter data to include only relevant samples
otu_filter <- dat2 %>% dplyr::select(taxonomy, all_of(groups$Sample))

# Parse taxonomy data into a hierarchical structure
myobj <- parse_tax_data(otu_filter,
                        class_cols = "taxonomy",
                        class_sep = ";",
                        class_regex = "^(.+)__(.+)$",
                        class_key = c(
                          tax_rank = "taxon_rank",
                          tax_name = "taxon_name"))

# Filter to include only bacterial taxa
myobj <- metacoder::filter_taxa(myobj, taxon_names == "Bacteria", subtaxa = TRUE)

# Remove taxa with no reads across all samples
has_no_reads <- rowSums(myobj$data$tax_data[, groups$Sample]) == 0
myobj_d_tax <- metacoder::filter_obs(myobj, "tax_data", !has_no_reads, drop_taxa = TRUE)

# Further filter taxa to include only those with valid names
obj_o <- myobj_d_tax %>%
  metacoder::filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>%
  metacoder::filter_taxa(taxon_ranks == "g", supertaxa = TRUE)

# Calculate taxon abundance for each group
obj_o$data$tax_abund <- calc_taxon_abund(obj_o, "tax_data", cols = groups$Sample)

# Compare groups and calculate statistical differences
obj_o$data$diff_table <- compare_groups(obj_o, 
                                        data = "tax_abund",
                                        cols = groups$Sample,
                                        groups = groups$Type,
                                        func = function(abund_1, abund_2) {
                                          log_ratio_median <- log2(median(abund_1) / median(abund_2))
                                          if (is.nan(log_ratio_median)) {
                                            log_ratio_median <- 0
                                          }
                                          
                                          log_ratio_mean <- log2(mean(abund_1) / mean(abund_2))
                                          if (is.nan(log_ratio_mean)) {
                                            log_ratio_mean <- 0
                                          }
                                          list(log2_median_ratio = log_ratio_median,
                                               log2_mean_ratio = log_ratio_mean,
                                               median_diff = median(abund_1) - median(abund_2),
                                               mean_diff = mean(abund_1) - mean(abund_2),
                                               wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)
                                        })

# Adjust p-values for multiple comparisons
obj_o <- mutate_obs(obj_o, "diff_table", wilcox_p_value_adj = p.adjust(wilcox_p_value, method = "BH"))

# Set insignificant log2 mean ratios to zero
obj_o$data$diff_table$log2_mean_ratio[obj_o$data$diff_table$wilcox_p_value_adj > 0.05] <- 0

# Generate a heat tree visualization
set.seed(5)
ROW_COLOR = "#d79e2e"
COL_COLOR = "#c47673"
heat_tree_matrix(obj_o,
                 data = "diff_table",
                 node_size = n_obs,
                 node_size_range = c(0.018, 0.04),
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio",
                 node_color_range = c(ROW_COLOR, "#DDDDDD", COL_COLOR),
                 row_label_color = ROW_COLOR,
                 col_label_color = COL_COLOR,
                 initial_layout = "large-graph",
                 layout = "reingold-tilford",
                 make_node_legend = TRUE,
                 output_file = "figures/Fig2c.pdf")
