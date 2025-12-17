# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read ARG carrier contig (ACC) taxonomy
acc_tax <- read.csv("data/ARG_carrying_contig_taxonomy.csv")

# Read project metadata
project_metadata <- read.csv("data/Metadata_metagenome.csv")

# Read ARG carrying contig quantification data (TPM)
acc_tpm <- read.csv("data/ARG_carrying_contig_TPM.csv")

# Read ARG annotation results (based on contigs)
arg_anno_res <- read.csv("data/ARG_annotation_results.csv")

# Read ARG and MGE carrying contigs
amcc_table <- read.csv("data/ARG_and_MGE_carrying_contig.csv")

# Filter samples by agricultural and non-agricultural
sample_agriculture <- (project_metadata %>% filter(Agricultural == "Agricultural"))$Sample
sample_non_agriculture <- (project_metadata %>% filter(Agricultural == "Non-agricultural"))$Sample

gene_cluster_df <- read.csv("data/Contig_gene_cluster_table.csv")

gene_cluster_nr_df <- gene_cluster_df %>% select(source) %>%
  left_join(acc_tax %>% select(contig_name, Species), by = c("source" = "contig_name")) %>% unique()

gene_cluster_nr_prop_df <- gene_cluster_nr_df %>% left_join(acc_tpm, by = c("source" = "Contig")) %>%
  mutate(total_abun_agri = rowSums(across(sample_agriculture), na.rm = TRUE)) %>% 
  mutate(mean_abun_agri = total_abun_agri / length(sample_agriculture), .groups = "drop") %>% 
  mutate(total_abun_non_agri = rowSums(across(sample_non_agriculture), na.rm = TRUE)) %>%
  mutate(mean_abun_non_agri = total_abun_non_agri / length(sample_non_agriculture), .groups = "drop") %>% 
  mutate(agri_prop = mean_abun_agri / (mean_abun_agri + mean_abun_non_agri)) %>% 
  mutate(non_agri_prop = mean_abun_non_agri / (mean_abun_non_agri + mean_abun_agri)) %>% 
  relocate(c("agri_prop", "non_agri_prop"), .after = Species) %>% 
  select(1:4)

for(i in 1:nrow(gene_cluster_nr_prop_df)) {
  tmp_df <- data.frame(
    category = c("agri", "non_agri"),
    value = c(gene_cluster_nr_prop_df[i,"agri_prop"], gene_cluster_nr_prop_df[i,"non_agri_prop"])
  )
  p <- ggplot(tmp_df, aes(x = "", y = value, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = c("agri" = "#e5d190", "non_agri" = "#c53a32")) +
    theme_void() +
    ggtitle(paste0(gene_cluster_nr_prop_df[i,]$source))
  ggsave(p, filename = paste0("figures/Fig2j/", gene_cluster_nr_prop_df[i,]$source, ".pdf"), family = "ArialMT", width=10)
}