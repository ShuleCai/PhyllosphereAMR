# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(dplyr)
library(stringr)
library(treemap)
library(hash)

# Read ARG carrier contig (ACC) taxonomy data
acc_tax <- read.csv("data/ARG_carrying_contig_taxonomy.csv") 

# Read project metadata
project_metadata <- read.csv("data/Metadata_metagenome.csv")

# Read ARG carrying contig quantification data (TPM)
acc_tpm <- read.csv("data/ARG_carrying_contig_TPM.csv")

# Read ARG and MGE carrying contigs
amcc_names <- read.csv("data/ARG_and_MGE_co-carrying_contig_names.csv")

# Filter samples for agricultural environments
sample_agriculture <- (project_metadata %>% filter(Agricultural == "Agricultural"))$Sample

# Process agricultural data
acc_agri_union_tpm <- data.frame(contig_name=acc_tpm$Contig,
                              with_MGE = acc_tpm$Contig %in% amcc_names$contig_name,
                              abun_avg=rowMeans(acc_tpm[, sample_agriculture])) %>%
  left_join(acc_tax, by = "contig_name") %>% 
  mutate(random_color = ifelse(with_MGE, "#ef8b70","#fae2ae"),
         Genus_plot = paste0(Family, "_", 1:nrow(acc_tpm)),
         abun_avg_sqrt = abun_avg) %>% 
  filter(!is.na(Order) & abun_avg > 0)

# Generate a treemap visualization
pdf("figures/Fig2a.pdf", family="ArialMT", width=8, height = 3)
treemap(acc_agri_union_tpm,
        index=c("Phylum", "Class", "Order", "Family", "Genus_plot"),
        vSize="abun_avg_sqrt",
        type="color",
        vColor="random_color",
        algorithm="squarified",
        fontsize.labels = c(0, 0, 0, 5, 0),
        fontface.labels = c("bold", "plain", "plain", "plain", "plain"),
        align.labels = list(c("left", "top"), c("right", "top"), c("right", "bottom"), c("left", "top"), c("center", "center")),
        fontcolor.labels = c("blue", "red", "orange", "black", "black"),
        inflate.labels = TRUE,
        bg.labels = 0,
        title = "",
        border.col = "grey50",
        lowerbound.cex.labels = 0,
        border.lwds = 0.5)
dev.off()

