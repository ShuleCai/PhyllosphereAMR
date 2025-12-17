# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(tidyverse)
library(umap)
library(cluster)
library(dplyr)

# Read the original data
origin_data <- read.csv("./data/MAG_bigtable.csv", check.names = F)

# Define ARG, MGE, and VFG column names
ARG <- c("Aminoglycoside", "Bacitracin", "Beta_lactam", "Chloramphenicol", 
         "Defensin", "Fosfomycin", "Macrolide-Lincosamide-Streptogramin", 
         "Multidrug", "Other_peptide_antibiotics", "Polymyxin", 
         "Quinolone", "Rifamycin", "Tetracycline")

MGE <- c("Integration/Excision", "Phage", "Replication/Recombination/Repair", 
         "Stability/Transfer/Defense", "Transfer")

VFG <- c("Adherence", "Antimicrobial Activity/Competitive Advantage", 
         "Biofilm", "Effector Delivery System", "Exoenzyme", 
         "Exotoxin", "Immune Modulation", "Invasion", "Motility", 
         "Nutritional/Metabolic Factor", "Others", "Regulation", 
         "Stress Survival")

# Calculate counts for each category
selected_data <- origin_data %>%
  mutate(ARG_count = rowSums(dplyr::select(., all_of(ARG))),
         MGE_count = rowSums(dplyr::select(., all_of(MGE))),
         VFG_count = rowSums(dplyr::select(., all_of(VFG)))) %>%
  dplyr::select(ARG_count, MGE_count, VFG_count)

# Perform UMAP dimensionality reduction
umap_result <- umap(selected_data)

# Determine the optimal number of clusters using silhouette scores
silhouette_scores <- numeric(10)
for (i in 3:10) {
  kmeans_result <- kmeans(selected_data, centers = i)
  silhouette_scores[i] <- silhouette(kmeans_result$cluster, dist(selected_data)) %>%
    summary() %>%
    .$avg.width
}

# Identify the best number of clusters
best_k <- which.max(silhouette_scores)

# Perform K-means clustering with the optimal number of clusters
kmeans_result <- kmeans(selected_data, centers = best_k)

# Add clustering results to UMAP output
umap_df <- as.data.frame(umap_result$layout) %>%
  mutate(Cluster = as.factor(kmeans_result$cluster))

# Plot UMAP results
library(ggthemes)
library(extrafont)
ggplot(umap_df, aes(x = V1, y = V2, color = Cluster)) +
  geom_point() +
  labs(title = NULL, x = "UMAP 1", y = "UMAP 2") +
  theme_base()
ggsave(filename = "./figures/Fig4a.pdf", family = "ArialMT", width = 5.8, height = 5, units = "in")

# Add cluster labels to the original data
selected_data_with_cluster <- selected_data
selected_data_with_cluster$Cluster <- as.factor(kmeans_result$cluster)
origin_data_s <- origin_data %>% mutate(ARG_count = selected_data$ARG_count,
                                        MGE_count = selected_data$MGE_count,
                                        VFG_count = selected_data$VFG_count,
                                        Cluster = selected_data_with_cluster$Cluster)

# Generate parallel coordinate plots
library(GGally)
pdf("./figures/Fig4b.pdf", family = "ArialMT", width = 6, height = 5)
ggparcoord(origin_data_s %>% mutate(Cluster = as.factor(Cluster)), columns = 43:45, groupColumn = "Cluster", scale = "uniminmax",
           showPoints = T, alphaLines = 0.6) +
  labs(y = "Normalized richness") +
  scale_color_manual(values = c("red", "#6f9af3", "#69b3a2")) +
  theme_bw() +
  theme()
dev.off()


