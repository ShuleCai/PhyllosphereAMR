# Load required libraries
library(tidyverse)
library(umap)
library(cluster)
library(GGally)
library(clusterSim)
library(fpc)
library(patchwork)
library(ggthemes)
library(extrafont)

# Set random seed for reproducibility
set.seed(123)

origin_data <- read.csv("/path/to/data/source_data.csv", check.names = FALSE)

# Define feature categories for analysis

# Antibiotic Resistance Genes (ARG) categories
ARG <- c(
  "Aminoglycoside", "Bacitracin", "Beta_lactam", "Chloramphenicol",
  "Defensin", "Fosfomycin", "Macrolide-Lincosamide-Streptogramin",
  "Multidrug", "Other_peptide_antibiotics", "Polymyxin",
  "Quinolone", "Rifamycin", "Tetracycline"
)

# Mobile Genetic Elements (MGE) categories
MGE <- c(
  "Integration/Excision", "Phage", "Replication/Recombination/Repair",
  "Stability/Transfer/Defense", "Transfer"
)

# Virulence Factor Genes (VFG) categories
VFG <- c(
  "Adherence", "Antimicrobial Activity/Competitive Advantage",
  "Biofilm", "Effector Delivery System", "Exoenzyme",
  "Exotoxin", "Immune Modulation", "Invasion", "Motility",
  "Nutritional/Metabolic Factor", "Others", "Regulation",
  "Stress Survival"
)

# Calculate counts for each category and select relevant columns
selected_data <- origin_data %>%
  mutate(
    ARG_count = rowSums(dplyr::select(., all_of(ARG))),
    MGE_count = rowSums(dplyr::select(., all_of(MGE))),
    VFG_count = rowSums(dplyr::select(., all_of(VFG)))
  ) %>%
  dplyr::select(ARG_count, MGE_count, VFG_count)

# Perform UMAP dimensionality reduction
umap_result <- umap(selected_data)

# Evaluate clustering performance across different k values
k_values <- 3:12
silhouette_scores <- numeric(length(k_values))
davies_bouldin_scores <- numeric(length(k_values))
calinski_harabasz_scores <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  k <- k_values[i]

  # Perform K-means clustering
  kmeans_result <- kmeans(selected_data, centers = k, nstart = 25)

  # Calculate Silhouette Score
  silhouette_avg <- mean(silhouette(
    kmeans_result$cluster,
    dist(as.matrix(selected_data))
  )[, 3])
  silhouette_scores[i] <- silhouette_avg

  # Calculate Davies-Bouldin Score
  davies_bouldin_scores[i] <- index.DB(as.matrix(selected_data),
    cl = kmeans_result$cluster
  )$DB

  # Calculate Calinski-Harabasz Score
  cluster_stats <- cluster.stats(
    d = dist(as.matrix(selected_data)),
    clustering = kmeans_result$cluster
  )
  calinski_harabasz_scores[i] <- cluster_stats$ch
}

# Create results dataframe
results <- data.frame(
  k = k_values,
  Silhouette_Score = silhouette_scores,
  Davies_Bouldin_Score = davies_bouldin_scores,
  Calinski_Harabasz_Score = calinski_harabasz_scores
)

# Create visualization of clustering evaluation metrics
font_size <- 10
axis_font_size <- 8
optimal_k <- 3 # Based on the evaluation above

# Silhouette Score plot
p1 <- ggplot(results, aes(x = k, y = Silhouette_Score)) +
  geom_line(color = "orange") +
  geom_point(color = "orange3") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "black") +
  labs(
    title = "Silhouette Score",
    x = "Number of Clusters (K)",
    y = "Silhouette Score"
  ) +
  theme_base() +
  theme(
    plot.title = element_text(size = font_size, hjust = 0.5),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_font_size)
  )

# Davies-Bouldin Score plot
p2 <- ggplot(results, aes(x = k, y = Davies_Bouldin_Score)) +
  geom_line(color = "orange") +
  geom_point(color = "orange3") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "black") +
  labs(
    title = "Davies-Bouldin Score",
    x = "Number of Clusters (K)",
    y = "Davies-Bouldin Score"
  ) +
  theme_base() +
  theme(
    plot.title = element_text(size = font_size, hjust = 0.5),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_font_size)
  )

# Combine and save evaluation plots
combined_plot <- p1 + p2 + plot_layout(nrow = 1)
ggsave(combined_plot,
  filename = "/path/to/figure/kmeans_score.pdf",
  family = "ArialMT",
  height = 4, width = 7
)

# Perform K-means clustering with optimal cluster number
kmeans_result <- kmeans(selected_data, centers = best_k)

# Add clustering results to UMAP coordinates
umap_df <- as.data.frame(umap_result$layout) %>%
  mutate(Cluster = as.factor(kmeans_result$cluster))

# Create and save UMAP visualization
ggplot(umap_df, aes(x = V1, y = V2, color = Cluster)) +
  geom_point() +
  labs(
    title = "UMAP Plot with K-means Clustering",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  ) +
  theme_base()

ggsave(
  filename = "/path/to/figure/clutser_umap.pdf",
  family = "ArialMT",
  width = 5.8, height = 5, units = "in"
)

# Add cluster labels to the original data
origin_data_s <- origin_data %>%
  mutate(
    ARG_count = selected_data$ARG_count,
    MGE_count = selected_data$MGE_count,
    VFG_count = selected_data$VFG_count,
    Cluster = as.factor(kmeans_result$cluster)
  )

# Create parallel coordinates plot
pdf("/path/to/figure/cluster_parallel_coordinates_plot.pdf",
  family = "ArialMT", width = 6, height = 5
)

ggparcoord(origin_data_s,
  columns = which(names(origin_data_s) %in% c("ARG_count", "MGE_count", "VFG_count")),
  groupColumn = "Cluster",
  scale = "uniminmax",
  showPoints = TRUE,
  alphaLines = 0.6
) +
  labs(y = "Normalized Richness") +
  scale_color_manual(values = c("red", "darkblue", "#69b3a2")) +
  theme_bw()

dev.off()
