library(ClusterR)
library(cluster)
library(fpc)

# Define the range of k values to test
k_values <- 2:16
n_k <- length(k_values)

# Initialize vectors to store results
wss_mbkm <- numeric(n_k) # WSS values for Elbow Method
silhouette_scores <- numeric(n_k) # Silhouette scores
davies_bouldin_scores <- numeric(n_k) # Davies-Bouldin index
calinski_harabasz_scores <- numeric(n_k) # Calinski-Harabasz index

final_result <- read.csv("/path/to/data/Final_1_risk_abun_pred_mean_sd.csv", row.names = 1)

# Extract the predicted values for clustering
y_values <- final_result$mean_predict
n <- length(y_values) # Use length() instead of nrow() as y_values is a vector

# Loop through k values to calculate the clustering metrics
for (i in seq_along(k_values)) {
  k <- k_values[i]
  cat("Processing k =", k, "...\n")

  # Perform MiniBatchKMeans clustering (convert input to matrix form)
  mbkm <- MiniBatchKmeans(matrix(y_values, ncol = 1),
    clusters = k,
    batch_size = 1000, max_iters = 800,
    initializer = "kmeans++"
  )

  # Get the cluster labels (in vector form)
  mbkm_clusters <- predict_MBatchKMeans(matrix(y_values, ncol = 1), mbkm$centroids)

  # 1. Calculate WSS (Within-cluster sum of squares) for the Elbow method - Optimized for vectors
  wss <- 0
  for (j in 1:k) {
    # Retrieve the points in each cluster
    cluster_points <- y_values[mbkm_clusters == j]
    if (length(cluster_points) > 0) {
      wss <- wss + sum((cluster_points - mbkm$centroids[j, 1])^2)
    }
  }
  wss_mbkm[i] <- wss

  # 2. Calculate the Silhouette Score (Sampling method to avoid large matrices)
  sample_size <- min(10000, n) # Sample up to a maximum of 10000 points
  set.seed(124)
  sampled_indices <- sample(1:n, size = sample_size)

  # Calculate the Silhouette Score for the sampled data
  sample_data <- y_values[sampled_indices]
  sample_clusters <- mbkm_clusters[sampled_indices]

  if (length(unique(sample_clusters)) == 1) {
    silhouette_scores[i] <- NA # Cannot compute if there's only one cluster
  } else {
    sample_dists <- dist(sample_data)
    sil <- silhouette(sample_clusters, sample_dists)
    silhouette_scores[i] <- mean(sil[, 3])
  }

  # 3. Calculate the Davies-Bouldin Index - Optimized for vectors
  # Calculate the average distance for each cluster
  s <- numeric(k)
  for (cluster_idx in 1:k) {
    cluster_points <- y_values[mbkm_clusters == cluster_idx]
    if (length(cluster_points) > 0) {
      s[cluster_idx] <- mean(abs(cluster_points - mbkm$centroids[cluster_idx, 1]))
    } else {
      s[cluster_idx] <- 0
    }
  }

  # Calculate the Davies-Bouldin Index
  db_sum <- 0
  for (i_clust in 1:k) {
    max_ratio <- -Inf
    for (j_clust in 1:k) {
      if (i_clust != j_clust) {
        # Calculate the distance between cluster centers
        centroid_dist <- abs(mbkm$centroids[i_clust, 1] - mbkm$centroids[j_clust, 1])
        if (centroid_dist > 0) {
          ratio <- (s[i_clust] + s[j_clust]) / centroid_dist
          if (ratio > max_ratio) max_ratio <- ratio
        }
      }
    }
    db_sum <- db_sum + max_ratio
  }
  davies_bouldin_scores[i] <- db_sum / k

  # 4. Calculate the Calinski-Harabasz Index - Optimized for vectors
  global_mean <- mean(y_values)
  tss <- sum((y_values - global_mean)^2) # Total sum of squares
  bss <- tss - wss # Between-cluster sum of squares
  ch <- (bss / (k - 1)) / (wss / (n - k))
  calinski_harabasz_scores[i] <- ch
}

# Create a dataframe to store the results
results <- data.frame(
  k = k_values,
  wss = wss_mbkm,
  silhouette = silhouette_scores,
  davies_bouldin = davies_bouldin_scores,
  calinski_harabasz = calinski_harabasz_scores
)

# Plot the Elbow method (WSS vs k values)
plot(results$k, results$wss,
  type = "b", pch = 19, col = "blue",
  main = "Elbow Method for MiniBatchKMeans",
  xlab = "Number of clusters (k)", ylab = "WSS"
)

# Libraries for visualization
library(patchwork)
library(ggthemes)
library(ggplot2)
library(extrafont)
font_size <- 10
axis_font_size <- 8
x_value <- 12
results_s <- results[2:14, ]

# Plot the Silhouette Score as a line plot
p1 <- ggplot(results_s, aes(x = k, y = silhouette)) +
  geom_line(color = "orange") +
  geom_point(color = "orange3") +
  geom_vline(xintercept = x_value, linetype = "dashed", color = "black") +
  labs(title = "Silhouette Score", x = "Number of Clusters (K)", y = "Silhouette Score") +
  theme_base() +
  theme(
    panel.background = element_blank(),
    plot.title = element_text(size = font_size, hjust = 0.5),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_font_size)
  )

# Plot the Davies-Bouldin Score as a line plot
p2 <- ggplot(results_s, aes(x = k, y = davies_bouldin)) +
  geom_line(color = "orange") +
  geom_point(color = "orange3") +
  geom_vline(xintercept = x_value, linetype = "dashed", color = "black") +
  labs(title = "Davies-Bouldin Score", x = "Number of Clusters (K)", y = "Davies-Bouldin Score") +
  theme_base() +
  theme(
    panel.background = element_blank(),
    plot.title = element_text(size = font_size, hjust = 0.5),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_font_size)
  )

# Plot the Calinski-Harabasz Score as a line plot
p3 <- ggplot(results_s, aes(x = k, y = calinski_harabasz)) +
  geom_line(color = "orange") +
  geom_point(color = "orange3") +
  geom_vline(xintercept = x_value, linetype = "dashed", color = "black") +
  labs(title = "Calinski-Harabasz Score", x = "Number of Clusters (K)", y = "Calinski - Harabasz Score") +
  theme_base() +
  theme(
    panel.background = element_blank(),
    plot.title = element_text(size = font_size, hjust = 0.5),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_font_size)
  )

# Use patchwork to combine the plots
combined_plot <- p1 + p2 + plot_layout(nrow = 1)
combined_plot
ggsave("/path/to/data/Final_2_pred_risk_bestk_minibatch.pdf",
  family = "ArialMT",
  height = 4, width = 8
)
results_s %>% write.csv("/path/to/data/Final_2_pred_risk_bestk_minibatch.csv")

# Consider the best k value based on silhouette score
best_k_silhouette <- 12

# Step 2: Perform k-means clustering with the best k value
set.seed(124)

# Perform MiniBatchKMeans clustering (input converted to matrix)
mbkm <- MiniBatchKmeans(matrix(y_values, ncol = 1),
  clusters = best_k_silhouette,
  batch_size = 1000, max_iters = 800,
  initializer = "kmeans++"
)
mbkm_clusters <- predict_MBatchKMeans(matrix(y_values, ncol = 1), mbkm$centroids)
# Save the trained model and centroids
saveRDS(mbkm, "/path/to/data/Final_2_minibatch_kmeans_model_best12.rds")
write.csv(mbkm$centroids, "/path/to/data/Final_2_minibatch_kmeans_model_best12_centroids.csv")

# Retrieve the cluster centroids and sort them in ascending order, then get the sorted indices
centroids <- mbkm$centroids[, 1]
sorted_indices <- order(centroids)

# Create a mapping of old cluster labels to new labels based on sorted centroids
label_mapping <- 1:length(centroids)
names(label_mapping) <- sorted_indices

# Update the cluster labels based on the new order of centroids
new_mbkm_clusters <- label_mapping[as.character(mbkm_clusters)]

# Step 3: Use t-SNE for 2D visualization of the clusters
# t-SNE input needs to be in matrix format
set.seed(123)
sampled_indices <- sample(1:length(y_values), size = 500)
tsne_result <- Rtsne(as.matrix(y_values[sampled_indices]), dims = 2, check_duplicates = FALSE)

# Create a new dataframe to store the t-SNE results
tsne_df <- data.frame(
  x = tsne_result$Y[, 1],
  y = tsne_result$Y[, 2],
  cluster = as.factor(new_mbkm_clusters[sampled_indices])
)

# Step 4: Plot the t-SNE 2D visualization
library(ggplot2)
ggplot(tsne_df, aes(x = x, y = y, color = cluster)) +
  geom_point() +
  labs(
    title = "t-SNE Visualization of k-means Clustering (500 Random Samples)",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2"
  ) +
  theme_minimal()

# Save the t-SNE plot to a file
ggsave(filename = "/path/to/data/Final_2_minibatch_kmeans_model_best12_tnse.pdf", width = 6, height = 6)

# Add the new cluster labels to the final result dataframe
final_result_risk <- final_result %>% mutate(Risk = new_mbkm_clusters)

# Save the updated dataframe with cluster labels
final_result_risk %>%
  write.csv("/path/to/data/Final_3_world_pixel_Risk.csv")

# Load the country boundaries dataset
library(rnaturalearth)
countries <- ne_countries(scale = "medium", returnclass = "sf")

# Load the river data
rivers <- ne_download(
  scale = "medium", type = "rivers_lake_centerlines",
  category = "physical", returnclass = "sf"
)

# Generate an extended 12-color Spectral color palette
library(RColorBrewer)
spectral_pal <- rev(brewer.pal(11, "Spectral")) # Reverse the Spectral color palette
extended_spectral <- colorRampPalette(spectral_pal)(12) # Extend to 12 colors

# Create the risk level map with the new color palette and geographical features
ggplot() +
  # Plot the risk level heatmap
  geom_raster(
    data = final_result_risk,
    aes(x = x, y = y, fill = as.factor(Risk)),
    na.rm = FALSE
  ) + # Keep NA values
  geom_sf(data = countries, color = "grey50", alpha = 0.2, size = 0.1, fill = NA) +
  geom_sf(data = rivers, color = "lightblue", alpha = 0.5, size = 0.05) +
  # Set the coordinate system (focus on land area)
  coord_sf(ylim = c(-60, 90), expand = FALSE, crs = "+proj=longlat") + # Use latitude-longitude coordinate system
  # Set the color scale
  scale_fill_manual(
    name = "Risk Level", # Legend title
    values = extended_spectral, # Use the extended color palette
    na.value = "white", # Set NA values to white
    guide = guide_legend( # Customize legend guide
      direction = "horizontal", # Position the legend horizontally
      title.position = "top", # Title on top
      title.hjust = 0.5, # Center the title
      label.position = "bottom", # Labels at the bottom
      nrow = 1, # Arrange legend items in a single row
      keywidth = unit(0.9, "cm"), # Legend key width
      keyheight = unit(0.3, "cm") # Legend key height
    )
  ) +
  # Set the theme for the plot
  theme_void() +
  theme(
    legend.position = "bottom", # Position legend at the bottom
    legend.justification = "center", # Center the legend
    legend.box = "horizontal", # Place legend items horizontally
    legend.key.spacing = unit(0, "cm"), # Remove spacing between legend keys
    legend.margin = margin(t = 0, r = 0, b = 5, l = 0), # Reduce margin at the bottom
    legend.text = element_text(size = 9, margin = margin(t = 2)), # Style the legend text
    legend.title = element_text(face = "bold", size = 10, margin = margin(b = 5)), # Style the legend title
    plot.margin = margin(10, 10, 5, 10) # Adjust plot margins
  )

# Save the risk level map to a file
ggsave("/path/to/data/Final_3_global_mapping_risk.pdf", width = 10, height = 8)
ggsave("/path/to/data/Final_3_global_mapping_risk.png")

# Plot x, y subplot for risk distribution
library(tidyverse)
library(cowplot)
library(scales)

# Prepare the data: final_result_risk contains x, y, Risk columns
# Ensure Risk is treated as a factor (levels 1-12)
final_result_risk <- final_result_risk %>%
  mutate(Risk = factor(Risk, levels = 1:12))

# Create the extended 12-color Spectral palette (consistent with the map)
spectral_pal <- rev(brewer.pal(11, "Spectral"))
extended_spectral <- colorRampPalette(spectral_pal)(12)

# Modify the longitude bar plot function (invert the y-axis)
create_x_barplot <- function(data) {
  x_data <- data %>%
    count(x, Risk) %>%
    complete(x = seq(-180, 180), Risk, fill = list(n = 0))

  ggplot(x_data, aes(x = x, y = n, fill = Risk)) +
    geom_bar(stat = "identity", width = 1, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = extended_spectral) +
    scale_x_continuous(limits = c(-180, 180), expand = expansion(0)) +
    scale_y_reverse(expand = expansion(c(0.05, 0))) + # Reverse y-axis to extend values downwards
    labs(x = NULL, y = NULL) +
    theme_minimal_hgrid() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
}

# Function to create latitude-based bar plots (y-axis direction)
create_y_barplot <- function(data) {
  # Calculate the risk distribution at each latitude (y)
  y_data <- data %>%
    count(y, Risk) %>%
    complete(y = seq(-60, 90), Risk, fill = list(n = 0))

  ggplot(y_data, aes(x = y, y = n, fill = Risk)) +
    geom_bar(stat = "identity", width = 1, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = extended_spectral) +
    scale_x_continuous(limits = c(-60, 90), expand = expansion(0)) +
    scale_y_continuous(expand = expansion(c(0, 0.05))) +
    labs(x = NULL, y = NULL) +
    coord_flip() +
    theme_minimal_hgrid() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
}

# Create the main map (slightly modified from the original code)
main_map <- ggplot() +
  geom_raster(
    data = final_result_risk,
    aes(x = x, y = y, fill = Risk),
    na.rm = FALSE
  ) +
  geom_sf(data = countries, color = "grey50", alpha = 0.2, size = 0.1, fill = NA) +
  geom_sf(data = rivers, color = "lightblue", alpha = 0.5, size = 0.05) +
  coord_sf(ylim = c(-60, 90), expand = FALSE) +
  scale_fill_manual(
    values = extended_spectral, # Use the extended color palette
    na.value = "white", # Set NA values to white
    guide = guide_legend( # Customize legend guide
      direction = "horizontal", # Position the legend horizontally
      title.position = "top", # Title on top
      title.hjust = 0.5, # Center the title
      label.position = "bottom", # Position the labels at the bottom
      nrow = 1, # Arrange legend items in a single row
      keywidth = unit(0.9, "cm"), # Legend key width
      keyheight = unit(0.3, "cm") # Legend key height
    )
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_blank(),
    legend.position = "none", # Hide the legend
    plot.margin = margin(0, 0, 0, 0) # Remove plot margins
  )

# Create the x and y bar plots for risk distribution
x_bar <- create_x_barplot(final_result_risk)
y_bar <- create_y_barplot(final_result_risk)

# Align the plots (longitude bar plot, main map, and latitude bar plot)
aligned_plots <- align_plots(
  x_bar, main_map,
  align = "hv", axis = "lrtb"
)

# Combine the plots into one final plot with the map in the center, x bar plot at the bottom, and y bar plot on the right side
pdf("/path/to/data/Final_3_global_mapping_risk_lonlat.pdf", family = "ArialMT")
ggdraw() +
  draw_plot(aligned_plots[[2]], x = 0, y = 0.1, width = 0.8, height = 0.8) + # Main map
  draw_plot(aligned_plots[[1]], x = 0, y = 0.2, width = 0.8, height = 0.1) + # x bar plot (bottom)
  draw_plot(y_bar, x = 0.78, y = 0.277, width = 0.08, height = 0.44) + # y bar plot (right side)
  draw_plot_label(
    c("Longitude (X)", "Latitude (Y)"),
    x = c(0.3, 0.9),
    y = c(0.2, 0.6),
    angle = c(0, -90),
    size = 10
  )
dev.off()

png("/path/to/data/Final_3_global_mapping_risk_lonlat.png", family = "ArialMT")
ggdraw() +
  draw_plot(aligned_plots[[2]], x = 0, y = 0.1, width = 0.8, height = 0.8) + # Main map
  draw_plot(aligned_plots[[1]], x = 0, y = 0.2, width = 0.8, height = 0.1) + # x bar plot (bottom)
  draw_plot(y_bar, x = 0.78, y = 0.277, width = 0.08, height = 0.44) + # y bar plot (right side)
  draw_plot_label(
    c("Longitude (X)", "Latitude (Y)"),
    x = c(0.3, 0.9),
    y = c(0.2, 0.6),
    angle = c(0, -90),
    size = 10
  )
dev.off()
