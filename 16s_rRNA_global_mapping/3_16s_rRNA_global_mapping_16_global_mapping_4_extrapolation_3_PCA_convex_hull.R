library(tidyverse)
library(FactoMineR)
library(geometry)
library(sf)

# Load training datasets for agricultural and non-agricultural data
training_agri <- read.csv("/path/to/data/agri_models_x_2_feature_selected_iqr.csv")
training_nonagri <- read.csv("/path/to/data/nonagri_models_x_2_feature_selected_iqr.csv")

# Load the imputed feature data
res_extract_imputed <- read.csv("/path/to/data/models_x_1_original.csv")

# Combine selected feature names from agricultural and non-agricultural datasets
res_bands <- res_extract_imputed %>% dplyr::select(union(names(training_agri), c()))

# 1. Data Preparation and PCA Analysis
# Perform PCA on the combined feature dataset (res_bands)
pca_result <- PCA(res_bands, scale.unit = TRUE, graph = FALSE, ncp = 30)

# Extract eigenvalues and calculate cumulative variance
eigenvalues <- pca_result$eig
cumulative_variance <- cumsum(eigenvalues[, 2]) / 100 # Cumulative variance explained
k <- which.max(cumulative_variance >= 0.85) # Select the number of components that explain 85% of the variance

# 2. Generate all pairs of principal components
pairs <- combn(1:k, 2, simplify = FALSE)
num_pairs <- length(pairs) # Number of pairs of principal components

# 3. Compute Convex Hull for each pair of principal components
convex_hulls <- lapply(pairs, function(idx) {
  pair_scores <- pca_result$ind$coord[, idx] # Get the PCA scores for the selected pair
  convhulln(pair_scores, "FA") # Compute the convex hull in FA space
})

# 4. Prepare prediction data based on the selected principal components
pca_pred <- predict(pca_result, preddf_combined_f[, names(res_bands)])$coord[, 1:k]

# 5. Check if each point falls within the convex hull for each pair of components
in_hull_matrix <- matrix(0, nrow = nrow(pca_pred), ncol = num_pairs)

# Loop through each pair of components
for (i in 1:num_pairs) {
  print(i)
  pair_idx <- pairs[[i]] # Get the pair of principal components
  conv_hull <- convex_hulls[[i]] # Get the convex hull for this pair

  # Loop through each point in the prediction data
  for (j in 1:nrow(pca_pred)) {
    point <- matrix(pca_pred[j, pair_idx], nrow = 1) # Extract the point's coordinates for the pair
    in_hull_matrix[j, i] <- as.integer(inhulln(conv_hull, point)) # Check if the point is inside the convex hull
  }
}

# Calculate the percentage of pairs that each point falls within the convex hull
in_hull_percentage <- rowSums(in_hull_matrix) / num_pairs * 100

# 6. Create a data frame with results
results_df <- data.frame(
  x = preddf_combined_f$x, # X-coordinate
  y = preddf_combined_f$y, # Y-coordinate
  in_hull_percentage = in_hull_percentage # Percentage of convex hulls covered
)

# 7. Create binned data for visualization
binned_data <- results_df %>%
  mutate(range_bin = cut(
    in_hull_percentage,
    breaks = seq(0, 100, by = 10), # Define breaks for bins
    include.lowest = TRUE,
    labels = c(
      "<10%", "10-20%", "20-30%", "30-40%", "40-50%",
      "50-60%", "60-70%", "70-80%", "80-90%", ">90%"
    ),
    right = FALSE
  ))

# 8. Calculate the number of pixels in each bin
result_counts <- binned_data %>%
  group_by(range_bin) %>%
  summarise(count = n(), .groups = "drop")

# 9. Calculate the percentage of pixels in each bin
result_counts <- result_counts %>%
  mutate(percentage = count / sum(count) * 100)

# 10. Ensure all bins are represented in the final result
all_bins <- tibble(range_bin = factor(
  c(
    "<10%", "10-20%", "20-30%", "30-40%", "40-50%",
    "50-60%", "60-70%", "70-80%", "80-90%", ">90%"
  ),
  levels = c(
    "<10%", "10-20%", "20-30%", "30-40%", "40-50%",
    "50-60%", "60-70%", "70-80%", "80-90%", ">90%"
  )
))

result_counts <- all_bins %>%
  left_join(result_counts, by = "range_bin") %>%
  replace_na(list(count = 0, percentage = 0)) # Replace NA values with 0

# 11. Create a horizontal bar plot of the binned data
ggplot(result_counts, aes(x = range_bin, y = percentage)) +
  geom_bar(stat = "identity", fill = "grey40", width = 0.8) +
  labs(
    subtitle = paste(
      k, "PCs (", round(cumulative_variance[k] * 100, 1), "% variance explained),",
      num_pairs, "pairs"
    ),
    x = "Percentage of convex hulls the pixels fall within",
    y = "Percentage of terrestrial pixels"
  ) +
  theme_linedraw() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.minor.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_flip() # Display the plot horizontally
ggsave("/path/to/figure/extrapolation_2_convexhull_percentage_boxplot.pdf", family = "ArialMT", width = 5, height = 6)
write.csv(result_counts, "/path/to/data/extrapolation_2_convexhull_percentage_boxplot.csv", row.names = F)

# 12. Generate a world map showing the percentage of convex hull coverage
ggplot() +
  geom_raster(
    data = results_df,
    aes(x = x, y = y, fill = in_hull_percentage),
    na.rm = TRUE
  ) +
  coord_sf(ylim = c(-60, 90), expand = FALSE, crs = "+proj=longlat") +
  scale_fill_viridis_c(
    name = "Percentage of convex hulls the pixels fall within",
    option = "D", # Use the viridis "inferno" color scheme
    direction = -1, # Reverse the color scale (high values = darker colors)
    na.value = "lightgray",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(8, "cm"),
      barheight = unit(0.4, "cm"),
      ticks.colour = "lightgray"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_text(size = 10, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(title = paste0("Proportion of PC pairs covered by convex hulls (", k, " PCs, ", num_pairs, " pairs)"))
ggsave("/path/to/figure/extrapolation_2_convexhull_percentage_worldmap.pdf", family = "ArialMT")
saveRDS(results_df, "/path/to/data/extrapolation_2_convexhull_percentage_worldmap_mapdata.rds")
