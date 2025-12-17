library(tidyverse)
library(FactoMineR)
library(geometry)
library(sf)

# 1. Data preparation and PCA analysis
pca_result <- PCA(res_bands, scale.unit = TRUE, graph = FALSE, ncp = 20)

# Extract eigenvalues and cumulative variance explained
eigenvalues <- pca_result$eig
cumulative_variance <- cumsum(eigenvalues[, 2]) / 100 # Cumulative variance ratio
k <- which.max(cumulative_variance >= 0.85) # Find the number of components explaining >= 85% variance

# 2. Generate all combinations of principal component pairs
pairs <- combn(1:k, 2, simplify = FALSE)
num_pairs <- length(pairs)
num_pairs # Number of PC pairs

# 3. Construct convex hulls for each pair of principal components
convex_hulls <- lapply(pairs, function(idx) {
  pair_scores <- pca_result$ind$coord[, idx]
  convhulln(pair_scores, "FA") # Calculate convex hulls
})

# 4. Prepare prediction data
pca_pred <- predict(pca_result, preddf_combined_f)$coord[, 1:k] # Project predictions onto PCA space

# 5. Check which points fall inside the convex hulls
in_hull_matrix <- matrix(0, nrow = nrow(pca_pred), ncol = num_pairs)

for (i in 1:num_pairs) {
  pair_idx <- pairs[[i]]
  conv_hull <- convex_hulls[[i]]

  for (j in 1:nrow(pca_pred)) {
    point <- matrix(pca_pred[j, pair_idx], nrow = 1)
    in_hull_matrix[j, i] <- as.integer(inhulln(conv_hull, point)) # Check if point is inside the convex hull
  }
}

# Calculate the proportion of points inside the convex hulls
in_hull_percentage <- rowSums(in_hull_matrix) / num_pairs * 100 # Proportion of hulls covered for each point

# 6. Create a result data frame
results_df <- data.frame(
  x = preddf_combined_f$x,
  y = preddf_combined_f$y,
  in_hull_percentage = in_hull_percentage # Store percentage in the data frame
)

# 7. Create binning data based on the calculated percentages
binned_data <- results_df %>%
  mutate(range_bin = cut(
    in_hull_percentage,
    breaks = seq(0, 100, by = 10),
    include.lowest = TRUE,
    labels = c(
      "<10%", "10-20%", "20-30%", "30-40%", "40-50%",
      "50-60%", "60-70%", "70-80%", "80-90%", ">90%"
    ),
    right = FALSE # Left-closed, right-open intervals
  ))

# 8. Calculate the count and percentage of pixels in each bin
result_counts <- binned_data %>%
  group_by(range_bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count) * 100)

# Ensure that all bins are represented in the data (even if count = 0)
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

# 9. Plot the histogram for the bin percentages (horizontal bar chart)
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
  coord_flip() # Horizontal display of bars
ggsave("/path/to/figure/extrapolation_2_convexhull_percentage_boxplot.pdf", family = "ArialMT")
write.csv(result_counts, "/path/to/data/extrapolation_2_convexhull_percentage_boxplot.csv", row.names = FALSE)

# 10. Plot the uncertainty map
ggplot() +
  geom_raster(
    data = results_df,
    aes(x = x, y = y, fill = in_hull_percentage),
    na.rm = TRUE
  ) +
  coord_sf(ylim = c(-60, 90), expand = FALSE, crs = "+proj=longlat") +
  scale_fill_viridis_c(
    name = "Percentage of convex hulls the pixels fall within",
    option = "D", # Using viridis' "inferno" scheme
    direction = -1, # Reverse color direction (high value = darker color)
    na.value = "lightgray",
    limits = c(0, 100),
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
    plot.margin = margin(10, 10, 10, 10),
    panel.background = element_rect(fill = "lightgray"),
  ) +
  labs(title = paste0("Proportion of PC pairs covered by convex hulls (", k, " PCs, ", num_pairs, " pairs)"))
ggsave("/path/to/figure/extrapolation_2_convexhull_percentage_worldmap.pdf", family = "ArialMT")
saveRDS(results_df, "/path/to/data/extrapolation_2_convexhull_percentage_worldmap_mapdata.rds")
