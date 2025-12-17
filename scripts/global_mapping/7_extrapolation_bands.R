library(tidyverse)
library(sf)
library(rnaturalearth)
library(cowplot)
library(scales)

# Read input data files for agricultural and non-agricultural models
training_agri <- read.csv("/path/to/data/agri_models_x_2_feature_selected.csv")
training_nonagri <- read.csv("/path/to/data/nonagri_models_x_2_feature_selected.csv")

# Select columns based on the union of the agricultural and non-agricultural models
res_bands <- res_extract_imputed %>% dplyr::select(union(names(training_agri), names(training_nonagri)))
sub_stack_combined <- tif_stack_crop[[names(res_bands)]] # 73 bands in the data

# Process the data and create a data frame from the raster stack
preddf_combined <- as.data.frame(sub_stack_combined, xy = TRUE, na.rm = TRUE)
preddf_combined_f <- xy_both_filtered2 %>% left_join(preddf_combined, by = c("x", "y"))

# 1. Calculate the valid range for each variable
feature_ranges <- sapply(res_bands, function(x) {
  c(
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
})

# 2. Check whether each value in the prediction data frame is within its corresponding variable's range
in_range_matrix <- matrix(NA,
  nrow = nrow(preddf_combined_f),
  ncol = ncol(res_bands)
)

for (i in 1:ncol(res_bands)) {
  feature_name <- names(res_bands)[i]
  min_val <- feature_ranges["min", feature_name]
  max_val <- feature_ranges["max", feature_name]

  in_range_matrix[, i] <- as.integer(
    preddf_combined_f[[feature_name]] >= min_val &
      preddf_combined_f[[feature_name]] <= max_val
  )
}

# 3. Calculate the proportion of variables within range for each pixel
pixel_percent_in_range <- rowMeans(in_range_matrix) * 100

# 4. Create custom binning for the data
# Define binning rules: one bin for values <50%, and 5% bins for values >50%
breaks_vector <- c(seq(0, 100, by = 10))
level_labels <- c(
  "<10%", "10-20%", "20-30%", "30-40%", "40-50%",
  "50-60%", "60-70%", "70-80%", "80-90%", ">90%"
)

binned_data <- tibble(percent_in_range = pixel_percent_in_range) %>%
  mutate(range_bin = cut(percent_in_range,
    breaks = breaks_vector,
    include.lowest = TRUE,
    right = FALSE, # Left-closed, right-open intervals
    labels = level_labels
  ))

# 5. Calculate the percentage of pixels in each bin
result <- binned_data %>%
  group_by(range_bin) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# Ensure all bins have corresponding rows (even if the count is 0)
all_bins <- tibble(range_bin = factor(level_labels, levels = level_labels))
result <- all_bins %>%
  left_join(result, by = "range_bin") %>%
  replace_na(list(count = 0, percentage = 0))

# 6. Create a bar plot of the percentage distribution
ggplot(result, aes(x = range_bin, y = percentage)) +
  geom_bar(stat = "identity", fill = "grey40", width = 0.8) +
  labs(
    x = "Percentage of bands within interpolation range",
    y = "Percentage of terrestrial pixels"
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 0),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_flip() # Display the bar chart horizontally
ggsave("/path/to/figure/extrapolation_1_bands_percentage_boxplot.pdf", family = "ArialMT")
write.csv(result, "/path/to/data/extrapolation_1_bands_percentage_boxplot.csv", row.names = FALSE)

# 7. Create a bands uncertainty map
mapdata_bands <- data.frame(x = preddf_combined_f$x, y = preddf_combined_f$y, in_range_pct = pixel_percent_in_range)

# Plot the uncertainty map
ggplot() +
  # Plot the uncertainty heatmap (continuous variable)
  geom_raster(
    data = mapdata_bands,
    aes(x = x, y = y, fill = in_range_pct),
    na.rm = FALSE
  ) +
  # Set the coordinate system (limit to land areas)
  coord_sf(ylim = c(-60, 90), expand = FALSE, crs = "+proj=longlat") +
  # Use viridis continuous color scale
  scale_fill_viridis_c(
    name = "Percentage of bands within interpolation range", # Title for the color legend
    option = "D", # Use viridis "D" option
    direction = -1, # Reverse the color scale (higher percentage = darker color)
    na.value = "lightgray", # Set NA values as light gray
    guide = guide_colorbar(
      direction = "horizontal", # Horizontal color legend
      title.position = "top", # Title position at the top
      title.hjust = 0.5, # Title alignment
      barwidth = unit(6, "cm"), # Legend bar width
      barheight = unit(0.4, "cm"), # Legend bar height
      ticks.colour = "lightgray" # Set tick color to light gray
    )
  ) +
  # Set the theme
  theme_void() +
  theme(
    legend.position = "bottom", # Position the legend at the bottom
    legend.justification = "center", # Center the legend
    legend.title = element_text(size = 10, vjust = 1), # Title font size for the legend
    plot.margin = margin(10, 10, 10, 10),
    panel.background = element_rect(fill = "lightgray"),
    plot.title = element_text(hjust = 0.5, size = 12)
  ) +
  labs(title = "Proportion of bands within interpolation range")
ggsave("/path/to/figure/extrapolation_1_bands_percentage_worldmap.pdf", family = "ArialMT")
saveRDS(mapdata_bands, "/path/to/data/extrapolation_1_bands_percentage_worldmap_mapdata.rds")
