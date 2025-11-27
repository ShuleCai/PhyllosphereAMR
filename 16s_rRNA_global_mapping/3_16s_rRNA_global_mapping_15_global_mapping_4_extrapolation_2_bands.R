###### Extrapolation - Range of bands - Combined ######
# 1. Combine prediction data and feature data
training_agri <- read.csv("/home/data/t010622/namespace/phyllo/16S/ML/agri_models_x_2_feature_selected_iqr.csv")
training_nonagri <- read.csv("/home/data/t010622/namespace/phyllo/16S/ML/nonagri_models_x_2_feature_selected_iqr.csv")
res_extract_imputed <- read.csv("/home/data/t010622/namespace/phyllo/16S/ML/models_x_1_original.csv")
res_bands <- res_extract_imputed %>% dplyr::select(union(names(training_agri), names(training_nonagri)))

# 2. Stack raster data for combined bands
tif_stack_crop_filled <- stack(filled_stack, rs_landcover, rs_agri2010)
sub_stack_combined <- tif_stack_crop_filled[[names(res_bands)]] # 73 bands

# 3. Convert raster stack to data frame and filter using xy_both_filtered2
preddf_combined <- as.data.frame(sub_stack_combined, xy = TRUE, na.rm = TRUE)
preddf_combined_f <- xy_both_filtered2 %>% left_join(preddf_combined, by = c("x", "y"))

# 4. Compute the valid range for each feature
feature_ranges <- sapply(res_bands, function(x) {
  c(min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
})

# 5. Check if the prediction values fall within the valid range
in_range_matrix <- matrix(NA, nrow = nrow(preddf_combined_f), ncol = ncol(res_bands))

for (i in 1:ncol(res_bands)) {
  feature_name <- names(res_bands)[i]
  min_val <- feature_ranges["min", feature_name]
  max_val <- feature_ranges["max", feature_name]

  in_range_matrix[, i] <- as.integer(
    preddf_combined_f[[feature_name]] >= min_val &
      preddf_combined_f[[feature_name]] <= max_val
  )
}

# 6. Calculate the percentage of variables within their valid range for each pixel
pixel_percent_in_range <- rowMeans(in_range_matrix) * 100

# 7. Bin the percentage values into custom categories
breaks_vector <- c(0, seq(50, 100, by = 5), 105)
level_labels <- c(
  "<50%", "50-55%", "55-60%", "60-65%", "65-70%",
  "70-75%", "75-80%", "80-85%", "85-90%", "90%-95%", "95-100%", "100%"
)
binned_data <- tibble(percent_in_range = pixel_percent_in_range) %>%
  mutate(range_bin = cut(percent_in_range,
    breaks = breaks_vector,
    include.lowest = TRUE, right = FALSE,
    labels = level_labels
  ))

# 8. Calculate percentage of pixels in each bin
result <- binned_data %>%
  group_by(range_bin) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# 9. Ensure all bins have corresponding rows (even if count is 0)
all_bins <- tibble(range_bin = factor(level_labels, levels = level_labels))
result <- all_bins %>%
  left_join(result, by = "range_bin") %>%
  replace_na(list(count = 0, percentage = 0))

# 10. Create bar plot to visualize the distribution of coefficient variation
ggplot(result, aes(x = range_bin, y = percentage)) +
  geom_bar(stat = "identity", fill = "grey40", width = 0.8) +
  labs(x = "Percentage of bands within interpolation range", y = "Percentage of terrestrial pixels") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 0), panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_flip()

ggsave("~/namespace/phyllo/16S/ML/final/extrapolation_1_bands_percentage_boxplot.pdf", family = "ArialMT", width = 5, height = 6)
write.csv(result, "~/namespace/phyllo/16S/ML/final/extrapolation_1_bands_percentage_boxplot.csv", row.names = F)

# 11. Create a map of the percentage of bands within the valid range
mapdata_bands <- data.frame(x = preddf_combined_f$x, y = preddf_combined_f$y, in_range_pct = pixel_percent_in_range)

ggplot() +
  geom_raster(data = mapdata_bands, aes(x = x, y = y, fill = in_range_pct), na.rm = FALSE) +
  coord_sf(ylim = c(-60, 90), expand = FALSE, crs = "+proj=longlat") +
  scale_fill_viridis_c(
    name = "Percentage of bands within interpolation range",
    option = "D", direction = -1, na.value = "lightgray",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(6, "cm"),
      barheight = unit(0.4, "cm")
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, vjust = 1),
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

ggsave("~/namespace/phyllo/16S/ML/final/extrapolation_1_bands_percentage_worldmap.pdf", family = "ArialMT")
saveRDS(mapdata_bands, "~/namespace/phyllo/16S/ML/final/extrapolation_1_bands_percentage_worldmap_mapdata.rds")
