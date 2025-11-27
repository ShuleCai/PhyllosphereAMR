library(tidyverse)
library(sf)
library(rnaturalearth)
library(cowplot)
library(scales)

# Calculate the proportion of pixels with variation > 0.1
sum(final_result$variation > 0.1) / nrow(final_result) # Output: 0.0006281026

# 4. Create a variation map
# Plot the map
ggplot() +
  # Plot the uncertainty heatmap (continuous variable)
  geom_raster(
    data = final_result,
    aes(x = x, y = y, fill = variation),
    na.rm = FALSE
  ) +
  # Set the coordinate system (limit to land areas)
  coord_sf(ylim = c(-60, 90), expand = FALSE, crs = "+proj=longlat") +
  # Use the "viridis" continuous color scale
  scale_fill_viridis_c(
    name = "Ensemble coefficient of variation", # Title for the color legend
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
  labs(title = "Ensemble coefficient of variation")
ggsave("/path/to/project/ML/final/extrapolation_0_coefficient_variation_worldmap.pdf", family = "ArialMT")

# 5. Create custom bins for the data
# Define binning rules: one bin for <50%, and then 5% bins for values >50%
breaks_vector <- c(seq(0, 0.5, by = 0.05), 1)
level_labels <- c(
  "<5%", "5-10%", "10-15%", "15-20%", "20-25%", "25-30%",
  "30-35%", "35-40%", "40-45%", "45-50%", ">50%"
)

# Binning the data based on the variation percentage
binned_data <- tibble(percent_in_range = final_result$variation) %>%
  mutate(range_bin = cut(percent_in_range,
    breaks = breaks_vector,
    include.lowest = TRUE,
    right = FALSE, # Left-closed, right-open intervals
    labels = level_labels
  ))

# 6. Calculate the percentage of pixels in each bin
result <- binned_data %>%
  group_by(range_bin) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# Ensure all bins have corresponding rows (even if the count is 0)
all_bins <- tibble(range_bin = factor(level_labels, levels = level_labels %>% rev()))
result <- all_bins %>%
  left_join(result, by = "range_bin") %>%
  replace_na(list(count = 0, percentage = 0))

# 7. Create a bar plot of the percentage distribution
ggplot(result, aes(x = range_bin, y = percentage)) +
  geom_bar(stat = "identity", fill = "grey40", width = 0.8) +
  labs(
    x = "Percentage of coefficient variation",
    y = "Percentage of terrestrial pixels"
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 0),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_flip() # Flip the bars to display horizontally
ggsave("/path/to/project/ML/final/extrapolation_0_coefficient_variation_boxplot.pdf", family = "ArialMT")
write.csv(result, "/path/to/project/ML/final/extrapolation_0_coefficient_variation_boxplot.csv", row.names = FALSE)
