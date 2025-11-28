# Load WDPA points data
wdpa_pts <- readRDS("/path/to/data/WDPA/lsl_processed/pts_protected_final.rds")

# Merge predictions with WDPA data
predictions_future_wdpa <- predictions_future %>% left_join(wdpa_pts, by = c("x", "y"))

# Combine predictions and compute weighted predictions based on agricultural land cover and protection status
combined_predictions_future <- predictions_future_wdpa %>%
  dplyr::select(-run) %>%
  left_join(
    preddf_agri_future_f %>%
      dplyr::select(x, y, AgriculturalLandCover_2010, AgriculturalLandCover_2050, landcover = Global_MCD12Q1_IGBP_500m_2023_resampled),
    by = c("x", "y")
  ) %>%
  mutate(
    # Business-as-usual scenario (100% agricultural expansion)
    pred_2050_agriexpan100 = pred_agri * AgriculturalLandCover_2050 + pred_nona * (1 - AgriculturalLandCover_2050),

    # Predictions considering protected areas and land cover types
    pred_2050_protect1 = ifelse(landcover %in% 1:5 | landcover == 11 | strictly_protected_area == TRUE,
      pred_agri * AgriculturalLandCover_2010 + pred_nona * (1 - AgriculturalLandCover_2010),
      ifelse(protected_area == TRUE,
        pred_agri * (AgriculturalLandCover_2010 + AgriculturalLandCover_2050) / 2 + pred_nona * (1 - (AgriculturalLandCover_2010 + AgriculturalLandCover_2050) / 2),
        pred_agri * AgriculturalLandCover_2050 + pred_nona * (1 - AgriculturalLandCover_2050)
      )
    ),

    # Predictions for non-agricultural expansion and protected areas
    pred_2050_protect2 = ifelse(landcover %in% 1:5 | landcover == 11 | strictly_protected_area == TRUE,
      pred_nona,
      ifelse(protected_area == TRUE,
        pred_agri * AgriculturalLandCover_2010 + pred_nona * (1 - AgriculturalLandCover_2010),
        pred_agri * (AgriculturalLandCover_2010 * 0.2 + AgriculturalLandCover_2050 * 0.8) + pred_nona * (1 - (AgriculturalLandCover_2010 * 0.2 + AgriculturalLandCover_2050 * 0.8))
      )
    )
  )

# Merge with 2010 predictions
combined_predictions_merge <- combined_predictions_future %>%
  inner_join(combined_predictions %>% dplyr::select(x, y, iteration, pred_2010),
    by = c("x", "y", "iteration")
  )

# Get continent information using world map (SF object)
library(sf)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create SF object for unique points in the predictions data
unique_points <- unique(combined_predictions_merge[, c("x", "y")])
points_sf <- sf::st_as_sf(unique_points, coords = c("x", "y"), crs = 4326)

# Join nearest continent information
joined <- sf::st_join(points_sf, world["continent"], join = st_nearest_feature) %>%
  mutate(x = st_coordinates(.)[, 1], y = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()

# Merge continent information with predictions
combined_predictions_merge_continent <- combined_predictions_merge %>%
  dplyr::left_join(joined, by = c("x", "y")) %>%
  dplyr::filter(!is.na(continent)) # Remove points without continent information

# Add a global data entry (global as an additional continent)
global_data <- combined_predictions_merge %>%
  mutate(continent = "Global") %>%
  bind_rows(combined_predictions_merge_continent)

# Calculate mean predictions for each continent/region per iteration
continent_means <- global_data %>%
  group_by(continent, iteration) %>%
  summarise(
    mean_pred_2010 = mean(pred_2010, na.rm = TRUE),
    mean_pred_agriexpan100 = mean(pred_2050_agriexpan100, na.rm = TRUE),
    mean_pred_protect1 = mean(pred_2050_protect1, na.rm = TRUE),
    mean_pred_protect2 = mean(pred_2050_protect2, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate growth rates for each continent/region and scenario
continent_growth_rates <- continent_means %>%
  mutate(
    growth_agriexpan100 = (mean_pred_agriexpan100 - mean_pred_2010) / mean_pred_2010,
    growth_protect1 = (mean_pred_protect1 - mean_pred_2010) / mean_pred_2010,
    growth_protect2 = (mean_pred_protect2 - mean_pred_2010) / mean_pred_2010
  ) %>%
  dplyr::select(continent, iteration, starts_with("growth"))

# Compute summary statistics (mean, CI) for the growth rates
stats_data <- continent_growth_rates %>%
  group_by(continent) %>%
  summarise(
    across(
      c(growth_agriexpan100, growth_protect1, growth_protect2),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        ci_low = ~ {
          clean_data <- na.omit(.x)
          mean_val <- mean(clean_data, na.rm = TRUE)
          n_valid <- length(clean_data)
          mean_val - 1.96 * sd(clean_data, na.rm = TRUE) / sqrt(n_valid)
        },
        ci_high = ~ {
          clean_data <- na.omit(.x)
          mean_val <- mean(clean_data, na.rm = TRUE)
          n_valid <- length(clean_data)
          mean_val + 1.96 * sd(clean_data, na.rm = TRUE) / sqrt(n_valid)
        }
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  pivot_longer(
    cols = -continent,
    names_to = c("scenario", ".value"),
    names_pattern = "growth_(.+)_(mean|ci_low|ci_high)"
  )

# Set the continent order with 'Global' at the end
continent_order <- c(setdiff(unique(stats_data$continent), "Global"), "Global")
stats_data$continent <- factor(stats_data$continent, levels = continent_order)

# Create labels for the scenarios
scenario_labels <- c(
  "agriexpan100" = "Business-as-usual",
  "protect1" = "Partially regulated",
  "protect2" = "Strongly regulated"
)

# Add scenario labels to the data
plot_data <- stats_data %>%
  mutate(scenario_label = factor(scenario_labels[scenario],
    levels = scenario_labels
  ))

# Custom colors for continents
manual_colors <- c(
  "Africa" = "#da4b50",
  "Asia" = "#f78f5b",
  "Europe" = "#fed489",
  "North America" = "#ffffc3",
  "Oceania" = "#d6eea1",
  "South America" = "#8acfa9",
  "Global" = "#4196b5"
)

# Create bar chart for projected change across continents
combined_plot <- ggplot(plot_data, aes(x = scenario_label, y = 100 * mean, fill = continent)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    width = 0.85,
    alpha = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = 100 * ci_low, ymax = 100 * ci_high),
    position = position_dodge(width = 0.9),
    width = 0.4,
    color = "black",
    linewidth = 0.5
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Projected change (%)",
    fill = NULL
  ) +
  scale_fill_manual(values = manual_colors) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  guides(fill = guide_legend(
    title.position = "top",
    title.hjust = 0.5,
    keywidth = 0.8,
    keyheight = 0.8
  ))

# Save the bar chart for projected changes across continents and scenarios
ggsave(combined_plot, filename = "/path/to/figure/Future2050_change_scenarios_barchart.pdf", family = "ArialMT", width = 7, height = 4)

# Export the data for further analysis or use
plot_data %>%
  as.data.frame() %>%
  write.csv("/path/to/data/Future2050_change_scenarios_barchart.csv")

# Adjust scenario labels and order for the plot
scenario_levels <- c("Business-as-usual", "Partially regulated", "Strongly regulated")
plot_data$scenario_label <- factor(plot_data$scenario_label, levels = scenario_levels)

# Adjust the data for plotting with custom transformations and calculations
plot_data_modified <- plot_data %>%
  group_by(continent) %>%
  mutate(
    protect1_value = mean[scenario == "protect1"],
    protect2_value = mean[scenario == "protect2"],
    mean_new = case_when(
      scenario == "agriexpan100" ~ mean - protect1_value,
      scenario == "protect1" ~ mean - protect2_value,
      scenario == "protect2" ~ mean
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-protect2_value, -protect1_value) %>% # Remove temporary columns
  mutate(continent = factor(continent, levels = c("Africa", "South America", "Asia", "Global", "North America", "Oceania", "Europe"))) %>%
  dplyr::select(continent, mean_new, scenario_label) %>%
  rbind(data.frame(continent = unique(plot_data$continent), mean_new = 0.1, scenario_label = "FILL"))

# Create label data for positioning text in the plot
label_data <- plot_data_modified %>%
  group_by(continent) %>%
  summarise(max_mean = sum(mean_new)) %>%
  mutate(label_y = max_mean * 1.05)

# Create a polar plot (rose chart) to visualize the data across continents by scenario
ggplot(plot_data_modified, aes(x = continent, y = mean_new, fill = scenario_label)) +
  geom_col(
    position = "stack", # Stack the bars
    width = 0.95, # Control the width of each bar
    color = "black", # Set the border color of bars
    linewidth = 0.3 # Control the thickness of the borders
  ) +
  coord_polar(theta = "x", start = 0) + # Create the circular (rose) chart
  scale_fill_manual(
    values = c(
      "Strongly regulated" = "#ffffd8", # Light yellow - Strong regulation
      "Partially regulated" = "#87b4c9", # Light blue - Partial regulation
      "Business-as-usual" = "#d78687" # Light red - Business as usual
    ),
    name = "Scenarios",
    guide = guide_legend(reverse = TRUE) # Reverse the legend order for clarity
  ) +
  geom_text(
    data = label_data,
    aes(x = continent, y = label_y, label = continent),
    inherit.aes = FALSE,
    size = 3.5,
    color = "black",
    fontface = "bold"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    plot.caption = element_text(hjust = 0.5, color = "gray50", margin = margin(t = 10)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.box.spacing = unit(0.2, "cm")
  )

# Save the rose chart (polar plot) as a PDF
ggsave("/path/to/figure/Future2050_change_scenarios_roses.pdf")
