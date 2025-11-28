# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(sf)
library(viridis)
library(purrr)

# Load final prediction results (combined from different runs)
predictions_all <- readRDS("/path/to/data/Final_1_risk_abun_pred_mean_sd_run_all.rds")

# Merge all prediction files into a single dataframe (arranging by MAE)
predictions_1k2k <- readRDS("/path/to/data/Final_1_risk_abun_pred_mean_sd_run1k_2k_ranks.rds")
predictions_3k4k <- readRDS("/path/to/data/Final_1_risk_abun_pred_mean_sd_run3k_4k_ranks.rds")
predictions_5k6k <- readRDS("/path/to/data/Final_1_risk_abun_pred_mean_sd_run5k_6k_ranks.rds")
predictions_6k7k <- readRDS("/path/to/data/Final_1_risk_abun_pred_mean_sd_run6k_7k_ranks.rds")
predictions_7k8k <- readRDS("/path/to/data/Final_1_risk_abun_pred_mean_sd_run7k_8k_ranks.rds")
predictions_8k9k <- readRDS("/path/to/data/Final_1_risk_abun_pred_mean_sd_run8k_9k_ranks.rds")
predictions_9k10k <- readRDS("/path/to/data/Final_1_risk_abun_pred_mean_sd_run9k_10k_ranks.rds")
predictions_10k11k <- readRDS("/path/to/data/Final_1_risk_abun_pred_mean_sd_run10k_11k_ranks.rds")

# Combine predictions and arrange them by MAE
predictions_all <- rbind(
  predictions_all,
  predictions_1k2k,
  predictions_3k4k,
  predictions_5k6k,
  predictions_6k7k,
  predictions_7k8k,
  predictions_8k9k,
  predictions_9k10k
) %>%
  arrange(mae)

# Save the combined predictions for future use
saveRDS(predictions_all, "/path/to/data/Final_1_risk_abun_pred_mean_sd_run_all.rds")

# 1. Compute robust statistics (mean, median, SD, IQR, MAD) for each pixel
final_result_used <- combined_predictions %>%
  group_by(x, y) %>%
  summarize(
    mean_predict_2010 = mean(pred_2010, na.rm = TRUE),
    sd_predict_2010 = sd(pred_2010, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Plot coefficient of variation (CV) map
plot_data <- final_result_used %>%
  dplyr::mutate(cv_sd_m = ifelse(sd_predict_2010 > 1, 1, sd_predict_2010))

# Create the variation heatmap
ggplot() +
  geom_raster(data = plot_data, aes(x = x, y = y, fill = cv_sd_m), na.rm = FALSE) +
  coord_sf(ylim = c(-60, 90), expand = FALSE, crs = "+proj=longlat") +
  scale_fill_viridis_c(
    name = "Ensemble Coefficient of Variation",
    option = "D",
    direction = -1,
    na.value = "lightgray",
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

ggsave("/path/to/figure/extrapolation_0_coefficient_variation_worldmap.pdf", family = "ArialMT")

# 3. Bin coefficient variation values into specific ranges
breaks_vector <- c(0, seq(20, 100, by = 10))
level_labels <- c("<20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", ">90%")

binned_data <- tibble(percent_in_range = 100 * plot_data$cv_sd_m) %>%
  mutate(range_bin = cut(percent_in_range,
    breaks = breaks_vector,
    include.lowest = TRUE, right = FALSE,
    labels = level_labels
  ))

# Calculate the percentage of pixels in each bin
result <- binned_data %>%
  group_by(range_bin) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# Ensure all bins are represented (even those with 0 pixels)
all_bins <- tibble(range_bin = factor(level_labels, levels = level_labels %>% rev()))
result <- all_bins %>%
  left_join(result, by = "range_bin") %>%
  replace_na(list(count = 0, percentage = 0))

# 4. Create bar plot for coefficient variation distribution
ggplot(result, aes(x = range_bin, y = percentage)) +
  geom_bar(stat = "identity", fill = "grey40", width = 0.8) +
  labs(x = "Percentage of Coefficient Variation", y = "Percentage of Terrestrial Pixels") +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 0),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_flip()

ggsave("/path/to/figure/extrapolation_0_coefficient_variation_boxplot.pdf", family = "ArialMT", width = 5, height = 6)

# Save the results of the binning
write.csv(result, "/path/to/data/extrapolation_0_coefficient_variation_boxplot.csv", row.names = F)
