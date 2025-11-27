library(h2o)
h2o.init(nthreads = 5, max_mem_size = "80G")

# Load the final models for non-agricultural and agricultural data
aml_nona_leader <- readRDS("/path/to/data/nonagri_final_model.rds")
aml_agri_leader <- readRDS("/path/to/data/agri_final_model.rds")

# Load the future ensemble stack and other relevant raster layers
future_stack <- rast("/path/to/data/future2050/future_ensembled_stack_data.tif")
rs_landcover <- raster("/path/to/data/Global_MCD12Q1_IGBP_500m_2023_resampled.tif")
rs_agri2010 <- raster("/path/to/data/AgriculturalLandCover_2010_resampled.tif")
rs_agri2050 <- raster("/path/to/data/AgriculturalLandCover_2050_resampled.tif")
rs_Koppen <- raster("/path/to/data/Bio_199_Koppen.tif")

# Combine the future stack with additional raster layers
tif_stack_crop_future <- c(future_stack, rast(rs_landcover), rast(rs_agri2010), rast(rs_agri2050), rast(rs_Koppen))

# Subset the raster stack for agricultural and non-agricultural features
sub_stack_agri_future <- tif_stack_crop_future[[c(
  names(read.csv("/path/to/data/agri_models_x_2_feature_selected_iqr.csv")),
  "AgriculturalLandCover_2010", "Bio_199_Koppen", "AgriculturalLandCover_2050", "Global_MCD12Q1_IGBP_500m_2023_resampled"
)]]
sub_stack_nona_future <- tif_stack_crop_future[[c(
  names(read.csv("/path/to/data/nonagri_models_x_2_feature_selected_iqr.csv")),
  "AgriculturalLandCover_2010", "Bio_199_Koppen", "AgriculturalLandCover_2050", "Global_MCD12Q1_IGBP_500m_2023_resampled"
)]]

# Convert the selected raster stack subsets into data frames for prediction
preddf_agri_future <- as.data.frame(sub_stack_agri_future, xy = TRUE, na.rm = TRUE)
preddf_nona_future <- as.data.frame(sub_stack_nona_future, xy = TRUE, na.rm = TRUE)

# Filter by climate zone and land cover to find overlapping coordinates between agricultural and non-agricultural data
xy_both_future <- inner_join(data.frame(x = preddf_agri_future$x, y = preddf_agri_future$y),
  data.frame(x = preddf_nona_future$x, y = preddf_nona_future$y),
  by = c("x", "y")
) %>%
  distinct(x, y) %>%
  left_join(preddf_agri_future %>% dplyr::select(x, y, Bio_199_Koppen, landcover = Global_MCD12Q1_IGBP_500m_2023_resampled))

# Calculate overlap percentage
nrow(xy_both_future) / nrow(preddf_agri_future) # 96.26%

# Check for missing Koppen climate zone values
coor_16S <- read.csv("/path/to/data/p8_16s_ML_source_data.csv") %>% dplyr::select(lon = longitude, lat = latitude)
coordinates(coor_16S) <- c("lon", "lat")
cbind(coor_16S %>% as.data.frame(), Koppen = raster::extract(rs_Koppen, coor_16S)) %>%
  filter(is.na(Koppen)) %>%
  View()

Koppen_existed <- c(raster::extract(rs_Koppen, coor_16S) %>% unique(), 1, 2) # NA values found are 1, 2
Koppen_existed

# First filtration: remove land cover values from specific categories (e.g., snow, barren, water bodies)
xy_both_future_filtered1 <- xy_both_future %>% filter(!landcover %in% c(15, 16, 17))
nrow(xy_both_future_filtered1) / nrow(xy_both_future) # 87.50%

# Second filtration: filter by Koppen Climate Zone
xy_both_future_filtered2 <- xy_both_future_filtered1 %>% filter(Bio_199_Koppen %in% Koppen_existed)
nrow(xy_both_future_filtered2) / nrow(xy_both_future) # 80.85%

# Join the filtered data with the original prediction data
preddf_agri_future_f <- xy_both_future_filtered2 %>% left_join(preddf_agri_future, by = c("x", "y"))
preddf_nona_future_f <- xy_both_future_filtered2 %>% left_join(preddf_nona_future, by = c("x", "y"))

# Convert the filtered data frames into H2O format
h2o_preddf_future_agri <- as.h2o(preddf_agri_future_f)
h2o_preddf_future_nona <- as.h2o(preddf_nona_future_f)

# Extract GBM model parameters (useful for analysis and reproduction of model results)
extract_gbm_params <- function(model) {
  keep_params <- c(
    "ntrees", "max_depth", "learn_rate", "sample_rate", "col_sample_rate",
    "col_sample_rate_per_tree", "distribution", "nfolds", "fold_assignment",
    "score_tree_interval", "min_rows", "stopping_metric", "stopping_tolerance",
    "min_split_improvement", "keep_cross_validation_predictions",
    "auto_rebalance", "gainslift_bins", "calibration_method"
  )

  params <- model@allparameters
  model_params <- list()

  for (param in keep_params) {
    if (!is.null(params[[param]])) {
      model_params[[param]] <- params[[param]]
    }
  }

  model_params[c("x", "y")] <-
    list(model@parameters$x, model@parameters$y, model@parameters$training_frame)

  return(model_params)
}

# Extract XGBoost parameters (useful for analyzing the training process)
extract_xgb_params <- function(model) {
  original_params <- model@params$actual
  original_params["seed"] <- NULL
  original_params["response_column"] <- NULL
  original_params["fold_assignment"] <- NULL
  original_params["y"] <- "Enterobacteriaceae"
  original_params[["x"]] <- model@allparameters$x
  model_params <- original_params[4:length(original_params)]
  return(model_params)
}

# Extract model parameters for agricultural and non-agricultural models
params_agri <- extract_xgb_params(aml_agri_leader)
params_nona <- extract_xgb_params(aml_nona_leader)

# Build train and test frames for model training
split_agri <- as.h2o(cbind(
  read.csv("/path/to/data/agri_models_x_2_feature_selected_iqr.csv"),
  data.frame(Enterobacteriaceae = read.csv("/path/to/data/agri_models_y_iqr.csv")$Enterobacteriaceae)
)) %>%
  h2o.splitFrame(ratios = 0.9, seed = 33)
train_agri <- split_agri[[1]] # 90% training set
test_agri <- split_agri[[2]] # 10% test set

split_nona <- as.h2o(cbind(
  read.csv("/path/to/data/nonagri_models_x_2_feature_selected_iqr.csv"),
  data.frame(Enterobacteriaceae = read.csv("/path/to/data/nonagri_models_y_iqr.csv")$Enterobacteriaceae)
)) %>%
  h2o.splitFrame(ratios = 0.9, seed = 9)
train_nona <- split_nona[[1]] # 90% training set
test_nona <- split_nona[[2]] # 10% test set

# Retraining and prediction function for future scenarios
retrain_and_predict_future <- function(seed) {
  print(seed)

  # Train models using the extracted parameters
  nona_model <- do.call(h2o.xgboost, c(params_nona, list(seed = seed, training_frame = train_nona)))
  agri_model <- do.call(h2o.xgboost, c(params_agri, list(seed = seed, training_frame = train_agri)))

  # Make predictions
  pred_nona <- as.data.frame(h2o.predict(nona_model, h2o_preddf_future_nona))$predict
  pred_agri <- as.data.frame(h2o.predict(agri_model, h2o_preddf_future_agri))$predict

  # Ensure no negative predictions by setting them to 0
  pred_nona <- ifelse(pred_nona < 0, 0, pred_nona)
  pred_agri <- ifelse(pred_agri < 0, 0, pred_agri)

  # Clean up memory by removing models
  h2o.rm(nona_model)
  gc()

  # Return prediction results as a data frame
  data.frame(
    x = preddf_nona_future_f$x, y = preddf_nona_future_f$y,
    pred_agri = pred_agri,
    pred_nona = pred_nona, iteration = seed
  )
}

# Execute retraining and prediction for 100 iterations
n_iterations <- 100
predictions_future <- purrr::map_df(1:n_iterations, retrain_and_predict_future, .id = "run")

# Combine predictions and compute weighted predictions
combined_predictions_future <- predictions_future %>%
  dplyr::select(-run) %>%
  left_join(preddf_agri_future_f %>% dplyr::select(x, y, AgriculturalLandCover_2010, AgriculturalLandCover_2050), by = c("x", "y")) %>%
  mutate(
    pred_2050_agriexpan0 = ifelse(AgriculturalLandCover_2010 < AgriculturalLandCover_2050,
      pred_agri * AgriculturalLandCover_2010 + pred_nona * (1 - AgriculturalLandCover_2010),
      pred_agri * AgriculturalLandCover_2050 + pred_nona * (1 - AgriculturalLandCover_2050)
    ),
    pred_2050_agriexpan100 = pred_agri * AgriculturalLandCover_2050 + pred_nona * (1 - AgriculturalLandCover_2050),
    pred_2050_agriexpan50 = ifelse(AgriculturalLandCover_2010 < AgriculturalLandCover_2050,
      pred_agri * (AgriculturalLandCover_2010 + AgriculturalLandCover_2050) / 2 + pred_nona * (1 - (AgriculturalLandCover_2010 + AgriculturalLandCover_2050) / 2),
      pred_agri * AgriculturalLandCover_2050 + pred_nona * (1 - AgriculturalLandCover_2050)
    )
  )

# Calculate mean and standard deviation for each location in the future dataset
final_result_future <- combined_predictions_future %>%
  group_by(x, y) %>%
  summarize(
    mean_predict_future_agri2010 = mean(pred_2050_agriexpan0, na.rm = TRUE),
    mean_predict_future_agri2050 = mean(pred_2050_agriexpan100, na.rm = TRUE),
    .groups = "drop"
  )

# Save the results to CSV files
combined_predictions_future %>% write.csv("/path/to/data/Future2050_1_risk_abun_pred_mean_sd_run100.csv")
final_result_future %>% write.csv("/path/to/data/Future2050_1_risk_abun_pred_mean_sd.csv")

# Plot the future prediction results using natural earth for country borders and rivers
library(rnaturalearth)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(cowplot)

# Define the color palette for the map
magma_pal <- magma(256, begin = 0, direction = -1)

# Create the main map for 2050 prediction (Agricultural expansion)
main_map <- ggplot() +
  geom_raster(
    data = final_result_future,
    aes(x = x, y = y, fill = mean_predict_future_agri2050),
    na.rm = FALSE
  ) +
  geom_sf(data = countries, color = "grey50", alpha = 0.2, size = 0.1, fill = NA) +
  geom_sf(data = rivers, color = "lightblue", alpha = 0.5, size = 0.04) +
  coord_sf(ylim = c(-60, 90), expand = FALSE) +
  scale_fill_gradientn(
    colors = magma_pal,
    na.value = "white",
    limits = c(0, 2),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(10, "cm"),
      barheight = unit(0.3, "cm")
    )
  ) +
  theme_void() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    plot.title = element_blank(),
    plot.text = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

# Density plots for the predictions
x_density <- ggplot(final_result_future, aes(x = x, y = mean_predict_future_agri2050)) +
  geom_pointdensity(
    aes(color = after_stat(density)),
    size = 0.8,
    shape = 16,
    adjust = 2,
    method = "kde2d"
  ) +
  scale_color_gradientn(
    colors = magma_pal,
    guide = guide_colorbar(
      direction = "vertical",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.3, "cm"),
      barheight = unit(1.5, "cm")
    )
  ) +
  labs(y = NULL, x = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

y_density <- ggplot(final_result_future, aes(x = mean_predict_future_agri2050, y = y)) +
  geom_pointdensity(
    aes(color = after_stat(density)),
    size = 0.8,
    shape = 16,
    adjust = 2,
    method = "kde2d"
  ) +
  scale_color_gradientn(
    colors = magma_pal %>% rev(),
    guide = guide_colorbar(
      direction = "vertical",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.3, "cm"),
      barheight = unit(1.5, "cm")
    )
  ) +
  labs(y = NULL, x = NULL) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

# Combine the main map, x and y density plots
library(patchwork)
design <- "
  11112
  11112
  33334
"

combined_plot <- wrap_plots(
  main_map, # Position 1
  y_density, # Position 2
  x_density, # Position 3
  design = design,
  heights = c(3, 3, 0.5), # Row height proportions
  widths = c(4, 0.8, 4) # Column width proportions
)

# Save the individual and combined plots
ggsave(y_density, filename = "/path/to/data/Future2050_1_worldmap_subplot_y.png", height = 5, width = 1, dpi = 600)
ggsave(x_density, filename = "/path/to/data/Future2050_1_worldmap_subplot_x.png", height = 1, width = 6.5, dpi = 600)
ggsave(main_map, filename = "/path/to/data/Future2050_1_worldmap.png", height = 7, width = 9, dpi = 600)

# Optional: Save the legend plot (if needed)
legend_plot <- ggplot() +
  geom_raster(
    data = final_result_future,
    aes(x = x, y = y, fill = mean_predict_future_agri2050),
    na.rm = FALSE
  ) +
  geom_sf(data = countries, color = "grey50", alpha = 0.2, size = 0.1, fill = NA) +
  geom_sf(data = rivers, color = "lightblue", alpha = 0.5, size = 0.04) +
  coord_sf(ylim = c(-60, 90), expand = FALSE) +
  scale_fill_gradientn(
    colors = magma_pal,
    na.value = "white",
    limits = c(0, 2),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(10, "cm"),
      barheight = unit(0.3, "cm")
    )
  ) +
  theme_void() +
  theme(
    plot.title = element_blank(),
    plot.text = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0)
  ) +
  labs(fill = "Enterobacteriaceae Abundance (%)")

ggsave(legend_plot, filename = "/path/to/data/Future2050_1_worldmap_legend.pdf", family = "ArialMT")

# Plot distribution density for both 2010 and 2050 data
ggplot(data = final_result_future, aes(x = mean_predict_future_agri2050)) +
  geom_histogram(
    aes(y = after_stat(density)), # Change y-axis to density
    binwidth = 0.03, # Adjust binwidth
    fill = "#fbf7c0", # Fill color
    color = "grey", # Border color
    alpha = 0.7 # Transparency
  ) +
  geom_density(
    color = "#b14474", # Red density curve
    linewidth = 0.6, # Line width
    alpha = 0.8, # Transparency
    adjust = 4
  ) +
  geom_hline(yintercept = 0, color = "#b14474", linewidth = 0.5) +
  coord_cartesian(xlim = c(0, 3)) +
  labs(
    title = NULL,
    x = "Abundance (%)",
    y = "Distribution Density (%)"
  ) +
  theme_base() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# Save the distribution density plot
ggsave("/path/to/data/Future2050_1_distribution_density.pdf", family = "ArialMT", width = 3.5, height = 3)

# Calculating summary statistics (mean, SD, standard error) for 2010 and 2050 predictions
library(dplyr)

# Calculate statistics for 2010 data
stats_2010 <- final_result %>%
  group_by(x) %>%
  summarise(
    mean_2010 = mean(mean_predict_2010, na.rm = TRUE),
    sd_2010 = sd(mean_predict_2010, na.rm = TRUE),
    n_2010 = n(),
    se_2010 = sd_2010 / sqrt(n_2010)
  ) %>%
  ungroup()

# Calculate statistics for 2050 data
stats_2050 <- final_result_future %>%
  group_by(x) %>%
  summarise(
    mean_2050 = mean(mean_predict_future_agri2050, na.rm = TRUE),
    sd_2050 = sd(mean_predict_future_agri2050, na.rm = TRUE),
    n_2050 = n(),
    se_2050 = sd_2050 / sqrt(n_2050)
  ) %>%
  ungroup()

# Combine both 2010 and 2050 statistics
combined_stats <- inner_join(stats_2010, stats_2050, by = "x") %>%
  mutate(
    lower_2010 = mean_2010 - 1.96 * se_2010, # 95% confidence interval for 2010
    upper_2010 = mean_2010 + 1.96 * se_2010,
    lower_2050 = mean_2050 - 1.96 * se_2050, # 95% confidence interval for 2050
    upper_2050 = mean_2050 + 1.96 * se_2050
  )

# Plot the comparison between 2010 and 2050 using line plots with confidence intervals
ggplot(combined_stats) +
  # 2010 data and confidence interval
  geom_ribbon(
    aes(x = x, ymin = lower_2010, ymax = upper_2010),
    fill = "#398aa4", alpha = 0.3
  ) +
  geom_line(
    aes(x = x, y = mean_2010),
    color = "#398aa4", size = 0.5, alpha = 1
  ) +
  # 2050 data and confidence interval
  geom_ribbon(
    aes(x = x, ymin = lower_2050, ymax = upper_2050),
    fill = "#da2c48", alpha = 0.3
  ) +
  geom_line(
    aes(x = x, y = mean_2050),
    color = "#da2c48", size = 0.5, alpha = 1
  ) +
  labs(
    x = NULL,
    y = "Abundance (%)",
    title = NULL,
    subtitle = NULL
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", size = 0.2)
  ) +
  guides(color = guide_legend(title = "Year"))

# Save the line plot comparing 2010 and 2050 predictions
ggsave("/path/to/data/Future2050_compare2010_lineplot.pdf", family = "ArialMT", width = 145, height = 20, units = "mm")
