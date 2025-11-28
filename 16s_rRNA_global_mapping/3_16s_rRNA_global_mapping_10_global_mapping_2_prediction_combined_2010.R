library(h2o)
h2o.init(nthreads = 5, max_mem_size = "80G")

# Load final models for non-agricultural and agricultural data
aml_nona_leader <- readRDS("/path/to/data/nonagri_final_model.rds")
aml_agri_leader <- readRDS("/path/to/data/agri_final_model.rds")

# Load the filled raster stack and other relevant raster layers
filled_stack <- readRDS("/path/to/data/tif_stack_data_leaf_filled.rds")
rs_landcover <- raster("/path/to/data/Global_MCD12Q1_IGBP_500m_2023_resampled.tif")
rs_agri2010 <- raster("/path/to/data/AgriculturalLandCover_2010_resampled.tif")
rs_Koppen <- raster("/path/to/data/Bio_199_Koppen.tif")

# Combine the filled raster stack with additional data layers
tif_stack_crop_filled <- stack(filled_stack, rs_landcover, rs_agri2010)

# Subset the raster stack for agricultural and non-agricultural features
sub_stack_agri_filled <- tif_stack_crop_filled[[c(
  names(read.csv("/path/to/data/agri_models_x_2_feature_selected_iqr.csv")),
  "AgriculturalLandCover_2010", "Bio_199_Koppen", "Global_MCD12Q1_IGBP_500m_2023_resampled"
)]]
sub_stack_nona_filled <- tif_stack_crop_filled[[c(
  names(read.csv("/path/to/data/nonagri_models_x_2_feature_selected_iqr.csv")),
  "AgriculturalLandCover_2010", "Bio_199_Koppen", "Global_MCD12Q1_IGBP_500m_2023_resampled"
)]]

# Prepare the data frames for prediction
preddf_agri <- as.data.frame(sub_stack_agri_filled, xy = TRUE, na.rm = TRUE)
preddf_nona <- as.data.frame(sub_stack_nona_filled, xy = TRUE, na.rm = TRUE)

# Filter by climate zone & land cover to find overlapping coordinates
xy_both <- inner_join(data.frame(x = preddf_agri$x, y = preddf_agri$y),
  data.frame(x = preddf_nona$x, y = preddf_nona$y),
  by = c("x", "y")
) %>%
  distinct(x, y) %>%
  left_join(preddf_agri %>% dplyr::select(x, y, Bio_199_Koppen, landcover = Global_MCD12Q1_IGBP_500m_2023_resampled))
nrow(xy_both) / nrow(preddf_agri) # 96.23%

# Handle extraction of Koppen climate zone
coor_16S <- read.csv("/path/to/data/p8_16s_ML_source_data.csv") %>% dplyr::select(lon = longitude, lat = latitude)
coordinates(coor_16S) <- c("lon", "lat")

# Extract Koppen values for the coordinates and filter out missing values
cbind(coor_16S %>% as.data.frame(), Koppen = raster::extract(rs_Koppen, coor_16S)) %>%
  filter(is.na(Koppen)) %>%
  View()
Koppen_existed <- c(raster::extract(rs_Koppen, coor_16S) %>% unique(), 1, 2) # NA values are 1, 2
Koppen_existed

# First filtration: remove land cover values from certain categories
xy_both_filtered1 <- xy_both %>% filter(!landcover %in% c(15, 16, 17))
nrow(xy_both_filtered1) / nrow(xy_both) # 87.49%

# Second filtration: filter by Koppen Climate Zone
xy_both_filtered2 <- xy_both_filtered1 %>% filter(Bio_199_Koppen %in% Koppen_existed)
nrow(xy_both_filtered2) / nrow(xy_both) # 80.85%

# Join the filtered data with original prediction data
preddf_agri_f <- xy_both_filtered2 %>% left_join(preddf_agri, by = c("x", "y"))
preddf_nona_f <- xy_both_filtered2 %>% left_join(preddf_nona, by = c("x", "y"))

# Convert data to H2O frame
h2o_preddf_agri <- as.h2o(preddf_agri_f)
h2o_preddf_nona <- as.h2o(preddf_nona_f)

# Function to extract GBM model parameters
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

# Function to extract XGBoost model parameters
extract_xgb_params <- function(model) {
  original_model <- model # Assuming aml is the trained AutoML model
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

# Function for retraining and prediction
retrain_and_predict <- function(seed) {
  print(seed)
  # Train models
  nona_model <- do.call(h2o.xgboost, c(params_nona, list(seed = seed, training_frame = train_nona)))
  agri_model <- do.call(h2o.xgboost, c(params_agri, list(seed = seed, training_frame = train_agri)))

  # Make predictions
  pred_nona <- as.data.frame(h2o.predict(nona_model, h2o_preddf_nona))$predict
  pred_agri <- as.data.frame(h2o.predict(agri_model, h2o_preddf_agri))$predict

  # Modify predictions: replace negative values with 0
  pred_nona <- ifelse(pred_nona < 0, 0, pred_nona)
  pred_agri <- ifelse(pred_agri < 0, 0, pred_agri)

  # Clean up memory
  h2o.rm(nona_model)
  gc()

  # Return the prediction results
  data.frame(x = preddf_nona_f$x, y = preddf_nona_f$y, pred_agri = pred_agri, pred_nona = pred_nona, iteration = seed)
}

# Execute retraining and prediction for 100 iterations
n_iterations <- 100
predictions <- purrr::map_df(1:n_iterations, retrain_and_predict, .id = "run")

# Combine predictions and compute weighted predictions
combined_predictions <- predictions %>%
  dplyr::select(-run) %>%
  left_join(preddf_agri_f %>% dplyr::select(x, y, AgriculturalLandCover_2010), by = c("x", "y")) %>%
  mutate(pred_2010 = pred_agri * AgriculturalLandCover_2010 + pred_nona * (1 - AgriculturalLandCover_2010))

# Calculate mean and standard deviation for each location
final_result <- combined_predictions %>%
  group_by(x, y) %>%
  summarize(
    mean_predict_2010 = mean(pred_2010, na.rm = TRUE),
    sd_predict_2010 = sd(pred_2010, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate variation (coefficient of variation)
final_result <- final_result %>% mutate(variation_2010 = sd_predict_2010 / mean_predict_2010)

# Save the final result (optional)
# final_result %>% write.csv("/path/to/data/Final_1_risk_abun_pred_mean_sd.csv")
final_result <- read.csv("/path/to/data/Final_1_risk_abun_pred_mean_sd.csv", row.names = 1)

# Libraries for plotting
library(rnaturalearth)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(cowplot)

# Get country borders (medium resolution)
countries <- ne_countries(scale = "medium", returnclass = "sf")

# Get major rivers (medium resolution)
rivers <- ne_download(
  scale = "medium", type = "rivers_lake_centerlines",
  category = "physical", returnclass = "sf"
)

# Define color palette (magma)
magma_pal <- magma(256, begin = 0, direction = -1)

# Create the main map
main_map <- ggplot() +
  geom_raster(
    data = final_result,
    aes(x = x, y = y, fill = mean_predict_2010),
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

# Density plot for x-axis (longitude)
x_density <- ggplot(final_result, aes(x = x, y = mean_predict_2010)) +
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
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1)
  )

# Density plot for y-axis (latitude)
y_density <- ggplot(final_result, aes(x = mean_predict_2010, y = y)) +
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
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1)
  )

# Combine the plots
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
  widths = c(4, 0.8, 2.5) # Column width proportions
)

# Save individual plots and combined plot
ggsave(y_density, filename = "/path/to/figure/Final_1_worldmap_subplot_y.png", height = 5, width = 1, dpi = 600)
ggsave(x_density, filename = "/path/to/figure/Final_1_worldmap_subplot_x.png", height = 1, width = 6.5, dpi = 600)
ggsave(main_map, filename = "/path/to/figure/Final_1_worldmap.png", height = 7, width = 9, dpi = 600)

# Add legend (optional)
legend_plot <- ggplot() +
  geom_raster(
    data = final_result,
    aes(x = x, y = y, fill = mean_predict_2010),
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

ggsave(legend_plot, filename = "/path/to/figure/Final_1_worldmap_legend.pdf", family = "ArialMT")

# Histogram and Density Plot
ggplot(data = final_result, aes(x = mean_predict_2010)) +
  geom_histogram(
    binwidth = 0.03, # Adjust binwidth based on data range
    fill = "#fbf7c0", # Fill color
    color = "#fbf7c0", # Border color
    alpha = 0.8 # Transparency
  ) +
  labs(
    title = "Frequency Distribution of mean_predict_2010",
    x = "Predicted Value",
    y = "Frequency"
  ) +
  theme_minimal() + # Clean theme
  theme(plot.title = element_text(hjust = 0.5)) # Center the title

ggsave("/path/to/figure/Final_1_distribution_density.pdf", family = "ArialMT", width = 3.5, height = 3)

# Density plot with overlayed density curve
ggplot(data = final_result, aes(x = mean_predict_2010)) +
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
    y = "Distribution density (%)"
  ) +
  theme_base() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave("/path/to/figure/Final_1_distribution_density.pdf", family = "ArialMT", width = 3.5, height = 3)

# Final result density plot for abundance
ggplot(data = final_result, aes(x = mean_predict_2010)) +
  geom_histogram(
    aes(y = after_stat(density)), # Change y-axis to density (percentage)
    binwidth = 0.03, # Adjust binwidth based on data range
    fill = "#fbf7c0", # Fill color
    color = "grey", # Border color
    alpha = 0.7 # Transparency
  ) +
  # Add a red density curve
  geom_density(
    color = "#b14474", # Red curve
    linewidth = 0.6, # Line width
    alpha = 0.8, # Transparency
    adjust = 4
  ) +
  # Add horizontal line at y = 0
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

# Save the density plot as a PDF
ggsave("/path/to/figure/Final_1_distribution_density.pdf", family = "ArialMT", width = 3.5, height = 3)

# Final Map Plot with Country Boundaries and Rivers, Color Gradient for Predicted Abundance
ggplot() +
  geom_raster(
    data = final_result,
    aes(x = x, y = y, fill = mean_predict_2010),
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

# Save the final world map with predicted abundance
ggsave("/path/to/figure/Final_1_worldmap.png", height = 7, width = 9, dpi = 600)

# Optional: Save additional maps or plots as needed
ggsave("/path/to/figure/Final_1_worldmap_legend.pdf", family = "ArialMT")
