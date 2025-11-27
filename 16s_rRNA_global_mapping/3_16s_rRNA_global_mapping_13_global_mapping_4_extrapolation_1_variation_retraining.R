# Load necessary packages
library(h2o) # H2O library for machine learning models
library(dplyr) # Data manipulation package

# Initialize H2O with multiple threads and memory allocation
h2o.init(nthreads = 5, max_mem_size = "80G", port = 54321)

# Load pre-trained models for agricultural and non-agricultural data
aml_nona_leader <- readRDS("/path/to/project/nonagri_final_model.rds")
aml_agri_leader <- readRDS("/path/to/project/agri_final_model.rds")

# Prepare training dataset for agricultural data
split_agri <- as.h2o(cbind(
  read.csv("/path/to/project/agri_models_x_2_feature_selected_iqr.csv"),
  data.frame(Enterobacteriaceae = read.csv("/path/to/project/agri_models_y_iqr.csv")$Enterobacteriaceae)
)) %>%
  h2o.splitFrame(ratios = 0.9, seed = 33)
train_agri <- split_agri[[1]] # 90% for training
test_agri <- split_agri[[2]] # 10% for testing

# Prepare training dataset for non-agricultural data
split_nona <- as.h2o(cbind(
  read.csv("/path/to/project/nonagri_models_x_2_feature_selected_iqr.csv"),
  data.frame(Enterobacteriaceae = read.csv("/path/to/project/nonagri_models_y_iqr.csv")$Enterobacteriaceae)
)) %>%
  h2o.splitFrame(ratios = 0.9, seed = 9)
train_nona <- split_nona[[1]] # 90% for training
test_nona <- split_nona[[2]] # 10% for testing

# Function to extract GBM (Gradient Boosting Machine) model parameters
extract_gbm_params <- function(model) {
  keep_params <- c(
    "ntrees", "max_depth", "learn_rate", "sample_rate", "col_sample_rate", "col_sample_rate_per_tree",
    "distribution", "nfolds", "fold_assignment", "score_tree_interval", "min_rows", "stopping_metric",
    "stopping_tolerance", "min_split_improvement", "keep_cross_validation_predictions",
    "keep_cross_validation_models", "auto_rebalance", "gainslift_bins", "calibration_method",
    "min_split_improvement", "col_sample_rate_per_tree", "col_sample_rate_change_per_level",
    "col_sample_rate", "huber_alpha", "tweedie_power", "quantile_alpha", "distribution",
    "learn_rate_annealing", "build_tree_one_node", "pred_noise_bandwidth"
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
  original_model <- model # Assume aml is a pre-trained AutoML object
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

# Load prediction data for non-agricultural and agricultural models
preddf_nona_f <- readRDS("/path/to/project/nonagri_preddf_f.rds")
preddf_agri_f <- readRDS("/path/to/project/agri_preddf_f.rds")

# Convert prediction data to H2O frames
h2o_preddf_agri <- as.h2o(preddf_agri_f)
h2o_preddf_nona <- as.h2o(preddf_nona_f)

# Function to retrain and make predictions
retrain_and_predict <- function(seed) {
  print(seed)
  # Train the models
  nona_model <- do.call(h2o.xgboost, c(params_nona, list(seed = seed, training_frame = train_nona)))
  agri_model <- do.call(h2o.xgboost, c(params_agri, list(seed = seed, training_frame = train_agri)))

  # Make predictions
  pred_nona <- as.data.frame(h2o.predict(nona_model, h2o_preddf_nona))$predict
  pred_agri <- as.data.frame(h2o.predict(agri_model, h2o_preddf_agri))$predict

  # Key modification: replace all negative predictions with 0
  pred_nona <- ifelse(pred_nona < 0, 0, pred_nona)
  pred_agri <- ifelse(pred_agri < 0, 0, pred_agri)

  # Clean up memory
  h2o.rm(nona_model)
  gc()

  # Return predictions as a data frame
  data.frame(x = preddf_nona_f$x, y = preddf_nona_f$y, pred_agri = pred_agri, pred_nona = pred_nona, iteration = seed)
}

# Run retraining and prediction for the best seeds
predictions <- purrr::map_df(best_seeds, retrain_and_predict, .id = "run")

# Combine results and compute weighted predictions
combined_predictions <- predictions %>%
  dplyr::select(-run) %>%
  left_join(preddf_agri_f %>% dplyr::select(x, y, AgriculturalLandCover_2010), by = c("x", "y")) %>%
  mutate(pred_2010 = pred_agri * AgriculturalLandCover_2010 + pred_nona * (1 - AgriculturalLandCover_2010))

# Save combined predictions
combined_predictions %>% saveRDS("/path/to/project/Final_1_risk_abun_pred_mean_sd_run_best100.rds")

# Calculate robust statistics for each location
final_result <- combined_predictions %>%
  group_by(x, y) %>%
  summarize(
    mean_predict_2010 = mean(pred_2010, na.rm = TRUE),
    median_predict_2010 = median(pred_2010, na.rm = TRUE),
    sd_predict_2010 = sd(pred_2010, na.rm = TRUE),
    iqr_predict_2010 = IQR(pred_2010, na.rm = TRUE), # Interquartile Range
    mad_predict_2010 = mad(pred_2010, na.rm = TRUE), # Median Absolute Deviation

    # Add requested variation metrics
    cv_iqr = iqr_predict_2010 / median_predict_2010, # Coefficient of Variation for IQR
    cv_mad = mad_predict_2010 / median_predict_2010, # Coefficient of Variation for MAD

    # Retain original coefficient of variation for reference
    cv_sd = sd_predict_2010 / mean_predict_2010,
    .groups = "drop"
  )

# Save final results with robust statistics
final_result %>% saveRDS("/path/to/project/Final_1_risk_abun_pred_mean_sd_run_best100_stats.rds")

# Display the first few rows of the final result
head(final_result)

# Calculate MSE (Mean Squared Error) based best seeds
best_seeds <- head(predictions_all %>% arrange(mse), 100)$iteration
predictions <- purrr::map_df(best_seeds, retrain_and_predict, .id = "run")

# Combine results and calculate weighted predictions
combined_predictions <- predictions %>%
  dplyr::select(-run) %>%
  left_join(preddf_agri_f %>% dplyr::select(x, y, AgriculturalLandCover_2010), by = c("x", "y")) %>%
  mutate(pred_2010 = pred_agri * AgriculturalLandCover_2010 + pred_nona * (1 - AgriculturalLandCover_2010))

# Save predictions based on MSE
combined_predictions %>% saveRDS("/path/to/project/Final_1_risk_abun_pred_mean_sd_run_best100_mse.rds")

# Calculate robust statistics for each location based on MSE
final_result <- combined_predictions %>%
  group_by(x, y) %>%
  summarize(
    mean_predict_2010 = mean(pred_2010, na.rm = TRUE),
    median_predict_2010 = median(pred_2010, na.rm = TRUE),
    sd_predict_2010 = sd(pred_2010, na.rm = TRUE),
    iqr_predict_2010 = IQR(pred_2010, na.rm = TRUE), # Interquartile Range
    mad_predict_2010 = mad(pred_2010, na.rm = TRUE), # Median Absolute Deviation

    # Add requested variation metrics
    cv_iqr = iqr_predict_2010 / median_predict_2010, # Coefficient of Variation for IQR
    cv_mad = mad_predict_2010 / median_predict_2010, # Coefficient of Variation for MAD

    # Retain original coefficient of variation for reference
    cv_sd = sd_predict_2010 / mean_predict_2010,
    .groups = "drop"
  )

# Save final results with robust statistics based on MSE
final_result %>% saveRDS("/path/to/project/Final_1_risk_abun_pred_mean_sd_run_best100_stats_mse.rds")
