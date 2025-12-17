library(h2o)

# Initialize H2O cluster with 5 threads and 80 GB of memory
h2o.init(nthreads = 5, max_mem_size = "80G")

aml <- readRDS("/path/to/data/automl_res_seed2_basicmodels_2000.rds")

# Get the leaderboard of AutoML models
lb <- h2o.get_leaderboard(object = aml, extra_columns = "ALL") %>% as.data.frame()
h2o.varimp_plot(aml@leader) # Plot variable importance of the best model

# Read feature matrix X and target variable y for non-agricultural data
X <- read.csv("/path/to/data/nonagri_models_x_2_feature_selected.csv")
y <- read.csv("/path/to/data/nonagri_models_y_1_original.csv")$Abun_trans

# Combine the feature matrix X and the target variable y into a complete dataset
seed <- 2
h2o_data <- as.h2o(cbind(X, data.frame(Abun_trans = y)))

# Split the data into training (90%) and testing (10%) sets
split <- h2o.splitFrame(
  h2o_data,
  ratios = 0.9,
  seed = seed
)
train <- split[[1]] # 90% training data
test <- split[[2]] # 10% testing data

# Create a larger dataset for performance testing (prediction time)
tt <- as.h2o(cbind(X, data.frame(Abun_trans = y)) %>% rbind(cbind(X, data.frame(Abun_trans = y)))
  %>% rbind(cbind(X, data.frame(Abun_trans = y)))
  %>% rbind(cbind(X, data.frame(Abun_trans = y)))
  %>% rbind(cbind(X, data.frame(Abun_trans = y)))
  %>% rbind(cbind(X, data.frame(Abun_trans = y)))
  %>% rbind(cbind(X, data.frame(Abun_trans = y)))
  %>% rbind(cbind(X, data.frame(Abun_trans = y)))
  %>% rbind(cbind(X, data.frame(Abun_trans = y)))
  %>% rbind(cbind(X, data.frame(Abun_trans = y))))

# GBM model evaluation
lb_gbm <- h2o.get_best_model(aml, algorithm = "GBM", criterion = "RMSE")
start_time <- Sys.time()
predict(lb_gbm, tt)
end_time <- Sys.time()
perf_lb_gbm <- lb %>%
  filter(model_id == lb_gbm@model_id) %>%
  mutate(
    r2 = h2o.r2(h2o.performance(lb_gbm, newdata = train)),
    r2_cv_sd = lb_gbm@model$cross_validation_metrics_summary$sd[6],
    predict_per_row = as.numeric((end_time - start_time) / nrow(tt))
  )

# DRF model evaluation
lb_drf <- h2o.get_best_model(aml, algorithm = "DRF", criterion = "RMSE")
start_time <- Sys.time()
predict(lb_drf, tt)
end_time <- Sys.time()
perf_lb_drf <- lb %>%
  filter(model_id == lb_drf@model_id) %>%
  mutate(
    r2 = h2o.r2(h2o.performance(lb_drf, newdata = train)),
    r2_cv_sd = lb_drf@model$cross_validation_metrics_summary$sd[6],
    predict_per_row = as.numeric((end_time - start_time) / nrow(tt))
  )

# GLM model evaluation
lb_glm <- h2o.get_best_model(aml, algorithm = "GLM", criterion = "RMSE")
start_time <- Sys.time()
predict(lb_glm, tt)
end_time <- Sys.time()
perf_lb_glm <- lb %>%
  filter(model_id == lb_glm@model_id) %>%
  mutate(
    r2 = h2o.r2(h2o.performance(lb_glm, newdata = train)),
    r2_cv_sd = lb_glm@model$cross_validation_metrics_summary$sd[7],
    predict_per_row = as.numeric((end_time - start_time) / nrow(tt))
  )

# XGBoost model evaluation
lb_xgb <- h2o.getModel("XGBoost_grid_1_AutoML_191_20250605_154609_model_1674")
start_time <- Sys.time()
predict(lb_xgb, tt)
end_time <- Sys.time()
perf_lb_xgb <- lb %>%
  filter(model_id == lb_xgb@model_id) %>%
  mutate(
    r2 = h2o.r2(h2o.performance(lb_xgb, newdata = train)),
    r2_cv_sd = lb_xgb@model$cross_validation_metrics_summary$sd[6],
    predict_per_row = as.numeric((end_time - start_time) / nrow(tt))
  )

# StackedEnsemble model evaluation
lb_ensem <- h2o.get_best_model(aml, algorithm = "StackedEnsemble", criterion = "RMSE")
start_time <- Sys.time()
predict(lb_ensem, tt)
end_time <- Sys.time()
perf_lb_ensem <- lb %>%
  filter(model_id == lb_ensem@model_id) %>%
  mutate(
    r2 = h2o.r2(h2o.performance(lb_ensem, newdata = train)),
    r2_cv_sd = lb_ensem@model$cross_validation_metrics_summary$sd[7],
    predict_per_row = as.numeric((end_time - start_time) / nrow(tt))
  )

# Combine all model performance metrics
perf_all <- rbind(perf_lb_gbm, perf_lb_drf, perf_lb_glm, perf_lb_xgb) %>%
  mutate(predict_per_row_ms = 1000 * predict_per_row) %>%
  dplyr::select(-predict_time_per_row_ms, -predict_per_row)

# Save performance data to CSV
perf_all %>% write.csv("/path/to/data/nonagri_1_model_radar_performance.csv")

# Inverse certain metrics (except r2)
reverse_cols <- c("rmse", "mse", "mae", "rmsle", "mean_residual_deviance", "training_time_ms", "predict_per_row_ms", "r2_cv_sd")
perf_all[reverse_cols] <- -perf_all[reverse_cols]

# Apply Min-Max normalization to the data
min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

scaled_data <- as.data.frame(lapply(perf_all[, c(reverse_cols, "r2")], min_max_scale))

# Prepare data for radar chart
radar_data <- cbind(
  model = perf_all$model_id,
  scaled_data
)[4:1, ]

# Add extreme value rows (all 1s and all 0s)
radar_data <- rbind(
  rep(1, ncol(scaled_data)) %>% setNames(names(scaled_data)),
  rep(0, ncol(scaled_data)) %>% setNames(names(scaled_data)),
  radar_data
)

# Plot radar chart
colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A") %>% rev()

library(fmsb)
library(extrafont)

pdf("/path/to/data/nonagri_1_model_radar.pdf", family = "ArialMT", width = 5.5, height = 5)
radarchart(
  radar_data %>% dplyr::select(-model, -mean_residual_deviance, -mse, -r2) %>%
    rename(
      RMSE = rmse, Robustness = r2_cv_sd, MAE = mae, RMSLE = rmsle,
      `Training time` = training_time_ms, `Computational efficiency` = predict_per_row_ms
    ),
  axistype = 1,
  pcol = colors,
  pfcol = alpha(colors, 0.15),
  plwd = 2,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  axislabcol = "grey50",
  vlcex = 0.8,
  title = "",
  calcex = 0.6,
  caxislabels = c("", "0.25", "0.5", "0.75", "1.0"),
  seg = 4
)

# Add legend
legend(
  "topright",
  legend = c("GBM", "RFR", "GLM", "XGB"),
  bty = "n",
  pch = 20,
  col = colors %>% rev(),
  text.col = "black",
  cex = 0.8
)
dev.off()

# Save radar chart data
radar_data %>% write.csv("/path/to/data/nonagri_1_model_radar.csv")

# Save leaderboard excluding specific models
lb %>%
  filter(!model_id %in% lb$model_id[1:6]) %>%
  write.csv("/path/to/data/nonagri_1_model_leaderboard.csv")

# Step 1: Extract original model parameters
original_model <- lb_gbm # Assuming lb_gbm is the trained AutoML model
original_params <- original_model@params$actual

# List of parameters to keep from the model
keep_params <- c(
  "ntrees", "max_depth", "learn_rate", "sample_rate",
  "col_sample_rate", "col_sample_rate_per_tree",
  "distribution", "nfolds", "fold_assignment",
  "score_tree_interval", "min_rows", "stopping_metric",
  "stopping_tolerance", "min_split_improvement",
  "keep_cross_validation_predictions", "keep_cross_validation_models"
)

# Create a list for the selected parameters
model_params <- list()
for (param in keep_params) {
  if (!is.null(original_params[[param]])) {
    model_params[[param]] <- original_params[[param]]
  }
}

# Step 2: Loop to train models with different seeds (100 seeds)
feature_importance_list <- list()

# Loop to train 100 models with different seeds and extract feature importance
for (seed in 1:100) { # Looping from seed 1 to 100
  # Update the parameters with the current seed
  current_params <- c(model_params, list(seed = seed))

  # Train the model with the updated parameters
  model <- do.call(h2o.gbm, c(
    list(
      x = names(X), # Feature column names
      y = "Abun_trans", # Target column name
      training_frame = train # Training data
    ),
    current_params
  ))

  # Extract feature importance for each model
  fi <- h2o.varimp(model) %>%
    mutate(seed = as.integer(seed)) %>% # Add the seed column
    dplyr::select(variable, scaled_importance, seed)

  # Store the feature importance data
  feature_importance_list[[seed]] <- fi

  # Clean up memory after each model
  h2o.rm(model)
}

# Step 3: Combine the feature importance data and sort by median importance
df <- bind_rows(feature_importance_list) %>% mutate(variable = factor(variable))

# Calculate the median feature importance for each variable and sort
median_importance <- df %>%
  group_by(variable) %>%
  summarise(median_imp = median(scaled_importance)) %>%
  arrange(desc(median_imp))

# Reorder the factor levels based on the median importance
df$variable <- factor(df$variable, levels = median_importance$variable)

# Step 4: Plot the top 10 features based on median importance using a boxplot
top10_vars <- median_importance %>%
  slice_head(n = 10) %>% # Select the top 10 features based on median importance
  pull(variable)

df_top10 <- df %>%
  filter(variable %in% top10_vars) %>% # Filter to keep only the top 10 features
  mutate(variable = factor(variable, levels = rev(top10_vars))) # Reverse order to have the most important at the top

# Load the ggplot2 library for creating plots
library(ggthemes)

# Create the boxplot for the top 10 most important features
ggplot(df_top10, aes(x = scaled_importance, y = variable)) + # Swap x and y mappings
  geom_boxplot(fill = "lightblue", alpha = 0.8, width = 0.6) +
  theme_bw() +
  labs(
    title = "Top 10 Feature Importance Stability (100 Seeds)",
    x = "Scaled Importance",
    y = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0)) + # Start importance axis from 0
  geom_vline(xintercept = 0, linetype = "dashed") # Add a dashed line at 0 for reference
ggsave(filename = "/path/to/figure/nonagri_2_feature_impr_100repeat.pdf", height = 5, width = 6)

# Save the feature importance data to CSV for further analysis
df %>%
  data.frame() %>%
  write.csv("/path/to/data/nonagri_2_feature_impr_100repeat.csv", row.names = FALSE)
df_top10 %>%
  as.data.frame() %>%
  write.csv("/path/to/data/nonagri_2_feature_impr_100repeat_top10feature.csv", row.names = FALSE)

# Step 5: Calculate relative feature importance and plot using color gradients
df_m <- df %>% mutate(relative_importance = 100 * (df %>% group_by(seed) %>% summarise(relative_importance = scaled_importance / sum(scaled_importance)))[, 2] %>% pull())

# Get the top 10 features by relative importance
top10_vars <- df_m %>%
  group_by(variable) %>%
  summarise(median_imp = median(relative_importance)) %>%
  arrange(desc(median_imp)) %>%
  slice_head(n = 10) %>% # Get the top 10 features
  pull(variable)

# Filter to keep only the top 10 features
df_top10 <- df_m %>%
  filter(variable %in% top10_vars) %>% # Filter for top 10 features
  mutate(variable = factor(variable, levels = rev(top10_vars))) # Reverse order for better visualization

# Create a custom color gradient (from dark red to light blue)
color_palette <- colorRampPalette(c("#d55740", "#fcf1b0", "#4b80b5"), bias = 0.8)

# Calculate median relative importance for each feature
df_median <- df_top10 %>%
  group_by(variable) %>%
  summarise(median_importance = median(relative_importance)) %>%
  ungroup()

# Merge the median importance data back into the original data
df_top10_mean <- df_top10 %>%
  left_join(df_median, by = "variable")

# Plot the top 10 features with a color gradient based on their median relative importance
ggplot(df_top10_mean, aes(x = relative_importance, y = variable)) +
  geom_boxplot(
    aes(fill = median_importance), # Use median importance for color mapping
    width = 0.6,
    alpha = 0.6
  ) +
  scale_fill_gradientn(
    colours = color_palette(100) %>% rev(), # Use a 100-step color gradient
    guide = "none" # Hide the legend
  ) +
  theme_bw() +
  labs(
    x = "Relative variable importance (%)",
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave(filename = "/path/to/data/nonagri_2_feature_impr_100repeat.pdf", height = 4, width = 6)

# Save the final feature importance data
df_top10_mean %>% write.csv("/path/to/data/nonagri_2_feature_impr_100repeat_top10feature.csv", row.names = FALSE)
