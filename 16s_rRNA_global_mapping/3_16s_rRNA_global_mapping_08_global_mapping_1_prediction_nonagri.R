library(h2o)
# Initialize the H2O instance with 5 threads and 80GB of memory
h2o.init(nthreads = 5, max_mem_size = "80G")

# Load the AutoML model from a saved RDS file (use a dummy path)
aml <- readRDS("/path/to/data/AutoML_nonagri/optimal_iqr/automl_res_seed37_MAE_iqr.rds")

# Get the leaderboard of models from the AutoML run
lb <- h2o.get_leaderboard(object = aml, extra_columns = "ALL") %>% as.data.frame()

# Plot variable importance for a specific model in the leaderboard
h2o.varimp_plot(h2o.getModel(lb$model_id[49]))

# Save the best model from the leaderboard
h2o.getModel(lb$model_id[49]) %>% saveRDS("/path/to/data/final/nonagri_final_model.rds")

# Read the feature matrix X and the target variable y
X <- read.csv("/path/to/data/nonagri_models_x_2_feature_selected_iqr.csv")
y <- read.csv("/path/to/data/nonagri_models_y_iqr.csv")$Enterobacteriaceae

# Combine X and y to create the full dataset
seed <- 37
h2o_data <- as.h2o(cbind(X, data.frame(Enterobacteriaceae = y)))

# Split the dataset into training and testing sets
split <- h2o.splitFrame(h2o_data, ratios = 0.9, seed = seed)
train <- split[[1]] # 80% training set
test <- split[[2]] # 20% test set

# Create a larger dataset for predicting performance time experiments
tt <- as.h2o(cbind(X, data.frame(Enterobacteriaceae = y)) %>%
  rbind(cbind(X, data.frame(Enterobacteriaceae = y))) %>%
  rbind(cbind(X, data.frame(Enterobacteriaceae = y))) %>%
  rbind(cbind(X, data.frame(Enterobacteriaceae = y))) %>%
  rbind(cbind(X, data.frame(Enterobacteriaceae = y))) %>%
  rbind(cbind(X, data.frame(Enterobacteriaceae = y))) %>%
  rbind(cbind(X, data.frame(Enterobacteriaceae = y))) %>%
  rbind(cbind(X, data.frame(Enterobacteriaceae = y))) %>%
  rbind(cbind(X, data.frame(Enterobacteriaceae = y))) %>%
  rbind(cbind(X, data.frame(Enterobacteriaceae = y))))

# Set performance metric for the models
metric <- "MAE"

## GBM model ##
lb_gbm <- h2o.getModel(lb$model_id[104])
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

## DRF model ##
lb_drf <- h2o.get_best_model(aml, algorithm = "DRF", criterion = metric)
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

## GLM model ##
lb_glm <- h2o.get_best_model(aml, algorithm = "GLM", criterion = metric)
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

## XGBoost model ##
lb_xgb <- h2o.getModel(lb$model_id[49])
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

# Combine all performance metrics
perf_all <- rbind(perf_lb_gbm, perf_lb_drf, perf_lb_glm, perf_lb_xgb) %>%
  mutate(predict_per_row_ms = 1000 * predict_per_row) %>%
  dplyr::select(-predict_time_per_row_ms, -predict_per_row)

# Save the performance metrics to a CSV file
perf_all %>% write.csv("/path/to/data/final/nonagri_1_model_radar_performance.csv")

# Reverse some performance metrics (all except R2)
reverse_cols <- c("rmse", "mse", "mae", "rmsle", "mean_residual_deviance", "training_time_ms", "predict_per_row_ms", "r2_cv_sd")
perf_all[reverse_cols] <- -perf_all[reverse_cols]

# Min-Max scaling for normalization
min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

scaled_data <- as.data.frame(lapply(perf_all[, c(reverse_cols, "r2")], min_max_scale))

# Prepare radar chart data
radar_data <- cbind(
  model = perf_all$model_id,
  scaled_data
)[4:1, ]

# Add extreme rows (all 1s and 0s)
radar_data <- rbind(
  rep(1, ncol(scaled_data)) %>% setNames(names(scaled_data)),
  rep(0, ncol(scaled_data)) %>% setNames(names(scaled_data)),
  radar_data
)

# Create radar chart
colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A") %>% rev()

library(fmsb)
library(extrafont)
pdf("/path/to/data/final/nonagri_1_model_radar.pdf", family = "ArialMT", width = 5.5, height = 5)
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

# Save the radar chart data
radar_data %>% write.csv("/path/to/data/final/nonagri_1_model_radar.csv")
lb[c(49:nrow(lb)), ] %>% write.csv("/path/to/data/final/nonagri_1_model_leaderload.csv")

# Extract model parameters for re-training with different seeds
original_model <- lb_xgb # Example model to extract parameters
original_params <- original_model@params$actual
original_params["seed"] <- NULL
original_params["response_column"] <- NULL
original_params["fold_assignment"] <- NULL
model_params <- original_params[4:length(original_params)]

# Loop to train models with different seeds (100 times)
feature_importance_list <- list()

for (seed in 1:100) { # Loop through seeds from 1 to 100
  current_params <- c(model_params, list(seed = seed))

  # Train the model
  model <- do.call(h2o.xgboost, c(
    list(
      x = X %>% names(), # Feature names
      y = "Enterobacteriaceae", # Target variable
      training_frame = train # Training data
    ),
    current_params
  ))

  # Extract feature importance for each seed
  fi <- h2o.varimp(model) %>%
    mutate(seed = as.integer(seed)) %>% # Add seed column
    dplyr::select(variable, scaled_importance, seed)

  feature_importance_list[[seed]] <- fi

  # Clear the model from memory
  h2o.rm(model)
}

# Step 3: Combine feature importance data and sort by median importance
df <- bind_rows(feature_importance_list) %>% mutate(variable = factor(variable))

# Sort by median importance
median_importance <- df %>%
  group_by(variable) %>%
  summarise(median_imp = median(scaled_importance)) %>%
  arrange(desc(median_imp))

# Set factor levels based on median importance
df$variable <- factor(df$variable, levels = median_importance$variable)

# Step 4: Create a boxplot for the top 10 most important features (flip axes)
top10_vars <- median_importance %>%
  slice_head(n = 10) %>% # Select top 10 features based on median importance
  pull(variable)

df_top10 <- df %>%
  filter(variable %in% top10_vars) %>% # Filter for top 10 features
  mutate(variable = factor(variable, levels = rev(top10_vars))) # Reverse order for boxplot

# Save the feature importance data
df %>%
  data.frame() %>%
  write.csv("/path/to/data/final/nonagri_2_feature_impr_100repeat.csv", row.names = F)
df_top10 %>%
  as.data.frame() %>%
  write.csv("/path/to/data/final/nonagri_2_feature_impr_100repeat_top10feature.csv", row.names = F)

# Step 5: Calculate relative importance for each feature based on the sum of importance across seeds
df_m <- df %>% mutate(relative_importance = 100 * (df %>% group_by(seed) %>%
  summarise(relative_importance = scaled_importance / sum(scaled_importance)))[, 2] %>% pull())

# Select the top 10 most important features based on relative importance
top10_vars <- df_m %>%
  data.frame() %>%
  group_by(variable) %>%
  summarise(median_imp = median(relative_importance)) %>%
  arrange(desc(median_imp)) %>%
  slice_head(n = 10) %>% # Select top 10 features
  pull(variable)

df_top10 <- df_m %>%
  data.frame() %>%
  filter(variable %in% top10_vars) %>% # Filter for top 10 features
  mutate(variable = factor(variable, levels = rev(top10_vars))) # Reverse order

# Step 6: Create a custom color palette for the boxplot
color_palette <- colorRampPalette(c("#d55740", "#fcf1b0", "#4b80b5"), bias = 0.8)

# Calculate median importance for each feature to use as fill color in the boxplot
df_median <- df_top10 %>%
  group_by(variable) %>%
  summarise(median_importance = median(relative_importance)) %>%
  ungroup()

# Merge the median importance back into the original data
df_top10_mean <- df_top10 %>%
  left_join(df_median, by = "variable")

# Step 7: Create and save the boxplot
ggplot(df_top10_mean, aes(x = relative_importance, y = variable)) +
  geom_boxplot(
    aes(fill = median_importance), # Map median importance to color
    width = 0.6,
    alpha = 0.6
  ) +
  scale_fill_gradientn(
    colours = color_palette(100) %>% rev(), # Use a 100-level color gradient
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
ggsave(filename = "/path/to/data/final/nonagri_2_feature_impr_100repeat.pdf", height = 4, width = 6)

# Create another version of the plot without the y-axis text
ggplot(df_top10_mean, aes(x = relative_importance, y = variable)) +
  geom_boxplot(
    aes(fill = median_importance), # Map median importance to color
    width = 0.6,
    alpha = 0.6
  ) +
  scale_fill_gradientn(
    colours = color_palette(100) %>% rev(), # Use a 100-level color gradient
    guide = "none" # Hide the legend
  ) +
  theme_bw() +
  labs(
    x = "Relative variable importance (%)",
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave(filename = "/path/to/data/final/nonagri_2_feature_impr_100repeat_wo_text.pdf", height = 4, width = 3.5)

# Save the data for relative importance and top 10 features
df_top10_mean %>% write.csv("/path/to/data/final/nonagri_2_feature_impr_100repeat_top10.csv", row.names = F)
