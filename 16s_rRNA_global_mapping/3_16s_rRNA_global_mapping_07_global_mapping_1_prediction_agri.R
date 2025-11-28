library(h2o)

# Initialize the H2O cluster with specified number of threads and memory size
h2o.init(nthreads = 5, max_mem_size = "80G")

# Load the AutoML model from a previously saved RDS file
aml <- readRDS("/path/to/data/AutoML_agri/optimal_iqr/automl_res_seed28_MAE_iqr.rds")

# Plot variable importance for the best model (leader model)
h2o.varimp_plot(aml@leader)

# Get the leaderboard and convert it to a data frame
lb <- h2o.get_leaderboard(object = aml, extra_columns = "ALL") %>% as.data.frame()

# Plot variable importance for the 12th model in the leaderboard
h2o.varimp_plot(h2o.getModel(lb$model_id[12]))

# Save the 12th model in the leaderboard as an RDS file
h2o.getModel(lb$model_id[12]) %>% saveRDS("/path/to/data/final/agri_final_model.rds")

# Read the feature matrix X and target variable y
X <- read.csv("/path/to/data/agri_models_x_2_feature_selected.csv") # Feature matrix
y <- read.csv("/path/to/data/agri_models_y_1_original.csv")$Enterobacteriaceae # Target variable

# Combine the features and target variable into a complete dataset
seed <- 28
h2o_data <- as.h2o(cbind(X, data.frame(Enterobacteriaceae = y)))

# Split the dataset into training (90%) and testing (10%) sets
split <- h2o.splitFrame(
  h2o_data,
  ratios = 0.9,
  seed = seed
)
train <- split[[1]] # 90% training set
test <- split[[2]] # 10% testing set

# Create a larger dataset for prediction time performance experiment
tt <- as.h2o(cbind(X, data.frame(Enterobacteriaceae = y)) %>% rbind(cbind(X, data.frame(Enterobacteriaceae = y)))
  %>% rbind(cbind(X, data.frame(Enterobacteriaceae = y)))
  %>% rbind(cbind(X, data.frame(Enterobacteriaceae = y)))
  %>% rbind(cbind(X, data.frame(Enterobacteriaceae = y)))
  %>% rbind(cbind(X, data.frame(Enterobacteriaceae = y)))
  %>% rbind(cbind(X, data.frame(Enterobacteriaceae = y)))
  %>% rbind(cbind(X, data.frame(Enterobacteriaceae = y)))
  %>% rbind(cbind(X, data.frame(Enterobacteriaceae = y)))
  %>% rbind(cbind(X, data.frame(Enterobacteriaceae = y))))

# Model performance metrics for MAE (Mean Absolute Error)
metric <- "MAE"

## GBM Model ##
# Get the best GBM model based on the specified metric
lb_gbm <- h2o.get_best_model(aml, algorithm = "GBM", criterion = metric)
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

## DRF Model ##
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

## GLM Model ##
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

## XGBoost Model ##
lb_xgb <- h2o.getModel(lb$model_id[12]) # Use the 12th model for XGBoost
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

# Combine all performance data into one data frame
perf_all <- rbind(perf_lb_gbm, perf_lb_drf, perf_lb_glm, perf_lb_xgb) %>%
  mutate(predict_per_row_ms = 1000 * predict_per_row) %>%
  dplyr::select(-predict_time_per_row_ms, -predict_per_row)

# Save the performance metrics to CSV
perf_all %>% write.csv("/path/to/data/agri_1_model_radar_performance.csv")

# Step 1: Reverse the values (except for R-squared)
reverse_cols <- c("rmse", "mse", "mae", "rmsle", "mean_residual_deviance", "training_time_ms", "predict_per_row_ms", "r2_cv_sd")
perf_all[reverse_cols] <- -perf_all[reverse_cols]

# Step 2: Min-Max normalization
min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

scaled_data <- as.data.frame(lapply(perf_all[, c(reverse_cols, "r2")], min_max_scale))

# Step 3: Prepare data for radar chart
radar_data <- cbind(
  model = perf_all$model_id,
  scaled_data
)[4:1, ]

# Add rows for maximum and minimum values (0 and 1)
radar_data <- rbind(
  rep(1, ncol(scaled_data)) %>% setNames(names(scaled_data)),
  rep(0, ncol(scaled_data)) %>% setNames(names(scaled_data)),
  radar_data
)

# Step 4: Create and save radar chart
colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A") %>% rev()
library(fmsb)
library(extrafont)
pdf("/path/to/data/agri_1_model_radar.pdf", family = "ArialMT", width = 5.5, height = 5)
radarchart(
  radar_data %>% dplyr::select(-model, -mean_residual_deviance, -mse, -r2) %>%
    rename(
      RMSE = rmse, Robustness = r2_cv_sd,
      MAE = mae, RMSLE = rmsle,
      `Training time` = training_time_ms,
      `Computational efficiency` = predict_per_row_ms
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

# Save the radar chart data to a CSV file
radar_data %>% write.csv("/path/to/data/agri_1_model_radar.csv")

# Save the leaderboard from the 12th model to a CSV file
lb[12:nrow(lb), ] %>%
  write.csv("/path/to/data/agri_1_model_leaderload.csv")

# Step 1: Extract original model parameters
# Assuming `aml` is a trained AutoML object or `lb_gbm` is a trained model
original_model <- lb_xgb # Use the XGBoost model from the leaderboard as the original model
original_params <- original_model@params$actual

# Remove unnecessary parameters (e.g., "seed", "response_column")
original_params["seed"] <- NULL
original_params["response_column"] <- NULL
original_params["fold_assignment"] <- NULL

# Create a list of model parameters to be used for training
model_params <- original_params[4:length(original_params)]

# Step 2: Loop to train models with different seeds (1 to 100)
feature_importance_list <- list()

for (seed in 1:100) { # Loop through seeds 1 to 100
  # Set the seed and combine parameters for training
  current_params <- c(model_params, list(seed = seed))

  # Train the model using the specified parameters
  model <- do.call(h2o.xgboost, c(
    list(
      x = X %>% names(), # Feature column names
      y = "Enterobacteriaceae", # Target column name
      training_frame = train # Original training dataset
    ),
    current_params
  ))

  # Extract feature importance
  fi <- h2o.varimp(model) %>%
    mutate(seed = as.integer(seed)) %>% # Add seed column
    dplyr::select(variable, scaled_importance, seed)

  # Store feature importance for the current seed
  feature_importance_list[[seed]] <- fi

  # Clean up the model to free memory
  h2o.rm(model)
}

# Step 3: Merge feature importance data and sort by average importance
df <- bind_rows(feature_importance_list) %>% mutate(variable = factor(variable))

# Calculate the median importance for each variable
median_importance <- df %>%
  group_by(variable) %>%
  summarise(median_imp = median(scaled_importance)) %>%
  arrange(desc(median_imp))

# Set the factor levels based on median importance
df$variable <- factor(df$variable, levels = median_importance$variable)

# Step 4: Plot the top 10 features based on median importance (boxplot with horizontal orientation)
top10_vars <- median_importance %>%
  slice_head(n = 10) %>% # Select the top 10 features with highest median importance
  pull(variable)

df_top10 <- df %>%
  filter(variable %in% top10_vars) %>% # Filter for the top 10 features
  mutate(variable = factor(variable, levels = rev(top10_vars))) # Reverse order for the plot

# Use custom color palette for the boxplot
library(ggthemes)
color_palette <- colorRampPalette(c("#d55740", "#fcf1b0", "#4b80b5"), bias = 0.8)

# Plot boxplot for the top 10 features
ggplot(df_top10, aes(x = scaled_importance, y = variable)) +
  geom_boxplot(
    aes(fill = median_imp), # Map the median importance to fill color
    width = 0.6,
    alpha = 0.6
  ) +
  scale_fill_gradientn(
    colours = color_palette(100) %>% rev(), # Use a 100-color gradient
    guide = "none" # Hide the color legend
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

# Save the boxplot to a PDF
ggsave(filename = "/path/to/data/agri_2_feature_impr_100repeat.pdf", height = 4, width = 6)

# Save the data for the boxplot to a CSV
df_top10 %>%
  as.data.frame() %>%
  write.csv("/path/to/data/agri_2_feature_impr_100repeat_top10feature.csv", row.names = F)

# Step 5: Calculate the relative importance of each feature across seeds
df_m <- df %>% mutate(relative_importance = 100 * (df %>% group_by(seed) %>%
  summarise(relative_importance = scaled_importance / sum(scaled_importance)))[, 2] %>% pull())

# Get the top 10 features based on relative importance
top10_vars <- df_m %>%
  data.frame() %>%
  group_by(variable) %>%
  summarise(median_imp = median(relative_importance)) %>%
  arrange(desc(median_imp)) %>%
  slice_head(n = 10) %>% # Select the top 10 features based on median relative importance
  pull(variable)

df_top10 <- df_m %>%
  data.frame() %>%
  filter(variable %in% top10_vars) %>% # Filter for the top 10 features
  mutate(variable = factor(variable, levels = rev(top10_vars))) # Reverse order for the plot

# Plot the relative importance for the top 10 features
ggplot(df_top10, aes(x = relative_importance, y = variable)) +
  geom_boxplot(
    aes(fill = median_imp), # Map the median importance to fill color
    width = 0.6,
    alpha = 0.6
  ) +
  scale_fill_gradientn(
    colours = color_palette(100) %>% rev(), # Use a 100-color gradient
    guide = "none" # Hide the color legend
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

# Save the boxplot without axis labels to a PDF
ggsave(filename = "/path/to/data/final/agri_2_feature_impr_100repeat_wo_text.pdf", height = 4, width = 3.5)
