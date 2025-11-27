library(h2o)

# Initialize H2O cluster with specified number of threads and memory size
h2o.init(nthreads = 10, max_mem_size = "80G")
# h2o.shutdown()  # Uncomment to shut down the H2O cluster after processing

# Read the feature matrix X and target variable y
X <- read.csv("/path/to/data/agri_models_x_2_feature_selected.csv") # Feature matrix
y <- read.csv("/path/to/data/agri_models_y_1_original.csv")$Enterobacteriaceae # Target variable

# Combine the features and target variable into a complete dataset
data <- cbind(X, data.frame(Enterobacteriaceae = y))

# Convert the data into an H2O frame (automatically handles categorical encoding)
h2o_data <- as.h2o(data)

# Loop to run the model training and evaluation for 50 different seeds
for (seed in 1:50) {
  # Split the data into training (90%) and testing (10%) sets
  split <- h2o.splitFrame(
    h2o_data,
    ratios = 0.9,
    seed = seed
  )
  train <- split[[1]] # 90% training set
  test <- split[[2]] # 10% testing set

  # Run AutoML with specified parameters
  aml <- h2o.automl(
    y = "Enterobacteriaceae", # Name of the target variable
    training_frame = train, # Training dataset
    max_models = 50, # Maximum number of models to run
    max_runtime_secs = 18000, # Maximum runtime (5 hours)
    stopping_metric = "RMSE", # Stopping criterion based on RMSE (for regression tasks)
    sort_metric = "RMSE", # Metric to sort the leaderboard (use RMSE)
    seed = seed,
    include_algos = c("DRF", "XGBOOST", "GBM"), # Algorithms to include in AutoML (Decision Trees, XGBoost, and Gradient Boosting)
    verbosity = "info", # Verbosity level of the output logs
    nfolds = 5, # 5-fold cross-validation
    # export_checkpoints_dir = "/path/to/data/MAG/figures/ML/AutoML_agri/AutoML_checkpoints"  # Uncomment if you want to export checkpoints
  )

  # Save the trained AutoML model results to a file
  saveRDS(aml, paste0("./seeds/automl_res_seed", seed, ".rds"))

  # Plot SHAP summary plot for feature importance
  png(paste0("./seeds/automl_res_seed", seed, "_shap.png"), family = "ArialMT")
  h2o.shap_summary_plot(aml@leader, newdata = test)
  dev.off()

  # Plot variable importance for the best model
  pdf(paste0("./seeds/automl_res_seed", seed, "_imp.pdf"), family = "ArialMT")
  h2o.varimp_plot(aml@leader)
  dev.off()
}
