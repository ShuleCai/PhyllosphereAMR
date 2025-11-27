library(h2o)

# Initialize the H2O cluster with 10 threads and 80 GB of memory
h2o.init(nthreads = 10, max_mem_size = "80G")
# h2o.shutdown()  # Uncomment to shut down H2O cluster when done

# Read the feature matrix X and target variable y
X <- read.csv("/path/to/data/MAG/figures/ML/models_x_2_feature_selected.csv")
y <- read.csv("/path/to/data/MAG/figures/ML/models_y_1_original.csv")$Abun_trans

# Combine the features and target variable into a complete dataset
data <- cbind(X, data.frame(Abun_trans = y))

# Convert the data to an H2O frame (automatically handles categorical variables encoding)
h2o_data <- as.h2o(data)

# Loop over seeds for model training
for (seed in 1:100) {
  # Split the data into 90% training and 10% test sets
  split <- h2o.splitFrame(
    h2o_data,
    ratios = 0.9,
    seed = seed
  )
  train <- split[[1]] # 90% training data
  test <- split[[2]] # 10% test data

  # Train an AutoML model with the specified settings
  aml <- h2o.automl(
    y = "Abun_trans", # Name of the target variable column
    training_frame = train, # Training data
    max_models = 20, # Maximum number of models to train
    max_runtime_secs = 18000, # Maximum runtime (5 hours)
    stopping_metric = "RMSE", # Early stopping metric (RMSE for regression)
    sort_metric = "RMSE", # Metric for sorting models
    seed = seed,
    # Specify the algorithms to include in the AutoML run
    include_algos = c("DRF", "XGBOOST", "GBM"),
    verbosity = "info",
    nfolds = 5, # 5-fold cross-validation
    # export_checkpoints_dir = "/path/to/project/AutoML_Full/AutoML_checkpoints"
  )

  # Save the AutoML model result for this seed
  saveRDS(aml, paste0("/path/to/project/AutoML_Full/seeds/automl_res_seed", seed, ".rds"))

  # Generate SHAP summary plot for model interpretability and save as PNG
  png(paste0("/path/to/project/AutoML_Full/seeds/automl_res_seed", seed, "_shap.png"), family = "ArialMT")
  h2o.shap_summary_plot(aml@leader, newdata = test)
  dev.off()

  # Generate variable importance plot for the best model and save as PDF
  pdf(paste0("/path/to/project/AutoML_Full/seeds/automl_res_seed", seed, "_imp.pdf"), family = "ArialMT")
  h2o.varimp_plot(aml@leader)
  dev.off()
}
