library(h2o)

# Initialize H2O cluster with 5 threads and 80 GB memory
h2o.init(nthreads = 5, max_mem_size = "80G")
# h2o.shutdown()  # Uncomment to shut down H2O cluster when done

# Read the feature matrix X and target variable y for Non-Agricultural data
X <- read.csv("/path/to/data/MAG/figures/ML/nonagri_models_x_2_feature_selected.csv")
y <- read.csv("/path/to/data/MAG/figures/ML/nonagri_models_y_1_original.csv")$Abun_trans

# Combine the feature matrix X and the target variable y into a complete dataset
data <- cbind(X, data.frame(Abun_trans = y))

# Convert the dataset to an H2O frame (automatically handles categorical variables encoding)
h2o_data <- as.h2o(data)

# Loop through specified seed values (in this case, seed = 2)
for (seed in c(2)) {
  # Split the data into training (90%) and testing (10%) sets
  split <- h2o.splitFrame(
    h2o_data,
    ratios = 0.9,
    seed = seed
  )

  # Training set (90% of data) and test set (10% of data)
  train <- split[[1]]
  test <- split[[2]]

  # Run H2O AutoML with specified settings
  aml <- h2o.automl(
    y = "Abun_trans", # Name of the target variable
    training_frame = train, # Training data
    max_models = 2000, # Maximum number of models to train
    max_runtime_secs = 18000, # Maximum runtime in seconds (5 hours)
    stopping_metric = "RMSE", # Metric for early stopping (RMSE for regression tasks)
    sort_metric = "RMSE", # Metric to sort models by
    seed = seed,
    # Specify the algorithms to include in the AutoML run
    include_algos = c("DRF", "XGBOOST", "GBM", "GLM", "StackedEnsemble"),
    verbosity = "info",
    nfolds = 5, # 5-fold cross-validation
    # export_checkpoints_dir = "/path/to/project/AutoML_nonagri/optimal/AutoML_checkpoints"
  )

  # Save the AutoML result for this seed
  saveRDS(aml, paste0("/path/to/project/AutoML_nonagri/optimal/automl_res_seed", seed, "_basicmodels_2000.rds"))

  # Uncomment the following lines to generate SHAP summary plot (interpretability)
  # png(paste0("/path/to/project/AutoML_nonagri/optimal/automl_res_seed", seed, "_shap.png"), family = "ArialMT")
  # h2o.shap_summary_plot(aml@leader, newdata = test)
  # dev.off()

  # Uncomment the following lines to generate variable importance plot for the best model
  # pdf(paste0("/path/to/project/AutoML_nonagri/optimal/automl_res_seed", seed, "_imp.pdf"), family = "ArialMT")
  # h2o.varimp_plot(aml@leader)
  # dev.off()
}
