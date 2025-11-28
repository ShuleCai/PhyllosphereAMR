# Read the metadata and sample data
metadata_df <- read.csv("/path/to/data/Metadata_phyllosphere.csv")
MAG_abun_A <- read.csv("/path/to/data/ML/models_y_1_original.csv", row.names = 1)

# Split samples into Agricultural and Non-Agricultural
sample_names_agricultural <- (bigtable %>% filter(Agricultural == "Agricultural"))$Sample
sample_names_non_agricultural <- (bigtable %>% filter(Agricultural != "Agricultural"))$Sample

##################################################################
#                Agricultural Data Processing                   #
##################################################################

# Filter data for Agricultural samples
data_y_table <- read.csv("/path/to/data/ML/models_y_1_original.csv") %>%
  filter(Sample %in% sample_names_agricultural)

# Read and preprocess X data for Agricultural samples
data_x_origin <- read.csv("/path/to/data/ML/models_x_1_original.csv", row.names = read.csv("/path/to/data/ML/models_y_1_original.csv")$Sample) %>%
  dplyr::select(-lon, -lat, -Bio_199_Koppen) %>%
  dplyr::filter(row.names(.) %in% sample_names_agricultural)

# Check if row names of X and Y match
identical(rownames(data_x_origin), data_y_table$Sample) # check

# Define the response variable
data_y <- data_y_table$Abun_trans

##################################################################
#                Step 1: Correlation-Based Feature Selection     #
##################################################################

# Initialize vectors to store p-values and correlation coefficients (r values)
p_values <- numeric(ncol(data_x_origin))
names(p_values) <- colnames(data_x_origin)
r_values <- numeric(ncol(data_x_origin))
names(r_values) <- colnames(data_x_origin)

# Calculate correlation coefficients and p-values for each feature
for (i in 1:ncol(data_x_origin)) {
  test_result <- cor.test(data_x_origin[[i]], data_y, method = "spearman") # Spearman correlation
  p_values[i] <- test_result$p.value
  r_values[i] <- test_result$estimate
}

# Save the correlation results to a CSV file
data.frame(p = p_values, r = r_values) %>% write.csv("/path/to/data/ML/agri_feature_selection_1_var_cor.csv")

# Sort the results by p-value
data.frame(p = p_values, r = r_values) %>% arrange(p)

# Set significance level for feature selection
sig_level <- 0.05 # Significance level

# Select features based on significance level (p-value)
selected_cor <- names(p_values)[which(p_values < sig_level)]
cat("Number of features selected by significance test:", length(selected_cor), "\n")
selected_cor

# Remove unselected features and convert integer columns to numeric
data_x_selected_1 <- data_x_origin[selected_cor] %>% mutate_if(is.integer, as.numeric)

##################################################################
#                Step 2: MRMR Feature Selection                  #
##################################################################

# Load mRMRe package for MRMR feature selection
library(mRMRe)

# Prepare data for MRMR
mrmr.data <- mRMR.data(data = data.frame(target = data_y, data_x_selected_1))

# Perform MRMR feature selection with a specified number of features
set.seed(3)
final_feature_number <- 60
mrmr <- mRMR.ensemble(
  data = mrmr.data, target_indices = 1,
  feature_count = final_feature_number, solution_count = 1
)

# Get the selected features and sort them
selected_features <- mrmr@feature_names[mrmr@filters[[1]]]
selected_features %>% sort()

# Create a new dataset with selected features
data_x_selected_2 <- data_x_selected_1[selected_features]
setdiff(names(data_x_selected_1), names(data_x_selected_2))

# Save the MRMR results and the selected features to CSV
saveRDS(mrmr.data, "/path/to/data/ML/agri_feature_selection_2_mrmr.rds")
write.csv(data_x_selected_2, "/path/to/data/ML/agri_models_x_2_feature_selected.csv", row.names = FALSE)

##################################################################
#                Non-Agricultural Data Processing                #
##################################################################

# Filter data for Non-Agricultural samples
data_y_table <- read.csv("/path/to/data/ML/models_y_1_original.csv") %>%
  filter(Sample %in% sample_names_non_agricultural)

# Read and preprocess X data for Non-Agricultural samples
data_x_origin <- read.csv("/path/to/data/ML/models_x_1_original.csv", row.names = read.csv("/path/to/data/ML/models_y_1_original.csv")$Sample) %>%
  dplyr::select(-lon, -lat, -Bio_199_Koppen) %>%
  dplyr::filter(row.names(.) %in% sample_names_non_agricultural)

# Remove specific columns for Non-Agricultural data
delete_cols <- c(
  "Bio_156_NitrogenFertilizerApplication", "Bio_157_NitrogenManureProduction",
  "Bio_158_PhosphorusFertilizerApplication", "Bio_159_PhosphorusManureProduction",
  "Bio_161_GlyphosateCorn", "Bio_162_AtrazineCorn", "Bio_163_AcetochlorCorn",
  "Bio_164_GlyphosateSoybean", "Bio_165_PendimethalinSoybean", "Bio_166_TrifluralinSoybean",
  "Bio_167_GlyphosateWheat", "Bio_168_Wheat2.4.d", "Bio_169_McpaWheat", "Bio_170_GlyphosateCotton",
  "Bio_171_DichloropropeneCotton", "Bio_172_TrifluralinCotton"
)

# Remove unwanted columns
data_x_origin <- data_x_origin %>% dplyr::select(-all_of(delete_cols))

##################################################################
#                Step 1: Correlation-Based Feature Selection     #
##################################################################

# Reinitialize p-values and correlation coefficients for Non-Agricultural data
p_values <- numeric(ncol(data_x_origin))
names(p_values) <- colnames(data_x_origin)
r_values <- numeric(ncol(data_x_origin))
names(r_values) <- colnames(data_x_origin)

# Calculate correlation coefficients and p-values for each feature
for (i in 1:ncol(data_x_origin)) {
  test_result <- cor.test(data_x_origin[[i]], data_y, method = "spearman") # Spearman correlation
  p_values[i] <- test_result$p.value
  r_values[i] <- test_result$estimate
}

# Save the correlation results to a CSV file
data.frame(p = p_values, r = r_values) %>% write.csv("/path/to/data/ML/nonagri_feature_selection_1_var_cor.csv")

# Sort the results by p-value
data.frame(p = p_values, r = r_values) %>% arrange(p)

# Set significance level for feature selection
sig_level <- 0.05 # Significance level

# Select features based on significance level (p-value)
selected_cor <- names(p_values)[which(p_values < sig_level)]
cat("Number of features selected by significance test:", length(selected_cor), "\n")
selected_cor

# Remove unselected features and convert integer columns to numeric
data_x_selected_1 <- data_x_origin[selected_cor] %>% mutate_if(is.integer, as.numeric)

##################################################################
#                Step 2: MRMR Feature Selection                  #
##################################################################

# Prepare data for MRMR
mrmr.data <- mRMR.data(data = data.frame(target = data_y, data_x_selected_1))

# Perform MRMR feature selection with a specified number of features
set.seed(3)
final_feature_number <- 60
mrmr <- mRMR.ensemble(
  data = mrmr.data, target_indices = 1,
  feature_count = final_feature_number, solution_count = 1
)

# Get the selected features and sort them
selected_features <- mrmr@feature_names[mrmr@filters[[1]]]
selected_features %>% sort()

# Create a new dataset with selected features
data_x_selected_2 <- data_x_selected_1[selected_features]
setdiff(names(data_x_selected_1), names(data_x_selected_2))

# Save the MRMR results and the selected features to CSV
saveRDS(mrmr.data, "/path/to/data/nonagri_feature_selection_2_mrmr.rds")
write.csv(data_x_selected_2, "/path/to/data/nonagri_models_x_2_feature_selected.csv", row.names = FALSE)
