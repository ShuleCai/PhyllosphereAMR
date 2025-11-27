##################################################################
#                Agricultural                    #
##################################################################
# Read the target variable (Enterobacteriaceae) and explanatory variables for agricultural samples
data_y_table <- read.csv("/path/to/data/16S/ML/agri_models_y_1_original.csv") %>% dplyr::select(-lon, -lat)
data_x_origin <- read.csv("/path/to/data/16S/ML/agri_models_x_1_original.csv") %>%
  tibble::column_to_rownames("sample") %>%
  dplyr::select(-lon, -lat, -Bio_199_Koppen)

# Check if the rownames of data_x_origin match the sample names in data_y_table
identical(rownames(data_x_origin), data_y_table$sample) # Check

# Assign the target variable for agricultural samples
data_y <- data_y_table$Enterobacteriaceae

##################################################################
#                Step 1: Correlation-based Feature Selection    #
##################################################################
# Initialize vectors to store p-values and correlation coefficients (r-values)
p_values <- numeric(ncol(data_x_origin))
names(p_values) <- colnames(data_x_origin)
r_values <- numeric(ncol(data_x_origin))
names(r_values) <- colnames(data_x_origin)

# Calculate the correlation coefficient and p-value for each feature
for (i in 1:ncol(data_x_origin)) {
  test_result <- cor.test(data_x_origin[[i]], data_y, method = "spearman") # Use Spearman correlation
  p_values[i] <- test_result$p.value
  r_values[i] <- test_result$estimate
}

# Display the p-values and r-values in ascending order
data.frame(p = p_values, r = r_values) %>% arrange(p)

# Save the correlation results to a CSV file
data.frame(p = p_values, r = r_values) %>% write.csv("/path/to/data/16S/ML/agri_feature_selection_1_var_cor.csv")

# Set significance level for feature selection
sig_level <- 0.05 # Significance level for p-value

# Feature selection based on p-value threshold (can also apply correlation threshold)
selected_cor <- names(p_values)[which(p_values < sig_level)]

cat("Number of features selected based on significance: ", length(selected_cor), "\n")
selected_cor

# Remove non-selected features
names(p_values) %>% setdiff(selected_cor)

# Select the features and convert integer columns to numeric
data_x_selected_1 <- data_x_origin[selected_cor] %>% mutate_if(is.integer, as.numeric)

##################################################################
#                      Step 2: MRMR Feature Selection           #
##################################################################
library(mRMRe)

# Create the mRMR data object (target is Enterobacteriaceae)
mrmr.data <- mRMR.data(data = data.frame(target = data_y, data_x_selected_1))

# Set seed for reproducibility and define the number of final features
set.seed(4)
final_feature_numer <- 41

# Run mRMR ensemble feature selection
mrmr <- mRMR.ensemble(
  data = mrmr.data, target_indices = 1,
  feature_count = final_feature_numer, solution_count = 1
)

# Extract and sort selected features
selected_features <- mrmr@feature_names[mrmr@filters[[1]]]
selected_features %>% sort()

# Specify columns to remove (if any) based on domain knowledge or other criteria
delete_cols <- c("Bio_091_ClayContent")
data_x_selected_2 <- data_x_selected_1[selected_features %>% setdiff(delete_cols)]

# Check that features were removed correctly
setdiff(names(data_x_selected_1), names(data_x_selected_2))

# Save the mRMR data object and the selected feature dataset
saveRDS(mrmr.data, "/path/to/data/16S/ML/agri_feature_selection_2_mrmr.rds")
write.csv(data_x_selected_2, "/path/to/data/16S/ML/agri_models_x_2_feature_selected.csv", row.names = F)

##################################################################
#                Non-Agricultural                    #
##################################################################
# Read the target variable (Enterobacteriaceae) and explanatory variables for non-agricultural samples
data_y_table <- read.csv("/path/to/data/16S/ML/nonagri_models_y_1_original.csv") %>% dplyr::select(-lon, -lat)
data_x_origin <- read.csv("/path/to/data/16S/ML/nonagri_models_x_1_original.csv") %>%
  tibble::column_to_rownames("sample") %>%
  dplyr::select(-lon, -lat, -Bio_199_Koppen)

# Check if the rownames of data_x_origin match the sample names in data_y_table
identical(rownames(data_x_origin), data_y_table$sample) # Check

# Assign the target variable for non-agricultural samples
data_y <- data_y_table$Enterobacteriaceae

##################################################################
#                Step 1: Correlation-based Feature Selection    #
##################################################################
# Initialize vectors to store p-values and correlation coefficients (r-values)
p_values <- numeric(ncol(data_x_origin))
names(p_values) <- colnames(data_x_origin)
r_values <- numeric(ncol(data_x_origin))
names(r_values) <- colnames(data_x_origin)

# Calculate the correlation coefficient and p-value for each feature
for (i in 1:ncol(data_x_origin)) {
  test_result <- cor.test(data_x_origin[[i]], data_y, method = "spearman") # Use Spearman correlation
  p_values[i] <- test_result$p.value
  r_values[i] <- test_result$estimate
}

# Display the p-values and r-values in ascending order
data.frame(p = p_values, r = r_values) %>% arrange(p)

# Save the correlation results to a CSV file
data.frame(p = p_values, r = r_values) %>% write.csv("/path/to/data/16S/ML/nonagri_feature_selection_1_var_cor.csv")

# Set significance level for feature selection
sig_level <- 0.05 # Significance level for p-value

# Feature selection based on p-value threshold (can also apply correlation threshold)
selected_cor <- names(p_values)[which(p_values < sig_level)]

cat("Number of features selected based on significance: ", length(selected_cor), "\n")
selected_cor

# Remove non-selected features
names(p_values) %>% setdiff(selected_cor)

# Select the features and convert integer columns to numeric
data_x_selected_1 <- data_x_origin[selected_cor] %>% mutate_if(is.integer, as.numeric)

##################################################################
#                      Step 2: MRMR Feature Selection           #
##################################################################
library(mRMRe)

# Create the mRMR data object (target is Enterobacteriaceae)
mrmr.data <- mRMR.data(data = data.frame(target = data_y, data_x_selected_1))

# Set seed for reproducibility and define the number of final features
set.seed(3)
final_feature_numer <- 40

# Run mRMR ensemble feature selection
mrmr <- mRMR.ensemble(
  data = mrmr.data, target_indices = 1,
  feature_count = final_feature_numer, solution_count = 1
)

# Extract and sort selected features
selected_features <- mrmr@feature_names[mrmr@filters[[1]]]
selected_features %>% sort()

# Select the features for the final dataset
data_x_selected_2 <- data_x_selected_1[selected_features]

# Check that the features were correctly selected
setdiff(names(data_x_selected_1), names(data_x_selected_2))

# Save the mRMR data object and the selected feature dataset
saveRDS(mrmr.data, "/path/to/data/16S/ML/nonagri_feature_selection_2_mrmr.rds")
write.csv(data_x_selected_2, "/path/to/data/16S/ML/nonagri_models_x_2_feature_selected.csv", row.names = F)
