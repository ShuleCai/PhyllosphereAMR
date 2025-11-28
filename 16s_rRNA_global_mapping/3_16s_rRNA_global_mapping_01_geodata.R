# Read 16S metadata including latitude and longitude from the dataset
metadata_df_16S_latlon <- read.csv("/path/to/data//p8_16s_ML_source_data.csv") %>%
  rename(lon = longitude, lat = latitude)

# Transform and log-transform Enterobacteriaceae abundance data in multiple ways for skewness and normality tests
y_values <- metadata_df_16S_latlon$Enterobacteriaceae
y_values <- log(100 * metadata_df_16S_latlon$Enterobacteriaceae + 1) # log-transform with small constant
y_values <- log(100 * metadata_df_16S_latlon$Enterobacteriaceae + 0.000001) # another log transformation with epsilon
y_values <- sqrt(metadata_df_16S_latlon$Enterobacteriaceae) # square root transformation
y_values <- asinh(5 * metadata_df_16S_latlon$Enterobacteriaceae) # inverse hyperbolic sine transformation
y_values <- asinh(sqrt(metadata_df_16S_latlon$Enterobacteriaceae)) # combination of sqrt and asinh transformation

# Compute skewness using the 'moments' package
moments::skewness(y_values)
moments::skewness(log(y_values + 1)) # log transformation with offset
moments::skewness(log(100 * y_values + 0.000001)) # another log transform with a small constant

# Perform Kolmogorov-Smirnov test for normality comparison
ks.test(y_values, "pnorm", mean = mean(y_values), sd = sd(y_values))

# Visualize the distribution of transformed values with a histogram and normal distribution overlay
ggplot(data.frame(y = y_values), aes(x = y)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "lightblue") + # Density histogram
  stat_function(fun = dnorm, args = list(mean = mean(y_values), sd = sd(y_values)), color = "red", size = 1.2) + # Normal curve
  geom_vline(xintercept = mean(y_values), linetype = "dashed", color = "blue") + # Mean line
  labs(x = "ARG abundance", y = "Distribution density") + # Axis labels
  theme_base() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "darkgray"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Load raster data and handle spatial data with sf and terra packages
library(raster)

tif_stack <- readRDS("/path/to/data/tif_stack_data_197all.rds") # Raster data stack
nlayers(tif_stack) # Check number of layers in the raster stack

# Assign geographic coordinates to the metadata (latitude and longitude)
lon_lat_16S <- metadata_df_16S_latlon
coor_16S <- lon_lat_16S
coordinates(coor_16S) <- c("lon", "lat")
coor_16S_sf <- st_as_sf(coor_16S)

# Load raster data using the 'terra' package
tif_terra <- rast(tif_stack)

# Extract values from the raster stack at the coordinates
# ext_16S <- raster::extract(tif_stack, coor_16S)  # Alternative extract method (deprecated)
ext_16S_terra <- terra::extract(tif_terra, coor_16S_sf, method = "neareast") # Efficient extraction with 'terra'

# Combine extracted raster values with the metadata
res_extract_16S <- cbind(lon_lat_16S, ext_16S_terra)

# Identify and count missing values (NA)
is.na(res_extract_16S)
is.na(res_extract_16S) %>% rowSums()

# Calculate the percentage of NA values for each column
na_percentage_16S <- sapply(res_extract_16S, function(x) mean(is.na(x)) * 100)

# Create a data frame to view columns with high NA percentages
result_16S <- data.frame(Column = names(na_percentage_16S), NA_Percentage = na_percentage_16S)

# Sort the result by NA percentage in descending order
result_16S <- result_16S[order(result_16S$NA_Percentage, decreasing = TRUE), ]
print(result_16S %>% arrange(NA_Percentage))

# Identify columns with NA percentage exceeding a threshold (e.g., 15%)
threshold <- 15
columns_to_remove <- result_16S$Column[result_16S$NA_Percentage > threshold]
print(paste("Columns with more than", threshold, "% missing values:"))
print(columns_to_remove)

# Identify columns with a single unique value (i.e., constant columns)
unique_counts <- sapply(res_extract_16S, function(x) length(unique(x)))
single_value_columns <- names(unique_counts[unique_counts == 1])
one_value_cols <- single_value_columns

# Identify columns with excessive missing data
toomuchna_cols <- c("Bio_113_BuffaloDensity", "Bio_119_DucksDensity", "Bio_121_FeedNitrogenPigs", ...)

# Identify columns that should be excluded due to poor performance in models
nogood_cols <- c("Bio_173_PropanilRice", "Bio_174_ThiobencarbRice", "Bio_175_GlyphosateRice", ...)

# Define a function to calculate the mode (for discrete columns)
calculate_mode <- function(x) {
  if (all(is.na(x))) {
    return(NA) # Handle all-NA case
  }
  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))] # Return the most frequent value
}

# Fill missing values in the dataset with appropriate methods (mean for continuous, mode for discrete)
res_extract_imputed_16S <- res_extract_16S %>%
  mutate(
    across(all_of(mean_cols), ~ if_else(is.na(.), mean(., na.rm = TRUE), .)), # Replace NAs with mean for continuous variables
    across(all_of(discrete_cols), ~ if_else(is.na(.), calculate_mode(.), .)) # Replace NAs with mode for discrete variables
  ) %>%
  dplyr::select(-all_of(columns_to_remove), -all_of(one_value_cols), -all_of(toomuchna_cols), -all_of(nogood_cols))

# Check that there are no missing values left for a specific column
all(!is.na(res_extract_imputed_16S$Bio_199_Koppen))

# Export the final dataset to CSV files
res_extract_imputed_16S %>%
  dplyr::mutate(Enterobacteriaceae = 100 * Enterobacteriaceae) %>%
  dplyr::select(Enterobacteriaceae) %>%
  write.csv("/path/to/data/models_y.csv", row.names = F)
res_extract_imputed_16S[, 22:151] %>%
  write.csv("/path/to/data/models_x_1_original.csv", row.names = F)

# Separate agricultural and non-agricultural samples and export results
sample_names_agricultural <- (metadata_df_16S_latlon %>% filter(Agricultural == "Agricultural"))$sample
sample_names_non_agricultural <- (metadata_df_16S_latlon %>% filter(Agricultural != "Agricultural"))$sample

res_extract_imputed_16S_agri <- res_extract_imputed_16S %>% filter(Agricultural == "Agricultural")
res_extract_imputed_16S_agri %>%
  dplyr::mutate(Enterobacteriaceae = 100 * Enterobacteriaceae) %>%
  dplyr::select(sample, lon, lat, Enterobacteriaceae) %>%
  write.csv("/path/to/data/agri_models_y.csv", row.names = F)
res_extract_imputed_16S_agri[, c(1, 3, 4, 22:151)] %>%
  write.csv("/path/to/data/agri_models_x_1_original.csv", row.names = F)

# Repeat for non-agricultural samples
res_extract_imputed_16S_nona <- res_extract_imputed_16S %>% filter(Agricultural == "Non-agricultural")
res_extract_imputed_16S_nona %>%
  dplyr::mutate(Enterobacteriaceae = 100 * Enterobacteriaceae) %>%
  dplyr::select(sample, lon, lat, Enterobacteriaceae) %>%
  write.csv("/path/to/data/nonagri_models_y.csv", row.names = F)
res_extract_imputed_16S_nona[, c(1, 3, 4, 22:151)] %>%
  write.csv("/path/to/data/nonagri_models_x_1_original.csv", row.names = F)

# Handle outliers with IQR and Z-score methods, and export results for both

# Read the imputed dataset and integrate it with the Agricultural classification data
res_extract_imputed_16S <- read.csv("/path/to/data/models_x_1_original.csv") %>%
  mutate(Agricultural = metadata_df_16S_latlon$Agricultural, Enterobacteriaceae = 100 * metadata_df_16S_latlon$Enterobacteriaceae)

# Separate the dataset into agricultural and non-agricultural samples
agri_data <- res_extract_imputed_16S %>% filter(Agricultural == "Agricultural")
nona_data <- res_extract_imputed_16S %>% filter(Agricultural == "Non-agricultural")

# Inspect outlier thresholds using quantiles (for non-agricultural data)
y_values <- nona_data$Enterobacteriaceae
quantile(y_values, 0.75) + 1.5 * (quantile(y_values, 0.75) - quantile(y_values, 0.25)) # Calculate upper bound for outliers using IQR

# Define a function to remove outliers based on IQR (Interquartile Range)
remove_outliers_iqr <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25) # First quartile
  Q3 <- quantile(data[[column]], 0.75) # Third quartile
  IQR_val <- Q3 - Q1 # Calculate IQR
  lower_bound <- Q1 - 1.5 * IQR_val # Lower bound for outlier detection
  upper_bound <- Q3 + 1.5 * IQR_val # Upper bound for outlier detection

  # Return data within the bounds
  return(data[data[[column]] >= lower_bound & data[[column]] <= upper_bound, ])
}

# Define a function to remove outliers based on Z-score (standard deviation)
remove_outliers_zscore <- function(data, column, threshold = 3) {
  mean_val <- mean(data[[column]]) # Mean of the column
  sd_val <- sd(data[[column]]) # Standard deviation of the column
  z_scores <- abs((data[[column]] - mean_val) / sd_val) # Calculate Z-scores

  # Return data within the Z-score threshold
  return(data[z_scores <= threshold, ])
}

# Process agricultural samples using IQR and Z-score outlier removal methods
agri_data <- res_extract_imputed_16S_agri %>%
  dplyr::mutate(Enterobacteriaceae = 100 * Enterobacteriaceae)

# Apply IQR-based outlier removal on agricultural data
agri_iqr <- remove_outliers_iqr(agri_data, "Enterobacteriaceae")
agri_iqr_y <- agri_iqr %>% dplyr::select(sample, lon, lat, Enterobacteriaceae) # Outcome values after outlier removal
agri_iqr_x <- agri_iqr[, c(1, 3, 4, 22:151)] # Explanatory variables after outlier removal

# Apply Z-score-based outlier removal on agricultural data
agri_zscore <- remove_outliers_zscore(agri_data, "Enterobacteriaceae")
agri_zscore_y <- agri_zscore %>% dplyr::select(sample, lon, lat, Enterobacteriaceae) # Outcome values after Z-score filtering
agri_zscore_x <- agri_zscore[, c(1, 3, 4, 22:151)] # Explanatory variables after Z-score filtering

# Summarize the results for the Z-score filtered agricultural data
agri_zscore_y$Enterobacteriaceae %>% summary() # Summary statistics for Enterobacteriaceae after Z-score filtering

# Calculate the proportion of remaining rows after outlier removal
nrow(agri_zscore) / nrow(agri_data) # Proportion of rows kept after Z-score filtering
nrow(agri_iqr) / nrow(agri_data) # Proportion of rows kept after IQR-based filtering

# Export the results for agricultural samples after outlier processing
write.csv(agri_iqr_y, "/path/to/data/agri_models_y_iqr.csv", row.names = F)
write.csv(agri_iqr_x, "/path/to/data/agri_models_x_iqr.csv", row.names = F)
write.csv(agri_zscore_y, "/path/to/data/agri_models_y_zscore.csv", row.names = F)
write.csv(agri_zscore_x, "/path/to/data/agri_models_x_zscore.csv", row.names = F)

# Process non-agricultural samples using the same outlier removal methods
nona_data <- res_extract_imputed_16S_nona %>%
  dplyr::mutate(Enterobacteriaceae = 100 * Enterobacteriaceae)

# Apply IQR-based outlier removal on non-agricultural data
nona_iqr <- remove_outliers_iqr(nona_data, "Enterobacteriaceae")
nona_iqr_y <- nona_iqr %>% dplyr::select(sample, lon, lat, Enterobacteriaceae) # Outcome values after IQR filtering
nona_iqr_x <- nona_iqr[, c(1, 3, 4, 22:151)] # Explanatory variables after IQR filtering

# Apply Z-score-based outlier removal on non-agricultural data
nona_zscore <- remove_outliers_zscore(nona_data, "Enterobacteriaceae")
nona_zscore_y <- nona_zscore %>% dplyr::select(sample, lon, lat, Enterobacteriaceae) # Outcome values after Z-score filtering
nona_zscore_x <- nona_zscore[, c(1, 3, 4, 22:151)] # Explanatory variables after Z-score filtering

# Export the results for non-agricultural samples after outlier processing
write.csv(nona_iqr_y, "/path/to/data/nonagri_models_y_iqr.csv", row.names = F)
write.csv(nona_iqr_x, "/path/to/data/nonagri_models_x_iqr.csv", row.names = F)
write.csv(nona_zscore_y, "/path/to/data/nonagri_models_y_zscore.csv", row.names = F)
write.csv(nona_zscore_x, "/path/to/data/nonagri_models_x_zscore.csv", row.names = F)

# Output the statistics on outlier removal for both agricultural and non-agricultural samples
cat("Agricultural sample outlier statistics:\n")
cat("Outliers removed by IQR method:", nrow(agri_data) - nrow(agri_iqr), "records\n")
cat("Outliers removed by Z-score method:", nrow(agri_data) - nrow(agri_zscore), "records\n\n")

cat("Non-agricultural sample outlier statistics:\n")
cat("Outliers removed by IQR method:", nrow(nona_data) - nrow(nona_iqr), "records\n")
cat("Outliers removed by Z-score method:", nrow(nona_data) - nrow(nona_zscore), "records\n")
