# Read the metadata and MAG abundance data
metadata_df <- read.csv("/path/to/data/Metadata_phyllosphere.csv")
MAG_abun_A <- read.csv("/path/to/data/models_y_1_original.csv", row.names = 1)

# Join metadata with MAG abundance and log-transformed sum values
metadata_risk <- metadata_df %>% left_join(data.frame(
  Sample = row.names(MAG_abun_A),
  risk_level = MAG_abun_A$cluster,
  Abun_trans = log(MAG_abun_A$sum + 1)
), by = "Sample")

# Load the raster package for spatial data
library(raster)
# tif_stack <- readRDS("/path/to/data/GIS_rasters/tif_stack_data_197all.rds")
tif_stack # Assuming tif_stack is already available in the environment
nlayers(tif_stack)

# Select longitude and latitude columns for spatial analysis
lon_lat <- metadata_risk %>% dplyr::select(lon = longitude, lat = latitude)
coor <- lon_lat
coordinates(coor) <- c("lon", "lat")

# Extract raster values at the selected coordinates
ext <- raster::extract(tif_stack, coor)
res_extract <- cbind(lon_lat, ext)
is.na(res_extract)

# Calculate the percentage of missing (NA) values for each column
na_percentage <- sapply(res_extract, function(x) mean(is.na(x)) * 100)

# Convert the results into a dataframe for easier viewing
result <- data.frame(
  Column = names(na_percentage),
  NA_Percentage = na_percentage
)

# Sort the columns based on NA percentage in descending order
result <- result[order(result$NA_Percentage, decreasing = TRUE), ]

# Display the NA percentage results
print(result)

# Select columns with NA percentage greater than a specific threshold (e.g., 10%)
threshold <- 10
columns_to_remove <- result$Column[result$NA_Percentage > threshold]
print(paste("Columns with NA percentage greater than", threshold, "%:"))
print(columns_to_remove)

# Plot a histogram for a specific raster column (e.g., Bio_177_HumanDevelopmentIndex)
hist(res_extract$Bio_177_HumanDevelopmentIndex, breaks = "Sturges", col = "lightblue", border = "white", main = "Histogram of Continuous Data", xlab = "Value")

# Define columns with zero values
zero_cols <- c(
  "Bio_126_WheatYield", "Bio_127_RiceYield", "Bio_128_MaizeYield",
  "Bio_129_BarleyYield", "Bio_130_SorghumYield", "Bio_131_PearlMilletYield",
  "Bio_132_SmallMilletYield", "Bio_133_SoybeanYield", "Bio_134_TeaYield",
  "Bio_135_PotatoYield", "Bio_136_SweetPotatoYield", "Bio_137_CassavaYield",
  "Bio_138_BeanYield", "Bio_139_TropicalFruitYield", "Bio_140_CocoaYield",
  "Bio_141_TobaccoYield", "Bio_142_SesameseedYield", "Bio_143_CottonYield"
)

# Define columns to be filled with the mean value
mean_cols <- c(
  "Bio_105_SoilMicrobialBiomassNitrogen", "Bio_106_CNRatio", "Bio_107_SoilMicrobialBiomassCarbon",
  "Bio_099_ApatitePhosphorus", "Bio_100_LabileInorganicP", "Bio_101_OccludedP",
  "Bio_102_OrganicP", "Bio_103_SecondaryMineralP", "Bio_104_TotalP",
  "Bio_177_HumanDevelopmentIndex", "Bio_090_BulkDensity", "Bio_091_ClayContent",
  "Bio_092_SandContent", "Bio_093_SiltContent", "Bio_094_CoarseFragments",
  "Bio_095_OrganicCarbonDensity", "Bio_096_pHWater", "Bio_097_Totalnitrogen",
  "Bio_098_SoilOrganicCarbon", "Bio_182_ProportionAntimicrobials", "Bio_124_FeedNitrogenDairyCattle",
  "Bio_113_BuffaloDensity"
)

# Define columns with a single value
one_value_cols <- c(
  "Bio_082_SugarCropRainfed", "Bio_069_C3Arctic", "Bio_060_NeedleleafDeciduousTreeBoreal",
  "Bio_061_BroadleafEvergreenTreeTropical", "Bio_062_BroadleafEvergreenTreeTemperate",
  "Bio_063_BroadleafDeciduousTreeTropical", "Bio_083_SugarCropIrrigated",
  "Bio_132_SmallMilletYield"
)

# Define columns with too many missing values
toomuchna_cols <- c(
  "Bio_113_BuffaloDensity", "Bio_119_DucksDensity", "Bio_121_.FeedNitrogenPigs",
  "Bio_122_FeedNitrogenBeefCattle", "Bio_123_FeedNitrogenChicken", "Bio_124_FeedNitrogenDairyCattle",
  "Bio_182_ProportionAntimicrobials", "Bio_126_WheatYield", "Bio_127_RiceYield",
  "Bio_128_MaizeYield", "Bio_129_BarleyYield", "Bio_130_SorghumYield",
  "Bio_131_PearlMilletYield", "Bio_132_SmallMilletYield", "Bio_133_SoybeanYield",
  "Bio_134_TeaYield", "Bio_135_PotatoYield", "Bio_136_SweetPotatoYield",
  "Bio_137_CassavaYield", "Bio_138_BeanYield", "Bio_139_TropicalFruitYield",
  "Bio_140_CocoaYield", "Bio_141_TobaccoYield", "Bio_142_SesameseedYield",
  "Bio_143_CottonYield", "Bio_044_Coefficient", "Bio_045_Contrast",
  "Bio_046_Correlation", "Bio_047_Dissimilarity", "Bio_048_Entropy",
  "Bio_049_Evenness", "Bio_050_Homogeneity", "Bio_051_Maximum",
  "Bio_052_Range", "Bio_053_Shannon", "Bio_054_Simpson",
  "Bio_055_StandardDeviation", "Bio_056_Uniformity", "Bio_057_Variance"
)

# Define columns to be filled with the mode value
discrete_cols <- c("Bio_144_DevelopmentThreatIndex", "Bio_199_Koppen")

# Function to calculate the mode of a vector
calculate_mode <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  } # Handle all NA case
  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))]
}

# Perform imputation (filling NA values) for various columns
res_extract_imputed <- res_extract %>%
  mutate(
    across(all_of(zero_cols), ~ coalesce(., 0)), # Fill NAs with 0 for zero columns
    across(all_of(mean_cols), ~ if_else( # Fill NAs with mean for mean columns
      is.na(.),
      mean(., na.rm = TRUE),
      .
    )),
    across(all_of(discrete_cols), ~ if_else( # Fill NAs with mode for discrete columns
      is.na(.),
      calculate_mode(.),
      .
    ))
  ) %>%
  dplyr::select(-all_of(columns_to_remove), -all_of(one_value_cols), -all_of(toomuchna_cols))

# Check that there are no more NAs for the specified columns
all(!is.na(res_extract_imputed$Bio_199_Koppen))

# Write the imputed data to a CSV file
res_extract_imputed %>% write.csv("/path/to/data/models_x_1_original.csv", row.names = FALSE)
metadata_risk %>% write.csv("/path/to/data/models_y_1_original.csv", row.names = FALSE)
