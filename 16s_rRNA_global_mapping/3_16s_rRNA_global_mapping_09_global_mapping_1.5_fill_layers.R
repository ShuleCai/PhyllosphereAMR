library(raster)
library(RANN)
library(tidyverse)

# Convert raster stack to data frame including coordinates, removing NA if specified
raster_stack <- as.data.frame(sub_stack_agri, xy = TRUE, na.rm = TRUE)

# Function: calculate the percentage of NA cells for each raster layer
calculate_na_percentage <- function(raster_stack) {
  # Initialize an empty data frame to store NA statistics
  result_df <- data.frame(
    Layer = character(),
    Total_Cells = numeric(),
    NA_Cells = numeric(),
    NA_Percentage = numeric(),
    stringsAsFactors = FALSE
  )

  # Iterate over each layer in the raster stack
  for (i in 1:nlayers(raster_stack)) {
    current_layer <- raster_stack[[i]] # Get the current layer

    total_cells <- ncell(current_layer) # Count total cells
    na_cells <- sum(is.na(getValues(current_layer))) # Count NA cells
    na_percentage <- (na_cells / total_cells) * 100 # Calculate NA percentage

    # Append results to the data frame
    result_df <- rbind(result_df, data.frame(
      Layer = names(current_layer),
      Total_Cells = total_cells,
      NA_Cells = na_cells,
      NA_Percentage = na_percentage
    ))
  }

  return(result_df)
}

# Calculate NA percentages for agricultural and non-agricultural raster stacks
na_stats <- calculate_na_percentage(sub_stack_agri)
na_stats_agri <- na_stats
na_stats_nona <- calculate_na_percentage(sub_stack_nona)

# Display results
print(na_stats)
na_stats %>% View()
plot(is.na(sub_stack_agri[["Bio_197_NDVI"]]))

# Identify column differences between CSV files from different sources
names(read.csv("/path/to/data/agri_models_x_2_feature_selected_iqr.csv")) %>%
  union(names(read.csv("/path/to/data/nonagri_models_x_2_feature_selected_iqr.csv"))) %>%
  setdiff(
    names(read.csv("/path/to/data/agri_models_x_2_feature_selected.csv")) %>%
      union(names(read.csv("/path/to/data/nonagri_models_x_2_feature_selected.csv")))
  ) %>%
  sort()

# ------------------------------------------------------------------------
# Nearest Neighbor Filling for Missing Values in Specific Raster Layers
# ------------------------------------------------------------------------

# Select reference layer (agricultural mask)
reference_layer <- subset(tif_stack, "Bio_011_MeanTemperatureOfColdestQuarter")

# Specify layers that require NA value filling
layers_to_fill <- c(
  "Bio_179_LeafNitrogenConcentration",
  "Bio_180_LeafPhosphorusConcentration",
  "Bio_178_SpecificLeafArea",
  "Bio_181_LeafDryMatterContent"
)

# Function: Fill NA values in raster layers using nearest neighbor interpolation within reference region
fill_na_with_nearest <- function(raster_stack, ref_layer, layers) {
  filled_stack <- raster_stack # Initialize filled raster stack

  # Extract coordinates of valid cells (non-NA) in reference layer
  ref_points <- xyFromCell(ref_layer, which(!is.na(getValues(ref_layer))))

  for (layer_name in layers) {
    current_layer <- subset(raster_stack, layer_name) # Get the current layer

    # Identify cells that need filling: NA cells within reference region
    target_points <- xyFromCell(current_layer, which(
      !is.na(getValues(ref_layer)) & # Agricultural region
        is.na(getValues(current_layer)) # NA in current layer
    ))

    if (nrow(target_points) > 0) {
      # Coordinates of known (non-NA) cells
      known_points <- xyFromCell(current_layer, which(!is.na(getValues(current_layer))))

      # Perform nearest neighbor search using RANN
      nn_idx <- nn2(known_points, target_points, k = 1)$nn.idx

      # Extract values from nearest neighbors
      known_values <- extract(current_layer, known_points)
      target_values <- known_values[nn_idx]

      # Update raster cells with filled values
      current_layer[cellFromXY(current_layer, target_points)] <- target_values
    }

    # Update filled raster stack
    filled_stack[[layer_name]] <- current_layer
  }

  return(filled_stack)
}

# Execute NA filling
filled_stack <- fill_na_with_nearest(tif_stack, reference_layer, layers_to_fill)

# Compare NA percentages before and after filling
original_na_stats <- calculate_na_percentage(sub_stack_agri)
filled_na_stats <- calculate_na_percentage(filled_stack)

comparison <- merge(original_na_stats, filled_na_stats,
  by = "Layer",
  suffixes = c("_original", "_filled")
)
comparison$NA_Reduction <- comparison$NA_Percentage_original - comparison$NA_Percentage_filled
