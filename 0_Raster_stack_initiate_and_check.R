library(raster)
library(dplyr)

# List all .tif files from the specified directory and subdirectories
file_path_all <- list.files("/path/to/data/GIS_rasters/", pattern = ".tif$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

# Extract the numeric values following "Bio_" in the filenames
nums <- as.numeric(sub(".*Bio_([0-9]+).*", "\\1", file_path_all))

# Sort the files based on the extracted numeric values
file_path_all_ordered <- file_path_all[order(nums)]
file_path_all_ordered

# Create a raster stack by stacking the ordered .tif files
tif_stack <- stack(file_path_all_ordered)

# Save the raster stack as an RDS file for later use (commented out to avoid overwriting)
# saveRDS(tif_stack, "/path/to/data/GIS_rasters/tif_stack_data_194all.rds")

# Alternatively, read the previously saved raster stack from an RDS file
tif_stack <- readRDS("/path/to/data/GIS_rasters/tif_stack_data_197all.rds")
tif_stack

# Check the number of layers in the raster stack
nlayers(tif_stack)

# Optionally, plot the entire raster stack
# plot(tif_stack)

# View the names of the raster layers in the stack
# names(tif_stack)

# Loop through the raster stack layers and perform checks for each one
i <- 1
tif_stack[[i]]

# Plot where the raster cells are equal to 0 (NA values can also be analyzed similarly)
plot(tif_stack[] == 0)

# View the raster stack in a tabular format (optional)
tif_stack %>% View()

# Calculate the minimum values for all layers in the stack
min_values <- cellStats(tif_stack, "min")

# Check raster stacks in different directories by reading all .tif files and stacking them
stack(list.files("/path/to/data/GIS_rasters/1_WorldClim/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/2_climond/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/3_CGIAR-CSI", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/4_USGS/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/5_EarthEnv/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/6_GCAM/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/7_SoilGrids/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/8_Earthdata/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/9_FAO_Animals/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/10_FAO_mapping_Google_Earth/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/11_CGIAR-CSI/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/12_Earthdata/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/13_Dryad/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))
stack(list.files("/path/to/data/GIS_rasters/14_Thomas_AMR_livestock/", pattern = ".tif$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE))

# Load a specific raster file from the ordered list by index
i <- 194
raster_tmp <- raster(file_path_all_ordered[i])
raster_tmp

# Check statistical properties of the raster (min, max, mean, standard deviation, and NA counts)
cellStats(raster_tmp, "min")
cellStats(raster_tmp, "max")
cellStats(raster_tmp, "mean")
cellStats(raster_tmp, "sd")
cellStats(raster_tmp, "countNA")

# Plot the raster showing NA values
plot(is.na(raster_tmp))

# Plot the raster
plot(raster_tmp)

# Plot where the raster values are equal to 128
plot(raster_tmp == 128)

# Modify the raster by setting values equal to 0 or greater than 100 to NA
raster_tmp_ad <- raster_tmp
raster_tmp_ad[raster_tmp_ad == 0] <- NA
raster_tmp_ad[raster_tmp_ad > 100] <- NA

# Plot the modified raster
plot(raster_tmp_ad)

# Check the count of NA values in the modified raster
cellStats(raster_tmp_ad, "countNA")

# Save the modified raster back to the original location, overwriting the previous version
writeRaster(raster_tmp_ad, file_path_all_ordered[i], overwrite = TRUE)
