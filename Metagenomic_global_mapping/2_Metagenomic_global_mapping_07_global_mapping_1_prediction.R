# Map predictions to visualize crop risk levels and save the results as a PDF and PNG
# Final predictions map for crop risk

library(ggplot2)
library(raster)

# Load the necessary country and river shapefiles (optional for context)
# countries <- ne_countries(scale = "medium", returnclass = "sf")
# rivers <- ne_download(scale = "medium", type = "rivers_lake_centerlines", category = "physical", returnclass = "sf")

# Visualize the final predictions using ggplot
ggplot() +
  geom_raster(data = final_result, aes(x = x, y = y, fill = mean_predict), na.rm = TRUE) +
  scale_fill_viridis_c(
    name = "Mean Prediction of Risk",
    option = "inferno",
    direction = -1,
    na.value = "lightgray",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(6, "cm"),
      barheight = unit(0.4, "cm"),
      ticks.colour = "lightgray"
    )
  ) +
  coord_sf(ylim = c(-60, 90), expand = FALSE, crs = "+proj=longlat") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_text(size = 10, vjust = 1),
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0.5, size = 12)
  ) +
  labs(title = "Predicted Crop Risk Levels (Transformed Abundance)")

# Save the plot as both PDF and PNG
ggsave("/home/data/t010622/namespace/phyllo/MAG/figures/ML/final/global_mapping_risk.pdf", width = 10, height = 8)
ggsave("/home/data/t010622/namespace/phyllo/MAG/figures/ML/final/global_mapping_risk.png", width = 10, height = 8)

# Optionally save the final prediction results to CSV for reporting
write.csv(final_result, "/home/data/t010622/namespace/phyllo/MAG/figures/ML/final/Final_1_risk_abun_pred_mean_sd.csv", row.names = FALSE)
