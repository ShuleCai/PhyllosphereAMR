# Map predictions to visualize crop risk levels and save the results as a PDF and PNG
# Final predictions map for crop risk

library(ggplot2)
library(raster)

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
ggsave("/path/to/figure/global_mapping_risk.pdf", width = 10, height = 8)
ggsave("/path/to/figure/global_mapping_risk.png", width = 10, height = 8)
