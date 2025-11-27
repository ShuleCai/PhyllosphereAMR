library(tidyverse)
library(sf)
library(rnaturalearth)
library(cowplot)
library(scales)

# Define risk groups
final_result_risk <- final_result_risk %>%
  mutate(
    risk_group = case_when(
      as.numeric(Risk) <= 3 ~ "Low risk", # Low risk group
      as.numeric(Risk) <= 6 ~ "Medium low risk", # Medium low risk group
      as.numeric(Risk) <= 9 ~ "Medium high risk", # Medium high risk group
      TRUE ~ "High risk" # High risk group
    ),
    risk_group = factor(
      risk_group,
      levels = c(
        "Low risk", "Medium low risk",
        "Medium high risk", "High risk"
      ) # Set risk group order
    )
  )

# Get continent information
world <- ne_countries(scale = "medium", returnclass = "sf")
points_sf <- st_as_sf(final_result_risk, coords = c("x", "y"), crs = 4326) # Convert points to sf object
joined <- st_join(points_sf, world["continent"], join = st_nearest_feature) %>%
  st_drop_geometry() # Remove geometry column

# Calculate proportion of risk groups per continent and per risk level
continent_risk_data <- joined %>%
  drop_na(continent) %>%
  count(continent, Risk, risk_group) %>%
  group_by(continent) %>%
  mutate(
    total_continent = sum(n), # Total points per continent
    prop_group = n / total_continent # Proportion of each risk group per continent
  ) %>%
  group_by(continent, risk_group) %>%
  mutate(
    total_group = sum(n), # Total points per risk group
    prop_risk_in_group = n / total_group # Proportion of each risk within the group
  ) %>%
  ungroup() %>%
  mutate(
    Risk = factor(Risk, levels = 1:12) # Ensure correct ordering of risk levels
  )

# Add global risk data for overall calculation
global_risk_data <- joined %>%
  mutate(continent = "Global") %>% # Create "Global" continent for global data
  count(continent, Risk, risk_group) %>%
  group_by(continent) %>%
  mutate(
    total_continent = sum(n), # Total points globally
    prop_group = n / total_continent # Proportion of each risk group globally
  ) %>%
  group_by(continent, risk_group) %>%
  mutate(
    total_group = sum(n), # Total points in each risk group globally
    prop_risk_in_group = n / total_group # Proportion of each risk level in the global group
  ) %>%
  ungroup() %>%
  mutate(
    Risk = factor(Risk, levels = 1:12) # Ensure correct ordering of risk levels
  )

# Combine global data with continent-level data
all_continent_data <- bind_rows(
  global_risk_data, # Include global risk data
  continent_risk_data # Add continent-specific risk data
)

# Define the list of continents, including "Global"
all_continents <- c("Global", "Africa", "Asia", "Europe", "North America", "South America", "Oceania")

# Determine the maximum y-axis limit for consistency across plots
risk_group_prop <- all_continent_data %>%
  group_by(continent, risk_group) %>%
  summarise(sum = sum(prop_group))
max_prop <- max(risk_group_prop$sum) * 1.00
y_limit <- c(0, max_prop) # Set y-axis limit

# Modify the plotting function (adjust text size for Global label)
create_continent_plot <- function(continent_name, data) {
  plot_data <- data %>%
    filter(continent == continent_name) %>%
    arrange(risk_group, desc(Risk)) %>%
    mutate(Risk_label = as.integer(Risk)) # Add integer label for risk level

  # Adjust text size for "Global"
  text_size <- 2.5

  ggplot(plot_data, aes(x = risk_group, y = prop_group, fill = Risk)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2, position = position_stack(reverse = TRUE)) +
    geom_text(
      aes(label = Risk_label),
      position = position_stack(vjust = 0.5, reverse = T), # Place text labels in the correct stack order
      size = text_size,
      color = "black"
    ) +
    scale_fill_manual(values = extended_spectral) + # Use extended color palette for risk groups
    scale_y_continuous(
      limits = y_limit, # Set y-axis limits
      labels = scales::percent_format(accuracy = 0.1), # Format as percentage
      expand = expansion(c(0, 0.05)) # Add extra space at the top of the bar plot
    ) +
    labs(
      title = ifelse(continent_name == "Global", "Global", continent_name), # Set title based on continent
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        size = 10, # Set title size
        hjust = 0.5,
        face = "bold", # Make title bold
        margin = margin(b = 5) # Add margin below title
      ),
      legend.position = "none", # Remove legend from the plot
      panel.grid.major.x = element_blank(), # Remove vertical grid lines
      panel.grid.major.y = element_line(color = "gray70", linetype = "dashed", size = 0.2), # Add dashed horizontal grid lines
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 45, hjust = 1, size = 8, # Rotate x-axis labels and adjust size
        margin = margin(t = 2, b = 2) # Adjust margin for x-axis labels
      ),
      axis.text.y = element_text(size = 8), # Adjust size of y-axis labels
      panel.border = element_rect(color = "black", fill = NA, size = 0.5), # Add border around panel
      plot.margin = margin(5, 5, 5, 5) # Add margin around the plot
    )
}

# Create plots for all continents (including Global)
all_continent_plots <- map(all_continents, create_continent_plot, data = all_continent_data)

# Combine all plots into a single row
plots_row <- plot_grid(
  plotlist = all_continent_plots,
  nrow = 1,
  align = "hv",
  axis = "lb"
)

# Create a legend plot
legend_plot <- create_legend_plot()
legend_grob <- get_legend(legend_plot)

# Combine the continent plots and legend into a final layout
final_plot <- plot_grid(
  ggdraw(), # Empty plot for positioning
  plots_row, # Continent plots
  legend_grob, # Legend
  ncol = 1,
  rel_heights = c(0.1, 1, 0.15) # Adjust relative heights for each section
)

# Add a unified y-axis label
y_label <- ggdraw() +
  draw_label("Risk Level Proportion", angle = 90, size = 10)

# Final layout combining the y-axis label with the plot
final_plot <- plot_grid(
  y_label,
  final_plot,
  ncol = 2,
  rel_widths = c(0.03, 0.97) # Adjust width ratio for y-axis label and plot
)

# Output the final plot to a PDF file
pdf("/path/to/project/final/Final_3_global_mapping_risk_continents_prop_boxplot.pdf",
  family = "ArialMT",
  width = 14, height = 3.2
)
final_plot
dev.off()
