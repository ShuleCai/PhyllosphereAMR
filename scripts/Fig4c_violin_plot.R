# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggthemes)
library(ggpubr)
library(agricolae)

# Read the data
data_s <- read.csv('./data/MAG_bigtable_w_Clusters.csv')

# Filter and select relevant columns
selected_data <- data_s %>%
  dplyr::select(ARG_count, MGE_count, VFG_count, Cluster) %>%
  na.omit()  # Remove rows with missing values

# Transform data to long format
long_data <- tidyr::gather(selected_data, key = "variable", value = "value", -Cluster)

# Define a function to draw violin plots with significance annotations
draw_violin_plot <- function(data_subset, var_name) {
  # Perform Kruskal-Wallis test and calculate significance groups
  kruskal_result <- kruskal(data_subset$value, data_subset$Cluster, group = TRUE)
  groups <- kruskal_result$groups
  groups$Cluster <- rownames(groups)
  groups$Cluster <- as.factor(groups$Cluster)
  
  # Create violin plot
  p <- ggplot(data_subset, aes(x = factor(Cluster), y = value, fill = factor(Cluster))) +
    geom_violin() +
    scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#00A087")) +
    theme_bw() +
    labs(y = var_name, x = ifelse(var_name == "VFG count", "Cluster", "")) +
    theme(legend.position = "none") +
    geom_text(data = groups, aes(x = Cluster, y = max(data_subset$value) * 1.05, label = groups), size = 4)
  
  return(p)
}

# Generate violin plots for each variable
plot_ARG_count <- draw_violin_plot(long_data %>% filter(variable == "ARG_count"), "ARG count")
plot_MGE_count <- draw_violin_plot(long_data %>% filter(variable == "MGE_count"), "MGE count")
plot_VFG_count <- draw_violin_plot(long_data %>% filter(variable == "VFG_count"), "VFG count")

# Extract legend
legend <- get_legend(ggplot(long_data, aes(x = factor(Cluster), y = value, fill = factor(Cluster))) +
                       geom_violin() +
                       scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#00A087")) +
                       theme_minimal() +
                       theme(legend.position = "bottom"))

# Combine plots and legend
pdf("./figures/Fig4c.pdf", family = "ArialMT")
ggarrange(plot_ARG_count, plot_MGE_count, plot_VFG_count, legend,
          ncol = 1, nrow = 4, heights = c(1, 1, 1, 0.2))
dev.off()
