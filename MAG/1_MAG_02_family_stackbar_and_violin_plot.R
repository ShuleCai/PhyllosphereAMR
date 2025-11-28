library(ggplot2)
library(dplyr)
library(ggpubr)
library(agricolae)

data_s <- read.csv("/path/to/data/source_data.csv") %>%
  mutate(Cluster = ifelse(Cluster == 3, "Class A",
    ifelse(Cluster == 2, "Class B", "Class C")
  ))

# ==============================================================================
# PART 1: Family Composition Stacked Bar Chart by Cluster
# ==============================================================================

# Calculate percentage of each family within each cluster
percentage_data <- data_s %>%
  group_by(Cluster, Family) %>%
  summarise(count = n()) %>%
  group_by(Cluster) %>%
  mutate(percentage = count / sum(count) * 100)

# Combine families with less than 2% abundance into "Other" category
tmp <- percentage_data %>% filter(percentage < 2)
percentage_data_other <- percentage_data %>%
  rbind(list(Cluster = 3, Family = "Other", count = sum(tmp$count), percentage = sum(tmp$percentage))) %>%
  filter(percentage >= 2)

# Define color palette for visualization
nature_colors <- c(
  "#B09C85", "#00A087", "#8491B4", "#F39B7F",
  "#3C5488", "#91D1C2", "#E64B35"
) %>% rev()

# Create and save stacked bar chart showing family composition by cluster
pdf("/path/to/figure/cluster_family_stackbar.pdf", family = "ArialMT", width = 5, height = 6)

ggplot(percentage_data_other, aes(x = as.factor(Cluster), y = percentage, fill = Family)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  coord_flip() +
  theme_par() +
  scale_fill_manual(values = nature_colors) +
  labs(x = "Cluster", y = "Percentage", fill = "Family") +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(color = "black"),
    legend.text = element_text(size = 10),
    legend.position = "bottom"
  )

dev.off()

# ==============================================================================
# PART 2: Violin Plots with Statistical Testing for ARG, MGE, VFG counts
# ==============================================================================

# Select relevant columns and remove missing values
selected_data <- data_s %>%
  dplyr::select(ARG_count, MGE_count, VFG_count, Cluster) %>%
  na.omit()

# Convert data to long format for plotting
long_data <- tidyr::gather(selected_data, key = "variable", value = "value", -Cluster)

# Custom function to create violin plots with significance annotations
draw_violin_plot <- function(data_subset, var_name) {
  # Perform Kruskal-Wallis test and post-hoc grouping
  kruskal_result <- kruskal(data_subset$value, data_subset$Cluster, group = TRUE)
  groups <- kruskal_result$groups
  groups$Cluster <- rownames(groups)
  groups$Cluster <- as.factor(groups$Cluster)

  # Create violin plot
  p <- ggplot(data_subset, aes(x = factor(Cluster), y = value, fill = factor(Cluster))) +
    geom_violin(trim = TRUE, scale = "width") +
    scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#00A087")) +
    theme_bw() +
    labs(y = var_name, x = ifelse(var_name == "VFG count", "Cluster", "")) +
    theme(legend.position = "none") +
    # Add significance annotations
    geom_text(
      data = groups,
      aes(x = Cluster, y = max(data_subset$value) * 1.05, label = groups),
      size = 4
    )

  return(p)
}

# Generate individual violin plots for each variable
plot_ARG_count <- draw_violin_plot(long_data %>% filter(variable == "ARG_count"), "ARG count")
plot_MGE_count <- draw_violin_plot(long_data %>% filter(variable == "MGE_count"), "MGE count")
plot_VFG_count <- draw_violin_plot(long_data %>% filter(variable == "VFG_count"), "VFG count")

# Extract common legend
legend <- get_legend(
  ggplot(long_data, aes(x = factor(Cluster), y = value, fill = factor(Cluster))) +
    geom_violin() +
    scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#00A087"), name = "Cluster") +
    theme_minimal() +
    theme(legend.position = "bottom")
)

# Combine plots and save as PDF
pdf("/path/to/figure/cluster_ARG_MGE_VFG_violin.pdf", family = "ArialMT")
ggarrange(plot_ARG_count, plot_MGE_count, plot_VFG_count, legend,
  ncol = 1, nrow = 4, heights = c(1, 1, 1, 0.2)
)
dev.off()
