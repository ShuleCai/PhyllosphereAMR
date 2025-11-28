library(tidyverse)
library(ggplot2)
library(broom)
library(ggpubr)

data <- read.csv("/path/to/data/source_data.csv", check.names = FALSE, row.names = 1)

# Define the columns to plot (Antibiotics)
columns_to_plot <- c(
  "Aminoglycoside", "Bacitracin", "Beta_lactam", "Chloramphenicol",
  "Defensin", "Fosfomycin", "Macrolide-Lincosamide-Streptogramin",
  "Multidrug", "Other_peptide_antibiotics", "Polymyxin",
  "Quinolone", "Rifamycin", "Tetracycline"
)

# Define MGE (Mobile Genetic Elements) column names
MGE <- c(
  "Integration/Excision", "Phage", "Replication/Recombination/Repair",
  "Stability/Transfer/Defense", "Transfer"
)

# Define VFG (Virulence Factors) column names
VFG <- c(
  "Adherence", "Antimicrobial Activity/Competitive Advantage",
  "Biofilm", "Effector Delivery System", "Exoenzyme",
  "Exotoxin", "Immune Modulation", "Invasion", "Motility",
  "Nutritional/Metabolic Factor", "Others", "Regulation",
  "Stress Survival"
)

# Loop through each VFG column to create plots
for (col in VFG) {
  # Perform ANOVA to check differences between groups
  anova_model <- aov(get(col) ~ CR, data = data)

  # Perform Tukey's HSD test for multiple comparisons
  tukey_result <- TukeyHSD(anova_model)

  # Convert the Tukey result to a tidy data frame for easier manipulation
  tukey_df <- tidy(tukey_result)

  # Use the multcompView package to compute significance letters for visual representation
  library(multcompView)
  resultLetters <- multcompLetters4(anova_model, tukey_result)

  # Extract significance letters
  sig_letters <- as.data.frame.list(resultLetters$CR)

  # Define colors for plot fill and points
  color <- c("#5CB85C", "#337AB7", "#F0AD4E")
  value <- c("#85B22E", "#5F80B4", "#E29827")

  # Create the boxplot with jitter points, adding significance markers and error bars
  p <- ggplot(data_s, aes(x = CR, y = get(col), group = CR, fill = CR, color = CR)) +
    geom_jitter(alpha = 0.3, size = 2.5, height = 0, width = 0.2) + # Jitter points for better visibility
    geom_boxplot(alpha = 0.5, size = 1.5, width = 0.7, fill = color) + # Boxplot with fill color
    scale_color_manual(
      limits = c("Class A", "Class B", "Class C"),
      values = value
    ) + # Custom color scale for CR groups
    labs(
      x = "Core resistome type", # x-axis label
      y = col
    ) + # y-axis label based on the current VFG column
    theme_bw() + # Use a clean black-and-white theme
    stat_compare_means(
      comparisons = list(
        c("Class A", "Class B"),
        c("Class A", "Class C"),
        c("Class B", "Class C")
      ),
      label = "p.signif", # Show p-value significance labels
      vjust = 0.5, # Vertical adjustment for label placement
      method = "t.test", # Use t-test for significance comparison
      hide.ns = TRUE, # Hide non-significant results
      bracket.size = 1, # Size of brackets for comparison
      tip.length = 0.0, # Tip length for the brackets
    ) +
    theme(legend.position = "none") # Hide the legend

  # Save the plot as a PDF file for each VFG column
  ggsave(paste0("/path/to/figure/CR_VFG_subtype_bar_sig/", col, "_barplot.pdf"), p, family = "ArialMT")
}
