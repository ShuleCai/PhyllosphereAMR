# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggthemes)

# Read the data
data_s <- read.csv('./data/MAG_bigtable_w_Clusters.csv')

# Calculate the percentage of each family within each cluster
percentage_data <- data_s %>%
  group_by(Cluster, Family) %>%
  summarise(count = n()) %>%
  group_by(Cluster) %>%
  mutate(percentage = count / sum(count) * 100)

# Combine families with less than 2% into "Other"
tmp <- percentage_data %>% filter(percentage < 2)
percentage_data_other <- percentage_data %>% rbind(list(Cluster = 3, Family = "Other", count = sum(tmp$count), percentage = sum(tmp$percentage))) %>% 
  filter(percentage >= 2)

# Define color palette
nature_colors <- c("#B09C85", "#00A087", "#8491B4", "#F39B7F", 
                   "#3C5488", "#91D1C2", "#E64B35") %>% rev

# Set plot resolution
pdf("./figures/Fig4d.pdf", family = "ArialMT", width = 5, height = 6)

# Create horizontal stacked bar plot
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
