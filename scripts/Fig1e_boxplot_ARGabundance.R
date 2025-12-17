# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(ggpubr)
library(dplyr)

# Read metadata for metagenome samples
bigtable <- read.csv("./data/Metadata_metagenome.csv")

# Read ARG abundance data
ARG_abun <- read.csv("./data/ARG_abundance_tpm.csv")

# Set row names as ARG identifiers
row.names(ARG_abun) <- ARG_abun$ARG

# Remove the ARG column from the abundance data
ARG_abun <- ARG_abun %>% dplyr::select(-ARG)

# Transpose the abundance data and convert to a data frame
ARG_abun_t <- ARG_abun %>% t() %>% as.data.frame()

# Calculate normalized abundance and detection count for each sample
ARG_abun_syth <- data.frame(Sample=rownames(ARG_abun_t),
                            ARG_normalized_abundance=rowSums(ARG_abun_t),
                            ARG_dection_count=(ARG_abun_t>0) %>% rowSums)

# Merge metadata with calculated abundance metrics
bigtable_w_ARG <- bigtable %>% left_join(ARG_abun_syth, by="Sample")

# Define the grouping variable and value to analyze
group="Agricultural"
value="ARG_normalized_abundance"

# Filter and prepare data for analysis
bigtable_w_ARG_g <- bigtable_w_ARG %>% mutate(Value=get(value), Group=get(group)) %>% 
  filter(Group != "Forbidden" & Group != "WholeLeaf")

# Perform ANOVA to test differences between groups
fit <- aov(Value ~ Group, data = bigtable_w_ARG_g)
summary(fit)

# Perform LSD test for pairwise comparisons between groups
var <- agricolae::LSD.test(fit, "Group", p.adj="BH")
var

# Extract p-value and determine significance symbols
p_value = summary(fit)[[1]][[5]]
signif_symbol = ifelse(p_value <= 0.001, "***", ifelse(p_value <= 0.01, "**", ifelse(p_value <= 0.05, "*", "ns")))

# Create a boxplot with jittered points and significance annotations
ggplot(bigtable_w_ARG_g, aes(Group, Value, group = Group, colour = Group)) +
  geom_jitter(alpha = 0.3, size = 2.5) +
  geom_boxplot(alpha= 0.5, size = 1.5, width = 0.7, fill = rev(c("#F0AD4E","#D9534F"))) + 
  scale_color_manual(limits = c("Agricultural", "Non-agricultural"), values = rev(c("#E29827","#922927"))) + 
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, colour = "black") +
  theme_bw() +
  stat_compare_means(
    comparisons = list(c("Agricultural", "Non-agricultural")),
    label = "p.signif",
    vjust = 0.5,
    method = "t.test",
    hide.ns = TRUE,
    bracket.size = 1,
    tip.length = 0.0,
    size = 5
  ) +   
  labs(title = NULL, x = "", y = "ARG abundance (tpm)") + 
  theme(legend.position = "none",
        plot.title = element_text(size = 15,
                                  hjust = 0.5,
                                  colour = "black",
                                  face = "bold"),
        axis.title.y = element_text(size = 12,
                                    color = "black",
                                    vjust = 1.9,
                                    hjust = 0.5,
                                    angle = 90),
        legend.title = element_text(size = 12,
                                    color = "black"),
        legend.text = element_text(size = 10,
                                   color = "black"),
        axis.text.x = element_text(size = 10,
                                   color = "black",
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 0),
        axis.text.y = element_text(size = 10,
                                   color = "black",
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 0))

# Save the plot as a PDF file
ggsave(filename = "./figures/Fig1e.pdf",
       family = "ArialMT", width = 4, height = 3.5)
