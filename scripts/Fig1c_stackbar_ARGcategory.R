# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Read AMR category abundance data
ARG_type <- read.csv("./data/AMRcategory_abundance_tpm.csv")

# Set row names as ARG categories and remove the category column
rownames(ARG_type) <- ARG_type$ARG_Category
ARG_type <- ARG_type %>% dplyr::select(-ARG_Category)

# Read metadata for metagenome samples
bigtable <- read.csv("./data/Metadata_metagenome.csv")

# Normalize ARG abundance data by dividing each value by the column sum
ARG_abun_normalized <- ARG_type %>% select_if(~sum(.)!=0) %>% mutate(across(everything(), ~ ./sum(.)))

# Filter and arrange metadata to match normalized ARG abundance data
bigtable_s <- bigtable %>% filter(Sample %in% names(ARG_abun_normalized)) %>% 
  arrange(match(Sample, names(ARG_abun_normalized))) %>% mutate(Group = Agricultural)

# Select top 9 ARG categories based on mean abundance
names_top <- names(ARG_abun_normalized %>% as.data.frame %>% apply(1, FUN=mean) %>% sort(decreasing = T))[1:9]

# Combine remaining categories into "Others"
names_diff <- setdiff(rownames(ARG_abun_normalized), names_top)
class_df <- rbind(ARG_abun_normalized[names_top, ], Others = apply(ARG_abun_normalized[names_top, ], 2, function(x){1 - sum(x)}))

# Define color palette for the plot
values <- rev(c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC6666","#9999CC","#66CC99", "#ADD1E5"))

# Add group information to the metadata
df_w_type <- bigtable_s %>% mutate(Type = Agricultural)

# Filter metadata to include only samples present in the abundance data
names <- intersect(df_w_type$Sample, colnames(class_df))
df_w_type <- df_w_type %>% filter(Sample %in% names)

# Subset abundance data to include only filtered samples
class_final <- class_df[, names]

# Reshape abundance data for plotting
class.2 <- class_final
class.2 <- cbind(ClassID = rownames(class.2), class.2) %>% data.frame()
class.2$ClassID <- factor(class.2$ClassID, levels = rev(class.2$ClassID))

# Melt the data for ggplot
library(reshape2)
class.gg <- melt(class.2, id.vars = "ClassID", variable.name = "SampleID", value.name = "Abundance")

# Merge melted data with metadata
class.gg.group <- class.gg %>% merge(df_w_type, by.x = "SampleID", by.y = "Sample", all.x = TRUE) %>%
  mutate(Abundance = as.numeric(Abundance))

# Create a stacked bar plot for ARG categories
ggplot(class.gg.group %>% dplyr::rename(`ARG category` = ClassID), aes(Type, 100 * Abundance, fill = `ARG category`)) +
  geom_bar(stat = "summary_bin", fun = mean, position = 'stack', 
           width = 0.5, size = 0.3, color = "black",
           show.legend = TRUE) +
  scale_y_continuous(expand = c(0,0))+
  labs(x = '', y = 'Relative abundance (%)')+
  scale_fill_manual(values = values) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10))

# Save the plot as a PDF file
ggsave("./figures/Fig1c.pdf", width = 4, height = 4)

