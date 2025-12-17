# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(readxl)
library(dplyr)
library(igraph)
library(viridis)
library(ggraph)
library(RColorBrewer)

# Read the SARG annotation table
SARG_table <- read.csv("./data/Contig_SARG_annotation.csv")

# Extract unique ARG subtype information from the SARG table
ARG_subtype_bigtable <- SARG_table[, 2:4] %>% unique

# Read ARG abundance data
ARG_subtype_abun <- read.csv("./data/ARG_abundance_tpm.csv")

# Match ARG abundance data with ARG subtype information
ARG_subtype_abun <- ARG_subtype_abun[match(ARG_subtype_bigtable$ARG, ARG_subtype_abun$ARG), ]

# Remove the ARG column and filter columns with non-zero sums
df_t <- ARG_subtype_abun %>% select(-ARG)
df_t <- df_t[, colSums(df_t) > 0]

# Calculate relative abundance for each ARG subtype
ARG_subtype_rela_abun <- df_t %>% apply(MARGIN = 2, function(x){x/sum(x)})

# Compute the mean relative abundance for each ARG subtype
df_mean <- ARG_subtype_rela_abun %>% apply(MARGIN = 1, mean)

# Create a data frame for nodes with ARG subtype information
df_sub <- data.frame(name = ARG_subtype_abun$ARG,
                     size = df_mean) %>% 
  left_join(SARG_table[, 2:3] %>% unique, by = c("name" = "ARG")) %>% 
  dplyr::rename(fill_color = ARG_Category)

# Summarize ARG categories and their sizes
df_cate <- SARG_table[, 2:3] %>% unique %>% mutate(size = df_mean) %>% 
  group_by(ARG_Category) %>% 
  dplyr::summarise(size = sum(size)) %>% 
  dplyr::rename(name = ARG_Category) %>% 
  mutate(fill_color = NA)

# Combine ARG subtype and category data into a single node data frame
df_nodes <- df_sub %>% rbind(df_cate)

# Create edges for the graph based on ARG subtype and category relationships
ARG_edges <- SARG_table[, 2:3] %>% unique
names(ARG_edges) <- c("from", "to")

# Set row names for the node data frame
rownames(df_nodes) <- df_nodes$name

# Create a graph object using the nodes and edges
mygraph <- graph_from_data_frame(ARG_edges, vertices=df_nodes %>% arrange(desc(size)))

# Generate a circular layout for the graph
data = create_layout(mygraph, layout = "circlepack", weight = size, sort.by = NULL, direction = "out", circular = T)

# Assign fill colors to nodes based on their categories
data$fill_color = as.character(data$fill_color)
data$fill_color[is.na(data$fill_color)] = "AA"
data$fill_color = factor(data$fill_color, levels = unique(data$fill_color))

# Define color palettes for the graph
mi = c ("#7CB1A4", "#99BFB2", "#5872A5", "#8aa2ca", "#97D6F3", "#527A72", "#CFE3ED", "#CEC7D9")
mi = c ("#ca6668", "#d23e51", "#5972a2", "#e46a56", "#fbae6a", "#5c5c5c", "#9e2337")
fil_colors <- mi[1 + (0:(nrow(df_cate) - 1)) %% length(mi)]
line_color = "grey50"

# Create a circular graph visualization
g_circular <- ggraph(data) +
  geom_node_circle(aes(fill = as.factor(fill_color), color = as.factor(depth)), size = 0.4) +
  scale_fill_manual(values=c("white", fil_colors)) +
  scale_color_manual(values=c("0" = line_color, "1" = line_color, "2" = line_color, "3" = line_color, "4"=line_color) ) +
  geom_node_text(aes(label=name), size = 3, repel = F) +
  theme_void() +
  theme(legend.position="FALSE")

# Save the circular graph as a PDF file
ggsave(g_circular, filename = "./figures/Fig1a.pdf", family = "ArialMT")


