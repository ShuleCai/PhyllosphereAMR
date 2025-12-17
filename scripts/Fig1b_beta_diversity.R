rm(list=ls()) # Clear the workspace by removing all objects

# Load necessary libraries
library(vegan)
library(ape)
library(ggplot2)
library(grid)
library(dplyr)
library(multcomp)
library(patchwork)
library(agricolae)

# Read ARG abundance data and set ARG as row names
ARG_abun <- read.csv("./data/ARG_abundance_tpm.csv") %>% tibble::column_to_rownames("ARG")

# Read metadata for metagenome samples
bigtable <- read.csv("./data/Metadata_metagenome.csv")

# Filter ARG abundance data to remove rows with zero sums
ARG_beta <- ARG_abun %>% t %>% data.frame() %>% filter(rowSums(.) > 0)

# Calculate Bray-Curtis dissimilarity matrix
ARG_dist <- vegdist(ARG_beta, method = "bray", binary = F)

# Perform Principal Coordinates Analysis (PCoA)
ARG_pcoa <- cmdscale(ARG_dist, k=2, eig=T)

# Convert PCoA results to a data frame
ARG_pcoa_points <- as.data.frame(ARG_pcoa$points)

# Calculate the percentage of variance explained by each axis
sum_eig <- sum(ARG_pcoa$eig)
eig_percent <- round(ARG_pcoa$eig/sum_eig*100,1)

# Rename PCoA axes
colnames(ARG_pcoa_points) <- paste0("PCoA", 1:2)

# Merge PCoA results with metadata
ARG_pcoa_result <- ARG_pcoa_points %>% mutate(Sample = rownames(.)) %>% left_join(bigtable, by = "Sample")

# Perform PERMANOVA to test differences between groups
set.seed(2)
ARG.div <- adonis2(ARG_beta ~ Agricultural, data = ARG_pcoa_result, permutations = 999, method="bray")
print(ARG.div)

# Format PERMANOVA results for display
ARG.div_adonis <- paste0("adonis R2: ",round(ARG.div$R2[1],2), "; P-value: ", ARG.div$`Pr(>F)`)
print(ARG.div_adonis)

# Define color palette for groups
cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999", "#ADD1E5")

# Group samples by agricultural type and filter out unwanted groups
group_name="Agricultural"
ARG_pcoa_result_g <- ARG_pcoa_result %>% mutate(Group = factor(get(group_name))) %>% filter(Group != "WholeLeaf")

# Calculate maximum values for PCoA axes for each group
yd1 <- ARG_pcoa_result_g %>% group_by(Group) %>% summarise(Max = max(PCoA1))
yd2 <- ARG_pcoa_result_g %>% group_by(Group) %>% summarise(Max = max(PCoA2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1

# Perform ANOVA and Tukey's HSD test for PCoA1
fit1 <- aov(PCoA1~Group,data = ARG_pcoa_result_g)
tuk1<-glht(fit1, linfct=mcp(Group="Tukey"))
res1 <- cld(tuk1,alpha=0.05, decreasing = T)

# Perform ANOVA and Tukey's HSD test for PCoA2
fit2 <- aov(PCoA2~Group,data = ARG_pcoa_result_g)
tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
res2 <- cld(tuk2,alpha=0.05, decreasing = T)

# Perform LSD test for PCoA axes
lsd1 = LSD.test(fit1, "Group", p.adj="BH")
lsd2 = LSD.test(fit2, "Group", p.adj="BH")

# Create a data frame for annotation text in plots
test <- data.frame(PCoA1 = res1$mcletters$Letters,PCoA2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,Group = yd1$Group) %>% 
  left_join(lsd2$groups %>% dplyr::select(PCoA_BH2 = groups) %>% mutate(Group=rownames(.)), by="Group") %>% 
  left_join(lsd1$groups %>% dplyr::select(PCoA_BH1 = groups) %>% mutate(Group=rownames(.)), by="Group")

# Create boxplot for PCoA1
p1 <- ggplot(ARG_pcoa_result_g,aes(Group, PCoA1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd1,label = PCoA_BH1),
            size = 5,color = "black",fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=10),
        axis.text.x=element_blank(),
        legend.position = "none")

# Create boxplot for PCoA2
p3 <- ggplot(ARG_pcoa_result_g,aes(Group,PCoA2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd2,label = PCoA_BH2),
            size = 5,color = "black",fontface = "bold") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=10,angle = 45,
                                 vjust = 1,hjust = 1),
        axis.text.y=element_blank(),
        legend.position = "none")

# Create scatter plot for PCoA1 vs PCoA2
p2 <- ggplot(ARG_pcoa_result_g, aes(PCoA1, PCoA2)) +
  geom_point(aes(fill=Group),size=4,pch=21, alpha=0.9)+
  scale_fill_manual(values=cbbPalette,name = "Group")+
  labs(x = paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y = paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
  ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
  theme(text=element_text(size=10))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=12),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12,vjust = 7, face="bold"),
        axis.title.y=element_text(colour='black', size=12,vjust = -2, face="bold"),
        axis.text=element_text(colour='black',size=10),
        legend.title=element_text(size = 10,face = "bold"),
        legend.text=element_text(size=10),
        legend.key=element_blank(),
        legend.position = c(0.63,0.15),
        legend.background = element_rect(colour = "black"),
        # legend.key.height=unit(1,"cm")
  ) +
  guides(fill = guide_legend(ncol = 1))

# Perform PERMANOVA for the combined plot
otu.adonis=adonis2(ARG_dist~Group,data = ARG_pcoa_result_g,distance = "bray")

# Create a plot for PERMANOVA results
p4 <- ggplot() +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:",  "\nR2 = ",round(otu.adonis[[3]][1],3),  "\np-value < ",otu.adonis[[5]][1],sep = "")),
            size = 3) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

# Combine all plots into a single layout
p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)

# Display the combined plot
p5

# Save the combined plot as a PDF file
ggsave(p5, filename = "./figures/Fig1b.pdf")
