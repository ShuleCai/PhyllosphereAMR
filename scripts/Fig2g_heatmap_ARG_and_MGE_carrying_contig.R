# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(dplyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(tidyr)
library(extrafont)
library(viridis)
library(tibble)

# Read ARG annotation results (based on contigs)
arg_anno_res <- read.csv("data/ARG_annotation_results.csv")

# Read ARG risk rank table
arg_riskrank <- read.csv("data/SARG_riskrank_table.txt", sep = '\t')

# Read ARG carrying contig (ACC) taxonomy data
acc_tax <- read.csv("data/ARG_carrying_contig_taxonomy.csv")

# Read ARG carrying contig quantification data (TPM)
acc_tpm <- read.csv("data/ARG_carrying_contig_TPM.csv")

# Read ARG and MGE carrying contig information
amcc_table <- read.csv("data/ARG_and_MGE_carrying_contig.csv")

# Read project metadata
project_metadata <- read.csv("data/Metadata_metagenome.csv")

# Merge ARG annotation results with risk rank data
arg_anno_res <- arg_anno_res %>% left_join(arg_riskrank %>% select(variant, rank), by = c("sseqid" = "variant"))

# Filter TPM data for ARG and MGE carrying contigs
amcc_tpm_df <- acc_tpm %>% filter(Contig %in% amcc_table$contig_name)

# Merge ARG annotation results with family taxonomy
arg_anno_res_s_family <- arg_anno_res %>% dplyr::select(contig_name, Type, Subtype, rank) %>% 
  left_join(acc_tax %>% select(contig_name, Family), by = "contig_name")

# Filter TPM data by type and calculate aggregated statistics
amcc_tpm_df_Type <- amcc_tpm_df %>% 
  inner_join(arg_anno_res_s_family %>% filter(!is.na(Family)) %>% select(contig_name, Type) %>% unique, by = c("Contig" = "contig_name"))
amcc_tpm_df_Type_ag <- amcc_tpm_df_Type %>% select(-Type, -Contig) %>% aggregate(by = list(amcc_tpm_df_Type$Type), FUN = sum)

# Calculate frequency and abundance statistics
amcc_tpm_df_Type_stats <- data.frame(Name = amcc_tpm_df_Type_ag$Group.1,
                                      Frequency = ((amcc_tpm_df_Type_ag %>% select(-Group.1)) > 0) %>% rowMeans,
                                      Abundance = amcc_tpm_df_Type_ag %>% select(-Group.1) %>% rowMeans) %>% 
  arrange(desc(Frequency), desc(Abundance))

# Reorder data based on calculated statistics
amcc_tpm_df_Type_ag <- amcc_tpm_df_Type_ag[match(amcc_tpm_df_Type_stats$Name, amcc_tpm_df_Type_ag$Group.1), ]
rownames(amcc_tpm_df_Type_ag) <- amcc_tpm_df_Type_ag$Group.1 %>% tools::toTitleCase(.)
amcc_tpm_df_Type_ag <- amcc_tpm_df_Type_ag %>% select(-Group.1)

# Calculate family-level abundance and rank statistics
amcc_tpm_df_Family <- arg_anno_res_s_family %>% select(contig_name, Family) %>% unique %>% right_join(amcc_tpm_df, by = c("contig_name" = "Contig"))
amcc_tpm_df_Family_abun <- amcc_tpm_df_Family %>% mutate(total_abundance = rowMeans(across(3:ncol(.)))) %>% 
  select(contig_name, Family, total_abundance) %>% 
  left_join(arg_anno_res_s_family %>% select(contig_name, Type) %>% unique, by = c("contig_name" = "contig_name")) %>% 
  filter(!is.na(Family)) %>% 
  select(-contig_name) %>% 
  aggregate(total_abundance ~ Family + Type, FUN = sum) %>% 
  tidyr::pivot_wider(names_from = Family, values_from = total_abundance, values_fill = 0) %>% 
  mutate(Type = tools::toTitleCase(Type)) %>% 
  tibble::column_to_rownames("Type") 
Family_order <- amcc_tpm_df_Family_abun %>% summarise(across(where(is.numeric), sum)) %>% t %>% data.frame %>% rename(sum = ".") %>% 
  tibble::rownames_to_column("Family") %>% arrange(desc(sum)) %>% head(10) %>% select(Family) %>% pull

amcc_Family_rank_abundance <- amcc_tpm_df_Family %>% mutate(total_abundance = rowMeans(across(3:ncol(.)))) %>%
  filter(!is.na(Family)) %>% 
  left_join(arg_anno_res_s_family %>% select(contig_name, rank) %>% unique, by = "contig_name") %>% 
  select(contig_name, Family, rank, total_abundance) %>% 
  aggregate(total_abundance ~ Family + rank, FUN = sum) %>% 
  pivot_wider(names_from = Family, values_from = total_abundance, values_fill = 0) %>% 
  filter(rank != "notassessed")

amcc_Family_rank_richness <- amcc_tpm_df_Family %>% select(contig_name, Family) %>% 
  filter(!is.na(Family)) %>% 
  left_join(arg_anno_res_s_family %>% select(contig_name, Subtype, rank) %>% unique, by = "contig_name") %>% 
  filter(rank != "notassessed") %>% 
  select(-contig_name) %>% 
  unique %>% 
  group_by(rank, Family) %>% 
  summarise(class_type_cnt = n(), .groups = "drop") %>% 
  tidyr::pivot_wider(names_from = Family, values_from = class_type_cnt, values_fill = 0) 
  
amcc_Family_rank_abundance_ordered <- (amcc_Family_rank_abundance %>% tibble::column_to_rownames("rank"))[4:1, Family_order]
amcc_Family_rank_richness_ordered <- (amcc_Family_rank_richness %>% tibble::column_to_rownames("rank"))[4:1, Family_order]
amcc_tpm_df_Family_abun_ordered <- amcc_tpm_df_Family_abun[amcc_tpm_df_Type_stats$Name %>% tools::toTitleCase(), Family_order]

# Generate heatmap annotations
p_heatmap_annotation <- HeatmapAnnotation(`Normalized\nabundance (tpm)` = anno_barplot(amcc_Family_rank_abundance_ordered %>% t %>% as.data.frame,
                                                                                gp = gpar(fill = c("#d87e44", "#2b97e2", "#68d05c", "#dd526c"), lwd = 0.5),
                                                                                height = unit(1.2, "cm"),
                                                                                outline = FALSE,
                                                                                axis = TRUE,
                                                                                bar_width = 0.5),
                                          `Richness` = anno_barplot(amcc_Family_rank_richness_ordered %>% t %>% as.data.frame,
                                                                    gp = gpar(fill = c("#d87e44", "#2b97e2", "#68d05c", "#dd526c"), lwd = 0.5),
                                                                    height = unit(1.2, "cm"),
                                                                    outline = FALSE,
                                                                    axis = T,
                                                                    bar_width=0.5),
                                          show_annotation_name = T,
                                          show_legend = T,
                                          annotation_name_gp = gpar(fontsize = 10)
)

# Prepare data for heatmap plotting
plot_df <- amcc_tpm_df_Family_abun_ordered
plot_df[plot_df==0] <- 0.005
plot_df <- log10(plot_df)

# Create heatmap
p1 <- ComplexHeatmap::pheatmap(data.matrix(plot_df),
                               display_numbers = F,
                               cellwidth = 20,
                               cellheight = 16,
                               scale = "none",
                               angle_col = c("45"),
                               cluster_rows = F,
                               color = rocket(n=101, begin = 1, end = 0.4, direction = 1),
                               breaks = seq(-2.3, 2, length.out = 101),
                               treeheight_col = 20,
                               show_colnames = TRUE,
                               border_color = "grey20",
                               show_rownames = T,
                               legend =  T,
                               legend_breaks = c(-2, -1, 0, 1, 2),
                               name = "Abundance (tpm)",
                               top_annotation = p_heatmap_annotation,
                               cluster_cols = T,
                               fontsize = 10
                               )

# Frequency annotation
col_anno_Fre <- rowAnnotation(Frequency = anno_barplot(amcc_tpm_df_Type_stats$Frequency,
                                                       gp = gpar(fill = c("steelblue")),
                                                       axis_param = list(side = "top")), 
                              border = F,
                              show_annotation_name = TRUE,
                              annotation_name_gp = gpar(fontsize = 10),
                              annotation_name_side = "bottom",
                              annotation_name_rot = 45
)

# Abundance boxplot annotation
plot_abundance_sample <- amcc_tpm_df_Type_ag[match(amcc_tpm_df_Type_stats$Name %>% tools::toTitleCase(), rownames(amcc_tpm_df_Type_ag)), ]
plot_abundance_sample[plot_abundance_sample==0] = 1e-6
plot_abundance_sample_log10 <- log10(plot_abundance_sample)
col_anno_Abun_box <- rowAnnotation(`Abundance (tpm)` = anno_boxplot(plot_abundance_sample_log10 %>% t %>% as.data.frame,
                                                                    outline = FALSE,
                                                                    gp = gpar(fill = c("grey")),
                                                                    axis_param = list(side = "top",
                                                                                      at = c(-6, -4, -2, 0, 2, 4),
                                                                                      labels = c("0", "10^-4", "10^-2", "10^0", "10^2", "10^4")),
                                                                    width = unit(12, "mm")), 
                                   border = F,
                                   show_annotation_name = TRUE,
                                   annotation_name_gp = gpar(fontsize = 10),
                                   annotation_name_side = "bottom",
                                   annotation_name_rot = 45
)

# Ratio plot annotation
amcc_tpm_df_Type_ag %>% names %>% identical(project_metadata$Sample) # check

df_agricultural <- amcc_tpm_df_Type_ag %>% t %>% as.data.frame %>% 
  mutate(Type = project_metadata$Agricultural) 
df_agricultural <- df_agricultural %>% aggregate(. ~ Type, data = ., FUN=mean)  %>%  
  column_to_rownames("Type") %>% t %>% as.data.frame
df_agricultural <- df_agricultural / rowSums(df_agricultural)

col_agi <- rowAnnotation("Agricultural/Non-agricultural" = anno_barplot(as.matrix(df_agricultural),
                                                                        gp = gpar(fill = c("#c02553","#f3c65e") ),
                                                                        axis_param = list(side = 'bottom',
                                                                                          direction = "normal",
                                                                                          at = c(0, 1),
                                                                                          labels = c("Agricultural", "Non-agricultural"),
                                                                                          labels_rot = 45,
                                                                                          gp = gpar(fontsize = 10))),
                         show_annotation_name = FALSE
)

# Combine all annotations and the heatmap
p_patch <- col_anno_Fre + col_anno_Abun_box + col_agi + p1

# Save the final plot to a PDF file
pdf("figures/Fig2g.pdf", family = "ArialMT", width=10)
p_patch
dev.off()
