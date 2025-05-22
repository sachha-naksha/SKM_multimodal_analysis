library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(escape)
library(BiocParallel)
library(rogme)
library(msigdbr)
library(effsize)  # for Cliff's delta
library(pheatmap)
library(RColorBrewer)
library(tibble)


#Load in the cells
###########################################################################
#combined <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_nonTfh_Tfh_combo.rds")
#combined <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_LN_nonTfh_Tfh_combo.rds")
combined <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/LN_Tfh_45_Th2_noNaive_noCD8_combo.rds")

###########################################################################

#Load msigDB pathway 
###########################################################################
GS.hallmark <- getGeneSets(library = "H")
#GS.PID <- getGeneSets(library = "C2", subcategory = "CP:PID", species = "Homo sapiens")
names(GS.hallmark)
###########################################################################

#Run analysis
###########################################################################
ssg_seurat = runEscape(combined, method = "ssGSEA", gene.sets = GS.hallmark, min_size=2,
                       new.assay.name = "escape.ssGSEA", normalize = TRUE, BPPARAM = SnowParam(workers = 2), alpha=0.75)

# Reorder the levels of the cell_type metadata
ssg_seurat$cell_type <- factor(ssg_seurat$cell_type, levels = c("Th2", "Tfh"))

#Optional
#ssg_seurat <- performNormalization(ssg_seurat, 
#                                 assay = "escape.ssGSEA", 
#                                 gene.sets = GS.hallmark)

saveRDS(ssg_seurat, file = "/ix/djishnu/Alisa/Tfh/ForPaper/scRNA/scGSEA_hallmark_pathway_analysis.RSD")
#saveRDS(ssg_seurat, file = "/ix/djishnu/Alisa/Tfh/ForPaper/scRNA/scGSEA_PID_pathway_analysis.RSD")
###########################################################################

#Plot heatmap top 20 based on average enrichment score difference
###########################################################################
# Transpose to have cells as rows
scores <- as.data.frame(t(ssg_seurat@assays$escape.ssGSEA@data))
scores$cell_type <- ssg_seurat@meta.data$cell_type  # adjust to your metadata column if needed

# Get average enrichment per group
avg_scores <- scores %>%
  group_by(cell_type) %>%
  summarise(across(where(is.numeric), mean))



# Extract just the enrichment values
group1 <- as.numeric(avg_scores[1, -1])
group2 <- as.numeric(avg_scores[2, -1])
names(group1) <- colnames(avg_scores)[-1]
names(group2) <- colnames(avg_scores)[-1]

# Compute log2FC as a named numeric vector
log2FC <- log2((group1 + 1e-5) / (group2 + 1e-5))

# Get top 20 by absolute value
top_pathways <- names(sort(abs(log2FC), decreasing = TRUE))[1:20]
print(top_pathways)

p <- heatmapEnrichment(ssg_seurat, 
                       group.by = "cell_type",
                       #gene.set.use = rownames(ssg_seurat@assays$escape.ssGSEA@data),
                       gene.set.use = top_pathways,
                       assay = "escape.ssGSEA", 
                       palette = "Spectral")

p + theme(
  axis.text.x = element_text(angle = 90,size=10),
  legend.text = element_text(angle = 90)  # optional: adjust legend text size
)
###########################################################################

#GS plot for 1 pathway of interest w/ out cliffs delta
###########################################################################
print(ssg_seurat@assays$escape.ssGSEA@data)
my_comparison <- list(c("Tfh", "Th2"))

gs <- geyserEnrichment(ssg_seurat, 
                       assay = "escape.ssGSEA",
                       group.by = "cell_type",
                       gene.set = "HALLMARK-IL6-JAK-STAT3-SIGNALING",
                       palette = 'Spectral')
gs + stat_compare_means(comparisons = my_comparison, method = "wilcox.test", label = "p.signif")


# Add significance test
gs <- gs + stat_compare_means(comparisons = my_comparison, 
                              method = "wilcox.test", 
                              label = "p.signif")
gs

ggsave(filename = paste0("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Gut_scRNA/Results/Non_tfh_Tfh_pathway_analysis/R01_final/LN_geyser_SMAD2.pdf"),
       plot = gs,
       device = "pdf",
       width = 7,
       height = 5)
###########################################################################


# Make violin plot
###########################################################################
p <- VlnPlot(ssg_seurat, features =  "HALLMARK-IL6-JAK-STAT3-SIGNALING",  assay = "escape.ssGSEA", group.by = "cell_type", pt.size = 0)

# Add stat bar with asterisk label
p + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
  stat_compare_means(comparisons = my_comparison, method = "wilcox.test", label = "p.signif")

###########################################################################

#Loop through all pathways and make a summary statistic and saves
########################################################################################################
# Extract enrichment matrix
enrichment_matrix <- as.data.frame(t(ssg_seurat[["escape.ssGSEA"]]@data))

# Add metadata for cell type classification
enrichment_matrix$cell_type <- ssg_seurat$cell_type

# Initialize results list
results <- list()

# Loop through each pathway
for (pathway in colnames(enrichment_matrix)[!colnames(enrichment_matrix) %in% "cell_type"]) {
  
  # Subset scores by cell type
  tfh_scores <- enrichment_matrix %>% filter(cell_type == "Tfh") %>% pull(pathway)
  th2_scores <- enrichment_matrix %>% filter(cell_type == "Th2") %>% pull(pathway)
  
  # Compute means
  mean_tfh <- mean(tfh_scores, na.rm = TRUE)
  mean_th2 <- mean(th2_scores, na.rm = TRUE)
  
  # Cliff's delta
  cliff <- tryCatch({
    effsize::cliff.delta(tfh_scores, th2_scores)$estimate
  }, error = function(e) NA)
  
  # Wilcoxon p-value
  pval <- tryCatch({
    wilcox.test(tfh_scores, th2_scores)$p.value
  }, error = function(e) NA)
  
  # Append to results
  results[[pathway]] <- data.frame(
    Pathway = pathway,
    Mean_Tfh = mean_tfh,
    Mean_Th2 = mean_th2,
    Cliffs_Delta = cliff,
    P_Value = pval
  )
}

# Combine and write CSV
summary_df <- do.call(rbind, results)

write.csv(summary_df, "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/LN_scRNA/Results/pid_pathways/pid_ssGSEA_Tfh_vs_Th2_summary.csv", row.names = FALSE)

table(enrichment_matrix$cell_type)

########################################################################################################

#Plot and save Top 20 pathways w/ cliffs delta
########################################################################################################

#summary_df <- read.csv("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/LN_scRNA/Results/hallmark_pathways/hallmark_ssGSEA_Tfh_vs_Th1_summary.csv")

top20 <- summary_df %>%
  mutate(abs_cliff = abs(Cliffs_Delta)) %>%
  arrange(desc(abs_cliff)) %>%
  slice(1:20)

# Define comparisons
my_comparison <- list(c("Tfh", "Th2"))

# Create folder for plots
#dir.create("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/LN_scRNA/Results/hallmark_pathways/plots", showWarnings = FALSE)

# Plot each top pathway
for (i in 1:nrow(top20)) {
  gs_name <- top20$Pathway[i]
  cliff_val <- round(top20$Cliffs_Delta[i], 3)
  p_val <- signif(top20$P_Value[i], 3)
  
  p <- geyserEnrichment(
    ssg_seurat,
    assay = "escape.ssGSEA",
    group.by = "cell_type",
    gene.set = gs_name,
    palette = "Spectral"
  ) +
    stat_compare_means(comparisons = my_comparison, method = "wilcox.test", label = "p.signif") +
    ggtitle(gs_name) +
    annotate("text", x = 1.5, y = Inf, label = paste0("Cliff's Δ = ", cliff_val, 
                                                      "\nP = ", p_val), vjust = 1.5, hjust = 0.5, size = 4)
  
  ggsave(
    filename = paste0("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/LN_scRNA/Results/hallmark_pathways/", gsub("[/: ]", "_", gs_name), "_TEST.pdf"),
    plot = p,
    width = 5,
    height = 4
  )
}

########################################################################################################


#Other type of heatmap...
###################################################################################

heatmap_df <- top20 %>%
  select(Pathway, Mean_Tfh, Mean_Th2)

# Use base R to set rownames
rownames(heatmap_df) <- heatmap_df$Pathway

# Optional: scale each row (pathway)
# heatmap_matrix <- t(scale(t(heatmap_matrix)))

# Drop the Pathway column
heatmap_matrix <- as.matrix(heatmap_df[, c("Mean_Tfh", "Mean_Th2")])


pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
  main = "Top 20 Hallmark Pathways (Mean ssGSEA Score)",
  angle_col = 45,
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 40,
  cellheight = 12,
  palette = "Spectral"
)
###################################################################################