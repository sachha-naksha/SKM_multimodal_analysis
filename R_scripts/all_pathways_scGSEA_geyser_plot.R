library(Seurat)
library(ggplot2)
library(harmony)
library(Azimuth)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(escape)
library(ggplot2)
library(BiocParallel)
library(rogme)

#Load in the cells
###########################################################################
#combined <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_nonTfh_Tfh_combo.rds")
#combined <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_LN_nonTfh_Tfh_combo.rds")
#combined <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/LN_Tfh_45_Th1_noNaive_noCD8_combo.rds")
combined <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/LN_Tfh_45_Th2_noNaive_noCD8_combo.rds")

#Use this to match your gene lists to make sure the gene is in your dataset before running enrichment
seurat_genes <- rownames(combined)
###########################################################################



#Load in the sets of interest - I had a bunch of different .csv, .gmt, .txt etc
###########################################################################

#SatishRF
SatishRF_list <- read.csv('/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Riberio_pathways/SatishRF.csv', 
                          header = FALSE, stringsAsFactors = FALSE)
SatishRF <- SatishRF_list[[1]]
Restriction_Factors <- intersect(unlist(SatishRF), seurat_genes)

#IL10
IL10_H1_list <- read.csv('/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Riberio_pathways/IL10_H12_FDR5_up.csv', 
                         header = TRUE, stringsAsFactors = FALSE)
IL10<-IL10_H1_list$SYMBOL
IL10 <- intersect(unlist(IL10), seurat_genes)

#IFNI
line_IFNI <- readLines("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Riberio_pathways/Interferome_I_genelist.txt")
all_genes <- lapply(line_IFNI, function(line) {
  fields <- strsplit(line, "\t")[[1]]
  genes <- fields[-c(1, 2)]  # remove first two columns
  genes[genes != ""]         # remove empty entries
})

IFNI <- all_genes
IFNI <- intersect(unlist(IFNI), seurat_genes)

#IFNII
line_IFNII <- readLines("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Riberio_pathways/Interferome_II_genelist.txt")
all_genes <- lapply(line_IFNII, function(line) {
  fields <- strsplit(line, "\t")[[1]]
  genes <- fields[-c(1, 2)]  # remove first two columns
  genes[genes != ""]         # remove empty entries
})

IFNII <- all_genes
IFNII <- intersect(unlist(IFNII), seurat_genes)

#IFNIII
line_IFNIII <- readLines("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Riberio_pathways/Interferome_III_genelist.txt")
all_genes <- lapply(line_IFNIII, function(line) {
  fields <- strsplit(line, "\t")[[1]]
  genes <- fields[-c(1, 2)]  # remove first two columns
  genes[genes != ""]         # remove empty entries
})

IFNIII <- all_genes
IFNIII <- intersect(unlist(IFNIII), seurat_genes)

#Elias 
line_Elias_Haddad <- readLines("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Riberio_pathways/TFH_EliasHaddad.gmt")

all_genes <- lapply(line_Elias_Haddad, function(line) {
  fields <- strsplit(line, "\t")[[1]]
  genes <- fields[-c(1, 2)]  # remove first two columns
  genes[genes != ""]         # remove empty entries
})

Elias_Haddad<- all_genes
Tfh_Markers <- intersect(unlist(Elias_Haddad), seurat_genes)
###########################################################################



#Make a list of all gene lists of interest
###########################################################################
pathway_names <- c("Restriction-Factors", "IL10", "IFNI", "IFNII","IFNIII", "Tfh-Markers")
gene_set <- list(
  Restriction_Factors = Restriction_Factors,
  IL10 = IL10,
  IFNI = IFNI,
  IFNII = IFNII,
  IFNIII = IFNIII,
  Tfh_Markers = Tfh_Markers)
###########################################################################


#Run ssGSEA
###########################################################################
ssg_seurat = runEscape(combined, method = "ssGSEA", gene.sets = gene_set, min_size=2,
                       new.assay.name = "escape.ssGSEA", normalize = TRUE, BPPARAM = SnowParam(workers = 2), alpha=0.75)


DefaultAssay(ssg_seurat)  <- "escape.ssGSEA"

# Reorder in the way you want it to show up on your plot
ssg_seurat$cell_type <- factor(ssg_seurat$cell_type, levels = c("Th2", "Tfh"))

#Gives the comparison for the p-value calculation in the plot
#my_comparison <- list(c("Tfh", "non-Tfh"))
my_comparison <- list(c("Th2", "Tfh"))


# Convert underscores to hyphens to match ssGSEA output (Anything with a hyphen in the pathway gets changed by ssGSEA from _ to -)
valid_gene_sets <- names(gene_set)
converted_names <- gsub("_", "-", valid_gene_sets)

# Check which ones exist in the ssGSEA data
existing_sets <- rownames(ssg_seurat[["escape.ssGSEA"]]@data)
valid_gene_sets <- valid_gene_sets[converted_names %in% existing_sets]
###########################################################################


#Plot GS no cliffs delta for 1 pathway
###########################################################################
gs <- geyserEnrichment(ssg_seurat, 
                       assay = "escape.ssGSEA",
                       group.by = "cell_type",
                       gene.set = "IFNI",
                       palette = 'Spectral')
gs + stat_compare_means(comparisons = my_comparison, method = "wilcox.test", label = "p.signif")
###########################################################################

#Plot GS with Cliffs delta for 1 pathway
###########################################################################


# Reorder the levels of the cell_type metadata
ssg_seurat$cell_type <- factor(ssg_seurat$cell_type, levels = c("Th1", "Tfh"))

# Then regenerate the plot
gs <- geyserEnrichment(ssg_seurat, 
                       assay = "escape.ssGSEA",
                       group.by = "cell_type",
                       gene.set = "IFNI",
                       palette = 'Spectral')

# Add statistical comparisons
gs + stat_compare_means(comparisons = my_comparison, method = "wilcox.test", label = "p.signif")


# Get enrichment matrix from the assay directly
ssgsea_mat <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")

# Combine with metadata
df <- data.frame(
  score = ssgsea_mat["Tfh-Markers", ],
  group = ssg_seurat$cell_type
)
df <- df[df$group %in% c("Th1", "Tfh"), ]
df$group <- factor(df$group, levels = c("Th1", "Tfh"))

cliff <- cidv2(
  df[df$group == "Th2", "score"],
  df[df$group == "Tfh", "score"],
  alpha = 0.05
)[c("d.hat", "p.value")]

# Format label
delta_label <- paste0("Cliff's Δ = ", round(cliff$d.hat, 2), 
                      ", p = ", signif(cliff$p.value, 2))


# Plot
gs <- geyserEnrichment(ssg_seurat, 
                       assay = "escape.ssGSEA",
                       group.by = "cell_type",
                       gene.set = "Tfh-Markers",
                       palette = 'Spectral')

gs + 
  stat_compare_means(comparisons = my_comparison, method = "wilcox.test", label = "p.signif") +
  annotate("text", x = 1.5, y = max(df$score), 
           label = delta_label, vjust = 0.5, size = 5)

###########################################################################

#Plot all pathways, no cliffs delta 
###########################################################################
for (gene_name in valid_gene_sets) {
  cat("Plotting:", gene_name, "\n")
  
  converted_name <- gsub("_", "-", gene_name)
  
  # Generate plot
  gs <- geyserEnrichment(ssg_seurat, 
                         assay = "escape.ssGSEA",
                         group.by = "cell_type",
                         gene.set = converted_name,
                         palette = 'Spectral')
  
  # Add significance test
  gs <- gs + stat_compare_means(comparisons = my_comparison, 
                                method = "wilcox.test", 
                                label = "p.signif")
  
  # Save PDF
  ggsave(filename = paste0("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Gut_scRNA/Results/Non_tfh_Tfh_pathway_analysis/R01_final/LN_geyser_", gene_name, ".pdf"),
         plot = gs,
         device = "pdf",
         width = 7,
         height = 5)
}
###########################################################################
#GS Plot and Save files with Cliffs delta for all pathways
###########################################################################
for (gene_name in valid_gene_sets) {
  cat("Plotting:", gene_name, "\n")
  
  converted_name <- gsub("_", "-", gene_name)
  
  # Extract enrichment matrix from assay
  ssgsea_mat <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")
  
  # Make data frame of scores and group labels
  df <- data.frame(
    score = ssgsea_mat[converted_name, ],
    group = ssg_seurat$cell_type
  )
  df <- df[df$group %in% c("Th1", "Tfh"), ]
  df$group <- factor(df$group, levels = c("Th1", "Tfh"))
  
  # Run Cliff's Delta with cidv2
  cliff <- cidv2(
    df[df$group == "Th1", "score"],
    df[df$group == "Tfh", "score"],
    alpha = 0.05
  )[c("d.hat", "p.value")]
  
  # Create annotation label
  delta_label <- paste0("Cliff's Δ = ", round(cliff$d.hat, 2), 
                        ", p = ", signif(cliff$p.value, 2))
  
  # Generate plot
  gs <- geyserEnrichment(ssg_seurat, 
                         assay = "escape.ssGSEA",
                         group.by = "cell_type",
                         gene.set = converted_name,
                         palette = 'Spectral')
  
  # Add significance test and Cliff's Delta annotation
  gs <- gs + 
    annotate("text", x = -1.5, y = max(df$score, na.rm = TRUE), 
             label = delta_label, vjust = 0.5, size = 5)
  
  # Save PDF
  ggsave(filename = paste0("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Gut_scRNA/Results/Non_tfh_Tfh_pathway_analysis/R01_final/LN_geyser_", gene_name, "cliffs_Th1_Tfh.pdf"),
         plot = gs,
         device = "pdf",
         width = 7,
         height = 5)
}
###########################################################################

#GS Plot w/ Cliffs delta for all pathways and Save summary file, final
###########################################################################
# Prepare results list
results_list <- list()

# Get ssGSEA enrichment matrix once
ssgsea_mat <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")

for (gene_name in valid_gene_sets) {
  cat("Plotting:", gene_name, "\n")
  
  converted_name <- gsub("_", "-", gene_name)
  
  # Check if gene set is present in assay
  if (!(converted_name %in% rownames(ssgsea_mat))) {
    cat("  Skipping - gene set not found\n")
    next
  }
  
  # Build data frame of enrichment scores and group labels
  df <- data.frame(
    score = ssgsea_mat[converted_name, ],
    group = ssg_seurat$cell_type
  )
  df <- df[df$group %in% c("Th2", "Tfh"), ]
  df$group <- factor(df$group, levels = c("Th2", "Tfh"))
  
  # Skip if invalid data
  if (all(is.na(df$score)) || length(unique(df$score)) <= 1) {
    cat("  Skipping - no variability\n")
    next
  }
  
  # Run Cliff's Delta
  cliff <- tryCatch({
    cidv2(
      df[df$group == "Th2", "score"],
      df[df$group == "Tfh", "score"],
      alpha = 0.05
    )[c("d.hat", "p.value")]
  }, error = function(e) {
    cat("  Skipping - cidv2 failed\n")
    return(NULL)
  })
  
  if (is.null(cliff)) next
  
  # Compute group means
  avg_th1 <- mean(df$score[df$group == "Th2"], na.rm = TRUE)
  avg_tfh <- mean(df$score[df$group == "Tfh"], na.rm = TRUE)
  
  # Save result for summary
  results_list[[gene_name]] <- data.frame(
    gene_set = gene_name,
    d.hat = cliff[["d.hat"]],
    p.value = cliff[["p.value"]],
    mean_Th1 = avg_th1,
    mean_Tfh = avg_tfh
  )
  
  # Plot label
  delta_label <- paste0("Cliff's Δ = ", round(cliff[["d.hat"]], 2), 
                        ", p = ", signif(cliff[["p.value"]], 2))
  
  # Generate geyser plot
  gs <- geyserEnrichment(ssg_seurat, 
                         assay = "escape.ssGSEA",
                         group.by = "cell_type",
                         gene.set = converted_name,
                         palette = 'Spectral')
  
  gs <- gs + stat_compare_means(comparisons = my_comparison, 
                                  method = "wilcox.test", 
                                  label = "p.signif") + annotate("text", x = 1.5, y = max(df$score, na.rm = TRUE), 
             label = delta_label, vjust = 0.5, size = 5)
  
  # Save plot
  ggsave(filename = paste0("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/LN_scRNA/Results/pathways_of_interest/Th2/LN_geyser_", gene_name, "_cliffs_Th2_Tfh.pdf"),
         plot = gs,
         device = "pdf",
         width = 7,
         height = 5)
}

# Write summary CSV
if (length(results_list) > 0) {
  summary_df <- do.call(rbind, results_list)
  write.csv(summary_df,
            "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/LN_scRNA/Results/pathways_of_interest/Th2/LN_cliffs_delta_summary_Th2.csv",
            row.names = FALSE)
  cat("Summary CSV saved.\n")
} else {
  cat("No valid results to write to summary.\n")
}

##################################################################################


