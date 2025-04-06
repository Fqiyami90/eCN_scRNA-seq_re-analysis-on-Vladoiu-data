
####################--------- Integration of E12-P14 Mouse data: re-analysis on Vladoiu's scRNA-seq Data ---------###

#All datasets downloaded on March 2nd, 2025, at 11:28 am (Winnipeg time)
#Farshid Ghiyamihoor
#Analysis performed Between March 2nd to April 5th, 2025

#####################################################################################################################################
############-------------Step 1. Install and load required packages----------###################

install.packages("Seurat")
install.packages("Matrix")  # Required for sparse matrix handling
install.packages('devtools')
devtools::install_github('immunogenomics/presto')
install.packages("dplyr")
install.packages("ggplot2")
install.packages("readxl")
install.packages("stringr")
install.packages("openxlsx")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("patchwork")
install.packages("BiocManager")
BiocManager::install("BiocParallel", type = "binary")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db", force = TRUE)  # for mouse
install.packages("enrichR")
BiocManager::install("ReactomePA")

library(Seurat)
library(Matrix)
library(devtools)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(BiocParallel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichR)
library(ReactomePA)
#####################################################################################################################################
############-------------Step 2: Load and Label Subsetted Seurat eCN Objects for E12, E16, P0, P7 and P14

# Set the working directory
setwd("C:/Users/ghiyamif/OneDrive - University of Manitoba/eCN_project/Datasets/Vladoiu/E12_P14")

e12 <- readRDS("S_Object_Vladoiu_E12_eCN_FG4.rds")
e16 <- readRDS("S_Object_Vladoiu_E16_eCN_FG4.rds")
p0 <- readRDS("S_Object_Vladoiu_P0_eCN_FG4.rds")
p7 <- readRDS("S_Object_Vladoiu_P7_eCN_FG4.rds")
p14 <- readRDS("S_Object_Vladoiu_P14_eCN_FG4.rds")

# Add metadata to distinguish the five datasets
e12$time <- "E12"
e16$time <- "E16"
p0$time <- "P0"
p7$time <- "P7"
p14$time <- "P14"

#####################################################################################################################################
############-------------Step 3: Prepare for Integration

# Create a list of objects to integrate
data_list <- list(E12 = e12, E16 = e16, P0 = p0, P7 = p7, P14 = p14)

# Select features for integration (common highly variable features)
integration_features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 2000)

# (Optional- I did this) Scale and run PCA on each object using the integration features
data_list <- lapply(data_list, function(x) {
  x <- ScaleData(x, features = integration_features, verbose = FALSE)
  x <- RunPCA(x, features = integration_features, verbose = FALSE)
})

#####################################################################################################################################
############-------------Step 4: Integrate Datasets Using reciprocal PCA (RPCA)

# Find integration anchors using the RPCA method
anchors <- FindIntegrationAnchors(object.list = data_list, 
                                  anchor.features = integration_features, 
                                  reduction = "rpca")

# Integrate the datasets to correct for batch effects
integrated_data <- IntegrateData(anchorset = anchors)

# Save the integrated Seurat object to a file
saveRDS(integrated_data, file = "S_Object_Vladoiu_integrated_E12_P14_FG6.rds")

#####################################################################################################################################
############-------------Step 5: Find variable features (optional to get a plot)

integrated_data <- readRDS("S_Object_Vladoiu_integrated_E12_P14_FG6.rds")

#Find variable features
integrated_data <- FindVariableFeatures(integrated_data, selection.method = "vst", nfeatures = 2000)

head(VariableFeatures(integrated_data))

# Plot variable features
plot1 <- VariableFeaturePlot(integrated_data)
# Label the top variable features on the plot (e.g., top 10 by variance)
top10 <- head(VariableFeatures(integrated_data), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#####################################################################################################################################
############-------------Step 6: Dimensionality Reduction and Clustering
#Perform scaling, PCA, and dimensionality reduction (tSNE and UMAP),
#followed by neighborhood graph construction and clustering of the integrated data.

#Switch to the integrated assay for downstream analysis (when you start with reading/reloading the integrated_data)
#Loading the object with readRDS() resets the session-specific settings like the default assay. 
DefaultAssay(integrated_data) <- "integrated"

# Check QC is not necessary as I work on filtered and normalized data for integration.

# Scale the integrated data
integrated_data <- ScaleData(integrated_data, verbose = FALSE)

# Run PCA on the integrated data
integrated_data <- RunPCA(integrated_data, npcs = 30, verbose = FALSE)

# Run tSNE using the top PCs (adjust dims if necessary)
integrated_data <- RunTSNE(integrated_data, reduction = "pca", dims = 1:20)

# If needed (I did not) Run UMAP using the top PCs (adjust dims if necessary)
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:20)

# Compute the nearest neighbor graph and perform clustering
integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.25)

#####################################################################################################################################
############-------------Step 7: Visualize Clusters and Time Points

# Visualize the integrated tSNE colored by time point
DimPlot(integrated_data, reduction = "tsne", group.by = "time") +
  ggtitle("Integrated tSNE: Time Point")

# (I did not) Or Visualize the integrated UMAP colored by time point
DimPlot(integrated_data, reduction = "umap", group.by = "time")

#####################################################################################################################################
############-------------Step 8: Differential Expression Analysis – Global Comparison

#Set Identity to Time Point for Differential Expression
#Set the Identity to the Time Point: This allows you to directly compare cells from P0 versus E16.
Idents(integrated_data) <- integrated_data$time

#Perform DEA: Use FindAllMarkers() to find markers across all time points, comparing each cluster/group to all others:

markers_all_timepoints <- FindAllMarkers(integrated_data, 
                                         test.use = "wilcox", 
                                         min.pct = 0.25, 
                                         logfc.threshold = 0.25)

# Save the markers (RESULTS AFTER DEA) to a separate file
saveRDS(markers_all_timepoints, file = "S_Object_Vladoiu_integrated_E12_P14_markers_all_timepoints_FG7.rds")
# Save the Seurat object (with PCA, tSNE, clusters, etc.)
saveRDS(integrated_data, file = "S_Object_Vladoiu_integrated_E12_P14_with_clusters_FG8.rds")

#When you need: Load the Marker Data
markers_all_timepoints <- readRDS("S_Object_Vladoiu_integrated_E12_P14_markers_all_timepoints_FG7.rds")
#Load the Seurat Object
integrated_data <- readRDS("S_Object_Vladoiu_integrated_E12_P14_with_clusters_FG8.rds")

#####################################################################################################################################
############-------------Step 9: Visualize Top Markers

##Visualize Top 10 Markers (Heatmap)
# For each time point, get the top 10 markers by avg_log2FC
top_markers <- markers_all_timepoints %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Plot a heatmap of the top markers
DoHeatmap(integrated_data, features = top_markers$gene) +
  theme(legend.position = "right") +                                # ensures the expression legend stays
  guides(color = guide_colorbar(), fill = guide_colorbar()) +
  theme(legend.box = "vertical", legend.title = element_text(size = 12)) +
  theme(legend.key.size = unit(1, "cm")) +
  theme(legend.text = element_text(size = 10)) +
  theme(legend.title.align = 1) +
  guides(group = "none")                #this removes the identity legend (the time point bar)

#####################################################################################################################################
############-------------Step 10: Differential Expression Analysis – Pairwise Comparison

#Set Identity to Time Point for Differential Expression (you may did already in previous step)
Idents(integrated_data) <- integrated_data$time

#Compare E12 vs E16:  Use the FindMarkers function to identify differentially expressed genes between the E12 and E16 groups:
markers_E12_vs_E16 <- FindMarkers(integrated_data, ident.1 = "E12", ident.2 = "E16",
                                  test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
#Compare E16 vs P0
markers_E16_vs_P0 <- FindMarkers(integrated_data, ident.1 = "E16", ident.2 = "P0",
                                 test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
#Compare P0 vs P7
markers_P0_vs_P7 <- FindMarkers(integrated_data, ident.1 = "P0", ident.2 = "P7",
                                test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
#Compare P7 vs P14
markers_P7_vs_P14 <- FindMarkers(integrated_data, ident.1 = "P7", ident.2 = "P14",
                                 test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)

#Wilcoxon rank sum test, Minimum 25% of cells expressing the gene, Log-fold change threshold

####--------Visualize Top Pairwise DE Genes  (Change time point for other pairs)

# Extract top 20 genes from each pairwise comparison results (e.g., E12 vs E16) 
top_genes_E12_vs_E16 <- markers_E12_vs_E16 %>%
  arrange(p_val_adj) %>%
  head(20) %>%
  rownames()

top_genes_E16_vs_P0 <- markers_E16_vs_P0 %>%
  arrange(p_val_adj) %>%
  head(20) %>%
  rownames()

top_genes_P0_vs_P7 <- markers_P0_vs_P7 %>%
  arrange(p_val_adj) %>%
  head(20) %>%
  rownames()

top_genes_P7_vs_P14 <- markers_P7_vs_P14 %>%
  arrange(p_val_adj) %>%
  head(20) %>%
  rownames()

# Subset only E12 and E16 cells (repeat for the rest of pairs)
subset_e12_e16 <- subset(integrated_data, idents = c("E12", "E16"))

subset_e16_p0 <- subset(integrated_data, idents = c("E16", "P0"))

subset_p0_p7 <- subset(integrated_data, idents = c("P0", "P7"))

subset_p7_p14 <- subset(integrated_data, idents = c("P7", "P14"))

# Create individual heatmaps
heatmap_E12_E16 <- DoHeatmap(subset_e12_e16, features = top_genes_E12_vs_E16, group.by = "time") +
  guides(color = guide_colorbar(), fill = guide_colorbar())

heatmap_E16_P0 <- DoHeatmap(subset_e16_p0, features = top_genes_E16_vs_P0, group.by = "time") +
  guides(color = guide_colorbar(), fill = guide_colorbar())

heatmap_P0_P7 <- DoHeatmap(subset_p0_p7, features = top_genes_P0_vs_P7, group.by = "time") +
  guides(color = guide_colorbar(), fill = guide_colorbar())

heatmap_P7_P14 <- DoHeatmap(subset_p7_p14, features = top_genes_P7_vs_P14, group.by = "time") +
  guides(color = guide_colorbar(), fill = guide_colorbar())

# Combine into 2x2 grid
combined_heatmaps <- (heatmap_E12_E16 | heatmap_E16_P0) / (heatmap_P0_P7 | heatmap_P7_P14)

# Print combined plot
combined_heatmaps

#####################################################################################################################################
############-------------Step 11: Find Up- and Down-regulated genes as paiwise comparison 

# When comparing two time point: we need to find Up- and Down-regulated genes in later time point only

# Genes upregulated in E16 (compared to E12)
up_E16 <- markers_E12_vs_E16 %>%
  filter(avg_log2FC < 0 & p_val_adj < 0.05)

# Select top 20 most significant upregulated genes
top_up_E16 <- up_E16 %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20) %>%
  rownames()

# Generate heatmap
DoHeatmap(subset_e12_e16, features = top_up_E16, group.by = "time") +
  theme(legend.position = "right")

# Save the upregulated genes to an Excel file
write.xlsx(up_E16, file = "upregulated_markers_E16_vs_12.xlsx")
#----------------------------------------------------------------------------------
# Genes upregulated in P0 (compared to E16)
up_P0 <- markers_E16_vs_P0 %>%
  filter(avg_log2FC < 0 & p_val_adj < 0.05)

# Select top 20 most significant upregulated genes
top_up_P0 <- up_P0 %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20) %>%
  rownames()

# Generate heatmap
DoHeatmap(subset_e16_p0, features = top_up_P0, group.by = "time") +
  theme(legend.position = "right")

write.xlsx(up_P0, file = "upregulated_markers_P0_vs_E16.xlsx")
#----------------------------------------------------------------------------------
# Genes upregulated in P7 (compared to P0)
up_P7 <- markers_P0_vs_P7 %>%
  filter(avg_log2FC < 0 & p_val_adj < 0.05)

# Select top 20 most significant upregulated genes
top_up_P7 <- up_P7 %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20) %>%
  rownames()

# Generate heatmap
DoHeatmap(subset_p0_p7, features = top_up_P7, group.by = "time") +
  theme(legend.position = "right")

write.xlsx(up_P7, file = "upregulated_markers_P7_vs_P0.xlsx")
#----------------------------------------------------------------------------------
# Genes upregulated in P14 (compared to P7)
up_P14 <- markers_P7_vs_P14 %>%
  filter(avg_log2FC < 0 & p_val_adj < 0.05)

# Select top 20 most significant upregulated genes
top_up_P14 <- up_P14 %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20) %>%
  rownames()

# Ensure 'time' is a factor with the desired order (P7 on the left, P14 on the right)
subset_p7_p14$time <- factor(subset_p7_p14$time, levels = c("P7", "P14"))
# Reorder the columns of the data to match the desired factor order
subset_p7_p14 <- subset_p7_p14[, order(subset_p7_p14$time)]

# Generate heatmap
DoHeatmap(subset_p7_p14, features = top_up_P14, group.by = "time") +
  theme(legend.position = "right")

write.xlsx(up_P14, file = "upregulated_markers_P14_vs_P7.xlsx")
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# Genes downregulated in E16 (compared to E12) (i.e., was up in E12)
down_E16 <- markers_E12_vs_E16 %>%
  filter(avg_log2FC > 0 & p_val_adj < 0.05)

# Select top 20
top_down_E16 <- down_E16 %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20) %>%
  rownames()

# Heatmap
DoHeatmap(subset_e12_e16, features = top_down_E16, group.by = "time") +
  theme(legend.position = "right")

write.xlsx(down_E16, file = "downregulated_markers_E16_vs_E12.xlsx")
#-------------------------------------------------------------------------------------
# Genes downregulated in P0 (compared to E16) (i.e., was up in E16)
down_P0 <- markers_E16_vs_P0 %>%
  filter(avg_log2FC > 0 & p_val_adj < 0.05)

# Select top 20
top_down_P0 <- down_P0 %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20) %>%
  rownames()

# Heatmap
DoHeatmap(subset_e16_p0, features = top_down_P0, group.by = "time") +
  theme(legend.position = "right")

write.xlsx(down_P0, file = "downregulated_markers_P0_vs_E16.xlsx")
#-------------------------------------------------------------------------------------
# Genes downregulated in P7 (compared to P0) (i.e., was up in P0)
down_P7 <- markers_P0_vs_P7 %>%
  filter(avg_log2FC > 0 & p_val_adj < 0.05)

# Select top 20
top_down_P7 <- down_P7 %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20) %>%
  rownames()

# Heatmap
DoHeatmap(subset_p0_p7, features = top_down_P7, group.by = "time") +
  theme(legend.position = "right")

write.xlsx(down_P7, file = "downregulated_markers_P7_vs_P0.xlsx")
#-------------------------------------------------------------------------------------
# Genes downregulated in P14 (compared to P7) (i.e., was up in P7)
down_P14 <- markers_P7_vs_P14 %>%
  filter(avg_log2FC > 0 & p_val_adj < 0.05)

# Select top 20
top_down_P14 <- down_P14 %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20) %>%
  rownames()

# Heatmap
DoHeatmap(subset_p7_p14, features = top_down_P14, group.by = "time") +
  theme(legend.position = "right")

write.xlsx(down_P14, file = "downregulated_markers_P14_vs_P7.xlsx")
#---------------#Visualize in combine way ---------------------------
### Combine All Heatmaps into a Grid (using patchwork)
# Re-run your individual heatmaps and assign them to objects

# E12 vs E16
plot_E16 <- DoHeatmap(subset_e12_e16, features = top_up_E16, group.by = "time") +
  ggtitle("Upregulated in E16") +
  theme(legend.position = "right")

# E16 vs P0
plot_P0 <- DoHeatmap(subset_e16_p0, features = top_up_P0, group.by = "time") +
  ggtitle("Upregulated in P0") +
  theme(legend.position = "right")

# P0 vs P7
plot_P7 <- DoHeatmap(subset_p0_p7, features = top_up_P7, group.by = "time") +
  ggtitle("Upregulated in P7") +
  theme(legend.position = "right")

# P7 vs P14
plot_P14 <- DoHeatmap(subset_p7_p14, features = top_up_P14, group.by = "time") +
  ggtitle("Upregulated in P14") +
  theme(legend.position = "right")

# Combine all into 2x2 grid using patchwork
library(patchwork)
plot_E16 + plot_P0 + plot_P7 + plot_P14 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

#####################################################################################################################################
############-------------Step 12: Run GO Enrichment for Each Upregulated Gene Set

#1: Get Entrez IDs from gene symbols
# Example: Convert top upregulated genes in E16 from gene symbols to Entrez IDs
genes_E16_symbols <- rownames(up_E16)

# Use bitr to convert to Entrez
genes_E16_entrez <- bitr(genes_E16_symbols, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Mm.eg.db)

#2: Biological Process: Run GO Enrichment for E16
ego_E16_BP <- enrichGO(gene = genes_E16_entrez$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",           # Biological Process
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE)
#3:Molecular Function (MF)
ego_E16_MF <- enrichGO(gene = genes_E16_entrez$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       ont = "MF",  # Molecular Function
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

#4: Cellular Component (CC)
ego_E16_CC <- enrichGO(gene = genes_E16_entrez$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       ont = "CC",  # Cellular Component
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

# Save enrichment results as CSV files
write.csv(as.data.frame(ego_E16_BP), file = "GO_BP_results_E16.csv")
write.csv(as.data.frame(ego_E16_MF), file = "GO_MF_results_E16.csv")
write.csv(as.data.frame(ego_E16_CC), file = "GO_CC_results_E16.csv")

# Plot GO enrichment (GO trms)
# Dotplot for Biological Process (BP)
dotplot(ego_E16_BP, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in E16_BP")

# Dotplot for Molecular Function (MF)
dotplot(ego_E16_MF, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in E16_MF")

# Dotplot for Cellular Component (CC)
dotplot(ego_E16_CC, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in E16_CC")


#-------------------------------------------------------------------------------------
#1: Get Entrez IDs from gene symbols
genes_P0_symbols <- rownames(up_P0)

# Use bitr to convert to Entrez
genes_P0_entrez <- bitr(genes_P0_symbols, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Mm.eg.db)

#2: Run GO Enrichment for P0
ego_P0_BP <- enrichGO(gene = genes_P0_entrez$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",           # Biological Process
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE)
#3:Molecular Function (MF)
ego_P0_MF <- enrichGO(gene = genes_P0_entrez$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       ont = "MF",  # Molecular Function
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.1,
                       readable = TRUE)

#4: Cellular Component (CC)
ego_P0_CC <- enrichGO(gene = genes_P0_entrez$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       ont = "CC",  # Cellular Component
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

# Save enrichment results as CSV files
write.csv(as.data.frame(ego_P0_BP), file = "GO_BP_results_P0.csv")
write.csv(as.data.frame(ego_P0_MF), file = "GO_MF_results_P0.csv")
write.csv(as.data.frame(ego_P0_CC), file = "GO_CC_results_P0.csv")

# Plot GO enrichment (GO trms)
# Dotplot for Biological Process (BP)
dotplot(ego_P0_BP, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in P0_BP")

# Dotplot for Molecular Function (MF)
dotplot(ego_P0_MF, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in P0_MF")

# Dotplot for Cellular Component (CC)
dotplot(ego_P0_CC, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in P0_CC")
#-------------------------------------------------------------------------------------
#1: Get Entrez IDs from gene symbols
genes_P7_symbols <- rownames(up_P7)

# Use bitr to convert to Entrez
genes_P7_entrez <- bitr(genes_P7_symbols, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Mm.eg.db)

#2: Run GO Enrichment for P7
ego_P7_BP <- enrichGO(gene = genes_P7_entrez$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",           # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   readable = TRUE)
#3:Molecular Function (MF)
ego_P7_MF <- enrichGO(gene = genes_P7_entrez$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "MF",  # Molecular Function
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.1,
                      readable = TRUE)

#4: Cellular Component (CC)
ego_P7_CC <- enrichGO(gene = genes_P7_entrez$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "CC",  # Cellular Component
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      readable = TRUE)

# Save enrichment results as CSV files
write.csv(as.data.frame(ego_P7_BP), file = "GO_BP_results_P7.csv")
write.csv(as.data.frame(ego_P7_MF), file = "GO_MF_results_P7.csv")
write.csv(as.data.frame(ego_P7_CC), file = "GO_CC_results_P7.csv")

# Plot GO enrichment (GO trms)
# Dotplot for Biological Process (BP)
dotplot(ego_P7_BP, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in P7_BP")

# Dotplot for Molecular Function (MF)
dotplot(ego_P7_MF, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in P7_MF")

# Dotplot for Cellular Component (CC)
dotplot(ego_P7_CC, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in P7_CC")
#-------------------------------------------------------------------------------------
#1: Get Entrez IDs from gene symbols
genes_P14_symbols <- rownames(up_P14)

# Use bitr to convert to Entrez
genes_P14_entrez <- bitr(genes_P14_symbols, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Mm.eg.db)

#2: Run GO Enrichment for E16
ego_P14_BP <- enrichGO(gene = genes_P14_entrez$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",           # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   readable = TRUE)
#3:Molecular Function (MF)
ego_P14_MF <- enrichGO(gene = genes_P14_entrez$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "MF",  # Molecular Function
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.1,
                      readable = TRUE)

#4: Cellular Component (CC)
ego_P14_CC <- enrichGO(gene = genes_P14_entrez$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "CC",  # Cellular Component
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      readable = TRUE)

# Save enrichment results as CSV files
write.csv(as.data.frame(ego_P14_BP), file = "GO_BP_results_P14.csv")
write.csv(as.data.frame(ego_P14_MF), file = "GO_MF_results_P14.csv")
write.csv(as.data.frame(ego_P14_CC), file = "GO_CC_results_P14.csv")

# Plot GO enrichment (GO trms)
# Dotplot for Biological Process (BP)
dotplot(ego_P14_BP, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in P14_BP")

# Dotplot for Molecular Function (MF)
dotplot(ego_P14_MF, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in P14_MF")

# Dotplot for Cellular Component (CC)
dotplot(ego_P14_CC, showCategory = 15) + ggtitle("GO Enrichment: Upregulated in P14_CC")
#-------------------------------------------------------------------------------------













