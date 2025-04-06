
###############------------Mouse P0: re-analysis on Vladoiu's scRNA-seq Data----------------------######################

##The single-cell RNA sequencing (scRNA-seq) dataset for mouse developmental time point embryonic day 16 (E16, N=9) has been taken from:
#a published article from European Genome-phenome Archive (EGA) repository with accession (GSE118068).
##Dataset includes the three following files which is typical for scRNA-seq data:
#1. Gene expression matrix (.mtx) – contains gene expression counts
#2.Genes/features file (.tsv or .txt) – contains gene names
#3.Barcodes/cell file (.tsv or .txt) – contains cell identifiers

#All datasets downloaded on March 2nd, 2025, at 11:28 am (Winnipeg time)
#Farshid Ghiyamihoor
#Analysis performed Between March 2nd to April 5th, 2025

#####################################################################################################################################
############-------------Step 1. Install and load required packages

install.packages("Seurat")
install.packages("Matrix")  # Required for sparse matrix handling
install.packages('devtools')
devtools::install_github('immunogenomics/presto')
install.packages("dplyr")
install.packages("ggplot2")
install.packages("readxl")
install.packages("stringr")
install.packages("openxlsx")

library(Seurat)
library(Matrix)
library(devtools)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(openxlsx)

#####################################################################################################################################
############-------------Step 2. Read the matrix file and assign barcodes to expression matrix

# Set the working directory
setwd("C:/Users/ghiyamif/OneDrive - University of Manitoba/eCN_project/Datasets/Vladoiu/P0")

# Load the sparse matrix
expression_matrix <- readMM("GSM3318004_P0_matrix.mtx")

# Load gene names (features)
features <- read.table("GSM3318004_P0_genes.tsv", 
                       header = FALSE, stringsAsFactors = FALSE)

# Load cell barcodes - cell identifiers
##### If your barcode folder contains multiple files (e.g., 7383 items) instead of a single .tsv or .txt file, 
##### it likely means that each file represents one barcode (cell), with the barcode itself being the filename.
##### 
##### Instead of reading the contents of the files (which are likely empty), you only need to extract the filenames 
##### to get the barcodes.
##### 
##### Solution: Extract Barcodes from Filenames
##### Use list.files() to get all filenames in the folder.
##### This will give you a vector of barcodes, which can be assigned as column names in your expression matrix.

barcode_folder <- "GSM3318004_P0_barcodes/"

# List all file names (without reading the files)
barcodes <- list.files(barcode_folder)

# Print the first few barcodes
head(barcodes)

###Assign barcodes to Expression Matrix
colnames(expression_matrix) <- barcodes  #Assign barcodes as column names

###Now ensure Data is Correctly Formatted. Before creating the Seurat object, check:
dim(expression_matrix)  # Should be genes x cells: [1] 27998  4809
length(barcodes)        # Should match the number of columns in expression_matrix: [1] 4809
nrow(features)          # Should match the number of rows in expression_matrix: [1] 27998


# Ensure gene names are unique
features$V2 <- make.unique(features$V2)
#assign gene names as row names
rownames(expression_matrix) <- features$V2

#####################################################################################################################################
############-------------Step 3: Create a Seurat Object (Seurat v5 uses CreateSeuratObject)

seurat_object <- CreateSeuratObject(counts = expression_matrix, 
                                    project = "Mouse_P0", 
                                    min.cells = 3, 
                                    min.features = 200)

head(rownames(seurat_object))
##It returns: [1] "Xkr4"   "Sox17"  "Mrpl15" "Lypla1" "Tcea1"  "Rgs20", which are correct official gene symbols for mouse.

#Save the Seurat Object
saveRDS(seurat_object, file = "S_Object_Vladoiu_P0_raw_FG1.rds")

class(seurat_object)

#####################################################################################################################################
############-------------Step 4: Quality ControL. Check:

#Number of detected genes (low values indicate empty droplets)
#Percentage of mitochondrial genes (high values suggest dying cells)

#You can work on the saved raw data in step 3 (which is a seurat object)
seurat_object <- readRDS("S_Object_Vladoiu_P0_raw_FG1.rds")

# Calculate percentage of mitochondrial genes
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")

# Visualize quality control metrics
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#####################################################################################################################################
############-------------Step 5: Filter Low-Quality Cells. You can filter cells based on mitochondrial percentage, gene count, and RNA count.

#For example: Now Filtering to set a condition to
#Remove cells with too few or too many detected genes.
#Filter out cells with high mitochondrial content (which may indicate dying cells).

seurat_object <- subset(seurat_object, 
                        subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & 
                          nCount_RNA < 20000 & percent.mt < 5)
# Visualize post-filtering
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Check available layers
Layers(seurat_object)

#Save the seurat object
saveRDS(seurat_object, file = "S_Object_Vladoiu_P0_filtered_FG2.rds")

class(seurat_object)

#####################################################################################################################################
############-------------Step 6: Normalize the data

#You can work on the saved filtered data in step 5 (which is a seurat object)
seurat_object <- readRDS("S_Object_Vladoiu_P0_filtered_FG2.rds")

seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize")

saveRDS(seurat_object, file = "S_Object_Vladoiu_P0_normalized_FG3.rds")

Layers(seurat_object)

#####################################################################################################################################
############-------------Step 7: Subset/Isolate excitatory (glutamatergic) neurons

#You can work on the saved normalized data in step 6 (which is a seurat object)
seurat_object <- readRDS("S_Object_Vladoiu_P0_normalized_FG3.rds")

# Subset glutamatergic neurons based on specific markers for glutamatergic specifications
eCN_candidates <- subset(seurat_object, subset = (Slc17a6 > 0 | Atoh1 > 0 | Tcf4 > 0))

# Isolate eCN (excitatory cerebellar nuclei neurons) from eCN_candidates based on specific markers for eCN and other excitatory neurons (GCs and UBCs)
eCN_pure <- subset(eCN_candidates,  
                   (Meis2 > 0.5 | Lhx9 > 0.5 | Tbr1 > 0.5 | Olig2 > 0.5
                    | Lmx1a > 0.5 | Nhlh2 > 0.5 | Sox11 > 0.5) &  
                     Eomes == 0 & Neurod1 == 0 & Pax6 == 0)
#specific markers for eCN or expressed in subset of eCN: Meis2, Lhx9, Tbr1, Olig2, Lmx1a, Nhlh2, Sox11
#I excluded granule cells (GCs) based on their specific markers: Neurod1 and Pax6
#I excluded unipolar brush cells (GCs) based on their specific marker: Eomes

dim(eCN_pure)  # Returned: [1] 17418  866

#Save subsetted  eCN_pure as P0_eCN
saveRDS(seurat_object, file = "S_Object_Vladoiu_P0_eCN_FG4.rds")

#####################################################################################################################################
############-------------Step 8: Find variable features and Scale the data

#You can work on the saved normalized data (which is a seurat object)
eCN_pure <- readRDS("S_Object_Vladoiu_P0_eCN_FG4.rds")

#Find variable features (so PCA can work):
eCN_pure <- FindVariableFeatures(eCN_pure, selection.method = "vst", nfeatures = 2000)

head(VariableFeatures(eCN_pure))

# Plot variable features
plot1 <- VariableFeaturePlot(eCN_pure)
# Label the top variable features on the plot (e.g., top 10 by variance)
top10 <- head(VariableFeatures(eCN_pure), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#Scale the data
eCN_pure <- ScaleData(eCN_pure, assay = "RNA")

Layers(eCN_pure)

#####################################################################################################################################
############-------------Step 9: Run PCA

eCN_pure <- RunPCA(eCN_pure, assay = "RNA")

#####################################################################################################################################
############-------------Step 10: Run UMAP or tSNE (one of them is needed for paper)

# Run tSNE after performing PCA:
eCN_pure <- RunTSNE(eCN_pure, reduction = "pca", dims = 1:10)

#Run UMAP after performing PCA:
eCN_pure <- RunUMAP(eCN_pure, reduction = "pca", dims = 1:10)

#####################################################################################################################################
############-------------Step 11: Clustering

#FindNeighbors(): It computed the nearest neighbor graph and shared nearest neighbor (SNN) for clustering.
eCN_pure <- FindNeighbors(eCN_pure, dims = 1:30)

eCN_pure <- FindClusters(eCN_pure, resolution = 0.25)  #Higher resolution: More clusters - finer splitting

#####################################################################################################################################
############-------------Step 12: DEA (Differential Expression Analysis) and Marker Selection

#Determine Distinguishing Markers for Each Cluster

cluster_markers <- FindAllMarkers(eCN_pure, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the clustered seurat object 
saveRDS(eCN_pure, file = "S_Object_Vladoiu_P0_eCN_clustered_FG5.rds")
# Read the clustered seurat object
eCN_pure <- readRDS("S_Object_Vladoiu_P0_eCN_clustered_FG5.rds")

# Print the top markers for each cluster - Marker Selection
top_5_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top_10_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top_20_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top_50_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

#####################################################################################################################################
############-------------Step 13: Cluster annotation

#Annotate the Clusters to generate annotated tSNE or umap

new.cluster.ids <- c("eCN1", "eCN2", "eCN3", "eCN4", "eCN5", "eCN6", "eCN7", "eCN8", "eCN9", "eCN10", "eCN11", "eCN12", "eCN13")
names(new.cluster.ids) <- levels(eCN_pure)
eCN_pure <- RenameIdents(eCN_pure, new.cluster.ids)

DimPlot(eCN_pure, reduction = "tsne", label = TRUE, label.size = 3) +
  labs(color = "Cluster")

#####################################################################################################################################
############-------------Step 14: Save top markers

# Create a list of data frames with sheet names
sheets5 <- list("All Markers" = cluster_markers, "Top Markers" = top_5_markers)
sheets10 <- list("All Markers" = cluster_markers, "Top Markers" = top_10_markers)
sheets20 <- list("All Markers" = cluster_markers, "Top Markers" = top_20_markers)
sheets50 <- list("All Markers" = cluster_markers, "Top Markers" = top_50_markers)

# Write the list to a single Excel file
write.xlsx(sheets5, file = "P0_markers5.xlsx")
write.xlsx(sheets10, file = "P0_markers10.xlsx")
write.xlsx(sheets20, file = "P0_markers20.xlsx")
write.xlsx(sheets50, file = "P0_markers50.xlsx")

#####################################################################################################################################
############-------------Step 15: Plot marker genes or DEG

#To create plot
# Select top 5 or markers for each cluster
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Rename cluster identities
names(new.cluster.ids) <- levels(eCN_pure)  # Assign labels to existing cluster levels
eCN_pure <- RenameIdents(eCN_pure, new.cluster.ids)

# Option 1: DotPlot using top 5 marker genes: Scales better than heatmap and Still shows both expression level + cell proportion
# I did this
DotPlot(eCN_pure, features = unique(top_5_markers$gene)) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1))


#Option 2: Heatmap
DoHeatmap(
  eCN_pure,
  features = top_markers$gene,  # Use the top marker genes
  group.by = "ident", # Group cells by clusters
  label = TRUE                  # Add cluster labels
) + 
  theme(
    axis.text.y = element_text(size = 11),  # Adjust y-axis text size for readability
    FontSize(x.text = 10, y.text = 8),
    axis.text.x = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  ) + 
  guides(color = guide_colorbar())  # Keep only the colorbar for expression



