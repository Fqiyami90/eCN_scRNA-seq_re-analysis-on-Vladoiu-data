scRNA-seq Analysis Pipeline for Snca+ Neuronal Population Profiling in the developing cerebellum

Description

This repository contains an R-based analysis workflow for single-cell RNA sequencing (scRNA-seq) data to profile neuronal populations, including:

Quality control and filtering

Normalization and scaling

Dimensionality reduction (PCA, UMAP)

Clustering and cell-type annotation

Differential gene expression analysis

Functional enrichment (GO and GSEA)


Tools Used

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
