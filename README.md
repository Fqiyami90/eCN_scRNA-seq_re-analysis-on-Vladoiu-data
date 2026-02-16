single-cell RNA sequencing (scRNA-seq) Analysis Pipeline for Snca+ Neuronal Population Profiling in the developing cerebellum

Description

This repository contains an R-based analysis workflow for scRNA-seq data from following published article to profile neuronal populations:

Childhood cerebellar tumours mirror conserved fetal transcriptional programs 2019
Maria C. Vladoiu, Ibrahim El-Hamamy, Laura K. Donovan, Hamza Farooq, Borja L. Holgado, Yogi Sundaravadanam, Vijay Ramaswamy, Liam D. Hendrikse, Sachin Kumar, Stephen C. Mack, John J. Y. Lee, Vernon Fong, Kyle Juraschka, David Przelicki, Antony Michealraj, Patryk Skowron, Betty Luu, Hiromichi Suzuki, A. Sorana Morrissy, Florence M. G. Cavalli, Livia Garzia, Craig Daniels, Xiaochong Wu, Maleeha A. Qazi, …Michael D. Taylor

The pipline includes:

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
