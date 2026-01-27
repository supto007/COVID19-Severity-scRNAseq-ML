################################################################################
# Project: Integrative Single-Cell Transcriptomics of COVID-19 Severity
# Author: Md. Shakil Mahmud Supto
# Data Source: NCBI GEO (GSE150728)
# Description: QC, SCTransform, Harmony Integration, and Clustering of PBMCs
################################################################################

# 1. Load Libraries
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)

# 2. Data Loading & Initial QC (1.5 Lakh Cells)
# Assuming 'merged_data' contains 13 donors
data.filtered <- subset(merged_data, subset = nFeature_RNA > 200 & 
                          nFeature_RNA < 4000 & 
                          percent.mt < 10)
# Result: ~1.36 Lakh cells retained

# 3. Subsetting for Computational Efficiency
# Selecting 50,000 representative cells (29,648 COVID / 20,352 Healthy)
set.seed(123)
data_50k <- data.filtered[, sample(colnames(data.filtered), 50000)]

# 4. Normalization and Variance Stabilization (SCTransform v2)
data_50k <- SCTransform(data_50k, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

# 5. Dimensionality Reduction (PCA)
data_50k <- RunPCA(data_50k, verbose = FALSE)

# 6. Batch Correction (Harmony Integration across 13 Donors)
data_50k <- RunHarmony(data_50k, group.by.vars = "donor", plot_convergence = TRUE)

# 7. Clustering and UMAP Visualization
data_50k <- RunUMAP(data_50k, reduction = "harmony", dims = 1:30)
data_50k <- FindNeighbors(data_50k, reduction = "harmony", dims = 1:30)
data_50k <- FindClusters(data_50k, resolution = 0.8) # Identified 23 clusters

# 8. Marker Identification
all_markers <- FindAllMarkers(data_50k, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top3_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

# 9. Visualization
# UMAP Plot
DimPlot(data_50k, reduction = "umap", label = TRUE) + NoLegend()

# Heatmap for Top 3 Genes per Cluster
DoHeatmap(data_50k, features = top3_markers$gene) + NoLegend()

# DotPlot for Canonical Markers
# Replace 'features' with your specific marker list
DotPlot(data_50k, features = c("CD3D", "CD14", "MS4A1", "FCGR3A", "GNLY")) + RotatedAxis()

# 10. Status: WORK IN PROGRESS
# Upcoming: Differential Gene Expression (DGE) between COVID-19 vs Healthy cohorts.