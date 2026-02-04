################################################################################
# Project: Integrative Single-Cell Transcriptomics of COVID-19 Severity
# Title: Discovery of S100A8/S100A9/ITGB2 Signaling Trio
# Author: Md. Shakil Mahmud Supto
# Data Source: NCBI GEO (GSE150728)
# Description: QC, SCTransform v2, Harmony Integration, and Marker Discovery
################################################################################

# 1. Load Libraries
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(magrittr)

# 2. Data Loading & Rigorous QC
# Initial cohort: ~1.5 Lakh cells across 13 donors
# Filtering: Genes per cell (200-4000), Mitochondrial content (<10%)
data.filtered <- subset(merged_data, subset = nFeature_RNA > 200 & 
                          nFeature_RNA < 4000 & 
                          percent.mt < 10)

# Result: ~1.36 Lakh high-quality cells retained

# 3. Computational Optimization (Downsampling)
# Selecting 50,000 representative cells to ensure balanced analysis
# (Approx. 29,648 COVID-19 / 20,352 Healthy)
set.seed(123)
data_50k <- data.filtered[, sample(colnames(data.filtered), 50000)]

# 4. Normalization & Variance Stabilization
# Using SCTransform v2 for superior batch effect handling
data_50k <- SCTransform(data_50k, 
                        method = "glmGamPoi", 
                        vars.to.regress = "percent.mt", 
                        verbose = FALSE)

# 5. Dimensionality Reduction (PCA)
data_50k <- RunPCA(data_50k, verbose = FALSE)

# 6. Multi-Donor Integration (Harmony)
# Correcting batch effects across 13 independent donors
data_50k <- RunHarmony(data_50k, group.by.vars = "donor", plot_convergence = TRUE)

# 7. Unsupervised Clustering & UMAP
# Using Harmony embeddings for downstream manifold learning
data_50k <- RunUMAP(data_50k, reduction = "harmony", dims = 1:30)
data_50k <- FindNeighbors(data_50k, reduction = "harmony", dims = 1:30)
data_50k <- FindClusters(data_50k, resolution = 0.8) # Identified 23 discrete clusters

# 8. Differential Expression & Marker Discovery
# Identifying top conserved markers for each cluster
all_markers <- FindAllMarkers(data_50k, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top3_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

# 9. Visualization of Clinical Biomarkers (The Trio)
# Checking the distribution of S100A8, S100A9, and ITGB2
VlnPlot(data_50k, features = c("S100A8", "S100A9", "ITGB2"), pt.size = 0, ncol = 3)
FeaturePlot(data_50k, features = c("S100A8", "S100A9", "ITGB2"))

# 10. Data Export for Python Validation
# The processed object is now ready for LIANA (Cell-Cell Communication) 
# and External Validation in Python/Scanpy.
saveRDS(data_50k, file = "COVID19_Discovery_Integrated.rds")

print("Discovery Phase Pipeline Completed Successfully.")