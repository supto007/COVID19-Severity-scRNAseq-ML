# Integrative Single-Cell Transcriptomics and Machine Learning Identify an S100A8/S100A9/ITGB2 Signaling Trio Across Neutrophil and Plasma Cell Lineages as Predictive Biomarkers of COVID-19 Severity

**Author:** **Md. Shakil Mahmud Supto** **Role:** DVM & Student
**Affiliation:** Sylhet Agricultural University, Sylhet

# Data Source & Study Cohort
**Internal Discovery Cohort (Training/Discovery)**
Retrieved from the NCBI Gene Expression Omnibus (GEO) accession number: **GSE150728**
Original Publication: Wilk et al. (2020), "A single-cell atlas of the peripheral immune response in severe COVID-19".
**Study Population:** A total of 13 donor samples were integrated.
**COVID-19 Patients:** 7 individuals 
**Healthy Controls:** 6 individuals.

**External Validation Details**
Source: CellXGene Census (Genome sequence Archive for Human: HRA001149).
Data Scale: 318,894 cells (Independent Meta-cohort).
Target Population: High-resolution Neutrophil and Plasma Cell lineages.
Clinical Groups: COVID-19 vs. Healthy Controls.
Significance: Validates the signature's performance across diverse platforms and global patient populations.

## Project Overview
This repository contains the bioinformatic pipeline and machine learning approach used to identify and validate a novel 3-gene signaling trio (S100A8, S100A9, ITGB2) as predictive biomarkers for COVID-19 severity. By integrating high-resolution single-cell RNA sequencing (scRNA-seq) data from both internal discovery and massive external validation cohorts, we demonstrate the robust diagnostic potential of these neutrophil and plasma cell-driven markers.

## Key Highlights
**Discovery Dataset:** Analyzed 136,236 cells (Wilk et al., 2020) after rigorous QC to identify severity-linked transcriptomic shifts.
**Signaling Trio:** Identified S100A8, S100A9, and ITGB2 as core biomarkers through Differential Gene Expression (DGE) and Ligand-Receptor interaction analysis (LIANA).
**Massive External Validation:** Validated the signature on an independent cohort of 318,894 neutrophils retrieved from the CellXGene database.
**Predictive Performance:** The signature achieved a statistically significant AUC of 0.66 (Internal) and AUC of 0.65 (External) with a P-value of 0.00e+00, confirming its generalizability across different platforms and patient populations.

## Bioinformatic & ML Pipeline
### 1. Discovery Phase (R / Seurat)
**Preprocessing:** Filtered cells based on mitochondrial content (<10%) and feature counts (200-4000).
**Normalization & Integration:** Used SCTransform v2 for variance stabilization and Harmony for multi-donor batch effect removal.
**Clustering:** Identified 23 distinct immune clusters representing the peripheral immune landscape.
**Discovery:** Identified the S100A8/S100A9/ITGB2 trio as key drivers of COVID-19 severity in myeloid and lymphoid lineages.

### 2. Validation Phase (Python / Scanpy)
**Stacking Logic:** Developed a Composite Signature Scoring model by averaging the Z-normalized expression of the three target genes.
**Statistical Analysis:** Applied the Mann-Whitney U test to confirm significant distribution shifts between COVID-19 and Healthy groups.
**Performance Metrics:** Generated ROC-AUC curves to evaluate the classification power of the signature on 300k+ unseen cells.

## 📊 Visual Results
The results demonstrate a clear molecular shift in COVID-19 patients. While individual cell variability exists, the composite score of the $S100A8/S100A9/ITGB2$ trio provides a statistically robust signal for disease identification. The near-identical performance in both internal (0.66) and external (0.65) cohorts highlights the stability of this biomarker trio.

## Tech Stack
* **Language:** R, Python
* **R Environment:** Seurat v4/v5, Harmony, LIANA, ggplot2, ComplexHeatmap
* **Python Environment:** Scanpy, Pandas, Scikit-learn, Seaborn, Matplotlib
