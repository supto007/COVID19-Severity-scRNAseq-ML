# Integrative Single-Cell Transcriptomics and Machine Learning Identify Neutrophil-Driven Signaling Signatures as Predictive Biomarkers of COVID-19 Severity

**Author:** **Md. Shakil Mahmud Supto** **Role:** DVM & Student
**Affiliation:** Sylhet Agricultural University, Sylhet

# Data Source & Study Cohort
 Retrieved from the NCBI Gene Expression Omnibus (GEO) accession number: **GSE150728**
Original Publication: Wilk et al. (2020), "A single-cell atlas of the peripheral immune response in severe COVID-19".
**Study Population:** A total of 13 donor samples were integrated.
**COVID-19 Patients:** 7 individuals 
**Healthy Controls:** 6 individuals.

### 🧬 Key Highlights
* **Data Scale:** Analyzed **136,236 cells** after rigorous QC from 13 donors.
* **Balanced Subset:** Downsampled to a representative **50,000 cells** (29,648 COVID-19 / 20,352 Healthy) for high-fidelity scaling.
* **Advanced Pipeline:** Employed **SCTransform v2** for variance stabilization and **Harmony** for multi-donor batch correction.
* **Clustering:** Identified **23 discrete cellular clusters** representing the diverse immune landscape.

## 🛠️ Bioinformatic Pipeline
1. **Preprocessing & QC:** Filtering cells based on mitochondrial content (<10%) and feature counts (200-4000).
2. **Normalization:** SCTransform (2,000 HVGs).
3. **Integration:** Harmony-based batch effect removal across 13 donors.
4. **Dimensionality Reduction:** PCA and UMAP visualization.
5. **Cluster Analysis:** Identification of 23 clusters and top marker discovery (Differential Gene Expression).
6. **Visualization:** Cluster-specific heatmaps (Top 3 genes) and canonical marker DotPlots.

## 📊 Results Summary
The analysis revealed significant transcriptomic shifts in neutrophil populations and other myeloid lineages. The predictive signatures identified here serve as potential biomarkers for clinical severity assessment in COVID-19 patients.

## 💻 Tech Stack
* **Language:** R
* **Main Package:** Seurat v4/v5
* **Integration:** Harmony
* **Visualization:** ggplot2, ComplexHeatmap
