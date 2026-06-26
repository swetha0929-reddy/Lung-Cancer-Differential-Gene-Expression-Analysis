# Lung Cancer Differential Gene Expression Analysis Using Microarray Dataset GSE19804

## Project Overview

This project presents a bioinformatics workflow for identifying differentially expressed genes (DEGs) associated with lung cancer using the publicly available microarray dataset **GSE19804** from the NCBI Gene Expression Omnibus (GEO).

The analysis was performed using the **R programming language** and the **Bioconductor limma package** to compare gene expression profiles between lung cancer and normal lung tissue samples. Statistical analysis, visualization, and result generation were carried out through a reproducible computational pipeline.

This project was completed as part of the **PGDP in Bioinformatics, Genomics and Data Science**.

---

# Objectives

* Download the lung cancer dataset (GSE19804) from GEO.
* Perform preprocessing of gene expression data.
* Construct experimental design and contrast matrices.
* Identify differentially expressed genes using the limma package.
* Classify significant genes into upregulated and downregulated groups.
* Visualize results using Volcano Plot and Heatmap.
* Export reproducible analysis outputs for further biological interpretation.

---

# Dataset

**Source:** NCBI Gene Expression Omnibus (GEO)

**Dataset ID:** GSE19804

**Platform:** GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array)

**Samples:** 120 lung tissue samples

---

# Tools and Technologies

* R (Version 4.6.0)
* RStudio
* GEOquery
* limma
* ggplot2
* pheatmap
* Bioconductor

---

# Analysis Workflow

1. Download GEO dataset
2. Extract expression matrix and sample metadata
3. Prepare experimental groups
4. Create design matrix
5. Construct contrast matrix
6. Perform differential expression analysis using limma
7. Apply empirical Bayes statistics
8. Identify significant genes
9. Generate Volcano Plot
10. Generate Heatmap of Top 50 DEGs
11. Export results and summary files

---

# Results Summary

| Parameter             | Value       |
| --------------------- | ----------- |
| Dataset               | GSE19804    |
| Total genes analyzed  | 54,675      |
| Significant genes     | 1,960       |
| Upregulated genes     | 591         |
| Downregulated genes   | 1,369       |
| Statistical Method    | limma       |
| Adjusted P-value      | < 0.05      |
| Fold Change Threshold | |logFC| > 1 |

---

# Repository Structure

```
Lung-Cancer-Differential-Gene-Expression-Analysis/

├── data/
│   ├── Sample_Information.csv
│   └── README.md
│
├── figures/
│   ├── Volcano_Plot.png
│   ├── Heatmap_Top50_DEGs.png
│   ├── Workflow_Diagram.png
│   └── README.md
│
├── results/
│   ├── Differential_Expression_Results.csv
│   ├── Significant_Genes.csv
│   ├── Analysis_Summary.txt
│   └── README.md
│
├── scripts/
│   ├── dge_analysis.R
│   └── README.md
│
├── report/
│   ├── Mini_Project_Report.pdf
│   └── README.md
│
├── presentation/
│   ├── Project_Presentation.pptx
│   └── README.md
│
└── README.md
```

---

# Key Outputs

* Differential Expression Results
* Significant Genes
* Volcano Plot
* Heatmap of Top 50 Differentially Expressed Genes
* Analysis Summary
* Complete R Script

---

# Applications

The workflow developed in this project can be applied to:

* Cancer biomarker discovery
* Differential gene expression analysis
* Transcriptomics research
* Precision medicine
* Computational biology
* Bioinformatics research

---

# Future Improvements

* Gene Ontology (GO) Enrichment Analysis
* KEGG Pathway Analysis
* Protein-Protein Interaction (PPI) Network Analysis
* Survival Analysis
* Validation using independent datasets
* RNA-Seq based differential expression analysis

---

# Author

**Thamanna Swetha Reddy**

PGDP in Bioinformatics, Genomics and Data Science

Bversity School of Bioscience

---

# License

This repository is intended for educational and research purposes.
