# Interactive and Modular RNA-Seq Analysis Pipeline
## Overview
This Repository contains a comprehensive, interactive R script designed for an **end-to-end RNA Sequencing data analysis workflow.** It allows users to manage data processing, normalization, quality control, Differential Gene Expression analysis, advanced gene set enrichment analysis and Pathway enrichment analysis through a command-line interface.

## Pipeline Features
| Module | Letter | Analysis Step | Key Functions |
| :--- | :--- | :--- | :--- |
| Data Import | A | Reads RSEM `.genes.results`files, merges them and generates unified **Count, TPM and FPKM** expression matrices. | `build_expression_matrices`|
| Normalization | B | Performs **Variance Stabilizing Transformation** and offers optional **Quantile Normalization** on the VST data. | `normalize_vst`, `normalize_qn`|
| Quality Control | C | Generates essential QC plots for selected expression data (TPM, FPKM, VST): **Density Plot, Sample-to-Sample Correlation Heatmap and Principal Component Analysis.** | `plot_density`, `correlation_plot`, `plot_pca`|
| DGE Analysis | D | Run DESeq2 analysis for user-defined comparisons. Outputs are raw DESeq2 results, filtered Differentially Expressed Genes (given a p-value treshold) and high-resolution Volcano Plots. | `differentially_expressed_genes`|
| Pathway Analysis | E | Modular pathway analysis: E1 - filter DEGs by KEGG pathway; E2 - visualize expression (Heatmap, PCA) for genes in a specified pathway; E3 - Over-Representation Analysis (ORA, KEGG, GO); E4 - Heatmaps of Log2FC for pathway DEGs; E5 - Gene Set Enrichment Analysis (GSEA using MSigDB or custom GMT files). | `filter_degs_by_pathways`, `enrichment_analysis_ora`, `gsea_analysis`|
| Gene Lists | F | Generation of shared and unique gene lists. | `unique_genes_list`, `shared_genes_list`|

---

## 🛠️ Installation and Dependencies

This script relies on several popular Bioconductor and CRAN packages. You must have R and RStudio installed.

### 📦 R Packages

Install the necessary packages in your R environment:

```R
# Install Bioconductor packages if not already present
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "preprocessCore", "EnhancedVolcano"))

# Install CRAN packages
install.packages(c("dplyr", "ggplot2", "tidyr", "KEGGREST", "corrplot", "tidyverse",
                   "pheatmap", "RColorBrewer", "reshape", "forcats"))
```
### 🧰 Required R Packages

```r
library(dplyr)
library(ggplot2)
library(tidyr)
library(KEGGREST)
library(corrplot)
library(DESeq2)
library(EnhancedVolcano)
library(preprocessCore)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(reshape)
library(forcats)
```
## 📂 Input Data Requirements
Depending on the module the script will prompt for the FULL path to the following required files:
1. Project Name - A single, descriptive name for the output folder that shouldn't include white spaces or dots.
2. RSEM files folder - The directory containing all individual `*.genes.results` files. E.g. moonG.genes.results
3. Sample Metadata Table - A CSV file with columns including:
   - First Column: Sample IDs (must match column names in expression matrices. Same order!)
   - Group: Primary factor for DESeq2 design (e.g. replicates in groups).
   - Condition: Additional factor for PCA plot grouping (e.g. control, treatment)
4. Ensembl Mapping Table - A csv file with columns: gene_id (Ensembl ID), and Gene.name (Gene Symbol)
5. DESeq2 Comparison Table (*Module D*) - A csv file defining contrasts with columns: `Factor` (should be Group), `Normal` (Reference Level), `Treatment` (Contrast Level)
6. Pathway IDs Table (*Module E*) - A csv file with KEGG IDs for filtering and visualization.
7. GMT File (*Module E5*) - FULL path to a downloaded MSigDB (c2, c5, etc.) or custom GMT file.


## 📝 Usage

**1. Save the Script:** Save the provided R code as an R file (e.g., `analysis_pipeline.R`)

**2. Open the Terminal and start R.**

**3. Run the Main Function:** Launch your R environment and execute the following:
```R
source("analysis_pipeline.R")
run_analysis()
```
**4. Follow the Prompts:** The script will enter an interactive mode, asking for file paths, the letters of the modules to run (e.g. ABCDE) and the base matrix type (TPM, FPKM, VST) for QC and visualization.
**5. Data Loading Flexibility:** If you skip modules A (Import) and B (Normalization), the script will automatically check the project output folder for previously generated files (e.g., `[ProjectName]_vst_names.csv` or prompt you to load external pre-calculated expression matrices.

## 📊 Output Structure

All results are automatically written to a structured output folder organized by project name.  
Default output path:  

```R
~/Desktop/Analysis_RNA_Seq/
└── [Project_Name]/
    ├── count_files/
    │   ├── [Project_Name]_count_genenames.csv
    │   ├── [Project_Name]_TPM_genenames.csv
    │   ├── [Project_Name]_vst_names.csv
    │   └── [Project_Name]_qn_vst_names.csv
    ├── Density_Plot_[DataType].png
    ├── Correlation_Heatmap.png
    ├── [GeneSet]_PCA.png
    ├── results_[Contrast].csv      (Full DESeq2 output)
    ├── DEGs_[Contrast].csv         (Filtered significant DEGs)
    ├── VolcanoPlot_[Contrast].png
    ├── E1_FilteredDEGs_[Contrast].csv
    ├── E3_KEGG_Enrichment_[Contrast].csv
    └── E5_GSEA_Results_[ProjectName].csv
```

Each major module in the script saves its outputs automatically:

| Analysis step | Output file(s) | Description |
|----------------|----------------|--------------|
| **RSEM import & merge** | `_count.csv`, `_TPM.csv`, `_FPKM.csv` | Combined count and expression matrices |
| **Normalization (VST/QN)** | `vst_names.csv`, `qn_vst_names.csv` | Normalized expression data |
| **Correlation & PCA** | `.pdf` or `.png` (plots) | Visual representation of sample relationships |
| **Unique/shared genes** | `<group>_genes_unique.csv`, `genes_shared.csv` | Lists of group-specific or shared genes |
| **Differential expression** | `DEGs_<contrast>.csv` | Results of DESeq2 comparisons with adjusted p-values |
| **Pathway analysis (optional)** | `.csv` tables or plots | Gene sets per KEGG pathway |

All outputs are timestamped and automatically organized by project for reproducibility.
