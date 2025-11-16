# Modular-RNA-Sequencing-Analysis
A modular R-based RNA-Seq analysis pipeline. Starting from a count file or rsem output, the pipeline can be used for merging RSEM outputs, normalization (VST/QN), visualization (PCA, heatmaps, boxplots), differential expression analysis (DESeq2)and pathway and gene-level exploration.

# 🧬 RNA-Seq Analysis Pipeline

## Project Overview

This R script provides a **complete, end-to-end pipeline for RNA-Sequencing data analysis**. It handles raw RSEM output files (or a combined count file), performs data normalization (VST and Quantile), conducts quality control (PCA and correlation heatmaps), identifies uniquely and commonly expressed genes, and performs **Differential Gene Expression (DGE) analysis using DESeq2**. It also includes functions for visualizing results and preparing data for pathway analysis (KEGG).

## 🚀 Key Features

* **Data Import & Quantification:** Reads in RSEM `.genes.results` files and generates unified count, TPM, and FPKM matrices.
* **Normalization:** Implements **VST (Variance Stabilizing Transformation)** via DESeq2 and optional **Quantile Normalization (QN)**.
* **Quality Control (QC):**
    * Generates a **Sample-to-Sample Correlation Heatmap** using VST data.
    * Performs **Principal Component Analysis (PCA)** for overall sample clustering and pathway-specific gene sets.
    * Generates **Gene Expression Distribution Boxplots** to assess sample similarity.
* **Co-expression Analysis:** Identifies **uniquely expressed** and **shared genes** across different experimental groups.
* **Differential Gene Expression (DGE):** Runs **DESeq2** for robust statistical comparison between user-defined conditions.
* **Pathway Tools:** Includes initial functions for fetching **KEGG Pathway gene lists** for downstream analysis.

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
## 🧰 Required R Packages

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

## 📝 Usage

### 1. Configuration (Command Line Arguments)

The script is designed to be highly configurable via a set of initial variables. **Before running the script**, you must set the following parameters at the beginning of the `Analysis_script.R` file:

| Variable | Description | Example |
| :--- | :--- | :--- |
| `folderPath_count_files` | Path to the directory containing RSEM `.genes.results` files. | `/path/to/rsem/output/` |
| `project_name` | A name for your project, used for output folder creation. | `My_Treatment_Study` |
| `sample_table` | A data frame containing sample metadata (e.g., sample ID, `group`, `condition`). | `read.csv("sample_metadata.csv")` |
| `list_of_comparisons` | A list defining the pairs for DESeq2 DGE analysis. | `list(c("group", "normal", "treatment"))` |
| `Pathways_of_interest` | A list of pathway names or a file path to a list of pathways for subsetting. | `c("Apoptosis", "p53 Signaling")` |

---

### 2. Running the Pipeline Steps

The following steps outline the execution order for the main parts of the analysis:

#### A. Data Processing and Matrix Generation

This step reads your raw RSEM files and generates count, TPM, and FPKM matrices.

```R
--- RNA-Seq Analysis Pipeline ---
Enter the FULL path to the folder containing RSEM (.genes.results) files:
Enter a Project Name without white spaces: 
Enter the FULL path to the Sample Metadata Table (.csv file with columns: sample, group, condition): 
Enter the FULL path to the Ensembl ID to Gene Symbol mapping table (.csv with 'gene_id' and 'Gene_name'):                                                                                                                   

Available Analysis Options:
 A: Data Import & Matrix Creation (Count, TPM, FPKM)
 B: Normalization (VST, optional Quantile)
 C: Quality Control Plots (Correlation Plot, Density, PCA)
 D: Differential Gene Expression (DESeq2)
 E: Gene Set/Pathway Analysis
   E1: Filter DEGs per pathway
   E2: Heatmap per Pathway
   E3: Enrichment Analysis & Over-Representation Analysis
   E4: DEGs Heatmaps
   E5: GSEA (Gene Set Enrichment Analysis)
 F: Shared/Unique Gene Lists

Enter the letters for analyses to run (e.g., ABCDEF):
Choose expression matrix for downstream analysis (TPM/FPKM/VST):


```


## 📊 Output Structure

All results are automatically written to a structured output folder organized by project name.  
Default output path:  

```R
~/Desktop/Analysis_RNA_Seq/
├── <project_name>/
│   ├── count_files/
│   │   ├── _TPM.csv
│   │   ├── _FPKM.csv
│   │   ├── _count.csv
│   │   └── vst_names.csv (Normalized VST data)
│   ├── Differential_Gene_Expression_Analysis/
│   │   ├── results_<comp>.csv (Full DESeq2 output for each comparison)
│   │   └── DEGs_<comp>_with_names.csv (Filtered list of significant DEGs)
│   └── Coexpression_Analysis/
│       ├── JR_genes_shared.csv (Genes expressed in all groups)
│       └── JR_<group>_genes_unique.csv (Genes unique to each group)
```

Each major module in the script saves its outputs automatically:

| Analysis step | Output file(s) | Description |
|----------------|----------------|--------------|
| **RSEM import & merge** | `_count.csv`, `_TPM.csv`, `_FPKM.csv` | Combined count and expression matrices |
| **Normalization (VST/QN)** | `vst_names.csv`, `qn_vst_names.csv` | Normalized expression data |
| **Correlation & PCA** | `.pdf` or `.png` (plots) | Visual representation of sample relationships |
| **Unique/shared genes** | `JR_<group>_genes_unique.csv`, `JR_genes_shared.csv` | Lists of group-specific or shared genes |
| **Differential expression** | `DEGs_<contrast>.csv` | Results of DESeq2 comparisons with adjusted p-values |
| **Pathway analysis (optional)** | `.csv` tables or plots | Gene sets per KEGG pathway |

All outputs are timestamped and automatically organized by project for reproducibility.
