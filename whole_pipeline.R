## ----------------------------------------------------------------------------
## 1. Library Setup and helper functions
## ----------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(preprocessCore)
library(pheatmap)
library(RColorBrewer)
library(KEGGREST)
library(corrplot)
library(matrixStats)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(DOSE)
library(data.table)
library(readr)

# Helper function for user input
get_user_input <- function(prompt) {
  return(readline(prompt = prompt))
}

## ----------------------------------------------------------------------------
## 2. Core data processing functions
## ----------------------------------------------------------------------------

build_expression_matrices <- function(files, output_dir, ensembl_table, project_name) {
  cat("Reading RSEM files and building matrices (Aggressive memory cleanup enabled)...\n")
  sample_names <- sub("\\.genes\\.results$", "", basename(files))
  
  # 1. Initialize Matrices using the first file
  first_file <- files[1]
  
  # Use data.table::fread for speed
  count_dt <- data.table::fread(first_file, select = c("gene_id", "expected_count"))
  TPM_dt <- data.table::fread(first_file, select = c("gene_id", "TPM"))
  FPKM_dt <- data.table::fread(first_file, select = c("gene_id", "FPKM"))
  
  # Rename columns to sample name
  data.table::setnames(count_dt, "expected_count", basename(first_file))
  data.table::setnames(TPM_dt, "TPM", basename(first_file))
  data.table::setnames(FPKM_dt, "FPKM", basename(first_file))
  
  # 2. Iterative, Aggressive Merge (Samples 2 through N)
  for (i in 2:length(files)) {
    file <- files[i]
    sample <- basename(file)
    cat(paste("Merging file", i, "/", length(files), ":", sample, "\n"))
    
    # Read data for the current file
    current_dt <- data.table::fread(file, select = c("gene_id", "expected_count", "TPM", "FPKM"))
    
    # --- COUNT Merge ---
    count_merge_dt <- current_dt[, c("gene_id", "expected_count"), with = FALSE]
    data.table::setnames(count_merge_dt, "expected_count", sample)
    count_dt <- merge(count_dt, count_merge_dt, by = "gene_id", all = TRUE)
    rm(count_merge_dt); gc() # Aggressive cleanup
    
    # --- TPM Merge ---
    TPM_merge_dt <- current_dt[, c("gene_id", "TPM"), with = FALSE]
    data.table::setnames(TPM_merge_dt, "TPM", sample)
    TPM_dt <- merge(TPM_dt, TPM_merge_dt, by = "gene_id", all = TRUE)
    rm(TPM_merge_dt); gc() 
    
    # --- FPKM Merge ---
    FPKM_merge_dt <- current_dt[, c("gene_id", "FPKM"), with = FALSE]
    data.table::setnames(FPKM_merge_dt, "FPKM", sample)
    FPKM_dt <- merge(FPKM_dt, FPKM_merge_dt, by = "gene_id", all = TRUE)
    rm(FPKM_merge_dt); gc() 
    
    # Explicitly remove the current data chunk after its use
    rm(current_dt); gc() 
  }
  
  cat("Merging complete. Finalizing matrices...\n")
  
  # 3. Final Conversion and Cleanup
  finalize_matrix <- function(dt, type) {
    # Convert to tibble/data.frame for compatibility with tidyverse
    df <- as_tibble(dt) 
    colnames(df)[-1] <- sample_names # Ensure correct sample names
    df[is.na(df)] <- 0
    df$gene_id <- sub("\\..*", "", df$gene_id) # Remove version numbers
    
    # Merge with gene symbols (assuming ensembl_table is small enough)
    df_merged <- df %>%
      left_join(ensembl_table, by = "gene_id") %>%
      relocate(Gene.name)
    
    return(df_merged)
  }
  
  count_df_merged <- finalize_matrix(count_dt, "count")
  rm(count_dt); gc()
  
  TPM_df_merged <- finalize_matrix(TPM_dt, "TPM")
  rm(TPM_dt); gc()
  
  FPKM_df_merged <- finalize_matrix(FPKM_dt, "FPKM")
  rm(FPKM_dt); gc()
  
  # --- Writing Files --- (The writing function remains the same as before)
  count_path <- file.path(output_dir, "count_files")
  dir.create(count_path, recursive = TRUE, showWarnings = FALSE)
  
  # Function to write (using readr::write_csv for speed)
  readr::write_csv(count_df_merged, file = file.path(count_path, paste0(project_name, "_count_genenames.csv")))
  readr::write_csv(TPM_df_merged, file = file.path(count_path, paste0(project_name, "_TPM_genenames.csv")))
  readr::write_csv(FPKM_df_merged, file = file.path(count_path, paste0(project_name, "_FPKM_genenames.csv")))
  
  cat("Expression matrices created and saved.\n")
  return(list(count_df = count_df_merged, TPM_df = TPM_df_merged, FPKM_df = FPKM_df_merged))
}

normalize_vst <- function(sample_table, count_file, output_dir, project_name) {
  cat("Performing VST Normalization. \n")
  
  #Prepare Count Matrix
  count_matrix <- count_file %>%
    filter(!duplicated(Gene.name)) %>%
    filter(!is.na(Gene.name)) %>%
    column_to_rownames("Gene.name") %>%
    dplyr::select(-gene_id) %>%
    filter(rowSums(.)>=10) %>%
    round()%>%
    as.matrix()
  
  #Prepare ColData
  sample_table_filtered <- sample_table %>%
    filter(sample_id %in% colnames(count_matrix)) %>%
    column_to_rownames("sample_id")
  
  #Create DESeq2 object and VST
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_table_filtered,
    design = ~Group
  )
  
vst_data <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_data)

#Convert back to data frame for consistent saving
vst_df <- as.data.frame(vst_matrix) %>%
  rownames_to_column("Gene.name")

#Safe to file
path_name <- file.path(output_dir, "count_files", paste0(project_name, "_vst_names.csv"))
write_csv(vst_df, path_name)

cat("VST normalized data saved to:", path_name, ". \n")
return(vst_df)
}

normalize_qn <- function(vst_data, output_dir, project_name) {
  cat("Performing Quantile Normalization. \n")
  #Prepare matrix
  vst_matrix <- vst_data %>%
    column_to_rownames("Gene.name") %>%
    as.matrix()
  
  #Quantile Normalization
  qn_vst_matrix <- normalize.quantiles(vst_matrix)
  rownames(qn_vst_matrix) <- rownames(vst_matrix)
  colnames(qn_vst_matrix) <- colnames(vst_matrix)
  
  #Convert back to data frame for consistent saving
  qn_vst_df <- as.data.frame(qn_vst_matrix) %>%
    rownames_to_column("Gene.name")
  
  #Safe to file
  path_name <- file.path(output_dir, "count_files", paste0(project_name, "_qn_vst_names.csv"))
  write_csv(qn_vst_df, path_name)
  
  cat("Quantile VST normalized data saved to:", path_name, "\n")
  return(qn_vst_df)
}

load_analysis_data <- function(project_name, data_type_choice, count_files_dir) {
  cat(paste0("\n--- Attempting to load existing ", data_type_choice, " data ---\n"))
  analysis_df <- NULL
  
  # Determine expected file name based on user choice
  if (data_type_choice == "VST") {
    file_name <- paste0(project_name, "_vst_names.csv")
  } else {
    file_name <- paste0(project_name, "_", tolower(data_type_choice), "_genenames.csv")
  }
  
  analysis_file_path <- file.path(count_files_dir, file_name)
  
  if (file.exists(analysis_file_path)) {
    analysis_df <- read_csv(analysis_file_path, show_col_types = FALSE)
    cat(paste("Successfully loaded file from project folder:", analysis_file_path, "\n"))
  } else {
    cat(paste("Project file not found: ", analysis_file_path, "\n"))
    
    # Prompt for external file
    external_path <- get_user_input(paste0("Enter the FULL path to your pre-calculated ", data_type_choice, " file: "))
    if (file.exists(external_path)) {
      analysis_df <- read_csv(external_path, show_col_types = FALSE)
      cat(paste("Successfully loaded external file:", external_path, "\n"))
    } else {
      cat("Error: External file not found or invalid path.\n")
    }
  }
  
  # Check if Gene.name column exists, as required by C/E/F
  if (!is.null(analysis_df) && !"Gene.name" %in% colnames(analysis_df)) {
    stop("Loaded analysis file is missing the required 'Gene.name' column. Cannot proceed with C, E, or F.")
  }
  
  return(analysis_df)
}

## ----------------------------------------------------------------------------
## 3. Core Analysis & Plotting Functions
## ----------------------------------------------------------------------------

plot_density <- function(expression_matrix, project_name, output_dir, data_type) {
  cat("Generating Gene Expression Density Plot. \n")
  #Convert matrix to long format for ggplot
  df_long <- expression_matrix %>%
    as.data.frame() %>%
    rownames_to_column("Gene.name") %>%
    pivot_longer(
      cols= -Gene.name,
      names_to = "Sample",
      values_to = "Expression"
    )
  #Plot density
  plot <- ggplot(df_long, aes(x=Expression, color = Sample)) + 
    geom_density(linewidth = 1) + 
    labs(
      title = paste("Gene Expression Density Plot (", project_name, "-", data_type, ")"),
      x = paste("Expression Value (", data_type, ")"),
      y = "Density"
    ) +
    theme_minimal()+
    theme(legend.position = "none") # Hide legend for clarity if many samples
  
  #Save plot
  file_name <- file.path(output_dir, paste0("Density_Plot_", data_type, ".png"))
  ggsave(file_name, plot, width = 8, height = 5)
  cat("Density Plot saved to:", file_name, "\n")
}

correlation_plot <- function(sample_table, vst_matrix, output_dir) {
  cat("Generating Sample-to-Sample Correlation Heatmap...\n")
  # Ensure matrix is properly subsetted to matching samples
  sample_ids <- intersect(colnames(vst_matrix), sample_table$sample_id)
  vst_matrix_subset <- vst_matrix[, sample_ids]
  
  # Calculate correlation matrix
  correlation_matrix <- cor(vst_matrix_subset, method = "pearson")
  
  # Prepare annotation data frame
  annotation_df <- sample_table %>%
    filter(sample_id %in% sample_ids) %>%
    dplyr::select(sample_id, Group) %>%
    column_to_rownames("sample_id") %>%
    rename(Group = Group)
  
  # Generate Heatmap
  pheatmap(
    mat = correlation_matrix,
    annotation_col = annotation_df,
    main = "Sample-to-Sample Correlation Heatmap (VST/TPM/FPKM)",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    border_color = "gray60",
    fontsize = 8,
    display_numbers = FALSE,
    filename = file.path(output_dir, "Correlation_Heatmap.png")
  )
}


plot_pca <- function(sample_table, norm_data, gene_subset, name_gene_set, output_dir) {
  cat(paste("Generating PCA plot for:", name_gene_set, "\n"))
  
  # 1. Gene Selection
  if (is.numeric(gene_subset)) { # If gene_subset is a number (ntop)
    ntop <- gene_subset
    rowvars <- rowVars(norm_data)
    select_genes <- names(rowvars)[order(rowvars, decreasing = TRUE)][1:min(ntop, length(rowvars))]
    title <- paste0("PCA (Top ", ntop, " most variable genes)")
  } else { # If gene_subset is a vector of gene names (pathway)
    select_genes <- intersect(gene_subset, rownames(norm_data))
    if (length(select_genes) < 5) {
      cat(paste("Skipping PCA for", name_gene_set, ": too few genes (", length(select_genes), ").\n"))
      return(NULL)
    }
    title <- paste0("PCA: ", name_gene_set)
  }
  # 2. PCA Calculation
  pca_res <- prcomp(t(norm_data[select_genes, ]), scale. = FALSE)
  percentVar <- pca_res$sdev ^ 2 / sum(pca_res$sdev ^ 2)
  
  # 3. Prepare Plotting Data Frame
  plot_data <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  ) %>%
    rownames_to_column("sample_id") %>%
    left_join(sample_table, by = "sample_id")
  # 4. Plotting
  plot <- ggplot(data = plot_data, aes(
    x = PC1,
    y = PC2,
    # Assumes 'group' and 'condition' columns exist in sample_table
    colour = Group,
    shape = condition
  )) +
    geom_point(size = 3) +
    geom_text(aes(label = sample_id), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    ggtitle(title) +
    coord_fixed() +
    theme_minimal() +
    theme(plot.title = element_text(size = 12))
  # Note: Color scales would be defined here if provided
  
  # Save Plot
  file_name <- paste0(gsub(" ", "_", name_gene_set), "_PCA.png")
  ggsave(file.path(output_dir, file_name), plot, width = 6, height = 5)
  
  return(plot)
}


differentially_expressed_genes <- function(sample_table, count_file, list_of_comparisons, padj_treshold, output_dir) {
  cat("Starting DESeq2 Analysis...\n")
  # 1. Prepare DESeq2 object
  count_matrix <- count_file %>%
    filter(!duplicated(Gene.name)) %>%
    filter(!is.na(Gene.name)) %>%
    column_to_rownames("Gene.name") %>%
    dplyr::select(-gene_id) %>%
    filter(rowSums(.) >= 10) %>%
    round() %>%
    as.matrix()
  
  sample_table_filtered <- sample_table %>%
    filter(sample_id %in% colnames(count_matrix)) %>%
    column_to_rownames("sample_id")
  
  # Use the first factor in the comparison table for the design
  main_factor <- list_of_comparisons$Factor[1]
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_table_filtered,
    design = ~ 1 
  )
  # Ensure factor column is a factor and the design is updated
  colData(dds)[[main_factor]] <- factor(colData(dds)[[main_factor]])
  design(dds) <- formula(paste("~", main_factor))
  
  dds <- DESeq(dds, quiet = TRUE)
  cat("DESeq2 object created and normalized.\n")
  
  # 2. Process Comparisons from the CSV table
  for (i in 1:nrow(list_of_comparisons)) {
    comp <- list_of_comparisons[i, ]
    factor_col <- comp$Factor
    level_1 <- comp$Normal # Reference
    level_2 <- comp$Treatment # Contrast
    
    contrast_name <- paste0(level_2, "_vs_", level_1)
    cat(paste0("\n--- Analyzing: ", contrast_name, " ---\n"))
    
    contrast_vector <- c(factor_col, level_2, level_1)
    
    # Extract results
    res <- results(dds, contrast = contrast_vector, alpha = padj_treshold)
    res_df <- as.data.frame(res) %>% rownames_to_column("Gene.name")
    
    # Save raw results
    output_deseq2_file <- file.path(output_dir, paste0("results_", contrast_name, ".csv"))
    write_csv(res_df, file = output_deseq2_file)
    
    # Filter for significant DEGs
    degs <- res_df %>%
      filter(padj < padj_treshold)
    
    if (nrow(degs) > 0) {
      output_deg_file <- file.path(output_dir, paste0("DEGs_", contrast_name, ".csv"))
      write_csv(degs, file = output_deg_file)
      cat(paste0("...Found ", nrow(degs), " DEGs. Saved to '", output_deg_file, "'\n"))
    } else {
      cat("...No differentially expressed genes found for this comparison.\n")
    }
    
    # Volcano Plot
    plot_volcano <- EnhancedVolcano(
      res_df,
      lab = res_df$Gene.name,
      x = 'log2FoldChange',
      y = 'pvalue',
      title = paste0('Volcano Plot: ', contrast_name),
      pCutoff = padj_treshold,
      FCcutoff = 1, # Example log2FC threshold
      labSize = 3.0
    )
    
    volcano_file <- file.path(output_dir, paste0("VolcanoPlot_", contrast_name, ".png"))
    ggsave(volcano_file, plot_volcano, width = 8, height = 10)
    cat("Volcano Plot saved.\n")
  }
}

filter_degs_by_pathway <- function(deg_df, pathway_ids, ensembl_table, output_dir, comparison_name) {
  cat("\n--- E1: Filtering DEGs by KEGG Pathway ---\n")
  
  # 1. Map Gene Names to Entrez IDs (needed for KEGG query)
  gene_names <- deg_df$Gene.name
  entrez_ids <- suppressMessages(bitr(
    geneID = gene_names,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ))
  
  # 2. Get genes belonging to input pathways
  genes_in_pathways <- c()
  for (path_id in pathway_ids) {
    kegg_path <- try(keggGet(path_id), silent = TRUE)
    if (!inherits(kegg_path, "try-error") && length(kegg_path[[1]]$GENE) > 0) {
      # Extract Entrez IDs (numeric IDs in the KEGG structure)
      path_genes <- as.character(kegg_path[[1]]$GENE[seq(1, length(kegg_path[[1]]$GENE), by = 2)])
      genes_in_pathways <- c(genes_in_pathways, path_genes)
    }
  }
  
  # 3. Filter DEGs
  genes_in_pathways <- unique(genes_in_pathways)
  filtered_degs <- deg_df %>%
    left_join(entrez_ids, by = c("Gene.name" = "SYMBOL")) %>%
    filter(ENTREZID %in% genes_in_pathways) %>%
    dplyr::select(-ENTREZID) # Remove temporary column
  
  file_name <- file.path(output_dir, paste0("E1_FilteredDEGs_", comparison_name, ".csv"))
  write_csv(filtered_degs, file_name)
  cat(paste("E1 Complete. Found", nrow(filtered_degs), "DEGs in specified pathways. Saved to:", file_name, "\n"))
  
  return(filtered_degs)
}

visualize_pathway_expression <- function(analysis_df, pathway_ids, output_dir, data_type, project_name, sample_table) {
  cat("\n--- E2: Pathway Expression Visualization ---\n")
  
  # 1. Get Gene Symbols from pathways
  gene_list_all <- pathway_gene_list_symbols(pathway_ids) # See new helper function below
  genes_for_plot <- unique(unlist(gene_list_all))
  
  # 2. Prepare expression matrix
  exp_matrix <- analysis_df %>%
    filter(Gene.name %in% genes_for_plot) %>%
    filter(!is.na(Gene.name), !duplicated(Gene.name)) %>%
    column_to_rownames("Gene.name") %>%
    #dplyr::select(-gene_id) %>%
    as.matrix()
  
  if (nrow(exp_matrix) < 5) {
    cat("Skipping visualization: Too few genes found in pathways.\n")
    return(NULL)
  }
  
  # 3. Generate PCA Plot (using the existing plot_pca function)
  cat("Generating PCA for pathway genes...\n")
  plot_pca(sample_table, exp_matrix, gene_subset = genes_for_plot, 
           name_gene_set = paste(project_name, "Pathway Expression"), 
           output_dir = output_dir)
  
  # 4. Generate Heatmap
  cat("Generating Heatmap for pathway genes...\n")
  sample_ids <- intersect(colnames(exp_matrix), sample_table$sample_id)
  annotation_df <- sample_table %>%
    filter(sample_id %in% sample_ids) %>%
    dplyr::select(sample_id, Group) %>%
    column_to_rownames("sample_id") %>%
    rename(Group = Group)
  
  pheatmap(
    mat = exp_matrix,
    annotation_col = annotation_df,
    main = paste("Pathway Gene Expression (", data_type, ")"),
    scale = "row", # Scale by row for gene contrast
    show_rownames = FALSE,
    filename = file.path(output_dir, paste0("E2_Pathway_Heatmap_", data_type, ".png"))
  )
  cat("E2 Complete. Heatmap and PCA saved.\n")
}

enrichment_analysis_ora <- function(deg_df, output_dir, comparison_name, padj_treshold) {
  cat("\n--- E3: Enrichment Analysis (KEGG & GO) ---\n")
  
  # 1. Map Gene Names to Entrez IDs (ORA requires Entrez IDs)
  gene_names <- deg_df$Gene.name
  entrez_ids <- suppressMessages(bitr(
    geneID = gene_names,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ))
  # Filter DEGs to those with valid Entrez IDs
  gene_vector <- entrez_ids$ENTREZID
  
  if (length(gene_vector) < 10) {
    cat("Skipping ORA: Too few DEGs with Entrez IDs for reliable enrichment.\n")
    return(NULL)
  }
  
  # 2. KEGG Enrichment
  kegg_res <- enrichKEGG(
    gene = gene_vector,
    organism = 'hsa', # Assuming Human; should be interactive
    pAdjustMethod = "BH",
    pvalueCutoff = padj_treshold
  )
  
  # 3. GO Enrichment (Biological Process - BP, Cellular Component - CC, Molecular Function - MF)
  go_res <- enrichGO(
    gene = gene_vector,
    OrgDb = org.Hs.eg.db,
    ont = "BP", # Focusing on Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff = padj_treshold,
    readable = TRUE
  )
  
  # 4. Save Results
  write_csv(as_tibble(kegg_res), file.path(output_dir, paste0("E3_KEGG_Enrichment_", comparison_name, ".csv")))
  write_csv(as_tibble(go_res), file.path(output_dir, paste0("E3_GO_Enrichment_BP_", comparison_name, ".csv")))
  
  # 5. Plotting (Dot plots and Bar plots)
  if (nrow(kegg_res) > 0) {
    kegg_dot <- dotplot(kegg_res, showCategory = 10, title = paste0("KEGG Enrichment: ", comparison_name))
    ggsave(file.path(output_dir, paste0("E3_KEGG_DotPlot_", comparison_name, ".png")), kegg_dot, width = 8, height = 5)
  }
  if (nrow(go_res) > 0) {
    go_bar <- barplot(go_res, showCategory = 10, title = paste0("GO (BP) Enrichment: ", comparison_name))
    ggsave(file.path(output_dir, paste0("E3_GO_BarPlot_", comparison_name, ".png")), go_bar, width = 8, height = 5)
  }
  cat("E3 Complete. Enrichment results and plots saved.\n")
}

heatmap_pathway_degs <- function(deg_df, res_df, pathway_ids, output_dir, comparison_name, padj_treshold) {
  cat("\n--- E4: Heatmaps of Pathway DEGs (Log2FC & padj) ---\n")
  
  # 1. Get Gene Symbols from pathways (using the new helper function)
  gene_list_all <- pathway_gene_list_symbols(pathway_ids)
  
  # Prepare the data matrix using LFC and padj from the full results
  plot_data <- res_df %>%
    filter(!is.na(padj)) %>% # Ensure we only plot genes with a result
    dplyr::select(Gene.name, log2FoldChange, padj) %>%
    # Filter for the relevant comparison and ensure no duplicates
    filter(!duplicated(Gene.name))
  
  if (nrow(plot_data) == 0) {
    cat("Skipping E4: No valid log2FoldChange data found.\n")
    return(NULL)
  }
  
  # 2. Iterate through each pathway for individual heatmap
  for (path_name in names(gene_list_all)) {
    path_genes <- gene_list_all[[path_name]]
    
    # Filter the full results data by genes in the current pathway
    pathway_lfc_data <- plot_data %>%
      filter(Gene.name %in% path_genes)
    
    if (nrow(pathway_lfc_data) < 3) {
      next
    }
    
    # 3. Create the LFC Matrix for the Heatmap Body
    lfc_matrix <- pathway_lfc_data %>%
      dplyr::select(Gene.name, log2FoldChange) %>%
      column_to_rownames("Gene.name") %>%
      as.matrix()
    print("hellooooo")
    # 4. Create Annotation for padj (Heatmap Annotation Column/Row)
    # We'll use padj to color code the rows (genes)
    padj_annotation <- pathway_lfc_data %>%
      dplyr::select(Gene.name, padj) %>%
      column_to_rownames("Gene.name") %>%
      mutate(Significance = case_when(
        padj < padj_treshold ~ paste0("Sig (<", padj_treshold, ")"),
        TRUE ~ "Not Sig"
      )) %>%
      dplyr::select(Significance)
    
    # 5. Generate Heatmap (pheatmap is column-centric, so LFC for one comparison is a column)
    # Transpose LFC matrix so genes are rows (standard view)
    pheatmap(
      mat = t(lfc_matrix), # Transposed to make genes rows, LFC a single column
      annotation_row = padj_annotation,
      main = paste0("LFC of DEGs in ", str_split_fixed(path_name, " - ", 2)[1]),
      cluster_rows = FALSE, # Single row, don't cluster
      cluster_cols = TRUE, # Cluster the genes
      show_colnames = FALSE, # Hide the single 'log2FoldChange' column label
      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), # Red/Blue for LFC
      breaks = seq(-max(abs(lfc_matrix)), max(abs(lfc_matrix)), length.out = 100), # Symmetric scale
      filename = file.path(output_dir, paste0("E4_Heatmap_LFC_", gsub(" - .*", "", path_name), "_", comparison_name, ".png"))
    )
  }
  cat("E4 Complete. LFC Heatmaps for pathway DEGs saved.\n")
}

gsea_analysis <- function(analysis_df, output_dir, project_name, padj_treshold) {
  cat("\n--- E5: Gene Set Enrichment Analysis (GSEA) ---\n")
  
  # 1. Load Gene Sets (MSigDB C2, C5, or Custom)
  gsea_choice <- get_user_input("Choose gene set source (MSigDB/Custom): ")
  gmt_file_path <- NULL
  
  if (tolower(gsea_choice) == "msigdb") {
    db_choice <- get_user_input("Choose MSigDB set (C2 for curated pathways, C5 for GO terms): ")
    if (tolower(db_choice) == "c2") {
      # In a real pipeline, you would download/load an actual C2.all.vX.gmt file
      # For simplicity, we'll prompt for the file path if it's not locally available
      gmt_file_path <- get_user_input("Enter FULL path to your downloaded C2.all.vX.gmt file: ")
    } else if (tolower(db_choice) == "c5") {
      gmt_file_path <- get_user_input("Enter FULL path to your downloaded C5.all.vX.gmt file: ")
    }
  } else if (tolower(gsea_choice) == "custom") {
    gmt_file_path <- get_user_input("Enter FULL path to your custom GMT file (.gmt): ")
  }
  
  if (!file.exists(gmt_file_path)) {
    stop("GSEA: Specified GMT file not found.")
  }
  
  # 2. Load Gene Sets
  gene_sets <- getGmt(gmt_file_path)
  
  # 3. Create Ranked Gene List (using median expression or log2FC if available)
  # Since we only have VST (normalized expression), we'll rank by variance (most variable first).
  # NOTE: GSEA typically uses log2FC, but for VST, variance is a good proxy for activity.
  vst_matrix <- analysis_df %>% 
    filter(!is.na(Gene.name), !duplicated(Gene.name)) %>%
    column_to_rownames("Gene.name") %>% 
    dplyr::select(-gene_id) %>% 
    as.matrix()
  
  gene_var <- rowVars(vst_matrix)
  names(gene_var) <- rownames(vst_matrix)
  
  # Rank by variance (high variance = high rank)
  ranked_genes <- sort(gene_var, decreasing = TRUE)
  
  # 4. GSEA Calculation (requires Entrez IDs)
  # Map Gene Symbols to Entrez IDs for the ranked list
  rank_df <- data.frame(SYMBOL = names(ranked_genes), RANK = ranked_genes)
  rank_entrez <- suppressMessages(bitr(
    geneID = names(ranked_genes),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ))
  
  # Create final ranked list using Entrez IDs
  final_rank_df <- rank_df %>%
    inner_join(rank_entrez, by = "SYMBOL") %>%
    # Rank by the original variance metric
    arrange(desc(RANK))
  
  final_ranked_list <- final_rank_df$RANK
  names(final_ranked_list) <- final_rank_df$ENTREZID
  
  cat("Running GSEA (may take a few minutes)...\n")
  gsea_res <- GSEA(
    geneList = final_ranked_list,
    geneSet = gene_sets,
    pvalueCutoff = padj_treshold
  )
  
  # 5. Save Results
  gsea_df <- as_tibble(gsea_res)
  write_csv(gsea_df, file.path(output_dir, paste0("E5_GSEA_Results_", project_name, ".csv")))
  
  cat("E5 Complete. GSEA results saved.\n")
}

pathway_gene_list_symbols <- function(pathway_ids, organism = "hsa") {
  pathway_genes <- list()
  for (path_id in pathway_ids) {
    pathway_info <- try(keggGet(path_id), silent = TRUE)
    if (!inherits(pathway_info, "try-error") && !is.null(pathway_info[[1]]$GENE)) {
      path_name <- pathway_info[[1]]$DESCRIPTION
      gene_info <- pathway_info[[1]]$GENE
      gene_descriptions <- gene_info[seq(2, length(gene_info), by = 2)]
      gene_symbols <- sub(";.*", "", gene_descriptions)
      pathway_genes[[path_name]] <- unique(gene_symbols)
    }
  }
  return(pathway_genes)
}

## ----------------------------------------------------------------------------
## 4. The interactive main function
## ----------------------------------------------------------------------------

#The main interactive function
run_analysis <- function() {
  # ---Step 1: Get Project Info ---
  cat("--- RNA-Seq Analysis Pipeline ---\n")
  
  folderPath <- get_user_input("Enter the FULL path to the folder containing RSEM (.genes.results) files: ")
  if (!dir.exists(folderPath)) {
    stop("Error: Folder path does not exist. PLease check the path.")
  }
  
  project_name <- get_user_input("Enter a Project Name without white spaces: ")
  
  #Load Sample Table
  sample_table_path <- get_user_input("Enter the FULL path to the Sample Metadata Table (.csv file with columns: sample, group, condition): ")
  if (!file.exists(sample_table_path)) {
    stop("Error: Sample table file not found.")
  }
  sample_table <- read_csv(sample_table_path, show_col_types = FALSE) %>%
    rename(sample_id = 1) #ensure the first column is named for consistency !!!NEEDED?!!!
  
  #Load Ensembl Mapping Table (for mapping gene.ids to gene symbols)
  ensembl_path <- get_user_input("Enter the FULL path to the Ensembl ID to Gene Symbol mapping table (.csv with 'gene_id' and 'Gene.name'): ")
  if (!file.exists(ensembl_path)) {
    stop("Error: Ensembl mapping file not found.")
  }
  ensembl_table <- read_csv(ensembl_path, show_col_types = FALSE)
  
  # --- Step 2: Select Analysis Options & Data Type ---
  cat("\nAvailable Analysis Options:\n")
  cat(" A: Data Import & Matrix Creation (Count, TPM, FPKM)\n")
  cat(" B: Normalization (VST, optional Quantile)\n")
  cat(" C: Quality Control Plots (Correlation Heatmap, Density, PCA)\n")
  cat(" D: Differential Gene Expression (DESeq2) - REQUIRES COUNT FILE\n")
  cat(" E: Gene Set/Pathway Analysis (Advanced)\n")
  cat(" F: Shared/Unique Gene Lists\n")
  
  selected_options <- get_user_input("Enter the letters for analyses to run (e.g., ABCDEF): ")
  selected_options <- toupper(selected_options)
  
  data_type_choice <- get_user_input("Choose expression matrix for downstream analysis (TPM/FPKM/VST): ")
  data_type_choice <- toupper(data_type_choice)
  if (!data_type_choice %in% c("TPM", "FPKM", "VST")) {
    stop("Invalid choice. Please select TPM, FPKM, or VST.")
  }
  
  # --- Step 3: Setup Output and Initialize Data Containers ---
  output_dir <- file.path("~", "Desktop", "Analysis_RNA_Seq", project_name)
  count_files_dir <- file.path(output_dir, "count_files")
  dir.create(count_files_dir, recursive = TRUE, showWarnings = FALSE)
  cat(paste("\nOutput directory created:", output_dir, "\n"))
  
  count_df_merged <- NULL # Required for DESeq2 (D) and Normalization (B)
  analysis_df <- NULL    # The matrix chosen for C, E, F (VST, TPM, or FPKM)
  
  
  # A: Data Import & Matrix Creation
  if (grepl("A", selected_options)) {
    cat("\n--- Running Option A: Data Import & Matrix Creation ---\n")
    files <- list.files(path = folderPath, pattern = "\\.genes\\.results$", full.names = TRUE)
    if (length(files) == 0) {
      stop("No RSEM files found in the specified folder.")
    }
    matrix_list <- build_expression_matrices(files, output_dir, ensembl_table, project_name)
    count_df_merged <- matrix_list$count_df
    
    if (data_type_choice == "TPM") analysis_df <- matrix_list$TPM_df
    if (data_type_choice == "FPKM") analysis_df <- matrix_list$FPKM_df
    
    cat("Option A complete.\n")
  }
  
  # --- Pre-Analysis Data Load (Flexible Loading for C, D, E, F) ---
  
  # 1. Load Analysis Data (VST/TPM/FPKM) if required and A was skipped
  if (!grepl("A", selected_options) && grepl("[CEF]", selected_options)) {
    analysis_df <- load_analysis_data(project_name, data_type_choice, count_files_dir)
  }
  
  # 2. Load Count File for D (DESeq2) and B (Normalization), if not run in A
  if (is.null(count_df_merged) && grepl("[BD]", selected_options)) {
    count_file_path <- file.path(count_files_dir, paste0(project_name, "_count_genenames.csv"))
    if (file.exists(count_file_path)) {
      count_df_merged <- read_csv(count_file_path, show_col_types = FALSE)
      cat("Loaded previous COUNT file for DESeq2/Normalization.\n")
    } else if (grepl("D", selected_options)) {
      stop("Cannot run Differential Expression (D). Raw COUNT file is required but not found.")
    }
  }
  
  
  # B: Normalization (VST is mandatory if VST is chosen or D is selected)
  if (grepl("B", selected_options) || (data_type_choice == "VST" && is.null(analysis_df))) {
    
    if (!is.null(count_df_merged)) {
      cat("\n--- Running Option B: Normalization (VST) ---\n")
      
      # Run VST normalization only if analysis_df (VST) is not already loaded
      if (data_type_choice == "VST" && !is.null(analysis_df)) {
        cat("VST data already loaded. Skipping VST calculation.\n")
        vst_df <- analysis_df # Use the loaded data for potential QN step
      } else {
        vst_df <- normalize_vst(sample_table, count_df_merged, output_dir, project_name)
      }
      
      # --- QN Logic ---
      qn_choice <- get_user_input("Perform Quantile Normalization (QN) after VST? (yes/no): ")
      if (tolower(qn_choice) == "yes") {
        # QN function receives the VST data (either freshly calculated or loaded)
        vst_df <- normalize_qn(vst_df, output_dir, project_name)
        cat("VST and Quantile Normalization complete. (File saved with '_qn_vst_names.csv')\n")
      } else {
        cat("VST Normalization complete.\n")
      }
      
      # If VST was the data type chosen, update the analysis_df with the final normalized version
      if (data_type_choice == "VST") analysis_df <- vst_df
      
    } else {
      cat("Skipping B: Cannot run normalization without a COUNT file.\n")
    }
  }
  
  
  # --- Final Check before Analysis Steps ---
  if(is.null(analysis_df) && grepl("[CEF]", selected_options)) {
    stop(paste("Cannot run C, E, or F. The required data type (", data_type_choice, ") is not available."))
  }
  
  
  # C: Quality Control Plots
  if (grepl("C", selected_options) && !is.null(analysis_df)) {
    cat("\n--- Running Option C: Quality Control Plots ---\n")
    
    qc_matrix <- analysis_df %>%
      filter(!is.na(Gene.name)) %>%
      column_to_rownames("Gene.name") %>%
      #dplyr::select(-gene_id) %>%
      as.matrix()
    
    plot_density(qc_matrix, project_name, output_dir, data_type_choice)
    correlation_plot(sample_table, qc_matrix, output_dir)
    plot_pca(sample_table, qc_matrix, gene_subset = 5000, name_gene_set = paste(project_name, "PCA Top 5000"), output_dir)
    cat("Option C complete. QC Plots saved.\n")
  }
  
  # D: Differential Gene Expression (DESeq2)
  if (grepl("D", selected_options) && !is.null(count_df_merged)) {
    cat("\n--- Running Option D: Differential Gene Expression ---\n")
    comp_path <- get_user_input("Enter the FULL path to the CSV file with DESeq2 comparisons: ")
    if (!file.exists(comp_path)) {
      stop("Error: Comparison file not found.")
    }
    list_of_comparisons <- read_csv(comp_path, show_col_types = FALSE)
    
    required_cols <- c("Factor", "Normal", "Treatment")
    if (!all(required_cols %in% colnames(list_of_comparisons))) {
      stop(paste("Comparison file must contain columns:", paste(required_cols, collapse = ", ")))
    }
    padj_treshold <- as.numeric(get_user_input("Enter p-adjust threshold (e.g., 0.05): "))
    
    differentially_expressed_genes(sample_table, count_df_merged, list_of_comparisons, padj_treshold, output_dir)
    cat("Option D complete. DESeq2 results and Volcano Plots saved.\n")
  }
  
  # E: Gene Set/Pathway Analysis
  if (grepl("E", selected_options) && !is.null(analysis_df)) {
    cat("\n--- Running Option E: Advanced Pathway Analysis ---\n")
    
    # 1. Get comparison information
    comparison_name <- get_user_input("Enter the name of the comparison you want to analyze (e.g., Treatment_vs_Normal): ")
    
    # Load FULL DESeq2 results (res_df) and DEGs (deg_df)
    res_file <- file.path(output_dir, paste0("results_", comparison_name, ".csv")) # Full results
    deg_file <- file.path(output_dir, paste0("DEGs_", comparison_name, ".csv")) # Filtered DEGs
    
    res_df <- NULL
    deg_df <- NULL
    
    if (file.exists(res_file)) {
      res_df <- read_csv(res_file, show_col_types = FALSE)
      cat(paste("Loaded full results file:", res_file, "\n"))
      if(file.exists(deg_file)) {
        deg_df <- read_csv(deg_file, show_col_types = FALSE)
        cat(paste("Loaded DEG file:", deg_file, "\n"))
      }
    } else {
      cat("WARNING: Full DESeq2 results file not found. Running DESeq2 (D) is required for E1, E3, and E4.\n")
    }
    
    # 2. Sub-options Menu
    pathway_options <- get_user_input("Which pathway analysis step(s) to run? (1:Filter DEGs, 2:Pathway Viz, 3:Enrichment ORA, 4:DEG Heatmaps, 5:GSEA): ")
    
    # 3. Interactive Pathway List Input (Required for E1, E2, E4)
    pathway_ids <- NULL
    if (grepl("[124]", pathway_options)) {
      pathway_ids_string <- get_user_input("Enter KEGG pathway IDs (e.g., hsa04010,hsa04110) separated by commas: ")
      pathway_ids <- str_split(pathway_ids_string, pattern = ",", simplify = TRUE) %>% str_trim()
    }
    
    # 4. Execute Selected Steps
    
    # E1: Filter DEGs by Pathway
    if (grepl("1", pathway_options) && !is.null(deg_df) && !is.null(pathway_ids)) {
      filter_degs_by_pathway(deg_df, pathway_ids, ensembl_table, output_dir, comparison_name)
    }
    
    # E2: Pathway Visualization (Heatmap & PCA using Normalized data)
    if (grepl("2", pathway_options) && !is.null(pathway_ids)) {
      visualize_pathway_expression(analysis_df, pathway_ids, output_dir, data_type_choice, project_name, sample_table)
    }
    
    # E3: Enrichment Analysis (ORA - KEGG and GO)
    if (grepl("3", pathway_options) && !is.null(deg_df)) {
      padj_treshold <- as.numeric(get_user_input("Enter p-adjust threshold for ORA (e.g., 0.05): "))
      enrichment_analysis_ora(deg_df, output_dir, comparison_name, padj_treshold)
    }
    
    # E4: Heatmaps for Pathway DEGs (Uses LFC & padj from res_df)
    if (grepl("4", pathway_options) && !is.null(res_df) && !is.null(pathway_ids)) {
      padj_treshold <- as.numeric(get_user_input("Enter p-adjust threshold for LFC annotation (e.g., 0.05): "))
      heatmap_pathway_degs(deg_df, res_df, pathway_ids, output_dir, comparison_name, padj_treshold)
    }
    
    # E5: Gene Set Enrichment Analysis (GSEA)
    if (grepl("5", pathway_options)) {
      padj_treshold <- as.numeric(get_user_input("Enter p-adjust threshold for GSEA (e.g., 0.05): "))
      gsea_analysis(analysis_df, output_dir, project_name, padj_treshold)
    }
    
    cat("Option E complete. Advanced Pathway Analysis finished.\n")
  }
  
  # F: Shared/Unique Gene Lists
  if (grepl("F", selected_options) && !is.null(analysis_df)) {
    cat("\n--- Running Option F: Shared/Unique Gene Lists ---\n")
    
    # 1. Get user input for groups and cutoff
    groups_string <- get_user_input("Enter the sample group names to compare (separated by comma, e.g., GroupA,GroupB,GroupC): ")
    groups_to_compare <- str_split(groups_string, pattern = ",", simplify = TRUE) %>% str_trim()
    
    expression_cutoff <- as.numeric(get_user_input(paste0("Enter expression cutoff (min ", data_type_choice, ") for a gene to be considered 'Expressed' in a group (e.g., 1 or 5): ")))
    
    if (length(groups_to_compare) < 2) {
      cat("Skipping F: Must specify at least two groups for comparison.\n")
    } else {
      # 2. Filter expression data to only include the samples in the specified groups
      # First, identify the samples belonging to the requested groups
      samples_in_groups <- sample_table %>%
        filter(group %in% groups_to_compare) %>%
        pull(sample_id)
      
      # Prepare data matrix, including only the samples needed
      exp_data <- analysis_df %>%
        filter(!is.na(Gene.name)) %>%
        dplyr::select(Gene.name, all_of(samples_in_groups))
      
      # 3. Determine 'Expressed' Genes per Group
      expressed_genes_list <- list()
      for (g in groups_to_compare) {
        # Get sample IDs for the current group
        current_samples <- sample_table %>%
          filter(group == g) %>%
          pull(sample_id)
        
        # Filter expression data for this group
        group_exp_data <- exp_data %>%
          dplyr::select(Gene.name, all_of(current_samples))
        
        # A gene is considered expressed if its mean expression across the group 
        # is greater than the cutoff
        expressed_genes_in_group <- group_exp_data %>%
          rowwise() %>%
          mutate(Mean_Expression = mean(c_across(all_of(current_samples)), na.rm = TRUE)) %>%
          ungroup() %>%
          filter(Mean_Expression >= expression_cutoff) %>%
          pull(Gene.name)
        
        expressed_genes_list[[g]] <- expressed_genes_in_group
        cat(paste0("Found ", length(expressed_genes_in_group), " genes expressed in ", g, ".\n"))
      }
      
      # 4. Find Shared Genes (Intersection)
      shared_genes <- Reduce(intersect, expressed_genes_list)
      
      # 5. Find Unique Genes (Set difference)
      unique_genes_df <- data.frame(Gene.name = character(), Group = character())
      for (g in groups_to_compare) {
        # Unique to this group, not shared with any other group
        other_groups_genes <- unlist(expressed_genes_list[setdiff(names(expressed_genes_list), g)])
        unique_to_group <- setdiff(expressed_genes_list[[g]], other_groups_genes)
        
        unique_genes_df <- unique_genes_df %>%
          bind_rows(data.frame(Gene.name = unique_to_group, Group = g))
      }
      
      # 6. Save Results
      shared_df <- data.frame(Gene.name = shared_genes, Type = "Shared")
      
      # Combine and save all
      all_results_df <- shared_df %>%
        bind_rows(unique_genes_df %>% rename(Type = Group)) %>%
        arrange(Type)
      
      file_name <- file.path(output_dir, paste0(project_name, "_Shared_Unique_Genes.csv"))
      write_csv(all_results_df, file = file.path(output_dir, file_name))
      
      cat(paste0("Found ", length(shared_genes), " shared genes.\n"))
      cat(paste0("Shared and unique gene lists saved to: ", file_name, "\n"))
    }
    
    cat("Option F complete.\n")
  }
}
