### Analysis script

## input: rsem files (.genes.results) or count file (.csv)

## libraries
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
library(corrplot)

## command line
folderPath_count_files <- xyz
project_name <- yjz
Groups_comparison <- jgj
Pathways_of_interest <- sds
sample_table <- xyx #should have sample, group, condition
list_of_comparisons <- command_line_3 #should have "column in sample table", "normal", "treatment"
##maybe which normalization method + scale colour function

#### for everything write down how things should be saved

### -------------- COMMAND LINE DO! read in files & create single count, TPM, FPKM files -----------------------------

read_rsem_files <-
  function(folder = folderPath) {
    files <- list.files(path = folderPath, pattern = "\\.genes\\.results$", full.names = TRUE)
    counts_list <- lapply(files, function(f) {
      read_tsv(f, show_col_types = FALSE)[, c("gene_id", "expected_count", "TPM", "FPKM")]
    })
    
    return(counts_list)
  }

build_TPM_FPKM_matrix <- function(files) {
  
  # Extract sample names from file names
  sample_names <- sub("\\.genes\\.results$", "", basename(files))
  
  # Read each file and extract gene_id + TPM or gene_id + FPKM
  TPM_list <- map(files, function(file) {
    read_tsv(file, show_col_types = FALSE) %>%
      dplyr::select(gene_id, TPM)
  })
  
  FPKM_list <- map(files, function(file) {
    read_tsv(file, show_col_types = FALSE) %>%
      dplyr::select(gene_id, FPKM)
  })
  
  counts_list <- map(files, function(file) {
    read_tsv(file, show_col_types = FALSE) %>%
      dplyr::select(gene_id, expected_count)
  })
  # Name each list element by sample
  names(TPM_list) <- sample_names
  names(FPKM_list) <- sample_names
  names(counts_list) <- sample_names
  
  # Rename "expected_count" column in each to match sample name
  TPM_list <- imap(TPM_list, function(df, sample_name) {
    colnames(df)[colnames(df) == "TPM"] <- sample_name
    df
  })  
  FPKM_list <- imap(FPKM_list, function(df, sample_name) {
    colnames(df)[colnames(df) == "FPKM"] <- sample_name
    df
  }) 
  counts_list <- imap(counts_list, function(df, sample_name) {
    colnames(df)[colnames(df) == "expected_count"] <- sample_name
    df
  })
  # Full join all data frames by gene_id
  TPM_df <- purrr::reduce(TPM_list, full_join, by = "gene_id")
  FPKM_df <- purrr::reduce(FPKM_list, full_join, by = "gene_id")
  count_df <- purrr::reduce(counts_list, full_join, by = "gene_id")
  
  # Replace NAs with 0
  TPM_df[is.na(TPM_df)] <- 0
  FPKM_df[is.na(FPKM_df)] <- 0
  count_df[is.na(count_df)] <- 0
  
  path_name <- paste0("~/Desktop/Analysis_RNA_Seq/", project_name, "/count_files/")
  write.csv(TPM_df, file = paste0(path_name, "_TPM.csv")) ### TO-DO: create folder or automatic path; man könnte nach dem Namen der Person/ Projekt fragen und dafür nen Ordner anlegen
  write.csv(FPKM_df, file = paste0(path_name, "_FPKM.csv"))
  colnames(count_df)[0] <- c("gene_id")
  count_df$gene_id <- sub("\\..*", "", count_df$gene_id)
  write.csv(count_df, file = paste0(path_name, "_count.csv"))
  count_df <- merge(x = count_df, y = ensembl_table, by = "gene_id", all.x = TRUE)
  write.csv(count_df, file=paste0(path_name, "_count_genenames.csv"))
  
}

files <- read_rsem_files(folder = Path) ### TO-DO: command line argument
build_TPM_FPKM_matrix(files)

### -------------- COMMAND LINE DO! (sample_table + qn or not) vst normalization -----------------------------
normalization_vst <- 
  function(sample_table, count_file, ensembl) {
    #count file with gene names (column: Gene.name)
    #remove duplicates due to removed version numbers
    count_file <- count_file[!duplicated(count_file$Gene.name), ]
    #make gene names rownnames
    rownames(count_file) <- count_file$Gene.name
    count_file <- select(count_file, -Gene.name)
    #remove rows with small expression
    count_file <- count_file[rowSums(count_file[, -1]) >= 10, ]
    #create DESeq2 object (used for normalization)
    dds <- DESeqDataSetFromMatrix(countData = round(count_file),
                                  colData = sample_table,
                                  design = ~Group)
    #normalize data
    vst_data <- vst(dds, blind = FALSE)
    vst_data <- as.matrix(assay(vst_data))
    
    #safe normalization table on computer
    path_name <- paste0("~/Desktop/Analysis_RNA_Seq/", project_, "/count_files/")
    write.csv(vst_data_names, paste0(path_name, "vst_names.csv"))
    return(vst_data_names)
  }

normalization_qn <-
  function(sample_table, vst_data, ensembl){
    #quantile Normalization - to remove effect of different techniques between samples
    qn_vstdata <- normalize.quantiles(vst_data)
    rownames(qn_vstdata) <- rownames(vst_data)
    colnames(qn_vstdata) <- colnames(vst_data)
    qn_vstdata <- as.data.frame(qn_vstdata)
    
    #safe normalization table on computer
    path_name <- paste0("~/Desktop/Analysis_RNA_Seq/", project_, "/count_files/")
    write.csv(qn_vst_data_names, paste0(path_name, "qn_vst_names.csv"))
  }
### -------------- !!!ADAPT TO DIFFERENT INPUTS! gene expression distribution boxplot -----------------------------

get_matrix <- 
  function(normalized_counts, selectgene){ #selectgene must be a vector
    #selectgene > 0 --> subset counts
    if (is.na(selectgene)) {
      matrix <- matrix[complete.cases(normalized_counts),]
      return(matrix)
    } else {
      matrix <- normalized_counts[selectgene, ]
      matrix <- matrix[complete.cases(matrix)]
      return(matrix)
    }
  }


# Function to plot boxplot
plot_boxplot <-
  function(matrix = get_matrix(normalized_counts, selectgenes),
           colD = sample_table) {
    
    #set rownames of coldata - needed
    rownames(colD) <- colD[, 1]
    
    #calculate euclidean distance between sample in a certain gene set
    euclid <-
      as.matrix(dist(t(matrix), method = "euclidean"))
    
    #subset only the distance to PHH and liver controls
    #comdist <- euclid[grep("PHH|Liver", rownames(euclid)), ]
    
    #transform the euclidean distance into distance-based similarity
    comsim <- (max(euclid) - euclid) / max(euclid)
    
    #remove the rownames, it will mess up downstream table transformation
    rownames(comsim) <- NULL
    
    #we change the table format from "wide" into "long"
    box_dist <-
      melt(comsim, varnames = c("Number", "Names"))
    box_dist <-
      box_dist[, -1] #we don't need the first column
    #adding author and cell type information
    com_box_dist <-
      merge(box_dist,
            colD,
            by.x = "Names",
            by.y = 0,
            all = TRUE)
    
    #combining author and cell type so replicates will be plotted in the same box
    com_box_dist <- com_box_dist %>%
      unite(
        com,
        Author,
        Type,
        Type2,
        sep = "_",
        na.rm = T,
        remove = F
      )
    
    plot <-
      ggplot(com_box_dist,
             aes(
               x = fct_reorder(com, value),
               y = value,
               fill = Type,
               colour = factor(Source)
             )) +

      geom_boxplot() +
      ylim(0, 1) +
      xlab("") +
      ylab(expression("Distance-based similarity (to PHH/Liver)")) +
      theme(axis.text.x = element_text(
        angle = 270,
        vjust = 0.5,
        hjust = 0
      )) +
      scale_colour_manual(values = c("User" = "red",
                                     "Training" = "black"),
                          name = "Data source") +
      scale_fill_manual(values = get_color_scale())
    return(plot)
  }

### -------------- COMMAND LINE DO! (sample_table) correlation Plot -----------------------------

correlation_plot <- 
  function(sample_table, vst_table) {
    
    correlation_matrix <- cor(vst_table, method = "pearson")
    annotation_df <- data.frame(Group=sample_table$group)
    rownames(annotation_df) <- colnames(vst_tables)
    
    correlation_matrix <- correlation_matrix %>%
      filter(rowSums(.) > 300) # Keep rows where the sum of FPKM values (excluding GeneID) is greater than 0
   
     pheatmap(
      mat = correlation_matrix,
      annotation_col = annotation_df,
      annotation_row = annotation_df,
      main = "Sample-to-Sample Correlation Heatmap (VST)",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      border_color = "gray60",
      fontsize = 8,
      display_numbers = FALSE
    )

  }

### -------------- COMMAND LINE DO! (sample_table, gene_subset, name gene set) PCA plot for pathways or top x-----------------------------
plot_pca <- function(sample_table, norm_data, gene_subset, name_gene_set) {
    #define the groups
    intgroup = colnames(sample_table)
    #calculate row variance and order them descending to select the highest variance genes or genes in pathway
    rowvar <- rowVars(as.matrix(norm_data[, ]))
    if (is.numeric(gene_subset)) {
      select <- order(rowvar, decreasing = TRUE)[seq_len(min(ntop, length(rowvar)))]
    } else {
      pathway_genes <- gene_subset$gene_symbol
      select <- intersect(pathway_genes, rownames(norm_data))
    }
    #create the PCA
    pca_res <- prcomp(t(norm_data[select, ]), scale. =FALSE)
    percentVar <- pca_res$sdev ^ 2 / sum(pca_res$sdev ^ 2)
    
    #prepare the label for the PCA plot
    intgroup_df <- as.data.frame(sample_table[, intgroup, drop = FALSE])

    #combine the PCA with the label
    d <- data.frame(
      PC1 = pca_res$x[, 1],
      PC2 = pca_res$x[, 2],
      intgroup_df,
      name = colnames(norm_data[, ])
    )
    
    if (is.numeric(gene_subset)) {
      #plot the PCA
      plot <- ggplot(data=d, aes_string(
        x = "PC1",
        y = "PC2",
        colour = "condition"
      )) +
        #theme_bw() +
        #theme(legend.text = element_text(size=8), legend.position = "bottom")+
        #theme(plot.title = "PCA ")
        geom_point(aes(shape = treatment), size = 2) +
        xlab(paste0("PC1: ", round(percentVar[1]*100), "% variance")) +
        ylab(paste0("PC2: ", round(percentVar[2]*100), "% variance")) +
        coord_fixed() +
        scale_colour_manual(values = get_color_scale())
      return(plot)
    } else {
      #plot the PCA
      plot <- ggplot(data=d, aes_string(
        x = "PC1",
        y = "PC2",
        colour = "condition"
        )) +
        #theme_bw() +
        #theme(legend.text = element_text(size=8), legend.position = "bottom")+
        ggtitle(name_gene_set)+
        theme(plot.title = element_text(size=10, vjust = 15))+
        geom_point(aes(shape = treatment), size = 2) +
        xlab(paste0("PC1: ", round(percentVar[1]*100), "% variance")) +
        ylab(paste0("PC2: ", round(percentVar[2]*100), "% variance")) +
        coord_fixed() +
        scale_colour_manual(values = get_color_scale())
      return(plot) 
      }
  }
 
### -------------- Create lists of uniquely expressed genes and shared genes -----------------------------

#include ven diagram when max 4 lists

unique_genes <- function(sample_table, count_table){
  cat("\n--- Uniquely Expressed Genes per Group ---\n")
  # identify expressed genes for each group
  expressed_genes_by_group <- list()
  
  for (group_name in unique(sample_table$group)) {
    # Get samples belonging to the current group
    current_group_samples <- sample_table %>%
      filter(group == group_name) %>%
      pull(sample)
    
    # Filter FPKM data for these samples
    group_count_data <- count_table %>%
      select(all_of(current_group_samples))
    
    # Identify genes expressed (vst)> threshold) in at least one sample in this group
    # We transpose, apply a check, and then select genes where at least one sample meets the criteria
    expressed_genes <- group_count_data %>%
      rownames_to_column(var="gene_id") %>%
      pivot_longer(
        cols = -gene_id,
        names_to = "Sample",
        values_to = "Value"
      ) %>%
      distinct(gene_id) %>% # Get unique GeneIDs
      pull(gene_id)
    
    # Strip version numbers for the set operations and Venn Diagram plotting
    expressed_genes_by_group[[group_name]] <- expressed_genes_no_version
  }
  
  # Create an empty list to store unique genes for each group
  unique_genes <- list()
  
  # Iterate through each group to find its unique genes
  for (i in seq_along(expressed_genes_by_group)) {
    current_group_name <- names(expressed_genes_by_group)[i]
    current_group_genes <- expressed_genes_by_group[[i]]
    
    # Get genes from all other groups
    other_groups_genes <- unlist(expressed_genes_by_group[-i])
    
    # Find genes unique to the current group
    unique_to_current_group <- setdiff(current_group_genes, other_groups_genes)
    
    # Filter the gene_id_mapping table to get original versioned IDs and names
    unique_to_current_group_with_names <- gene_id_mapping %>%
      filter(gene_id_withoutversion %in% unique_to_current_group_unversioned) %>%
      select(gene_id_version, Gene.name) %>% # Select both original ID and name
      distinct(gene_id_version, .keep_all = TRUE) %>% # Handle potential duplicates
      arrange(gene_id_version) # Order for consistent output
    
    unique_genes_with_names_list[[current_group_name]] <- unique_to_current_group_with_names
    
    cat(paste0("Genes unique to ", current_group_name, " (", nrow(unique_to_current_group_with_names), " genes):\n"))
    path_name <- paste0("~/Desktop/S3/Analysis_Jette/Gene_Expression_Analysis/Coexpression_Analysis/JR_", current_group_name,"_genes_unique.csv")
    write.csv(unique_to_current_group_with_names, path_name)
  }
  return(unique_genes_with_names_list)
}

shared_genes <- function(sample_table, count_table){
  if (num_groups > 1) {
    # Use Reduce with intersect to find common elements across all vectors in the list
    genes_shared_by_all_unversioned_IDs <- Reduce(intersect, expressed_genes_by_group)
    if (length(genes_shared_by_all_unversioned_IDs) > 0) {
      #merge with gene_names_df to get gene names
      genes_shared_by_all_with_names <- gene_id_mapping %>%
        filter(gene_id_withoutversion %in% genes_shared_by_all_unversioned_IDs) %>%
        select(gene_id_version, Gene.name) %>%
        arrange(gene_id_version)
      cat(paste0("Total genes shared by all groups: ", nrow(genes_shared_by_all_with_names), "\n"))
    } else {
      cat("No genes are shared by all groups.\n")
    }
  } else {
    cat("Cannot compute shared genes for less than 2 groups.\n")
  }
  write.csv(genes_shared_by_all_with_names, file="~/Desktop/S3/Analysis_Jette/Gene_Expression_Analysis/Coexpression_Analysis/JR_genes_shared.csv")
  return(genes_shared_by_all_with_names)
}
### -------------- !!!COMMAND LINE: GROUPS! Differential gene expression analysis -----------------------------
differentially_expressed_genes <- 
  function(sample_table, count_file, list_of_comparisons, padj_treshold) {

    #remove duplicates due to removed version numbers
    count_file <- count_file[!duplicated(count_file$Gene.name), ]
    #make gene names rownnames
    rownames(count_file) <- count_file$Gene.name
    count_file <- select(count_file, -Gene.name)
    #remove rows with small expression
    count_file <- count_file[rowSums(count_file[, -1]) >= 10, ]
    #create DESeq2 object (used for normalization)
    dds <- DESeqDataSetFromMatrix(countData = round(count_file),
                                  colData = sample_table,
                                  design = ~group)
    cat("\nRunning DESeq2 analysis...\n")
    dds <- DESeq(dds)
    cat("DESeq2 analysis complete.\n")

    
    for (comp in list_of_comparisons) {
      contrast_name <- paste0(comp[2], "_vs_", comp[3])
      cat(paste0("\n--- Analyzing: ", contrast_name, " ---\n"))
      ## define significance tresholds
      padj_treshold <- 0.05
      
      # Extract results
      res <- results(dds, contrast, comp)
      res_df <- as.data.frame(res)
      output_deseq2_file <- paste0("~/Desktop/S3/Analysis_Jette/Differential_Gene_Expression_Analysis/results_", contrast_name, ".csv")
      write.csv(res, file = output_deseq2_file, row.names = FALSE)
      
      cat("Raw DESeq2 results (first 6 rows):\n")
      print(head(res_df))
      
      # filter for significant differentially expressed genes
      degs <- res_df %>%
        filter(padj<padj_treshold)
      cat(paste0("\nNumber of DEGs (padj < ", padj_threshold, "): ", nrow(degs), "\n"))
      output_deg_file <- paste0("~/Desktop/S3/Analysis_Jette/Differential_Gene_Expression_Analysis/DEGs_", contrast_name, ".csv")
      write.csv(degs, file = output_deg_file, row.names = FALSE)
      
      # if no names in yet then merge with ensemble again... (TRY NOT TO)
      degs_with_names <- merge(x = degs, y = ensembl_table, by = "gene_id", all.x = TRUE)
      
      
    }
        
        # Save DEGs to a CSV file
        output_deg_file_name <- paste0("~/Desktop/S3/Analysis_Jette/Differential_Gene_Expression_Analysis/DEGs_", contrast_name, "_with_names.csv")
        write.csv(degs_with_names, file = output_deg_file_name, row.names = FALSE)
        cat(paste0("\nDifferentially expressed genes saved to '", output_deg_file_name, "'\n"))
       else {
        cat("No differentially expressed genes found for this comparison at the specified thresholds.\n")
      }
    
  }
### -------------- COMMAND LINE: PATHWAYS! Read interesting pathways in -----------------------------

### command line argument with list of pathways or NULL
### command line argument with organism like mouse or human
pathway_gene_list <- 
  function(pathways = list_pathways, organism = "") {
    pathways <- read.csv("/Users/new_jette/Desktop/S3/Analysis_Costina/pathways.csv", header=TRUE)
    if (organism == "mouse") {
      pathways <- keggList("pathway", "mm")
    } else {
      pathways <- keggList("pathway", "hsa")
    }
    
    pathway_ids <- pathways$Pathway_ID
    pathway_genes <- list()
    for (path in pathway_ids) {
      print(path)
      pathway_info <- keggGet(path)
      genes_data <- pathway_info[[1]]$GENE
      gene_descriptions <- genes_data[seq(2, length(genes_data), by = 2)]
      gene_symbols <- sub(";.*", "", gene_descriptions)
      pathway_genes[[pathways[path]]] <- gene_symbols
    }
    
    return(pathway_genes)
  }


### -------------- Pathway analysis based on gene expression -----------------------------
## PCA per pathway
## heatmap for pathways with less than 50 genes
## ssGSEA?

### -------------- Pathway analysis based on differential expressed genes -----------------------------
## table of how many significant DEGs in a pathway
## PCA per pathway
## heatmap for pathways with less than 50 genes
## GSEA





## Possibilities: Volcano-Plot noch hinzufügen
