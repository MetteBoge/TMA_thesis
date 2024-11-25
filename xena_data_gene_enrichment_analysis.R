#!/zhome/94/f/147417/miniconda3/bin/R

rm(list = ls())
set.seed(123)  # Use any consistent value

## Load libraries and data
library(tidyverse)
#library(GSVA)
#library(Hmisc) # rcorr
#library(reticulate) #using python
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(kableExtra)
library(Cairo)
library(rmarkdown)
library(webshot2) # save_kable
library(ggridges)
packageVersion("msigdbr")

# Control operations: 
load_tpm <- FALSE
plot_mds <- FALSE

cancer_abbr_list <- c("LAML", "CHOL", "CESC", "DLBC", "HNSC", "BRCA", "PAAD", "SKCM", "LUAD", "LUSC")


file_path_xenadata <- "/zhome/94/f/147417/tma_thesis/data/xena/"
path_save_plots <- "/zhome/94/f/147417/tma_thesis/data/GSEA_plots/"

file_name_pheno1 <- "TCGA_PANCAN_sampleType_primaryDisease"
file_name_pheno2 <- "TCGA_PANCAN_curated_clinical_data"


target_list <- c("CD274", "PDCD1LG2", "IDO1", "TDO2", "ARG1", "ARG2", "TGFB1", "CCL22", "CCL2", "LGALS3","LRRC32", "IL10", "SIGLEC15", "BIRC7")

# Phenotype table for simple sample filtering: 
df_st_pd <- read.csv(paste(file_path_xenadata, 
                           file_name_pheno1, ".tsv", sep = ""),sep = "\t")

df_clin_pd <- read.csv(paste(file_path_xenadata, 
                             file_name_pheno2, ".tsv", sep = ""),sep = "\t")

df_pheno <- df_clin_pd %>% 
  full_join(df_st_pd, by = "sample")

### GO, gene set


# Load GO Biological process data: msigdbr_collections() to see category
# Using package version '7.5.1'
# Using GO:BP
GO_database <- msigdbr(species = "Homo sapiens",
                       category = "C5",
                       subcategory = "GO:BP")





### Filter samples:

run_gsea <- function(c_type, sub_ctype){
  # Collect sample IDs in lists
  sample_id_of_interest <- list()

  #for (c_type in cancer_abbr_list){
  df_pheno_c_type <- df_pheno %>% 
    filter(cancer.type.abbreviation == c_type)

  # Filter sample ids to only include those of tumor tissue. For LAML: only PBMC. 
  # For Melanoma, both metastatic and tumor. 
  if (c_type == "LAML"){
    # Load data tables to list: 
    phenotype_ctype <- list()
    file_name_pheno <- paste0(file_path_xenadata, "TCGA_", c_type, "_phenotype.tsv")
    phenotype_ctype[[c_type]] <- read.csv(file_name_pheno, sep = "\t") %>% 
      dplyr::rename(sample = sampleID)
    
  # Choose those who do not have atra_exposure
  no_atra_exp <- phenotype_ctype[[c_type]] %>% 
    filter(atra_exposure != "YES") %>% 
    dplyr::select(sample) %>% 
    unlist(use.names = F)
  
  sample_id_of_interest[[c_type]] <- df_pheno_c_type %>% 
    filter(sample_type == "Primary Blood Derived Cancer - Peripheral Blood",
            sample %in% no_atra_exp) %>%
    dplyr::select(sample)
  
  } else if(c_type == "SKCM"){
    
    if (sub_ctype == "SKCM_tumor"){
      sample_id_of_interest[["SKCM"]] <- df_pheno_c_type %>% 
        filter(sample_type == "Primary Tumor") %>%
        dplyr::select(sample)
      
    } else if (sub_ctype == "SKCM_metastatic"){
      sample_id_of_interest[["SKCM"]] <- df_pheno_c_type %>% 
        filter(sample_type == "Metastatic") %>%
        dplyr::select(sample)
    }

    # Filter out Male if BRCA (as MOrten told us to)
  } else if(c_type == "BRCA"){
      sample_id_of_interest[[c_type]] <- df_pheno_c_type %>% 
      filter(sample_type == "Primary Tumor", 
            gender == "FEMALE") %>%
      dplyr::select(sample) 
    
  } else {
    sample_id_of_interest[[c_type]] <- df_pheno_c_type %>% 
      filter(sample_type == "Primary Tumor") %>%
      dplyr::select(sample) 
  }


  ### Load gene expression data


  # Load TPM gene exp data: 

  if (!exists("df_tpm_data")){
    df_tpm_data <- list()
    
    df_tpm_data[[c_type]] <- read_tsv(paste0(file_path_xenadata,
                                            "TCGA_pancan_gene_expression_tpm_",
                                            c_type, ".tsv"), 
                                  show_col_types = FALSE)
    
    # samples to select (with phenotype, gene exp data and correct sample type)
    samples_to_select <- intersect(unlist(sample_id_of_interest[[c_type]]),
                                  colnames(df_tpm_data[[c_type]]))
    
    # Filter samples out which are not of interest: 
    df_tpm_data[[c_type]] <- df_tpm_data[[c_type]] %>% 
      dplyr::select(sample, all_of(samples_to_select))
    
    # Change gene IDs to gene names: 
    geneName_ensg <- read_csv("/zhome/94/f/147417/tma_thesis/data/targets_ensg.csv")
    geneName_ensg <- geneName_ensg %>% 
      dplyr::rename(gene_ensg_shrt = converted_alias) %>% 
      dplyr::select(name, gene_ensg_shrt)
    
    df_tpm_data[[c_type]] <- df_tpm_data[[c_type]] %>% 
      mutate(gene_ensg_shrt = gsub("\\..*", "", sample)) %>% 
      dplyr::select(-sample) %>% 
      relocate(gene_ensg_shrt) %>% 
      left_join(geneName_ensg, by = "gene_ensg_shrt") %>%
      dplyr::rename(sample = name) %>% 
      relocate(sample) %>% 
      mutate(sample = ifelse(is.na(sample), gene_ensg_shrt, sample)) %>% 
      dplyr::select(-gene_ensg_shrt)
  }


  ## Correlate genes

  #Make pairwise correlation of our target gene vs all other genes (one target gene at a time)


  # Create an empty list to store correlation results
  correlation_results <- list()

  gene_expression_matrix <- df_tpm_data[[c_type]] %>% 
    column_to_rownames("sample") %>% 
    as.matrix()

  ### Filter out genes that have only one unique gene exp value for all samples: 
  # Identify genes with more than one unique value across samples
  variable_genes <- apply(gene_expression_matrix, 1, function(x) length(unique(x)) > 1)
  # Filter the gene expression matrix
  gene_expression_matrix <- gene_expression_matrix[variable_genes, ]

  # Pass the R object to Python
  #py$gene_expression_matrix = gene_expression_matrix



  ### One gene: correlation, GO



  R_element_name <- ifelse(c_type == "SKCM", 
                          paste0(file_path_xenadata, 'GSEA_Rdata/gene_correlation_dataframes_', sub_ctype, '.RData'),
                            paste0(file_path_xenadata, 'GSEA_Rdata/gene_correlation_dataframes_', c_type, '.RData'))

  if (!file.exists(R_element_name)){

    # List to save for each target gene
    correlation_stats_2 <- list()
    n_all_genes <- list()
    n_sig_genes <- list()
    n_sig01_genes <- list()
    n_sig001_genes <- list()
    
    # ChatGPT suggest using FALSE for large sample sizes, hence it will assume a distribution and be more computationally efficient. For smaller sample size using TRUE will give more reliable p values. ChatGPT suggested a cutoff of 129. 
    if (ncol(gene_expression_matrix) > 129){
      compute_exact_p_val <- FALSE
    } else {
      compute_exact_p_val <- TRUE
    }
    
    
    for (gene_x in target_list){
      print(paste("Running correlation analysis and adding GO terms for", gene_x))
      print(paste("Using compute_exact_p_val =", compute_exact_p_val, ", and suppressing warnings related to cor.test"))
      
      op <- options(warn = (-1)) # suppress warnings
      # Compute p-values for correlations (using cor.test for each gene)
      correlation_stats <- apply(gene_expression_matrix, 1, function(gene) {
        test <- cor.test(gene_expression_matrix[gene_x, ], 
                        gene, 
                        method = "spearman",
                        exact = compute_exact_p_val)
        c(correlation = test$estimate, p.value = test$p.value)
      })
      options(op) # reset the default value, if you want

      p_values <- correlation_stats["p.value",]
      
      # Adjust p-values using Benjamini-Hochberg (FDR control)
      adjusted_p_values <- p.adjust(p_values, method = "BH")
      
      correlation_stats_2[[gene_x]] <- rbind(correlation_stats, adjusted_p_values) %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate(abs_rho = abs(correlation.rho)) %>% 
        rownames_to_column("gene")
        
      # Target genes only
      correlation_stats_2_TG <- correlation_stats_2[[gene_x]] %>% 
        filter(!grepl("^ENSG", gene)) %>% 
        mutate(gene_symbol = gene)
      
      # All genes
      correlation_stats_2_AG <- correlation_stats_2[[gene_x]] %>% 
        filter(grepl("^ENSG", gene)) %>% 
        mutate(ensembl_gene = gene)
      
      # Add GO information:
      correlation_stats_2_TG <- correlation_stats_2_TG %>% 
        left_join(GO_database, 
                  by="gene_symbol")
      correlation_stats_2_AG <- correlation_stats_2_AG %>% 
        left_join(dplyr::select(GO_database, -gene_symbol), 
                  by="ensembl_gene") 
      
      # Merge for AG and TG again, remove genes without GO info (e.g. ENSG00000109790) and trim columns
      correlation_stats_2[[gene_x]] <- bind_rows(correlation_stats_2_TG,
                correlation_stats_2_AG) %>%
        arrange(-abs_rho) %>% 
        filter(!is.na(gs_cat)) %>% # Filter out genes without GO info
        # Remove all columns with only same unique value for all rows: 
        dplyr::select(where(~ n_distinct(.) > 1 )) %>% 
        dplyr::select(-starts_with("human")) %>% 
        # Remove GOBP in front of GS name: 
        mutate(gs_name = gsub("GOBP_", "", gs_name))
    }

    saveRDS(correlation_stats_2, file = R_element_name)
    
  } else {
    
    correlation_stats_2 <- readRDS(R_element_name)
    
  }




  ### GSEA



  R_element_name_2 <- ifelse(c_type == "SKCM", 
                          paste0(file_path_xenadata, 'GSEA_Rdata/gsea_results_', sub_ctype, '.RData'),
                            paste0(file_path_xenadata, 'GSEA_Rdata/gsea_results_', c_type, '.RData'))

  if (!file.exists(R_element_name_2)){
    
    gsea_results <- list()
    
    for (gene_x in target_list){
      print(paste("Running GSEA analysis for", gene_x))
      
      # Remove gene itself from correlation: 
      cor_df <- correlation_stats_2[[gene_x]] %>% 
        filter(gene != gene_x)
      
      # Prepare ranked gene list
      ranked_genes <- cor_df %>%
        filter(adjusted_p_values < 0.05) %>% 
        arrange(desc(correlation.rho)) %>%  # Sort by correlation
        select(gene, correlation.rho) %>% unique() %>% 
        deframe() 
      
      # If doesnt correlate with any other genes (significantly)
      if (length(ranked_genes) < 1){
        warning(paste(gene_x, "for", c_type, "does not correlate with any other genes significantly (adj.p.val). No GSEA analysis to be carried out."))
        
      } else {
        # Perform GSEA
        gsea_results[[gene_x]] <- GSEA(
          geneList = ranked_genes, 
          TERM2GENE = cor_df %>% select(gs_name, gene),
          pvalueCutoff = 0.05
        )
      }
    }
    
    saveRDS(gsea_results, file = R_element_name_2)
    
  } else {
    
    gsea_results <- readRDS(R_element_name_2)
  }




  ### GSE plots and tables


  # Function to add a newline after the middle word if nchar > 20
  add_line_break <- function(string) {
    if (nchar(string) > 20) {
      # Split string into words
      words <- strsplit(string, " ")[[1]]
      # Calculate the middle index
      mid_idx <- ceiling(length(words) / 2)
      # Insert a newline after the middle word
      paste(c(words[1:mid_idx], "\n", words[(mid_idx + 1):length(words)]), collapse = " ")
    } else {
      string  # Return string as-is if nchar <= 20
    }
  }

  ####
  # Function making first letter capital: 
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }



  #gene_x <- target_list[11]

  summary_table <- list()

  for (gene_x in names(gsea_results)){
    print(paste("Plotting GSEA results for", gene_x))
    
    gsea_res_gene <- gsea_results[[gene_x]]
    gsea_df <- as.data.frame(gsea_res_gene)
    
    if (nrow(gsea_res_gene) > 0){
      # Visualize GSEA results
      
      ##############
      ## DOT PLOT ##
      # # # # # # ##
      gsea_plot <- dotplot(gsea_res_gene, showCategory = 10) +
        scale_y_discrete(labels = function(x) tools::toTitleCase(tolower(gsub("_", " ", x)))) + # Convert y-axis labels to Title Case 
        ggtitle(paste("GSEA Results for", gene_x, "in", c_type))  # Add a title here
      
      ggsave(filename = ifelse(c_type == "SKCM",
                      paste0(path_save_plots, sub_ctype,
                      "_dotplot_", gene_x,
                      ".pdf"),
                      paste0(path_save_plots, 
                      c_type, "_dotplot_", gene_x,
                      ".pdf")),
            width = 10,  # Adjust width as needed
            height = 6,  # Adjust height as needed
            units = "in",
            device = cairo_pdf)
      
      
      
      ####################
      ## ENRICHMENT MAP ##
      # # # # # # # # ## #
      gsea_res_gene_tmp <- gsea_res_gene
      
      gsea_res_gene_tmp@result$Description <- gsea_res_gene_tmp@result$Description %>% 
        tolower() %>%  
        tools::toTitleCase() %>% 
        gsub("_", " ", .) %>% 
        sapply(add_line_break)
      
      # Compute pairwise similarity of enriched terms
      term_sim <- pairwise_termsim(gsea_res_gene_tmp)
      
      emapplot(term_sim, showCategory = 10, min_edge = 0.1) + 
        ggtitle(paste("GSEA Results for", gene_x, "in", c_type))  # Add a title here
      
      ggsave(filename = ifelse(c_type == "SKCM",
                    paste0(path_save_plots, sub_ctype,
                    "_enrichment_map_", gene_x,
                    ".pdf"),
                    paste0(path_save_plots, 
                    c_type, "_enrichment_map_", gene_x,
                    ".pdf")),
          width = 10,  # Adjust width as needed
          height = 6,  # Adjust height as needed
          units = "in",
          device = cairo_pdf)
      
    
      ###############
      ## RIDGEPLOT ##
      # # # # # # # #
      ridgeplot(gsea_res_gene) + labs(x = "enrichment distribution") +
      scale_y_discrete(labels = function(x) x %>% 
                        gsub("_", " ", .) %>% 
                        tolower() %>% 
                        tools::toTitleCase() %>% 
                        sapply(add_line_break)) + 
        ggtitle(paste("GSEA Results for", gene_x, "in", c_type))  # Add a title here
      
      ggsave(filename = ifelse(c_type == "SKCM",
                    paste0(path_save_plots, sub_ctype,
                    "_ridgeplot_", gene_x,
                    ".pdf"),
                    paste0(path_save_plots, 
                    c_type, "_ridgeplot_", gene_x,
                    ".pdf")),
          width = 10,  # Adjust width as needed
          height = 6,  # Adjust height as needed
          units = "in",
          device = cairo_pdf)
      
      
      ###################
      ## SUMMARY TABLE ##
      # # # # # # # # # #
      # gsea_df already only contains significant genesets only ! (adj.p.val < 0.05)
      cor_df <- correlation_stats_2[[gene_x]] %>% 
        filter(gene != gene_x)
      
      gsea_df_2 <- gsea_df %>% 
        mutate(Description = Description %>% 
                tolower() %>% 
                gsub("_", " ", .) %>% 
                tools::toTitleCase() %>% 
                sapply(firstup)) 
      
      # Calculate summary metrics
      summary_stats <- list(
        total_gene_sets_tested = length(unique(cor_df$gs_name)),  # Total gene sets analyzed
        significant_gene_sets = nrow(gsea_df_2),  # Significant gene sets
        top_NES_pathways = gsea_df_2 %>%
          arrange(-abs(NES)) %>%
          slice(1:10) %>%
          select(Description, NES, p.adjust),  # Top 10 pathways
        top_sig_pathways = gsea_df_2 %>%
          arrange(p.adjust) %>%
          slice(1:10) %>%
          select(Description, NES, p.adjust),
        genes_per_set = summary(gsea_df_2$setSize),  # Average size of significant sets
        total_genes_tested = cor_df$gene %>% unique() %>% length(),  # Total genes tested in correlation
        significantly_correlated_genes = cor_df %>% 
          filter(adjusted_p_values < 0.05) %>% 
          dplyr::select(gene) %>% unique() %>% nrow() # Significant genes
      )
      
      # Top NES pathways (converted to a long string)
      top_NES_pathways_str <- summary_stats$top_NES_pathways %>%
        mutate(entry = paste0(Description, " (NES: ", round(NES, 2), 
                              ", p.adjust: ", signif(p.adjust, 2), ")")) %>%
        pull(entry) %>%
        paste(collapse = "\n ")
      
      # Top significant pathways (converted to a long string)
      top_sig_pathways_str <- summary_stats$top_sig_pathways %>%
        mutate(entry = paste0(Description, " (NES: ", round(NES, 2), 
                              ", p.adjust: ", signif(p.adjust, 2), ")")) %>%
        pull(entry) %>%
        paste(collapse = "\n ")
          
          
      # Combine all summary information into one data frame
      tmp_table <- tibble(
        Metric = c(
          "Total Gene Sets Tested",
          "Significant Gene Sets",
          "Total Genes Tested",
          "Significantly Correlated Genes",
          "Genes per Set (Set Size)",
          "Top 10 Pathways by NES",
          "Top 10 Most Significant Pathways"
        ),
        Value = c(
          # Main metrics
          summary_stats$total_gene_sets_tested,
          summary_stats$significant_gene_sets,
          summary_stats$total_genes_tested,
          summary_stats$significantly_correlated_genes,
          paste0("Min: ", summary_stats$genes_per_set[["Min."]],
                ", Q1: ", summary_stats$genes_per_set[["1st Qu."]],
                ", Median: ", summary_stats$genes_per_set[["Median"]],
                ", Mean: ", summary_stats$genes_per_set[["Mean"]],
                ", Q3: ", summary_stats$genes_per_set[["3rd Qu."]],
                ", Max: ", summary_stats$genes_per_set[["Max."]]),
          top_NES_pathways_str,
          top_sig_pathways_str
        )
      )
      
        
      tmp_table <- tmp_table %>% 
        as.matrix() %>% 
        t() %>% 
        as.data.frame() 
      
      colnames(tmp_table) <- tmp_table[1,]
      
      summary_table[[gene_x]] <- tmp_table[-1, ] %>% 
        mutate(Gene = gene_x) %>% 
        relocate(Gene)
      
      rownames(summary_table[[gene_x]]) <- NULL

    } else {
      warning(paste(gene_x, "has no significantly enriched gene sets!"))
    }
  }

  full_summary_table <- summary_table %>% 
    bind_rows()


  # Make pretty table
  table <- full_summary_table %>% 
    kable(caption = paste("Summary of GSEA analysis for", 
                          c_type, "of", gene_x),
          escape = F,
          format = "markdown") %>%   
    kable_styling(full_width = T, position = "center", fixed_thead = T,
                  latex_options = c("striped", "hold_position", "longtable")) %>%
    column_spec(1, width = "5cm") %>%
    column_spec(2, width = "25cm", extra_css = "white-space: pre-line;") 

  table_file_name <- ifelse(c_type == "SKCM",
                            paste0(path_save_plots,
                                  sub_ctype,
                                  "_summaryTable.html"),
                            paste0(path_save_plots,
                                  c_type,
                                  "_summaryTable.html"))

  table %>% 
    save_kable(file = table_file_name)

}






for (c_type in cancer_abbr_list){
  print(c_type)
  if (c_type == "SKCM"){
    sub_ctype_list <- c("SKCM_tumor", "SKCM_metastatic")
    for (sub_ctype in sub_ctype_list){

      run_gsea(c_type, sub_ctype)
      }
  } else {
    sub_ctype <- "SKCM_tumor" # not used
    run_gsea(c_type, sub_ctype)
  }

}
