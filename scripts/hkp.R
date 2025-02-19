# Mai Nguyen 2/2/2025
# Human housekeeping proteome (HKP) analysis 

# Import required lib
library(data.table)
library(tidyr)
library(dplyr)
library(umap)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(tibble)
library(purrr)
library(clusterProfiler)
library(org.Hs.eg.db)  
library(DOSE) 
# Set seed
set.seed(30)

# Functions
# Function to read TSV files 
read_data = function(filepath) {return(read.csv(filepath, sep = '\t', header = TRUE))}

# Function to reshape dataset into wide format
reshape_wide = function(data, sample_column) {data = as.data.table(data)
  # Rename 'Gene.name' to 'Gene_name'
  setnames(data, old = 'Gene.name', new = 'Gene_name')
  # Reshape the dataset
  wide_format = dcast(data, Gene + Gene_name ~ get(sample_column), value.var = 'nTPM', fill = 0)
  return(wide_format)}

# Function to reshape dataset into long format
reshape_long = function(data) {
  # Gene column is unique (ENSG-ID) instead of Gene_name
  data_no_gene = dplyr::select(data, -Gene_name)  
  # Transpose data to long format
  long_format = as.data.frame(t(data_no_gene))  
  # Set first row as column names (Genes)
  colnames(long_format) = long_format[1, ]
  long_format = long_format[-1, ] 
  # Convert all values to numeric
  long_format[] = lapply(long_format, function(x) as.numeric(as.character(x)))
  return(long_format)}

# Function to plot UMAP and return xy coordinates
compute_umap = function(data, dataset_name) {
  data = as.data.frame(data)
  umapxy = umap(data, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')
  umap_xy = data.frame(UMAP1 = umapxy$layout[, 1], UMAP2 = umapxy$layout[, 2])
  plot_title = paste('UMAP', dataset_name)
  umap_plot = ggplot(umap_xy, aes(x = UMAP1, y = UMAP2)) +  
    geom_point(size = 2, alpha = 0.4, shape = 21, 
               fill = 'lightblue', color = 'dodgerblue3', stroke = 0.2) +  
    labs(title = plot_title, x = 'UMAP1', y = 'UMAP2') + 
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
          axis.title = element_text(size = 14))
  return(list(umap_xy = umap_xy, umap_plot = umap_plot))}

# Function to plot PCA and return PCA object
compute_pca = function(data, dataset_name) {
  data = as.data.frame(data)
  pca_result = PCA(data, graph = FALSE)
  plot_title = paste('PCA', dataset_name)
  pca_plot = fviz_pca_ind(pca_result, geom.ind = 'point', pointsize = 2, alpha.ind = 0.7, 
                          col.ind = 'lightblue') +
    labs(title = plot_title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
          axis.title = element_text(size = 14))
  return(list(pca_result = pca_result, pca_plot = pca_plot))}

# Function to calculate expression percentage
calculate_expression_percentage = function(data, ntpm_threshold = 1, percentage_threshold = 100) {
  num_samples = ncol(data) - 2  
  binary_matrix = if (ntpm_threshold == 0) data[, -c(1,2)] > 0 else data[, -c(1,2)] >= ntpm_threshold
  expressed_counts = rowSums(binary_matrix)
  expressed_percentage = (expressed_counts / num_samples) * 100
  result = data.frame(Gene = data[,1], Gene_name = data[,2], Expressed_percentage = expressed_percentage)
  always_expressed_genes = result %>%
    filter(Expressed_percentage >= percentage_threshold) %>%
    select(Gene, Gene_name)
  return(list(expression_data = result, always_expressed = always_expressed_genes))}

# Function to retrieve gene counts dynamically
get_gene_counts = function(dataset_name, ntpm, perc_threshold) {
  key_ntpm = paste0('ntpm_', ntpm)
  key_perc = paste0('ntpm_', ntpm, '_perc_', perc_threshold)
  results_list[[key_ntpm]][[key_perc]][[dataset_name]]$always_expressed %>% nrow()}

# Function to compute overlaps or gene shared between 4 data sets
compute_gene_overlap = function(..., dataset_names = NULL) {
  datasets = list(...)
  # Extract only Gene IDs for each dataset
  gene_lists = lapply(datasets, function(df) df$Gene)
  # Create a unique list of all genes across datasets
  all_genes = data.frame(Gene = unique(unlist(gene_lists)))
  # Add boolean columns for presence in each dataset
  for (i in seq_along(datasets)) {
    dataset_name = ifelse(is.null(dataset_names), paste0('Dataset_', i), dataset_names[i])
    all_genes[[dataset_name]] = all_genes$Gene %in% gene_lists[[i]]}
  # Count the number of datasets each gene appears in
  all_genes$Overlap_Count = rowSums(all_genes[, -1])  
  # Categorize genes based on overlap level
  all_genes$Overlap_Category = case_when(
    all_genes$Overlap_Count == length(datasets) ~ 'Shared by All',
    all_genes$Overlap_Count == length(datasets) - 1 ~ 'Shared by 3',
    all_genes$Overlap_Count == length(datasets) - 2 ~ 'Shared by 2',
    TRUE ~ 'Unique to One')
  shared_by_all = all_genes %>%
    filter(Overlap_Category == 'Shared by All') %>%
    select(Gene)
  return(list(overlap_summary = all_genes, shared_genes = shared_by_all))}

# Function to plot a stacked bar of overlaps
plot_gene_overlap = function(overlap_data) {
  overlap_summary = overlap_data %>%
    select(-Gene) %>%
    pivot_longer(cols = -c(Overlap_Count, Overlap_Category), 
                 names_to = 'Dataset', values_to = 'Present') %>%
    filter(Present) %>%
    count(Dataset, Overlap_Category)
  overlap_summary$Overlap_Category = factor(overlap_summary$Overlap_Category, 
                                            levels = c('Unique to One', 'Shared by 2', 'Shared by 3', 'Shared by All'))
    ggplot(overlap_summary, aes(x = Dataset, y = n, fill = Overlap_Category)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = c('Shared by All' = 'darkblue', 
                                 'Shared by 3' = 'blue', 
                                 'Shared by 2' = 'lightblue', 
                                 'Unique to One' = 'gray')) +
    labs(title = 'Gene Overlap Across Datasets',
         x = 'Dataset', y = 'Number of Genes', fill = 'Overlap Level') +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = 'bold'))}

# Function to filter dataset by a list and specified column
filter_by_list = function(data, filter_list, filter_column) {
  return(data %>% filter(!!sym(filter_column) %in% filter_list))}

# Function to prepare dataset for UMAP (set Gene as row names, remove Gene_name)
prepare_for_umap = function(data, row_name_column = 'Gene') {
  data = data %>% column_to_rownames(var = row_name_column)
  if('Gene_name' %in% colnames(data)) {  
    data = data %>% select(-Gene_name)  }
  return(data)}

# Data analysis
# Pre-processing 
# Read data
hpa_rna_single_cell = read_data('data/rna_single_cell_type.tsv')
hpa_rna_brain = read_data('data/rna_single_nuclei_cluster_type.tsv')
hpa_rna_immune = read_data('data/rna_immune_cell.tsv')
hpa_rna_cell_line = read_data('data/rna_celline.tsv')

# Convert data to wide format
hpa_rna_single_cell_wide = reshape_wide(hpa_rna_single_cell, 'Cell.type')
hpa_rna_brain_wide = reshape_wide(hpa_rna_brain, 'Cluster.type')
hpa_rna_immune_wide = reshape_wide(hpa_rna_immune, 'Immune.cell')
hpa_rna_cell_line_wide = reshape_wide(hpa_rna_cell_line, 'Cell.line')

# Convert data to long format
hpa_rna_single_cell_long = reshape_long(hpa_rna_single_cell_wide)
hpa_rna_brain_long = reshape_long(hpa_rna_brain_wide)
hpa_rna_immune_long = reshape_long(hpa_rna_immune_wide)
hpa_rna_cell_line_long = reshape_long(hpa_rna_cell_line_wide)

# Overall visualization
# UMAP: all seemed normal 
# Single cell data
umap_single_results = compute_umap(hpa_rna_single_cell_long, 'Single cell data')
umap_single_xy = umap_single_results$umap_xy 
umap_single_results$umap_plot

# Brain data
umap_brain_results = compute_umap(hpa_rna_brain_long, 'Brain data')
umap_brain_xy = umap_brain_results$umap_xy 
umap_brain_results$umap_plot

# Immune data
umap_immune_results = compute_umap(hpa_rna_immune_long, 'Immune data')
umap_immune_xy = umap_immune_results$umap_xy 
umap_immune_results$umap_plot

# Cell line data 
umap_celline_results = compute_umap(hpa_rna_cell_line_long, 'Cell line data')
umap_celline_xy = umap_celline_results$umap_xy 
umap_celline_results$umap_plot

# PCA
# Single cell
pca_single_results = compute_pca(hpa_rna_single_cell_long, 'Single cell data')
pca_single_results$pca_plot
# Brain data
pca_brain_results = compute_pca(hpa_rna_brain_long, 'Brain data')
pca_brain_results$pca_plot
# Immune data
pca_immune_results = compute_pca(hpa_rna_immune_long, 'Immune data')
pca_immune_results$pca_plot  
# Cell line
pca_celline_results = compute_pca(hpa_rna_cell_line_long, 'Cell line data')
pca_celline_results$pca_plot  

# Always expressed analysis: find the list of genes (RNA) that found expressed in all samples
datasets = list(single = hpa_rna_single_cell_wide, brain = hpa_rna_brain_wide,
  immune = hpa_rna_immune_wide, cell_line = hpa_rna_cell_line_wide)

# Threshold
ntpm_thresholds = c(0, 0.2, 0.5, 1)
percentage_thresholds = c(100, 98, 95)
results_list = set_names(map(ntpm_thresholds, function(ntpm) {
    set_names(map(percentage_thresholds, function(perc) {
      map(datasets, ~ calculate_expression_percentage(.x, ntpm, perc))
    }), paste0('ntpm_', ntpm, '_perc_', percentage_thresholds))}), paste0('ntpm_', ntpm_thresholds))

# Generate summary tables and plots
plots = lapply(ntpm_thresholds, function(ntpm) {
  summary_table = map_dfr(names(datasets), function(dataset) {
    gene_counts = map_int(percentage_thresholds, ~ get_gene_counts(dataset, ntpm, .x))
    names(gene_counts) = c('T100', 'T98', 'T95')
    data.frame(Dataset = dataset, T100 = gene_counts['T100'], T98 = gene_counts['T98'],
               T95 = gene_counts['T95'], Total_Genes = nrow(datasets[[dataset]]))}) %>%
    mutate(Extra_100 = T100, Extra_98 = T98 - T100,
      Extra_95 = T95 - T98, Extra_Less95 = Total_Genes - T95) %>%
    select(Dataset, Extra_100, Extra_98, Extra_95, Extra_Less95, Total_Genes)
  # Transform for stacked bar plot
  stacked_data = summary_table %>%
    pivot_longer(cols = starts_with('Extra'), names_to = 'Threshold', values_to = 'Gene_Count') %>%
    mutate(Total_Genes = rep(summary_table$Total_Genes, each = 4),
           Threshold = factor(Threshold, levels = c('Extra_Less95', 'Extra_95', 'Extra_98', 'Extra_100'),
                              labels = c('<95%', '95%', '98%', '100%')), Percentage = (Gene_Count / Total_Genes) * 100)
  ntpm_label = ifelse(ntpm == 0, 'nTPM > 0', paste0('nTPM ≥ ', ntpm))
  ggplot(stacked_data, aes(x = Dataset, y = Gene_Count, fill = Threshold)) +
    geom_bar(stat = 'identity', color = 'black') +
    geom_text(aes(label = paste0(round(Percentage, 1), '%')),
              position = position_stack(vjust = 0.5), size = 4, color = 'white') +
    scale_fill_manual(name = 'Expression Threshold', 
                      values = c('100%' = 'darkblue', '98%' = 'blue', '95%' = 'lightblue', '<95%' = 'gray')) +
    labs(title = paste0('Always Expressed Genes Across Datasets (', ntpm_label, ')'),
         x = 'Dataset', y = 'Number of Genes') + theme_minimal() +
    theme(plot.title = element_text(hjust = 0.25, size = 12, face = 'bold'))})

# Display plots
invisible(lapply(plots, print))

# Extract always-expressed gene lists for each dataset 
# We used the strict threshold: nTPM >= 1 in 100%
always_expressed_single = results_list[['ntpm_1']][['ntpm_1_perc_100']][['single']]$always_expressed
always_expressed_brain = results_list[['ntpm_1']][['ntpm_1_perc_100']][['brain']]$always_expressed
always_expressed_immune = results_list[['ntpm_1']][['ntpm_1_perc_100']][['immune']]$always_expressed
always_expressed_celline = results_list[['ntpm_1']][['ntpm_1_perc_100']][['cell_line']]$always_expressed

# Overlaps of always expressed analysis
overlap_results = compute_gene_overlap(always_expressed_single, always_expressed_brain, 
                                       always_expressed_immune, always_expressed_celline, 
                                       dataset_names = c('Single_Cell', 'Brain', 'Immune', 'Cell_Line'))
overlap_summary = overlap_results$overlap_summary
shared_genes = overlap_results$shared_genes  

# Stacked bar plot
plot_gene_overlap(overlap_summary)

# Filter datasets for overlapping always-expressed genes
overlaps_single = filter_by_list(hpa_rna_single_cell_wide, shared_genes$Gene, 'Gene')
overlaps_brain = filter_by_list(hpa_rna_brain_wide, shared_genes$Gene, 'Gene')
overlaps_immune = filter_by_list(hpa_rna_immune_wide, shared_genes$Gene, 'Gene')
overlaps_celline = filter_by_list(hpa_rna_cell_line_wide, shared_genes$Gene, 'Gene')

# Prepare datasets for UMAP
overlap_always_expressed_single = prepare_for_umap(overlaps_single)
overlap_always_expressed_brain = prepare_for_umap(overlaps_brain)
overlap_always_expressed_immune = prepare_for_umap(overlaps_immune)
overlap_always_expressed_celline = prepare_for_umap(overlaps_celline)

# Run UMAP using the existing function
umap_single_overlaps = compute_umap(overlap_always_expressed_single, 'Single cell data - Overlap genes expression')
umap_brain_overlaps = compute_umap(overlap_always_expressed_brain, 'Brain data - Overlap genes expression')
umap_immune_overlaps = compute_umap(overlap_always_expressed_immune, 'Immune data - Overlap genes expression')
umap_celline_overlaps = compute_umap(overlap_always_expressed_celline, 'Cell line data - Overlap genes expression')

# Display UMAP plots
print(umap_single_overlaps$umap_plot)
print(umap_brain_overlaps$umap_plot)
print(umap_immune_overlaps$umap_plot)
print(umap_celline_overlaps$umap_plot)

# Only single-cell data: we observed 2 distinct clusters
# Extract UMAP coordinates from the overlap UMAP analysis
umap_overlaps_single_xy = umap_single_overlaps$umap_xy
umap_overlaps_single_xy$Gene = rownames(overlap_always_expressed_single)

# Split genes into two clusters based on UMAP2 threshold
cluster_high = umap_overlaps_single_xy %>% filter(UMAP2 > 10) %>% pull(Gene)
cluster_low = umap_overlaps_single_xy %>% filter(UMAP2 <= 10) %>% pull(Gene)

# Convert to Entrez IDs for each cluster
entrez_high = bitr(cluster_high, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
entrez_low = bitr(cluster_low, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

# Run GO enrichment analysis for each cluster
go_high = enrichGO(gene = entrez_high$ENTREZID, OrgDb = org.Hs.eg.db, 
                   keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = 'BH', 
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05)

go_low = enrichGO(gene = entrez_low$ENTREZID, OrgDb = org.Hs.eg.db, 
                  keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = 'BH', 
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)

# Generate enrichment dot plots for both clusters
dotplot(go_high, showCategory = 10, title = 'Top 10 GO terms - Single cell data - Alwasy expressed genes (1)')
dotplot(go_low, showCategory = 10, title = 'Top 10 GO erms - Single cell data - Alwasy expressed genes (2)')

























