# Mai Nguyen 2/2/2025
# Human housekeeping proteome (HKP) analysis 
# Function file

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

