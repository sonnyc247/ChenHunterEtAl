# this code figures out the cell type mappings for the SZBD multicohort data

# this dataset uses different cell type labels, so we need to figure out how to map them onto the labels we used for the aging paper

## load in pseudobulk data from stored object
scz_snrnaseq = readRDS('/external/rprshnas01/netdata_kcni/stlab/SZBDMulticohort/pseudobulk_mean_logcounts.rds')

# rename cell type labels so they're a bit prettier
names(scz_snrnaseq@assays@data) = names(scz_snrnaseq@assays@data) %>% str_remove(., '-')

# calculate means per each cell type in external snRNAseq data
scz_data_cell_type_means = lapply(names(scz_snrnaseq@assays@data), function(cell_type_name){
  print(cell_type_name)
  temp_df = scz_snrnaseq@assays@data[[cell_type_name]]
  
  ret_df = rowMeans(temp_df) %>% as.data.frame()
  colnames(ret_df) = cell_type_name
  return(ret_df)
}) %>% bind_cols()


# sce is the allen single cell object from before
# calculate means per each cell type in external snRNAseq data
aibs_cell_type_averages = lapply(all_cell_type, function(cell_type_name){
  cell_sample_ids = aibs_small_metadata %>% filter(aging_paper_label == cell_type_name) %>% pull(sample_name)
  cell_averages = rowMeans(sce[, cell_sample_ids]) %>% as.data.frame()
  colnames(cell_averages) = cell_type_name
  return(cell_averages)
}) %>% bind_cols()


# find common genes between allen and external snRNAseq dataset
commongenes = intersect(rownames(aibs_cell_type_averages), rownames(scz_data_cell_type_means))

# this calculates the correlation matrix using all genes (could probably use variable genes)
cor_mat = cor(aibs_cell_type_averages[commongenes, ], scz_data_cell_type_means[commongenes, ])

# from that correlation matrix, i looked at it manually and then came up with the mappings shown below

# cell meta now has a field called aging_paper_label which matches the cell types we used in the aging paper
cell_meta = cell_meta %>% mutate(aging_paper_label = 
                                   case_when(startsWith(as.character(Celltype), 'Ex') ~ 'Excitatory',
                                             startsWith(as.character(Celltype), 'In-Rosehip') ~ 'LAMP5',
                                             Celltype %in% c('In-Reelin', 'In-VIP') ~ 'VIP', 
                                             startsWith(as.character(Celltype), 'In-PV') ~ 'PVALB',
                                             Celltype == 'In-SST' ~ 'SST', 
                                             Celltype == 'Oli' ~ 'Oligodendrocyte', 
                                             Celltype == 'OPC' ~ 'OPC',
                                             Celltype == 'Ast' ~ 'Astrocyte',
                                             Celltype == 'Mic' ~ 'Microglia',
                                             Celltype == 'Endo' ~ 'Endothelial',
                                             Celltype == 'Pericytes' ~ 'Pericyte',
                                             TRUE ~ 'Other')
)