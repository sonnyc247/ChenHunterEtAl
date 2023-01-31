# create snCTP data frames from stored metadata on scc

library(tidyverse)
library(data.table)
library(cowplot)


# load in cell metadata from single cell objects from mt sinai and maclean snRNAseq dataset
actionnet_summary = readRDS('/external/rprshnas01/netdata_kcni/stlab/SZBDMulticohort/ACTIONet_summary.rds')

# medadata from individuals from dataset
indiv_meta = read_tsv('/external/rprshnas01/netdata_kcni/stlab/SZBDMulticohort/individual_metadata.tsv')

# redefine cell metadata object
cell_meta = actionnet_summary$metadata


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



### calculate snCTPs using the harmonized cell type labels that match to those used in the aging paper
# now we count up cells per type per individual and calculate proportions
total_cell_counts_per_cell_type_per_indiv = cell_meta %>% 
  group_by(ID, aging_paper_label, .drop = FALSE) %>% 
  tally() %>%
  ungroup %>%
  complete(ID, aging_paper_label, fill = list(n = 0))

total_cell_counts_per_indiv = cell_meta %>% group_by(ID, .drop = FALSE) %>% tally(name = 'total_counts')

cell_prop_df = left_join(total_cell_counts_per_cell_type_per_indiv, total_cell_counts_per_indiv) %>% 
  mutate(cell_prop = n / total_counts)



indiv_meta_long = left_join(indiv_meta, cell_prop_df)

indiv_meta_wide = indiv_meta_long %>% pivot_wider(id_cols = ID, names_from = "aging_paper_label", values_from = "cell_prop")

indiv_meta_wide = left_join(indiv_meta, indiv_meta_wide)

indiv_meta_wide = left_join(indiv_meta_wide, cell_prop_df %>% select(ID, total_counts))

# indiv_meta_wide = indiv_meta_wide %>% mutate(CMC_ID_full = paste0('CMC_', CMC_ID))

scale_this <- function(x) as.vector(scale(x))


indiv_meta_wide_output = indiv_meta_wide %>% 
  filter(Phenotype == 'CON') %>% 
  dplyr::rename(dataset_subject_id = ID, 
                           ageOfDeath = Age, 
                           PMI = PMI) %>%
  mutate(Sex = case_when(Gender == 'Male' ~ 'M',
                         Gender == 'Female' ~ 'F',
                         TRUE ~ 'Other'), 
         Institution = case_when(Cohort == 'MtSinai' ~ 'MSSM',
                                 Cohort == 'McLean' ~ 'McLean',
                                  TRUE ~ 'Other'), 
         dataset = case_when(Cohort == 'MtSinai' ~ 'MtSinai_snRNAseq',
                             Cohort == 'McLean' ~ 'McLean_snRNAseq',
                                 TRUE ~ 'Other'), 
         bulk_sample_id = paste0('MSSM_RNA_PFC_', parse_number(CMC_ID)),
         subject_id = case_when(Cohort == 'MtSinai' ~ bulk_sample_id,
                                  Cohort == 'McLean' ~ dataset_subject_id,
                                  TRUE ~ 'Other'), 
         total_cell_counts = total_counts
         
  ) %>% 
  select(dataset_subject_id, Astrocyte:VIP, ageOfDeath, Sex, PMI, Institution, dataset, subject_id, bulk_sample_id, total_cell_counts) %>%
  group_by(Institution) %>%
  mutate(age_scaled = scale(ageOfDeath, center = T, scale = T) %>% as.vector(), 
         pmi_scaled = scale(PMI, center = T, scale = T) %>% as.vector()) %>% 
  ungroup() %>% 
  distinct(subject_id, .keep_all = T) %>%
  select(subject_id, everything())

write.csv(indiv_meta_wide_output, 'code/collab/data/snCTP_results_harmonized.csv')


### calculate snCTPs using the original cell type labels used in the SZBDMulticohort dataset

# now we count up cells per type per individual and calculate proportions
total_cell_counts_per_cell_type_per_indiv = cell_meta %>% 
  group_by(ID, Celltype, .drop = FALSE) %>% 
  tally() %>%
  ungroup %>%
  complete(ID, Celltype, fill = list(n = 0))

total_cell_counts_per_indiv = cell_meta %>% group_by(ID, .drop = FALSE) %>% tally(name = 'total_counts')

cell_prop_df = left_join(total_cell_counts_per_cell_type_per_indiv, total_cell_counts_per_indiv) %>% 
  mutate(cell_prop = n / total_counts)



indiv_meta_long = left_join(indiv_meta, cell_prop_df)

indiv_meta_wide = indiv_meta_long %>% pivot_wider(id_cols = ID, names_from = "Celltype", values_from = "cell_prop")

indiv_meta_wide = left_join(indiv_meta, indiv_meta_wide)

indiv_meta_wide = left_join(indiv_meta_wide, cell_prop_df %>% select(ID, total_counts))

# indiv_meta_wide = indiv_meta_wide %>% mutate(CMC_ID_full = paste0('CMC_', CMC_ID))

scale_this <- function(x) as.vector(scale(x))


indiv_meta_wide_output = indiv_meta_wide %>% 
  filter(Phenotype == 'CON') %>% 
  dplyr::rename(dataset_subject_id = ID, 
                ageOfDeath = Age, 
                PMI = PMI) %>%
  mutate(Sex = case_when(Gender == 'Male' ~ 'M',
                         Gender == 'Female' ~ 'F',
                         TRUE ~ 'Other'), 
         Institution = case_when(Cohort == 'MtSinai' ~ 'MSSM',
                                 Cohort == 'McLean' ~ 'McLean',
                                 TRUE ~ 'Other'), 
         dataset = case_when(Cohort == 'MtSinai' ~ 'MtSinai_snRNAseq',
                             Cohort == 'McLean' ~ 'McLean_snRNAseq',
                             TRUE ~ 'Other'), 
         bulk_sample_id = paste0('MSSM_RNA_PFC_', parse_number(CMC_ID)),
         subject_id = case_when(Cohort == 'MtSinai' ~ bulk_sample_id,
                                Cohort == 'McLean' ~ dataset_subject_id,
                                TRUE ~ 'Other'), 
         total_cell_counts = total_counts
         
  ) %>% 
  select(dataset_subject_id, `Ex-NRGN`:`Pericytes`, ageOfDeath, Sex, PMI, Institution, dataset, subject_id, bulk_sample_id, total_cell_counts) %>%
  group_by(Institution) %>%
  mutate(age_scaled = scale(ageOfDeath, center = T, scale = T) %>% as.vector(), 
         pmi_scaled = scale(PMI, center = T, scale = T) %>% as.vector()) %>% 
  ungroup() %>%
  distinct(subject_id, .keep_all = T) %>%
  select(subject_id, everything())

write.csv(indiv_meta_wide_output, 'code/collab/data/snCTP_results_fine_grained.csv')


