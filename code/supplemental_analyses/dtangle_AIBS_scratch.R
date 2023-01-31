library(tidyverse)
library(Seurat)
library(edgeR)
library(dtangle)
library(limma)



### LOAD IN AND NORMALIZE CMC DATA - this is used for bulk tissue deconvolution

# get together the common mind metadata

# cmc has two types of metadata, 
# 1. clinical metadata, which is about the characteristics of each subject
# 2. and rnaseq metadata, which is more about the rnaseq sample quality and characteristics

# get clinical meta
cmc_clinical_path = '/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/Clinical/'
cmc_clinical_csv = paste0(cmc_clinical_path, 'CMC_Human_clinical_metadata.csv')

cmc_clinical_meta = read_csv(cmc_clinical_csv)
colnames(cmc_clinical_meta) = make.names(colnames(cmc_clinical_meta), unique = T)

# get rnaseq meta
cmc_rnaseq_path = '/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/RNAseq/Release4/Metadata/'
cmc_rnaseq_meta_path = paste0(cmc_rnaseq_path, 'CMC_Human_rnaSeq_metadata.csv')

cmc_rnaseq_meta = read_csv(cmc_rnaseq_meta_path)
colnames(cmc_rnaseq_meta) = make.names(colnames(cmc_rnaseq_meta), unique = T)

# create a new df called cmc_meta that merges the clinical and rnaseq meta - note that this duplicates data from some subjects
cmc_meta = full_join(cmc_clinical_meta, cmc_rnaseq_meta)

# parse and massage and format subject data fields 

cmc_meta$Dx = factor(cmc_meta$Dx, levels = c('Control', 'SCZ', 'BP', 'AFF', 'undetermined'))
cmc_meta$Reported.Gender = factor(cmc_meta$Reported.Gender, levels = c('Male', 'Female'))
cmc_meta$Ethnicity = factor(cmc_meta$Ethnicity, levels = c('Caucasian', 'African-American', 'Asian', 'Hispanic', '(Multiracial)'))

# convert ages to numbers
cmc_meta = cmc_meta %>% mutate(Age_norm = case_when(Age.of.Death == '90+' ~ 90, 
                                                    TRUE ~ as.numeric(Age.of.Death)))
# get numeric subject id - this is required for mapping to etienne's data from Pitt cohort
cmc_meta = cmc_meta %>% mutate(Subject = parse_number(SampleID))


### Now read in gene expression data from CommonMind subjects
gene_counts_ob = readRDS('/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/RNAseq/Release4/QuantitatedExpression/geneCountsMerged.RDS')

gene_counts_ob %>% dim

# this is the expression data from Pitt Penn and MSSM
cmc_expr_1 = read_tsv('/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/RNAseq/Release4/QuantitatedExpression/MSSM.Penn.Pitt_DLPFC.featureCount.tsv.gz')
# 
cmc_expr_2 = read_tsv('/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/RNAseq/Release4/QuantitatedExpression/NIMH.HBCC.featureCount.tsv.gz')
# 
cmc_expr = cbind(cmc_expr_1, cmc_expr_2[, 7:387]) 

# cmc_expr = cmc_expr_1 # just using data from main cmc dataset (no nimh-hbcc)

# cmc_expr_counts_only = cmc_expr[, 7:630] %>% as.data.frame()

cmc_expr_counts_only = cmc_expr[, 7:1011] %>% as.data.frame()
rownames(cmc_expr_counts_only) = cmc_expr$Geneid %>% substr(., 1, 15) %>% make.names(., unique = T)

### perform normalization for gene expression data from CMC

cmc_expr_cpm = edgeR::cpm(cmc_expr_counts_only, prior.count = 1, log = T)
cmc_expr_cpm = cmc_expr_cpm %>% as.data.frame() %>% tibble::rownames_to_column(var = "ensembl_id")



# just use samples from MSSM and remove samples that don't have valid sample ids (not sure why those samples don't)
use_samples = cmc_meta %>% filter(Institution == 'MSSM') %>% distinct(Individual.ID, .keep_all = T) %>% pull(SampleID)
use_samples = intersect(use_samples, colnames(cmc_expr_cpm))

cmc_expr_use = cmc_expr_cpm[, c('ensembl_id', use_samples)]

cmc_meta_use = cmc_meta %>% filter(SampleID %in% use_samples)

cmc_expr_cpm = cmc_expr_use

cmc_meta = cmc_meta_use

hgnc_mapping = read_tsv('~/mathys_analysis/hgnc_complete_set.txt')  %>% distinct(ensembl_gene_id, .keep_all = T)

# rownames(cmc_expr_counts_only) = cmc_expr$Geneid %>% substr(., 1, 15) %>% make.names(., unique = T)
# cmc_expr_cpm = edgeR::cpm(cmc_expr_counts_only, prior.count = 0.1, log = T)
# cmc_expr_cpm = cmc_expr_cpm %>% as.data.frame() %>% tibble::rownames_to_column(var = "ensembl_id")
rownames(cmc_expr_cpm) = cmc_expr_cpm$ensembl_id

hgnc_mapping %>% dplyr::select(ensembl_gene_id, symbol) %>% distinct(ensembl_gene_id, .keep_all = T)

matching_ids = intersect(hgnc_mapping$ensembl_gene_id, rownames(cmc_expr_cpm))
matching_gene_names = hgnc_mapping[hgnc_mapping$ensembl_gene_id %in% matching_ids, 'symbol']
cmc_expr_matching_gene_names = cmc_expr_cpm[matching_ids,  ]

rownames(cmc_expr_matching_gene_names) = make.names(matching_gene_names$symbol)

me = cmc_expr_matching_gene_names[-1]
# me is the final bulk expression dataset that you want to deconvolve


### now we'll load in the single cell gene expression data for generating markers and using for reference deconvolution

## try to load in aibs single cell seurat object (maybe the one used for aging paper?) and then use for dtangle

# sonny - this is the most recent Seurat object I could find
aibs_seurat_path = '/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_AIBS_obj_update_07JUN21.rds'
aibs_seurat = readRDS(aibs_seurat_path)

# pull metadata for aibs cells like we did in the aging paper - we'll likely need to change this!
aibs_small_metadata = aibs_seurat@meta.data %>% filter(region_label == 'CgG' | NeuN == 'Non-neuronal') 

# generate a new cell type label that matches roughly what we used in the aging paper
aibs_small_metadata = aibs_small_metadata %>% 
  mutate(aging_paper_label = subclass_label_collapsed)

# map use class label for glutamateric cells and call these Excitatory
aibs_small_metadata$aging_paper_label[aibs_small_metadata$class_label == 'Glutamatergic'] = 'Excitatory' 

# pull expression data for cells defined by metadata above
aibs_seurat_small = aibs_seurat@assays$RNA@counts[, aibs_small_metadata$sample_name]

# 
# write_rds(aibs_seurat_small, 'dtangle_test/aibs_seurat_small.rds')
# write_rds(aibs_small_metadata, 'dtangle_test/aibs_small_metadata.rds')
# 
# 
# aibs_seurat_small = read_rds('~/commonmind_analyses/dtangle_test/aibs_seurat_small.rds')
# aibs_small_metadata = read_rds('~/commonmind_analyses/dtangle_test/aibs_small_metadata.rds')

sce = aibs_seurat_small %>% cpm(prior.count = 1, log = T) # unclear what to use for prior count, but the output needs to be log transformed i think

### sonny - we need to talk about what set of genes to use for this - as this affects the marker gene identification done by dtangle below
# my vote is to use the harmonized list of genes that we ended up using in the bulk analysis (e.g., genes available across all/most datasets)

commongenes <- intersect (rownames(me), rownames(sce))

sce <- sce[pmatch(commongenes, rownames(sce)), ]

# me is defined above by the target bulk dataset we want to deconvlve
me <- me[pmatch(commongenes, rownames(me)), ]


y <- cbind(sce, me)
y <- normalizeBetweenArrays(y) # this step takes a while given how large the single cell dataset is
y <- t(y)

all_cell_type <- unique(aibs_small_metadata$aging_paper_label)

pure_samples <- lapply(1:length(all_cell_type), function(i) {
  which(aibs_small_metadata$aging_paper_label == all_cell_type[i])
})

names(pure_samples) = all_cell_type

# you don't actually need bulk data here, you find markers with just the single-cell reference data
marker_list = find_markers(y,pure_samples=pure_samples,data_type="rna-seq",marker_method='ratio')


num_use_markers = 50 # number of top markers to use per cell type
n_markers = sapply(1:K,function(i){num_use_markers})

# you can also have dtangle pick how many markers per cell type to use for you
# q = .1
# quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
# K = length(pure_samples)
# n_markers = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})
# n_markers
                   

# based on the expresion dataset you're trying to deconvolve, you'd change data type to rna-microarray i believe
marks = marker_list$L
dc <- dtangle(y, pure_samples=pure_samples, n_markers=n_markers, data_type = 'rna-seq', markers = marks)

final_est <- dc$estimates[(dim(sce)[2]+1):dim(y)[1],]
colnames(final_est) <-  all_cell_type

# final_est has the markers 

head(final_est)

write.csv(final_est, 'code/collab/data/mssm_dtangle_estimates.csv')
