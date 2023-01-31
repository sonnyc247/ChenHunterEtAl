library(tidyverse)
library(magrittr)
library(cowplot)
library(edgeR)
library(markerGeneProfile)
library(readxl)


# get sample metadata from gtex website
gtex_samples = read_delim(url("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), delim = '\t')

gtex_data_path = '/external/rprshnas01/netdata_kcni/stlab/Public/GTEx/'

# read in full junction dataset into workspace - careful, this file is fucking huge (14 GB)
reads_data_file_name = 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct'
gtex_reads = read_tsv(paste0(gtex_data_path, reads_data_file_name), skip=2)


# get samples from brain regions from cortex

use_cortex_regions = c("Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)")

use_samples = gtex_samples %>% filter(SMTSD %in% use_cortex_regions) %>% pull(SAMPID)

# add a column called SUBJID to metadata that enables linking with subject metadata
gtex_samples = gtex_samples %>% 
  separate(SAMPID, into = c('first', 'second'), remove = F, sep = '-') %>% 
  unite(col = "SUBJID", first:second, sep = '-', remove = F) 

get_sample_col_names = intersect(c('Name', 'Description', use_samples), colnames(gtex_reads))
brain_gtex_reads = gtex_reads[,get_sample_col_names]

brain_gtex_metadata = gtex_samples %>% filter(SAMPID %in% use_samples)


# subject metadata
gtex_subject_meta = read_delim(url('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt'), delim = '\t')


### perform normalization for gene expression data from GTEX

gtex_brain_expr_cpm = edgeR::cpm(brain_gtex_reads[, 3:387], prior.count = 1)
rownames(gtex_brain_expr_cpm) = str_extract(brain_gtex_reads$Name, 'ENSG\\d+')
gtex_brain_expr_cpm = gtex_brain_expr_cpm %>% as.data.frame() %>% tibble::rownames_to_column(var = "ensembl_id")


### GET MARKERS FOR MGP ANALYSIS
# update this to get sonnys markers from emmas paper
sonny_markers = read_csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'))
colnames(sonny_markers) = colnames(sonny_markers) %>% make.names() %>% tolower()

# I find it helpful to map some gene symbols to ensembl ids manually using mappings from hgnc, you can get those from here: 
# http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt

hgnc_mapping = read_tsv('~/mathys_analysis/hgnc_complete_set.txt')

# now, this is the list of sonnys markers with entrez ids and ensembl ids where possible
sonny_hgnc_merged_markers = left_join(sonny_markers %>% dplyr::rename(entrez_id = entrez.gene.id), 
                                      hgnc_mapping %>% distinct(entrez_id, .keep_all = T)%>% 
                                        dplyr::select(entrez_id, ensembl_gene_id) %>% 
                                        dplyr::rename(ensembl_id = ensembl_gene_id)) %>% 
  dplyr::select(x1, gene, entrez_id, ensembl_id, -ensembl.gene.id, everything()) %>% 
  group_by(subclass) %>% 
  arrange(subclass, -average.log.fold.change) %>% 
  ungroup()

# get ensembl list of markers
new_markers = sonny_hgnc_merged_markers %>% filter(used.in.mgp == "TRUE")
new_cell_types = new_markers %>% filter(!is.na(subclass)) %>% pull(subclass) %>% unique
new_marker_list_ensembl <- lapply(new_cell_types, function(cell_type){
  return(new_markers %>% filter(subclass == cell_type, 
                                ensembl_id %in% gtex_brain_expr_cpm$ensembl_id,
  ) %>% pull(ensembl_id))
})
names(new_marker_list_ensembl) <- new_cell_types
print(new_cell_types)

# new_marker_list_ensembl reflects the markers that we'll use in the MGP calculation


#### PERFORM MGP ANALYSIS

# just use samples from Pitt and remove samples that don't have valid sample ids (not sure why those samples don't)

use_samples = intersect(use_samples, colnames(gtex_brain_expr_cpm))


# perform MGP calculations, I set seekConsensus and removeMinority to FALSE as we're using Micaela's consensus markers
estimations <-  mgpEstimate(
  exprData= gtex_brain_expr_cpm[, c('ensembl_id', use_samples)], # filters samples for those defined in use_samples
  genes=new_marker_list_ensembl,
  geneColName='ensembl_id',
  outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
  geneTransform = NULL, # this is the default option for geneTransform
  groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus=FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority=FALSE) 

# rm(cmc_expr_1, cmc_expr_cpm, cmc_expr, cmc_expr_counts_only, gene_counts_ob)

# calculate qc metrics as well
pc1_exp = lapply(new_cell_types, function(cell_type){
  s = estimations$trimmedPCAs[[cell_type]] %>% summary()
  return(s$importance[3, 1])
}) %>% unlist()

mgp_qc_df = data.frame(subclass = new_cell_types, removedMarkerRatios = estimations$removedMarkerRatios, pc1_exp = pc1_exp )


mgp_df_big = estimations$estimates %>% as.data.frame() %>% tibble::rownames_to_column(var = "SAMPID") 

# merge mgp df with metadata df
brain_gtex_mgps_wide = left_join(mgp_df_big , brain_gtex_metadata) # this is the final mgp df 

brain_gtex_mgps_wide = left_join(brain_gtex_mgps_wide , gtex_subject_meta) # this is the final mgp df 

# make a plot
brain_gtex_mgps_wide %>% ggplot(aes(x = AGE, y = SST)) + 
  geom_boxplot(outlier.size = 0) + 
  geom_jitter(alpha = .5) + 
  theme_cowplot() + 
  facet_wrap(~SMTSD) 
