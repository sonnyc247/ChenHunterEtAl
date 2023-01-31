#### Packages ####

library(dplyr)
library(markerGeneProfile)
library(ggplot2)
library(ggpubr)

#### Marker list importing and processing ####

markers_df <- read.csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/All_hodge_regions/new_ALLReg_results_SubclassWithExc.csv'))

cell_types = markers_df$subclass %>% unique()
marker_list = lapply(cell_types, function(cell_type){
  return(markers_df %>% filter(subclass == cell_type) %>% pull(gene) %>% unlist())
})
names(marker_list) = cell_types

### Run the code for each dataset once before continuing with the code below, to get a list of which markers ends up being used/are available in each dataset
## Use "removeMinority = T" for mgpEstimate, store used markers in markers_used_total_df

### Finalizing markers (after first run)

#markers_used_total_df_Init <- markers_used_total_df
marker_temp <- as.data.frame(table(markers_used_total_df$rowname))
marker_temp <- marker_temp[marker_temp$Freq >= 7,]

marker_temp <- marker_temp[marker_temp$Var1 %in% CMC_genes,]
marker_temp <- marker_temp[marker_temp$Var1 %in% NABEC_genes,]
marker_temp <- marker_temp[marker_temp$Var1 %in% UKBEC_genes,]
marker_temp <- marker_temp[marker_temp$Var1 %in% BA11_genes,]
marker_temp <- marker_temp[marker_temp$Var1 %in% BA47_genes,]

markers_df <- markers_df[markers_df$gene %in% marker_temp$Var1,]

cell_types = markers_df$subclass %>% unique()
marker_list = lapply(cell_types, function(cell_type){
  return(markers_df %>% filter(subclass == cell_type) %>% pull(gene) %>% unlist())
})
names(marker_list) = cell_types

remove(cell_types)
remove(markers_df)
remove(marker_temp)

### annotate and export reader-friendly marker list

markers_df_original <- read.csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/All_hodge_regions/new_ALLReg_results_SubclassWithExc.csv'))
markers_df_original$used <- markers_df_original$gene %in% unname(unlist(marker_list))
markers_df_original <- markers_df_original[,-1]
markers_df_original <- markers_df_original[,-8]
names(markers_df_original) <- c("Gene symbol",
                                "Entrez gene ID",
                                "Ensembl gene ID",
                                "Cell type",
                                "Prevalence % in cell type ('pct.1')",
                                "Prevalence % in all other cells ('pct.2')",
                                "Average log 2 fold change",
                                "AUC from ROC test",
                                "Power from ROC test",
                                "Mast test p",
                                "Mast test FDR",
                                "Used in MGP")

write.csv(markers_df_original, "./data/output/tables/Supplementary_Table_1_AllRegions.csv")

#
#### CMC (RNAseq data) ####

### select subjects
METADATA <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/cmc_data/METADATA.rds")
METADATA <- METADATA[METADATA$Dx == "Control",]
METADATA <- METADATA[METADATA$ageOfDeath >= 15,]
METADATA$is_duplicated <- duplicated(METADATA$IndividualID) # ID replicates
replicateDF <- METADATA[METADATA$IndividualID %in% METADATA[METADATA$is_duplicated == T, "IndividualID"], c("SampleID", "IndividualID")]

### count matrix preprocessing
geneCountsMerged <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/cmc_data/geneCountsMerged.rds")
geneCountsMerged <- geneCountsMerged[,METADATA$SampleID]
geneCountsMerged <- geneCountsMerged[rowSums(geneCountsMerged) >= 15, ]

# sum replicates and adjust metadata

for (subject in unique(replicateDF$IndividualID)) {
  
  geneCountsMerged$temp <- rowSums(geneCountsMerged[, METADATA[METADATA$IndividualID == subject, "SampleID"]])
  colnames(geneCountsMerged)[ncol(geneCountsMerged)] <- paste0("PittRNAseq_", readr::parse_number(METADATA[METADATA$IndividualID == subject, "SampleID"])[1])
  
  newidx <- nrow(METADATA)+1
  METADATA[newidx,] <- METADATA[METADATA$is_duplicated == T & METADATA$IndividualID == subject,]
  METADATA[newidx,"SampleID"] <- paste0("PittRNAseq_", readr::parse_number(METADATA[METADATA$IndividualID == subject, "SampleID"])[1])
  
}

row.names(METADATA)[462:471] <- METADATA$SampleID[462:471]
METADATA <- METADATA[!(row.names(METADATA) %in% replicateDF$SampleID),]
geneCountsMerged <- geneCountsMerged[, !names(geneCountsMerged) %in% replicateDF$SampleID]
geneCountsMerged <- edgeR::cpm(geneCountsMerged, log = TRUE, prior.count = 1) # prior count 1!!!

### ensembl id mapping
#hgnc_mapping <- readr::read_tsv("./data/hgnc_complete_set_2021-06-01.txt") # retrieved from 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2021-06-01.txt'
hgnc_mapping <- read.delim("/external/rprshnas01/kcni/ychen/git/Aging_MGP_private/data/hgnc_complete_set_2021-06-01.txt")
hgnc_mapping <- hgnc_mapping[,c("symbol", "ensembl_gene_id")]

geneCountsMerged <- tibble::rownames_to_column(as.data.frame(geneCountsMerged))
geneCountsMerged$rowname <- gsub("\\..*","",geneCountsMerged$rowname) # remove version id from ensembl IDs
geneCountsMerged <- merge(hgnc_mapping, geneCountsMerged, by.x = "ensembl_gene_id", by.y = "rowname", all = F) # ENSG00000254876 and ENSG00000230417 duplicated; none of their gene symbols are markers
row.names(geneCountsMerged) <- make.names(geneCountsMerged$symbol)
geneCountsMerged <- geneCountsMerged[,3:ncol(geneCountsMerged)]

### mgp run

mgp_results_total_df <- temp_estimates_df[0,] # run one instance of the loop below (comment out appropriate rbind line), then run this
markers_used_total_df <- marker_accounting_df[0,] # first run only; run one instance of the loop below (comment out appropriate rbind line), then run this; to compile list of markers used

for (institution in unique(METADATA$Institution)){
  
  ## filter
  temp_metadata <- METADATA[METADATA$Institution == institution,]
  temp_genecounts <- geneCountsMerged[,colnames(geneCountsMerged) %in% temp_metadata$SampleID]
  temp_genecounts <- tibble::rownames_to_column(as.data.frame(temp_genecounts))
  
  ## run
  temp_estimates <- mgpEstimate(temp_genecounts, 
                               genes = marker_list, 
                               geneColName = 'rowname', 
                               outlierSampleRemove = F,
                               geneTransform = NULL,
                               groups = NULL, 
                               seekConsensus = F,
                               removeMinority = F) # T on first run
  
  ## process and scale
  temp_estimates_df <- temp_estimates$estimates %>% as.data.frame()
  temp_estimates_df <- scale(temp_estimates_df) %>% as.data.frame()
  temp_estimates_df <- tibble::rownames_to_column(temp_estimates_df)
  
  ## append
  metadata_extract <- temp_metadata[,c(1,8,6,9,3)]
  temp_estimates_df <- merge(temp_estimates_df, metadata_extract, by.x = "rowname", by.y = "SampleID")
  mgp_results_total_df <- rbind(mgp_results_total_df, temp_estimates_df)
  
  ## check markers
  marker_accounting_df <- tibble::rownames_to_column(bind_rows(temp_estimates$usedMarkerExpression, .id = "column_label"))
  print(institution)
  print(table(marker_accounting_df$column_label))
  marker_accounting_df$Dataset <- institution
  marker_accounting_df <- marker_accounting_df[,c("rowname", "column_label", "Dataset")]
  #markers_used_total_df <- rbind(markers_used_total_df, marker_accounting_df) ## as appropriate (use on first run of analysis to get used markers)

}

#CMC_genes <- temp_genecounts$rowname #first run only

remove(geneCountsMerged)
remove(hgnc_mapping)
remove(marker_accounting_df)
remove(METADATA)
remove(temp_genecounts)
remove(temp_metadata)
remove(institution)

#
#### NABEC (GSE36192; already log-transformed Microarray data) #### 

### preprocessing

NABEC_microarray <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/BEC_microarray/GSE36192_out.rds")
NABEC_microarray <- NABEC_microarray[NABEC_microarray$region == "frontal cortex",]
NABEC_microarray <- NABEC_microarray[NABEC_microarray$age_years >= 15,] 

## get and process gene expression data only
NABEC_micro_dataonly <- NABEC_microarray[,c(1,6:19618)] 
row.names(NABEC_micro_dataonly) <- NABEC_micro_dataonly$geo_accession
NABEC_micro_dataonly <- NABEC_micro_dataonly[,-1]
NABEC_micro_dataonly <- as.data.frame(t(NABEC_micro_dataonly))

## quartile filtering
NABEC_micro_dataonly$gene_median_exp <- NABEC_micro_dataonly %>% as.matrix() %>% matrixStats::rowMedians()
quantile_filter <- quantile(NABEC_micro_dataonly$gene_median_exp)[2] %>% unname()
NABEC_micro_filtereddata <- filter(NABEC_micro_dataonly, gene_median_exp >= quantile_filter)
NABEC_micro_filtereddata <- NABEC_micro_filtereddata[,-430]
NABEC_micro_filtereddata <- tibble::rownames_to_column(NABEC_micro_filtereddata)

### Main mgp run
temp_estimates <- mgpEstimate(NABEC_micro_filtereddata, 
                              genes = marker_list, 
                              geneColName = 'rowname', 
                              outlierSampleRemove = F,
                              geneTransform = NULL,
                              groups = NULL, 
                              seekConsensus = F,
                              removeMinority = F)

## process
temp_estimates_df <- temp_estimates$estimates %>% as.data.frame()
temp_estimates_df <- scale(temp_estimates_df) %>% as.data.frame()
temp_estimates_df <- tibble::rownames_to_column(temp_estimates_df)

## check markers
marker_accounting_df <- tibble::rownames_to_column(bind_rows(temp_estimates$usedMarkerExpression, .id = "column_label"))
table(marker_accounting_df$column_label) 
marker_accounting_df$Dataset <- "NABEC"
marker_accounting_df <- marker_accounting_df[,c("rowname", "column_label", "Dataset")]
markers_used_total_df <- rbind(markers_used_total_df, marker_accounting_df) ## as appropriate (use on first run of analysis to get used markers)


## append
metadata_extract <- NABEC_microarray[,c(1:3,5)]
metadata_extract$Institution <- "NABEC"
temp_estimates_df <- merge(temp_estimates_df, metadata_extract, by.x = "rowname", by.y = "geo_accession")
colnames(temp_estimates_df)[15:18] <- colnames(mgp_results_total_df)[15:18]
mgp_results_total_df <- rbind(mgp_results_total_df, temp_estimates_df)

#NABEC_genes <- NABEC_micro_filtereddata$rowname # first run only

remove(NABEC_micro_dataonly)
remove(NABEC_micro_filtereddata)
remove(NABEC_microarray) 

#
#### UKBEC (GSE60862; already log-transformed Microarray data) #### 

### preprocessing

UKBEC_microarray <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/BEC_microarray/GSE60862_out.rds")
UKBEC_microarray <- UKBEC_microarray[UKBEC_microarray$region == "frontal cortex",]

## get and process gene expression only
UKBEC_micro_dataonly <- UKBEC_microarray[,c(1,9:20820)]
row.names(UKBEC_micro_dataonly) <- UKBEC_micro_dataonly$geo_accession
UKBEC_micro_dataonly <- UKBEC_micro_dataonly[,-1]
UKBEC_micro_dataonly <- as.data.frame(t(UKBEC_micro_dataonly))

## quartile filtering
UKBEC_micro_dataonly$gene_median_exp <- UKBEC_micro_dataonly %>% as.matrix() %>% matrixStats::rowMedians()
quantile_filter <- quantile(UKBEC_micro_dataonly$gene_median_exp)[2] %>% unname()
UKBEC_micro_filtereddata <- filter(UKBEC_micro_dataonly, gene_median_exp >= quantile_filter)
UKBEC_micro_filtereddata <- UKBEC_micro_filtereddata[,-128]
UKBEC_micro_filtereddata <- tibble::rownames_to_column(UKBEC_micro_filtereddata)

### Main run
temp_estimates <- mgpEstimate(UKBEC_micro_filtereddata, 
                              genes = marker_list, 
                              geneColName = 'rowname', 
                              outlierSampleRemove = F,
                              geneTransform = NULL,
                              groups = NULL, 
                              seekConsensus = F,
                              removeMinority = F) # T on first run

## process
temp_estimates_df <- temp_estimates$estimates %>% as.data.frame()
temp_estimates_df <- scale(temp_estimates_df) %>% as.data.frame()
temp_estimates_df <- tibble::rownames_to_column(temp_estimates_df)

## check markers
marker_accounting_df <- tibble::rownames_to_column(bind_rows(temp_estimates$usedMarkerExpression, .id = "column_label"))
table(marker_accounting_df$column_label) 
marker_accounting_df$Dataset <- "UKBEC"
marker_accounting_df <- marker_accounting_df[,c("rowname", "column_label", "Dataset")]
markers_used_total_df <- rbind(markers_used_total_df, marker_accounting_df) ## as appropriate (use on first run of analysis to get used markers)

## append
metadata_extract <- UKBEC_microarray[,c(1:3,6)]
metadata_extract$Institution <- "UKBEC"
temp_estimates_df <- merge(temp_estimates_df, metadata_extract, by.x = "rowname", by.y = "geo_accession")
colnames(temp_estimates_df)[15:18] <- colnames(mgp_results_total_df)[15:18]
mgp_results_total_df <- rbind(mgp_results_total_df, temp_estimates_df)

#UKBEC_genes <- UKBEC_micro_filtereddata$rowname # first run only

remove(UKBEC_micro_dataonly)
remove(UKBEC_micro_filtereddata)
remove(UKBEC_microarray) 
remove(quantile_filter)

#
#### BA11 (already quartile-filtered Microarray data) ####

### Pre processing
BA11_expr_filtered <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/pitt_microarray/BA11_expr_filtered.rds")
row.names(BA11_expr_filtered) <- BA11_expr_filtered$Gene_Symbol
BA11_expr_filtered <- BA11_expr_filtered[,-1]
BA11_expr_filtered <- log2(BA11_expr_filtered)
BA11_expr_filtered <- tibble::rownames_to_column(BA11_expr_filtered)

### Main run
temp_estimates <- mgpEstimate(BA11_expr_filtered, 
                              genes = marker_list, 
                              geneColName = 'rowname', 
                              outlierSampleRemove = F,
                              geneTransform = NULL,
                              groups = NULL, 
                              seekConsensus = F,
                              removeMinority = F) # T for first run

## process
temp_estimates_df <- temp_estimates$estimates %>% as.data.frame()
temp_estimates_df <- scale(temp_estimates_df) %>% as.data.frame()
temp_estimates_df <- tibble::rownames_to_column(temp_estimates_df)

## check markers
marker_accounting_df <- tibble::rownames_to_column(bind_rows(temp_estimates$usedMarkerExpression, .id = "column_label"))
table(marker_accounting_df$column_label) 
marker_accounting_df$Dataset <- "BA11"
marker_accounting_df <- marker_accounting_df[,c("rowname", "column_label", "Dataset")]
markers_used_total_df <- rbind(markers_used_total_df, marker_accounting_df) ## as appropriate (use on first run of analysis to get used markers)

## append
sample_meta <- read.csv('/external/rprshnas01/netdata_kcni/stlab/bulk_expression/pittsburgh_aging_cohort/pitt_microarray_cleaned.csv') 
sample_meta$region_id <- paste0(sample_meta$ID, "_", sample_meta$Region)

metadata_extract <- sample_meta[,c(9,7,5,2)]
metadata_extract$Institution <- "Pitt"
temp_estimates_df <- merge(temp_estimates_df, metadata_extract, by.x = "rowname", by.y = "region_id")
colnames(temp_estimates_df)[15:18] <- colnames(mgp_results_total_df)[15:18]
mgp_results_total_df <- rbind(mgp_results_total_df, temp_estimates_df)

#BA11_genes <- BA11_expr_filtered$rowname # first run only

## from here on, add new dataset column, which is needed to distinguish pitt RNAseq and microarray data

mgp_results_total_df$dataset <- mgp_results_total_df$Institution
mgp_results_total_df[endsWith(mgp_results_total_df$rowname, "_BA11"),"dataset"] <- "BA11"

remove(BA11_expr_filtered)

#
#### BA47 (already quartile-filtered Microarray data) ####

### Pre processing
BA47_expr_filtered <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/pitt_microarray/BA47_expr_filtered.rds")
row.names(BA47_expr_filtered) <- BA47_expr_filtered$Gene_Symbol
BA47_expr_filtered <- BA47_expr_filtered[,-1]
BA47_expr_filtered <- log2(BA47_expr_filtered) 
BA47_expr_filtered <- tibble::rownames_to_column(BA47_expr_filtered)

### Main run
temp_estimates <- mgpEstimate(BA47_expr_filtered, 
                              genes = marker_list, 
                              geneColName = 'rowname', 
                              outlierSampleRemove = F,
                              geneTransform = NULL,
                              groups = NULL, 
                              seekConsensus = F,
                              removeMinority = F) # T for first run

## process
temp_estimates_df <- temp_estimates$estimates %>% as.data.frame()
temp_estimates_df <- scale(temp_estimates_df) %>% as.data.frame()
temp_estimates_df <- tibble::rownames_to_column(temp_estimates_df)

## check markers
marker_accounting_df <- tibble::rownames_to_column(bind_rows(temp_estimates$usedMarkerExpression, .id = "column_label"))
table(marker_accounting_df$column_label) 
marker_accounting_df$Dataset <- "BA47"
marker_accounting_df <- marker_accounting_df[,c("rowname", "column_label", "Dataset")]
markers_used_total_df <- rbind(markers_used_total_df, marker_accounting_df) ## as appropriate (use on first run of analysis to get used markers)

## append (Using same metadata_extract as produced for BA11)

temp_estimates_df <- merge(temp_estimates_df, metadata_extract, by.x = "rowname", by.y = "region_id")
colnames(temp_estimates_df)[15:18] <- colnames(mgp_results_total_df)[15:18]
temp_estimates_df$dataset <- "BA47"
mgp_results_total_df <- rbind(mgp_results_total_df, temp_estimates_df)

#BA47_genes <- BA47_expr_filtered$rowname # first run only

remove(BA47_expr_filtered)
remove(marker_accounting_df)
remove(metadata_extract)
remove(sample_meta)
remove(temp_estimates_df)
remove(temp_estimates)

#
#### Other post-processing ####

### Scaling
mgp_results_total_df$age_scaled <- scale(mgp_results_total_df$ageOfDeath)
mgp_results_total_df$pmi_scaled <- scale(mgp_results_total_df$PMI)

### Sex variable harmonization

unique(mgp_results_total_df$Sex)
mgp_results_total_df[mgp_results_total_df$Sex == "Male", "Sex"] <- "M"
mgp_results_total_df[mgp_results_total_df$Sex == "male", "Sex"] <- "M"
mgp_results_total_df[mgp_results_total_df$Sex == "XY", "Sex"] <- "M"
mgp_results_total_df[mgp_results_total_df$Sex == "Female", "Sex"] <- "F"
mgp_results_total_df[mgp_results_total_df$Sex == "female", "Sex"] <- "F"
mgp_results_total_df[mgp_results_total_df$Sex == "XX", "Sex"] <- "F"

### clarify subject IDs (mainly for subjects from Pitt)

mgp_results_total_df$subject_id <- mgp_results_total_df$rowname
mgp_results_total_df[mgp_results_total_df$Institution == "Pitt","subject_id"] <- readr::parse_number(mgp_results_total_df[mgp_results_total_df$Institution == "Pitt","subject_id"])
colnames(mgp_results_total_df)[1] <- "dataset_subject_id"

### formalize some values
unique(mgp_results_total_df$dataset)

mgp_results_total_df[mgp_results_total_df$Institution == "NIMH-HBCC", "Institution"] <- "NIMH"
mgp_results_total_df$dataset <- paste0(mgp_results_total_df$Institution, " RNAseq")

mgp_results_total_df[endsWith(mgp_results_total_df$dataset_subject_id, "_BA11"),"dataset"] <- "Pitt BA11 Microarray"
mgp_results_total_df[endsWith(mgp_results_total_df$dataset_subject_id, "_BA47"),"dataset"] <- "Pitt BA47 Microarray"
mgp_results_total_df[mgp_results_total_df$dataset == "NABEC RNAseq","dataset"] <- "NABEC Microarray"
mgp_results_total_df[mgp_results_total_df$dataset == "UKBEC RNAseq","dataset"] <- "UKBEC Microarray"

### Save/export

write.csv(mgp_results_total_df, "./data/output/tables/Allregion_MGP_Results.csv")
