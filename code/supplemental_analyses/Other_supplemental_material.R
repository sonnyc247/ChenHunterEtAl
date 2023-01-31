### Run after MGP analysis, modeling analysis, and graphing code

#### Packages and markers ####

library(dplyr)
library(ggplot2)

### List of interested genes

markers_df <- read.csv("./data/output/tables/Supplementary_Table_1_AllRegions.csv", row.names = 1)
markers_df <- markers_df[markers_df$Used.in.MGP == T,]
markers_df <- markers_df[markers_df$Cell.type %in% c("SST", "Excitatory", "Astrocyte"), ]
markers_df <- markers_df %>% group_by(Cell.type) %>% top_n(20, Average.log.2.fold.change)

#### CMC ####

METADATA <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/cmc_data/METADATA.rds")
METADATA <- METADATA[METADATA$Dx == "Control",]
METADATA <- METADATA[METADATA$ageOfDeath >= 15,]
METADATA$is_duplicated <- duplicated(METADATA$IndividualID)

### count matrix preprocessing

geneCountsMerged <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/cmc_data/geneCountsMerged.rds")
geneCountsMerged <- geneCountsMerged[,METADATA$SampleID]
geneCountsMerged <- geneCountsMerged[rowSums(geneCountsMerged) >= 15, ]

for (subject in METADATA[METADATA$is_duplicated == T, "IndividualID"]) {
  
  geneCountsMerged$temp <- rowSums(geneCountsMerged[, METADATA[METADATA$IndividualID == subject, "SampleID"]])
  colnames(geneCountsMerged)[ncol(geneCountsMerged)] <- paste0("PittRNAseq_", readr::parse_number(METADATA[METADATA$IndividualID == subject, "SampleID"])[1])
  
}

geneCountsMerged <- geneCountsMerged[, !(names(geneCountsMerged) %in% METADATA[METADATA$IndividualID %in% METADATA[METADATA$is_duplicated == T, "IndividualID"], "SampleID"])]
geneCountsMerged <- edgeR::cpm(geneCountsMerged, log = TRUE, prior.count = 1) 

### ensembl id mapping

hgnc_mapping <- readr::read_tsv("./data/hgnc_complete_set_2021-06-01.txt") # retrieved from 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2021-06-01.txt'
hgnc_mapping <- hgnc_mapping[,c("symbol", "ensembl_gene_id")]

geneCountsMerged <- tibble::rownames_to_column(as.data.frame(geneCountsMerged))
geneCountsMerged$rowname <- gsub("\\..*","",geneCountsMerged$rowname) # remove version id from ensembl IDs
geneCountsMerged <- merge(hgnc_mapping, geneCountsMerged, by.x = "ensembl_gene_id", by.y = "rowname", all = F)
row.names(geneCountsMerged) <- make.names(geneCountsMerged$symbol)
geneCountsMerged <- geneCountsMerged[,3:ncol(geneCountsMerged)]
geneCountsMerged <- t(as.matrix(geneCountsMerged[markers_df$Gene.symbol,]))

### merge metadata (from mgp results)

mgp_results_total_df <- read.csv("./data/output/tables/Allregion_MGP_Results.csv", row.names = 1)
geneCountsMerged <- merge(as.data.frame(geneCountsMerged), mgp_results_total_df[,c("dataset_subject_id", "age_scaled", "Sex", "pmi_scaled", "dataset")], by.x = "row.names", by.y = "dataset_subject_id")

### get bet values

gene_beta_results = data.frame(term = character(), 
                           estimate = numeric(),
                           std.error = numeric(),
                           p.value = numeric(), 
                           gene = character(),
                           dataset = character())

for (gene in markers_df$Gene.symbol){
  for (dataset in unique(geneCountsMerged$dataset)){
    test_df <- geneCountsMerged[,c(gene, "age_scaled", "Sex", "pmi_scaled", "dataset")]
    test_df <- test_df[test_df$dataset == dataset,]
    formula = paste(gene, '~ age_scaled + Sex + pmi_scaled', sep = ' ')
    model <- lm(formula, test_df)
    model_df <- broom::tidy(model)
    model_df <- model_df[c(2:4),] 
    model_df <- model_df[,-4]
    model_df$gene <- gene
    model_df$dataset <- dataset
    gene_beta_results <- rbind(gene_beta_results, model_df)
  }
  
}

#### NABEC ####

### preprocessing

NABEC_microarray <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/BEC_microarray/GSE36192_out.rds")
NABEC_microarray <- NABEC_microarray[NABEC_microarray$region == "frontal cortex",]
NABEC_microarray <- NABEC_microarray[NABEC_microarray$age_years >= 15,] 
row.names(NABEC_microarray) <- NABEC_microarray$geo_accession
NABEC_microarray <- NABEC_microarray[,markers_df$Gene.symbol]
NABEC_microarray <- merge(NABEC_microarray, mgp_results_total_df[,c("dataset_subject_id", "age_scaled", "Sex", "pmi_scaled", "dataset")], by.x = "row.names", by.y = "dataset_subject_id")

for (gene in markers_df$Gene.symbol){

    test_df <- NABEC_microarray[,c(gene, "age_scaled", "Sex", "pmi_scaled", "dataset")]
    formula = paste(gene, '~ age_scaled + Sex + pmi_scaled', sep = ' ')
    model <- lm(formula, test_df)
    model_df <- broom::tidy(model)
    model_df <- model_df[c(2:4),] 
    model_df <- model_df[,-4]
    model_df$gene <- gene
    model_df$dataset <- unique(test_df$dataset)
    gene_beta_results <- rbind(gene_beta_results, model_df)
  
}

#### UKBEC ####

### preprocessing

UKBEC_microarray <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/BEC_microarray/GSE60862_out.rds")
UKBEC_microarray <- UKBEC_microarray[UKBEC_microarray$region == "frontal cortex",]
row.names(UKBEC_microarray) <- UKBEC_microarray$geo_accession
UKBEC_microarray <- UKBEC_microarray[,markers_df$Gene.symbol]
UKBEC_microarray <- merge(UKBEC_microarray, mgp_results_total_df[,c("dataset_subject_id", "age_scaled", "Sex", "pmi_scaled", "dataset")], by.x = "row.names", by.y = "dataset_subject_id")

for (gene in markers_df$Gene.symbol){
  
  test_df <- UKBEC_microarray[,c(gene, "age_scaled", "Sex", "pmi_scaled", "dataset")]
  formula = paste(gene, '~ age_scaled + Sex + pmi_scaled', sep = ' ')
  model <- lm(formula, test_df)
  model_df <- broom::tidy(model)
  model_df <- model_df[c(2:4),] 
  model_df <- model_df[,-4]
  model_df$gene <- gene
  model_df$dataset <- unique(test_df$dataset)
  gene_beta_results <- rbind(gene_beta_results, model_df)
  
}

#### BA11 #### 

### Pre processing
BA11_expr_filtered <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/pitt_microarray/BA11_expr_filtered.rds")
row.names(BA11_expr_filtered) <- BA11_expr_filtered$Gene_Symbol
BA11_expr_filtered <- BA11_expr_filtered[markers_df$Gene.symbol,-1]
BA11_expr_filtered <- t(log2(BA11_expr_filtered))
BA11_expr_filtered <- merge(as.data.frame(BA11_expr_filtered), mgp_results_total_df[,c("dataset_subject_id", "age_scaled", "Sex", "pmi_scaled", "dataset")], by.x = "row.names", by.y = "dataset_subject_id")

for (gene in markers_df$Gene.symbol){
  
  test_df <- BA11_expr_filtered[,c(gene, "age_scaled", "Sex", "pmi_scaled", "dataset")]
  formula = paste(gene, '~ age_scaled + Sex + pmi_scaled', sep = ' ')
  model <- lm(formula, test_df)
  model_df <- broom::tidy(model)
  model_df <- model_df[c(2:4),] 
  model_df <- model_df[,-4]
  model_df$gene <- gene
  model_df$dataset <- unique(test_df$dataset)
  gene_beta_results <- rbind(gene_beta_results, model_df)
  
}

#### BA47 ####

### Pre processing
BA47_expr_filtered <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/pitt_microarray/BA47_expr_filtered.rds")
row.names(BA47_expr_filtered) <- BA47_expr_filtered$Gene_Symbol
BA47_expr_filtered <- BA47_expr_filtered[markers_df$Gene.symbol,-1]
BA47_expr_filtered <- t(log2(BA47_expr_filtered))

BA47_expr_filtered <- merge(as.data.frame(BA47_expr_filtered), mgp_results_total_df[,c("dataset_subject_id", "age_scaled", "Sex", "pmi_scaled", "dataset")], by.x = "row.names", by.y = "dataset_subject_id")

for (gene in markers_df$Gene.symbol){
  
  test_df <- BA47_expr_filtered[,c(gene, "age_scaled", "Sex", "pmi_scaled", "dataset")]
  formula = paste(gene, '~ age_scaled + Sex + pmi_scaled', sep = ' ')
  model <- lm(formula, test_df)
  model_df <- broom::tidy(model)
  model_df <- model_df[c(2:4),] 
  model_df <- model_df[,-4]
  model_df$gene <- gene
  model_df$dataset <- unique(test_df$dataset)
  gene_beta_results <- rbind(gene_beta_results, model_df)
  
}

#### Plot ####

gene_beta_results$dataset <- factor(gene_beta_results$dataset, 
                                             levels = c("MSSM RNAseq",
                                                        "Penn RNAseq",
                                                        "Pitt RNAseq",
                                                        "NIMH RNAseq",
                                                        "Pitt BA11 Microarray",
                                                        "Pitt BA47 Microarray",
                                                        "NABEC Microarray",
                                                        "UKBEC Microarray"))

gene_beta_results <- merge(gene_beta_results, unique(markers_df[,c("Gene.symbol", "Cell.type")]), by.x = "gene", by.y = "Gene.symbol")

gene_beta_results$Cell.type <- factor(gene_beta_results$Cell.type, 
                                              levels = c("SST",
                                                         "Excitatory",
                                                         "Astrocyte"))

sup_figure_df <- gene_beta_results[gene_beta_results$term == "age_scaled",-2]
#markers_df <- read.csv("./data/output/tables/Supplementary_Table_1.csv")
markers_df <- markers_df[markers_df$Gene.symbol %in% sup_figure_df$gene, c("Gene.symbol", "Average.log.2.fold.change")]
sup_figure_df <- merge(sup_figure_df, markers_df, by.x = "gene", by.y = "Gene.symbol")
sup_figure_df$gene <- factor(sup_figure_df$gene)
sup_figure_df$gene <- reorder(sup_figure_df$gene, sup_figure_df$Average.log.2.fold.change, decreasing = T)

ggplot(sup_figure_df) +
  geom_tile(aes(x = dataset, y = gene, fill = estimate)) +
  scale_fill_gradientn(limits = c(-1,1),
                       colours = c("blue", "white", "red"),
                       breaks = c(-1,0,1),
                       labels = format(c(-1,0,1))) +
  facet_wrap(~Cell.type, ncol = 1, scales = "free_y") +
  labs(x = "Dataset", y = "Marker gene", fill = "Std. beta") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 20, vjust = 0.95, hjust = 0.99))

ggsave(width = 140,
       height = 275,
       dpi = 300, 
       units = "mm",
       bg = "white",
       path = "./data/output/figures/",
       filename = "SupFig1AllBrainMarkers.pdf",
       device = "pdf")
