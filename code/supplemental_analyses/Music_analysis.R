#devtools::install_github('xuranw/MuSiC') # as needed

library(MuSiC)
library(Seurat)
library(Biobase) # needed to make annotated dfs
library(tidyverse)
library(ggpubr)

#
#### prep reference ####

Seu_AIBS_obj <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_AIBS_obj_update_07JUN21.rds")

# make desired taxonomy (groups excitatory cells together)
Seu_AIBS_obj$ExClass_OtherSubclass <- Seu_AIBS_obj$subclass_label # make custom metadata column, starting from subclasses
Seu_AIBS_obj$ExClass_OtherSubclass[Seu_AIBS_obj$ExClass_OtherSubclass %in% c("IT",
                                                                             "L4 IT",
                                                                             "L5 ET",
                                                                             "L5/6 IT Car3",
                                                                             "L5/6 NP",
                                                                             "L6 CT",
                                                                             "L6b")] <- "Excitatory" #set to excitatory
table(Seu_AIBS_obj$ExClass_OtherSubclass, useNA = "always") #check our data composition

# subsample data to facilitate use in MuSiC
Seu_AIBS_obj$filter_criteria <- paste0(Seu_AIBS_obj$external_donor_name_label, "_",
                                       Seu_AIBS_obj$region_label, "_", 
                                       Seu_AIBS_obj$subclass_label) # for every donor, region, and subclass
sort(table(Seu_AIBS_obj$filter_criteria, useNA = "always"))

Idents(Seu_AIBS_obj) <- "filter_criteria"
Seu_AIBS_filtered <- subset(Seu_AIBS_obj, downsample = 500)
sort(table(Seu_AIBS_filtered$filter_criteria, useNA = "always"))

rm(Seu_AIBS_obj)

# convert to ExpressionSet

pdata <- new("AnnotatedDataFrame", data=Seu_AIBS_filtered@meta.data)
AIBS_ES <- ExpressionSet(assayData = as.matrix(Seu_AIBS_filtered@assays$RNA@counts), 
                         phenoData = pdata)

saveRDS(AIBS_ES, "~/git/Ex_Env_Storage/Aging_MGP_private/AIBS_ES.rds")
rm(Seu_AIBS_filtered, Seu_AIBS_obj, pdata)
gc()

#### CMC (RNAseq data) MuSiC analysis ####

### Data preprocessing

# select subjects
METADATA <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/cmc_data/METADATA.rds")
METADATA <- METADATA[METADATA$Dx == "Control",]
METADATA <- METADATA[METADATA$ageOfDeath >= 15,]

# count matrix preprocessing (qc)
geneCountsMerged <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Aging_MGP_private/cmc_data/geneCountsMerged.rds")
geneCountsMerged <- geneCountsMerged[,METADATA$SampleID]
geneCountsMerged <- geneCountsMerged[rowSums(geneCountsMerged) >= 15, ]

# ensembl id mapping
hgnc_mapping <- readr::read_tsv("./data/hgnc_complete_set_2021-06-01.txt") # retrieved from 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2021-06-01.txt'
hgnc_mapping <- hgnc_mapping[,c("symbol", "ensembl_gene_id")]

geneCountsMerged <- rownames_to_column(as.data.frame(geneCountsMerged))
geneCountsMerged$rowname <- gsub("\\..*","",geneCountsMerged$rowname) # remove version id from ensembl IDs
geneCountsMerged <- merge(hgnc_mapping, geneCountsMerged, by.x = "ensembl_gene_id", by.y = "rowname", all = F)
row.names(geneCountsMerged) <- make.names(geneCountsMerged$symbol)
geneCountsMerged <- geneCountsMerged[,3:ncol(geneCountsMerged)]

### MuSiC run

# make data into ES
pdata <- new("AnnotatedDataFrame", data=METADATA)
Bulk_ES <- ExpressionSet(assayData = as.matrix(geneCountsMerged), 
                         phenoData = pdata)

# MuSiC run
Result <- music_prop(bulk.eset = Bulk_ES, 
                     sc.eset = AIBS_ES, 
                     clusters = 'ExClass_OtherSubclass',
                     samples = 'external_donor_name_label')
setdiff(unique(AIBS_ES$ExClass_OtherSubclass), colnames(Result$Est.prop.weighted)) # all cell types estimated

saveRDS(Result, "./data/output/result_objects/MuSiC_CMC.rds")

# QC checks
hist(Result$r.squared.full)
max(Result$Var.prop)
min(Result$Var.prop)
min(Music_result[Music_result$music_prop > 0, "music_prop"])

# process and format results
Music_result <- Result$Est.prop.weighted
Music_result <- rownames_to_column(as.data.frame(Music_result))
Music_result <- merge(Music_result, 
                      METADATA[,c("SampleID", "ageOfDeath", "Sex", "PMI", "Institution")], 
                      by.x = "rowname", 
                      by.y = "SampleID", 
                      all.x = T, 
                      all.y = F)
Music_result$Dataset <- Music_result$Institution

# wrap up

write.csv(Music_result, "./data/output/tables/MuSiC_CMC_table.csv")
all_music_results <- Music_result

rm(Bulk_ES, geneCountsMerged, hgnc_mapping, METADATA, Music_result, pdata, Result)
gc()

#
#### Do age effect estimates with CMC music results ####

#library(lmerTest)

mgp_results_MuSiC <- read.csv("./data/output/tables/MuSiC_CMC_table.csv", row.names = 1)

mgp_results_MuSiC$Sex <- as.factor(mgp_results_MuSiC$Sex)
celltypes <- colnames(mgp_results_MuSiC)[2:14]

mega_results_MuSiC <- data.frame(term = character(), 
                               estimate = numeric(),
                               std.error = numeric(),
                               p.value = numeric(), 
                               celltype = character())

mgp_results_MuSiC$age_scaled <- scale(mgp_results_MuSiC$ageOfDeath)
mgp_results_MuSiC$pmi_scaled <- scale(mgp_results_MuSiC$PMI)

for (cell in celltypes) {
  formula <- paste(cell, '~ 1 + age_scaled + pmi_scaled  + Sex + Institution', sep = ' ') # all rownames/subject ID unique, can't use
  model <- lm(data = mgp_results_MuSiC, formula = formula) # thus no random effect, thus use lm
  model_df <- broom::tidy(model)
  model_df <- model_df[c(2:4),] 
  model_df <- model_df[,-4]
  model_df$celltype <- cell
  mega_results_MuSiC <- rbind(mega_results_MuSiC, model_df)
}

mega_results_MuSiC$p.value.adjusted <- stats::p.adjust(mega_results_MuSiC$p.value, method = "BH")
mega_results_MuSiC_age <- mega_results_MuSiC[mega_results_MuSiC$term == "age_scaled",c(5,2:4,6)]
mega_results_MuSiC_sex <- mega_results_MuSiC[mega_results_MuSiC$term == "SexXY",c(5,2:4,6)]

# save results

names(mega_results_MuSiC_age) <- c("Cell type",
                             "Std. beta coefficient",
                             "Std. error",
                             "p value",
                             "FDR")

write.csv(mega_results_MuSiC_age, "./data/output/tables/MUSIC_Mega_analysis.csv")

### individual analysis 

indiv_results_MuSiC <- data.frame(term = character(), 
                                  estimate = numeric(),
                                  std.error = numeric(),
                                  p.value = numeric(), 
                                  celltype = character(),
                                  dataset = character())

datasets <- unique(mgp_results_MuSiC$Dataset)

for (cell in celltypes){
  for (dataset in datasets){
    test_df <- mgp_results_MuSiC[,c(cell, "age_scaled", "Sex", "pmi_scaled", "Dataset")]
    test_df <- test_df[test_df$Dataset == dataset,]
    formula = paste(cell, '~ age_scaled + Sex + pmi_scaled', sep = ' ')
    model <- lm(formula, test_df)
    model_df <- broom::tidy(model)
    model_df <- model_df[c(2:4),] 
    model_df <- model_df[,-4]
    model_df$celltype <- cell
    model_df$dataset <- dataset
    indiv_results_MuSiC <- rbind(indiv_results_MuSiC, model_df)
  }
}

indiv_age_results <- indiv_results_MuSiC[indiv_results_MuSiC$term == "age_scaled",]
indiv_age_results$p.value.adjusted <- stats::p.adjust(indiv_age_results$p.value, method = "BH")
indiv_age_results <- indiv_age_results[,c(6,5,2:4,7)]

### save results

names(indiv_age_results) <- c("Dataset",
                              "Cell type",
                              "Std. beta coefficient",
                              "Std. error",
                              "p value",
                              "FDR")

write.csv(indiv_results_MuSiC, "./data/output/tables/MUSIC_Indiv_analysis.csv")

#### Preliminary plots of MuSiC Results ####

figure_2_df <- mgp_results_MuSiC

#V2 plot:

# rename

names(figure_2_df)[19] <- "Dataset"

figure_2_df$Dataset <- factor(figure_2_df$Dataset, levels = c("MSSM RNAseq",
                                                              "Penn RNAseq",
                                                              "Pitt RNAseq",
                                                              "NIMH RNAseq",
                                                              "Pitt BA11 Microarray",
                                                              "Pitt BA47 Microarray",
                                                              "NABEC Microarray",
                                                              "UKBEC Microarray"))

figure_2_df$Dataset <- factor(figure_2_df$Dataset, levels = c("MSSM",
                                                              "Penn",
                                                              "Pitt",
                                                              "NIMH-HBCC",
                                                              "BA11",
                                                              "BA47",
                                                              "NABEC",
                                                              "UKBEC"))

gathered_figure_2_df <- tidyr::gather(figure_2_df, key = "Celltype", value = "Scaled_MGP", Excitatory, SST, Oligodendrocyte)
gathered_figure_2_df <- tidyr::gather(figure_2_df, key = "Celltype", value = "MuSiC_estimate", 2:14)

gathered_figure_2_df$Celltype <- factor(gathered_figure_2_df$Celltype, levels = c("SST",
                                                                                  "Excitatory",
                                                                                  "Oligodendrocyte"))

ggplot(data = gathered_figure_2_df, aes(x = ageOfDeath, y = MuSiC_estimate)) + 
  geom_point(size = 0.35, alpha = 0.3) + 
  geom_smooth(method = 'lm', se = F, size = 0.45, color = "red") +
  geom_smooth(method = "loess", se = F, size = 0.45) +
  facet_grid(Celltype ~ Dataset, switch = "y", scale="free_y") +
  scale_y_continuous(position = "right") +
  xlab('Age of death') +
  ylab('Standardized rCTP') +
  stat_cor(aes(label = ..r.label..),  label.x = 3, size = 3) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black', size = 8),
        panel.spacing = unit(0.3, "lines"),
        text = element_text(size = 8), 
        axis.text = element_text(size = 6.5),
        axis.title.y.right =  element_text(margin = margin(t = 0, r = 0, b = 0, l = 5), size = 12),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0), size = 12),
        strip.text.x = element_text(size = 12))

ggsave(width = 180,
       height = 500,
       dpi = 300, 
       units = "mm",
       bg = "white",
       path = "/external/rprshnas01/kcni/ychen/git/Aging_MGP_private/data/output/figures/",
       filename = "MuSiC_scatter.pdf",
       device = "pdf")

#### Compare with snCTPs ####

mgp_results_MuSiC <- read.csv("./data/output/tables/MuSiC_CMC_table.csv", row.names = 1)
mgp_results_MuSiC <- mgp_results_MuSiC[mgp_results_MuSiC$Dataset == "MSSM",]
mgp_results_MuSiC <- mgp_results_MuSiC[,1:14]
mgp_results_MuSiC <- gather(mgp_results_MuSiC, 2:14, key = "CellType", value = "MuSiC_Estimate")
mgp_results_MuSiC$mergecol <- paste0(mgp_results_MuSiC$rowname, mgp_results_MuSiC$CellType)

snCTP_results_harmonized <- read.csv('code/collab/data/snCTP_results_harmonized.csv', row.names = 1)
snCTP_results_harmonized <- snCTP_results_harmonized[,c(1,3:13)]
snCTP_results_harmonized <- gather(snCTP_results_harmonized, 2:12, key = "CellType", value = "sn_Estimate")
snCTP_results_harmonized$mergecol <- paste0(snCTP_results_harmonized$subject_id, snCTP_results_harmonized$CellType)

setdiff(unique(mgp_results_MuSiC$CellType), unique(snCTP_results_harmonized$CellType)) # no pax6 and vlmc in sn
length(intersect(mgp_results_MuSiC$mergecol, snCTP_results_harmonized$mergecol))

MuSiCvSN <- merge(mgp_results_MuSiC, snCTP_results_harmonized[,c("sn_Estimate", "mergecol")], by = "mergecol")
MuSiCvSN <- MuSiCvSN[,-1]

ggplot(MuSiCvSN, aes(x = sn_Estimate, y = MuSiC_Estimate)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~CellType, scales = "free") +
  theme_classic()

ggsave(width = 300,
       height = 300,
       dpi = 300, 
       units = "mm",
       bg = "white",
       path = "/external/rprshnas01/kcni/ychen/git/Aging_MGP_private/data/output/figures/",
       filename = "MuSiCvSN_scatter.pdf",
       device = "pdf")
