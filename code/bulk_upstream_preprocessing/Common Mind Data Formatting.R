source('Functions.R')


# import commonmind data
commonmind_expression_raw = read_tsv(file = '/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/RNAseq/Release4/QuantitatedExpression/MSSM.Penn.Pitt_DLPFC.featureCount.tsv.gz')
rownames(commonmind_expression_raw) = commonmind_expression_raw$Geneid

commonmind_meta = read_csv(file = '/external/rprshnas01/netdata_kcni/stlab/bulk_expression/commonmind/CMC_MSSM-Penn-Pitt_Clinical.csv')
commonmind_meta$DLPFC_RNA_Sequencing_Sample_ID = gsub('\n', '', commonmind_meta$DLPFC_RNA_Sequencing_Sample_ID)
commonmind_meta$DLPFC_RNA_Sequencing_Sample_ID = gsub('\\\\', '', commonmind_meta$DLPFC_RNA_Sequencing_Sample_ID)
commonmind_meta[commonmind_meta$Age_of_Death == "90+", 'Age_of_Death'] = '90'
commonmind_meta$Age_of_Death = commonmind_meta$Age_of_Death %>% as.numeric()

cm_MSSM_meta = filter(commonmind_meta, Institution == 'MSSM')
cm_pitt_meta = filter(commonmind_meta, Institution == 'Pitt')
cm_penn_meta = filter(commonmind_meta, Institution == 'Penn')

#mutate the Geneid column to remove the suffix on the end of the ensembl ids                                                                                                                          
cm_ensemble_ids = vector()
for (id in commonmind_expression_raw$Geneid){
  new_id = substr(id, 1, 15)
  cm_ensemble_ids = c(cm_ensemble_ids, new_id)
}

commonmind_expression_raw$ensembl_id = cm_ensemble_ids

# convert ensembl ids to gene symbol
symbols <- mapIds(org.Hs.eg.db, keys = commonmind_expression_raw$ensembl_id, keytype = "ENSEMBL", column="SYMBOL")
symbols_df = symbols %>% as.data.frame()
symbols_df$ensembl_id = commonmind_expression_raw$ensembl_id
colnames(symbols_df) = c('gene_symbol', 'ensembl_id')



#add gene_symbols to expression df
commonmind_expression_raw = merge(symbols_df, commonmind_expression_raw, by = 'ensembl_id')
commonmind_expression_raw = filter(commonmind_expression_raw, !(is.na(gene_symbol)))

#remove duplicated genes
commonmind_gene_duplicates = duplicated(commonmind_expression_raw$gene_symbol)

commonmind_expression_raw = commonmind_expression_raw[!commonmind_gene_duplicates, ]
rownames(commonmind_expression_raw) = commonmind_expression_raw$gene_symbol

#seperate into the different cohorts
MSSM_controls_meta_df = filter(cm_MSSM_meta, Dx == 'Control')
MSSM_ids = MSSM_controls_meta_df$DLPFC_RNA_Sequencing_Sample_ID
MSSM_ids = intersect(MSSM_ids, colnames(commonmind_expression_raw))

pitt_controls_meta_df = filter(cm_pitt_meta, Dx == 'Control')
pitt_ids = pitt_controls_meta_df$DLPFC_RNA_Sequencing_Sample_ID
pitt_ids = intersect(pitt_ids, colnames(commonmind_expression_raw))

penn_controls_meta_df = filter(cm_penn_meta, Dx == 'Control')
penn_ids = penn_controls_meta_df$DLPFC_RNA_Sequencing_Sample_ID
penn_ids = intersect(penn_ids, colnames(commonmind_expression_raw))

#separate into seperate cohorts
MSSM_expression = commonmind_expression_raw[, MSSM_ids]
pitt_expression = commonmind_expression_raw[, pitt_ids]
penn_expression = commonmind_expression_raw[, penn_ids]


#filter out genes with med expression under 15 reads
MSSM_expression$med = rowMedians(MSSM_expression %>% as.matrix())
MSSM_expression_filtered = filter(MSSM_expression, med >= 15)

pitt_expression$med = rowMedians(pitt_expression %>% as.matrix())
pitt_expression_filtered = filter(pitt_expression, med >= 15)

penn_expression$med = rowMedians(penn_expression %>% as.matrix())
penn_expression_filtered = filter(penn_expression, med >= 15)

#filter and scale data - prep for MGP calculation ----
MSSM_expression_filtered$gene_symbol = rownames(MSSM_expression_filtered)
MSSM_expression_filtered = edgeR::cpm(MSSM_expression_filtered[-166], log = T, prior.count = 1) %>% as.data.frame()
MSSM_expression_filtered$gene_symbol = rownames(MSSM_expression_filtered)
MSSM_expression_filtered = dplyr::select(MSSM_expression_filtered, gene_symbol, everything())

pitt_expression_filtered$gene_symbol = rownames(pitt_expression_filtered)
pitt_expression_filtered = edgeR::cpm(pitt_expression_filtered[-86], log = T, prior.count = 1) %>% as.data.frame()
pitt_expression_filtered$gene_symbol = rownames(pitt_expression_filtered)
pitt_expression_filtered = dplyr::select(pitt_expression_filtered, gene_symbol, everything())

penn_expression_filtered$gene_symbol = rownames(penn_expression_filtered)
penn_expression_filtered = edgeR::cpm(penn_expression_filtered[-40], log = T, prior.count = 1) %>% as.data.frame()
penn_expression_filtered$gene_symbol = rownames(penn_expression_filtered)
penn_expression_filtered = dplyr::select(penn_expression_filtered, gene_symbol, everything())
