source('Functions.R')
library(dplyr)
library(readxl)
#save file path in variable
expression_file = '/external/rprshnas01/netdata_kcni/stlab/bulk_expression/pittsburgh_aging_cohort/Age-BA11-47-Expr-Stats.xlsx'

#read in and format metadata dfs ----
sample_meta_colnames = read_excel(expression_file, range = 'A2:A8', col_names = F) %>% unlist() %>% unname()

sample_meta_ba47 = read_excel(expression_file, range = 'J3:HL9', col_names = F)
sample_meta_df_ba47 = sample_meta_ba47 %>% t() %>% as.data.frame()
colnames(sample_meta_df_ba47) = sample_meta_colnames
sample_meta_df_ba47$Region = 'BA47'
sample_meta_df_ba47 = mutate(sample_meta_df_ba47, ID = paste(ID, Region, sep = '_'))

sample_meta_ba11 = read_excel(expression_file, range = 'HM3:PO9', col_names = F)
sample_meta_df_ba11 = sample_meta_ba11 %>% t() %>% as.data.frame()
colnames(sample_meta_df_ba11) = sample_meta_colnames
sample_meta_df_ba11$Region = 'BA11'
sample_meta_df_ba11 = mutate(sample_meta_df_ba11, ID = paste(ID, Region, sep = '_'))

sample_ids = c(sample_meta_df_ba47$ID, sample_meta_df_ba11$ID)

#read in and format expression data ----
gene_exp_matrix = read_excel(expression_file, range = 'J10:PO33306', col_names = F)
colnames(gene_exp_matrix) = sample_ids

gene_names = read_excel(expression_file, range = 'A10:C33306', col_names = F)
colnames(gene_names) = c('Probeset_ID', 'Gene_Ascension', 'Gene_Symbol')

gene_exp_matrix$Probeset_ID = gene_names$Probeset_ID
gene_exp_matrix$Gene_Symbol = gene_names$Gene_Symbol

gene_exp_BA11 = dplyr::select(gene_exp_matrix, Probeset_ID, Gene_Symbol, sample_meta_df_ba11$ID) %>% as.data.frame()
gene_exp_BA47 = dplyr::select(gene_exp_matrix, Probeset_ID, Gene_Symbol, sample_meta_df_ba47$ID) %>% as.data.frame()
# filter for most variable probeset ----

#get median expr of the two datasets ## quantile first?
BA11_25_quart = gene_exp_BA11[, -1:-2] %>% as.matrix() %>% rowMedians() %>% quantile() %>% unname() %>% {.[[2]]}
BA47_25_quart = gene_exp_BA47[, -1:-2]  %>% as.matrix() %>% rowMedians() %>% quantile() %>% unname() %>% {.[[2]]}


#we filter these datasets differently from the ukbec and nabec datasets because
#these datasets contain multiple probesets for the same gene, and we want to
#reduce it down to just one probeset
BA11_expr_filtered = mostVariable(gene_exp_BA11[, -1], genes = 'Gene_Symbol', 
                             threshold = BA11_25_quart, threshFun = median)

BA47_expr_filtered = mostVariable(gene_exp_BA47[, -1], genes = 'Gene_Symbol', 
                             threshold = BA47_25_quart, threshFun = median)

