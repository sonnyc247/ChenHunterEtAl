source('Aging-Expression-Analysis/Functions.R')

# read in data
GSE36192 = read_rds('~/Final_Aging_Analysis/Aging-Expression-Analysis/GSE36192_out.rds')

#seperate into different brain regions, and then into metadata and expression

GSE36192_fc = filter(GSE36192, region == 'frontal cortex')
GSE36192_fc = filter(GSE36192_fc, age_years >= 15)
GSE36192_fc_meta = GSE36192_fc[, 1:5]
GSE36192_fc_expr = GSE36192_fc[, 6:19618]
rownames(GSE36192_fc_expr) = GSE36192_fc_meta$geo_accession

#filtering----

#transpose expression dataframe, make row of median expr per gene
GSE36192_fc_expr = t(GSE36192_fc_expr) %>% as.data.frame()
GSE36192_fc_expr$med = GSE36192_fc_expr %>% as.matrix() %>% rowMedians()

GSE36192_expr_25_quart = quantile(GSE36192_fc_expr$med)[2] %>% unname()

GSE36192_expr_filtered = filter(GSE36192_fc_expr, med >= GSE36192_expr_25_quart)
GSE36192_expr_filtered$gene = rownames(GSE36192_expr_filtered)
GSE36192_expr_filtered = dplyr::select(GSE36192_expr_filtered, gene, everything())
GSE36192_expr_filtered = GSE36192_expr_filtered[, 1:430] # removing med column

