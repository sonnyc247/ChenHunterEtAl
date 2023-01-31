source('Aging-Expression-Analysis/Functions.R')

# read in data
GSE60862 = read_rds('~/Final_Aging_Analysis/Aging-Expression-Analysis/GSE60862_out.rds')

#filter for frontal cortex data, and seperate into meta and expr data
GSE60862_fc = filter(GSE60862, region == 'frontal cortex')
GSE60862_fc_meta = GSE60862_fc[, 1:8]
GSE60862_fc_expr = GSE60862_fc[, 9:20820]
rownames(GSE60862_fc_meta) = GSE60862_fc_meta$geo_accession

#filtering ----
GSE60862_expr = GSE60862_fc[, -1:-8] %>% t() %>% as.data.frame()
rownames(GSE60862_expr) = colnames(GSE60862_fc[, -1:-8])
colnames(GSE60862_expr) = GSE60862_fc$geo_accession

GSE60862_expr$med = GSE60862_expr %>% as.matrix() %>% rowMedians()
GSE60862_expr_25_quart = quantile(GSE60862_expr$med)[2] %>% unname() #quartiles by default?

GSE60862_expr_filtered = filter(GSE60862_expr, med >= GSE60862_expr_25_quart)
#remove med column
GSE60862_expr_filtered = GSE60862_expr_filtered[, 1:127] # remove med column
GSE60862_expr_filtered$gene = rownames(GSE60862_expr_filtered)
GSE60862_expr_filtered = dplyr::select(GSE60862_expr_filtered, gene, everything())
