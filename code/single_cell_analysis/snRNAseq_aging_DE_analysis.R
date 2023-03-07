library(edgeR)
library(data.table)
library(plyr)
library(lme4)
library(ggrepel)
library(ggplot2)

setwd("aging_genexp_analysis/")

directory_name <- "All_Cohorts"

# Load SZBD pseudobulk data
SZBD_PBcounts_filename <- "/external/rprshnas01/netdata_kcni/stlab/SZBDMulticohort/pseudobulk_fine_grained.counts.tsv"
SZBD_PBmetadata_filename <- "/external/rprshnas01/netdata_kcni/stlab/SZBDMulticohort/pseudobulk_fine_grained.metadata.tsv"

counts = fread(SZBD_PBcounts_filename)
metadata = fread(SZBD_PBmetadata_filename)

# Only keep controls
counts_con = counts[metadata$Phenotype=="CON",]
metadata_con = metadata[metadata$Phenotype=="CON",]

# Group cell classes of interest
counts_Excitatory <- counts_con[startsWith(as.character(Celltype), 'Ex'),]
counts_LAMP5 <- counts_con[startsWith(as.character(Celltype), 'In-Rosehip'),]
counts_VIP <- counts_con[Celltype %in% c('In-Reelin', 'In-VIP'),]
counts_PVALB <- counts_con[startsWith(as.character(Celltype), 'In-PV'),]
counts_SST <- counts_con[Celltype == 'In-SST',]
counts_OPC <- counts_con[Celltype == 'OPC',]
counts_Oligodendrocyte <- counts_con[Celltype == 'Oli',]
counts_Astrocyte <- counts_con[Celltype == 'Ast',]

metadata_Excitatory <- metadata_con[startsWith(as.character(Celltype), 'Ex'),]
metadata_LAMP5 <- metadata_con[startsWith(as.character(Celltype), 'In-Rosehip'),]
metadata_VIP <- metadata_con[Celltype %in% c('In-Reelin', 'In-VIP'),]
metadata_PVALB <- metadata_con[startsWith(as.character(Celltype), 'In-PV'),]
metadata_SST <- metadata_con[Celltype == 'In-SST',]
metadata_OPC <- metadata_con[Celltype == 'OPC',]
metadata_Oligodendrocyte <- metadata_con[Celltype == 'Oli',]
metadata_Astrocyte <- metadata_con[Celltype == 'Ast',]

# Merge pseudobulk counts with common unique_donor_ID
sum_cell_type_counts_of_same_donor <- function(cdf) {
  cdf_m <- c()
  for (xx in unique(cdf$unique_donor_ID)){
    t <- unname(colSums(cdf[cdf$unique_donor_ID %in% xx,3:ncol(cdf)]))
    cdf_m <- rbind(cdf_m,t)
  }
  row_names <- unique(cdf$unique_donor_ID)
  rownames(cdf_m) <- row_names
  col_names <- names(cdf)[3:ncol(cdf)]
  colnames(cdf_m) <- col_names
  cdf <- data.frame(cdf_m)
  return(cdf)
}

counts_Excitatory <- sum_cell_type_counts_of_same_donor(counts_Excitatory)
counts_LAMP5 <- sum_cell_type_counts_of_same_donor(counts_LAMP5)
counts_VIP <- sum_cell_type_counts_of_same_donor(counts_VIP)
counts_PVALB <- sum_cell_type_counts_of_same_donor(counts_PVALB)
counts_SST <- sum_cell_type_counts_of_same_donor(counts_SST)
counts_OPC <- sum_cell_type_counts_of_same_donor(counts_OPC)
counts_Oligodendrocyte <- sum_cell_type_counts_of_same_donor(counts_Oligodendrocyte)
counts_Astrocyte <- sum_cell_type_counts_of_same_donor(counts_Astrocyte)

# Merge pseudobulk metadata with common unique_donor_ID
reshape_metadata_and_sum_cell_counts <- function(metadat) {
  metadat <- subset(metadat, select = -c(Celltype)) # drop Celltype column
  
  df <- data.frame(matrix(ncol = ncol(metadat), nrow = 0))
  colnames(df) <- colnames(metadat)
  for (xx in unique(metadat$unique_donor_ID)){
    unique_metadat <- metadat[metadat$unique_donor_ID %in% xx,1:ncol(metadat)]
    dropped_metadata <- unique_metadat[c(1),] # drop duplicate rows
    dropped_metadata$num_cells <- sum(unique_metadat$num_cells) # update cell count with summed count across merged cell types
    df[nrow(df) + 1,] = dropped_metadata
  }
  return(df)
}

metadata_Excitatory <- reshape_metadata_and_sum_cell_counts(metadata_Excitatory)
metadata_LAMP5 <- reshape_metadata_and_sum_cell_counts(metadata_LAMP5)
metadata_VIP <- reshape_metadata_and_sum_cell_counts(metadata_VIP)
metadata_PVALB <- reshape_metadata_and_sum_cell_counts(metadata_PVALB)
metadata_SST <- reshape_metadata_and_sum_cell_counts(metadata_SST)
metadata_OPC <- reshape_metadata_and_sum_cell_counts(metadata_OPC)
metadata_Oligodendrocyte <- reshape_metadata_and_sum_cell_counts(metadata_Oligodendrocyte)
metadata_Astrocyte <- reshape_metadata_and_sum_cell_counts(metadata_Astrocyte)

# Remove people with fewer than 10 cells of that type
counts_Excitatory_qc1 <- counts_Excitatory[(metadata_Excitatory$num_cells > 10),]
counts_LAMP5_qc1 <- counts_LAMP5[(metadata_LAMP5$num_cells > 10),]
counts_VIP_qc1 <- counts_VIP[(metadata_VIP$num_cells > 10),]
counts_PVALB_qc1 <- counts_PVALB[(metadata_PVALB$num_cells > 10),]
counts_SST_qc1 <- counts_SST[(metadata_SST$num_cells > 10),]
counts_OPC_qc1 <- counts_OPC[(metadata_OPC$num_cells > 10),]
counts_Oligodendrocyte_qc1 <- counts_Oligodendrocyte[(metadata_Oligodendrocyte$num_cells > 10),]
counts_Astrocyte_qc1 <- counts_Astrocyte[(metadata_Astrocyte$num_cells > 10),]

metadata_Excitatory_qc1 <- metadata_Excitatory[(metadata_Excitatory$num_cells > 10),]
metadata_LAMP5_qc1 <- metadata_LAMP5[(metadata_LAMP5$num_cells > 10),]
metadata_VIP_qc1 <- metadata_VIP[(metadata_VIP$num_cells > 10),]
metadata_PVALB_qc1 <- metadata_PVALB[(metadata_PVALB$num_cells > 10),]
metadata_SST_qc1 <- metadata_SST[(metadata_SST$num_cells > 10),]
metadata_OPC_qc1 <- metadata_OPC[(metadata_OPC$num_cells > 10),]
metadata_Oligodendrocyte_qc1 <- metadata_Oligodendrocyte[(metadata_Oligodendrocyte$num_cells > 10),]
metadata_Astrocyte_qc1 <- metadata_Astrocyte[(metadata_Astrocyte$num_cells > 10),]

# Only keep genes with >1 pseudobulk count in more than 80% of subjects for a given cell type
counts_Excitatory_qc2 <- counts_Excitatory_qc1[, (colSums(counts_Excitatory_qc1 > 1)/nrow(counts_Excitatory_qc1) > 0.8)]
counts_LAMP5_qc2 <- counts_LAMP5_qc1[, (colSums(counts_LAMP5_qc1 > 1)/nrow(counts_LAMP5_qc1) > 0.8)]
counts_VIP_qc2 <- counts_VIP_qc1[, (colSums(counts_VIP_qc1 > 1)/nrow(counts_VIP_qc1) > 0.8)]
counts_PVALB_qc2 <- counts_PVALB_qc1[, (colSums(counts_PVALB_qc1 > 1)/nrow(counts_PVALB_qc1) > 0.8)]
counts_SST_qc2 <- counts_SST_qc1[, (colSums(counts_SST_qc1 > 1)/nrow(counts_SST_qc1) > 0.8)]
counts_OPC_qc2 <- counts_OPC_qc1[, (colSums(counts_OPC_qc1 > 1)/nrow(counts_OPC_qc1) > 0.8)]
counts_Oligodendrocyte_qc2 <- counts_Oligodendrocyte_qc1[, (colSums(counts_Oligodendrocyte_qc1 > 1)/nrow(counts_Oligodendrocyte_qc1) > 0.8)]
counts_Astrocyte_qc2 <- counts_Astrocyte_qc1[, (colSums(counts_Astrocyte_qc1 > 1)/nrow(counts_Astrocyte_qc1) > 0.8)]

# Assign covariates, create design matrix, and normalize count matrix
dmat_Excitatory <- model.matrix(~1+Age+log10(num_cells)+PMI+Cohort+Gender, metadata_Excitatory_qc1)
dmat_LAMP5 <- model.matrix(~1+Age+log10(num_cells)+PMI+Cohort+Gender, metadata_LAMP5_qc1)
dmat_VIP <- model.matrix(~1+Age+log10(num_cells)+PMI+Cohort+Gender, metadata_VIP_qc1)
dmat_PVALB <- model.matrix(~1+Age+log10(num_cells)+PMI+Cohort+Gender, metadata_PVALB_qc1)
dmat_SST <- model.matrix(~1+Age+log10(num_cells)+PMI+Cohort+Gender, metadata_SST_qc1)
dmat_OPC <- model.matrix(~1+Age+log10(num_cells)+PMI+Cohort+Gender, metadata_OPC_qc1)
dmat_Oligodendrocyte <- model.matrix(~1+Age+log10(num_cells)+PMI+Cohort+Gender, metadata_Oligodendrocyte_qc1)
dmat_Astrocyte <- model.matrix(~1+Age+log10(num_cells)+PMI+Cohort+Gender, metadata_Astrocyte_qc1)

Norm_Excitatory <- calcNormFactors(DGEList(t(counts_Excitatory_qc2)))
Norm_LAMP5 <- calcNormFactors(DGEList(t(counts_LAMP5_qc2)))
Norm_VIP <- calcNormFactors(DGEList(t(counts_VIP_qc2)))
Norm_PVALB <- calcNormFactors(DGEList(t(counts_PVALB_qc2)))
Norm_SST <- calcNormFactors(DGEList(t(counts_SST_qc2)))
Norm_OPC <- calcNormFactors(DGEList(t(counts_OPC_qc2)))
Norm_Oligodendrocyte <- calcNormFactors(DGEList(t(counts_Oligodendrocyte_qc2)))
Norm_Astrocyte <- calcNormFactors(DGEList(t(counts_Astrocyte_qc2)))

y_Excitatory <- voom(Norm_Excitatory, design=dmat_Excitatory, plot = T)
y_LAMP5 <- voom(Norm_LAMP5, design=dmat_LAMP5, plot = T)
y_VIP <- voom(Norm_VIP, design=dmat_VIP, plot = T)
y_PVALB <- voom(Norm_PVALB, design=dmat_PVALB, plot = T)
y_SST <- voom(Norm_SST, design=dmat_SST, plot = T)
y_OPC <- voom(Norm_OPC, design=dmat_OPC, plot = T)
y_Oligodendrocyte <- voom(Norm_Oligodendrocyte, design=dmat_Oligodendrocyte, plot = T)
y_Astrocyte <- voom(Norm_Astrocyte, design=dmat_Astrocyte, plot = T)

fit_Excitatory <- lmFit(y_Excitatory, dmat_Excitatory)
fit_LAMP5 <- lmFit(y_LAMP5, dmat_LAMP5)
fit_VIP <- lmFit(y_VIP, dmat_VIP)
fit_PVALB <- lmFit(y_PVALB, dmat_PVALB)
fit_SST <- lmFit(y_SST, dmat_SST)
fit_OPC <- lmFit(y_OPC, dmat_OPC)
fit_Oligodendrocyte <- lmFit(y_Oligodendrocyte, dmat_Oligodendrocyte)
fit_Astrocyte <- lmFit(y_Astrocyte, dmat_Astrocyte)

tmp_Excitatory <- contrasts.fit(fit_Excitatory,coef=2)
tmp_LAMP5 <- contrasts.fit(fit_LAMP5,coef=2)
tmp_VIP <- contrasts.fit(fit_VIP,coef=2)
tmp_PVALB <- contrasts.fit(fit_PVALB,coef=2)
tmp_SST <- contrasts.fit(fit_SST,coef=2)
tmp_OPC <- contrasts.fit(fit_OPC,coef=2)
tmp_Oligodendrocyte <- contrasts.fit(fit_Oligodendrocyte,coef=2)
tmp_Astrocyte <- contrasts.fit(fit_Astrocyte,coef=2)

tmp_Excitatory <- eBayes(tmp_Excitatory)
tmp_LAMP5 <- eBayes(tmp_LAMP5)
tmp_VIP <- eBayes(tmp_VIP)
tmp_PVALB <- eBayes(tmp_PVALB)
tmp_SST <- eBayes(tmp_SST)
tmp_OPC <- eBayes(tmp_OPC)
tmp_Oligodendrocyte <- eBayes(tmp_Oligodendrocyte)
tmp_Astrocyte <- eBayes(tmp_Astrocyte)

top_Excitatory.table <- topTable(tmp_Excitatory, sort.by = "P", n = Inf, adjust.method = "BH")
top_LAMP5.table <- topTable(tmp_LAMP5, sort.by = "P", n = Inf, adjust.method = "BH")
top_VIP.table <- topTable(tmp_VIP, sort.by = "P", n = Inf, adjust.method = "BH")
top_PVALB.table <- topTable(tmp_PVALB, sort.by = "P", n = Inf, adjust.method = "BH")
top_SST.table <- topTable(tmp_SST, sort.by = "P", n = Inf, adjust.method = "BH")
top_OPC.table <- topTable(tmp_OPC, sort.by = "P", n = Inf, adjust.method = "BH")
top_Oligodendrocyte.table <- topTable(tmp_Oligodendrocyte, sort.by = "P", n = Inf, adjust.method = "BH")
top_Astrocyte.table <- topTable(tmp_Astrocyte, sort.by = "P", n = Inf, adjust.method = "BH")

rownames(top_Excitatory.table) <- paste("Excitatory",rownames(top_Excitatory.table))
rownames(top_LAMP5.table) <- paste("LAMP5",rownames(top_LAMP5.table))
rownames(top_VIP.table) <- paste("VIP",rownames(top_VIP.table))
rownames(top_PVALB.table) <- paste("PVALB",rownames(top_PVALB.table))
rownames(top_SST.table) <- paste("SST",rownames(top_SST.table))
rownames(top_OPC.table) <- paste("OPC",rownames(top_OPC.table))
rownames(top_Oligodendrocyte.table) <- paste("Oligodendrocyte",rownames(top_Oligodendrocyte.table))
rownames(top_Astrocyte.table) <- paste("Astrocyte",rownames(top_Astrocyte.table))

#top_Excitatory.table <- cbind(top_Excitatory.table,Celltype=rep("Excitatory",nrow(top_Excitatory.table)))
#top_LAMP5.table <- cbind(top_LAMP5.table,Celltype=rep("LAMP5",nrow(top_LAMP5.table)))
#top_VIP.table <- cbind(top_VIP.table,Celltype=rep("VIP",nrow(top_VIP.table)))
#top_PVALB.table <- cbind(top_PVALB.table,Celltype=rep("PVALB",nrow(top_PVALB.table)))
#top_SST.table <- cbind(top_SST.table,Celltype=rep("SST",nrow(top_SST.table)))
#top_OPC.table <- cbind(top_OPC.table,Celltype=rep("OPC",nrow(top_OPC.table)))
#top_Oligodendrocyte.table <- cbind(top_Oligodendrocyte.table,Celltype=rep("Oligodendrocyte",nrow(top_Oligodendrocyte.table)))
#top_Astrocyte.table <- cbind(top_Astrocyte.table,Celltype=rep("Astrocyte",nrow(top_Astrocyte.table)))

combined_table <- rbind(top_Excitatory.table,
                        top_LAMP5.table,
                        top_VIP.table,
                        top_PVALB.table,
                        top_SST.table,
                        top_OPC.table,
                        top_Oligodendrocyte.table,
                        top_Astrocyte.table)

write.csv(combined_table, paste(directory_name,"/","SonnyPaper_Results.csv", sep=""), row.names=TRUE)

