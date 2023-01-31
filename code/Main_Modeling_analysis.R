#### Packages ####

library(ggplot2)
library(lmerTest)
library(dplyr)

#### Mega analysis ####

mgp_results_total_df <- read.csv("./data/output/tables/Allregion_MGP_Results.csv", row.names = 1)
mgp_results_total_df$Sex <- as.factor(mgp_results_total_df$Sex)
mgp_results_total_df <- mgp_results_total_df[-5]
celltypes <- colnames(mgp_results_total_df)[2:13]

mega_results = data.frame(term = character(), 
                          estimate = numeric(),
                          std.error = numeric(),
                          p.value = numeric(), 
                          celltype = character())

for (cell in celltypes) {
  formula <- paste(cell, '~ 1 + age_scaled + pmi_scaled  + Sex + (1|subject_id) + Institution', sep = ' ')
  model <- lmer(data = mgp_results_total_df, formula = formula, REML = T)
  model_df <- broom.mixed::tidy(model)
  model_df <- model_df[2:4,c(3:5,8)] 
  model_df$celltype <- cell
  mega_results <- rbind(mega_results, model_df)
}

mega_results$p.value.adjusted <- stats::p.adjust(mega_results$p.value, method = "BH")
mega_results_age <- mega_results[mega_results$term == "age_scaled",c(5,2:4,6)]
mega_results_sex <- mega_results[mega_results$term == "SexM",c(5,2:4,6)]

### save results

names(mega_results_age) <- c("Cell type",
                             "Std. beta coefficient",
                             "Std. error",
                             "p value",
                             "FDR")

write.csv(mega_results_age, "./data/output/tables/Supplementary_Table_2_AllBrainMarkers.csv", row.names = F)

names(mega_results_sex) <- c("Cell type",
                             "Std. beta coefficient",
                             "Std. error",
                             "p value",
                             "FDR")

write.csv(mega_results_sex, "./data/output/tables/Supplementary_Table_3_AllBrainMarkers.csv", row.names = F)

#
#### Individual dataset analysis #### 

#mgp_results_total_df$Sex <- as.factor(mgp_results_total_df$Sex) ## as needed
#celltypes <- colnames(mgp_results_total_df)[2:14] ## as needed
datasets <- unique(mgp_results_total_df$dataset) 
indiv_results <- data.frame(term = character(), 
                           estimate = numeric(),
                           std.error = numeric(),
                           p.value = numeric(), 
                           celltype = character(),
                           dataset = character())

for (cell in celltypes){
  for (dataset in datasets){
    test_df <- mgp_results_total_df[,c(cell, "age_scaled", "Sex", "pmi_scaled", "dataset")]
    test_df <- test_df[test_df$dataset == dataset,]
    formula = paste(cell, '~ age_scaled + Sex + pmi_scaled', sep = ' ')
    model <- lm(formula, test_df)
    model_df <- broom::tidy(model)
    model_df <- model_df[c(2:4),] 
    model_df <- model_df[,-4]
    model_df$celltype <- cell
    model_df$dataset <- dataset
    indiv_results <- rbind(indiv_results, model_df)
  }
  
}

indiv_age_results <- indiv_results[indiv_results$term == "age_scaled",]
indiv_age_results$p.value.adjusted <- stats::p.adjust(indiv_age_results$p.value, method = "BH")
indiv_age_results <- indiv_age_results[,c(6,5,2:4,7)]

### save results

names(indiv_age_results) <- c("Dataset",
                              "Cell type",
                              "Std. beta coefficient",
                              "Std. error",
                              "p value",
                              "FDR")

write.csv(indiv_age_results, "./data/output/tables/Supplementary_Table_4_AllBrainMarkers.csv", row.names = F)

