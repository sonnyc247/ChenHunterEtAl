#### code required to make the graphs in the paper ####

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggpubr)

#### figure 1 age distribution histogram ####

fig1_df <- mgp_results_total_df[,c("ageOfDeath", "dataset")]
fig1_df$dataset <- as.character(fig1_df$dataset)
fig1_df <- fig1_df[fig1_df$dataset != "Pitt BA47 Microarray",]
fig1_df[fig1_df$dataset == "Pitt BA11 Microarray", "dataset"] <- "Pitt BA11 and BA47 Microarray"
fig1_df_total <- fig1_df
fig1_df_total$dataset <- "Combined datasets"
fig1_df <- rbind(fig1_df, fig1_df_total)
remove(fig1_df_total)

fig1_df$dataset <- factor(fig1_df$dataset, levels = c("MSSM RNAseq",
                                                      "Penn RNAseq",
                                                      "Pitt RNAseq",
                                                      "NIMH RNAseq",
                                                      "Pitt BA11 and BA47 Microarray",
                                                      "NABEC Microarray",
                                                      "UKBEC Microarray",
                                                      "Combined datasets"))

ggplot(fig1_df, aes(x = ageOfDeath, fill = dataset)) + 
          geom_histogram(binwidth = 3) + 
          geom_hline(yintercept=-0.1, size = 1.5) + 
          scale_y_continuous(expand = c(0,0)) + 
          facet_wrap(~dataset, ncol = 1, scales = "free_y") +
          scale_fill_manual(values = c("#F4BA4C", "#F2875D", "#E66497", "#7846D9", "blue", "#45ADB4", "#57C555", "gray")) +
          xlab('Age of subjects') + 
          ylab('Number of subjects') + 
          theme(legend.position = "none", 
                text = element_text(size=14), 
                axis.text = element_text(size = 10), 
                strip.background = element_rect(fill= "white"),
                strip.text = element_text(colour = 'black'))

ggsave(width = 180,
       height = 200,
       dpi = 300, 
       units = "mm", 
       bg = "white",
       path = "/external/rprshnas01/kcni/ychen/git/Aging_MGP_private/output/",
       filename = "Figure_1.pdf",
       device = "pdf")

#### figure 2 illustrative MGP scatter plots ####

mgp_results_total_df <- read.csv("./data/output/tables/Allregion_MGP_Results.csv", row.names = 1)

figure_2_df <- mgp_results_total_df[,c("Astrocyte", "SST", "Excitatory", "ageOfDeath", "dataset")]
#figure_2_df <- mgp_results_total_df

#V2 plot:

# rename
names(figure_2_df)[5] <- "Dataset"
#names(figure_2_df)[19] <- "Dataset"

figure_2_df$Dataset <- factor(figure_2_df$Dataset, levels = c("MSSM RNAseq",
                                                              "Penn RNAseq",
                                                              "Pitt RNAseq",
                                                              "NIMH RNAseq",
                                                              "Pitt BA11 Microarray",
                                                              "Pitt BA47 Microarray",
                                                              "NABEC Microarray",
                                                              "UKBEC Microarray"))

gathered_figure_2_df <- tidyr::gather(figure_2_df, key = "Celltype", value = "Scaled_MGP", Excitatory, SST, Astrocyte)
#gathered_figure_2_df <- tidyr::gather(figure_2_df, key = "Celltype", value = "Scaled_MGP", 2:14)

gathered_figure_2_df$Celltype <- factor(gathered_figure_2_df$Celltype, levels = c("SST",
                                                                                  "Excitatory",
                                                                                  "Astrocyte"))

ggplot(data = gathered_figure_2_df, aes(x = ageOfDeath, y = Scaled_MGP)) + 
  geom_point(size = 0.35, alpha = 0.3) + 
  geom_smooth(method = 'lm', se = F, size = 0.45) +
  #geom_smooth(method = "loess", se = F, size = 0.45) +
  facet_grid(Dataset ~ Celltype, switch = "y", scale="free_y") +
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
       height = 250,
       dpi = 300, 
       units = "mm",
       bg = "white",
       path = "./data/output/figures/",
       filename = "Figure_2_AllBrainMarkers.pdf",
       device = "pdf")

# conf interval values
cor.test(figure_2_df[figure_2_df$Dataset == "UKBEC Microarray", "SST"],
         figure_2_df[figure_2_df$Dataset == "UKBEC Microarray", "ageOfDeath"])
cor.test(figure_2_df[figure_2_df$Dataset == "NIMH RNAseq", "SST"],
         figure_2_df[figure_2_df$Dataset == "NIMH RNAseq", "ageOfDeath"])

cor.test(figure_2_df[figure_2_df$Dataset == "UKBEC Microarray", "Excitatory"],
         figure_2_df[figure_2_df$Dataset == "UKBEC Microarray", "ageOfDeath"])
cor.test(figure_2_df[figure_2_df$Dataset == "NIMH RNAseq", "Excitatory"],
         figure_2_df[figure_2_df$Dataset == "NIMH RNAseq", "ageOfDeath"])

cor.test(figure_2_df[figure_2_df$Dataset == "Pitt RNAseq", "Astrocyte"],
         figure_2_df[figure_2_df$Dataset == "Pitt RNAseq", "ageOfDeath"])
cor.test(figure_2_df[figure_2_df$Dataset == "NABEC Microarray", "Astrocyte"],
         figure_2_df[figure_2_df$Dataset == "NABEC Microarray", "ageOfDeath"])

#### figure 3 mega analysis plot ####

term_of_interest <- "age_scaled"
#term_of_interest <- "sexM" ## if desired for plot
#term_of_interest <- "pmi_scaled" ## if desired for plot

mega_results$class <- "Neuronal"
mega_results[mega_results$celltype %in% c("OPC",
                                          "Pericyte",
                                          "Endothelial",
                                          "Oligodendrocyte",
                                          "VLMC",
                                          "Astrocyte") , "class"] <- "Non-Neuronal"

#plotting beta values from mega analysis

mega_results$celltype <- factor(mega_results$celltype, levels = c("SST",
                                                                  "VIP",
                                                                  "PAX6",
                                                                  "LAMP5",
                                                                  "PVALB",
                                                                  "Excitatory",
                                                                  "OPC",
                                                                  "Oligodendrocyte",
                                                                  "Endothelial",
                                                                  "Astrocyte",
                                                                  "VLMC",
                                                                  "Pericyte"))

ggplot(mega_results[mega_results$term == term_of_interest,], aes(x = celltype)) + 
  geom_bar(stat = 'identity', width = 0.8, aes(y = estimate, fill = estimate), colour = "black") +
  geom_errorbar((aes(ymin = estimate - std.error, ymax = estimate + std.error))) +
  scale_fill_gradient2(low = "blue", high = "red") +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  facet_grid(~ class, scales = "free_x", space = "free") +
  xlab('Cell types') + 
  ylab('Mega analysis std. beta') +
  theme_classic() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 20, hjust = 0.95, vjust = 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 11),
        text = element_text(size = 10), 
        axis.text = element_text(size = 7.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) 

ggsave(width = 180, 
       dpi = 300, 
       units = "mm",
       bg = "white",
       path = "./data/output/figures/",
       filename = "Figure_3_AllBrainMarkers.pdf",
       device = "pdf")

#### figure 4 individual analyses plot ####

# pre-process (add dummy group for spacing)
fig4_df <- indiv_age_results
#fig4_df[nrow(fig4_df)+1,] <- fig4_df[1,]
#fig4_df[nrow(fig4_df),3:6] <- 0

#fig4_df$`Cell type` <- as.character(fig4_df$`Cell type`)
#fig4_df[nrow(fig4_df),"Cell type"] <- "remove"

#set order of datasets and cell types
fig4_df$Dataset <- factor(fig4_df$Dataset, levels = c("MSSM RNAseq",
                                                      "Penn RNAseq",
                                                      "Pitt RNAseq",
                                                      "NIMH RNAseq",
                                                      "Pitt BA11 Microarray",
                                                      "Pitt BA47 Microarray",
                                                      "NABEC Microarray",
                                                      "UKBEC Microarray"))

fig4_df$`Cell type` <- factor(fig4_df$`Cell type`, levels = c("SST",
                                                              "VIP",
                                                              "PAX6",
                                                              "LAMP5",
                                                              "PVALB",
                                                              "Excitatory",
                                                              "OPC",
                                                              "Oligodendrocyte",
                                                              "Endothelial",
                                                              "Astrocyte",
                                                              "VLMC",
                                                              "Pericyte"))

colnames(fig4_df)[1:4] <- c("dataset", "celltype", "estimate", "std.error")

ggplot(fig4_df, aes(x = dataset)) + 
  geom_bar(stat = 'identity', width = 0.8, aes(y = estimate, fill = dataset)) + 
  facet_wrap(~celltype, nrow = 2, scales = "free_x") + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.3) +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  ylab('Age std. beta') +
  labs(fill = "Dataset:") +
  scale_fill_manual(values = c("#F4BA4C", "#F2875D", "#E66497", "#7846D9", "#4578BF", "#4E52F2", "#45ADB4", "#57C555")) +
  theme(text = element_text(size=8), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        #legend.margin = margin(t = -15, r = 0, b = 0, l = -25, unit = "pt"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.25),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(size = 8.45),
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black', size = 8.45),
        panel.spacing = unit(0, "lines"),
        legend.title=element_blank(),
        legend.text=element_text(size=8.45)) + 
  guides(fill = guide_legend(nrow = 2, 
                             byrow = TRUE))

ggsave(width = 150, 
       dpi = 300, 
       units = "mm", 
       bg = "white",
       path = "./data/output/figures/",
       filename = "Figure_4_AllBrainMarkers.pdf",
       device = "pdf")

#
#### supplemental figures ####

### separate and combined

supplemental_data_for_plot <- readRDS("~/collabgit/Aging-Expression-Analysis/Figures: Code and Data/gene_betas_Ex.rds")
supplemental_data_for_plot <- readRDS("~/collabgit/Aging-Expression-Analysis/Figures: Code and Data/gene_betas_SST.rds")
supplemental_data_for_plot <- readRDS("~/collabgit/Aging-Expression-Analysis/Figures: Code and Data/gene_betas_Oligo.rds")
ref <- readRDS("~/collabgit/Aging-Expression-Analysis/Figures: Code and Data/gene_betas_Oligo_old.rds")
ref <- unique(ref$gene)

supplemental_data_for_plot[supplemental_data_for_plot$dataset == "pitt_microarray_BA11", "dataset"] <- "Pitt BA11 Microarray"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "pitt_microarray_BA47", "dataset"] <- "Pitt BA47 Microarray"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "MSSM", "dataset"] <- "MSSM RNAseq"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "penn", "dataset"] <- "Penn RNAseq"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "pitt_RNAseq", "dataset"] <- "Pitt RNAseq"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "UKBEC", "dataset"] <- "UKBEC Microarray"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "NABEC", "dataset"] <- "NABEC Microarray"

supplemental_data_for_plot$dataset <- factor(supplemental_data_for_plot$dataset, 
                                             levels = c("MSSM RNAseq",
                                                        "Penn RNAseq",
                                                        "Pitt RNAseq",
                                                        "Pitt BA11 Microarray",
                                                        "Pitt BA47 Microarray",
                                                        "NABEC Microarray",
                                                        "UKBEC Microarray"))

supplemental_data_for_plot <- supplemental_data_for_plot[supplemental_data_for_plot$gene %in% ref,]

oligo_plot <- ggplot(data = supplemental_data_for_plot) +
  geom_tile(aes(x = dataset, y = gene, fill = beta)) +
  scale_fill_gradientn(limits = c(-0.046,0.046),
                       colours = c("blue", "white", "red"),
                       breaks = c(-0.04,0,0.04),
                       labels = format(c(-0.04,0,0.04))) +
  #scale_fill_gradient(low = 'blue', high = 'red') +
  labs(x = "Dataset", y = "Marker gene", fill = "Std. beta") +
  theme_classic() +
  #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.background = element_rect(fill = 'grey'))
  theme(axis.text.x = element_text(angle = 20, vjust = 0.95, hjust = 0.99), panel.background = element_rect(fill = 'grey'))


ggarrange(ex_plot,
          sst_plot,
          oligo_plot,
          nrow = 3,
          align = "v",
          labels = c("a)", "b)", "c)"),
          font.label = list(size = 12, face = "plain"))

ggsave(width = 140,
       height = 280,
       dpi = 300, 
       units = "mm",
       bg = "white",
       path = "/external/rprshnas01/kcni/ychen/collabgit/Aging-Expression-Analysis/Figures: Code and Data/",
       filename = "SupFig1_Sep.pdf",
       device = "pdf")

### all in one

supplemental_data_for_plot <- readRDS("~/collabgit/Aging-Expression-Analysis/Figures: Code and Data/gene_betas_Ex.rds")
SST_data <- readRDS("~/collabgit/Aging-Expression-Analysis/Figures: Code and Data/gene_betas_SST.rds")
Oligo_data <- readRDS("~/collabgit/Aging-Expression-Analysis/Figures: Code and Data/gene_betas_Oligo.rds")

supplemental_data_for_plot$celltype <- "Excitatory"
SST_data$celltype <- "SST"
Oligo_data$celltype <- "Oligodendrocyte"

supplemental_data_for_plot <- rbind(supplemental_data_for_plot, SST_data, Oligo_data)
remove(SST_data)
remove(Oligo_data)

supplemental_data_for_plot[supplemental_data_for_plot$dataset == "pitt_microarray_BA11", "dataset"] <- "Pitt BA11 Microarray"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "pitt_microarray_BA47", "dataset"] <- "Pitt BA47 Microarray"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "MSSM", "dataset"] <- "MSSM RNAseq"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "penn", "dataset"] <- "Penn RNAseq"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "pitt_RNAseq", "dataset"] <- "Pitt RNAseq"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "UKBEC", "dataset"] <- "UKBEC Microarray"
supplemental_data_for_plot[supplemental_data_for_plot$dataset == "NABEC", "dataset"] <- "NABEC Microarray"

supplemental_data_for_plot$dataset <- factor(supplemental_data_for_plot$dataset, 
                                             levels = c("MSSM RNAseq",
                                                        "Penn RNAseq",
                                                        "Pitt RNAseq",
                                                        "Pitt BA11 Microarray",
                                                        "Pitt BA47 Microarray",
                                                        "NABEC Microarray",
                                                        "UKBEC Microarray"))

supplemental_data_for_plot$celltype <- factor(supplemental_data_for_plot$celltype, 
                                              levels = c("Excitatory",
                                                         "SST",
                                                         "Oligodendrocyte"))

ggplot(data = supplemental_data_for_plot) +
  geom_tile(aes(x = dataset, y = gene, fill = beta)) +
  scale_fill_gradient(low = 'red', high = 'blue') +
  facet_wrap(~celltype, ncol = 1, scales = "free") +
  labs(x = "Dataset", y = "Marker gene", fill = "Std. beta") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1))

ggsave(width = 100,
       height = 400,
       dpi = 300, 
       units = "mm",
       bg = "white",
       path = "/external/rprshnas01/kcni/ychen/collabgit/Aging-Expression-Analysis/Figures: Code and Data/",
       filename = "SupFig1.pdf",
       device = "pdf")
