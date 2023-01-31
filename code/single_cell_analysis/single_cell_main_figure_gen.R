library(broom)
library(stats)
library(ggrepel)
library(cowplot)

# compare snCTPs with rCTPs

snCTP_results_harmonized = read_csv('code/collab/data/snCTP_results_harmonized.csv')

# load in mgp results from sonny's github or local
mgp_results = read_csv('data/Misc/MGP_Results.csv')

# load in effect sizes of beta coeffs from sonny's prior analysis
bulk_tissue_age_effects = read_csv('data/output/tables/Supplementary_Table_2.csv')
bulk_tissue_age_effects = bulk_tissue_age_effects %>% dplyr::rename(cell_type = 'Cell type', bulk_age_beta = `Std. beta coefficient`)


# do stats modelling to estimate std. beta coeffs based on snCTPs - note that this is doing a mega-analysis of the 2 cohorts

snCTPs_results_harmonized_long = snCTP_results_harmonized %>% pivot_longer(cols = Astrocyte:VIP, names_to = 'cell_type', values_to = 'snCTP')

cell_type_list = unique(snCTPs_results_harmonized_long$cell_type)

snCTP_model_effects <- snCTPs_results_harmonized_long %>% 
  
  filter(total_cell_counts > 500) %>% 
  distinct(cell_type, Institution, subject_id, .keep_all = T) %>%
  
  # group aging_paper_label stacked data by cell_type
  group_by(cell_type) %>%
  
  # fit all the cell_type_prop data accorting to the model 
  # using the broom package to tidy the results 
  do(tidy(lm(scale(snCTP) ~ scale(ageOfDeath)  +  
               scale(PMI) + 
               Sex + Institution,  data = .))) %>%
  
  # unstack the data and adjust for multiple comparisons using the Benjamini-Hochberg method
  ungroup() %>%
  mutate(padj = p.adjust(`p.value`, method = 'BH')) %>%
  
  # clean up the names the the term column
  mutate(term = recode(term,
                       `(Intercept)` = "Intercept",
                       `SexM` = "gender:Male",
                       `DxSCZ` = "SCZ",
                       `DxBP` = "BP",
                       `DxMDD` = "MDD",
                       `EthnicityAfrican-American` = "ethnicity:Black", ## Q in the model ???
                       `scale(ageOfDeath)` = "age",
                       `scale(PMI)` = "PMI"))


# generate bar graph with error bars showing std. beta coeffs of age effects based on snCTP mega-analysis
beta_coeff_plot = snCTP_model_effects  %>% filter(term == "age") %>% ggplot(aes(x = cell_type, y = estimate)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) + 
  ylab('std. Beta coeff') + theme_cowplot() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) 


# generate examples showing age associations in snCTPs in both cohorts
snCTP_examples = snCTPs_results_harmonized_long %>% filter(cell_type %in% c('SST', 'Excitatory', 'Astrocyte')) %>% 
  mutate(cell_type = factor(cell_type, levels = c('SST', 'Excitatory', 'Astrocyte')), 
         Institution = factor(Institution, levels = c('MSSM', 'McLean'))) %>%
  distinct(subject_id, cell_type, .keep_all = T) %>% 
  ggplot(aes(x = ageOfDeath, y = snCTP * 100), alpha = .5) + 
  theme_cowplot() + 
  geom_smooth(method = "lm", se = F, color = 'grey') + 
  geom_point(alpha = 0.75) + 
  ylab('Cell type proportion (snCTP, %)') + 
  xlab('Age of death (years)') + 
  facet_grid(rows = vars(cell_type), cols = vars(Institution), scales = "free_y")


# look at consistency of age-effects between bulk and snCTP analyses
single_nucleus_age_effects =  snCTP_model_effects %>% filter(term == "age") %>% dplyr::rename(single_nucleus_age_beta = estimate)
joined_bulk_sn_effects = inner_join(bulk_tissue_age_effects, single_nucleus_age_effects, by = 'cell_type')


bulk_vs_sn_effects_plot = joined_bulk_sn_effects %>% 
  ggplot(aes(x = bulk_age_beta, y = single_nucleus_age_beta, label = cell_type)) + 
  geom_vline(xintercept = 0, alpha = .25) + 
  geom_hline(yintercept = 0, alpha = .25) + 
  geom_abline(slope = 1, intercept = 0, alpha = .25) + 
  geom_point() + 
  geom_text_repel() + 
  xlab('Age effects, bulk tissue (std. Beta)') + 
  ylab('Age effects, single-nucleus (std. Beta)') + 
  theme_cowplot()


# generate final plots 
top_plot_half = plot_grid(snCTP_examples, bulk_vs_sn_effects_plot, ncol = 2, rel_widths = c(1.25, 1))


final_snCTP_plot = plot_grid(top_plot_half, beta_coeff_plot, nrow = 2, rel_heights = c(1, 1))
ggsave(final_snCTP_plot, filename = 'code/collab/graphs/final_snCTP_plot.pdf', width = 8, height = 8)

### figure out consistency between rCTPs (using MGP), abCTPs (using dtangle) and snCTPs for subjects with both in the MSSM cohort



rCTPs_long = mgp_results %>% filter(Institution == 'MSSM') %>% select(subject_id, Oligodendrocyte:PVALB) %>% pivot_longer(cols = Oligodendrocyte:PVALB, names_to = 'cell_type', values_to = 'rCTP')

mgp_snctp_joined = inner_join(snCTPs_results_harmonized_long %>% filter(total_cell_counts > 500, Institution == 'MSSM'), 
                              rCTPs_long, by = c('subject_id', 'cell_type')) 

# these are the mssm based estimates we made in dtangle_aibs_scratch.R
dtangle_estimates = read_csv('code/collab/data/mssm_dtangle_estimates.csv')

dtangle_estimates_long = dtangle_estimates %>% dplyr::rename(subject_id = X1) %>% 
  pivot_longer(cols = Oligodendrocyte:PAX6, names_to = 'cell_type', values_to = 'abCTP')

mgp_snctp_joined = inner_join(mgp_snctp_joined, 
                              dtangle_estimates_long, by = c('subject_id', 'cell_type')) 
  

mgp_vs_snctp_consistency_plot = mgp_snctp_joined %>% filter(cell_type %in% c('SST', 'Excitatory', 'Astrocyte')) %>% 
  mutate(cell_type = factor(cell_type, levels = c('SST', 'Excitatory', 'Astrocyte')), 
         Institution = factor(Institution, levels = c('MSSM', 'McLean'))) %>%
  distinct(subject_id, cell_type, .keep_all = T) %>% 
  ggplot(aes(y = rCTP, x = snCTP * 100), alpha = .5) + 
  theme_cowplot() + 
  geom_smooth(method = "lm", se = F, color = 'blue') + 
  geom_point(alpha = 0.75) + 
  xlab('Cell type proportion (single-nucleus, %)') + 
  ylab('Rel. proportion (MGP, AU)') + 
  facet_wrap(~cell_type, scales = "free")

abctp_vs_snctp_consistency_plot = mgp_snctp_joined %>% filter(cell_type %in% c('SST', 'Excitatory', 'Astrocyte')) %>% 
  mutate(cell_type = factor(cell_type, levels = c('SST', 'Excitatory', 'Astrocyte')), 
         Institution = factor(Institution, levels = c('MSSM', 'McLean'))) %>%
  distinct(subject_id, cell_type, .keep_all = T) %>% 
  ggplot(aes(y = abCTP* 100, x = snCTP * 100), alpha = .5) + 
  theme_cowplot() + 
  geom_abline(slope = 1, intercept = 0, color = 'grey') + 
  geom_smooth(method = "lm", se = F, color = 'blue') + 
  geom_point(alpha = 0.75) + 
  xlab('Cell type proportion (single-nucleus, %)') + 
  ylab('Abs. proportion (dtangle, %)') + 
  facet_wrap(~cell_type, scales = "free")


consistency_corr_df = mgp_snctp_joined %>% group_by(cell_type) %>% summarize(rctp_snctp_corr = cor(snCTP, rCTP, method = 'spearman'), 
                                                                             abctp_snctp_corr = cor(snCTP, abCTP, method = 'spearman'))

consistency_corr_df_long = consistency_corr_df %>% 
  pivot_longer(cols = rctp_snctp_corr:abctp_snctp_corr, names_to = 'corr_type', values_to = 'corr_value') %>%
  mutate(corr_algo = case_when(corr_type == 'rctp_snctp_corr' ~ 'MGP', 
                               corr_type == 'abctp_snctp_corr' ~ 'dtangle'), 
         corr_algo = factor(corr_algo, levels = c('MGP', 'dtangle')))

mgp_snctp_bar_graph = consistency_corr_df_long %>% ggplot(aes(x = cell_type, y = corr_value, fill = corr_algo)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_cowplot() + 
  ylab('bulk vs single-nucleus CTPs (Spearman)') + 
  xlab('cell type') 

final_consistency_plot = plot_grid(mgp_vs_snctp_consistency_plot, abctp_vs_snctp_consistency_plot, mgp_snctp_bar_graph, nrow = 3)
ggsave(final_consistency_plot, filename = 'code/collab/graphs/final_consistency_plot.pdf', width = 8, height = 12)


