

NK_all <- check_genesNK(1,4, df_expanded, ex_all, res_df)

samples_all <- clin_all %>% 
  filter(GEO_ID %in% colnames(ex_all)) %>% 
  select(GEO_ID, MET_SITE)


# All 
pv_all <- data.frame(Gene = unique(NK_all$final_symb),
                     p_value = rep(NA,length(unique(NK_all$final_symb))))
sig_genes_all <- sign_genes_all(clin_all, ex_all, pv_all)
sig_df_all <- sig_genes_all[[2]]
sig_genes_all <- sig_genes_all[[1]]

sig_genes_all <- sort(sig_genes_all)