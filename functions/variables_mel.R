

NK_mel <- check_genesNK(1,4, df_expanded, ex_mel, res_df)

samples_mel <- clin_mel %>% 
  filter(GEO_ID %in% colnames(ex_mel)) %>% 
  select(GEO_ID, MET_SITE)

# Melanoma
pv_mel <- data.frame(Gene = unique(NK_mel$final_symb),
                     p_value = rep(NA,length(unique(NK_mel$final_symb))))
sig_genes_mel <- sign_genes(clin_mel, ex_mel, pv_mel, 'Melanoma')
sig_df_mel <- sig_genes_mel[[2]]
sig_genes_mel <- sig_genes_mel[[1]]