

NK_pancan <- check_genesNK(1,4, df_expanded, ex_pancan, res_df)

samples_pancan <- clin_pancan %>% 
  filter(SAMPLE_ID %in% colnames(ex_pancan)) %>% 
  select(SAMPLE_ID, TUMOR_SITE)

pv_pancan <- data.frame(Gene = unique(NK_pancan$final_symb),
                        p_value = rep(NA,length(unique(NK_pancan$final_symb))))
sig_genes_pancan <- sign_genes(clin_pancan, ex_pancan, pv_pancan, 'Pancreas')
sig_df_pancan <- sig_genes_pancan[[2]]
sig_genes_pancan <- sig_genes_pancan[[1]]
