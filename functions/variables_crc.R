

NK_crc <- check_genesNK(1,4, df_expanded, ex_crc, res_df)

samples_crc <- clin_crc %>% 
  filter(GEO_ID %in% colnames(ex_crc)) %>% 
  select(GEO_ID, MET_SITE)

# CRC
pv_crc <- data.frame(Gene = unique(NK_crc$final_symb),
                     p_value = rep(NA,length(unique(NK_crc$final_symb))))
sig_genes_crc <- sign_genes(clin_crc, ex_crc, pv_crc, 'CRC')
sig_df_crc <- sig_genes_crc[[2]]
sig_genes_crc <- sig_genes_crc[[1]]

sig_genes_crc <- sort(sig_genes_crc)