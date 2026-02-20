### 
###
### Supplementary Figure 5
###
### 


load('data/datasets_lung_soil.RData')
source("functions.R")
source("variables.R")
source("variables_pancan.R")
source("variables_all.R")
source("variables_crc.R")
source("variables_mel.R")


# Filter expression data only with the NK genes
ex_filt_crc <- as.data.frame(ex_crc[rownames(ex_crc) %in% NK_crc$final_symb,])

ex_filt_mel <- as.data.frame(ex_mel[rownames(ex_mel) %in% NK_mel$final_symb,
                                    colnames(ex_mel) %in% samples_mel$GEO_ID])

ex_filt_pancan <- as.data.frame(ex_pancan[rownames(ex_pancan) %in% NK_pancan$final_symb,
                                          colnames(ex_pancan) %in% samples_pancan$SAMPLE_ID])

ex_filt_all <- as.data.frame(ex_all[rownames(ex_all) %in% NK_all$final_symb,
                                    colnames(ex_all) %in% samples_all$GEO_ID])


# Change some genes by their synonym in the pancreas dataset to have these ones
# in common with the other datasets
## NECTIN2 ---> PVRL2
## NECTIN3 ---> PVRL3
## NECTIN1 ---> PVRL1

for(i in seq_along(rownames(ex_filt_pancan))){
  if(rownames(ex_filt_pancan)[i] == 'NECTIN2'){
    rownames(ex_filt_pancan)[i] <- 'PVRL2'
  } else if(rownames(ex_filt_pancan)[i] == 'NECTIN3'){
    rownames(ex_filt_pancan)[i] <- 'PVRL3'
  } else if (rownames(ex_filt_pancan)[i] == 'NECTIN1'){
    rownames(ex_filt_pancan)[i] <- 'PVRL1'
  }
}

# Also in significant gene vector in pancan
for (i in seq_along(sig_genes_pancan)){
  if(sig_genes_pancan[i] == 'NECTIN2'){
    sig_genes_pancan[i] <- 'PVRL2'
  } else if(sig_genes_pancan[i] == 'NECTIN3'){
    sig_genes_pancan[i] <- 'PVRL3'
  } else if (sig_genes_pancan[i] == 'NECTIN1'){
    sig_genes_pancan[i] <- 'PVRL1'
  }
}

# Now we have to add not shared genes with NA value
all_genes <- union(rownames(ex_filt_all), rownames(ex_filt_crc))
all_genes <- union(all_genes,rownames(ex_filt_pancan))
all_genes <- union(all_genes,rownames(ex_filt_mel))

# Add the missed ones in each expression matrix
for(i in seq_along(all_genes)){
  g <- all_genes[i]
  
  if(!g %in% rownames(ex_filt_all)){
    new_row <- as.list(rep(NA, ncol(ex_filt_all)))
    ex_filt_all <- rbind(ex_filt_all, new_row)
    rownames(ex_filt_all)[nrow(ex_filt_all)] <- g
  }
  
  if(!g %in% rownames(ex_filt_crc)){
    new_row <- as.list(rep(NA, ncol(ex_filt_crc)))
    ex_filt_crc <- rbind(ex_filt_crc, new_row)
    rownames(ex_filt_crc)[nrow(ex_filt_crc)] <- g
  }
  
  if(!g %in% rownames(ex_filt_mel)){
    new_row <- as.list(rep(NA, ncol(ex_filt_mel)))
    ex_filt_mel <- rbind(ex_filt_mel, new_row)
    rownames(ex_filt_mel)[nrow(ex_filt_mel)] <- g
  }
  
  if(!g %in% rownames(ex_filt_pancan)){
    new_row <- as.list(rep(NA, ncol(ex_filt_pancan)))
    ex_filt_pancan <- rbind(ex_filt_pancan, new_row)
    rownames(ex_filt_pancan)[nrow(ex_filt_pancan)] <- g
  }
}

# Now they have the same genes (N = 120)
dim(ex_filt_all);
dim(ex_filt_mel);
dim(ex_filt_crc);
dim(ex_filt_pancan)


# Keep only the significant genes
# Condition: if there is any gene which isn't significant in any of the 
# dataset, remove it

total_sig_genes <- union(sig_genes_all, sig_genes_crc)
total_sig_genes <- union(total_sig_genes, sig_genes_pancan)
total_sig_genes <- union(total_sig_genes, sig_genes_mel)

# Filter expression data with these genes
ex_filt_all.sign <- ex_filt_all[rownames(ex_filt_all) %in% total_sig_genes,]
ex_filt_crc.sign <- ex_filt_crc[rownames(ex_filt_crc) %in% total_sig_genes,]
ex_filt_pancan.sign <- ex_filt_pancan[rownames(ex_filt_pancan) %in% total_sig_genes,]
ex_filt_mel.sign <- ex_filt_mel[rownames(ex_filt_mel) %in% total_sig_genes,]

dim(ex_filt_all.sign);
dim(ex_filt_mel.sign);
dim(ex_filt_crc.sign);
dim(ex_filt_pancan.sign)

# Recalculate NK variables by using the complete expression matrix
NK_melan_filt <- check_genesNK(1,4,df_expanded, ex_filt_mel.sign, res_df)
NK_crc_filt <- check_genesNK(1,4,df_expanded, ex_filt_crc.sign, res_df)
NK_pancan_filt <- check_genesNK(1,4,df_expanded, ex_filt_pancan.sign, res_df)
NK_all_filt <- check_genesNK(1,4,df_expanded, ex_filt_all.sign, res_df)


NK_all_filt$geneset <- as.factor(NK_all_filt$geneset)
split_all_filt <- split(NK_all_filt, NK_all_filt$geneset)

NK_melan_filt$geneset <- as.factor(NK_melan_filt$geneset)
split_melan_filt <- split(NK_melan_filt, NK_melan_filt$geneset)

NK_pancan_filt$geneset <- as.factor(NK_pancan_filt$geneset)
split_pancan_filt <- split(NK_pancan_filt, NK_pancan_filt$geneset)

NK_crc_filt$geneset <- as.factor(NK_crc_filt$geneset)
split_crc_filt <- split(NK_crc_filt, NK_crc_filt$geneset)


# Obtain scaled z-score for each dataset by using the original
# (which contain all genes)

## CRC
geneset_names <- names(split_crc_filt)
scaled_crc <- lapply(geneset_names, function(g) median_zs_data(ex_filt_crc.sign, 
                                                               samples_crc, 
                                                               'CRC',
                                                               sig_genes_crc,
                                                               g,
                                                               split_crc_filt))
scaled_crc <- setNames(scaled_crc, geneset_names)

## Melanoma
geneset_names <- names(split_melan_filt)
scaled_mel <- lapply(geneset_names, function(g) median_zs_data(ex_filt_mel.sign, 
                                                               samples_mel,
                                                               'Melanoma',
                                                               sig_genes_mel, 
                                                               g, 
                                                               split_melan_filt))
scaled_mel <- setNames(scaled_mel, geneset_names)

## All
geneset_names <- names(split_all_filt)
scaled_all <- lapply(geneset_names, function(g) median_zs_data(ex_filt_all.sign, 
                                                               samples_all,
                                                               'All',
                                                               sig_genes_all, 
                                                               g,
                                                               split_all_filt))
scaled_all <- setNames(scaled_all, geneset_names)

## Pancreas
geneset_names <- names(split_pancan_filt)
scaled_pancan <- lapply(geneset_names, 
                        function(g) median_zs_data(ex_filt_pancan.sign,
                                                   samples_pancan,
                                                   'Pancreas',
                                                   sig_genes_pancan, 
                                                   g,
                                                   split_pancan_filt))
scaled_pancan <- setNames(scaled_pancan, geneset_names)


# Prepare data
all_datav1 <- lapply(geneset_names, function(g) prepare_data(scaled_all,
                                                             scaled_crc,
                                                             scaled_pancan,
                                                             scaled_mel, g))
all_datav1 <- setNames(all_datav1, geneset_names)

# Set color palette (valid for all gene sets)
library(RColorBrewer)
color_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

# Set breaks range for each group
breaks_range <- lapply(geneset_names, function(g) define_breaks(all_datav1,g))
breaks_range <- setNames(breaks_range, geneset_names)


# Heatmaps

# All
heatmap_all <- lapply(geneset_names, function(g) heatmap_zscore_sign(ex_filt_all.sign, 
                                                                     samples_all,
                                                                     sig_genes_all,
                                                                     g,
                                                                     split_all_filt, 'All',
                                                                     color_palette,
                                                                     breaks_range, FALSE))
heatmap_all <- setNames(heatmap_all, geneset_names)

# CRC
heatmap_crc <- lapply(geneset_names, function(g) heatmap_zscore_sign(ex_filt_crc.sign, 
                                                                     samples_crc,
                                                                     sig_genes_crc,
                                                                     g,
                                                                     split_crc_filt, 'CRC',
                                                                     color_palette,
                                                                     breaks_range, FALSE))
heatmap_crc <- setNames(heatmap_crc, geneset_names)


# Melanoma
heatmap_melan <- lapply(geneset_names, function(g) heatmap_zscore_sign(ex_filt_mel.sign, 
                                                                       samples_mel,
                                                                       sig_genes_mel,
                                                                       g,
                                                                       split_melan_filt, 
                                                                       'Melanoma',
                                                                       color_palette,
                                                                       breaks_range, TRUE))
heatmap_melan <- setNames(heatmap_melan, geneset_names)


# Pancreas
heatmap_pancan <- lapply(geneset_names, function(g) heatmap_zscore_sign(ex_filt_pancan.sign,
                                                                        samples_pancan,
                                                                        sig_genes_pancan,
                                                                        g,
                                                                        split_pancan_filt,
                                                                        'Pancreas',
                                                                        color_palette,
                                                                        breaks_range, FALSE))
heatmap_pancan <- setNames(heatmap_pancan, geneset_names)


# Generate figure
res.suppl.fig5 <- supp.figure5(heatmap_all, heatmap_crc, heatmap_pancan, heatmap_melan, geneset_names)