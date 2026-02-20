### 
###
### Figure 1
###
### 

load('data/datasets_lung_soil.RData')
source("functions.R")
source("variables.R")
source("variables_all.R")

## Figure 1A

boxplot_axis <- function(clin_data, exp_data, gene, y_axis = FALSE, y_lim = NULL){
  t <- as.data.frame(t(exp_data))
  t$GEO_ID <- rownames(t)
  
  merged_df <- merge(clin_data, t, by = "GEO_ID")
  merged_df$MET_SITE <- factor(merged_df$MET_SITE) # factor type
  
  count_df <- merged_df %>% 
    group_by(MET_SITE) %>% 
    summarise(count = n())
  
  # Keep only samples with N > 1
  filtered_sites <- count_df$count > 1
  filtered_sites <- count_df$MET_SITE[filtered_sites]
  
  # Filter original dataframe
  merged_df <- merged_df %>% 
    filter(MET_SITE %in% filtered_sites)
  
  # Clean data
  merged_df <- merged_df %>%
    filter(!is.na(merged_df[, gene])) %>%
    filter(is.finite(merged_df[, gene]))  # Remove any NA, Inf or -Inf
  
  if(any(rownames(exp_data) == gene)) {
    pos <- which(colnames(merged_df) == gene)
    data_gene <- as.numeric(exp_data[gene, ])
    
    # Normal distribution evaluation
    km <- ks.test(data_gene, "pnorm", mean = mean(data_gene), sd = sd(data_gene))
    
    if(km$p.value >= 0.05) {
      # ANOVA test
      anova_result <- aov(merged_df[, pos] ~ MET_SITE, data = merged_df)
      pv <- summary(anova_result)[[1]]$`Pr(>F)`[1]
      

      g <- ggplot(merged_df, aes(x = MET_SITE,
                                 y = .data[[gene]],
                                 color = MET_SITE)) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = "white", size = 0.8) +
        
        scale_color_manual(values = c(
          "Lung" = "green3",
          "Brain" = "maroon2",
          "Bone" = "mediumpurple3",
          "Liver" = "darkorange")) +
        
        labs(#title = gene,
          y = if (y_axis) "Enrichment score" else NULL,
          x = gene,
          color = NULL) +
        
        theme_classic(base_size = 10) +
        
        theme(
          legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_line(size = 0.3)) +
        
        annotate("text", x = Inf, y = Inf, 
                 label = paste('p = ', ifelse(pv < 1e-5, format(pv, scientific = TRUE, digits = 4),
                                              round(pv, digits = 5))),
                 hjust = 1.1, vjust = 1.3, size = 3.5)
      
      if (!is.null(y_lim)) {
        g <- g + coord_cartesian(ylim = y_lim)
      }
      
      if (!y_axis) {
        g <- g + theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank()
        )
      } else {
        g <- g + theme(
          axis.text.y = element_text(size = 8)
        )
      }
      return(g)
    } else {
      # Kruskal Wallis
      kr_result <- kruskal.test(merged_df[, pos] ~ MET_SITE, data = merged_df)
      pv <- kr_result$p.value
      
      g <- ggplot(merged_df, aes(x = MET_SITE,
                                 y = .data[[gene]],
                                 color = MET_SITE)) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = "white", size = 0.8) +
        
        scale_color_manual(values = c(
          "Lung" = "green3",
          "Brain" = "maroon2",
          "Bone" = "mediumpurple3",
          "Liver" = "darkorange")) +
        
        labs(#title = gene,
          y = if (y_axis) "Enrichment score" else NULL,
          x = gene,
          color = NULL) +
        
        theme_classic(base_size = 10) +
        
        theme(
          legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_line(size = 0.3)) +
        
        annotate("text", x = Inf, y = Inf, 
                 label = paste('p = ', ifelse(pv < 1e-5, format(pv, scientific = TRUE, digits = 4),
                                              round(pv, digits = 5))),
                 hjust = 1.1, vjust = 1.3, size = 3.5) 
      
      if (!is.null(y_lim)) {
        g <- g + coord_cartesian(ylim = y_lim)
      }
      
      if (!y_axis) {
        g <- g + theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank()
        )
      } else {
        g <- g + theme(
          axis.text.y = element_text(size = 8)
        )
      }
      
      return(g)
    }
  }
}

## Deconvolution with ConsensusTME
ConsensusTME::cancerAll
consensus <- consensusTMEAnalysis(bulkExp=as.matrix(ex_all), 
                                  cancerType = "COAD", 
                                  statMethod = "ssgsea",
                                  immuneScore=T)
rownames(consensus) <- gsub("_", " ", rownames(consensus))
rownames(consensus)

t_cells <- rownames(consensus)[which(grepl('T ', rownames(consensus)))]
ncol <- 4
t_plots <- list()
y_lim <- c(-0.5, 0.5)

for (i in seq_along(t_cells)) {
  y_axis <- ((i - 1) %% ncol == 0)
  t_plots[[i]] <- boxplot_axis(samples_all, consensus, t_cells[i], 
                     y_axis = y_axis, y_lim = y_lim)
}

library(patchwork)
combined <- wrap_plots(t_plots, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ggsave(
#   "prubas/tcell_figure.tiff",
#   combined,
#   device = "tiff",
#   width = 180,
#   height = 130,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw"
# )

ggsave(
  "figures_def2026/Fig1_A.tiff", combined, device = "tiff",
  width = 180, height = 130, units = "mm", dpi = 600, compression = "lzw"
)


## Figure 1B

nk <- boxplot_axis(samples_all, consensus,'NK cells', y_axis = TRUE)
grid.newpage()
grid.draw(nk)

# ggsave(
#   "prubas/NKcell_figure.tiff",
#   nk,
#   device = "tiff",
#   width = 88,
#   height = 70,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw"
# )

ggsave(
  "figures_def2026/Fig1_B.tiff", nk, device = "tiff",
  width = 88, height = 70, units = "mm", dpi = 600, compression = "lzw"
)



## Figure 1C

# Significative genes

pvalue_gene_all <- function(clin_data, exp_data, gene){
  
  t <- as.data.frame(t(exp_data))
  t$GEO_ID <- rownames(t)
  
  merged_df <- merge(clin_data, t, by = "GEO_ID")
  merged_df$MET_SITE <- factor(merged_df$MET_SITE) # factor type
  
  count_df <- merged_df %>% 
    group_by(MET_SITE) %>% 
    summarise(count = n())
  
  # Keep these samples with N > 1
  filtered_sites <- count_df$count > 1
  filtered_sites <- count_df$MET_SITE[filtered_sites]
  
  # Filter original dataframe
  merged_df <- merged_df %>% 
    filter(MET_SITE %in% filtered_sites)
  
  ex_filtered <- exp_data %>% select(which(colnames(exp_data) %in% merged_df$GEO_ID))
  
  if(any(rownames(ex_filtered)==gene)){
    pos <- which(colnames(merged_df)==gene)
    data_gene <- as.numeric(ex_filtered[gene, ])
    # Normal distribution evaluation
    km <- ks.test(data_gene,"pnorm", mean=mean(data_gene), sd=sd(data_gene))
    if(km$p.value >= 0.05){
      # ANOVA test
      anova_result <- aov(merged_df[,pos] ~ MET_SITE, data = merged_df)
      pv <- summary(anova_result)[[1]]$`Pr(>F)`[1]
      return(pv)
    } else{
      # Kruskal Wallis test
      kr_result <- kruskal.test(merged_df[,pos] ~ MET_SITE, data = merged_df)
      pv <- kr_result$p.value
      return(pv)
    }
  } else{
    res <- paste0(gene, ' is not in the expression matrix')
    print(res)
    return(NA)
  }
}

sign_genes_all <- function(clin_data, exp_data, df){
  for (i in seq_along(df$Gene)){
    g <- df$Gene[i]
    res <- pvalue_gene_all(clin_data, exp_data, g)
    df$p_value[i] <- res
    if(is.na(res)){
      df$is_significant[i] <- NA
    } else if(res <= 0.05){
      df$is_significant[i] <- TRUE
    } else{
      df$is_significant[i] <- FALSE
    }
  }
  sig_genes <- df %>% 
    filter(is_significant == TRUE) %>% 
    select(Gene)
  
  sig_genes <- sig_genes$Gene
  return(list(sig_genes, df))
}

# Calculate p-values to select the significative genes 
pv_all <- data.frame(Gene = unique(NK_all$final_symb),
                     p_value = rep(NA,length(unique(NK_all$final_symb))))
sig_genes_all <- sign_genes_all(clin_all, ex_all, pv_all)
sig_df_all <- sig_genes_all[[2]]
sig_genes_all <- sig_genes_all[[1]]

# Selected identificators
samples_keep <- samples_all$GEO_ID

# New expression matrix with these samples
ex_filtered <- ex_all[sig_genes_all, samples_keep ]

ex.st <- scale(t(ex_filtered), center=T, scale=T)

# GENESETS: Calculate distance between genesets and cluster distance using Ward linkage

d <- dist(ex.st, method = "euclidean") # distance matrix
hcl1 <- hclust(d, method="ward.D2") 
column_dend = as.dendrogram(hcl1)

orden_orig <- order.dendrogram(column_dend)
hcl1_reversed <- rotate(hcl1, rev(orden_orig))
column_dend = as.dendrogram(hcl1_reversed)

d <- dist(t(ex.st), method = "euclidean") # distance matrix
hcl2 <- hclust(d, method="ward.D2") 
row_dend = as.dendrogram(hcl2)

paletteLength <- 100

library(RColorBrewer)
myColor <- colorRampPalette(brewer.pal(11, "RdBu"))(paletteLength)
myColor <- rev(myColor)

myBreaks <- seq(-2, 2, length.out = paletteLength+1)
myBreaks


# Define clusters
clusters <- cutree(hcl1_reversed, k = 3)
cluster_df <- data.frame(GEO_ID = names(clusters), id_cluster = clusters)

for (i in seq_along(samples_all$GEO_ID)){
  for (k in seq_along(cluster_df$GEO_ID)){
    if(samples_all$GEO_ID[i] == cluster_df$GEO_ID[k]){
      cluster_df$MET_SITE[k] <- samples_all$MET_SITE[i]
    }
  }
}

for (k in seq_along(cluster_df$GEO_ID)){
  if(cluster_df$id_cluster[k] == 1){
    cluster_df$cluster[k] <- 'Low'
  } else if(cluster_df$id_cluster[k] == 2){
    cluster_df$cluster[k] <- 'Medium'
  } else if(cluster_df$id_cluster[k] == 3){
    cluster_df$cluster[k] <- 'High'
  }
}

annotation <- data.frame(Metastatic_site = cluster_df$MET_SITE
)
row.names(annotation) <- samples_all$GEO_ID
annotation_colors <- list(Metastatic_site = c(Liver = "darkorange", 
                                              Lung = "green3",
                                              Brain = 'maroon2' ,
                                              Bone = 'mediumpurple3')
)


p <- pheatmap(
  mat= t(ex.st),
  color = myColor,
  breaks = myBreaks,
  cluster_cols = as.hclust(hcl1_reversed),
  cluster_rows = hcl2,
  cutree_cols = 3,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 11,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  legend = TRUE,
  border_color= NA
)

ggsave(
  "figures_def2026/Fig1_C.tiff", p, device = "tiff",
  width = 180, height = 130, units = "mm", dpi = 600, compression = "lzw", bg ='white'
)


## Figure 1D

cluster_df$cluster <- factor(
  cluster_df$cluster,
  levels = c("High", "Medium", "Low")
)

s_staked <- ggplot(cluster_df, aes(x = cluster, fill = MET_SITE)) +
  geom_bar(position = "fill") + 
  scale_fill_manual(values = c('Liver' = "darkorange", 
                               'Lung' = "green3",
                               'Brain' = 'maroon2' ,
                               'Bone' = 'mediumpurple3')) +
  scale_y_continuous(labels = function(x) x*100) +
  labs(fill = "Metastatic_site", y = 'Frequency (%)') +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7.5),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 6))

s_staked

# ggsave(
#   "prubas/Fig1_D.tiff", s_staked, device = "tiff",
#   width = 88, height = 70, units = "mm", dpi = 600, compression = "lzw"
# )

ggsave(
  "figures_def2026/Fig1_D.tiff", s_staked, device = "tiff",
  width = 88, height = 70, units = "mm", dpi = 600, compression = "lzw"
)
