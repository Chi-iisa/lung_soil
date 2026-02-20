library(tidyverse)
library(ConsensusTME)
library(gplots)
library(reshape)
library(pheatmap)
library(grid)
library(gridExtra)
library(dendextend)
library(patchwork)

check_genesNK <- function(i, j, genes_df, ex_mt, df_res){
  gene_symb <- genes_df[i,j]
  genes_ex_mt <- rownames(ex_mt)
  if(i <= nrow(genes_df) & j >= 4 & j <= 17){
    if(is.na(gene_symb)){
      df_res$final_symb[i] <- df_res$Hugo_Symb[i]
      df_res$is_synonym_used[i] <- NA
      # Restart j value for consulting the next gene
      j <- 4 
      return(check_genesNK(i+1, j, genes_df, ex_mt, df_res))
    } else {
      if(any(genes_ex_mt == gene_symb)){
        df_res$final_symb[i] <- gene_symb
        hugo <- df_res$Hugo_Symb[i]
        if(gene_symb == hugo){
          df_res$is_synonym_used[i] <- FALSE
        } else{
          df_res$is_synonym_used[i] <- TRUE
        }
        return(check_genesNK(i+1, j, genes_df, ex_mt, df_res))
      } else {
        if (j == 17){
          df_res$final_symb[i] <- df_res$Hugo_Symb[i]
          df_res$is_synonym_used[i] <- NA 
          j <- 4 # Restart j value for consulting the next gene
          return(check_genesNK(i+1, j, genes_df, ex_mt, df_res))
        } else {
          return(check_genesNK(i, j+1, genes_df, ex_mt, df_res))
        }
      }
    }
  }
  
  return(df_res)
}



pvalue_gene_lung_liver <- function(clin_data, exp_data, gene, dataset){
  
  t <- as.data.frame(t(exp_data))
  
  if (dataset %in% c('All', 'Melanoma', 'CRC')){
    t$GEO_ID <- rownames(t)
    
    merged_df <- merge(clin_data, t, by = "GEO_ID")
    # Consider only the distant metastasis sites
    distant <- c('Liver', 'Lung')
    merged_df <- merged_df[merged_df$MET_SITE %in% distant,]
    
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
    
    ex_filtered <- as.data.frame(exp_data) %>% 
      select(which(colnames(as.data.frame(exp_data)) %in% merged_df$GEO_ID))
  } else if(dataset == 'Pancreas'){
    t$SAMPLE_ID <- rownames(t)
    
    merged_df <- merge(clin_data, t, by = "SAMPLE_ID")
    # Consider only the distant metastasis sites
    distant <- c('Liver', 'Lung')
    merged_df <- merged_df[merged_df$TUMOR_SITE %in% distant,]
    
    merged_df$TUMOR_SITE <- factor(merged_df$TUMOR_SITE) # factor type
    
    count_df <- merged_df %>% 
      group_by(TUMOR_SITE) %>% 
      summarise(count = n())
    
    # Keep these samples with N > 1
    filtered_sites <- count_df$count > 1
    filtered_sites <- count_df$TUMOR_SITE[filtered_sites]
    
    # Filter original dataframe
    merged_df <- merged_df %>% 
      filter(TUMOR_SITE %in% filtered_sites)
    
    ex_filtered <- as.data.frame(exp_data) %>% 
      select(which(colnames(as.data.frame(exp_data)) %in% merged_df$SAMPLE_ID))
  }
  
  if(any(rownames(ex_filtered)==gene)){
    pos <- which(colnames(merged_df)==gene)
    data_gene <- as.numeric(c(ex_filtered[gene,]))
    # Verify if it follows a normal distribution
    km <- ks.test(data_gene,"pnorm", mean=mean(data_gene), sd=sd(data_gene))
    if(km$p.value >= 0.05){
      # t-test test
      if(dataset %in% c('All', 'Melanoma', 'CRC')){
        x <- merged_df %>% 
          select(GEO_ID, MET_SITE, `gene`) %>% 
          filter(MET_SITE == 'Lung')
        y <- merged_df %>% 
          select(GEO_ID, MET_SITE, `gene`) %>% 
          filter(MET_SITE == 'Liver')
      } else if(dataset == 'Pancreas'){
        x <- merged_df %>% 
          select(SAMPLE_ID, TUMOR_SITE, `gene`) %>% 
          filter(TUMOR_SITE == 'Lung')
        y <- merged_df %>% 
          select(SAMPLE_ID, TUMOR_SITE, `gene`) %>% 
          filter(TUMOR_SITE == 'Liver')
      }
      
      x[,gene] <- as.numeric(x[,gene])
      y[,gene] <- as.numeric(y[,gene])
      
      x_data <- x[,gene]
      y_data <- y[,gene]
      
      if (nrow(x) == 0 || nrow(y) == 0) {
        stop("One of those is empty, it isn't possible to do the t-test")
      }
      if (!is.numeric(x_data) || !is.numeric(y_data)) {
        stop("Non numeric data is used for t-test calculation")
      }
      
      t_result <- t.test(x_data, y_data)
      pv <- t_result$p.value
      return(pv)
    } else{
      # Wilcoxon test
      if(dataset %in% c('All', 'Melanoma', 'CRC')){
        w_result <- wilcox.test(merged_df[,pos] ~ MET_SITE, data = merged_df)
      } else if(dataset == 'Pancreas'){
        w_result <- wilcox.test(merged_df[,pos] ~ TUMOR_SITE, data = merged_df)
      }
      
      pv <- w_result$p.value
      return(pv)
    }
  } else{
    res <- paste0('The gene ', gene, ' is not in the expression matrix')
    print(res)
    return(NA)
  }
}


sign_genes <- function(clin_data, exp_data, df, dataset){
  for (i in seq_along(df$Gene)){
    g <- df$Gene[i]
    res <- pvalue_gene_lung_liver(clin_data, exp_data, g, dataset)
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
  return(list(sig_genes,df))
}


## > 2 groups comparison pvalue / significative genes

pvalue_gene_all <- function(clin_data, exp_data, gene){
  
  t <- as.data.frame(t(exp_data))
  t$GEO_ID <- rownames(t)
  
  merged_df <- merge(clin_data, t, by = "GEO_ID")
  # Consider only the distant metastasis sites
  distant <- c('Bone','Brain', 'Liver', 'Lung')
  merged_df <- merged_df[merged_df$MET_SITE %in% distant,]
  
  merged_df <- merged_df %>% 
    filter(TISSUE == "metastasis")
  
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
    # Verify if it follows a normal distribution
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
    res <- paste0('The gene ', gene, ' is not in the expression matrix')
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


gs.sign.boxplot <- function(clin_data, exp_data, sign.gene.vector, sign.gs.vector, geneset){
  pos.consulta <- which(sign.gs.vector == geneset)
  genes.gs <- sign.gene.vector[pos.consulta]
  
  res <- lapply(genes.gs, function(g) boxplot_gene_livlun(clin_data, exp_data, g))
  ncol.grid <- length(res)
  res.patch <- wrap_plots(res, ncol = 6)
  
  p <- res.patch + plot_annotation(
    title = geneset,
    theme = theme(
      plot.title = element_text(face = "bold", size = 16
                                # , hjust = 0.5
      )
    )
  )
  return(p)
}


median_zs_data <- function(ex_filtered_ini, samples_keep, dataset, 
                           genes_sign, geneset, split_df){
  if(nrow(split_df[[`geneset`]]) != 0) {
    
    data_geneset <- as.data.frame(split_df[[`geneset`]])
    
    filtered_geneset <- data_geneset %>% 
      filter(!is.na(is_synonym_used)) %>% 
      select(final_symb)
    
    filtered_geneset_vector <- filtered_geneset$final_symb
    rownames(ex_filtered_ini) <- gsub("\\.", "-", rownames(ex_filtered_ini))
    ex_filtered <- ex_filtered_ini[filtered_geneset_vector,]
    
    expression_data_transposed <- as.data.frame(t(ex_filtered))
    
    if(dataset %in% c('All','CRC','Melanoma')){
      expression_data_transposed$GEO_ID <- rownames(expression_data_transposed)
      
      # Combine expression data with clinical data (met_sites)
      expression_data_with_sites <- merge(expression_data_transposed, samples_keep, 
                                          by = 'GEO_ID')
      
      expression_data_with_sites <- expression_data_with_sites[, -which(names(expression_data_with_sites) == 'GEO_ID')]
      
      # Median calculation by met_site group
      median_expression <- expression_data_with_sites %>%
        group_by(MET_SITE) %>%
        summarise(across(everything(), median, na.rm = TRUE))
      
    } else if(dataset == 'Pancreas'){
      expression_data_transposed$SAMPLE_ID <- rownames(expression_data_transposed)
      
      # Combine expression data with clinical data (met_sites)
      expression_data_with_sites <- merge(expression_data_transposed, samples_keep, 
                                          by = 'SAMPLE_ID')
      
      expression_data_with_sites <- expression_data_with_sites[, -which(names(expression_data_with_sites) == 'SAMPLE_ID')]
      
      # Median calculation by tumor_site group
      median_expression <- expression_data_with_sites %>%
        group_by(TUMOR_SITE) %>%
        summarise(across(everything(), median, na.rm = TRUE))
    }
    
    # Scale median values into z-score
    # Remove met_site info column
    median_expression_scaled <- as.data.frame(scale(median_expression[, -1]))
    return(median_expression_scaled)
  } else{
    print(paste0('There is not any gene from ', geneset, ' geneset'))
  }
}


## Supplementary Fig 5
prepare_data <- function(data1, data2, data3, data4, geneset){
  d1 <- data1[[`geneset`]]
  d2 <- data2[[`geneset`]]
  d3 <- data3[[`geneset`]]
  d4 <- data4[[`geneset`]]
  
  my_list <- c(d1,d2,d3,d4)
  pos_data <- !sapply(my_list, is.character)
  available_data <- my_list[pos_data]
  res_data <- c(unlist(available_data))
  
  return (res_data)
}


define_breaks <- function(data, geneset){
  my_data <- data[[`geneset`]]
  
  breaks_range <- seq(min(my_data, na.rm = TRUE), max(my_data, na.rm = TRUE),
                      length.out = 101)
  return(breaks_range)
}


heatmap_zscore_sign <- function(ex_filtered_ini, samples_keep, 
                                genes_sign, geneset, split_df, dataset,
                                color_palette, break_range, show_legend){
  print(geneset)
  if(nrow(split_df[[`geneset`]]) != 0) {
    
    data_geneset <- as.data.frame(split_df[[`geneset`]])
    
    filtered_geneset <- data_geneset %>% 
      filter(!is.na(is_synonym_used)) %>% 
      select(final_symb)
    
    filtered_geneset_vector <- filtered_geneset$final_symb
    rownames(ex_filtered_ini) <- gsub("\\.", "-", rownames(ex_filtered_ini))
    ex_filtered <- ex_filtered_ini[filtered_geneset_vector,]
    
    expression_data_transposed <- as.data.frame(t(ex_filtered))
    
    # Combine expression data with clinical data (met_sites)
    if(dataset %in% c('All','Melanoma','CRC')){
      expression_data_transposed$GEO_ID <- rownames(expression_data_transposed)
      expression_data_with_sites <- merge(expression_data_transposed, 
                                          samples_keep, by = 'GEO_ID')
      expression_data_with_sites <- expression_data_with_sites[, -which(names(expression_data_with_sites) == 'GEO_ID')]
      
      # Median calculation by met_site group
      median_expression <- expression_data_with_sites %>%
        group_by(MET_SITE) %>%
        summarise(across(everything(), median, na.rm = TRUE))
      
      # Scale median values into z-score
      # Remove met_site info column
      median_expression_scaled <- as.data.frame(scale(median_expression[, -1]))  
      median_expression$MET_SITE <- str_to_title(median_expression$MET_SITE)
      
      # Define met_sites as rownames
      rownames(median_expression_scaled) <- median_expression$MET_SITE
      
    } else if(dataset == 'Pancreas'){
      expression_data_transposed$SAMPLE_ID <- rownames(expression_data_transposed)
      expression_data_with_sites <- merge(expression_data_transposed, 
                                          samples_keep, by = 'SAMPLE_ID')
      expression_data_with_sites <- expression_data_with_sites[, -which(names(expression_data_with_sites) == 'SAMPLE_ID')]
      
      # Median calculation by TUMOR_SITE group
      median_expression <- expression_data_with_sites %>%
        group_by(TUMOR_SITE) %>%
        summarise(across(everything(), median, na.rm = TRUE))
      
      
      # Scale median values into z-score
      # Remove tumor_site info column
      median_expression_scaled <- as.data.frame(scale(median_expression[, -1]))  
      median_expression$TUMOR_SITE <- str_to_title(median_expression$TUMOR_SITE)
      
      # Define tumor_site as row names
      rownames(median_expression_scaled) <- median_expression$TUMOR_SITE
    }
    
    median_expression_scaled <- t(median_expression_scaled)
    rownames(median_expression_scaled) <- gsub("\\.", "-", rownames(median_expression_scaled))
    
    # Define asterisk gene data
    asterisks_labels <- matrix("", nrow = nrow(median_expression_scaled), 
                               ncol = ncol(median_expression_scaled),
                               dimnames = list(rownames(median_expression_scaled), 
                                               colnames(median_expression_scaled)))
    # For each significant gene
    for (gene in genes_sign) {
      if(gene %in% rownames(median_expression_scaled)){
        # Get z-scores
        gene_zscores <- median_expression_scaled[gene, ]
        # Find the maximum one
        max_zscore_index <- which.max(gene_zscores)
        asterisks_labels[gene, max_zscore_index] <- "*"
      } 
    }
    
    # Prepare heatmap input matrix only with the genes of interest
    if(nrow(median_expression_scaled)==1){
      mt_filtered <- median_expression_scaled
    } else{
      mt_filtered <- median_expression_scaled[filtered_geneset_vector,]
    }

    
    if(any(duplicated(rownames(mt_filtered)))){
      rownames(asterisks_labels) <- gsub("\\-", ".", rownames(asterisks_labels))
      mt_filtered <- median_expression_scaled[filtered_geneset_vector,]
      # Remove any duplicate
      mt_filtered <- mt_filtered[!duplicated(rownames(mt_filtered)), ]
      # asterisks_labels must to have the same dimensions
      asterisks_labels <- asterisks_labels[rownames(mt_filtered),]
    }
    
    n_genes <- nrow(mt_filtered)
    
    # Generate heatmap
    pheatmap(mt_filtered,
             display_numbers = asterisks_labels,
             fontsize_number = 12,
             cluster_rows = FALSE,  
             cluster_cols = FALSE,
             angle_col = 0,
             main = dataset,
             fontsize = 6,
             fontsize_row = 8,
             fontsize_col = 8,
             color = color_palette,
             breaks = break_range[[`geneset`]],
             silent = TRUE,
             legend = show_legend,
             show_rownames = show_legend,
             na_col = 'white',
             cellheight = 20,
             cellwidth = 25
    )
  } else{
    print(paste0('There is not any gene from ', geneset, ' geneset'))
  }
}

get_grob_size_mm <- function(gt) {
  width_mm  <- convertWidth(sum(gt$widths), "mm", valueOnly = TRUE)
  height_mm <- convertHeight(sum(gt$heights), "mm", valueOnly = TRUE)
  
  list(width = width_mm, height = height_mm)
}


supp.figure5 <- function(heat1, heat2, heat3, heat4, g_names){
  layout_list <- list()
  
  max_cols <- 5
  for (i in seq_along(g_names)){
    group_gene <- g_names[i]
    h1 <- heat1[[`group_gene`]]
    h2 <- heat2[[`group_gene`]]
    h3 <- heat3[[`group_gene`]]
    h4 <- heat4[[`group_gene`]]
    my_list <- list(h1, h2, h3, h4)
    
    group_text <- ifelse(grepl("Ligands", group_gene),
                         paste0(toupper(substr(group_gene,1,1)), 
                                tolower(substr(group_gene,2,nchar(group_gene)))),
                         group_gene)
    group_label <- textGrob(group_text, rot = 90, 
                            x = 0.8, 
                            gp = gpar(fontsize = 10, fontface = "bold"))
    
    # Keep only these datasets with available geneset
    pos_data <- !sapply(my_list, is.character)
    available_data <- my_list[pos_data]
    gtable_list <- lapply(available_data, function(x) x$gtable)
    
    # Complete the layout with white spaces
    if (length(gtable_list) < (max_cols - 1)) {
      gtable_list <- c(gtable_list, replicate((max_cols - 1) - length(gtable_list), 
                                              nullGrob(), simplify = FALSE))
    }
    
    # Generate a new row with title and heatmaps
    row_layout <- arrangeGrob(
      grobs = c(list(group_label), gtable_list),
      widths = c(0.02, rep(0.24, length(gtable_list))),
      ncol = max_cols
    )
    
    # h <- onvertHeight(sum(row_layout$heights), "mm", TRUE)
    
    sizes <- lapply(gtable_list, get_grob_size_mm)
    heights <- sapply(sizes, `[[`, "height")
    
    h <- max(heights) + 20
    
    ggsave(paste0("figures_def2026/Suppl-Fig5-",i,'_',gsub(" ", "_", group_gene),".tiff"),
           row_layout, device = "tiff", width = 180, height = h, units = "mm", 
           dpi = 600, compression = "lzw", bg = 'white')
    
    print(paste0('Saved: ', group_gene))
  }
}


gs.sign.boxplot <- function(clin_data, exp_data, sign.gene.vector, sign.gs.vector, geneset){
  pos.consulta <- which(sign.gs.vector == geneset)
  genes.gs <- sign.gene.vector[pos.consulta]
  
  res <- lapply(genes.gs, function(g) boxplot_gene_livlun(clin_data, exp_data, g))
  ncol.grid <- length(res)
  res.patch <- wrap_plots(res, ncol = case_when(
    (geneset %in% c('Activating Ligands', 'HLA-II')) ~ 4,
    geneset =='Activating NK receptors' ~ 5,
    TRUE ~ 6)
  )
  
  group_text <- ifelse(grepl("Ligands", geneset),
                       paste0(toupper(substr(geneset,1,1)), 
                              tolower(substr(geneset,2,nchar(geneset)))),
                       geneset)
  p <- res.patch + plot_annotation(
    title = group_text,
    theme = theme(
      plot.title = element_text(face = "bold", size = 10))
  )
  return(p)
}

gs.sign.boxplot.all <- function(clin_data, exp_data, sign.gene.vector, sign.gs.vector, geneset){
  pos.consulta <- which(sign.gs.vector == geneset)
  genes.gs <- sign.gene.vector[pos.consulta]
  res <- lapply(genes.gs, function(g) boxplot_gene_axis(clin_data, exp_data, g))
  ncol.grid <- length(res)
  res.patch <- wrap_plots(res, ncol = ifelse(geneset == 'Activating NK receptors', 3, 4))
  
  group_text <- ifelse(grepl("Ligands", geneset),
                       paste0(toupper(substr(geneset,1,1)), 
                              tolower(substr(geneset,2,nchar(geneset)))),
                       geneset)
  
  p <- res.patch + plot_annotation(
    title = group_text,
    theme = theme(
      plot.title = element_text(face = "bold", size = 10)
    )
  )
  return(p)
}


plot.tis.ips <- function(clin_data, file.data, dataset, score, y_axis = FALSE, y_lim = NULL){
  
  if (score =='TIS'){
    title <- 'TIS (T-Cell inflammatory signature)'
    t <- read.table(file.data)
    data_score <- as.numeric(t[,1])
  } else if (score =='IPS') {
    title <- 'Immunophenoscore'
    t <- as.data.frame(t(read.table(file.data)))
    data_score <- t$IPS
  }
  
  if(dataset %in% c('Melanoma', 'CRC', 'All')){
    t$GEO_ID <- rownames(t)
    merged_df <- merge(clin_data, t, by = "GEO_ID")
    merged_df$MET_SITE <- factor(merged_df$MET_SITE) # factor type
    
    count_df <- merged_df %>% 
      group_by(MET_SITE) %>% 
      summarise(count = n())
    
    # Keep samples with N > 1
    filtered_sites <- count_df$count > 1
    filtered_sites <- count_df$MET_SITE[filtered_sites]
    
    # Filter original dataframe
    merged_df <- merged_df %>% 
      filter(MET_SITE %in% filtered_sites)
    
    
    # Verify if it follows a normal distribution 
    print(paste("Normality test for type:", score))
    km <- ks.test(data_score, "pnorm", mean=mean(data_score), sd=sd(data_score))
    print(km)
    
    if(km$p.value >= 0.05 & mean(data_score) != 0){
      # t-test
      x <- merged_df %>% 
        select(GEO_ID, MET_SITE, `score`) %>% 
        filter(MET_SITE == 'Lung')
      y <- merged_df %>% 
        select(GEO_ID, MET_SITE, `score`) %>% 
        filter(MET_SITE == 'Liver')
      
      x[,score] <- as.numeric(x[,score])
      y[,score] <- as.numeric(y[,score])
      
      x_data <- x[,score]
      y_data <- y[,score]
      
      if (nrow(x) == 0 || nrow(y) == 0) {
        stop("One of those is empty, it isn't possible to do the t-test")
      }
      if (!is.numeric(x_data) || !is.numeric(y_data)) {
        stop("Non numeric data is used for t-test calculation")
      }
      
      t_result <- t.test(x_data, y_data)
      pv <- t_result$p.value
      g <- ggplot(merged_df, aes(x = MET_SITE,
                                 y = .data[[score]],
                                 color = MET_SITE)) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = "white", size = 0.8) +
        
        scale_color_manual(values = c(
          "Lung" = "green3",
          "Liver" = "darkorange")) +
        
        labs(y = if (y_axis) "Enrichment score" else NULL,
             x = title,
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
    } else if(km$p.value < 0.05 & mean(data_score) != 0){
      # Wilcoxon test
      w_result <- wilcox.test(merged_df[,score] ~ MET_SITE, data = merged_df)
      pv <- w_result$p.value
      g <- ggplot(merged_df, aes(x = MET_SITE,
                                 y = .data[[score]],
                                 color = MET_SITE)) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = "white", size = 0.8) +
        
        scale_color_manual(values = c(
          "Lung" = "green3",
          "Liver" = "darkorange")) +
        
        labs(y = if (y_axis) "Enrichment score" else NULL,
             x = title,
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
  } else if (dataset =='Pancreas'){
    rownames(t) <- gsub("\\.", "-", rownames(t))
    t$SAMPLE_ID <- rownames(t)
    merged_df <- merge(clin_data, t, by = "SAMPLE_ID")
    merged_df$TUMOR_SITE <- factor(merged_df$TUMOR_SITE) # factor type
    
    count_df <- merged_df %>% 
      group_by(TUMOR_SITE) %>% 
      summarise(count = n())
    
    # Keep samples with N > 1
    filtered_sites <- count_df$count > 1
    filtered_sites <- count_df$TUMOR_SITE[filtered_sites]
    
    # Filter original dataframe
    merged_df <- merged_df %>% 
      filter(TUMOR_SITE %in% filtered_sites)
    
    data_score <- as.numeric(t[,1])
    # Verify if it follows a normal distribution 
    print(paste("Normality test for type:", score))
    km <- ks.test(data_score, "pnorm", mean=mean(data_score), sd=sd(data_score))
    print(km)
    
    if(km$p.value >= 0.05 & mean(data_score) != 0){
      # t-test
      x <- merged_df %>% 
        select(SAMPLE_ID, TUMOR_SITE, `score`) %>% 
        filter(TUMOR_SITE == 'Lung')
      y <- merged_df %>% 
        select(SAMPLE_ID, TUMOR_SITE, `score`) %>% 
        filter(TUMOR_SITE == 'Liver')
      
      x[,score] <- as.numeric(x[,score])
      y[,score] <- as.numeric(y[,score])
      
      x_data <- x[,score]
      y_data <- y[,score]
      
      if (nrow(x) == 0 || nrow(y) == 0) {
        stop("One of those is empty, it isn't possible to do the t-test")
      }
      if (!is.numeric(x_data) || !is.numeric(y_data)) {
        stop("Non numeric data is used for t-test calculation")
      }
      
      t_result <- t.test(x_data, y_data)
      pv <- t_result$p.value
      g <- ggplot(merged_df, aes(x = TUMOR_SITE,
                                 y = .data[[score]],
                                 color = TUMOR_SITE)) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = "white", size = 0.8) +
        
        scale_color_manual(values = c(
          "Lung" = "green3",
          "Liver" = "darkorange")) +
        
        labs(y = if (y_axis) "Enrichment score" else NULL,
             x = title,
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
    } else if(km$p.value < 0.05 & mean(data_score) != 0){
      # Wilcoxon test
      w_result <- wilcox.test(merged_df[,score] ~ MET_SITE, data = merged_df)
      pv <- w_result$p.value
      g <- ggplot(merged_df, aes(x = TUMOR_SITE,
                                 y = .data[[score]],
                                 color = TUMOR_SITE)) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = "white", size = 0.8) +
        
        scale_color_manual(values = c(
          "Lung" = "green3",
          "Liver" = "darkorange")) +
        
        labs(y = if (y_axis) "Enrichment score" else NULL,
             x = title,
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
