### 
###
### Supplementary Figure 2
###
### 


load('data/datasets_lung_soil.RData')
source("functions.R")
source("variables.R")
source("variables_crc.R")

NK_crc <- NK_crc %>%
  group_by(final_symb) %>%
  slice(1) %>% 
  arrange(final_symb)


sig_df_crc.edit <- sig_df_crc %>%
  arrange(Gene) %>% 
  mutate(Geneset = ifelse(Gene == NK_crc$final_symb, 
                          NK_crc$geneset, NA))

geneset_sig.crc <- sig_df_crc.edit %>% 
  filter(is_significant ==TRUE) %>% 
  select(Geneset)
geneset_sig.crc <- geneset_sig.crc$Geneset


boxplot_gene_livlun <- function(clin_data, exp_data, gene, y_axis = TRUE, y_lim = NULL){
  
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
  
  if(any(rownames(exp_data)==gene)){
    pos <- which(colnames(merged_df)==gene)
    data_gene <- as.numeric(exp_data[gene, ])
    # Verify if it follows a normal distribution 
    print(paste("Normality test for gene:", gene))
    km <- ks.test(data_gene, "pnorm", mean=mean(data_gene), sd=sd(data_gene))
    print(km)
    
    if(km$p.value >= 0.05 & mean(data_gene) != 0){
      # t-test (bc there are only two independent groups)
      x <- merged_df %>% 
        select(GEO_ID, MET_SITE, `gene`) %>% 
        filter(MET_SITE == 'Lung')
      y <- merged_df %>% 
        select(GEO_ID, MET_SITE, `gene`) %>% 
        filter(MET_SITE == 'Liver')
      
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
      g <- ggplot(data = merged_df, aes(x = MET_SITE, 
                                        y = merged_df[,pos],
                                        color = MET_SITE)) +
        stat_boxplot(geom = "errorbar", # Boxplot with error bars 
                     width = 0.2) +
        geom_boxplot(fill = 'white',
                     size = 0.5) +
        scale_color_manual(values = c("Lung" = "green3",
                                      'Liver' = 'darkorange')) +
        annotate("text", x = Inf, y = Inf, 
                 label = paste('p = ', ifelse(pv < 1e-5, format(pv, scientific = TRUE, digits = 3),
                                              round(pv, digits = 5))),
                 hjust = 1.1, vjust = 1.3, size = 2) +
        ggtitle(gene) + 
        labs(y = if (y_axis) "Log2 expression" else NULL) +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(size = 8.5, hjust = 0.5), 
              axis.title.x = element_blank(),
              axis.title.y = if (y_axis) element_text(size = 6.5) else element_blank(),
              axis.text.x = element_text(size = 6.5),
              axis.line = element_line(size = 0.3),
              legend.position = "none",  
              panel.background = element_blank(),  
              plot.background = element_blank(),  
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.y = if (y_axis) element_text(size = 6.5) else element_blank(),
              # axis.ticks.y = if (y_axis) element_line() else element_blank(),
              axis.line.y = if (y_axis) element_line(size = 0.3) else element_blank())
      
      # Adjust limits
      if (!is.null(y_lim)) {
        g <- g + coord_cartesian(ylim = y_lim)
      }
      print(paste0(gene, ' | t-test pv: ', pv))
      return(ggplotGrob(g+guides(fill = 'none'))) # Remove legend
    } else if(km$p.value < 0.05 & mean(data_gene) != 0){
      # Wilcoxon test
      w_result <- wilcox.test(merged_df[,pos] ~ MET_SITE, data = merged_df)
      pv <- w_result$p.value
      g <- ggplot(data = merged_df, aes(x = MET_SITE, 
                                        y = merged_df[,pos],
                                        color = MET_SITE)) +
        stat_boxplot(geom = "errorbar",
                     width = 0.2) +
        geom_boxplot(fill = 'white',
                     size = 0.5) +
        scale_color_manual(values = c("Lung" = "green3",
                                      'Liver' = 'darkorange')) +
        annotate("text", x = Inf, y = Inf, 
                 label = paste('p = ', ifelse(pv < 1e-5, format(pv, scientific = TRUE, digits = 3),
                                              round(pv, digits = 5))),
                 hjust = 1.1, vjust = 1.3, size = 2) +
        ggtitle(gene) + 
        labs(y = if (y_axis) "Log2 expression" else NULL) +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(size = 8.5, hjust = 0.5), 
              axis.title.x = element_blank(),
              axis.title.y = if (y_axis) element_text(size = 6.5) else element_blank(),
              axis.text.x = element_text(size = 6.5),
              axis.line = element_line(size = 0.3),
              legend.position = "none",  
              panel.background = element_blank(),  
              plot.background = element_blank(),  
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.y = if (y_axis) element_text(size = 6.5) else element_blank(),
              # axis.ticks.y = if (y_axis) element_line() else element_blank(),
              axis.line.y = if (y_axis) element_line(size = 0.3) else element_blank())
      
      # Adjust limits
      if (!is.null(y_lim)) {
        g <- g + coord_cartesian(ylim = y_lim)
      }
      print(paste0(gene, ' | Wilcoxon test pv: ', pv))
      return(ggplotGrob(g+guides(fill = 'none'))) # Remove legend
    } else{
      res <- paste0('The gene ', gene, ' has null values in the expression matrix')
      print(res)
    }
  } else{
    res <- paste0('The gene ', gene, ' is not in the expression matrix')
    print(res)
  }
}

gset.plots <- lapply(unique(geneset_sig.crc),
                     function(g) gs.sign.boxplot(samples_crc, ex_crc,sig_genes_crc, 
                                                 geneset_sig.crc,g))



library(patchwork)
library(cowplot)
library(ggplotify)

row1 <- plot_grid(
  gset.plots[[2]], gset.plots[[7]], gset.plots[[13]],
  ncol = 3, nrow = 1,
  rel_widths = c(2,1,0.66)
)
row1

row2 <- plot_grid(
  gset.plots[[4]], gset.plots[[6]], gset.plots[[5]],
  ncol = 3, nrow = 1,
  rel_widths = c(1.8,0.6,1.1)
)
row2

row3 <- plot_grid(
  gset.plots[[3]], gset.plots[[8]],
  ncol = 2, nrow = 1,
  rel_widths = c(2,0.66)
)
row3

row4 <- plot_grid(
  gset.plots[[1]],gset.plots[[9]], 
  ncol = 2, nrow = 1,
  rel_widths = c(1,1)
)
row4

row5 <- plot_grid(
  gset.plots[[12]], gset.plots[[11]],gset.plots[[10]], 
  ncol = 3, nrow = 1,
  rel_widths = c(0.8,0.8,1.6)
)
row5

suppl.fig2 <- plot_grid(
  row1, row2, row3, row4, row5,
  ncol = 1, nrow = 5,
  rel_heights = c(1,1,1,1,1)
)
suppl.fig2

suppl.fig2_backg <- suppl.fig2 +
  theme(plot.background = element_rect(fill = "white", color = NA))

save_plot(
  "figures_def2026/Suppl-Fig2_v1.tiff", 
  suppl.fig2_backg,
  base_width = 385 / 25.4,
  base_height = 220 / 25.4, 
  dpi = 600       
)