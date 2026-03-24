### 
###
### Supplementary Figure 1
###
### 


load('data/datasets_lung_soil.RData')
source("functions.R")
source("variables.R")
source("variables_all.R")

NK_all <- NK_all %>%
  group_by(final_symb) %>%
  slice(1) %>% 
  arrange(final_symb)

sig_df_all.gs <- sig_df_all %>%
  arrange(Gene) %>% 
  mutate(Geneset = ifelse(Gene == NK_all$final_symb, 
                          NK_all$geneset, NA))

geneset_sig.all <- sig_df_all.gs %>% 
  filter(is_significant ==TRUE) %>% 
  select(Geneset)
geneset_sig.all <- geneset_sig.all$Geneset

boxplot_gene_axis <- function(clin_data, exp_data, gene, 
                              y_axis = TRUE, y_lim = NULL) {
  
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
      
      g <- ggplot(data = merged_df, aes(x = MET_SITE, 
                                        y = merged_df[, gene],
                                        color = MET_SITE)) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = 'white', size = 0.8) +
        scale_color_manual(values = c("Lung" = "green3", 
                                      'Brain' = 'maroon2', 
                                      'Bone' = 'mediumpurple3', 
                                      'Liver' = 'darkorange')) +
        annotate("text", x = Inf, y = Inf, 
                 label = paste('p = ', ifelse(pv < 1e-5, format(pv, scientific = TRUE, digits = 4),
                                              round(pv, digits = 5))),
                 hjust = 1.1, vjust = 1.3, size = 3) +
        ggtitle(gene) + 
        labs(y = if (y_axis) "Log2 expression" else NULL) +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(size = 8.5, hjust = 0.5), 
              axis.title.x = element_blank(),
              axis.title.y = if (y_axis) element_text(size = 7.5) else element_blank(),
              axis.text.x = element_text(size = 7.5),
              axis.line = element_line(size = 0.3),
              legend.position = "none",  
              panel.background = element_blank(),  
              plot.background = element_blank(),  
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.y = if (y_axis) element_text(size = 7.5) else element_blank(),
              axis.line.y = if (y_axis) element_line(size = 0.3) else element_blank())
      
      # Adjust limits
      if (!is.null(y_lim)) {
        g <- g + coord_cartesian(ylim = y_lim)
      }
      print(paste0(gene, ' pv: ', pv))
      return(ggplotGrob(g + guides(fill = 'none'))) # Remove legend
    } else {
      # Kruskal Wallis
      kr_result <- kruskal.test(merged_df[, pos] ~ MET_SITE, data = merged_df)
      pv <- kr_result$p.value
      g <- ggplot(data = merged_df, aes(x = MET_SITE, 
                                        y = merged_df[, gene],
                                        color = MET_SITE)) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = 'white', size = 0.8) +
        scale_color_manual(values = c("Lung" = "green3", 
                                      'Brain' = 'maroon2', 
                                      'Bone' = 'mediumpurple3', 
                                      'Liver' = 'darkorange')) +
        annotate("text", x = Inf, y = Inf, 
                 label = paste('p = ', ifelse(pv < 1e-5, format(pv, scientific = TRUE, digits = 4),
                                              round(pv, digits = 5))),
                 hjust = 1.1, vjust = 1.3, size = 3) +
        ggtitle(gene) + 
        labs(y = if (y_axis) "Log2 expression" else NULL) +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(size = 8.5, hjust = 0.5), 
              axis.title.x = element_blank(),
              axis.title.y = if (y_axis) element_text(size = 7.5) else element_blank(),
              axis.text.x = element_text(size = 7.5),
              axis.line = element_line(size = 0.3),
              legend.position = "none",  
              panel.background = element_blank(),  
              plot.background = element_blank(),  
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.y = if (y_axis) element_text(size = 7.5) else element_blank(),
              axis.line.y = if (y_axis) element_line(size = 0.3) else element_blank())
      
      # Adjust limits
      if (!is.null(y_lim)) {
        g <- g + coord_cartesian(ylim = y_lim)
      }
      print(paste0(gene, ' pv: ', pv))
      return(ggplotGrob(g + guides(fill = 'none'))) # Remove legend
    }
  }
}

gset.plots <- lapply(unique(geneset_sig.all),
                     function(g) gs.sign.boxplot.all(samples_all, ex_all, sig_genes_all,
                                                     geneset_sig.all,g))



library(patchwork)
library(cowplot)
library(ggplotify)

row1 <- plot_grid(
  gset.plots[[2]], gset.plots[[4]], gset.plots[[5]],gset.plots[[6]],
  ncol = 4, nrow = 1,
  rel_widths = c(2,1,0.6,0.6)
)
row1

row2 <- plot_grid(
  gset.plots[[3]], gset.plots[[9]], gset.plots[[10]],
  ncol = 3, nrow = 1,
  rel_widths = c(3.5,0.65,0.65)
)
row2

row3b <- plot_grid(
  gset.plots[[8]], gset.plots[[7]],
  ncol = 2, nrow = 1,
  rel_widths = c(1,2.4)
)
row3b

row3_slamf <- plot_grid(
  gset.plots[[11]],
  ncol = 1, nrow = 1
)
row3_slamf

row3a <- plot_grid(
  gset.plots[[1]], row3b,
  ncol = 1, nrow = 2,
  rel_widths = c(1,1)
)
row3a

row3 <- plot_grid(
  row3a, row3_slamf,
  ncol = 2, nrow = 1,
  rel_widths = c(1,0.15)
)
row3

suppl.fig1 <- plot_grid(
  row1,row2,row3,
  nrow = 3, ncol = 1,
  rel_heights = c(1,1,2)
)
suppl.fig1

suppl.fig1_backg <- suppl.fig1 +
  theme(plot.background = element_rect(fill = "white", color = NA))

save_plot(
  "figures_def2026/Suppl-Fig1_v8.tiff", 
  suppl.fig1_backg,
  base_width = 385 / 25.4,
  base_height = 220 / 25.4, 
  dpi = 600       
)
