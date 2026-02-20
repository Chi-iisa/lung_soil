### 
###
### Figure 2
###
### 

library(ConsensusTME)
library(ggplot2)
library(gplots)
library(reshape)
library(dplyr)
library(pheatmap)
library(grid)
library(gridExtra)
library(dendextend)

load('data/datasets_lung_soil.RData')
source("functions.R")
source("variables.R")
source("variables_crc.R")


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
      # t-test
      x <- merged_df %>% 
        select(GEO_ID, MET_SITE, `g`) %>% 
        filter(MET_SITE == 'Lung')
      y <- merged_df %>% 
        select(GEO_ID, MET_SITE, `g`) %>% 
        filter(MET_SITE == 'Liver')
      
      x[,g] <- as.numeric(x[,g])
      y[,g] <- as.numeric(y[,g])
      
      x_data <- x[,g]
      y_data <- y[,g]
      
      if (nrow(x) == 0 || nrow(y) == 0) {
        stop("One of those is empty, it isn't possible to do the t-test")
      }
      if (!is.numeric(x_data) || !is.numeric(y_data)) {
        stop("Non numeric data is used for t-test calculation")
      }
      
      t_result <- t.test(x_data, y_data)
      pv <- t_result$p.value
      g <- ggplot(merged_df, aes(x = MET_SITE,
                                 y = .data[[gene]],
                                 color = MET_SITE)) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = "white", size = 0.8) +
        scale_color_manual(values = c(
          "Lung" = "green3",
          "Liver" = "darkorange")) +
        labs(y = if (y_axis) "Enrichment score" else NULL,
             x = gene, color = NULL) +
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
      # Wilcoxon test
      w_result <- wilcox.test(merged_df[,pos] ~ MET_SITE, data = merged_df)
      pv <- w_result$p.value
      
      g <- ggplot(merged_df, aes(x = MET_SITE,
                                 y = .data[[gene]],
                                 color = MET_SITE)) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = "white", size = 0.8) +
        
        scale_color_manual(values = c(
          "Lung" = "green3",
          "Liver" = "darkorange")) +
      labs(y = if (y_axis) "Enrichment score" else NULL,
             x = gene, color = NULL) +
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

# Deconvolution with ConsensusTME
ConsensusTME::cancerAll
consensus <- consensusTMEAnalysis(bulkExp=as.matrix(ex_crc), 
                                  cancerType = "COAD", 
                                  statMethod = "ssgsea",
                                  immuneScore=T)
rownames(consensus) <- gsub("_", " ", rownames(consensus))
rownames(consensus)

t_cells <- rownames(consensus)[which(grepl('T ', rownames(consensus)))]
ncol <- 4
t_plots <- list()
y_lim <- range(consensus[t_cells,], na.rm = TRUE)

for (i in seq_along(t_cells)) {
  y_axis <- ((i - 1) %% ncol == 0)
  t_plots[[i]] <- boxplot_axis(samples_crc, consensus, t_cells[i], 
                               y_axis = y_axis, y_lim = y_lim)
}

library(patchwork)
combined <- wrap_plots(t_plots, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ggsave(
#   "prubas/tcell_figure2.tiff",
#   combined,
#   device = "tiff",
#   width = 180,
#   height = 130,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw"
# )

ggsave(
  "figures_def2026/Fig2_A.tiff", combined, device = "tiff",
  width = 180, height = 130, units = "mm", dpi = 600, compression = "lzw"
)

## Figure 2b
nk <- boxplot_axis(samples_crc, consensus,'NK cells', y_axis = TRUE)
grid.newpage()
grid.draw(nk)
ggsave(
  "figures_def2026/Fig2_B.tiff", nk, device = "tiff",
  width = 88, height = 70, units = "mm", dpi = 600, compression = "lzw"
)


## Figure 2C
# Prepare data
ex_crc_sign <- ex_crc[rownames(ex_crc) %in% sig_genes_crc,]

ex.sign.st <- scale(t(ex_crc_sign), center=T, scale=T)
ex.sign.st <- na.omit(ex.sign.st)

# GENESETS: Calculate distance between genesets and cluster distance using Ward linkage
d <- dist(ex.sign.st, method = "euclidean") # distance matrix
hcl1 <- hclust(d, method="ward.D2") 
column_dend = as.dendrogram(hcl1)

d <- dist(t(ex.sign.st), method = "euclidean") # distance matrix
hcl2 <- hclust(d, method="ward.D2") 
row_dend = as.dendrogram(hcl2)

paletteLength <- 100

library(RColorBrewer)
myColor <- colorRampPalette(brewer.pal(11, "RdBu"))(paletteLength)
myColor <- rev(myColor)

myBreaks <- seq(-2, 2, length.out = paletteLength + 1)
myBreaks


annotation <- data.frame(Metastatic_site = clin_crc$MET_SITE)
row.names(annotation) <- clin_crc$GEO_ID
annotation_colors <- list(Metastatic_site = c(Liver = "darkorange", 
                                              Lung = "green3"))

p <- pheatmap(
  mat= t(ex.sign.st),                  
  color = myColor,
  breaks = myBreaks,
  cluster_cols = hcl1,
  cluster_rows = hcl2,     
  show_rownames = TRUE,                    
  show_colnames = FALSE,                   
  fontsize_row = 5.75,                    
  annotation_col = annotation,             
  annotation_colors = annotation_colors,   
  legend = TRUE,
  border_color= NA,
)

ggsave(
  "figures_def2026/Fig2_C.tiff", p, device = "tiff",
  width = 180, height = 130, units = "mm", dpi = 600, compression = "lzw", bg = 'white'
)


## Figure 2d - TIS enrichment
# genes selected
TIS.gset <- list('TIS' = c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"))

library(GSVA)
params <- gsvaParam(as.matrix(ex_crc), TIS.gset)
GSVA <- gsva(params, verbose = FALSE)
GSVA <- as.data.frame(t(GSVA))

# Save data
write.table(GSVA, file="results/crc_tis_signature.txt", sep="\t", col.names = T, row.names = T, quote = F)

# Plot
file.tis.crc <-'results/crc_tis_signature.txt'
tis.crc <- plot.tis.ips(clin_crc, file.tis.crc, 'CRC', 'TIS', y_axis = TRUE)
grid.draw(tis.crc)
ggsave(
  "figures_def2026/Fig2_D.tiff", tis.crc, device = "tiff",
  width = 88, height = 70, units = "mm", dpi = 600, compression = "lzw"
)

## Figure 2E - IPS
library(devtools)
library(MCPcounter); library(estimate)

# Run Immunophenoscore
ipsmap <- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=2)
    }
  }
  return(ips)
}

# expression data
gene_expression <- as.data.frame(ex_crc)
sample_names <- colnames(gene_expression)

# Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
IPSG <- read.table("data/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
unique_ips_genes <-as.vector(unique(IPSG$NAME))

# add gene values for genes not in the expression matrix
mean_values <- apply(gene_expression, 2, mean)
genes <- IPSG$GENE[!(IPSG$GENE %in% rownames(gene_expression))]

for (gene in genes) {
  if (!(gene %in% rownames(gene_expression))) {
    new_row <- as.data.frame(t(mean_values))
    rownames(new_row) <- gene
    gene_expression <- rbind(gene_expression, new_row)
  }
}
gene_expression <- gene_expression[order(rownames(gene_expression)),]

IPS<-NULL
MHC<-NULL
CP<-NULL
EC<-NULL
SC<-NULL
AZ<-NULL

for (i in 1:length(sample_names)) {
  GE<-gene_expression[[i]]
  mGE<-mean(GE)
  sGE<-sd(GE)
  Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
  W1<-IPSG$WEIGHT
  WEIGHT<-NULL
  MIG<-NULL
  k<-1
  for (gen in unique_ips_genes) {
    MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)])
    WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k<-k+1
  }
  WG<-MIG*WEIGHT
  MHC[i]<-mean(WG[1:10])
  CP[i]<-mean(WG[11:20])
  EC[i]<-mean(WG[21:24])
  SC[i]<-mean(WG[25:26])
  AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
  IPS[i]<-ipsmap(AZ[i])
}
DF<-data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
rownames(DF) <- DF$SAMPLE
DF$SAMPLE <- NULL
ips <- t(DF)

# Save data
write.table(ips,file="results/IPS_scores_crc.txt",
            row.names=TRUE, col.names=TRUE, quote=FALSE,sep="\t")
# Plot
file.ips.crc <-"results/IPS_scores_crc.txt"
ips.crc <- plot.tis.ips(clin_crc, file.ips.crc, 'CRC', 'IPS', y_axis = TRUE)
grid.draw(ips.crc)
ggsave(
  "figures_def2026/Fig2_E.tiff", ips.crc, device = "tiff",
  width = 88, height = 70, units = "mm", dpi = 600, compression = "lzw"
)
