# Script to analyze RNAseq data from BMEC STAR Rerun
# Uses DE table from AL to perform GO Analysis
# Uses DESeq2 object to generate clustered heatmap

# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.4.1

# 1) Setup -------------------------------------------------------------
library(data.table) # 1.16.4
library(tidyverse) # 2.0.0
library(ggplot2) # 3.5.1
library(patchwork) # 1.3.0
library(pheatmap) # 1.0.12
library(tidyr) # 1.3.1
library(enrichplot) # 1.26.6
library(UpSetR) # 1.4.0
library(DESeq2) # 1.46.0
library(clusterProfiler) # 4.14.4
library(org.Mm.eg.db) # 3.20.0

# Clear R's brain
rm(list = ls())

# 2) Data Cleanup and Formatting -------------------------------------

# Set working Directory
setwd("/Users/tylerlewy/Desktop/BMEC Star")

# Pull Data from DESeq Table
data <- as_tibble(fread("Objects/gene_clusters_top500_genes.tsv"))

# Cluster 1 Genes
# Pull genes names 
gene.ids.1 <- data%>%
  filter(cluster == 1)%>%
  pull(gene)

# Perform enrichment using mouse genome and GO:BP ontology
ego.c1 <- enrichGO(gene = gene.ids.1,
                    OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL",
                    ont = "BP",  
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)  

# Need to add a measurement for similarity between nodes for emapplot to work
to.emap.c1 <- pairwise_termsim(ego.c1)

emap.plot.c1 <- emapplot(to.emap.c1, 
                          showCategory = 15, 
                          group = T, 
                          node_label = "category",
                          layout = "fr", 
                          group_style = "ggforce",
                          label_group_style = "ggforce", 
                          nCluster = 3,
                          nWords = 5, 
                          min_edge = 0.4, 
                          label_format = 10)

enrich.export.c1 <- to.emap.c1@result%>%
  filter(p.adjust < 0.05)%>%
  arrange(desc(zScore))%>%
  dplyr::select(c("ID", "Description", "FoldEnrichment", "zScore", "p.adjust", "geneID"))

write_csv(enrich.export.c1,
          "/Users/tylerlewy/Desktop/BMEC Star/Objects/enrichment_c1.csv")

ggsave(emap.plot.c1,
       file=paste0("/Users/tylerlewy/Desktop/BMEC Star/Plots/emap.c1.pdf"),
       dpi = 300, units = c("cm"),width = 30, height = 30)

# Cluster 2 Genes
# Pull genes names 
gene.ids.2 <- data%>%
  filter(cluster == 2)%>%
  pull(gene)

# Perform enrichment using mouse genome and GO:BP ontology
ego.c2 <- enrichGO(gene = gene.ids.2,
                   OrgDb = org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",  
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)  

# Need to add a measurement for similarity between nodes for emapplot to work
to.emap.c2 <- pairwise_termsim(ego.c2)

emap.plot.c2 <- emapplot(to.emap.c2, 
                         showCategory = 15, 
                         group = T, 
                         node_label = "category",
                         layout = "fr", 
                         group_style = "ggforce",
                         label_group_style = "ggforce", 
                         nCluster = 3,
                         nWords = 5, 
                         min_edge = 0.4, 
                         label_format = 10)

enrich.export.c2 <- to.emap.c2@result%>%
  filter(p.adjust < 0.05)%>%
  arrange(desc(zScore))%>%
  dplyr::select(c("ID", "Description", "FoldEnrichment", "zScore", "p.adjust", "geneID"))

write_csv(enrich.export.c2,
          "/Users/tylerlewy/Desktop/BMEC Star/Objects/enrichment_c2.csv")

ggsave(emap.plot.c2,
       file=paste0("/Users/tylerlewy/Desktop/BMEC Star/Plots/emap.c2.pdf"),
       dpi = 300, units = c("cm"),width = 30, height = 30)

# Cluster 3 Genes
# Pull genes names 
gene.ids.3 <- data%>%
  filter(cluster == 3)%>%
  pull(gene)

# Perform enrichment using mouse genome and GO:BP ontology
ego.c3 <- enrichGO(gene = gene.ids.3,
                   OrgDb = org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",  
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)  

# Need to add a measurement for similarity between nodes for emapplot to work
to.emap.c3 <- pairwise_termsim(ego.c3)

emap.plot.c3 <- emapplot(to.emap.c3, 
                         showCategory = 15, 
                         group = T, 
                         node_label = "category",
                         layout = "fr", 
                         group_style = "ggforce",
                         label_group_style = "ggforce", 
                         nCluster = 3,
                         nWords = 5, 
                         min_edge = 0.4, 
                         label_format = 10)

enrich.export.c3 <- to.emap.c3@result%>%
  filter(p.adjust < 0.05)%>%
  arrange(desc(zScore))%>%
  dplyr::select(c("ID", "Description", "FoldEnrichment", "zScore", "p.adjust", "geneID"))

write_csv(enrich.export.c3,
          "/Users/tylerlewy/Desktop/BMEC Star/Objects/enrichment_c3.csv")

ggsave(emap.plot.c3,
       file=paste0("/Users/tylerlewy/Desktop/BMEC Star/Plots/emap.c3.pdf"),
       dpi = 300, units = c("cm"),width = 30, height = 30)
