# Script to analyze processed sequencing data for XR259
# Uses table of sig genes provided by AL to perform GO Analysis and emap plotting

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
setwd("/Users/tylerlewy/Desktop/XR259/")

# Pull Data from significance table
data <- as_tibble(fread("Objects/XR259_ISG_heatmap_individual.tsv"))

# 3) GO Enrichment -------------------------------------------------

# Uses table provided by AL

# Pull genes that were significantly upregulated following poly(I:C) and IFN
gene.ids.c1 <- data%>%
  filter(cluster == 1)%>%
  pull(gene_id)

# Perform enrichment using mouse genome and GO:BP ontology
ego.c1 <- enrichGO(gene = gene.ids.c1,
                    OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL",
                    ont = "BP",  
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)  

# Create a table of significant results ordered by z score
enrich.export.c1 <- ego.c1@result%>%
  filter(p.adjust < 0.05)%>%
  arrange(desc(zScore))%>%
  dplyr::select(c("ID", "Description", "FoldEnrichment", "zScore", "p.adjust"))

# Export results as csv
write_csv(enrich.export.c1,
          "/Users/tylerlewy/Desktop/XR259/Objects/go_enrichment_cluster1.csv")

# Repeat for cluster 2 (genes upregulated by poly(I:C) only)
gene.ids.c2 <- data%>%
  filter(cluster == 2)%>%
  pull(gene_id)

# Perform enrichment using mouse genome and GO:BP ontology
ego.c2 <- enrichGO(gene = gene.ids.c2,
                    OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL",
                    ont = "BP",  
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)  

# Create a table of significant results ordered by z score
enrich.export.c2 <- ego.c2@result%>%
  filter(p.adjust < 0.05)%>%
  arrange(desc(zScore))%>%
  dplyr::select(c("ID", "Description", "FoldEnrichment", "zScore", "p.adjust"))

# Export results as csv
write_csv(enrich.export.c2,
          "/Users/tylerlewy/Desktop/XR259/Objects/go_enrichment_cluster2.csv")

# 5) Plotting - Unique Enrichment Plots ---------------------------------------

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

ggsave(emap.plot.c1,
       file=paste0("/Users/tylerlewy/Desktop/XR259/Plots/emap.cluster1.pdf"),
       dpi = 300, units = c("cm"),width = 30, height = 30)


# Repeat for cluster 2
to.emap.c2 <- pairwise_termsim(ego.c2)

emap.plot.c2 <- emapplot(to.emap.c2, 
                         showCategory = 15, 
                         group = T, 
                         node_label = "category",
                         layout = "fr", 
                         group_style = "ggforce",
                         label_group_style = "ggforce", 
                         nCluster = 4,
                         nWords = 5, 
                         min_edge = 0.3, 
                         label_format = 10)
emap.plot.c2

ggsave(emap.plot.c2,
       file=paste0("/Users/tylerlewy/Desktop/XR259/Plots/emap.cluster2.pdf"),
       dpi = 300, units = c("cm"),width = 30, height = 30)
