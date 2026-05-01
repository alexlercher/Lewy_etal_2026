# Load libraries
library(data.table)
library(ggplot2)
library(biomaRt)
library(ggrepel)

# Load data
data <- fread("~/Desktop/Paper Revisions/1_BMEC_IFNa_vs_BMEC_PBS_DEG.tsv")

# Load mart
mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

# Filter on GO terms
go_terms <- c("GO:0005125", # Cytokine Activity
              "GO:0008009" # Chemokine Activity
              )

go_genes <- getBM(
  attributes = c(
    "external_gene_name",
    "ensembl_gene_id",
    "go_id"
  ),
  filters = "go",
  values = go_terms,
  mart = mart
)

cytokine.genes <- unique(go_genes$external_gene_name)

toPlot <- data%>%
  filter(gene_id %in% cytokine.genes)%>%
  mutate(sig = if_else(padj < 0.05 & abs(log2FoldChange) > 1, "Y", "N"))

## 1) Compute symmetric x limits
x_lim <- max(abs(toPlot$log2FoldChange), na.rm = TRUE)

## 2) Define y cutoff
padj_cutoff <- -log10(0.05)

ggplot(
  toPlot,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    colour = sig
  )
) +
  ## points
  geom_point(size = 1.5, alpha = 0.8) +
  
  ## labels: significant genes only
  geom_label_repel(
    data = filter(toPlot, sig == "Y"),
    aes(label = gene_id),
    colour = "black",
    label.colour = NA,
    size = 3,
    max.overlaps = Inf
  ) +
  
  ## 3) Cutoff lines
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.4,
    colour = "black"
  ) +
  geom_hline(
    yintercept = padj_cutoff,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "black"
  ) +
  
  ## 4) Symmetric x-axis
  scale_x_continuous(
    limits = c(-x_lim, x_lim)
  ) +
  
  ## color mapping
  scale_colour_manual(
    values = c(
      "N" = "grey70",
      "Y"  = "red"
    )
  ) +
  
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  ) +
  labs(
    x = "log2 Fold Change",
    y = expression(-log[10](adjusted~p))
  )+ 
  ggtitle(label = "BMEC Cytokine Response to IFNa",
          subtitle = "GO:0005125 - Cytokine Activity \nGO:0008009 - Chemokine Activity")


 ggsave("~/Desktop/Paper Revisions/BMEC_IFN_Cytokine_Induction_3.pdf")
