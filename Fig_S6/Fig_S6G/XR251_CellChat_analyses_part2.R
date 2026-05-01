# Script to perform comparative analysis of previously generated cellchat objects
# Used cellchat pipeline to create individual cellchat objects for snRNAseq brain 
# samples from animals harvested 02 or 16h post IV IFNa or vehicle controls
# 
#
# 1) Setup ----------------------------------------------------------------
# Load required libraries
library(CellChat) # 2.2.0
library(patchwork) # 1.3.1
library(ComplexHeatmap) # 2.26.0

rm(list = ls())

setwd("/Users/lewytyg/Desktop/IFN_snRNAseq/plots/")

# 2) Data Import -----------------------------------------------------------
# Load the individual cellchat objects
cellchat.00 <- readRDS("/Users/lewytyg/Desktop/IFN_snRNAseq/objects/cellchat_00H.rds")
cellchat.02 <- readRDS("/Users/lewytyg/Desktop/IFN_snRNAseq/objects/cellchat_02H.rds")
cellchat.16 <- readRDS("/Users/lewytyg/Desktop/IFN_snRNAseq/objects/cellchat_16H.rds")

# Combine those objects into a list then use it to merge using mergeCellChat() 
cellchat.list <- list("00h" = cellchat.00, 
                      "02h" = cellchat.02, 
                      "16h" = cellchat.16)
cellchat.merged <- mergeCellChat(cellchat.list, 
                                 add.names = names(cellchat.list))
cellchat.merged

# Circle Plot to display outgoing microglia signals
par(mfrow = c(1,1), 
    xpd=TRUE,
    oma = c(0, 0, 3, 0))
netVisual_diffInteraction(cellchat.merged, 
                          weight.scale = T, 
                          comparison = c(1,3),
                          sources.use = 8)
mtext(
  "Outgoing signals microglia 16h vs. ctrl",
  side = 3,
  outer = TRUE,
  line = 1,
  cex = 1.4,
  font = 2
)

# Alternatively we can compare communication through DE analysis
# Define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "16h"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat.merged <- identifyOverExpressedGenes(cellchat.merged, 
                                              group.dataset = "datasets", 
                                              pos.dataset = pos.dataset, 
                                              features.name = features.name, 
                                              only.pos = FALSE, 
                                              thresh.pc = 0.1, 
                                              thresh.fc = 0.1, 
                                              thresh.p = 0.05)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat.merged, 
                     features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat.merged, 
                              net = net, 
                              datasets = "16h",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat.merged, 
                                net = net, 
                                datasets = "00h",
                                ligand.logFC = -0.1, 
                                receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat.merged)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat.merged)

pairLR.use.up = net.up[, "interaction_name", drop = F]
# Differential Microglia --> Neuron Signaling
gg1 <- netVisual_bubble(cellchat.merged, 
                        pairLR.use = pairLR.use.up, 
                        sources.use = 8, 
                        targets.use = c(3,6,7), 
                        comparison = c(1, 3),  
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Up-regulated Microglia --> Neuron signaling at ", 
                                            names(cellchat.list)[3])) + 
  coord_flip()
gg1
# Microglia --> Excitatory + Inhibitory + Interneurons - 16H vs 00H
gg1 <- rankNet(cellchat.merged, 
               mode = "comparison", 
               comparison = c(1,3),
               sources.use = c(8),
               signaling = c("NRG", "Glutamate", "NRXN", "NEGR", "PTPR", "ADGRL", "GABA-A"),
               targets.use = c(3,6,7),
               stacked = T, 
               do.stat = TRUE)
gg2 <- rankNet(cellchat.merged, 
               mode = "comparison", 
               comparison = c(1,3),
               sources.use = c(8),
               signaling = c("NRG", "Glutamate", "NRXN", "NEGR", "PTPR", "ADGRL", "GABA-A"),
               targets.use = c(3,6,7),
               stacked = F, 
               do.stat = TRUE)
gg1 + gg2

