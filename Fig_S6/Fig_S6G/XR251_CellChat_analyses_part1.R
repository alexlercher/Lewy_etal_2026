# Script to analyze snRNAseq dataset for L-R interactions
# Using vignette from https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
# as a guide

# 1) Setup -------------------------------------------------------------
## Load the required libraries
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)


# 2) Data Preprocessing -------------------------------------------------
# Read in the Seurat object
seurat_obj <- readRDS("C:\\Users\\Tyler\\Desktop\\snRNAseq\\seu_obj_v12192025_AL.rds")

# Make sure the default assay is RNA. Cellchat cannot use the SVT normalized values
DefaultAssay(seurat_obj) <- "RNA"

# Cellchat needs normalized counts in the RNA slot otherwise it will fail
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)

# Now we want to subset our object into the three individual conditions
# Pull out the condition names
conditions <- unique(seurat_obj$condition)

# Pull out each condition into its own list element
seurat_list <- lapply(conditions, function(cond) {
  subset(seurat_obj, subset = condition == cond)
})

# Label the list
names(seurat_list) <- conditions

# Remove the original object for now to save on memory
rm(seurat_obj)

# 3) Cellchat Processing ---------------------------------------

# For each condition we will perform the same pipeline
# We pull out data layer of the Seurat object
data.input <- GetAssayData(
  seurat_list[[1]],
  assay = "RNA",
  layer = "data"
)

# We pull out the metadata (in this case cell type and condition)
meta <- data.frame(
  cell_type = seurat_list[[1]]$cell_type,
  condition = seurat_list[[1]]$condition,
  row.names = colnames(data.input)
)

# We create the cellchat object and indicate cell_type as the variable of interest
cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "cell_type"
)

# Load in the mouse database
cellchat@DB <- CellChatDB.mouse

# Filters expression data to include only genes present in the dataset
cellchat <- subsetData(cellchat)

# Identify genes that are overexpressed in some cell types
cellchat <- identifyOverExpressedGenes(cellchat)

# Filters database to include only L-R pairs that are present in the dataset
cellchat <- identifyOverExpressedInteractions(cellchat)

# Calculates the probability of each L-R pair between each indicated cell type
cellchat <- computeCommunProb(cellchat,
                              type = "thresholdedMean")

# Removes probabilities for relationships with 10 or fewer cells
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Combines individual interactions into pathway level score
cellchat <- computeCommunProbPathway(cellchat)

# Combines all pair probabilities into pathway level scores for each cell type pair
cellchat <- aggregateNet(cellchat)

# Add the results to a list of cellchat results
cellchat_list <- list()
cellchat_list[[1]] <- cellchat

# Save the cellchat object
saveRDS(cellchat_list[[1]], file = "C:\\Users\\Tyler\\Desktop\\snRNAseq\\outputs\\cellchat_00H.rds")

# Remove the temporary object and repeat for the other two conditions
rm(cellchat)

# RUN FOR 02H CONDITION
# For each condition we will perform the same pipeline
# We pull out data layer of the Seurat object
data.input <- GetAssayData(
  seurat_list[[2]],
  assay = "RNA",
  layer = "data"
)

# We pull out the metadata (in this case cell type and condition)
meta <- data.frame(
  cell_type = seurat_list[[2]]$cell_type,
  condition = seurat_list[[2]]$condition,
  row.names = colnames(data.input)
)

# We create the cellchat object and indicate cell_type as the variable of interest
cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "cell_type"
)

# Load in the mouse database
cellchat@DB <- CellChatDB.mouse

# Filters expression data to include only genes present in the dataset
cellchat <- subsetData(cellchat)

# Identify genes that are overexpressed in some cell types
cellchat <- identifyOverExpressedGenes(cellchat)

# Filters database to include only L-R pairs that are present in the dataset
cellchat <- identifyOverExpressedInteractions(cellchat)

# Calculates the probability of each L-R pair between each indicated cell type
cellchat <- computeCommunProb(cellchat,
                              type = "thresholdedMean")

# Removes probabilities for relationships with 10 or fewer cells
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Combines individual interactions into pathway level score
cellchat <- computeCommunProbPathway(cellchat)

# Combines all pair probabilities into pathway level scores for each cell type pair
cellchat <- aggregateNet(cellchat)

# Add the results to a list of cellchat results
cellchat_list[[2]] <- cellchat

# Save the cellchat object
saveRDS(cellchat_list[[2]], file = "C:\\Users\\Tyler\\Desktop\\snRNAseq\\outputs\\cellchat_02H.rds")


# Remove the temporary object and repeat for last condition
rm(cellchat)


# ANALYZE 16H DATASET
# For each condition we will perform the same pipeline
# We pull out data layer of the Seurat object
data.input <- GetAssayData(
  seurat_list[[3]],
  assay = "RNA",
  layer = "data"
)

# We pull out the metadata (in this case cell type and condition)
meta <- data.frame(
  cell_type = seurat_list[[3]]$cell_type,
  condition = seurat_list[[3]]$condition,
  row.names = colnames(data.input)
)

# We create the cellchat object and indicate cell_type as the variable of interest
cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "cell_type"
)

# Load in the mouse database
cellchat@DB <- CellChatDB.mouse

# Filters expression data to include only genes present in the dataset
cellchat <- subsetData(cellchat)

# Identify genes that are overexpressed in some cell types
cellchat <- identifyOverExpressedGenes(cellchat)

# Filters database to include only L-R pairs that are present in the dataset
cellchat <- identifyOverExpressedInteractions(cellchat)

# Calculates the probability of each L-R pair between each indicated cell type
cellchat <- computeCommunProb(cellchat,
                              type = "thresholdedMean")

# Removes probabilities for relationships with 10 or fewer cells
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Combines individual interactions into pathway level score
cellchat <- computeCommunProbPathway(cellchat)

# Combines all pair probabilities into pathway level scores for each cell type pair
cellchat <- aggregateNet(cellchat)

# Add the results to a list of cellchat results
cellchat_list[[3]] <- cellchat

# Remove the temporary object
rm(cellchat)

# Save the cellchat object
saveRDS(cellchat_list[[3]], file = "C:\\Users\\Tyler\\Desktop\\snRNAseq\\outputs\\cellchat_16H.rds")