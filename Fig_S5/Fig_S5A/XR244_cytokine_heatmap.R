# Script to take raw serum cytokine values from XR244 create a plot depicting PI Score

# Load in required libraries
library(data.table)
library(tidyverse)
library(patchwork)
library(tidyheatmaps)
library(UpSetR)

# Load in the serum cytokine data

# Data that is imported is a table of Luminex cytokine concentrations from XR244

# Sample Vehicle 3 was omitted because there was some error with sample.

concs <- fread("/Users/tylerlewy/Desktop/XR244/XR244_1_serum_data_formatted_TL.csv")

# Create a trimmed matrix to make calculations easier
concs.trimmed <- concs[,!c(1:3)]

# Create a separate table to store mean data and various stat tests
# Start with cytokine means
results <- tibble(
  Cytokine = colnames(concs.trimmed),
  `Vehicle Mean` = colMeans(subset(concs, Treatment == "Vehicle")[,!c(1:3)]), 
  `Poly(I:C) Mean` = colMeans(subset(concs, Treatment == "Poly(I:C)")[,!c(1:3)]), 
  `DPI1 Mean` = colMeans(subset(concs, Treatment == "WNV 1DPI")[,!c(1:3)]),
  `DPI2 Mean` = colMeans(subset(concs, Treatment == "WNV 2DPI")[,!c(1:3)])
)

# Calculate the fold changes
results <- results%>%
  mutate(`Poly(I:C) FC` = `Poly(I:C) Mean`/`Vehicle Mean`)%>%
  mutate(`DPI1 FC` = `DPI1 Mean`/`Vehicle Mean`)%>%
  mutate(`DPI2 FC` = `DPI2 Mean`/`Vehicle Mean`)

# Calculate the p values using Student's T Test
# Need to convert the tibbles to dataframes in order to use in the loop
conc.vehicle <- as.data.frame(subset(concs, Treatment == "Vehicle"))
conc.poly <- as.data.frame(subset(concs, Treatment == "Poly(I:C)"))
conc.dpi1 <- as.data.frame(subset(concs, Treatment == "WNV 1DPI"))
conc.dpi2 <- as.data.frame(subset(concs, Treatment == "WNV 2DPI"))

# Create an empty list to write into
t.test.pic <- list()
t.test.dpi1 <- list()
t.test.dpi2 <- list()

# Loop to run the t test function for each of the three conditions
for (i in 4:dim(concs)[2]) {
  cytokine = colnames(concs)[i]
  temp.pic <- filter(concs, Treatment %in% c("Vehicle", "Poly(I:C)"))
  t.test.pic[[i-3]] <- t.test(x = conc.vehicle[,i], y = conc.poly[,i])
  t.test.dpi1[[i-3]] <- t.test(x = conc.vehicle[,i], y = conc.dpi1[,i])
  t.test.dpi2[[i-3]] <- t.test(x = conc.vehicle[,i], y = conc.dpi2[,i])
}

# Pull the p values out as a named vector and adjust for multiple t tests
results <- results%>%
  mutate(`Poly(I:C) pval` = sapply(t.test.pic, function(res) res$p.value))%>%
  mutate(`DPI1 pval` = sapply(t.test.dpi1, function(res) res$p.value))%>%
  mutate(`DPI2 pval` = sapply(t.test.dpi2, function(res) res$p.value))%>%
  mutate(`Poly(I:C) padj` = p.adjust(`Poly(I:C) pval`, method = "BH"))%>%
  mutate(`DPI1 padj` = p.adjust(`DPI1 pval`, method = "BH"))%>%
  mutate(`DPI2 padj` = p.adjust(`DPI2 pval`, method = "BH"))

# Calculate a PI Score [-log10(padj) * log2(FC)]
results <- results%>%
  mutate(`Poly(I:C) PI` = -log10(results$`Poly(I:C) padj`) * log2(results$`Poly(I:C) FC`))%>%
  mutate(`DPI1 PI` = -log10(results$`DPI1 padj`) * log2(results$`DPI1 FC`))%>%
  mutate(`DPI2 PI` = -log10(results$`DPI2 padj`) * log2(results$`DPI2 FC`))

# Plot the results
# Need to format the table into a longer format
toPlot <- results %>%
  dplyr::select(c("Cytokine", "Poly(I:C) PI", "DPI1 PI", "DPI2 PI"))%>%
  pivot_longer(cols = !Cytokine, names_to = "Treatment", values_to = "PI Score")%>%
  mutate(Treatment = str_sub(Treatment, start = 1, end = -4))%>%
  mutate(Treatment = fct_relevel(Treatment, "Poly(I:C)", "DPI1", "DPI2"))%>%
  arrange(desc(`PI Score`))

toPlot.wnv <- results %>%
  dplyr::select(c("Cytokine", "DPI1 PI", "DPI2 PI"))%>%
  pivot_longer(cols = !Cytokine, names_to = "Treatment", values_to = "PI Score")%>%
  mutate(Treatment = str_sub(Treatment, start = 1, end = -4))%>%
  mutate(Treatment = fct_relevel(Treatment, "DPI1", "DPI2"))%>%
  arrange(desc(`PI Score`))

# I will create a heatmap and extract the clustered row names
serum.all <- tidy_heatmap(toPlot, 
                          rows = Treatment, 
                          columns = Cytokine, 
                          values = `PI Score`, 
                          cluster_cols = F,
                          clustering_distance_cols = "euclidean",
                          cellwidth = 9,
                          cellheight = 12, 
                          border_color = "black", 
                          angle_col = 45,
                          legend = TRUE,
                          # scale = "row", 
                          colors = c("#F2F2F2", "red")
)

ggsave("/Users/tylerlewy/Desktop/XR244/Plots/serum_all_piscores.svg", 
       serum.all)
