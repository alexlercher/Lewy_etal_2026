# Script to take raw serum cytokine values from TL338 and TL359 create a plot depicting PI Score

# Load in required libraries
library(data.table)
library(tidyverse)
library(patchwork)
library(tidyheatmaps)
library(UpSetR)

# Load in the serum cytokine data

# Data that is imported is a merged table of Luminex cytokine concentrations from TL338 and TL359

# Column containing 6CKine data was omitted because nearly all samples were above the limit of detection

# Sample Vehicle 7 was omitted because there was error with sample. 

# Data for IFNɑ was determined not by Luminex but by ELISA using a separate cohort of mice. 

concs <- fread("/Users/tylerlewy/Desktop/Dissertation/Combined Figures/Luminex/Serum_Cytokine_Combined_Conc.csv")

# Create a trimmed matrix to make calculations easier
concs.trimmed <- concs[,!c(1:3)]

# Create a separate table to store mean data and various stat tests
# Start with cytokine means
results <- tibble(
  Cytokine = colnames(concs.trimmed),
  `Vehicle Mean` = colMeans(subset(concs, Treatment == "Vehicle")[,!c(1:3)]), 
  `Poly(I:C) Mean` = colMeans(subset(concs, Treatment == "Poly(I:C)")[,!c(1:3)]), 
  `LPS Mean` = colMeans(subset(concs, Treatment == "LPS")[,!c(1:3)]),
  `ODN Mean` = colMeans(subset(concs, Treatment == "ODN")[,!c(1:3)])
)

# Calculate the fold changes
results <- results%>%
  mutate(`Poly(I:C) FC` = `Poly(I:C) Mean`/`Vehicle Mean`)%>%
  mutate(`LPS FC` = `LPS Mean`/`Vehicle Mean`)%>%
  mutate(`ODN FC` = `ODN Mean`/`Vehicle Mean`)

# Calculate the p values using Student's T Test
# Need to convert the tibbles to dataframes in order to use in the loop
conc.vehicle <- as.data.frame(subset(concs, Treatment == "Vehicle"))
conc.poly <- as.data.frame(subset(concs, Treatment == "Poly(I:C)"))
conc.lps <- as.data.frame(subset(concs, Treatment == "LPS"))
conc.odn <- as.data.frame(subset(concs, Treatment == "ODN"))

# Create an empty list to write into
t.test.pic <- list()
t.test.lps <- list()
t.test.odn <- list()

# Loop to run the t test function for each of the three conditions
for (i in 4:dim(concs)[2]) {
  cytokine = colnames(concs)[i]
  temp.pic <- filter(concs, Treatment %in% c("Vehicle", "Poly(I:C)"))
  t.test.pic[[i-3]] <- t.test(x = conc.vehicle[,i], y = conc.poly[,i])
  t.test.lps[[i-3]] <- t.test(x = conc.vehicle[,i], y = conc.lps[,i])
  t.test.odn[[i-3]] <- t.test(x = conc.vehicle[,i], y = conc.odn[,i])
}

# Pull the p values out as a named vector and adjust for multiple t tests
results <- results%>%
  mutate(`Poly(I:C) pval` = sapply(t.test.pic, function(res) res$p.value))%>%
  mutate(`LPS pval` = sapply(t.test.lps, function(res) res$p.value))%>%
  mutate(`ODN pval` = sapply(t.test.odn, function(res) res$p.value))%>%
  mutate(`Poly(I:C) padj` = p.adjust(`Poly(I:C) pval`, method = "BH"))%>%
  mutate(`LPS padj` = p.adjust(`LPS pval`, method = "BH"))%>%
  mutate(`ODN padj` = p.adjust(`ODN pval`, method = "BH"))

# Calculate a PI Score [-log10(padj) * log2(FC)]
results <- results%>%
  mutate(`Poly(I:C) PI` = -log10(results$`Poly(I:C) padj`) * log2(results$`Poly(I:C) FC`))%>%
  mutate(`LPS PI` = -log10(results$`LPS padj`) * log2(results$`LPS FC`))%>%
  mutate(`ODN PI` = -log10(results$`ODN padj`) * log2(results$`ODN FC`))

# Plot the results
# Need to format the table into a longer format
toPlot <- results %>%
  dplyr::select(c("Cytokine", "Poly(I:C) PI", "LPS PI", "ODN PI"))%>%
  pivot_longer(cols = !Cytokine, names_to = "Treatment", values_to = "PI Score")%>%
  mutate(Treatment = str_sub(Treatment, start = 1, end = -4))%>%
  mutate(Treatment = fct_relevel(Treatment, "Poly(I:C)", "LPS", "ODN"))%>%
  arrange(desc(`PI Score`))
  
# I will create a heatmap and extract the clustered row names
serum.hm <- tidy_heatmap(toPlot, 
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
ggsave("/Users/tylerlewy/Desktop/Dissertation/Combined Figures/Luminex/serum_piscores.svg", 
       serum.hm)

# Plot an upset plot to quantify shared vs unique responses
# Modify the results tibble to include a PI Threshold binary
toPlot.upset <- results %>%
  select(c(Cytokine, `Poly(I:C) PI`, `LPS PI`, `ODN PI`))%>%
  mutate(`Poly(I:C)` = if_else(`Poly(I:C) PI` > 4, 1, 0))%>%
  mutate(LPS = if_else(`LPS PI` > 4, 1, 0))%>%
  mutate(ODN = if_else(`ODN PI` > 4, 1, 0))%>%
  select(c(Cytokine, `Poly(I:C)`, LPS, ODN))

toPlot.upset <- as.data.frame(toPlot.upset)

upset.Serum <- upset(toPlot.upset, 
      nsets = 3,
      sets = c("ODN", "LPS", "Poly(I:C)"),
      keep.order = TRUE,
      order.by = "freq", 
      empty.intersections = TRUE,
      text.scale = 3,
      sets.bar.color = c("black","black","black"), 
      point.size = 6,
      set_size.scale_max = 15
      #ggtitle("Serum Cytokine Induction", subtitle = "PI Score > 4")
      )

svg("/Users/tylerlewy/Desktop/Dissertation/Combined Figures/Luminex/upset_serum.svg")
upset.Serum
dev.off()
