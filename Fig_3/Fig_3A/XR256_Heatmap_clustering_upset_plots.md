XR256_Heatmap_clustering_Upset_plot
================
Alexander Lercher
2026-04-29

## LOAD PACKAGES

``` r
library(readr)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ purrr     1.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ ggplot2   4.0.1     ✔ tibble    3.3.0
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.2
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(dplyr)
library(pheatmap)
library(ComplexUpset)
library(ggplot2)
library(grid)
```

## DATA IMPORT AND CLEANUP

``` r
list_of_files <- list.files(path = "../XR256_DESeq2/output",
                            recursive = FALSE,
                            pattern = "DEG.tsv",
                            full.names = TRUE)

data <- readr::read_tsv(list_of_files, id = NULL)
```

    ## Rows: 114059 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): gene_id, comparison_name
    ## dbl (7): comparison, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
data_fpkm <- read.delim("../XR256_DESeq2/output/XR256_FPKM.tsv")

# what are the comparisons made
data_summary <- data %>%
  group_by(comparison) %>%  
  distinct(comparison_name)
```

## FILTER FOR AND EXTRACT DEGs

``` r
# define log2FC and adjPval cutoffs
log2FC_cutoff = 1
adjpval_cutoff = 0.05

# subset data by cutoff criteria
significant_DEG <- data %>%
  group_by(comparison) %>%
  filter(abs(log2FoldChange) > log2FC_cutoff & padj < adjpval_cutoff)

# count significant genes per condition
count_DEG <-  significant_DEG %>%
  group_by(comparison, comparison_name) %>%
  count(nrow(comparison)) %>%
  dplyr::rename(DEG_count = n) %>%
  print()
```

    ## # A tibble: 6 × 3
    ## # Groups:   comparison, comparison_name [6]
    ##   comparison comparison_name DEG_count
    ##        <dbl> <chr>               <int>
    ## 1          1 LPS_vs_PBS            454
    ## 2          2 ODN_vs_PBS            110
    ## 3          3 PIC_vs_PBS            389
    ## 4          4 PIC_vs_LPS            227
    ## 5          5 PIC_vs_ODN            299
    ## 6          6 LPS_vs_ODN            341

``` r
# put DEG of individual comparisons into list
list_DEG_comparisons <- list()
for(i in 1:length(count_DEG$comparison)){
  subset_comparison <- subset(significant_DEG, significant_DEG$comparison == i)
  comparison_name <- unique(subset_comparison$comparison_name)
  list_DEG_comparisons[[paste0(i,"_",comparison_name)]] <- subset_comparison
  rm(subset_comparison,comparison_name)
}
```

## CURATE FPKM DATA

``` r
# add means for all conditions
data_fpkm <- data_fpkm %>%
  mutate(
    mean_PBS  = rowMeans(dplyr::select(., starts_with("PBS"))),
    mean_LPS = rowMeans(dplyr::select(., starts_with("LPS"))),
    mean_ODN  = rowMeans(dplyr::select(., starts_with("ODN"))),
    mean_PIC = rowMeans(dplyr::select(., starts_with("PIC")))
  )

# rename column gene to gene_id column
data_fpkm <- data_fpkm %>%
  dplyr::rename(gene_id = gene)
```

## LOAD ISG DATA

``` r
# import ISG list from this study, generated in XR189_ISG_analyses
ISG_list  <- read.delim("../../XR249_TL396_iv_IFN-a_brain_RNAseq/XR249_RNAseq_analyses_AL/XR249_ISG_extraction/output/XR249_ISG_list_brain.tsv")

# add a TRUE column for all ISGs
ISG_list$ISG = T
```

## FPKM FILTER

``` r
# FPKM values only for filtering
fpkm_only <- data_fpkm %>%
  dplyr::select(!contains(c("mean")))

# calculate average FPKM level across all samples and define FPKM cutoff
fpkm_only <- fpkm_only %>%
  mutate(average_fpkm = rowSums(fpkm_only[2:length(fpkm_only)]/(length(fpkm_only)-1)))

FPKM_cutoff <- 1

fpkm_only %>%
  ggplot(aes(x=average_fpkm)) +
    geom_density(fill="#8a4655", color="#e9ecef", alpha=0.8) +
    scale_x_continuous(trans='log10') +
    ggtitle("Average FPKM before cutoff") +
    geom_vline(xintercept = FPKM_cutoff) +
    theme_bw()
```

![](XR256_Heatmap_clustering_upset_plots_files/figure-gfm/FPKM%20filter-1.png)<!-- -->

``` r
# filter according to FPKM cutoff
fpkm_only <- fpkm_only %>%
  filter(average_fpkm > FPKM_cutoff)

fpkm_only %>%
  ggplot(aes(x=average_fpkm)) +
  geom_density(fill="#8a4655", color="#e9ecef", alpha=0.8) +
  scale_x_continuous(trans='log10') +
  ggtitle("Average FPKM after cutoff") +
  geom_vline(xintercept = FPKM_cutoff) +
  theme_bw()
```

![](XR256_Heatmap_clustering_upset_plots_files/figure-gfm/FPKM%20filter-2.png)<!-- -->

``` r
# drop average FPKM column
fpkm_only <- fpkm_only %>%
  dplyr::select(!average_fpkm)

# loop over comparisons of interest
# 1) cluster genes
# 2) plot heatmaps
# 3) add cluster info to gene list

DEG_comparisons_cluster_list <- list()
DEG_comparisons_cluster_DF <- data.frame()

for(i in 1:length(list_DEG_comparisons)){
# select condition of interest
conditionOI <- list_DEG_comparisons[[i]]

# extract DEG from condition of interest
genesOI <- data.frame(gene_id = conditionOI$gene_id)

# get FPKM for DEG of interest
fpkmOI <- left_join(genesOI,fpkm_only)
fpkmOI <- na.omit(fpkmOI)
row.names(fpkmOI) <- fpkmOI$gene_id
fpkmOI$gene_id <- NULL

#--------------------------------------------------------------------
# PLOT HEATMAP with ROW Z SCORES
#--------------------------------------------------------------------
# define formula to calculate row Z score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# calculate row Z score for FPKMs for each row (gene) of interest
fpkmOI_norm <- t(apply(fpkmOI, 1, cal_z_score))

# annotate sample columns
my_sample_col <- data.frame(sample = c(rep("PBS",4),
                                       rep("LPS",4),
                                       rep("ODN",4),
                                       rep("PIC",4)),
                            condition = c(rep("control",4),
                                          rep("treated",12)))

row.names(my_sample_col) <- colnames(fpkmOI)

# calculate dendrogram
hc <- hclust(dist(fpkmOI_norm), method = "complete")
#as.dendrogram(hc) %>%  #no need to print dendrogram
#  plot(horiz = T)

# obtain gene names as per dendrogram
#genes_order <- rev(row.names(fpkmOI_norm)[hc$order]) #no need to print gene names

# add cluster information to heatmap
numberofclusters = 6
my_gene_col <- cutree(hc, k = numberofclusters)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "1",
                                           ifelse(test = my_gene_col == 2, yes = "2",
                                                  ifelse(test = my_gene_col == 3, yes = "3",
                                                         ifelse(test = my_gene_col == 4, yes = "4",
                                                                ifelse(test = my_gene_col == 5, yes = "5", no = "6"))))))

# merge ISG_list heatmap row/gene names to add ISG info
my_gene_col$gene_id <- rownames(my_gene_col)
my_gene_col <- my_gene_col %>%
  left_join(ISG_list) %>%
  replace(is.na(.), FALSE)
my_gene_col$ISG <-as.integer(my_gene_col$ISG)
rownames(my_gene_col) <- my_gene_col$gene_id
my_gene_col$gene_id <- NULL

pheatmap(fpkmOI_norm, annotation_col = my_sample_col, annotation_row = my_gene_col, fontsize_row = 2,
         show_rownames = T, main = paste0(i,":",as.character(data_summary[i,2])),
         #color = colorRampPalette(c("midnightblue", "ivory", "firebrick3"))(50),
         cutree_rows = numberofclusters,
         filename = paste0("output/clustering/",i,"_",as.character(data_summary[i,2]),".pdf"))

# add cluster and ISG information to DEG list
my_gene_col$gene_id <- rownames(my_gene_col) 
conditionOI_cluster <- conditionOI %>%
  left_join(my_gene_col) %>%
  arrange(cluster) %>%
  na.omit()
DEG_comparisons_cluster_list[[paste0(i,"_",as.character(data_summary[i,2]))]] <- conditionOI_cluster
DEG_comparisons_cluster_DF <- rbind(DEG_comparisons_cluster_DF,conditionOI_cluster)

write_tsv(conditionOI_cluster, paste0("output/clustering/",i,"_",as.character(data_summary[i,2]),"_cluster.tsv"))

}
```

    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`

## CLUSTER GENES SIGNIFICANT IN AT LEAST ONE CONDITION

``` r
# get all significant DEG for  comparisons
sig_DEG <- significant_DEG %>%
  group_by(comparison_name)

# reduce to unique significant DEG for comparisons
unique_sig_DEG <- data.frame(gene_id = unique(sig_DEG$gene_id))

# get FPKM for DEG of interest
fpkmOI <- left_join(unique_sig_DEG,fpkm_only)
```

    ## Joining with `by = join_by(gene_id)`

``` r
fpkmOI <- na.omit(fpkmOI)
row.names(fpkmOI) <- fpkmOI$gene_id
fpkmOI$gene_id <- NULL

# define formula to calculate row Z score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# calculate row Z score for FPKMs for each row (gene) of interest
fpkmOI_norm <- t(apply(fpkmOI, 1, cal_z_score))

# annotate sample columns
my_sample_col <- data.frame(sample = c(rep("PBS",4),
                                       rep("LPS",4),
                                       rep("ODN",4),
                                       rep("PIC",4)),
                            condition = c(rep("control",4),
                                          rep("treated",12)))

row.names(my_sample_col) <- colnames(fpkmOI)

# calculate dendrogram
hc <- hclust(dist(fpkmOI_norm), method = "complete")

# add cluster information to heatmap
numberofclusters = 6
my_gene_col <-cutree(hc, k = numberofclusters)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "1",
                                           ifelse(test = my_gene_col == 2, yes = "2",
                                                  ifelse(test = my_gene_col == 3, yes = "3",
                                                         ifelse(test = my_gene_col == 4, yes = "4",
                                                                ifelse(test = my_gene_col == 5, yes = "5", no = "6"))))))


# merge ISG_list heatmap row/gene names to add ISG info
my_gene_col$gene_id <- rownames(my_gene_col)
my_gene_col <- my_gene_col %>%
  left_join(ISG_list) %>%
  replace(is.na(.), FALSE)
```

    ## Joining with `by = join_by(gene_id)`

``` r
my_gene_col$ISG <-as.integer(my_gene_col$ISG)
rownames(my_gene_col) <- my_gene_col$gene_id
my_gene_col$gene_id <- NULL

pheatmap::pheatmap(fpkmOI_norm, annotation_col = my_sample_col, annotation_row = my_gene_col, 
           show_rownames = F, main = "significant in at least one comparison",
           #color = colorRampPalette(c("midnightblue", "ivory", "firebrick3"))(50),
           cutree_rows = numberofclusters)
```

![](XR256_Heatmap_clustering_upset_plots_files/figure-gfm/cluster%20genes%20significant%20in%20at%20least%20one%20condition-1.png)<!-- -->

``` r
# add cluster and ISG information to DEG list
my_gene_col$gene_id <- rownames(my_gene_col) 
sig_DEG_cluster <- sig_DEG %>%
  left_join(my_gene_col) %>%
  arrange(cluster) %>%
  na.omit()
```

    ## Joining with `by = join_by(gene_id)`

``` r
write_tsv(sig_DEG_cluster, "output/clustering/all_DEG_cluster.tsv")
```

## GENERATE UPSET PLOTS

``` r
# filter for conditions of interest
conditionsOI <- c("PBS")
deg_list <- DEG_comparisons_cluster_list[str_detect(names(list_DEG_comparisons), conditionsOI)]

# Extract gene_id vectors
gene_sets <- map(deg_list, ~ .x$gene_id)

# Build long to wide presence/absence matrix
df_wide <- imap_dfr(
  gene_sets,
  ~ tibble(gene_id = .x, set = .y)
) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = set,
    values_from = value,
    values_fill = 0
  )

# Plot Upset Plots
p <- ComplexUpset::upset(
  df_wide,
  intersect = names(df_wide)[-1],
  base_annotations = list(
    'Intersection size' = intersection_size(counts = TRUE)
  ),
  themes = list(
    default = theme_minimal(),
    annotation = theme_minimal(),
    intersection_matrix = theme_minimal(),
    intersections = theme_minimal(),
    overall_sizes = theme_minimal()
  )
) +
  ggtitle("Shared DEGs Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5))
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## ℹ The deprecated feature was likely used in the ComplexUpset package.
    ##   Please report the issue at
    ##   <https://github.com/krassowski/complex-upset/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
print(p)
```

![](XR256_Heatmap_clustering_upset_plots_files/figure-gfm/generate%20upset%20plots-1.png)<!-- -->

## APPENDIX

``` r
# which R packages and versions?
if ("devtools" %in% installed.packages()) devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.2 (2025-10-31)
    ##  os       macOS Sonoma 14.7.3
    ##  system   aarch64, darwin20
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       America/New_York
    ##  date     2026-04-29
    ##  pandoc   3.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
    ##  quarto   1.5.57 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package      * version date (UTC) lib source
    ##  bit            4.6.0   2025-03-06 [1] CRAN (R 4.5.0)
    ##  bit64          4.6.0-1 2025-01-16 [1] CRAN (R 4.5.0)
    ##  cachem         1.1.0   2024-05-16 [1] CRAN (R 4.5.0)
    ##  cli            3.6.5   2025-04-23 [1] CRAN (R 4.5.0)
    ##  colorspace     2.1-2   2025-09-22 [1] CRAN (R 4.5.0)
    ##  ComplexUpset * 1.3.3   2021-12-11 [1] CRAN (R 4.5.0)
    ##  crayon         1.5.3   2024-06-20 [1] CRAN (R 4.5.0)
    ##  devtools       2.4.6   2025-10-03 [1] CRAN (R 4.5.0)
    ##  digest         0.6.39  2025-11-19 [1] CRAN (R 4.5.2)
    ##  dplyr        * 1.1.4   2023-11-17 [1] CRAN (R 4.5.0)
    ##  ellipsis       0.3.2   2021-04-29 [1] CRAN (R 4.5.0)
    ##  evaluate       1.0.5   2025-08-27 [1] CRAN (R 4.5.0)
    ##  farver         2.1.2   2024-05-13 [1] CRAN (R 4.5.0)
    ##  fastmap        1.2.0   2024-05-15 [1] CRAN (R 4.5.0)
    ##  forcats      * 1.0.1   2025-09-25 [1] CRAN (R 4.5.0)
    ##  fs             1.6.6   2025-04-12 [1] CRAN (R 4.5.0)
    ##  generics       0.1.4   2025-05-09 [1] CRAN (R 4.5.0)
    ##  ggplot2      * 4.0.1   2025-11-14 [1] CRAN (R 4.5.2)
    ##  glue           1.8.0   2024-09-30 [1] CRAN (R 4.5.0)
    ##  gtable         0.3.6   2024-10-25 [1] CRAN (R 4.5.0)
    ##  hms            1.1.4   2025-10-17 [1] CRAN (R 4.5.0)
    ##  htmltools      0.5.9   2025-12-04 [1] CRAN (R 4.5.2)
    ##  knitr          1.51    2025-12-20 [1] CRAN (R 4.5.2)
    ##  labeling       0.4.3   2023-08-29 [1] CRAN (R 4.5.0)
    ##  lifecycle      1.0.4   2023-11-07 [1] CRAN (R 4.5.0)
    ##  lubridate    * 1.9.4   2024-12-08 [1] CRAN (R 4.5.0)
    ##  magrittr       2.0.4   2025-09-12 [1] CRAN (R 4.5.0)
    ##  memoise        2.0.1   2021-11-26 [1] CRAN (R 4.5.0)
    ##  otel           0.2.0   2025-08-29 [1] CRAN (R 4.5.0)
    ##  patchwork      1.3.2   2025-08-25 [1] CRAN (R 4.5.0)
    ##  pheatmap     * 1.0.13  2025-06-05 [1] CRAN (R 4.5.0)
    ##  pillar         1.11.1  2025-09-17 [1] CRAN (R 4.5.0)
    ##  pkgbuild       1.4.8   2025-05-26 [1] CRAN (R 4.5.0)
    ##  pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.5.0)
    ##  pkgload        1.4.1   2025-09-23 [1] CRAN (R 4.5.0)
    ##  purrr        * 1.2.0   2025-11-04 [1] CRAN (R 4.5.0)
    ##  R6             2.6.1   2025-02-15 [1] CRAN (R 4.5.0)
    ##  RColorBrewer   1.1-3   2022-04-03 [1] CRAN (R 4.5.0)
    ##  readr        * 2.1.6   2025-11-14 [1] CRAN (R 4.5.2)
    ##  remotes        2.5.0   2024-03-17 [1] CRAN (R 4.5.0)
    ##  rlang          1.1.6   2025-04-11 [1] CRAN (R 4.5.0)
    ##  rmarkdown      2.30    2025-09-28 [1] CRAN (R 4.5.0)
    ##  rstudioapi     0.17.1  2024-10-22 [1] CRAN (R 4.5.0)
    ##  S7             0.2.1   2025-11-14 [1] CRAN (R 4.5.2)
    ##  scales         1.4.0   2025-04-24 [1] CRAN (R 4.5.0)
    ##  sessioninfo    1.2.3   2025-02-05 [1] CRAN (R 4.5.0)
    ##  stringi        1.8.7   2025-03-27 [1] CRAN (R 4.5.0)
    ##  stringr      * 1.6.0   2025-11-04 [1] CRAN (R 4.5.0)
    ##  tibble       * 3.3.0   2025-06-08 [1] CRAN (R 4.5.0)
    ##  tidyr        * 1.3.2   2025-12-19 [1] CRAN (R 4.5.2)
    ##  tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.5.0)
    ##  tidyverse    * 2.0.0   2023-02-22 [1] CRAN (R 4.5.0)
    ##  timechange     0.3.0   2024-01-18 [1] CRAN (R 4.5.0)
    ##  tzdb           0.5.0   2025-03-15 [1] CRAN (R 4.5.0)
    ##  usethis        3.2.1   2025-09-06 [1] CRAN (R 4.5.0)
    ##  utf8           1.2.6   2025-06-08 [1] CRAN (R 4.5.0)
    ##  vctrs          0.6.5   2023-12-01 [1] CRAN (R 4.5.0)
    ##  vroom          1.6.7   2025-11-28 [1] CRAN (R 4.5.2)
    ##  withr          3.0.2   2024-10-28 [1] CRAN (R 4.5.0)
    ##  xfun           0.55    2025-12-16 [1] CRAN (R 4.5.2)
    ##  yaml           2.3.12  2025-12-10 [1] CRAN (R 4.5.2)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
