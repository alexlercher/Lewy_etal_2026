XR261_DESeq2_top500
================
Alexander Lercher
2026-04-29

## LOAD PACKAGES

``` r
library(tidyr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(readr)
library(stringr)
library(ggplot2)
library(gplots)
```

    ## 
    ## ---------------------
    ## gplots 3.3.0 loaded:
    ##   * Use citation('gplots') for citation info.
    ##   * Homepage: https://talgalili.github.io/gplots/
    ##   * Report issues: https://github.com/talgalili/gplots/issues
    ##   * Ask questions: https://stackoverflow.com/questions/tagged/gplots
    ##   * Suppress this message with: suppressPackageStartupMessages(library(gplots))
    ## ---------------------

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
library(RColorBrewer)
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     explain

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:gplots':
    ## 
    ##     space

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: GenomicRanges

    ## Loading required package: Seqinfo

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(pheatmap)
library(edgeR)
```

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:DESeq2':
    ## 
    ##     plotMA

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

## DATA IMPORT AND CLEANUP

``` r
# import table with reads of all samples
data <- read_csv("input/raw_counts.csv")
```

    ## Rows: 41079 Columns: 20
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (1): Gene.name
    ## dbl (19): Length, brain.PBS.1, brain.PBS.2, brain.PBS.3, brain.PBS.4, brain....
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# remove columns not needed
data_clean <- data %>%
  dplyr::select(-Length)

# Rename columng Gene.name to gene, make it first column and remove NA genes
data_clean <- data_clean %>%
  dplyr::rename(gene = Gene.name) %>%
  dplyr::select(gene, everything()) %>%
  filter(!is.na(gene))

# remove duplicates, keep the one with highest counts
# Group by gene name, calculate the sum of counts across samples, and keep the row with the highest sum
data_clean <- data_clean %>%
  rowwise() %>%
  mutate(sum_except_one = sum(c_across(-gene))) %>%
  group_by(gene) %>%
  arrange(desc(sum_except_one)) %>%
  slice_head() %>%
  ungroup() %>%
  dplyr::select(!sum_except_one)

# define data_clean as data_raw
data_raw <- data_clean

# make gene names row names and remove gene column
data_raw <- as.data.frame(data_raw)
rownames(data_raw) <- as.vector(data_raw$gene)
data_raw$gene <- NULL
rm(data,data_clean)
```

## IMPORT METADATA AND COMPARISONS OF INTEREST

``` r
# import table with metadata of all samples
metadata <- read.table("input/metadata.txt", header = T)

# make sure metadata row order is same as column names of data
# get desired order of row names
desired_order <- names(data_raw)

# reorder the rows based on the desired order
metadata <- metadata %>%
  arrange(match(variable, desired_order))

# import table with groups you want to compare with
# add column with comparison name
# add column with comparison number
DE_groups <- read.table("input/group_comparisons.txt", header = T)

comparison_name <- list()
comparison_number <- list()
for(i in 1:nrow(DE_groups)){
  x <- paste0(DE_groups[i,1],"_vs_",DE_groups[i,2])
  comparison_name[[i]] <- x
  comparison_number[[i]] <- i
}
DE_groups$Comparison <- comparison_name
DE_groups$ComparisonNumber <- comparison_number
rm(comparison_name, comparison_number)
```

## RUN DESEQ2

``` r
# DESeq analyses and create DESeq object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(data_raw),
                              colData = metadata,
                              design = ~ description)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

## GET TOP 500 VARIABLE GENES

``` r
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)
gene_var <- apply(vsd_mat, 1, var)

top_n_genes <- 500
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n_genes]
vsd_top <- vsd_mat[top_genes, ]
```

## PCA FOR TOP 500 VARIABLE GENES

``` r
p <- plotPCA(vsd, intgroup = "description", ntop = 500) +
  scale_color_brewer(palette = "Set1") + 
  theme_bw() + 
  ggtitle("PCA top 500 variable genes") +
  theme(aspect.ratio = 1) +
  geom_text(aes(label = name), vjust = -0.5)
```

    ## using ntop=500 top features by variance

``` r
print(p)
```

![](XR261_DESeq2_top500_files/figure-gfm/PCA%20for%20top%20500%20variable%20genes-1.png)<!-- -->

## HEATMAP TOP 500 GENES

``` r
cluster_count <- 3
# Compute correlation distance for rows (genes) and clustering
dist_rows <- as.dist(1 - cor(t(vsd_top), method = "pearson"))
hc_rows <- hclust(dist_rows, method = "average")

# cut into 3 clusters
gene_clusters <- cutree(hc_rows, k = cluster_count)

# calculate cpm
cpm_mat <- cpm(counts(dds), log = FALSE)
cpm_top <- cpm_mat[top_genes, ]

# save to tsv
cluster_df <- data.frame(
  gene      = names(gene_clusters),
  cluster   = gene_clusters
)

cluster_df <- cbind(cluster_df, as.data.frame(cpm_top[names(gene_clusters), ]))

write.table(cluster_df,
            file = paste0("output/gene_clusters_top",as.character(top_n_genes),"_genes.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# prepare row annotation to use my calculated clusters
row_anno <- data.frame(Cluster = factor(gene_clusters[top_genes]))
rownames(row_anno) <- top_genes

# plot heatmap
p <- pheatmap::pheatmap(vsd_top,
             scale = "row",
             show_rownames = FALSE,
             annotation_col = as.data.frame(colData(dds)[, "description", drop = FALSE]),
             annotation_row = row_anno,
             cluster_rows = hc_rows,
             main = paste0("Clustering of Top ",as.character(top_n_genes)," genes"),
             cluster_cols = TRUE,
             fontsize_col = 10)
```

![](XR261_DESeq2_top500_files/figure-gfm/heatmap%20of%20top%20500%20genes-1.png)<!-- -->

``` r
grid::grid.newpage()
grid::grid.draw(p$gtable)
```

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
    ##  package              * version  date (UTC) lib source
    ##  abind                  1.4-8    2024-09-12 [1] CRAN (R 4.5.0)
    ##  Biobase              * 2.70.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocGenerics         * 0.56.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocParallel           1.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  bit                    4.6.0    2025-03-06 [1] CRAN (R 4.5.0)
    ##  bit64                  4.6.0-1  2025-01-16 [1] CRAN (R 4.5.0)
    ##  bitops                 1.0-9    2024-10-03 [1] CRAN (R 4.5.0)
    ##  cachem                 1.1.0    2024-05-16 [1] CRAN (R 4.5.0)
    ##  caTools                1.18.3   2024-09-04 [1] CRAN (R 4.5.0)
    ##  cli                    3.6.5    2025-04-23 [1] CRAN (R 4.5.0)
    ##  codetools              0.2-20   2024-03-31 [1] CRAN (R 4.5.2)
    ##  crayon                 1.5.3    2024-06-20 [1] CRAN (R 4.5.0)
    ##  DelayedArray           0.36.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  DESeq2               * 1.50.2   2025-11-17 [1] Bioconductor 3.22 (R 4.5.2)
    ##  devtools               2.4.6    2025-10-03 [1] CRAN (R 4.5.0)
    ##  digest                 0.6.39   2025-11-19 [1] CRAN (R 4.5.2)
    ##  dplyr                * 1.1.4    2023-11-17 [1] CRAN (R 4.5.0)
    ##  edgeR                * 4.8.2    2025-12-23 [1] https://bioc-release.r-universe.dev (R 4.5.2)
    ##  ellipsis               0.3.2    2021-04-29 [1] CRAN (R 4.5.0)
    ##  evaluate               1.0.5    2025-08-27 [1] CRAN (R 4.5.0)
    ##  farver                 2.1.2    2024-05-13 [1] CRAN (R 4.5.0)
    ##  fastmap                1.2.0    2024-05-15 [1] CRAN (R 4.5.0)
    ##  fs                     1.6.6    2025-04-12 [1] CRAN (R 4.5.0)
    ##  generics             * 0.1.4    2025-05-09 [1] CRAN (R 4.5.0)
    ##  GenomicRanges        * 1.62.1   2025-12-08 [1] Bioconductor 3.22 (R 4.5.2)
    ##  ggplot2              * 4.0.1    2025-11-14 [1] CRAN (R 4.5.2)
    ##  glue                   1.8.0    2024-09-30 [1] CRAN (R 4.5.0)
    ##  gplots               * 3.3.0    2025-11-30 [1] CRAN (R 4.5.2)
    ##  gtable                 0.3.6    2024-10-25 [1] CRAN (R 4.5.0)
    ##  gtools                 3.9.5    2023-11-20 [1] CRAN (R 4.5.0)
    ##  hms                    1.1.4    2025-10-17 [1] CRAN (R 4.5.0)
    ##  htmltools              0.5.9    2025-12-04 [1] CRAN (R 4.5.2)
    ##  IRanges              * 2.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  KernSmooth             2.23-26  2025-01-01 [1] CRAN (R 4.5.2)
    ##  knitr                  1.51     2025-12-20 [1] CRAN (R 4.5.2)
    ##  labeling               0.4.3    2023-08-29 [1] CRAN (R 4.5.0)
    ##  lattice                0.22-7   2025-04-02 [1] CRAN (R 4.5.2)
    ##  lifecycle              1.0.4    2023-11-07 [1] CRAN (R 4.5.0)
    ##  limma                * 3.66.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  locfit                 1.5-9.12 2025-03-05 [1] CRAN (R 4.5.0)
    ##  magrittr               2.0.4    2025-09-12 [1] CRAN (R 4.5.0)
    ##  Matrix                 1.7-4    2025-08-28 [1] CRAN (R 4.5.0)
    ##  MatrixGenerics       * 1.22.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  matrixStats          * 1.5.0    2025-01-07 [1] CRAN (R 4.5.0)
    ##  memoise                2.0.1    2021-11-26 [1] CRAN (R 4.5.0)
    ##  otel                   0.2.0    2025-08-29 [1] CRAN (R 4.5.0)
    ##  pheatmap             * 1.0.13   2025-06-05 [1] CRAN (R 4.5.0)
    ##  pillar                 1.11.1   2025-09-17 [1] CRAN (R 4.5.0)
    ##  pkgbuild               1.4.8    2025-05-26 [1] CRAN (R 4.5.0)
    ##  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.5.0)
    ##  pkgload                1.4.1    2025-09-23 [1] CRAN (R 4.5.0)
    ##  purrr                  1.2.0    2025-11-04 [1] CRAN (R 4.5.0)
    ##  R6                     2.6.1    2025-02-15 [1] CRAN (R 4.5.0)
    ##  RColorBrewer         * 1.1-3    2022-04-03 [1] CRAN (R 4.5.0)
    ##  Rcpp                   1.1.0    2025-07-02 [1] CRAN (R 4.5.0)
    ##  readr                * 2.1.6    2025-11-14 [1] CRAN (R 4.5.2)
    ##  remotes                2.5.0    2024-03-17 [1] CRAN (R 4.5.0)
    ##  rlang                  1.1.6    2025-04-11 [1] CRAN (R 4.5.0)
    ##  rmarkdown              2.30     2025-09-28 [1] CRAN (R 4.5.0)
    ##  rstudioapi             0.17.1   2024-10-22 [1] CRAN (R 4.5.0)
    ##  S4Arrays               1.10.1   2025-12-01 [1] Bioconductor 3.22 (R 4.5.0)
    ##  S4Vectors            * 0.48.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  S7                     0.2.1    2025-11-14 [1] CRAN (R 4.5.2)
    ##  scales                 1.4.0    2025-04-24 [1] CRAN (R 4.5.0)
    ##  Seqinfo              * 1.0.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  sessioninfo            1.2.3    2025-02-05 [1] CRAN (R 4.5.0)
    ##  SparseArray            1.10.8   2025-12-18 [1] Bioconductor 3.22 (R 4.5.2)
    ##  statmod                1.5.1    2025-10-09 [1] CRAN (R 4.5.0)
    ##  stringi                1.8.7    2025-03-27 [1] CRAN (R 4.5.0)
    ##  stringr              * 1.6.0    2025-11-04 [1] CRAN (R 4.5.0)
    ##  SummarizedExperiment * 1.40.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  tibble                 3.3.0    2025-06-08 [1] CRAN (R 4.5.0)
    ##  tidyr                * 1.3.2    2025-12-19 [1] CRAN (R 4.5.2)
    ##  tidyselect             1.2.1    2024-03-11 [1] CRAN (R 4.5.0)
    ##  tzdb                   0.5.0    2025-03-15 [1] CRAN (R 4.5.0)
    ##  usethis                3.2.1    2025-09-06 [1] CRAN (R 4.5.0)
    ##  vctrs                  0.6.5    2023-12-01 [1] CRAN (R 4.5.0)
    ##  vroom                  1.6.7    2025-11-28 [1] CRAN (R 4.5.2)
    ##  withr                  3.0.2    2024-10-28 [1] CRAN (R 4.5.0)
    ##  xfun                   0.55     2025-12-16 [1] CRAN (R 4.5.2)
    ##  XVector                0.50.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  yaml                   2.3.12   2025-12-10 [1] CRAN (R 4.5.2)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
