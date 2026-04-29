XR251_TF_enrichment
================
Alexander Lercher
2026-04-29

## LOAD PACKAGES

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.6
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ ggplot2   4.0.1     ✔ tibble    3.3.0
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.2
    ## ✔ purrr     1.2.0     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(dplyr)
library(ggplot2)
library(dorothea)
library(viper)
```

    ## Loading required package: Biobase
    ## Loading required package: BiocGenerics
    ## Loading required package: generics
    ## 
    ## Attaching package: 'generics'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     as.difftime
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     explain
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union
    ## 
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min
    ## 
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

``` r
library(decoupleR)
```

    ## 
    ## Attaching package: 'decoupleR'
    ## 
    ## The following object is masked from 'package:dorothea':
    ## 
    ##     run_viper

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors
    ## Loading required package: stats4
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand
    ## 
    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname
    ## 
    ## Loading required package: IRanges
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce
    ## 
    ## Loading required package: GenomicRanges
    ## Loading required package: Seqinfo
    ## Loading required package: SummarizedExperiment
    ## Loading required package: MatrixGenerics
    ## Loading required package: matrixStats
    ## 
    ## Attaching package: 'matrixStats'
    ## 
    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count
    ## 
    ## 
    ## Attaching package: 'MatrixGenerics'
    ## 
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
    ## 
    ## The following object is masked from 'package:Biobase':
    ## 
    ##     rowMedians

``` r
library(pheatmap)
library(Seurat)
```

    ## Loading required package: SeuratObject
    ## Loading required package: sp
    ## 
    ## Attaching package: 'sp'
    ## 
    ## The following object is masked from 'package:IRanges':
    ## 
    ##     %over%
    ## 
    ## 
    ## Attaching package: 'SeuratObject'
    ## 
    ## The following object is masked from 'package:SummarizedExperiment':
    ## 
    ##     Assays
    ## 
    ## The following object is masked from 'package:GenomicRanges':
    ## 
    ##     intersect
    ## 
    ## The following object is masked from 'package:Seqinfo':
    ## 
    ##     intersect
    ## 
    ## The following object is masked from 'package:IRanges':
    ## 
    ##     intersect
    ## 
    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     intersect
    ## 
    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     intersect
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, t
    ## 
    ## 
    ## Attaching package: 'Seurat'
    ## 
    ## The following object is masked from 'package:SummarizedExperiment':
    ## 
    ##     Assays

## LOAD in SEURAT RDS file modified by me

``` r
data <- readRDS(file = "input/seu_obj_v12192025_AL.rds")
```

## EXTRACT PSEUDOBULK DATA

``` r
pseudobulk <- AggregateExpression(
  data,
  group.by = c("cell_type","condition"),
  assays = "RNA",
  slot = "counts"
)$RNA
```

    ## Names of identity class contain underscores ('_'), replacing with dashes ('-')
    ## This message is displayed once every 8 hours.

## CREATE METADATA

``` r
# get list of ISG names and save to file
meta <- data.frame(
  sample = colnames(pseudobulk)
)
meta$celltype <- sub("_.*", "", meta$sample)
meta$sample_id <- sub(".*_", "", meta$sample)

meta$condition <- ifelse(meta$sample_id == c("00-hours"),
                         "control", "treated")
rownames(meta) <- meta$sample

# safety net
stopifnot(all(colnames(pseudobulk) == rownames(meta)))
stopifnot(is.data.frame(meta))
stopifnot("condition" %in% colnames(meta))
```

## RUN DESEQ2

``` r
dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk,
  colData = meta,
  design = ~condition
)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
vsd <- vst(dds, blind = FALSE)
expr_matrix <- assay(vsd)
```

## RUN TF ACTIVITY ENRICHMENT

``` r
# get TF networks
data("dorothea_mm")
regulon <- dorothea_mm %>%
  filter(confidence %in% c("A","B","C"))

# run viper
tf_activities <- run_viper(
  mat = expr_matrix,
  network = regulon,
  .source = "tf",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

# change dashes to underscores
tf_activities <- tf_activities %>%
  dplyr::mutate(
    condition = gsub("-","_", condition),
    condition = gsub(" ","_", condition)
  )

# make wide table and write to file
tf_activities_wide <- tf_activities %>%
  dplyr::select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  rownames_to_column(var = "TF") %>%
  write_tsv("output/TF_activities/XR251_TF_activity_enrichment_per_sample.tsv")

# add cell type and time info to activities table
tf_activities$celltype <- sub("_(\\d+)_hours$", "", tf_activities$condition)
tf_activities$time <- stringr::str_extract(tf_activities$condition, "\\d+_hours")

# extract all cell type names and time points
all_celltypes <- unique(tf_activities$celltype)
all_timepoints <- unique(tf_activities$time)

# subset by cell type of interest
for(i in 1:length(all_celltypes)) {
tf_activities_OI <- tf_activities%>%
  dplyr::filter(celltype == all_celltypes[[i]])

# make wide table
tf_matrix_OI <- tf_activities_OI %>%
  dplyr::select(source, time, score) %>%
  pivot_wider(names_from = time, values_from = score) %>%
  column_to_rownames("source")

# identify most variable TFs
tf_sd <- apply(tf_matrix_OI, 1, sd, na.rm = TRUE)

# get most variable TFs
number_of_tfs <- 50
top_tfs <- names(sort(tf_sd, decreasing = TRUE))[1:number_of_tfs]

# subset the matrix
tf_matrix_top <- tf_matrix_OI[top_tfs, ]

# plot heatmap and save to file
my_palette <- colorRampPalette(c("#6BAED6", "#FDFDFD", "#FB6A4A"))(50)
p <- pheatmap::pheatmap(
      tf_matrix_top,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      scale = "row",
      color = my_palette,
      height = 12,
      width = 8,
      main = paste0("Top ",as.character(number_of_tfs)," variable TFs in ",gsub("_"," ",all_celltypes[i])," (row z norm)"))
grid::grid.newpage()
grid::grid.draw(p$gtable)
}
```

![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-1.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-2.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-3.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-4.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-5.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-6.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-7.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-8.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-9.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-10.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-11.png)<!-- -->![](XR251_TF_enrichment_files/figure-gfm/run%20TF%20activity%20enrichment-12.png)<!-- -->

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
    ##  bcellViper             1.46.0   2025-10-30 [1] Bioconductor 3.22 (R 4.5.0)
    ##  Biobase              * 2.70.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocGenerics         * 0.56.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocParallel           1.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  bit                    4.6.0    2025-03-06 [1] CRAN (R 4.5.0)
    ##  bit64                  4.6.0-1  2025-01-16 [1] CRAN (R 4.5.0)
    ##  cachem                 1.1.0    2024-05-16 [1] CRAN (R 4.5.0)
    ##  class                  7.3-23   2025-01-01 [1] CRAN (R 4.5.2)
    ##  cli                    3.6.5    2025-04-23 [1] CRAN (R 4.5.0)
    ##  cluster                2.1.8.1  2025-03-12 [1] CRAN (R 4.5.2)
    ##  codetools              0.2-20   2024-03-31 [1] CRAN (R 4.5.2)
    ##  cowplot                1.2.0    2025-07-07 [1] CRAN (R 4.5.0)
    ##  crayon                 1.5.3    2024-06-20 [1] CRAN (R 4.5.0)
    ##  data.table             1.18.0   2025-12-24 [1] CRAN (R 4.5.2)
    ##  decoupleR            * 2.16.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.0)
    ##  DelayedArray           0.36.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  deldir                 2.0-4    2024-02-28 [1] CRAN (R 4.5.0)
    ##  DESeq2               * 1.50.2   2025-11-17 [1] Bioconductor 3.22 (R 4.5.2)
    ##  devtools               2.4.6    2025-10-03 [1] CRAN (R 4.5.0)
    ##  digest                 0.6.39   2025-11-19 [1] CRAN (R 4.5.2)
    ##  dorothea             * 1.22.0   2025-11-04 [1] Bioconductor 3.22 (R 4.5.0)
    ##  dotCall64              1.2      2024-10-04 [1] CRAN (R 4.5.0)
    ##  dplyr                * 1.1.4    2023-11-17 [1] CRAN (R 4.5.0)
    ##  e1071                  1.7-17   2025-12-18 [1] CRAN (R 4.5.2)
    ##  ellipsis               0.3.2    2021-04-29 [1] CRAN (R 4.5.0)
    ##  evaluate               1.0.5    2025-08-27 [1] CRAN (R 4.5.0)
    ##  farver                 2.1.2    2024-05-13 [1] CRAN (R 4.5.0)
    ##  fastDummies            1.7.5    2025-01-20 [1] CRAN (R 4.5.0)
    ##  fastmap                1.2.0    2024-05-15 [1] CRAN (R 4.5.0)
    ##  fitdistrplus           1.2-4    2025-07-03 [1] CRAN (R 4.5.0)
    ##  forcats              * 1.0.1    2025-09-25 [1] CRAN (R 4.5.0)
    ##  fs                     1.6.6    2025-04-12 [1] CRAN (R 4.5.0)
    ##  future                 1.68.0   2025-11-17 [1] CRAN (R 4.5.2)
    ##  future.apply           1.20.1   2025-12-09 [1] CRAN (R 4.5.0)
    ##  generics             * 0.1.4    2025-05-09 [1] CRAN (R 4.5.0)
    ##  GenomicRanges        * 1.62.1   2025-12-08 [1] Bioconductor 3.22 (R 4.5.2)
    ##  ggplot2              * 4.0.1    2025-11-14 [1] CRAN (R 4.5.2)
    ##  ggrepel                0.9.6    2024-09-07 [1] CRAN (R 4.5.0)
    ##  ggridges               0.5.7    2025-08-27 [1] CRAN (R 4.5.0)
    ##  globals                0.18.0   2025-05-08 [1] CRAN (R 4.5.0)
    ##  glue                   1.8.0    2024-09-30 [1] CRAN (R 4.5.0)
    ##  goftest                1.2-3    2021-10-07 [1] CRAN (R 4.5.0)
    ##  gridExtra              2.3      2017-09-09 [1] CRAN (R 4.5.0)
    ##  gtable                 0.3.6    2024-10-25 [1] CRAN (R 4.5.0)
    ##  hms                    1.1.4    2025-10-17 [1] CRAN (R 4.5.0)
    ##  htmltools              0.5.9    2025-12-04 [1] CRAN (R 4.5.2)
    ##  htmlwidgets            1.6.4    2023-12-06 [1] CRAN (R 4.5.0)
    ##  httpuv                 1.6.16   2025-04-16 [1] CRAN (R 4.5.0)
    ##  httr                   1.4.7    2023-08-15 [1] CRAN (R 4.5.0)
    ##  ica                    1.0-3    2022-07-08 [1] CRAN (R 4.5.0)
    ##  igraph                 2.2.1    2025-10-27 [1] CRAN (R 4.5.0)
    ##  IRanges              * 2.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  irlba                  2.3.5.1  2022-10-03 [1] CRAN (R 4.5.0)
    ##  jsonlite               2.0.0    2025-03-27 [1] CRAN (R 4.5.0)
    ##  kernlab                0.9-33   2024-08-13 [1] CRAN (R 4.5.0)
    ##  KernSmooth             2.23-26  2025-01-01 [1] CRAN (R 4.5.2)
    ##  knitr                  1.51     2025-12-20 [1] CRAN (R 4.5.2)
    ##  later                  1.4.4    2025-08-27 [1] CRAN (R 4.5.0)
    ##  lattice                0.22-7   2025-04-02 [1] CRAN (R 4.5.2)
    ##  lazyeval               0.2.2    2019-03-15 [1] CRAN (R 4.5.0)
    ##  lifecycle              1.0.4    2023-11-07 [1] CRAN (R 4.5.0)
    ##  listenv                0.10.0   2025-11-02 [1] CRAN (R 4.5.0)
    ##  lmtest                 0.9-40   2022-03-21 [1] CRAN (R 4.5.0)
    ##  locfit                 1.5-9.12 2025-03-05 [1] CRAN (R 4.5.0)
    ##  lubridate            * 1.9.4    2024-12-08 [1] CRAN (R 4.5.0)
    ##  magrittr               2.0.4    2025-09-12 [1] CRAN (R 4.5.0)
    ##  MASS                   7.3-65   2025-02-28 [1] CRAN (R 4.5.2)
    ##  Matrix                 1.7-4    2025-08-28 [1] CRAN (R 4.5.0)
    ##  MatrixGenerics       * 1.22.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  matrixStats          * 1.5.0    2025-01-07 [1] CRAN (R 4.5.0)
    ##  memoise                2.0.1    2021-11-26 [1] CRAN (R 4.5.0)
    ##  mime                   0.13     2025-03-17 [1] CRAN (R 4.5.0)
    ##  miniUI                 0.1.2    2025-04-17 [1] CRAN (R 4.5.0)
    ##  mixtools               2.0.0.1  2025-03-08 [1] CRAN (R 4.5.0)
    ##  nlme                   3.1-168  2025-03-31 [1] CRAN (R 4.5.2)
    ##  otel                   0.2.0    2025-08-29 [1] CRAN (R 4.5.0)
    ##  parallelly             1.46.0   2025-12-12 [1] CRAN (R 4.5.2)
    ##  patchwork              1.3.2    2025-08-25 [1] CRAN (R 4.5.0)
    ##  pbapply                1.7-4    2025-07-20 [1] CRAN (R 4.5.0)
    ##  pheatmap             * 1.0.13   2025-06-05 [1] CRAN (R 4.5.0)
    ##  pillar                 1.11.1   2025-09-17 [1] CRAN (R 4.5.0)
    ##  pkgbuild               1.4.8    2025-05-26 [1] CRAN (R 4.5.0)
    ##  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.5.0)
    ##  pkgload                1.4.1    2025-09-23 [1] CRAN (R 4.5.0)
    ##  plotly                 4.11.0   2025-06-19 [1] CRAN (R 4.5.0)
    ##  plyr                   1.8.9    2023-10-02 [1] CRAN (R 4.5.0)
    ##  png                    0.1-8    2022-11-29 [1] CRAN (R 4.5.0)
    ##  polyclip               1.10-7   2024-07-23 [1] CRAN (R 4.5.0)
    ##  progressr              0.18.0   2025-11-06 [1] CRAN (R 4.5.0)
    ##  promises               1.5.0    2025-11-01 [1] CRAN (R 4.5.0)
    ##  proxy                  0.4-29   2025-12-29 [1] CRAN (R 4.5.2)
    ##  purrr                * 1.2.0    2025-11-04 [1] CRAN (R 4.5.0)
    ##  R6                     2.6.1    2025-02-15 [1] CRAN (R 4.5.0)
    ##  RANN                   2.6.2    2024-08-25 [1] CRAN (R 4.5.0)
    ##  RColorBrewer           1.1-3    2022-04-03 [1] CRAN (R 4.5.0)
    ##  Rcpp                   1.1.0    2025-07-02 [1] CRAN (R 4.5.0)
    ##  RcppAnnoy              0.0.22   2024-01-23 [1] CRAN (R 4.5.0)
    ##  RcppHNSW               0.6.0    2024-02-04 [1] CRAN (R 4.5.0)
    ##  readr                * 2.1.6    2025-11-14 [1] CRAN (R 4.5.2)
    ##  remotes                2.5.0    2024-03-17 [1] CRAN (R 4.5.0)
    ##  reshape2               1.4.5    2025-11-12 [1] CRAN (R 4.5.0)
    ##  reticulate             1.44.1   2025-11-14 [1] CRAN (R 4.5.2)
    ##  rlang                  1.1.6    2025-04-11 [1] CRAN (R 4.5.0)
    ##  rmarkdown              2.30     2025-09-28 [1] CRAN (R 4.5.0)
    ##  ROCR                   1.0-11   2020-05-02 [1] CRAN (R 4.5.0)
    ##  RSpectra               0.16-2   2024-07-18 [1] CRAN (R 4.5.0)
    ##  rstudioapi             0.17.1   2024-10-22 [1] CRAN (R 4.5.0)
    ##  Rtsne                  0.17     2023-12-07 [1] CRAN (R 4.5.0)
    ##  S4Arrays               1.10.1   2025-12-01 [1] Bioconductor 3.22 (R 4.5.0)
    ##  S4Vectors            * 0.48.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  S7                     0.2.1    2025-11-14 [1] CRAN (R 4.5.2)
    ##  scales                 1.4.0    2025-04-24 [1] CRAN (R 4.5.0)
    ##  scattermore            1.2      2023-06-12 [1] CRAN (R 4.5.0)
    ##  sctransform            0.4.2    2025-04-30 [1] CRAN (R 4.5.0)
    ##  segmented              2.1-4    2025-02-28 [1] CRAN (R 4.5.0)
    ##  Seqinfo              * 1.0.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  sessioninfo            1.2.3    2025-02-05 [1] CRAN (R 4.5.0)
    ##  Seurat               * 5.4.0    2025-12-14 [1] CRAN (R 4.5.2)
    ##  SeuratObject         * 5.3.0    2025-12-12 [1] CRAN (R 4.5.2)
    ##  shiny                  1.12.1   2025-12-09 [1] CRAN (R 4.5.0)
    ##  sp                   * 2.2-0    2025-02-01 [1] CRAN (R 4.5.0)
    ##  spam                   2.11-1   2025-01-20 [1] CRAN (R 4.5.0)
    ##  SparseArray            1.10.8   2025-12-18 [1] Bioconductor 3.22 (R 4.5.2)
    ##  spatstat.data          3.1-9    2025-10-18 [1] CRAN (R 4.5.0)
    ##  spatstat.explore       3.6-0    2025-11-22 [1] CRAN (R 4.5.2)
    ##  spatstat.geom          3.6-1    2025-11-20 [1] CRAN (R 4.5.2)
    ##  spatstat.random        3.4-3    2025-11-21 [1] CRAN (R 4.5.2)
    ##  spatstat.sparse        3.1-0    2024-06-21 [1] CRAN (R 4.5.0)
    ##  spatstat.univar        3.1-5    2025-11-17 [1] CRAN (R 4.5.2)
    ##  spatstat.utils         3.2-0    2025-09-20 [1] CRAN (R 4.5.0)
    ##  stringi                1.8.7    2025-03-27 [1] CRAN (R 4.5.0)
    ##  stringr              * 1.6.0    2025-11-04 [1] CRAN (R 4.5.0)
    ##  SummarizedExperiment * 1.40.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  survival               3.8-3    2024-12-17 [1] CRAN (R 4.5.2)
    ##  tensor                 1.5.1    2025-06-17 [1] CRAN (R 4.5.0)
    ##  tibble               * 3.3.0    2025-06-08 [1] CRAN (R 4.5.0)
    ##  tidyr                * 1.3.2    2025-12-19 [1] CRAN (R 4.5.2)
    ##  tidyselect             1.2.1    2024-03-11 [1] CRAN (R 4.5.0)
    ##  tidyverse            * 2.0.0    2023-02-22 [1] CRAN (R 4.5.0)
    ##  timechange             0.3.0    2024-01-18 [1] CRAN (R 4.5.0)
    ##  tzdb                   0.5.0    2025-03-15 [1] CRAN (R 4.5.0)
    ##  usethis                3.2.1    2025-09-06 [1] CRAN (R 4.5.0)
    ##  uwot                   0.2.4    2025-11-10 [1] CRAN (R 4.5.0)
    ##  vctrs                  0.6.5    2023-12-01 [1] CRAN (R 4.5.0)
    ##  viper                * 1.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  viridisLite            0.4.2    2023-05-02 [1] CRAN (R 4.5.0)
    ##  vroom                  1.6.7    2025-11-28 [1] CRAN (R 4.5.2)
    ##  withr                  3.0.2    2024-10-28 [1] CRAN (R 4.5.0)
    ##  xfun                   0.55     2025-12-16 [1] CRAN (R 4.5.2)
    ##  xtable                 1.8-4    2019-04-21 [1] CRAN (R 4.5.0)
    ##  XVector                0.50.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  yaml                   2.3.12   2025-12-10 [1] CRAN (R 4.5.2)
    ##  zoo                    1.8-15   2025-12-15 [1] CRAN (R 4.5.2)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
