XR261_XR252_TF_analyses
================
Alexander Lercher
2026-04-29

## LOAD PACKAGES

``` r
library(dorothea)
library(viper)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

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

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomicRanges

    ## Loading required package: Seqinfo

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

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

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     rowMedians

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
    ## ✖ lubridate::%within%() masks IRanges::%within%()
    ## ✖ dplyr::collapse()     masks IRanges::collapse()
    ## ✖ dplyr::combine()      masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::count()        masks matrixStats::count()
    ## ✖ dplyr::desc()         masks IRanges::desc()
    ## ✖ tidyr::expand()       masks S4Vectors::expand()
    ## ✖ dplyr::filter()       masks stats::filter()
    ## ✖ dplyr::first()        masks S4Vectors::first()
    ## ✖ dplyr::lag()          masks stats::lag()
    ## ✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()       masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ dplyr::rename()       masks S4Vectors::rename()
    ## ✖ lubridate::second()   masks S4Vectors::second()
    ## ✖ lubridate::second<-() masks S4Vectors::second<-()
    ## ✖ dplyr::slice()        masks IRanges::slice()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(limma)
```

    ## 
    ## Attaching package: 'limma'
    ## 
    ## The following object is masked from 'package:DESeq2':
    ## 
    ##     plotMA
    ## 
    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(grid)
```

## XR261 DATA - BMEC: polyIC vs ctrl

``` r
#--------------------------------------------------------------------
# LOAD CLUSTERING RESULTS FILTER for CLUSTER OF INTEREST
#--------------------------------------------------------------------
clusterOI = 1
data <- readr::read_tsv("output/XR261_ISG_heatmap_individual.tsv")
```

    ## Rows: 1766 Columns: 17
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): gene_id
    ## dbl (16): BMEC_PBS_1, BMEC_PBS_2, BMEC_PBS_3, BMEC_PIC_1, BMEC_PIC_2, BMEC_P...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
genesOI <- data %>%
  dplyr::filter(cluster == clusterOI) %>%
  dplyr::select(gene_id)

#--------------------------------------------------------------------
# DESeq2 OBJECT IMPORT, FILTER SAMPLES and get EXPRESSION MATRIX of GENES OI
#--------------------------------------------------------------------
# load DESeq2 object
dds <- readRDS("../XR261_DESeq2/output/XR261_dds_processed.rds")


# extract sample metadata
coldata <- as.data.frame(colData(dds))
coldata$condition <- factor(coldata$condition)

# define samples of interest
samplesOI <- coldata %>%
  dplyr::filter(condition %in% c("BMEC")) %>%
  rownames()

# subset DESeq2 object
ddsOI <- dds[, samplesOI]

# extract sample metadata of interest
coldataOI <- as.data.frame(colData(ddsOI))
coldataOI$condition <- factor(coldataOI$condition)

# extract expression matrix
rld <- rlog(ddsOI, blind = FALSE)
expr_matrix <- assay(rld)

# filter to genes of interest
expr_matrixOI <- expr_matrix[
  rownames(expr_matrix) %in% genesOI$gene_id,
]

#--------------------------------------------------------------------
# TF ENRICHMENT ANALYSES
#--------------------------------------------------------------------
# Get high confidence regulons
data("dorothea_mm")
regulon <- dorothea_mm %>%
  dplyr::filter(confidence %in% c("A","B","C"))

# run VIPER
tf_activities <- run_viper(
  mat = expr_matrixOI,
  network = regulon,
  .source = "tf",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

# write to file
write_tsv(tf_activities, "output/TF_analyses/XR261_TF_activity_enrichment_per_sample.tsv")

#--------------------------------------------------------------------
# LIMMA STATISTICAL ANALYSES ACROSS CONDITIONS
#--------------------------------------------------------------------
# convert to wide results table
tf_matrix <- tf_activities %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source")

# align columns of matrix with coldata (conditions)
tf_matrix <- tf_matrix[, rownames(coldataOI)]

# design matrix
design <- model.matrix(~ treatment, data = coldataOI)

# run limma
fit <- lmFit(tf_matrix, design)
fit <- eBayes(fit)
tf_results <- topTable(fit, number = Inf)
```

    ## Removing intercept from test coefficients

``` r
# write to file
to_file <- tf_results %>%
  rownames_to_column("TF")
write_tsv(to_file, "output/TF_analyses/XR261_TF_activity_differential_analyses.tsv")

#--------------------------------------------------------------------
# PLOT RESULTS AS VOLCANO PLOT
#--------------------------------------------------------------------
# define cutoffs
logFC_cutoff <- 1
pval_cutoff <- 0.05

# add significance column
tf_results$TF <- rownames(tf_results)
tf_results$Significant <- with(tf_results, adj.P.Val < pval_cutoff & abs(logFC) > logFC_cutoff)

# Volcano plot
p <- ggplot(tf_results, aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(aes(color = Significant), size = 5, alpha = 0.5) +
      scale_color_manual(values = c("grey", "darkred")) +
      geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "grey20") +
      geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey20") +
      geom_text_repel(
        data = tf_results[tf_results$Significant, ],
        aes(label = TF),
        max.overlaps = 20
      ) +
      ylim(0,20) +
      xlim(-3.5,3.5) +
      theme_bw() +
      theme(aspect.ratio = 1) +
      labs(
        title = "Differential TF activity BMEC PIC vs ctrl",
        x = "log2 Fold Change (TF activity)",
        y = "-log10(adj.P.Val)",
        color = "Significant"
      )
print(p)
```

![](XR261_XR252_TF_analyses_files/figure-gfm/XR261%20data%20-%20BMEC:%20polyIC%20vs%20ctrl-1.png)<!-- -->

``` r
#--------------------------------------------------------------------
# PLOT RESULTS AS HEATMAP with INDIVIDUAL VALUES
#--------------------------------------------------------------------
# filter significant TFs
sig_TFs <- tf_results %>% 
  dplyr::filter(adj.P.Val < pval_cutoff) %>%
  pull(TF)

# subset wide TF matrix
tf_sig_matrix <- tf_matrix[sig_TFs, ]

# row-wise z-score normalization
tf_z <- t(scale(t(tf_sig_matrix)))

# plot heatmap
p <- pheatmap::pheatmap(tf_z, 
             cluster_rows = TRUE, 
             cluster_cols = TRUE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             annotation_col = coldata["treatment"],
             main = "Differential TF activity BMEC PIC vs ctrl")
```

![](XR261_XR252_TF_analyses_files/figure-gfm/XR261%20data%20-%20BMEC:%20polyIC%20vs%20ctrl-2.png)<!-- -->

``` r
grid::grid.newpage()
grid::grid.draw(p$gtable)

rm(list = ls())
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
    ##  bcellViper             1.46.0   2025-10-30 [1] Bioconductor 3.22 (R 4.5.0)
    ##  Biobase              * 2.70.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocGenerics         * 0.56.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocParallel           1.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  bit                    4.6.0    2025-03-06 [1] CRAN (R 4.5.0)
    ##  bit64                  4.6.0-1  2025-01-16 [1] CRAN (R 4.5.0)
    ##  cachem                 1.1.0    2024-05-16 [1] CRAN (R 4.5.0)
    ##  class                  7.3-23   2025-01-01 [1] CRAN (R 4.5.2)
    ##  cli                    3.6.5    2025-04-23 [1] CRAN (R 4.5.0)
    ##  codetools              0.2-20   2024-03-31 [1] CRAN (R 4.5.2)
    ##  crayon                 1.5.3    2024-06-20 [1] CRAN (R 4.5.0)
    ##  data.table             1.18.0   2025-12-24 [1] CRAN (R 4.5.2)
    ##  decoupleR            * 2.16.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.0)
    ##  DelayedArray           0.36.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  DESeq2               * 1.50.2   2025-11-17 [1] Bioconductor 3.22 (R 4.5.2)
    ##  devtools               2.4.6    2025-10-03 [1] CRAN (R 4.5.0)
    ##  digest                 0.6.39   2025-11-19 [1] CRAN (R 4.5.2)
    ##  dorothea             * 1.22.0   2025-11-04 [1] Bioconductor 3.22 (R 4.5.0)
    ##  dplyr                * 1.1.4    2023-11-17 [1] CRAN (R 4.5.0)
    ##  e1071                  1.7-17   2025-12-18 [1] CRAN (R 4.5.2)
    ##  ellipsis               0.3.2    2021-04-29 [1] CRAN (R 4.5.0)
    ##  evaluate               1.0.5    2025-08-27 [1] CRAN (R 4.5.0)
    ##  farver                 2.1.2    2024-05-13 [1] CRAN (R 4.5.0)
    ##  fastmap                1.2.0    2024-05-15 [1] CRAN (R 4.5.0)
    ##  forcats              * 1.0.1    2025-09-25 [1] CRAN (R 4.5.0)
    ##  fs                     1.6.6    2025-04-12 [1] CRAN (R 4.5.0)
    ##  generics             * 0.1.4    2025-05-09 [1] CRAN (R 4.5.0)
    ##  GenomicRanges        * 1.62.1   2025-12-08 [1] Bioconductor 3.22 (R 4.5.2)
    ##  ggplot2              * 4.0.1    2025-11-14 [1] CRAN (R 4.5.2)
    ##  ggrepel              * 0.9.6    2024-09-07 [1] CRAN (R 4.5.0)
    ##  glue                   1.8.0    2024-09-30 [1] CRAN (R 4.5.0)
    ##  gtable                 0.3.6    2024-10-25 [1] CRAN (R 4.5.0)
    ##  hms                    1.1.4    2025-10-17 [1] CRAN (R 4.5.0)
    ##  htmltools              0.5.9    2025-12-04 [1] CRAN (R 4.5.2)
    ##  htmlwidgets            1.6.4    2023-12-06 [1] CRAN (R 4.5.0)
    ##  httr                   1.4.7    2023-08-15 [1] CRAN (R 4.5.0)
    ##  IRanges              * 2.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  jsonlite               2.0.0    2025-03-27 [1] CRAN (R 4.5.0)
    ##  kernlab                0.9-33   2024-08-13 [1] CRAN (R 4.5.0)
    ##  KernSmooth             2.23-26  2025-01-01 [1] CRAN (R 4.5.2)
    ##  knitr                  1.51     2025-12-20 [1] CRAN (R 4.5.2)
    ##  labeling               0.4.3    2023-08-29 [1] CRAN (R 4.5.0)
    ##  lattice                0.22-7   2025-04-02 [1] CRAN (R 4.5.2)
    ##  lazyeval               0.2.2    2019-03-15 [1] CRAN (R 4.5.0)
    ##  lifecycle              1.0.4    2023-11-07 [1] CRAN (R 4.5.0)
    ##  limma                * 3.66.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  locfit                 1.5-9.12 2025-03-05 [1] CRAN (R 4.5.0)
    ##  lubridate            * 1.9.4    2024-12-08 [1] CRAN (R 4.5.0)
    ##  magrittr               2.0.4    2025-09-12 [1] CRAN (R 4.5.0)
    ##  MASS                   7.3-65   2025-02-28 [1] CRAN (R 4.5.2)
    ##  Matrix                 1.7-4    2025-08-28 [1] CRAN (R 4.5.0)
    ##  MatrixGenerics       * 1.22.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  matrixStats          * 1.5.0    2025-01-07 [1] CRAN (R 4.5.0)
    ##  memoise                2.0.1    2021-11-26 [1] CRAN (R 4.5.0)
    ##  mixtools               2.0.0.1  2025-03-08 [1] CRAN (R 4.5.0)
    ##  nlme                   3.1-168  2025-03-31 [1] CRAN (R 4.5.2)
    ##  otel                   0.2.0    2025-08-29 [1] CRAN (R 4.5.0)
    ##  parallelly             1.46.0   2025-12-12 [1] CRAN (R 4.5.2)
    ##  pheatmap             * 1.0.13   2025-06-05 [1] CRAN (R 4.5.0)
    ##  pillar                 1.11.1   2025-09-17 [1] CRAN (R 4.5.0)
    ##  pkgbuild               1.4.8    2025-05-26 [1] CRAN (R 4.5.0)
    ##  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.5.0)
    ##  pkgload                1.4.1    2025-09-23 [1] CRAN (R 4.5.0)
    ##  plotly                 4.11.0   2025-06-19 [1] CRAN (R 4.5.0)
    ##  proxy                  0.4-29   2025-12-29 [1] CRAN (R 4.5.2)
    ##  purrr                * 1.2.0    2025-11-04 [1] CRAN (R 4.5.0)
    ##  R6                     2.6.1    2025-02-15 [1] CRAN (R 4.5.0)
    ##  RColorBrewer           1.1-3    2022-04-03 [1] CRAN (R 4.5.0)
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
    ##  segmented              2.1-4    2025-02-28 [1] CRAN (R 4.5.0)
    ##  Seqinfo              * 1.0.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  sessioninfo            1.2.3    2025-02-05 [1] CRAN (R 4.5.0)
    ##  SparseArray            1.10.8   2025-12-18 [1] Bioconductor 3.22 (R 4.5.2)
    ##  statmod                1.5.1    2025-10-09 [1] CRAN (R 4.5.0)
    ##  stringi                1.8.7    2025-03-27 [1] CRAN (R 4.5.0)
    ##  stringr              * 1.6.0    2025-11-04 [1] CRAN (R 4.5.0)
    ##  SummarizedExperiment * 1.40.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  survival               3.8-3    2024-12-17 [1] CRAN (R 4.5.2)
    ##  tibble               * 3.3.0    2025-06-08 [1] CRAN (R 4.5.0)
    ##  tidyr                * 1.3.2    2025-12-19 [1] CRAN (R 4.5.2)
    ##  tidyselect             1.2.1    2024-03-11 [1] CRAN (R 4.5.0)
    ##  tidyverse            * 2.0.0    2023-02-22 [1] CRAN (R 4.5.0)
    ##  timechange             0.3.0    2024-01-18 [1] CRAN (R 4.5.0)
    ##  tzdb                   0.5.0    2025-03-15 [1] CRAN (R 4.5.0)
    ##  usethis                3.2.1    2025-09-06 [1] CRAN (R 4.5.0)
    ##  vctrs                  0.6.5    2023-12-01 [1] CRAN (R 4.5.0)
    ##  viper                * 1.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  viridisLite            0.4.2    2023-05-02 [1] CRAN (R 4.5.0)
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
