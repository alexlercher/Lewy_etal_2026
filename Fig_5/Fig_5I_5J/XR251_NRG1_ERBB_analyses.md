XR251_NRG1_ERBB_analyses
================
Alexander Lercher
2026-04-29

## LOAD PACKAGES

``` r
library(Seurat)
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, t

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
library(pheatmap)
library(circlize)
```

    ## ========================================
    ## circlize version 0.4.17
    ## CRAN page: https://cran.r-project.org/package=circlize
    ## Github page: https://github.com/jokergoo/circlize
    ## Documentation: https://jokergoo.github.io/circlize_book/book/
    ## 
    ## If you use it in published research, please cite:
    ## Gu, Z. circlize implements and enhances circular visualization
    ##   in R. Bioinformatics 2014.
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(circlize))
    ## ========================================

``` r
library(scales)
```

    ## 
    ## Attaching package: 'scales'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard
    ## 
    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

``` r
library(ggplot2)
library(grid)
```

## LOAD in SEURAT RDS file modified by me

``` r
data_seu <- readRDS(file = "input/seu_obj_v12192025_AL.rds")
```

## LOAD TF ENRICHMENT DATA

``` r
data_tf <- read_tsv("output/TF_activities/XR251_TF_activity_enrichment_per_sample.tsv")
```

    ## Rows: 273 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): TF
    ## dbl (36): Astrocytes_00_hours, Astrocytes_02_hours, Astrocytes_16_hours, End...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

## PLOT and SAVE HEATMAP OF TFs and CELLS of interest

``` r
# filter for neurons and TFs
neurons <- c("Excitatory_neurons","Inhibitory_neurons","Interneurons")
TFs <- c("Elk1","Stat3","Stat5a","Stat5b","Creb1")

for(i in 1:length(neurons)){
data_tf_filtered <- data_tf %>%
  dplyr::filter(TF %in% TFs) %>%
  dplyr::select(TF,
                starts_with(neurons[i])
  )

data_tf_filtered_matrix <- data_tf_filtered %>%
  column_to_rownames("TF") %>%
  as.matrix()

# plot heatmap and save to file
my_palette <- colorRampPalette(c("#759CA1", "#FDFDFD", "#DEA07E"))(50)
p <- pheatmap::pheatmap(
    data_tf_filtered_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "row",
    color = my_palette,
    height = 5,
    width = 3,
    main = paste0("ErbB downstream TFs in ",gsub("_"," ",neurons[i])),
    fontsize = 6)
grid::grid.newpage()
grid::grid.draw(p$gtable)
}
```

![](XR251_NRG1_ERBB_analyses_files/figure-gfm/plot%20and%20save%20heatmap%20of%20TFs%20and%20cells%20of%20interest-1.png)<!-- -->![](XR251_NRG1_ERBB_analyses_files/figure-gfm/plot%20and%20save%20heatmap%20of%20TFs%20and%20cells%20of%20interest-2.png)<!-- -->![](XR251_NRG1_ERBB_analyses_files/figure-gfm/plot%20and%20save%20heatmap%20of%20TFs%20and%20cells%20of%20interest-3.png)<!-- -->

## PLOT and SAVE DOTPLOT OF LIGANDS of interest in MICROGLIA

``` r
ligands <- c("Nrg1","Nrg2","Nrg3","Nrg4")

subset_data_seu <- subset(data_seu, idents = "Microglia")

p <- DotPlot(subset_data_seu, features = ligands, group.by = "condition") +
    scale_size(
      range = c(4, 20),
      breaks = c(1, 5, 10, 25, 50, 75),
      limits = c(0, 100)) +
    scale_color_gradient(low = "lightgrey", high = "firebrick") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size =12),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          aspect.ratio = 1) +
    ggtitle(paste0("Nrg ligand expression in Microglia")) +
    coord_flip()
```

    ## Warning: Scaling data with a low number of groups may produce misleading
    ## results

    ## Scale for size is already present.
    ## Adding another scale for size, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

``` r
print(p)
```

![](XR251_NRG1_ERBB_analyses_files/figure-gfm/plot%20and%20save%20dotplot%20of%20ligands%20of%20interest%20in%20microglia-1.png)<!-- -->

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
    ##  package          * version date (UTC) lib source
    ##  abind              1.4-8   2024-09-12 [1] CRAN (R 4.5.0)
    ##  bit                4.6.0   2025-03-06 [1] CRAN (R 4.5.0)
    ##  bit64              4.6.0-1 2025-01-16 [1] CRAN (R 4.5.0)
    ##  cachem             1.1.0   2024-05-16 [1] CRAN (R 4.5.0)
    ##  circlize         * 0.4.17  2025-12-08 [1] CRAN (R 4.5.2)
    ##  cli                3.6.5   2025-04-23 [1] CRAN (R 4.5.0)
    ##  cluster            2.1.8.1 2025-03-12 [1] CRAN (R 4.5.2)
    ##  codetools          0.2-20  2024-03-31 [1] CRAN (R 4.5.2)
    ##  colorspace         2.1-2   2025-09-22 [1] CRAN (R 4.5.0)
    ##  cowplot            1.2.0   2025-07-07 [1] CRAN (R 4.5.0)
    ##  crayon             1.5.3   2024-06-20 [1] CRAN (R 4.5.0)
    ##  data.table         1.18.0  2025-12-24 [1] CRAN (R 4.5.2)
    ##  deldir             2.0-4   2024-02-28 [1] CRAN (R 4.5.0)
    ##  devtools           2.4.6   2025-10-03 [1] CRAN (R 4.5.0)
    ##  digest             0.6.39  2025-11-19 [1] CRAN (R 4.5.2)
    ##  dotCall64          1.2     2024-10-04 [1] CRAN (R 4.5.0)
    ##  dplyr            * 1.1.4   2023-11-17 [1] CRAN (R 4.5.0)
    ##  ellipsis           0.3.2   2021-04-29 [1] CRAN (R 4.5.0)
    ##  evaluate           1.0.5   2025-08-27 [1] CRAN (R 4.5.0)
    ##  farver             2.1.2   2024-05-13 [1] CRAN (R 4.5.0)
    ##  fastDummies        1.7.5   2025-01-20 [1] CRAN (R 4.5.0)
    ##  fastmap            1.2.0   2024-05-15 [1] CRAN (R 4.5.0)
    ##  fitdistrplus       1.2-4   2025-07-03 [1] CRAN (R 4.5.0)
    ##  forcats          * 1.0.1   2025-09-25 [1] CRAN (R 4.5.0)
    ##  fs                 1.6.6   2025-04-12 [1] CRAN (R 4.5.0)
    ##  future             1.68.0  2025-11-17 [1] CRAN (R 4.5.2)
    ##  future.apply       1.20.1  2025-12-09 [1] CRAN (R 4.5.0)
    ##  generics           0.1.4   2025-05-09 [1] CRAN (R 4.5.0)
    ##  ggplot2          * 4.0.1   2025-11-14 [1] CRAN (R 4.5.2)
    ##  ggrepel            0.9.6   2024-09-07 [1] CRAN (R 4.5.0)
    ##  ggridges           0.5.7   2025-08-27 [1] CRAN (R 4.5.0)
    ##  GlobalOptions      0.1.3   2025-11-28 [1] CRAN (R 4.5.2)
    ##  globals            0.18.0  2025-05-08 [1] CRAN (R 4.5.0)
    ##  glue               1.8.0   2024-09-30 [1] CRAN (R 4.5.0)
    ##  goftest            1.2-3   2021-10-07 [1] CRAN (R 4.5.0)
    ##  gridExtra          2.3     2017-09-09 [1] CRAN (R 4.5.0)
    ##  gtable             0.3.6   2024-10-25 [1] CRAN (R 4.5.0)
    ##  hms                1.1.4   2025-10-17 [1] CRAN (R 4.5.0)
    ##  htmltools          0.5.9   2025-12-04 [1] CRAN (R 4.5.2)
    ##  htmlwidgets        1.6.4   2023-12-06 [1] CRAN (R 4.5.0)
    ##  httpuv             1.6.16  2025-04-16 [1] CRAN (R 4.5.0)
    ##  httr               1.4.7   2023-08-15 [1] CRAN (R 4.5.0)
    ##  ica                1.0-3   2022-07-08 [1] CRAN (R 4.5.0)
    ##  igraph             2.2.1   2025-10-27 [1] CRAN (R 4.5.0)
    ##  irlba              2.3.5.1 2022-10-03 [1] CRAN (R 4.5.0)
    ##  jsonlite           2.0.0   2025-03-27 [1] CRAN (R 4.5.0)
    ##  KernSmooth         2.23-26 2025-01-01 [1] CRAN (R 4.5.2)
    ##  knitr              1.51    2025-12-20 [1] CRAN (R 4.5.2)
    ##  labeling           0.4.3   2023-08-29 [1] CRAN (R 4.5.0)
    ##  later              1.4.4   2025-08-27 [1] CRAN (R 4.5.0)
    ##  lattice            0.22-7  2025-04-02 [1] CRAN (R 4.5.2)
    ##  lazyeval           0.2.2   2019-03-15 [1] CRAN (R 4.5.0)
    ##  lifecycle          1.0.4   2023-11-07 [1] CRAN (R 4.5.0)
    ##  listenv            0.10.0  2025-11-02 [1] CRAN (R 4.5.0)
    ##  lmtest             0.9-40  2022-03-21 [1] CRAN (R 4.5.0)
    ##  lubridate        * 1.9.4   2024-12-08 [1] CRAN (R 4.5.0)
    ##  magrittr           2.0.4   2025-09-12 [1] CRAN (R 4.5.0)
    ##  MASS               7.3-65  2025-02-28 [1] CRAN (R 4.5.2)
    ##  Matrix             1.7-4   2025-08-28 [1] CRAN (R 4.5.0)
    ##  matrixStats        1.5.0   2025-01-07 [1] CRAN (R 4.5.0)
    ##  memoise            2.0.1   2021-11-26 [1] CRAN (R 4.5.0)
    ##  mime               0.13    2025-03-17 [1] CRAN (R 4.5.0)
    ##  miniUI             0.1.2   2025-04-17 [1] CRAN (R 4.5.0)
    ##  nlme               3.1-168 2025-03-31 [1] CRAN (R 4.5.2)
    ##  otel               0.2.0   2025-08-29 [1] CRAN (R 4.5.0)
    ##  parallelly         1.46.0  2025-12-12 [1] CRAN (R 4.5.2)
    ##  patchwork          1.3.2   2025-08-25 [1] CRAN (R 4.5.0)
    ##  pbapply            1.7-4   2025-07-20 [1] CRAN (R 4.5.0)
    ##  pheatmap         * 1.0.13  2025-06-05 [1] CRAN (R 4.5.0)
    ##  pillar             1.11.1  2025-09-17 [1] CRAN (R 4.5.0)
    ##  pkgbuild           1.4.8   2025-05-26 [1] CRAN (R 4.5.0)
    ##  pkgconfig          2.0.3   2019-09-22 [1] CRAN (R 4.5.0)
    ##  pkgload            1.4.1   2025-09-23 [1] CRAN (R 4.5.0)
    ##  plotly             4.11.0  2025-06-19 [1] CRAN (R 4.5.0)
    ##  plyr               1.8.9   2023-10-02 [1] CRAN (R 4.5.0)
    ##  png                0.1-8   2022-11-29 [1] CRAN (R 4.5.0)
    ##  polyclip           1.10-7  2024-07-23 [1] CRAN (R 4.5.0)
    ##  progressr          0.18.0  2025-11-06 [1] CRAN (R 4.5.0)
    ##  promises           1.5.0   2025-11-01 [1] CRAN (R 4.5.0)
    ##  purrr            * 1.2.0   2025-11-04 [1] CRAN (R 4.5.0)
    ##  R6                 2.6.1   2025-02-15 [1] CRAN (R 4.5.0)
    ##  RANN               2.6.2   2024-08-25 [1] CRAN (R 4.5.0)
    ##  RColorBrewer       1.1-3   2022-04-03 [1] CRAN (R 4.5.0)
    ##  Rcpp               1.1.0   2025-07-02 [1] CRAN (R 4.5.0)
    ##  RcppAnnoy          0.0.22  2024-01-23 [1] CRAN (R 4.5.0)
    ##  RcppHNSW           0.6.0   2024-02-04 [1] CRAN (R 4.5.0)
    ##  readr            * 2.1.6   2025-11-14 [1] CRAN (R 4.5.2)
    ##  remotes            2.5.0   2024-03-17 [1] CRAN (R 4.5.0)
    ##  reshape2           1.4.5   2025-11-12 [1] CRAN (R 4.5.0)
    ##  reticulate         1.44.1  2025-11-14 [1] CRAN (R 4.5.2)
    ##  rlang              1.1.6   2025-04-11 [1] CRAN (R 4.5.0)
    ##  rmarkdown          2.30    2025-09-28 [1] CRAN (R 4.5.0)
    ##  ROCR               1.0-11  2020-05-02 [1] CRAN (R 4.5.0)
    ##  RSpectra           0.16-2  2024-07-18 [1] CRAN (R 4.5.0)
    ##  rstudioapi         0.17.1  2024-10-22 [1] CRAN (R 4.5.0)
    ##  Rtsne              0.17    2023-12-07 [1] CRAN (R 4.5.0)
    ##  S7                 0.2.1   2025-11-14 [1] CRAN (R 4.5.2)
    ##  scales           * 1.4.0   2025-04-24 [1] CRAN (R 4.5.0)
    ##  scattermore        1.2     2023-06-12 [1] CRAN (R 4.5.0)
    ##  sctransform        0.4.2   2025-04-30 [1] CRAN (R 4.5.0)
    ##  sessioninfo        1.2.3   2025-02-05 [1] CRAN (R 4.5.0)
    ##  Seurat           * 5.4.0   2025-12-14 [1] CRAN (R 4.5.2)
    ##  SeuratObject     * 5.3.0   2025-12-12 [1] CRAN (R 4.5.2)
    ##  shape              1.4.6.1 2024-02-23 [1] CRAN (R 4.5.0)
    ##  shiny              1.12.1  2025-12-09 [1] CRAN (R 4.5.0)
    ##  sp               * 2.2-0   2025-02-01 [1] CRAN (R 4.5.0)
    ##  spam               2.11-1  2025-01-20 [1] CRAN (R 4.5.0)
    ##  spatstat.data      3.1-9   2025-10-18 [1] CRAN (R 4.5.0)
    ##  spatstat.explore   3.6-0   2025-11-22 [1] CRAN (R 4.5.2)
    ##  spatstat.geom      3.6-1   2025-11-20 [1] CRAN (R 4.5.2)
    ##  spatstat.random    3.4-3   2025-11-21 [1] CRAN (R 4.5.2)
    ##  spatstat.sparse    3.1-0   2024-06-21 [1] CRAN (R 4.5.0)
    ##  spatstat.univar    3.1-5   2025-11-17 [1] CRAN (R 4.5.2)
    ##  spatstat.utils     3.2-0   2025-09-20 [1] CRAN (R 4.5.0)
    ##  stringi            1.8.7   2025-03-27 [1] CRAN (R 4.5.0)
    ##  stringr          * 1.6.0   2025-11-04 [1] CRAN (R 4.5.0)
    ##  survival           3.8-3   2024-12-17 [1] CRAN (R 4.5.2)
    ##  tensor             1.5.1   2025-06-17 [1] CRAN (R 4.5.0)
    ##  tibble           * 3.3.0   2025-06-08 [1] CRAN (R 4.5.0)
    ##  tidyr            * 1.3.2   2025-12-19 [1] CRAN (R 4.5.2)
    ##  tidyselect         1.2.1   2024-03-11 [1] CRAN (R 4.5.0)
    ##  tidyverse        * 2.0.0   2023-02-22 [1] CRAN (R 4.5.0)
    ##  timechange         0.3.0   2024-01-18 [1] CRAN (R 4.5.0)
    ##  tzdb               0.5.0   2025-03-15 [1] CRAN (R 4.5.0)
    ##  usethis            3.2.1   2025-09-06 [1] CRAN (R 4.5.0)
    ##  uwot               0.2.4   2025-11-10 [1] CRAN (R 4.5.0)
    ##  vctrs              0.6.5   2023-12-01 [1] CRAN (R 4.5.0)
    ##  viridisLite        0.4.2   2023-05-02 [1] CRAN (R 4.5.0)
    ##  vroom              1.6.7   2025-11-28 [1] CRAN (R 4.5.2)
    ##  withr              3.0.2   2024-10-28 [1] CRAN (R 4.5.0)
    ##  xfun               0.55    2025-12-16 [1] CRAN (R 4.5.2)
    ##  xtable             1.8-4   2019-04-21 [1] CRAN (R 4.5.0)
    ##  yaml               2.3.12  2025-12-10 [1] CRAN (R 4.5.2)
    ##  zoo                1.8-15  2025-12-15 [1] CRAN (R 4.5.2)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
