XR251_GO_score
================
Alexander Lercher
2026-04-29

## LOAD PACKAGES

``` r
library(GO.db)
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

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

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 

``` r
library(org.Mm.eg.db)
```

    ## 

``` r
library(AnnotationDbi)
library(Seurat)
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 
    ## Attaching package: 'sp'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     %over%

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     intersect

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     intersect

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     intersect

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
    ## ✖ lubridate::%within%() masks IRanges::%within%()
    ## ✖ dplyr::collapse()     masks IRanges::collapse()
    ## ✖ dplyr::combine()      masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::desc()         masks IRanges::desc()
    ## ✖ tidyr::expand()       masks S4Vectors::expand()
    ## ✖ dplyr::filter()       masks stats::filter()
    ## ✖ dplyr::first()        masks S4Vectors::first()
    ## ✖ dplyr::lag()          masks stats::lag()
    ## ✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()       masks IRanges::reduce()
    ## ✖ dplyr::rename()       masks S4Vectors::rename()
    ## ✖ lubridate::second()   masks S4Vectors::second()
    ## ✖ lubridate::second<-() masks S4Vectors::second<-()
    ## ✖ dplyr::select()       masks AnnotationDbi::select()
    ## ✖ dplyr::slice()        masks IRanges::slice()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(patchwork)
library(ComplexHeatmap)
```

    ## Loading required package: grid
    ## ========================================
    ## ComplexHeatmap version 2.26.0
    ## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    ## Github page: https://github.com/jokergoo/ComplexHeatmap
    ## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    ## 
    ## If you use it in published research, please cite either one:
    ## - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
    ## - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    ##     genomic data. Bioinformatics 2016.
    ## 
    ## 
    ## The new InteractiveComplexHeatmap package can directly export static 
    ## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(ComplexHeatmap))
    ## ========================================

``` r
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

## LOAD in SEURAT RDS file modified by me

``` r
data <- readRDS(file = "input/seu_obj_v12192025_AL.rds")
```

## LOAD in ISG data from BRAIN and BMEC

``` r
ISG_brain <- read_tsv("input/XR249_ISG_expression_brain.tsv") %>%
  dplyr::pull(gene_id) %>%
  stringr::str_trim()
```

    ## Rows: 112 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): gene_id
    ## dbl (10): PBS_1, PBS_2, PBS_3, IFNa_1, IFNa_2, IFNa_3, log2FoldChange, padj,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ISG_BMEC <- read_tsv("input/XR252_ISG_expression_BMEC.tsv") %>%
  dplyr::pull(gene_id) %>%
  stringr::str_trim()
```

    ## Rows: 951 Columns: 13
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): gene_id
    ## dbl (12): BMEC_PBS_1, BMEC_PBS_2, BMEC_PBS_3, BMEC_PBS_4, BMEC_IFNa_1, BMEC_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Protective_ISG_signature <- read_csv("input/XR256_heatmap_cluster1_genes.csv") %>%
  dplyr::pull(gene_id) %>%
  stringr::str_trim()
```

    ## Rows: 68 Columns: 1
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (1): gene_id
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ISG <- list(ISG_brain = ISG_brain,
            ISG_BMEC = ISG_BMEC,
            protective_ISG = Protective_ISG_signature)
```

## ADD MODULE SCORES

``` r
for(i in 1:length(ISG)){
data <- AddModuleScore(
  object = data,
  features = list(ISG[[i]]),
  name = paste0(names(ISG[i]),"_score")
)
}
```

    ## Warning: The following features are not present in the object: Iigp1c, not
    ## searching for symbol synonyms

    ## Warning: The following features are not present in the object: 2610024D14Rik,
    ## Apol7c, Apon, Atp5if1, Ctps1, Dhh, Duoxa2, Gm20556, Gm45521, Gm8579, Gvin-ps2,
    ## Gzma, Gzmb, Iigp1c, Kars1, LOC100041034, Ly6i, Mir6385, Mir6386, Mir7679, Nkg7,
    ## Or1o4, Phxr4, Rars1, Rd3l, Ren1, Ripply3, Scgb3a1, Serpina3i, Snord12, Snord13,
    ## Snord53, Trim34b, Vars1, Zng1, not searching for symbol synonyms

## PLOT MODULE SCORES

``` r
for(i in 1:length(ISG)){
  colname <- paste0(names(ISG[i]), "_score1")
  
# VIOLIN PLOT
# get expression values
  df <- data@meta.data %>%
    dplyr::select(cell_type, condition, all_of(colname))

# plot
  p <- ggplot(df, aes(x = condition, y = .data[[colname]], fill = condition)) +
        geom_violin() +
        facet_wrap(~cell_type, ncol = 3, scales = "fixed") +  # 3 panels per row
        scale_fill_manual(values = c("00_hours"="#E69F00", "02_hours"="#56B4E9", "16_hours"="#009E73")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              aspect.ratio = 1) +
        ylab("ISG Module Score (log1p)") +
        xlab("Time Point") +
        scale_y_continuous(trans="log1p")
  print(p)
  
# DOT PLOT
# calculate percent expressing and avg score
  df_summary <- data@meta.data %>%
    dplyr::select(cell_type, condition, all_of(colname)) %>%
    group_by(cell_type, condition) %>%
    summarize(
      avg_score = mean(.data[[colname]]),
      pct_expressing = mean(.data[[colname]] > 0) * 100,  # fraction of cells with positive score
      .groups = "drop"
    ) %>%
    group_by(cell_type) %>%
    mutate(avg_score_total = mean(avg_score)) %>%
    ungroup()

# order by highest expressing cell type
  df_summary$cell_type <- factor(df_summary$cell_type,
                                 levels = df_summary %>%
                                   distinct(cell_type, avg_score_total) %>%
                                   arrange(avg_score_total) %>%
                                   pull(cell_type))
  
# save table to file
  df_summary_save <- df_summary %>%
    dplyr::select(!avg_score_total) %>%
    write_tsv(paste0("output/ISG_score/",names(ISG[i]),"_module_scores.tsv"))

# plot
  p <- ggplot(df_summary, aes(x = condition, y = cell_type)) +
    geom_point(aes(size = pct_expressing, color = avg_score)) +
    scale_color_gradient(low = "lightgrey", high = "firebrick") +
    scale_size(range = c(3, 10)) +  # adjust point size range
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Time Point") +
    ylab("Cell Type") +
    ggtitle("ISG Brain Module Score")
  print(p)

  # HEATMAP
  # get matrix of module set expression
  mat_avg <- df_summary %>%
    dplyr::select(condition, cell_type, avg_score) %>%
    tidyr::pivot_wider(
      names_from = cell_type,
      values_from = avg_score
    ) %>%
    tibble::column_to_rownames("condition") %>%
    as.matrix()
  
  # get annotation dataframe of average pct_expressing across all time points
  pct_df <- df_summary %>%
    mutate(cell_type = as.character(cell_type)) %>%
    group_by(cell_type) %>%
    summarise(pct_expressing = mean(pct_expressing), .groups = "drop") %>%
    tibble::column_to_rownames("cell_type")
  
  # enforce exact column order
  pct_df <- pct_df[colnames(mat_avg), , drop = FALSE]
  
  # safety check
  stopifnot(!any(is.na(pct_df$pct_expressing)))
  
  # define annotation (prepare to plot)
  ha <- HeatmapAnnotation(
    df = pct_df,
    col = list(
      pct_expressing = colorRamp2(
        c(0,100),
        c("#f7fcf5", "#238b45")
      )
    ),
    show_annotation_name = FALSE
  )
  
  # define color scale for module scores
  # lower end of signal usually around 0.03 in Seurat module scores
  low_signal <- (abs(max(mat_avg) - min(mat_avg)))/9.5 + min(mat_avg)
  middle_signal <- (abs(max(mat_avg) - min(mat_avg)))/6 + min(mat_avg)
  
  # define color range
  col_scale <- colorRamp2(
    c(min(mat_avg, na.rm= TRUE), low_signal, middle_signal, max(mat_avg, na.rm = TRUE)),
    c("#8c6bb1", "#cce6ff", "#fff5b1", "#fb8072")
  )
  
  # plot and save heatmap
  ht <- ComplexHeatmap::Heatmap(
      mat_avg,
      name = paste0(names(ISG[i])," module score"),
      col = col_scale,
      top_annotation = ha,
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      row_names_gp = gpar(fontsize = 10),
      column_names_gp = gpar(fontsize = 9),
      show_row_dend = FALSE,
      show_column_dend = TRUE,
      heatmap_legend_param = list(title = "Module score")
    )
  
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
    
  # UMAP per TIME POINT
  # compute global min/max for feature
  global_vals <-  data@meta.data[[colname]]
  min_val <- min(global_vals, na.rm = TRUE)
  max_val <- max(global_vals, na.rm = TRUE)
  middle_val <- 0.2
  
  # split data by time points
  objs <- SplitObject(data, split.by = "condition")
  
  p_list <- lapply(names(objs), function(nm) {
    FeaturePlot(
      objs[[nm]],
      features = colname,
      reduction = "umap",
      order = TRUE,
      pt.size = 1
    ) +
      scale_color_viridis_c(option = "viridis",
                            limits = c(min_val, max_val),
                            values= rescale(c(min_val,middle_val, max_val))) +
      ggtitle(nm) +               # set title to condition name
      theme(aspect.ratio = 1,
            axis.title = element_text(size = 26),
            axis.text = element_text(size = 26),
            plot.title = element_text(size = 24),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 26))
  })
  print(p_list)
}
```

![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-1.png)<!-- -->![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-2.png)<!-- -->![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-3.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

    ## [[1]]

![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-4.png)<!-- -->

    ## 
    ## [[2]]

![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-5.png)<!-- -->

    ## 
    ## [[3]]

![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-6.png)<!-- -->
![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-7.png)<!-- -->![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-8.png)<!-- -->![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-9.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

    ## [[1]]

![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-10.png)<!-- -->

    ## 
    ## [[2]]

![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-11.png)<!-- -->

    ## 
    ## [[3]]

![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-12.png)<!-- -->
![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-13.png)<!-- -->![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-14.png)<!-- -->![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-15.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

    ## [[1]]

![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-16.png)<!-- -->

    ## 
    ## [[2]]

![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-17.png)<!-- -->

    ## 
    ## [[3]]

![](XR251_ISG_protective_gene_score_files/figure-gfm/plot%20module%20scores-18.png)<!-- -->

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
    ##  AnnotationDbi    * 1.72.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  Biobase          * 2.70.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocGenerics     * 0.56.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  Biostrings         2.78.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  bit                4.6.0   2025-03-06 [1] CRAN (R 4.5.0)
    ##  bit64              4.6.0-1 2025-01-16 [1] CRAN (R 4.5.0)
    ##  blob               1.2.4   2023-03-17 [1] CRAN (R 4.5.0)
    ##  cachem             1.1.0   2024-05-16 [1] CRAN (R 4.5.0)
    ##  circlize         * 0.4.17  2025-12-08 [1] CRAN (R 4.5.2)
    ##  cli                3.6.5   2025-04-23 [1] CRAN (R 4.5.0)
    ##  clue               0.3-66  2024-11-13 [1] CRAN (R 4.5.0)
    ##  cluster            2.1.8.1 2025-03-12 [1] CRAN (R 4.5.2)
    ##  codetools          0.2-20  2024-03-31 [1] CRAN (R 4.5.2)
    ##  colorspace         2.1-2   2025-09-22 [1] CRAN (R 4.5.0)
    ##  ComplexHeatmap   * 2.26.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  cowplot            1.2.0   2025-07-07 [1] CRAN (R 4.5.0)
    ##  crayon             1.5.3   2024-06-20 [1] CRAN (R 4.5.0)
    ##  data.table         1.18.0  2025-12-24 [1] CRAN (R 4.5.2)
    ##  DBI                1.2.3   2024-06-02 [1] CRAN (R 4.5.0)
    ##  deldir             2.0-4   2024-02-28 [1] CRAN (R 4.5.0)
    ##  devtools           2.4.6   2025-10-03 [1] CRAN (R 4.5.0)
    ##  digest             0.6.39  2025-11-19 [1] CRAN (R 4.5.2)
    ##  doParallel         1.0.17  2022-02-07 [1] CRAN (R 4.5.0)
    ##  dotCall64          1.2     2024-10-04 [1] CRAN (R 4.5.0)
    ##  dplyr            * 1.1.4   2023-11-17 [1] CRAN (R 4.5.0)
    ##  ellipsis           0.3.2   2021-04-29 [1] CRAN (R 4.5.0)
    ##  evaluate           1.0.5   2025-08-27 [1] CRAN (R 4.5.0)
    ##  farver             2.1.2   2024-05-13 [1] CRAN (R 4.5.0)
    ##  fastDummies        1.7.5   2025-01-20 [1] CRAN (R 4.5.0)
    ##  fastmap            1.2.0   2024-05-15 [1] CRAN (R 4.5.0)
    ##  fitdistrplus       1.2-4   2025-07-03 [1] CRAN (R 4.5.0)
    ##  forcats          * 1.0.1   2025-09-25 [1] CRAN (R 4.5.0)
    ##  foreach            1.5.2   2022-02-02 [1] CRAN (R 4.5.0)
    ##  fs                 1.6.6   2025-04-12 [1] CRAN (R 4.5.0)
    ##  future             1.68.0  2025-11-17 [1] CRAN (R 4.5.2)
    ##  future.apply       1.20.1  2025-12-09 [1] CRAN (R 4.5.0)
    ##  generics         * 0.1.4   2025-05-09 [1] CRAN (R 4.5.0)
    ##  GetoptLong         1.1.0   2025-11-28 [1] CRAN (R 4.5.2)
    ##  ggplot2          * 4.0.1   2025-11-14 [1] CRAN (R 4.5.2)
    ##  ggrepel            0.9.6   2024-09-07 [1] CRAN (R 4.5.0)
    ##  ggridges           0.5.7   2025-08-27 [1] CRAN (R 4.5.0)
    ##  GlobalOptions      0.1.3   2025-11-28 [1] CRAN (R 4.5.2)
    ##  globals            0.18.0  2025-05-08 [1] CRAN (R 4.5.0)
    ##  glue               1.8.0   2024-09-30 [1] CRAN (R 4.5.0)
    ##  GO.db            * 3.22.0  2025-12-09 [1] Bioconductor
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
    ##  IRanges          * 2.44.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  irlba              2.3.5.1 2022-10-03 [1] CRAN (R 4.5.0)
    ##  iterators          1.0.14  2022-02-05 [1] CRAN (R 4.5.0)
    ##  jsonlite           2.0.0   2025-03-27 [1] CRAN (R 4.5.0)
    ##  KEGGREST           1.50.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
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
    ##  magick             2.9.0   2025-09-08 [1] CRAN (R 4.5.0)
    ##  magrittr           2.0.4   2025-09-12 [1] CRAN (R 4.5.0)
    ##  MASS               7.3-65  2025-02-28 [1] CRAN (R 4.5.2)
    ##  Matrix             1.7-4   2025-08-28 [1] CRAN (R 4.5.0)
    ##  matrixStats        1.5.0   2025-01-07 [1] CRAN (R 4.5.0)
    ##  memoise            2.0.1   2021-11-26 [1] CRAN (R 4.5.0)
    ##  mime               0.13    2025-03-17 [1] CRAN (R 4.5.0)
    ##  miniUI             0.1.2   2025-04-17 [1] CRAN (R 4.5.0)
    ##  nlme               3.1-168 2025-03-31 [1] CRAN (R 4.5.2)
    ##  org.Mm.eg.db     * 3.22.0  2025-12-09 [1] Bioconductor
    ##  otel               0.2.0   2025-08-29 [1] CRAN (R 4.5.0)
    ##  parallelly         1.46.0  2025-12-12 [1] CRAN (R 4.5.2)
    ##  patchwork        * 1.3.2   2025-08-25 [1] CRAN (R 4.5.0)
    ##  pbapply            1.7-4   2025-07-20 [1] CRAN (R 4.5.0)
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
    ##  rjson              0.2.23  2024-09-16 [1] CRAN (R 4.5.0)
    ##  rlang              1.1.6   2025-04-11 [1] CRAN (R 4.5.0)
    ##  rmarkdown          2.30    2025-09-28 [1] CRAN (R 4.5.0)
    ##  ROCR               1.0-11  2020-05-02 [1] CRAN (R 4.5.0)
    ##  RSpectra           0.16-2  2024-07-18 [1] CRAN (R 4.5.0)
    ##  RSQLite            2.4.5   2025-11-30 [1] CRAN (R 4.5.2)
    ##  rstudioapi         0.17.1  2024-10-22 [1] CRAN (R 4.5.0)
    ##  Rtsne              0.17    2023-12-07 [1] CRAN (R 4.5.0)
    ##  S4Vectors        * 0.48.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  S7                 0.2.1   2025-11-14 [1] CRAN (R 4.5.2)
    ##  scales           * 1.4.0   2025-04-24 [1] CRAN (R 4.5.0)
    ##  scattermore        1.2     2023-06-12 [1] CRAN (R 4.5.0)
    ##  sctransform        0.4.2   2025-04-30 [1] CRAN (R 4.5.0)
    ##  Seqinfo            1.0.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
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
    ##  XVector            0.50.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  yaml               2.3.12  2025-12-10 [1] CRAN (R 4.5.2)
    ##  zoo                1.8-15  2025-12-15 [1] CRAN (R 4.5.2)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
