XR251_DEG_plot
================
Alexander Lercher
2026-04-29

## LOAD PACKAGES

``` r
library(ggplot2)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.6
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ lubridate 1.9.4     ✔ tibble    3.3.0
    ## ✔ purrr     1.2.0     ✔ tidyr     1.3.2
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(knitr)
library(kableExtra)
```

    ## 
    ## Attaching package: 'kableExtra'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

## LOAD in DEG files

``` r
list_of_files <- list.files(
  path = "output/DEG",
  pattern = "\\.tsv$",
  full.names = TRUE
)

# read in files and put them in a list
data_list <- lapply(list_of_files, read_tsv)
```

    ## Rows: 29957 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): gene, cluster, comparison
    ## dbl (5): p_val, avg_log2FC, pct.1, pct.2, p_val_adj
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 20498 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): gene, cluster, comparison
    ## dbl (5): p_val, avg_log2FC, pct.1, pct.2, p_val_adj
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 24355 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): gene, cluster, comparison
    ## dbl (5): p_val, avg_log2FC, pct.1, pct.2, p_val_adj
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# rename list acc to input file
names(data_list) <- tools::file_path_sans_ext(basename(list_of_files))
```

## IDENTIFY NUMBER OF DEGs

``` r
# set DEG cutoffs
log2FC_cutoff <- 0.5
padj_cutoff <- 0.05

# filter for DEGs
data_list_sig <- lapply(data_list,
                        dplyr::filter,
                        abs(avg_log2FC) >= log2FC_cutoff,
                        p_val_adj < padj_cutoff)

# count number of DEG
deg_counts_list <- lapply(
  data_list_sig,
  function(df) {
    df %>%
      mutate(direction = if_else(avg_log2FC > 0, "up", "down"),
             direction = factor(direction, levels = c("up", "down"))) %>%
      count(cluster, direction) %>%
      mutate(n = if_else(direction == "down", -n, n))
  }
)

# print number of DEG
lapply(names(deg_counts_list), function(nm) {
  cat("DEG counts for:", nm, "\n")
  
  deg_counts_list[[nm]] %>%
    kable(
      format = "pipe",
      digits = 0,
      caption = paste("DEG counts for", nm)
    )
})
```

    ## DEG counts for: DEG_02_vs_00 
    ## DEG counts for: DEG_16_vs_00 
    ## DEG counts for: DEG_16_vs_02

    ## [[1]]
    ## 
    ## 
    ## Table: DEG counts for DEG_02_vs_00
    ## 
    ## |cluster                      |direction |    n|
    ## |:----------------------------|:---------|----:|
    ## |Astrocytes                   |up        |  276|
    ## |Astrocytes                   |down      | -451|
    ## |Endothelial cells            |up        |  366|
    ## |Endothelial cells            |down      | -241|
    ## |Excitatory neurons           |up        |  667|
    ## |Excitatory neurons           |down      | -178|
    ## |Fibroblast                   |up        |  233|
    ## |Fibroblast                   |down      |  -91|
    ## |Fibroblast-like              |up        |  221|
    ## |Fibroblast-like              |down      | -187|
    ## |Inhibitory neurons           |up        |  305|
    ## |Inhibitory neurons           |down      | -437|
    ## |Interneurons                 |up        |  138|
    ## |Interneurons                 |down      |  -23|
    ## |Microglia                    |up        |  762|
    ## |Microglia                    |down      | -251|
    ## |OPCs                         |up        |  133|
    ## |OPCs                         |down      | -151|
    ## |Oligodendrocytes             |up        |  185|
    ## |Oligodendrocytes             |down      | -183|
    ## |Pericytes                    |up        |   90|
    ## |Pericytes                    |down      |  -48|
    ## |Ventricular Epithelial cells |up        |  200|
    ## |Ventricular Epithelial cells |down      | -296|
    ## 
    ## [[2]]
    ## 
    ## 
    ## Table: DEG counts for DEG_16_vs_00
    ## 
    ## |cluster                      |direction |    n|
    ## |:----------------------------|:---------|----:|
    ## |Astrocytes                   |up        |   86|
    ## |Astrocytes                   |down      | -160|
    ## |Endothelial cells            |up        |   26|
    ## |Endothelial cells            |down      |  -28|
    ## |Excitatory neurons           |up        |   64|
    ## |Excitatory neurons           |down      |  -32|
    ## |Fibroblast                   |up        |   12|
    ## |Fibroblast                   |down      |  -13|
    ## |Fibroblast-like              |up        |   18|
    ## |Fibroblast-like              |down      |  -28|
    ## |Inhibitory neurons           |up        |   27|
    ## |Inhibitory neurons           |down      | -167|
    ## |Interneurons                 |up        |   12|
    ## |Interneurons                 |down      |   -1|
    ## |Microglia                    |up        | 1050|
    ## |Microglia                    |down      | -345|
    ## |OPCs                         |up        |   32|
    ## |OPCs                         |down      |  -52|
    ## |Oligodendrocytes             |up        |   35|
    ## |Oligodendrocytes             |down      |  -57|
    ## |Pericytes                    |up        |   23|
    ## |Pericytes                    |down      |  -21|
    ## |Ventricular Epithelial cells |up        |    2|
    ## |Ventricular Epithelial cells |down      |  -21|
    ## 
    ## [[3]]
    ## 
    ## 
    ## Table: DEG counts for DEG_16_vs_02
    ## 
    ## |cluster                      |direction |    n|
    ## |:----------------------------|:---------|----:|
    ## |Astrocytes                   |up        |  171|
    ## |Astrocytes                   |down      | -228|
    ## |Endothelial cells            |up        |  160|
    ## |Endothelial cells            |down      | -305|
    ## |Excitatory neurons           |up        |  121|
    ## |Excitatory neurons           |down      | -293|
    ## |Fibroblast                   |up        |   63|
    ## |Fibroblast                   |down      | -207|
    ## |Fibroblast-like              |up        |   26|
    ## |Fibroblast-like              |down      |  -88|
    ## |Inhibitory neurons           |up        |  182|
    ## |Inhibitory neurons           |down      | -292|
    ## |Interneurons                 |up        |    5|
    ## |Interneurons                 |down      |  -25|
    ## |Microglia                    |up        |  303|
    ## |Microglia                    |down      | -277|
    ## |OPCs                         |up        |   17|
    ## |OPCs                         |down      | -112|
    ## |Oligodendrocytes             |up        |  149|
    ## |Oligodendrocytes             |down      | -195|
    ## |Pericytes                    |up        |   44|
    ## |Pericytes                    |down      |  -46|
    ## |Ventricular Epithelial cells |up        |  230|
    ## |Ventricular Epithelial cells |down      | -340|

``` r
# save number of DEGs
all_deg_counts <- dplyr::bind_rows(
  deg_counts_list,
  .id = "comparison"
) %>%
  write_csv("output/DEG/DEG_counts.csv")
```

## PLOT NUMBER OF DEGs

``` r
# reorder to match the dotplot generated previously
new_order <- rev(c("Excitatory neurons",
                   "Inhibitory neurons",
                   "Interneurons",
                   "Oligodendrocytes",
                   "OPCs",
                   "Astrocytes",
                   "Microglia",
                   "Endothelial cells",
                   "Pericytes",
                   "Fibroblast-like",
                   "Fibroblast",
                   "Ventricular Epithelial cells"))

deg_counts_list <- lapply(
  deg_counts_list,
  function(df) {
    df %>%
      mutate(cluster = factor(cluster, levels = new_order))
  }
)

# PLOT BAR GRAPH
# maximum changed genes for any cluster and any comparison
max_abs_all <- max(
  sapply(deg_counts_list, function(df) max(abs(df$n)))
)

# generate plots, save in list
plot_list <- setNames(
  lapply(
  names(deg_counts_list),
  function(nm) {
    ggplot(deg_counts_list[[nm]],
           aes(x = factor(cluster), y = n, fill = direction)) +
      geom_col() +
      scale_fill_brewer(palette = "Set1") +
      labs(
        title = nm,
        x = "Cluster",
        y = "Number of DEGs"
      ) +
      ylim(-max_abs_all,max_abs_all) +
      theme_bw() +
      theme(aspect.ratio = 1) +
      coord_flip()
  }
),
names(deg_counts_list))
plot_list
```

    ## $DEG_02_vs_00

![](XR251_DEG_plots_files/figure-gfm/plot%20number%20of%20DEGs-1.png)<!-- -->

    ## 
    ## $DEG_16_vs_00

![](XR251_DEG_plots_files/figure-gfm/plot%20number%20of%20DEGs-2.png)<!-- -->

    ## 
    ## $DEG_16_vs_02

![](XR251_DEG_plots_files/figure-gfm/plot%20number%20of%20DEGs-3.png)<!-- -->

``` r
print(plot_list)
```

    ## $DEG_02_vs_00

![](XR251_DEG_plots_files/figure-gfm/plot%20number%20of%20DEGs-4.png)<!-- -->

    ## 
    ## $DEG_16_vs_00

![](XR251_DEG_plots_files/figure-gfm/plot%20number%20of%20DEGs-5.png)<!-- -->

    ## 
    ## $DEG_16_vs_02

![](XR251_DEG_plots_files/figure-gfm/plot%20number%20of%20DEGs-6.png)<!-- -->

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
    ##  kableExtra   * 1.4.0   2024-01-24 [1] CRAN (R 4.5.0)
    ##  knitr        * 1.51    2025-12-20 [1] CRAN (R 4.5.2)
    ##  labeling       0.4.3   2023-08-29 [1] CRAN (R 4.5.0)
    ##  lifecycle      1.0.4   2023-11-07 [1] CRAN (R 4.5.0)
    ##  lubridate    * 1.9.4   2024-12-08 [1] CRAN (R 4.5.0)
    ##  magrittr       2.0.4   2025-09-12 [1] CRAN (R 4.5.0)
    ##  memoise        2.0.1   2021-11-26 [1] CRAN (R 4.5.0)
    ##  otel           0.2.0   2025-08-29 [1] CRAN (R 4.5.0)
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
    ##  svglite        2.2.2   2025-10-21 [1] CRAN (R 4.5.0)
    ##  systemfonts    1.3.1   2025-10-01 [1] CRAN (R 4.5.0)
    ##  textshaping    1.0.4   2025-10-10 [1] CRAN (R 4.5.0)
    ##  tibble       * 3.3.0   2025-06-08 [1] CRAN (R 4.5.0)
    ##  tidyr        * 1.3.2   2025-12-19 [1] CRAN (R 4.5.2)
    ##  tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.5.0)
    ##  tidyverse    * 2.0.0   2023-02-22 [1] CRAN (R 4.5.0)
    ##  timechange     0.3.0   2024-01-18 [1] CRAN (R 4.5.0)
    ##  tzdb           0.5.0   2025-03-15 [1] CRAN (R 4.5.0)
    ##  usethis        3.2.1   2025-09-06 [1] CRAN (R 4.5.0)
    ##  vctrs          0.6.5   2023-12-01 [1] CRAN (R 4.5.0)
    ##  viridisLite    0.4.2   2023-05-02 [1] CRAN (R 4.5.0)
    ##  vroom          1.6.7   2025-11-28 [1] CRAN (R 4.5.2)
    ##  withr          3.0.2   2024-10-28 [1] CRAN (R 4.5.0)
    ##  xfun           0.55    2025-12-16 [1] CRAN (R 4.5.2)
    ##  xml2           1.5.1   2025-12-01 [1] CRAN (R 4.5.2)
    ##  yaml           2.3.12  2025-12-10 [1] CRAN (R 4.5.2)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
