XR251_GSEA_Enrichment
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
library(clusterProfiler)
```

    ## 
    ## clusterProfiler v4.18.4 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## G Yu. Thirteen years of clusterProfiler. The Innovation. 2024,
    ## 5(6):100722
    ## 
    ## Attaching package: 'clusterProfiler'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(org.Mm.eg.db)
```

    ## Loading required package: AnnotationDbi
    ## Loading required package: stats4
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
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## Loading required package: IRanges
    ## Loading required package: S4Vectors
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     rename
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
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     slice
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
    ## 
    ## Attaching package: 'AnnotationDbi'
    ## 
    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     select
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
library(dplyr)
library(ggplot2)
library(viridis)
```

    ## Loading required package: viridisLite

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
# put them into one big dataframe
data <- bind_rows(data_list)
```

## EXTRACT SIGNIFICANT GENES IN CELL TYPE AND COMPARISON OF INTEREST

``` r
# set significance criteria
pval_cutoff <- 0.05
log2_cutoff <- 0.5

# genes significant in condition of interest and cluster of interest
sig_COI <- data %>%
  dplyr::filter(
    p_val_adj < pval_cutoff,
    abs(avg_log2FC) >= log2_cutoff
  )
```

## ENRICHMENT ANALYSES

``` r
# extract clusters and comparisons
cluster_name <- unique(data$cluster)
comparison_name <- unique(data$comparison)

# loop over each comparison and within it each cluster to get DEGs
# only interested in 16h vs 0h
i <- 2
#for(i in 1:length(comparison_name)){
  genes_of_interest <- sig_COI %>%
    dplyr::filter(comparison == comparison_name[[i]])

# only interested in microglia
y <- 8
  #for(y in 1:length(cluster_name)){
    genes_of_interest_per_cluster <- genes_of_interest %>%
      dplyr::filter(cluster == cluster_name[[y]])
```

## GSEA / gseGO enrichment

``` r
gse <- tryCatch(
  {
    # create ranked gene vector
    gene_list <- genes_of_interest_per_cluster %>%
      dplyr::mutate(logFC_signed = avg_log2FC * -log10(p_val_adj + 1e-10)) %>%
      dplyr::arrange(desc(logFC_signed)) %>%
      dplyr::select(gene, logFC_signed)
    
    gene_list_vect <- gene_list$logFC_signed
    names(gene_list_vect) <- gene_list$gene
    
    # run gseGO
    gseGO(
      geneList     = gene_list_vect,
      OrgDb        = org.Mm.eg.db,
      ont          = "BP",
      keyType      = "SYMBOL",
      minGSSize    = 10,
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      nPermSimple  = 10000,
      verbose      = FALSE
    )
  },
  error = function(e) NULL
)
```

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (0.07% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some of the pathways the P-values were likely overestimated. For
    ## such pathways log2err is set to NA.

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some pathways, in reality P-values are less than 1e-10. You can
    ## set the `eps` argument to zero for better estimation.

``` r
if(is.null(gse)){
  message("GSEA returned NULL for cluster ", cluster_name[[y]], " | comparison ", comparison_name[[i]])
} else if(nrow(as.data.frame(gse)) == 0){
  message("No enriched GO BP terms found for cluster ", cluster_name[[y]], " | comparison ", comparison_name[[i]])
} else {
  
  # save full GSEA table
  gse %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(p.adjust)) %>%
    arrange(desc(NES)) %>%
    write_tsv(
      file = paste0(
        "output/GSEA_enrichments/",
        comparison_name[[i]], "/",
        cluster_name[[y]], "_GSEA.tsv"
    )
  )
  
  # simplify GO terms
  gse_simpl <- simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min)
  gse_df <- as.data.frame(gse_simpl) %>% filter(!is.na(p.adjust))
  
  #---------------------------------------------
  # Positive NES (top 10)
  #---------------------------------------------
  top_pos <- gse_df %>%
    filter(NES > 0) %>%
    slice_max(NES, n = 10, with_ties = FALSE) %>%
    mutate(Description = forcats::fct_reorder(Description, NES))
  
  if(nrow(top_pos) > 0){
    p_pos <- ggplot(top_pos, aes(x = NES, y = Description)) +
      geom_point(aes(size = setSize, color = p.adjust)) +
      scale_size_continuous(range = c(8,20)) +
      scale_color_viridis_c(option = "viridis", direction = 1, name = "Adjusted p-value") +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
      ) +
      labs(
        title = paste0("Top 10 positively enriched GO BP terms: ", cluster_name[[y]], " | ", comparison_name[[i]]),
        x = "NES (upregulated)", y = ""
      )
    
    print(p_pos)
  }
  
  }
```

![](XR251_GSEA_enrichment_files/figure-gfm/GSEA%20gseGO%20enrichment-1.png)<!-- -->

``` r
print(paste0(cluster_name[[y]]," done!"," We are in comparison ",comparison_name[[i]],". ",length(cluster_name)-y," of ",length(cluster_name),
             " clusters in this comparison to go. ",length(comparison_name)-i," comparisons to go in total."))
```

    ## [1] "Microglia done! We are in comparison 16_vs_00. 4 of 12 clusters in this comparison to go. 1 comparisons to go in total."

``` r
#  }
#}
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
    ##  package           * version date (UTC) lib source
    ##  AnnotationDbi     * 1.72.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  ape                 5.8-1   2024-12-16 [1] CRAN (R 4.5.0)
    ##  aplot               0.2.9   2025-09-12 [1] CRAN (R 4.5.0)
    ##  Biobase           * 2.70.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocGenerics      * 0.56.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocParallel        1.44.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  Biostrings          2.78.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  bit                 4.6.0   2025-03-06 [1] CRAN (R 4.5.0)
    ##  bit64               4.6.0-1 2025-01-16 [1] CRAN (R 4.5.0)
    ##  blob                1.2.4   2023-03-17 [1] CRAN (R 4.5.0)
    ##  cachem              1.1.0   2024-05-16 [1] CRAN (R 4.5.0)
    ##  cli                 3.6.5   2025-04-23 [1] CRAN (R 4.5.0)
    ##  cluster             2.1.8.1 2025-03-12 [1] CRAN (R 4.5.2)
    ##  clusterProfiler   * 4.18.4  2025-12-15 [1] Bioconductor 3.22 (R 4.5.2)
    ##  codetools           0.2-20  2024-03-31 [1] CRAN (R 4.5.2)
    ##  cowplot             1.2.0   2025-07-07 [1] CRAN (R 4.5.0)
    ##  crayon              1.5.3   2024-06-20 [1] CRAN (R 4.5.0)
    ##  data.table          1.18.0  2025-12-24 [1] CRAN (R 4.5.2)
    ##  DBI                 1.2.3   2024-06-02 [1] CRAN (R 4.5.0)
    ##  devtools            2.4.6   2025-10-03 [1] CRAN (R 4.5.0)
    ##  digest              0.6.39  2025-11-19 [1] CRAN (R 4.5.2)
    ##  DOSE                4.4.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  dplyr             * 1.1.4   2023-11-17 [1] CRAN (R 4.5.0)
    ##  ellipsis            0.3.2   2021-04-29 [1] CRAN (R 4.5.0)
    ##  enrichplot          1.30.4  2025-12-01 [1] Bioconductor 3.22 (R 4.5.0)
    ##  evaluate            1.0.5   2025-08-27 [1] CRAN (R 4.5.0)
    ##  farver              2.1.2   2024-05-13 [1] CRAN (R 4.5.0)
    ##  fastmap             1.2.0   2024-05-15 [1] CRAN (R 4.5.0)
    ##  fastmatch           1.1-6   2024-12-23 [1] CRAN (R 4.5.0)
    ##  fgsea               1.36.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  fontBitstreamVera   0.1.1   2017-02-01 [1] CRAN (R 4.5.0)
    ##  fontLiberation      0.1.0   2016-10-15 [1] CRAN (R 4.5.0)
    ##  fontquiver          0.2.1   2017-02-01 [1] CRAN (R 4.5.0)
    ##  forcats           * 1.0.1   2025-09-25 [1] CRAN (R 4.5.0)
    ##  fs                  1.6.6   2025-04-12 [1] CRAN (R 4.5.0)
    ##  gdtools             0.4.4   2025-10-06 [1] CRAN (R 4.5.0)
    ##  generics          * 0.1.4   2025-05-09 [1] CRAN (R 4.5.0)
    ##  ggforce             0.5.0   2025-06-18 [1] CRAN (R 4.5.0)
    ##  ggfun               0.2.0   2025-07-15 [1] CRAN (R 4.5.0)
    ##  ggiraph             0.9.2   2025-10-07 [1] CRAN (R 4.5.0)
    ##  ggnewscale          0.5.2   2025-06-20 [1] CRAN (R 4.5.0)
    ##  ggplot2           * 4.0.1   2025-11-14 [1] CRAN (R 4.5.2)
    ##  ggplotify           0.1.3   2025-09-20 [1] CRAN (R 4.5.0)
    ##  ggrepel             0.9.6   2024-09-07 [1] CRAN (R 4.5.0)
    ##  ggtangle            0.0.9   2025-11-30 [1] CRAN (R 4.5.2)
    ##  ggtree              4.0.3   2025-12-18 [1] Bioconductor 3.22 (R 4.5.2)
    ##  glue                1.8.0   2024-09-30 [1] CRAN (R 4.5.0)
    ##  GO.db               3.22.0  2025-12-09 [1] Bioconductor
    ##  GOSemSim            2.36.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  gridExtra           2.3     2017-09-09 [1] CRAN (R 4.5.0)
    ##  gridGraphics        0.5-1   2020-12-13 [1] CRAN (R 4.5.0)
    ##  gson                0.1.0   2023-03-07 [1] CRAN (R 4.5.0)
    ##  gtable              0.3.6   2024-10-25 [1] CRAN (R 4.5.0)
    ##  hms                 1.1.4   2025-10-17 [1] CRAN (R 4.5.0)
    ##  htmltools           0.5.9   2025-12-04 [1] CRAN (R 4.5.2)
    ##  htmlwidgets         1.6.4   2023-12-06 [1] CRAN (R 4.5.0)
    ##  httr                1.4.7   2023-08-15 [1] CRAN (R 4.5.0)
    ##  igraph              2.2.1   2025-10-27 [1] CRAN (R 4.5.0)
    ##  IRanges           * 2.44.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  jsonlite            2.0.0   2025-03-27 [1] CRAN (R 4.5.0)
    ##  KEGGREST            1.50.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  knitr               1.51    2025-12-20 [1] CRAN (R 4.5.2)
    ##  labeling            0.4.3   2023-08-29 [1] CRAN (R 4.5.0)
    ##  lattice             0.22-7  2025-04-02 [1] CRAN (R 4.5.2)
    ##  lazyeval            0.2.2   2019-03-15 [1] CRAN (R 4.5.0)
    ##  lifecycle           1.0.4   2023-11-07 [1] CRAN (R 4.5.0)
    ##  lubridate         * 1.9.4   2024-12-08 [1] CRAN (R 4.5.0)
    ##  magrittr            2.0.4   2025-09-12 [1] CRAN (R 4.5.0)
    ##  MASS                7.3-65  2025-02-28 [1] CRAN (R 4.5.2)
    ##  Matrix              1.7-4   2025-08-28 [1] CRAN (R 4.5.0)
    ##  memoise             2.0.1   2021-11-26 [1] CRAN (R 4.5.0)
    ##  nlme                3.1-168 2025-03-31 [1] CRAN (R 4.5.2)
    ##  org.Mm.eg.db      * 3.22.0  2025-12-09 [1] Bioconductor
    ##  otel                0.2.0   2025-08-29 [1] CRAN (R 4.5.0)
    ##  patchwork           1.3.2   2025-08-25 [1] CRAN (R 4.5.0)
    ##  pillar              1.11.1  2025-09-17 [1] CRAN (R 4.5.0)
    ##  pkgbuild            1.4.8   2025-05-26 [1] CRAN (R 4.5.0)
    ##  pkgconfig           2.0.3   2019-09-22 [1] CRAN (R 4.5.0)
    ##  pkgload             1.4.1   2025-09-23 [1] CRAN (R 4.5.0)
    ##  plyr                1.8.9   2023-10-02 [1] CRAN (R 4.5.0)
    ##  png                 0.1-8   2022-11-29 [1] CRAN (R 4.5.0)
    ##  polyclip            1.10-7  2024-07-23 [1] CRAN (R 4.5.0)
    ##  purrr             * 1.2.0   2025-11-04 [1] CRAN (R 4.5.0)
    ##  qvalue              2.42.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  R.methodsS3         1.8.2   2022-06-13 [1] CRAN (R 4.5.0)
    ##  R.oo                1.27.1  2025-05-02 [1] CRAN (R 4.5.0)
    ##  R.utils             2.13.0  2025-02-24 [1] CRAN (R 4.5.0)
    ##  R6                  2.6.1   2025-02-15 [1] CRAN (R 4.5.0)
    ##  rappdirs            0.3.3   2021-01-31 [1] CRAN (R 4.5.0)
    ##  RColorBrewer        1.1-3   2022-04-03 [1] CRAN (R 4.5.0)
    ##  Rcpp                1.1.0   2025-07-02 [1] CRAN (R 4.5.0)
    ##  readr             * 2.1.6   2025-11-14 [1] CRAN (R 4.5.2)
    ##  remotes             2.5.0   2024-03-17 [1] CRAN (R 4.5.0)
    ##  reshape2            1.4.5   2025-11-12 [1] CRAN (R 4.5.0)
    ##  rlang               1.1.6   2025-04-11 [1] CRAN (R 4.5.0)
    ##  rmarkdown           2.30    2025-09-28 [1] CRAN (R 4.5.0)
    ##  RSQLite             2.4.5   2025-11-30 [1] CRAN (R 4.5.2)
    ##  rstudioapi          0.17.1  2024-10-22 [1] CRAN (R 4.5.0)
    ##  S4Vectors         * 0.48.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  S7                  0.2.1   2025-11-14 [1] CRAN (R 4.5.2)
    ##  scales              1.4.0   2025-04-24 [1] CRAN (R 4.5.0)
    ##  scatterpie          0.2.6   2025-09-12 [1] CRAN (R 4.5.0)
    ##  Seqinfo             1.0.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  sessioninfo         1.2.3   2025-02-05 [1] CRAN (R 4.5.0)
    ##  stringi             1.8.7   2025-03-27 [1] CRAN (R 4.5.0)
    ##  stringr           * 1.6.0   2025-11-04 [1] CRAN (R 4.5.0)
    ##  systemfonts         1.3.1   2025-10-01 [1] CRAN (R 4.5.0)
    ##  tibble            * 3.3.0   2025-06-08 [1] CRAN (R 4.5.0)
    ##  tidydr              0.0.6   2025-07-25 [1] CRAN (R 4.5.0)
    ##  tidyr             * 1.3.2   2025-12-19 [1] CRAN (R 4.5.2)
    ##  tidyselect          1.2.1   2024-03-11 [1] CRAN (R 4.5.0)
    ##  tidytree            0.4.6   2023-12-12 [1] CRAN (R 4.5.0)
    ##  tidyverse         * 2.0.0   2023-02-22 [1] CRAN (R 4.5.0)
    ##  timechange          0.3.0   2024-01-18 [1] CRAN (R 4.5.0)
    ##  treeio              1.34.0  2025-10-30 [1] Bioconductor 3.22 (R 4.5.1)
    ##  tweenr              2.0.3   2024-02-26 [1] CRAN (R 4.5.0)
    ##  tzdb                0.5.0   2025-03-15 [1] CRAN (R 4.5.0)
    ##  usethis             3.2.1   2025-09-06 [1] CRAN (R 4.5.0)
    ##  vctrs               0.6.5   2023-12-01 [1] CRAN (R 4.5.0)
    ##  viridis           * 0.6.5   2024-01-29 [1] CRAN (R 4.5.0)
    ##  viridisLite       * 0.4.2   2023-05-02 [1] CRAN (R 4.5.0)
    ##  vroom               1.6.7   2025-11-28 [1] CRAN (R 4.5.2)
    ##  withr               3.0.2   2024-10-28 [1] CRAN (R 4.5.0)
    ##  xfun                0.55    2025-12-16 [1] CRAN (R 4.5.2)
    ##  XVector             0.50.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  yaml                2.3.12  2025-12-10 [1] CRAN (R 4.5.2)
    ##  yulab.utils         0.2.3   2025-12-15 [1] CRAN (R 4.5.2)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
