XR244_1_brain_cytokines
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
library(grid)
```

## DATA IMPORT

``` r
# import
experiment <- c("XR244_1")
data <- read_csv(paste0("input/",experiment,"_brain_data_formatted.csv"))
```

    ## Rows: 29 Columns: 20
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (1): Cytokine
    ## dbl (19): PBS_WNV_1, PBS_WNV_2, PBS_WNV_3, PBS_WNV_4, PBS_WNV_5, PIC_WNV_1, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# define cytokines as row names
data <- data %>%
  column_to_rownames(var = "Cytokine")
```

## SELECT SIGNIFICANT CONDITIONS

``` r
# filter for samples that are significant in at least one comparison
# add TRUE/FALSE columns for significance
# if sum of all TRUE/FALSE columns is >1, keep row
data_sig <- data %>%
  mutate(sig_pic_pbs = if_else(pval_PBS_PIC < 0.05,T,F),
         sig_pic_nai = if_else(pval_PIC_NAI < 0.05,T,F),
         sig_pbs_nai = if_else(pval_PBS_NAI < 0.05,T,F)) %>%
  mutate_at(vars(sig_pic_pbs,
                 sig_pic_nai,
                 sig_pbs_nai), ~replace_na(.,0)) %>%
  mutate(sum_pval = rowSums(select(., starts_with("sig_")))) %>%
  filter(sum_pval > 0) %>%
  select(!sum_pval)
```

## PLOT HEATMAP with ROW Z SCORES

``` r
# define formula to calculate row Z score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
# select individual values (all replicates) to plot
data_sig_values_only <- data_sig[1:13]

# select mean values (per condition) to plot
data_sig_values_only_mean <- data_sig[14:16]

# combine in list
data_sig_list <- list("individual values" = data_sig_values_only, 
                          "means" = data_sig_values_only_mean)

# calculate row Z score for each row
data_norm <- list()
for(i in 1:length(data_sig_list)){
  n <- t(apply(data_sig_list[[i]], 1, cal_z_score))
  data_norm[[names(data_sig_list)[[i]]]] <- n
}

# annotate sample columns for individual samples
my_sample_col_individual <- data.frame(sample = c(rep("pbs_wnv",5),
                                       rep("pic_wnv",5),
                                       rep("nai",3)
                                       ))
row.names(my_sample_col_individual) <- colnames(data_sig_values_only)

# annotate sample columns for means
my_sample_col_mean <- data.frame(sample = c("pbs_wnv",
                                            "pic_wnv",
                                            "nai")
)
row.names(my_sample_col_mean) <- colnames(data_sig_values_only_mean)

my_sample_col <- list("individual values" = my_sample_col_individual,
                      "means" = my_sample_col_mean)

# calculate dendrogram
hc <- lapply(X = data_norm, FUN = function(x){
  hclust(dist(x), method = "complete")
})

# add cluster information to heatmap
numberofclusters = 5
my_cytokine_col <- list()
for(i in 1:length(hc)){
  a <- cutree(hc[[i]], k = numberofclusters)
  b <- data.frame(cluster = ifelse(test = a == 1, yes = "1",
                                   ifelse(test = a == 2, yes = "2",
                                          ifelse(test = a == 3, yes = "3",
                                                ifelse(test = a == 4, yes = "4", no = "5")))))
  my_cytokine_col[[names(hc)[[i]]]] <- b
}

# add column cytokine to dataframes
for(i in 1:length(my_cytokine_col)){
  my_cytokine_col[[i]]$cytokine <- rownames(my_cytokine_col[[i]])
}

# info about significant p-value and condition to heatmap annotation
data_sig_info <- data_sig %>%
  select(sig_pic_pbs,
         sig_pic_nai,
         sig_pbs_nai)
data_sig_info$cytokine <- rownames(data_sig_info)

# merge information about significance to cluster
for(i in 1:length(my_cytokine_col)){
  my_cytokine_col[[i]] <- left_join(my_cytokine_col[[i]],data_sig_info)
  my_cytokine_col[[i]]$sig_pic_pbs <- as.integer(my_cytokine_col[[i]]$sig_pic_pbs)
  my_cytokine_col[[i]]$sig_pic_nai <- as.integer(my_cytokine_col[[i]]$sig_pic_nai)
  my_cytokine_col[[i]]$sig_pbs_nai <- as.integer(my_cytokine_col[[i]]$sig_pbs_nai)
  rownames(my_cytokine_col[[i]]) <- my_cytokine_col[[i]]$cytokine
  my_cytokine_col[[i]]$cytokine <- NULL
}
```

    ## Joining with `by = join_by(cytokine)`
    ## Joining with `by = join_by(cytokine)`

``` r
# generate and save heatmap
for(i in 1:length(data_norm)){
  p <- pheatmap::pheatmap(data_norm[[i]], annotation_col = my_sample_col[[i]], annotation_row = my_cytokine_col[[i]],
           show_rownames = T, main = paste0(experiment," - ","cytokine profile of ",names(data_norm[i])),
           #color = colorRampPalette(c("midnightblue", "ivory", "firebrick3"))(50),
           cutree_rows = numberofclusters)
grid::grid.newpage()
grid::grid.draw(p$gtable)
}
```

![](XR244_1_brain_cytokines_files/figure-gfm/plot%20heatmap%20with%20row%20z%20scores-1.png)<!-- -->![](XR244_1_brain_cytokines_files/figure-gfm/plot%20heatmap%20with%20row%20z%20scores-2.png)<!-- -->

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
    ##  knitr          1.51    2025-12-20 [1] CRAN (R 4.5.2)
    ##  lifecycle      1.0.4   2023-11-07 [1] CRAN (R 4.5.0)
    ##  lubridate    * 1.9.4   2024-12-08 [1] CRAN (R 4.5.0)
    ##  magrittr       2.0.4   2025-09-12 [1] CRAN (R 4.5.0)
    ##  memoise        2.0.1   2021-11-26 [1] CRAN (R 4.5.0)
    ##  otel           0.2.0   2025-08-29 [1] CRAN (R 4.5.0)
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
