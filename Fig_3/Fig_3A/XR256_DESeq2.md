XR256_DESeq2
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

## DATA IMPORT AND CLEANUP

``` r
# import table with reads of all samples
data <- read_csv("input/raw_counts.csv")
```

    ## Rows: 41079 Columns: 18
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (1): Gene.name
    ## dbl (17): Length, LPS.1, LPS.2, LPS.3, LPS.4, ODN.1, ODN.2, ODN.3, ODN.4, PB...
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

# reorder columns to PBS, LPS, ODN and PIC
data_clean <- data_clean %>%
  dplyr::select(gene,
                starts_with("PBS"),
                starts_with("LPS"),
                starts_with("ODN"),
                starts_with("PIC"))

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
data_raw <- data_clean %>%
  dplyr::filter(gene != "a")

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

``` r
# save DESeq2 object for separate analyes
saveRDS(dds, file = "output/XR256_dds_processed.rds")
```

## FILTER OUT GENES WITH LESS THAN THRESHOLD COUNT ACROSS ALL CONDITIONS

``` r
threshold_count <- 10
keep <- rowSums(counts(dds)) >= threshold_count
head(keep)
```

    ## 0610005C13Rik 0610006L08Rik 0610009E02Rik 0610009L18Rik 0610010K14Rik 
    ##          TRUE          TRUE          TRUE          TRUE          TRUE 
    ## 0610025J13Rik 
    ##          TRUE

``` r
dds <- dds[keep,]
rm(threshold_count, keep)
```

## EXAMINE DESEQ2 RESULTS

``` r
# calculated correction factors and dispersion plot
as.data.frame(sizeFactors(dds))
```

    ##       sizeFactors(dds)
    ## PBS.1        1.1234538
    ## PBS.2        1.1148898
    ## PBS.3        1.0335702
    ## PBS.4        1.0186495
    ## LPS.1        1.0005899
    ## LPS.2        1.0236272
    ## LPS.3        1.0614183
    ## LPS.4        1.0766814
    ## ODN.1        0.5127088
    ## ODN.2        0.9654827
    ## ODN.3        1.0082904
    ## ODN.4        1.0864279
    ## PIC.1        0.8978934
    ## PIC.2        0.9578021
    ## PIC.3        1.1930349
    ## PIC.4        1.2269818

``` r
plotDispEsts(dds,
             genecol = "darkgray",
             fitcol = "black",
             finalcol = "tomato3")
```

![](XR256_DESeq2_files/figure-gfm/examine%20DESeq2%20results-1.png)<!-- -->

``` r
#--------------------------------------------------------------------
# SAMPLE DISTANCE PLOTS
#--------------------------------------------------------------------
# do sample distance plot for all samples
rld <- rlog(dds)
sampleDist <- dist(t(assay(rld)))
as.matrix(sampleDist)[1:3,1:3]
```

    ##          PBS.1    PBS.2    PBS.3
    ## PBS.1  0.00000 10.81160 10.45822
    ## PBS.2 10.81160  0.00000 10.72226
    ## PBS.3 10.45822 10.72226  0.00000

``` r
sampleDistMatrix <- as.matrix(sampleDist)
rownames(sampleDistMatrix) <- paste(rld$description, rld$replicate, sep="_")
colnames(sampleDistMatrix) <- paste(rld$description, rld$replicate, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Purples")))

heatmap.2(sampleDistMatrix, trace = "none", col = colors, cexRow = 0.5, cexCol = 0.5, main = "Sample to sample distances")
```

![](XR256_DESeq2_files/figure-gfm/examine%20DESeq2%20results-2.png)<!-- -->

``` r
# clean up
rm(rld, sampleDist,sampleDistMatrix)

#--------------------------------------------------------------------
# PCA PLOTS
#--------------------------------------------------------------------
# do PCA plot for all samples
vsdata <- vst(dds, blind = F)
p <- plotPCA(vsdata, intgroup = "description") + 
    scale_color_brewer(palette = "Set1") + 
    theme_bw() + 
    ggtitle("PCA all samples") +
    theme(aspect.ratio = 1) + 
    geom_text(aes(label=name))
```

    ## using ntop=500 top features by variance

``` r
print(p)
```

![](XR256_DESeq2_files/figure-gfm/examine%20DESeq2%20results-3.png)<!-- -->

## GET DESeq2 RESULTS FOR COMPARISONS OF INTEREST

``` r
# generate function to generate output files for comparisons of interest
deseq_all_comparisons <- function(deseq_data,pval_cutoff) {
  datalist = list()
  
  for (i in 1:nrow(DE_groups)){
    numerator <- DE_groups[[i,1]]
    denominator <- DE_groups[[i,2]]
    
    #Get results from specific contrasts
    results_contrast <- results(deseq_data, contrast = c("description", numerator, denominator))
    results_contrast_wo_na=results_contrast[!is.na(results_contrast$pvalue),]
    results_contrast_wo_na$gene_id=rownames(results_contrast_wo_na)
    results_contrast_sign=results_contrast_wo_na[results_contrast_wo_na$pvalue<=pval_cutoff,]
    results_contrast_sign$comparison=i
    datalist[[i]]=results_contrast_sign
  }
  big_data = do.call(rbind, datalist)
  return(big_data)
}

# generate file with comparisons of interest
pval_cutoff = 1
results <- deseq_all_comparisons(dds,pval_cutoff)

head(results(dds, tidy=TRUE))
```

    ##             row    baseMean log2FoldChange      lfcSE       stat    pvalue
    ## 1 0610005C13Rik  13.0274663     0.74249251 0.54721410  1.3568592 0.1748259
    ## 2 0610006L08Rik   0.7052578    -2.09316700 2.88661700 -0.7251281 0.4683735
    ## 3 0610009E02Rik  12.2613953    -0.49246239 0.57583245 -0.8552182 0.3924304
    ## 4 0610009L18Rik  43.9845264    -0.04704137 0.35322418 -0.1331771 0.8940533
    ## 5 0610010K14Rik 764.6751201    -0.08401000 0.08382175 -1.0022458 0.3162249
    ## 6 0610025J13Rik   1.7695448    -1.10150226 1.59009701 -0.6927265 0.4884812
    ##        padj
    ## 1 0.7685347
    ## 2        NA
    ## 3 0.9133184
    ## 4 0.9948922
    ## 5 0.8800925
    ## 6        NA

``` r
summary(results)
```

    ## 
    ## out of 155862 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 5371, 3.4%
    ## LFC < 0 (down)     : 4090, 2.6%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 41803, 27%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
# fact check that number of comparisons in DESeq2 output is same as comparisons of interest in input file
length(unique(data.frame(results)$comparison)) == nrow(DE_groups)
```

    ## [1] TRUE

``` r
# generate table with comparison names
mycomparison_names <- DE_groups %>%
  dplyr::select(ComparisonNumber, Comparison) %>%
  dplyr::rename(comparison = ComparisonNumber,
                comparison_name = Comparison) %>%
  print()
```

    ##   comparison comparison_name
    ## 1          1      LPS_vs_PBS
    ## 2          2      ODN_vs_PBS
    ## 3          3      PIC_vs_PBS
    ## 4          4      PIC_vs_LPS
    ## 5          5      PIC_vs_ODN
    ## 6          6      LPS_vs_ODN

## MA PLOTS FOR COMPARISONS

``` r
# do MA plots for all comparisons
for(i in 1:length(unique(results$comparison))){
  selected_res <- subset(results, comparison == i)
  plotMA(selected_res, colSig = "tomato3", 
         main = paste0("MA plot of condition ",i,": ",mycomparison_names$comparison_name[[i]]))
}
```

![](XR256_DESeq2_files/figure-gfm/MA%20plots%20for%20comparisons-1.png)<!-- -->![](XR256_DESeq2_files/figure-gfm/MA%20plots%20for%20comparisons-2.png)<!-- -->![](XR256_DESeq2_files/figure-gfm/MA%20plots%20for%20comparisons-3.png)<!-- -->![](XR256_DESeq2_files/figure-gfm/MA%20plots%20for%20comparisons-4.png)<!-- -->![](XR256_DESeq2_files/figure-gfm/MA%20plots%20for%20comparisons-5.png)<!-- -->![](XR256_DESeq2_files/figure-gfm/MA%20plots%20for%20comparisons-6.png)<!-- -->

## WRITE DESeq2 RESULTS TO FILE

``` r
# filter for genes with p-value cutoff of choice
padj_cutoff = 1
results_cutoff <- subset(results, results$padj <= padj_cutoff)

# save these results as data frame and add information on comparisons
df_results <- data.frame(results_cutoff)
df_results <-  merge(df_results,mycomparison_names, by=c("comparison"))

# save each comparison as separate data frames
comparisons_list <- list()
for(i in 1:length(unique(df_results$comparison))){
  x <- subset(df_results,df_results$comparison == i)
  y <- filter(mycomparison_names, comparison == i)
  y <- y$comparison_name
  x$comparison_name <- NULL
  x$comparison_name <- as.character(y)
  comparisons_list[[paste0(i,"_",y)]] <- x
  rm(x,y)
}

for(i in 1:length(comparisons_list)){
  write_tsv(comparisons_list[[i]], paste0("output/",names(comparisons_list)[[i]],"_DEG.tsv"))
}
```

## GET FPM/CPM, CLEAN UP TABLE AND WRITE TO FILE

``` r
# get fpm/cpm from dds
data_cpm <- fpm(dds)
data_cpm <- as.data.frame(data_cpm)
data_cpm$gene <- rownames(data_cpm)
data_cpm  <- data_cpm %>%
  dplyr::select(gene, everything())

# import new col names for data
colnames_data <- read.delim("input/col_names_data.tsv")

# check whether old col name order is correct
colnames(data_cpm) == colnames_data$colnames_old
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [16] TRUE TRUE

``` r
# rename data with new colnames
colnames(data_cpm) <- colnames_data$colnames_new

# order columns
# order by sample number suffix 1,2,3,4 as "1$"
# order by sample prefix PBS,.. as "^PBS"
data_cpm <- data_cpm %>%
  dplyr::select(
    gene,
    matches("1$"),
    matches("2$"),
    matches("3$"),
    matches("4$")) %>% 
  dplyr::select(
    gene,
    matches("^PBS"),
    matches("^LPS"),
    matches("^ODN"),
    matches("PIC")
  )

# save cpm to file
write_tsv(data_cpm, "output/XR256_CPM.tsv")
```

## GET FPKM FROM DDS, CLEAN UP TABLE AND WRITE TO FILE

``` r
# get gene order from dds object
dds_gene_order <- data.frame(SYMBOL = rownames(dds))

# get transcript lengths (actual transcript, not genomic length)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
```

    ## Loading required package: GenomicFeatures

    ## Loading required package: AnnotationDbi

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
transcriptLengths <- transcriptLengths(TxDb.Mmusculus.UCSC.mm10.knownGene)
# rename gene_id column to ENTREZID
transcriptLengths <- transcriptLengths %>%
  dplyr::rename(ENTREZID = gene_id)

# get ENTREZID and SYMBOLS from all genes
library(org.Mm.eg.db)
```

    ## 

``` r
mm <- org.Mm.eg.db
my.geneID <- unique(transcriptLengths$ENTREZID)
gene_list_geneID <- AnnotationDbi::select(mm, 
                                          keys = my.geneID,
                                          columns = c("ENTREZID", "SYMBOL"),
                                          keytype = "ENTREZID")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
# merge all transcript lengths to all genes
transcriptLengths <- left_join(gene_list_geneID, transcriptLengths)
```

    ## Joining with `by = join_by(ENTREZID)`

``` r
# only keep longest transcript in case some genes have multiple transcripts
longestTranscripts <- transcriptLengths %>%
  na.omit() %>%
  group_by(ENTREZID) %>%
  arrange(desc(tx_len)) %>%
  slice_head()

# merge transcript lengths with genes from dds object (order will be kept)
dds_gene_order_length <- left_join(dds_gene_order,longestTranscripts)
```

    ## Joining with `by = join_by(SYMBOL)`

``` r
# check whether order of SYMBOLS of new table is the same as original DDS object
summary(dds_gene_order$SYMBOL == dds_gene_order_length$SYMBOL)
```

    ##    Mode    TRUE 
    ## logical   25980

``` r
# define basepairs vector
basepairs <- as.vector(dds_gene_order_length$tx_len)

# add basepairs vector to DDS object
mcols(dds)$basepairs <- basepairs

# calculate FPKMs and omit NA (genes with no annotated transcript length)
data_fpkm <- fpkm(dds)
data_fpkm <- na.omit(data_fpkm)
data_fpkm <- as.data.frame(data_fpkm)
data_fpkm$gene <- rownames(data_fpkm)
data_fpkm  <- data_fpkm %>%
  dplyr::select(gene, everything())

# import new col names for data
colnames_data <- read.delim("input/col_names_data.tsv")

# check whether old col name order is correct
colnames(data_fpkm) == colnames_data$colnames_old
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [16] TRUE TRUE

``` r
# rename data with new colnames
colnames(data_fpkm) <- colnames_data$colnames_new

# order columns
# order by sample number suffix 1,2,3,4 as "1$"
# order by sample prefix PBS,.. as "^PBS"
data_fpkm <- data_fpkm %>%
  dplyr::select(
    gene,
    matches("1$"),
    matches("2$"),
    matches("3$"),
    matches("4$")) %>% 
  dplyr::select(
    gene,
    matches("^PBS"),
    matches("^LPS"),
    matches("^ODN"),
    matches("PIC")
  )

# save FPKM to file
write_tsv(data_fpkm, "output/XR256_FPKM.tsv")
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
    ##  package                            * version   date (UTC) lib source
    ##  abind                                1.4-8     2024-09-12 [1] CRAN (R 4.5.0)
    ##  AnnotationDbi                      * 1.72.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  Biobase                            * 2.70.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocGenerics                       * 0.56.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocIO                               1.20.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  BiocParallel                         1.44.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  Biostrings                           2.78.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  bit                                  4.6.0     2025-03-06 [1] CRAN (R 4.5.0)
    ##  bit64                                4.6.0-1   2025-01-16 [1] CRAN (R 4.5.0)
    ##  bitops                               1.0-9     2024-10-03 [1] CRAN (R 4.5.0)
    ##  blob                                 1.2.4     2023-03-17 [1] CRAN (R 4.5.0)
    ##  cachem                               1.1.0     2024-05-16 [1] CRAN (R 4.5.0)
    ##  caTools                              1.18.3    2024-09-04 [1] CRAN (R 4.5.0)
    ##  cigarillo                            1.0.0     2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  cli                                  3.6.5     2025-04-23 [1] CRAN (R 4.5.0)
    ##  codetools                            0.2-20    2024-03-31 [1] CRAN (R 4.5.2)
    ##  crayon                               1.5.3     2024-06-20 [1] CRAN (R 4.5.0)
    ##  curl                                 7.0.0     2025-08-19 [1] CRAN (R 4.5.0)
    ##  DBI                                  1.2.3     2024-06-02 [1] CRAN (R 4.5.0)
    ##  DelayedArray                         0.36.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  DESeq2                             * 1.50.2    2025-11-17 [1] Bioconductor 3.22 (R 4.5.2)
    ##  devtools                             2.4.6     2025-10-03 [1] CRAN (R 4.5.0)
    ##  digest                               0.6.39    2025-11-19 [1] CRAN (R 4.5.2)
    ##  dplyr                              * 1.1.4     2023-11-17 [1] CRAN (R 4.5.0)
    ##  ellipsis                             0.3.2     2021-04-29 [1] CRAN (R 4.5.0)
    ##  evaluate                             1.0.5     2025-08-27 [1] CRAN (R 4.5.0)
    ##  farver                               2.1.2     2024-05-13 [1] CRAN (R 4.5.0)
    ##  fastmap                              1.2.0     2024-05-15 [1] CRAN (R 4.5.0)
    ##  fs                                   1.6.6     2025-04-12 [1] CRAN (R 4.5.0)
    ##  generics                           * 0.1.4     2025-05-09 [1] CRAN (R 4.5.0)
    ##  GenomicAlignments                    1.46.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  GenomicFeatures                    * 1.62.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  GenomicRanges                      * 1.62.1    2025-12-08 [1] Bioconductor 3.22 (R 4.5.2)
    ##  ggplot2                            * 4.0.1     2025-11-14 [1] CRAN (R 4.5.2)
    ##  glue                                 1.8.0     2024-09-30 [1] CRAN (R 4.5.0)
    ##  gplots                             * 3.3.0     2025-11-30 [1] CRAN (R 4.5.2)
    ##  gtable                               0.3.6     2024-10-25 [1] CRAN (R 4.5.0)
    ##  gtools                               3.9.5     2023-11-20 [1] CRAN (R 4.5.0)
    ##  hms                                  1.1.4     2025-10-17 [1] CRAN (R 4.5.0)
    ##  htmltools                            0.5.9     2025-12-04 [1] CRAN (R 4.5.2)
    ##  httr                                 1.4.7     2023-08-15 [1] CRAN (R 4.5.0)
    ##  IRanges                            * 2.44.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  KEGGREST                             1.50.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  KernSmooth                           2.23-26   2025-01-01 [1] CRAN (R 4.5.2)
    ##  knitr                                1.51      2025-12-20 [1] CRAN (R 4.5.2)
    ##  labeling                             0.4.3     2023-08-29 [1] CRAN (R 4.5.0)
    ##  lattice                              0.22-7    2025-04-02 [1] CRAN (R 4.5.2)
    ##  lifecycle                            1.0.4     2023-11-07 [1] CRAN (R 4.5.0)
    ##  locfit                               1.5-9.12  2025-03-05 [1] CRAN (R 4.5.0)
    ##  magrittr                             2.0.4     2025-09-12 [1] CRAN (R 4.5.0)
    ##  Matrix                               1.7-4     2025-08-28 [1] CRAN (R 4.5.0)
    ##  MatrixGenerics                     * 1.22.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  matrixStats                        * 1.5.0     2025-01-07 [1] CRAN (R 4.5.0)
    ##  memoise                              2.0.1     2021-11-26 [1] CRAN (R 4.5.0)
    ##  org.Mm.eg.db                       * 3.22.0    2025-12-09 [1] Bioconductor
    ##  otel                                 0.2.0     2025-08-29 [1] CRAN (R 4.5.0)
    ##  pillar                               1.11.1    2025-09-17 [1] CRAN (R 4.5.0)
    ##  pkgbuild                             1.4.8     2025-05-26 [1] CRAN (R 4.5.0)
    ##  pkgconfig                            2.0.3     2019-09-22 [1] CRAN (R 4.5.0)
    ##  pkgload                              1.4.1     2025-09-23 [1] CRAN (R 4.5.0)
    ##  png                                  0.1-8     2022-11-29 [1] CRAN (R 4.5.0)
    ##  purrr                                1.2.0     2025-11-04 [1] CRAN (R 4.5.0)
    ##  R6                                   2.6.1     2025-02-15 [1] CRAN (R 4.5.0)
    ##  RColorBrewer                       * 1.1-3     2022-04-03 [1] CRAN (R 4.5.0)
    ##  Rcpp                                 1.1.0     2025-07-02 [1] CRAN (R 4.5.0)
    ##  RCurl                                1.98-1.17 2025-03-22 [1] CRAN (R 4.5.0)
    ##  readr                              * 2.1.6     2025-11-14 [1] CRAN (R 4.5.2)
    ##  remotes                              2.5.0     2024-03-17 [1] CRAN (R 4.5.0)
    ##  restfulr                             0.0.16    2025-06-27 [1] CRAN (R 4.5.0)
    ##  rjson                                0.2.23    2024-09-16 [1] CRAN (R 4.5.0)
    ##  rlang                                1.1.6     2025-04-11 [1] CRAN (R 4.5.0)
    ##  rmarkdown                            2.30      2025-09-28 [1] CRAN (R 4.5.0)
    ##  Rsamtools                            2.26.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  RSQLite                              2.4.5     2025-11-30 [1] CRAN (R 4.5.2)
    ##  rstudioapi                           0.17.1    2024-10-22 [1] CRAN (R 4.5.0)
    ##  rtracklayer                          1.70.1    2025-12-20 [1] https://bioc-release.r-universe.dev (R 4.5.2)
    ##  S4Arrays                             1.10.1    2025-12-01 [1] Bioconductor 3.22 (R 4.5.0)
    ##  S4Vectors                          * 0.48.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  S7                                   0.2.1     2025-11-14 [1] CRAN (R 4.5.2)
    ##  scales                               1.4.0     2025-04-24 [1] CRAN (R 4.5.0)
    ##  Seqinfo                            * 1.0.0     2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  sessioninfo                          1.2.3     2025-02-05 [1] CRAN (R 4.5.0)
    ##  SparseArray                          1.10.8    2025-12-18 [1] Bioconductor 3.22 (R 4.5.2)
    ##  stringi                              1.8.7     2025-03-27 [1] CRAN (R 4.5.0)
    ##  stringr                            * 1.6.0     2025-11-04 [1] CRAN (R 4.5.0)
    ##  SummarizedExperiment               * 1.40.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  tibble                               3.3.0     2025-06-08 [1] CRAN (R 4.5.0)
    ##  tidyr                              * 1.3.2     2025-12-19 [1] CRAN (R 4.5.2)
    ##  tidyselect                           1.2.1     2024-03-11 [1] CRAN (R 4.5.0)
    ##  TxDb.Mmusculus.UCSC.mm10.knownGene * 3.10.0    2025-06-06 [1] Bioconductor
    ##  tzdb                                 0.5.0     2025-03-15 [1] CRAN (R 4.5.0)
    ##  usethis                              3.2.1     2025-09-06 [1] CRAN (R 4.5.0)
    ##  vctrs                                0.6.5     2023-12-01 [1] CRAN (R 4.5.0)
    ##  vroom                                1.6.7     2025-11-28 [1] CRAN (R 4.5.2)
    ##  withr                                3.0.2     2024-10-28 [1] CRAN (R 4.5.0)
    ##  xfun                                 0.55      2025-12-16 [1] CRAN (R 4.5.2)
    ##  XML                                  3.99-0.20 2025-11-08 [1] CRAN (R 4.5.0)
    ##  XVector                              0.50.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
    ##  yaml                                 2.3.12    2025-12-10 [1] CRAN (R 4.5.2)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
