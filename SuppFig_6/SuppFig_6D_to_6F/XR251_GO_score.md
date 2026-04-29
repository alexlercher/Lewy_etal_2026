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

## GET GENES FOR GO of INTEREST

``` r
go_id <- c("GO:0006955","GO:0006954","GO:0034340")

# Get Entrez IDs for the GO term
go_list <- list()
for(i in 1:length(go_id)){
  genes <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = go_id[i],
    keytype = "GOALL",               # include all evidence codes
    columns = c("ENTREZID", "SYMBOL")
  )
  # remove duplicates
  genes <- genes %>% distinct(ENTREZID, SYMBOL)
  # keep unique symbols, remove NA
  genes_symbol <- na.omit(unique(genes$SYMBOL))
  go_list[[i]] <- genes_symbol

go_term <- Term(GOTERM[[go_id[i]]])
go_term <- gsub(" ", "_", go_term)

names(go_list)[i] <- paste0(go_id[i],"_",go_term)

print(paste0(go_id[[i]],": ",gsub("_"," ",go_term)," done! ",length(go_id)-i," GO terms of ",length(go_id)," to go!"))
}
```

    ## 'select()' returned 1:many mapping between keys and columns

    ## [1] "GO:0006955: immune response done! 2 GO terms of 3 to go!"

    ## 'select()' returned 1:many mapping between keys and columns

    ## [1] "GO:0006954: inflammatory response done! 1 GO terms of 3 to go!"

    ## 'select()' returned 1:many mapping between keys and columns

    ## [1] "GO:0034340: response to type I interferon done! 0 GO terms of 3 to go!"

## ADD MODULE SCORES

``` r
for(i in 1:length(go_list)){
  data <- AddModuleScore(
    object = data,
    features = list(go_list[[i]]),
    name = paste0(names(go_list[i]),"_score")
  )
  print(paste0(i," module scores done, ",length(go_id)-i," to go!"))
}
```

    ## Warning: The following features are not present in the object: Cfd, Aicda,
    ## Alox15, Ang3, Ang2, Apoa4, Asgr2, Azgp1, Btn1a1, Cd19, Cd1d2, Cd3d, Cd3e, Cd3g,
    ## Cd8a, Cdh17, Cxcr3, Ccr1l1, Ccr3, Ccr8, Dmbt1, Csf2, Csf3, Cst9, Ctsg, Defb2,
    ## Defa1, Defa29, Defa-rs10, Defa-rs12, Defa-rs2, Defa31, Defa10, Defa11, Defa12,
    ## Defa13, Defa14, Defa15, Defa16, Defa25, Defa3, Defa4, Defa5, Defa6, Defa7,
    ## Defa8, Defa9, Ear2, F2, F2rl1, Fcer1a, Ms4a2, Fga, Fpr-rs3, Fpr-rs4, Gba1,
    ## Gcsam, Gfi1, Lilrb4b, Gpr33, Gzmb, Gzmc, H2-L, H2-M10.1, H2-M2, H2-M9, H2-Oa,
    ## H2-Q10, H2-T3, Ifi44l, Hprt1, Ifna1, Ifna11, Ifna2, Ifna4, Ifna5, Ifna6, Ifna7,
    ## Ifna9, Ifnab, Ifnb1, Ifng, Cd79b, Ighg2b, Ighg1, Igh-J, Igh-VJ558, Igkv1-117,
    ## Iglv1, Il12b, Il13, Il17a, Il18rap, Il1rn, Il3, Il9, Il9r, Ins1, Acod1, Itgb2l,
    ## Ivl, Klrc1, Klrc2, Krt16, Krt6a, Lep, Lgals6, Xcl1, Ltf, Klrb1c, Ncr1, Marco,
    ## Mbl1, Mmp12, Mmp3, Mmp7, Clec4d, Nkx2-3, Nppb, Reg3b, Pck1, Pdcd1, Prf1, Phb1,
    ## Pira2, Lilra6, Ccl21a, Bpifa1, Prg2, Rab17, Raet1a, Raet1b, Raet1c, Rag1, Rag2,
    ## Reg1, Reg2, Reg3a, Reg3g, Trim10, S100a8, S100a9, Apcs, Ccl11, Ccl20, Cxcl15,
    ## Cxcl2, Sftpd, Serpinb9b, Serpinb9f, Serpinb9e, Serpinb9d, Spn, Sprr2a1, Znrf4,
    ## Vamp7, Tcra, Trav5-4, Trgv1, Trgv2, Trgv4, Trgv7, Otop1, Trex1, Umod, Vpreb1b,
    ## Vpreb3, Wap, Defa17, Fyb1, Klk7, Ubd, Ifi202b, Psg19, Sh2d1b1, Pla2g2f, Reg3d,
    ## Elane, Cfhr1, Tnfsf14, Clec4f, Prg3, Ear3, Ear4, Rnase2b, Sit1, Il36a, Mefv,
    ## Cxcl13, Pfpl, Ccl24, Ripk3, Raet1d, Clec4e, Clec4n, Clec7a, Pf4, Psg23, Ppbp,
    ## U90926, Pglyrp2, Klrc3, Il20, Sectm1b, Trem1, Trem3, Il21, Fcamr, Wfdc21,
    ## S100a14, Gkn2, Mptx1, Defa21, Hamp2, Ppp1r14bl, Ceacam11, Ceacam14, Gp2,
    ## Bpifa5, Ceacam12, Defa20, Wfdc15a, H2-Q8, Klk5, Cd177, Trim15, Cd209f, Cd209b,
    ## Fcmr, Il36b, Ceacam13, Clec4b1, Cd209g, Zdhhc11, Trim29, Psg21, Nkg7, Spink5,
    ## Exosc6, Ceacam5, Lrrc15, Defb29, Defb41, Defb12, Spag11a, Klrb1b, Fcrl2, Tlr9,
    ## Gimap3, Sucnr1, Kars1, Il24, Clec2i, Ear6, Ear10, Rnase2a, Serpinb9g, Fgg,
    ## Oas1d, Ceacam15, Slc22a13, Eprs1, Igkv3-1, Igkv3-7, Igkv6-15, Igkv6-20, Pgc,
    ## Fgb, Trbv26, H2-Q9, Iglv2, H2-M10.3, Iglc1, Iglc2, Iglc3, Psg25, Psg28, Il25,
    ## Cd209c, Cd209d, Cd209e, Cd209a, Btnl10, Wfdc12, Wfdc15b, Rfpl4, Igh-VX24,
    ## Ighv9-1, Trav16, Wfdc5, Sectm1a, Nlrp4b, Lypd11, Trav3-4, Ighv1-54, Trim52,
    ## Trav7-6, Ighv1-78, Igkv4-79, Trim50, Il36g, Il1f10, Cd300lb, Cd300ld, Ang4,
    ## Trav13-5, H2-M10.4, H2-M11, H2-M1, H2-M10.5, Cxcr1, Bpifb1, Ifna13, Ifna16,
    ## Ifne, Oas1e, Igkv8-19, Klrb1f, Klrh1, Lypd10, Nlrp9a, Ffar2, Mrgprb1, Trim60,
    ## Igkv4-57, Igkv4-74, Btnl9, Ighv2-3, Ighv14-3, Ighv6-6, Igha, Gapt, Slc30a8,
    ## Colec10, Iigp1c, Tnfsf18, Ifna15, Ifna12, Oas1f, Igkv1-135, Igkv1-132,
    ## Igkv11-125, Igkv9-124, Igkv17-127, Igkv14-100, Igkv4-86, Igkv4-81, Igkv7-33,
    ## Psg22, Nlrp9b, Nlrp4a, Mrgprx2, Mrgprb2, Defb8, Triml1, Apol10a, Kir3dl1,
    ## Defb15, Defb35, Cd207, Defb19, Oas1h, Il17f, Leap2, Trim61, Defb36, Mill1,
    ## Wfdc16, Vsig4, Trav9n-4, Ifnz, Defb20, Klri2, Skint4, Trav7d-6, Trdv2-2,
    ## Fpr-rs6, Fpr-rs7, Tnfsf15, Skint7, Il19, Bpi, Skint9, Cxcl3, Nlrp9c, Ifnl2,
    ## Gsdmc2, Lcn10, Trim75, H2-M10.2, Ifnl3, Defb37, Defb38, Defb40, Ighe, Ighg2a,
    ## Ighg, Ighg3, Ighv10-1, Ighv10-3, Ighv1-64, Ighv1-66, H2-Eb2, Igkv2-112,
    ## Igkv1-110, Igkv16-104, Igkv4-50, Igkv5-43, Igkv6-25, Clec4b2, Igkv4-58,
    ## Ceacam23, Defa39, Defa22, Trav7d-4, Cd300ld3, Ighv1-19, Ighv1-56, Ighv1-52,
    ## Igkv4-92, Igkv12-89, Igkv5-37, Igkv8-30, Igkv4-57-1, Igkv4-68, Pglyrp4,
    ## Trav6d-7, Igkv4-72, Igkv4-70, Igkv4-55, Igkv4-54, Trav7-4, Mir21a, Mir155,
    ## Ifnk, Cd300c, H2-M10.6, Defb21, Ifna14, Iglv3, Wfdc13, Ube2frt, Ighv6-7,
    ## Ighv1-4, Ighv1-59, Ighv8-9, Defb45, Igkv1-122, Igkv9-120, Igkv1-99, Igkv10-95,
    ## Igkv4-91, Igkv4-90, Igkv13-85, Igkv6-32, Igkv8-28, Igkv8-27, Trim34b, Psg20,
    ## Igkv4-71, Ighv1-84, Ighv12-3, Ighv1-67, Fcrlb, Gm13271, Igkv12-98, Gpr31b,
    ## Trav14-1, Trav15d-1-dv6d-1, Trav9-2, Trav15-2-dv6-2, Defb22, Nlrp4e, Defa23,
    ## Defa24, Ang5, Ear14, Ccl26, Ighv5-9, Serpinb9h, Sh2d1b2, Gm5849, Gm13283,
    ## Gm13272, Gm13276, Gm13277, Gm13275, Igkv14-111, Igkv4-80, Igkv4-78, Igkv5-45,
    ## Igkv12-44, Psg27, Igkv4-53, Igkv4-61, Trav13n-3, Trav16d-dv11, Trav14-3,
    ## Trav7-3, Trav15n-2, Trav19, Btnl2, Psg26, Igkv4-51, Ighv1-69, Igkv5-48,
    ## Ighv1-71, Ighv1-72, Igkv12-41, Ighv1-77, Igkv5-39, Igkv12-38, Igkv18-36,
    ## Igkv1-35, Igkv8-34, Igkv6-29, Igkv8-26, Trav6n-6, Igkv8-21, Trbv1, Ighv1-75,
    ## Gm6377, Btnl6, Usp17le, Igkv3-4, Trdv5, Igkv3-2, Apol11a, Defa28, Defa26,
    ## Trav7d-5, Trav9-1, Trav6d-3, Trav17, Igkv1-133, Igkv1-131, Igkv14-130,
    ## Ighv5-12-4, Igkv14-126, Igkv9-123, Igkv2-109, Ighv6-4, Igkv4-69, Igkv4-63,
    ## Ighv1-34, Defb23, Wfdc10, Ighv16-1, Ighv9-2, Ighv7-3, Ighv14-4, Ighv3-4,
    ## Ighv13-2, Ighv6-3, Ighv8-2, Ighv1-11, Ighv1-12, Ighv1-15, Ighv1-16, Ighv1-18,
    ## Ighv1-26, Ighv1-31, Ighv1-36, Ighv1-42, Ighv1-43, Ighv1-47, Ighv8-4, Ighv1-49,
    ## Ighv8-6, Trav7-5, Trav9-4, Ighv2-6, Ang6, Gm12250, Rnf7l, Btnl4, Ighv7-4,
    ## Ighv3-5, Ighv5-16, Trav5-1, Trav7-2, Trav14d-3-dv8, Defa38, Trav7n-4, Ighv9-4,
    ## Ighv8-11, Igkv4-62, Igkv6-23, Trav13-2, Ighv2-5, Ighv6-5, Trav13d-4, Skint1,
    ## Ighv2-9, Trav11d, Igkv8-16, Ighv8-5, Trav13d-1, Ighv5-9-1, Defb26, Defb43,
    ## Defb25, Defb18, Trav12-3, Defa42, C1rb, Igkv17-121, Trav3-1, Trav6-1, Trav8d-1,
    ## Trav6d-5, Trav12d-3, Igkv10-94, Trav14-2, Trav9d-4, Trav5d-4, Trav13d-3,
    ## Trav10n, Trav7n-5, Trav12n-3, Igkv4-59, Trav10, Igkv8-18, Igkv6-17, Igkv6-14,
    ## Igkv6-13, Igkv3-12, Igkv3-10, Igkv3-9, Igkv3-5, Igkv3-3, H2-T26, Ighv2-4,
    ## Ighv5-12, Ighv14-2, Ighv3-3, Ighv1-5, Ighv1-7, Ighv1-9, Ighv1-20, Ighv1-37,
    ## Ighv1-62-3, Ighv1-80, Ighv1-81, Ighv1-85, H60c, Trav7-1, Trav9d-1, Ighv1-61,
    ## Trav8-1, Fcrl6, Igkv8-24, Igkv13-84, Igkv1-88, Igkv19-93, Igkv10-96,
    ## Igkv15-103, Igkv9-129, Igkv2-137, Igkv12-46, Mir301, Mir326, Mir181b-1, Mir324,
    ## Mir181b-2, Igkv20-101-2, Ighv5-2, Ighv2-2, Ighv5-4, Ighv5-6, Ighv5-15, Ighv2-7,
    ## Ighv5-17, Ighv7-1, Ighv7-2, Ighv14-1, Ighv4-1, Ighv3-1, Ighv11-1, Ighv11-2,
    ## Ighv9-3, Ighv3-6, Ighv3-8, Ighv15-2, Ighv1-22, Ighv1-23, Ighv1-24, Ighv1-39,
    ## Ighv1-50, Ighv1-53, Ighv1-55, Ighv8-8, Ighv1-58, Ighv1-63, Ighv8-12, Ighv2-9-1,
    ## Ighv2-6-8, Wfdc17, Trav9d-3, Btnl1, Pira12, Defa30, Gapdhrt2, Trbv3, Pira13,
    ## Defa35, Defa41, Defa40, Defa27, Defa37, Defa34, Trav7d-2, Trav7d-3, Gapdhrt,
    ## Trav6d-6, Trav4d-3, Ccl21b, Trav13n-1, Trav4n-3, Trav5n-4, Cd300ld4, Trav9n-3,
    ## Trav14n-3, Tigit, Trav12-2, Cd200l2, Trav10d, Trav15d-2-dv6d-2, Trdv2-1,
    ## Trav13-4-dv7, Trav6-7-dv9, Trav4-4-dv10, Trav23, Trdv4, Mir873a, Mir511, Trbv2,
    ## Trbv4, Trbv5, Trbv12-1, Trbv13-1, Trbv12-2, Trbv13-2, Trbv13-3, Trbv14, Trbv15,
    ## Trbv16, Trbv17, Trbv19, Trbv20, Trbv21, Trbv23, Trbv29, Trbv30, Trav6n-5,
    ## Trav6n-7, Trav7n-6, Trav8n-2, Trav16n, Trav13n-4, Trav4-2, Trav6-5, Trav6-6,
    ## Trav11, Trav13-1, Trav4-3, Defa2, Zfp683, AY761185, Mptx2, H2-Ea, H2-T27,
    ## Gm17416, Ighv8-13, Ighv1-74, Ighv1-76, Ighv1-82, Pvrig, not searching for
    ## symbol synonyms

    ## [1] "1 module scores done, 2 to go!"

    ## Warning: The following features are not present in the object: Adipoq, Ahsg,
    ## Alox15, Fabp4, Bdkrb1, Chi3l1, Cxcr3, Ccr1l1, Ccr3, Csrp3, Elf3, F2, F2rl1,
    ## Fabp1, Fcer1a, Fpr-rs3, Fpr-rs4, Gba1, Gnat2, Gpr33, Gzma, Gzmb, Gzmc, Hp,
    ## Ifnb1, Ifng, Ighg2b, Ighg1, Il12b, Il13, Il17a, Il18rap, Il1rn, Ins1, Acod1,
    ## Itgb2l, Klk1b1, Krt16, Lep, Xcl1, Cma1, Tpsb2, Mmp8, Orm1, Orm2, Orm3, Reg3b,
    ## Prf1, Pla2g2a, Ccl21a, Proc, Ptgir, Reg3a, Reg3g, S100a8, S100a9, Saa1, Saa2,
    ## Saa4, Apcs, Ccl11, Ccl20, Cxcl15, Cxcl2, Serpina1a, Spn, Vamp7, Tff2, Trex1,
    ## Umod, Scgb1a1, Agr2, Nlrp5, Ifi202b, Pla2g2e, Elane, Naip7, Il36a, Mefv,
    ## Cxcl13, Ccl24, Clec7a, Pf4, Ppbp, Pglyrp2, Gsdma, Gkn2, Mptx1, Duoxa2, Odam,
    ## Il36b, Nkg7, Gsdma2, Nlrp14, Tlr9, Sucnr1, Kars1, Chil4, Lgals2, Il25, Ghsr,
    ## Nlrp4b, Il36g, Il1f10, Hrh4, Nlrp9a, Ffar2, Mrgpra3, Il22ra2, Tnfsf18, Nlrp9b,
    ## Nlrp4a, Il17f, Gsdmc3, H2bc1, Cd200r1l, Fpr-rs6, Fpr-rs7, Cxcl3, Nlrp9c,
    ## Gsdmc2, Lcn10, Ighe, Ighg2a, Ighg3, Mir21a, Mir147, Mir155, Ugt1a1, Spink7,
    ## Kprp, Gpr31b, Nlrp4e, Gsdma3, Ccl26, Gm5849, Gm12250, Trav7-2, Mir301, Mir324,
    ## Msmp, Trav7d-2, Ccl21b, Cd200l2, Mir883b, Mptx2, Mir7116, Mir7578, not
    ## searching for symbol synonyms

    ## [1] "2 module scores done, 1 to go!"

    ## Warning: The following features are not present in the object: Ifna1, Ifna11,
    ## Ifna2, Ifna4, Ifna5, Ifna6, Ifna7, Ifna9, Ifnab, Ifnb1, Mmp12, Trex1, Oas1d,
    ## Ifna13, Ifna16, Ifne, Oas1e, Ifna15, Ifna12, Oas1f, Oas1h, Ifnz, Ifnk, Ifna14,
    ## Gm13271, Gm13283, Gm13272, Gm13276, Gm13277, Gm13275, not searching for symbol
    ## synonyms

    ## [1] "3 module scores done, 0 to go!"

## PLOT MODULE SCORES

``` r
for(i in 1:length(go_list)){
  colname <- paste0(names(go_list[i]), "_score1")
  # sanitize GO name
  go_name <- gsub(":", "_", names(go_list)[i])
  # make output directory
  go_dir <- file.path("output", "GO_score", go_name)
  
  if (!dir.exists(go_dir)) {
    dir.create(go_dir, recursive = TRUE)
  }
  
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
        ylab("go_list Module Score (log1p)") +
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
    write_tsv(paste0(go_dir,"/",go_name,"_module_scores.tsv"))
  
  # plot
  p <- ggplot(df_summary, aes(x = condition, y = cell_type)) +
        geom_point(aes(size = pct_expressing, color = avg_score)) +
        scale_color_gradient(low = "lightgrey", high = "firebrick") +
        scale_size(range = c(1, 5)) +  # adjust point size range
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("Time Point") +
        ylab("Cell Type") +
        ggtitle("GO Module Score")
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
  # min and max are clear but need to fraction range of measurement to get low and middle numbers
  low_signal <- (abs(max(mat_avg) - min(mat_avg)))/9.5 + min(mat_avg)
  middle_signal <- (abs(max(mat_avg) - min(mat_avg)))/6 + min(mat_avg)
  
  # define color range
  col_scale <- colorRamp2(
    c(min(mat_avg, na.rm= TRUE), low_signal, middle_signal, max(mat_avg, na.rm = TRUE)),
    c("#8c6bb1", "#cce6ff", "#fff5b1", "#fb8072")
  )
  
  # plot and save heatmap
  ht <- ComplexHeatmap:: Heatmap(
    mat_avg,
    name = paste0(gsub("_"," ",go_name)," module score"),
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
      ggtitle(nm) +
      theme(aspect.ratio = 1)
  })
  p_combined <- wrap_plots(p_list, ncol =3)
  print(p_combined)

print(paste0(gsub("_"," ",go_name)," done! ",as.character(length(go_id)-i)," of ",length(go_id)," remaining"))
}
```

![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-1.png)<!-- -->![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-2.png)<!-- -->![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-3.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-4.png)<!-- -->

    ## [1] "GO 0006955 immune response done! 2 of 3 remaining"

![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-5.png)<!-- -->![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-6.png)<!-- -->![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-7.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-8.png)<!-- -->

    ## [1] "GO 0006954 inflammatory response done! 1 of 3 remaining"

![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-9.png)<!-- -->![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-10.png)<!-- -->![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-11.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](XR251_GO_score_files/figure-gfm/plot%20module%20scores-12.png)<!-- -->

    ## [1] "GO 0034340 response to type I interferon done! 0 of 3 remaining"

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
