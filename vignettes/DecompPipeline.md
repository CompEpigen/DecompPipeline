---
title: 'DecompPipeline: Preprocessing of DNA Methylation data for MeDeCom'
author: "Michael Scherer, Pavlo Lutsik"
date: '2019-09-27'
output:
  html_document:
    fig_height: 5
    fig_width: 5
    keep_md: yes
    mathjax: default
    number_sections: no
    toc: yes
  pdf_document:
    toc: yes
bibliography: biblio.bib
vignette: >
  %\VignetteIndexEntry{MeDeCom}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# Introduction

*DecompPipeline* is an R package created for preprocessing DNA Methylation data for deconvolution of complex tissue samples using one of the methods [MeDeCom](http://public.genetik.uni-sb.de/medecom/), [RefFreeEWAS](https://cran.r-project.org/web/packages/RefFreeEWAS/index.html) or [EDec](https://github.com/BRL-BCM/EDec). Briefly, non-negative matrix factorization of complex methylation data sets is performed to discover latent methylation components (LMCs). Those components can, for instance, represent cell types, but are not limited to that. The *DecompPipeline* uses the [RnBeads](https://rnbeads.org) package for handling input DNA methylation data and is able to handle both BeadArray (27k, 450k, EPIC) and bisulfite sequencing (RRBS, WGBS) data. Necessary steps include stringent filtering and selecting the correct subset of CpG sites for downstream MeDeCom analysis. Also check out [FactorViz](https://github.com/lutsik/FactorViz) for the visualization of deconvolution results.

# Installation

You can install the *DecompPipeline* through GitHub using *devtools*:


```r
install.packages("devtools",repos="https://cran.uni-muenster.de/")
devtools::install_github("lutsik/DecompPipeline")
```

# Using DecompPipeline
## Default Pipeline

The main workhorse for *DecompPipeline* is ```start_decomp_pipeline```. It only requires few inputs, including DNA methylation data either in the form of a *matrix\data.frame* or as an *RnBSet* object. We will now discuss the rich functionalities that are available for having and *RnBSet* as input. A short introduction on how to import DNA methylation data using *RnBeads* is given in the last section. For details on the many options available within the *DecompPipeline*, see the function's documentation. We will further explain the options in this vignette.



## CpG Filtering

This filtering step involves removing potentially unreliable and/or problematic CpGs from further analysis and has dedicated functions both for bisulfite sequencing and array based data sets. We will discuss the two data types separately:

### Array based data sets

Filtering CpG sites of array based data sets (27k, 450k, EPIC) involves setting the minimum number of required beads on the chip (```MIN_N_BEADS```). Furthermore, low and high intensity outliers can be removed by a quantile approach, which removes the highest (```MAX_INT_QUANT```) and lowest quantile (```MIN_INT_QUANT```). In addition, all sites containing any missing value (```FILTER_NA```), outside of CpG context (```FILTER_CONTEXT```), mapping to an annotated Single Nucleotide Polymorphism (```FILTER_SNP```, ```snp.list```) and on the sex chromosomes (```FILTER_SOMATIC```) can be omitted. Further options are available and described in the function's documentation. The function also provides options to normalize data using the methods available in the *RnBeads* R package and will return the processed data set and further information on the steps executed.


```r
data("small.RnBeadSet")
data.prep <- prepare_data(RNB_SET = rnb.set.example,
                          NORMALIZATION = "wm.dasen",
                          MIN_N_BEADS = 5,
                          MIN_INT_QUANT = 0.05,
                          MAX_INT_QUANT = 0.95,
                          FILTER_NA = T,
                          FILTER_SNP = T,
                          FILTER_CONTEXT = FALSE,
                          FILTER_SOMATIC = FALSE)
names(data.prep)
```

```
## [1] "quality.filter"   "annot.filter"     "total.filter"    
## [4] "rnb.set.filtered" "info"
```

### Bisulfite sequencing based data sets

For bisulfite sequencing data sets, different filtering criteria apply. First, a absolute coverage threshold can be specified with ```MIN_COVERAGE``` to remove all sites with lower coverage. Similar to array-based data sets, upper and lower quantile of coverage can be omitted using ```MIN_COVG_QUANT``` and ```MAX_COVG_QUANT```. In complete accordance with array-based data sets, sites having missing values, located at annotated SNPs and on sex chromosomes can be removed.


```r
rnb.set <- load.rnb.set(system.file("extdata/small_rnbSet.zip",package="DecompPipeline"))
data.prep.bs <- prepare_data_BS(RNB_SET = rnb.set,
                                MIN_COVERAGE = 5,
                                MIN_COVG_QUANT = 0.1,
                                MAX_COVG_QUANT = 0.9,
                                FILTER_NA = T,
                                FILTER_SOMATIC = F,
                                FILTER_SNP = F)
```

```
## opening ff /tmp/Rtmpw8TrQZ/ff6f623a453210.ff
```

```
## opening ff /tmp/Rtmpw8TrQZ/ff6f627581c0dd.ff
```

```r
names(data.prep.bs)
```

```
## [1] "quality.filter"   "annot.filter"     "total.filter"    
## [4] "rnb.set.filtered"
```

## Selecting subsets of CpGs

Since performing MeDeCom on complete 450k/EPIC or BS datasets is still computationally infeasible, it is crucial to select sites for subsequent analysis that might define LMCs. The *DecompPipeline* provides multiple options to preselect those CpGs and we will briefly introduce each of them:

* **pheno** This option selects the markers that define sample identity for the grouping information given through *rnb.sample.groups*. Briefly, the limma method is used to define differentially methylated sites between the two groups and uses those markers for subsequent analysis.
* **houseman2012** This option selects 50,000 sites that were found to be cell-type specific using the  @houseman_refbased method on blood cell types. The Houseman method was employed on the @reiniusRef dataset and is thus only applicable to blood data sets.
* **houseman2014** Here, sites are selected as cell-type specific by [RefFreeEWAS](https://cran.rstudio.com/web/packages/RefFreeEWAS/index.html). For further information, see @Houseman2014. This method is applicable to any kind of data set.
* **jaffe2014** Another list of supposedly cell-type specific CpG sites reported in @Jaffe2014.
* **rowFstat** If reference methylation profiles are provided through ```REF_DATA_SET``` and ```REF_PHENO_COLUMN```, sites are selected as those being linked to the reference cell types using an F-test.
* **random** Sites are randomly selected from all possible sites.
* **pca** This option selects the sites that have most influence on the principal components. The number of principal components calculated is determined by ```N_PRIN_COMP```.
* **var** The most variable sites across all samples are selected and used for subsequent analysis (DEFAULT option).
* **hybrid** Selects half of the sites randomly and the other half as the most variable.
* **range** This options selects the sites that have a difference between the minimum and maximum value across all samples higher than ```RANGE_DIFF```.
* **custom** The sites to be used are provided by the user with a file containing row indices. The file needs to be provided in ```CUSTOM_MARKER_FILE```.
* **all** Using all sites available in the input.
* **pcadapt** Uses principal component analysis as implemented in the *bigstats* R package to determine sites that are significantly linked to the potential cell types. This requires specifying K a priori (argument ```K.prior```). We thank Florian Prive and Sophie Achard for providing the idea and parts of the codes.
* **edec_stage0** Employs EDec's stage 0 to infer cell-type specific markers. By default EDec's example reference data is provided.

For most of the options (except for **houseman2012**, **jaffe2014**, and **range**) the number of selected sites can be specified using the parameter ```N_MARKERS```. In contrast to CpG filtering, subset selection is independent of the data type (array-based and BS). The function returns a list, with each entry containing row indices of the selected sites:


```r
cg_subsets <- prepare_CG_subsets(rnb.set=data.prep$rnb.set.filtered,
                                 MARKER_SELECTION = c("houseman2012","var"),
                                 N_MARKERS = 500
)
lengths(cg_subsets)
```

```
## houseman2012          var 
##          110          500
```

## Starting MeDeCom

After these preprocessing steps, you are ready to perfom the actual MeDeCom analysis using the ```start_medecom_analysis``` function. To store output in a format that is later on readable by FactorViz, you need to set the flag ```factorviz.outputs```. Further parameters are described in detail in the reference manual.


```r
md.res <- start_medecom_analysis(rnb.set=data.prep$rnb.set.filtered,
                                 cg_groups = cg_subsets,
                                 Ks=2:5,
                                 LAMBDA_GRID = c(0.01,0.001),
                                 factorviz.outputs = T)
```

```
## [1] "Did not write the variable dump: should only be executed from an environment with all the variables set"
## [2019-09-27 14:42:04, Main:] checking inputs
## [2019-09-27 14:42:04, Main:] preparing data
## [2019-09-27 14:42:04, Main:] preparing jobs
## [2019-09-27 14:42:04, Main:] 176 factorization runs in total
## [2019-09-27 14:48:59, Main:] 100 runs complete
## [2019-09-27 15:01:21, Main:] finished all jobs. Creating the object
```

## Executing DecompPipeline

You can also peform all the steps above, by just calling a single function:


```r
md.res <- start_decomp_pipeline(rnb.set=rnb.set,
                                Ks=2:5,
                                lambda.grid = c(0.01,0.001),
                                factorviz.outputs = T,
                                marker.selection = c("houseman2012","var"),
                                n.markers = 50,
                                min.n.beads = 5,
                                min.int.quant = 0.05,
                                max.int.quant = 0.95,
                                filter.na = T,
                                filter.snp  = T,
                                filter.context = FALSE,
                                filter.somatic = FALSE,
                                normalization="wm.dasen")
```

```
## [1] "Did not write the variable dump: should only be executed from an environment with all the variables set"
## [2019-09-27 15:02:03, Main:] checking inputs
## [2019-09-27 15:02:03, Main:] preparing data
## [2019-09-27 15:02:03, Main:] preparing jobs
## [2019-09-27 15:02:03, Main:] 176 factorization runs in total
## [2019-09-27 15:02:58, Main:] 100 runs complete
## [2019-09-27 15:03:30, Main:] finished all jobs. Creating the object
```

# Data Import through RnBeads

We recommend to use the [RnBeads](https://rnbeads.org) package to provide methylation data to *DecompPipeline*. *RnBeads* can handle BS and array-based datasets and provides an extensive toolset. BS data can be loaded directly from BED-files generated through methylation data mapping and calling software, such as bismark or bsmap. For array-based datasets, IDAT-files can be directly loaded or GEO accession numbers provided to download data from the repository, among other import options. We refer to the [RnBeads vignette](http://bioconductor.org/packages/release/bioc/vignettes/RnBeads/inst/doc/RnBeads.pdf) for further descriptions on the data import options. Data Import is the only module to be exectued ahead of using *DecompPipeline*.


```r
idat.dir <- "~/idats"
sample.annotation <- "~/sample_annotation.csv"
rnb.set <- rnb.execute.import(data.source = c(idat.dir,sample.annotation),data.type = "infinium.idat.dir")
```

# R session
Here is the output of `sessionInfo()` on the system on which this document was compiled:

```
## R version 3.5.3 (2019-03-11)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 8 (jessie)
## 
## Matrix products: default
## BLAS: /usr/lib/atlas-base/atlas/libblas.so.3.0
## LAPACK: /usr/lib/atlas-base/atlas/liblapack.so.3.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] RnBeads.hg19_1.14.0                    
##  [2] DecompPipeline_0.2.5                   
##  [3] R.utils_2.7.0                          
##  [4] R.oo_1.22.0                            
##  [5] R.methodsS3_1.7.1                      
##  [6] MeDeCom_0.2.2                          
##  [7] RnBeads_2.3.1                          
##  [8] plyr_1.8.4                             
##  [9] methylumi_2.28.0                       
## [10] minfi_1.28.0                           
## [11] bumphunter_1.24.5                      
## [12] locfit_1.5-9.1                         
## [13] iterators_1.0.10                       
## [14] foreach_1.4.4                          
## [15] Biostrings_2.50.1                      
## [16] XVector_0.22.0                         
## [17] SummarizedExperiment_1.12.0            
## [18] DelayedArray_0.8.0                     
## [19] BiocParallel_1.16.2                    
## [20] FDb.InfiniumMethylation.hg19_2.2.0     
## [21] org.Hs.eg.db_3.7.0                     
## [22] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
## [23] GenomicFeatures_1.34.1                 
## [24] AnnotationDbi_1.44.0                   
## [25] reshape2_1.4.3                         
## [26] scales_1.0.0                           
## [27] Biobase_2.42.0                         
## [28] illuminaio_0.24.0                      
## [29] matrixStats_0.54.0                     
## [30] limma_3.38.2                           
## [31] gridExtra_2.3                          
## [32] ggplot2_3.1.0                          
## [33] fields_9.6                             
## [34] maps_3.3.0                             
## [35] spam_2.2-0                             
## [36] dotCall64_1.0-0                        
## [37] ff_2.2-14                              
## [38] bit_1.1-14                             
## [39] cluster_2.0.7-1                        
## [40] MASS_7.3-51.1                          
## [41] GenomicRanges_1.34.0                   
## [42] GenomeInfoDb_1.18.1                    
## [43] IRanges_2.16.0                         
## [44] S4Vectors_0.20.1                       
## [45] BiocGenerics_0.28.0                    
## [46] RUnit_0.4.32                           
## [47] gplots_3.0.1                           
## [48] gtools_3.8.1                           
## [49] pracma_2.2.2                           
## [50] Rcpp_1.0.0                             
## [51] usethis_1.4.0                          
## [52] devtools_2.0.1                         
## [53] knitr_1.21                             
## 
## loaded via a namespace (and not attached):
##  [1] backports_1.1.2          lazyeval_0.2.1          
##  [3] splines_3.5.3            digest_0.6.18           
##  [5] gdata_2.18.0             magrittr_1.5            
##  [7] memoise_1.1.0            remotes_2.0.2           
##  [9] readr_1.2.1              annotate_1.60.0         
## [11] siggenes_1.56.0          prettyunits_1.0.2       
## [13] colorspace_1.3-2         blob_1.1.1              
## [15] xfun_0.4                 dplyr_0.7.8             
## [17] callr_3.0.0              crayon_1.3.4            
## [19] RCurl_1.95-4.11          genefilter_1.64.0       
## [21] bindr_0.1.1              GEOquery_2.50.0         
## [23] survival_2.43-3          glue_1.3.0              
## [25] registry_0.5             gtable_0.2.0            
## [27] zlibbioc_1.28.0          pkgbuild_1.0.2          
## [29] Rhdf5lib_1.4.1           HDF5Array_1.10.0        
## [31] DBI_1.0.0                rngtools_1.3.1          
## [33] bibtex_0.4.2             xtable_1.8-3            
## [35] progress_1.2.0           mclust_5.4.2            
## [37] preprocessCore_1.44.0    httr_1.3.1              
## [39] RColorBrewer_1.1-2       pkgconfig_2.0.2         
## [41] reshape_0.8.8            XML_3.98-1.16           
## [43] tidyselect_0.2.5         rlang_0.3.1             
## [45] munsell_0.5.0            tools_3.5.3             
## [47] cli_1.0.1                RSQLite_2.1.1           
## [49] evaluate_0.12            stringr_1.3.1           
## [51] processx_3.2.0           bit64_0.9-7             
## [53] fs_1.2.6                 beanplot_1.2            
## [55] caTools_1.17.1.1         purrr_0.2.5             
## [57] bindrcpp_0.2.2           nlme_3.1-137            
## [59] doRNG_1.7.1              nor1mix_1.2-3           
## [61] xml2_1.2.0               biomaRt_2.38.0          
## [63] compiler_3.5.3           tibble_1.4.2            
## [65] stringi_1.2.4            ps_1.2.1                
## [67] desc_1.2.0               lattice_0.20-38         
## [69] Matrix_1.2-16            multtest_2.38.0         
## [71] pillar_1.3.0             data.table_1.11.8       
## [73] bitops_1.0-6             rtracklayer_1.42.1      
## [75] R6_2.3.0                 KernSmooth_2.23-15      
## [77] sessioninfo_1.1.1        codetools_0.2-15        
## [79] assertthat_0.2.0         pkgload_1.0.2           
## [81] rhdf5_2.26.0             openssl_1.1             
## [83] pkgmaker_0.27            rprojroot_1.3-2         
## [85] withr_2.1.2              GenomicAlignments_1.18.0
## [87] Rsamtools_1.34.0         GenomeInfoDbData_1.2.0  
## [89] hms_0.4.2                quadprog_1.5-5          
## [91] tidyr_0.8.2              base64_2.0              
## [93] DelayedMatrixStats_1.4.0 base64enc_0.1-3
```

# References
