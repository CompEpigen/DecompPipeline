---
title: "DNA Methylation Deconvolution Protocol"
author: "Michael Scherer, Pavlo Lutsik, Petr Nazarov, Tony Kamoa"
date: "July 31, 2019"
output:
  pdf_document: default
    toc: true
    toc_depth: 3
    df_print: kable
    highlight: tango
    documentclass: scrartcl
    mainfont: ebgaramond
bibliography: bibliography.bib
---



# Introduction

This protocol aims at guiding reasearcher how to employ deconvolution of methylomes obtained from complex tissue. It will start with data retrieval from a public resource, but is equally applicable to in-house generated data. We will furthermore focus on the Illumina BeadChip series as a data source, although the protocol is also compatible with bisulfite sequencing that provides single base pair resolution.
Deconvolution here refers to creating two matrices (proportion matrix A and methylation pattern matrix T) from a single matrix of input DNA methylation data (dimension CpGs x samples). Non-negative matrix factorization can be employed for this task, and we will discuss some of the advantages and caveats of the methods.

# Protocol

## Data Retrival

### Obtaining data from a public resource (duration ~5h)

We focus on DNA methylation data from cancer patients that has been generated in The Cancer Genome Atlas (TCGA) project. Since lung cancer has been shown to be a premier candidate for DNA methylation based deconvolution, we selected the lung adenocarcinoma dataset from the TCGA website (dataset TCGA-LUAD, https://portal.gdc.cancer.gov/legacy-archive/search/f). The dataset was generated using the Illumina Infinum 450k BeadChip and comprises 461 samples. The clinical metadata of the samples is available at https://portal.gdc.cancer.gov/projects/TCGA-LUAD and lists 585 samples. The discrepancy between the number comes from recent progress within TCGA. We used the Genomic Data Commons (GDC) data download tool (https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) to download the intensity data (IDAT) files listed in the manifest file and its associated metadata. This metadata also includes the mapping between each of the samples and the IDAT files. To create a final mapping and to prepare the files for downstream analysis, the following code was employed.


```r
clinical.data <- read.table("annotation/clinical.tsv",sep="\t",header=T)
idat.files <- list.files("idat",full.names = T)
meta.files <- list.files(idat.files[1],full.names = T)
untar(meta.files[3],exdir = idat.files[1])
meta.files <- untar(meta.files[3],list=T)
meta.info <- read.table(file.path(idat.files[1],meta.files[5]),sep="\t",header=T)
meta.info <- meta.info[match(unique(meta.info$Comment..TCGA.Barcode.),meta.info$Comment..TCGA.Barcode.),]
match.meta.clin <- match(clinical.data$submitter_id,substr(meta.info$Comment..TCGA.Barcode.,1,12))
anno.frame <- na.omit(data.frame(clinical.data,meta.info[match.meta.clin,]))
anno.frame$barcode <- unlist(lapply(lapply(as.character(anno.frame$Array.Data.File),function(x)strsplit(x,"_")),function(x)paste(x[[1]][1],x[[1]][2],sep="_")))
anno.frame$Sentrix_ID <- unlist(lapply(lapply(as.character(anno.frame$Array.Data.File),function(x)strsplit(x,"_")),function(x)paste(x[[1]][1])))
anno.frame$Sentrix_Position <- unlist(lapply(lapply(as.character(anno.frame$Array.Data.File),function(x)strsplit(x,"_")),function(x)paste(x[[1]][2])))
anno.frame$healthy_cancer <- ifelse(grepl("11A",anno.frame$Comment..TCGA.Barcode.),"healthy","cancer")
write.table(anno.frame,"annotation/sample_annotation.tsv",quote=F,row.names = F,sep="\t")
anno.frame <- read.table("annotation/sample_annotation.tsv",quote=F,row.names = F,sep="\t")

#' write idat files to parent directory
lapply(idat.files,function(x){
  is.idat <- list.files(x,pattern = ".idat",full.names = T)
  file.copy(is.idat,"idat/")
  unlink(x,recursive = T)
})
```

## Data Processing

### Data Import and Quality Control in RnBeads (~3h)

After downloading the data, it has to be processed into a format that can be used by downstream software. We used RnBeads to convert the files into a data object and performed basic quality control steps on the dataset. Most notably, analysis options need to be specified for RnBeads, either through an XML file or in the command line. We will follow the latter strategy here, and deactivate the preprocessing, exploratory, covariate inference and differential methylation modules. In the next step, we specify the input to RnBeads: the created sample annotation sheet, the folder in which the IDAT files are stored and a folder to which the HTML report is to be saved. We additionally recommend to specify a temporary directory for the analysis. Then we start the RnBeads analysis.


```r
suppressPackageStartupMessages(library(RnBeads))
rnb.options(
  assembly="hg19",
  identifiers.column="submitter_id",
  import=T,
  import.default.data.type="idat.dir",
  import.table.separator="\t",
  import.sex.prediction=T,
  qc=T,
  preprocessing=F,
  exploratory=F,
  inference=F,
  differential=F,
  export.to.bed=F,
  export.to.trackhub=NULL,
  export.to.csv=F
)
sample.anno <- "annotation/sample_annotation.tsv"
idat.folder <- "idat/"
dir.report <- paste0("report",Sys.Date(),"/")
temp.dir <- "/DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp"
options(fftempdir=temp.dir)
rnb.set <- rnb.run.analysis(dir.reports = dir.report, sample.sheet = sample.anno, data.dir = idat.folder)
```

```
## 2019-08-31 09:30:43     1.1  STATUS STARTED RnBeads Pipeline
## 2019-08-31 09:30:43     1.1    INFO     Initialized report index and saved to index.html
## 2019-08-31 09:30:44     1.1  STATUS     STARTED Loading Data
## 2019-08-31 09:30:44     1.1    INFO         Number of cores: 1
## 2019-08-31 09:30:44     1.1    INFO         Loading data of type "idat.dir"
## 2019-08-31 09:30:44     1.1  STATUS         STARTED Loading Data from IDAT Files
## 2019-08-31 09:30:46     1.1    INFO             Detected platform: HumanMethylation450
## 2019-08-31 09:57:14     1.5  STATUS         COMPLETED Loading Data from IDAT Files
## 2019-08-31 11:00:07     2.0  STATUS         Loaded data from idat/
## 2019-08-31 11:01:10     7.6  STATUS         Predicted sex for the loaded samples
## 2019-08-31 11:01:41     7.1  STATUS         Added data loading section to the report
## 2019-08-31 11:01:41     7.1  STATUS         Loaded 461 samples and 485577 sites
## 2019-08-31 11:01:41     7.1    INFO         Output object is of type RnBeadRawSet
## 2019-08-31 11:01:41     7.1  STATUS     COMPLETED Loading Data
## 2019-08-31 11:15:50     7.1    INFO     Initialized report index and saved to index.html
## 2019-08-31 11:15:51     7.1  STATUS     STARTED Quality Control
## 2019-08-31 11:15:51     7.1    INFO         Number of cores: 1
## 2019-08-31 11:15:51     7.1  STATUS         STARTED Quality Control Section
```

```
## 2019-08-31 11:16:21     2.0  STATUS             Added quality control box plots
```

```
## 2019-08-31 11:21:14     2.0  STATUS             Added quality control bar plots
```

```
## 2019-08-31 11:21:44     2.0  STATUS             Added negative control boxplots
## 2019-08-31 11:21:44     2.0  STATUS         COMPLETED Quality Control Section
## 2019-08-31 11:21:44     2.0  STATUS         STARTED Visualizing SNP Probe Data
## 2019-08-31 11:21:44     2.0  STATUS             STARTED Mixups Visualization Section
```

```
## 2019-08-31 11:22:20     5.4  STATUS                 Added SNP Heatmap
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
```

```
## Also defined by 'spam'
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
```

```
## Also defined by 'spam'
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
```

```
## Also defined by 'spam'
```

```
## 2019-08-31 11:22:20     5.4  STATUS                 Calculated Manhattan distances between samples based on SNP probes
```

```
## 2019-08-31 11:22:24     5.4  STATUS                 Added SNP-based Distances
## 2019-08-31 11:22:24     5.4  STATUS             COMPLETED Mixups Visualization Section
## 2019-08-31 11:22:24     5.4  STATUS         COMPLETED Visualizing SNP Probe Data
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ffaf03c64acb5.ff
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ffaf033c86bb6.ff
```

```
## 2019-08-31 11:23:21     7.4  STATUS     COMPLETED Quality Control
## 2019-08-31 11:23:21     7.4    INFO     Initialized report index and saved to index.html
## 2019-08-31 11:23:21     7.4  STATUS     STARTED Saving RData
## 2019-08-31 11:23:21     7.4  STATUS     COMPLETED Saving RData
## 2019-08-31 11:23:21     7.4  STATUS COMPLETED RnBeads Pipeline
```

RnBeads creates an interactive HTML report, specifying the steps performed and the associated results. Data was of good quality such that is can be used for further analysis. (Include two screenshots from the RnBeads report)

### Preprocessing and Filtering

For further analysis, we use the DecompPipeline package (https://github.com/lutsik/DecompPipeline), which provides a comprehensive workflow including crucial data preparation steps for methylome deconvolution experiments. The options are provided through the individual function parameters. We follow a stringent filtering strategy. First, all samples having fewer than 3 beads covered are filtered, as well as those probes that are in the 0.05 and 0.95 overall intensity quantiles, respectively. We then remove all probes containing missing values, outside of CpG context, that overlap with annotated SNPs, on the sex chromosomes and probes that have been shown to be cross-reactive on the chip. Then, BMIQ normalization [@bmiq] is employed to account for the chip's design bias.
Accounting for potential confounding factor is crucial in epigenomic studies. Especially, the influence of donor sex and age on the DNA methylation pattern is well-studied and strong. Furthermore, genetic differences between groups of individuals due to different origins may influence the DNA methylation paterrn. We used Independent Component Analysis (ICA) to account for DNA methylation differences that are due to these confounding factors. ICA detects components in the data accounting for most of the variance similar to PCA, but does not require orthogonality of the components but statistical independence. We used an external library (http://sablab.net/scripts/LibICA.r) for performing ICA to adjust for sex, age, race and ethnicity.

```r
suppressPackageStartupMessages(library(DecompPipeline))
data.prep <- prepare_data(RNB_SET = rnb.set,
                          analysis.name = "TCGA_LUAD2",
                          NORMALIZATION = "bmiq",
                          FILTER_BEADS = T,
                          MIN_N_BEADS = 3,
                          FILTER_INTENSITY = T,
                          MIN_INT_QUANT = 0.001,
                          MAX_INT_QUANT = 0.999,
                          FILTER_NA = T,
                          FILTER_CONTEXT = T,
                          FILTER_SNP = T,
                          FILTER_SOMATIC = T,
                          FILTER_CROSS_REACTIVE = T,
                          execute.lump=T,
                          remove.ICA=T,
                          conf.fact.ICA=c("age_at_diagnosis","race","gender","ethnicity"),
                          ica.setting=c("alpha.fact"=1e-5,"save.report"=T,"ntry"=10,"nmax"=50,"ncores"=10))
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ffaf0605f749e.ff
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ffaf067eb4ae1.ff
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ffaf03e38f42.ff
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ffaf03b2f13f9.ff
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ffaf0372921bb.ff
```

```
## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion

## Warning in getComponentNumber(rnb.set, conf.factor, nmin = nmin, nmax =
## nmax, : NAs introduced by coercion
```

```
## Loading required package: fastICA
```

```
## Loading required package: doMC
```

```
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 10 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 13:59:37 
## *** Done! 
## *** End time: 2019-08-31 17:02:18 
## Time difference of 3.044657 hours
## List of 10
##  $ : num [1:230223, 1:10] -0.4004 -0.415 -0.1183 -0.0512 -0.1045 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:10] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:10] 0.5778 0.238 0.0955 -0.0368 0.0677 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:10] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:10] -0.4041 -0.4215 -0.1232 -0.0562 -0.1238 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:10] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:10] -0.966 -0.678 -1.085 1.314 0.958 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:10] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:10] -0.943 -0.646 -1.094 1.305 0.942 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:10] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:10] 0.35013 0.34455 0.20516 -0.1234 0.00404 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:10] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:10] -0.4142 -0.4391 -0.1111 -0.0435 -0.0971 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:10] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:10] 0.5682 0.2421 0.0971 -0.0371 0.0741 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:10] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:10] -0.3985 -0.4167 -0.1202 -0.0529 -0.1158 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:10] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:10] 0.408 0.4243 0.0969 0.0423 0.0851 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:10] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 11 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 17:03:29 
## *** Done! 
## *** End time: 2019-08-31 17:38:19 
## Time difference of 34.83294 mins
## List of 10
##  $ : num [1:230223, 1:11] 0.4069 0.2793 0.2051 -0.1103 0.0771 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:11] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:11] -0.959 -0.696 -1.087 1.313 0.971 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:11] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:11] -0.958 -0.676 -1.088 1.308 0.975 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:11] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:11] -0.3926 -0.4252 -0.1247 -0.0569 -0.1091 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:11] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:11] -0.939 -0.629 -1.073 1.314 0.946 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:11] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:11] -0.6852 -0.0883 -0.1206 0.0163 -0.1997 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:11] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:11] -0.3116 0.0709 0.1591 0.0676 -0.1193 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:11] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:11] -0.948 -0.7 -1.087 1.315 0.969 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:11] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:11] -0.942 -0.692 -1.088 1.315 0.961 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:11] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:11] -0.3332 -0.0605 0.1395 0.1412 -0.1176 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:11] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 12 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 17:40:04 
## *** Done! 
## *** End time: 2019-08-31 17:49:45 
## Time difference of 9.679189 mins
## List of 10
##  $ : num [1:230223, 1:12] 0.4327 0.0409 0.1508 0.038 0.4025 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:12] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:12] 0.41958 0.43849 0.12737 0.02707 -0.00168 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:12] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:12] 0.826 0.499 1.051 -1.35 -0.881 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:12] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:12] -0.411194 -0.424342 -0.126414 -0.021715 0.000994 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:12] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:12] -0.4385 -0.0335 -0.1513 -0.0409 -0.4256 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:12] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:12] 0.5036 0.0265 0.1055 0.0558 0.4052 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:12] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:12] 0.97 0.69 1.086 -1.313 -0.982 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:12] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:12] -0.3677 0.0476 -0.0901 -0.1104 -0.4792 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:12] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:12] 0.4938 -0.0882 0.0744 0.1234 0.4833 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:12] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:12] -0.45006 -0.00549 -0.10429 -0.08245 -0.45764 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:12] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 13 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 17:51:51 
## *** Done! 
## *** End time: 2019-08-31 18:01:23 
## Time difference of 9.530311 mins
## List of 10
##  $ : num [1:230223, 1:13] -0.0179 0.2083 -0.0241 -0.1523 -0.3926 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:13] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:13] 0.777 0.414 1.12 -1.385 -0.795 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:13] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:13] 0.986 0.673 1.09 -1.313 -0.971 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:13] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:13] 0.3725 0.4743 0.0975 0.0282 -0.0224 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:13] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:13] 0.3957 0.4808 0.1207 0.0287 -0.0301 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:13] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:13] -0.983 -0.689 -1.108 1.31 0.948 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:13] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:13] 0.989 0.681 1.097 -1.311 -0.96 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:13] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:13] 0.37785 0.45148 0.1196 0.03112 -0.00767 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:13] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:13] 0.3146 -0.0343 0.083 0.1394 0.4459 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:13] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:13] 0.838 0.516 1.048 -1.346 -0.898 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:13] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 14 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 18:03:02 
## *** Done! 
## *** End time: 2019-08-31 18:11:05 
## Time difference of 8.062823 mins
## List of 10
##  $ : num [1:230223, 1:14] 0.4143 0.4865 0.0954 0.0331 -0.0319 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:14] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:14] -0.4156 -0.4811 -0.1053 -0.0404 0.0288 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:14] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:14] -0.4503 -0.0892 -0.1243 -0.1026 -0.348 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:14] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:14] 0.446 0.077 0.119 0.112 0.35 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:14] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:14] 0.995 0.68 1.089 -1.311 -0.978 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:14] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:14] 1.013 0.714 1.115 -1.291 -0.961 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:14] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:14] -0.993 -0.689 -1.09 1.312 0.963 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:14] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:14] 0.227 -0.116 0.042 0.195 0.425 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:14] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:14] -0.813 -0.506 -1.05 1.352 0.896 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:14] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:14] -0.435 -0.074 -0.114 -0.123 -0.366 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:14] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 15 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 18:12:33 
## *** Done! 
## *** End time: 2019-08-31 18:23:36 
## Time difference of 11.0486 mins
## List of 10
##  $ : num [1:230223, 1:15] -0.894 -0.569 -1.056 1.347 0.945 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:15] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:15] 0.34853 0.24177 0.08601 -0.00763 0.09155 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:15] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:15] 0.35532 0.23312 0.09579 -0.00543 0.09606 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:15] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:15] 0.35639 0.23036 0.08624 -0.00735 0.09522 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:15] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:15] 0.79 0.488 1.057 -1.357 -0.872 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:15] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:15] -0.948 -0.626 -1.058 1.343 0.955 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:15] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:15] -0.38416 -0.34689 -0.08244 -0.00451 -0.03354 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:15] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:15] 0.35744 0.22195 0.08831 -0.00572 0.10849 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:15] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:15] -0.961 -0.623 -1.064 1.333 0.972 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:15] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:15] 0.883 0.558 1.057 -1.348 -0.927 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:15] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 16 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 18:25:30 
## *** Done! 
## *** End time: 2019-08-31 18:32:39 
## Time difference of 7.161136 mins
## List of 10
##  $ : num [1:230223, 1:16] 1.007 0.677 1.096 -1.31 -0.983 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:16] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:16] -0.4077 -0.2838 -0.1121 0.0158 -0.0347 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:16] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:16] 0.4027 0.2798 0.1126 -0.0162 0.0264 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:16] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:16] -0.409 -0.276 -0.1103 0.0165 -0.0387 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:16] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:16] 0.927 0.592 1.052 -1.355 -0.972 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:16] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:16] -0.2192 0.5817 0.1795 -0.0841 -0.0326 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:16] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:16] -1 -0.666 -1.089 1.316 0.967 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:16] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:16] -0.421 -0.374 -0.1122 0.0127 0.0159 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:16] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:16] 0.4191 -0.0664 0.0832 0.0445 0.2762 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:16] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:16] 0.387 -0.1806 0.0379 0.0793 0.3577 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:16] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 17 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 18:34:29 
## *** Done! 
## *** End time: 2019-08-31 18:38:51 
## Time difference of 4.359877 mins
## List of 10
##  $ : num [1:230223, 1:17] -0.923 -0.596 -1.07 1.346 0.971 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:17] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:17] 0.3567 -0.0339 0.0856 0.1185 0.2358 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:17] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:17] -1.021 -0.645 -1.098 1.301 0.94 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:17] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:17] 0.3693 -0.0575 0.0993 0.2203 0.2077 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:17] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:17] 0.3605 -0.0319 0.0891 0.1163 0.2336 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:17] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:17] -0.4273 0.0292 -0.0931 -0.1984 -0.1858 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:17] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:17] 0.2276 0.4768 0.2003 0.0603 -0.147 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:17] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:17] 0.99 0.666 1.091 -1.322 -0.962 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:17] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:17] -0.38494 -0.00653 -0.10689 -0.09005 -0.2231 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:17] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:17] -0.2565 -0.4989 -0.1944 -0.0594 0.1575 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:17] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 18 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 18:40:26 
## *** Done! 
## *** End time: 2019-08-31 18:45:14 
## Time difference of 4.795503 mins
## List of 10
##  $ : num [1:230223, 1:18] -0.363 0.132 -0.089 -0.193 -0.253 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:18] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:18] 0.3466 -0.0269 0.0924 0.1021 0.2088 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:18] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:18] 0.989 0.683 1.098 -1.284 -0.976 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:18] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:18] 0.974 0.685 1.084 -1.292 -1.005 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:18] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:18] 0.3834 0.0122 0.1083 0.0901 0.2335 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:18] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:18] 0.987 0.688 1.089 -1.293 -0.981 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:18] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:18] -0.962 -0.671 -1.084 1.3 0.984 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:18] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:18] 0.241 0.1053 0.07 0.1949 -0.0965 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:18] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:18] -0.3853 -0.0124 -0.1072 -0.0902 -0.2247 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:18] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:18] -0.986 -0.691 -1.094 1.288 0.98 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:18] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 19 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 18:46:41 
## *** Done! 
## *** End time: 2019-08-31 18:51:34 
## Time difference of 4.879213 mins
## List of 10
##  $ : num [1:230223, 1:19] 0.963 0.673 1.092 -1.3 -0.974 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:19] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:19] 0.1958 0.5147 0.1973 0.0941 -0.1098 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:19] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:19] -0.3405 0.0351 -0.0826 -0.1081 -0.2107 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:19] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:19] 0.3328 -0.0646 0.0756 0.1147 0.2125 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:19] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:19] 0.3519 -0.0504 0.0704 0.1139 0.2348 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:19] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:19] 0.849 0.631 1.061 -1.324 -1.009 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:19] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:19] -0.955 -0.66 -1.085 1.305 0.983 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:19] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:19] 0.24552 0.00104 0.05827 0.17526 0.03279 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:19] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:19] 0.21527 -0.00509 0.05321 0.15985 0.02273 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:19] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:19] 0.3719 -0.0275 0.0787 0.1059 0.2036 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:19] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 20 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 18:52:17 
## *** Done! 
## *** End time: 2019-08-31 18:58:00 
## Time difference of 5.716855 mins
## List of 10
##  $ : num [1:230223, 1:20] -0.3778 0.0491 -0.05 -0.1027 -0.27 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:20] -0.2275 -0.5376 -0.1722 -0.0829 0.0584 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:20] -0.1727 -0.0346 -0.1028 -0.1816 0.0974 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:20] 0.2057 0.5021 0.176 0.0873 -0.0567 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:20] -0.2136 -0.5353 -0.1773 -0.0897 0.0726 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:20] -0.2155 -0.5349 -0.1786 -0.0859 0.0685 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:20] 0.2302 0.0291 0.1019 0.1819 -0.0445 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:20] -0.5 0.333 0.165 -0.194 0.239 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:20] 0.2319 0.5509 0.1833 0.0733 -0.0773 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:20] 0.985 0.67 1.092 -1.297 -0.955 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 21 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 18:58:47 
## *** Done! 
## *** End time: 2019-08-31 19:05:07 
## Time difference of 6.341652 mins
## List of 10
##  $ : num [1:230223, 1:21] 0.2438 -0.0525 -0.0649 -0.2111 0.3099 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:21] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:21] 0.991 0.68 1.094 -1.293 -0.951 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:21] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:21] 0.2111 0.0217 0.1048 0.1775 -0.06 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:21] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:21] -0.2116 -0.5594 -0.1816 -0.0713 0.072 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:21] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:21] 0.2633 -0.0265 0.0909 0.0556 0.2088 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:21] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:21] -0.1729 -0.5538 -0.1965 -0.079 0.0803 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:21] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:21] -0.2659 0.0322 -0.0837 -0.0514 -0.1914 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:21] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:21] -0.4037 -0.2895 -0.0776 0.208 -0.152 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:21] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:21] -0.2495 0.0579 -0.0728 -0.0658 -0.1977 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:21] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:21] -0.224 -0.6006 -0.1691 -0.0762 0.086 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:21] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 22 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 19:05:54 
## *** Done! 
## *** End time: 2019-08-31 19:12:02 
## Time difference of 6.143035 mins
## List of 10
##  $ : num [1:230223, 1:22] -0.35861 -0.00655 -0.11405 -0.07301 -0.1925 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] -0.35971 -0.00611 -0.06698 -0.06948 -0.2027 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] 0.34333 0.00921 0.07797 0.05995 0.22354 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] 0.34723 -0.00277 0.0705 0.06798 0.20097 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] 1.088 0.735 1.082 -1.278 -0.947 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] 0.821 0.641 1.051 -1.328 -0.976 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] 1.059 0.726 1.071 -1.293 -0.96 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] -0.337 -0.0276 -0.0813 -0.065 -0.2136 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] -0.1678 -0.5266 -0.1893 -0.0796 0.076 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] -0.0973 -0.4795 -0.1865 -0.0548 0.081 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 23 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 19:12:50 
## *** Done! 
## *** End time: 2019-08-31 19:19:50 
## Time difference of 6.996149 mins
## List of 10
##  $ : num [1:230223, 1:23] 0.3438 0.0232 0.1178 0.0668 0.2121 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:23] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:23] -0.3337 -0.0209 -0.0854 -0.0573 -0.2178 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:23] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:23] 0.34183 0.00713 0.07137 0.07441 0.21593 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:23] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:23] -0.0599 -0.5554 -0.1993 -0.0788 0.0811 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:23] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:23] -0.3805 -0.0726 -0.075 -0.1655 0.0376 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:23] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:23] -1.011 -0.71 -1.058 1.308 0.975 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:23] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:23] -0.0877 -0.5882 -0.198 -0.0835 0.0803 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:23] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:23] 0.0775 0.5619 0.1994 0.0728 -0.0942 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:23] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:23] -0.3887 -0.0862 -0.0805 -0.176 0.0407 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:23] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:23] -0.596 -0.115 -0.172 -0.344 -0.187 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:23] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 24 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 19:20:39 
## *** Done! 
## *** End time: 2019-08-31 19:28:32 
## Time difference of 7.876361 mins
## List of 10
##  $ : num [1:230223, 1:24] 0.988 0.698 1.06 -1.31 -0.988 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:24] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:24] 1.035 0.712 1.082 -1.292 -0.973 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:24] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:24] 0.4026 0.0303 0.0811 0.0734 0.2344 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:24] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:24] 0.0698 0.5452 0.2028 0.0852 -0.0736 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:24] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:24] -0.3857 -0.0324 -0.0942 -0.0812 -0.1938 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:24] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:24] 0.0442 0.5836 0.192 0.0809 -0.1082 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:24] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:24] 0.3801 0.0243 0.0903 0.072 0.2364 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:24] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:24] 0.3713 0.0162 0.0834 0.0739 0.2361 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:24] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:24] 0.819 0.643 1.053 -1.327 -0.957 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:24] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:24] -0.0641 -0.566 -0.1988 -0.0787 0.1012 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:24] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 25 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 19:29:41 
## *** Done! 
## *** End time: 2019-08-31 19:41:51 
## Time difference of 12.16293 mins
## List of 10
##  $ : num [1:230223, 1:25] -0.0468 -0.5968 -0.1865 -0.0712 0.1167 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:25] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:25] -0.213 -0.243 -0.205 -0.239 0.228 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:25] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:25] -0.0711 -0.5674 -0.1964 -0.0896 0.0789 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:25] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:25] -0.0374 -0.5982 -0.1923 -0.0765 0.1167 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:25] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:25] -0.148 -0.236 -0.155 -0.246 0.244 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:25] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:25] -0.3382 -0.0942 -0.0748 -0.1226 -0.2434 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:25] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:25] -0.95 -0.72 -1.045 1.31 0.993 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:25] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:25] 0.3818 0.1047 0.0705 0.1072 0.2223 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:25] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:25] -0.3634 -0.0961 -0.0718 -0.1084 -0.2093 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:25] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:25] 0.139045 0.460301 0.202012 0.073194 0.000956 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:25] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 26 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 19:44:05 
## *** Done! 
## *** End time: 2019-08-31 19:53:51 
## Time difference of 9.768758 mins
## List of 10
##  $ : num [1:230223, 1:26] 0.454 0.1589 0.0511 0.1052 0.3504 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:26] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:26] -0.1484 -0.4498 -0.1987 -0.0707 -0.0087 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:26] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:26] 0.0472 0.5818 0.1907 0.0711 -0.1326 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:26] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:26] 0.0131 0.5966 0.1906 0.063 -0.1636 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:26] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:26] -0.967 -0.731 -1.046 1.304 0.969 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:26] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:26] 0.457 0.164 0.041 0.109 0.346 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:26] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:26] 0.0149 0.6054 0.1947 0.0642 -0.1696 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:26] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:26] 0.4192 0.1543 0.0516 0.1209 0.3348 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:26] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:26] 0.956 0.724 1.044 -1.307 -0.975 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:26] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:26] 0.0406 0.5795 0.1957 0.0792 -0.1285 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:26] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 27 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 19:54:52 
## *** Done! 
## *** End time: 2019-08-31 20:05:50 
## Time difference of 10.96851 mins
## List of 10
##  $ : num [1:230223, 1:27] -0.4013 -0.1687 -0.0319 -0.1199 -0.3208 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:27] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:27] -0.403 -0.182 -0.031 -0.121 -0.346 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:27] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:27] 0.0484 0.5779 0.2046 0.0763 -0.1732 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:27] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:27] -0.0448 -0.5769 -0.2103 -0.066 0.168 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:27] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:27] 0.0482 0.5545 0.2166 0.0746 -0.1587 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:27] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:27] -0.0312 -0.5872 -0.2071 -0.0643 0.1808 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:27] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:27] 0.461 0.138 -0.015 0.188 0.166 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:27] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:27] -0.0517 -0.5814 -0.2093 -0.0656 0.1699 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:27] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:27] -0.21 -0.269 -0.29 -0.217 0.426 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:27] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:27] -0.0256 -0.5901 -0.2055 -0.0627 0.1825 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:27] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 28 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 20:07:13 
## *** Done! 
## *** End time: 2019-08-31 20:18:22 
## Time difference of 11.15398 mins
## List of 10
##  $ : num [1:230223, 1:28] 0.2983 0.0187 0.3688 0.3626 -0.1618 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:28] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:28] 1.111 0.743 1.09 -1.277 -0.942 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:28] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:28] 0.218 0.217 0.25 0.207 -0.359 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:28] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:28] 0.985 0.724 1.046 -1.307 -0.961 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:28] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:28] 0.0413 0.59 0.2039 0.069 -0.1711 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:28] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:28] 0.247 0.259 0.269 0.216 -0.426 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:28] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:28] 0.258 0.279 0.279 0.207 -0.385 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:28] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:28] -0.1147 -0.5193 -0.2127 -0.0797 0.0968 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:28] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:28] -0.252 -0.288 -0.29 -0.204 0.385 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:28] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:28] 1.137 0.753 1.107 -1.274 -0.912 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:28] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 29 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 20:19:19 
## *** Done! 
## *** End time: 2019-08-31 20:31:51 
## Time difference of 12.53569 mins
## List of 10
##  $ : num [1:230223, 1:29] 0.4148 0.1703 0.0446 0.1253 0.316 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:29] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:29] -0.0855 -0.5565 -0.2591 -0.079 0.0409 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:29] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:29] 0.257 0.274 0.286 0.202 -0.359 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:29] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:29] 0.0697 0.5767 0.245 0.0697 -0.0741 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:29] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:29] -0.373 -0.1753 -0.0517 -0.1272 -0.3151 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:29] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:29] -0.395 -0.1712 -0.0432 -0.1305 -0.3161 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:29] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:29] -0.0872 -0.5554 -0.2507 -0.0804 0.0503 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:29] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:29] -0.1362 -0.5553 -0.24 -0.0765 0.0578 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:29] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:29] -0.3887 -0.1846 -0.0487 -0.1242 -0.3168 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:29] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:29] -0.3754 -0.1793 -0.0524 -0.1204 -0.3127 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:29] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 30 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 20:32:57 
## *** Done! 
## *** End time: 2019-08-31 20:47:30 
## Time difference of 14.55301 mins
## List of 10
##  $ : num [1:230223, 1:30] -0.287 -0.26 -0.284 -0.193 0.336 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:30] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:30] 0.915 0.699 1.038 -1.309 -0.961 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:30] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:30] 0.44036 0.11754 -0.00502 0.17102 0.16064 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:30] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:30] 0.28 0.268 0.29 0.191 -0.347 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:30] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:30] 0.393 0.185 0.062 0.113 0.323 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:30] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:30] -1.131 -0.747 -1.111 1.279 0.905 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:30] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:30] -0.0961 -0.5166 -0.2459 -0.0708 0.0109 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:30] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:30] -0.0886 -0.517 -0.2439 -0.0681 -0.0124 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:30] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:30] -0.144 -0.5077 -0.255 -0.0858 0.0098 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:30] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:30] -0.44369 -0.12292 -0.00352 -0.17692 -0.13843 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:30] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 31 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 20:48:33 
## *** Done! 
## *** End time: 2019-08-31 21:03:12 
## Time difference of 14.65043 mins
## List of 10
##  $ : num [1:230223, 1:31] 0.97 0.7 1.03 -1.326 -0.944 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:31] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:31] 0.1478 0.5065 0.2267 0.0361 0.0398 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:31] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:31] 0.953 0.692 1.026 -1.328 -0.953 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:31] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:31] -1.115 -0.739 -1.083 1.303 0.896 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:31] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:31] -0.977 -0.714 -1.034 1.319 0.939 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:31] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:31] -0.991 -0.714 -1.036 1.323 0.931 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:31] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:31] 0.1197 0.4848 0.2345 0.0543 0.031 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:31] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:31] -0.3892 -0.1642 -0.0427 -0.1227 -0.3118 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:31] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:31] -0.976 -0.717 -1.04 1.323 0.951 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:31] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:31] 0.4468 0.1352 0.0162 0.1973 0.1487 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:31] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 32 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 21:04:15 
## *** Done! 
## *** End time: 2019-08-31 21:20:17 
## Time difference of 16.04605 mins
## List of 10
##  $ : num [1:230223, 1:32] 0.4037 0.1665 0.0427 0.1151 0.3301 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:32] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:32] 1.003 0.713 1.035 -1.323 -0.923 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:32] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:32] -0.225 -0.25 -0.289 -0.209 0.372 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:32] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:32] -0.08226 -0.49116 -0.22365 -0.05028 0.00763 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:32] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:32] -0.3893 -0.1825 -0.0512 -0.1215 -0.3263 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:32] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:32] -0.26 -0.257 -0.292 -0.193 0.359 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:32] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:32] 0.0223 0.5148 0.2212 0.054 -0.026 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:32] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:32] 0.3919 0.1735 0.0533 0.1056 0.3439 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:32] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:32] -0.3625 -0.1768 -0.0443 -0.1267 -0.2721 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:32] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:32] -0.3935 -0.1694 -0.0408 -0.1179 -0.3046 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:32] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 33 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 21:21:20 
## *** Done! 
## *** End time: 2019-08-31 21:39:07 
## Time difference of 17.78439 mins
## List of 10
##  $ : num [1:230223, 1:33] -0.3274 -0.1775 -0.0358 -0.1172 -0.2939 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:33] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:33] 0.4278 0.1699 -0.0236 0.1892 0.1428 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:33] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:33] 0.4134 0.1643 -0.0104 0.186 0.1344 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:33] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:33] -0.1125 -0.0852 -0.39 -0.3919 0.2194 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:33] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:33] 0.4455 0.1762 -0.0107 0.1908 0.1667 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:33] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:33] -0.3344 -0.214 -0.0464 -0.0973 -0.2922 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:33] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:33] -0.2817 -0.2011 -0.0464 -0.1253 -0.293 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:33] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:33] -0.3557 -0.1549 -0.0389 -0.1129 -0.3219 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:33] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:33] -0.0487 -0.504 -0.2261 -0.0496 0.014 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:33] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:33] -0.3401 -0.1877 -0.0377 -0.115 -0.2813 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:33] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 34 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 21:40:11 
## *** Done! 
## *** End time: 2019-08-31 21:59:08 
## Time difference of 18.94929 mins
## List of 10
##  $ : num [1:230223, 1:34] 0.372 0.274 0.277 0.144 -0.26 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:34] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:34] -0.4138 -0.217 -0.0431 -0.0608 -0.3475 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:34] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:34] 0.2173 0.5313 0.2232 -0.0059 0.1017 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:34] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:34] -0.20388 -0.53008 -0.22349 0.00418 -0.08948 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:34] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:34] -0.24165 -0.51264 -0.23356 0.00161 -0.13002 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:34] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:34] 1.006 0.707 1.032 -1.327 -0.921 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:34] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:34] 0.352 0.218 0.262 0.14 -0.221 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:34] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:34] 0.921 0.685 1.023 -1.327 -0.939 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:34] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:34] 0.25555 0.51141 0.22872 -0.00419 0.13314 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:34] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:34] -0.18629 -0.52934 -0.22449 0.00107 -0.08359 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:34] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 35 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 22:00:14 
## *** Done! 
## *** End time: 2019-08-31 22:20:46 
## Time difference of 20.53046 mins
## List of 10
##  $ : num [1:230223, 1:35] -0.4262 -0.1986 -0.0432 -0.0702 -0.3681 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:35] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:35] -0.206441 -0.514672 -0.22838 -0.000518 -0.097577 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:35] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:35] 0.4206 0.2023 0.0474 0.0618 0.3562 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:35] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:35] 0.4284 0.2035 0.0444 0.0516 0.3647 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:35] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:35] 0.978 0.718 1.024 -1.327 -0.917 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:35] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:35] 0.2175 0.5184 0.22 -0.0156 0.0403 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:35] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:35] 0.339 0.256 0.277 0.155 -0.26 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:35] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:35] 0.4276 0.1729 0.0383 0.0656 0.3727 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:35] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:35] -0.20139 -0.51523 -0.2226 0.00476 -0.09501 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:35] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:35] -0.254168 -0.50173 -0.236758 0.000353 -0.125986 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:35] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 36 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 22:21:58 
## *** Done! 
## *** End time: 2019-08-31 22:43:23 
## Time difference of 21.41267 mins
## List of 10
##  $ : num [1:230223, 1:36] -0.29007 -0.52097 -0.26514 -0.00908 -0.15016 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:36] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:36] 0.4573 0.1814 0.0671 0.0771 0.3904 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:36] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:36] 0.27466 0.51895 0.26113 0.00937 0.14509 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:36] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:36] -0.2825 -0.2382 0.024 -0.2232 -0.0377 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:36] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:36] 0.24536 0.541262 0.252433 -0.000537 0.120486 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:36] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:36] 0.456 0.1865 0.0677 0.0819 0.389 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:36] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:36] -0.294 -0.268 -0.225 -0.12 0.284 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:36] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:36] -0.4433 -0.1993 -0.0713 -0.0797 -0.3797 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:36] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:36] -0.1063 -0.3681 -0.2707 -0.0339 0.2454 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:36] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:36] 0.4636 0.2276 0.0762 0.0729 0.3731 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:36] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 37 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 22:44:35 
## *** Done! 
## *** End time: 2019-08-31 23:07:05 
## Time difference of 22.49572 mins
## List of 10
##  $ : num [1:230223, 1:37] -0.2697 -0.4812 -0.2586 -0.0215 -0.1361 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:37] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:37] -0.2956 -0.1397 0.0531 -0.2048 -0.0689 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:37] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:37] -0.263 -0.508 -0.249 -0.012 -0.109 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:37] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:37] 0.1143 -0.1209 -0.0572 -0.2626 0.3663 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:37] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:37] 0.265 0.3547 0.2298 0.0941 -0.2595 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:37] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:37] -0.992 -0.725 -1.033 1.324 0.909 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:37] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:37] -0.4699 -0.2007 -0.0698 -0.0699 -0.3648 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:37] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:37] -0.953 -0.705 -1.025 1.325 0.939 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:37] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:37] -0.3023 -0.4831 -0.2559 -0.0139 -0.1502 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:37] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:37] 0.25616 0.50937 0.24623 0.00667 0.11533 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:37] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 38 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 23:08:17 
## *** Done! 
## *** End time: 2019-08-31 23:32:29 
## Time difference of 24.20173 mins
## List of 10
##  $ : num [1:230223, 1:38] -0.4542 -0.197 -0.0714 -0.0898 -0.3611 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:38] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:38] 0.2902 0.1465 -0.0349 0.2086 0.0392 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:38] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:38] -0.298 -0.2759 -0.2012 -0.0835 0.1426 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:38] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:38] 0.2833 0.4975 0.2664 0.0196 0.1216 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:38] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:38] 0.4582 0.2047 0.0684 0.0824 0.3601 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:38] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:38] -0.4595 -0.2226 -0.0806 -0.0709 -0.3425 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:38] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:38] 0.29 0.1581 -0.0374 0.2102 0.0397 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:38] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:38] -0.2906 -0.1527 0.0495 -0.207 -0.0383 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:38] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:38] -0.0952 0.1398 0.2743 -0.2438 -0.1439 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:38] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:38] 0.2442 0.4969 0.2668 0.0324 0.1193 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:38] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 39 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-08-31 23:33:41 
## *** Done! 
## *** End time: 2019-08-31 23:58:47 
## Time difference of 25.1137 mins
## List of 10
##  $ : num [1:230223, 1:39] 1.068 0.738 1.071 -1.315 -0.896 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:39] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:39] 0.1307 0.1917 0.2659 0.0422 -0.1875 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:39] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:39] -0.32393 -0.52116 -0.29004 -0.00327 -0.15676 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:39] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:39] -0.986 -0.73 -1.038 1.324 0.93 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:39] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:39] -0.4284 -0.2182 -0.0771 -0.0924 -0.3365 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:39] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:39] 0.32797 0.52107 0.28561 -0.00291 0.13858 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:39] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:39] -0.4613 -0.2029 -0.0884 -0.069 -0.3912 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:39] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:39] 1.001 0.729 1.04 -1.322 -0.927 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:39] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:39] 0.1866 0.1245 -0.1048 0.253 -0.0302 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:39] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:39] -0.1549 -0.1397 0.1238 -0.2037 -0.0749 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:39] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 40 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 00:00:03 
## *** Done! 
## *** End time: 2019-09-01 00:30:28 
## Time difference of 30.41466 mins
## List of 10
##  $ : num [1:230223, 1:40] -0.462 -0.2186 -0.0861 -0.0593 -0.3842 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:40] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:40] -0.0971 -0.2119 -0.2744 -0.0318 0.0674 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:40] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:40] 0.36207 0.50792 0.28863 0.00972 0.13294 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:40] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:40] -0.3396 -0.5235 -0.2792 -0.0106 -0.1071 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:40] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:40] -0.2685 -0.2664 -0.2053 -0.0583 0.0682 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:40] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:40] 0.4565 0.221 0.0865 0.0572 0.3862 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:40] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:40] 0.466 0.206 0.074 0.074 0.389 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:40] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:40] -0.4563 -0.2124 -0.0867 -0.065 -0.3948 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:40] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:40] 0.4588 0.2185 0.0833 0.0633 0.382 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:40] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:40] 0.4601 0.1994 0.076 0.0723 0.3969 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:40] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 41 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 00:33:08 
## *** Done! 
## *** End time: 2019-09-01 01:04:49 
## Time difference of 31.68772 mins
## List of 10
##  $ : num [1:230223, 1:41] 0.972 0.728 1.035 -1.319 -0.957 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:41] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:41] 0.3594 0.5056 0.2833 0.0134 0.0567 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:41] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:41] -0.976 -0.732 -1.029 1.326 0.941 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:41] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:41] -0.1679 -0.1523 0.0964 -0.2425 -0.137 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:41] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:41] -0.4065 -0.2423 -0.088 -0.0822 -0.4149 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:41] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:41] 0.3459 0.5016 0.2814 0.0121 0.0477 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:41] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:41] -0.1726 -0.1605 0.0794 -0.2571 -0.1087 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:41] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:41] 0.29696 0.51709 0.26353 0.00787 0.01515 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:41] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:41] 0.3365 0.4951 0.2836 0.0145 0.054 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:41] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:41] -0.2399 -0.2445 -0.1902 -0.0593 0.1404 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:41] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 42 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 01:08:03 
## *** Done! 
## *** End time: 2019-09-01 01:42:34 
## Time difference of 34.51695 mins
## List of 10
##  $ : num [1:230223, 1:42] -0.21954 -0.18068 -0.1707 -0.00677 0.15092 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:42] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:42] -0.3237 -0.4746 -0.2516 0.0222 -0.0371 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:42] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:42] 0.4938 0.2293 0.0786 0.0744 0.4083 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:42] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:42] -0.4752 -0.2529 -0.0994 -0.09 -0.3838 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:42] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:42] -0.2352 -0.1865 -0.16222 -0.00775 0.12945 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:42] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:42] 0.3154 0.4635 0.248 -0.0164 0.0355 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:42] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:42] 0.352 0.4477 0.2593 -0.0218 0.11 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:42] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:42] -0.4744 -0.2315 -0.0931 -0.0752 -0.4053 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:42] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:42] -0.3136 -0.4747 -0.2482 0.0239 -0.0243 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:42] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:42] -0.3144 -0.473 -0.2375 0.0208 -0.0217 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:42] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 43 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 01:46:42 
## *** Done! 
## *** End time: 2019-09-01 02:22:39 
## Time difference of 35.95354 mins
## List of 10
##  $ : num [1:230223, 1:43] 0.964 0.718 1.02 -1.328 -0.937 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:43] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:43] -2.52e-01 -1.89e-01 -1.65e-01 9.05e-05 1.21e-01 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:43] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:43] 0.4557 0.2745 0.1174 0.0779 0.3588 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:43] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:43] -0.454 -0.308 -0.1183 -0.0896 -0.3149 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:43] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:43] 0.379 0.012 0.355 0.308 0.092 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:43] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:43] 0.2309 0.17942 0.16302 0.00762 -0.10973 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:43] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:43] -0.466 -0.2815 -0.1095 -0.0801 -0.3463 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:43] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:43] 0.3022 0.488 0.2542 -0.0216 -0.012 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:43] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:43] -0.3772 -0.4754 -0.2869 0.0118 -0.0227 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:43] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:43] 0.236109 0.189235 0.163614 -0.000282 -0.115437 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:43] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 44 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 02:27:00 
## *** Done! 
## *** End time: 2019-09-01 03:01:47 
## Time difference of 34.7858 mins
## List of 10
##  $ : num [1:230223, 1:44] -0.977 -0.712 -1.031 1.317 0.93 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:44] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:44] -0.3107 -0.4418 -0.2553 0.0401 -0.0796 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:44] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:44] 0.3963 0.4184 0.2818 -0.0464 0.1246 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:44] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:44] 0.3727 0.4133 0.2777 -0.0508 0.1338 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:44] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:44] 0.1762 0.1756 -0.0924 0.282 0.0874 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:44] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:44] -0.3566 -0.4395 -0.2546 0.0493 -0.0858 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:44] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:44] 0.465 0.299 0.115 0.118 0.27 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:44] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:44] 0.455 0.324 0.122 0.116 0.265 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:44] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:44] 0.3295 0.445 0.2493 -0.0477 0.074 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:44] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:44] 0.3394 0.4469 0.2515 -0.0521 0.1 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:44] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 45 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 03:05:44 
## *** Done! 
## *** End time: 2019-09-01 03:48:00 
## Time difference of 42.2618 mins
## List of 10
##  $ : num [1:230223, 1:45] 0.1 0.2135 0.2835 0.0894 -0.0264 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:45] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:45] 0.3117 0.4395 0.2554 -0.0403 0.0969 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:45] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:45] -0.42 -0.303 -0.139 -0.126 -0.312 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:45] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:45] 0.386 0.301 0.141 0.106 0.323 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:45] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:45] -0.383 -0.419 -0.279 0.048 -0.129 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:45] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:45] 0.2978 0.0189 0.358 0.3031 0.0939 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:45] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:45] 0.387 0.307 0.135 0.122 0.336 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:45] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:45] -0.331 -0.44 -0.255 0.053 -0.106 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:45] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:45] 0.3209 0.4401 0.2659 -0.0378 0.1025 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:45] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:45] 0.2101 0.1524 0.1714 -0.0233 0.0211 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:45] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 46 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 03:52:03 
## *** Done! 
## *** End time: 2019-09-01 04:57:01 
## Time difference of 1.082895 hours
## List of 10
##  $ : num [1:230223, 1:46] 0.1861 0.1754 0.1667 -0.0552 0.0427 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:46] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:46] -0.1652 -0.1614 -0.1727 0.0541 -0.0477 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:46] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:46] -0.1812 -0.1783 -0.1695 0.0531 -0.0348 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:46] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:46] -0.3479 -0.422 -0.2668 0.0637 -0.0968 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:46] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:46] -0.1773 -0.2009 0.0572 -0.2438 -0.0747 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:46] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:46] -0.409 -0.288 -0.151 -0.148 -0.286 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:46] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:46] 0.401 0.282 0.149 0.141 0.316 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:46] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:46] -0.983 -0.709 -1.034 1.312 0.93 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:46] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:46] 0.173 0.1788 -0.0644 0.2712 0.1457 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:46] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:46] -0.3198 -0.4447 -0.2622 0.0505 -0.0859 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:46] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 47 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 05:00:51 
## *** Done! 
## *** End time: 2019-09-01 05:53:45 
## Time difference of 52.89993 mins
## List of 10
##  $ : num [1:230223, 1:47] 0.399 0.275 0.159 0.149 0.293 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:47] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:47] 0.458 0.295 0.144 0.152 0.24 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:47] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:47] 0.404 0.269 0.147 0.169 0.276 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:47] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:47] 0.1276 0.1598 0.161 -0.0462 0.0338 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:47] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:47] 0.341 0.278 0.149 0.117 0.301 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:47] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:47] -0.2114 -0.0587 -0.1833 0.3023 0.0845 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:47] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:47] -0.2038 -0.1863 0.0607 -0.2491 -0.1491 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:47] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:47] -0.2135 -0.296 0.204 0.0734 -0.6711 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:47] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:47] -0.435 -0.274 -0.145 -0.151 -0.301 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:47] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:47] 0.956 0.707 1.034 -1.308 -0.937 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:47] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 48 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 05:56:53 
## *** Done! 
## *** End time: 2019-09-01 06:37:48 
## Time difference of 40.90832 mins
## List of 10
##  $ : num [1:230223, 1:48] -0.1535 -0.1804 -0.3095 -0.1046 0.0634 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:48] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:48] 0.2823 0.4528 0.2421 -0.0507 0.1079 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:48] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:48] -0.1916 -0.1604 0.0657 -0.2621 -0.1442 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:48] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:48] 0.285 0.4426 0.2435 -0.0568 0.109 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:48] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:48] 0.398 0.281 0.143 0.134 0.318 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:48] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:48] 0.386 0.277 0.141 0.126 0.33 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:48] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:48] 0.358 0.264 0.14 0.135 0.316 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:48] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:48] -0.396 -0.291 -0.139 -0.131 -0.296 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:48] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:48] 0.399 0.287 0.139 0.137 0.304 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:48] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:48] 0.3046 0.4368 0.2465 -0.0555 0.12 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:48] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 49 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 06:41:13 
## *** Done! 
## *** End time: 2019-09-01 07:34:04 
## Time difference of 52.84757 mins
## List of 10
##  $ : num [1:230223, 1:49] 0.2938 0.4851 0.2202 -0.064 -0.0172 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:49] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:49] -0.376 -0.26 -0.147 -0.136 -0.306 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:49] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:49] 0.14 0.267 -0.171 -0.113 0.659 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:49] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:49] -0.30009 -0.48113 -0.23563 0.06234 0.00907 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:49] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:49] 0.3138 0.462 0.2594 -0.0499 0.0251 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:49] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:49] 0.2962 0.4844 0.2307 -0.0638 -0.0192 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:49] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:49] -0.2776 -0.4561 -0.2253 0.0228 -0.058 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:49] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:49] -0.3397 -0.4562 -0.2567 0.0549 -0.0136 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:49] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:49] 0.325 0.4687 0.2507 -0.0654 0.0283 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:49] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:49] -0.399 -0.275 -0.13 -0.149 -0.312 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:49] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
## *** Starting parallel calculation on 10 core(s)...
## *** System: unix 
## *** 50 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 07:37:21 
## *** Done! 
## *** End time: 2019-09-01 08:29:20 
## Time difference of 51.98905 mins
## List of 10
##  $ : num [1:230223, 1:50] -0.968 -0.688 -1.039 1.313 0.856 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:50] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:50] 0.1141 0.1952 0.1467 -0.0504 -0.0641 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:50] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:50] 0.974 0.671 1.033 -1.333 -0.83 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:50] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:50] 0.3437 0.4466 0.2527 -0.0678 0.0555 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:50] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:50] -0.331 -0.245 -0.164 -0.135 -0.344 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:50] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:50] -0.354 -0.255 -0.162 -0.126 -0.317 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:50] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:50] -0.2041 -0.21 0.0733 -0.2768 -0.0366 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:50] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:50] 0.2398 0.4473 0.2407 -0.0505 0.0718 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:50] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:50] 0.11538 0.20096 0.14874 -0.05776 -0.00564 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:50] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:50] -0.319 -0.25 -0.164 -0.131 -0.347 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:50] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
```

```
## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion
```

```
## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion

## Warning in removeFactor(rnb.set, fact = conf.factor, ncomp = ncomp, ntry =
## ntry, : NAs introduced by coercion
```

```
## Warning in if (!fact %in% names(Var)) {: the condition has length > 1 and
## only the first element will be used
```

```
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 22 components, 10 runs, 230223 features, 461 samples.
## *** Start time: 2019-09-01 08:33:27 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 10 
## try # 2 of 10 
## try # 3 of 10 
## try # 4 of 10 
## try # 5 of 10 
## try # 6 of 10 
## try # 7 of 10 
## try # 8 of 10 
## try # 9 of 10 
## try # 10 of 10 
## *** Done! 
## *** End time: 2019-09-01 09:08:47 
## Time difference of 35.32122 mins
## List of 10
##  $ : num [1:230223, 1:22] -0.3499 -0.0128 -0.0793 -0.0784 -0.2112 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:230223] "cg24669183" "cg15560884" "cg17505339" "cg23803172" ...
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] -0.184 -0.5453 -0.1877 -0.075 0.0772 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:230223] "cg24669183" "cg15560884" "cg17505339" "cg23803172" ...
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] -0.3727 -0.0823 -0.1833 -0.3498 -0.2355 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:230223] "cg24669183" "cg15560884" "cg17505339" "cg23803172" ...
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] 1.09 0.728 1.082 -1.282 -0.948 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:230223] "cg24669183" "cg15560884" "cg17505339" "cg23803172" ...
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] -0.1876 -0.5532 -0.188 -0.0745 0.0799 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:230223] "cg24669183" "cg15560884" "cg17505339" "cg23803172" ...
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] 0.3393 0.0141 0.0818 0.0763 0.2242 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:230223] "cg24669183" "cg15560884" "cg17505339" "cg23803172" ...
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] -1.083 -0.734 -1.085 1.282 0.943 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:230223] "cg24669183" "cg15560884" "cg17505339" "cg23803172" ...
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] -0.19 -0.548 -0.188 -0.078 0.075 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:230223] "cg24669183" "cg15560884" "cg17505339" "cg23803172" ...
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] -0.1838 -0.5455 -0.1873 -0.0785 0.0764 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:230223] "cg24669183" "cg15560884" "cg17505339" "cg23803172" ...
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
##  $ : num [1:230223, 1:22] 0.3485 0.0259 0.0819 0.0589 0.2188 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:230223] "cg24669183" "cg15560884" "cg17505339" "cg23803172" ...
##   .. ..$ : chr [1:22] "ic.1" "ic.2" "ic.3" "ic.4" ...
## NULL
## Calculate ||X-SxM|| and r2 between component weights
## Correlate rows of S between tries
## Build consensus ICA
## Analyse stability
##                  ic.1       ic.2        ic.3       ic.4        ic.5
## cg24669183 -0.1547295 -0.4320545 -0.36317044 -1.0517664 -0.10622529
## cg15560884 -0.2803537 -0.2340899 -0.07078829 -0.7298743 -0.49905708
## cg17505339 -0.2240179 -0.1341435 -0.08727780 -1.0987710 -0.13401748
## cg23803172 -0.3056196 -0.1289506 -0.04617918  1.3020026 -0.15580314
## cg03128332 -0.1398510  0.1130155 -0.18720620  0.9283063  0.08714486
## cg18147296 -1.6905931 -0.5789715 -0.78820957 -0.5807884 -0.60709497
##                  ic.6        ic.7        ic.8        ic.9       ic.10
## cg24669183 -0.1783935  0.06509611  0.14438946 -0.29217952 -0.34331712
## cg15560884  0.4404296 -0.13667746  0.35753764 -0.39096648 -0.45617591
## cg17505339  0.2722434 -0.26222428  0.21347301 -0.04882235 -0.18721375
## cg23803172  0.1678416 -0.05989675 -0.13081278  0.13322699 -0.03437058
## cg03128332 -0.5379390  0.05229467 -0.03190098 -0.14342388 -0.07800132
## cg18147296  2.3135109 -1.01537300  0.64742470  0.45559358  0.50159555
##                  ic.11       ic.12       ic.13       ic.14        ic.15
## cg24669183  0.23692564 -0.43425379  0.64225592 -0.38245636 -0.115451367
## cg15560884 -0.06653967 -0.16225547 -0.35547960 -0.12517927 -0.197966854
## cg17505339 -0.02677486 -0.09945771 -0.01999728  0.18098922 -0.007959121
## cg23803172  0.14977805 -0.19536945 -0.06586736  0.03921036 -0.285005272
## cg03128332 -0.07236817 -0.20297714 -0.14952277 -0.25101716 -0.293246754
## cg18147296  0.41464793 -0.85554186  0.95992154  0.57862940 -0.003342885
##                 ic.16       ic.17       ic.18       ic.19      ic.20
## cg24669183 -0.4726277 -0.05541974  0.78027031  0.90933223 -0.5235399
## cg15560884 -1.7030178 -0.64826833 -0.07684828 -0.22861381 -0.2344407
## cg17505339  0.1191328 -0.04534413  0.19519776 -0.05027251  0.8514392
## cg23803172 -0.4812601 -0.11390001 -0.13434776  0.05859110  0.9527040
## cg03128332  1.6964959  0.45482826  0.12991987  0.27766816  1.1670942
## cg18147296 -0.4497652 -0.32860006 -1.38632343 -1.38219229 -1.2054599
##                  ic.21      ic.22
## cg24669183 -0.07010370  0.7665255
## cg15560884 -1.21474502 -0.3943747
## cg17505339  0.09684685 -1.1371824
## cg23803172  0.02908491  0.8685530
## cg03128332  0.46050049  0.8473540
## cg18147296  0.65591034  1.0148547
```

```
## Loading required package: pheatmap
```

```
## Loading required package: vioplot
```

```
## Loading required package: sm
```

```
## Warning in fun(libname, pkgname): couldn't connect to display "localhost:
## 10.0"
```

```
## Package 'sm', version 2.2-5.6: type help(sm) for summary information
```

```
## 
## Attaching package: 'sm'
```

```
## The following object is masked from 'package:MASS':
## 
##     muscle
```

```
## Loading required package: zoo
```

```
## 
## Attaching package: 'zoo'
```

```
## The following objects are masked from 'package:base':
## 
##     as.Date, as.Date.numeric
```

```
## Working with component # 1
```

```
## [1] "ANOVA"
```

```
## Working with component # 2
```

```
## [1] "ANOVA"
```

```
## Working with component # 3
```

```
## [1] "ANOVA"
```

```
## Working with component # 4
```

```
## [1] "ANOVA"
```

```
## Working with component # 5
```

```
## [1] "ANOVA"
```

```
## Working with component # 6
```

```
## [1] "ANOVA"
```

```
## Working with component # 7
```

```
## [1] "ANOVA"
```

```
## Working with component # 8
```

```
## [1] "ANOVA"
```

```
## Working with component # 9
```

```
## [1] "ANOVA"
```

```
## Working with component # 10
```

```
## [1] "ANOVA"
```

```
## Working with component # 11
```

```
## [1] "ANOVA"
```

```
## Working with component # 12
```

```
## [1] "ANOVA"
```

```
## Working with component # 13
```

```
## [1] "ANOVA"
```

```
## Working with component # 14
```

```
## [1] "ANOVA"
```

```
## Working with component # 15
```

```
## [1] "ANOVA"
```

```
## Working with component # 16
```

```
## [1] "ANOVA"
```

```
## Working with component # 17
```

```
## [1] "ANOVA"
```

```
## Working with component # 18
```

```
## [1] "ANOVA"
```

```
## Working with component # 19
```

```
## [1] "ANOVA"
```

```
## Working with component # 20
```

```
## [1] "ANOVA"
```

```
## Working with component # 21
```

```
## [1] "ANOVA"
```

```
## Working with component # 22
```

```
## [1] "ANOVA"
```

![plot of chunk data_preparation](figure/data_preparation-1.png)

```
## No id variables; using all as measure variables
```

```
## Saving 8.3 x 11.7 in image
```

```r
names(data.prep)
```

```
## [1] "quality.filter"   "annot.filter"     "total.filter"    
## [4] "rnb.set.filtered" "info"
```

### Selecting informative features (CpGs)

The next, crucial, step is selecting a subset of sites that are informative about the cell type composition of your sample. This can be done in various ways, and DecompPipeline provides a list of them through the ```prepare_CG_subsets``` function. However, we focus on a single option, which is typically employed in epigenomic studies: selecting the most variable sites across the samples. Since many sites are constant for all samples, focusing on the ones that show the highest variablity across the samples is sensible. We assume that we do not lose information by not considering those sites that do not vary at all. Here, we focus on the 5,000 most variable sites.


```r
cg_subset <- prepare_CG_subsets(rnb.set=data.prep$rnb.set.filtered,
                                MARKER_SELECTION = "var",
                                N_MARKERS = 5000)
names(cg_subset)
```

```
## [1] "var"
```

## Methylome Deconvolution

### Performing Deconvolution

In this step, the actual deconvolution experiment is performed. There are different approaches, which are conceptually similar, yet different in their performance, running time and robustness. Among others, EDec, RefFreeCellMix from the RefFreeEWAS package and MeDeCom can be used to execute non-negative matrix factorization on your data. This will lead to two matrices, the proportions matrix of potential cell types (here referred to as LMCs) and the matrix of those pure profiles. We here focus on MeDeCom as the Deconvolution tool, although DecompPipeline also supports RefFreeCellMix and EDec.


```r
md.res <- start_medecom_analysis(
  rnb.set=data.prep$rnb.set.filtered,
  cg_groups = cg_subset,
  Ks=2:15,
  LAMBDA_GRID = c(0,10^-(2:5)),
  factorviz.outputs = T,
  analysis.name = "TCGA_LUAD2",
  cores = 15
)
```

```
## Loading required package: Rcpp
```

```
## Loading required package: pracma
```

```
## 
## Attaching package: 'pracma'
```

```
## The following object is masked from 'package:sm':
## 
##     nile
```

```
## The following object is masked from 'package:ff':
## 
##     quad
```

```
## The following object is masked from 'package:bit':
## 
##     is.sorted
```

```
## Loading required package: gtools
```

```
## 
## Attaching package: 'gtools'
```

```
## The following object is masked from 'package:pracma':
## 
##     logit
```

```
## The following object is masked from 'package:R.utils':
## 
##     capture
```

```
## Loading required package: RUnit
```

```
## Warning: replacing previous import 'gtools::logit' by 'pracma::logit' when
## loading 'MeDeCom'
```

```
## [1] "Did not write the variable dump: should only be executed from an environment with all the variables set"
## [2019-09-01 09:47:26, Main:] checking inputs
## [2019-09-01 09:47:26, Main:] preparing data
## [2019-09-01 09:47:26, Main:] preparing jobs
## [2019-09-01 09:47:26, Main:] 3570 factorization runs in total
## [2019-09-03 15:29:49, Main:] finished all jobs. Creating the object
```

## Downstream analysis

After performing deconvolution, results need to be visualized and interpreted. Most notably, the contribution matrix can be linked to phenotypic information about the samples to indicate different cellular compositions of the groups and the LMC matrix can be used to determine what the components represent. For visualization and downstream analysis, we use FactorViz. LOLA or GO enrichment analysis can be employed on sites that are specifically methylated/unmethylated in one of the LMCs.


```r
suppressPackageStartupMessages(library(FactorViz))
startFactorViz(file.path(getwd(),"TCGA_LUAD2","FactorViz_outputs"))
```

# References
