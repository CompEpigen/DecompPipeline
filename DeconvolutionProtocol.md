---
title: "DNA Methylation Deconvolution Protocol"
author: "Michael Scherer, Pavlo Lutsik, Petr Nazarov, Tony Kamoa"
date: "May 12, 2019"
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
## 2019-06-07 15:12:27     1.1  STATUS STARTED RnBeads Pipeline
## 2019-06-07 15:12:27     1.1    INFO     Initialized report index and saved to index.html
## 2019-06-07 15:12:28     1.1  STATUS     STARTED Loading Data
## 2019-06-07 15:12:28     1.1    INFO         Number of cores: 1
## 2019-06-07 15:12:28     1.1    INFO         Loading data of type "idat.dir"
## 2019-06-07 15:12:28     1.1  STATUS         STARTED Loading Data from IDAT Files
## 2019-06-07 15:12:30     1.1    INFO             Detected platform: HumanMethylation450
## 2019-06-07 15:34:36     1.5  STATUS         COMPLETED Loading Data from IDAT Files
## 2019-06-07 16:20:50     2.0  STATUS         Loaded data from idat/
## 2019-06-07 16:21:44     7.6  STATUS         Predicted sex for the loaded samples
## 2019-06-07 16:22:11     7.1  STATUS         Added data loading section to the report
## 2019-06-07 16:22:11     7.1  STATUS         Loaded 461 samples and 485577 sites
## 2019-06-07 16:22:11     7.1    INFO         Output object is of type RnBeadRawSet
## 2019-06-07 16:22:11     7.1  STATUS     COMPLETED Loading Data
## 2019-06-07 16:36:14     7.1    INFO     Initialized report index and saved to index.html
## 2019-06-07 16:36:14     7.1  STATUS     STARTED Quality Control
## 2019-06-07 16:36:14     7.1    INFO         Number of cores: 1
## 2019-06-07 16:36:14     7.1  STATUS         STARTED Quality Control Section
```

```
## 2019-06-07 16:36:49     2.0  STATUS             Added quality control box plots
```

```
## 2019-06-07 16:41:56     2.0  STATUS             Added quality control bar plots
```

```
## 2019-06-07 16:42:22     2.0  STATUS             Added negative control boxplots
## 2019-06-07 16:42:22     2.0  STATUS         COMPLETED Quality Control Section
## 2019-06-07 16:42:22     2.0  STATUS         STARTED Visualizing SNP Probe Data
## 2019-06-07 16:42:22     2.0  STATUS             STARTED Mixups Visualization Section
```

```
## 2019-06-07 16:42:56     5.4  STATUS                 Added SNP Heatmap
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
## 2019-06-07 16:42:56     5.4  STATUS                 Calculated Manhattan distances between samples based on SNP probes
```

```
## 2019-06-07 16:43:00     5.4  STATUS                 Added SNP-based Distances
## 2019-06-07 16:43:00     5.4  STATUS             COMPLETED Mixups Visualization Section
## 2019-06-07 16:43:00     5.4  STATUS         COMPLETED Visualizing SNP Probe Data
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ff7dc16ee9ef5b.ff
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ff7dc14516b2a.ff
```

```
## 2019-06-07 16:43:53     7.7  STATUS     COMPLETED Quality Control
## 2019-06-07 16:43:53     7.7    INFO     Initialized report index and saved to index.html
## 2019-06-07 16:43:53     7.7  STATUS     STARTED Saving RData
## 2019-06-07 16:43:53     7.7  STATUS     COMPLETED Saving RData
## 2019-06-07 16:43:53     7.7  STATUS COMPLETED RnBeads Pipeline
```

RnBeads creates an interactive HTML report, specifying the steps performed and the associated results. Data was of good quality such that is can be used for further analysis. (Include two screenshots from the RnBeads report)

### Preprocessing and Filtering

For further analysis, we use the DecompPipeline package (https://github.com/lutsik/DecompPipeline), which provides a comprehensive workflow including crucial data preparation steps for methylome deconvolution experiments. The options are provided through the individual function parameters. We follow a stringent filtering strategy. First, all samples having fewer than 3 beads covered are filtered, as well as those probes that are in the 0.05 and 0.95 overall intensity quantiles, respectively. We then remove all probes containing missing values, outside of CpG context, that overlap with annotated SNPs, on the sex chromosomes and probes that have been shown to be cross-reactive on the chip. Then, BMIQ normalization [@bmiq] is employed to account for the chip's design bias.
Accounting for potential confounding factor is crucial in epigenomic studies. Especially, the influence of donor sex on the DNA methylation pattern is well-studied and strong. We used Independent Component Analysis (ICA) to account for DNA methylation differences that are due to sex. ICA detects components in the data accounting for most of the variance similar to PCA, but does not require orthogonality of the components but statistical independence. We used an external library (http://sablab.net/scripts/LibICA.r) for performing ICA to adjust for sex.

```r
suppressPackageStartupMessages(library(DecompPipeline))
```

```
## Warning: replacing previous import 'gtools::logit' by 'pracma::logit' when
## loading 'MeDeCom'
```

```r
data.prep <- prepare_data(RNB_SET = rnb.set,
                          analysis.name = "TCGA_LUAD",
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
                          conf.fact.ICA="gender",
                          ica.setting=c("alpha.fact"=1e-5))
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ff7dc194821f0.ff
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ff7dc145b1d579.ff
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ff7dc159d07925.ff
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ff7dc1511ee463.ff
```

```
## opening ff /DEEP_fhgfs/projects/mscherer/data/450K/TCGA/LUAD/tmp/ff7dc14fe6b7f9.ff
```

```
## Loading required package: fastICA
```

```
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 10 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:35:03 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:36:19 
## Time difference of 1.257379 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 3.914845e-06
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 11 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:36:21 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:37:32 
## Time difference of 1.185299 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 4.618214e-07
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 12 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:37:34 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:38:56 
## Time difference of 1.353956 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 8.07913e-05
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 13 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:38:58 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:40:22 
## Time difference of 1.393692 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 8.253542e-05
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 14 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:40:24 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:41:49 
## Time difference of 1.416987 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 4.797991e-06
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 15 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:41:52 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:43:28 
## Time difference of 1.598963 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 9.294011e-07
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 16 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:43:31 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:45:17 
## Time difference of 1.776674 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 0.0004529832
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 17 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:45:20 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:47:20 
## Time difference of 2.007774 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 4.989633e-06
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 18 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:47:24 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:49:28 
## Time difference of 2.079608 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 2.363297e-07
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 19 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:49:32 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:51:43 
## Time difference of 2.190718 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 9.885586e-06
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 20 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:51:47 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:54:08 
## Time difference of 2.354854 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 2.200403e-05
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 21 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:54:11 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:57:04 
## Time difference of 2.873583 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 1.354768e-05
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 22 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 18:57:07 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 18:59:59 
## Time difference of 2.859485 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 9.786937e-08
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 23 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 19:00:03 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 19:03:23 
## Time difference of 3.345457 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 3.355455e-08
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 24 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 19:03:27 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 19:07:04 
## Time difference of 3.623443 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 5.450879e-07
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 25 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 19:07:09 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 19:09:55 
## Time difference of 2.773285 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 2.173957e-08
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 26 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 19:09:59 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 19:13:31 
## Time difference of 3.526244 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 9.057925e-09
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 27 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 19:13:35 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 19:17:44 
## Time difference of 4.145353 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 1.158838e-07
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 28 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 19:17:48 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 19:22:25 
## Time difference of 4.607401 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 1.894004e-08
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 29 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 19:22:29 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 19:26:46 
## Time difference of 4.288464 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 5.409245e-09
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 30 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 19:26:51 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 19:31:16 
## Time difference of 4.410147 mins
## Calculate ||X-SxM|| and r2 between component weights
## [1] 3.561972e-07
## *** Starting  calculation on 1 core(s)...
## *** System: unix 
## *** 29 components, 1 runs, 230223 features, 461 samples.
## *** Start time: 2019-06-07 19:31:44 
## Execute one-core analysis, showing progress every 1 run(s)
## try # 1 of 1 
## *** Done! 
## *** End time: 2019-06-07 19:36:05 
## Time difference of 4.33801 mins
## Calculate ||X-SxM|| and r2 between component weights
```

```
## No id variables; using all as measure variables
```

```
## Saving 7 x 7 in image
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
  analysis.name = "TCGA_LUAD",
  cores = 15
)
```

```
## [1] "Did not write the variable dump: should only be executed from an environment with all the variables set"
## [2019-06-07 20:08:43, Main:] checking inputs
## [2019-06-07 20:08:43, Main:] preparing data
## [2019-06-07 20:08:43, Main:] preparing jobs
## [2019-06-07 20:08:43, Main:] 3570 factorization runs in total
## [2019-06-09 09:42:13, Main:] finished all jobs. Creating the object
```

## Downstream analysis

After performing deconvolution, results need to be visualized and interpreted. Most notably, the contribution matrix can be linked to phenotypic information about the samples to indicate different cellular compositions of the groups and the LMC matrix can be used to determine what the components represent. For visualization and downstream analysis, we use FactorViz. LOLA or GO enrichment analysis can be employed on sites that are specifically methylated/unmethylated in one of the LMCs.


```r
suppressPackageStartupMessages(library(FactorViz))
startFactorViz(file.path(getwd(),"TCGA_LUAD","FactorViz_outputs"))
```

# References
