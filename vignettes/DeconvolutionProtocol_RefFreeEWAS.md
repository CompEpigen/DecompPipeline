---
title: "DNA Methylation Deconvolution Protocol - RefFreeEWAS"
author: Michael Scherer, Petr Nazarov, Reka Toth, Shashwat Sahay, Tony Kamoa, Valentin
  Maurer, Nikita Vedenev, Christoph Plass, Thomas Lengauer, Joern Walter, and Pavlo Lutsik
date: "May 03, 2020"
output:
  html_document: default
  pdf_document: default
bibliography: bibliography.bib
---

# Installation and Data Retrieval
Pleas follow the instructions provided in the [protocol](protocol.html) using MeDeCom.

# Data import
Analogously to the general protocol, we use [RnBeads](https://rnbeads.org) for data import.


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
temp.dir <- tempdir()
options(fftempdir=temp.dir)
rnb.set <- rnb.run.analysis(dir.reports = dir.report,
                            sample.sheet = sample.anno,
                            data.dir = idat.folder)
```

# Preprocessing and Filtering
The preprocessing and filtering is also performed in line with the general [protocol](protocol.html) using [*DecompPipeline*](https://github.com/CompEpigen/DecompPipeline).


```r
suppressPackageStartupMessages(library(DecompPipeline))
data.prep <- prepare.data(rnb.set = rnb.set,
                          analysis.name = "TCGA_RefFree",
                          normalization = "bmiq",
                          filter.beads = TRUE,
                          min.n.beads = 3,
                          filter.intensity = TRUE,
                          min.int.quant = 0.001,
                          max.int.quant = 0.999,
                          filter.na = TRUE,
                          filter.context = TRUE,
                          filter.snp = TRUE,
                          filter.sex.chromosomes = TRUE,
                          filter.cross.reactive = TRUE,
                          execute.lump = TRUE)
```

# Confounding factor adjustment
We use Independent Component Analysis (ICA) to adjust for the potential confounding factors age, race, gender, and ethnicity.


```r
data.prep <- prepare.data(rnb.set = data.prep$rnb.set.filtered,
			  analysis.name = "TCGA_RefFree",
			  normalization = "none",
		  	filter.beads = FALSE,
			  filter.intensity = F,
        filter.na = FALSE,
        filter.context = FALSE,
        filter.snp = FALSE,
        filter.sex.chromosomes = FALSE,
        filter.cross.reactive = FALSE,
        remove.ICA = TRUE,
        conf.fact.ICA = c("age_at_diagnosis","race","gender","ethnicity"),
        ica.setting = c("alpha.fact"=1e-5,"save.report"=TRUE,
                        "ntry"=10,"nmin"=20,"nmax"=50,"ncores"=10))
```

# Feature selection
We select the 5,000 most variable CpGs across the samples for downstream analysis.


```r
cg_subset <- prepare.CG.subsets(rnb.set=data.prep$rnb.set.filtered,
                                marker.selection = "var",
                                n.markers = 5000)
```

# Deconvolution using RefFreeCellMix

We use the *RefFreeCellMix* function from the [RefFreeEWAS](https://cran.r-project.org/web/packages/RefFreeEWAS/index.html) package for deconvolution analysis. In contrast to *MeDeCom*, *RefFreeCellMix* does not use a regularization for the entries of the matrix to be at the extreme values 0 or 1. Thus, we just specify the number of components to be tested and the preprocessed object, along with the selected features.


```r
md.res <- start.refreeewas.analysis(
  rnb.set=data.prep$rnb.set.filtered,
  cg_groups = cg_subset,
  Ks=2:15,
  factorviz.outputs=TRUE,
  work.dir = "TCGA_RefFree"
)
```

# Downstream analysis

The resulting object is an object of type *MeDeComSet*, since the output from *RefFreeCellMix* has internally been transformed to the *MeDeCom* infrastructure. This enables loading the data set directly into the visualization tool [*FactorViz*](https://github.com/CompEpigen/FactorViz)


```r
suppressPackageStartupMessages(library(FactorViz))
startFactorViz(file.path(getwd(),"TCGA_RefFree","FactorViz_outputs"))
```
