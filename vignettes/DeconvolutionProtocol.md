---
title: "DNA Methylation Deconvolution Protocol"
author: Michael Scherer, Petr Nazarov, Reka Toth, Shashwat Sahay, Tony Kamoa, Valentin
  Maurer, Nikita Vedenev, Christoph Plass, Thomas Lengauer, Joern Walter, and Pavlo Lutsik
date: "May 03, 2020"
output:
  html_document: default
  pdf_document: default
bibliography: bibliography.bib
---



# Introduction

This protocol aims at guiding reasearcher how to employ deconvolution of methylomes obtained from complex tissues. It will start with data retrieval from a public resource, but is equally applicable to in-house generated data. We will furthermore focus on the Illumina BeadChip series as a data source, although the protocol is also compatible with bisulfite sequencing, or any technology yielding single base pair resolution.
Deconvolution here refers to creating two matrices (proportion matrix A and methylation pattern matrix T) from a single matrix of input DNA methylation data (dimension CpGs x samples). Non-negative matrix factorization is employed for this task, and we will discuss some of the advantages and caveats of the methods.

# Installation (duration up to 2 h)

1. If R is not yet installed, follow the instructions at https://cran.r-project.org/. Create a working directory on a filesystem partition with sufficient storage capacity. Throughout the analysis, ~20-30 Gb of free storage will be required. Be sure that all the files are downloaded into this working directory and that the code is executed within the directory.

2a. **Linux** To execute the protocol described here, several R-packages (*RnBeads*, *DecompPipeline*, *MeDeCom*, and *FactorViz*) are required. First, you should have a recent R version installed and all your packages updated. The installation of the software is appliable for most *Linux distributions* directly through R, while for *MacOS* the binary release of *MeDeCom* https://github.com/lutsik/MeDeCom/releases/download/v0.3.0/MeDeCom_0.3.0.tgz should be used. The protocol is also available as a Docker container https://hub.docker.com/r/mscherer/medecom, which is currently the only option for Windows operating systems. For the remainder of this protocol, we assume a common Linux distribution as the operating system.


```r
install.packages(c("devtools","BiocManager"))
BiocManager::install(c("RnBeads","RnBeads.hg19","RnBeads.mm10","RnBeads.hg38"),dependencies=TRUE)
devtools::install_github(c("lutsik/MeDeCom","CompEpigen/DecompPipeline","CompEpigen/FactorViz"))
library(DecompPipeline)
```

2b. **macOS** In a first step, you will have to retrieve the binary version of *MeDeCom* from GitHub. Afterwards, you can install the remaining packages similar to the *Linux* case.


```r
install.packages(c("devtools","BiocManager"))
BiocManager::install(c("RnBeads","RnBeads.hg19","RnBeads.mm10","RnBeads.hg38"),dependencies=TRUE)
install.packages("https://github.com/lutsik/MeDeCom/releases/download/v1.0.0/MeDeCom_1.0.0.tgz",repos=NULL,type="binary")
devtools::install_github(c("CompEpigen/DecompPipeline","CompEpigen/FactorViz"))
library(DecompPipeline)
```

2c. **Windows** A dedicated branch in the GitHub repository handling installation on Windows machines is available and can be installed similar to the *Linux* case. 


```r
BiocManager::install(c("RnBeads",
        "RnBeads.hg19",
        "RnBeads.mm10",
        "RnBeads.hg38"),
   dependencies=TRUE)
library(devtools)
devtools::install_github("lutsik/MeDeCom",ref="windows")
devtools::install_github(
    c("CompEpigen/DecompPipeline","CompEpigen/FactorViz")
)
library(DecompPipeline)
```

2 d. **Cross-platform using Docker** If not yet available, install Docker on your machine. The [Docker website](https://docs.docker.com) provides detailed installation instructions for all major operating systems and computational environments, including Linux, Windows and MacOS. Specifically, on Windows follow [these](https://docs.docker.com/docker-for-windows/install/) instructions. After successful installation, start Docker desktop, be sure that a X Server is running, open a new PowerShell window and type:


```bash
# Windows
docker run  -e DISPLAY=<YOUR_IP>:0.0 -it mscherer/medecom

# Linux
xhost +"local:docker@"
docker run --env DISPLAY=$DISPLAY --privileged --volume $XAUTH:/root/.Xauthority --network=host --volume /tmp/.X11-unix:/tmp/.X11-unix --rm --rm -it mscherer/medecom
```

An interactive R-session will start which can be used to run this protocol.

# Protocol

## Data Retrival

### Obtaining data from a public resource (duration ~5 h)

3. Use the Genomic Data Commons (GDC) data download tool (https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) to download the IDAT files with the TCGA-LUAD manifest file (http://epigenomics.dkfz.de/DecompProtocol/data/gdc_manifest.2019-01-23.txt) and its associated metadata by running the download tool in a terminal session initiated in the working directory:


```bash
gdc-client download –m gdc_manifest.2019-01-23.txt
```

Store the resulting files in a new directory *idat*.

4. Retrieve the clinical metadata from https://portal.gdc.cancer.gov/projects/TCGA-LUAD by clicking on the “Clinical” button and unpacking the downloaded *.tar.gz* file into a new subdirectory annotation within your working directory. The remaining parsing is performed through an R session initiated in the same directory:


```r
library(tidyverse)
clinical.data <- read.table(
	"annotation/clinical.tsv",
	sep="\t", header = TRUE)
idat.files <- list.files("idat",full.names = TRUE)
meta.files <- list.files(idat.files[1],full.names = TRUE)
untar(meta.files[grepl(".tar.gz",meta.files)],exdir = idat.files[1])
meta.files <- untar(meta.files[grepl(".tar.gz",meta.files)], list = TRUE)

meta.info <- read.table(
	file.path(idat.files[1],meta.files[grepl(".sdrf.txt",meta.files)]),
	sep="\t",header=TRUE)

meta.info <- meta.info[
	match(unique(meta.info$Comment..TCGA.Barcode.),
	        meta.info$Comment..TCGA.Barcode.
		),
]
anno.frame <- meta.info %>% 
	mutate("submitter_id"=substr(Comment..TCGA.Barcode., 1, 12)) %>%
	inner_join(x=clinical.data, by=c("submitter_id"="submitter_id")) %>%
	distinct(submitter_id, submitter_id.4, .keep_all = TRUE) %>%
	mutate(barcode=gsub("_Grn.idat|_Red.idat", "", Array.Data.File),
		Sentrix_ID=sub("_.*$","", barcode),
		Sentrix_Position=sub("^[^_]+_","",barcode),
		healthy_cancer=ifelse(grepl("11A",Comment..TCGA.Barcode.),
			"healthy","cancer"))

write.table(anno.frame,
	"annotation/sample_annotation.tsv",
	quote = FALSE, row.names = FALSE, sep= "\t")
```

5. Copy the IDAT files into a single directory idat for downstream analysis.


```r
#' write idat files to parent directory
lapply(idat.files,function(x){
  is.idat <- list.files(x,pattern = ".idat", full.names = TRUE)
  file.copy(is.idat,"idat/")
  unlink(x,recursive = TRUE)
})
```

## Data Processing

### Data import (~2 h)

6. *RnBeads* converts the files into a data object and performs basic quality control. Analysis options have to be specified for *RnBeads*, either through an XML file, or through the command line. Deactivate the preprocessing, exploratory, covariate inference, export and differential methylation modules, such that *RnBeads* only performs data import and quality control.


```r
suppressPackageStartupMessages(library(RnBeads))
rnb.options(
  assembly = "hg19",
  identifiers.column = "submitter_id",
  import = TRUE,
  import.default.data.type = "idat.dir",
  import.table.separator = "\t",
  import.sex.prediction = TRUE,
  qc = TRUE,
  preprocessing = FALSE,
  exploratory = FALSE,
  inference = FALSE,
  differential = FALSE,
  export.to.bed = FALSE,
  export.to.trackhub = NULL,
  export.to.csv = FALSE
)
```

7. Specify the input to *RnBeads*: the sample annotation sheet created at the data retrieval step, the
directory in which the IDAT files are stored and a directory to which the HTML report is to be saved.
Additionally, specify a temporary directory and start the RnBeads analysis.


```r
sample.anno <- "annotation/sample_annotation.tsv"
idat.folder <- "idat/"
dir.report <- paste0("report",Sys.Date(),"/")
temp.dir <- "/tmp"
options(fftempdir=temp.dir)
rnb.set <- rnb.run.analysis(dir.reports = dir.report, 
	sample.sheet = sample.anno,
	data.dir = idat.folder)
```

**PAUSE POINT** We provide an example report containing the processed RnBSet object for further analysis. It can be obtained from [http://epigenomics.dkfz.de/downloads/DecompProtocol/RnBeads_Report_TCGA_LUAD/index.html](http://epigenomics.dkfz.de/downloads/DecompProtocol/RnBeads_Report_TCGA_LUAD/index.html) or directly loaded via:


```r
download.file("http://epigenomics.dkfz.de/downloads/DecompProtocol/rnbSet_unnormalized.zip", destfile=paste0(tempdir(),"/RnBSet"))
rnb.set <- load.rnb.set(paste0(tempdir(),"/RnBSet"))
```

8. *RnBeads* creates an interactive HTML report, specifying the steps performed and the associated results. Data was of good quality such that it can be used for further analysis. 

**PAUSE POINT** For reference, we provide a complete RnBeads report on the supplementary website
http://epigenomics.dkfz.de/downloads/DecompProtocol/RnBeads_Report_TCGA_LUAD/.

### Preprocessing and Filtering (22 h)

9. For further data preparation and analysis steps use the *DecompPipeline* package. Processing options are provided through individual function parameters. Follow a stringent filtering strategy: (i) Filter CpGs covered by less than 3 beads, and probes that are in the 0.001 and 0.999 overall intensity quantiles (low and high intensity outliers). (ii) Remove all probes containing missing values in any of the samples. (iii) Discard sites outside of CpG context, overlapping annotated SNPs, located on the sex chromosomes and potentially cross-reactive probes. Finally, apply BMIQ normalization to account for the bias introduced by the two Infinium probe designs.


```r
suppressPackageStartupMessages(library(DecompPipeline))
data.prep <- prepare.data(rnb.set = rnb.set,
	analysis.name = "TCGA_LUAD",
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
	execute.lump=TRUE)
```

**PAUSE POINT** We provide a list of CpGs that are selected for the downstream analysis at http://epigenomics.dkfz.de/downloads/DecompProtocol/sites_passing_complete_filtering.csv. You can directly select those sites from you RnBSet object.


```r
rem.sites <- rep(TRUE,nrow(meth(rnb.set)))
sel.sites <-
read.csv("http://epigenomics.dkfz.de/downloads/DecompProtocol/sites_passing_co
mplete_filtering.csv")
rem.sites[row.names(annotation(rnb.set))%in%as.character(sel.sites[,1])] <-
FALSE
rnb.set <- remove.sites(rnb.set,rem.sites)
data.prep <- list(rnb.set.filtered=rnb.set)
```

10. Adjustment for potential confounding factors is crucial in epigenomic studies and the influences
of, for instance, donor sex, age, and genotype on the DNA methylation pattern are well-studied 28,29 .
Use Independent Component Analysis (ICA, see Materials) to account for DNA methylation differences
that are due to sex, age, race, and ethnicity. Confounding factor adjustment alters the overall data distribution. Go to the analysis directory generated by DecompPipeline (*TCGA_LUAD*) and investigate the associations between Independent Components and the specified confounding factors. In case the reported p-value is lower than the parameter alpha.fact , ICA will automatically adjust for the corresponding confounding factor. Furthermore, inspect if the distribution of DNA methylation is preserved. Confounding factor adjustment is currently only supported for Illumina BeadArray datasets and has not yet been extended to bisulfite sequencing data.


```r
suppressPackageStartupMessages(library(DecompPipeline))
data.prep <- prepare.data(rnb.set = data.prep$rnb.set.filtered,
	normalization = "none",
	analysis.name = "TCGA_LUAD",
	filter.beads=FALSE,
	filter.intensity=FALSE,
	remove.ICA=TRUE,
	filter.na = FALSE,
	filter.context = FALSE,
	filter.snp = FALSE,
	filter.sex.chromosomes = FALSE,
	filter.cross.reactive = FALSE,
	conf.fact.ICA=c("age_at_diagnosis","race","gender","ethnicity"),
	ica.setting=c("alpha.fact"=1e-5,"save.report"=TRUE,
		"ntry"=10,"nmin"=20,"nmax"=50,"ncores"=10))
```
### Selection of CpG subsets (1 min)

11. Select a subset of sites to be used for deconvolution. DecompPipeline provides different options (see documentation of *prepare_CG_subsets*) through the *prepare_CG_subsets* function. Select the 5,000 most highly variable CpGs across the samples.


```r
cg_subset <- prepare.CG.subsets(rnb.set = data.prep$rnb.set.filtered,
                                marker.selection = "var",
                                n.markers = 5000)
names(cg_subset)
```

## Methylome Deconvolution

### Performing Deconvolution (54 h)

12. Perform the deconvolution experiment. Use MeDeCom with a grid of values for the number of components (K) ranging from 2 to 15, which covers homogeneous to heterogeneous samples. Also, specify a grid for the regularization parameter (λ) from strong (0.01) to no regularization (0). We here focus on *MeDeCom* as the Deconvolution tool, although *DecompPipeline* also supports *RefFreeCellMix* and *EDec*.


```r
md.res <- start_medecom_analysis(
  rnb.set=data.prep$rnb.set.filtered,
  cg_groups = cg_subset,
  Ks=2:15,
  lambda.grid = c(0,10^-(2:5)),
  factorviz.outputs = T,
  analysis.name = "TCGA_May",
  cores = 15
)
```
13. (Optional) Decide on the Ddeconvolution tool to be used of the DNA methylation data matrix is at
the core of this protocol, and different methods have been proposed. In addition to our own
method MeDeCom, this protocol supports *RefFreeCellMix* and *EDec*.

**PAUSE POINT** The final MeDeCom result is stored in a format that can be directly imported with FactorViz. We provide the final result for exploration at http://epigenomics.dkfz.de/downloads/DecompProtocol/FactorViz_outputs.tar.gz.


```r
download.file("http://epigenomics.dkfz.de/downloads/DecompProtocol/FactorViz_outputs.tar.gz",destfile=paste0(tempdir(),"/FactorViz_outputs.tar.gz"))
untar(paste0(tempdir(),"/FactorViz_outputs.tar.gz"))
```

## Downstream analysis (3 h)

14. Start the *FactorViz* application to visualize and interactively explore the deconvolution results.


```r
suppressPackageStartupMessages(library(FactorViz))
startFactorViz(file.path(getwd(),"TCGA_LUAD","FactorViz_outputs"))
```
15. (Optional). Load the MeDeComSet object in an additional R session started in the working
directory, where deconvolution has been performed. The object is stored in the FactorViz_outputs directory, and can be used for additional analysis.

16. Determine the number of LMCs (K) and the regularization parameter (λ) based on the cross-validation error. First, go to panel “K selection” to plot cross-validation error for the range of Ks specified earlier. We recommend selecting K values as a trade-off between too high complexity, i.e. fitting the noise in the data (higher K values) and insufficient degrees of freedom (lower K values). This typically corresponds to the saddle point where cross-validation error starts to level out.

17. For a fixed K, select a value for the regularization parameter (λ) by proceeding to panel “Lambda selection”. An optimal value for λ will often correspond to a local minimum of the cross-validation error. On the other hand, noticeable changes in other statistics, such as objective value or root mean squared error (RMSE), can point at a different λ value.

18. In the “Proportions” panel of FactorViz, visualize LMC proportions using heatmaps. Proportion heatmaps can be visually annotated with available qualitative and quantitative traits (see dropdown “Color samples by”). Further visualization options include stacked barplot and lineplots.

19. In the “Meta-analysis” panel, associate LMC proportions with quantitative and qualitative traits using correlation- and t-tests. Further sample annotations, such as mutational load (https://www.cbioportal.org/study/clinicalData?id=luad_tcga_pan_can_atlas_2018) and tumor purity scores can be obtained from public repositories and related publications (https://static-content.springer.com/esm/art%3A10.1038%2Fncomms3612/MediaObjects/41467_2013_BFncomms3612_MOESM489_ESM.xlsx). In addition, scripts to conduct such analysis are available on the Supplementary Website (http://epigenomics.dkfz.de/DecompProtocol/).

20. Explore the LMCs in the panel “LMCs”. Several visualization options are available such as Multidimensional Scaling and histograms. Reference profiles can be used for joint visualization in case those are available.

21. Determine DNA methylation sites that are specifically hypo- and hypermethylated in an LMC by comparing the methylation values in the LMC matrix for each LMC to the median of the remaining LMCs in the “Meta-analysis” panel. LMC-specific sites can be used for either a gene-centric Gene Ontology analysis (select ”Enrichments” and “GO Enrichments” in dropdown ”Analysis” and “Output type”) or for region-based enrichment analysis with the LOLA package (select “LOLA Enrichments” in dropdown ”Output type”). The differential sites can be exported for further analysis, and the resulting plots can be stored by clicking on the “PDF” button.

**PAUSE POINT** Additional built-in options for visualization of LMCs and proportions are available via the *MeDeCom* functions ```plotLMCs``` and ```plotProportions``` in R. Finally, the LMC and proportions matrix can be extracted from the MeDeComSet object, and the genomic annotation of CpGs (```ann.C```) can be obtained from the *FactorViz* output for a custom downstream analysis:


```r
LMCs <- getLMCs(medecom.set, K=7, lambda=0.001)
LMC.proportions <- getProportions(medecom.set, K=7, lambda=0.001)
load("FactorViz_outputs/ann_C.RData")
head(ann.C)
```
# References
