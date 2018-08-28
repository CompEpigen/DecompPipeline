suppressPackageStartupMessages(library(RnBeads))
library(MeDeCom)
options(fftempdir="/scratch/divanshu/dataset1/temp",disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)
rnb.options(disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)

rnb.set <- load.rnb.set("./2018-07-08_integrated_rnb.clean.RDS.zip")

source("/scratch/divanshu/newset/filtering.R")
res<-prepare_data(
  RNB_SET=rnb.set, 
  WORK_DIR="/scratch/divanshu/dataset1/temp",
  DATASET="adgbbSorted",
  DATA_SUBSET="frontal",
  SAMPLE_SELECTION_COL=NA,
  SAMPLE_SELECTION_GREP=NA,
  PHENO_COLUMNS=NA,
  ID_COLUMN=NA,
  NORMALIZATION="none",
  REF_CT_COLUMN=NA,
  REF_RNB_SET=NA,
  REF_RNB_CT_COLUMN=NA,
  PREPARE_TRUE_PROPORTIONS=FALSE,
  TRUE_A_TOKEN=NA,
  HOUSEMAN_A_TOKEN=NA,
  ESTIMATE_HOUSEMAN_PROP=FALSE,
  FILTER_BEADS=!is.null(rnb.set@covg.sites),
  FILTER_INTENSITY=inherits(rnb.set, "RnBeadRawSet"),
  FILTER_NA=TRUE,
  FILTER_CONTEXT=TRUE,
  FILTER_SNP=TRUE,
  FILTER_SOMATIC=TRUE,
  MATLAB_EXPORT=FALSE
)
rnb.set.filtered <- res$rnb.set.filtered

meth.data <- meth(rnb.set.filtered)
sds<-apply(meth.data, 1, sd)
sortedsdsrownumb <-order(sds, decreasing=TRUE)
sortedsds <- sds[sortedsdsrownumb]

pdf('finalplot.pdf')
plot(sortedsds[1:60000])
dev.off()

cg_subsets <- sortedsdsrownumb[1:20000]
row20k <- list(cg_subsets)
source("/scratch/divanshu/newset/start_analysis.R")
md.res<-start_medecom_analysis(
  rnb.set=rnb.set.filtered,
  WORK_DIR="/scratch/divanshu/newset/medecomrun",
  cg_groups=row20k,
  Ks=2:15,
  LAMBDA_GRID=c(0,10^(-5:-1)),
  SAMPLE_SUBSET=NULL,
  K_FIXED=NULL,
  WRITE_FILES=FALSE,
  startT=NULL,
  startA=NULL,
  CLUSTER_SUBMIT=FALSE,
  CLUSTER_RDIR=NA,
  CLUSTER_HOSTLIST="*",
  CLUSTER_MEMLIMIT="5G",
  #		MAX_JOBS=1000,
  #		WAIT_TIME="30m",
  #		PORTIONS=FALSE,
  #		JOB_FILE=NA,
  CLEANUP=FALSE
)

saveRDS(md.res,"medecomoutput.rds")
library("FactorViz")
startFactorViz(md.res,input.data=rnb.set.filtered)
