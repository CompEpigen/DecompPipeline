suppressPackageStartupMessages(library(RnBeads))
options(fftempdir="/scratch/divanshu/newset/temp",disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)
rnb.options(disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)

###############    extracting sample sheet  ################################################################

rnb.set<-rnb.execute.import(data.source="GSE42861", data.type="GEO")

#########      add barcode   
x <- pheno(rnb.set)
barcode <- x$supplementary_file.1
barcode <- gsub("_Red.idat.gz","",barcode)
barcode <- sub(".+?_","",barcode)
x <- cbind(x,barcode)

###########  rename file in datasets as per the barcode 
#####   for file in `ls -1 | grep idat`; do new_file_name=`echo $file | sed -e 's/^GSM[0-9]\+_//'`; echo $new_file_name; mv $file $new_file_name; done

############################   creating rnb.set #####################################################

rnb.set<-rnb.execute.import(data.source=list("/scratch/divanshu/newset/datasets", x), data.type="idat.dir")
x <- rnb.set@pheno
##############    select pheno data to be retained 
#############  x <- x[c(1,,,,,,,,,,,,,,,,,,,)]
#############  rnb.set@pheno <- x
saveRDS(rnb.set,"/scratch/divanshu/newset/rnb.set.rds")

######################    select healthy 100 sample from 684 samples 

x <- rnb.set@pheno
y <- x$subject
temp <- -1
for(i in 1:689){
  if(y[i] == y[1]){ temp <- cbind(temp,i)} }
temp <- temp[-1]
x <- c(1:689)
y <- setdiff(x,temp)
a <- length(y)
y <- y[100:a]
y <- c(temp,y)
rnb.set <- remove.samples(rnb.set, y)

#####################     filter snp probes    
##########          prepare data function by pavlo (https://github.com/lutsik/DecompPipeline/blob/master/R/data_preparation.R)

rnb.set <- readRDS("/scratch/divanshu/newset/healthyrnb.set.rds")

source("/scratch/divanshu/newset/filtering.R")
res<-prepare_data(
  RNB_SET=rnb.set, 
  WORK_DIR="/scratch/divanshu/newset/temp",
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
saveRDS(rnb.set.filtered,"filtered_rnb.set.rds")
########################    selecting top rows based on sd ##################################

meth.data <- meth(rnb.set.filtered)
sds<-apply(meth.data, 1, sd)
sortedsdsrownumb <-order(sds, decreasing=TRUE)

sortedsds <- sds[sortedsdsrownumb]

pdf('plot.pdf')
plot(sortedsds)
dev.off()

pdf('finalplot.pdf')
plot(sortedsds[1:60000])
dev.off()

############    send plot to remoter desktop ################################
######  scp ./plot.pdf c010-generic@c010-ngs:/home/c010-generic/Documents/divanshu/684sampledataset
######  scp ./finalplot.pdf c010-generic@c010-ngs:/home/c010-generic/Documents/divanshu/684sampledataset


###################       running medecom for selected top 25K

cg_subsets <- sortedsdsrownumb[1:25000]
row25k <- list(cg_subsets)
#################### RUN medecom 
#########################     (https://github.com/lutsik/DecompPipeline/blob/master/R/start_analysis.R)
## scp c010-generic@c010-ngs:/home/c010-generic/Documents/divanshu/684sampledataset/runmedecomanalysis.R /scratch/divanshu/newset/start_analysis.R

source("/scratch/divanshu/newset/start_analysis.R")
md.res<-start_medecom_analysis(
  rnb.set=rnb.set.filtered,
  WORK_DIR="/scratch/divanshu/newset/medecomrun",
  cg_groups=row25k,
  Ks=11:15,
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

###########################  Run factorViz
#########                        (https://github.com/lutsik/FactorViz)
##### scp -r c010-generic@c010-ngs:/home/c010-generic/Documents/divanshu/684sampledataset/FactorViz-master/R /scratch/divanshu/newset/FactorViz

library(shiny)
library("shinyURL")
library("FactorViz")
startFactorViz(md.res,input.data=rnb.set.filtered)

########333  port forwarding                 ssh guptad@c010-srv1 -L 4972:localhost:4972 -fN

###############################################################################################################################3
#######################     Getting TREF for comparing withthe other datasets 

md.res <- readRDS("./medecomoutput.RDS")

######    check the one with minimum error
###################      x$cve


x <- md.res@output[[1]]
a <- x[[1]]
b <- a$T
c <- b[9,4]
Tmedecom <- c[[1]]

library(pheatmap)
png(' heatmap of medecom resultant T .png')
pheatmap(as.matrix(Tmedecom))
dev.off()



