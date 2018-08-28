suppressPackageStartupMessages(library(RnBeads))
rnb.options(qc = FALSE,identifiers.column = "sample_id",exploratory.region.profiles=  vector('character'))

rnb.set <- readRDS("healthyrnbset.rds")

file <- pheno(rnb.set)
x <- file$title
y <- file$celltype
z <- file$geo_accession
q <- file$geo_id
q <- paste0(x,"_",y,"_",z,"_",q)
file$sample_id <- q
rnb.set@pheno <- file

system("rm -r /scratch/divanshu/finalsets/analysisreports/temp")
x <- rnb.run.analysis(dir.reports = "/scratch/divanshu/finalsets/analysisreports/temp", data.source = rnb.set,data.type = "rnb.set")


############################ ploting heatmap #################################################


completehealtyset <- readRDS("completehealtyset.rds")
library(genefilter)

######   fac = factor(floor(runif(ncol(Dmatrix))*10))
ph <- pheno(completehealtyset)
Dmatrix <- meth(completehealtyset)
x <- rowFtests(Dmatrix,ph$celltypegeneral)
temp <- x
sortedrownumb <- order(temp$statistic , decreasing=TRUE) 
###data <- temp$statistic[sortedrownumb]
library(pheatmap)
png(' completehealtyset with cell general type.png')
pheatmap(as.matrix(Dmatrix[sortedrownumb[1:20000],]))
dev.off()


###################              making TREF from the dataset to compare with the other data

temprnb <- readRDS("/scratch/divanshu/newset/final_medecomoutput.rds")
cg_sub <- temprnb@parameters$cg_subsets[[1]]


ph <- pheno(completehealtyset)
Dmatrix <- meth(completehealtyset)
Dmatrix <- Dmatrix[cg_sub,]

a <- unique(ph$celltypegeneral)

Tref <-  matrix(data = 0, nrow = 25000, ncol = length(a))

for(i in 1:length(a)){
  b <- which(ph$celltypegeneral == a[i])
  Tref[,i] <- rowMeans(Dmatrix[,b]) 
}


library(pheatmap)
png(' heatmap of Tref by taking averages .png')
pheatmap(as.matrix(Tref))
dev.off()

###############################################################################




