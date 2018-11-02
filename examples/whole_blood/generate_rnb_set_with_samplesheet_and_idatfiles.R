# cd /home/guptad/divanshu
# mkdir exp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cd exp%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# mkdir temp
# mkdir datasets
# cd datasets
# wget '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
# mv '%%%%%%%%%%%%%%%%%%%%%%%%%5' &&&&&&&&&&.tar
# tar -xvf %%%%%%%.tar
# gunzip -d *.gz
# for file in `ls -1 | grep idat`; do new_file_name=`echo $file | sed -e 's/^GSM[0-9]\+_//'`; echo $new_file_name; mv $file $new_file_name; done
##  scp c010-generic@c010-ngs:/home/c010-generic/Documents/divanshu/pipeline/sheet.csv /
# scp -r /home/guptad/divanshu/rnbsets c010-generic@c010-ngs:/home/c010-generic/Documents/divanshu/pipeline

suppressPackageStartupMessages(library(RnBeads))
options(fftempdir="/scratch/divanshu/newset/temp",disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)
rnb.options(disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)

rnb.set<-rnb.execute.import(data.source="GSE42861", data.type="GEO")

x <- pheno(rnb.set)
barcode <- x$supplementary_file.1
barcode <- gsub("_Red.idat.gz","",barcode)
barcode <- sub(".+?_","",barcode)
x <- cbind(x,barcode)

rnb.set<-rnb.execute.import(data.source=list("/scratch/divanshu/newset/datasets", x), data.type="idat.dir")
saveRDS(rnb.set,"/scratch/divanshu/newset/rnb.set.rds")

save.image("rnb.setimage.RData")

source("/home/guptad/divanshu/pipeline/code.R")
res<-prepare_data(
  RNB_SET=rnb.set, 
  WORK_DIR=file.path("/home/guptad/divanshu/exp2"),
  DATASET="adgbbSorted",
  DATA_SUBSET="frontal"
)

source("/home/guptad/divanshu/pipeline/cpgsubset.R")

cg_subsets<-prepare_CG_subsets(
  res$rnb.set.filtered,
  MARKER_SELECTION=c("var5k", "var10k")
)
save.image("final.RData")
