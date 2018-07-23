suppressPackageStartupMessages(library(RnBeads))
source("/home/guptad/divanshu/scripts/00_raw_signals.R")


options(fftempdir="/home/guptad/divanshu/exp4/temp",disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)
rnb.options(disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)

pd <- data.frame(SampleID = paste0("GSM",c(861635:861694)))

intensity_list<-read.raw.signals.1("datasets/GSE59250_signal_intensities.txt", cnames = c(" Unmethylated Signal", " Methylated Signal", " Detection Pval"),
                                   probe.names = NULL, verbose = TRUE)



suppressPackageStartupMessages(library(RnBeads))
options(fftempdir="/home/guptad/divanshu/exp4/temp",disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)
rnb.options(disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)

rnb.set.raw<-RnBeadRawSet(pheno=a, M=intensity_list$M, U=intensity_list$U, p.values=intensity_list$pvalues)

saveRDS(rnb.set.raw,"/scratch/divanshu/finalsets/exp4rnb.set.raw.rds")


save.image("rnb.setimage.RData")