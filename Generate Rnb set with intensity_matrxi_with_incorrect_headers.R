suppressPackageStartupMessages(library(RnBeads))
source("/home/guptad/divanshu/scripts/00_raw_signals.R")

options(fftempdir="/home/guptad/divanshu/exp3/temp",disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)
rnb.options(disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)

fname <- "/home/guptad/divanshu/exp3/datasets/GSE59250_signal_intensities.txt"
content.lines <- scan(fname, "", sep = "\n", quiet = TRUE)
content.lines <- content.lines[!grepl("^\\s*$", content.lines)]
content.lines <- strsplit(sub("\\s+$", "", content.lines[-1]), "\t", fixed = TRUE)


temp <- strsplit(sub("\\s+$", "", content.lines[-1]), "\t", fixed = TRUE)

N = 482000
i.ids <- seq(2, 1303, by = 3L)
init.matrix <- function(v) {
  matrix(v, nrow = 482000, ncol = length(i.ids))
}
mu <- init.matrix(as.integer(NA))
mm <- init.matrix(as.integer(NA))
mp <- init.matrix(as.double(NA))

for (i in 1:482000) {
  txt <- content.lines[[i]]
  mu[i, ] <- as.integer(txt[i.ids])
  mm[i, ] <- as.integer(txt[i.ids + 1L])
  mp[i, ] <- as.double(txt[i.ids + 2L])
}
intensity_list <- list("U" = mu, "M" = mm, "pvalues" = mp)

Probes = "temp"
for (i in 1:482000) {
  Probes <- rbind(Probes,content.lines[[i]][1])
}
Probes <- Probes[-1]

pd <- read.table("/home/guptad/divanshu/exp3/samplesheet.csv")
y <- gsub(".*\\.","",x)
pd <- cbind(pd,y)

suppressPackageStartupMessages(library(RnBeads))
options(fftempdir="/home/guptad/divanshu/exp3/temp",disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)
rnb.options(disk.dump.big.matrices = FALSE,disk.dump.bigff = FALSE)

rnb.set.raw<-RnBeadRawSet(pheno=a,probes = Probes,M=intensity_list$M, U=intensity_list$U, p.values=intensity_list$pvalues)
saveRDS(rnb.set.raw,"/scratch/divanshu/finalsets/exp3rnb.set.raw.rds")


save.image("rnb.setimage.RData")



