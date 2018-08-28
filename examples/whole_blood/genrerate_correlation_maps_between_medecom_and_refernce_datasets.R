suppressPackageStartupMessages(library(RnBeads))
library(MeDeCom)
library(pheatmap)

##########  working directory /sctrach/divanshu/Correlation
filteredrnb.set <- readRDS("./filteredrnb.set.rds")
md.res <- readRDS("./medecomoutk_5_15.rds")
a <- md.res@outputs[[1]]
Tmed <- a$T
meth.data <- meth(filteredrnb.set,row.names = TRUE)
sds<-apply(meth.data, 1, sd)
sortedsdsrownumb <-order(sds, decreasing=TRUE)
selectedmeth.data <- meth.data[sortedsdsrownumb[1:25000],]
rnames <- rownames(selectedmeth.data)


generatedrnb.set <- readRDS("./completehealtyset.rds")
generatedmeth.data <- meth(generatedrnb.set,row.names = TRUE)
selectedfromgenerated <- generatedmeth.data[rnames,]

ph <- pheno(generatedrnb.set)
a <- unique(ph$celltype)
Trefofhealthy <-  matrix(data = 0, nrow = 25000, ncol = length(a))
colnames(Trefofhealthy) <- a
for(i in 1:length(a)){
  b <- which(ph$celltype == a[i])
  if(length(b) >= 2)
    Trefofhealthy[,i] <- rowMeans(selectedfromgenerated[,b],na.rm = TRUE) 
  else
    Trefofhealthy[,i] <- selectedfromgenerated[,b]
}

MeDeCom:::components.heatmap(Tmed[[numb]],Trefofhealthy,centered = TRUE)



numb <- 44
cormatbw_Tmedecom_and_Trefofhealthy <- cor(Trefofhealthy,Tmed[[numb]], use = "complete.obs")
png('try till success part 2 .png')
pheatmap(as.matrix(cormatbw_Tmedecom_and_Trefofhealthy))
dev.off()

