###############################################################################
## These scripts take RnBeads data from disk and perform ICA in order to:
##  - get optimal number of components, getComponentNumber()
##  - detect features linked to several factors, getFeatures()
##  - remove effect of a single confounding factor, removeFactor()
## P.Nazarov in collaboration with M.Scherer, T.Kaoma, P.Lutsik
## ver. 2019-03-24
###############################################################################

##=============================================================================
## Scripts of LibICA.r that could be used / added to the package
##    runICA() - runs consensus ICA with parallelization;
##    getGenesICA() - gets features contributing to each component
##    getGO() - assigns contributing features to Gene Ontologies (functional annotation)
##    saveGO() - saves features and enriched GO terms
##    saveReport() - saves report of ICA including top features, GO-terms and clinical factors.
##=============================================================================

## New scripts

##=============================================================================
## getComponentNumber(file,fact,nmin,nmax)
## Searches for the number of components required to identify factor effect
## with maximized significance
##    file - file path to RnBeads dataset
##    fact - vector of factor names to be investigated
##    nmin - minimal number of components to try
##    nmax - maximal number of components to try
##    ntry - number of ICA runs (1 is fastest, but leads to noisy log-p-value profiles)
##    thr.sd - minimal st.dev. of the beta value of a methylation site to be included into the analysis
##             set it to 0 to avoid filtering
## Returns: optimal number of components (scalar)
##-----------------------------------------------------------------------------
getComponentNumber = function(rnb.set, fact, nmin=3, nmax=20, ntry=1, ncores = 1, thr.sd=0.05, assume.numbers = TRUE){
  # get methylation data & annotations as data frames
  meth.data = meth(rnb.set)
  pheno.data = pheno(rnb.set)
  anno.data = annotation(rnb.set)
  ikeep = apply(meth.data,1,sd)>thr.sd
  logger.info(sprintf("Only sites with SD > %g were kept: %d of %d\n",thr.sd,sum(ikeep),length(ikeep)))
  
  ## check whether factors are, in fact, numbers
  if (assume.numbers){
    icol=1
    for (icol in 1:ncol(pheno.data)){
      if (class(pheno.data[[icol]])=="factor"){
        if (nlevels(pheno.data[[icol]]) > nrow(pheno.data)*0.25) {## if nlevels > 25% of nsamples
          v = pheno.data[[icol]]
          v =as.numeric( as.character(v))
          if (sum(!is.na(v)) > nrow(pheno.data)*0.25) {
            pheno.data[[icol]] = v
            logger.info(paste("Assuming numeric data for pheno column `",colnames(pheno.data)[icol]))
          }
        }
      }
    }
  }
  
  Var = df2factor(pheno.data,maxlev=20)
  if(sum(fact %in% names(Var))!= length(fact)) {
    logger.warning(paste("Error in getComponentNumber: poor factor or wrong factor name `",fact[! fact %in% names(Var)],"factor might have more than 20 levels, but cannot be treated as numeric."))
    return(NA)
  }
  
  if (nmin>=nmax) {
    logger.error("Error in getComponentNumber: nmin >= nmax\n")
    return(NA)
  }
  
  ncomp = nmin:nmax
  logPV = data.frame(matrix(nrow=length(ncomp),ncol=length(fact)))
  rownames(logPV) = paste0("nc.",ncomp)
  colnames(logPV) = fact
  nc = ncomp[1]
  for (nc in ncomp){
    logger.info(paste("getMinCompNumber: working with",nc,"components"))
    IC = runICA(meth.data,ncomp=nc, ntry = ntry, ncores = ncores)
    logPV[paste0("nc.",nc),] = 1
    for (fct in fact){
      for (ic in 1:nc){
        mod = aov(IC$M[ic,]~Var[[fct]])
        logPV[paste0("nc.",nc),fct] = min(logPV[paste0("nc.",nc),fct],summary(mod)[[1]][1,"Pr(>F)"])
      }
    }
    logger.info(paste(logPV[paste0("nc.",nc),],collapse="; "))
    logPV[paste0("nc.",nc),] = -log10(logPV[paste0("nc.",nc),])
  }
  
  
  return(ncomp[which.max(apply(logPV,1,min))])
  ##ToDo: later we could perform shape analysis - fitting with sigmoid
  ## e.g.
  #sigreg = function(x,p){return(p[1]+p[2]/(1+exp(-p[4]^2*(x-p[3]))))}
  #x = ncomp
  #y = logPV
  #p0 = list(a=min(logPV),b=max(logPV)-min(logPV),c=mean(ncomp),d=1)
  #res = try(nls( y ~ a+b/(1+exp(-d^2*(x-c))),start=p0,trace=F,control=nls.control(warnOnly=T),algorithm="port",weights=y),silent = T)
  #xest = seq(min(x),max(x),(max(x)-min(x))*0.001)
  #yest = sigreg (xest,summary(res)$coefficients[,1])
  #plot(x,y,pch=19)
  #lines(xest,yest,col=2)
}

##=============================================================================
## getFeatures(file,fact,ncomp,alpha)
## Generates feature indexes that are important for the top component linked to the factor
##    file - file path to RnBeads dataset
##    fact - factor name(s) to be investigated
##    ncomp - number of components in ICA deconvolution
##    ntry - number of ICA runs
##    alpha.fact - significance level linking component to a factor (< 1e-5 is recommended)
##    alpha.feat - significance level detecting features contributing to component
## Returns: logical vector showing whether features are linked to any of `fact` factors
##-----------------------------------------------------------------------------
getFeatures = function(rnb.set, fact, ncomp=3, ntry = 1, alpha.fact =1e-20, alpha.feat=0.01, assume.numbers = TRUE){
  ## get methylation data & annotations as data frames
  meth.data = meth(rnb.set,row.names=T)
  pheno.data = pheno(rnb.set)
  anno.data = annotation(rnb.set)
  
  ## check whether factors are, in fact, numbers
  if (assume.numbers){
    icol=1
    for (icol in 1:ncol(pheno.data)){
      if (class(pheno.data[[icol]])=="factor"){
        if (nlevels(pheno.data[[icol]]) > nrow(pheno.data)*0.25) {## if nlevels > 25% of nsamples
          v = pheno.data[[icol]]
          v =as.numeric( as.character(v))
          if (sum(!is.na(v)) > nrow(pheno.data)*0.25) {
            pheno.data[[icol]] = v
            logger.info(paste("Assuming numeric data for pheno column `",colnames(pheno.data)[icol]))
          }
        }
      }
    }
  }
  
  
  ## variables = factors
  Var = df2factor(pheno.data,maxlev=20)
  if(!all(fact %in% names(Var))) {
    logger.error(paste("Error in getFeatures: poor factor or wrong factor name `",fact))
    return(NA)
  }
  ## ICA
  IC = runICA(meth.data,ncomp=ncomp, ntry = ntry)
  ## assign components to factors
  PV = matrix(nrow=length(fact),ncol=ncomp)
  rownames(PV)=fact
  for (fa in fact){ 
    for (ic in 1:ncomp){
      mod = aov(IC$M[ic,]~Var[[fa]])
      PV[fa,ic] = summary(mod)[[1]][1,"Pr(>F)"]
    }
  }
  ## feature names are required to select and combine them
  #features = paste("f",1:nrow(meth.data),sep=".")
  features = rownames(meth.data)
  rownames(IC$X) = features
  rownames(IC$S) = features
  Genes=getGenesICA(IC,alpha.feat)
  sign.features = NULL
  for (fa in fact){ 
    for (ic in 1:ncomp){
      if (PV[fa,ic]<alpha.fact) {
        logger.info(paste("for",fa,"we take component",ic))
        sign.features = unique(c(sign.features,c(rownames(Genes[[ic]]$pos),rownames(Genes[[ic]]$neg))))
      }
    }
  }
  return(features %in% sign.features)
}

##=============================================================================
## removeFactors(file,facts,ncomp,alpha)
## Generates features that are important for the top component linked to MULTIPLE factors
##    file - file path to RnBeads dataset
##    fact - the names of several factors to be removed
##    ncomp - number of components in ICA deconvolution
##    ntry - number of ICA runs
##    alpha.fact - significance level linking component to a factor
##    qthr - quantile at which boundaries are put on the corrected signal (to make it [0,1])
## Returns: matrix of corrected beta values
##-----------------------------------------------------------------------------
removeFactor = function(rnb.set, fact, ncomp=3, ntry = 1, alpha.fact =1e-20, qthr=0.01, assume.numbers = TRUE, save.report=NULL){
  ## get methylation data & annotations as data frames
  meth.data = meth(rnb.set,row.names=T)
  pheno.data = pheno(rnb.set)
  anno.data = annotation(rnb.set)
  
  ## check whether factors are, in fact, numbers
  if (assume.numbers){
    icol=1
    for (icol in 1:ncol(pheno.data)){
      if (class(pheno.data[[icol]])=="factor"){
        if (nlevels(pheno.data[[icol]]) > nrow(pheno.data)*0.25) {## if nlevels > 25% of nsamples
          v = pheno.data[[icol]]
          v =as.numeric( as.character(v))
          if (sum(!is.na(v)) > nrow(pheno.data)*0.25) {
            pheno.data[[icol]] = v
            logger.info(paste("Assuming numeric data for pheno column `",colnames(pheno.data)[icol]))
          }
        }
      }
    }
  }
  
  
  ## variables = factors
  Var = df2factor(pheno.data,maxlev=20)
  if(!fact %in% names(Var)) {
    logger.info(paste("Error in getFeatures: poor factor or wrong factor name `",fact))
    return(NA)
  }
  ## ICA
  IC = runICA(meth.data,ncomp=ncomp, ntry = ntry)
  if(!is.null(save.report)){
    print(head(IC$S))
    saveReport(IC=IC,Var=Var[,fact],file=save.report) 
    #dev.off()
  }
  ## assign components to removed factor
  pv = double(ncomp)+1
  for (fa in fact){ ## loop for factors
    for (ic in 1:ncomp){  ## loop for components
      mod = aov(IC$M[ic,]~Var[[fa]])
      pv1 = summary(mod)[[1]][1,"Pr(>F)"]
      if (pv1 < alpha.fact){
        logger.info(paste("Component",ic,"is linked to",fa,"factor, p-value=",pv1))
      }
      pv[ic] = min(pv[ic],pv1)
    }
  }
  
  ## check whether any component is linked to confounding factor(s)
  if (sum(pv < alpha.fact)==0){
    logger.info(paste("Cannot find a component linked to",paste(fact,collapse=","),"factor(s). Min p-value=",min(pv),",while alpha=",alpha.fact))
    return(meth.data)
  }
  
  ##correction
  M1 = IC$M
  M1[pv < alpha.fact,] = 0
  meth.data1 = IC$S %*% M1
  meth.data1 = (meth.data1 - quantile(meth.data1,qthr)) / (quantile(meth.data1,1-qthr) - quantile(meth.data1,qthr))
  meth.data1[meth.data1<0] = 0
  meth.data1[meth.data1>1] = 1
  return(meth.data1)
  ## ToDo: compare with original distribution of `meth.data`!
}

#' run.rnb.ICA
#' 
#' Performs ICA to remove the effect of a confounding factor from the methylation matrix
#' 
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation information
#' @param conf.factor A vector of column names in the sample annotation sheet of \code{rnb.set} representing confounding factors to be 
#'                removed.
#' @param ica.setting Optional argument. If specified a named vector of arguments to be used. Can be one of the following. 
#' @param nmin Minimum number of components to be used
#' @param nmax Maximum number of components to be used
#' @param ntry Further argument for runICA
#' @param thr.sd Threshold for the standard deviation across samples. Only sites with a standard deviation larger than this
#'                threshold are kept. 0 deactivates filtering.
#' @param alpha.fact Significance level for the factor
#' @param type Analysis type to be performed. Can be either \code{"remove"} to remove the effect of the confounding
#'              factor or \code{"keep"} to export the sites that are linked to the confounding factor. Note that for
#'              \code{"keep"}, the factor effect is not removed, but only reported.
#' @param out.folder Folder to store ICA's output and diagnostic plots
#' @return The modified \code{rnb.set} object with updated methylation values. The effect of the confouding factor is removed
#'      using independent component analysis (ICA).
#' @author Michael Scherer. ICA code was generated by Peter Nazarov and Tony Kamoa
#' @export
run.rnb.ICA <- function(rnb.set,conf.factor,ica.setting=NULL,nmin=10,nmax=30,ntry=1,thr.sd=0,
                        alpha.fact=1e-10,type="remove",ncores=1,out.folder=NULL,save.report=F,alpha.feat=0.01){
  if(!inherits(rnb.set,"RnBSet")){
    logger.error("Invalid value for rnb.set, needs to be RnBSet object")
  }
  if(!is.null(ica.setting)){
    if(!is.na(ica.setting["nmin"])){
      nmin <- as.numeric(ica.setting["nmin"])
    }
    if(!is.na(ica.setting["nmax"])){
      nmax <- as.numeric(ica.setting["nmax"])
    }
    if(!is.na(ica.setting["ntry"])){
      ntry <- as.numeric(ica.setting["ntry"])
    }
    if(!is.na(ica.setting["thr.sd"])){
      thr.sd <- as.numeric(ica.setting["thr.sd"])
    }
    if(!is.na(ica.setting["alpha.fact"])){
      alpha.fact <- as.numeric(ica.setting["alpha.fact"])
    }
    if(!is.na(ica.setting["save.report"])){
      save.report <- ica.setting["save.report"]
    }
    if(!is.na(ica.setting["ncores"])){
      ncores <- as.numeric(ica.setting["ncores"])
    }
    if(!is.na(ica.setting["type"])){
      type <- ica.setting["type"]
    }
    if(!is.na(ica.setting["alpha.feat"])){
      alpha.feat <- as.numeric(ica.setting["alpha.feat"])
    }
  }
  if(nmax>length(samples(rnb.set))){
    logger.warning(paste("Number of components cannot be larger than number of samples, setting to",length(samples(rnb.set))))
    nmax <- length(samples(rnb.set))
  }
  if(nmin>=nmax){
    logger.warning(paste("nmin cannot be smaller than nmax, set to",nmax-1))
    nmin <- nmax-1
  }
  if(!isImputed(rnb.set)){
    rnb.set <- rnb.execute.imputation(rnb.set)
  }
  source("http://sablab.net/scripts/LibICA.r")
  logger.start("Determining number of components")
  ncomp <- getComponentNumber(rnb.set,conf.factor,nmin=nmin,nmax=nmax,ntry=ntry,thr.sd=thr.sd,ncores = ncores)
  logger.completed()
  if(type=="remove"){
    if(!is.null(out.folder)){
      old.meth <- meth(rnb.set)  
    }
    logger.start("Removing factor effect")
    if(save.report){
      if(!requireNamespace("sm")) install.packages("sm",repos="https://cloud.r-project.org/")
      if(!requireNamespace("vioplot")) install.packages("vioplot",repos="https://cloud.r-project.org/")
      if(!requireNamespace("pheatmap")) install.packages("pheatmap",repos="https://cloud.r-project.org/")
      new.meth <- removeFactor(rnb.set,fact = conf.factor, ncomp = ncomp,ntry = ntry, alpha.fact = alpha.fact,save.report=file.path(out.folder,"violin_report_ICA.pdf"))
    }else{
      new.meth <- removeFactor(rnb.set,fact = conf.factor, ncomp = ncomp,ntry = ntry, alpha.fact = alpha.fact)
    }
    logger.completed()
    rnb.set <- updateMethylationSites(rnb.set,new.meth)
    if(!is.null(out.folder)){
      old.meth <- melt(old.meth)
      new.meth <- melt(new.meth)
      to.plot <- data.frame("Before"=old.meth$value,"After"=new.meth$value)
      to.plot <- melt(to.plot)
      colnames(to.plot) <- c("Transformation","Beta")
      plot <- ggplot(to.plot,aes(x=Beta,y=..density..,color=Transformation))+geom_density()+theme_bw()
      ggsave(file.path(out.folder,"ICA_beta_comparison.pdf"),plot,device = "pdf")
    }
  }else{
    logger.start("Determine CpGs linked to factor")
    sel.features <- which(getFeatures(rnb.set,fact = conf.factor,ncomp = ncomp, ntry = ntry, alpha.fact = alpha.fact, alpha.feat = alpha.feat))
    if(is.null(out.folder)){
      out.folder <- file.path(getwd(),"Decomp_output")
      if(!dir.exists(out.folder)){
        dir.create(out.folder)
      }
    }
    write.csv(sel.features,file.path(out.folder,paste("sites_linked_to",paste(conf.factor,collapse = "_"),"data.csv",sep="_")))
    logger.completed()
  }
  return(rnb.set)
}
