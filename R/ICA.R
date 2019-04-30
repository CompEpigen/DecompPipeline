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
##    fact - factor name to be investigated
##    nmin - minimal number of components to try
##    nmax - maximal number of components to try
##    ntry - number of ICA runs (1 is fastest, but leads to noisy log-p-value profiles)
##    thr.sd - minimal st.dev. of the beta value of a methylation site to be included into the analysis
##             set it to 0 to avoid filtering
## Returns: optimal number of components (scalar)
##-----------------------------------------------------------------------------
getComponentNumber = function(rnb.set, fact, nmin=3, nmax=20, ntry=1, thr.sd=0.05){
  # get methylation data & annotations as data frames
  meth.data = meth(rnb.set)
  pheno.data = pheno(rnb.set)
  anno.data = annotation(rnb.set)
  ikeep = apply(meth.data,1,sd)>thr.sd
  logger.info(sprintf("Only sites with SD > %g were kept: %d of %d",thr.sd,sum(ikeep),length(ikeep)))
  
  Var = df2factor(pheno.data,maxlev=20)
  if(!fact %in% names(Var)) {
    logger.error(paste("Error in getComponentNumber: poor factor or wrong factor name `",fact))
    return(NA)
  }
  
  if (nmin>=nmax) {
    logger.error(paste("Error in getComponentNumber: nmin >= nmax"))
    return(NA)
  }
  
  ncomp = nmin:nmax
  logPV = double(length(ncomp))+NA
  for (nc in ncomp){
    logger.start(paste("getMinCompNumber: working with",nc,"components"))
    IC = runICA(meth.data,ncomp=nc, ntry = ntry)
    pv=1
    for (ic in 1:nc){
      mod = aov(IC$M[ic,]~Var[[fact]])
      pv = min(pv,summary(mod)[[1]][1,"Pr(>F)"])
    }
    logPV[which(is.na(logPV))[1]] = -log10(pv)
    write(print(pv),file=logger.getfiles(),append = T)
    logger.completed()
  }
  return(ncomp[which.max(logPV)])
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
getFeatures = function(rnb.set, fact, ncomp=3, ntry = 1, alpha.fact =1e-20, alpha.feat=0.01){
  rnb.set <- rnb.execute.imputation(rnb.set)
  ## get methylation data & annotations as data frames
  meth.data = meth(rnb.set)
  pheno.data = pheno(rnb.set)
  anno.data = annotation(rnb.set)
  ## variables = factors
  Var = df2factor(pheno.data,maxlev=20)
  if(!all(fact %in% names(Var))) {cat("Error in getFeatures: poor factor or wrong factor name `",fact,"`\n");return(NA)}
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
  features = paste("f",1:nrow(meth.data),sep=".")
  rownames(IC$X) = features
  rownames(IC$S) = features
  Genes=getGenesICA(IC,alpha.feat)
  sign.features = NULL
  for (fa in fact){ 
    for (ic in 1:ncomp){
      if (PV[fa,ic]<alpha.fact) {
        cat("for",fa,"we take component",ic,"\n")
        sign.features = unique(c(sign.features,c(rownames(Genes[[ic]]$pos),rownames(Genes[[ic]]$neg))))
      }
    }
  }
  return(features %in% sign.features)
}

##=============================================================================
## remFactor(file,fact,ncomp,alpha)
## Generates features that are important for the top component linked to the factor
##    file - file path to RnBeads dataset
##    fact - the name of the factor to be removed
##    ncomp - number of components in ICA deconvolution
##    ntry - number of ICA runs
##    alpha.fact - significance level linking component to a factor
##    qthr - quantile at which boundaries are put on the corrected signal (to make it [0,1])
## Returns: matrix of corrected beta values
##-----------------------------------------------------------------------------
removeFactor = function(rnb.set, fact, ncomp=3, ntry = 1, alpha.fact =1e-20, qthr=0.01){
  ## get methylation data & annotations as data frames
  meth.data = meth(rnb.set)
  pheno.data = pheno(rnb.set)
  anno.data = annotation(rnb.set)
  ## variables = factors
  Var = df2factor(pheno.data,maxlev=20)
  if(!fact %in% names(Var)){
    logger.error(paste("Error in getFeatures: poor factor or wrong factor name `",fact))
    return(NA)
  }
  ## ICA
  IC = runICA(meth.data,ncomp=ncomp, ntry = ntry)
  ## assign components to removed factor
  pv = double(ncomp)
  for (ic in 1:ncomp){
    mod = aov(IC$M[ic,]~Var[[fact]])
    pv[ic] = summary(mod)[[1]][1,"Pr(>F)"]
    if (pv[ic] < alpha.fact){
      logger.info(paste("Component",ic,"is linked to",fact,"factor, p-value=",pv[ic]))
    } 
  }
  
  ## check whether any component is linked to confounding factor
  if (sum(pv < alpha.fact)==0){
    logger.info(paste("Cannot find a component linked to",fact,"factor. Min p-value=",min(pv),",while alpha=",alpha.fact))
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
#' @param conf.factor A column name in the sample annotation sheet of \code{rnb.set} representing the confounding factor to be 
#'                removed.
#' @param ica.setting Optional argument. If specified a named vector of arguments to be used. Can be one of the following. 
#' @param nmin Minimum number of components to be used
#' @param nmax Maximum number of components to be used
#' @param ntry Further argument for runICA
#' @param thr.sd Threshold for the standard deviation across samples. Only sites with a standard deviation larger than this
#'                threshold are kept. 0 deactivates filtering.
#' @param alpha.fact Significance level for the factor
#' @return The modified \code{rnb.set} object with updated methylation values. The effect of the confouding factor is removed
#'      using independent component analysis (ICA).
#' @author Michael Scherer. ICA code was generated by Peter Nazarov and Tony Kamoa
#' @export
run.rnb.ICA <- function(rnb.set,conf.factor,ica.setting=NULL,nmin=10,nmax=30,ntry=1,thr.sd=0,alpha.fact=1e-10,out.folder=NULL){
  if(!inherits(rnb.set,"RnBSet")){
    logger.error("Invalid value for rnb.set, needs to be RnBSet object")
  }
  if(!is.null(ica.setting)){
    if(!is.na(ica.setting["nmin"])){
      nmin <- ica.setting["nmin"]
    }
    if(!is.na(ica.setting["nmax"])){
      nmax <- ica.setting["nmax"]
    }
    if(!is.na(ica.setting["ntry"])){
      nrty <- ica.setting["ntry"]
    }
    if(!is.na(ica.setting["thr.sd"])){
      thr.sd <- ica.setting["thr.sd"]
    }
    if(!is.na(ica.setting["alpha.fact"])){
      alpha.fact <- ica.setting["alpha.fact"]
    }
  }
  if(nmin>nmax){
    logger.warning("nmin cannot be smaller than nmax, set to nmax")
    nmin <- nmax
  }
  if(nmax>length(samples(rnb.set))){
    logger.warning(paste("Number of components cannot be larger than number of samples, setting to",length(samples(rnb.set))))
    nmax <- length(samples(rnb.set))
  }
  if(!isImputed(rnb.set)){
    rnb.set <- rnb.execute.imputation(rnb.set)
  }
  source("http://sablab.net/scripts/LibICA.r")
  logger.start("Determining number of components")
  ncomp <- getComponentNumber(rnb.set,conf.factor,nmin=nmin,nmax=nmax,ntry=ntry,thr.sd=thr.sd)
  logger.completed()
  if(!is.null(out.folder)){
    old.meth <- meth(rnb.set)  
  }
  logger.start("Removing factor effect")
  new.meth <- removeFactor(rnb.set,fact = conf.factor, ncomp = ncomp,ntry = ntry, alpha.fact = alpha.fact)
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
  return(rnb.set)
}
