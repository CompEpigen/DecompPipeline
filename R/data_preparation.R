#
# Data preparation script for the parameter tuning run
#

#################################################################################################################################
#' GLOBALS
#################################################################################################################################

#################################################################################################################################
#' FUNCTIONS
#################################################################################################################################

#' prepare_data
#' 
#' This functions prepares Illumina BeadChip data for a MeDeCom run.
#' 
#' @param RNB_SET An object of type \code{\link[RnBeads]{RnBSet-class}} for which analysis is to be performed.
#' @param WORK_DIR A path to a existing directory, in which the results are to be stored
#' @param analysis.name A string representing the dataset for which analysis is to be performed. Only used to create a folder with a 
#'                 descriptive name of the analysis.
#' @param SAMPLE_SELECTION_COL A column name in the phenotypic table of \code{RNB_SET} used to selected a subset of samples for
#'                 analysis that contain the string given in \code{SAMPLE_SELECTION_GREP}.
#' @param SAMPLE_SELECTION_GREP A string used for selecting samples in the column \code{SAMPLE_SELECTION_COL}.
#' @param PHENO_COLUMNS Vector of column names in the phenotypic table of \code{RNB_SET} that is kept and exported for further 
#'                 exploration.
#' @param ID_COLUMN Sample-specific ID column name in \code{RNB_SET}
#' @param NORMALIZATION Normalization method to be performed before employing MeDeCom. Can be one of \code{"none","dasen","illumina","noob"}.
#' @param REF_CT_COLUMN Column name in \code{RNB_SET} used to extract methylation information on the reference cell types.
#' @param REF_RNB_SET An object of type \code{\link[RnBeads]{RnBSet-class}} containing methylation information on reference cell types.
#' @param REF_RNB_CT_COLUMN Column name in \code{REF_RNB_SET} used to extract methylation information on the reference cell types.
#' @param PREPARE_TRUE_PROPORTIONS Flag indicating if true proportions are either available in \code{RNB_SET} or to be estimated 
#'                          with Houseman's reference-based deconvolution approach.
#' @param TRUE_A_TOKEN String present in the column names of \code{RNB_SET} used for selecting the true proportions of the corresponding
#'                      cell types.
#' @param HOUSEMAN_A_TOKEN Similar to \code{TRUE_A_TOKEN}, but not containing the true proportions, rather the estimated proportions
#'                      by Houseman's method.
#' @param ESTIMATE_HOUSEMAN_PROP If neither \code{TRUE_A_TOKEN} nor \code{HOUSEMAN_A_TOKEN} are given, the proportions of the reference
#'                      cell type are estimated with Houseman's approach.
#' @param FILTER_BEADS Flag indicating, if site-filtering based on the number of beads available is to be conducted.
#' @param MIN_N_BEADS Minimum number of beads required in each sample for the site to be considered for adding to MeDeCom.
#' @param FILTER_INTENSITY  Flag indicating if sites should be removed according to the signal intensities (the lowest and highest quantiles
#'                      given by \code{MIN_INT_QUANT} and \code{MAX_INT_QUANT}). Note that all sites are removed that have a value outside of
#'                      the provided quantile range in either of the channels and in any of the samples. 
#' @param MIN_INT_QUANT Lower quantile of intensities which is to be removed.
#' @param MAX_INT_QUANT Upper quantile of intensities which is to be removed.
#' @param FILTER_NA Flag indicating if sites with any missing values are to be removed or not.
#' @param FILTER_CONTEXT Flag indicating if only CG probes are to be kept.
#' @param FILTER_SNP Flag indicating if annotated SNPs are to be removed from the list of sites according to RnBeads' SNP list. Or as the sites
#'                  specified in \code{snp.list}.
#' @param dist.snps Flag indicating if SNPs are to removed by determining if the pairwise differences between the CpGs in the samples are trimodally
#'                  distributed as it is frequently found around SNPs.                  
#' @param snp.list Path to a file containing CpG IDs of known SNPs to be removed from the analysis, if \code{FILTER_SNP} is \code{TRUE}.
#' @param FILTER_SOMATIC Flag indicating if only somatic probes are to be kept.
#' @param FILTER_CROSS_REACTIVE Flag indicating if sites showing cross reactivity on the array are to be removed.
#' @param remove.ICA Flag indicating if independent component analysis is to be executed to remove potential confounding factor.
#'             If \code{TRUE},conf.fact.ICA needs to be specified.
#' @param conf.fact.ICA A vector of column names in the sample annotation sheet representing potential confounding factors.
#' @param ica.setting Named vector of settings passed to run.rnb.ica. Options are \code{nmin, nmax, thres.sd, alpha.fact, save.report, alpha.feat, type, ncores}. See 
#'             \code{\link{run.rnb.ica}} for further details. NULL indicates the default setting.
#' @param execute.lump Flag indicating if the LUMP algorithm is to be used for estimating the amount of immune cells in a particular sample.
#' @return A list with four elements: \itemize{
#'           \item quality.filter The indices of the sites that survived quality filtering
#' }
#' @export
prepare_data<-function(
		RNB_SET, 
		WORK_DIR=getwd(),
		analysis.name="analysis",
		SAMPLE_SELECTION_COL=NA,
		SAMPLE_SELECTION_GREP=NA,
		PHENO_COLUMNS=NA,
		ID_COLUMN=rnb.getOption("identifiers.column"),
		NORMALIZATION="none",
		REF_CT_COLUMN=NA,
		REF_RNB_SET=NULL,
		REF_RNB_CT_COLUMN=NA,
		PREPARE_TRUE_PROPORTIONS=FALSE,
		TRUE_A_TOKEN=NA,
		HOUSEMAN_A_TOKEN=NA,
		ESTIMATE_HOUSEMAN_PROP=FALSE,
		FILTER_BEADS=!is.null(RNB_SET@covg.sites),
		MIN_N_BEADS=3,
		FILTER_INTENSITY=inherits(RNB_SET, "RnBeadRawSet"),
		MIN_INT_QUANT = 0.01,
		MAX_INT_QUANT = 0.99, 
		FILTER_NA=TRUE,
		FILTER_CONTEXT=TRUE,
		FILTER_SNP=TRUE,
		FILTER_SOMATIC=TRUE,
		FILTER_CROSS_REACTIVE=T,
		remove.ICA=F,
		conf.fact.ICA=NULL,
		ica.setting=NULL,
		snp.list=NULL,
		execute.lump=FALSE,
		dist.snps=FALSE
){
	suppressPackageStartupMessages(require(RnBeads))

	OUTPUTDIR <- file.path(WORK_DIR, analysis.name)
	if(!file.exists(OUTPUTDIR)){
	  dir.create(OUTPUTDIR)
	}
	log.file <- file.path(OUTPUTDIR,"analysis.log")
	if(!file.exists(log.file)){
	  if(logger.isinitialized()){
	    logger.close()
	    logger.start(fname=log.file)
	  }else{
	    logger.start(fname=log.file)
	  }
	}
	
	################################# PREPARE THE PARAMETER TUNING RUN ############################################
	
	if(is.character(RNB_SET)){
		rnb.set<-load.rnb.set(RNB_SET)
	}else if(inherits(RNB_SET,"RnBSet")){
		rnb.set<-RNB_SET
	}
	
	############################### SELECTION OF SAMPLES
	
	#meth.ref<-meth(rnb.set.comb, row.names=TRUE)[,is.na(rnb.set.comb@pheno$diseaseState)]
	
	#pd.ref<-pheno(rnb.set.comb)[is.na(rnb.set.comb@pheno$diseaseState),]
	
	
	if(!is.na(SAMPLE_SELECTION_COL)){
		
		SAMP_SUBS<-grep(SAMPLE_SELECTION_GREP, pheno(rnb.set)[[SAMPLE_SELECTION_COL]])
		
		rms<-setdiff(1:length(samples(rnb.set)), SAMP_SUBS)
		
		rnb.set<-remove.samples(rnb.set, rms)
	}
	
	############################### NORMALIZATION
	if(!(NORMALIZATION %in% c("none"))){
		if(NORMALIZATION=="illumina"){
			rnb.set<-rnb.execute.normalization(rnb.set, method="illumina", bgcorr.method="none")
		}else if(NORMALIZATION=="dasen"){
			rnb.set<-rnb.execute.normalization(rnb.set, method="wm.dasen", bgcorr.method="none")
		}else if(NORMALIZATION=="noob"){
			rnb.set<-rnb.execute.normalization(rnb.set, method="none", bgcorr.method="methylumi.noob")
		}
	}
	############################### PREPARATION OF THE REFERENCE DATA
	
	meth.rnb<-meth(rnb.set)
	pd<-pheno(rnb.set)
	if(!is.na(REF_CT_COLUMN)){
		subs<-is.na(pd[[REF_CT_COLUMN]])
	}else{
		subs<-1:nrow(pd)
	}
	if(!is.na(PHENO_COLUMNS)){
		pheno.data<-pd[subs,PHENO_COLUMNS,drop=FALSE]
		save(pheno.data, file=sprintf("%s/pheno.RData", OUTPUTDIR))
	}
	
	if(!is.null(ID_COLUMN)){
		sample_ids<-pd[subs,ID_COLUMN]
		saveRDS(sample_ids, file=sprintf("%s/sample_ids.RDS", OUTPUTDIR))	
	}
	
	if(!is.null(REF_RNB_SET) && !is.na(REF_RNB_CT_COLUMN)){
		
		rnb.set.ref<-load.rnb.set(REF_RNB_SET)
		meth.rnb.ref<-meth(rnb.set.ref)
		pd.ref<-pheno(rnb.set.ref)
		
		ct<-pd.ref[[REF_RNB_CT_COLUMN]]
		nnas<-!is.na(ct)
		ct<-ct[nnas]
		
		meth.ref.ct<-lapply(unique(ct), function(cn) rowMeans(meth.rnb.ref[,nnas][,ct==cn]))
		trueT<-do.call("cbind", meth.ref.ct)
		colnames(trueT)<-unique(ct)
		
		ref.set.save <- file.path(OUTPUTDIR,"RefSet")
		if(!file.exists(ref.set.save)) dir.create(ref.set.save)
		save(trueT, file=sprintf("%s/data.set.RData", ref.set.save))
		save(pd.ref, file=sprintf("%s/pheno.RData", ref.set.save))
		
	}else if(!is.na(REF_CT_COLUMN)){
		
		ct<-pd[[REF_CT_COLUMN]]
		nnas<-!is.na(ct)
		ct<-ct[nnas]
		
		meth.ref.ct<-lapply(unique(ct), function(cn) rowMeans(meth.rnb[,nnas][,ct==cn]))
		trueT<-do.call("cbind", meth.ref.ct)
		colnames(trueT)<-unique(ct)
		
		save(trueT, file=sprintf("%s/trueT.RData", OUTPUTDIR))
	}
	
	############################### TRUE PROPORTIONS 
	if(PREPARE_TRUE_PROPORTIONS){
		if(!is.na(TRUE_A_TOKEN)){
			
			trueA<-t(na.omit(pd[subs,grep(TRUE_A_TOKEN, colnames(pd))]))
			trueA <- apply(trueA,c(1,2),as.numeric)
			rownames(trueA)<-gsub(TRUE_A_TOKEN, "", colnames(pd)[grep(TRUE_A_TOKEN, colnames(pd))])
			
			save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
			
		}
		
		if(!is.na(HOUSEMAN_A_TOKEN)){
			
			trueA<-t(na.omit(pd[,grep(HOUSEMAN_A_TOKEN, colnames(pd))]))
			trueA <- apply(trueA,c(1,2),as.numeric)
			rownames(trueA)<-gsub(HOUSEMAN_A_TOKEN, "", colnames(pd)[grep(HOUSEMAN_A_TOKEN, colnames(pd))])
			
			save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
			
		}else if(ESTIMATE_HOUSEMAN_PROP){
			
			if(!is.na(REF_RNB_CT_COLUMN)){
				rnb.set.ref<-remove.samples(rnb.set.ref, which(is.na(pheno(rnb.set.ref)[,REF_RNB_CT_COLUMN])))
				rnb.set<-combine(rnb.set, rnb.set.ref)
				REF_CT_COLUMN<-REF_RNB_CT_COLUMN
			}
			print("Estimating proportions using the Houseman et al, 2012 method")
			res<-estimateProportionsCP(rnb.set, REF_CT_COLUMN, NA, 2000, full.output = TRUE)
			
			trueA<-t(res$contributions.nonneg)
			trueA[trueA<1e-5]<-0
			
			trueA<-sweep(trueA, 2, colSums(trueA),"/")
			
			save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
		}
	}
	
	if(!is.na(REF_CT_COLUMN)){
		meth.data<-meth.rnb[,subs]
	}else{
		meth.data<-meth.rnb
	}
	colnames(meth.data)<-NULL
	
	save(meth.data, file=sprintf("%s/data.set.RData", OUTPUTDIR))
	
	if(execute.lump){
	  lump.est <- rnb.execute.lump(rnb.set)
	  rnb.set <- addPheno(rnb.set,as.vector(lump.est),"LUMP_estimate")
	}
	
	####################### FILTERING
	################################# QUALITY FILTERING ######################################
	FILTER_QUALITY<- FILTER_BEADS || FILTER_INTENSITY
	
	if(FILTER_QUALITY){
		M.raw<-RnBeads:::M(rnb.set, row.names=TRUE)
		U.raw<-RnBeads:::U(rnb.set, row.names=TRUE)
		b.raw<-RnBeads:::covg(rnb.set, row.names=TRUE)
		
		if(!is.na(REF_CT_COLUMN)){
			
			M.raw<-M.raw[,subs,drop=FALSE]
			U.raw<-U.raw[,subs,drop=FALSE]
			b.raw<-b.raw[,subs, drop=FALSE]
			
		}
		
		saveRDS(M.raw, file.path(OUTPUTDIR, "Mint.RDS"))
		saveRDS(U.raw, file.path(OUTPUTDIR, "Uint.RDS"))
		saveRDS(b.raw, file.path(OUTPUTDIR, "Nbeads.RDS"))
		
		qual.filter<-filter.quality(rnb.set,beads=FILTER_BEADS, min.beads=MIN_N_BEADS, intensity = FILTER_INTENSITY, 
		                            min.int.quant=MIN_INT_QUANT, max.int.quant=MAX_INT_QUANT,subs = subs)
		
		save(qual.filter, file=sprintf("%s/quality.filter.RData", OUTPUTDIR))
		
	}else{
		qual.filter<-1:nrow(rnb.set@meth.sites)
	}
	if(FILTER_NA){
	  qual.filter <- filter.nas(rnb.set,subs=subs,qual.filter)
	}
	########################################## ANNOTATION FILTERING ###################################################
	FILTER_ANNOTATION<-FILTER_CONTEXT || FILTER_SNP || FILTER_SOMATIC
	
	if(FILTER_ANNOTATION){
		
		annot.filter<-filter.annotation(rnb.set, context = FILTER_CONTEXT, snp = FILTER_SNP, snp.list = snp.list,
		                                somatic = FILTER_SOMATIC, qual.filter = qual.filter, dist.snps = dist.snps)
	
		save(annot.filter, file=sprintf("%s/annotation.filter.RData", OUTPUTDIR))
	
		}else{
		annot.filter<-1:nrow(rnb.set@meth.sites)
	}
	
	total.filter<-intersect(qual.filter, annot.filter)
	logger.info(paste("Removing",nsites(rnb.set)-length(total.filter),"sites, retaining ",length(total.filter)))
	rnb.set.f<-remove.sites(rnb.set, setdiff(1:nrow(rnb.set@meth.sites), total.filter))
	
	if(FILTER_CROSS_REACTIVE && inherits(rnb.set.f, "RnBeadSet")){
	  cross.reactive.filter <- rnb.execute.cross.reactive.removal(rnb.set.f)
	  logger.info(paste(length(cross.reactive.filter$filtered),"sites removed in cross-reactive filtering"))
	  rnb.set.f <- cross.reactive.filter$dataset
	}
	
	if(remove.ICA){
	  if(!inherits(rnb.set.f,"RnBSet")){
	    logger.error("ICA only applicable to RnBSet objects.")
	  }
	  logger.start("Removing confounding factors using ICA")
	  rnb.set.f <- run.rnb.ICA(rnb.set.f,conf.fact.ICA,out.folder=OUTPUTDIR,ica.setting=ica.setting)
	  logger.completed()
	}
	
	analysis_info<-list()
	
	analysis_info$QUALITY_FILTERING <- sprintf("%s%s%s%s%s%s",
	                                           ifelse(FILTER_BEADS, "Beads", ""),
	                                           ifelse(FILTER_INTENSITY, "Intensity", ""),
	                                           ifelse(FILTER_NA, "Missing", ""),
	                                           ifelse(FILTER_CONTEXT, "Context", ""),
	                                           ifelse(FILTER_SNP, "SNP", ""),
	                                           ifelse(FILTER_SOMATIC, "Somatic", "")
	)
	
	analysis_info$NORMALIZATION <- NORMALIZATION
	
	res <- list(quality.filter=qual.filter, annot.filter=annot.filter, total.filter=total.filter, rnb.set.filtered=rnb.set.f, info=analysis_info)
	if(exists("trueT")){
	  res$RefMeth <- trueT[total.filter,]
	}
	if(exists("trueA")){
	  res$RefProps <- trueA
	}
	return(res)
}

#' filter.quality
#' 
#' This functions filters the CpG sites in the given rnb.set for quality criteria specified in the arguments.
#' 
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} containing the CpG sites used for filtering, all well as intensity and
#'                 coverage informatio, if provided.
#' @param beads Flag indicating if sites not having more than \code{min.beads} number of beads in all samples are to be removed.
#' @param min.beads Integer specifying the minimum number of beads required for including the site in the analysis.
#' @param intensity Flag indicating if sites are to be filtered according to their intensity information.
#' @param min.int.quant Lower quantile of intensities which is to be removed.
#' @param max.int.quant Upper quantile of intensities which is to be removed.
#' @param subs Argument specifying the subset of samples to be used
#' @return The indices of the sites that are to be kept.
#' @details If \code{intensity} is set, those probes with lower/higher intensity in one of the channels than \code{min.int.quant}/
#'            \code{max.int.quant} are removed.
filter.quality<-function(
  rnb.set,
  beads=TRUE,
  min.beads,
  intensity=TRUE,
  min.int.quant=0.1,
  max.int.quant=0.95,
  subs
){
  
  annot <- annotation(rnb.set)
  
  qf<-1:nrow(annot)
  
  M.raw <- M(rnb.set, row.names=TRUE)
  U.raw <- U(rnb.set, row.names=TRUE)
  b.raw <- covg(rnb.set, row.names=TRUE)
  
  if(!is.null(subs)){
    M.raw<-M.raw[,subs,drop=FALSE]
    U.raw<-U.raw[,subs,drop=FALSE]
    b.raw<-b.raw[,subs, drop=FALSE]
  }
  
  if(beads){
    b.raw<-RnBeads:::covg(rnb.set, row.names=TRUE)
    qf.b<-which(rowSums(b.raw>=min.beads)==ncol(b.raw))
    logger.info(paste(length(setdiff(qf,qf.b)),"sites removed in bead count filtering."))
    qf<-intersect(qf, qf.b)
  }
  
  if(intensity){
    
    MplusU<-M.raw+U.raw
    
    hm450_ann <- annotation(rnb.set)
    
    MplusU.I<-MplusU[hm450_ann$Design=="I",]
    MplusU.II<-MplusU[hm450_ann$Design=="II",]
    
    MU.q001.I<-sort(as.numeric(MplusU.I))[ceiling(min.int.quant*nrow(MplusU.I)*ncol(MplusU.I))]
    MU.q099.I<-sort(as.numeric(MplusU.I))[ceiling(max.int.quant*nrow(MplusU.I)*ncol(MplusU.I))]
    
    MU.q001.II<-sort(as.numeric(MplusU.II))[ceiling(min.int.quant*nrow(MplusU.II)*ncol(MplusU.II))]
    MU.q099.II<-sort(as.numeric(MplusU.II))[ceiling(max.int.quant*nrow(MplusU.II)*ncol(MplusU.II))]
    
    MplusU.f<-matrix(FALSE, nrow=nrow(MplusU), ncol=ncol(MplusU))
    
    MplusU.f[hm450_ann$Design=="I",]<-MplusU.I>MU.q001.I & MplusU.I<MU.q099.I
    MplusU.f[hm450_ann$Design=="II",]<-MplusU.II>MU.q001.II & MplusU.II<MU.q099.II
    
    qf.MU<-which(rowSums(MplusU.f)==ncol(MplusU.f))
    logger.info(paste(length(setdiff(qf,qf.MU)),"sites removed in intensity filtering."))
    
    qf<-intersect(qf, qf.MU)
  }
  return(qf)
}

#' bigFF.row.apply
#' 
#' This routine applies a function to chunks of the dataset that is stored in disk
#' either with `ff` or `BigFfMat`.
#' 
#' @param mat A disk-based matrix of type \code{\link{ff}} or \code{\link{BigFfMat}}.
#' @param FUN The function to be applied to each chunk of the matrix. Should be a function
#'             that computes matrix statistics, such as \code{rowMeans}.
#' @param iter.count Number of chunks to be created. The larger this number, the slower
#'             the computation, but the lower the disk usage.
#' @param ... Further arguments passed to FUN.
#' @return A vector summarizing the results of the function. Final structure is determined
#'          by FUN.
#' @author Michael Scherer
#' @noRd
bigFF.row.apply <- function(mat,FUN,iter.count=1000,...){
  chunk.size <- floor(nrow(mat)/iter.count)
  iter <- 1
  res <- c()
  while(iter+chunk.size < nrow(mat)){
    chunk <- mat[iter:(iter+chunk.size-1),]
    res <- c(res,FUN(chunk,...))
    iter <- iter+chunk.size
#    print(paste(round((iter/nrow(mat))*100,2)," percent completed"))
  }
  chunk <- mat[iter:nrow(mat),]
  res <- c(res,FUN(chunk,...))
  return(res)
}

#' filter.quality.covg
#' 
#' This functions filters the CpG sites in the given rnb.set for quality criteria specified in the arguments.
#' 
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} containing the CpG sites used for filtering, all well as intensity and
#'                 coverage informatio, if provided.
#' @param min.covg Integer specifying the minimum number of beads required for including the site in the analysis.
#' @param min.covg.quant Lower quantile of coverage which is to be removed.
#' @param max.covg.quant Upper quantile of coverage which is to be removed.
#' @return The indices of the sites that are to be kept.
#' @details The probes with lower/higher coverage than \code{min.int.quant}/
#'            \code{max.int.quant} are removed.
filter.quality.covg <- function(
  rnb.set,
  min.covg,
  min.covg.quant=0.05,
  max.covg.quant=0.95
){
  
  qf <- 1:nsites(rnb.set)
  covg.data <- rnb.set@covg.sites  

  mins <- bigFF.row.apply(covg.data,rowMins,iter.count = 100,na.rm=T)
  qf.b <- which(mins >= min.covg)
  logger.info(paste(length(setdiff(qf,qf.b)),"sites removed in absolute coverage filtering."))
  qf <- intersect(qf, qf.b)
  quants <- bigFF.row.apply(covg.data,quantile,iter.count=100,probs=c(min.covg.quant,max.covg.quant),na.rm=T)
  lower.covg <- mean(quants[seq(1,length(quants)-1,by=2)])
  upper.covg <- mean(quants[seq(2,length(quants),by=2)])
  
  maxs <- bigFF.row.apply(covg.data,rowMaxs,iter.count = 100,na.rm=T)         
  qf.covg <- which((mins>lower.covg & maxs<upper.covg))
  logger.info(paste(length(setdiff(qf,qf.covg)),"sites removed in quantile coverage filtering."))
    
  qf<-intersect(qf, qf.covg)
  return(qf)
}

#' filter.nas
#' 
#' This function removes any site that contains a missing methylation value in a sample.
#' 
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} containing methylation information.
#' @param qf.qual Vector of indices that survided quality filtering
#' @param subs Optional argument specifying the subset of samples to be used
#' @return Vector of indices that survived NA filtering.
filter.nas <- function(rnb.set,
                       subs,
                       qf.qual){
  meth.data <- meth(rnb.set)[,subs,drop=FALSE]
  na.filter<-which(rowSums(is.na(meth.data))<1)
  logger.info(paste(length(setdiff(qf.qual,na.filter)),"sites removed in NA filtering"))
  qf.qual<-intersect(qf.qual, na.filter)
  return(qf.qual)
}

#' filter.nas.biseq
#' 
#' This function removes any site that contains a missing methylation value in a sample.
#' 
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} containing methylation information.
#' @param qf.qual Vector of indices that survided quality filtering
#' @param subs Optional argument specifying the subset of samples to be used
#' @return Vector of indices that survived NA filtering.
filter.nas.biseq <- function(rnb.set,
                       subs,
                       qf.qual){
  meth.data <- rnb.set@meth.sites
  na.filter <- which(bigFF.row.apply(meth.data,function(x,s.subset){
    rowSums(is.na(x[,s.subset]))<1
  },iter.count=100,s.subset=subs))
  logger.info(paste(length(setdiff(qf.qual,na.filter)),"sites removed in NA filtering"))
  qf.qual<-intersect(qf.qual, na.filter)
  return(qf.qual)
}


#' filter.annotation
#' 
#' This function removes sites in SNP, sex-chromosomal, or non-CG context.
#' 
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} containg annotation information.
#' @param snp Flag indicating if snps are to be removed. Either SNPs annotated in the RnBeads annotation object of specified as an
#'            additional file with one SNP identifier per row in \code{snp.list}.
#' @param snp.list Path to a file containing known SNPs. One SNP identifier should be present per row.
#' @param somatic Flag indicating if only somatic probes are to be kept.
#' @param context Flag indicating if probes in non-CpG context are to be removed.
#' @param qual.filter Vector of indices removed during quality filtering.
#' @param dist.snps Flag indicating of potential SNPs are to be determined by selecting those sites that have trimodally
#'           distributed (differences 0, 0.5 and 1.0) pairwise differences across the samples.
#' @return A vector of indices of sites surviving the annotation filter criteria.
filter.annotation<-function(
  rnb.set,
  snp=TRUE,
  snp.list,
  somatic=TRUE,
  context=TRUE,
  qual.filter=NULL,
  dist.snps=FALSE)
{
  annot<-annotation(rnb.set)
  
  if(!is.null(qual.filter)){
    probe.ind.filtered <- qual.filter
  }else{
    probe.ind.filtered <- 1:nrow(annot)
  }
  
  if(snp){
    snp.filter <- probe.ind.filtered
    if(is.null(snp.list)){
      snp.filter<-which(!annot$`SNPs 3` & !annot$`SNPs 5` & !annot$`SNPs Full`)
    }else{
      snps <- readLines(snp.list)
      snp.filter <- which(!(row.names(annot) %in% snps))
    }
    if(dist.snps){
      require("LaplacesDemon")
      meth.data <- meth(rnb.set)
      meth.data <- meth(rnb.set)
      pair.dist.rand <- apply(meth.data,1,function(snp){
        pair.snp <- c()	
        for(snp2 in snp){
          diff <- abs(snp-snp2)
          pair.snp <- c(pair.snp,diff)
        }
        pair.snp
      })
      modes.rand <- unlist(lapply(pair.dist.rand,function(snp)is.trimodal(unlist(snp))))
      logger.info(paste("Filtered",length(setdiff(snp.filter,which(!modes.rand))),"CpGs in distribution SNP filtering"))
      snp.filter <- intersect(snp.filter,which(!modes.rand))
      rm(meth.data)
    }
    logger.info(paste(length(setdiff(probe.ind.filtered,snp.filter)),"sites removed in SNP filtering"))
    probe.ind.filtered<-intersect(probe.ind.filtered, snp.filter)
  }
  
  if(somatic){
    somatic.filter<-which(!annot$Chromosome %in% c("chrX", "chrY"))
    logger.info(paste(length(setdiff(probe.ind.filtered,somatic.filter)),"sites removed in somatic sites filtering"))
    probe.ind.filtered<-intersect(probe.ind.filtered, somatic.filter)
  }
  
  if(context && inherits(rnb.set, "RnBeadSet")){
    context.filter<-grep("cg", rownames(annot))
    logger.info(paste(length(setdiff(probe.ind.filtered,context.filter)),"sites removed in CG context filtering"))
    probe.ind.filtered<-intersect(probe.ind.filtered, context.filter)
  }
  return(probe.ind.filtered)
}

#' filter.annotation.biseq
#' 
#' This function removes sites in SNP or sex-chromosomal context.
#' 
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} containg annotation information.
#' @param snp Flag indicating if snps are to be removed. Either SNPs annotated in the RnBeads annotation object of specified as an
#'            additional file with one SNP identifier per row in \code{snp.list}.
#' @param snp.list Path to a file containing known SNPs. One SNP identifier should be present per row.
#' @param somatic Flag indicating if only somatic probes are to be kept.
#' @param qual.filter Vector of indices removed during quality filtering.
#' @return A vector of indices of sites surviving the annotation filter criteria.
#' @details This functions uses information on the sites on the chip and transfers this knowledge to BS data sets.
filter.annotation.biseq<-function(
  rnb.set,
  snp=TRUE,
  snp.list,
  somatic=TRUE,
  qual.filter=NULL)
{
  annot <- annotation(rnb.set)

  if(!is.null(qual.filter)){
    probe.ind.filtered <- qual.filter
  }else{
    probe.ind.filtered <- 1:nrow(annot)
  }
  
  if(snp){
    if(is.null(snp.list)){
      snp.filter<-which(is.na(annot$SNPs))
    }else{
      snps <- read.table(snp.list,sep="\t")
	  snps <- GRanges(Rle(snps[,1]),IRanges(start=as.numeric(snps[,2]),end=as.numeric(snps[,2])+1))
      anno.granges <- GRanges(Rle(annot$Chromosome),IRanges(start=annot$Start,end=annot$End))
	  op <- findOverlaps(anno.granges,snps)
      snp.filter <- queryHits(op)
    }
    logger.info(paste(length(setdiff(probe.ind.filtered,snp.filter)),"sites removed in SNP filtering"))
    probe.ind.filtered<-intersect(probe.ind.filtered, snp.filter)
  }
  
  if(somatic){
    somatic.filter<-which(!annot$Chromosome %in% c("chrX", "chrY"))
    logger.info(paste(length(setdiff(probe.ind.filtered,somatic.filter)),"sites removed in somatic sites filtering"))
    probe.ind.filtered<-intersect(probe.ind.filtered, somatic.filter)
  }  
  return(probe.ind.filtered)
}

#' prepare_data_BS
#' 
#' This functions prepares sequencing data sets for a MeDeCom run.
#' 
#' @param RNB_SET An object of type \code{\link[RnBeads]{RnBiseqSet-class}} for which analysis is to be performed.
#' @param WORK_DIR A path to a existing directory, in which the results are to be stored
#' @param analysis.name A string representing the dataset for which analysis is to be performed. Only used to create a folder with a 
#'                 descriptive name of the analysis.
#' @param SAMPLE_SELECTION_COL A column name in the phenotypic table of \code{RNB_SET} used to selected a subset of samples for
#'                 analysis that contain the string given in \code{SAMPLE_SELECTION_GREP}.
#' @param SAMPLE_SELECTION_GREP A string used for selecting samples in the column \code{SAMPLE_SELECTION_COL}.
#' @param REF_CT_COLUMN Column name in \code{RNB_SET} used to extract methylation information on the reference cell types.
#' @param PHENO_COLUMNS Vector of column names in the phenotypic table of \code{RNB_SET} that is kept and exported for further 
#'                 exploration.
#' @param PREPARE_TRUE_PROPORTIONS Flag indicating if true proportions are either available in \code{RNB_SET} or to be estimated 
#'                          with Houseman's reference-based deconvolution approach.
#' @param TRUE_A_TOKEN String present in the column names of \code{RNB_SET} used for selecting the true proportions of the corresponding
#'                      cell types.
#' @param HOUSEMAN_A_TOKEN Similar to \code{TRUE_A_TOKEN}, but not containing the true proportions, rather the estimated proportions
#'                      by Houseman's method.
#' @param ID_COLUMN Sample-specific ID column name in \code{RNB_SET}
#' @param FILTER_COVERAGE Flag indicating, if site-filtering based on coverage is to be conducted.
#' @param MIN_COVERAGE Minimum number of reads required in each sample for the site to be considered for adding to MeDeCom.
#' @param MIN_COVG_QUANT Lower quantile of coverages. Values lower than this value will be ignored for analysis.
#' @param MAX_COVG_QUANT Upper quantile of coverages. Values higher than this value will be ignored for analysis.
#' @param FILTER_NA Flag indicating if sites with any missing values are to be removed or not.
#' @param FILTER_SNP Flag indicating if annotated SNPs are to be removed from the list of sites according to RnBeads' SNP list.
#' @param snp.list Path to a file containing positions of known SNPs to be removed from the analysis, if \code{FILTER_SNP} is \code{TRUE}. The coordinates must be the 
#'             provided in the same genome assembly as RNB_SET. The file must be a tab-separated value (tsv) file with only one header line an the following meaning of 
#'             the rows: 1st row: chromosome, 2nd row: position of the SNP on the chromosome
#' @param FILTER_SOMATIC Flag indicating if only somatic probes are to be kept.
#' @return A list with four elements: \itemize{
#'           \item quality.filter The indices of the sites that survived quality filtering
#' }
#' @export
prepare_data_BS <- function(
		RNB_SET, 
		WORK_DIR=getwd(),
		analysis.name="analysis",
		SAMPLE_SELECTION_COL=NA,
		SAMPLE_SELECTION_GREP=NA,
		REF_CT_COLUMN=NA,
		PHENO_COLUMNS=NA,
		PREPARE_TRUE_PROPORTIONS=FALSE,
		TRUE_A_TOKEN=NA,
		HOUSEMAN_A_TOKEN=NA,
		ID_COLUMN=rnb.getOption("identifiers.column"),
		FILTER_COVERAGE = hasCovg(RNB_SET),
		MIN_COVERAGE=5,
		MIN_COVG_QUANT=0.05,
		MAX_COVG_QUANT=0.95,
		FILTER_NA=TRUE,
		FILTER_SNP=TRUE,
		snp.list=NULL,
		FILTER_SOMATIC=TRUE
){
	suppressPackageStartupMessages(require(RnBeads))

	OUTPUTDIR <- file.path(WORK_DIR, analysis.name)
	if(!file.exists(OUTPUTDIR)){
	  dir.create(OUTPUTDIR)
	}
	log.file <- file.path(OUTPUTDIR,"analysis.log")
	if(!file.exists(log.file)){
	  if(logger.isinitialized()){
	    logger.close()
	    logger.start(fname=log.file)
	  }else{
	    logger.start(fname=log.file)
	  }
	}
	if(is.character(RNB_SET)){
		rnb.set<-load.rnb.set(RNB_SET)
	}else if(inherits(RNB_SET,"RnBSet")){
		rnb.set<-RNB_SET
	}	
	if(!is.na(SAMPLE_SELECTION_COL)){
		
		SAMP_SUBS<-grep(SAMPLE_SELECTION_GREP, pheno(rnb.set)[[SAMPLE_SELECTION_COL]])
		
		rms<-setdiff(1:length(samples(rnb.set)), SAMP_SUBS)
		
		rnb.set<-remove.samples(rnb.set, rms)
	}
		
	meth.rnb <- rnb.set@meth.sites
	pd<-pheno(rnb.set)
	if(!is.na(REF_CT_COLUMN)){
	  subs<-is.na(pd[[REF_CT_COLUMN]])
	}else{
	  subs<-1:nrow(pd)
	}
	if(!is.na(PHENO_COLUMNS)){
		pheno.data<-pd[,PHENO_COLUMNS,drop=FALSE]
		save(pheno.data, file=sprintf("%s/pheno.RData", OUTPUTDIR))
	}
	
	if(!is.null(ID_COLUMN)){
		sample_ids<-pd[,ID_COLUMN]
		saveRDS(sample_ids, file=sprintf("%s/sample_ids.RDS", OUTPUTDIR))	
	}
	if(PREPARE_TRUE_PROPORTIONS){
	  if(!is.na(TRUE_A_TOKEN)){
	    
	    trueA<-t(na.omit(pd[subs,grep(TRUE_A_TOKEN, colnames(pd))]))
	    trueA <- apply(trueA,c(1,2),as.numeric)
	    rownames(trueA)<-gsub(TRUE_A_TOKEN, "", colnames(pd)[grep(TRUE_A_TOKEN, colnames(pd))])
	    
	    save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
	    
	  }
	  
	  if(!is.na(HOUSEMAN_A_TOKEN)){
	    
	    trueA<-t(na.omit(pd[,grep(HOUSEMAN_A_TOKEN, colnames(pd))]))
	    trueA <- apply(trueA,c(1,2),as.numeric)
	    rownames(trueA)<-gsub(HOUSEMAN_A_TOKEN, "", colnames(pd)[grep(HOUSEMAN_A_TOKEN, colnames(pd))])
	    
	    save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
	    
	  }
	}
	if(!is.na(REF_CT_COLUMN)){
	  ct<-pd[[REF_CT_COLUMN]]
	  nnas<-!is.na(ct)
	  ct<-ct[nnas]
	  
	  meth.ref.ct<-lapply(unique(ct), function(cn) rowMeans(meth.rnb[,nnas][,ct==cn]))
	  trueT<-do.call("cbind", meth.ref.ct)
	  colnames(trueT)<-unique(ct)
	  
	  save(trueT, file=sprintf("%s/trueT.RData", OUTPUTDIR))
	}
	if(!is.na(REF_CT_COLUMN)){
	  meth.data<-meth.rnb[,subs]
	}else{
	  meth.data<-meth.rnb
	}
	colnames(meth.data) <- NULL	
	save(meth.data, file=sprintf("%s/data.set.RData", OUTPUTDIR))
	if(FILTER_COVERAGE){
		qual.filter <- filter.quality.covg(rnb.set, min.covg = MIN_COVERAGE,  
		                            min.covg.quant = MIN_COVG_QUANT, max.covg.quant=MAX_COVG_QUANT)		
		save(qual.filter, file=sprintf("%s/quality.filter.RData", OUTPUTDIR))		
	}else{
		qual.filter<-1:nrow(rnb.set@meth.sites)
	}
	if(FILTER_NA){
	  qual.filter <- filter.nas.biseq(rnb.set,subs=1:length(samples(rnb.set)),qual.filter)
	}
	FILTER_ANNOTATION <- FILTER_SNP || FILTER_SOMATIC
	
	if(FILTER_ANNOTATION){
		
		annot.filter <- filter.annotation.biseq(rnb.set, snp = FILTER_SNP, snp.list = snp.list,
		                                somatic = FILTER_SOMATIC, qual.filter = qual.filter)
	
		save(annot.filter, file=sprintf("%s/annotation.filter.RData", OUTPUTDIR))
	
		}else{
		annot.filter<-1:nrow(rnb.set@meth.sites)
	}
	
	total.filter<-intersect(qual.filter, annot.filter)
	logger.info(paste("Removing",nsites(rnb.set)-length(total.filter),"sites, retaining ",length(total.filter)))
	rnb.set.f<-remove.sites(rnb.set, setdiff(1:nrow(rnb.set@meth.sites), total.filter))
	
	res <- list(quality.filter=qual.filter, annot.filter=annot.filter, total.filter=total.filter, rnb.set.filtered=rnb.set.f)
	if(exists("trueT")){
	  res$RefMeth <- trueT[total.filter,]
	}
	if(exists("trueA")){
	  res$RefProps <- trueA
	}
	
	return(res)
}
