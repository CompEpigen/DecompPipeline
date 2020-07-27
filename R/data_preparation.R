#
# Data preparation functions for Illumina BeadArray and Bisulfite sequencing data
# Authors: Pavlo Lutsik and Michael Scherer
#

#################################################################################################################################
#' FUNCTIONS
#################################################################################################################################

#' prepare.data
#' 
#' This functions prepares Illumina BeadChip data for a MeDeCom/EDec/RefFreeCellMix run.
#' 
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} for which analysis is to be performed.
#' @param work.dir A path to a existing directory, in which the results are to be stored
#' @param analysis.name A string representing the dataset for which analysis is to be performed. Only used to create a folder with a 
#'                 descriptive name of the analysis.
#' @param sample.selection.col A column name in the phenotypic table of \code{rnb.set} used to selected a subset of samples for
#'                 analysis that contain the string given in \code{sample.selection.col}.
#' @param sample.selection.grep A string used for selecting samples in the column \code{sample.selection.grep}.
#' @param pheno.columns Vector of column names in the phenotypic table of \code{rnb.set} that is kept and exported for further 
#'                 exploration.
#' @param id.column Sample-specific ID column name in \code{rnb.set}
#' @param normalization Normalization method to be performed before employing MeDeCom. Can be one of \code{"none", "dasen", "illumina", "noob", "bmiq"}.
#' @param ref.ct.column Column name in \code{rnb.set} used to extract methylation information on the reference cell types.
#' @param ref.rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} containing methylation information on reference cell types.
#' @param ref.rnb.ct.column Column name in \code{ref.rnb.set} used to extract methylation information on the reference cell types.
#' @param prepare.true.proportions Flag indicating if true proportions are either available in \code{rnb.set} or to be estimated 
#'                          with Houseman's reference-based deconvolution approach.
#' @param true.A.token String present in the column names of \code{rnb.set} used for selecting the true proportions of the corresponding
#'                      cell types.
#' @param houseman.A.token Similar to \code{true.A.token}, but not containing the true proportions, rather the estimated proportions
#'                      by Houseman's method.
#' @param estimate.houseman.prop If neither \code{true.A.token} nor \code{houseman.A.token} are given, the proportions of the reference
#'                      cell type are estimated with Houseman's approach.
#' @param filter.beads Flag indicating, if site-filtering based on the number of beads available is to be conducted.
#' @param min.n.beads Minimum number of beads required in each sample for the site to be considered for adding to MeDeCom.
#' @param filter.intensity  Flag indicating if sites should be removed according to the signal intensities (the lowest and highest quantiles
#'                      given by \code{min.int.quant} and \code{max.int.quant}). Note that all sites are removed that have a value outside of
#'                      the provided quantile range in either of the channels and in any of the samples. 
#' @param min.int.quant Lower quantile of intensities which is to be removed.
#' @param max.int.quant Upper quantile of intensities which is to be removed.
#' @param filter.na Flag indicating if sites with any missing values are to be removed or not.
#' @param filter.context Flag indicating if only CG probes are to be kept.
#' @param filter.snp Flag indicating if annotated SNPs are to be removed from the list of sites according to RnBeads' SNP list. Or as the sites
#'                  specified in \code{snp.list}.
#' @param dist.snps Flag indicating if SNPs are to removed by determining if the pairwise differences between the CpGs in the samples are trimodally
#'                  distributed as it is frequently found around SNPs.                  
#' @param snp.list Path to a file containing CpG IDs of known SNPs to be removed from the analysis, if \code{filter.snp} is \code{TRUE}.
#' @param filter.sex.chromosomes Flag indicating if only somatic probes are to be kept.
#' @param filter.cross.reactive Flag indicating if sites showing cross reactivity on the array are to be removed.
#' @param remove.ICA Flag indicating if independent component analysis is to be executed to remove potential confounding factor.
#'             If \code{TRUE},conf.fact.ICA needs to be specified.
#' @param conf.fact.ICA A vector of column names in the sample annotation sheet representing potential confounding factors.
#' @param ica.setting Named vector of settings passed to run.rnb.ica. Options are \code{nmin, nmax, thres.sd, alpha.fact, save.report, alpha.feat, type, ncores}. See 
#'             \code{\link{run.rnb.ICA}} for further details. NULL indicates the default setting.
#' @param execute.lump Flag indicating if the LUMP algorithm is to be used for estimating the amount of immune cells in a particular sample.
#' @return A list with four elements: \itemize{
#'           \item quality.filter The indices of the sites that survived quality filtering
#' }
#' @export
#' @author Michael Scherer, Pavlo Lutsik
#' @import RnBeads
#' @import R.utils
prepare.data<-function(
		rnb.set, 
		work.dir=getwd(),
		analysis.name="analysis",
		sample.selection.col=NA,
		sample.selection.grep=NA,
		pheno.columns=NA,
		id.column=rnb.getOption("identifiers.column"),
		normalization="none",
		ref.ct.column=NA,
		ref.rnb.set=NULL,
		ref.rnb.ct.column=NA,
		prepare.true.proportions=FALSE,
		true.A.token=NA,
		houseman.A.token=NA,
		estimate.houseman.prop=FALSE,
		filter.beads=!is.null(rnb.set@covg.sites),
		min.n.beads=3,
		filter.intensity=inherits(rnb.set, "RnBeadRawSet"),
		min.int.quant = 0.001,
		max.int.quant = 0.999, 
		filter.na=TRUE,
		filter.context=TRUE,
		filter.snp=TRUE,
		filter.sex.chromosomes=TRUE,
		filter.cross.reactive=T,
		remove.ICA=F,
		conf.fact.ICA=NULL,
		ica.setting=NULL,
		snp.list=NULL,
		execute.lump=FALSE,
		dist.snps=FALSE
){

	OUTPUTDIR <- file.path(work.dir, analysis.name)
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
	
	if(is.character(rnb.set)){
		rnb.set<-load.rnb.set(rnb.set)
	}else if(inherits(rnb.set,"RnBSet")){
		rnb.set<-rnb.set
	}
	
	if(!is.null(conf.fact.ICA)){
		if(!all(conf.fact.ICA%in%colnames(pheno(rnb.set)))){
			missing.names <- conf.fact.ICA[!conf.fact.ICA%in%colnames(pheno(rnb.set))]
			conf.fact.ICA <- conf.fact.ICA[conf.fact.ICA%in%colnames(pheno(rnb.set))]
			logger.warning(paste("Missing confounding factor names",missing.names,", only using",conf.fact.ICA))
		}
	}

	############################### SELECTION OF SAMPLES
	
	#meth.ref<-meth(rnb.set.comb, row.names=TRUE)[,is.na(rnb.set.comb@pheno$diseaseState)]
	
	#pd.ref<-pheno(rnb.set.comb)[is.na(rnb.set.comb@pheno$diseaseState),]
	
	
	if(!is.na(sample.selection.col)){
		
		SAMP_SUBS<-grep(sample.selection.grep, pheno(rnb.set)[[sample.selection.col]])
		
		rms<-setdiff(1:length(samples(rnb.set)), SAMP_SUBS)
		
		rnb.set<-remove.samples(rnb.set, rms)
	}
	
	############################### normalization
	if(!(normalization %in% c("none"))){
		if(normalization=="illumina"){
			rnb.set<-rnb.execute.normalization(rnb.set, method="illumina", bgcorr.method="none")
		}else if(normalization=="dasen"){
			rnb.set<-rnb.execute.normalization(rnb.set, method="wm.dasen", bgcorr.method="none")
		}else if(normalization=="noob"){
			rnb.set<-rnb.execute.normalization(rnb.set, method="none", bgcorr.method="methylumi.noob")
		}else if(normalization=="bmiq"){
			rnb.set<-rnb.execute.normalization(rnb.set, method="bmiq", bgcorr.method="none")
		}else{
			logger.warning("No valid normalization method specified (illumina, dasen, noob, bmiq allowed), skipping normalization")	
		}
	}
	############################### PREPARATION OF THE REFERENCE DATA
	
	meth.rnb<-meth(rnb.set)
	pd<-pheno(rnb.set)
	if(!is.na(ref.ct.column)){
		subs<-is.na(pd[[ref.ct.column]])
	}else{
		subs<-1:nrow(pd)
	}
	if(!is.na(pheno.columns)){
		pheno.data<-pd[subs,pheno.columns,drop=FALSE]
		save(pheno.data, file=sprintf("%s/pheno.RData", OUTPUTDIR))
	}
	
	if(!is.null(id.column)){
		sample_ids<-pd[subs,id.column]
		saveRDS(sample_ids, file=sprintf("%s/sample_ids.RDS", OUTPUTDIR))	
	}
	
	if(!is.null(ref.rnb.set) && !is.na(ref.rnb.ct.column)){
		
		rnb.set.ref<-load.rnb.set(ref.rnb.set)
		meth.rnb.ref<-meth(rnb.set.ref)
		pd.ref<-pheno(rnb.set.ref)
		
		ct<-pd.ref[[ref.rnb.ct.column]]
		nnas<-!is.na(ct)
		ct<-ct[nnas]
		
		meth.ref.ct<-lapply(unique(ct), function(cn) rowMeans(meth.rnb.ref[,nnas][,ct==cn]))
		trueT<-do.call("cbind", meth.ref.ct)
		colnames(trueT)<-unique(ct)
		
		ref.set.save <- file.path(OUTPUTDIR,"RefSet")
		if(!file.exists(ref.set.save)) dir.create(ref.set.save)
		save(trueT, file=sprintf("%s/data.set.RData", ref.set.save))
		save(pd.ref, file=sprintf("%s/pheno.RData", ref.set.save))
		
	}else if(!is.na(ref.ct.column)){
		
		ct<-pd[[ref.ct.column]]
		nnas<-!is.na(ct)
		ct<-ct[nnas]
		
		meth.ref.ct<-lapply(unique(ct), function(cn) rowMeans(meth.rnb[,nnas][,ct==cn]))
		trueT<-do.call("cbind", meth.ref.ct)
		colnames(trueT)<-unique(ct)
		
		save(trueT, file=sprintf("%s/trueT.RData", OUTPUTDIR))
	}
	
	############################### TRUE PROPORTIONS 
	if(prepare.true.proportions){
		if(!is.na(true.A.token)){
			
			trueA<-t(na.omit(pd[subs,grep(true.A.token, colnames(pd))]))
			trueA <- apply(trueA,c(1,2),as.numeric)
			rownames(trueA)<-gsub(true.A.token, "", colnames(pd)[grep(true.A.token, colnames(pd))])
			
			save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
			
		}
		
		if(!is.na(houseman.A.token)){
			
			trueA<-t(na.omit(pd[,grep(houseman.A.token, colnames(pd))]))
			trueA <- apply(trueA,c(1,2),as.numeric)
			rownames(trueA)<-gsub(houseman.A.token, "", colnames(pd)[grep(houseman.A.token, colnames(pd))])
			
			save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
			
		}else if(estimate.houseman.prop){
			
			if(!is.na(ref.rnb.ct.column)){
				rnb.set.ref<-remove.samples(rnb.set.ref, which(is.na(pheno(rnb.set.ref)[,ref.rnb.ct.column])))
				rnb.set<-combine(rnb.set, rnb.set.ref)
				ref.ct.column<-ref.rnb.ct.column
			}
			print("Estimating proportions using the Houseman et al, 2012 method")
			res<-estimateProportionsCP(rnb.set, ref.ct.column, NA, 2000, full.output = TRUE)
			
			trueA<-t(res$contributions.nonneg)
			trueA[trueA<1e-5]<-0
			
			trueA<-sweep(trueA, 2, colSums(trueA),"/")
			
			save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
		}
	}
	
	if(!is.na(ref.ct.column)){
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
	FILTER_QUALITY<- filter.beads || filter.intensity
	
	if(FILTER_QUALITY){
		M.raw<-RnBeads::M(rnb.set, row.names=TRUE)
		U.raw<-RnBeads::U(rnb.set, row.names=TRUE)
		b.raw<-RnBeads::covg(rnb.set, row.names=TRUE)
		
		if(!is.na(ref.ct.column)){
			
			M.raw<-M.raw[,subs,drop=FALSE]
			U.raw<-U.raw[,subs,drop=FALSE]
			b.raw<-b.raw[,subs, drop=FALSE]
			
		}
		
		saveRDS(M.raw, file.path(OUTPUTDIR, "Mint.RDS"))
		saveRDS(U.raw, file.path(OUTPUTDIR, "Uint.RDS"))
		saveRDS(b.raw, file.path(OUTPUTDIR, "Nbeads.RDS"))
		
		qual.filter<-filter.quality(rnb.set,beads=filter.beads, min.beads=min.n.beads, intensity = filter.intensity, 
		                            min.int.quant=min.int.quant, max.int.quant=max.int.quant,subs = subs)
		
		save(qual.filter, file=sprintf("%s/quality.filter.RData", OUTPUTDIR))
		
	}else{
		qual.filter<-1:nrow(rnb.set@meth.sites)
	}
	if(filter.na){
	  qual.filter <- filter.nas(rnb.set,subs=subs,qual.filter)
	}
	########################################## ANNOTATION FILTERING ###################################################
	FILTER_ANNOTATION<-filter.context || filter.snp || filter.sex.chromosomes
	
	if(FILTER_ANNOTATION){
		
		annot.filter<-filter.annotation(rnb.set, context = filter.context, snp = filter.snp, snp.list = snp.list,
		                                somatic = filter.sex.chromosomes, qual.filter = qual.filter, dist.snps = dist.snps)
	
		save(annot.filter, file=sprintf("%s/annotation.filter.RData", OUTPUTDIR))
	
		}else{
		annot.filter<-1:nrow(rnb.set@meth.sites)
	}
	
	total.filter<-intersect(qual.filter, annot.filter)
	logger.info(paste("Removing",nsites(rnb.set)-length(total.filter),"sites, retaining ",length(total.filter)))
	rnb.set.f<-remove.sites(rnb.set, setdiff(1:nrow(rnb.set@meth.sites), total.filter))
	
	if(filter.cross.reactive && inherits(rnb.set.f, "RnBeadSet")){
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
	                                           ifelse(filter.beads, "Beads", ""),
	                                           ifelse(filter.intensity, "Intensity", ""),
	                                           ifelse(filter.na, "Missing", ""),
	                                           ifelse(filter.context, "Context", ""),
	                                           ifelse(filter.snp, "SNP", ""),
	                                           ifelse(filter.sex.chromosomes, "Somatic", "")
	)
	
	analysis_info$NORMALIZATION <- normalization
	
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
#' @noRd
#' @author Michael Scherer, Pavlo Lutsik
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
    b.raw<-RnBeads::covg(rnb.set, row.names=TRUE)
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
#' @noRd
#' @author Michael Scherer, Pavlo Lutsik
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
#' @noRd
#' @author Pavlo Lutsik, Michael Scherer
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
#' @noRd
#' @author Pavlo Lutsik, Michael Scherer
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
#' @noRd
#' @author Michael Scherer, Pavlo Lutsik
#' @import LaplacesDemon
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
      if(!requireNamespace('LaplacesDemon')){
        stop("Missing required package 'LaplacesDemon'. Please install it.")
      }
      meth.data <- meth(rnb.set)
      logger.start("SNP computation")
      pair.dist <- apply(meth.data,1,function(snp){
        this.cpg <- lapply(snp,function(snp2){
          abs(snp-snp2)
        })
        is.trimodal(unlist(this.cpg))
      })
      logger.info(paste("Filtered",length(setdiff(snp.filter,which(!pair.dist))),"CpGs in distribution SNP filtering"))
      snp.filter <- intersect(snp.filter,which(!pair.dist))
      rm(meth.data)
      logger.completed()
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
#' @noRd
#' @author Michael Scherer, Pavlo Lutsik
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

#' prepare.data.BS
#' 
#' This functions prepares sequencing data sets for a MeDeCom run.
#' 
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBiseqSet-class}} for which analysis is to be performed.
#' @param work.dir A path to a existing directory, in which the results are to be stored
#' @param analysis.name A string representing the dataset for which analysis is to be performed. Only used to create a folder with a 
#'                 descriptive name of the analysis.
#' @param sample.selection.col A column name in the phenotypic table of \code{rnb.set} used to selected a subset of samples for
#'                 analysis that contain the string given in \code{sample.selection.col}.
#' @param sample.selection.grep A string used for selecting samples in the column \code{sample.selection.grep}.
#' @param ref.ct.column Column name in \code{rnb.set} used to extract methylation information on the reference cell types.
#' @param pheno.columns Vector of column names in the phenotypic table of \code{rnb.set} that is kept and exported for further 
#'                 exploration.
#' @param prepare.true.proportions Flag indicating if true proportions are either available in \code{rnb.set} or to be estimated 
#'                          with Houseman's reference-based deconvolution approach.
#' @param true.A.token String present in the column names of \code{rnb.set} used for selecting the true proportions of the corresponding
#'                      cell types.
#' @param houseman.A.token Similar to \code{true.A.token}, but not containing the true proportions, rather the estimated proportions
#'                      by Houseman's method.
#' @param id.column Sample-specific ID column name in \code{rnb.set}
#' @param filter.coverage Flag indicating, if site-filtering based on coverage is to be conducted.
#' @param min.coverage Minimum number of reads required in each sample for the site to be considered for adding to MeDeCom.
#' @param min.covg.quant Lower quantile of coverages. Values lower than this value will be ignored for analysis.
#' @param max.covg.quant Upper quantile of coverages. Values higher than this value will be ignored for analysis.
#' @param filter.na Flag indicating if sites with any missing values are to be removed or not.
#' @param filter.snp Flag indicating if annotated SNPs are to be removed from the list of sites according to RnBeads' SNP list.
#' @param snp.list Path to a file containing positions of known SNPs to be removed from the analysis, if \code{filter.snp} is \code{TRUE}. The coordinates must be the 
#'             provided in the same genome assembly as rnb.set. The file must be a tab-separated value (tsv) file with only one header line an the following meaning of 
#'             the rows: 1st row: chromosome, 2nd row: position of the SNP on the chromosome
#' @param filter.sex.chromosomes Flag indicating if only somatic probes are to be kept.
#' @param execute.lump Flag indicating if the LUMP algorithm is to be used for estimating the amount of immune cells in a particular sample.
#' @return A list with four elements: \itemize{
#'           \item quality.filter The indices of the sites that survived quality filtering
#' }
#' @export
#' @author Michael Scherer, Pavlo Lutsik
#' @import RnBeads
#' @import R.utils
prepare.data.BS <- function(
		rnb.set, 
		work.dir=getwd(),
		analysis.name="analysis",
		sample.selection.col=NA,
		sample.selection.grep=NA,
		ref.ct.column=NA,
		pheno.columns=NA,
		prepare.true.proportions=FALSE,
		true.A.token=NA,
		houseman.A.token=NA,
		id.column=rnb.getOption("identifiers.column"),
		filter.coverage = hasCovg(rnb.set),
		min.coverage=5,
		min.covg.quant=0.05,
		max.covg.quant=0.95,
		filter.na=TRUE,
		filter.snp=TRUE,
		snp.list=NULL,
		filter.sex.chromosomes=TRUE,
		execute.lump=FALSE
){

	OUTPUTDIR <- file.path(work.dir, analysis.name)
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
	if(is.character(rnb.set)){
		rnb.set<-load.rnb.set(rnb.set)
	}else if(inherits(rnb.set,"RnBSet")){
		rnb.set<-rnb.set
	}	
	if(!is.na(sample.selection.col)){
		
		SAMP_SUBS<-grep(sample.selection.grep, pheno(rnb.set)[[sample.selection.col]])
		
		rms<-setdiff(1:length(samples(rnb.set)), SAMP_SUBS)
		
		rnb.set<-remove.samples(rnb.set, rms)
	}
		
	meth.rnb <- rnb.set@meth.sites
	pd<-pheno(rnb.set)
	if(!is.na(ref.ct.column)){
	  subs<-is.na(pd[[ref.ct.column]])
	}else{
	  subs<-1:nrow(pd)
	}
	if(!is.na(pheno.columns)){
		pheno.data<-pd[,pheno.columns,drop=FALSE]
		save(pheno.data, file=sprintf("%s/pheno.RData", OUTPUTDIR))
	}
	
	if(!is.null(id.column)){
		sample_ids<-pd[,id.column]
		saveRDS(sample_ids, file=sprintf("%s/sample_ids.RDS", OUTPUTDIR))	
	}
	if(prepare.true.proportions){
	  if(!is.na(true.A.token)){
	    
	    trueA<-t(na.omit(pd[subs,grep(true.A.token, colnames(pd))]))
	    trueA <- apply(trueA,c(1,2),as.numeric)
	    rownames(trueA)<-gsub(true.A.token, "", colnames(pd)[grep(true.A.token, colnames(pd))])
	    
	    save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
	    
	  }
	  
	  if(!is.na(houseman.A.token)){
	    
	    trueA<-t(na.omit(pd[,grep(houseman.A.token, colnames(pd))]))
	    trueA <- apply(trueA,c(1,2),as.numeric)
	    rownames(trueA)<-gsub(houseman.A.token, "", colnames(pd)[grep(houseman.A.token, colnames(pd))])
	    
	    save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
	    
	  }
	}
	if(!is.na(ref.ct.column)){
	  ct<-pd[[ref.ct.column]]
	  nnas<-!is.na(ct)
	  ct<-ct[nnas]
	  
	  meth.ref.ct<-lapply(unique(ct), function(cn) rowMeans(meth.rnb[,nnas][,ct==cn]))
	  trueT<-do.call("cbind", meth.ref.ct)
	  colnames(trueT)<-unique(ct)
	  
	  save(trueT, file=sprintf("%s/trueT.RData", OUTPUTDIR))
	}
	if(!is.na(ref.ct.column)){
	  meth.data<-meth.rnb[,subs]
	}else{
	  meth.data<-meth.rnb
	}
	colnames(meth.data) <- NULL	
	save(meth.data, file=sprintf("%s/data.set.RData", OUTPUTDIR))
	if(filter.coverage){
		qual.filter <- filter.quality.covg(rnb.set, min.covg = min.coverage,  
		                            min.covg.quant = min.covg.quant, max.covg.quant=max.covg.quant)		
		save(qual.filter, file=sprintf("%s/quality.filter.RData", OUTPUTDIR))		
	}else{
		qual.filter<-1:nrow(rnb.set@meth.sites)
	}
	if(filter.na){
	  qual.filter <- filter.nas.biseq(rnb.set,subs=1:length(samples(rnb.set)),qual.filter)
	}
	FILTER_ANNOTATION <- filter.snp || filter.sex.chromosomes
	
	if(FILTER_ANNOTATION){
		
		annot.filter <- filter.annotation.biseq(rnb.set, snp = filter.snp, snp.list = snp.list,
		                                somatic = filter.sex.chromosomes, qual.filter = qual.filter)
	
		save(annot.filter, file=sprintf("%s/annotation.filter.RData", OUTPUTDIR))
	
		}else{
		annot.filter<-1:nrow(rnb.set@meth.sites)
	}
	if(execute.lump){
	  lump.est <- rnb.execute.lump(rnb.set)
	  rnb.set <- addPheno(rnb.set,as.vector(lump.est),"LUMP_estimate")
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
