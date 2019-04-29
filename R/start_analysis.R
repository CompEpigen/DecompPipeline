#' start.analysis
#' 
#' Wrapper function to start one of the deconvolution algorithms \code{MeDeCom}, \code{RefFreeEWAS} or \code{EDec}
#' 
#' @param meth.data A \code{matrix} or \code{data.frame} containing methylation information. If NULL, methylation information needs to be provided
#'                   through \code{rnb.set}
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation and sample meta information.
#' @param work.dir Working directory for the analysis.
#' @param cg_groups List of CpG indices used for the analysis. Can be computed by \code{\link{prepare_CG_subsets}}.
#' @param Ks Vector of integers used as components in MeDeCom.
#' @param factorviz.outputs Flag indicating, if outputs should be stored to be compatible with FactorViz for data exploration
#' @param method The method to be used for deconvolution. Can be one of \code{MeDeCom}, \code{RefFreeCellMix} or \code{EDec}.
#' @author Michael Scherer
#' @export
start.analysis <- function(meth.data=NULL,
                           rnb.set=NULL,
                           cg_groups,
                           Ks,
                           work.dir,
                           factorviz.outputs=FALSE,
                           method="MeDeCom",
                           ...){
  all.methods <- c("MeDeCom","RefFreeCellMix","EDec")
  if(!method %in% all.methods){
    stop(paste0("Invalid value for method. Needs to be one of ",all.methods))
  }
  if(method == "MeDeCom"){
    md.res <- start_medecom_analysis(meth.data = meth.data,
                                     rnb.set = rnb.set,
                                     cg_groups = cg_groups,
                                     Ks = Ks,
                                     WORK_DIR = work.dir,
                                     factorviz.outputs = factorviz.outputs,
                                     ...)
  }else if(method == "RefFreeCellMix"){
    md.res <- start.refreeewas.analysis(meth.data=meth.data,
                                                    rnb.set=rnb.set,
                                                    cg_groups=cg_groups,
                                                    Ks=Ks,
                                                    work.dir=work.dir,
                                                    factorviz.outputs=factorviz.outputs)
  }else if(method == "EDec"){
    md.res <- start.edec.analysis(meth.data=meth.data,
                                  rnb.set=rnb.set,
                                  cg_groups=cg_groups,
                                  Ks=Ks,
                                  work.dir = work.dir,
                                  factorviz.outputs = factorviz.outputs)
  }
  return(md.res)
}

#' start.edec.analysis
#' 
#' This function executes EDec for the specified CpGs and the number of cell types K.
#' 
#' @param meth.data A \code{matrix} or \code{data.frame} containing methylation information. If NULL, methylation information needs to be provided
#'                   through \code{rnb.set}
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation and sample meta information.
#' @param cg_groups List of CpG indices used for the analysis. Can be computed by \code{\link{prepare_CG_subsets}}.
#' @param Ks The number of cell types to be tested. Can be a single numeric value or an array of numbers.
#' @param work.dir The working directory to be used.
#' @param factorviz.outputs Flag indicating, if outputs should be stored to be compatible with FactorViz for data exploration
#' @return An object of type \code{\link{MeDeComSet}} containing the results of the EDec experiment.
#' @author Michael Scherer
#' @export
start.edec.analysis <- function(meth.data=NULL,
                                rnb.set=NULL,
                                cg_groups,
                                Ks,
                                work.dir=getwd(),
                                factorviz.outputs=FALSE){
  if(is.null(meth.data) && is.null(rnb.set)){
    logger.error("No input methylation data provided")
  }
  if(is.null(meth.data)){
    if(inherits(rnb.set,"RnBSet")){
      meth.data <- meth(rnb.set,row.names=T)
    }else{
      logger.error("Invalid value for rnb.set")
    }
  }else if(!(grepl("cg",row.names(meth.data)))){
    stop("Rownames of methylation data need to be provided for EDec")
  }
  require("EDec")
  T.all <- list()
  A.all <- list()
  rss.all <- list()
  res.all <- list()
  for(i.group in 1:length(cg_groups)){
    logger.start(paste("Processing group:",i.group))
    group <- cg_groups[[i.group]]
    group <- row.names(meth.data)[group]
    rss.vec <- c()
    T.all <- list()
    A.all <- list()
    for(j.K in 1:length(Ks)){
      K <- Ks[j.K]
      logger.start(paste("Processing K:",K))
      edec.res <- run_edec_stage_1(meth.data,
                                   informative_loci = group,
                                   num_cell_types = K)
      rss.vec <- c(rss.vec,edec.res$res.sum.squares)
      T.all[[j.K]] <- edec.res$methylation
      A.all[[j.K]] <- t(edec.res$proportions)
      logger.completed()
    }
    rss.all[[i.group]] <- rss.vec
    res.all[[i.group]] <- list(T=T.all,A=A.all)
    logger.completed()
  }
  result <- as.MeDeComSet(res.all,cg_subsets=1:length(cg_groups),Ks=Ks,rss=rss.all,m.orig=nrow(meth.data),n.orig=ncol(meth.data))
  result@parameters$GROUP_LISTS <- cg_groups
  if(factorviz.outputs){
    store.path <- file.path(work.dir,"FactorViz_outputs")
    if(!file.exists(store.path)){
      dir.create(store.path)
    }
    if(!is.null(rnb.set)){
      result@parameters$ASSEMBLY <- assembly(rnb.set)
      ann.C <- annotation(rnb.set)
      ann.S <- pheno(rnb.set)
      save(ann.C,file=file.path(store.path,"ann_C.RData"))
      save(ann.S,file=file.path(store.path,"ann_S.RData"))
    }
    medecom.set <- result
    save(medecom.set,file=file.path(store.path,"medecom_set.RData"))
    save(meth.data,file=file.path(store.path,"meth_data.RData"))
  }
  return(result)
}

#' start.refreeewas.analysis
#' 
#' This function executes RefFreeCellMix for the specified CpGs and the number of cell types K.
#' 
#' @param meth.data A \code{matrix} or \code{data.frame} containing methylation information. If NULL, methylation information needs to be provided
#'                   through \code{rnb.set}
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation and sample meta information.
#' @param cg_groups List of CpG indices used for the analysis. Can be computed by \code{\link{prepare_CG_subsets}}.
#' @param Ks The number of cell types to be tested. Can be a single numeric value or an array of numbers.
#' @param work.dir The working directory to be used.
#' @param factorviz.outputs Flag indicating, if outputs should be stored to be compatible with FactorViz for data exploration
#' @return An object of type \code{\link{MeDeComSet}} containing the results of the RefFreeCellMix experiment.
#' @author Michael Scherer
#' @export
start.refreeewas.analysis <- function(meth.data=NULL,
                                      rnb.set=NULL,
                                      cg_groups,
                                      Ks,
                                      work.dir=getwd(),
                                      factorviz.outputs=FALSE){
  if(is.null(meth.data) && is.null(rnb.set)){
    logger.error("No input methylation data provided")
  }
  if(is.null(meth.data)){
    if(inherits(rnb.set,"RnBSet")){
      meth.data <- meth(rnb.set)
    }else{
      logger.error("Invalid value for rnb.set")
    }
  }
  res.all <- list()
  devis.all <- list()
  for(i.group in 1:length(cg_groups)){
    logger.start(paste("Processing group:",i.group))
    group <- cg_groups[[i.group]]
    meth.sset <- meth.data[group,]
    res.sset <- RefFreeCellMixArray(meth.sset,Klist=Ks)
    devis <- tryCatch(RefFreeCellMixArrayDevianceBoots(res.sset,Y=meth.sset),error=function(e)e)
    if(inherits(devis,"error")){
      devis <- rep(NA,length(Ks))
    }else{
      devis <- colMeans(devis)
    }
    devis.all[[i.group]] <- devis
    res.all[[i.group]] <- res.sset
    logger.completed()
  }
  result <- as.MeDeComSet(res.all,cg_subsets=1:length(cg_groups),Ks=Ks,deviances=devis.all,m.orig=nrow(meth.data),n.orig=ncol(meth.data))
  result@parameters$GROUP_LISTS <- cg_groups
  if(factorviz.outputs){
    store.path <- file.path(work.dir,"FactorViz_outputs")
    if(!file.exists(store.path)){
      dir.create(store.path)
    }
    if(!is.null(rnb.set)){
      result@parameters$ASSEMBLY <- assembly(rnb.set)
      ann.C <- annotation(rnb.set)
      ann.S <- pheno(rnb.set)
      save(ann.C,file=file.path(store.path,"ann_C.RData"))
      save(ann.S,file=file.path(store.path,"ann_S.RData"))
    }
    medecom.set <- result
    save(medecom.set,file=file.path(store.path,"medecom_set.RData"))
    save(meth.data,file=file.path(store.path,"meth_data.RData"))
  }
  return(result)
}

#' start_medecom_analysis
#' 
#' Wrapper for runMeDeCom, for data preprocessed through the DecompPipeline
#' 
#' @param meth.data A \code{matrix} or \code{data.frame} containing methylation information. If NULL, methylation information needs to be provided
#'                   through \code{rnb.set}
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation and sample meta information.
#' @param WORK_DIR Working directory for the analysis.
#' @param cg_groups List of CpG indices used for the analysis. Can be computed by \code{\link{prepare_CG_subsets}}.
#' @param Ks Vector of integers used as components in MeDeCom.
#' @param LAMBDA_GRID Vector of doubles representing the regularization parameter in MeDeCom.
#' @param SAMPLE_SUBSET Vector of indices of samples to be included in the analysis. If \code{NULL}, all samples are included.
#' @param K_FIXED Columns in the T matrix that should be fixed. If \code{NULL}, no columns are fixed.
#' @param WRITE_FILES Flag indicating if intermediate results are to be stored.
#' @param factorviz.outputs Flag indicating, if outputs should be stored to be compatible with FactorViz for data exploration
#' @param opt.method Optimization method to be used. Either MeDeCom.quadPen or MeDeCom.cppTAfact (default).
#' @param startT Inital matrix for T.
#' @param startA Initial matrix for A.
#' @param trueT True value for the T matrix.
#' @param trueA True value for the A matrix.
#' @param analysis.name Name of the analysis.
#' @param folds Integer representing the number of folds used in the analysis.
#' @param cores Integer representing the number of cores to be used in the analysis.
#' @param itermax Maximum number of iterations
#' @param ninit Number if initialtions.
#' @param CLUSTER_SUBMIT Flag indicating, if the jobs are to be submitted to a scientific compute cluster (only SGE supported).
#' @param CLUSTER_RDIR Path to an executable version of R.
#' @param CLUSTER_HOSTLIST Regular expression, on which basis hosts are selected in the cluster environment.
#' @param CLUSTER_MEMLIMIT the \code{memlimit} resource value of the cluster submission.
#' @param CLEANUP Flag indicating if temprary files are to be deleted.
#' @param analysis_info Information to be saved about the analysis. Just stored as info.
#' @param LAMBDA_GRID_TYPE String represent the lambda grid that was chosen. Just stored as info.
#' @param ANALYSIS_TOKEN String specifying the type of analysis that was conducted. Just stored as info.
#' @return An object of type \code{\link{MeDeComSet}} containing the results of the MeDeCom experiment.
#' @export
start_medecom_analysis<-function(
    meth.data=NULL,
		rnb.set=NULL,
		WORK_DIR=getwd(),
		cg_groups,
		Ks,
		LAMBDA_GRID,
		SAMPLE_SUBSET=NULL,
		K_FIXED=NULL,
		WRITE_FILES=TRUE,
		factorviz.outputs=F,
		opt.method = "MeDeCom.cppTAfact",
		startT=NULL,
		startA=NULL,
		trueT = NULL,
		trueA = NULL,
		analysis.name="MeDeComRun",
		folds=10,
		cores=1,
		itermax=1000,
		ninit=100,
		CLUSTER_SUBMIT=FALSE,
		CLUSTER_RDIR=NA,
		CLUSTER_HOSTLIST="*",
		CLUSTER_MEMLIMIT="5G",
#		MAX_JOBS=1000,
#		WAIT_TIME="30m",
#		PORTIONS=FALSE,
#		JOB_FILE=NA,
		CLEANUP=FALSE,
		analysis_info=NULL,
		LAMBDA_GRID_TYPE="standard",
		ANALYSIS_TOKEN="customAnalysis"
){
	library(MeDeCom)
	library(R.utils)
	
	#RDIR="/TL/deep-share/archive00/software/bin"
	#.libPaths(sprintf("%s/Rlib_test", WORK_DIR))
  
  if(!is.null(analysis_info)){
    ANALYSIS_ID<-paste(
      analysis_info$DATASET, 
      analysis_info$DATA_SUBSET,
      analysis_info$NORMALIZATION,
      analysis_info$QUALITY_FILTERING, 
      analysis_info$MARKER_SELECTION,
      LAMBDA_GRID_TYPE,
      ANALYSIS_TOKEN,
      sep="_")
  }else{
    ANALYSIS_ID<-"customAnalysis"
    analysis_info<-list()
  }
  analysis_info$ANALYSIS<-ANALYSIS_ID
  
  WORK_DIR <- file.path(WORK_DIR,analysis.name)
  if(!file.exists(WORK_DIR)){
    dir.create(WORK_DIR)
  }
  log.file <- file.path(WORK_DIR,"analysis.log")
  if(!file.exists(log.file)){
    if(logger.isinitialized()){
      logger.close()
      logger.start(fname=log.file)
    }else{
      logger.start(fname=log.file)
    }
  }
	
	if(is.na(CLUSTER_RDIR)){
		if(any(grepl("deep", R.utils:::System$getHostname()))){
			#RDIR="/TL/deep-share/archive00/software/packages/R/R-3.1.1/bin"
			CLUSTER_RDIR="/TL/deep-share/archive00/software/packages/R/R-devel_20160307/bin"
		}else if(any(grepl("t7600", R.utils:::System$getHostname()))){
			CLUSTER_RDIR="/usr/bin/"
		}else{
			CLUSTER_RDIR="/usr/bin/"
		}
	}
	
#	SRCDIR=sprintf("%s/projects/parameter_tuning/src", WORK_DIR)
#	DATA_DIR=sprintf("%s/projects/parameter_tuning/data", WORK_DIR)
#	GLOBAL_DD<-sprintf("%s/projects/parameter_tuning/data/common", WORK_DIR)
	
	if(TRUE){
		print("Did not write the variable dump: should only be executed from an environment with all the variables set")	
	}else{
		var_list<-c(ANALYSIS)
		system(sprintf("mkdir %s", WORK_DIR))
		dump(ls()[var_list], file=file.path(WORK_DIR, "analysis_settings.RDump"))
	}
	
  if(is.null(meth.data)){
  	if(is.null(rnb.set)){
  		load(sprintf("%s/data.set.RData", WORK_DIR))
  	}else{
  		meth.data<-meth(rnb.set)
  	}
  }
	
#	if(!DATASET %in% c("sim", "simReal")){
#		load(sprintf("%s/indices.RData", GLOBAL_DD))
#		load(sprintf("%s/snp.filter.RData", GLOBAL_DD))
#		load(sprintf("%s/somatic.filter.RData", GLOBAL_DD))
#	}
	
#	if(file.exists(sprintf("%s/quality.filter.RData", DD))){
#		load(sprintf("%s/quality.filter.RData", DD))
#	}
	
#	if(file.exists(sprintf("%s/Mint.RDS", DD))){
#		M.raw<-readRDS(sprintf("%s/Mint.RDS", DD))
#	}else{
#		M.raw<-NULL
#	}
	
#	if(file.exists(sprintf("%s/Uint.RDS", DD))){
#		U.raw<-readRDS(sprintf("%s/Uint.RDS", DD))
#	}else{
#		U.raw<-NULL
#	}
	
#	if(file.exists(sprintf("%s/Nbeads.RDS", DD))){
#		b.raw<-readRDS(sprintf("%s/Nbeads.RDS", DD))
#	}else{
#		b.raw<-NULL
#	}
	
	if(WRITE_FILES){
		saveRDS(LAMBDA_GRID, file=sprintf("%s/lambda_grid.RDS", WORK_DIR))
		if(!is.null(SAMPLE_SUBSET)){
			saveRDS(SAMPLE_SUBSET, file=sprintf("%s/sample_subset.RDS", WORK_DIR))
		}else{
			SAMPLE_SUBSET<-1:ncol(meth.data)
		}
	}
	
	groups<-1:length(cg_groups)
	
#	if(QUALITY_FILTERING %in% c("standard","MplusU") && !is.null(M.raw) && !is.null(U.raw)){
#		
#		if(!"MIN_INT_QUANT" %in% ls()){
#			MIN_INT_QUANT<-0.1
#		}
#		if(!"MAX_INT_QUANT" %in% ls()){
#			MAX_INT_QUANT<-0.95
#		}
#		if(!"MIN_N_BEADS" %in% ls()){
#			MIN_N_BEADS<-3
#		}
#		
#		MplusU<-M.raw+U.raw
#		
#		hm450_ann<-readRDS(file.path(GLOBAL_DD, "hm450_cg_annot.RDS"))
#		
#		MplusU.I<-MplusU[hm450_ann$Design=="I",]
#		MplusU.II<-MplusU[hm450_ann$Design=="II",]
#		
#		MU.q001.I<-sort(as.numeric(MplusU.I))[ceiling(MIN_INT_QUANT*nrow(MplusU.I)*ncol(MplusU.I))]
#		MU.q099.I<-sort(as.numeric(MplusU.I))[ceiling(MAX_INT_QUANT*nrow(MplusU.I)*ncol(MplusU.I))]
#		
#		MU.q001.II<-sort(as.numeric(MplusU.II))[ceiling(MIN_INT_QUANT*nrow(MplusU.II)*ncol(MplusU.II))]
#		MU.q099.II<-sort(as.numeric(MplusU.II))[ceiling(MAX_INT_QUANT*nrow(MplusU.II)*ncol(MplusU.II))]
#		
#		MplusU.f<-matrix(FALSE, nrow=nrow(MplusU), ncol=ncol(MplusU))
#		
#		MplusU.f[hm450_ann$Design=="I",]<-MplusU.I>MU.q001.I & MplusU.I<MU.q099.I
#		MplusU.f[hm450_ann$Design=="II",]<-MplusU.II>MU.q001.II & MplusU.II<MU.q099.II
#		
#		qf.MU<-which(rowSums(MplusU.f)==ncol(MplusU.f))
#		
#		if(!is.null(b.raw)){
#			qf.b<-which(rowSums(b.raw>=MIN_N_BEADS)==ncol(b.raw))
#			qf.MU.beads<-intersect(qf.MU, qf.b)
#		}else{
#			qf.MU.beads<-qf.MU
#		}
#		
#	}
	
#	for(group in groups){
#		
#		# LOAD THE RESULTS
#		#if(DATASET %in% c("sim", "simReal")){
#		ind<-1:nrow(meth.data)	
#		#}else{
#		#	ind<-sort(unlist(indices[GROUP_LISTS[[group]]]))
#		#}
#		
#		if(QUALITY_FILTERING=="standard"){
#			ind<-intersect(ind, qf.MU.beads)
#			ind<-intersect(ind, snp.filter)
#			ind<-intersect(ind, somatic.filter)
#		}else if(QUALITY_FILTERING=="snpANDsomatic"){
#			ind<-intersect(ind, snp.filter)
#			ind<-intersect(ind, somatic.filter)
#		}else if(QUALITY_FILTERING=="somatic"){
#			ind<-intersect(ind, somatic.filter)
#		}else if(QUALITY_FILTERING=="snp"){
#			ind<-intersect(ind, snp.filter)
#		}else if(QUALITY_FILTERING=="MplusU"){
#			ind<-intersect(ind, qf.MU.beads)
#		}
#		
#	}
	if(!is.null(K_FIXED)){
		saveRDS(K_FIXED, file=sprintf("%s/fixed_T_cols.RDS", WORK_DIR))
	}else{
		K_FIXED<-NULL
	}
	
#	if("START" %in% ls() && file.exists(START)){
#		system(sprintf("cp %s %s/start.RData", START, WORK_DIR))
#		load(file.path(WORK_DIR, START))
#		startT=result$T
#		startA=result$A
#	}else{
#		startT=NULL
#		startA=NULL
#	}
	
	if(!is.null(trueT)){
	  if(is.character(trueT)){
	    if(!file.exists(trueT)){
	      logger.error(paste("File for trueT",trueT,"does not exist."))
	    }else{
		    load(trueT)
	    }
	  }else if(!is.matrix(trueT)){
	    logger.error("Invalid value for trueT")
	  }
	}
	
	if(!is.null(trueA)){
	  if(is.character(trueA)){
	    if(!file.exists(trueA)){
	      logger.error(paste("File for trueA",trueA,"does not exist."))
	    }else{
	      load(trueA)
	    }
	  }else if(!is.matrix(trueA)){
	    logger.error("Invalid value for trueA")
	  }
	}
#	if(PORTIONS && is.na(JOB_FILE)){
#		JOB_FILE<-"/tmp/job_file"
#	}
#	cluster_submit<-function(qsub_string, portions=PORTIONS, job_script_file=JOB_FILE){
#		
#		if(portions){
#			cat(qsub_string, file=JOB_FILE, append=TRUE, sep="\n")
#		}else{
#			system(qsub_string)
#			#print(qsub_string)
#		}
#		
#	}
	
	if(CLUSTER_SUBMIT){
		cluster.settings=list(R_bin_dir=CLUSTER_RDIR, host_pattern=CLUSTER_HOSTLIST, mem_limit=CLUSTER_MEMLIMIT)
	}else{
		cluster.settings=NULL
	}
	
	result<-runMeDeCom(data=meth.data, 
			Ks=Ks,
			lambdas=LAMBDA_GRID,
			opt.method=opt.method,
			cg_subsets=cg_groups,
			sample_subset=SAMPLE_SUBSET,
			startT=startT,
			startA=startA,
			trueT=trueT,
			trueA=trueA,
			fixed_T_cols=K_FIXED,
			NINIT=ninit, 
			ITERMAX=itermax, 
			NFOLDS=folds,
			N_COMP_LAMBDA=4,
			NCORES=cores,
			analysis.name=analysis.name,
			use.ff=FALSE,
			cluster.settings=cluster.settings,
			temp.dir=WORK_DIR,
			cleanup=CLEANUP,
			verbosity=1L,
			time.stamps=TRUE
			#,random.seed=RANDOM_SEED
	)
	
	### TODO: Fix this once
	
	result@parameters$ANALYSIS <- analysis.name
	result@parameters$GROUP_LISTS <- cg_groups
	# result@parameters$cg_subsets <- c(1:length(cg_groups))
	# result@parameters$ANALYSIS <- ANALYSIS_ID
	# result@parameters$NORMALIZATION <- NORMALIZATION
	result@parameters$ITERMAX<-itermax
	#    result@parameters$MARKER_SELECTION<- MARKER_SELECTION
	result@parameters$NFOLDS<-folds
	#    result@parameters$ANALYSIS_TOKEN<-""
	result@parameters$NINIT<-ninit
	#    result@parameters$DATASET<-DATASET 
	#    result@parameters$DATA_SUBSET<-DATA_SUBSET 
	
	analysis_info$GROUP_LISTS<-cg_groups
	analysis_info$ANALYSIS_DATE<-date()
	result@dataset_info<-c(result@dataset_info, analysis_info)
	
	if(WRITE_FILES){
		saveRDS(result, file=file.path(WORK_DIR, "collected.result.RDS"))
	}
	
	if(factorviz.outputs){
	  store.path <- file.path(WORK_DIR,"FactorViz_outputs")
	  if(!file.exists(store.path)){
	    dir.create(store.path)
	  }
	  if(!is.null(rnb.set)){
	    result@parameters$ASSEMBLY <- assembly(rnb.set)
	    ann.C <- annotation(rnb.set)
	    ann.S <- pheno(rnb.set)
	    save(ann.C,file=file.path(store.path,"ann_C.RData"))
	    save(ann.S,file=file.path(store.path,"ann_S.RData"))
	  }
	  medecom.set <- result
	  save(medecom.set,file=file.path(store.path,"medecom_set.RData"))
	  save(meth.data,file=file.path(store.path,"meth_data.RData"))
	  if(!is.null(trueT)){
	    ref.meth <- trueT
	    save(ref.meth,file=file.path(store.path,"ref_meth.RData"))
	  }
	  if(!is.null(trueA)){
	    ref.props <- trueA
	    save(ref.props,file=file.path(store.path,"ref_props.RData"))
	  }
	}
	return(result)
}

#' start_decomp_pipeline
#' 
#' CPG FILTERING (BeadChip)
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} for which analysis is to be performed.
#' @param Ks Vector of integers used as components in MeDeCom.
#' @param lambda.grid Vector of doubles representing the regularization parameter in MeDeCom.
#' @param work.dir A path to a existing directory, in which the results are to be stored
#' @param analysis.name A string representing the dataset for which analysis is to be performed. Only used to create a folder with a 
#'                 descriptive name of the analysis.
#' @param sample.selection.col A column name in the phenotypic table of \code{RNB_SET} used to selected a subset of samples for
#'                 analysis that contain the string given in \code{SAMPLE_SELECTION_GREP}.
#' @param sample.selection.grep A string used for selecting samples in the column \code{SAMPLE_SELECTION_COL}.
#' @param pheno.cols Vector of column names in the phenotypic table of \code{RNB_SET} that is kept and exported for further 
#'                 exploration.
#' @param id.column Sample-specific ID column name in \code{RNB_SET}
#' @param normalization Normalization method to be performed before employing MeDeCom. Can be one of \code{"none","dasen","illumina","noob"} (BeadChip only).
#' @param ref.ct.column Column name in \code{RNB_SET} used to extract methylation information on the reference cell types.
#' @param ref.rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} containing methylation information on reference cell types (BeadChip only).
#' @param ref.rnb.ct.column Column name in \code{REF_RNB_SET} used to extract methylation information on the reference cell types (BeadChip only).
#' @param prepare.true.proportions Flag indicating if true proportions are either available in \code{RNB_SET} or to be estimated 
#'                          with Houseman's reference-based deconvolution approach (BeadChip only).
#' @param true.A.token String present in the column names of \code{RNB_SET} used for selecting the true proportions of the corresponding
#'                      cell types.
#' @param houseman.A.token Similar to \code{TRUE_A_TOKEN}, but not containing the true proportions, rather the estimated proportions
#'                      by Houseman's method (BeadChip only).
#' @param estimate.houseman.prop If neither \code{TRUE_A_TOKEN} nor \code{HOUSEMAN_A_TOKEN} are given, the proportions of the reference
#'                      cell type are estimated with Houseman's approach (BeadChip only).
#' @param filter.beads Flag indicating, if site-filtering based on the number of beads available is to be conducted (BeadChip only).
#' @param min.n.beads Minimum number of beads required in each sample for the site to be considered for adding to MeDeCom (BeadChip only).
#' @param filter.intensity  Flag indicating if sites should be removed according to the signal intensities (the lowest and highest quantiles
#'                      given by \code{MIN_INT_QUANT} and \code{MAX_INT_QUANT}) (BeadChip only).
#' @param min.int.quant Lower quantile of intensities which is to be removed (BeadChip only).
#' @param max.int.quant Upper quantile of intensities which is to be removed (BeadChip only).
#' @param filter.na Flag indicating if sites with any missing values are to be removed or not.
#' @param filter.context Flag indicating if only CG probes are to be kept (BeadChip only).
#' @param filter.cross.reactive Flag indicating if sites showing cross reactivity on the array are to be removed.
#' @param execute.lump Flag indicating if the LUMP algorithm is to be used for estimating the amount of immune cells in a particular sample.
#' @param remove.ICA Flag indicating if independent component analysis is to be executed to remove potential confounding factor.
#'             If \code{TRUE},conf.fact.ICA needs to be specified.
#' @param conf.fact.ICA Column name in the sample annotation sheet representing a potential confounding factor.
#' @param ica.setting Optional argument setting up ICA.
#' @param filter.snp Flag indicating if annotated SNPs are to be removed from the list of sites according to RnBeads' SNP list. (@TODO: we
#'                     could provide an addititional list of SNPs, similar to RnBeads blacklist for filtering)
#' @param snp.list Path to a file containing CpG IDs of known SNPs to be removed from the analysis, if \code{FILTER_SNP} is \code{TRUE}.
#' @param filter.somatic Flag indicating if only somatic probes are to be kept.
#' CPG FILTERING (BS)
#' @param filter.coverage Flag indicating, if site-filtering based on coverage is to be conducted (BS only).
#' @param min.coverage Minimum number of reads required in each sample for the site to be considered for adding to MeDeCom (BS only).
#' @param min.covg.quant Lower quantile of coverages. Values lower than this value will be ignored for analysis (BS only).
#' @param max.covg.quant Upper quantile of coverages. Values higher than this value will be ignored for analysis (BS only).
#' CG_SUBSET SELECTION
#' @param marker.selection A vector of strings representing marker selection methods. Available method are \itemize{
#'                                  \item{"\code{all}"} Using all sites available in the input.
#'                                  \item{"\code{pheno}"} Selected are the top \code{N_MARKERS} site that differ between the phenotypic
#'                                         groups defined in data preparation or by \code{\link{rnb.sample.groups}}. Those are
#'                                         selected by employing limma on the methylation matrix.
#'                                  \item{"\code{houseman2012}"} The 50k sites reported as cell-type specific in the Houseman's reference-
#'                                         based deconvolution. See Houseman et.al. 2012.
#'                                  \item{"\code{houseman2014}"} Selects the sites said to be linked to cell type composition by \code{RefFreeEWAS},
#'                                         which is similar to surrogate variable analysis. See Houseman et.al. 2014.
#'                                  \item{"\code{jaffe2014}"} The sites stated as related to cell-type composition Jaffe et.al. 2014.
#'                                  \item{"\code{rowFstat}"} Markers are selected as those found to be associated to the reference cell
#'                                         types with F-statistics. If this option is selected, \code{REF_DATA_SET} and \code{REF_PHENO_COLUMN}
#'                                         need to be specified.
#'                                  \item{"\code{random}"} Sites are randomly selected.
#'                                  \item{"\code{pca}"} Sites are selected as those with most influence on the principal components.
#'                                  \item{"\code{var}"} Selects the most variable sites.
#'                                  \item{"\code{hybrid}"} Selects (N_MARKERS/2) most variable and (N_MARKERS/2) random sites.
#'                                  \item{"\code{range}"} Selects the sites with the largest difference between minimum and maximum
#'                                       across samples.
#'                                  \item{"\code{pcadapt}"} Uses principal component analysis as implemented in the \code{"bigstats"}
#'                                       R package to determine sites that are significantly linked to the potential cell types. This
#'                                       requires specifying K a priori (argument \code{K.prior}). We thank Florian Prive and Sophie
#'                                       Achard for providing the idea and parts of the codes.
#'                                  \item{"\code{edec_stage0}} Employs EDec's stage 0 to infer cell-type specific markers. By default
#'                                       EDec's example reference data is provided. If a specific data set is to be provided, it needs
#'                                       to be done through \code{REF_DATA_SET}.
#'                                  \item{"\code{custom}"} Specifying a custom file with indices.
#'                         }
#' @param n.markers The number of sites to be selected. Defaults to 5000.
#' @param write.files Flag indicating if the selected sites are to be stored on disk.
#' @param n.prin.comp Optional argument deteriming the number of prinicipal components used for selecting the most important sites.
#' @param range.diff Optional argument specifying the difference between maximum and minimum required.
#' @param custom.marker.file Optional argument containing a file that specifies the indices used for employing MeDeCom.
#' @param store.heatmaps Flag indicating if a heatmap of the selected input sites is to be create from the input methylation matrix.
#'                       The files are then stored in the 'heatmaps' folder in WD.
#' @param heatmap.sample.col Column name in the phenotypic table of \code{rnb.set}, used for creating a color scheme in the heatmap.
#' @param sample.subset Vector of indices of samples to be included in the analysis. If \code{NULL}, all samples are included.
#' @param k.fixed Columns in the T matrix that should be fixed. If \code{NULL}, no columns are fixed.
#' @param K.prior K determined from visual inspection. Only has an influence, if \code{MARKER_SELECTION="pcadapt"}.
#' @param factorviz.outputs Flag indicating, if outputs should be stored to be compatible with FactorViz for data exploration
#' @param opt.method Optimization method to be used. Either MeDeCom.quadPen or MeDeCom.cppTAfact (default).
#' @param startT Inital matrix for T.
#' @param startA Initial matrix for A.
#' @param folds Integer representing the number of folds used in the analysis.
#' @param cores Integer representing the number of cores to be used in the analysis.
#' @param itermax Maximum number of iterations
#' @param ninit Number if initialtions.
#' @param cluster.submit Flag indicating, if the jobs are to be submitted to a scientific compute cluster (only SGE supported).
#' @param cluster.Rdir Path to an executable version of R.
#' @param cluster.hostlist Regular expression, on which basis hosts are selected in the cluster environment.
#' @param cluster.memlimit the \code{memlimit} resource value of the cluster submission.
#' @param cleanup Flag indicating if temprary files are to be deleted.
#' @return An object of type \code{\link{MeDeComSet}} containing the results of the MeDeCom experiment.
#' @export
#' @author Michael Scherer

start_decomp_pipeline <- function(rnb.set,
                                  Ks,
                                  lambda.grid,
                                  work.dir=getwd(),
                                  factorviz.outputs=F,
                                  analysis.name="Analysis",
                                  sample.selection.col=NA,
                                  sample.selection.grep=NA,
                                  pheno.cols=NA,
                                  id.column=rnb.getOption("identifiers.column"),
                                  normalization="none",
                                  ref.ct.column=NA,
                                  ref.rnb.set=NULL,
                                  ref.rnb.ct.column=NA,
                                  prepare.true.proportions=F,
                                  true.A.token=NA,
                                  houseman.A.token=NA,
                                  estimate.houseman.prop=F,
                                  filter.beads=!is.null(rnb.set@covg.sites),
                                  min.n.beads=3,
                                  filter.intensity=inherits(rnb.set, "RnBeadRawSet"),
                                  min.int.quant = 0.1,
                                  max.int.quant = 0.95, 
                                  filter.na=TRUE,
                                  filter.context=TRUE,
                                  filter.cross.reactive=TRUE,
                                  execute.lump=FALSE,
                                  remove.ICA=FALSE,
                                  conf.fact.ICA=FALSE,
                                  ica.setting=NULL,
                                  filter.snp=TRUE,
                                  filter.somatic=TRUE,
                                  snp.list=NULL,
                                  filter.coverage = hasCovg(rnb.set),
                                  min.coverage=5,
                                  min.covg.quant=0.05,
                                  max.covg.quant=0.95,
                                  marker.selection="var",
                                  n.markers=5000,
                                  write.files=FALSE,
                                  n.prin.comp=10,
                                  range.diff=0.05,
                                  custom.marker.file="",
                                  store.heatmaps=F,
                                  heatmap.sample.col=NULL,
                                  sample.subset=NULL,
                                  k.fixed=NULL,
                                  K.prior=NULL,
                                  opt.method = "MeDeCom.cppTAfact",
                                  startT=NULL,
                                  startA=NULL,
                                  folds=10,
                                  cores=1,
                                  itermax=1000,
                                  ninit=100,
                                  cluster.submit=FALSE,
                                  cluster.Rdir=NA,
                                  cluster.hostlist="*",
                                  cluster.memlimit="5G",
                                  cleanup=FALSE
                                  ){
  if(inherits(rnb.set,"RnBeadSet")){
    data.prep <- prepare_data(RNB_SET=rnb.set,
                              WORK_DIR=work.dir,
                              analysis.name=analysis.name,
                              SAMPLE_SELECTION_COL=sample.selection.col,
                              SAMPLE_SELECTION_GREP=sample.selection.grep,
                              PHENO_COLUMNS=pheno.cols,
                              ID_COLUMN=id.column,
                              NORMALIZATION=normalization,
                              REF_CT_COLUMN=ref.ct.column,
                              REF_RNB_SET=ref.rnb.set,
                              REF_RNB_CT_COLUMN=ref.rnb.ct.column,
                              PREPARE_TRUE_PROPORTIONS=prepare.true.proportions,
                              TRUE_A_TOKEN=true.A.token,
                              HOUSEMAN_A_TOKEN=houseman.A.token,
                              ESTIMATE_HOUSEMAN_PROP=estimate.houseman.prop,
                              FILTER_BEADS=filter.beads,
                              MIN_N_BEADS=min.n.beads,
                              FILTER_INTENSITY=filter.intensity,
                              MIN_INT_QUANT = min.int.quant,
                              MAX_INT_QUANT = max.int.quant, 
                              FILTER_NA=filter.na,
                              FILTER_CONTEXT=filter.context,
                              FILTER_CROSS_REACTIVE=filter.cross.reactive,
                              execute.lump=execute.lump,
                              FILTER_SNP=filter.snp,
                              FILTER_SOMATIC=filter.somatic,
                              snp.list=snp.list,
                              remove.ICA=remove.ICA,
                              conf.fact.ICA=conf.fact.ICA,
                              ica.setting=ica.setting
    )
  }else if(inherits(rnb.set,"RnBiseqSet")){
    data.prep <- prepare_data_BS(RNB_SET = rnb.set,
                                WORK_DIR = work.dir,
                                analysis.name = analysis.name,
                                SAMPLE_SELECTION_COL = sample.selection.col,
                                SAMPLE_SELECTION_GREP = sample.selection.grep,
                                REF_CT_COLUMN=ref.ct.column,
                                PHENO_COLUMNS=pheno.cols,
                                PREPARE_TRUE_PROPORTIONS=prepare.true.proportions,
                                TRUE_A_TOKEN=true.A.token,
                                HOUSEMAN_A_TOKEN=houseman.A.token,
                                ID_COLUMN=id.column,
                                FILTER_COVERAGE = filter.coverage,
                                MIN_COVERAGE=min.coverage,
                                MIN_COVG_QUANT=min.covg.quant,
                                MAX_COVG_QUANT=max.covg.quant,
                                FILTER_NA=filter.na,
                                FILTER_SNP=filter.snp,
                                snp.list=snp.list,
                                FILTER_SOMATIC=filter.somatic
      
    )
  }
  cg_subsets <- prepare_CG_subsets(rnb.set=data.prep$rnb.set.filtered,
                                     MARKER_SELECTION=marker.selection,
                                     N_MARKERS=n.markers,
                                     WRITE_FILES=write.files,
                                     WD=work.dir,
                                     REF_DATA_SET=ref.rnb.set,
                                     REF_PHENO_COLUMN=ref.rnb.ct.column,
                                     N_PRIN_COMP=n.prin.comp,
                                     RANGE_DIFF=range.diff,
                                     CUSTOM_MARKER_FILE=custom.marker.file,
                                     store.heatmaps=store.heatmaps,
                                     heatmap.sample.col=heatmap.sample.col,
                                     K.prior = K.prior
  )
  if("RefMeth" %in% names(data.prep)){
    trueT <- data.prep$RefMeth
  }else{
    trueT <- NULL
  }
  if("RefProps" %in% names(data.prep)){
    trueA <- data.prep$RefProps
  }else{
    trueA <- NULL
  }
  medecom.result <- start_medecom_analysis(rnb.set=data.prep$rnb.set.filtered,
                                             WORK_DIR=work.dir,
                                             cg_groups=cg_subsets,
                                             Ks=Ks,
                                             LAMBDA_GRID=lambda.grid,
                                             SAMPLE_SUBSET=sample.subset,
                                             K_FIXED=k.fixed,
                                             WRITE_FILES=write.files,
                                             factorviz.outputs=factorviz.outputs,
                                             opt.method = opt.method,
                                             startT = startT,
                                             startA = startA,
                                             trueT = trueT,
                                             trueA = trueA,
                                             analysis.name=analysis.name,
                                             folds=folds,
                                             cores=cores,
                                             itermax=itermax,
                                             ninit=ninit,
                                             CLUSTER_SUBMIT=cluster.submit,
                                             CLUSTER_RDIR=cluster.Rdir,
                                             CLUSTER_HOSTLIST=cluster.hostlist,
                                             CLUSTER_MEMLIMIT=cluster.memlimit,
                                             CLEANUP=cleanup
  )
  return(medecom.result)
}

#' A small RnBeads object used to run the examples.
#' @name rnb.set.example
#' @docType data
#' @author Michael Scherer
NULL
