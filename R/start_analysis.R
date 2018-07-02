#' start_medecom_analysis
#' 
#' Wrapper for runMeDeCom, for data preprocessed through the DecombPipeline
#' 
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation and sample meta information.
#' @param WORK_DIR Working directory for the analysis.
#' @param cg_groups List of CpG indices used for the analysis. Can be computed by \code{\link{prepare_CG_subsets}}.
#' @param Ks Vector of integers used as components in MeDeCom.
#' @param LAMBDA_GRID Vector of doubles representing the regularization parameter in MeDeCom.
#' @param SAMPLE_SUBSET Vector of indices of samples to be included in the analysis. If \code{NULL}, all samples are included.
#' @param K_FIXED Columns in the T matrix that should be fixed. If \code{NULL}, no columns are fixed.
#' @param WRITE_FILES Flag indicating if intermediate results are to be stored.
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
#' @return An object of type \code{\link{MeDeComSet}} containing the results of the MeDeCom experiment.
#' @export
start_medecom_analysis<-function(
		rnb.set=NULL,
		WORK_DIR,
		cg_groups,
		Ks,
		LAMBDA_GRID,
		SAMPLE_SUBSET=NULL,
		K_FIXED=NULL,
		WRITE_FILES=TRUE,
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
		CLEANUP=FALSE
){
	library(MeDeCom)
	library(R.utils)
	
	#RDIR="/TL/deep-share/archive00/software/bin"
	#.libPaths(sprintf("%s/Rlib_test", WORK_DIR))
  
  WORK_DIR <- file.path(WORK_DIR,analysis.name)
	
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
	
	if(is.null(rnb.set)){
		load(sprintf("%s/data.set.RData", WORK_DIR))
	}else{
		meth.data<-meth(rnb.set)
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
	# if(!is.null(A_LOWER) && is.null(A_UPPER)){
	# 	saveRDS(A_LOWER, file=sprintf("%s/A_lower.RDS", WORK_DIR))
	# 	saveRDS(A_UPPER, file=sprintf("%s/A_upper.RDS", WORK_DIR))
	# }
	
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
	result@parameters$cg_subsets <- c(1:length(cg_groups))
	result@parameters$ITERMAX<-itermax
	result@parameters$NFOLDS<- folds
	result@parameters$NINIT<-ninit
	help<- NULL
	for ( i in 1:length(cg_groups)){
		help <- append(help, cg_groups[[i]])
	}
	result@parameters$ORIGINAL_GROUP_LISTS<- help
	
	if(WRITE_FILES){
		saveRDS(result, file=file.path(WORK_DIR, "collected.result.RDS"))
	}
	
	return(result)
}
