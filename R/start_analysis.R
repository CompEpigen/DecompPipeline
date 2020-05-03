#' start.analysis
#' 
#' Wrapper function to start one of the deconvolution algorithms \code{MeDeCom}, \code{RefFreeEWAS} or \code{EDec}
#' 
#' @param meth.data A \code{matrix} or \code{data.frame} containing methylation information. If NULL, methylation information needs to be provided
#'                   through \code{rnb.set}
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation and sample meta information.
#' @param work.dir Working directory for the analysis.
#' @param cg.groups List of CpG indices used for the analysis. Can be computed by \code{\link{prepare.CG.subsets}}.
#' @param Ks Vector of integers used as components in MeDeCom.
#' @param factorviz.outputs Flag indicating, if outputs should be stored to be compatible with FactorViz for data exploration
#' @param method The method to be used for deconvolution. Can be one of \code{MeDeCom}, \code{RefFreeCellMix} or \code{EDec}.
#' @author Michael Scherer
#' @export
start.analysis <- function(meth.data=NULL,
                           rnb.set=NULL,
                           cg.groups,
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
                                     cg.groups = cg.groups,
                                     Ks = Ks,
                                     work.dir = work.dir,
                                     factorviz.outputs = factorviz.outputs,
                                     ...)
  }else if(method == "RefFreeCellMix"){
    md.res <- start.refreeewas.analysis(meth.data=meth.data,
                                                    rnb.set=rnb.set,
                                                    cg.groups=cg.groups,
                                                    Ks=Ks,
                                                    work.dir=work.dir,
                                                    factorviz.outputs=factorviz.outputs)
  }else if(method == "EDec"){
    md.res <- start.edec.analysis(meth.data=meth.data,
                                  rnb.set=rnb.set,
                                  cg.groups=cg.groups,
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
#' @param cg.groups List of CpG indices used for the analysis. Can be computed by \code{\link{prepare.CG.subsets}}.
#' @param Ks The number of cell types to be tested. Can be a single numeric value or an array of numbers.
#' @param work.dir The working directory to be used.
#' @param factorviz.outputs Flag indicating, if outputs should be stored to be compatible with FactorViz for data exploration
#' @return An object of type \code{\link{MeDeComSet}} containing the results of the EDec experiment.
#' @author Michael Scherer
#' @export
start.edec.analysis <- function(meth.data=NULL,
                                rnb.set=NULL,
                                cg.groups,
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
  for(i.group in 1:length(cg.groups)){
    logger.start(paste("Processing group:",i.group))
    group <- cg.groups[[i.group]]
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
  result <- as.MeDeComSet(res.all,cg_subsets=1:length(cg.groups),Ks=Ks,rss=rss.all,m.orig=nrow(meth.data),n.orig=ncol(meth.data))
  result@parameters$GROUP_LISTS <- cg.groups
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
#' @param cg.groups List of CpG indices used for the analysis. Can be computed by \code{\link{prepare.CG.subsets}}.
#' @param Ks The number of cell types to be tested. Can be a single numeric value or an array of numbers.
#' @param work.dir The working directory to be used.
#' @param factorviz.outputs Flag indicating, if outputs should be stored to be compatible with FactorViz for data exploration
#' @return An object of type \code{\link{MeDeComSet}} containing the results of the RefFreeCellMix experiment.
#' @author Michael Scherer
#' @export
start.refreeewas.analysis <- function(meth.data=NULL,
                                      rnb.set=NULL,
                                      cg.groups,
                                      Ks,
                                      work.dir=getwd(),
                                      factorviz.outputs=FALSE){
  if(!requireNamespace("RefFreeEWAS")){
    stop("Please install ReFreeEWAS")
  }else{
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
    for(i.group in 1:length(cg.groups)){
      logger.start(paste("Processing group:",i.group))
      group <- cg_groups[[i.group]]
      meth.sset <- meth.data[group,]
      res.sset <- RefFreeEWAS::RefFreeCellMixArray(meth.sset,Klist=Ks)
      devis <- tryCatch(RefFreeEWAS::RefFreeCellMixArrayDevianceBoots(res.sset,Y=meth.sset),error=function(e)e)
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
}

#' start.medecom.analysis
#' 
#' Wrapper for runMeDeCom, for data preprocessed through the DecompPipeline
#' 
#' @param meth.data A \code{matrix} or \code{data.frame} containing methylation information. If NULL, methylation information needs to be provided
#'                   through \code{rnb.set}
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation and sample meta information.
#' @param work.dir Working directory for the analysis.
#' @param cg.groups List of CpG indices used for the analysis. Can be computed by \code{\link{prepare.CG.subsets}}.
#' @param Ks Vector of integers used as components in MeDeCom.
#' @param lambda.grid Vector of doubles representing the regularization parameter in MeDeCom.
#' @param sample.subset Vector of indices of samples to be included in the analysis. If \code{NULL}, all samples are included.
#' @param k.fixed Columns in the T matrix that should be fixed. If \code{NULL}, no columns are fixed.
#' @param write.files Flag indicating if intermediate results are to be stored.
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
#' @param cluster.submit Flag indicating, if the jobs are to be submitted to a scientific compute cluster (only SGE supported).
#' @param cluster.R.dir Path to an executable version of R.
#' @param cluster.hostlist Regular expression, on which basis hosts are selected in the cluster environment.
#' @param cluster.memlimit the \code{memlimit} resource value of the cluster submission.
#' @param cleanup Flag indicating if temprary files are to be deleted.
#' @param analysis.info Information to be saved about the analysis. Just stored as info.
#' @param lambda.grid_TYPE String represent the lambda grid that was chosen. Just stored as info.
#' @param analysis.token String specifying the type of analysis that was conducted. Just stored as info.
#' @return An object of type \code{\link{MeDeComSet}} containing the results of the MeDeCom experiment.
#' @export
#' @author Pavlo Lutsik, Michael Scherer
start.medecom.analysis<-function(
    meth.data=NULL,
		rnb.set=NULL,
		work.dir=getwd(),
		cg.groups,
		Ks,
		lambda.grid,
		sample.subset=NULL,
		k.fixed=NULL,
		write.files=TRUE,
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
		cluster.submit=FALSE,
		cluster.R.dir=NA,
		cluster.hostlist="*",
		cluster.memlimit="5G",
		cleanup=FALSE,
		analysis.info=NULL,
		lambda.grid_TYPE="standard",
		analysis.token="customAnalysis"
){
	require(MeDeCom)
	library(R.utils)
	
  
  if(!is.null(analysis.info)){
    ANALYSIS_ID<-paste(
      analysis.info$DATASET, 
      analysis.info$DATA_SUBSET,
      analysis.info$NORMALIZATION,
      analysis.info$QUALITY_FILTERING, 
      analysis.info$MARKER_SELECTION,
      lambda.grid_TYPE,
      analysis.token,
      sep="_")
  }else{
    ANALYSIS_ID<-"customAnalysis"
    analysis.info<-list()
  }
  analysis.info$ANALYSIS<-ANALYSIS_ID
  
  work.dir <- file.path(work.dir,analysis.name)
  if(!file.exists(work.dir)){
    dir.create(work.dir)
  }
  log.file <- file.path(work.dir,"analysis.log")
  if(!file.exists(log.file)){
    if(logger.isinitialized()){
      logger.close()
      logger.start(fname=log.file)
    }else{
      logger.start(fname=log.file)
    }
  }
	
	if(is.na(cluster.R.dir)){
		cluster.R.dir="/usr/bin/"
	}
			
  if(is.null(meth.data)){
  	if(is.null(rnb.set)){
  		load(sprintf("%s/data.set.RData", work.dir))
  	}else{
  		meth.data<-meth(rnb.set)
  	}
  }
		
	if(write.files){
		saveRDS(lambda.grid, file=sprintf("%s/lambda.grid.RDS", work.dir))
		if(!is.null(sample.subset)){
			saveRDS(sample.subset, file=sprintf("%s/sample.subset.RDS", work.dir))
		}else{
			sample.subset<-1:ncol(meth.data)
		}
	}
	
	groups<-1:length(cg.groups)
	
	if(!is.null(k.fixed)){
		saveRDS(k.fixed, file=sprintf("%s/fixed_T_cols.RDS", work.dir))
	}else{
		k.fixed<-NULL
	}
		
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
	if(cluster.submit){
		cluster.settings=list(R_bin_dir=cluster.R.dir, host_pattern=cluster.hostlist, mem_limit=cluster.memlimit)
	}else{
		cluster.settings=NULL
	}
	
	result<-runMeDeCom(data=meth.data, 
			Ks=Ks,
			lambdas=lambda.grid,
			opt.method=opt.method,
			cg_subsets=cg.groups,
			sample.subset=sample.subset,
			startT=startT,
			startA=startA,
			trueT=trueT,
			trueA=trueA,
			fixed_T_cols=k.fixed,
			NINIT=ninit, 
			ITERMAX=itermax, 
			NFOLDS=folds,
			N_COMP_LAMBDA=4,
			NCORES=cores,
			analysis.name=analysis.name,
			use.ff=FALSE,
			cluster.settings=cluster.settings,
			temp.dir=work.dir,
			cleanup=cleanup,
			verbosity=1L,
			time.stamps=TRUE
	)
	
	result@parameters$ANALYSIS <- analysis.name
	result@parameters$GROUP_LISTS <- cg.groups
	result@parameters$ITERMAX<-itermax
	result@parameters$NFOLDS<-folds
	result@parameters$NINIT<-ninit
	
	analysis.info$GROUP_LISTS<-cg.groups
	analysis.info$ANALYSIS_DATE<-date()
	result@dataset_info<-c(result@dataset_info, analysis.info)
	
	if(write.files){
		saveRDS(result, file=file.path(work.dir, "collected.result.RDS"))
	}
	
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

#' start.decomp.pipeline
#' 
#' Main workhorse of the DecompPipeline R-package. Performs preprocessing (\code{\link{prepare.data}} or \code{\link{prepare.data.BS}}),
#' CpG subset selection (\code{\link{prepare.CG.subsets}}) and deconvolution (\code{\link{start.medecom.analysis}}, \code{\link{start.refreeewas.analysis}}, \code{\link{start.edec.analysis}})
#' 
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} for which analysis is to be performed.
#' @param Ks Vector of integers used as components in MeDeCom.
#' @param lambda.grid Vector of doubles representing the regularization parameter in MeDeCom.
#' @param work.dir A path to a existing directory, in which the results are to be stored
#' @param analysis.name A string representing the dataset for which analysis is to be performed. Only used to create a folder with a 
#'                 descriptive name of the analysis.
#' @param sample.selection.col A column name in the phenotypic table of \code{rnb.set} used to selected a subset of samples for
#'                 analysis that contain the string given in \code{sample.selection.col}.
#' @param sample.selection.grep A string used for selecting samples in the column \code{sample.selection.col}.
#' @param pheno.cols Vector of column names in the phenotypic table of \code{rnb.set} that is kept and exported for further 
#'                 exploration.
#' @param id.column Sample-specific ID column name in \code{rnb.set}
#' @param normalization normalization method to be performed before employing MeDeCom. Can be one of \code{"none","dasen","illumina","noob"} (BeadChip only).
#' @param ref.ct.column Column name in \code{rnb.set} used to extract methylation information on the reference cell types.
#' @param ref.rnb.set An object of type \code{\link[RnBeads]{RnBSet-class}} containing methylation information on reference cell types (BeadChip only).
#' @param ref.rnb.ct.column Column name in \code{ref.rnb.set} used to extract methylation information on the reference cell types (BeadChip only).
#' @param prepare.true.proportions Flag indicating if true proportions are either available in \code{rnb.set} or to be estimated 
#'                          with Houseman's reference-based deconvolution approach (BeadChip only).
#' @param true.A.token String present in the column names of \code{rnb.set} used for selecting the true proportions of the corresponding
#'                      cell types.
#' @param houseman.A.token Similar to \code{true.A.token}, but not containing the true proportions, rather the estimated proportions
#'                      by Houseman's method (BeadChip only).
#' @param estimate.houseman.prop If neither \code{true.A.token} nor \code{houseman.A.token} are given, the proportions of the reference
#'                      cell type are estimated with Houseman's approach (BeadChip only).
#' @param filter.beads Flag indicating, if site-filtering based on the number of beads available is to be conducted (BeadChip only).
#' @param min.n.beads Minimum number of beads required in each sample for the site to be considered for adding to MeDeCom (BeadChip only).
#' @param filter.intensity  Flag indicating if sites should be removed according to the signal intensities (the lowest and highest quantiles
#'                      given by \code{min.int.quant} and \code{max.int.quant}) (BeadChip only).
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
#' @param filter.snp Flag indicating if annotated SNPs are to be removed from the list of sites according to RnBeads' SNP list. 
#' @param snp.list Path to a file containing CpG IDs of known SNPs to be removed from the analysis, if \code{filter.snp} is \code{TRUE}.
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
#' @param remove.correlated Flag indicating if highly correlated features are to be removed.
#' @param cor.threshold Numeric indicating a threshold above which sites are not to be considered in the feature selection.
#'          If \code{"quantile"}, sites correlated higher than the 95th quantile are removed.
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

start.decomp.pipeline <- function(rnb.set,
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
                                  min.int.quant = 0.001,
                                  max.int.quant = 0.999, 
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
                                  remove.correlated=FALSE,
                                  cor.threshold="quantile",
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
    data.prep <- prepare.data(rnb.set=rnb.set,
                              work.dir=work.dir,
                              analysis.name=analysis.name,
                              sample.selection.col=sample.selection.col,
                              sample.selection.col=sample.selection.grep,
                              pheno.columns=pheno.cols,
                              id.column=id.column,
                              normalization=normalization,
                              ref.ct.column=ref.ct.column,
                              ref.rnb.set=ref.rnb.set,
                              ref.rnb.ct.column=ref.rnb.ct.column,
                              prepare.true.proportions=prepare.true.proportions,
                              true.A.token=true.A.token,
                              houseman.A.token=houseman.A.token,
                              estimate.houseman.prop=estimate.houseman.prop,
                              filter.beads=filter.beads,
                              min.n.beads=min.n.beads,
                              filter.intensity=filter.intensity,
                              min.int.quant = min.int.quant,
                              max.int.quant = max.int.quant, 
                              filter.na=filter.na,
                              filter.context=filter.context,
                              filter.cross.reactive=filter.cross.reactive,
                              execute.lump=execute.lump,
                              filter.snp=filter.snp,
                              filter.somatic=filter.somatic,
                              snp.list=snp.list,
                              remove.ICA=remove.ICA,
                              conf.fact.ICA=conf.fact.ICA,
                              ica.setting=ica.setting
    )
  }else if(inherits(rnb.set,"RnBiseqSet")){
    data.prep <- prepare.data.BS(rnb.set = rnb.set,
                                work.dir = work.dir,
                                analysis.name = analysis.name,
                                sample.selection.col = sample.selection.col,
                                sample.selection.col = sample.selection.grep,
                                ref.ct.column=ref.ct.column,
                                pheno.columns=pheno.cols,
                                prepare.true.proportions=prepare.true.proportions,
                                true.A.token=true.A.token,
                                houseman.A.token=houseman.A.token,
                                id.column=id.column,
                                filter.coverage = filter.coverage,
                                min.coverage=min.coverage,
                                min.covg.quant=min.covg.quant,
                                max.covg.quant=max.covg.quant,
                                filter.na=filter.na,
                                filter.snp=filter.snp,
                                snp.list=snp.list,
                                filter.somatic=filter.somatic,
                                execute.lump=execute.lump
      
    )
  }
  cg_subsets <- prepare.CG.subsets(rnb.set=data.prep$rnb.set.filtered,
                                     marker.selection=marker.selection,
                                     n.markers=n.markers,
                                     remove.correlated = remove.correlated,
                                     cor.threshold = cor.threshold,
                                     write.files=write.files,
                                     work.dir=work.dir,
                                     ref.rnb.set=ref.rnb.set,
                                     ref.pheno.column=ref.rnb.ct.column,
                                     n.prin.comp=n.prin.comp,
                                     range.diff=range.diff,
                                     custom.marker.file=custom.marker.file,
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
  medecom.result <- start.medecom.analysis(rnb.set=data.prep$rnb.set.filtered,
                                             work.dir=work.dir,
                                             cg.groups=cg_subsets,
                                             Ks=Ks,
                                             lambda.grid=lambda.grid,
                                             sample.subset=sample.subset,
                                             k.fixed=k.fixed,
                                             write.files=write.files,
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
                                             cluster.submit=cluster.submit,
                                             cluster.R.dir=cluster.Rdir,
                                             cluster.hostlist=cluster.hostlist,
                                             cluster.memlimit=cluster.memlimit,
                                             cleanup=cleanup
  )
  return(medecom.result)
}

