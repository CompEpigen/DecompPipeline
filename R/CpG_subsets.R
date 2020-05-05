#' prepare.CG.subsets
#' 
#' This routine selects a subset of CpGs sites used for MeDeCom analysis. Different selection methods are supported.
#' 
#' @name prepare_CG_subsets
#' 
#' @param meth.data A \code{matrix} or \code{data.frame} containing methylation information. If NULL, methylation information needs to be provided
#'                   through \code{rnb.set}
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation, sample and optional coverage information.
#' @param marker.selection A vector of strings representing marker selection methods. Available method are \itemize{
#'                                  \item{"\code{all}"} Using all sites available in the input.
#'                                  \item{"\code{pheno}"} Selected are the top \code{n.markers} site that differ between the phenotypic
#'                                         groups defined in data preparation or by \code{\link{rnb.sample.groups}}. Those are
#'                                         selected by employing limma on the methylation matrix.
#'                                  \item{"\code{houseman2012}"} 50k sites determined to be cell-type specific for blood cell types using the Houseman's reference-
#'                                         based deconvolution and the Reinius et al. reference data set. See Houseman et.al. 2012 and Reinis et.al. 2012.
#'                                         NOTE: This option should only be used for whole blood data generated using the 450k array.
#'                                  \item{"\code{houseman2014}"} Selects the sites said to be linked to cell type composition by \code{RefFreeEWAS},
#'                                         which is similar to surrogate variable analysis. See Houseman et.al. 2014.
#'                                  \item{"\code{jaffe2014}"} The 600 sites stated as related to cell-type composition Jaffe et.al. 2014.
#'                                         NOTE: This option should only be used for whole blood data generated using the 450k array.
#'                                  \item{"\code{rowFstat}"} Markers are selected as those found to be associated to the reference cell
#'                                         types with F-statistics. If this option is selected, \code{ref.rnb.set} and \code{ref.pheno.column}
#'                                         need to be specified.
#'                                  \item{"\code{random}"} Sites are randomly selected.
#'                                  \item{"\code{pca}"} Sites are selected as those with most influence on the principal components.
#'                                  \item{"\code{var}"} Selects the most variable sites.
#'                                  \item{"\code{hybrid}"} Selects (n.markers/2) most variable and (n.markers/2) random sites.
#'                                  \item{"\code{range}"} Selects the sites with the largest difference between minimum and maximum
#'                                       across samples.
#'                                  \item{"\code{pcadapt}"} Uses principal component analysis as implemented in the \code{"bigstats"}
#'                                       R package to determine sites that are significantly linked to the potential cell types. This
#'                                       requires specifying K a priori (argument \code{K.prior}). We thank Florian Prive and Sophie
#'                                       Achard for providing the idea and parts of the codes.
#'                                  \item{"\code{edec_stage0}} Employs EDec's stage 0 to infer cell-type specific markers. By default
#'                                       EDec's example reference data is provided. If a specific data set is to be provided, it needs
#'                                       to be done through \code{ref.rnb.set}.
#'                                  \item{"\code{custom}"} Specifying a custom file with indices.
#'                         }
#' @param n.markers The number of sites to be selected. Defaults to 5000.
#' @param remove.correlated Flag indicating if highly correlated sites are to be removed
#' @param cor.threshold Numeric indicating a threshold above which sites are not to be considered in the feature selection.
#'          If \code{"quantile"}, sites correlated higher than the 95th quantile are removed.
#' @param write.files Flag indicating if the selected sites are to be stored on disk.
#' @param out.dir Path to the working directory used for analyis, or data preparation.
#' @param ref.rnb.set An object of type \code{\link{RnBSet-class}} or a path to such an object stored on disk, 
#'                      if \code{rowFstat} is selected.
#' @param ref.pheno.column Optional argument stating the column name of the phenotypic table of \code{ref.rnb.set} with
#'                      the reference cell type.
#' @param n.prin.comp Optional argument deteriming the number of prinicipal components used for selecting the most important sites.
#' @param range.diff Optional argument specifying the difference between maximum and minimum required.
#' @param custom.marker.file Optional argument containing an absolute path to a file that specifies the indices used for employing MeDeCom. Can be provided
#'                   either as an \code{RDS} file containing a vector of indices to select or as a \code{txt, csv, tsv} file containing each index
#'                   to be selected as a single row.
#' @param store.heatmaps Flag indicating if a heatmap of the selected input sites is to be create from the input methylation matrix.
#'                       The files are then stored in the 'heatmaps' folder in out.dir.
#' @param heatmap.sample.col Column name in the phenotypic table of \code{rnb.set}, used for creating a color scheme in the heatmap.
#' @param K.prior K determined from visual inspection. Only has an influence, if \code{marker.selection="pcadapt"}.
#' @return List of indices, one entry for each marker selection method specified by \code{marker.selection}. The indices correspond
#'          to the sites that should be used in \code{rnb.set}.
#' @details For methods "\code{houseman2012}" and "\code{jaffe2014}", a predefined set of markers is used. Since those correspond to
#'          absolute indices on the chip, the provided \code{rnb.set} must not be preprocessed and therefore still contain all sites.
#'          For the other metods, you may used \code{\link{prepare.data}} to filter sites for quality and context.
#' @references \itemize{
#'             \item{1.} Houseman, E. A., Accomando, W. P., Koestler, D. C., Christensen, B. C., Marsit, C. J., Nelson, H. H., ..., Kelsey, K.
#'                 T. (2012). DNA methylation arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics, 13. 
#'             \item{2.} Reinius, L. E., Acevedo, N., Joerink, M., Pershagen, G., Dahlen, S. E., Greco, D., ..., A., & Kere, J.
#'                    (2012). Differential DNA methylation in purified human blood cells: Implications for cell lineage and studies on disease susceptibility.
#'                    PLoS ONE, 7(7). https://doi.org/10.1371/journal.pone.0041361
#'             \item{3.} Houseman, E. A., Molitor, J., & Marsit, C. J. (2014). Reference-free cell mixture adjustments in analysis of 
#'                 DNA methylation data. Bioinformatics, 30(10), 1431-1439. https://doi.org/10.1093/bioinformatics/btu029
#'             \item{4.} Jaffe, A. E., & Irizarry, R. A. (2014). Accounting for cellular heterogeneity is critical in epigenome-wide 
#'                 association studies. Genome Biology, 15(2), R31. https://doi.org/10.1186/gb-2014-15-2-r31
#'                  
#' }
#' @export
#' @author Michael Scherer, Pavlo Lutsik
prepare.CG.subsets<-function(
    meth.data=NULL,
		rnb.set=NULL,
		marker.selection,
		n.markers=5000,
		remove.correlated=FALSE,
		cor.threshold="quantile",
		write.files=FALSE,
		work.dir=NA,
		ref.rnb.set=NULL,
		ref.pheno.column=NULL,
		n.prin.comp=10,
		range.diff=0.05,
		custom.marker.file="",
		store.heatmaps=F,
		heatmap.sample.col=NULL,
		K.prior=NULL
		)
{
  require("RnBeads")
	cg.groups<-list()
	
	groups<-1:length(marker.selection)
	
  if(is.null(meth.data)){
    if(inherits(rnb.set,"RnBSet")){
	    meth.data<-meth(rnb.set)
    }else{
      stop("Invalid value for rnb.set, needs to be RnBSet")
    }
  }
	
	if(!is.data.frame(meth.data)&&!is.matrix(meth.data)){
	  stop("Invalid value for meth.data")
	}
	  
	
	if(store.heatmaps){
	  if(!is.null(heatmap.sample.col)){
	     if(!heatmap.sample.col %in% colnames(pheno(rnb.set))){
	       logger.error("heatmap.sample.col not a column name of the phenotypic table")
	    }
	    trait <- as.factor(pheno(rnb.set)[,heatmap.sample.col])
	    palette <- rainbow(length(levels(trait)))
	    sample.cols <- palette[as.integer(trait)]
	  }else{
	    sample.cols <- NULL
	  }
	  if(!dir.exists(file.path(out.dir,"heatmaps"))){
	    dir.create(file.path(out.dir,"heatmaps"))
	  }
	}
	
	if(n.markers>nrow(meth.data)){
	  logger.error("You cannot select more markers than you have sites")
	}
	
	ind.all <- 1:nrow(meth.data)
	if(remove.correlated){
	  if(inherits(rnb.set,'RnBSet')){
  	  logger.start("Removing highly correlated features")
  	  if(!isImputed(rnb.set)){
  	    meth.data <- rnb.execute.imputation(meth.data)
  	  }
  	  logger.start("Correlation matrix computation for chromosomes separately")
  	  anno.set <- annotation(rnb.set)
  	  all.cor <- rep(NA,nrow(anno.set))
  	  for(chr in unique(anno.set$Chromosome)){
  	    to.cor <- meth.data[anno.set$Chromosome %in% chr,]
  	    if(any(is.na(to.cor))){
  	      logger.error(paste("Missing values in correlation for chromosome",chr))
  	    }
    	  cor.mat <- abs(cor(t(to.cor)))
    	  diag(cor.mat) <- NA
    	  max.cor <- apply(cor.mat,1,mean,na.rm=T)
    	  all.cor[anno.set$Chromosome %in% chr] <- max.cor
  	  }
  	  logger.completed()
  	  logger.start("removing features")
  	  if(!is.numeric(cor.threshold)){
  	    if(cor.threshold != "quantile"){
  	      logger.error("cor.threshold only allows numeric values and 'quantile'")
  	    }
  	    cor.threshold <- quantile(all.cor,0.95,na.rm = T)
  	  }
  	  logger.completed()
  	  if(cor.threshold < 0 | cor.threshold > 1){
  	    logger.error("cor.threshold needs to be between 0 and 1")
  	  }
  	  logger.info(paste("cor.threshold is",cor.threshold))
  	  logger.info(paste("Removed ",sum(all.cor>cor.threshold),"sites in correlation filtering"))
  	  ind.all <- intersect(ind.all,which(all.cor<=cor.threshold))
  	  logger.completed()
    }else{
	    logger.warning("Removing highly correlated features only applicable for 'RnBSet' input")
	  }
	}
	
	for(group in groups){
		
		ind<-ind.all
		
		if(marker.selection[group]=="all"){
		  #Do nothing
		  ind <- ind
		}
		
		if(marker.selection[group]=="pheno"){
			
		  if(file.exists(sprintf("%s/pheno.RData",out.dir))){
			  load(sprintf("%s/pheno.RData",out.dir))
		    sel.columns <- gsub(" ","_",colnames(pheno.data))
		  }else{
		    pheno.data <- pheno(rnb.set)
		    sel.columns <- gsub(" ","_",names(rnb.sample.groups(rnb.set)))
		  }
			
			X <- meth.data
			ncgs<-nrow(X)
			X <- na.omit(X)

			colnames(pheno.data) <- gsub(" ","_",colnames(pheno.data))
			na.pheno <- lapply(sel.columns,function(x)is.na(pheno.data[,x]))
			rem.samples <- rep(FALSE,ncol(X))
			for(entry in na.pheno){
			  rem.samples <- rem.samples | entry
			}
			X <- X[,!rem.samples]
			formula.text <- paste0("~0+", paste(sel.columns,collapse="+"))
			design <- model.matrix(as.formula(formula.text), data=pheno.data)
			fit <- limma::lmFit(X,design)
			fit <- limma::eBayes(fit)
			tstatDeltaAll<-abs(fit$t)
			
			ranks<-apply(tstatDeltaAll, 2, rank)
			maxRank<-apply(ranks, 1, max)
		
			ind<-ind[order(maxRank)[1:min(n.markers, length(ind))]]
		}
		
		if(marker.selection[group]=="houseman2012" ){
		  logger.warning("The houseman 2012 method should only be used for blood data sets generated on the 450k array.")
		  loc <- system.file(file.path("extdata","houseman.50k.markers.RDS"),package="DecompPipeline")
		  if(file.exists(loc)){
			  houseman.50k.markers<-readRDS(loc)
			  ind<-intersect(ind, houseman.50k.markers)
		  }
		}
		
		if(marker.selection[group]=="houseman2014"){
			require(RefFreeEWAS)
			
		  if(file.exists(sprintf("%s/pheno.RData",out.dir))){
		    load(sprintf("%s/pheno.RData",out.dir))
		    sel.columns <- gsub(" ","_",colnames(pheno.data))
		  }else{
		    pheno.data <- pheno(rnb.set)
		    sel.columns <- gsub(" ","_",names(rnb.sample.groups(rnb.set)))
		  }
		  
			X <- meth.data
			ncgs<-nrow(X)
			X <- na.omit(X)

			colnames(pheno.data) <- gsub(" ","_",colnames(pheno.data))
			na.pheno <- lapply(sel.columns,function(x)is.na(pheno.data[,x]))
			rem.samples <- rep(FALSE,ncol(X))
			for(entry in na.pheno){
			  rem.samples <- rem.samples | entry
			}
			X <- X[,!rem.samples]
			pheno.data <- pheno.data[!rem.samples,]
			level.problem <- unlist(lapply(sel.columns,function(x,mat){
			  length(unique(mat[,x])) == length(levels(mat[,x]))
			},pheno.data))
			sel.columns <- sel.columns[level.problem]
			formula.text <- paste0("~0+", paste(sel.columns,collapse="+"))
			design <- model.matrix(as.formula(formula.text), data=pheno.data)
			tmpBstar <- (X %*% design %*% solve(t(design)%*%design))
				
			## rescaling the residuals, and addition by Andres
			R <- X-tmpBstar %*% t(design)
			rescale.residual<-TRUE
			if(rescale.residual){
				R <- t(scale(t(R)))
			}
			d<-EstDimRMT(R)$dim
				
			logger.info(c("Estimated number of latent components is", d))
			rf <- RefFreeEwasModel(X, design, d)
				
			logger.status("Fitted the RefFreeEWAS model")
				
			Delta <- rf$Bstar-rf$Beta
				
			nboot<-100
			#### compute se delta from a bootstrap
			rfBoot <- BootRefFreeEwasModel(rf,nboot)
			logger.status("Pefrormed the bootstrap")
				
			seDelta <- apply(rfBoot[,,"B*",]-rfBoot[,,"B",], 1:2, sd)
			tstatDelta <- -abs(Delta)/seDelta
			ncgs <- nrow(meth.data)
			nnas <- apply(is.na(meth.data),1,sum)
			notna.rows <- which(nnas==0)
			if(is.null(dim(tstatDelta))){
					tstatDeltaAll<-rep(NA, ncgs)
					tstatDeltaAll[notna.rows]<-tstatDelta
				}else{
					tstatDeltaAll<-matrix(NA_real_, ncol=ncol(tstatDelta), nrow=ncgs)
					tstatDeltaAll[notna.rows,]<-tstatDelta
				}
			if(is.null(dim(tstatDelta))){
				ind<-ind[order(-tstatDeltaAll)[1:min(10000, length(ind))]]
			}else{
				ranks<-apply(-tstatDeltaAll, 2, rank)
				maxRank<-apply(ranks, 1, max)
				ind<-ind[order(maxRank)[1:min(n.markers, length(ind))]]
			}
		}
		
		if(marker.selection[group]=="jaffe2014"){
		  logger.warning("The jaffe 2014 method should only be used for blood data sets generated on the 450k array.")
			jaffe.markers <- readRDS(system.file(file.path("extdata","jaffe.irrizzary.markers.600.RDS"),package="DecompPipeline"))
			ind<-intersect(ind, jaffe.markers)
		}
		
		if(marker.selection[group]==("rowFstat")){
		  ####### This needs to be checked: it doesn't seem to do what it should
  		if(!(is.null(ref.rnb.set) || is.null(ref.pheno.column))){
  			
  			require(genefilter)
  			
  			load.env<-new.env(parent=emptyenv())
  			
  			if(inherits(ref.rnb.set,"RnBSet")){
  			  rnb.ref.set <- ref.rnb.set
  			}else{
  			  rnb.ref.set <- load.rnb.set(ref.rnb.set)
  			}
  			meth.data.ref <- meth(rnb.ref.set)
  			pheno.data.ref <- pheno(rnb.ref.set)
  			
  			marker.fstat<-rowFtests(meth.data.ref[ind,], as.factor(pheno.data.ref[[ref.pheno.column]]))
  			
  			subset<-which(rank(-marker.fstat$statistic)<=n.markers)
  			ind<-ind[subset]
  		}else{
  		  logger.error("ref.rnb.set and ref.pheno.column need to be specified, if rowFstat is selected.")
  		}
		}
		
		if(marker.selection[group]=="random"){
			subset<-sample.int(length(ind), min(n.markers, length(ind)))
			ind<-ind[subset]
		}
		
		if(marker.selection[group]=="pca"){
			
			meth.data.subs<-meth.data[ind,]
			pca<-prcomp(t(meth.data.subs))
			rot<-pca$rotation[,1:min(n.prin.comp,ncol(pca$rotation))]
			
			pca.ind<-integer()
			add.sites <- ceiling(n.markers/n.prin.comp)
			for(cix in 1:min(n.prin.comp,ncol(pca$rotation))){
				pca.ind<-union(pca.ind,order(abs(rot[,cix]), decreasing = TRUE)[1:add.sites])
			}
			ind<-ind[sort(pca.ind)]
		}
		
		if(marker.selection[group] == "var"){
			
			sds<-apply(meth.data[ind,], 1, sd)
			ind<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),n.markers)]]
		}
		
		if(marker.selection[group] == "hybrid"){
			
			sds<-apply(meth.data[ind,], 1, sd)
			
			var.set<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),ceiling(n.markers/2))]]
			random.set<-sample(setdiff(ind,var.set), floor(n.markers/2))
			
			ind<-sort(c(var.set, random.set))
			rm(var.set)
			rm(random.set)
		}
		
		if(marker.selection[group] == "range"){
			ranges<-apply(meth.data[ind,], 1, range)
			ranges<-ranges[2,]-ranges[1,]
			ind<-ind[ranges>range.diff]
			rm(ranges)
		}
		
		if(marker.selection[group] == "pcadapt"){
		  require("bigstatsr")
		  require("robust")
		  if(is.null(K.prior)){
		    stop("K.prior needs to be specific for pcadapat marker selection method")
		  }
		  meth.t <- t(meth.data[ind,])
		  meth.t <- as_FBM(meth.t)
		  svd <- big_SVD(meth.t,fun.scaling = big_scale(),k=K.prior)
		  singular.vectors <- svd$u
		  z.scores <- sapply(cols_along(singular.vectors), function(k) {
		    big_univLinReg(meth.t, singular.vectors[, k])$score
		  })
		  dist.math <- robust::covRob(z.scores, estim = "pairwiseGK")$dist
		  lpval <- pchisq(dist.math, df = K.prior, lower.tail = FALSE, log.p = TRUE) / log(10)
		  ind <- ind[order(lpval,decreasing = F)[1:n.markers]]
		}
		
		if(marker.selection[group]=="edec_stage0"){
		  require("EDec")
		  require("EDecExampleData")
		  if(is.null(ref.rnb.set)){
		    markers <- run_edec_stage_0(reference_meth = EDecExampleData::reference_meth,
		                                reference_classes = EDecExampleData::reference_meth_class,
		                                max_p_value = 1e-5,
		                                num_markers = n.markers)
		  }else{
		    if(inherits(ref.rnb.set,"RnBSet")){
		      rnb.ref.set <- ref.rnb.set
		    }else{
		      rnb.ref.set <- load.rnb.set(ref.rnb.set)
		    }
		    if(!ref.pheno.column %in% colnames(pheno(rnb.ref.set))){
		      stop("Supplied ref.pheno.column not in phenotypic information")
		    }
		    markers <- run_edec_stage_0(reference_meth = meth(rnb.ref.set),
		                                reference_classes = pheno(rnb.ref.set)[,ref.pheno.column],
		                                max_p_value = 1e-5,
		                                num_markers = n.markers)
		  }
		  if(is.null(rnb.set)){
		    if(!grep("cg",row.names(meth.data))){
		      stop("Row names of meth.data (cg-identifiers) need to be provided for EDec")
		    }
		    ind <- which(row.names(meth.data) %in% markers)
		  }else{
		    ind <- which(row.names(annotation(rnb.set)) %in% markers)
		  }
		}
		
		if(marker.selection[group]=="custom"){
			if(file.exists(custom.marker.file)){
			  if(grepl("RDS",custom.marker.file)){
				  custom_filter<-readRDS(custom.marker.file)
			  }else if(grepl("txt|csv|tsv",custom.marker.file)){
			    custom_filter <- as.numeric(readLines(custom.marker.file))
			  }else{
			    stop("Invalid value for custom.marker.file")
			  }
			  ind<-intersect(ind, custom_filter)
			}
		}
		
		if(store.heatmaps){
		  create.heatmap(meth.data[ind,],trait,sample.cols,palette,out.dir,marker.selection[group])
		}
		
		if(write.files){
			saveRDS(ind, file=sprintf("%s/cg_group_%d.RDS", out.dir, group))
		}
		
		cg.groups[[group]]<-ind
		
	}
	names(cg.groups) <- marker.selection

	return(cg.groups=cg.groups)
}


