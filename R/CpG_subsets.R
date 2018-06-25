#' prepare_CG_subsets
#' 
#' This routine selects a subset of CpGs sites used for MeDeCom analysis. Different selection methods are supported.
#' 
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation, sample and optional coverage information.
#' @param MARKER_SELECTION A vector of strings representing marker selection methods. Available method are \itemize{
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
#'                                  \item{"\code{custom}"} Specifying a custom file with indices.
#'                         }
#' @param N_MARKERS The number of sites to be selected. Defaults to 5000.
#' @param WRITE_FILES Flag indicating if the selected sites are to be stored on disk.
#' @param WD Path to the working directory used for analyis, or data preparation.
#' @param REF_DATA_SET Name of the reference data set in \code{WD}, if \code{rowFstat} is selected.
#' @param REF_PHENO_COLUMN Optional argument stating the column name of the phenotypic table of \code{REF_DATA_SET} with
#'                      the reference cell type.
#' @param N_PRIN_COMP Optional argument deteriming the number of prinicipal components used for selecting the most important sites.
#' @param RANGE_DIFF Optional argument specifying the difference between maximum and minimum required.
#' @param CUSTOM_MARKER_FILE Optional argument containing a file that specifies the indices used for employing MeDeCom.
#' @return List of indices, one entry for each marker selection method specified by \code{MARKER_SELECTION}. The indices correspond
#'          to the sites that should be used in \code{rnb.set}.
#' @details For methods "\code{houseman2012}" and "\code{jaffe2014}", a predefined set of markers is used. Since those correspond to
#'          absolute indices on the chip, the provided \code{rnb.set} must not be preprocessed and therefore still contain all sites.
#'          For the other metods, you may used \code{\link{prepare_data}} to filter sites for quality and context.
#' @references \itemize{
#'             \item{1.} Houseman, E. A., Accomando, W. P., Koestler, D. C., Christensen, B. C., Marsit, C. J., Nelson, H. H., ..., Kelsey, K.
#'                 T. (2012). DNA methylation arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics, 13. 
#'             \item{2.} Houseman, E. A., Molitor, J., & Marsit, C. J. (2014). Reference-free cell mixture adjustments in analysis of 
#'                 DNA methylation data. Bioinformatics, 30(10), 1431-1439. https://doi.org/10.1093/bioinformatics/btu029
#'             \item{3.} Jaffe, A. E., & Irizarry, R. A. (2014). Accounting for cellular heterogeneity is critical in epigenome-wide 
#'                 association studies. Genome Biology, 15(2), R31. https://doi.org/10.1186/gb-2014-15-2-r31
#' }
#' @export
prepare_CG_subsets<-function(
		rnb.set,
		MARKER_SELECTION,
		N_MARKERS=5000,
		WRITE_FILES=FALSE,
		WD=NA,
		REF_DATA_SET=NULL,
		REF_PHENO_COLUMN=NULL,
		N_PRIN_COMP=10,
		RANGE_DIFF=0.05,
		CUSTOM_MARKER_FILE=""
		)
{
  require("RnBeads")
	cg_groups<-list()
	
	groups<-1:length(MARKER_SELECTION)
	
	meth.data<-meth(rnb.set)
	
	for(group in groups){
		
		ind<-1:nrow(meth.data)	
		
		if(MARKER_SELECTION[group]=="pheno"){
			
		  if(file.exists(sprintf("%s/pheno.RData",WD))){
			  load(sprintf("%s/pheno.RData",WD))
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
		
			ind<-ind[order(maxRank)[1:min(N_MARKERS, length(ind))]]
		}
		
		if(MARKER_SELECTION[group]=="houseman2012" ){
		  loc <- system.file(file.path("extdata","houseman.50k.markers.RDS"),package="DecompPipeline")
		  if(file.exists(loc)){
			  houseman.50k.markers<-readRDS(loc)
			  ind<-intersect(ind, houseman.50k.markers)
		  }
		}
		
		if(MARKER_SELECTION[group]=="houseman2014"){
			require(RefFreeEWAS)
			
		  if(file.exists(sprintf("%s/pheno.RData",WD))){
		    load(sprintf("%s/pheno.RData",WD))
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
			if(FALSE){
				for(id in 1:10){
					system(sprintf("qsub -cwd -j y -o %s/rf_boot_%d.log -b y -V -N rfb_%d_%s -l h='%s' -l mem_free=%s %s/Rscript %s/rfEwasBoot.R %s %d", WD, id, id, ANALYSIS, HOSTLIST, MEMLIMIT, ANALYSIS, RDIR, SRCDIR, WD, id))
				}
				# load all bootstraps
				load("rfBoot-1.RData")
				rfBootGrandSum <- rfBoot.sum
				rfBootGrandSum2 <- rfBoot.sum2 
				for(i in 2:10){ # Now combine
					load(file.path(WD, "rfBoot-",i,".RData",sep=""))
					rfBootGrandSum <- rfBootGrandSum + rfBoot.sum
					rfBootGrandSum2 <- rfBootGrandSum2 + rfBoot.sum2 
				}
				# Calculate SE's from summary data
				rfBootGrandMean <- rfBootGrandSum/100
				rfBootGrandMean2 <- rfBootGrandSum2/100
				rfBootSE <- sqrt((100/99)*(rfBootGrandMean2-rfBootGrandMean*rfBootGrandMean))
		  }
			
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
				ind<-ind[order(maxRank)[1:min(N_MARKERS, length(ind))]]
			}
		}
		
		if(MARKER_SELECTION[group]=="jaffe2014"){
			jaffe.markers <- readRDS(system.file(file.path("extdata","jaffe.irrizzary.markers.600.RDS"),package="DecompPipeline"))
			ind<-intersect(ind, jaffe.markers)
		}
		
		if(MARKER_SELECTION[group]==("rowFstat")){
		  ####### This needs to be checked: it doesn't seem to do what it should
  		if(!(is.null(REF_DATA_SET) || is.null(REF_PHENO_COLUMN))){
  			
  			require(genefilter)
  			
  			load.env<-new.env(parent=emptyenv())
  			
  			load(file.path(WD, REF_DATA_SET, "data.set.RData"), envir = load.env)
  			load(file.path(WD, REF_DATA_SET, "pheno.RData"), envir = load.env)
  			
  			meth.data.ref<-get("trueT", envir = load.env)
  			pheno.data.ref<-get("pd.ref", envir = load.env)
  			
  			marker.fstat<-rowFtests(meth.data.ref[ind,], unique(as.factor(pheno.data.ref[[REF_PHENO_COLUMN]])))
  			
  			subset<-which(rank(-marker.fstat$statistic)<=N_MARKERS)
  			ind<-ind[subset]
  		}else{
  		  logger.error("REF_DATA_SET and REF_PHENO_COLUMN need to be specified, if rowFstat is selected.")
  		}
		}
		
		if(MARKER_SELECTION[group]=="random"){
			subset<-sample.int(length(ind), min(N_MARKERS, length(ind)))
			ind<-ind[subset]
		}
		
		if(MARKER_SELECTION[group]=="pca"){
			
			meth.data.subs<-meth.data[ind,]
			pca<-prcomp(t(meth.data.subs))
			rot<-pca$rotation[,1:min(N_PRIN_COMP,ncol(pca$rotation))]
			
			pca.ind<-integer()
			add.sites <- ceiling(N_MARKERS/N_PRIN_COMP)
			for(cix in 1:min(N_PRIN_COMP,ncol(pca$rotation))){
				pca.ind<-union(pca.ind,order(abs(rot[,cix]), decreasing = TRUE)[1:add.sites])
			}
			ind<-ind[sort(pca.ind)]	
		}
		
		if(MARKER_SELECTION[group] == "var"){
			
			sds<-apply(meth.data[ind,], 1, sd)
			ind<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),N_MARKERS)]]
		}
		
		if(MARKER_SELECTION[group] == "hybrid"){
			
			sds<-apply(meth.data[ind,], 1, sd)
			
			var.set<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),ceiling(N_MARKERS/2))]]
			random.set<-sample(setdiff(ind,var.set), floor(N_MARKERS/2))
			
			ind<-sort(c(var.set, random.set))
			rm(var.set)
			rm(random.set)
		}
		
		if(MARKER_SELECTION[group] == "range"){
			ranges<-apply(meth.data[ind,], 1, range)
			ranges<-ranges[2,]-ranges[1,]
			ind<-ind[ranges>RANGE_DIFF]
			rm(ranges)
		}
		
		if(MARKER_SELECTION[group]=="custom"){
			if(file.exists(sprintf("%s/%s",WD, CUSTOM_MARKER_FILE))){
				custom_filter<-readRDS(sprintf("%s/%s",WD,CUSTOM_MARKER_FILE))
				ind<-intersect(ind, custom_filter)
			}
		}
		
		if(WRITE_FILES){
			saveRDS(ind, file=sprintf("%s/cg_group_%d.RDS", WD, group))
		}
		
		cg_groups[[group]]<-ind
		
	}
	return(cg_groups)
	
}