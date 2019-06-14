#' prepare_CG_subsets
#' 
#' This routine selects a subset of CpGs sites used for MeDeCom analysis. Different selection methods are supported.
#' 
#' @param meth.data A \code{matrix} or \code{data.frame} containing methylation information. If NULL, methylation information needs to be provided
#'                   through \code{rnb.set}
#' @param rnb.set An object of type \code{\link{RnBSet-class}} containing methylation, sample and optional coverage information.
#' @param MARKER_SELECTION A vector of strings representing marker selection methods. Available method are \itemize{
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
#' @param N_MARKERS The number of sites to be selected. Defaults to 5000.
#' @param WRITE_FILES Flag indicating if the selected sites are to be stored on disk.
#' @param WD Path to the working directory used for analyis, or data preparation.
#' @param REF_DATA_SET An object of type \code{\link{RnBSet-class}} or a path to such an object stored on disk, 
#'                      if \code{rowFstat} is selected.
#' @param REF_PHENO_COLUMN Optional argument stating the column name of the phenotypic table of \code{REF_DATA_SET} with
#'                      the reference cell type.
#' @param N_PRIN_COMP Optional argument deteriming the number of prinicipal components used for selecting the most important sites.
#' @param RANGE_DIFF Optional argument specifying the difference between maximum and minimum required.
#' @param CUSTOM_MARKER_FILE Optional argument containing an absolute path to a file that specifies the indices used for employing MeDeCom. Can be provided
#'                   either as an \code{RDS} file containing a vector of indices to select or as a \code{txt, csv, tsv} file containing each index
#'                   to be selected as a single row.
#' @param store.heatmaps Flag indicating if a heatmap of the selected input sites is to be create from the input methylation matrix.
#'                       The files are then stored in the 'heatmaps' folder in WD.
#' @param heatmap.sample.col Column name in the phenotypic table of \code{rnb.set}, used for creating a color scheme in the heatmap.
#' @param K.prior K determined from visual inspection. Only has an influence, if \code{MARKER_SELECTION="pcadapt"}.
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
    meth.data=NULL,
		rnb.set=NULL,
		MARKER_SELECTION,
		N_MARKERS=5000,
		WRITE_FILES=FALSE,
		WD=NA,
		REF_DATA_SET=NULL,
		REF_PHENO_COLUMN=NULL,
		N_PRIN_COMP=10,
		RANGE_DIFF=0.05,
		CUSTOM_MARKER_FILE="",
		store.heatmaps=F,
		heatmap.sample.col=NULL,
		K.prior=NULL
		)
{
  require("RnBeads")
	cg_groups<-list()
	
	groups<-1:length(MARKER_SELECTION)
	
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
	  if(!dir.exists(file.path(WD,"heatmaps"))){
	    dir.create(file.path(WD,"heatmaps"))
	  }
	}
	
	if(N_MARKERS>nrow(meth.data)){
	  logger.error("You cannot select more markers than you have sites")
	}
	
	for(group in groups){
		
		ind<-1:nrow(meth.data)
		
		if(MARKER_SELECTION[group]=="all"){
		  #Do nothing
		  ind <- ind
		}
		
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
  			
  			if(inherits(REF_DATA_SET,"RnBSet")){
  			  rnb.ref.set <- REF_DATA_SET
  			}else{
  			  rnb.ref.set <- load.rnb.set(REF_DATA_SET)
  			}
  			meth.data.ref <- meth(rnb.ref.set)
  			pheno.data.ref <- pheno(rnb.ref.set)
  			
  			marker.fstat<-rowFtests(meth.data.ref[ind,], as.factor(pheno.data.ref[[REF_PHENO_COLUMN]]))
  			
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
		
		if(MARKER_SELECTION[group] == "pcadapt"){
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
		  ind <- ind[order(lpval,decreasing = F)[1:N_MARKERS]]
		}
		
		if(MARKER_SELECTION[group]=="edec_stage0"){
		  require("EDec")
		  require("EDecExampleData")
		  if(is.null(REF_DATA_SET)){
		    markers <- run_edec_stage_0(reference_meth = EDecExampleData::reference_meth,
		                                reference_classes = EDecExampleData::reference_meth_class,
		                                max_p_value = 1e-5,
		                                num_markers = N_MARKERS)
		  }else{
		    if(inherits(REF_DATA_SET,"RnBSet")){
		      rnb.ref.set <- REF_DATA_SET
		    }else{
		      rnb.ref.set <- load.rnb.set(REF_DATA_SET)
		    }
		    if(!REF_PHENO_COLUMN %in% colnames(pheno(rnb.ref.set))){
		      stop("Supplied REF_PHENO_COLUMN not in phenotypic information")
		    }
		    markers <- run_edec_stage_0(reference_meth = meth(rnb.ref.set),
		                                reference_classes = pheno(rnb.ref.set)[,REF_PHENO_COLUMN],
		                                max_p_value = 1e-5,
		                                num_markers = N_MARKERS)
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
		
		if(MARKER_SELECTION[group]=="custom"){
			if(file.exists(CUSTOM_MARKER_FILE)){
			  if(grepl("RDS",CUSTOM_MARKER_FILE)){
				  custom_filter<-readRDS(CUSTOM_MARKER_FILE)
			  }else if(grepl("txt|csv|tsv",CUSTOM_MARKER_FILE)){
			    custom_filter <- as.numeric(readLines(CUSTOM_MARKER_FILE))
			  }else{
			    stop("Invalid value for CUSTOM_MARKER_FILE")
			  }
			  ind<-intersect(ind, custom_filter)
			}
		}
		
		if(store.heatmaps){
		  create.heatmap(meth.data[ind,],trait,sample.cols,palette,WD,MARKER_SELECTION[group])
		}
		
		if(WRITE_FILES){
			saveRDS(ind, file=sprintf("%s/cg_group_%d.RDS", WD, group))
		}
		
		cg_groups[[group]]<-ind
		
	}
	names(cg_groups) <- MARKER_SELECTION

	return(cg_groups=cg_groups)
}

#' create.heatmap
#' 
#' This function create a heatmap of methylation values with the colors specified according to the trait.
#' 
#' @param meth.data The methylation data to be plotted
#' @param trait A column in the phenotypic table specifying a grouping.
#' @param sample.cols A vector of colors specifying the grouping in \code{trait}.
#' @param WD The working directory, in which the plots are to be added in a subfolder "heatmaps".
#' @param palette The color palette from which \code{sample.cols} was obtained.
#' @param sel.type The method used for selecting subsets of sites.
#' @author Michael Scherer
#' @noRd
create.heatmap <- function(meth.data,
                           trait,
                           sample.cols,
                           palette,
                           WD,
                           sel.type){
  if(!is.null(sample.cols)){
    png(file.path(WD,"heatmaps",paste0(sel.type,".png")))
    heatmap.2(meth.data,
      trace="none",
      ColSideColors=sample.cols
    )
    legend(x=0,y=1,levels(trait),col=palette,pch=15)
    dev.off()
  }else{
    png(file.path(WD,"heatmaps",paste0(sel.type,".png")))
    heatmap.2(meth.data,
              trace="none"
    )
    dev.off()
  }
}
