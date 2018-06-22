#' prepare_CG_subsets
#' 
#' This routine selects a subset of CpGs sites used for MeDeCom analysis. Different selection methods are supported.
#' 
#' @param rnb.set An object of type \code{\link{RnBSet}} containing methylation, sample and optional coverage information.
#' @param MARKER_SELECTION A vector of strings representing marker selection methods. Available method are \itemize{
#'                                  \item pheno Selected are the top \code{N_MARKERS} site that differ between the phenotypic
#'                                         groups defined in data preparation or by \code{\link{rnb.sample.groups}}. Those are
#'                                         selected by employing limma on the methylation matrix.
#'                                  \item houseman2012 The 50k sites reported as cell-type specific in the Houseman's reference-
#'                                         based deconvolution.
#'                         }
prepare_CG_subsets<-function(
		rnb.set,
		MARKER_SELECTION,
		N_MARKERS=5000,
		WRITE_FILES=FALSE,
		WD=NA
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
		  if(file.exists(sprintf("%s/houseman.50k.markers.RDS", WD))){
			  houseman.50k.markers<-readRDS(sprintf("%s/houseman.50k.markers.RDS", WD))
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
			houseman.50k.markers<-readRDS(sprintf("%s/jaffe.irrizzary.markers.600.RDS", DD))
			ind<-intersect(ind, houseman.50k.markers)
		}
		
		if(all(grepl("rowFstat", MARKER_SELECTION[group])) && "REF_DATA_SET" %in% ls() && "REF_PHENO_COLUMN" %in% ls()){
			
			require(genefilter)
			
			load.env<-new.env(parent=emptyenv())
			
			load(file.path(DATA_DIR, REF_DATA_SET, "data.set.RData"), envir = load.env)
			load(file.path(DATA_DIR, REF_DATA_SET, "pheno.RData"), envir = load.env)
			
			meth.data.ref<-get("meth.data", envir = load.env)
			pheno.data.ref<-get("pheno.data", envir = load.env)
			
			marker.fstat<-rowFtests(meth.data.ref[ind,], as.factor(pheno.data.ref[[REF_PHENO_COLUMN]]))
			
			if(MARKER_SELECTION[group]=="rowFstat1k"){
				NTOP<-1000
			}else if(MARKER_SELECTION[group]=="rowFstat5k"){
				NTOP<-5000
			}else if(MARKER_SELECTION[group]=="rowFstat10k"){
				NTOP<-10000
			}else if(MARKER_SELECTION[group]=="rowFstat15k"){
				NTOP<-15000
			}else if(MARKER_SELECTION[group]=="rowFstat20k"){
				NTOP<-20000
			}
			
			subset<-which(rank(-marker.fstat$statistic)<=NTOP)
			ind<-ind[subset]						
		}
		
		if(MARKER_SELECTION[group]=="random1k"){
			subset<-sample.int(length(ind), min(1000, length(ind)))
			ind<-ind[subset]
		}
		
		if(MARKER_SELECTION[group]=="random5k"){
			subset<-sample.int(length(ind), min(5000, length(ind)))
			ind<-ind[subset]
		}
		
		if(MARKER_SELECTION[group]=="random10k"){
			subset<-sample.int(length(ind), min(10000, length(ind)))
			ind<-ind[subset]
		}
		
		if(MARKER_SELECTION[group]=="random25k"){
			subset<-sample.int(length(ind), min(25000, length(ind)))
			ind<-ind[subset]
		}
		
		if(MARKER_SELECTION[group]=="random50k"){
			subset<-sample.int(length(ind), min(50000, length(ind)))
			ind<-ind[subset]
		}
		
		if(!"N_PRIN_COMP" %in% ls()){
			N_PRIN_COMP<-10
		}
		
		if(MARKER_SELECTION[group]=="pca500"){
			
			meth.data.subs<-meth.data[ind,]
			pca<-prcomp(t(meth.data.subs))
			rot<-pca$rotation[,1:min(N_PRIN_COMP,ncol(pca$rotation))]
			
			pca.ind<-integer()
			for(cix in 1:min(N_PRIN_COMP,ncol(pca$rotation))){
				pca.ind<-union(pca.ind,order(abs(rot[,cix]), decreasing = TRUE)[1:500])
			}
			ind<-ind[sort(pca.ind)]	
		}
		
		
		if(MARKER_SELECTION[group]=="pca1k"){
			meth.data.subs<-meth.data[ind,]
			pca<-prcomp(t(meth.data.subs))
			rot<-pca$rotation[,1:min(N_PRIN_COMP,ncol(pca$rotation))]
			
			pca.ind<-integer()
			for(cix in 1:min(N_PRIN_COMP,ncol(pca$rotation))){
				pca.ind<-union(pca.ind,order(abs(rot[,cix]), decreasing = TRUE)[1:1000])
			}
			ind<-ind[sort(pca.ind)]	
		}
		
		if(MARKER_SELECTION[group]=="pca5k"){
			meth.data.subs<-meth.data[ind,]
			pca<-prcomp(t(meth.data.subs))
			rot<-pca$rotation[,1:min(N_PRIN_COMP,ncol(pca$rotation))]
			
			pca.ind<-integer()
			for(cix in 1:min(N_PRIN_COMP,ncol(pca$rotation))){
				pca.ind<-union(pca.ind,order(abs(rot[,cix]), decreasing = TRUE)[1:5000])
			}
			ind<-ind[sort(pca.ind)]	
		}
		
		if(MARKER_SELECTION[group]=="pca10k"){
			meth.data.subs<-meth.data[ind,]
			pca<-prcomp(t(meth.data.subs))
			rot<-pca$rotation[,1:min(N_PRIN_COMP,ncol(pca$rotation))]
			
			pca.ind<-integer()
			for(cix in 1:min(N_PRIN_COMP,ncol(pca$rotation))){
				pca.ind<-union(pca.ind,order(abs(rot[,cix]), decreasing = TRUE)[1:10000])
			}
			ind<-ind[sort(pca.ind)]	
		}
		
		if(grepl("var", MARKER_SELECTION[group])){
			
			sds<-apply(meth.data[ind,], 1, sd)
			
			if(MARKER_SELECTION[group] == "var1k"){
				ind<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),1000)]]
			}else if(MARKER_SELECTION[group] == "var5k"){
				ind<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),5000)]]
			}else if(MARKER_SELECTION[group] == "var10k"){
				ind<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),10000)]]
			}else if(MARKER_SELECTION[group] == "var50k"){
				ind<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),50000)]]
			}else if(MARKER_SELECTION[group] == "var100k"){
				ind<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),100000)]]
			}
		}
		
		if(grepl("hybrid", MARKER_SELECTION[group])){
			
			sds<-apply(meth.data[ind,], 1, sd)
			
			if(MARKER_SELECTION[group] == "hybrid1k"){
				var.set<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),500)]]
				random.set<-sample(setdiff(ind,var.set), 500)
			}else if(MARKER_SELECTION[group] == "hybrid5k"){
				var.set<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),2500)]]
				random.set<-sample(setdiff(ind,var.set), 2500)
			}else if(MARKER_SELECTION[group] == "hybrid10k"){
				var.set<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),5000)]]
				random.set<-sample(setdiff(ind,var.set), 5000)
			}else if(MARKER_SELECTION[group] == "hybrid50k"){
				var.set<-ind[order(sds, decreasing=TRUE)[1:min(length(ind),25000)]]
				random.set<-sample(setdiff(ind,var.set), 25000)
			}
			
			ind<-sort(c(var.set, random.set))
			rm(var.set)
			rm(random.set)
		}
		
		if(grepl("var", MARKER_SELECTION[group])){
			
			ranges<-apply(meth.data[ind,], 1, range)
			
			ranges<-ranges[2,]-ranges[1,]
			
			if(MARKER_SELECTION[group] == "range005"){
				ind<-ind[ranges>0.05]
			}else if(MARKER_SELECTION[group] == "range01"){
				ind<-ind[ranges>0.1]
			}else if(MARKER_SELECTION[group] == "range015"){
				ind<-ind[ranges>0.15]
			}else if(MARKER_SELECTION[group] == "range02"){
				ind<-ind[ranges>0.2]
			}
			rm(ranges)
		}
		
		if(MARKER_SELECTION[group]=="custom"){
			if(file.exists(sprintf("%s/%s",DD, CUSTOM_MARKER_FILE))){
				custom_filter<-readRDS(sprintf("%s/%s",DD,CUSTOM_MARKER_FILE))
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