#
# Data preparation script for the parameter tuning run
#

prepare_data<-function(
		RNB_SET, 
		WORK_DIR,
		DATASET,
		DATA_SUBSET,
		SAMPLE_SELECTION_COL=NA,
		SAMPLE_SELECTION_GREP=NA,
		PHENO_COLUMNS=NA,
		ID_COLUMN=NA,
		NORMALIZATION="none",
		REF_CT_COLUMN=NA,
		REF_RNB_SET=NA,
		REF_RNB_CT_COLUMN=NA,
		PREPARE_TRUE_PROPORTIONS=FALSE,
		TRUE_A_TOKEN=NA,
		HOUSEMAN_A_TOKEN=NA,
		ESTIMATE_HOUSEMAN_PROP=FALSE,
		FILTER_BEADS=!is.null(RNB_SET@covg.sites),
		MIN_N_BEADS=3,
		FILTER_INTENSITY=inherits(RNB_SET, "RnBeadRawSet"),
		MIN_INT_QUANT=0.1,
		MAX_INT_QUANT=0.95,
		FILTER_NA=TRUE,
		FILTER_CONTEXT=TRUE,
		FILTER_SNP=TRUE,
		FILTER_SOMATIC=TRUE,
		MATLAB_EXPORT=TRUE
){
	#suppressPackageStartupMessages(library(MeDeCom))
	suppressPackageStartupMessages(require(RnBeads))
	suppressPackageStartupMessages(require(R.matlab))
	
	DATADIR<-file.path(WORK_DIR, "data")
	#DATADIR<-"/home/lutsik/Documents/science/projects/heterogeneity/data/parameter_tuning/test_output"
	#DATADIR<-"/mnt/deepfhfs/projects/parameter_tuning/data/"
	#ANALYSIS<-paste(DATASET, DATA_SUBSET, NORMALIZATION, sep="_")
	#
	
	OUTPUTDIR<-sprintf("%s/%s_%s_%s", DATADIR, DATASET, DATA_SUBSET, NORMALIZATION)
	system(sprintf("mkdir %s", OUTPUTDIR))
	
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
	if(!(NORMALIZATION %in% c("none","custom"))){
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
	
	if(!is.na(ID_COLUMN)){
		sample_ids<-pd[subs,ID_COLUMN]
		saveRDS(sample_ids, file=sprintf("%s/sample_ids.RDS", OUTPUTDIR))	
	}
	
	if(!is.na(REF_RNB_SET) && !is.na(REF_RNB_CT_COLUMN)){
		
		rnb.set.ref<-load.rnb.set(REF_RNB_SET)
		meth.rnb.ref<-meth(rnb.set.ref)
		pd.ref<-pheno(rnb.set.ref)
		
		ct<-pd.ref[[REF_RNB_CT_COLUMN]]
		nnas<-!is.na(ct)
		ct<-ct[nnas]
		
		meth.ref.ct<-lapply(unique(ct), function(cn) rowMeans(meth.rnb.ref[,nnas][,ct==cn]))
		trueT<-do.call("cbind", meth.ref.ct)
		colnames(trueT)<-unique(ct)
		
		save(trueT, file=sprintf("%s/trueT.RData", OUTPUTDIR))
		
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
		#meth.rheuma<-meth(rnb.set.comb,row.names=TRUE)[,!is.na(rnb.set.comb@pheno$diseaseState)]
		#pd.rheuma<-pheno(rnb.set.comb)[!is.na(rnb.set.comb@pheno$diseaseState),]
		
		if(!is.na(TRUE_A_TOKEN)){
			
			trueA<-t(na.omit(pd[subs,grep(TRUE_A_TOKEN, colnames(pd))]))
			rownames(trueA)<-gsub(TRUE_A_TOKEN, "", colnames(pd)[grep(TRUE_A_TOKEN, colnames(pd))])
			
			save(trueA, file=sprintf("%s/trueA.RData", OUTPUTDIR))
			
		}
		
		if(!is.na(HOUSEMAN_A_TOKEN)){
			
			Ahouseman2012<-t(na.omit(pd[,grep(HOUSEMAN_A_TOKEN, colnames(pd))]))
			rownames(Ahouseman2012)<-NULL
			
			save(Ahouseman2012, file=sprintf("%s/Ahouseman2012.RData", OUTPUTDIR))
			
		}else if(ESTIMATE_HOUSEMAN_PROP){
			
			if("REF_RNB_SET" %in% ls() && !is.na(REF_RNB_CT_COLUMN)){
				rnb.set.ref<-remove.samples(rnb.set.ref, which(is.na(pheno(rnb.set.ref)[,REF_RNB_CT_COLUMN])))
				rnb.set<-combine(rnb.set, rnb.set.ref)
				REF_CT_COLUMN<-REF_RNB_CT_COLUMN
			}
			print("Estimating proportions using the Houseman et al, 2012 method")
			res<-estimateProportionsCP(rnb.set, REF_CT_COLUMN, NA, 2000, full.output = TRUE)
			
			Ahouseman2012<-t(res$contributions.nonneg)
			Ahouseman2012[Ahouseman2012<1e-5]<-0
			
			Ahouseman2012<-sweep(Ahouseman2012, 2, colSums(Ahouseman2012),"/")
			
			save(Ahouseman2012, file=sprintf("%s/Ahouseman2012.RData", OUTPUTDIR))
		}
	}
	
	if(!is.na(REF_CT_COLUMN)){
		meth.data<-meth.rnb[,subs]
	}else{
		meth.data<-meth.rnb
	}
	colnames(meth.data)<-NULL
	
	save(meth.data, file=sprintf("%s/data.set.RData", OUTPUTDIR))
	
	
	####################### FILTERING
	
	get.probe.ind<-function(rnb.set){
		
		ann<-annotation(rnb.set)
		full.ann<-rnb.annotation2data.frame(rnb.get.annotation(rnb.set@target))
		ind<-match(ann$ID, full.ann$ID)
		return(ind)
		
	}
	
	filter.annotation<-function(
			rnb.set,
			snp=TRUE,
			somatic=TRUE,
			context=TRUE)
	{
		#probe.ind<-get.probe.ind(rnb.set)
		
		annot<-annotation(rnb.set)
		
		probe.ind.filtered<-1:nrow(annot)
		
		if(snp){
			snp.filter<-which(!annot$`SNPs 3` & !annot$`SNPs 5` & !annot$`SNPs Full`)
			probe.ind.filtered<-intersect(probe.ind.filtered, snp.filter)
		}
		
		if(somatic){
			somatic.filter<-which(!annot$Chromosome %in% c("chrX", "chrY"))
			probe.ind.filtered<-intersect(probe.ind.filtered, somatic.filter)
		}
		
		if(context && inherits(rnb.set, "RnBeadSet")){
			context.filter<-grep("cg", rownames(annot))
			probe.ind.filtered<-intersect(probe.ind.filtered, context.filter)
		}
		return(probe.ind.filtered)
	}
	
	
	filter.quality<-function(
			rnb.set,
			probe.ind,
			beads=TRUE,
			min.beads=3,
			intensity=TRUE,
			na=TRUE
			){
		
		annot<-rnb.annnotation2data.frame(annotation(rnb.set))
		
		qf<-1:nrow(annot)
		
		if(beads){
			b.raw<-RnBeads:::covg(rnb.set, row.names=TRUE)
			qf.b<-which(rowSums(b.raw>=min.beads)==ncol(b.raw))
			qf<-intersect(qf, qf.b)
		}
		
		if(intensity){
			
			MplusU<-M.raw+U.raw
			
			hm450_ann<-readRDS(file.path(GLOBAL_DD, "hm450_cg_annot.RDS"))
			
			MplusU.I<-MplusU[hm450_ann$Design=="I",]
			MplusU.II<-MplusU[hm450_ann$Design=="II",]
			
			MU.q001.I<-sort(as.numeric(MplusU.I))[ceiling(MIN_INT_QUANT*nrow(MplusU.I)*ncol(MplusU.I))]
			MU.q099.I<-sort(as.numeric(MplusU.I))[ceiling(MAX_INT_QUANT*nrow(MplusU.I)*ncol(MplusU.I))]
			
			MU.q001.II<-sort(as.numeric(MplusU.II))[ceiling(MIN_INT_QUANT*nrow(MplusU.II)*ncol(MplusU.II))]
			MU.q099.II<-sort(as.numeric(MplusU.II))[ceiling(MAX_INT_QUANT*nrow(MplusU.II)*ncol(MplusU.II))]
			
			MplusU.f<-matrix(FALSE, nrow=nrow(MplusU), ncol=ncol(MplusU))
			
			MplusU.f[hm450_ann$Design=="I",]<-MplusU.I>MU.q001.I & MplusU.I<MU.q099.I
			MplusU.f[hm450_ann$Design=="II",]<-MplusU.II>MU.q001.II & MplusU.II<MU.q099.II
			
			qf.MU<-which(rowSums(MplusU.f)==ncol(MplusU.f))
			
			
			
			qf<-intersect(qf, qf.MU)
			
		}
		
		if(na){
			
			##### filter for missing and non-CG probes
			
			na.filter<-which(rowSums(is.na(meth.data))<1)
			qf<-intersect(qf, na.filter)
		}
		
		return(qf)
	}
	
	################################# QUALITY FILTERING ######################################
	FILTER_QUALITY<- FILTER_BEADS && FILTER_INTENSITY && FILTER_NA
	
	if(FILTER_QUALITY){
		M.raw<-RnBeads:::M(rnb.set, row.names=TRUE)
		U.raw<-RnBeads:::U(rnb.set, row.names=TRUE)
		
		if("REF_CT_COLUMN" %in% ls()){
			
			M.raw<-M.raw[,subs,drop=FALSE]
			U.raw<-U.raw[,subs,drop=FALSE]
			
		}
		
		saveRDS(M.raw, file.path(OUTPUTDIR, "Mint.RDS"))
		saveRDS(U.raw, file.path(OUTPUTDIR, "Uint.RDS"))
		
		b.raw<-RnBeads:::covg(rnb.set, row.names=TRUE)
		
		if("REF_CT_COLUMN" %in% ls()){
			b.raw<-b.raw[,subs, drop=FALSE]
		}
		
		saveRDS(b.raw, file.path(OUTPUTDIR, "Nbeads.RDS"))
		
		
		#M<-M[rownames(meth.data),]
		#U<-U[rownames(meth.data),]
		
		if(MATLAB_EXPORT){
			#writeMat(con=sprintf("%s/M.mat", OUTPUTDIR), M=M.raw)
			#writeMat(con=sprintf("%s/U.mat", OUTPUTDIR), U=U.raw)
			writeMat(con=sprintf("%s/MandU.mat", OUTPUTDIR), M=M.raw, U=U.raw)
		}
		
		qual.filter<-filter.quality(rnb.set,beads=FILTER_BEADS, min.beads=MIN_N_BEADS, intensity = FILTER_INTENSITY, na=FILTER_NA )
		
		save(qual.filter, file=sprintf("%s/quality.filter.RData", OUTPUTDIR))
		
		if(MATLAB_EXPORT){
			writeMat(con=sprintf("%s/quality.filter.mat", OUTPUTDIR), qualityFilterMUbeads=qf.MU.beads)
		}
		
	}else{
		qual.filter<-1:nrow(rnb.set@meth.sites)
	}
	########################################## ANNOTATION FILTERING ###################################################
	FILTER_ANNOTATION<-FILTER_CONTEXT && FILTER_SNP && FILTER_SOMATIC
	
	if(FILTER_ANNOTATION){
		
		annot.filter<-filter.annotation(rnb.set, context = FILTER_CONTEXT, snp = FILTER_SNP, somatic = FILTER_SOMATIC)
	
		save(annot.filter, file=sprintf("%s/annotation.filter.RData", OUTPUTDIR))
	
		if(MATLAB_EXPORT){
			writeMat(con=sprintf("%s/annotation.filter.mat", OUTPUTDIR), annotFilter=annot.filter)
		}
	}else{
		annot.filter<-1:nrow(rnb.set@meth.sites)
	}
	
	total.filter<-intersect(qual.filter, annot.filter)
	rnb.set.f<-remove.sites(rnb.set, setdiff(1:nrow(rnb.set@meth.sites), intersect(qual.filter, annot.filter)))
	
	analysis_info<-list()
	
	analysis_info$DATASET<-DATASET 
	analysis_info$DATA_SUBSET<-DATA_SUBSET
	analysis_info$QUALITY_FILTERING <- sprintf("%s%s%s%s%s%s",
			ifelse(FILTER_BEADS, "Beads", ""),
			ifelse(FILTER_INTENSITY, "Intensity", ""),
			ifelse(FILTER_NA, "Missing", ""),
			ifelse(FILTER_NA, "Context", ""),
			ifelse(FILTER_NA, "SNP", ""),
			ifelse(FILTER_NA, "Somatic", "")
		)
	
	analysis_info$NORMALIZATION <- NORMALIZATION
	
	return(list(quality.filter=qual.filter, annot.filter=annot.filter, total.filter=total.filter, rnb.set.filtered=rnb.set.f, info=analysis_info))
}


