#
#          Neuronal data from Gasparoni et al., 2017         
#
.libPaths("/DEEP_fhgfs/projects/plutsik/Rlib_clean_RnBeads/")
PROJECT.DIR="/DEEP_fhgfs/projects/plutsik/projects/neuron/"
data<-read.table(file.path(PROJECT.DIR, "20170419_GasSort_GuiSort_manuscriptData.txt"))
pd<-read.table(file.path(PROJECT.DIR, "20170419_GasSort_GuiSort_manuscriptSampleSheet.txt"), sep='\t', header=TRUE)


probe.list<-rownames(data)

library(RnBeads)
rnb.options(disk.dump.big.matrices=FALSE)
rnb.set<-RnBeadSet(pd, probe.list, as.matrix(data))

#### neuron

res<-prepare_data(
		RNB_SET=rnb.set, 
		WORK_DIR=file.path(PROJECT.DIR),
		DATASET="adgbbSorted",
		DATA_SUBSET="neuron",
		SAMPLE_SELECTION_COL="Celltype",
		SAMPLE_SELECTION_GREP="Neuron",
		PHENO_COLUMNS=NA,
		ID_COLUMN=NA,
		NORMALIZATION="custom",
		REF_CT_COLUMN=NA,
		REF_RNB_SET=NA,
		REF_RNB_CT_COLUMN=NA,
		PREPARE_TRUE_PROPORTIONS=FALSE,
		TRUE_A_TOKEN=NA,
		HOUSEMAN_A_TOKEN=NA,
		ESTIMATE_HOUSEMAN_PROP=FALSE,
		FILTER_BEADS=!is.null(rnb.set@covg.sites),
		FILTER_INTENSITY=inherits(rnb.set, "RnBeadRawSet"),
		FILTER_NA=TRUE,
		FILTER_CONTEXT=TRUE,
		FILTER_SNP=TRUE,
		FILTER_SOMATIC=TRUE,
		MATLAB_EXPORT=TRUE
)


res_selection<-prepare_CG_subsets(
		res$rnb.set.filtered,
		MARKER_SELECTION=c("var10k", "var50k"),
		WD=file.path(PROJECT.DIR),
		analysis_info=res$info
)

md.res<-start_medecom_analysis(
		rnb.set=res$rnb.set.filtered,
		WORK_DIR=file.path(PROJECT.DIR),
		cg_groups=res_selection$cg_groups,
		Ks=2:10,
		LAMBDA_GRID=c(0,10^(-5:-1)),
		SAMPLE_SUBSET=NULL,
		K_FIXED=NA,
		WRITE_FILES=TRUE,
		startT=NULL,
		startA=NULL,
		NINIT=10, 
		ITERMAX=300, 
		NFOLDS=9,
		N_COMP_LAMBDA=4,
		NCORES=10,
		OPT_METHOD="cppTAfact",
		CLUSTER_SUBMIT=FALSE,
		CLUSTER_RDIR=NA,
		CLUSTER_HOSTLIST="*",
		CLUSTER_MEMLIMIT="5G",
#		MAX_JOBS=1000,
#		WAIT_TIME="30m",
#		PORTIONS=FALSE,
#		JOB_FILE=NA,
		CLEANUP=FALSE,
		analysis_info=res_selection$info
)


saveRDS(md.res, file="/DEEP_fhgfs/projects/plutsik/projects/neuron/md.res.neuron.var10k.var50k_new.RDS")


