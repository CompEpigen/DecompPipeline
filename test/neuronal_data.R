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



res<-prepare_data(
		RNB_SET=rnb.set, 
		WORK_DIR=file.path(PROJECT.DIR),
		DATASET="adgbbSorted",
		DATA_SUBSET="frontal",
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
		FILTER_BEADS=!is.null(rnb.set@covg.sites),
		FILTER_INTENSITY=inherits(rnb.set, "RnBeadRawSet"),
		FILTER_NA=TRUE,
		FILTER_CONTEXT=TRUE,
		FILTER_SNP=TRUE,
		FILTER_SOMATIC=TRUE,
		MATLAB_EXPORT=TRUE
)


cg_subsets<-prepare_CG_subsets(
		res$rnb.set.filtered,
		MARKER_SELECTION=c("var5k", "var10k")
)

md.res<-start_medecom_analysis(
		rnb.set=res$rnb.set.filtered,
		WORK_DIR=file.path(PROJECT.DIR),
		cg_groups=cg_subsets,
		Ks=2:10,
		LAMBDA_GRID=c(0,10^(-5:-1)),
		SAMPLE_SUBSET=NULL,
		K_FIXED=NULL,
		WRITE_FILES=TRUE,
		startT=NULL,
		startA=NULL,
		CLUSTER_SUBMIT=FALSE,
		CLUSTER_RDIR=NA,
		CLUSTER_HOSTLIST="*",
		CLUSTER_MEMLIMIT="5G",
#		MAX_JOBS=1000,
#		WAIT_TIME="30m",
#		PORTIONS=FALSE,
#		JOB_FILE=NA,
		CLEANUP=FALSE
)

