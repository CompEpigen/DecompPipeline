#
#          Neuronal data from Gasparoni et al., 2017         
#
#.libPaths("/DEEP_fhgfs/projects/plutsik/Rlib_clean_RnBeads/")
DATA.DIR="/DEEP_fhgfs/projects/plutsik/projects/neuron/"
PROJECT.DIR = "/TL/deep/projects/work/mscherer/projects/MeDeCom/test/"
#data<-read.table(file.path(DATA.DIR, "20170419_GasSort_GuiSort_manuscriptData.txt"))
#pd<-read.table(file.path(DATA.DIR, "20170419_GasSort_GuiSort_manuscriptSampleSheet.txt"), sep='\t', header=TRUE)


#probe.list<-rownames(data)

library(DecompPipeline)
library(RnBeads)
#rnb.options(disk.dump.big.matrices=FALSE)
#rnb.set<-RnBeadSet(pd, probe.list, as.matrix(data))

rnb.set <- load.rnb.set("/TL/deep/projects/nobackup/mage/data/publicationData/processed/TCGA_OV___AH/rnbeads_report/rnbSet_unnormalized/")

res<-prepare_data(
		RNB_SET=rnb.set, 
		WORK_DIR=file.path(PROJECT.DIR),
		SAMPLE_SELECTION_COL=NA,
		SAMPLE_SELECTION_GREP=NA,
		PHENO_COLUMNS="bcr",
		ID_COLUMN=NULL,
		NORMALIZATION="none",
		REF_CT_COLUMN=NA,
		REF_RNB_SET="/DEEP_fhgfs/projects/mscherer/data/450K/Reinius_Blood_Reference_unnormalized.zip",
		REF_RNB_CT_COLUMN="tissue/cell type",
		PREPARE_TRUE_PROPORTIONS=FALSE,
		TRUE_A_TOKEN=NA,
		HOUSEMAN_A_TOKEN=NA,
		ESTIMATE_HOUSEMAN_PROP=FALSE,
		FILTER_BEADS=!is.null(rnb.set@covg.sites),
		FILTER_INTENSITY=inherits(rnb.set, "RnBeadRawSet"),
		FILTER_NA=T,
		FILTER_CONTEXT=TRUE,
		FILTER_SNP=TRUE,
		FILTER_SOMATIC=TRUE,
		FILTER_CROSS_REACTIVE=T,
		snp.list="/DEEP_fhgfs/projects/mscherer/data/EPIC/Radar_Genetik/commonSNPs137.txt"
)


cg_subsets<-prepare_CG_subsets(
		rnb.set=res$rnb.set.filtered,
		MARKER_SELECTION=c("pheno","houseman2012","houseman2014","jaffe2014","rowFstat","random","pca","var","hybrid","range","custom","pcadapt","all"),
		WD=file.path(PROJECT.DIR,"data","foo_foo_none"),
		N_MARKERS = 4242,
		REF_DATA_SET = "/DEEP_fhgfs/projects/mscherer/data/450K/Reinius_Blood_Reference_unnormalized.zip",
		REF_PHENO_COLUMN = "tissue/cell type",
		N_PRIN_COMP = 2,
		RANGE_DIFF = 0.1,
		CUSTOM_MARKER_FILE = "/DEEP_fhgfs/projects/mscherer/data/Aussois/MeDeCom/marker_file_Sophie.txt",
		K.prior=2
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

