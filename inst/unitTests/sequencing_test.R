library(devtools)
devtools::install_local("/TL/deep/projects/work/mscherer/projects/MeDeCom/DecompPipeline/",lib="/TL/deep/projects/work/mscherer/projects/RnBeads/R/cleaner_lib/")
  library(DecompPipeline)

options(fftempdir="/DEEP_fhgfs/projects/mscherer/deep/tmp")  
  rnb.set <- load.rnb.set("/DEEP_fhgfs/projects/mscherer/data/RRBS/mouse/brahamMouse2017/report/rnbSet_unnormalized/")
  work.dir <- "/TL/deep/projects/work/mscherer/projects/MeDeCom/test/prepare_BS"
  analysis.name <- "TestBS"
  
  data.prep <- prepare_data_BS(rnb.set,
                               WORK_DIR=work.dir,
                               analysis.name=analysis.name)

subsets <- prepare_CG_subsets(
	data.prep$rnb.set.filtered,
	MARKER_SELECTION=c("pheno","random","pca","var","hybrid","range"),
		WD=file.path(work.dir,"data","foo_foo_none"),
		N_MARKERS = 4242,
		N_PRIN_COMP = 2,
		RANGE_DIFF = 0.1
)

