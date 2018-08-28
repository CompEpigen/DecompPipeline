start_medecom_analysis<-function(
        meth.data=NULL,
        rnb.set=NULL,
        WORK_DIR,
        cg_groups,
        Ks,
        LAMBDA_GRID,
        SAMPLE_SUBSET=NULL,
        K_FIXED=NA,
        WRITE_FILES=TRUE,
        startT=NULL,
        startA=NULL,
        A_LOWER=NULL,
        A_UPPER=NULL,
        NINIT=1, 
        ITERMAX=2, 
        NFOLDS=2,
        N_COMP_LAMBDA=4,
        NCORES=1,
        OPT_METHOD="cppTAfact",
        CLUSTER_SUBMIT=FALSE,
        CLUSTER_RDIR=NA,
        CLUSTER_HOSTLIST="*",
        CLUSTER_MEMLIMIT="5G",
#        MAX_JOBS=1000,
#        WAIT_TIME="30m",
#        PORTIONS=FALSE,
#        JOB_FILE=NA,
        CLEANUP=FALSE,
        analysis_info=NULL,
        LAMBDA_GRID_TYPE="standard",
        ANALYSIS_TOKEN="customAnalysis"
){
    library(MeDeCom)
    library(R.utils)
    
    if(!is.null(analysis_info)){
        ANALYSIS_ID<-paste(
            analysis_info$DATASET, 
            analysis_info$DATA_SUBSET,
            analysis_info$NORMALIZATION,
            analysis_info$QUALITY_FILTERING, 
            analysis_info$MARKER_SELECTION,
            LAMBDA_GRID_TYPE,
            ANALYSIS_TOKEN,
            sep="_")
    }else{
        ANALYSIS_ID<-"customAnalysis"
        analysis_info<-list()
    }
    analysis_info$ANALYSIS<-ANALYSIS_ID
    
    #RDIR="/TL/deep-share/archive00/software/bin"
    .libPaths(sprintf("%s/Rlib_test", WORK_DIR))
    
    if(is.na(CLUSTER_RDIR)){
        if(any(grepl("deep", R.utils:::System$getHostname()))){
            #RDIR="/TL/deep-share/archive00/software/packages/R/R-3.1.1/bin"
            CLUSTER_RDIR="/TL/deep-share/archive00/software/packages/R/R-devel_20160307/bin"
        }else if(any(grepl("t7600", R.utils:::System$getHostname()))){
            CLUSTER_RDIR="/usr/bin/"
        }else{
            CLUSTER_RDIR="/usr/bin/"
        }
    }
    
#    SRCDIR=sprintf("%s/projects/parameter_tuning/src", WORK_DIR)
#    DATA_DIR=sprintf("%s/projects/parameter_tuning/data", WORK_DIR)
#    GLOBAL_DD<-sprintf("%s/projects/parameter_tuning/data/common", WORK_DIR)
    
    DD<-file.path(WORK_DIR, "data")
    
    if(TRUE){
        print("Did not write the variable dump: should only be executed from an environment with all the variables set")    
    }else{
        var_list<-c(ANALYSIS_ID)
        system(sprintf("mkdir %s", WORK_DIR))
        dump(ls()[var_list], file=file.path(WORK_DIR, "analysis_settings.RDump"))
    }
    
    if(is.null(meth.data)){
        if(is.null(rnb.set)){
            load(sprintf("%s/data.set.RData", DD))
        }else{
            meth.data<-meth(rnb.set)
        }
    }
#    if(!DATASET %in% c("sim", "simReal")){
#        load(sprintf("%s/indices.RData", GLOBAL_DD))
#        load(sprintf("%s/snp.filter.RData", GLOBAL_DD))
#        load(sprintf("%s/somatic.filter.RData", GLOBAL_DD))
#    }
    
#    if(file.exists(sprintf("%s/quality.filter.RData", DD))){
#        load(sprintf("%s/quality.filter.RData", DD))
#    }
    
#    if(file.exists(sprintf("%s/Mint.RDS", DD))){
#        M.raw<-readRDS(sprintf("%s/Mint.RDS", DD))
#    }else{
#        M.raw<-NULL
#    }
    
#    if(file.exists(sprintf("%s/Uint.RDS", DD))){
#        U.raw<-readRDS(sprintf("%s/Uint.RDS", DD))
#    }else{
#        U.raw<-NULL
#    }
    
#    if(file.exists(sprintf("%s/Nbeads.RDS", DD))){
#        b.raw<-readRDS(sprintf("%s/Nbeads.RDS", DD))
#    }else{
#        b.raw<-NULL
#    }
    
    if(WRITE_FILES){
        saveRDS(LAMBDA_GRID, file=sprintf("%s/lambda_grid.RDS", WORK_DIR))
        if(!is.null(SAMPLE_SUBSET)){
            saveRDS(SAMPLE_SUBSET, file=sprintf("%s/sample_subset.RDS", WORK_DIR))
        }else{
            SAMPLE_SUBSET<-1:ncol(meth.data)
        }
    }
    
    groups<-1:length(cg_groups)
    
#    if(QUALITY_FILTERING %in% c("standard","MplusU") && !is.null(M.raw) && !is.null(U.raw)){
#        
#        if(!"MIN_INT_QUANT" %in% ls()){
#            MIN_INT_QUANT<-0.1
#        }
#        if(!"MAX_INT_QUANT" %in% ls()){
#            MAX_INT_QUANT<-0.95
#        }
#        if(!"MIN_N_BEADS" %in% ls()){
#            MIN_N_BEADS<-3
#        }
#        
#        MplusU<-M.raw+U.raw
#        
#        hm450_ann<-readRDS(file.path(GLOBAL_DD, "hm450_cg_annot.RDS"))
#        
#        MplusU.I<-MplusU[hm450_ann$Design=="I",]
#        MplusU.II<-MplusU[hm450_ann$Design=="II",]
#        
#        MU.q001.I<-sort(as.numeric(MplusU.I))[ceiling(MIN_INT_QUANT*nrow(MplusU.I)*ncol(MplusU.I))]
#        MU.q099.I<-sort(as.numeric(MplusU.I))[ceiling(MAX_INT_QUANT*nrow(MplusU.I)*ncol(MplusU.I))]
#        
#        MU.q001.II<-sort(as.numeric(MplusU.II))[ceiling(MIN_INT_QUANT*nrow(MplusU.II)*ncol(MplusU.II))]
#        MU.q099.II<-sort(as.numeric(MplusU.II))[ceiling(MAX_INT_QUANT*nrow(MplusU.II)*ncol(MplusU.II))]
#        
#        MplusU.f<-matrix(FALSE, nrow=nrow(MplusU), ncol=ncol(MplusU))
#        
#        MplusU.f[hm450_ann$Design=="I",]<-MplusU.I>MU.q001.I & MplusU.I<MU.q099.I
#        MplusU.f[hm450_ann$Design=="II",]<-MplusU.II>MU.q001.II & MplusU.II<MU.q099.II
#        
#        qf.MU<-which(rowSums(MplusU.f)==ncol(MplusU.f))
#        
#        if(!is.null(b.raw)){
#            qf.b<-which(rowSums(b.raw>=MIN_N_BEADS)==ncol(b.raw))
#            qf.MU.beads<-intersect(qf.MU, qf.b)
#        }else{
#            qf.MU.beads<-qf.MU
#        }
#        
#    }
    
#    for(group in groups){
#        
#        # LOAD THE RESULTS
#        #if(DATASET %in% c("sim", "simReal")){
#        ind<-1:nrow(meth.data)    
#        #}else{
#        #    ind<-sort(unlist(indices[GROUP_LISTS[[group]]]))
#        #}
#        
#        if(QUALITY_FILTERING=="standard"){
#            ind<-intersect(ind, qf.MU.beads)
#            ind<-intersect(ind, snp.filter)
#            ind<-intersect(ind, somatic.filter)
#        }else if(QUALITY_FILTERING=="snpANDsomatic"){
#            ind<-intersect(ind, snp.filter)
#            ind<-intersect(ind, somatic.filter)
#        }else if(QUALITY_FILTERING=="somatic"){
#            ind<-intersect(ind, somatic.filter)
#        }else if(QUALITY_FILTERING=="snp"){
#            ind<-intersect(ind, snp.filter)
#        }else if(QUALITY_FILTERING=="MplusU"){
#            ind<-intersect(ind, qf.MU.beads)
#        }
#        
#    }
    if(!is.null(A_LOWER) && is.null(A_UPPER)){
        saveRDS(A_LOWER, file=sprintf("%s/A_lower.RDS", WORK_DIR))
        saveRDS(A_UPPER, file=sprintf("%s/A_upper.RDS", WORK_DIR))
    }
    
    if(!is.na(K_FIXED)){
        saveRDS(K_FIXED, file=sprintf("%s/fixed_T_cols.RDS", WORK_DIR))
    }else{
        K_FIXED<-NULL
    }
    
#    if("START" %in% ls() && file.exists(START)){
#        system(sprintf("cp %s %s/start.RData", START, WORK_DIR))
#        load(file.path(WORK_DIR, START))
#        startT=result$T
#        startA=result$A
#    }else{
#        startT=NULL
#        startA=NULL
#    }
    
    if(file.exists(sprintf("%s/trueT.RData", DD))){
        load(sprintf("%s/trueT.RData", DD))
    }else{
        trueT=NULL
    }
    if(file.exists(sprintf("%s/trueA.RData", DD))){
        load(sprintf("%s/trueA.RData", DD))
    }else{
        trueA=NULL
    }
    
#    if(PORTIONS && is.na(JOB_FILE)){
#        JOB_FILE<-"/tmp/job_file"
#    }
#    cluster_submit<-function(qsub_string, portions=PORTIONS, job_script_file=JOB_FILE){
#        
#        if(portions){
#            cat(qsub_string, file=JOB_FILE, append=TRUE, sep="\n")
#        }else{
#            system(qsub_string)
#            #print(qsub_string)
#        }
#        
#    }
    
    if(CLUSTER_SUBMIT){
        cluster.settings=list(R_bin_dir=CLUSTER_RDIR, host_pattern=CLUSTER_HOSTLIST, mem_limit=CLUSTER_MEMLIMIT)
    }else{
        cluster.settings=NULL
    }
    
    result<-runMeDeCom(data=meth.data, 
            Ks=Ks,
            lambdas=LAMBDA_GRID,
            opt.method=sprintf("MeDeCom.%s",OPT_METHOD),
            cg_subsets=cg_groups,
            sample_subset=SAMPLE_SUBSET,
            startT=startT,
            startA=startA,
            trueT=trueT,
            trueA=trueA,
            fixed_T_cols=K_FIXED,
            NINIT=NINIT, 
            ITERMAX=ITERMAX, 
            NFOLDS=NFOLDS,
            N_COMP_LAMBDA=4,
            NCORES=NCORES,
            analysis.name=ANALYSIS_ID,
            use.ff=FALSE,
            cluster.settings=cluster.settings,
            temp.dir=WORK_DIR,
            cleanup=CLEANUP,
            verbosity=1L,
            time.stamps=TRUE
            #,random.seed=RANDOM_SEED
    )
    
    ### TODO: Fix this once
    
#    result@parameters$ANALYSIS <- ANALYSIS_ID
#    result@parameters$NORMALIZATION <- NORMALIZATION
    result@parameters$ITERMAX<-ITERMAX
#    result@parameters$MARKER_SELECTION<- MARKER_SELECTION
    result@parameters$NFOLDS<-NFOLDS
#    result@parameters$ANALYSIS_TOKEN<-""
    result@parameters$NINIT<-NINIT
#    result@parameters$DATASET<-DATASET 
#    result@parameters$DATA_SUBSET<-DATA_SUBSET 
    
    result@parameters$cg_subsets <- c(1:length(cg_groups))
    
    analysis_info$GROUP_LISTS<-cg_groups
    analysis_info$ANALYSIS_DATE<-date()
    result@dataset_info<-c(result@dataset_info, analysis_info)
    
    if(WRITE_FILES){
        saveRDS(result, file=file.path(WORK_DIR, "collected.result.RDS"))
    }
    
    return(result)
}




