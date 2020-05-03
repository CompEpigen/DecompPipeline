#' bigFF.row.apply
#' 
#' This routine applies a function to chunks of the dataset that is stored in disk
#' either with `ff` or `BigFfMat`.
#' 
#' @param mat A disk-based matrix of type \code{\link{ff}} or \code{\link{BigFfMat}}.
#' @param FUN The function to be applied to each chunk of the matrix. Should be a function
#'             that computes matrix statistics, such as \code{rowMeans}.
#' @param iter.count Number of chunks to be created. The larger this number, the slower
#'             the computation, but the lower the disk usage.
#' @param ... Further arguments passed to FUN.
#' @return A vector summarizing the results of the function. Final structure is determined
#'          by FUN.
#' @author Michael Scherer
#' @noRd
bigFF.row.apply <- function(mat,FUN,iter.count=1000,...){
  chunk.size <- floor(nrow(mat)/iter.count)
  iter <- 1
  res <- c()
  while(iter+chunk.size < nrow(mat)){
    chunk <- mat[iter:(iter+chunk.size-1),]
    res <- c(res,FUN(chunk,...))
    iter <- iter+chunk.size
#    print(paste(round((iter/nrow(mat))*100,2)," percent completed"))
  }
  chunk <- mat[iter:nrow(mat),]
  res <- c(res,FUN(chunk,...))
  return(res)
}

#' create.heatmap
#' 
#' This function create a heatmap of methylation values with the colors specified according to the trait.
#' 
#' @param meth.data The methylation data to be plotted
#' @param trait A column in the phenotypic table specifying a grouping.
#' @param sample.cols A vector of colors specifying the grouping in \code{trait}.
#' @param out.dir The working directory, in which the plots are to be added in a subfolder "heatmaps".
#' @param palette The color palette from which \code{sample.cols} was obtained.
#' @param sel.type The method used for selecting subsets of sites.
#' @author Michael Scherer
#' @noRd
create.heatmap <- function(meth.data,
                           trait,
                           sample.cols,
                           palette,
                           out.dir,
                           sel.type){
  if(!is.null(sample.cols)){
    png(file.path(out.dir,"heatmaps",paste0(sel.type,".png")))
    heatmap.2(meth.data,
      trace="none",
      ColSideColors=sample.cols
    )
    legend(x=0,y=1,levels(trait),col=palette,pch=15)
    dev.off()
  }else{
    png(file.path(out.dir,"heatmaps",paste0(sel.type,".png")))
    heatmap.2(meth.data,
              trace="none"
    )
    dev.off()
  }
}

#' A small RnBeads object used to run the examples.
#' @name rnb.set.example
#' @docType data
#' @author Michael Scherer
NULL
