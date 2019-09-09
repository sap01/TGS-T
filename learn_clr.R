#######################################################################################################
## Goal: Learn CLR net.
## Ref: https://github.com/cran/minet/blob/662cd27b4a7bf0067237a0c698451754de623419/R/netinf.R
##
LearnClr <- function(mim, clr.algo) {
  var.id<-NULL
  if(is.data.frame(mim)) {
    var.id <- names(mim)
    mim <- as.matrix(mim)
  }
  else if( is.matrix(mim) ) 
    var.id <- names(as.data.frame(mim))
  else stop("Supply a matrix-like argument")      
  if(ncol(mim)!=nrow(mim))
    stop("Argument matrix must be square")
  
  if (clr.algo == 'CLR1.2') {
    # res <- .Call( "clr", mim, nrow(mim), DUP=FALSE,PACKAGE="minet" )
    res <- .Call( "clr", mim, nrow(mim), DUP=FALSE)
  } else if (clr.algo == 'CLR1.3') {
    # res <- .Call( "clr1dot3", mim, nrow(mim), DUP=FALSE,PACKAGE="minet" )
    res <- .Call( "clr1dot3", mim, nrow(mim), DUP=FALSE)
  }
  
  dim(res) <- dim(mim)
  res <- as.matrix(res)
  rownames(res) <- var.id
  colnames(res) <- var.id
  res              
}
#######################################################################################################