## Goal: Compute Conditional Mutual Info between
## random variables 'v1' and 'v2' given
## random vector 'vcs'.
## Therefore, when 'vcs' is NULL, this 
## function returns the mutual info between
## random variables 'v1' and 'v2'.
## The implementation of this R function
## is adapted from the mutual info estimator
## function used in the PCA-CMI algo [1].
## The original funciton, namely 'cmi()' [2], 
## is written in MATLAB.
## 
##------------------------------------------------------------
## Begin: References
##------------------------------------------------------------
## [1] Zhang X, Zhao X M, He K, et al. Inferring gene 
## regulatory networks from gene expression data by 
## path consistency algorithm based on conditional 
## mutual information[J]. Bioinformatics, 2012, 28(1): 98-104.
## Companion website:
## https://sites.google.com/site/xiujunzhangcsb/software/pca-cmi
##
## [2] Function declaration: 'function cmiv=cmi(v1,v2,vcs)', 
## Source code file: http://www.comp-sysbio.org/grn/pca_cmi.m .
##------------------------------------------------------------
## End: References
##------------------------------------------------------------
##
## Compute Conditional Mutual Infortion (CMI) the way it is 
## done in the implementaion of the PCA-CMI algo
ComputeCmiPcaCmi <- function(v1,v2,vcs) {
  cmiv <- 0
  
  if(nargs() == 2) {
    c1 <-  var(v1)
    c2 <-  var(v2)
    v <-cbind(v1,v2)
    c3 <- det(cov(v))
    
    if (c3 != 0) {
      cmiv <- 0.5*log(c1*c2/c3) 
    } 
    # else {
    #   ## Perfectly correlated
    #   cmiv <- .Machine$double.xmax
    # }
    
  } else if(nargs() == 3) {   
    v <-cbind(v1,vcs)
    c1 <- det(cov(v))
    v <-cbind(v2,vcs)
    c2 <- det(cov(v))
    
    if (ncol(vcs) > 1) {
      c3 <- det(cov(vcs))
    } else {
      c3  <- var(vcs)
    }
    
    v <- cbind(v1,v2,vcs)
    c4 <- det(cov(v))
    
    if ((c3*c4) != 0) {
      cmiv <- 0.5 * log(abs((c1 * c2) / (c3 * c4)))
    }
  }
  return (cmiv)
}