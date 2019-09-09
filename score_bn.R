## This function is a wrapper function for calling
## 'all_fam_log_marg_lik_v2()' in 'src/score_bn.c'
##
ScoreBn <- function(local.dbn.input.data, node.sizes, scoring.func, 
                    tgt.node.idx, src.node.idx, num.src.nodes, 
                    init.path) {
  
  num.nodes <- ncol(local.dbn.input.data)
  levels <- rep(0, num.nodes)
  q.data <- quantize.matrix(local.dbn.input.data, levels)
  rm(num.nodes, levels)
  
  ## Since C index starts from zero
  tgt.node.idx <- (tgt.node.idx - 1)
  
  ## Since C index starts from zero.
  ## One would be deducted from each elt.
  src.node.idx <- (src.node.idx - 1)
  
  ## just to be sure
  storage.mode(node.sizes) <- "integer" 
  storage.mode(scoring.func) <- "integer"
  storage.mode(tgt.node.idx) <- "integer"
  storage.mode(src.node.idx) <- "integer"
  storage.mode(num.src.nodes) <- "integer"
  
  result <- .Call('all_fam_log_marg_lik_v2', q.data, node.sizes, 
                  scoring.func, tgt.node.idx, src.node.idx, 
                  num.src.nodes)
  return(result)
}

## This function is copied from
## https://github.com/cran/bnstruct/blob/663451121fd432fa30fe7cfbf6ca392f4dd4c8c2/R/util.R
##
quantize.matrix <- function(data, levels) {
  nr <- nrow(data)
  nc <- ncol(data)
  
  quant <- matrix(0,nr,nc)
  
  # print(levels)
  
  for( i in 1:nc )
  {
    if( levels[i] == 0 )  # already discrete
      quant[,i] <- as.matrix(data[,i],nr,1)
    else
    {
      quantiles <- quantile( data[,i], probs = (0:levels[i])/levels[i], na.rm = TRUE )
      # cut the range using the quantiles as break points.
      quant[,i] <- as.matrix( cut( data[,i], quantiles, labels=FALSE, include.lowest=TRUE),nr,1 )
    }
  }
  
  storage.mode(quant) <- "integer"
  # print(sapply(1:nc,function(x)max(quant[,x])))
  colnames(quant) <- colnames(data)
  return(quant)
}