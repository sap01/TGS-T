## Goal: Learn local Dynamic Bayesian Network (DBN)
##
#######################################################################################
## Goal: Learn local DBN using R package 'bnstruct'.
## This function allows nodes with less than two discrete levels.
#######################################################################################
LearnLocalDbnBnstruct <- function(local.dbn.input.data, 
                                  node.sizes, 
                                  scoring.func, 
                                  init.path) {
  
  ## Number of shortlisted source nodes.
  ## The last col is for the target node.
  num.sl.src.nodes <- (ncol(local.dbn.input.data) - 1)
  
  ########################################################################
  ## Begin: Find out the subset of source nodes with the best score
  ########################################################################
  
  ## Begin with the empty subset
  # curr.subset <- c()
  curr.subset.str <- as.list(rep(FALSE, num.sl.src.nodes))
  
  # print(local.dbn.input.data)
  # print(str(local.dbn.input.data))
  
  ## Target node index
  tgt.node.idx <- ncol(local.dbn.input.data)
  
  ########################################################################  
  ## Begin: Calc score of the empty subset
  ########################################################################
  
  ## Source node indices.
  ## Initialized with the empty subset.
  src.node.idx <- c()
  
  ## No. of source nodes
  num.src.nodes <- 0
  
  ## Initialize the current score
  curr.score <- NULL
  
  if (scoring.func == 'BIC') {
    
    ## score_bn.R
    ## 'scoring.func': 0 = BDeu, 1 = AIC, 2 = BIC
    ##
    curr.score <- ScoreBn(local.dbn.input.data, node.sizes, 
                          scoring.func = 2, 
                          tgt.node.idx, 
                          src.node.idx, num.src.nodes)
  } else {
    stop('Only supports BIC scoring function till now')
  }
  ########################################################################
  ## End: Calc score of the empty subset
  ########################################################################
  
  best.score <- curr.score
  best.subset.str <- curr.subset.str
  
  ########################################################################  
  ## Begin: Calc scores of the non-empty subsets
  ########################################################################
  
  num.subsets <- (2^num.sl.src.nodes)
  
  for (subset.idx in 2:num.subsets) {
    
    ## 'find.next.subset.str()' is defined in this script
    curr.subset.str <- find.next.subset.str(curr.subset.str)
    ##> 'curr.subset.str' = list(TRUE, FALSE, TRUE)
    
    src.node.idx <- which(curr.subset.str == TRUE, 
                          arr.ind = TRUE, 
                          useNames = FALSE)
    
    num.src.nodes <- length(src.node.idx)
    
    if (scoring.func == 'BIC') {
      
      ## score_bn.R
      ## 'scoring.func': 0 = BDeu, 1 = AIC, 2 = BIC
      ##
      curr.score <- ScoreBn(local.dbn.input.data, node.sizes, 
                            scoring.func = 2, 
                            tgt.node.idx, 
                            src.node.idx, num.src.nodes)
    } else {
      stop('Only supports BIC scoring function till now')
    }
    
    ## The higher the score, the fitter the model.
    ## Ref: Lines 38-42, defn of fbp() in 'src/smfast.c' of
    ## R package 'bnstruct' version 1.0.2.
    if (curr.score > best.score) {
      best.score <- curr.score
      best.subset.str <- curr.subset.str
    }
    
    # print(tgt.node.idx)
    # print(src.node.idx)
    # print(best.score)
    # print(curr.score)
    
  }
  
  ## Remove all large objects that are no more required
  rm(subset.idx, num.subsets, curr.score, curr.subset.str, best.score, 
     src.node.idx, tgt.node.idx, node.sizes)
  ########################################################################
  ## End: Calc scores of the non-empty subsets
  ########################################################################    
  
  ########################################################################
  ## End: Find out the subset of source nodes with the best score
  ########################################################################
  
  ## Shortlisted source node names
  ## E.g., [1] "v1_t1" "v2_t1" "v3_t1"
  sl.src.node.names <- colnames(local.dbn.input.data)[1:num.sl.src.nodes]
  
  ## Predicted source node names  
  pred.src.node.names <- sl.src.node.names[unlist(best.subset.str)]
  # print(pred.src.node.names)
  rm(sl.src.node.names, best.subset.str)
  pred.src.node.names <- as.list(pred.src.node.names)
  
  ## Return the list of predicted source node names    
  return(pred.src.node.names)  
}

#######################################################################################
## Goal: Learn local DBN using R package 'bnlearn'.
## This function requires that each node has at least two discrete levels.
#######################################################################################
LearnLocalDbnBnlearn <- function(local.dbn.input.data, 
                               scoring.func) {
  
  ## Number of nodes in the local DBN
  num.nodes.local.dbn <- ncol(local.dbn.input.data)
  
  ## E.g., [1] "v1_t2" 
  tgt.node.name <- colnames(local.dbn.input.data)[ncol(local.dbn.input.data)]
  
  ## Number of shortlisted source nodes.
  ## The last col is for the target node.
  num.sl.src.nodes <- (ncol(local.dbn.input.data) - 1) 
  
  ## E.g., [1] "v1_t1" "v2_t1" "v3_t1"
  sl.src.node.names <- colnames(local.dbn.input.data)[1:num.sl.src.nodes]
  
  ########################################################################
  ## Begin: Find out the subset of source nodes with the best score
  ########################################################################
  
  ## Begin with the empty subset
  # curr.subset <- c()
  curr.subset.str <- as.list(rep(FALSE, num.sl.src.nodes))
  
    ########################################################################
    ## Begin: Compose model string for the source nodes
    ########################################################################
    src.model.str <- c()
    
    for (node.name in sl.src.node.names) {
      node.to.add <- paste('[', 
                           node.name, 
                           ']', 
                           sep = '')
      
      src.model.str <- paste(src.model.str, 
                              node.to.add, 
                              sep = '')
    }
    rm(node.name)
    ##> 'src.model.str' = '[v1_t1][v2_t1][v3_t1]'
    ########################################################################
    ## End: Compose model string for the source nodes
    ########################################################################
  
  ## Model string for the target node
  tgt.model.str <- paste('[', 
                         tgt.node.name, 
                         ']', 
                         sep = '')
  ##> 'tgt.model.str' = '[v1_t2]'
  
  ## Current model string
  curr.model.str <- paste(src.model.str, 
                          tgt.model.str, 
                          sep = '')
  ##> 'curr.model.str' = '[v1_t1][v2_t1][v3_t1][v1_t2]'
  
  ## Current network model
  curr.net <- bnlearn::model2network(curr.model.str)
  
  ## Initialize the current score
  curr.score <- NULL
  
  local.dbn.input.data <- as.data.frame(local.dbn.input.data)
  
  ## Convert each column from int to factor
  for (col.idx in 1:ncol(local.dbn.input.data)) {
    local.dbn.input.data[, col.idx] <- as.factor(local.dbn.input.data[, col.idx])
  }
  rm(col.idx)
  
  print(local.dbn.input.data)
  print(str(local.dbn.input.data))
  
    ########################################################################  
    ## Begin: Calc score of the empty subset
    ########################################################################
    if (scoring.func == 'BIC') {
      curr.score <- bnlearn::score(curr.net, 
                                   local.dbn.input.data, 
                                   type = 'bic')
    } else {
      stop('Only supports BIC scoring function till now')
    }
    ########################################################################
    ## End: Calc score of the empty subset
    ########################################################################
  
  best.score <- curr.score
  best.subset.str <- curr.subset.str
  
    ########################################################################  
    ## Begin: Calc scores of the non-empty subsets
    ########################################################################
  
    num.subsets <- (2^num.sl.src.nodes)
  
    for (subset.idx in 2:num.subsets) {
      
      ## 'find.next.subset.str()' is defined in this script
      curr.subset.str <- find.next.subset.str(curr.subset.str)
      ##> 'curr.subset.str' = list(TRUE, FALSE, TRUE)
      
      ## Model string for the target node
      tgt.model.str <- paste('[', 
                             tgt.node.name, 
                             '|', 
                             sep = '')
      ##> 'tgt.model.str' = '[v1_t2|'
      
      is.first <- TRUE
      
      for (node.idx in length(curr.subset.str)) {
        if (curr.subset.str[[node.idx]]) {
          
          if (is.first) {
            is.first <- FALSE
            
          } else {
            tgt.model.str <- paste(tgt.model.str, 
                                   ':', 
                                   sep = '')
          }
          
          tgt.model.str <- paste(tgt.model.str, 
                                 sl.src.node.names[node.idx], 
                                 sep = '')
        }
      }
      rm(node.idx)
      
      tgt.model.str <- paste(tgt.model.str, 
                             ']',
                             sep = '')
      ##> tgt.model.str = '[v1_t2|v1_t1:v1_t3]'
      
      ## Current model string
      curr.model.str <- paste(src.model.str, 
                              tgt.model.str, 
                              sep = '')
      ##> 'curr.model.str' = '[v1_t1][v2_t1][v3_t1][v1_t2|V1_t1:v1_t3]'
      
      ## Current network model
      curr.net <- bnlearn::model2network(curr.model.str)
      
      if (scoring.func == 'BIC') {
        curr.score <- bnlearn::score(curr.net, 
                                     local.dbn.input.data, 
                                     type = 'bic')
      } else {
        stop('Only supports BIC scoring function till now')
      }
      
      ## The higher the score, the fitter the model.
      ## Ref: Section 'Note' of function 'score' in the manual of 
      ## R package 'bnlearn' version 4.3.
      if (curr.score > best.score) {
        best.score <- curr.score
        best.subset.str <- curr.subset.str
      }
      
    }
    rm(subset.idx, num.subsets)
    ########################################################################
    ## End: Calc scores of the non-empty subsets
    ########################################################################    
  
  ########################################################################
  ## End: Find out the subset of source nodes with the best score
  ########################################################################
  
  ## Predicted source node names  
  pred.src.node.names <- sl.src.node.names[unlist(best.subset.str)]
  rm(sl.src.node.names)
  pred.src.node.names <- as.list(pred.src.node.names)
  
  ## Return the list of predicted source node names    
  return(pred.src.node.names)  
}

#######################################################################################
## Goal: Find next subset string, given the previous subset string.
## Example 1:
## If 'prev.subset.str' = list(FALSE, FALSE)
## then it returns 'list(FALSE, TRUE)'.
##
## Example 2:
## If 'prev.subset.str' = list(FALSE, TRUE)
## then it returns 'list(TRUE, FALSE)'.
##
## So this function basically adds a TRUE to the last element of the
## list 'prev.subset.str'.
#######################################################################################
find.next.subset.str <- function(prev.subset.str) {
  
  num.nodes.local.dbn <- length(prev.subset.str)
  
  ## Initialize the carry bit
  carry.bit <- TRUE
  
  for (node.idx in num.nodes.local.dbn:1) {
    
    ## Cuurent bit
    curr.bit <- prev.subset.str[[node.idx]]
    
    if (!curr.bit & !carry.bit) {
      prev.subset.str[[node.idx]] <- FALSE
      carry.bit <- FALSE
    } else if (curr.bit & carry.bit) {
      prev.subset.str[[node.idx]] <- FALSE
      carry.bit <- TRUE
    } else {
      prev.subset.str[[node.idx]] <- TRUE
      carry.bit <- FALSE
    }
  }
  rm(node.idx)
  
  ## It is an in-place replacement i.e.
  ## the content of the previous subset
  ## string is modified to represent
  ## the next subset string
  return(prev.subset.str)
}