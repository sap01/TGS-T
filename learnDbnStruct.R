# Infer Dynamic Bayesian Network (DBN) parents of each node
# If cardinality of 'nbgh' is 'm', then the complexity of this psudocode
# is O(n * 2^m), as mentioned in Sec. 'Analysis of LBN computational complexity'
# of the LBN paper.

#---------------------------------
# Begin: Loading the Packages
#---------------------------------
library(bnstruct)
library(ggm)
#---------------------------------
# End: Loading the Packages
#---------------------------------

#---------------------------------
# Begin: Step 2: Decomposing the network
#---------------------------------
# Input parameters:
## G_MI: Adjacency matrix of the mutual information network. Rownames and colnames should be node names.
## discr.data: Discretized data. ToDo: Specify required discretization levels. Rows are time points and columns are nodes.

learnDbnStruct <- function(discr.input.data, mi.net.adj.matrix)
{
  
  # load('mi.net.adj.matrix.Rdata')
  # load('discr.input.Rdata')
  
  # In 'di.net.adj.matrix', rows are src nodes and cols are tgt nodes
  di.net.adj.matrix <- mi.net.adj.matrix
  di.net.adj.matrix[1:nrow(di.net.adj.matrix), 1:ncol(di.net.adj.matrix)] <- 0
  
  num.nodes <- ncol(discr.input.data) # Num of nodes
  num.timepts <- nrow(discr.input.data) # Num of time points
  
  for(centralNodeIdx in 1:num.nodes) # for each central node
  {
    centralNodeIdx <- 1
    
    centralNodeName <- rownames(mi.net.adj.matrix)[centralNodeIdx]
    
    # if central node does not have any neighbour in mutual information net
    if (sum(mi.net.adj.matrix[, centralNodeIdx]) == 0)
    {
      next
    }
    
    # List names of the central node's neighbours in G_MI 
    nbghNames <- rownames(mi.net.adj.matrix[which(mi.net.adj.matrix[, centralNodeIdx] == 1),])
    
    local.net.node.names <- c(centralNodeName, nbghNames)
    
    local.discr.data <- discr.input.data[, local.net.node.names]
    
    
    # #---------------------------------
    # # Begin: Local BN struct learning
    # #---------------------------------
    # local.BN.input.data <- bnstruct::BNDataset(local.discr.data,
    #                                         variables = colnames(local.discr.data),
    #                                         discreteness = rep(TRUE, ncol(local.discr.data)),
    #                                         node.sizes = rep(2, ncol(local.discr.data)),
    #                                         starts.from = 0)
    # 
    # # Default params: scoring.func = "BDeu"
    # local.BN <- bnstruct::learn.network(local.discr.data, algo = 'sm')
    # plot(local.BN)
    # 
    # # Rows are src nodes and cols are tgt nodes
    # local.BN.adj.matrix <- bnstruct::dag(local.BN)
    # local.BN.adj.matrix <- matrix(local.BN.adj.matrix, nrow = length(local.BN@variables), 
    #                               ncol = length(local.BN@variables),
    #                               dimnames = c(list(local.BN@variables), list(local.BN@variables)))
    # 
    # local.BN.adj.matrix[, centralNodeName]
    # #---------------------------------
    # # End: Local BN struct learning
    # #---------------------------------
    
    
    #---------------------------------
    # Begin: Local Unrolled DBN struct learning
    #---------------------------------
    local.DBN.input.data <- local.discr.data[1, ]
    for (rowIdx in 2:nrow(local.discr.data))
    {
      local.DBN.input.data <- cbind(local.DBN.input.data, local.discr.data[rowIdx, ])
    }
    
    local.DBN.input.data.var.names <- c()
    for (time.pt in 1:nrow(local.discr.data)) {
      for (var.name in colnames(local.discr.data)) {
        local.DBN.input.data.var.names <- c(local.DBN.input.data.var.names,
                                            paste(var.name, as.character(time.pt), sep = "_t"))
      }
    }
    colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
    
    local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
                                                          discreteness = rep(TRUE, ncol(local.DBN.input.data)),                                            
                                                          variables = colnames(local.DBN.input.data),
                                                          node.sizes = rep(2, ncol(local.DBN.input.data)),
                                                          starts.from = 0,
                                                          num.time.steps = nrow(local.discr.data))
    
    # algo = "mmhc", scoring.func = "BDeu"
    local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset, 
                                                           num.time.steps = 
                                                             num.time.steps(local.DBN.input.data.BNDataset))
    
    
    # plot(local.DBN)
    
    # Extracting the adjacency matrix of the local DBN
    local.unrolled.DBN.adj.matrix <- bnstruct::dag(local.unrolled.DBN)
    local.unrolled.DBN.adj.matrix <- matrix(local.unrolled.DBN.adj.matrix, 
                                            nrow = length(local.unrolled.DBN@variables), 
                                            ncol = length(local.unrolled.DBN@variables),
                                            dimnames = c(list(local.unrolled.DBN@variables), list(local.unrolled.DBN@variables)))
    #---------------------------------
    # End: Local Unrolled DBN struct learning
    #---------------------------------
    
    #---------------------------------
    # Begin: Roll up local unrolled DBN struct
    #---------------------------------
    local.unrolled.DBN.centralNodeNames <- c()
    for (time.pt in 1:nrow(local.discr.data)) {
      local.unrolled.DBN.centralNodeNames <- c(local.unrolled.DBN.centralNodeNames,
                                               paste(centralNodeName, as.character(time.pt), sep = "_t"))
    }
    
    local.unrolled.DBN.adj.matrix.tgt.centralNode <- local.unrolled.DBN.adj.matrix[, local.unrolled.DBN.centralNodeNames]
    
    # If the value corr. to a row in 'local.unrolled.DBN.adj.matrix.tgt.centralNode.single.col'  
    # is greater than zero, then the node corr. to the row name is a parent of the central node
    local.unrolled.DBN.adj.matrix.tgt.centralNode.single.col <- 
      matrix(rowSums(local.unrolled.DBN.adj.matrix.tgt.centralNode), 
             nrow = nrow(local.unrolled.DBN.adj.matrix.tgt.centralNode), ncol = 1,
             dimnames = c(list(rownames(local.unrolled.DBN.adj.matrix.tgt.centralNode)),
                          centralNodeName))
    
    local.rolled.DBN.parents.Idx <- c()
    for (rowIdx in 1:nrow(local.unrolled.DBN.adj.matrix.tgt.centralNode.single.col))
    {
      if (local.unrolled.DBN.adj.matrix.tgt.centralNode.single.col[rowIdx, ] > 0
          & !((rowIdx %% length(local.net.node.names)) %in% local.rolled.DBN.parents.Idx))
      {
        local.rolled.DBN.parents.Idx <- c(local.rolled.DBN.parents.Idx, rowIdx)
      }
    }
    local.rolled.DBN.parents.Idx <- sort(local.rolled.DBN.parents.Idx)
    
    local.rolled.DBN.parents.names <- local.net.node.names[local.rolled.DBN.parents.Idx]
    
    # Remove self loop
    local.rolled.DBN.parents.names <- local.rolled.DBN.parents.names[
      !(local.rolled.DBN.parents.names == centralNodeName)]
    
    di.net.adj.matrix[local.rolled.DBN.parents.names, centralNodeName] <- 1
    #---------------------------------
    # End: Roll up local unrolled DBN struct
    #---------------------------------
    
  }
  
  # Return the final directed net of the current iteration
  return (di.net.adj.matrix)
}