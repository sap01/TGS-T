# Goal: Infer Dynamic Bayesian Network (DBN)

#---------------------------------
# Begin: Loading the Packages
#---------------------------------
library(bnstruct)
library(ggm)
library(foreach)
library(doParallel)
#---------------------------------
# End: Loading the Packages
#---------------------------------

#---------------------------------
# Begin: Step 2: Decomposing the network
#---------------------------------
# Input parameters:
## input.data.3D: Dimensions {1 = time points, 2 = variables, 3 = samples under the same time point}. 
## mi.net.adj.matrix: Adjacency matrix of the mutual information network. Rownames and colnames should be node names.
## num.discr.levels: If input data is discretized, then number of discrete levels for each variable.
## Else if input data is continuous, then number of levels in which data needs to be discretized
## for performing the DBN structure learning.
## 
test_learnDbnStruct3dPar <- function(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels, num.nodes, num.timepts)
{
  # load('dream3.yeast1.size10.trajectory.3D.Rdata')
  # input.data.3D <- dream3.yeast1.size10.trajectory.3D.data  
  # load('mi.net.adj.matrix.Rdata')

  # In 'di.net.adj.matrix', rows are src nodes and cols are tgt nodes
  di.net.adj.matrix <- mi.net.adj.matrix
  di.net.adj.matrix[1:nrow(di.net.adj.matrix), 1:ncol(di.net.adj.matrix)] <- 0
  
  no_cores <- num.nodes
  cl <- parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(cl)

  # For each central node, learn local DBN struct in parallel.
  # Package names, specified by '.packages' must be copied to each worker core.
  # Current environmental variables that are referenced inside the foreach loop are copied
  # to each worker automatically.
  di.net.adj.matrix <- foreach::foreach(centralNodeIdx = 1:num.nodes, .combine = 'cbind', .packages = 'bnstruct')  %dopar%
  {
  #   # centralNodeIdx <- 1
  #   
    centralNodeName <- rownames(mi.net.adj.matrix)[centralNodeIdx]
  #   
  #   # if central node does not have any neighbour in mutual information net
  #   if (sum(mi.net.adj.matrix[, centralNodeIdx]) == 0)
  #   {
  #     next
  #   }
  #   
  #   # List names of the central node's neighbours in mi.net.adj.matrix
  #   nbghNames <- c()
  #   if (sum(mi.net.adj.matrix[, centralNodeIdx]) == 1) # Just 1 neighbour
  #   {
  #     for (nbrIdx in 1:n)
  #     {
  #       if (mi.net.adj.matrix[nbrIdx, centralNodeIdx] == 1)
  #       {
  #         nbghNames <- rownames(mi.net.adj.matrix)[nbrIdx]
  #         break
  #       }
  #     }
  #   }
  #   else if (sum(mi.net.adj.matrix[, centralNodeIdx]) > 1) # Multiple neighbours
  #   {
  #     nbghNames <- rownames(mi.net.adj.matrix[which(mi.net.adj.matrix[, centralNodeIdx] == 1),])
  #   }
  #   
  #   local.net.node.names <- c(centralNodeName, nbghNames)
  #   
  #   local.input.data.3D <- input.data.discr.3D[, local.net.node.names, ]
  #   
  #   
  #   # #---------------------------------
  #   # # Begin: Local BN struct learning
  #   # #---------------------------------
  #   # local.BN.input.data <- bnstruct::BNDataset(local.discr.data,
  #   #                                         variables = colnames(local.discr.data),
  #   #                                         discreteness = rep(TRUE, ncol(local.discr.data)),
  #   #                                         node.sizes = rep(2, ncol(local.discr.data)),
  #   #                                         starts.from = 0)
  #   # 
  #   # # Default params: scoring.func = "BDeu"
  #   # local.BN <- bnstruct::learn.network(local.discr.data, algo = 'sm')
  #   # plot(local.BN)
  #   # 
  #   # # Rows are src nodes and cols are tgt nodes
  #   # local.BN.adj.matrix <- bnstruct::dag(local.BN)
  #   # local.BN.adj.matrix <- matrix(local.BN.adj.matrix, nrow = length(local.BN@variables), 
  #   #                               ncol = length(local.BN@variables),
  #   #                               dimnames = c(list(local.BN@variables), list(local.BN@variables)))
  #   # 
  #   # local.BN.adj.matrix[, centralNodeName]
  #   # #---------------------------------
  #   # # End: Local BN struct learning
  #   # #---------------------------------
  #   
  #   
  #   #---------------------------------
  #   # Begin: Local Unrolled DBN struct learning
  #   #---------------------------------
  #   
  #   sampleIdx <- 1
  #   
  #   local.DBN.input.data <- matrix(local.input.data.3D[1, , sampleIdx], nrow = 1, 
  #                                  ncol = dim(local.input.data.3D)[2])
  #   
  #   for (time.pt in 2:dim(local.input.data.3D)[1])
  #   {
  #     data.to.combine <- matrix(local.input.data.3D[time.pt, , sampleIdx], nrow = 1, 
  #                               ncol = dim(local.input.data.3D)[2])
  #     local.DBN.input.data <- cbind(local.DBN.input.data, data.to.combine)
  #   }
  #   
  #   local.DBN.input.data.var.names <- c()
  #   for (time.pt in 1:dim(local.input.data.3D)[1]) {
  #     for (var.name in dimnames(local.input.data.3D)[2]) {
  #       local.DBN.input.data.var.names <- c(local.DBN.input.data.var.names,
  #                                           paste(var.name, as.character(time.pt), sep = "_t"))
  #     }
  #   }
  #   colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
  #   
  #   if (dim(local.input.data.3D)[3] > 1) # If there are multiple samples per time pt
  #   {
  #     for (sampleIdx in 2:dim(local.input.data.3D)[3])
  #     {
  #       sample.to.combine <- matrix(local.input.data.3D[1, , sampleIdx], nrow = 1, 
  #                                   ncol = dim(local.input.data.3D)[2])
  #       
  #       for (time.pt in 2:dim(local.input.data.3D)[1])
  #       {
  #         data.to.combine <- matrix(local.input.data.3D[time.pt, , sampleIdx], nrow = 1, 
  #                                   ncol = dim(local.input.data.3D)[2])
  #         sample.to.combine <- cbind(sample.to.combine, data.to.combine)
  #       }
  #       
  #       # local.DBN.input.data.var.names <- c()
  #       # for (time.pt in 1:dim(local.input.data.3D)[1]) {
  #       #   for (var.name in dimnames(local.input.data.3D)[2]) {
  #       #     local.DBN.input.data.var.names <- c(local.DBN.input.data.var.names,
  #       #                                         paste(var.name, as.character(time.pt), sep = "_t"))
  #       #   }
  #       # }
  #       # colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
  #       
  #       local.DBN.input.data <- rbind(local.DBN.input.data, sample.to.combine)
  #     }
  #   }
  #   
  #   rownames(local.DBN.input.data) <- as.vector(unlist(dimnames(local.input.data.3D)[3]))
  #   
  #   # local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
  #   #                                         discreteness = rep(TRUE, ncol(local.DBN.input.data)),                                            
  #   #                                         variables = colnames(local.DBN.input.data),
  #   #                                         node.sizes = rep(num.discr.levels, ncol(local.DBN.input.data)),
  #   #                                         starts.from = 0,
  #   #                                         num.time.steps = nrow(local.input.data.3D))
  #   
  #   local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
  #                                                         discreteness = rep(TRUE, ncol(local.DBN.input.data)),                                            
  #                                                         variables = colnames(local.DBN.input.data),
  #                                                         node.sizes = rep(num.discr.levels, ncol(local.DBN.input.data)),
  #                                                         num.time.steps = nrow(local.input.data.3D))
  #   
  #   # algo = "mmhc", scoring.func = "BDeu"
  #   local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset, 
  #                                                          num.time.steps = 
  #                                                            bnstruct::num.time.steps(local.DBN.input.data.BNDataset))
  #   
  #   # # The following four lines must be executed at the same time
  #   # save.plot.to.filename = paste(paste('LocalUnrolledDbn', centralNodeName, sep = '_'), '.jpg', sep = '')
  #   # jpeg(file = paste('LocalUnrolledDbn_', centralNodeName, '.jpg', sep = ''))
  #   # plot(local.unrolled.DBN)
  #   # dev.off()
  #   
  #   # Extracting the adjacency matrix of the local DBN
  #   local.unrolled.DBN.adj.matrix <- bnstruct::dag(local.unrolled.DBN)
  #   local.unrolled.DBN.adj.matrix <- matrix(local.unrolled.DBN.adj.matrix, 
  #                                           nrow = length(local.unrolled.DBN@variables), 
  #                                           ncol = length(local.unrolled.DBN@variables),
  #                                           dimnames = c(list(local.unrolled.DBN@variables), list(local.unrolled.DBN@variables)))
  #   #---------------------------------
  #   # End: Local Unrolled DBN struct learning
  #   #---------------------------------
  #   
  #   #---------------------------------
  #   # Begin: Roll up local unrolled DBN struct
  #   #---------------------------------
  #   local.unrolled.DBN.centralNodeNames <- c()
  #   for (time.pt in 1:nrow(local.input.data.3D)) {
  #     local.unrolled.DBN.centralNodeNames <- c(local.unrolled.DBN.centralNodeNames,
  #                                              paste(centralNodeName, as.character(time.pt), sep = "_t"))
  #   }
  #   
  #   local.unrolled.DBN.adj.matrix.tgt.centralNode <- local.unrolled.DBN.adj.matrix[, local.unrolled.DBN.centralNodeNames]
  #   
  #   # If the value corr. to a row in 'local.unrolled.DBN.adj.matrix.tgt.centralNode.single.col'  
  #   # is greater than zero, then the node corr. to the row name is a parent of the central node
  #   local.unrolled.DBN.adj.matrix.tgt.centralNode.single.col <- 
  #     matrix(rowSums(local.unrolled.DBN.adj.matrix.tgt.centralNode), 
  #            nrow = nrow(local.unrolled.DBN.adj.matrix.tgt.centralNode), ncol = 1,
  #            dimnames = c(list(rownames(local.unrolled.DBN.adj.matrix.tgt.centralNode)),
  #                         centralNodeName))
  #   
  #   local.rolled.DBN.parents.Idx <- c()
  #   for (rowIdx in 1:nrow(local.unrolled.DBN.adj.matrix.tgt.centralNode.single.col))
  #   {
  #     if (local.unrolled.DBN.adj.matrix.tgt.centralNode.single.col[rowIdx, ] > 0
  #         & !((rowIdx %% length(local.net.node.names)) %in% local.rolled.DBN.parents.Idx))
  #     {
  #       local.rolled.DBN.parents.Idx <- c(local.rolled.DBN.parents.Idx, (rowIdx %% length(local.net.node.names)))
  #     }
  #   }
  #   local.rolled.DBN.parents.Idx <- sort(local.rolled.DBN.parents.Idx)
  #   
  #   local.rolled.DBN.parents.names <- local.net.node.names[local.rolled.DBN.parents.Idx]
  #   
  #   # Remove self loop
  #   local.rolled.DBN.parents.names <- local.rolled.DBN.parents.names[
  #     !(local.rolled.DBN.parents.names == centralNodeName)]
    
    di.net.adj.matrix[1:centralNodeIdx, centralNodeIdx] <- 1
    local.rolled.DBN.parents.names <- colnames(mi.net.adj.matrix)[1:centralNodeIdx]
    writeLines('\n')
    print(paste('centralNodeName', centralNodeName,
                'has parents', paste(local.rolled.DBN.parents.names, collapse = ',')))
    writeLines('\n')
    #---------------------------------
    # End: Roll up local unrolled DBN struct
    #---------------------------------
    
    # Force each core to sleep for 300 secs
    Sys.sleep(300)
    
    # Return value for each 'foreach' iteration 
    di.net.adj.matrix[, centralNodeName]
  }
  
  parallel::stopCluster(cl)
  
  dimnames(di.net.adj.matrix) <- dimnames(mi.net.adj.matrix)
  
  # Return the final directed net of the current iteration
  return (di.net.adj.matrix)
}