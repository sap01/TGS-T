# Goal: Infer Dynamic Bayesian Network (DBN)

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
## input.data.3D: Dimensions {1 = time points, 2 = variables, 3 = samples under the same time point}. 
## mi.net.adj.matrix: Adjacency matrix of the mutual information network. Rownames and colnames should be node names.
## num.discr.levels: If input data is discretized, then number of discrete levels for each variable.
## Else if input data is continuous, then number of levels in which data needs to be discretized
## for performing the DBN structure learning.
## 
learnDbnStruct3dSer <- function(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels, num.nodes, num.timepts)
{
  # load('dream3.yeast1.size10.trajectory.3D.Rdata')
  # input.data.3D <- dream3.yeast1.size10.trajectory.3D.data  
  # load('mi.net.adj.matrix.Rdata')

  # In 'di.net.adj.matrix', rows are src nodes and cols are tgt nodes
  di.net.adj.matrix <- mi.net.adj.matrix
  di.net.adj.matrix[1:nrow(di.net.adj.matrix), 1:ncol(di.net.adj.matrix)] <- 0
  
  for(centralNodeIdx in 1:num.nodes) # for each central node
  {
    # centralNodeIdx <- 1
    
    centralNodeName <- rownames(mi.net.adj.matrix)[centralNodeIdx]
    
    # if central node does not have any neighbour in mutual information net
    if (sum(mi.net.adj.matrix[, centralNodeIdx]) == 0)
    {
      next
    }
    
    # List names of the central node's neighbours in mi.net.adj.matrix
    nbghNames <- c()
    if (sum(mi.net.adj.matrix[, centralNodeIdx]) == 1) # Just 1 neighbour
    {
      for (nbrIdx in 1:n)
      {
        if (mi.net.adj.matrix[nbrIdx, centralNodeIdx] == 1)
        {
          nbghNames <- rownames(mi.net.adj.matrix)[nbrIdx]
          break
        }
      }
    }
    else if (sum(mi.net.adj.matrix[, centralNodeIdx]) > 1) # Multiple neighbours
    {
      nbghNames <- rownames(mi.net.adj.matrix[which(mi.net.adj.matrix[, centralNodeIdx] == 1),])
    }
    
    local.net.node.names <- c(centralNodeName, nbghNames)
    
    local.input.data.3D <- input.data.discr.3D[, local.net.node.names, ]
    
    
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
    
    sampleIdx <- 1
    
    local.DBN.input.data <- matrix(local.input.data.3D[1, , sampleIdx], nrow = 1, 
                                   ncol = dim(local.input.data.3D)[2])
    
    for (time.pt in 2:dim(local.input.data.3D)[1])
    {
      data.to.combine <- matrix(local.input.data.3D[time.pt, , sampleIdx], nrow = 1, 
                                ncol = dim(local.input.data.3D)[2])
      local.DBN.input.data <- cbind(local.DBN.input.data, data.to.combine)
    }
    
    local.DBN.input.data.var.names <- c()
    for (time.pt in 1:dim(local.input.data.3D)[1]) {
      for (var.name in dimnames(local.input.data.3D)[2]) {
        local.DBN.input.data.var.names <- c(local.DBN.input.data.var.names,
                                            paste(var.name, as.character(time.pt), sep = "_t"))
      }
    }
    colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
    
    if (dim(local.input.data.3D)[3] > 1) # If there are multiple samples per time pt
    {
      for (sampleIdx in 2:dim(local.input.data.3D)[3])
      {
        sample.to.combine <- matrix(local.input.data.3D[1, , sampleIdx], nrow = 1, 
                                    ncol = dim(local.input.data.3D)[2])
        
        for (time.pt in 2:dim(local.input.data.3D)[1])
        {
          data.to.combine <- matrix(local.input.data.3D[time.pt, , sampleIdx], nrow = 1, 
                                    ncol = dim(local.input.data.3D)[2])
          sample.to.combine <- cbind(sample.to.combine, data.to.combine)
        }
        
        # local.DBN.input.data.var.names <- c()
        # for (time.pt in 1:dim(local.input.data.3D)[1]) {
        #   for (var.name in dimnames(local.input.data.3D)[2]) {
        #     local.DBN.input.data.var.names <- c(local.DBN.input.data.var.names,
        #                                         paste(var.name, as.character(time.pt), sep = "_t"))
        #   }
        # }
        # colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
        
        local.DBN.input.data <- rbind(local.DBN.input.data, sample.to.combine)
      }
    }
    
    rownames(local.DBN.input.data) <- as.vector(unlist(dimnames(local.input.data.3D)[3]))
    
    # local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
    #                                         discreteness = rep(TRUE, ncol(local.DBN.input.data)),                                            
    #                                         variables = colnames(local.DBN.input.data),
    #                                         node.sizes = rep(num.discr.levels, ncol(local.DBN.input.data)),
    #                                         starts.from = 0,
    #                                         num.time.steps = nrow(local.input.data.3D))
    
    local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
                                                          discreteness = rep(TRUE, ncol(local.DBN.input.data)),                                            
                                                          variables = colnames(local.DBN.input.data),
                                                          node.sizes = rep(num.discr.levels, ncol(local.DBN.input.data)),
                                                          num.time.steps = nrow(local.input.data.3D))
    
    # algo = "mmhc", scoring.func = "BDeu"
    local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset, 
                                                           num.time.steps = 
                                                             num.time.steps(local.DBN.input.data.BNDataset))
    
    # # The following four lines must be executed at the same time
    # save.plot.to.filename = paste(paste('LocalUnrolledDbn', centralNodeName, sep = '_'), '.jpg', sep = '')
    # jpeg(file = paste('LocalUnrolledDbn_', centralNodeName, '.jpg', sep = ''))
    # plot(local.unrolled.DBN)
    # dev.off()
    
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
    for (time.pt in 1:nrow(local.input.data.3D)) {
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
        local.rolled.DBN.parents.Idx <- c(local.rolled.DBN.parents.Idx, (rowIdx %% length(local.net.node.names)))
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

###############################################################################################################################
# Input parameters:
## input.data.3D: Dimensions {1 = time points, 2 = variables, 3 = samples under the same time point}. 
## mi.net.adj.matrix: Adjacency matrix of the mutual information network. Rownames and colnames should be node names.
## num.discr.levels: If input data is discretized, then number of discrete levels for each variable.
## Else if input data is continuous, then number of levels in which data needs to be discretized
## for performing the DBN structure learning.
## Uses layering.
## 
learnDbnStructLayer3dSer <- function(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels, num.nodes, num.timepts)
{
  #---------------------------------
  # Begin: Loading the Packages
  #---------------------------------
  library(bnstruct)
  # library(ggm)
  library(foreach)
  # library(doParallel)
  #---------------------------------
  # End: Loading the Packages
  #---------------------------------
  
  num.time.trans <- (num.timepts - 1)
  
  # '%do%' implies serial computing.
  # '.verbose = TRUE' is used for debugging.
  # '.packages' is nor required with '%do%'. If given, it is ignored. So it is not given.
  # 'when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0' means when central node does not have any neighbour in the mutual info net
  local.unrolled.DBN.adj.matrix.list <- 
    foreach::foreach(centralNodeIdx = 1:num.nodes, .verbose = TRUE) %:%
    foreach::foreach(time.trans.idx = 1:num.time.trans, .verbose = TRUE) %:%
    when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0) %do%
    {
      centralNodeName <- rownames(mi.net.adj.matrix)[centralNodeIdx]
      
      # List names of the central node's neighbours in mi.net.adj.matrix
      nbghNames <- c()
      if (sum(mi.net.adj.matrix[, centralNodeIdx]) == 1) # Just 1 neighbour
      {
        for (nbrIdx in 1:n)
        {
          if (mi.net.adj.matrix[nbrIdx, centralNodeIdx] == 1)
          {
            nbghNames <- rownames(mi.net.adj.matrix)[nbrIdx]
            break
          }
        }
      }
      else if (sum(mi.net.adj.matrix[, centralNodeIdx]) > 1) # Multiple neighbours
      {
        nbghNames <- rownames(mi.net.adj.matrix[which(mi.net.adj.matrix[, centralNodeIdx] == 1),])
      }
      
      local.net.node.names <- c(centralNodeName, nbghNames)
      
      # Select subdataset corr. to the nodes in 'local.net.node.names' and time points 'c(time.trans.idx, (time.trans.idx + 1))'  
      local.input.data.3D <- input.data.discr.3D[c(time.trans.idx, (time.trans.idx + 1)), local.net.node.names, ]
      
      #---------------------------------
      # Begin: Local Unrolled DBN struct learning
      #---------------------------------
      
      ## Begin: Generate 2D 'local.DBN.input.data' from 'local.input.data.3D'
      # Say, the local nodes are {v1, v2, v3} and the time points are {t1, t2}.
      # Then 'local.DBN.input.data.var.names' contains {v1_t1, v2_t1, v3_t1, v1_t2, v2_t2, v3_t2}.
      # In that case, 'local.DBN.input.data' will contain columns corr. to elements in 'local.DBN.input.data.var.names' and
      # rows corr. to different samples.
      
      ## Begin: Generate first row of 'local.DBN.input.data'
      sampleIdx <- 1
      local.DBN.input.data <- matrix(local.input.data.3D[1, , sampleIdx], nrow = 1, 
                                     ncol = dim(local.input.data.3D)[2])
      data.to.combine <- matrix(local.input.data.3D[2, , sampleIdx], nrow = 1, 
                                ncol = dim(local.input.data.3D)[2])
      local.DBN.input.data <- cbind(local.DBN.input.data, data.to.combine)
      rm(data.to.combine)
      
      #     for (time.pt in 2:dim(local.input.data.3D)[1])
      #     {
      #       data.to.combine <- matrix(local.input.data.3D[time.pt, , sampleIdx], nrow = 1, 
      #                                 ncol = dim(local.input.data.3D)[2])
      #       local.DBN.input.data <- cbind(local.DBN.input.data, data.to.combine)
      #     }
      
      local.DBN.input.data.var.names <- c()
      for (time.pt in time.trans.idx:(time.trans.idx + 1)) {
        for (var.name in dimnames(local.input.data.3D)[2]) {
          local.DBN.input.data.var.names <- c(local.DBN.input.data.var.names,
                                              paste(var.name, as.character(time.pt), sep = "_t"))
        }
      }
      colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
      ## End: Generate first row of 'local.DBN.input.data'
      
      ## Begin: Generate rest of the rows of 'local.DBN.input.data'
      if (dim(local.input.data.3D)[3] > 1) # If there are multiple samples per time pt
      {
        for (sampleIdx in 2:dim(local.input.data.3D)[3])
        {
          sample.to.combine <- matrix(local.input.data.3D[1, , sampleIdx], nrow = 1, 
                                      ncol = dim(local.input.data.3D)[2])
          data.to.combine <- matrix(local.input.data.3D[2, , sampleIdx], nrow = 1, 
                                    ncol = dim(local.input.data.3D)[2])
          sample.to.combine <- cbind(sample.to.combine, data.to.combine)
          # rm(data.to.combine)
          
          # for (time.pt in 2:dim(local.input.data.3D)[1])
          # {
          #   data.to.combine <- matrix(local.input.data.3D[time.pt, , sampleIdx], nrow = 1, 
          #                             ncol = dim(local.input.data.3D)[2])
          #   sample.to.combine <- cbind(sample.to.combine, data.to.combine)
          # }
          
          # local.DBN.input.data.var.names <- c()
          # for (time.pt in 1:dim(local.input.data.3D)[1]) {
          #   for (var.name in dimnames(local.input.data.3D)[2]) {
          #     local.DBN.input.data.var.names <- c(local.DBN.input.data.var.names,
          #                                         paste(var.name, as.character(time.pt), sep = "_t"))
          #   }
          # }
          # colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
          
          local.DBN.input.data <- rbind(local.DBN.input.data, sample.to.combine)
        }
      }
      
      rownames(local.DBN.input.data) <- as.vector(unlist(dimnames(local.input.data.3D)[3]))
      ## End: Generate rest of the rows of 'local.DBN.input.data'
      ## End: Generate 2D 'local.DBN.input.data' from 'local.input.data.3D'
      
      # local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
      #                                         discreteness = rep(TRUE, ncol(local.DBN.input.data)),                                            
      #                                         variables = colnames(local.DBN.input.data),
      #                                         node.sizes = rep(num.discr.levels, ncol(local.DBN.input.data)),
      #                                         starts.from = 0,
      #                                         num.time.steps = nrow(local.input.data.3D))
      
      
      # Begin: Uncomment this section after testing
      local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
                                                            discreteness = rep(TRUE, ncol(local.DBN.input.data)),
                                                            variables = colnames(local.DBN.input.data),
                                                            node.sizes = rep(num.discr.levels, ncol(local.DBN.input.data)),
                                                            num.time.steps = nrow(local.input.data.3D))
      
      # algo = "mmhc", scoring.func = "BDeu", layering = c()
      # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
      #                                                        num.time.steps =
      #                                                          bnstruct::num.time.steps(local.DBN.input.data.BNDataset))
      
      
      # algo = "sm", scoring.func = "BDeu", layering = c()
      # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
      #                                                      algo = 'sm',
      #                                                      num.time.steps =
      #                                                        bnstruct::num.time.steps(local.DBN.input.data.BNDataset))
      
      # algo = "sm", scoring.func = "BDeu", with layering
      # There are two time pts. represented by c(time.trans.idx, time.trans.idx + 1).
      # The nodes belonging to these time pts. are labeled with layer idx 1 and 2, resp.
      # A node with layer idx j can have parents from layer idx i such that i =< j.
      # layers <- c(rep(1, dim(local.input.data.3D)[2]), rep(2, dim(local.input.data.3D)[2]))
      # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
      #                                                        algo = 'sm',
      #                                                        num.time.steps =
      #                                                          bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
      #                                                        layering = layers)
      
      # algo = "sm", scoring.func = "BIC", with layering
      # There are two time pts. represented by c(time.trans.idx, time.trans.idx + 1).
      # The nodes belonging to these time pts. are labeled with layer idx 1 and 2, resp.
      # A node with layer idx j can have parents from layer idx i such that i =< j.
      layers <- c(rep(1, dim(local.input.data.3D)[2]), rep(2, dim(local.input.data.3D)[2]))
      local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
                                                             algo = 'sm',
                                                             scoring.func = 'BIC',
                                                             num.time.steps =
                                                               bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
                                                             layering = layers)
      
      # algo = "mmhc", scoring.func = "BDeu", with layering, initial.network = 'random.chain', seed = 12345
      # There are two time pts. represented by c(time.trans.idx, time.trans.idx + 1).
      # The nodes belonging to these time pts. are labeled with layer idx 1 and 2, resp.
      # A node with layer idx j can have parents from layer idx i such that i =< j.
      # layers <- c(rep(1, dim(local.input.data.3D)[2]), rep(2, dim(local.input.data.3D)[2]))
      # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
      #                                                        initial.network = 'random.chain',
      #                                                        seed = 12345,
      #                                                        algo = 'mmhc',
      #                                                        scoring.func = 'BDeu',
      #                                                        num.time.steps =
      #                                                          bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
      #                                                        layering = layers)
      
      # algo = "mmhc", scoring.func = "BDeu", with layering, initial.network = 'random.chain', seed = 12345
      # There are two time pts. represented by c(time.trans.idx, time.trans.idx + 1).
      # The nodes belonging to these time pts. are labeled with layer idx 1 and 2, resp.
      # A node with layer idx j can have parents from layer idx i such that i =< j.
      # layers <- c(rep(1, dim(local.input.data.3D)[2]), rep(2, dim(local.input.data.3D)[2]))
      # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
      #                                                        num.time.steps =
      #                                                          bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
      #                                                        layering = layers)
      
      
      # # The following four lines must be executed at the same time
      # save.plot.to.filename = paste(paste('LocalUnrolledDbn', centralNodeName, sep = '_'), '.jpg', sep = '')
      # jpeg(file = paste('LocalUnrolledDbn_', centralNodeName, '.jpg', sep = ''))
      # plot(local.unrolled.DBN)
      # dev.off()
      
      # Extracting the adjacency matrix of the local DBN
      local.unrolled.DBN.adj.matrix <- bnstruct::dag(local.unrolled.DBN)
      local.unrolled.DBN.adj.matrix <- matrix(local.unrolled.DBN.adj.matrix,
                                              nrow = length(local.unrolled.DBN@variables),
                                              ncol = length(local.unrolled.DBN@variables),
                                              dimnames = c(list(local.unrolled.DBN@variables), list(local.unrolled.DBN@variables)))
      # End: Uncomment this section after testing
      
      # # Begin: This section is for testing    
      # local.unrolled.DBN.adj.matrix <- matrix(0, 
      #                                         nrow = ncol(local.DBN.input.data),
      #                                         ncol = ncol(local.DBN.input.data),
      #                                         dimnames = list(colnames(local.DBN.input.data), colnames(local.DBN.input.data))) 
      # # End: This section is for testing
      
      # 'fixed = FALSE' represents that the given pattern is a regular expression.
      # local.unrolled.DBN.central.node.indices <- grep(paste('^', centralNodeName, sep = ''),
      #                                                 colnames(local.unrolled.DBN.adj.matrix),
      #                                                 fixed = FALSE)
      # Assuming there are > 1 time points. Otherwise, 'local.unrolled.DBN.adj.submatrix' would become a vector.
      # local.unrolled.DBN.adj.submatrix <- local.unrolled.DBN.adj.matrix[, local.unrolled.DBN.central.node.indices]
      
      local.unrolled.DBN.src.node.names <- local.unrolled.DBN@variables[1:dim(local.input.data.3D)[2]]
      local.unrolled.DBN.tgt.node.name <- paste(centralNodeName, as.character((time.trans.idx + 1)), sep = "_t")
      
      
      # Assuming there are > 1 time points. Otherwise, 'local.unrolled.DBN.adj.submatrix' would become a vector.
      local.unrolled.DBN.adj.submatrix <- matrix(
        local.unrolled.DBN.adj.matrix[local.unrolled.DBN.src.node.names, local.unrolled.DBN.tgt.node.name],
        nrow = length(local.unrolled.DBN.src.node.names),
        ncol = 1,
        dimnames = c(list(local.unrolled.DBN.src.node.names), list(local.unrolled.DBN.tgt.node.name)))
      
      # print(local.unrolled.DBN.adj.submatrix)
      #---------------------------------
      # End: Local Unrolled DBN struct learning
      #---------------------------------
      
      # Return value for each 'foreach' iteration 
      local.unrolled.DBN.adj.submatrix
    }
  
  # Begin: Unrolled DBN struct learning
  
  # Begin: Initialize 'unrolled.DBN.adj.matrix' as a zero matrix
  unrolled.DBN.adj.matrix.node.names <- c()
  for (time.pt in 1:num.timepts) {
    for (node.name in dimnames(input.data.discr.3D)[2]) {
      unrolled.DBN.adj.matrix.node.names <- c(unrolled.DBN.adj.matrix.node.names,
                                              paste(node.name, as.character(time.pt), sep = "_t"))
    }
  }
  
  unrolled.DBN.adj.matrix <- matrix(0, nrow = length(unrolled.DBN.adj.matrix.node.names), 
                                    ncol = length(unrolled.DBN.adj.matrix.node.names),
                                    dimnames = c(list(unrolled.DBN.adj.matrix.node.names),
                                                 list(unrolled.DBN.adj.matrix.node.names)))
  
  # End: Initialize 'unrolled.DBN.adj.matrix' as a zero matrix
  
  # 'local.unrolled.DBN.adj.matrix.list' is a list of lists of matrices.
  # The outer list contains 'local.unrolled.DBN.adj.matrix.list' inner lists.
  # Each inner list contains 'num.time.trans' matrices.
  # length(local.unrolled.DBN.adj.matrix.list) = num.nodes
  for (outer.list.idx in 1:length(local.unrolled.DBN.adj.matrix.list))
  {
    # length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]) = num.time.trans
    for (inner.list.idx in 1:length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]))
    {
      # 'submatrix.to.combine'
      submatrix.to.combine <- local.unrolled.DBN.adj.matrix.list[[outer.list.idx]][[inner.list.idx]]
      unrolled.DBN.adj.matrix[rownames(submatrix.to.combine), colnames(submatrix.to.combine)] <- submatrix.to.combine 
    }
  }
  
  # End: Unrolled DBN struct learning
  
  # Return the final directed net of the current iteration
  return (unrolled.DBN.adj.matrix)
}