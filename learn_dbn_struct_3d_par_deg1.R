## Goal: Infer Dynamic Bayesian Network (DBN)

#---------------------------------
# Begin: Step 2: Decomposing the network
#---------------------------------
## Goal: Unrolled DBN structure learning with all possible Markov Orders.
## Candidate parents: The target node itself and its CLR net neighbours at any previous and current time pt.
# Input parameters:
## input.data.3D: Dimensions {1 = time points, 2 = variables, 3 = samples under the same time point}. 
## mi.net.adj.matrix: Adjacency matrix of the mutual information network. Rownames and colnames should be node names.
## num.discr.levels: If input data is discretized, then number of discrete levels for each variable.
## Else if input data is continuous, then number of levels in which data needs to be discretized
## for performing the DBN structure learning.
## Does not use layering.
## 
learnDbnStruct3dParDeg1 <- function(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels, num.nodes, num.timepts)
{
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
  
  # load('dream3.yeast1.size10.trajectory.3D.Rdata')
  # input.data.3D <- dream3.yeast1.size10.trajectory.3D.data  
  # load('mi.net.adj.matrix.Rdata')

  # In 'di.net.adj.matrix', rows are src nodes and cols are tgt nodes
  di.net.adj.matrix <- mi.net.adj.matrix
  di.net.adj.matrix[1:nrow(di.net.adj.matrix), 1:ncol(di.net.adj.matrix)] <- 0
  
  no_cores <- num.nodes
  cl <- parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(cl)
  
  print('going inside foreach')
  
  # For each central node, learn local DBN struct in parallel.
  # Package names, specified by '.packages' must be copied to each worker core.
  # Current environmental variables that are referenced inside the foreach loop are copied
  # to each worker automatically.
  # local.unrolled.DBN.adj.matrix.list <- foreach::foreach(centralNodeIdx = 1:num.nodes, .packages = 'bnstruct')  %dopar% 
  # local.unrolled.DBN.adj.matrix.list <- as.list(1:10)
  
  # local.unrolled.DBN.adj.matrix.list <- foreach::foreach(centralNodeIdx = 1:num.nodes, .packages = 'bnstruct') %:% when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0)  %dopar%
  # {
  #   print('now inside foreach 1')
  #   
  #   centralNodeName <- rownames(mi.net.adj.matrix)[centralNodeIdx]
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
  #       local.DBN.input.data <- rbind(local.DBN.input.data, sample.to.combine)
  #     }
  #   }
  #   
  #   rownames(local.DBN.input.data) <- as.vector(unlist(dimnames(local.input.data.3D)[3]))
  #   
  #   local.unrolled.DBN.adj.matrix <- matrix(0, nrow = ncol(local.DBN.input.data), ncol = ncol(local.DBN.input.data),
  #                                           dimnames = list(colnames(local.DBN.input.data), colnames(local.DBN.input.data))) 
  #   
  #   # 'fixed = FALSE' represents that the given pattern is a regular expression.
  #   local.unrolled.DBN.central.node.indices <- grep(paste('^', centralNodeName, sep = ''),
  #                                                   colnames(local.unrolled.DBN.adj.matrix),
  #                                                   fixed = FALSE)
  #   
  #   # Assuming there are > 1 time points. Otherwise, 'local.unrolled.DBN.adj.submatrix' would become a vector.
  #   local.unrolled.DBN.adj.submatrix <- local.unrolled.DBN.adj.matrix[, local.unrolled.DBN.central.node.indices]
  #   
  #   #---------------------------------
  #   # End: Local Unrolled DBN struct learning
  #   #---------------------------------
  #   
  #   local.unrolled.DBN.adj.submatrix
  # }
  # 
  # print(local.unrolled.DBN.adj.matrix.list)
  
  local.unrolled.DBN.adj.matrix.list <- foreach::foreach(centralNodeIdx = 1:num.nodes, .packages = 'bnstruct') %:% when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0)  %dopar%
  {
    
    print('now inside foreach')
    
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
    
    
    # Begin: Uncomment this section after testing
    local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
                                                          discreteness = rep(TRUE, ncol(local.DBN.input.data)),
                                                          variables = colnames(local.DBN.input.data),
                                                          node.sizes = rep(num.discr.levels, ncol(local.DBN.input.data)),
                                                          num.time.steps = nrow(local.input.data.3D))

    # algo = "mmhc", scoring.func = "BDeu", no layering
    # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
    #                                                        num.time.steps =
    #                                                          bnstruct::num.time.steps(local.DBN.input.data.BNDataset))
    
    # algo = "mmhc", scoring.func = "BDeu", with layering.
    # A node with layer idx j can have parents from layer idx i such that i =< j.
    layers <- c()
    for (time.pt in 1:num.timepts)
    {
      layers <- c(layers, rep(time.pt, dim(local.input.data.3D)[2]))
    }
    local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
initial.network = 'random.chain', seed = 12345,                                                           
algo = 'mmhc',
                                                           scoring.func = 'BDeu',
                                                           num.time.steps =
                                                             bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
                                                           layering = layers)

    
    # algo = "sm", scoring.func = "BDeu", no layering
    # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
    #                                                        algo = 'sm',
    #                                                        num.time.steps =
    #                                                          bnstruct::num.time.steps(local.DBN.input.data.BNDataset))
    
    
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
    local.unrolled.DBN.central.node.indices <- grep(paste('^', centralNodeName, sep = ''),
                                          colnames(local.unrolled.DBN.adj.matrix),
                                          fixed = FALSE)
    
    # Assuming there are > 1 time points. Otherwise, 'local.unrolled.DBN.adj.submatrix' would become a vector.
    local.unrolled.DBN.adj.submatrix <- local.unrolled.DBN.adj.matrix[, local.unrolled.DBN.central.node.indices]
    
    #---------------------------------
    # End: Local Unrolled DBN struct learning
    #---------------------------------
    
    # Return value for each 'foreach' iteration 
    local.unrolled.DBN.adj.submatrix
  }
  
  print('foreach ended')
  
  # Shut down the cluster
  parallel::stopCluster(cl)
  
  print('after stopCluster')
  
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
  
  for (list.idx in 1:length(local.unrolled.DBN.adj.matrix.list))
  {
    # 'submatrix.to.combine'
    submatrix.to.combine <- local.unrolled.DBN.adj.matrix.list[[list.idx]]
    
	unrolled.DBN.adj.matrix[rownames(submatrix.to.combine), colnames(submatrix.to.combine)] <- submatrix.to.combine   
  }
  
  # End: Unrolled DBN struct learning
  
  # Return the final directed net of the current iteration
  return (unrolled.DBN.adj.matrix)
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
learnDbnStructLayer3dParDeg1 <- function(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels, num.nodes, num.timepts)
{
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
  
  num.time.trans <- (num.timepts - 1)
  
  no_cores <- num.nodes
  cl <- parallel::makeCluster(no_cores, outfile = paste(getwd(), 'asset/outfile.txt', sep = '/' ))
  doParallel::registerDoParallel(cl)
  
  # '.verbose = TRUE' is used for debugging
  # 'when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0' means when central node does not have any neighbour in the mutual info net
  local.unrolled.DBN.adj.matrix.list <- 
    foreach::foreach(centralNodeIdx = 1:num.nodes, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:% 
    # when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0) %:%
    foreach::foreach(time.trans.idx = 1:num.time.trans, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:%
    when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0) %do% # %dopar%
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
    # layers <- c(rep(1, dim(local.input.data.3D)[2]), rep(2, dim(local.input.data.3D)[2]))
    # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
    #                                                        algo = 'sm',
    #                                                        scoring.func = 'BIC',
    #                                                        num.time.steps =
    #                                                          bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
    #                                                        layering = layers)
    
    # algo = "mmhc", scoring.func = "BDeu", with layering, initial.network = 'random.chain', seed = 12345
    # There are two time pts. represented by c(time.trans.idx, time.trans.idx + 1).
    # The nodes belonging to these time pts. are labeled with layer idx 1 and 2, resp.
    # A node with layer idx j can have parents from layer idx i such that i =< j.
    layers <- c(rep(1, dim(local.input.data.3D)[2]), rep(2, dim(local.input.data.3D)[2]))
    local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
                                                           initial.network = 'random.chain',
                                                           seed = 12345,
                                                           algo = 'mmhc',
                                                           scoring.func = 'BDeu',
                                                           num.time.steps =
                                                             bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
                                                           layering = layers)
    
    # algo = "mmhc", scoring.func = "BDeu", with layering
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
  
  # Shut down the cluster
  parallel::stopCluster(cl)
  
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

###############################################################################################################################
## Goal: Unrolled DBN structure learning with Markov Order 0 and 1.
## Candidate parents: The target node itself and its CLR net neighbours at immediately previous and current time pt.
# Input parameters:
## input.data.3D: Dimensions {1 = time points, 2 = variables, 3 = samples under the same time point}. 
## mi.net.adj.matrix: Adjacency matrix of the mutual information network. Rownames and colnames should be node names.
## num.discr.levels: If input data is discretized, then number of discrete levels for each variable.
## Else if input data is continuous, then number of levels in which data needs to be discretized
## for performing the DBN structure learning.
## Uses layering.
## 
learnDbnStructLayer3dParDeg1 <- function(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels, num.nodes, num.timepts)
{
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
  
  num.time.trans <- (num.timepts - 1)
  
  no_cores <- num.nodes
  cl <- parallel::makeCluster(no_cores, outfile = paste(getwd(), 'asset/outfile.txt', sep = '/' ))
  doParallel::registerDoParallel(cl)
  
  # '.verbose = TRUE' is used for debugging
  # 'when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0' means when central node does not have any neighbour in the mutual info net
  local.unrolled.DBN.adj.matrix.list <- 
    foreach::foreach(centralNodeIdx = 1:num.nodes, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:% 
    # when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0) %:%
    foreach::foreach(time.trans.idx = 1:num.time.trans, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:%
    when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0) %do% # %dopar%
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
      # layers <- c(rep(1, dim(local.input.data.3D)[2]), rep(2, dim(local.input.data.3D)[2]))
      # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
      #                                                        algo = 'sm',
      #                                                        scoring.func = 'BIC',
      #                                                        num.time.steps =
      #                                                          bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
      #                                                        layering = layers)
      
      # algo = "mmhc", scoring.func = "BDeu", with layering, initial.network = 'random.chain', seed = 12345
      # There are two time pts. represented by c(time.trans.idx, time.trans.idx + 1).
      # The nodes belonging to these time pts. are labeled with layer idx 1 and 2, resp.
      # A node with layer idx j can have parents from layer idx i such that i =< j.
      layers <- c(rep(1, dim(local.input.data.3D)[2]), rep(2, dim(local.input.data.3D)[2]))
      local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
                                                             initial.network = 'random.chain',
                                                             seed = 12345,
                                                             algo = 'mmhc',
                                                             scoring.func = 'BDeu',
                                                             num.time.steps =
                                                               bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
                                                             layering = layers)
      
      # algo = "mmhc", scoring.func = "BDeu", with layering
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
  
  # Shut down the cluster
  parallel::stopCluster(cl)
  
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

###############################################################################################################################
## Goal: Unrolled DBN structure learning with Markov Order 1.
## Candidate parents: The target node itself and its CLR net neighbours at immediately previous time pt.
# Input parameters:
## input.data.3D: Dimensions {1 = time points, 2 = variables, 3 = samples under the same time point}. 
## mi.net.adj.matrix: Adjacency matrix of the mutual information network. Rownames and colnames should be node names.
## num.discr.levels: If input data is discretized, then number of discrete levels for each variable.
## Else if input data is continuous, then number of levels in which data needs to be discretized
## for performing the DBN structure learning.
## Uses layering.
## 
learnDbnStructMo1Layer3dParDeg1 <- function(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels, num.nodes, num.timepts, max.fanin)
{
  #---------------------------------
  # Begin: Loading the Packages
  #---------------------------------
  library(bnstruct)
  # library(ggm)
  library(foreach)
  library(doParallel)
  #---------------------------------
  # End: Loading the Packages
  #---------------------------------
  
  num.time.trans <- (num.timepts - 1)
  
  ## Start and register a parallel backend for parallel computing
  ## 10 cores to be used for grni server 
  # no_cores <- min(10, num.nodes, (parallel::detectCores() - 1))
  # cl <- parallel::makeCluster(no_cores, outfile = paste(getwd(), 'asset/outfile.txt', sep = '/' ))
  # doParallel::registerDoParallel(cl)
  
  # '.verbose = TRUE' is used for debugging
  # 'when(sum(mi.net.adj.matrix[, centralNodeIdx]) > 0' means when central node has at least one neighbour in the mutual info net
  # Use %do% amd %dopar% for serial and parallel computing, resp.
  local.unrolled.DBN.adj.matrix.list <- 
    foreach::foreach(centralNodeIdx = 1:num.nodes, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:% 
    # when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0) %:%
    foreach::foreach(time.trans.idx = 1:num.time.trans, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:%
    when(sum(mi.net.adj.matrix[, centralNodeIdx]) > 0) %do%
    {
      centralNodeName <- rownames(mi.net.adj.matrix)[centralNodeIdx]
      
      # List names of the central node's neighbours in mi.net.adj.matrix
      nbghNames <- c()
      if (sum(mi.net.adj.matrix[, centralNodeIdx]) == 1) # Just 1 neighbour
      {
        for (nbrIdx in 1:nrow(mi.net.adj.matrix))
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
      
      #---------------------------------
      # Begin: Local Unrolled DBN struct learning
      #---------------------------------
      
      ## Begin: Generate 2D 'local.DBN.input.data' from 'local.input.data.3D'
      # Say, the central node is v1 and the local nodes are {v1, v2, v3}, and the time points are {t1, t2}.
      # Then 'local.DBN.input.data.var.names' contains {v1_t1, v2_t1, v3_t1, v1_t2}.
      # In that case, 'local.DBN.input.data' will contain columns corr. to elements in 'local.DBN.input.data.var.names' and
      # rows corr. to different samples.
            
      num.samples <- dim(input.data.discr.3D)[3]
      local.DBN.input.data <- matrix(0, nrow = num.samples, ncol = (length(local.net.node.names) + 1))
      local.unrolled.DBN.src.node.names <- c()
      for (var.name in local.net.node.names) {
        local.unrolled.DBN.src.node.names <- c(local.unrolled.DBN.src.node.names,
                                            paste(var.name, as.character(time.trans.idx), sep = "_t"))
      }
      local.unrolled.DBN.tgt.node.name <- paste(centralNodeName, as.character(time.trans.idx + 1), sep = "_t")
      local.DBN.input.data.var.names <- c(local.unrolled.DBN.src.node.names, local.unrolled.DBN.tgt.node.name)
      colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
      for (node.name in local.net.node.names)
      {
        col.name <- paste(node.name, as.character(time.trans.idx), sep = "_t")
        local.DBN.input.data[, col.name] <- input.data.discr.3D[time.trans.idx, node.name, ]
      }
      col.name <- paste(centralNodeName, as.character(time.trans.idx + 1), sep = "_t")
      local.DBN.input.data[, col.name] <- input.data.discr.3D[(time.trans.idx + 1), centralNodeName, ]
      rm(col.name)
      
      ## End: Generate 2D 'local.DBN.input.data' from 'local.input.data.3D'
      
      local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
                                                            discreteness = rep(TRUE, ncol(local.DBN.input.data)),
                                                            variables = colnames(local.DBN.input.data),
                                                            node.sizes = rep(num.discr.levels, ncol(local.DBN.input.data)))
      
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
      # layers <- c(rep(1, (ncol(local.DBN.input.data) - 1)), 2)
      # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
      #                                                        algo = 'sm',
      #                                                        num.time.steps =
      #                                                          bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
      #                                                        layering = layers)
      # rm(layers)
            
      ## algo = "sm", scoring.func = "BIC", with layering
      ## There are two time pts. represented by c(time.trans.idx, time.trans.idx + 1).
      ## The nodes belonging to these time pts. are labeled with layer idx 1 and 2, resp.
      ## A node with layer idx j can have parents from layer idx i such that i =< j.
      layers <- c(rep(1, (ncol(local.DBN.input.data) - 1)), 2)
      local.unrolled.DBN <-  bnstruct::learn.network(local.DBN.input.data.BNDataset,
                                                             algo = 'sm',
                                                             scoring.func = 'BIC',
                                                             layering = layers)
      rm(layers)
      
      ## algo = "sm", scoring.func = "BIC", with layering, with max.fanin
      ## There are two time pts. represented by c(time.trans.idx, time.trans.idx + 1).
      ## The nodes belonging to these time pts. are labeled with layer idx 1 and 2, resp.
      ## A node with layer idx j can have parents from layer idx i such that i =< j.
      # layers <- c(rep(1, (ncol(local.DBN.input.data) - 1)), 2)
      # local.unrolled.DBN <-  bnstruct::learn.network(local.DBN.input.data.BNDataset,
      #                                                algo = 'sm',
      #                                                scoring.func = 'BIC',
      #                                                max.fanin = max.fanin,
      #                                                layering = layers)
      # rm(layers)
      
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
      
      # algo = "mmhc", scoring.func = "BDeu", with layering
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
      
      # This for loop checks whether parents are learnt only for 'local.unrolled.DBN.tgt.node.name'. If
      # so, then nothing is printed. Otherwise, prints the column(s) corr. to the undesired tgt node(s).
      for (col.idx in 1:(ncol(local.unrolled.DBN.adj.matrix) - 1))
      {
        if (sum(local.unrolled.DBN.adj.matrix[, col.idx]) > 0)
        {
          print('Erroneous column')
          print(local.unrolled.DBN.adj.matrix[, col.idx])
        }
      }
      
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
  
  print('End of foreach loops')
  
  ## Shut down the cluster
  # parallel::stopCluster(cl)
    
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
  
  ## Debugging
  # writeLines('length(local.unrolled.DBN.adj.matrix.list) = ', length(local.unrolled.DBN.adj.matrix.list), '\n')
  # print(length(local.unrolled.DBN.adj.matrix.list))
  
  
  # 'local.unrolled.DBN.adj.matrix.list' is a list of lists of matrices.
  # The outer list contains 'local.unrolled.DBN.adj.matrix.list' inner lists.
  # Each inner list contains 'num.time.trans' matrices.
  # length(local.unrolled.DBN.adj.matrix.list) = num.nodes
  for (outer.list.idx in 1:length(local.unrolled.DBN.adj.matrix.list))
  {
    ## Debugging
    # writeLines('outer.list.idx = ', outer.list.idx, '\n')
    # print(outer.list.idx)
    # writeLines('length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]) = ', 
    #            length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]), '\n')
    # print(length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]))
    # print(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]])
    
    # if 'local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]' is not an empty list
    if (length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]) > 0)
    {
      # length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]) = num.time.trans
      for (inner.list.idx in 1:length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]))
      {
        ## Debugging
        # writeLines('inner.list.idx = ', inner.list.idx, '\n')
        # print(inner.list.idx)
        
        # 'submatrix.to.combine'
        submatrix.to.combine <- local.unrolled.DBN.adj.matrix.list[[outer.list.idx]][[inner.list.idx]]
        unrolled.DBN.adj.matrix[rownames(submatrix.to.combine), colnames(submatrix.to.combine)] <- submatrix.to.combine 
      }
    }
  }
  
  # End: Unrolled DBN struct learning
  
  # Return the final directed net of the current iteration
  return (unrolled.DBN.adj.matrix)
}

###############################################################################################################################
## Goal: Unrolled DBN structure learning with Markov Order 1.
## This function is the newer version of learnDbnStructMo1Layer3dParDeg1().
## The only difference is in the size of unrolled DBN adjacency matrix. In earlier version,
## the size is ((V \times T) \times (v \times T)) where V = number of nodes and
## T = number of time points. But the size is too large when V is very large. For
## e.g., when V = 4028, the function can not execute successfully even with 32 GB main
## memory in grni server. In this version, the size is reduced by storing the unrolled DBN
## in an adjacency list of length (T - 1). The t^{th} element in the list is
## the predicted network adjacency matrix at the t^{th} time interval. Therefore, each list element
## is a binary matrix of dimension (V \times V). Hence, the total size of the unrolled DBN
## adjacency list is ((T - 1) \times (V \times V)).
##
## Candidate parents: The target node itself and its CLR net neighbours at immediately previous time pt.
# Input parameters:
## input.data.3D: Dimensions {1 = time points, 2 = variables, 3 = samples under the same time point}. 
## mi.net.adj.matrix: Adjacency matrix of the mutual information network. Rownames and colnames should be node names.
## num.discr.levels: If input data is discretized, then number of discrete levels for each variable.
## Else if input data is continuous, then number of levels in which data needs to be discretized
## for performing the DBN structure learning.
## Uses layering.
## 
LearnDbnStructMo1Layer3dParDeg1_v2 <- function(input.data.discr.3D, 
                                               mi.net.adj.matrix, 
                                               num.discr.levels, 
                                               num.nodes, 
                                               num.timepts, 
                                               max.fanin, 
                                               node.names, 
                                               clr.algo) {
  #---------------------------------
  # Begin: Loading the Packages
  #---------------------------------
  library(bnstruct)
  # library(ggm)
  library(foreach)
  library(doParallel)
  #---------------------------------
  # End: Loading the Packages
  #---------------------------------
  
  num.time.trans <- (num.timepts - 1)
  
  ## Start and register a parallel backend for parallel computing
  ## 10 cores to be used for grni server 
  # no_cores <- min(10, num.nodes, (parallel::detectCores() - 1))
  # cl <- parallel::makeCluster(no_cores, outfile = paste(getwd(), 'asset/outfile.txt', sep = '/' ))
  # doParallel::registerDoParallel(cl)
  
  # '.verbose = TRUE' is used for debugging
  # 'when(sum(mi.net.adj.matrix[, centralNodeIdx]) > 0' means when central node has at least one neighbour in the mutual info net
  # Use %do% amd %dopar% for serial and parallel computing, resp.
  local.unrolled.DBN.adj.matrix.list <- 
    foreach::foreach(centralNodeIdx = 1:num.nodes, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:% 
    # when(sum(mi.net.adj.matrix[, centralNodeIdx]) != 0) %:%
    foreach::foreach(time.trans.idx = 1:num.time.trans, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:%
    when(sum(mi.net.adj.matrix[, centralNodeIdx]) > 0) %do%
    {
      centralNodeName <- rownames(mi.net.adj.matrix)[centralNodeIdx]
      
      # List names of the central node's neighbours in mi.net.adj.matrix
      nbghNames <- c()
      if (sum(mi.net.adj.matrix[, centralNodeIdx]) == 1) # Just 1 neighbour
      {
        for (nbrIdx in 1:nrow(mi.net.adj.matrix))
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
      
      local.net.node.names <- c()
      if (clr.algo == 'CLR2.1') {
        local.net.node.names <- nbghNames
      } else if (clr.algo %in% c('CLR', 'CLR1.2', 'CLR2', 'CLR10')) {
        local.net.node.names <- c(centralNodeName, nbghNames)
      }
      
      #---------------------------------
      # Begin: Local Unrolled DBN struct learning
      #---------------------------------
      
      ## Begin: Generate 2D 'local.DBN.input.data' from 'local.input.data.3D'
      # Say, the central node is v1 and the local nodes are {v1, v2, v3}, and the time points are {t1, t2}.
      # Then 'local.DBN.input.data.var.names' contains {v1_t1, v2_t1, v3_t1, v1_t2}.
      # In that case, 'local.DBN.input.data' will contain columns corr. to elements in 'local.DBN.input.data.var.names' and
      # rows corr. to different samples.
      
      num.samples <- dim(input.data.discr.3D)[3]
      
      local.DBN.input.data <- matrix(0, nrow = num.samples, ncol = (length(local.net.node.names) + 1))
      
      local.unrolled.DBN.src.node.names <- c()
      for (var.name in local.net.node.names) {
        local.unrolled.DBN.src.node.names <- c(local.unrolled.DBN.src.node.names,
                                               paste(var.name, as.character(time.trans.idx), sep = "_t"))
      }
      rm(var.name)
      
      local.unrolled.DBN.tgt.node.name <- paste(centralNodeName, as.character(time.trans.idx + 1), sep = "_t")
      
      local.DBN.input.data.var.names <- c(local.unrolled.DBN.src.node.names, local.unrolled.DBN.tgt.node.name)
      
      colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
      
      for (node.name in local.net.node.names) {
        col.name <- paste(node.name, as.character(time.trans.idx), sep = "_t")
        local.DBN.input.data[, col.name] <- input.data.discr.3D[time.trans.idx, node.name, ]
      }
      rm(node.name)
      
      col.name <- paste(centralNodeName, as.character(time.trans.idx + 1), sep = "_t")
      local.DBN.input.data[, col.name] <- input.data.discr.3D[(time.trans.idx + 1), centralNodeName, ]
      rm(col.name)
      
      ## End: Generate 2D 'local.DBN.input.data' from 'local.input.data.3D'
      # save(local.DBN.input.data, file = 'local.DBN.input.data.RData')
      local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
                                                            discreteness = rep(TRUE, ncol(local.DBN.input.data)),
                                                            variables = colnames(local.DBN.input.data),
                                                            node.sizes = rep(num.discr.levels, 
                                                                             ncol(local.DBN.input.data)))
      
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
      # layers <- c(rep(1, (ncol(local.DBN.input.data) - 1)), 2)
      # local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset,
      #                                                        algo = 'sm',
      #                                                        num.time.steps =
      #                                                          bnstruct::num.time.steps(local.DBN.input.data.BNDataset),
      #                                                        layering = layers)
      # rm(layers)
      
      ## algo = "sm", scoring.func = "BIC", with layering
      ## There are two time pts. represented by c(time.trans.idx, time.trans.idx + 1).
      ## The nodes belonging to these time pts. are labeled with layer idx 1 and 2, resp.
      ## A node with layer idx j can have parents from layer idx i such that i =< j.
      layers <- c(rep(1, (ncol(local.DBN.input.data) - 1)), 2)
      local.unrolled.DBN <-  bnstruct::learn.network(local.DBN.input.data.BNDataset,
                                                     algo = 'sm',
                                                     scoring.func = 'BIC',
                                                     layering = layers)
      rm(layers)
      
      ## algo = "sm", scoring.func = "BIC", with layering, with max.fanin
      ## There are two time pts. represented by c(time.trans.idx, time.trans.idx + 1).
      ## The nodes belonging to these time pts. are labeled with layer idx 1 and 2, resp.
      ## A node with layer idx j can have parents from layer idx i such that i =< j.
      # layers <- c(rep(1, (ncol(local.DBN.input.data) - 1)), 2)
      # local.unrolled.DBN <-  bnstruct::learn.network(local.DBN.input.data.BNDataset,
      #                                                algo = 'sm',
      #                                                scoring.func = 'BIC',
      #                                                max.fanin = max.fanin,
      #                                                layering = layers)
      # rm(layers)
      
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
      
      # algo = "mmhc", scoring.func = "BDeu", with layering
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
      
      # This for loop checks whether parents are learnt only for 'local.unrolled.DBN.tgt.node.name'. If
      # so, then nothing is printed. Otherwise, prints the column(s) corr. to the undesired tgt node(s).
      for (col.idx in 1:(ncol(local.unrolled.DBN.adj.matrix) - 1))
      {
        if (sum(local.unrolled.DBN.adj.matrix[, col.idx]) > 0)
        {
          print('Erroneous column')
          print(local.unrolled.DBN.adj.matrix[, col.idx])
        }
      }
      
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
      
      ## Assuming there are > 1 time points. Otherwise, 'local.unrolled.DBN.adj.submatrix' would become a vector.
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
  print('End of foreach loops')
  
  # save(local.unrolled.DBN.adj.matrix.list, file = paste(getwd(), 'asset/local.unrolled.DBN.adj.matrix.list.RData', sep = '/'))
  
  ## Shut down the cluster
  # parallel::stopCluster(cl)
  
  # Begin: Unrolled DBN struct learning
  
  # Initialize the unrolled DBN adjacency list.
  ## It is a list of length (T - 1) where T = number of time points. The t^{th} element in the list is
  ## the predicted network adjacency matrix at the t^{th} time interval. Therefore, each list element
  ## is a binary matrix of dimension (V \times V) where V = number of nodes. Hence, the total size 
  ## of the unrolled DBN adjacency list is ((T - 1) \times (V \times V)).
  ## Each adjacency matrix is initialized with a zero matrix of dimension (V \times V). Its rows
  ## corr. to soruce nodes and the columns corr. to target nodes.
  unrolled.DBN.adj.matrix.list <- list()
  for (list.idx in 1:num.time.trans)
  {
    unrolled.DBN.adj.matrix.list[[list.idx]] <- matrix(0, nrow = num.nodes,
                                                       ncol = num.nodes,
                                                       dimnames = 
                                                         c(list(node.names), 
                                                           list(node.names)))
  }
  rm(list.idx)
  
  ## 'local.unrolled.DBN.adj.matrix.list' is a list of lists of matrices.
  ## The outer list contains length(local.unrolled.DBN.adj.matrix.list) number of inner lists.
  ## The length is \ge 1 and \le num.nodes. It is < num.nodes when there exists some node
  ## without any neighbour in the mutual info net. 
  ##
  ## Each inner list contains 'num.time.trans' matrices.
  for (outer.list.idx in 1:length(local.unrolled.DBN.adj.matrix.list))
  {
    ## Debugging
    # writeLines('outer.list.idx = ', outer.list.idx, '\n')
    # print(outer.list.idx)
    # writeLines('length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]) = ', 
    #            length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]), '\n')
    # print(length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]))
    # print(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]])
    
    ## if 'local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]' is not an empty list
    if (length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]) > 0)
    {
      ## length(local.unrolled.DBN.adj.matrix.list[[outer.list.idx]]) = num.time.trans
      ## for any value of outer.list.idx
      for (inner.list.idx in 1:num.time.trans)
      {
        ## Debugging
        # writeLines('inner.list.idx = ', inner.list.idx, '\n')
        # print(inner.list.idx)
        
        ## 'submatrix.to.combine'
        submatrix.to.combine <- local.unrolled.DBN.adj.matrix.list[[outer.list.idx]][[inner.list.idx]]
        
        ## Begin: remove the timestamps from the row and col names of submatrix.to.combine.
        ## E.g., a row or col name 'G2_t34' will be converted to just 'G2'.
        submatrix.to.combine.rownames <- c()
        for (row.idx in 1:nrow(submatrix.to.combine))
        {
          old.rowname <- rownames(submatrix.to.combine)[row.idx]
          
          ## Substitute '_t[0-9]+' pattern with empty string in old.rowname
          new.rowname <- sub('_t[0-9]+', '', old.rowname)
          
          submatrix.to.combine.rownames <- c(submatrix.to.combine.rownames, new.rowname)
        }
        rm(row.idx)
        
        rownames(submatrix.to.combine) <- submatrix.to.combine.rownames
        rm(submatrix.to.combine.rownames)
        
        submatrix.to.combine.colnames <- c()
        for (col.idx in 1:ncol(submatrix.to.combine))
        {
          old.colname <- colnames(submatrix.to.combine)[col.idx]
          
          ## Replace '_t[0-9]+' pattern with empty string in old.rowname
          new.colname <- sub('_t[0-9]+', '', old.colname)
          
          submatrix.to.combine.colnames <- c(submatrix.to.combine.colnames, new.colname)
        }
        rm(col.idx)
        
        colnames(submatrix.to.combine) <- submatrix.to.combine.colnames
        rm(submatrix.to.combine.colnames)
        ## End: remove the timestamps from the row and col names of submatrix.to.combine.
        
        
        # print('submatrix.to.combine')
        # print(submatrix.to.combine)
        # 
        # print('inner.list.idx')
        # print(inner.list.idx)
        # print('unrolled.DBN.adj.matrix.list[[inner.list.idx]]')
        # print(unrolled.DBN.adj.matrix.list[[inner.list.idx]])
        
        unrolled.DBN.adj.matrix.list[[inner.list.idx]][rownames(submatrix.to.combine), 
                                                       colnames(submatrix.to.combine)] <- submatrix.to.combine
        
      }
    }
  }
  rm(outer.list.idx)
  
  # End: Unrolled DBN struct learning
  
  return (unrolled.DBN.adj.matrix.list)
}

############################################################################################

############################################################################################
## Goal: Learn DBN structure of Markov order 1 where 
## candidate parents are selected using the CLR3 algo.
## It is a serial algorithmic implementation.
## 
LearnDbnStructMo1Clr3Ser <- function(input.data.discr.3D, mi.net.adj.matrix.list.filename, 
                                     num.discr.levels, num.nodes, num.timepts, max.fanin, 
                                     node.names, unrolled.DBN.adj.matrix.list)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  library(bnstruct)
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Load the adjacency matrices of time-interval-specific 
    ## CLR nets from the given file
    mi.net.adj.matrix.list <- NULL
    load(mi.net.adj.matrix.list.filename)
    
    ## Extract the CLR net specific to the current time interval
    mi.net.adj.matrix <- mi.net.adj.matrix.list[[time.pt.idx]]
    rm(mi.net.adj.matrix.list.filename)
    
    for (tgt.node.idx in 1:num.nodes) {
      
      ## If the current target node has no candidate parent,
      ## then skip to the next target node
      if (sum(mi.net.adj.matrix[, tgt.node.idx]) == 0) {
        next
      }
      
      ## E.g., 'v1_t2'.
      ## Do not use rownames.
      tgt.node.name <- colnames(mi.net.adj.matrix)[tgt.node.idx]
      
      ## E.g., 'v1'
      tgt.node.base.name <- node.names[tgt.node.idx]
      
      ## Not required for the CLR3 algo
      # tgt.node.name.as.src <- paste(tgt.node.base.name, as.character(time.pt.idx), sep = "_t")
      
      ## List names of the target node's neighbours in 'mi.net.adj.matrix'
      ngbr.names <- c()
      if (sum(mi.net.adj.matrix[, tgt.node.idx]) == 1) {
        ## Just one neighbour
        
        for (ngbr.idx in 1:nrow(mi.net.adj.matrix)) {
          if (mi.net.adj.matrix[ngbr.idx, tgt.node.idx] == 1) {
            
            ## Do not use colnames
            ngbr.names <- rownames(mi.net.adj.matrix)[ngbr.idx]
            
            break
          }
        }
        rm(ngbr.idx)
        
      } else if (sum(mi.net.adj.matrix[, tgt.node.idx]) > 1) {
        ## Multiple neighbours
        
        ## Do not use colnames
        ngbr.names <- rownames(mi.net.adj.matrix[which(mi.net.adj.matrix[, tgt.node.idx] == 1),])
      }
      
      ## Do not add 'tgt.node.name.as.src' to 'local.net.node.names'.
      ## Because, in case of CLR3 algo, if a target node (say, v7_t2)
      ## has v7_t1 as a candidate regulator, then v7-t1 already 
      ## belongs to 'ngbr.names'.
      # local.net.node.names <- c(tgt.node.name.as.src, ngbr.names)
      local.net.node.names <- ngbr.names
      
      ##------------------------------------------------------------
      ## Begin: Local Unrolled DBN struct learning
      ##------------------------------------------------------------
      
      num.samples.per.timept <- dim(input.data.discr.3D)[3]
      
      local.DBN.input.data <- matrix(0, nrow = num.samples.per.timept, 
                                     ncol = (length(local.net.node.names) + 1))
      
      
      
      local.unrolled.DBN.src.node.names <- c()
      
      ## Do not add 'tgt.node.name.as.src' to 'local.unrolled.DBN.src.node.names'.
      ## Because, in case of CLR3 algo, if a target node (say, v7_t2)
      ## has v7_t1 as a candidate regulator, then v7-t1 already 
      ## belongs to 'ngbr.names'.
      # local.unrolled.DBN.src.node.names <- c(ngbr.names, tgt.node.name.as.src)
      local.unrolled.DBN.src.node.names <- ngbr.names
      
      local.unrolled.DBN.src.node.names <- sort(local.unrolled.DBN.src.node.names)
      
      # rm(tgt.node.name.as.src)
      
      local.unrolled.DBN.tgt.node.name <- tgt.node.name
      
      local.DBN.input.data.var.names <- c(local.unrolled.DBN.src.node.names, local.unrolled.DBN.tgt.node.name)
      colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
      rm(local.DBN.input.data.var.names)
      
      for (node.name in local.net.node.names) {
        
        ## Substitute '_t[0-9]+' pattern with empty string in 'node.name'.
        ## If 'node.name' = 'v1_t1', then
        ## 'node.base.name' = 'v1'.
        node.base.name <- sub('_t[0-9]+', '', node.name)
        
        local.DBN.input.data[, node.name] <- input.data.discr.3D[time.pt.idx, node.base.name, ]
        
      }
      rm(node.name)
      
      local.DBN.input.data[, tgt.node.name] <- input.data.discr.3D[(time.pt.idx + 1), tgt.node.base.name, ]
      
      local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
                                                            discreteness = rep(TRUE, ncol(local.DBN.input.data)),
                                                            variables = colnames(local.DBN.input.data),
                                                            node.sizes = rep(num.discr.levels, 
                                                                             ncol(local.DBN.input.data)))
      
      layers <- c(rep(1, (ncol(local.DBN.input.data) - 1)), 2)
      local.unrolled.DBN <-  bnstruct::learn.network(local.DBN.input.data.BNDataset,
                                                     algo = 'sm',
                                                     scoring.func = 'BIC',
                                                     layering = layers)
      rm(layers)
      
      ## Extracting the adjacency matrix of the local DBN
      local.unrolled.DBN.adj.matrix <- bnstruct::dag(local.unrolled.DBN)
      local.unrolled.DBN.adj.matrix <- matrix(local.unrolled.DBN.adj.matrix,
                                              nrow = length(local.unrolled.DBN@variables),
                                              ncol = length(local.unrolled.DBN@variables),
                                              dimnames = c(list(local.unrolled.DBN@variables), 
                                                           list(local.unrolled.DBN@variables)))
      
      # This for loop checks whether parents are learnt only for 'local.unrolled.DBN.tgt.node.name'. If
      # so, then nothing is printed. Otherwise, prints the column(s) corr. to the undesired tgt node(s).
      for (col.idx in 1:(ncol(local.unrolled.DBN.adj.matrix) - 1)) {
        if (sum(local.unrolled.DBN.adj.matrix[, col.idx]) > 0) {
          print('Erroneous column')
          print(local.unrolled.DBN.adj.matrix[, col.idx])
        }
      }
      rm(col.idx)
      
      ## Assuming there are > 1 time points. 
      ## Otherwise, 'submatrix.to.combine' would become a vector.
      submatrix.to.combine <- matrix(
        local.unrolled.DBN.adj.matrix[local.unrolled.DBN.src.node.names, local.unrolled.DBN.tgt.node.name],
        nrow = length(local.unrolled.DBN.src.node.names),
        ncol = 1,
        dimnames = c(list(local.unrolled.DBN.src.node.names), list(local.unrolled.DBN.tgt.node.name)))
      
      ##------------------------------------------------------------
      ## End: Local Unrolled DBN struct learning
      ##------------------------------------------------------------
      
      ## Begin: remove the timestamps from the row names of 
      ## 'submatrix.to.combine'.
      ## E.g., a row name 'G2_t34' will be converted to just 'G2'.
      submatrix.to.combine.rownames <- c()
      for (row.idx in 1:nrow(submatrix.to.combine)) {
        old.rowname <- rownames(submatrix.to.combine)[row.idx]
        
        ## Substitute '_t[0-9]+' pattern with empty string in old.rowname
        new.rowname <- sub('_t[0-9]+', '', old.rowname)
        
        submatrix.to.combine.rownames <- c(submatrix.to.combine.rownames, new.rowname)
      }
      rm(row.idx)
      
      rownames(submatrix.to.combine) <- submatrix.to.combine.rownames
      rm(submatrix.to.combine.rownames)
      
      ## End: remove the timestamps from the row names of submatrix.to.combine.
      
      unrolled.DBN.adj.matrix.list[[time.pt.idx]][rownames(submatrix.to.combine), 
                                                  tgt.node.base.name] <- submatrix.to.combine
    }
    rm(tgt.node.idx)
  }
  rm(time.pt.idx)
  
  return(unrolled.DBN.adj.matrix.list)
}
############################################################################################