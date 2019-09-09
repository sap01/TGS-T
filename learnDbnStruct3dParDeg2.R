# Goal: Infer Dynamic Bayesian Network (DBN)

#---------------------------------
# Begin: Loading the Packages
#---------------------------------
library(bnstruct)
library(ggm)
library(parallel)
library(snow)
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
learnDbnStruct3dParDeg2 <- function(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels, num.nodes, num.timepts)
{
  # load('dream3.yeast1.size10.trajectory.3D.Rdata')
  # input.data.3D <- dream3.yeast1.size10.trajectory.3D.data  
  # load('mi.net.adj.matrix.Rdata')

  # In 'di.net.adj.matrix', rows are src nodes and cols are tgt nodes
  di.net.adj.matrix <- mi.net.adj.matrix
  di.net.adj.matrix[1:nrow(di.net.adj.matrix), 1:ncol(di.net.adj.matrix)] <- 0
  
  num.time.trans <- num.timepts - 1 # Num of time transitions 
  
  num.clnodes <- num.nodes # Num of cluster nodes
  num.cores.per.clnode <- num.time.trans # Num of cores per cluster node
  
  ## Start a multi-node cluster
  
  ## Start a FORK cluster 
  # Type of cluster to use: FORK or SOCK? 
  # Ref: https://www.r-bloggers.com/how-to-go-parallel-in-r-basics-tips/
  # Since 'grni' is used to test this program, which is a single Ubuntu Linux node
  # with 24 logical CPUs, 'SOCK' is only able to create (1 + 10) = 11 workers.
  # Whereas, 'FORK' is able to create 22 workers.
  # Uncomment the following line to start a FORK cluster
  # cl <- parallel::makeCluster(num.clnodes, type = "FORK")
  ## End(Start a FORK cluster)

  ## Start an MPI cluster
  cl <- snow::makeCluster(num.clnodes, type = 'MPI', outfile = '/home/p.saptarshi/R-3.3.2/projects/outfile.txt')
  ## End(Start an MPI cluster)
  
  ## End(Start a multi-node cluster)
  
  ## Fork multiple nodes per cluster node
  snow::clusterCall(cl, function() {
    library(doSNOW)
    doSNOW::registerDoSNOW(parallel::makeCluster(num.cores.per.clnode, type = "FORK"))
    NULL
  })
  registerDoSNOW(cl)
  ## End(Fork multiple nodes per cluster node)
  
  # For each central node, learn local DBN struct in parallel.
  # Package names, specified by '.packages' must be copied to each worker core.
  # Current environmental variables that are referenced inside the foreach loop are copied
  # to each worker automatically.
  local.unrolled.DBN.adj.matrix.list <- foreach::foreach(centralNodeIdx = 1:num.nodes, .combine = 'list', .packages = c('foreach', 'bnstruct'))  %:% 
    foreach::foreach(time.trans.idx = 1:num.time.trans, .combine = 'list', .packages = 'bnstruct')  %:% when (sum(mi.net.adj.matrix[, centralNodeIdx]) >= 1) %dopar% 
    {
      centralNodeName <- rownames(mi.net.adj.matrix)[centralNodeIdx]
      
      # List names of the central node's neighbours in mi.net.adj.matrix
      nbghNames <- c()
      if (sum(mi.net.adj.matrix[, centralNodeIdx]) == 1) # Just 1 neighbour
      {
        for (nbrIdx in 1:num.nodes)
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
      num.local.nodes <- length(local.net.node.names)
      
      # Input data[current and next time point, central node and its candidate parent nodes, all samples at each time point]
      local.input.data.3D <- input.data.discr.3D[time.trans.idx:(time.trans.idx + 1) , local.net.node.names, ]
      
      ## Local Unrolled DBN struct learning
      
      
      ## Convert 'local.DBN.input.data' from 3D array[time pt., variables, samples per time pt.]
      ## to a 2D matrix[samples per time pt., (time pt. * variables)]
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
      
      ## End(Convert 'local.DBN.input.data' from 3D array[time pt., variables, samples per time pt.]
      ## to a 2D matrix[samples per time pt., (time pt. * variables)])
      
      # local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
      #                                         discreteness = rep(TRUE, ncol(local.DBN.input.data)),                                            
      #                                         variables = colnames(local.DBN.input.data),
      #                                         node.sizes = rep(num.discr.levels, ncol(local.DBN.input.data)),
      #                                         starts.from = 0,
      #                                         num.time.steps = nrow(local.input.data.3D))
      
      # Only 2 time steps are considered i.e. the current time pt. and the immediate next time pt.
      local.DBN.input.data.BNDataset <- bnstruct::BNDataset(local.DBN.input.data,
                                                            discreteness = rep(TRUE, ncol(local.DBN.input.data)),                                            
                                                            variables = colnames(local.DBN.input.data),
                                                            node.sizes = rep(num.discr.levels, ncol(local.DBN.input.data)),
                                                            num.time.steps = 2)
      
      # algo = "mmhc", scoring.func = "BDeu"
      local.unrolled.DBN <-  bnstruct::learn.dynamic.network(local.DBN.input.data.BNDataset, 
                                                             num.time.steps = 
                                                               bnstruct::num.time.steps(local.DBN.input.data.BNDataset))
      
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
      
      ## End(Local Unrolled DBN struct learning)
      
      local.unrolled.DBN.adj.submatrix <- matrix(local.unrolled.DBN.adj.matrix[1:num.local.nodes, (num.local.nodes + 1)], 
                                              nrow = length(1:num.local.nodes), 
                                              ncol = 1,
                                              dimnames = c(list((local.unrolled.DBN@variables)[1:num.local.nodes]), 
                                                           list((local.unrolled.DBN@variables)[num.local.nodes + 1])))
      
      # Return value for the inner foreach loop
      local.unrolled.DBN.adj.submatrix
      
    }
  
  ## Shut down the parallel backend
  snow::stopCluster(cl)
  ## End(Shut down the parallel backend)
  
  ## Generate Unrolled DBN adjacency matrix
  
  ## Initialize 'unrolled.DBN.adj.matrix' as a zero matrix
  unrolled.DBN.adj.matrix.node.names <- c()
  for (time.pt in 1:num.timepts) {
    for (node.name in dimnames(local.input.data.3D)[2]) {
      unrolled.DBN.adj.matrix.node.names <- c(unrolled.DBN.adj.matrix.node.names,
                                          paste(node.name, as.character(time.pt), sep = "_t"))
    }
  }
  
  unrolled.DBN.adj.matrix <- matrix(0, nrow = length(unrolled.DBN.adj.matrix.node.names), 
                                    ncol = length(unrolled.DBN.adj.matrix.node.names),
                                    dimnames = c(list(unrolled.DBN.adj.matrix.node.names),
                                                 list(unrolled.DBN.adj.matrix.node.names)))
  
  ## End(Initialize 'unrolled.DBN.adj.matrix' as a zero matrix)
  
  for (centralNodeIdx in 1:num.nodes)
  {
    for (time.trans.idx in 1:num.time.trans)
    {
      # 'submatrix.to.combine' is a single column matrix
      submatrix.to.combine <- local.unrolled.DBN.adj.matrix.list[[centralNodeIdx]][[time.trans.idx]]
      
      for (row.idx in 1:nrow(submatrix.to.combine))
      {
        unrolled.DBN.adj.matrix[rownames(submatrix.to.combine)[row.idx], colnames(submatrix.to.combine)[1]] <- submatrix.to.combine[row.idx, 1]
      }
      
    }
  }
  ## End(Generate Unrolled DBN adjacency matrix)
  
  ## Return unrolled DBN adj matrix
  return(unrolled.DBN.adj.matrix)
}