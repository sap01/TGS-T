## Goal: Infer Dynamic Bayesian Network (DBN)

###############################################################################################################################
## Goal: Unrolled DBN structure learning with Markov Order 1.
## This version is the Lite version of learn_dbn_struct_3d_par_deg1.R::LearnDbnStructMo1Layer3dParDeg1_v2().
## It is a serial implemention.
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
LearnDbnStructMo1Layer3dSerDeg1_v2_Lite <- function(input.data.discr.3D, 
                                                    mi.net.adj.matrix, 
                                                    num.discr.levels, 
                                                    num.nodes, 
                                                    num.timepts, 
                                                    node.names, 
                                                    clr.algo, 
                                                    auto.reg.order, 
                                                    scoring.func) {
  #---------------------------------
  # Begin: Loading the Packages
  #---------------------------------
  library(bnlearn)
  # library(ggm)
  library(foreach)
  #---------------------------------
  # End: Loading the Packages
  #---------------------------------
  
  num.time.trans <- (num.timepts - 1)
  
  ## '.verbose = TRUE' is used for debugging
  ## 'when(sum(mi.net.adj.matrix[, tgt.node.idx]) > 0' means when central node has at least one neighbour in the mutual info net
  ## Use %do% amd %dopar% for serial and parallel computing, resp.
  local.unrolled.DBN.adj.matrix.list <- 
    foreach::foreach(tgt.node.idx = 1:num.nodes, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:% 
    # when(sum(mi.net.adj.matrix[, tgt.node.idx]) != 0) %:%
    foreach::foreach(time.trans.idx = 1:num.time.trans, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:%
    when(sum(mi.net.adj.matrix[, tgt.node.idx]) > 0) %do% {
      
      ## Already done in TGS.R
      # source(paste(init.path, 'learn_local_dbn.R', sep = '/'))
      # source(paste(init.path, 'score_bn.R', sep = '/'))
      # 
      # if (.Platform$OS.type == 'windows') {
      #   dyn.load(paste(init.path, 'src/windows/score_bn.dll', sep = '/'))
      # } else if (.Platform$OS.type == 'unix') {
      #   dyn.load(paste(init.path, 'src/unix/score_bn.so', sep = '/'))
      # }
      
      tgt.node.name <- rownames(mi.net.adj.matrix)[tgt.node.idx]
      
      # List names of the target node's neighbours in mi.net.adj.matrix
      nbr.names <- c()
      
      ## Just one neighbour
      if (sum(mi.net.adj.matrix[, tgt.node.idx]) == 1) {
        for (nbr.idx in 1:nrow(mi.net.adj.matrix)) {
          if (mi.net.adj.matrix[nbr.idx, tgt.node.idx] == 1) {
            nbr.names <- rownames(mi.net.adj.matrix)[nbr.idx]
            break
          }
        }
      } else if (sum(mi.net.adj.matrix[, tgt.node.idx]) > 1) {
        ## Multiple neighbours
        
        nbr.names <- rownames(mi.net.adj.matrix[which(mi.net.adj.matrix[, tgt.node.idx] == 1),])
      }
      
      local.net.node.names <- c()
      if (clr.algo == 'CLR2.1') {
        local.net.node.names <- nbr.names
      } else if (clr.algo %in% c('CLR', 'CLR1.2', 'CLR2', 'CLR10')) {
        local.net.node.names <- c(tgt.node.name, nbr.names)
      }
      
      #---------------------------------
      # Begin: Local Unrolled DBN struct learning
      #---------------------------------
      
      ## Begin: Generate 2D 'local.DBN.input.data' from 'local.input.data.3D'
      # Say, the target node is v1 and the local nodes are {v1, v2, v3}, 
      # and 'time.trans.idx' = 1 and 'auto.reg.order' = 1.
      # Then 'local.DBN.input.data.var.names' contains {v1_t1, v2_t1, v3_t1, v1_t2}.
      # In that case, 'local.DBN.input.data' is a matrix which
      # will contain columns corr. to elements 
      # in 'local.DBN.input.data.var.names' and
      # rows corr. to different samples.
      
      num.samples <- dim(input.data.discr.3D)[3]
      
      local.DBN.input.data <- matrix(0, nrow = num.samples, ncol = (length(local.net.node.names) + 1))
      
      local.unrolled.DBN.src.node.names <- c()
      for (var.name in local.net.node.names) {
        local.unrolled.DBN.src.node.names <- c(local.unrolled.DBN.src.node.names,
                                               paste(var.name, as.character(time.trans.idx), 
                                                     sep = "_t"))
      }
      rm(var.name)
      
      local.unrolled.DBN.tgt.node.name <- paste(tgt.node.name, 
                                                as.character(time.trans.idx + auto.reg.order), 
                                                sep = "_t")
      
      local.DBN.input.data.var.names <- c(local.unrolled.DBN.src.node.names, 
                                          local.unrolled.DBN.tgt.node.name)
      
      colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
      
      ## Fill in the data of the candidate source nodes
      ## in 'local.DBN.input.data'
      for (node.name in local.net.node.names) {
        col.name <- paste(node.name, as.character(time.trans.idx), sep = "_t")
        local.DBN.input.data[, col.name] <- input.data.discr.3D[time.trans.idx, node.name, ]
      }
      rm(node.name)
      
      ## Fill in the data of the target node
      ## in 'local.DBN.input.data'
      col.name <- paste(tgt.node.name, as.character(time.trans.idx + auto.reg.order), sep = "_t")
      local.DBN.input.data[, col.name] <- input.data.discr.3D[(time.trans.idx + auto.reg.order), 
                                                              tgt.node.name, ]
      rm(col.name)
      ## End: Generate 2D 'local.DBN.input.data' from 'local.input.data.3D'
      
      # print('before data frame')
      # local.DBN.input.data <- as.data.frame(local.DBN.input.data)
      
      ## source(paste(init.path, 'learn_local_dbn.R', sep = '/'))
      ## Returns the list of predicted source nodes for the local DBN.
      local.dbn.pred.src.nodes <- LearnLocalDbnBnstruct(local.DBN.input.data, 
                                                        rep(num.discr.levels, ncol(local.DBN.input.data)), 
                                                        scoring.func)
      
      ## Assuming there are > 1 time points. Otherwise, 'local.unrolled.DBN.adj.submatrix' would become a vector.
      ## Initialize the sub-matrix with zeroes.
      local.unrolled.DBN.adj.submatrix <- matrix(0, 
                                                 nrow = length(local.unrolled.DBN.src.node.names), 
                                                 ncol = 1, 
                                                 dimnames = c(list(local.unrolled.DBN.src.node.names), 
                                                              list(local.unrolled.DBN.tgt.node.name)))
      
      # print(local.unrolled.DBN.adj.submatrix)
      
      ## Only if the predicted source node set is non-empty
      if (length(local.dbn.pred.src.nodes) > 0) {
        for (src.node.name in local.dbn.pred.src.nodes) {
          local.unrolled.DBN.adj.submatrix[src.node.name, local.unrolled.DBN.tgt.node.name] <- 1
        }
        rm(src.node.name)
      }
      
      #---------------------------------
      # End: Local Unrolled DBN struct learning
      #---------------------------------
      
      print(local.unrolled.DBN.adj.submatrix)
      
      # Return value for each 'foreach' iteration 
      local.unrolled.DBN.adj.submatrix
    }
  print('End of foreach loops')
  
  # save(local.unrolled.DBN.adj.matrix.list, file = paste(getwd(), 'asset/local.unrolled.DBN.adj.matrix.list.RData', sep = '/'))
  
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

###############################################################################################################################
## Goal: Unrolled DBN structure learning with Markov Order 1.
## This function represents a multi-core implementation of function
## LearnDbnStructMo1Layer3dSerDeg1_v2_Lite(), also defined in this R script.
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
LearnDbnStructMo1Layer3dParDeg1_v2_Lite <- function(input.data.discr.3D, 
                                               mi.net.adj.matrix, 
                                               num.discr.levels, 
                                               num.nodes, 
                                               num.timepts, 
                                               node.names, 
                                               clr.algo, 
                                               auto.reg.order, 
                                               scoring.func, 
                                               max.num.cores, 
                                               output.filename, 
                                               init.path) {
  #---------------------------------
  # Begin: Loading the Packages
  #---------------------------------
  library(bnlearn)
  # library(ggm)
  library(foreach)
  
  library(doParallel)
  #---------------------------------
  # End: Loading the Packages
  #---------------------------------
  
  num.time.trans <- (num.timepts - 1)

  ## Start and register a parallel backend for parallel computing
  
  ## Available number of cores in the current system
  avl.num.cores <- (parallel::detectCores() - 1)
  
  ## 10 cores to be used for grni server 
  num.cores <- min(max.num.cores, num.nodes, avl.num.cores)
  cl <- parallel::makeCluster(num.cores, outfile = output.filename)
  rm(num.cores)
  doParallel::registerDoParallel(cl)
  
  ## '.verbose = TRUE' is used for debugging
  ## 'when(sum(mi.net.adj.matrix[, tgt.node.idx]) > 0' means when central node has at least one neighbour in the mutual info net
  ## Use %do% amd %dopar% for serial and parallel computing, resp.
  local.unrolled.DBN.adj.matrix.list <- 
    foreach::foreach(tgt.node.idx = 1:num.nodes, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:% 
    # when(sum(mi.net.adj.matrix[, tgt.node.idx]) != 0) %:%
    foreach::foreach(time.trans.idx = 1:num.time.trans, .packages = c('foreach', 'bnstruct'), .verbose = TRUE) %:%
    when(sum(mi.net.adj.matrix[, tgt.node.idx]) > 0) %dopar% {
      
      ## Already done in TGS.R
      # source(paste(init.path, 'learn_local_dbn.R', sep = '/'))
      # source(paste(init.path, 'score_bn.R', sep = '/'))
      # 
      # if (.Platform$OS.type == 'windows') {
      #   dyn.load(paste(init.path, 'src/windows/score_bn.dll', sep = '/'))
      # } else if (.Platform$OS.type == 'unix') {
      #   dyn.load(paste(init.path, 'src/unix/score_bn.so', sep = '/'))
      # }
      
      tgt.node.name <- rownames(mi.net.adj.matrix)[tgt.node.idx]
      
      # List names of the target node's neighbours in mi.net.adj.matrix
      nbr.names <- c()
      
      ## Just one neighbour
      if (sum(mi.net.adj.matrix[, tgt.node.idx]) == 1) {
        for (nbr.idx in 1:nrow(mi.net.adj.matrix)) {
          if (mi.net.adj.matrix[nbr.idx, tgt.node.idx] == 1) {
            nbr.names <- rownames(mi.net.adj.matrix)[nbr.idx]
            break
          }
        }
      } else if (sum(mi.net.adj.matrix[, tgt.node.idx]) > 1) {
        ## Multiple neighbours
        
        nbr.names <- rownames(mi.net.adj.matrix[which(mi.net.adj.matrix[, tgt.node.idx] == 1),])
      }
      
      local.net.node.names <- c()
      if (clr.algo == 'CLR2.1') {
        local.net.node.names <- nbr.names
      } else if (clr.algo %in% c('CLR', 'CLR1.2', 'CLR2', 'CLR10')) {
        local.net.node.names <- c(tgt.node.name, nbr.names)
      }
      
      #---------------------------------
      # Begin: Local Unrolled DBN struct learning
      #---------------------------------
      
      ## Begin: Generate 2D 'local.DBN.input.data' from 'local.input.data.3D'
      # Say, the target node is v1 and the local nodes are {v1, v2, v3}, 
      # and 'time.trans.idx' = 1 and 'auto.reg.order' = 1.
      # Then 'local.DBN.input.data.var.names' contains {v1_t1, v2_t1, v3_t1, v1_t2}.
      # In that case, 'local.DBN.input.data' is a matrix which
      # will contain columns corr. to elements 
      # in 'local.DBN.input.data.var.names' and
      # rows corr. to different samples.
      
      num.samples <- dim(input.data.discr.3D)[3]
      
      local.DBN.input.data <- matrix(0, nrow = num.samples, ncol = (length(local.net.node.names) + 1))
      
      local.unrolled.DBN.src.node.names <- c()
      for (var.name in local.net.node.names) {
        local.unrolled.DBN.src.node.names <- c(local.unrolled.DBN.src.node.names,
                                               paste(var.name, as.character(time.trans.idx), 
                                                     sep = "_t"))
      }
      rm(var.name)
      
      local.unrolled.DBN.tgt.node.name <- paste(tgt.node.name, 
                                                as.character(time.trans.idx + auto.reg.order), 
                                                sep = "_t")
      
      local.DBN.input.data.var.names <- c(local.unrolled.DBN.src.node.names, 
                                          local.unrolled.DBN.tgt.node.name)
      
      colnames(local.DBN.input.data) <- local.DBN.input.data.var.names
      
      ## Fill in the data of the candidate source nodes
      ## in 'local.DBN.input.data'
      for (node.name in local.net.node.names) {
        col.name <- paste(node.name, as.character(time.trans.idx), sep = "_t")
        local.DBN.input.data[, col.name] <- input.data.discr.3D[time.trans.idx, node.name, ]
      }
      rm(node.name)
      
      ## Fill in the data of the target node
      ## in 'local.DBN.input.data'
      col.name <- paste(tgt.node.name, as.character(time.trans.idx + auto.reg.order), sep = "_t")
      local.DBN.input.data[, col.name] <- input.data.discr.3D[(time.trans.idx + auto.reg.order), 
                                                              tgt.node.name, ]
      rm(col.name)
      ## End: Generate 2D 'local.DBN.input.data' from 'local.input.data.3D'
      
      # print('before data frame')
      # local.DBN.input.data <- as.data.frame(local.DBN.input.data)
      
      source(paste(init.path, 'learn_local_dbn.R', sep = '/'))
      
      ## source(paste(init.path, 'learn_local_dbn.R', sep = '/'))
      ## Returns the list of predicted source nodes for the local DBN.
      local.dbn.pred.src.nodes <- LearnLocalDbnBnstruct(local.DBN.input.data, 
                                                        rep(num.discr.levels, ncol(local.DBN.input.data)), 
                                                        scoring.func)
      
      ## Assuming there are > 1 time points. Otherwise, 'local.unrolled.DBN.adj.submatrix' would become a vector.
      ## Initialize the sub-matrix with zeroes.
      local.unrolled.DBN.adj.submatrix <- matrix(0, 
                                                 nrow = length(local.unrolled.DBN.src.node.names), 
                                                 ncol = 1, 
                                                 dimnames = c(list(local.unrolled.DBN.src.node.names), 
                                                              list(local.unrolled.DBN.tgt.node.name)))
      
      # print(local.unrolled.DBN.adj.submatrix)
      
      ## Only if the predicted source node set is non-empty
      if (length(local.dbn.pred.src.nodes) > 0) {
        for (src.node.name in local.dbn.pred.src.nodes) {
          local.unrolled.DBN.adj.submatrix[src.node.name, local.unrolled.DBN.tgt.node.name] <- 1
        }
        rm(src.node.name)
      }
      
      #---------------------------------
      # End: Local Unrolled DBN struct learning
      #---------------------------------
      
      print(local.unrolled.DBN.adj.submatrix)
      
      # Return value for each 'foreach' iteration 
      local.unrolled.DBN.adj.submatrix
    }
  print('End of foreach loops')
  
  # save(local.unrolled.DBN.adj.matrix.list, file = paste(getwd(), 'asset/local.unrolled.DBN.adj.matrix.list.RData', sep = '/'))
  
  ## Shut down the cluster
  parallel::stopCluster(cl)

  
  
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