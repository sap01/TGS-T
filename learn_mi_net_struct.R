ComputEntropy <- function(input.data)
{
  n <- ncol(input.data)
  entropy.matrix <- matrix(0, nrow = 1, ncol = n)

  for(i in 1:n)
  {
     c1 <-  var(input.data[,i])
     entropy.matrix[i] <- .5*log((2*pi*exp(1))*c1)
  }
  return(entropy.matrix) 
}

############################################################################################

LearnMiNetStructZstat <- function(mut.info.matrix, mi.net.adj.matrix, entropy.matrix, alpha)
{
  n <- ncol(mut.info.matrix)
  
  threshold <- qnorm(1-(alpha/2))
  
  for(i in 1:(n-1))
  {
    for(j in (i+1):n)
    { 
      value <- mut.info.matrix[i,j]/(entropy.matrix[i]+entropy.matrix[j])
       fisher_transform <- .5*log((1+value)/(1-value))
        value1 <- sqrt(n - 3) * abs(fisher_transform)
        
        if(value1 > threshold)
        {
          mi.net.adj.matrix[i,j] <- 1
          mi.net.adj.matrix[j,i] <- 1
        }
    }
  }
  return(mi.net.adj.matrix)
}

############################################################################################

LearnMiNetStructRowMedian <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes)
{
  for (rowIdx in 1:num.nodes)
  {
    threshold <- median(mut.info.matrix[rowIdx, -rowIdx])

    for (colIdx in 1:num.nodes)
    {
      if ((colIdx != rowIdx) & (mut.info.matrix[rowIdx, colIdx] >= threshold))
      {
        mi.net.adj.matrix[rowIdx, colIdx] <- 1          
      }
    }
  }
  
  return(mi.net.adj.matrix)
}

############################################################################################
## Goal: Learn CLR net. Replace all non-zero edge weights with 1.

library(minet)

LearnMiNetStructClr <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes)
{
  mi.net.adj.matrix.wt <- minet::clr(mut.info.matrix) # weighted adj matrix
  
  # Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
  mi.net.adj.matrix.wt[is.nan(mi.net.adj.matrix.wt)] <- 0
  
  # writeLines('\n mi.net.adj.matrix.wt = \n')
  # print(mi.net.adj.matrix.wt)
  save(mi.net.adj.matrix.wt, file = paste(getwd(), 'asset/mi.net.adj.matrix.wt.RData', sep = '/'))
  
  for (rowIdx in 1:num.nodes)
  {
    for (colIdx in 1:num.nodes)
    {
      if (mi.net.adj.matrix.wt[rowIdx, colIdx] != 0)
      {
        mi.net.adj.matrix[rowIdx, colIdx] <- 1
      }
    }
  }
  
  return(mi.net.adj.matrix)
}

############################################################################################
## Goal: Learn CLR net from a given discretized dataset. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.

LearnClrNetFromDiscrData <- function(input.data.discr, num.nodes, node.names)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  library(minet)
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
    
  # Initialize mutual information matrix with zeroes
  mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, dimnames = c(list(node.names), list(node.names)))
  
  ## Build mutual information matrix
  for (col.idx in 1:(num.nodes - 1)) {
    for (col.idx.2 in (colIdx + 1):num.nodes) {
      mut.info <- computeCmi(input.data.discr[, col.idx], input.data.discr[, col.idx.2])
      mut.info.matrix[col.idx, col.idx.2] <- mut.info
      mut.info.matrix[col.idx.2, col.idx] <- mut.info
    }
    rm(col.idx.2)
  }
  rm(col.idx)
    
  ## Estimate weighted adj matrix of the CLR net
  ## from the mutual info matrix
  mi.net.adj.matrix.wt <- minet::clr(mut.info.matrix)
    
    ## Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
    mi.net.adj.matrix.wt.curr.series[is.nan(mi.net.adj.matrix.wt.curr.series)] <- 0
    
    ## 'mi.net.adj.matrix.wt.curr.series' has non-neg values.
    ## Summing up all time-series-specific 'mi.net.adj.matrix.wt.curr.series' matrices
    ## produces a ('num.nodes' by 'num.nodes') matrix of non-neg
    ## values, namely 'mi.net.adj.matrix.wt'. Later 'mi.net.adj.matrix.wt'
    ## will be divided by 'num.time.series', thus producing the 
    ## arthmetic mean of all time-series-specific 
    ## 'mi.net.adj.matrix.wt.curr.series' matrices.
    mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt + mi.net.adj.matrix.wt.curr.series)
  
  ## Arthmetic mean of all time-series-specific 
  ## 'mi.net.adj.matrix.wt' matrices.
  mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt / num.time.series)
  
  save(mi.net.adj.matrix.wt, file = paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))
  
  ##############################################################
  ## Begin:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################
  
  ## For each target node
  for (col.idx in 1:num.nodes) {
    
    ## Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]
    
    ## Count number of neighbours having positive edge weight
    num.nbrs <- length(edge.wts[edge.wts > 0])
    
    if (num.nbrs >= max.fanin) {
      
      ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      ## Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
      
      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
      
    } else if (num.nbrs < max.fanin) {
      
      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  rm(col.idx)
  
  ##############################################################
  ## End:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################
  
  return(mi.net.adj.matrix)
  
}

#######################################################################################################
## Goal: Learn CLR net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
##
LearnClrNetMfi <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes, max.fanin, output.dirname) {
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  library(minet)
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  mi.net.adj.matrix.wt <- minet::clr(mut.info.matrix) # weighted adj matrix
  
  # Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
  mi.net.adj.matrix.wt[is.nan(mi.net.adj.matrix.wt)] <- 0
  
  # writeLines('\n mi.net.adj.matrix.wt = \n')
  # print(mi.net.adj.matrix.wt)
  save(mi.net.adj.matrix.wt, file = paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))

  # For each target node
  for (col.idx in 1:num.nodes) {
    # Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]
    
    # Count number of neighbours having positive edge weight
    num.nbrs <- length(edge.wts[edge.wts > 0])
    
    if (num.nbrs >= max.fanin) {
      # Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      # Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
      
      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
      
      ## The following line is not required since 'mi.net.adj.matrix' is initialized
      ## with all zeroes
      # mi.net.adj.matrix[-(valid.nbrs), col.idx] <- 0
    } else if (num.nbrs < max.fanin) {
      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  
  return(mi.net.adj.matrix)
}
############################################################################################

#######################################################################################################
## Goal: Learn CLR1.1 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
##
LearnClr1dot1NetMfi <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes, max.fanin, output.dirname) {
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  library(minet)
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  mi.net.adj.matrix.wt <- minet::clr(mut.info.matrix) # weighted adj matrix
  
  # Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
  mi.net.adj.matrix.wt[is.nan(mi.net.adj.matrix.wt)] <- 0
  
  # writeLines('\n mi.net.adj.matrix.wt = \n')
  # print(mi.net.adj.matrix.wt)
  save(mi.net.adj.matrix.wt, file = paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))
  
  # For each target node
  for (col.idx in 1:num.nodes) {
    # Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]
    
    # Count number of neighbours having positive edge weight
    num.nbrs <- length(edge.wts[edge.wts > 0])
    
    if (num.nbrs >= max.fanin) {
      # Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      # Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
      
      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
      
      ## The following line is not required since 'mi.net.adj.matrix' is initialized
      ## with all zeroes
      # mi.net.adj.matrix[-(valid.nbrs), col.idx] <- 0
    } else if (num.nbrs < max.fanin) {
      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  
  return(mi.net.adj.matrix)
}
############################################################################################

#######################################################################################################
## Goal: Learn CLR1.2 net. This is exactly same as CLR net.
## Therefore, this function is equivalent to function LearnClrNetMfi().
## The only difference is that in LearnClrNetMfi(), 
## the weighted CLR net adjacency matrix is estimated using the minet
## package's clr() function, whereas in this function, the estimation
## is implemented without depending on any external package.
##
LearnClr1dot2NetMfi <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes, max.fanin, output.dirname, init.path) {
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  # none
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Load external functions
  source(paste(init.path, 'learn_clr.R', sep = '/'))
  
  # dyn.load(paste(init.path, 'src/clr.so', sep = '/'))
  if (.Platform$OS.type == 'windows') {
    ## TODO(sap)
    # dyn.load(paste(init.path, 'src/windows/clr.dll', sep = '/'))
  } else if (.Platform$OS.type == 'unix') {
    dyn.load(paste(init.path, 'src/unix/clr.so', sep = '/'))
  }
  
  ## source(paste(init.path, 'learn_clr.R', sep = '/'))
  mi.net.adj.matrix.wt <- LearnClr(mut.info.matrix, 'CLR1.2') # weighted adj matrix
  
  # Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
  mi.net.adj.matrix.wt[is.nan(mi.net.adj.matrix.wt)] <- 0
  
  # writeLines('\n mi.net.adj.matrix.wt = \n')
  # print(mi.net.adj.matrix.wt)
  save(mi.net.adj.matrix.wt, file = paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))
  
  # For each target node
  for (col.idx in 1:num.nodes) {
    # Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]
    
    # Count number of neighbours having positive edge weight
    num.nbrs <- length(edge.wts[edge.wts > 0])
    
    if (num.nbrs >= max.fanin) {
      # Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      # Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
      
      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
      
      ## The following line is not required since 'mi.net.adj.matrix' is initialized
      ## with all zeroes
      # mi.net.adj.matrix[-(valid.nbrs), col.idx] <- 0
    } else if (num.nbrs < max.fanin) {
      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  
  return(mi.net.adj.matrix)
}
############################################################################################

#######################################################################################################
## Goal: Learn CLR1.3 net. It is same as learning the CLR net except
## when at least one among zi or zj is non-positive, the CLR1.3 edge weight is zero.
##
LearnClr1dot3NetMfi <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes, max.fanin, 
                                output.dirname, init.path) {
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  # none
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Load external functions
  source(paste(init.path, 'learn_clr.R', sep = '/'))
  
  # dyn.load(paste(init.path, 'src/clr.so', sep = '/'))
  if (.Platform$OS.type == 'windows') {
    ## TODO(sap)
    # dyn.load(paste(init.path, 'src/windows/clr.dll', sep = '/'))
  } else if (.Platform$OS.type == 'unix') {
    dyn.load(paste(init.path, 'src/unix/clr.so', sep = '/'))
  }
  
  ## source(paste(init.path, 'learn_clr.R', sep = '/'))
  mi.net.adj.matrix.wt <- LearnClr(mut.info.matrix, 'CLR1.3') # weighted adj matrix
  
  # Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
  mi.net.adj.matrix.wt[is.nan(mi.net.adj.matrix.wt)] <- 0
  
  # writeLines('\n mi.net.adj.matrix.wt = \n')
  # print(mi.net.adj.matrix.wt)
  save(mi.net.adj.matrix.wt, file = paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))
  
  # For each target node
  for (col.idx in 1:num.nodes) {
    # Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]
    
    # Count number of neighbours having positive edge weight
    num.nbrs <- length(edge.wts[edge.wts > 0])
    
    if (num.nbrs >= max.fanin) {
      # Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      # Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
      
      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
      
      ## The following line is not required since 'mi.net.adj.matrix' is initialized
      ## with all zeroes
      # mi.net.adj.matrix[-(valid.nbrs), col.idx] <- 0
    } else if (num.nbrs < max.fanin) {
      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  
  return(mi.net.adj.matrix)
}
############################################################################################

############################################################################################
## Goal: Learn CLR2 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
##
LearnClr2NetMfi <- function(input.data.discr, num.nodes, node.names, num.timepts, 
                            max.fanin, output.dirname, mi.net.adj.matrix)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  library(minet)
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Initialize weighted adjacency matrix of the mutual information network
  mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(node.names), list(node.names)))
  
  ## Total number of time series 
  ## = number of measurements (replicates) per time pt
  ## = (total number of measurements / number of time pts).
  num.time.series <- (nrow(input.data.discr) / num.timepts)
  
  for (time.series.idx in 1:num.time.series) {
    
    ## First time point of the current time series
    first.time.pt.curr.series <- (((time.series.idx - 1) * num.timepts) + 1)
    
    ## Last time point of the current time series
    last.time.pt.curr.series <- (time.series.idx * num.timepts)
    
    ## Discretized data of the current time series
    input.data.discr.curr.series <- input.data.discr[first.time.pt.curr.series:last.time.pt.curr.series, ]
    rm(first.time.pt.curr.series, last.time.pt.curr.series)
    
    # Initialize mutual information matrix with zeroes
    mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, dimnames = c(list(node.names), list(node.names)))
    
    ## Build mutual information matrix
    for (col.idx in 1:(num.nodes - 1)) {
      for (col.idx.2 in (col.idx + 1):num.nodes) {
        
        ## compute_cmi.R
        mut.info <- ComputeCmi(input.data.discr.curr.series[, col.idx], input.data.discr.curr.series[, col.idx.2])
        
        mut.info.matrix[col.idx, col.idx.2] <- mut.info
        mut.info.matrix[col.idx.2, col.idx] <- mut.info
      }
      rm(col.idx.2)
    }
    rm(col.idx)
    
    ## Estimate weighted adj matrix of the CLR net
    ## corr. to the current time series
    mi.net.adj.matrix.wt.curr.series <- minet::clr(mut.info.matrix)
    
    ## Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
    mi.net.adj.matrix.wt.curr.series[is.nan(mi.net.adj.matrix.wt.curr.series)] <- 0
    
    ## 'mi.net.adj.matrix.wt.curr.series' has non-neg values.
    ## Summing up all time-series-specific 'mi.net.adj.matrix.wt.curr.series' matrices
    ## produces a ('num.nodes' by 'num.nodes') matrix of non-neg
    ## values, namely 'mi.net.adj.matrix.wt'. Later 'mi.net.adj.matrix.wt'
    ## will be divided by 'num.time.series', thus producing the 
    ## arthmetic mean of all time-series-specific 
    ## 'mi.net.adj.matrix.wt.curr.series' matrices.
    mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt + mi.net.adj.matrix.wt.curr.series)
    
  }
  rm(time.series.idx)
  
  ## Arthmetic mean of all time-series-specific 
  ## 'mi.net.adj.matrix.wt' matrices.
  mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt / num.time.series)
  
  save(mi.net.adj.matrix.wt, file = paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))
  
  ##############################################################
  ## Begin:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################
  
  ## For each target node
  for (col.idx in 1:num.nodes) {
    
    ## Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]
    
    ## Count number of neighbours having positive edge weight
    num.nbrs <- length(edge.wts[edge.wts > 0])
    
    if (num.nbrs >= max.fanin) {
      
      ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      ## Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
      
      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
      
    } else if (num.nbrs < max.fanin) {
      
      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  rm(col.idx)
  
  ##############################################################
  ## End:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################
  
  return(mi.net.adj.matrix)

}

############################################################################################

############################################################################################
## Goal: Learn CLR2.1 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
##
LearnClrNetMfiVer2.1 <- function(input.data.discr, num.nodes, node.names, num.timepts, 
                                 max.fanin, output.dirname, mi.net.adj.matrix)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Initialize weighted adjacency matrix of the mutual information network
  mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                 dimnames = c(list(node.names), list(node.names)))
  
  ## Total number of time series 
  ## = number of measurements (replicates) per time pt
  ## = (total number of measurements / number of time pts).
  num.time.series <- (nrow(input.data.discr) / num.timepts)
  
  for (time.series.idx in 1:num.time.series) {
    
    ## First time point of the current time series
    first.time.pt.curr.series <- (((time.series.idx - 1) * num.timepts) + 1)
    
    ## Last time point of the current time series
    last.time.pt.curr.series <- (time.series.idx * num.timepts)
    
    ## Discretized data of the current time series
    input.data.discr.curr.series <- input.data.discr[first.time.pt.curr.series:last.time.pt.curr.series, ]
    rm(first.time.pt.curr.series, last.time.pt.curr.series)
    
    # Initialize mutual information matrix with zeroes
    mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, dimnames = c(list(node.names), list(node.names)))
    
    ## Build assymetric mutual information matrix
    for (row.idx in 1:num.nodes) {
      for (col.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        mut.info <- ComputeCmi(input.data.discr.curr.series[1:(num.timepts-1), row.idx], 
                               input.data.discr.curr.series[2:num.timepts, col.idx])
        
        mut.info.matrix[row.idx, col.idx] <- mut.info
      }
      rm(col.idx)
    }
    rm(row.idx)
    
    # print('mut.info.matrix')
    # print(time.series.idx)
    # print(mut.info.matrix)
    
    ## Initialize weighted adj matrix of the CLR net
    ## corr. to the current time series
    mi.net.adj.matrix.wt.curr.series <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                               dimnames = c(list(node.names), list(node.names)))
    ## Compute 'mi.net.adj.matrix.wt.curr.series[src.node.idx, tgt.node.idx]'
    for (src.node.idx in 1:num.nodes) {
      for (tgt.node.idx in 1:num.nodes) {
        
        if (mut.info.matrix[src.node.idx, tgt.node.idx] == .Machine$double.xmax) {
          ## Perfectly correlated
          
          mi.net.adj.matrix.wt.curr.series[src.node.idx, tgt.node.idx] <- .Machine$double.xmax
          
        } else {
          
          src.clr.mean <- mean(mut.info.matrix[src.node.idx, ])
          # if (src.clr.mean > .Machine$double.xmax) {
          #   src.clr.mean <- .Machine$double.xmax
          # }
          
          src.clr.sd <- sd(mut.info.matrix[src.node.idx, ])
          # if (src.clr.sd > .Machine$double.xmax) {
          #   src.clr.sd <- .Machine$double.xmax
          # }
          
          tgt.clr.mean <- mean(mut.info.matrix[, tgt.node.idx])
          # if (tgt.clr.mean > .Machine$double.xmax) {
          #   tgt.clr.mean <- .Machine$double.xmax
          # }
          
          tgt.clr.sd <- sd(mut.info.matrix[, tgt.node.idx])
          # if (tgt.clr.sd > .Machine$double.xmax) {
          #   tgt.clr.sd <- .Machine$double.xmax
          # }
          
          if ((src.clr.mean >= .Machine$double.xmax) | (src.clr.sd >= .Machine$double.xmax) | 
              (tgt.clr.mean >= .Machine$double.xmax) | (tgt.clr.sd >= .Machine$double.xmax)) {
            
            ## There exists far better candidate source node(s)
            mi.net.adj.matrix.wt.curr.series[src.node.idx, tgt.node.idx] <- 0
            
          } else {
            
            tmp <- 0
            if (src.clr.sd != 0) {
              tmp <- (mut.info.matrix[src.node.idx, tgt.node.idx] - src.clr.mean) / src.clr.sd
            }
            z.src <- max(0, tmp)
            
            ## Re-initilize in case 'tgt.clr.sd == 0'
            tmp <- 0
            if (tgt.clr.sd != 0) {
              tmp <- (mut.info.matrix[src.node.idx, tgt.node.idx] - tgt.clr.mean) / tgt.clr.sd
            }
            z.tgt <- max(0, tmp)
            
            rm(tmp)
            
            clr.edge.wt <- ((z.tgt)^2 + (z.src)^2)
            clr.edge.wt <- sqrt(clr.edge.wt)
            
            mi.net.adj.matrix.wt.curr.series[src.node.idx, tgt.node.idx] <- clr.edge.wt
          }
        }
        
      }
      rm(tgt.node.idx)
    }
    rm(src.node.idx)
    
    ## 'mi.net.adj.matrix.wt.curr.series' has non-neg values.
    ## Summing up all time-series-specific 'mi.net.adj.matrix.wt.curr.series' matrices
    ## produces a ('num.nodes' by 'num.nodes') matrix of non-neg
    ## values, namely 'mi.net.adj.matrix.wt'. Later 'mi.net.adj.matrix.wt'
    ## will be divided by 'num.time.series', thus producing the 
    ## arthmetic mean of all time-series-specific 
    ## 'mi.net.adj.matrix.wt.curr.series' matrices.
    mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt + mi.net.adj.matrix.wt.curr.series)
    
    ## '.Machine$double.xmax + .Machine$double.xmax = Inf'
    ## Prevent 'Inf'
    mi.net.adj.matrix.wt[mi.net.adj.matrix.wt > .Machine$double.xmax] <- .Machine$double.xmax
    
  }
  rm(time.series.idx)
  
  ## Arthmetic mean of all time-series-specific 
  ## 'mi.net.adj.matrix.wt' matrices.
  mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt / num.time.series)
  
  save(mi.net.adj.matrix.wt, file = paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))
  
  ##############################################################
  ## Begin:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################
  
  ## For each target node
  for (col.idx in 1:num.nodes) {
    
    ## Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]
    
    ## Count number of neighbours having positive edge weight
    num.nbrs <- length(edge.wts[edge.wts > 0])
    
    if (num.nbrs >= max.fanin) {
      
      ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      ## Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
      
      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
      
    } else if (num.nbrs < max.fanin) {
      
      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  rm(col.idx)
  
  ##############################################################
  ## End:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################
  
  return(mi.net.adj.matrix)
  
}

############################################################################################

############################################################################################
## Goal: Learn CLR3 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
## CLR3 edge weight = (z.parent^2 + z.tgt^2)^0.5 where
## z.parent = max(0, (MI(candidate.parent, tgt) - mean(MI(candidate.parent, .)) / sd(MI(candidate.parent, .)))) and
## z.tgt = max(0, (MI(candidate.parent, tgt) - mean(MI(., tgt)) / sd(MI(., tgt)))).
## If the tgt belongs to time point t_p, then the candidate parent must belong to time point t_(p + 1).
##
LearnClr3NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                            max.fanin, mi.net.adj.matrix.list) {
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]
    
    candidate.parent.node.names <- c()
    candidate.tgt.node.names <- c()
    for (curr.node.name in node.names) {
      parent.full.name <- paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- c(candidate.parent.node.names, parent.full.name)
      rm(parent.full.name)
      
      tgt.full.name <- paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- c(candidate.tgt.node.names, tgt.full.name)
      rm(tgt.full.name)
    }
    rm(curr.node.name)
    
    ## Initialize mutual information matrix with zeroes
    mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(candidate.parent.node.names), 
                                           list(candidate.tgt.node.names)))
    
    #################################################################################
    # ## Build mutual information matrix
    # for (col.idx in 1:(num.nodes - 1)) {
    #   for (col.idx.2 in (col.idx + 1):num.nodes) {
    #     
    #     ## compute_cmi.R
    #     ## (dim1 == 1) => time.pt.idx
    #     ## (dim1 == 2) => (time.pt.idx + 1)
    #     mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, col.idx, ], 
    #                            input.data.discr.3D.curr.ival[2, col.idx.2, ])
    #     
    #     mut.info.matrix[col.idx, col.idx.2] <- mut.info
    #     mut.info.matrix[col.idx.2, col.idx] <- mut.info
    #   }
    #   rm(col.idx.2)
    # }
    # rm(col.idx)
    #################################################################################
    ## Build mutual information matrix
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[2, tgt.idx, ])
        
        ## Mutual info matrix is asymmetric in this case
        mut.info.matrix[parent.idx, tgt.idx] <- mut.info
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    #################################################################################
    
    ## Initialize weighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(candidate.parent.node.names), 
                                           list(candidate.tgt.node.names)))
    
    candidate.parent.mean.sd <- matrix(0, nrow = num.nodes, ncol = 2)
    rownames(candidate.parent.mean.sd) <- candidate.parent.node.names
    colnames(candidate.parent.mean.sd) <- c('clr.mean', 'clr.sd')
    
    ## Calculate sample mean and sample standard deviation of the given nodes.
    ## It is calculated acc. to the logic in function 'clr()' in
    ## R package 'minet' (version 3.36.0).
    ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
    ## the calculation. Here, the same logic is re-implemented in R.
    ## Since the in-built 'mean()' and 'sd()' functions in 
    ## R version 3.3.2 follows the exact same logic, therefore, the
    ## re-implementation is straight-forward.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.mean.sd[parent.name, 'clr.mean'] <- mean(mut.info.matrix[parent.name, ])
      
      ## var <- 0
      ## for (each sample) {
      ##  sd <- (mean - sample val)
      ##  var <- var + (sd^2)
      ## }
      ## var <- var / (n -1) ## where n = number of samples
      ## sd <- sqrt(var)
      candidate.parent.mean.sd[parent.name, 'clr.sd'] <- sd(mut.info.matrix[parent.name, ])
    }
    rm(parent.name)
    
    for (tgt.node.name in candidate.tgt.node.names) {
      
      tgt.clr.mean <- mean(mut.info.matrix[, tgt.node.name])
      tgt.clr.sd <- sd(mut.info.matrix[, tgt.node.name])
      
      ## Edge weights of the CLR net are calculated acc. to
      ## 'minet::clr()'
      for (candidate.parent.name in candidate.parent.node.names) {
        
        tmp <- 0
        if (tgt.clr.sd != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - tgt.clr.mean) / tgt.clr.sd
        }
        z.tgt <- max(0, tmp)
        
        tmp <- 0
        if (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - 
                    candidate.parent.mean.sd[candidate.parent.name, 'clr.mean']) / 
            candidate.parent.mean.sd[candidate.parent.name, 'clr.sd']
        }
        z.parent <- max(0, tmp)

        rm(tmp)
        
        clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
        clr.edge.wt <- sqrt(clr.edge.wt)
        rm(z.tgt, z.parent)
        
        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      rm(candidate.parent.name)
      
    }
    rm(tgt.node.name)
    
    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                   dimnames = c(list(candidate.parent.node.names), 
                                                list(candidate.tgt.node.names)))
    
    ## For each target node
    for (col.idx in 1:num.nodes) {
      
      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]
      
      ## Count number of neighbours having positive edge weight
      num.nbrs <- length(edge.wts[edge.wts > 0])
      
      if (num.nbrs >= max.fanin) {
        
        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
        
        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
        
      } else if (num.nbrs < max.fanin) {
        
        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    rm(col.idx)
    
    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    
    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  rm(time.pt.idx)
  
  return(mi.net.adj.matrix.list)
}

############################################################################################

############################################################################################
## Goal: Learn CLR4 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
## If (z.parent and z.tgt both are greater than zero) then
## CLR4 edge weight = (z.parent^2 + z.tgt^2)^0.5 where
## z.parent = max(0, (MI(candidate.parent, tgt) - mean(MI(candidate.parent, .)) / sd(MI(candidate.parent, .)))) and
## z.tgt = max(0, (MI(candidate.parent, tgt) - mean(MI(., tgt)) / sd(MI(., tgt)))).
## Note that the tgt belongs to time point t_p and the candidate parent belongs to time point t_(p + 1).
## Else if (not both z.parent and z.tgt are greater than zero) then
## CLR4 edge weight = zero.
##
LearnClr4NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                            max.fanin, mi.net.adj.matrix.list)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]
    
    candidate.parent.node.names <- c()
    candidate.tgt.node.names <- c()
    for (curr.node.name in node.names) {
      parent.full.name <- paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- c(candidate.parent.node.names, parent.full.name)
      rm(parent.full.name)
      
      tgt.full.name <- paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- c(candidate.tgt.node.names, tgt.full.name)
      rm(tgt.full.name)
    }
    rm(curr.node.name)
    
    ## Initialize mutual information matrix with zeroes
    mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(candidate.parent.node.names), 
                                           list(candidate.tgt.node.names)))
    
    #################################################################################
    # ## Build mutual information matrix
    # for (col.idx in 1:(num.nodes - 1)) {
    #   for (col.idx.2 in (col.idx + 1):num.nodes) {
    #     
    #     ## compute_cmi.R
    #     ## (dim1 == 1) => time.pt.idx
    #     ## (dim1 == 2) => (time.pt.idx + 1)
    #     mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, col.idx, ], 
    #                            input.data.discr.3D.curr.ival[2, col.idx.2, ])
    #     
    #     mut.info.matrix[col.idx, col.idx.2] <- mut.info
    #     mut.info.matrix[col.idx.2, col.idx] <- mut.info
    #   }
    #   rm(col.idx.2)
    # }
    # rm(col.idx)
    #################################################################################
    ## Build mutual information matrix
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[2, tgt.idx, ])
        
        ## Mutual info matrix is asymmetric in this case
        mut.info.matrix[parent.idx, tgt.idx] <- mut.info
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    #################################################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                   dimnames = c(list(candidate.parent.node.names), 
                                                list(candidate.tgt.node.names)))
    
    candidate.parent.mean.sd <- matrix(0, nrow = num.nodes, ncol = 2)
    rownames(candidate.parent.mean.sd) <- candidate.parent.node.names
    colnames(candidate.parent.mean.sd) <- c('clr.mean', 'clr.sd')
    
    ## Calculate sample mean and sample standard deviation of the given nodes.
    ## It is calculated acc. to the logic in function 'clr()' in
    ## R package 'minet' (version 3.36.0).
    ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
    ## the calculation. Here, the same logic is re-implemented in R.
    ## Since the in-built 'mean()' and 'sd()' functions in 
    ## R version 3.3.2 follows the exact same logic, therefore, the
    ## re-implementation is straight-forward.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.mean.sd[parent.name, 'clr.mean'] <- mean(mut.info.matrix[parent.name, ])
      
      ## var <- 0
      ## for (each sample) {
      ##  sd <- (mean - sample val)
      ##  var <- var + (sd^2)
      ## }
      ## var <- var / (n -1) ## where n = number of samples
      ## sd <- sqrt(var)
      candidate.parent.mean.sd[parent.name, 'clr.sd'] <- sd(mut.info.matrix[parent.name, ])
    }
    rm(parent.name)
    
    for (tgt.node.name in candidate.tgt.node.names) {
      
      tgt.clr.mean <- mean(mut.info.matrix[, tgt.node.name])
      tgt.clr.sd <- sd(mut.info.matrix[, tgt.node.name])
      
      ## Edge weights of the CLR net are calculated acc. to
      ## 'minet::clr()'
      for (candidate.parent.name in candidate.parent.node.names) {
        
        tmp <- 0
        if (tgt.clr.sd != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - tgt.clr.mean) / tgt.clr.sd
        }
        z.tgt <- max(0, tmp)
        
        tmp <- 0
        if (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - 
                    candidate.parent.mean.sd[candidate.parent.name, 'clr.mean']) / 
            candidate.parent.mean.sd[candidate.parent.name, 'clr.sd']
        }
        z.parent <- max(0, tmp)
        
        rm(tmp)
        
        clr.edge.wt <- 0
        if ((z.tgt > 0) && (z.parent > 0)) {
          clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
          clr.edge.wt <- sqrt(clr.edge.wt)
        }
        
        rm(z.tgt, z.parent)
        
        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      rm(candidate.parent.name)
      
    }
    rm(tgt.node.name)
    
    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                dimnames = c(list(candidate.parent.node.names), 
                                             list(candidate.tgt.node.names)))
    
    ## For each target node
    for (col.idx in 1:num.nodes) {
      
      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]
      
      ## Count number of neighbours having positive edge weight
      num.nbrs <- length(edge.wts[edge.wts > 0])
      
      if (num.nbrs >= max.fanin) {
        
        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
        
        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
        
      } else if (num.nbrs < max.fanin) {
        
        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    rm(col.idx)
    
    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    
    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  rm(time.pt.idx)
  
  return(mi.net.adj.matrix.list)
}

############################################################################################

############################################################################################
## Goal: Learn CLR4.1 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
## If (z.parent and z.tgt both are greater than zero) then
## CLR4 edge weight = (z.parent^2 + z.tgt^2)^0.5 where
## z.parent = max(0, (MI(candidate.parent, tgt) - mean(MI(candidate.parent, .)) / sd(MI(candidate.parent, .))) - sd.thres) and
## z.tgt = max(0, (MI(candidate.parent, tgt) - mean(MI(., tgt)) / sd(MI(., tgt))) - sd.thres).
## Here, sd.thres is a user-defined parameter.
## Note that the tgt belongs to time point t_p and the candidate parent belongs to time point t_(p + 1).
## Else if (not both z.parent and z.tgt are greater than zero) then
## CLR4.1 edge weight = zero.
##
LearnClr4dot1NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                            max.fanin, mi.net.adj.matrix.list, sd.thres)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]
    
    candidate.parent.node.names <- c()
    candidate.tgt.node.names <- c()
    for (curr.node.name in node.names) {
      parent.full.name <- paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- c(candidate.parent.node.names, parent.full.name)
      rm(parent.full.name)
      
      tgt.full.name <- paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- c(candidate.tgt.node.names, tgt.full.name)
      rm(tgt.full.name)
    }
    rm(curr.node.name)
    
    ## Initialize mutual information matrix with zeroes
    mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(candidate.parent.node.names), 
                                           list(candidate.tgt.node.names)))
    
    #################################################################################
    # ## Build mutual information matrix
    # for (col.idx in 1:(num.nodes - 1)) {
    #   for (col.idx.2 in (col.idx + 1):num.nodes) {
    #     
    #     ## compute_cmi.R
    #     ## (dim1 == 1) => time.pt.idx
    #     ## (dim1 == 2) => (time.pt.idx + 1)
    #     mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, col.idx, ], 
    #                            input.data.discr.3D.curr.ival[2, col.idx.2, ])
    #     
    #     mut.info.matrix[col.idx, col.idx.2] <- mut.info
    #     mut.info.matrix[col.idx.2, col.idx] <- mut.info
    #   }
    #   rm(col.idx.2)
    # }
    # rm(col.idx)
    #################################################################################
    ## Build mutual information matrix
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[2, tgt.idx, ])
        
        ## Mutual info matrix is asymmetric in this case
        mut.info.matrix[parent.idx, tgt.idx] <- mut.info
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    #################################################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                   dimnames = c(list(candidate.parent.node.names), 
                                                list(candidate.tgt.node.names)))
    
    candidate.parent.mean.sd <- matrix(0, nrow = num.nodes, ncol = 2)
    rownames(candidate.parent.mean.sd) <- candidate.parent.node.names
    colnames(candidate.parent.mean.sd) <- c('clr.mean', 'clr.sd')
    
    ## Calculate sample mean and sample standard deviation of the given nodes.
    ## It is calculated acc. to the logic in function 'clr()' in
    ## R package 'minet' (version 3.36.0).
    ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
    ## the calculation. Here, the same logic is re-implemented in R.
    ## Since the in-built 'mean()' and 'sd()' functions in 
    ## R version 3.3.2 follows the exact same logic, therefore, the
    ## re-implementation is straight-forward.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.mean.sd[parent.name, 'clr.mean'] <- mean(mut.info.matrix[parent.name, ])
      
      ## var <- 0
      ## for (each sample) {
      ##  sd <- (mean - sample val)
      ##  var <- var + (sd^2)
      ## }
      ## var <- var / (n -1) ## where n = number of samples
      ## sd <- sqrt(var)
      candidate.parent.mean.sd[parent.name, 'clr.sd'] <- sd(mut.info.matrix[parent.name, ])
    }
    rm(parent.name)
    
    for (tgt.node.name in candidate.tgt.node.names) {
      
      tgt.clr.mean <- mean(mut.info.matrix[, tgt.node.name])
      tgt.clr.sd <- sd(mut.info.matrix[, tgt.node.name])
      
      ## Edge weights of the CLR net are calculated acc. to
      ## 'minet::clr()'
      for (candidate.parent.name in candidate.parent.node.names) {
        
        tmp <- 0
        if (tgt.clr.sd != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - tgt.clr.mean) / tgt.clr.sd
          tmp <- (tmp - sd.thres)
        }
        z.tgt <- max(0, tmp)
        
        tmp <- 0
        if (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - 
                    candidate.parent.mean.sd[candidate.parent.name, 'clr.mean']) / 
            candidate.parent.mean.sd[candidate.parent.name, 'clr.sd']
          tmp <- (tmp - sd.thres)
        }
        z.parent <- max(0, tmp)
        
        rm(tmp)
        
        clr.edge.wt <- 0
        if ((z.tgt > 0) && (z.parent > 0)) {
          clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
          clr.edge.wt <- sqrt(clr.edge.wt)
        }
        
        rm(z.tgt, z.parent)
        
        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      rm(candidate.parent.name)
      
    }
    rm(tgt.node.name)
    
    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                dimnames = c(list(candidate.parent.node.names), 
                                             list(candidate.tgt.node.names)))
    
    ## For each target node
    for (col.idx in 1:num.nodes) {
      
      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]
      
      ## Count number of neighbours having positive edge weight
      num.nbrs <- length(edge.wts[edge.wts > 0])
      
      if (num.nbrs >= max.fanin) {
        
        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
        
        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
        
      } else if (num.nbrs < max.fanin) {
        
        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    rm(col.idx)
    
    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    
    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  rm(time.pt.idx)
  
  return(mi.net.adj.matrix.list)
}

############################################################################################

############################################################################################
## Goal: Learn CLR4.2 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
## If (z.parent and z.tgt both are greater than zero) then
## CLR4 edge weight = (z.parent^2 + z.tgt^2)^0.5 where
## z.parent = max(0, (MI(candidate.parent, tgt) - mean(MI(candidate.parent, .)) / sd(MI(candidate.parent, .))) - tgt.clr.sd.thres);
## here,  
## tgt.clr.sd.thres = floor(max deviation / sd)
## = floor((max(MI(candidate.parent, .)) - mean(MI(candidate.parent, .))) / sd(MI(candidate.parent, .))).
## Similarly,
## z.tgt = max(0, (MI(candidate.parent, tgt) - mean(MI(., tgt)) / sd(MI(., tgt))) - candidate.parent.clr.sd.thres);
## here,
## candidate.parent.clr.sd.thres = floor(max deviation / sd)
## = floor((max(MI(., tgt)) - mean(MI(., tgt))) / sd(MI(., tgt))).
## 
## Note that the tgt belongs to time point t_p and the candidate parent belong to time point t_(p + 1).
## Else if (not both z.parent and z.tgt are greater than zero) then
## CLR4.2 edge weight = zero.
##
LearnClr4dot2NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                                max.fanin, mi.net.adj.matrix.list)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]
    
    candidate.parent.node.names <- c()
    candidate.tgt.node.names <- c()
    for (curr.node.name in node.names) {
      parent.full.name <- paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- c(candidate.parent.node.names, parent.full.name)
      rm(parent.full.name)
      
      tgt.full.name <- paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- c(candidate.tgt.node.names, tgt.full.name)
      rm(tgt.full.name)
    }
    rm(curr.node.name)
    
    ## Initialize mutual information matrix with zeroes
    mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(candidate.parent.node.names), 
                                           list(candidate.tgt.node.names)))
    
    #################################################################################
    # ## Build mutual information matrix
    # for (col.idx in 1:(num.nodes - 1)) {
    #   for (col.idx.2 in (col.idx + 1):num.nodes) {
    #     
    #     ## compute_cmi.R
    #     ## (dim1 == 1) => time.pt.idx
    #     ## (dim1 == 2) => (time.pt.idx + 1)
    #     mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, col.idx, ], 
    #                            input.data.discr.3D.curr.ival[2, col.idx.2, ])
    #     
    #     mut.info.matrix[col.idx, col.idx.2] <- mut.info
    #     mut.info.matrix[col.idx.2, col.idx] <- mut.info
    #   }
    #   rm(col.idx.2)
    # }
    # rm(col.idx)
    #################################################################################
    ## Build mutual information matrix
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[2, tgt.idx, ])
        
        ## Mutual info matrix is asymmetric in this case
        mut.info.matrix[parent.idx, tgt.idx] <- mut.info
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    #################################################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                   dimnames = c(list(candidate.parent.node.names), 
                                                list(candidate.tgt.node.names)))
    
    candidate.parent.mean.sd <- matrix(0, nrow = num.nodes, ncol = 3)
    rownames(candidate.parent.mean.sd) <- candidate.parent.node.names
    colnames(candidate.parent.mean.sd) <- c('clr.mean', 'clr.sd', 'clr.sd.thres')
    
    ## Calculate sample mean and sample standard deviation of the given nodes.
    ## It is calculated acc. to the logic in function 'clr()' in
    ## R package 'minet' (version 3.36.0).
    ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
    ## the calculation. Here, the same logic is re-implemented in R.
    ## Since the in-built 'mean()' and 'sd()' functions in 
    ## R version 3.3.2 follows the exact same logic, therefore, the
    ## re-implementation is straight-forward.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.mean.sd[parent.name, 'clr.mean'] <- mean(mut.info.matrix[parent.name, ])
      
      ## var <- 0
      ## for (each sample) {
      ##  sd <- (mean - sample val)
      ##  var <- var + (sd^2)
      ## }
      ## var <- var / (n -1) ## where n = number of samples
      ## sd <- sqrt(var)
      candidate.parent.mean.sd[parent.name, 'clr.sd'] <- sd(mut.info.matrix[parent.name, ])
      
      ## Calc candidate.parent.mean.sd[parent.name, 'clr.sd.thres'] 
      ########################################################################
      if (candidate.parent.mean.sd[parent.name, 'clr.sd'] != 0) {
        tmp <- max((mut.info.matrix[parent.name, ]))
        tmp <- (tmp - candidate.parent.mean.sd[parent.name, 'clr.mean'])
        tmp <- (tmp / candidate.parent.mean.sd[parent.name, 'clr.sd'])
        candidate.parent.mean.sd[parent.name, 'clr.sd.thres'] <- floor(tmp)
        rm(tmp)
        
      } else {
        candidate.parent.mean.sd[parent.name, 'clr.sd.thres'] <- 0
      }
      ########################################################################
      
      
    }
    rm(parent.name)
    
    for (tgt.node.name in candidate.tgt.node.names) {
      
      tgt.clr.mean <- mean(mut.info.matrix[, tgt.node.name])
      tgt.clr.sd <- sd(mut.info.matrix[, tgt.node.name])
      tgt.clr.max <- max(mut.info.matrix[, tgt.node.name])
      
      tgt.clr.sd.thres <- floor((tgt.clr.max - tgt.clr.mean) / tgt.clr.sd)
      
      ## Edge weights of the CLR net are calculated acc. to
      ## 'minet::clr()'
      for (candidate.parent.name in candidate.parent.node.names) {
        
        tmp <- 0
        if (tgt.clr.sd != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - tgt.clr.mean) / tgt.clr.sd
          tmp <- (tmp - tgt.clr.sd.thres)
        }
        z.tgt <- max(0, tmp)
        
        tmp <- 0
        if (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - 
                    candidate.parent.mean.sd[candidate.parent.name, 'clr.mean']) / 
            candidate.parent.mean.sd[candidate.parent.name, 'clr.sd']
          tmp <- (tmp - candidate.parent.mean.sd[candidate.parent.name, 'clr.sd.thres'])
        }
        z.parent <- max(0, tmp)
        
        rm(tmp)
        
        clr.edge.wt <- 0
        if ((z.tgt > 0) && (z.parent > 0)) {
          clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
          clr.edge.wt <- sqrt(clr.edge.wt)
        }
        
        rm(z.tgt, z.parent)
        
        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      rm(candidate.parent.name)
      
    }
    rm(tgt.node.name)
    
    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                dimnames = c(list(candidate.parent.node.names), 
                                             list(candidate.tgt.node.names)))
    
    ## For each target node
    for (col.idx in 1:num.nodes) {
      
      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]
      
      ## Count number of neighbours having positive edge weight
      num.nbrs <- length(edge.wts[edge.wts > 0])
      
      if (num.nbrs >= max.fanin) {
        
        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
        
        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
        
      } else if (num.nbrs < max.fanin) {
        
        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    rm(col.idx)
    
    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    
    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  rm(time.pt.idx)
  
  return(mi.net.adj.matrix.list)
}

############################################################################################

############################################################################################
## Goal: Learn CLR5 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
## If (z.parent and z.tgt both are greater than zero) then
## CLR5 edge weight = z.tgt where
## z.parent = max(0, (MI(candidate.parent, tgt) - mean(MI(candidate.parent, .)) / sd(MI(candidate.parent, .)))) and
## z.tgt = max(0, (MI(candidate.parent, tgt) - mean(MI(., tgt)) / sd(MI(., tgt)))).
## If the tgt belongs to time point t_p, then the candidate parent must belong to time point t_(p + 1).
## Else if (not both z.parent and z.tgt are greater than zero) then
## CLR5 edge weight = zero.
##
LearnClr5NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                            max.fanin, mi.net.adj.matrix.list)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]
    
    candidate.parent.node.names <- c()
    candidate.tgt.node.names <- c()
    for (curr.node.name in node.names) {
      parent.full.name <- paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- c(candidate.parent.node.names, parent.full.name)
      rm(parent.full.name)
      
      tgt.full.name <- paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- c(candidate.tgt.node.names, tgt.full.name)
      rm(tgt.full.name)
    }
    rm(curr.node.name)
    
    ## Initialize mutual information matrix with zeroes
    mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(candidate.parent.node.names), 
                                           list(candidate.tgt.node.names)))
    
    ## Build mutual information matrix
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[2, tgt.idx, ])
        
        ## Mutual info matrix is asymmetric in this case
        mut.info.matrix[parent.idx, tgt.idx] <- mut.info
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    #################################################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                   dimnames = c(list(candidate.parent.node.names), 
                                                list(candidate.tgt.node.names)))
    
    candidate.parent.mean.sd <- matrix(0, nrow = num.nodes, ncol = 2)
    rownames(candidate.parent.mean.sd) <- candidate.parent.node.names
    colnames(candidate.parent.mean.sd) <- c('clr.mean', 'clr.sd')
    
    ## Calculate sample mean and sample standard deviation of the given nodes.
    ## It is calculated acc. to the logic in function 'clr()' in
    ## R package 'minet' (version 3.36.0).
    ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
    ## the calculation. Here, the same logic is re-implemented in R.
    ## Since the in-built 'mean()' and 'sd()' functions in 
    ## R version 3.3.2 follows the exact same logic, therefore, the
    ## re-implementation is straight-forward.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.mean.sd[parent.name, 'clr.mean'] <- mean(mut.info.matrix[parent.name, ])
      
      ## var <- 0
      ## for (each sample) {
      ##  sd <- (mean - sample val)
      ##  var <- var + (sd^2)
      ## }
      ## var <- var / (n -1) ## where n = number of samples
      ## sd <- sqrt(var)
      candidate.parent.mean.sd[parent.name, 'clr.sd'] <- sd(mut.info.matrix[parent.name, ])
    }
    rm(parent.name)
    
    for (tgt.node.name in candidate.tgt.node.names) {
      
      tgt.clr.mean <- mean(mut.info.matrix[, tgt.node.name])
      tgt.clr.sd <- sd(mut.info.matrix[, tgt.node.name])
      
      ## Edge weights of the CLR net are calculated acc. to
      ## 'minet::clr()'
      for (candidate.parent.name in candidate.parent.node.names) {
        
        tmp <- 0
        if (tgt.clr.sd != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - tgt.clr.mean) / tgt.clr.sd
        }
        z.tgt <- max(0, tmp)
        
        tmp <- 0
        if (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - 
                    candidate.parent.mean.sd[candidate.parent.name, 'clr.mean']) / 
            candidate.parent.mean.sd[candidate.parent.name, 'clr.sd']
        }
        z.parent <- max(0, tmp)
        
        rm(tmp)
        
        clr.edge.wt <- 0
        if ((z.tgt > 0) && (z.parent > 0)) {
          # clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
          # clr.edge.wt <- sqrt(clr.edge.wt)
          clr.edge.wt <- z.tgt
        }
        
        rm(z.tgt, z.parent)
        
        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      rm(candidate.parent.name)
      
    }
    rm(tgt.node.name)
    
    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                dimnames = c(list(candidate.parent.node.names), 
                                             list(candidate.tgt.node.names)))
    
    ## For each target node
    for (col.idx in 1:num.nodes) {
      
      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]
      
      ## Count number of neighbours having positive edge weight
      num.nbrs <- length(edge.wts[edge.wts > 0])
      
      if (num.nbrs >= max.fanin) {
        
        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
        
        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
        
      } else if (num.nbrs < max.fanin) {
        
        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    rm(col.idx)
    
    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    
    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  rm(time.pt.idx)
  
  return(mi.net.adj.matrix.list)
}

############################################################################################
############################################################################################
## Goal: Learn CLR6 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
## If (z.parent and z.tgt both are greater than zero) then
## CLR6 edge weight = z.parent where
## z.parent = max(0, (MI(candidate.parent, tgt) - mean(MI(candidate.parent, .)) / sd(MI(candidate.parent, .)))) and
## z.tgt = max(0, (MI(candidate.parent, tgt) - mean(MI(., tgt)) / sd(MI(., tgt)))).
## If the tgt belongs to time point t_p, then the candidate parent must belong to time point t_(p + 1).
## Else if (not both z.parent and z.tgt are greater than zero) then
## CLR6 edge weight = zero.
##
LearnClr6NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                            max.fanin, mi.net.adj.matrix.list)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]
    
    candidate.parent.node.names <- c()
    candidate.tgt.node.names <- c()
    for (curr.node.name in node.names) {
      parent.full.name <- paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- c(candidate.parent.node.names, parent.full.name)
      rm(parent.full.name)
      
      tgt.full.name <- paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- c(candidate.tgt.node.names, tgt.full.name)
      rm(tgt.full.name)
    }
    rm(curr.node.name)
    
    ## Initialize mutual information matrix with zeroes
    mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(candidate.parent.node.names), 
                                           list(candidate.tgt.node.names)))
    
    ## Build mutual information matrix
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[2, tgt.idx, ])
        
        ## Mutual info matrix is asymmetric in this case
        mut.info.matrix[parent.idx, tgt.idx] <- mut.info
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    #################################################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                   dimnames = c(list(candidate.parent.node.names), 
                                                list(candidate.tgt.node.names)))
    
    candidate.parent.mean.sd <- matrix(0, nrow = num.nodes, ncol = 2)
    rownames(candidate.parent.mean.sd) <- candidate.parent.node.names
    colnames(candidate.parent.mean.sd) <- c('clr.mean', 'clr.sd')
    
    ## Calculate sample mean and sample standard deviation of the given nodes.
    ## It is calculated acc. to the logic in function 'clr()' in
    ## R package 'minet' (version 3.36.0).
    ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
    ## the calculation. Here, the same logic is re-implemented in R.
    ## Since the in-built 'mean()' and 'sd()' functions in 
    ## R version 3.3.2 follows the exact same logic, therefore, the
    ## re-implementation is straight-forward.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.mean.sd[parent.name, 'clr.mean'] <- mean(mut.info.matrix[parent.name, ])
      
      ## var <- 0
      ## for (each sample) {
      ##  sd <- (mean - sample val)
      ##  var <- var + (sd^2)
      ## }
      ## var <- var / (n -1) ## where n = number of samples
      ## sd <- sqrt(var)
      candidate.parent.mean.sd[parent.name, 'clr.sd'] <- sd(mut.info.matrix[parent.name, ])
    }
    rm(parent.name)
    
    for (tgt.node.name in candidate.tgt.node.names) {
      
      tgt.clr.mean <- mean(mut.info.matrix[, tgt.node.name])
      tgt.clr.sd <- sd(mut.info.matrix[, tgt.node.name])
      
      ## Edge weights of the CLR net are calculated acc. to
      ## 'minet::clr()'
      for (candidate.parent.name in candidate.parent.node.names) {
        
        tmp <- 0
        if (tgt.clr.sd != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - tgt.clr.mean) / tgt.clr.sd
        }
        z.tgt <- max(0, tmp)
        
        tmp <- 0
        if (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - 
                    candidate.parent.mean.sd[candidate.parent.name, 'clr.mean']) / 
            candidate.parent.mean.sd[candidate.parent.name, 'clr.sd']
        }
        z.parent <- max(0, tmp)
        
        rm(tmp)
        
        clr.edge.wt <- 0
        if ((z.tgt > 0) && (z.parent > 0)) {
          # clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
          # clr.edge.wt <- sqrt(clr.edge.wt)
          clr.edge.wt <- z.parent
        }
        
        rm(z.tgt, z.parent)
        
        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      rm(candidate.parent.name)
      
    }
    rm(tgt.node.name)
    
    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                dimnames = c(list(candidate.parent.node.names), 
                                             list(candidate.tgt.node.names)))
    
    ## For each target node
    for (col.idx in 1:num.nodes) {
      
      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]
      
      ## Count number of neighbours having positive edge weight
      num.nbrs <- length(edge.wts[edge.wts > 0])
      
      if (num.nbrs >= max.fanin) {
        
        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
        
        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
        
      } else if (num.nbrs < max.fanin) {
        
        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    rm(col.idx)
    
    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    
    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  rm(time.pt.idx)
  
  return(mi.net.adj.matrix.list)
}

############################################################################################
############################################################################################
## Goal: Learn CLR7 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
## If (z.parent and z.tgt both are greater than zero) then
## CLR7 edge weight = (z.parent^2 + z.tgt^2)^0.5 where
## z.parent = max(0, (MI(candidate.parent, tgt) - median(MI(candidate.parent, .))) and
## z.tgt = max(0, (MI(candidate.parent, tgt) - median(MI(., tgt)))).
## If the tgt belongs to time point t_p, then the candidate parent must belong to time point t_(p + 1).
## Else if (not both z.parent and z.tgt are greater than zero) then
## CLR7 edge weight = zero.
##
LearnClr7NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                            max.fanin, mi.net.adj.matrix.list)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]
    
    candidate.parent.node.names <- c()
    candidate.tgt.node.names <- c()
    for (curr.node.name in node.names) {
      parent.full.name <- paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- c(candidate.parent.node.names, parent.full.name)
      rm(parent.full.name)
      
      tgt.full.name <- paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- c(candidate.tgt.node.names, tgt.full.name)
      rm(tgt.full.name)
    }
    rm(curr.node.name)
    
    ## Initialize mutual information matrix with zeroes
    mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(candidate.parent.node.names), 
                                           list(candidate.tgt.node.names)))
    
    ## Build mutual information matrix
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[2, tgt.idx, ])
        
        ## Mutual info matrix is asymmetric in this case
        mut.info.matrix[parent.idx, tgt.idx] <- mut.info
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    #################################################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                   dimnames = c(list(candidate.parent.node.names), 
                                                list(candidate.tgt.node.names)))
    
    candidate.parent.median <- matrix(0, nrow = num.nodes, ncol = 1)
    rownames(candidate.parent.median) <- candidate.parent.node.names
    colnames(candidate.parent.median) <- c('clr.median')
    
    ## Calculate sample median of the given nodes.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.median[parent.name, 'clr.median'] <- median(mut.info.matrix[parent.name, ])
    }
    rm(parent.name)
    
    for (tgt.node.name in candidate.tgt.node.names) {
      
      tgt.clr.median <- median(mut.info.matrix[, tgt.node.name])
      
      for (candidate.parent.name in candidate.parent.node.names) {
        
        tmp <- 0
        tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - tgt.clr.median)
        z.tgt <- max(0, tmp)
        
        tmp <- 0
        tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - 
                  candidate.parent.median[candidate.parent.name, 'clr.median'])
        z.parent <- max(0, tmp)
        
        rm(tmp)
        
        clr.edge.wt <- 0
        if ((z.tgt > 0) && (z.parent > 0)) {
          clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
          clr.edge.wt <- sqrt(clr.edge.wt)
        }
        
        rm(z.tgt, z.parent)
        
        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      rm(candidate.parent.name)
      
    }
    rm(tgt.node.name)
    rm(candidate.parent.median)
    
    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                dimnames = c(list(candidate.parent.node.names), 
                                             list(candidate.tgt.node.names)))
    
    ## For each target node
    for (col.idx in 1:num.nodes) {
      
      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]
      
      ## Count number of neighbours having positive edge weight
      num.nbrs <- length(edge.wts[edge.wts > 0])
      
      if (num.nbrs >= max.fanin) {
        
        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
        
        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
        
      } else if (num.nbrs < max.fanin) {
        
        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    rm(col.idx)
    
    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    
    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  rm(time.pt.idx)
  
  return(mi.net.adj.matrix.list)
}

############################################################################################
############################################################################################
## Goal: Learn CLR8 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
## If (z.parent and z.tgt both are greater than zero) then
## CLR8 edge weight = z.parent where
## z.parent = max(0, (MI(candidate.parent, tgt) - mean(MI(candidate.parent, .)) / sd(MI(candidate.parent, .)))) and
## z.tgt = max(0, (MI(candidate.parent, tgt) - mean(MI(., tgt)) / sd(MI(., tgt)))).
## If the tgt belongs to time point t_p, then the candidate parent must belong to time point t_(p + 1).
## Else if (not both z.parent and z.tgt are greater than zero) then
## CLR8 edge weight = zero.
##
LearnClr8NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                            max.fanin, mi.net.adj.matrix.list, p.thres)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]
    
    candidate.parent.node.names <- c()
    candidate.tgt.node.names <- c()
    for (curr.node.name in node.names) {
      parent.full.name <- paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- c(candidate.parent.node.names, parent.full.name)
      rm(parent.full.name)
      
      tgt.full.name <- paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- c(candidate.tgt.node.names, tgt.full.name)
      rm(tgt.full.name)
    }
    rm(curr.node.name)
    
    ## Initialize mutual information matrix with zeroes
    mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(candidate.parent.node.names), 
                                           list(candidate.tgt.node.names)))
    
    ## Build mutual information matrix
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[2, tgt.idx, ])
        
        ## Mutual info matrix is asymmetric in this case
        mut.info.matrix[parent.idx, tgt.idx] <- mut.info
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    #################################################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                   dimnames = c(list(candidate.parent.node.names), 
                                                list(candidate.tgt.node.names)))
    
    candidate.parent.mean.sd <- matrix(0, nrow = num.nodes, ncol = 2)
    rownames(candidate.parent.mean.sd) <- candidate.parent.node.names
    colnames(candidate.parent.mean.sd) <- c('clr.mean', 'clr.sd')
    
    ## Calculate sample mean and sample standard deviation of the given nodes.
    ## It is calculated acc. to the logic in function 'clr()' in
    ## R package 'minet' (version 3.36.0).
    ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
    ## the calculation. Here, the same logic is re-implemented in R.
    ## Since the in-built 'mean()' and 'sd()' functions in 
    ## R version 3.3.2 follows the exact same logic, therefore, the
    ## re-implementation is straight-forward.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.mean.sd[parent.name, 'clr.mean'] <- mean(mut.info.matrix[parent.name, ])
      
      ## var <- 0
      ## for (each sample) {
      ##  sd <- (mean - sample val)
      ##  var <- var + (sd^2)
      ## }
      ## var <- var / (n -1) ## where n = number of samples
      ## sd <- sqrt(var)
      candidate.parent.mean.sd[parent.name, 'clr.sd'] <- sd(mut.info.matrix[parent.name, ])
    }
    rm(parent.name)
    
    for (tgt.node.name in candidate.tgt.node.names) {
      
      tgt.clr.mean <- mean(mut.info.matrix[, tgt.node.name])
      tgt.clr.sd <- sd(mut.info.matrix[, tgt.node.name])
      
      ## Edge weights of the CLR net are calculated acc. to
      ## 'minet::clr()'
      for (candidate.parent.name in candidate.parent.node.names) {
        
        num.samples.per.timept <- dim(input.data.discr.3D.curr.ival)[3]
        
        ## Calculate p-value.
        ## Ref::
        ## Code: https://www.cyclismo.org/tutorial/R/pValues.html
        ## Theory: https://www.dummies.com/education/math/statistics/how-to-determine-a-p-value-when-testing-a-null-hypothesis/
        p.val <- 1
        z.score <- 0
        if (tgt.clr.sd != 0) {
          z.score <- ((tgt.clr.mean - mut.info.matrix[candidate.parent.name, tgt.node.name]) / 
            (tgt.clr.sd / sqrt(num.samples.per.timept)))
          p.val <- 2 * pnorm(-abs(z.score))
          p.val <- (1 - p.val)
        }
        rm(z.score)
        
        z.tgt <- 0
        if (p.val <= p.thres) {
          z.tgt <- p.val
        }
        
        p.val <- 1
        z.score <- 0
        if (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] != 0) {
          z.score <- ((candidate.parent.mean.sd[candidate.parent.name, 'clr.mean'] - 
                         mut.info.matrix[candidate.parent.name, tgt.node.name]) / 
            (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] / sqrt(num.samples.per.timept)))
          p.val <- 2 * pnorm(-abs(z.score))
          p.val <- (1 - p.val)
        }
        rm(z.score)
        
        z.parent <- 0
        if (p.val <= p.thres) {
          z.parent <- p.val
        }
        rm(p.val)
        
        clr.edge.wt <- 0
        if ((z.tgt > 0) || (z.parent > 0)) {
          clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
          clr.edge.wt <- sqrt(clr.edge.wt)
        }
        
        rm(z.tgt, z.parent)
        
        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      rm(candidate.parent.name)
      rm(num.samples.per.timept)
      
    }
    rm(tgt.node.name)
    
    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                dimnames = c(list(candidate.parent.node.names), 
                                             list(candidate.tgt.node.names)))
    
    ## For each target node
    for (col.idx in 1:num.nodes) {
      
      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]
      
      ## Count number of neighbours having positive edge weight
      num.nbrs <- length(edge.wts[edge.wts > 0])
      
      if (num.nbrs >= max.fanin) {
        
        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
        
        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
        
      } else if (num.nbrs < max.fanin) {
        
        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    rm(col.idx)
    
    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    
    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  rm(time.pt.idx)
  
  return(mi.net.adj.matrix.list)
}

############################################################################################

############################################################################################
## Goal: Learn CLR9 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
## CLR9 edge weight = (z.parent^2 + z.tgt^2)^0.5 where
## z.parent = max(0, (MI(candidate.parent, tgt) - mean(MI(candidate.parent, .)) / sd(MI(candidate.parent, .)))) and
## z.tgt = max(0, (MI(candidate.parent, tgt) - mean(MI(., tgt)) / sd(MI(., tgt)))).
## If the tgt belongs to time point t_p, then the candidate parent must belong to time point t_(p + 1).
## The MI(*, *) values come from a 'refined' mutual info matrix, where the rows represent 
## the candidate parents and the cols represent the tgts. The refined matrix is derived
## from the raw mutual information matrix between candidate parents and tgts, where 
## the (i, j)^th elt represents the mutual information value between the i^th candidate
## parent and the j^th tgt. The derivation from the raw to refined matrix is done as follows:
## for the i^th candidate parent and the j^th tgt,
## if there exists another candidate parent k (!= i) such that
## MI(i, j) in the raw matrix is less than MI(k, j) in the raw matrix as well as
## less than MI(i, k), which is the mutual info value between the i^th and k^th parent,
## then MI(i, j) is reduced to zero in the refined matrix.
## Otherwise, MI(i, j) is directly copied from the raw matrix to the refined one.
##
LearnClr9NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                            max.fanin, mi.net.adj.matrix.list)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]
    
    candidate.parent.node.names <- c()
    candidate.tgt.node.names <- c()
    for (curr.node.name in node.names) {
      parent.full.name <- paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- c(candidate.parent.node.names, parent.full.name)
      rm(parent.full.name)
      
      tgt.full.name <- paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- c(candidate.tgt.node.names, tgt.full.name)
      rm(tgt.full.name)
    }
    rm(curr.node.name)
    
    ## Initialize raw mutual information matrix between parents and targets with zeroes
    mut.info.matrix.parent.tgt.raw <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                             dimnames = c(list(candidate.parent.node.names), 
                                                          list(candidate.tgt.node.names)))
    
    #################################################################################
    # ## Build mutual information matrix
    # for (col.idx in 1:(num.nodes - 1)) {
    #   for (col.idx.2 in (col.idx + 1):num.nodes) {
    #     
    #     ## compute_cmi.R
    #     ## (dim1 == 1) => time.pt.idx
    #     ## (dim1 == 2) => (time.pt.idx + 1)
    #     mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, col.idx, ], 
    #                            input.data.discr.3D.curr.ival[2, col.idx.2, ])
    #     
    #     mut.info.matrix.parent.tgt[col.idx, col.idx.2] <- mut.info
    #     mut.info.matrix.parent.tgt[col.idx.2, col.idx] <- mut.info
    #   }
    #   rm(col.idx.2)
    # }
    # rm(col.idx)
    #################################################################################
    ## Build mutual information matrix between parents and targets
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[2, tgt.idx, ])
        
        ## Mutual info matrix is asymmetric in this case
        mut.info.matrix.parent.tgt.raw[parent.idx, tgt.idx] <- mut.info
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    #################################################################################
    
    #################################################################################
    ## Build mutual information matrix between parents and parents
    
    ## Initialize mutual information matrix between parents and parents with zeroes
    mut.info.matrix.parent.parent <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                            dimnames = c(list(candidate.parent.node.names), 
                                                         list(candidate.parent.node.names)))
    
    for (parent.idx in 1:(num.nodes - 1)) {
      for (parent.idx.2 in (parent.idx + 1):num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[1, parent.idx.2, ])
        
        ## Mutual info matrix is symmetric in this case
        mut.info.matrix.parent.parent[parent.idx, parent.idx.2] <- mut.info
        mut.info.matrix.parent.parent[parent.idx.2, parent.idx] <- mut.info
      }
      rm(parent.idx.2)
    }
    rm(parent.idx)
    #################################################################################    
    
    #################################################################################
    ## Refine the raw parent-tgt mutual information matrix with the help of 
    ## parent-parent mutual info matrix
    
    ## Initialize the refined matrix with the raw one
    mut.info.matrix.parent.tgt.refined <- mut.info.matrix.parent.tgt.raw
    
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        parent.tgt.mut.info.raw <- mut.info.matrix.parent.tgt.raw[parent.idx, tgt.idx]
        for (parent.idx.2 in num.nodes) {
          if (parent.idx.2 != parent.idx) {
            parent.tgt.mut.info.raw.2 <- mut.info.matrix.parent.tgt.raw[parent.idx.2, tgt.idx]
            
            parent.parent.mut.info <- mut.info.matrix.parent.parent[parent.idx, parent.idx.2]
            
            if (parent.tgt.mut.info.raw < min(parent.tgt.mut.info.raw.2, parent.parent.mut.info)) {
              mut.info.matrix.parent.tgt.refined[parent.idx, tgt.idx] <- 0
            }
          }
        }
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    rm(mut.info.matrix.parent.parent, mut.info.matrix.parent.tgt.raw)
    
    #################################################################################
    
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                   dimnames = c(list(candidate.parent.node.names), 
                                                list(candidate.tgt.node.names)))
    
    candidate.parent.mean.sd <- matrix(0, nrow = num.nodes, ncol = 2)
    rownames(candidate.parent.mean.sd) <- candidate.parent.node.names
    colnames(candidate.parent.mean.sd) <- c('clr.mean', 'clr.sd')
    
    ## Calculate sample mean and sample standard deviation of the given nodes.
    ## It is calculated acc. to the logic in function 'clr()' in
    ## R package 'minet' (version 3.36.0).
    ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
    ## the calculation. Here, the same logic is re-implemented in R.
    ## Since the in-built 'mean()' and 'sd()' functions in 
    ## R version 3.3.2 follows the exact same logic, therefore, the
    ## re-implementation is straight-forward.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.mean.sd[parent.name, 'clr.mean'] <- mean(mut.info.matrix.parent.tgt.refined[parent.name, ])
      
      ## var <- 0
      ## for (each sample) {
      ##  sd <- (mean - sample val)
      ##  var <- var + (sd^2)
      ## }
      ## var <- var / (n -1) ## where n = number of samples
      ## sd <- sqrt(var)
      candidate.parent.mean.sd[parent.name, 'clr.sd'] <- sd(mut.info.matrix.parent.tgt.refined[parent.name, ])
    }
    rm(parent.name)
    
    for (tgt.node.name in candidate.tgt.node.names) {
      
      tgt.clr.mean <- mean(mut.info.matrix.parent.tgt.refined[, tgt.node.name])
      tgt.clr.sd <- sd(mut.info.matrix.parent.tgt.refined[, tgt.node.name])
      
      ## Edge weights of the CLR net are calculated acc. to
      ## 'minet::clr()'
      for (candidate.parent.name in candidate.parent.node.names) {
        
        tmp <- 0
        if (tgt.clr.sd != 0) {
          tmp <- (mut.info.matrix.parent.tgt.refined[candidate.parent.name, tgt.node.name] - tgt.clr.mean) / tgt.clr.sd
        }
        z.tgt <- max(0, tmp)
        
        tmp <- 0
        if (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] != 0) {
          tmp <- (mut.info.matrix.parent.tgt.refined[candidate.parent.name, tgt.node.name] - 
                    candidate.parent.mean.sd[candidate.parent.name, 'clr.mean']) / 
            candidate.parent.mean.sd[candidate.parent.name, 'clr.sd']
        }
        z.parent <- max(0, tmp)
        
        rm(tmp)
        
        clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
        clr.edge.wt <- sqrt(clr.edge.wt)
        rm(z.tgt, z.parent)
        
        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      rm(candidate.parent.name)
      
    }
    rm(tgt.node.name)
    
    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                dimnames = c(list(candidate.parent.node.names), 
                                             list(candidate.tgt.node.names)))
    
    ## For each target node
    for (col.idx in 1:num.nodes) {
      
      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]
      
      ## Count number of neighbours having positive edge weight
      num.nbrs <- length(edge.wts[edge.wts > 0])
      
      if (num.nbrs >= max.fanin) {
        
        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
        
        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
        
      } else if (num.nbrs < max.fanin) {
        
        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    rm(col.idx)
    
    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    
    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  rm(time.pt.idx)
  
  return(mi.net.adj.matrix.list)
}

############################################################################################

############################################################################################
## Goal: Learn CLR10 net. 
## TODO(sap)
## For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
##
LearnClr10NetMfi <- function(mut.info.matrix, mi.net.adj.matrix, 
                             num.nodes, max.fanin, grade.range, 
                             output.dirname) {
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  library(minet)
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Estimate CLR net adjacency matrix.
  ## It is a weighted undirected matrix.
  mi.net.adj.matrix.wt <- minet::clr(mut.info.matrix)
  
  # Replace 'NaN's with zeroes. 
  ## 'NaN' is produced when a corr. variable has zero variance.
  mi.net.adj.matrix.wt[is.nan(mi.net.adj.matrix.wt)] <- 0
  
  # writeLines('\n mi.net.adj.matrix.wt = \n')
  # print(mi.net.adj.matrix.wt)
  save(mi.net.adj.matrix.wt, file = paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))
  
  #####################################################################
  ## Begin: Apply grade threshold
  #####################################################################
  
  ## Backup
  mi.net.adj.matrix.wt.bak <- mi.net.adj.matrix.wt
  
  for (tgt.node.idx in 1:num.nodes) {
    tgt.mean <- mean(mi.net.adj.matrix.wt.bak[-(tgt.node.idx), tgt.node.idx])
    tgt.sd <- sd(mi.net.adj.matrix.wt.bak[-(tgt.node.idx), tgt.node.idx])
    tgt.max <- max(mi.net.adj.matrix.wt.bak[-(tgt.node.idx), tgt.node.idx])
    
    ## Max s.d. away = how many s.d. away is max from mean
    tgt.max.sda <- floor((tgt.max - tgt.mean) / tgt.sd)
    
    ## All nodes are candidate parents of the tgt node
    ## except the tgt node itself
    for (cand.parent.idx in (1:num.nodes)[-(tgt.node.idx)]) {
      
      cand.clr.wt <- mi.net.adj.matrix.wt.bak[cand.parent.idx, tgt.node.idx]
      
      if (cand.clr.wt == 0) {
        next
      }
      
      cand.sda <- floor((cand.clr.wt - tgt.mean) / tgt.sd)
      
      ## Candidate parent gets the top grade if
      ## it is same s.d. away from mean as that of max.
      ## In other words, candidate parent gets grade 1 if
      ## its s.d.a. is same as max's s.d.a.
      ## In general, candidate parent gets grade '(x+1)' if
      ## its s.d.a is 'x' less than that of max.
      cand.grade <- ((tgt.max.sda - cand.sda) + 1)
      
      ## If the candidate parent's grade is not within
      ## the user defined grade range, then
      ## reject the candidate.
      if (cand.grade > grade.range) {
        mi.net.adj.matrix.wt[cand.parent.idx, tgt.node.idx] <- 0
      }
    }
    rm(cand.parent.idx)
  }
  rm(tgt.node.idx)
  
  ## Remove the backup matrix
  rm(mi.net.adj.matrix.wt.bak)
  #####################################################################
  ## End: Apply grade threshold
  #####################################################################
  
  
  #####################################################################
  ## Begin: Apply max fan-in
  #####################################################################
  # For each target node
  for (col.idx in 1:num.nodes) {
    # Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]
    
    # Count number of neighbours having positive edge weight
    num.nbrs <- length(edge.wts[edge.wts > 0])
    
    if (num.nbrs >= max.fanin) {
      # Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      # Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
      
      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
      
      ## The following line is not required since 'mi.net.adj.matrix' is initialized
      ## with all zeroes
      # mi.net.adj.matrix[-(valid.nbrs), col.idx] <- 0
    } else if (num.nbrs < max.fanin) {
      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  #####################################################################
  ## End: Apply max fan-in
  #####################################################################
  
  return(mi.net.adj.matrix)
}
############################################################################################