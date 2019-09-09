############################################################################################
## Goal: Learn Correlation net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
##
LearnCorrNetMfi <- function(input.data.discr, num.nodes, node.names, num.timepts, 
                            max.fanin, output.dirname, corr.net.adj.matrix)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Initialize weighted adjacency matrix of the correlation network
  corr.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
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
    
    # Initialize correlation matrix with zeroes
    corr.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, dimnames = c(list(node.names), list(node.names)))
    
    ## Build assymetric correlation matrix
    for (row.idx in 1:num.nodes) {
      for (col.idx in 1:num.nodes) {
        
        corr <- cor(input.data.discr.curr.series[1:(num.timepts-1), row.idx], 
                    input.data.discr.curr.series[2:num.timepts, col.idx], 
                    method = 'spearman')
        
        corr.matrix[row.idx, col.idx] <- corr
      }
      rm(col.idx)
    }
    rm(row.idx)
    
    ## Estimate weighted adj matrix of the correlation net
    ## corr. to the current time series
    corr.net.adj.matrix.wt.curr.series <- corr.matrix
    for (row.idx in 1:num.nodes) {
      for (col.idx in 1:num.nodes) {
        corr.net.adj.matrix.wt.curr.series[row.idx, col.idx] <- 
          abs(corr.net.adj.matrix.wt.curr.series[row.idx, col.idx])
      }
      rm(col.idx)
    }
    rm(row.idx)
    
    ## 'corr.net.adj.matrix.wt.curr.series' has non-neg values.
    ## Summing up all time-series-specific 'corr.net.adj.matrix.wt.curr.series' matrices
    ## produces a ('num.nodes' by 'num.nodes') matrix of non-neg
    ## values, namely 'corr.net.adj.matrix.wt'. Later 'corr.net.adj.matrix.wt'
    ## will be divided by 'num.time.series', thus producing the 
    ## arthmetic mean of all time-series-specific 
    ## 'corr.net.adj.matrix.wt.curr.series' matrices.
    corr.net.adj.matrix.wt <- (corr.net.adj.matrix.wt + corr.net.adj.matrix.wt.curr.series)
    
  }
  rm(time.series.idx)
  
  ## Arthmetic mean of all time-series-specific 
  ## 'corr.net.adj.matrix.wt' matrices.
  corr.net.adj.matrix.wt <- (corr.net.adj.matrix.wt / num.time.series)
  
  save(corr.net.adj.matrix.wt, file = paste(output.dirname, 'corr.net.adj.matrix.wt.RData', sep = '/'))
  
  ##############################################################
  ## Begin:
  ## Estimate the unweighted adjacency matrix 'corr.net.adj.matrix'
  ## using the weighted adjacency matrix 'corr.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################
  
  ## For each target node
  for (col.idx in 1:num.nodes) {
    
    ## Weights of the edges with the target node
    edge.wts <- corr.net.adj.matrix.wt[, col.idx]
    
    ## Count number of neighbours having positive edge weight
    num.nbrs <- length(edge.wts[edge.wts > 0])
    
    if (num.nbrs >= max.fanin) {
      
      ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      ## Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
      
      corr.net.adj.matrix[valid.nbrs, col.idx] <- 1
      
    } else if (num.nbrs < max.fanin) {
      
      # Retain all the neighbours
      corr.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  rm(col.idx)
  
  ##############################################################
  ## End:
  ## Estimate the unweighted adjacency matrix 'corr.net.adj.matrix'
  ## using the weighted adjacency matrix 'corr.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################
  
  return(corr.net.adj.matrix)
  
}

############################################################################################