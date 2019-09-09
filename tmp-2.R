#######################################################################################################
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