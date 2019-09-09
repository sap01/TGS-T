## Goal: Count the number of feed-forward edges in a given undirected network.

CountFeedFwdEdgesUndi <- function(undi.net.adj.matrix) {
  
  num.nodes <- nrow(undi.net.adj.matrix)
  
  count <- 0
  
  for (node.idx.1 in 3:num.nodes) {
    for (node.idx.2 in 2:(node.idx.1 - 1)) {
      for (node.idx.3 in 1:(node.idx.2 - 1)) {
        
        if ((undi.net.adj.matrix[node.idx.1, node.idx.2] == 1) & 
          (undi.net.adj.matrix[node.idx.2, node.idx.3] == 1) & 
          (undi.net.adj.matrix[node.idx.1, node.idx.3] == 1)) {
          
          count <- (count + 1)
        }
      }
      rm(node.idx.3)
    }
    rm(node.idx.2)
  }
  rm(node.idx.1)
  
  print(paste('No. of feed-forward edge is: ', count, sep = ''))
  
}