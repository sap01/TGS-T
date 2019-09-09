## This R file contains multiple graph related functions.

####################################################################################################3
## Goal: Given a directed network adjacency matrix and a node name, returns names of all the nodes
## reachable from the given node. In the given directed network adjacency matrix,
## rows = source nodes, cols = tgt nodes. (i, j)th cell = 1 implies a directed edge from i
## to j; = 0 implies no directed edge from i to j. Row names of the given matrix must be
## the names of the nodes.
##
reachable.nodes <- function(di.net.adj.matrix, src.node.name)
{
  reachable.node.names <- c()
  
  ## If the given node has at least one target node
  if (length(which(di.net.adj.matrix[src.node.name, ] == 1)) > 0)
  {
    reachable.node.names <- names(which(di.net.adj.matrix[src.node.name, ] == 1))
    
    ## Avoid self loop
    to.traverse.node.names <- setdiff(reachable.node.names, src.node.name)
    
    while (length(to.traverse.node.names) != 0)
    {
      curr.src.node.name <- to.traverse.node.names[1]
      
      if (length(which(di.net.adj.matrix[curr.src.node.name, ] == 1)) > 0)
      {
        new.reachable.node.names <- names(which(di.net.adj.matrix[curr.src.node.name, ] == 1))
        
        ## Avoid directed cycle
        new.reachable.node.names <- setdiff(new.reachable.node.names, reachable.node.names)
        
        reachable.node.names <- union(reachable.node.names, new.reachable.node.names)
        
        to.traverse.node.names <- union(to.traverse.node.names, new.reachable.node.names)
      }
      
      to.traverse.node.names <- setdiff(to.traverse.node.names, curr.src.node.name) 
    }
  }
  
  return(reachable.node.names)
}