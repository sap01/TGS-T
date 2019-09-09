## Goal: Given a directed network, convert it into an undirected network

ConvertDinetToUndinet <- function(di.net) {
  
  ## Initilize undirected network with zeroes
  undi.net <- di.net
  undi.net[undi.net] <- 0
  
  for (row.idx in 1:nrow(di.net)) {
    for (col.idx in 1:ncol(di.net)) {
      if (di.net[row.idx, col.idx] == 1) {
        undi.net[row.idx, col.idx] <- 1
        undi.net[col.idx, row.idx] <- 1
      }
    }
    rm(col.idx)
  }
  rm(row.idx)
  
return(undi.net)  
}