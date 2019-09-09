## Goal: Generate synthetic networks and synthetic datasets using 
## R package 'EDISON' version 1.1.1
## 
##################################################################################################

##################################################################################################
## Goal: Generate synthetic time-varying gene regulatory networks
## using R package 'EDISON' version 1.1.1
##
## Input params:
## See the input params of function 'generateNetwork()' 
## in R package 'EDISON' version 1.1.1.
##
## Output params:
## See the output params of function 'generateNetwork()' 
## in R package 'EDISON' version 1.1.1.
## 'k_bar' must be at least 1.
##################################################################################################
GenEdiNet <- function(true.net.filename, true.net.adj.mx.filename, 
                      lambda_2 = 0.45, q = 10, min_phase_length = 1, k_bar = 5, l = 10, 
                      lambda_3 = 2, spacing = 1, gauss_weights = TRUE, same = FALSE, 
                      change_method ='sequential', cps = NULL) {
  
  ##------------------------------------------------------------
  ## Begin: Load the required packages
  ##------------------------------------------------------------
  library(EDISON)
  ##------------------------------------------------------------
  ## End: Load the required packa                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ges
  ##------------------------------------------------------------
  
  ## Generate synthetic networks with the given parameters.
  ##  Param 'fixed' must always be TRUE.
  edi.net <- EDISON::generateNetwork(lambda_2 = lambda_2,
                                     q = q,
                                     min_phase_length = min_phase_length,
                                     k_bar = k_bar,
                                     l = l,
                                     lambda_3 = lambda_3,
                                     spacing = spacing,
                                     gauss_weights = gauss_weights,
                                     same = same,
                                     change_method = change_method,
                                     fixed = TRUE,
                                     cps = cps)

  save(edi.net, file = true.net.filename)
  
  true.net.adj.matrix <- EdiNetToAdjMx(edi.net, q)
  
  # ## Distinct adjacency matrix or matrices of the true net
  # d.true.net.adj.matrix <- edi.net$network
  # 
  # node.names <- paste('G', 1:q, sep='')
  # 
  # ## Num of change pts
  # num.ch.pts <- edi.net$k
  # num.nets <- (num.ch.pts + 1)
  # 
  # rm(num.ch.pts)
  # 
  # for (net.idx in 1:num.nets) {
  #   rownames(d.true.net.adj.matrix[[net.idx]]) <- node.names
  #   colnames(d.true.net.adj.matrix[[net.idx]]) <- node.names
  # }
  # rm(net.idx)
  # 
  # rm(num.nets)
  # 
  # ## Locations of the change pts.
  # ## The i-th elt of vector 'edi.net$epsilon' 
  # ## represents the ending time pt of the j-th 
  # ## time segment.
  # end.pts <- edi.net$epsilon
  # 
  # true.net.adj.matrix <- vector(mode = 'list', length = l)
  # 
  # ## First time point of the first time segment
  # time.seg.start <- 1
  # 
  # ## For each distinct time segment
  # for (time.seg.idx in 1:length(d.true.net.adj.matrix)) {
  #   
  #   ## First time pt of the current time segment
  #   if (time.seg.idx == 1) {
  #     time.seg.start <- 1
  #   } else {
  #     time.seg.start <- (end.pts[(time.seg.idx - 1)] + 1)
  #   }
  #   
  #   ## Last time pt of the current time segment
  #   time.seg.end <- end.pts[time.seg.idx]
  #   
  #   ## The following statement does not work
  #   # true.net.adj.matrix[time.seg.start:time.seg.end] <- d.true.net.adj.matrix[[time.seg.idx]]
  #   
  #   for (time.pt in time.seg.start:time.seg.end) {
  #     true.net.adj.matrix[[time.pt]] <- d.true.net.adj.matrix[[time.seg.idx]]
  #   }
  #   rm(time.pt)
  #   
  # }
  # rm(time.seg.idx)
  
  save(true.net.adj.matrix, file = true.net.adj.mx.filename)
  
  # return(edi.net)
}
##################################################################################################

##################################################################################################
## Goal: Generate synthetic time-series gene expression data time-varying gene regulatory networks
## using R package 'EDISON' version 1.1.1
##
## Input params: 
## 'edi.net': A set of time-varying gene regulatory networks.
## See input param 'net' of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
## 'noise': See input param 'noise' of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
## 'num.time.series': Number of time series to be generated.
##
## Output params:
## See the output params of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
##################################################################################################
GenEdiData <- function(edi.net, noise = 0, num.time.series, output.filename) {
  
  ##------------------------------------------------------------
  ## Begin: Load the required packages
  ##------------------------------------------------------------
  library(EDISON)
  ##------------------------------------------------------------
  ## End: Load the required packa                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ges
  ##------------------------------------------------------------
  
  output.file.conn <- file(output.filename, open = 'wt')
  
  num.nodes <- nrow(edi.net$network[[1]])
  node.names <- paste('G', 1:num.nodes, sep='')
  node.names <- c('', node.names)
  node.names <- matrix(data = node.names, nrow = 1, ncol = length(node.names))
  
  write.table(node.names, file = output.file.conn, append = FALSE, 
              quote = FALSE, sep = "\t", eol = "\n", na = "NA", 
              dec = ".", row.names = FALSE,
              col.names = FALSE)
  
  rm(num.nodes, node.names)
  
  for (ts.idx in 1:num.time.series) {
    
    ## Simulate data using given networks.
    edi.data <- EDISON::simulateNetwork(noise = noise, net = edi.net) 
    
    ## A num.nodes by num.time.pts matrix
    edi.data <- edi.data$sim_data
    
    ## A num.time.pts by num.nodes matrix
    edi.data <- t(edi.data)
    
    if (ts.idx > 1) {
      write.table('\n', file = output.file.conn, append = TRUE, 
                  quote = FALSE, sep = "\t", eol = "", na = "NA", 
                  dec = ".", row.names = FALSE,
                  col.names = FALSE)
    }
    
    write.table(edi.data, file = output.file.conn, append = TRUE, 
                quote = FALSE, sep = "\t", eol = "\n", na = "NA", 
                dec = ".", row.names = TRUE,
                col.names = FALSE)
  }
  rm(ts.idx)
  
  close(output.file.conn)
  
}
##################################################################################################

##################################################################################################
## Goal: Generate synthetic time-series gene expression data time-varying gene regulatory networks
## using R package 'EDISON' version 1.1.1
##
## Input params: 
## 'edi.net': A set of time-varying gene regulatory networks.
## See input param 'net' of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
## 'noise': See input param 'noise' of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
## 'num.time.series': Number of time series to be generated.
##
## Output params:
## See the output params of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
##################################################################################################
GenLargerEdiNet <- function(edi.net, first.net, num.inter.seg.changes, true.net.filename, true.net.adj.mx.filename) {
  
  ## The obj 'edi.net' is a list containing the following five elements.
  ## network: 
  ## epsilon
  ## k
  ## changes
  ## l
  
  num.nodes <- ncol(first.net)
  
  num.cells <- (num.nodes)^2
  
  ## Remove row and col names if any
  first.net <- matrix(first.net, nrow = num.nodes, ncol = num.nodes)
  
  edi.net$network[[1]] <- first.net
  
  num.nets <- length(edi.net$network)
  
  if (num.nets > 1) {
    for (net.idx in 2:num.nets) {
      
      curr.net <- edi.net$network[[(net.idx - 1)]]
      
      cells.to.flip <- sample(1:num.cells, num.inter.seg.changes, replace = FALSE)
      
      for (cell.idx in cells.to.flip) {
        curr.val <- curr.net[cell.idx]
        
        if (curr.val == 0) {
          curr.net[cell.idx] <- 1
        } else if (curr.val == 1) {
          curr.net[cell.idx] <- 0
        }
      }
      rm(cell.idx)
      
      edi.net$network[[net.idx]] <- curr.net
    }
    rm(net.idx)
  }
  
  edi.net$changes[1:num.nets] <- num.inter.seg.changes
  
  save(edi.net, file = true.net.filename)
  
  ## Function 'EdiNetToAdjMx()' is defined in this file only
  true.net.adj.matrix <- EdiNetToAdjMx(edi.net, num.nodes)

  save(true.net.adj.matrix, file = true.net.adj.mx.filename)
  
  # return(edi.net)
}
##################################################################################################

##################################################################################################
## Goal: Generate synthetic time-series gene expression data time-varying gene regulatory networks
## using R package 'EDISON' version 1.1.1
##
## Input params: 
## 'edi.net': A set of time-varying gene regulatory networks.
## See input param 'net' of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
## 'noise': See input param 'noise' of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
## 'num.time.series': Number of time series to be generated.
##
## Output params:
## See the output params of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
##################################################################################################
EdiNetToAdjMx <- function(edi.net, num.nodes) {
  ## Distinct adjacency matrix or matrices of the true net
  d.true.net.adj.matrix <- edi.net$network
  
  node.names <- paste('G', 1:num.nodes, sep='')
  
  ## Num of change pts
  num.ch.pts <- edi.net$k
  num.nets <- (num.ch.pts + 1)
  
  rm(num.ch.pts)
  
  for (net.idx in 1:num.nets) {
    rownames(d.true.net.adj.matrix[[net.idx]]) <- node.names
    colnames(d.true.net.adj.matrix[[net.idx]]) <- node.names
  }
  rm(net.idx)
  
  rm(num.nets)
  
  ## Locations of the change pts.
  ## The i-th elt of vector 'edi.net$epsilon' 
  ## represents the ending time pt of the j-th 
  ## time segment.
  end.pts <- edi.net$epsilon
  
  true.net.adj.matrix <- vector(mode = 'list', length = edi.net$l)
  
  ## First time point of the first time segment
  time.seg.start <- 1
  
  ## For each distinct time segment
  for (time.seg.idx in 1:length(d.true.net.adj.matrix)) {
    
    ## First time pt of the current time segment
    if (time.seg.idx == 1) {
      time.seg.start <- 1
    } else {
      time.seg.start <- (end.pts[(time.seg.idx - 1)] + 1)
    }
    
    ## Last time pt of the current time segment
    time.seg.end <- end.pts[time.seg.idx]
    
    ## The following statement does not work
    # true.net.adj.matrix[time.seg.start:time.seg.end] <- d.true.net.adj.matrix[[time.seg.idx]]
    
    for (time.pt in time.seg.start:time.seg.end) {
      true.net.adj.matrix[[time.pt]] <- d.true.net.adj.matrix[[time.seg.idx]]
    }
    rm(time.pt)
    
  }
  rm(time.seg.idx)
  
  return(true.net.adj.matrix)
}
##################################################################################################