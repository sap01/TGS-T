# Goal: Discretize input data into the following two levels:
# Less than wild type => level 1
# Greater than equal to wild type => level 2
discretizeData.2L.wt.l <- function(input.data, input.wt.data.filename)
{
  input.wt.data <- read.table(input.wt.data.filename, header = TRUE, sep="\t")
  
  wt.data <- input.wt.data[1, -1]
  # wt.data[i] contains the wild type expression value of the i^{th} gene
  
  input.data.discr <- input.data
  # input.data.discr[,] <- 0
  
  for (colIdx in 1:ncol(input.data))
  {
    for (rowIdx in 1:nrow(input.data))
    {
      if (input.data.discr[rowIdx, colIdx] < wt.data[colIdx])
      {
        input.data.discr[rowIdx, colIdx] <- 1
      }
      else if (input.data.discr[rowIdx, colIdx] >= wt.data[colIdx])
      {
        input.data.discr[rowIdx, colIdx] <- 2
      }
    }
  }
  
  return(input.data.discr)
}

##########################################################################################################################
# Function's goal: Discretize input data into the following two levels:
# Less than or equal to wild type => level 1
# Greater than wild type => level 2
discretizeData.2L.wt.le <- function(input.data, input.wt.data.filename)
{
  input.wt.data <- read.table(input.wt.data.filename, header = TRUE, sep="\t")
  
  wt.data <- input.wt.data[1, -1]
  # wt.data[i] contains the wild type expression value of the i^{th} gene
  
  input.data.discr <- input.data
  # input.data.discr[,] <- 0
  
  for (colIdx in 1:ncol(input.data))
  {
    for (rowIdx in 1:nrow(input.data))
    {
      if (input.data.discr[rowIdx, colIdx] <= wt.data[colIdx])
      {
        input.data.discr[rowIdx, colIdx] <- 1
      }
      else if (input.data.discr[rowIdx, colIdx] >= wt.data[colIdx])
      {
        input.data.discr[rowIdx, colIdx] <- 2
      }
    }
  }
  
  return(input.data.discr)
}


##########################################################################################################################
## Function's goal: Discretize input data into the following three levels, given a tolerance:
## Let, expression value of gene i in a particular sample is x_i and
## wild type expression value of i is wt(i).
## Level 1: 0 =< x_i < (wt(i) - tolerance). Gene i is down regulated or knocked out.
## Level 2: (wt(i) - tolerance) =< x_i =< (wt(i) + tolerance). Expression of gene i is at a steady-state.
## Level 3: (wt(i) + tolerance) < x_i =< 1. Gene i is up regulated.
##
discretizeData.3L.wt <- function(input.data, input.wt.data.filename, tolerance, num.discr.levels)
{
  input.wt.data <- read.table(input.wt.data.filename, header = TRUE, sep="\t")
  
  wt.data <- input.wt.data[1, -1]
  # wt.data[i] contains the wild type expression value of the i^{th} gene
  
  input.data.discr <- input.data
  # input.data.discr[,] <- 0
  
  for (colIdx in 1:ncol(input.data))
  {
    for (rowIdx in 1:nrow(input.data))
    {
      if ((0 <= input.data.discr[rowIdx, colIdx]) && 
          (input.data.discr[rowIdx, colIdx] < (wt.data[colIdx] - tolerance)))
      {
        input.data.discr[rowIdx, colIdx] <- 1 # Level 1
        
      } else if (((wt.data[colIdx] - tolerance) <= input.data.discr[rowIdx, colIdx]) && 
                 (input.data.discr[rowIdx, colIdx] <= (wt.data[colIdx] + tolerance)))
      {
        input.data.discr[rowIdx, colIdx] <- 2 # Level 2
        
      } else if (((wt.data[colIdx] + tolerance) < input.data.discr[rowIdx, colIdx]) && 
                 (input.data.discr[rowIdx, colIdx] <= 1))
      {
        input.data.discr[rowIdx, colIdx] <- 3 # Level 3
        
      }     
    }
  }
  
  ## bnstruct::BNDataset() requires the discrete levels of a node to be in order starting from level 1. 
  ## Each level index must be in [1, number of discrete levels of that node].
  nodes.discr.sizes <- c()
  for (node.idx in 1:ncol(input.data.discr))
  {
    if (length(unique(input.data.discr[, node.idx])) > 0)
    {
      discr.levels <- sort(unique(input.data.discr[, node.idx]))
      
      for (discr.level.idx in 1:length(discr.levels))
      {
        input.data.discr[input.data.discr[, node.idx] == discr.levels[discr.level.idx], node.idx] <- discr.level.idx
      }
      
      node.discr.size <- length(discr.levels)
      nodes.discr.sizes <- c(nodes.discr.sizes, node.discr.size)
    }
  }
  
  return(list(input.data.discr, nodes.discr.sizes))
}

##########################################################################################################################
## Function's goal: Discretize input data into the following five levels:
## Let, expression value of gene i in a particular sample is x_i and
## wild type expression value of i is wt(i).
## Level 1: x_i = 0. Gene i is knocked out.
## Level 2: 0 < x_i < wt(i). Gene i is down regulated but not knocked out.
## Level 3: x_i = wt(i). Expression of gene i is at a steady-state.
## Level 4: wt(i) < x_i < 1. Gene i is up regulated but not maximally activated.
## Level 5: x_i = 1. Gene i is maximally activated.
discretizeData.5L.wt <- function(input.data, input.wt.data.filename, num.discr.levels)
{
  input.wt.data <- read.table(input.wt.data.filename, header = TRUE, sep="\t")
  
  wt.data <- input.wt.data[1, -1]
  # wt.data[i] contains the wild type expression value of the i^{th} gene
  
  input.data.discr <- input.data
  # input.data.discr[,] <- 0
  
  for (colIdx in 1:ncol(input.data))
  {
    for (rowIdx in 1:nrow(input.data))
    {
      if (input.data.discr[rowIdx, colIdx] == 0)
      {
        input.data.discr[rowIdx, colIdx] <- 1 # Level 1
        
      } else if ((input.data.discr[rowIdx, colIdx] > 0) && (input.data.discr[rowIdx, colIdx] < wt.data[colIdx]))
      {
        input.data.discr[rowIdx, colIdx] <- 2 # Level 2
        
      } else if (input.data.discr[rowIdx, colIdx] == wt.data[colIdx])
      {
        input.data.discr[rowIdx, colIdx] <- 3 # Level 3
        
      } else if ((input.data.discr[rowIdx, colIdx] > wt.data[colIdx]) && (input.data.discr[rowIdx, colIdx] < 1))
      {
        input.data.discr[rowIdx, colIdx] <- 4 # Level 4
        
      } else if (input.data.discr[rowIdx, colIdx] == 1)
      {
        input.data.discr[rowIdx, colIdx] <- 5 # Level 5
        
      }      
    }
  }
  
  ## bnstruct::BNDataset() requires the discrete levels of a node to be in order starting from level 1. 
  ## Each level index must be in [1, number of discrete levels of that node].
  for (node.idx in 1:ncol(input.data.discr))
  {
    if ((length(unique(input.data.discr[, node.idx])) > 0) &&
        (length(unique(input.data.discr[, node.idx])) != num.discr.levels))
    {
      discr.levels <- sort(unique(input.data.discr[, node.idx]))
      
      for (discr.level.idx in 1:length(discr.levels))
      {
        input.data.discr[input.data.discr[, node.idx] == discr.levels[discr.level.idx], node.idx] <- discr.level.idx
      }
    }
  }
  
  return(input.data.discr)
}

##########################################################################################################################
## Function's goal: Discretize input data into the following two levels as done in [Supporting Information p. 2, 1].
## For each gene, the expression values are first sorted; then the top 2
## extreme values in either end of the sorted list are discarded; last,
## the median of the remaining values is used as the threshold above
## which the value is binarized as 'Level 2' and 'Level 1' otherwise. Here, 2 means
## the expression of a gene is up-regulated, and 1 means down-regulated.
## input.data: rows = samples, cols = variables.
## Reference:
## 1. Ahmed, Amr, and Eric P. Xing. "Recovering time-varying networks of dependencies in social and biological studies." 
## Proceedings of the National Academy of Sciences 106.29 (2009): 11878-11883.]
##
discretizeData.2L.Tesla <- function(input.data)
{
  ## For each variable
  for (var.idx in 1:ncol(input.data))
  {
    ## Sorted in asc order
    conts.vals.sorted <- sort(input.data[, var.idx])
    
    if (length(conts.vals.sorted) > 4)
    {
      ## Discard lowest two values
      conts.vals.sorted <- conts.vals.sorted[3:length(conts.vals.sorted)]
      
      ## Discard highest two values
      conts.vals.sorted <- conts.vals.sorted[1:(length(conts.vals.sorted) - 2)]
    }
    
    ## Compute median
    discr.threshold <- median(conts.vals.sorted)
    
    ## For each sample value of the current variable
    for (sample.idx in 1:nrow(input.data))
    {
      if (input.data[sample.idx, var.idx] > discr.threshold)
      {
        input.data[sample.idx, var.idx] <- 2
      }
      else
      {
        input.data[sample.idx, var.idx] <- 1
      }
    }
  }
  
  return(input.data)
}
##########################################################################################################

##########################################################################################################
## Goal: Discretize time-series input data into two levels {1, 2} 
## using the following strategy. 
## The input data may or may not have multiple time series.
## For each gene, its expression value(s) at the first time point(s)
## are assinged level 1. 
## Level 1: If it belongs to the first time point or if it is lower than that 
## of the previous time point in the same time series.
## Level 2: If it is not lower than that of the previous time point in the same time series.
##########################################################################################################
discretizeData.2L.1 <- function(input.data, num.timepts) {
  
  ## Num of time series
  num.ts <- (nrow(input.data) / num.timepts)
  
  for (ts.idx in 1:num.ts) {
    
    ## Last and first rows of the current time series
    ## in the 'input.data'
    last.row <- (num.timepts * ts.idx)
    first.row <- (last.row - num.timepts + 1)
    
    ## Current time series  
    curr.ts <- input.data[first.row:last.row, ]
    
    ## Initialize discretized version of the 
    ## current time series
    curr.ts.d <- curr.ts
    
    ## Values at the first time point are assigned 
    ## level 1
    curr.ts.d[1, ] <- 1
    
    for (row.idx in 2:nrow(curr.ts)) {
      for (col.idx in 1:ncol(curr.ts)) {
        if (curr.ts[row.idx, col.idx] < curr.ts[(row.idx - 1), col.idx]) {
          curr.ts.d[row.idx, col.idx] <- 1
          
        } else {
          curr.ts.d[row.idx, col.idx] <- 2
        }
      }
      rm(col.idx)
    }
    rm(row.idx)
    
    ## Replace continuous data with discretized data
    input.data[first.row:last.row, ] <- curr.ts.d
  }
  rm(ts.idx)
  
  ## Return discretized version of the input data
  return(input.data)
}
########################################################################################################## 

##########################################################################################################
## Goal: Discretize time-series input data into three levels {1, 2, 3} 
## using the following strategy. 
## The input data may or may not have multiple time series.
## For each gene, its expression value(s) at the first time point(s)
## are assinged level 2. 
## Level 2: If it belongs to the first time point or if it is same as that 
## of the previous time point in the same time series.
## Level 3: If it is higher than that of the previous time point in the same time series.
## Level 1: If it is lower than that of the previous time point in the same time series.
##########################################################################################################
discretizeData.3L.1 <- function(input.data, num.timepts) {
  
  ## Num of time series
  num.ts <- (nrow(input.data) / num.timepts)
  
  for (ts.idx in 1:num.ts) {
    
    ## Last and first rows of the current time series
    ## in the 'input.data'
    last.row <- (num.timepts * ts.idx)
    first.row <- (last.row - num.timepts + 1)
    
    ## Current time series  
    curr.ts <- input.data[first.row:last.row, ]
    
    ## Initialize discretized version of the 
    ## current time series
    curr.ts.d <- curr.ts
    
    ## Values at the first time point are assigned 
    ## level 2
    curr.ts.d[1, ] <- 2
    
    for (row.idx in 2:nrow(curr.ts)) {
      for (col.idx in 1:ncol(curr.ts)) {
        if (curr.ts[row.idx, col.idx] < curr.ts[(row.idx - 1), col.idx]) {
          curr.ts.d[row.idx, col.idx] <- 1
          
        } else if (curr.ts[row.idx, col.idx] == curr.ts[(row.idx - 1), col.idx]) {
          curr.ts.d[row.idx, col.idx] <- 2
          
        } else if (curr.ts[row.idx, col.idx] > curr.ts[(row.idx - 1), col.idx]) {
          
          curr.ts.d[row.idx, col.idx] <- 3
        }
      }
      rm(col.idx)
    }
    rm(row.idx)
    
    ## Replace continuous data with discretized data
    input.data[first.row:last.row, ] <- curr.ts.d
  }
  rm(ts.idx)
  
  ## Return discretized version of the input data
  return(input.data)
}
########################################################################################################## 

##########################################################################################################
## Goal: Discretize time-series input data into three levels {1, 2, 3} 
## using the following strategy. 
## The input data may or may not have multiple time series.
## For each gene, its expression value(s) at the first time point(s)
## are assinged level 2. 
## Level 2: If it belongs to the first time point or if it is equal to 
## (e - e * d) where e = expression at the prev time pt in the same time series and 
## d = rate of mRNA decay corr. to the concerned gene.
## Level 3: If it is higher than (e - e * d).
## Level 1: If it is lower than (e - e * d).
##########################################################################################################
discretizeData.3L.2 <- function(input.data, num.timepts) {
  
  ## Num of time series
  num.ts <- (nrow(input.data) / num.timepts)
  
  
  num.genes <- ncol(input.data)
  
  for (gene.idx in 1:num.genes) {
    
    ## Determine the rate of mRNA decay for every gene
    ###########################################################
    ## Initialize the minimum mRNA decay rate of the given gene
    min.mrna.decay.rate <- 0
    
    for (ts.idx in 1:num.ts) {
      
      ## Last and first rows of the current time series
      ## in the 'input.data'
      last.row <- (num.timepts * ts.idx)
      first.row <- (last.row - num.timepts + 1)
      
      ## Current time series  for the given gene.
      ## It sd be a vector.
      curr.ts <- input.data[first.row:last.row, gene.idx]
      rm(first.row, last.row)
      
      for (curr.time.pt in 2:num.timepts) {
        
        curr.mrna.decay.rate <- 0
        
        if (curr.ts[curr.time.pt] < curr.ts[(curr.time.pt - 1)]) {
          curr.mrna.decay.rate <- ((curr.ts[(curr.time.pt - 1)] - curr.ts[curr.time.pt]) / curr.ts[(curr.time.pt - 1)])
        }
        
        if (curr.mrna.decay.rate < min.mrna.decay.rate) {
          min.mrna.decay.rate <- curr.mrna.decay.rate
        }
      }
      rm(curr.time.pt)
      
    }
    rm(ts.idx)
    ###########################################################
    
    ## Discretize data for the given gene
    ###########################################################
    for (ts.idx in 1:num.ts) {
      
      ## Last and first rows of the current time series
      ## in the 'input.data'
      last.row <- (num.timepts * ts.idx)
      first.row <- (last.row - num.timepts + 1)
      
      ## Current time series  
      curr.ts <- input.data[first.row:last.row, gene.idx]
      rm(first.row, last.row)
      
      ## Initialize discretized version of the 
      ## current time series for the given gene
      curr.ts.d <- curr.ts
      
      ## Value at the first time point are assigned 
      ## level 2 for the given gene
      curr.ts.d[1] <- 2
      
      for (curr.time.pt in 2:num.timepts) {
        
        ## Expected expression of the given gene in the current time pt
        ## assuming no up- or down-regulation has happened
        expr.no.reg <- (curr.ts[(curr.time.pt - 1)] * min.mrna.decay.rate)
        expr.no.reg <- (curr.ts[(curr.time.pt - 1)] - expr.no.reg)
        
        if (curr.ts[curr.time.pt] < expr.no.reg) {
          ## Down-regulated
          curr.ts.d[curr.time.pt] <- 1
        } else if (curr.ts[curr.time.pt] == expr.no.reg) {
          ## Neither down-regulated nor up-regulated
          curr.ts.d[curr.time.pt] <- 2
        } else if (curr.ts[curr.time.pt] > expr.no.reg) {
          ## ## Up-regulated
          curr.ts.d[curr.time.pt] <- 3
        }
      }
      rm(time.pt)
      
      ## Replace continuous data with discretized data for the given gene
      input.data[first.row:last.row, gene.idx] <- curr.ts.d
      
    }
    rm(ts.idx)
    ###########################################################
    
  }
  rm(gene.idx)
  
  rm(num.genes)
  
  ## Return discretized version of the input data
  return(input.data)
}
########################################################################################################## 