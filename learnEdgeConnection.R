computEntropy <- function(input.data)
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


learnEdgeConnectionZstat <- function(mut.info.matrix, mi.net.adj.matrix, entropy.matrix, alpha)
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

learnEdgeConnectionRowMedian <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes)
{
  for (rowIdx in 1:num.nodes)
  {
    threshold <- median(mut.info.matrix[colIdx, -colIdx])

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