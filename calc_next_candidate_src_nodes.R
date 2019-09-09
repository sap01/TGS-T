## Calculate next candidate source nodes


# x <- c('a', 'b', 'c')
# y <- c(TRUE, FALSE, TRUE)
# z <- x[y]
# print(z)
# length(z)

next.candidate.src.node.string <- calc.next.candidate.parent.string(curr.candidate.src.node.string, num.candidate.src.nodes)
next.candidate.src.node.set <- candidate.src.node.names[next.candidate.src.node.string]

calc.next.candidate.parent.string <- function(curr.candidate.src.node.string, num.candidate.src.nodes) {
  
  next.candidate.src.node.string <- vector(n)
  
  least.sig.bit <- curr.candidate.src.node.string[num.candidate.src.nodes]
  
  carry.bit <- FALSE
  
  if (least.sig.bit) {
    next.candidate.src.node.string[num.candidate.src.nodes] <- FALSE
    carry.bit <- TRUE
  } else {
    next.candidate.src.node.string[num.candidate.src.nodes] <- TRUE
  }
  
  for (loop.idx in (num.candidate.src.nodes - 1):1) {
    curr.bit <- curr.candidate.src.node.string[loop.idx]
    
    if (curr.bit & carry.bit) {
      next.candidate.src.node.string[loop.idx] <- FALSE
      carry.bit <- TRUE
    } else if (!curr.bit & !carry.bit) {
      next.candidate.src.node.string[loop.idx] <- FALSE
      carry.bit <- FALSE
    } else {
      next.candidate.src.node.string[loop.idx] <- TRUE
      carry.bit <- FALSE
    }
    
  }
  rm(loop.idx)
}

