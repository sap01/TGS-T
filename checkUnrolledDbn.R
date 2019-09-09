## Goal: Check whether the given unrolled DBN follows 1st Markov order or not
checkUnrolledDbn <- function(unrolled.DBN.adj.matrix)
{
  for(row.idx in 1:nrow(unrolled.DBN.adj.matrix))
  {
    for(col.idx in 1:ncol(unrolled.DBN.adj.matrix))
    {
      if (unrolled.DBN.adj.matrix[row.idx, col.idx] == 1)
      {
        src.node.name <- rownames(unrolled.DBN.adj.matrix)[row.idx]
        tgt.node.name <- colnames(unrolled.DBN.adj.matrix)[col.idx]
        
        src.time.pt.name <- strsplit(src.node.name, '_', fixed = TRUE)[[1]][2]
        src.time.pt <- substr(src.time.pt.name, 2, nchar(src.time.pt.name))  
        src.time.pt <- as.integer(src.time.pt)
        
        tgt.time.pt.name <- strsplit(tgt.node.name, '_', fixed = TRUE)[[1]][2]
        tgt.time.pt <- substr(tgt.time.pt.name, 2, nchar(tgt.time.pt.name))  
        tgt.time.pt <- as.integer(tgt.time.pt)

        if (src.time.pt != (tgt.time.pt - 1))
        {
          cat(src.node.name,'\t',tgt.node.name,'\n')
        }
      }
    }
  }
}