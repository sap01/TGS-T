## Goal: Given {an unrolled DBN adjacency matrix, a rolling method, the true rolled DBN adjacency matrix}, 
## roll the unrolled DBN adjacency matrix into a rolled DBN adjacency matrix and calculate its accuracy
## w.r.t. the given true rolled DBN adjacency matrix.
## calcPerfUnrolledDbn(10, c('G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10'), 21, unrolled.DBN.adj.matrix, 'any', FALSE)
# todo: add 'true.net.filename' or 'targetNet' as an input arg
calcPerfUnrolledDbn <- function(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, roll.method, allow.self.loop)
{
  setwd('/home/saptarshi/R/R-3.3.2/projects/repoprojldbn')
  
  # Input file for true network
  true.net.filename <- paste(getwd(), 'asset/LbnDataset10.RData', sep = '/')
  load(file = true.net.filename) # loads object 'dataset10TrueNet'
  targetNet <- dataset10TrueNet
  
  source('rollDbn.R')
  rolled.DBN.adj.matrix <- 
    rollDbn(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, roll.method, allow.self.loop)
  
  di.net.adj.matrix <- rolled.DBN.adj.matrix
  inferredNet <- di.net.adj.matrix
  
  source('calcperfdinet.R')
  ResultVsTrue <- calcperfdinet(inferredNet, targetNet, Result, num.nodes)
  writeLines('Result LDBN vs True = \n')
  print(ResultVsTrue)
  
}