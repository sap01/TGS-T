## Goal: Generate True net adjacency matrix and save as an R object

## Begin: Specify input params
# input.dir <- '/home/cse/Dropbox/Projects/ProjLbnExtn/ProjData/Dream3InSilicoChallengeData/Size10/DREAM3 gold standards'
# input.dir <- '/home/cse/Dropbox/Projects/ProjLbnExtn/ProjData/Dream3InSilicoChallengeData/Size50/DREAM3 gold standards'
input.dir <- '/home/cse/Dropbox/Projects/ProjLbnExtn/ProjData/Dream3InSilicoChallengeData/Size100/DREAM3 gold standards'

# input.filename <- 'DREAM3GoldStandard_InSilicoSize10_Yeast1.txt'
# input.filename <- 'DREAM3GoldStandard_InSilicoSize50_Yeast1.txt'
input.filename <- 'DREAM3GoldStandard_InSilicoSize100_Yeast1.txt'

input.file <- paste(input.dir, input.filename, sep = '/')

# output.filename <- 'DREAM3GoldStandard_InSilicoSize10_Yeast1_TrueNet.RData'
# output.filename <- 'DREAM3GoldStandard_InSilicoSize50_Yeast1_TrueNet.RData'
output.filename <- 'DREAM3GoldStandard_InSilicoSize100_Yeast1_TrueNet.RData'

output.file <- paste(input.dir, output.filename, sep = '/')
  
# num.nodes <- 10
# num.nodes <- 50
num.nodes <- 100

node.names <- c()
for (node.idx in 1:num.nodes)
{
  new.node.name <- paste('G', as.character(node.idx), sep = '')
  node.names <- c(node.names, new.node.name)
}
## End: Specify input params

# Param 'colClasses' must be given. Otherwise, the node names get converted into factors.
true.net.adj.list <- read.table(input.file, header = FALSE, sep = '\t', 
                                col.names = c('src.node', 'tgt.node', 'isEdge'), 
                                colClasses = c('character', 'character', 'integer'))

## Begin: Convert adjacency list to adjacency matrix
true.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = list(node.names, node.names))
for (row.idx in 1:nrow(true.net.adj.list))
{
  if (true.net.adj.list[row.idx, 'isEdge'] == 1)
  {
    src.node <- true.net.adj.list[row.idx, 'src.node']
    tgt.node <- true.net.adj.list[row.idx, 'tgt.node']
    true.net.adj.matrix[src.node, tgt.node] <- 1
  }
}

## Check number of edges
# length(true.net.adj.matrix[true.net.adj.matrix == 1])

## End: Conver adjacency list to adjacency matrix

save(true.net.adj.matrix, file = output.file)
