# Goal of this program: Implement the Local Dynamic Bayesian Network (LDBN) algorithm

# Input: A (n * N) numerical matrix where
# n = Number of genes
# N = Number of samples
# (i, j)^{th} cell contains a real value representing
# expression of the i^{th} gene in the j^{th} sample.

# Output: A single gene regulatory network on n genes.

# ToDo: Clarify the following points from the LBN authors.
# Q. How is the data discretized?
# Q. How is mutual information matrix calculated from the discretized data?
# Q. What search strategy is used to infer the local Bayesian structure?
# Q. Which estimator is used to calculate conditional mutual info?
# Q. Can You provide the output Networks from your implmentation ?

# ToDo: Correct the indentation of each line with the help of RStudio
# ToDo: Use the same variable names as in LBN paper
# ToDo: Use more self-explanatory variable names, like - 'cmi' instead of 'I'
# ToDo: Use comment generously so that the code is readable to your reviewers
# and readers who has never met you

# Remove all workspace objects
rm(list = ls())

#------------------------------------------------------------
# Begin: Load the Required Libraries
#------------------------------------------------------------
library(minet)
library(Rgraphviz)
library(bnlearn)
library(infotheo)
library(catnet)
library(igraph)
require(RoughSets)
library(gRain)
library(randomForest) # Required for GENIE3
#------------------------------------------------------------
# End: Load the Required Libraries
#------------------------------------------------------------

#------------------------------------------------------------
# Begin: Load the Required User-defined Functions
#------------------------------------------------------------
source('discretizeData.R')
source('computeCmi.R')
source('learnMiNetStruct.R') # 'learnEdgeConnection.R' is renamed to 'learnMiNetStruct.R'

# source('learnDbnStruct3dSer.R') # Serial programming
# source('learnDbnStruct3dPar.R') # Parallel programming
source('test_learnDbnStruct3dPar.R') # Parallel programming

source('calcperfdinet.R')

source('learnCmiNetStruct.R')
source('CompareNet.R')
#------------------------------------------------------------
# End: Load the Required User-defined Functions
#------------------------------------------------------------

#------------------------------------------------------------
# Begin: Input User Defined Params
#------------------------------------------------------------
# Input file for time-course gene expression data
input.data.filename <- paste(getwd(), 'asset/InSilicoSize10-Yeast1-trajectories.tsv', sep = '/')

# Number of time pts in the time-course data 
num.timepts <- 21

# Input file for wild type expression of genes 
input.wt.data.filename <- paste(getwd(), 'asset/InSilicoSize10-Yeast1-null-mutants.tsv', sep = '/')

# Input file for true network
true.net.filename <- paste(getwd(), 'asset/LbnDataset10.RData', sep = '/')

# If input data is discretized, then number of discrete levels for each variable.
num.discr.levels <- 2

# significance level for CMI
alpha <- 0.05

# save output in a file
# "asset/output<YYYY><MM><DD><hh><mm><ss>.txt". Eg. "asset/output20170419153235.txt"
output.filename <- paste('asset/output', format(Sys.time(), "%Y%m%d%H%M%S"), '.txt', sep = '')
output.filename <- paste(getwd(), output.filename, sep = '/')
#------------------------------------------------------------
# End: Input User Defined Params
#------------------------------------------------------------

#------------------------------------------------------------
# Begin: Main Program
#------------------------------------------------------------

output.file.conn <- file(output.filename, open = "wt")
sink(output.file.conn)

#------------------------------------------------------------
# Begin: Read True Net and LBN Original Implementation's  Output Net 
#------------------------------------------------------------
load(file = true.net.filename)

# plot(as(dataset10TrueNet ,"graphNEL"))
# plot(as(dataset10LbnNet ,"graphNEL"))
#------------------------------------------------------------
# End: Read True Net and LBN Original Implementation's  Output Net 
#------------------------------------------------------------

#------------------------------------------------------------
# Begin: Read Drosophila Life Cycle Data File
#------------------------------------------------------------
# node.names <- read.csv('genenames.txt', header = FALSE)
# node.names <- t(node.names)
 
# # Cols are nodes, rows are time points
# input.data <- read.table('drosophila.data', header = FALSE, col.names = node.names)
# num.nodes <- ncol(input.data)
#------------------------------------------------------------
# End: Read Drosophila Life Cycle Data File
#------------------------------------------------------------

#------------------------------------------------------------
# Begin: Read DREAM3 In Silico Challenge Data File
#------------------------------------------------------------
input.data <- read.table(input.data.filename, header = TRUE, sep="\t")

timepts.names <- input.data[1:num.timepts, 1]

# Remove first col i.e. time point names
input.data <- input.data[, -1]

node.names <- colnames(input.data)

num.nodes <- ncol(input.data)

num.samples.per.timept <- (nrow(input.data) / num.timepts)

input.data.discr <- discretizeData(input.data, input.wt.data.filename)

#------------------------------------------------------------
# End: Read DREAM3 In Silico Challenge Data File
#------------------------------------------------------------

#------------------------------------------------------------
# Begin: Learn MI Net Structure
#------------------------------------------------------------

# Initialize mutual information matrix with zeroes
mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, dimnames = c(list(node.names), list(node.names)))

# Build mutual information matrix
for (colIdx in 1:(num.nodes - 1))
{
  for (colIdx2 in (colIdx + 1):num.nodes)
  {
    mut.info <- computeCmi(input.data.discr[, colIdx], input.data.discr[, colIdx2])
    mut.info.matrix[colIdx, colIdx2] <- mut.info
    mut.info.matrix[colIdx2, colIdx] <- mut.info
  }
}
writeLines('mut.info.matrix = \n')
print(mut.info.matrix)

start.time <- proc.time() # start the timer

# Initialize mutual information network
mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, dimnames = c(list(node.names), list(node.names)))


# entropy.matrix <- computEntropy(input.data.discr) #----Verify the name
# mi.net.adj.matrix <- learnMiNetStructZstat(mut.info.matrix, mi.net.adj.matrix, entropy.matrix, alpha)

# mi.net.adj.matrix <- learnMiNetStructClr(mut.info.matrix, mi.net.adj.matrix, num.nodes)

mi.net.adj.matrix <- learnMiNetStructClr(mut.info.matrix, mi.net.adj.matrix, num.nodes)

writeLines('\n mi.net.adj.matrix = \n')
print(mi.net.adj.matrix)
save(mi.net.adj.matrix, file = paste(getwd(), 'asset/mi.net.adj.matrix.RData', sep = '/'))

# stop('MI net struct learning completed.')

# Rgraphviz::plot(as(mi.net.adj.matrix, 'graphNEL'))
#------------------------------------------------------------
# End: Learn MI Net Structure
#------------------------------------------------------------

#------------------------------------------------------------
# Begin: Learn Directed Network Structure
#------------------------------------------------------------

# load('dream3.yeast1.size10.trajectory.discr2.3D.RData') # loads 'DBN.input.data'

input.data.discr.matrix <- data.matrix(input.data.discr)

input.data.discr.3D <- array(NA, c(num.timepts, num.nodes, num.samples.per.timept),
                         dimnames = c(list(timepts.names), list(node.names),
                                      list(1:num.samples.per.timept)))

for (sampleIdx in 1:num.samples.per.timept)
{
  start.rowIdx <- (1 + (num.timepts * (sampleIdx - 1)))
  end.rowIIdx <- (num.timepts * sampleIdx)
  input.data.discr.3D[ , , sampleIdx] <- input.data.discr.matrix[start.rowIdx:end.rowIIdx, ]
}

# dream3.yeast1.size10.trajectory.3D.data <- input.data.3D
# save(dream3.yeast1.size10.trajectory.3D.data, file = 'dream3.yeast1.size10.trajectory.3D.RData')
# load('dream3.yeast1.size10.trajectory.3D.RData')
# input.data.3D <- dream3.yeast1.size10.trajectory.3D.data

# di.net.adj.matrix <- learnDbnStruct3dSer(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels,
#                                       num.nodes, num.timepts)
di.net.adj.matrix <- test_learnDbnStruct3dPar(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels,
                                         num.nodes, num.timepts)

writeLines('\n di.net.adj.matrix = \n')
print(di.net.adj.matrix)
save(di.net.adj.matrix, file = paste(getwd(), 'asset/di.net.adj.matrix.RData', sep = '/'))
#------------------------------------------------------------
# End: Learn LDBN Structures
#------------------------------------------------------------

writeLines('\n Iteration 1: After DBN = \n')
inferredNet <- di.net.adj.matrix
targetNet <- dataset10TrueNet
ResultVsTrue <- calcperfdinet(inferredNet, targetNet, Result, num.nodes)
writeLines('Result LDBN vs True = \n')
print(ResultVsTrue)

elapsed.time <- (proc.time() - start.time) # Stop the timer
writeLines('elapsed.time = \n')
print(elapsed.time)

sink()
stop('First iteration is over.')

#------------------------------------------------------------
# Begin: Learn CMI Net Structure
#------------------------------------------------------------
beta <- 0.03

cmi.net.adj.matrix <- learnCmiNetStruct(di.net.adj.matrix, input.data, node.names, beta)

#TEMP--------------cmi.net.adj.matrix.prev <- cmi.net.adj.matrix
#------------------------------------------------------------
# End: Learn CMI Net Structure
#------------------------------------------------------------
node.nbr.matrix <- mi.net.adj.matrix
# #------------------------------
# # Begin: Second Iteration Till Convergence 
# #------------------------------
 repeat
 {
   comp <- 0
   comp <- CompareNet(di.net.adj.matrix,cmi.net.adj.matrix,num.node)
   if(comp == 1)
	{
      break
	}
	else
	{ 
		node.nbr.matrix[1:num.nodes,1:num.nodes] <- 0 # will be use for storing the neighbour of the node.	  
		for (rowIdx in 1:(num.nodes - 1))
			{
			for (colIdx in (rowIdx + 1):num.nodes)
				{
				 if(cmi.net.adj.matrix[rowIdx,colIdx] == 1) 
					{
						node.nbr.matrix[rowIdx,colIdx] <- 1
						node.nbr.matrix[colIdx,rowIdx] <- 1
					}
				}
			}		
	#-----------kNN for k = 2--------------*
    node.2NN.matrix <- node.nbr.matrix
   
	for(i in 1:num.nodes) # for each node 'i'
		{
		for(j in 1:num.nodes) # for each node 'j'
			{
			if(node.nbr.matrix[i,j] == 1)
				{
				for(k in 1:num.nodes) # for each node 'k'
					{
					# If 'k' is a 2nd nearest neighbour of 'i'
					if(node.nbr.matrix[j,k] == 1 && k !=i)
						{
							node.2NN.matrix[i,k] <- 1;
						}
					}			
				}
			}
		}
	    di.net.adj.matrix <- learnDbnStruct3D(discr.input.data, node.2NN.matrix)
        cmi.net.adj.matrix <- learnCmiNetStruct(di.net.adj.matrix, input.data, node.names, beta)		
	}
 }
#------------------------------
# End: Second Iteration Till Convergence 
#------------------------------

plot(as(final_net, "graphNEL"))
targetNet <- #_______________________________True network to assign by the user.
Result <- calcperfdinet(cmi.net.adj.matrix, targetNet, Result, node.names) 
print(Result)

#------------------------------------------------------------
# Begin: Save R Session Info in a File
#------------------------------------------------------------
sink("sessionInfo.txt")
sessionInfo()
sink()
#------------------------------------------------------------
# End: Save R Session Info in a File
#------------------------------------------------------------
#------------------------------------------------------------
# End: Main Program
#------------------------------------------------------------

#------------------------------------------------------------
# Begin: References
#------------------------------------------------------------
# 1. Liu, Fei, et al. "Inference of Gene Regulatory Network Based on Local Bayesian Networks." PLoS Comput Biol 12.8 (2016): e1005024.
#------------------------------------------------------------
# End: References
#------------------------------------------------------------
