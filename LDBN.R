# Goal of this program: Implement the Local Dynamic Bayesian Network (LDBN) algorithm
# Run in 'grni' (172.16.112.130) server using the following commands:
# $ cd /home/saptarshi/R/R-3.3.2/projects/repoTGS/  
# $ nohup time /home/saptarshi/R/R-3.3.2/bin/Rscript /home/saptarshi/R/R-3.3.2/projects/repoTGS/LDBN.R &

## Begin: Input data format
## Input data format for cross-sectional dataset:
## A (N * V) numerical matrix where
## cols: V = Number of genes.
## rows: N = Number of samples; N > 1.
## (n, v)^{th} cell contains a real value representing
## expression of the v^{th} gene in the n^{th} sample.

## Input data format for time-series dataset:
## A (T * V * S) numerical matrix where
## 1st dim: T = Number of time points per time series; T > 1. Each time series must have
## equal no. of time points. todo: relax this restriction.
## 2nd dim: V = Number of genes.
## 3rd dim: S = Number of time series; T > 1.
## The first T row must corr. to the first time series,
## The second T row must corr. to the second time series, and so on.
## End: Input data format

# Output: A single gene regulatory network on the given genes.

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
setwd('/home/saptarshi/R/R-3.3.2/projects/repoprojldbn')

source('discretizeData.R')
source('computeCmi.R')
source('learnMiNetStruct.R') # 'learnEdgeConnection.R' is renamed to 'learnMiNetStruct.R'

# source('learnDbnStruct3dSer.R') # Serial programming
source('learnDbnStruct3dParDeg1.R') # Parallel programming with degree of parallelism 1
# source('learnDbnStruct3dParDeg2.R') # Parallel programming with degree of parallelism 2

source('calcPerfDiNet.R')

source('learnCmiNetStruct.R')
source('CompareNet.R')

source('rollDbn.R')

source('RDataToCytoscape.R')
#------------------------------------------------------------
# End: Load the Required User-defined Functions
#------------------------------------------------------------

#------------------------------------------------------------
# Begin: Input User Defined Params
#------------------------------------------------------------
## Input file for time-course gene expression data
## Yeast 1, 10 genes, 21 time pts, 4 sampels per gene per time pt, 10 di edges, noisy.
# input.data.filename <- paste(getwd(), 'asset/InSilicoSize10-Yeast1-trajectories.tsv', sep = '/')
## Yeast 1, 10 genes, 21 time pts, 4 sampels per gene per time pt, 10 di edges, noiseless.
# input.data.filename <- paste(getwd(), 'asset/InSilicoSize10-Yeast1-nonoise-trajectories.tsv', sep = '/')
## Yeast 1, 50 genes, 21 time pts, 23 samples per gene per time pt, 77 di edges, noisy.
# input.data.filename <- paste(getwd(), 'asset/InSilicoSize50-Yeast1-trajectories.tsv', sep = '/')
## Yeast 1, 50 genes, 21 time pts, 23 samples per gene per time pt, 77 di edges, noiseless.
# input.data.filename <- paste(getwd(), 'asset/InSilicoSize50-Yeast1-nonoise-trajectories.tsv', sep = '/')
## Yeast 1, 100 genes, 21 time pts, 46 samples per gene per time pt, 166 di edges, noisy.
# input.data.filename <- paste(getwd(), 'asset/InSilicoSize100-Yeast1-trajectories.tsv', sep = '/')
## Yeast 1, 100 genes, 21 time pts, 46 samples per gene per time pt, 166 di edges, noiseless.
# input.data.filename <- paste(getwd(), 'asset/InSilicoSize100-Yeast1-nonoise-trajectories.tsv', sep = '/')
## Drosophila melanogaster Life cycle data (DmLc), 4028 genes, 22 time pts, 3 sampels per gene per time pt, discretized.
# input.data.filename <- paste(getwd(), 'asset/DmLc.RData', sep = '/')
# input.data.filename <- paste(getwd(), 'asset/DmLcE.RData', sep = '/')
# input.data.filename <- paste(getwd(), 'asset/DmLcL.RData', sep = '/')
# input.data.filename <- paste(getwd(), 'asset/DmLcP.RData', sep = '/')
# input.data.filename <- paste(getwd(), 'asset/DmLcA.RData', sep = '/')
# input.data.filename <- paste(getwd(), 'asset/DmLc3E.RData', sep = '/')
# input.data.filename <- paste(getwd(), 'asset/DmLc3L.RData', sep = '/')
# input.data.filename <- paste(getwd(), 'asset/DmLc3P.RData', sep = '/')
input.data.filename <- paste(getwd(), 'asset/DmLc3A.RData', sep = '/')

## Number of time pts in the time-course data 
num.timepts <- 2
# num.timepts <- 3
# num.timepts <- 6
# num.timepts <- 10
# num.timepts <- 15
# num.timepts <- 21
# num.timepts <- 22

## Input file for wild type expression of genes 
# input.wt.data.filename <- paste(getwd(), 'asset/InSilicoSize10-Yeast1-null-mutants.tsv', sep = '/')
# input.wt.data.filename <- paste(getwd(), 'asset/InSilicoSize10-Yeast1-nonoise-null-mutants.tsv', sep = '/')
# input.wt.data.filename <- paste(getwd(), 'asset/InSilicoSize50-Yeast1-null-mutants.tsv', sep = '/')
# input.wt.data.filename <- paste(getwd(), 'asset/InSilicoSize50-Yeast1-nonoise-null-mutants.tsv', sep = '/')
# input.wt.data.filename <- paste(getwd(), 'asset/InSilicoSize100-Yeast1-null-mutants.tsv', sep = '/')
# input.wt.data.filename <- paste(getwd(), 'asset/InSilicoSize100-Yeast1-nonoise-null-mutants.tsv', sep = '/')

## Input file for true network. Loads R obj 'true.net.adj.matrix'
# true.net.filename <- paste(getwd(), 'asset/DREAM3GoldStandard_InSilicoSize10_Yeast1_TrueNet.RData', sep = '/')
# true.net.filename <- paste(getwd(), 'asset/DREAM3GoldStandard_InSilicoSize50_Yeast1_TrueNet.RData', sep = '/')
# true.net.filename <- paste(getwd(), 'asset/DREAM3GoldStandard_InSilicoSize100_Yeast1_TrueNet.RData', sep = '/')

# If input data is discretized, then number of discrete levels for each variable.
num.discr.levels <- 2

max.fanin <- 14

## significance level for CMI
# alpha <- 0.05

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
# load(file = true.net.filename)

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

load(input.data.filename) ## Loads DmLc
input.data <- DmLc
# input.data <- DmLc[, 1:10] ## Run the program for only first few genes in the input dataset
rm(DmLc)
timepts.names <- 1:num.timepts
#------------------------------------------------------------
# End: Read Drosophila Life Cycle Data File
#------------------------------------------------------------

#------------------------------------------------------------
# Begin: Read DREAM3 In Silico Challenge Data File
#------------------------------------------------------------
# input.data <- read.table(input.data.filename, header = TRUE, sep="\t")

# timepts.names <- input.data[1:num.timepts, 1]

## Remove first col i.e. time point names
# input.data <- input.data[, -1]

## Begin: Replace original node names with {v1, v2, ..., vV}
## V = total number of node in the input data.
## The replacement is crucial as unexpected characters in the node names
## may generate an error in further computation.

## Save the original node names in case they are required in future
orig.node.names <- colnames(input.data)

node.names <- c()
for (col.idx in 1:ncol(input.data))
{
  new.node.name <- paste('v', as.character(col.idx), sep = '')
  node.names <- c(node.names, new.node.name)
}
rm(col.idx)
colnames(input.data) <- node.names
## End: Replace original node names with {v1, v2, ..., vV}
## V = total number of node in the input data.


num.nodes <- ncol(input.data)

## Max fanin must be restricted to 14. Because, it is empirically observed that bnstruct::learn.network()
## function can learn a BN with upto 15 nodes without segmentation fault, given a 32 GB main memory. A max 
## fanin restriction of 14 limits the number of nodes in the to-be-learnt BN to 15 (1 target node and a 
## maximum of 14 candidate parents).
max.fanin <- min(num.nodes, 14)

num.samples.per.timept <- (nrow(input.data) / num.timepts)

input.data.discr <- input.data
# input.data.discr <- discretizeData.2L.wt.l(input.data, input.wt.data.filename)
# input.data.discr <- discretizeData.2L.Tesla(input.data)
# save(input.data.discr, file = paste(getwd(), 'asset/input.data.discr.RData', sep = '/'))
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
# writeLines('mut.info.matrix = \n')
# print(mut.info.matrix)
save(mut.info.matrix, file = paste(getwd(), 'asset/mut.info.matrix.RData', sep = '/'))

start.time <- proc.time() # start the timer

# Initialize mutual information network
mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, dimnames = c(list(node.names), list(node.names)))


# entropy.matrix <- computEntropy(input.data.discr) #----Verify the name
# mi.net.adj.matrix <- learnMiNetStructZstat(mut.info.matrix, mi.net.adj.matrix, entropy.matrix, alpha)
# mi.net.adj.matrix <- learnMiNetStructClr(mut.info.matrix, mi.net.adj.matrix, num.nodes)
mi.net.adj.matrix <- learnMiNetStructClrMfi(mut.info.matrix, mi.net.adj.matrix, num.nodes, max.fanin)

elapsed.time <- (proc.time() - start.time) # Check time taken by CLR
writeLines('elapsed.time just after CLR step= \n')
print(elapsed.time)

writeLines('\n mi.net.adj.matrix = \n')
# print(mi.net.adj.matrix)
save(mi.net.adj.matrix, file = paste(getwd(), 'asset/mi.net.adj.matrix.RData', sep = '/'))
rm(mut.info.matrix)
## Check max number of nbrs for a node in the MI net
# writeLines('\n Max number of nbrs = \n')
# print(max(colSums(mi.net.adj.matrix)))

## Identify which nodes have the max number of nbrs
# writeLines('\n Nodes having max number of nbrs = \n')
# print(colnames(mi.net.adj.matrix[, colSums(mi.net.adj.matrix) == max(colSums(mi.net.adj.matrix))]))

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

# Learn the unrolled DBN adj matrix
# unrolled.DBN.adj.matrix <- learnDbnStructLayer3dSer(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels,
#                                                         num.nodes, num.timepts)
# unrolled.DBN.adj.matrix <- learnDbnStruct3dParDeg1(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels,
#                                                     num.nodes, num.timepts)
# unrolled.DBN.adj.matrix <- learnDbnStructLayer3dParDeg1(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels,
#                                                       num.nodes, num.timepts)
# unrolled.DBN.adj.matrix <- learnDbnStructMo1Layer3dParDeg1(input.data.discr.3D, mi.net.adj.matrix, 
#                                                            num.discr.levels, num.nodes, 
#                                                            num.timepts, max.fanin)
unrolled.DBN.adj.matrix.list <- learnDbnStructMo1Layer3dParDeg1_v2(input.data.discr.3D, mi.net.adj.matrix, 
                                               num.discr.levels, num.nodes, num.timepts, max.fanin, 
                                               node.names)
# unrolled.DBN.adj.matrix <- learnDbnStruct3dParDeg2(input.data.discr.3D, mi.net.adj.matrix, num.discr.levels,
#                                          num.nodes, num.timepts)

# After unrolled.DBN.adj.matrix.RData is saved, the following command can be used to save it in .tsv file format
# write.table(unrolled.DBN.adj.matrix.RData, file = 'unrolled.DBN.adj.matrix.tsv', sep = '\t', row.names = TRUE, col.names = TRUE)
# save(unrolled.DBN.adj.matrix, file = paste(getwd(), 'asset/unrolled.DBN.adj.matrix.RData', sep = '/'))
save(unrolled.DBN.adj.matrix.list, file = paste(getwd(), 'asset/unrolled.DBN.adj.matrix.list.RData', sep = '/'))
rm(input.data.discr.3D, mi.net.adj.matrix)


# Learn the rolled DBN adj matrix
# rolled.DBN.adj.matrix <- rollDbn(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, roll.method, allow.self.loop)
# rolled.DBN.adj.matrix <- rollDbn(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, 'any', FALSE)
rolled.DBN.adj.matrix <- rollDbn_v2(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix.list, 
                       'any', FALSE)
di.net.adj.matrix <- rolled.DBN.adj.matrix
rm(rolled.DBN.adj.matrix)

# writeLines('\n di.net.adj.matrix = \n')
# print(di.net.adj.matrix)
## Change the node names back to the original node names
rownames(di.net.adj.matrix) <- orig.node.names
colnames(di.net.adj.matrix) <- orig.node.names
save(di.net.adj.matrix, file = paste(getwd(), 'asset/di.net.adj.matrix.RData', sep = '/'))

## Create an '.sif' file equivalent to the directed net adjacency matrix
## that is readable in Cytoscape.
adjmxToSif(di.net.adj.matrix)
# rm(unrolled.DBN.adj.matrix)
rm(unrolled.DBN.adj.matrix.list)
#------------------------------------------------------------
# End: Learn LDBN Structures
#------------------------------------------------------------

# writeLines('\n Iteration 1: After DBN = \n')

# inferredNet <- di.net.adj.matrix
rm(di.net.adj.matrix)

# targetNet <- true.net.adj.matrix

## Begin: Create the format for result
# Result <- matrix(0, nrow = 1, ncol = 8)
# colnames(Result) <- list('TPR', 'FPR', 'FDR', 'PPV', 'ACC', 'MCC',  'F', 'AUC')

# Result <- matrix(0, nrow = 1, ncol = 11)
# colnames(Result) <- list('TP', 'TN', 'FP', 'FN', 'TPR', 'FPR', 'FDR', 'PPV', 'ACC', 'MCC',  'F')
# ## End: Create the format for result

# ResultVsTrue <- calcPerfDiNet(inferredNet, targetNet, Result, num.nodes)
# writeLines('Result LDBN vs True = \n')
# print(ResultVsTrue)

elapsed.time <- (proc.time() - start.time) # Stop the timer
writeLines('elapsed.time = \n')
print(elapsed.time)

sink()
# stop('First iteration is over.')

# # Begin: Uncomment this section to include multiple iterations
# #------------------------------------------------------------
# # Begin: Learn CMI Net Structure
# #------------------------------------------------------------
# beta <- 0.03
# 
# cmi.net.adj.matrix <- learnCmiNetStruct(di.net.adj.matrix, input.data, node.names, beta)
# 
# #TEMP--------------cmi.net.adj.matrix.prev <- cmi.net.adj.matrix
# #------------------------------------------------------------
# # End: Learn CMI Net Structure
# #------------------------------------------------------------
# node.nbr.matrix <- mi.net.adj.matrix
# # #------------------------------
# # # Begin: Second Iteration Till Convergence 
# # #------------------------------
#  repeat
#  {
#    comp <- 0
#    comp <- CompareNet(di.net.adj.matrix,cmi.net.adj.matrix,num.node)
#    if(comp == 1)
# 	{
#       break
# 	}
# 	else
# 	{ 
# 		node.nbr.matrix[1:num.nodes,1:num.nodes] <- 0 # will be use for storing the neighbour of the node.	  
# 		for (rowIdx in 1:(num.nodes - 1))
# 			{
# 			for (colIdx in (rowIdx + 1):num.nodes)
# 				{
# 				 if(cmi.net.adj.matrix[rowIdx,colIdx] == 1) 
# 					{
# 						node.nbr.matrix[rowIdx,colIdx] <- 1
# 						node.nbr.matrix[colIdx,rowIdx] <- 1
# 					}
# 				}
# 			}		
# 	#-----------kNN for k = 2--------------*
#     node.2NN.matrix <- node.nbr.matrix
#    
# 	for(i in 1:num.nodes) # for each node 'i'
# 		{
# 		for(j in 1:num.nodes) # for each node 'j'
# 			{
# 			if(node.nbr.matrix[i,j] == 1)
# 				{
# 				for(k in 1:num.nodes) # for each node 'k'
# 					{
# 					# If 'k' is a 2nd nearest neighbour of 'i'
# 					if(node.nbr.matrix[j,k] == 1 && k !=i)
# 						{
# 							node.2NN.matrix[i,k] <- 1;
# 						}
# 					}			
# 				}
# 			}
# 		}
# 	    di.net.adj.matrix <- learnDbnStruct3D(discr.input.data, node.2NN.matrix)
#         cmi.net.adj.matrix <- learnCmiNetStruct(di.net.adj.matrix, input.data, node.names, beta)		
# 	}
#  }
# #------------------------------
# # End: Second Iteration Till Convergence 
# #------------------------------
# 
# plot(as(final_net, "graphNEL"))
# targetNet <- #_______________________________True network to assign by the user.
# Result <- calcperfdinet(cmi.net.adj.matrix, targetNet, Result, node.names) 
# print(Result)
# # End: Uncomment this section to include multiple iterations


#------------------------------------------------------------
# Begin: Save R Session Info in a File
#------------------------------------------------------------
sink("asset/sessionInfo.txt")
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
