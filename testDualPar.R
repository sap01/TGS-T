# Goal: Test Level 2 Nested Parallelism

#---------------------------------
# Begin: Loading the Packages
#---------------------------------
library(foreach)
library(parallel)
library(doSNOW)
#---------------------------------
# End: Loading the Packages
#---------------------------------

#---------------------------------
# Begin: Step 2: Decomposing the network
#---------------------------------
num.nodes <- 10
num.timepts <- 21

num.clnodes <- num.nodes
num.cores.per.clnode <- num.timepts

# Type of cluster to use: FORK or SOCK? 
# Ref: https://www.r-bloggers.com/how-to-go-parallel-in-r-basics-tips/
# Since 'grni' is used to test this program, which is a single Ubuntu Linux node
# with 24 logical CPUs, 'SOCK' is only able to create (1 + 10) = 11 workers.
# Whereas, 'FORK' is able to create 22 workers.
cl <- parallel::makeCluster(num.clnodes, type = "FORK")
parallel::clusterCall(cl, function() {
  library(doSNOW)
  doSNOW::registerDoSNOW(parallel::makeCluster(num.cores.per.clnode, type = "FORK"))
  NULL
})
registerDoSNOW(cl)

writeLines('\n nodename \n')
parallel::clusterEvalQ(cl, Sys.info()['nodename'])

outer.matrix <- foreach::foreach(node.idx = 1:num.nodes, .combine = 'cbind')  %:% 
  foreach::foreach(time.idx = 1:num.timepts, .combine = 'c')  %dopar% 
  {
    # Sys.sleep(2)
    time.idx + num.timepts * (node.idx -1)
  }

# outer.matrix <- foreach::foreach(node.idx = 1:num.nodes, .combine = 'cbind')  %dopar% 
# {
#   inner.matrix <- foreach::foreach(time.idx = 1:num.timepts, .combine = 'c')  %dopar% 
#   {
#     Sys.sleep(60)
#     time.idx + num.timepts * (node.idx -1)
#   }
#   
#   inner.matrix
# }

parallel::stopCluster(cl)

print(outer.matrix)
