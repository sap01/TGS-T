# Test nested foreach loop along with 'when' clause
# $ cd /home/saptarshi/R/R-3.3.2/projects/repoprojldbn/  
# $ nohup /home/saptarshi/R/R-3.3.2/bin/Rscript /home/saptarshi/R/R-3.3.2/projects/repoprojldbn/testnestedforeachwhen.R &

library(foreach)
library(doParallel)

no_cores <- 4
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

ipj <- 
foreach::foreach(i = 1:4, .packages = c('foreach')) %:% 
when(i != 3) %:%
foreach::foreach(j = 21:22, .packages = c('foreach','bnstruct')) %dopar%
{
i+j
}
save(ipj, file = 'ipj.RData')
