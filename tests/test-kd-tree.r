source("R/NN_kdtree.R")

spatial.dim=2 # number of spatial dimensions
n=43  # number of observed locs
m=4

# simulate locations
set.seed(10)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
}
ord = order_maxmin(locs)

locs = locs[ord,]

NNarray = findOrderedNN_kdtree2(locs, 4)
