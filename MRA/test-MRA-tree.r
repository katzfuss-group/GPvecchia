rm(list=ls())

source("~/GPvecchia/R/ordering_functions.R")
source("~/GPvecchia/tests/hierarchy.r")
source("~/GPvecchia/tests/hierarchy-plotting.r")

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

ind.tree = findOrderedHierarchyJ4(locs)
knot.tree = buildKnotTree(ind.tree, 3)
plot.locs.tree(ind.tree, locs, knots=knot.tree)