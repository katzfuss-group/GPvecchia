rm(list=ls())

setwd("/home/marcin/GPvecchia")
source("R/ordering_functions.R")
source("MRA/domain-tree.r")
source("MRA/knot-tree.r")
source("MRA/tree-plotting-methods.r")


spatial.dim=2 # number of spatial dimensions
n=20  # number of observed locs
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




ind.tree = domain.tree.J2(locs, m)
knot.tree = knot.tree(ind.tree, 2)
plot.locs.tree(ind.tree, locs, knots=knot.tree)

getNNmatrix(knot.tree)
