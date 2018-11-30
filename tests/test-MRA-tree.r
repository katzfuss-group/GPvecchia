rm(list=ls())

setwd("/home/marcin/GPvecchia")
source("R/ordering_functions.R")
source("MRA/domain-tree.r")
source("MRA/knot-tree.r")
source("MRA/tree-plotting-methods.r")
source("MRA/mraNN.r")

exactCase = FALSE

#exact case
if( exactCase ) {
  spatial.dim = 1
  exact = TRUE
  n=3; m=2; r=1
} else {
  spatial.dim=1 # number of spatial dimensions
  n=8 # number of observed locs
  m=3; r=2
}

# simulate locations
set.seed(10)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
  if( exactCase ) {
    ord = order_coordinate((locs))
    locs = matrix(locs[ord])
  }
} else {
  locs = cbind(runif(n),runif(n))
  ord = order_maxmin(locs)
  locs = locs[ord,]
}

#ind.tree = domain.tree.J2(locs, m)
#knt.tree = knot.tree(ind.tree, r, dim=1)
#plot.locs.tree(ind.tree, locs, knots=knt.tree)
#mra.tree = ord.knot.tree(knt.tree)

#print(mra.tree$knot.tree)
#print(mra.tree$ord)



#Narray = getNNmatrix(mra.tree$knot.tree)
#print(NNarray)
NNarray = findOrderedNN_mra(locs, J=4, r=1)
print(NNarray)
