#setwd("/home/marcin/GPvecchia")
#source("R/ordering_functions.R")
source("MRA/domain-tree.r")
source("MRA/knot-tree.r")
source("MRA/tree-plotting-methods.r")


findOrderedNN_mra = function(locs, J=4, r=2){

  ind.tree = domain.tree.J4(locs)
  knot.tree = knot.tree(ind.tree, 2)
  #plot.locs.tree(ind.tree, locs, knots=knot.tree)

  mat = getNNmatrix(knot.tree)

  return(mat)
}
