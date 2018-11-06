#setwd("/home/marcin/GPvecchia")
#source("R/ordering_functions.R")
source("MRA/domain-tree.r")
source("MRA/knot-tree.r")
source("MRA/tree-plotting-methods.r")


findOrderedNN_mra = function(locs, J, r=2){

  if( J==2 ) ind.tree = domain.tree.J2(locs, 3)
  else if( J==4 ) ind.tree = domain.tree.J4(locs, 4)
  else stop(paste(c("ERROR: partitioning domain into J=", J, " subregions not supported yet"), collapse=""))

  knt.tree = knot.tree(ind.tree, 1, dim=ncol(locs))
  #plot.locs.tree(ind.tree, locs, knots=knt.tree)

  mra.tree = ord.knot.tree(knot.tree)
  mat = getNNmatrix(mra.tree$knot.tree)

  return(mat)
}
