source("MRA/domain-tree.r")
source("MRA/knot-tree.r")

mra.tree = function(locs, J, m, plot=TRUE){
  
  ord_aux = order_coordinate(locs)
  reverse_map = seq(length(locs))[ord_aux]
  locs_aux = matrix(locs[ord_aux,], ncol=ncol(locs))
  
  if( J==2 ) ind.tree = domain.tree.J2(locs_aux, m)
  else if( J==4 ) ind.tree = domain.tree.J4(locs_aux, m)
  else stop(paste(c("ERROR: partitioning domain into J=", J, " subregions not supported yet"), collapse=""))
  
  knt.tree = knot.tree(ind.tree, 1, dim=ncol(locs))
  if( plot ) plot.locs.tree(ind.tree, locs_aux, knots=knt.tree)

  mra = ord.knot.tree(knt.tree)
  ord = reverse_map[mra$ord]
  
  return(list(tree=mra$knot.tree, ord=ord))
}