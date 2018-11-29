#setwd("/home/marcin/GPvecchia")
#source("R/ordering_functions.R")
source("MRA/domain-tree.r")
source("MRA/knot-tree.r")
source("MRA/tree-plotting-methods.r")


findOrderedNN_mra = function(locs, mra.params){

  n = length(locs)/ncol(locs)

  if( mra.params[['M']]==0 ) ind.tree = domain.tree.indep(locs, mra.params)
  else {
    if((length(mra.params$r)==2 && mra.params$r[2]==0) || (mra.params$J!=2 && mra.params$J!=4)){
      if( mra.params$M>1) warning("When J is neither 2 nor 4 we always set M to 1 and use the Full scale approximation")
      ind.tree = domain.tree.FSA(locs, mra.params$r[1])
    } else if( mra.params[['J']]==2 ) ind.tree = domain.tree.J2(locs, mra.params)
    else if( mra.params[['J']]==4 ) ind.tree = domain.tree.J4(locs, mra.params)
  }

  knt.tree = knot.tree(ind.tree, mra.params[['r']], dim=ncol(locs))
  plot.locs.tree(ind.tree, locs, knots=knt.tree)

  mat = getNNmatrix(knt.tree)

  return(mat)
}
