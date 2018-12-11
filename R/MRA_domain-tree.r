get.inds.from.children = function(ind, tree){
  inds = c()
  for(node in names(tree)){
    if(parent(node)==ind){
      inds=c(inds,tree[[node]])
    }
  }
  return(inds)
}


domain.tree = function( locs, mra.options ){

  n = length(locs)/ncol(locs)
  points = seq(n)

  r = mra.options[['r']]
  J = mra.options[['J']]
  M = mra.options[['M']]
  grid.tree = list(r=points)
  inds = genInds(M,J=J)
  for( ind in inds ) {
    if( child.id(ind)==1 ){
      par.inds = grid.tree[[parent(ind)]]
      par.locs = locs[par.inds,]
      if(ncol(locs)==1) par.locs = matrix(par.locs, ncol=1)
      j.ind = if( length(J)>=res(ind) ) res(ind) else 1
      clusters = cluster.equal(par.locs, K=J[j.ind], dim.start=res(ind)%%2+1)
      subregs = lapply(seq(J[j.ind]), function(i) par.inds[which(clusters==i)])
    }
    grid.tree[[ind]] = subregs[[child.id(ind)]]
  }

  return(grid.tree)
}
