findOrderedHierarchyJ4 = function( locs, r=2 ){
  
  r = 2; J = 4; n = length(locs)/ncol(locs)
  points = seq(n)
  
  M = floor(log((n/r)*(J-1) + 1)/log(J))-1

  addOnFirstRes = m - r*(1-J^M)/(1-J)
  grid.tree = list(r=points)
  inds = genInds(M)
  
  for( ind in inds){
    
    if( child.id(ind)==1 ){
      par.inds = grid.tree[[parent(ind)]]
      par.locs = locs[par.inds,]
      
      x_split = quantile(par.locs[,1], 0.5)
      y_split = quantile(par.locs[,2], 0.5)
      
      reg1 = par.inds[which(par.locs[,1]<=x_split & par.locs[,2]<=y_split)]
      reg2 = par.inds[which(par.locs[,1]> x_split & par.locs[,2]<=y_split)]
      reg3 = par.inds[which(par.locs[,1]<=x_split & par.locs[,2]> y_split)]
      reg4 = par.inds[which(par.locs[,1]> x_split & par.locs[,2]> y_split)]
      
      subregs = list(reg1, reg2, reg3, reg4)
    }

    grid.tree[[ind]] = subregs[[child.id(ind)]]
  }
  return(grid.tree)
}


buildKnotTree = function(locs.tree, r){
  
  M = get.M(locs.tree)
  
  knots = list()
  remaining = list(r=locs.tree[["r"]])
  for(ind in names(locs.tree)){
    node.locs = locs.tree[[ind]]
    available = intersect(node.locs, remaining[[parent(ind)]])
    if( res(ind)==M ) knots[[ind]] = available
    else knots[[ind]] = available[1:min(r, length(available))]
    remaining[[ind]] = setdiff(remaining[[parent(ind)]], knots[[ind]])

  }
  return(knots)
}