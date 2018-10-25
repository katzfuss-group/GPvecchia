source('MRA/tree-methods.r')
source('MRA/utility-functions.r')

# findOrderedHierarchy = function(locs, r=NULL, m=NULL, M=NULL, J=NULL){
#
#   if( is.null(r) && is.null(m) && is.null(M)) print("partitioning the domain requires m, r, J or M")
#
#   if(J==4) {
#
#   }
#
# }


domain.tree.J4 = function( locs, r=2 ){

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



domain.tree.J2 = function( locs, m ){

  n = length(locs)
  points = seq(n)
  M = choose.M( n, m )
  grid.tree = list(r=points)
  inds = genInds(M)
  for( ind in inds ) {
    if( child.id(ind)==1 ){
      par.inds = grid.tree[[parent(ind)]]
      par.locs = locs[par.inds,]

      if( (res %% 2) == 1 ){
        x_split = quantile(par.locs[,1],0.5)
        reg1 = par.inds[which(par.locs[,1]<=x_split)]
        reg2 = par.inds[which(par.locs[,1]>x_split)]
      } else {
        y_split = quantile(par.locs[,2],0.5)
        reg1 = par.inds[which(par.locs[,2]<=y_split)]
        reg2 = par.inds[which(par.locs[,2]>y_split)]
      }

      subregs = list(reg1, reg2)
    }

    grid.tree[[ind]] = subregs[[child.id(ind)]]
  }
}

# findOrderedHierarchyFSA = function( locs, J ){
#
#   n = length(locs)
#   points = seq(n)
#   M = 1
#   grid.tree = list(r=points)
#   inds =
#
#
# }
