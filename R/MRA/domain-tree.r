domain.tree.FSA = function(locs, r){
  domain.tree=list(r=seq(nrow(locs)))
  for( idx in 1:(nrow(locs)-r) ){
    locs.idx = idx + r
    id = paste("r",idx,sep="_")
    domain.tree[[id]] = locs.idx
  }

  D = fields::rdist(locs[-(1:r),], locs[1:r,])
  knot.regions = apply(D, 2, which.min)
  for( knot.idx in seq(r,1)) {
    region = knot.regions[knot.idx]
    region.id = paste("r", region, sep="_")
    domain.tree[[region.id]] = c(knot.idx, domain.tree[[region.id]])

  }
  return(domain.tree)
}



domain.tree.J4 = function( locs, mra.options ){

  r = mra.options[['r']]; J = 4; n = length(locs)/ncol(locs)
  points = seq(n)

  M = floor(log((n/r)*(J-1) + 1)/log(J)) # not sure if -1 should be left in or out
  #M = floor(log((n/r)*(J-1) + 1)/log(J)) - 1
  if( M==0 ) stop(paste(c('ERROR: n=', n, ' points is not enough for J=4 and r=', r,' basis functions'), collapse=""))
  addOnFirstRes = m - r*(1-J^M)/(1-J)
  grid.tree = list(r=points)
  inds = genInds(M)

  for( ind in inds){

    if( child.id(ind)==1 ){

        par.inds = grid.tree[[parent(ind)]]
        par.locs = locs[par.inds,]

        if( ncol(locs)==2 ) {
          x_split = quantile(par.locs[,1], 0.5)
          y_split = quantile(par.locs[,2], 0.5)

          reg1 = par.inds[which(par.locs[,1]<=x_split & par.locs[,2]<=y_split)]
          reg2 = par.inds[which(par.locs[,1]> x_split & par.locs[,2]<=y_split)]
          reg3 = par.inds[which(par.locs[,1]<=x_split & par.locs[,2]> y_split)]
          reg4 = par.inds[which(par.locs[,1]> x_split & par.locs[,2]> y_split)]

          subregs = list(reg1, reg2, reg3, reg4)

      } else if( ncol(locs)==1 ) {

        x_split_1 = quantile(par.locs, 0.25)
        x_split_2 = quantile(par.locs, 0.5)
        x_split_3 = quantile(par.locs, 0.75)

        reg1 = par.inds[which(par.locs<= x_split_1)]
        reg2 = par.inds[which(par.locs>x_split_1 & par.locs<=x_split_2)]
        reg3 = par.inds[which(par.locs>x_split_2 & par.locs<=x_split_3)]
        reg4 = par.inds[which(par.locs>x_split_3)]

        subregs = list(reg1, reg2, reg3, reg4)
      }
    }

    grid.tree[[ind]] = subregs[[child.id(ind)]]
  }
  return(grid.tree)
}



domain.tree.J2 = function( locs, mra.options ){

  n = length(locs)/ncol(locs)
  points = seq(n)

  r = mra.options[['r']]
  J = mra.options[['J']]
  M = mra.options[['M']]
  grid.tree = list(r=points)
  inds = genInds(M,J=c(2))

  for( ind in inds ) {
    if( child.id(ind)==1 ){
      par.inds = grid.tree[[parent(ind)]]
      par.locs = locs[par.inds,]

      if( ncol(locs)==2 ){
        if( (res(ind) %% 2) == 1 ) {
          x_split = quantile(par.locs[,1],0.5)
          reg1 = par.inds[which(par.locs[,1]<=x_split)]
          reg2 = par.inds[which(par.locs[,1]>x_split)]
        } else {
          y_split = quantile(par.locs[,2],0.5)
          reg1 = par.inds[which(par.locs[,2]<=y_split)]
          reg2 = par.inds[which(par.locs[,2]>y_split)]
        }
      } else if(ncol(locs)==1 ) {
        x_split = quantile(par.locs, 0.5)
        reg1 = par.inds[which(par.locs<x_split)]
        reg2 = par.inds[which(par.locs>=x_split)]
      }
      subregs = list(reg1, reg2)
    }

    grid.tree[[ind]] = subregs[[child.id(ind)]]
  }
  return(grid.tree)
}

# domain.tree.FSA = function( locs, mra.options ){
#
#   J = mra.options[['J']]
#   n = length(locs)/ncol(locs)
#   points = seq(n)
#   M = 1
#   grid.tree = list(r=points)
#   centers = locs[1:J,]
#
#   D = fields::rdist(locs[-seq(J),], centers)
#   cents = apply(D, 1, which.min)
#   for( j in 1:J){
#     ind = paste("r", j, sep="")
#     grid.tree[[ind]] = c(j, J + which(cents==j))
#   }
#
#   return(grid.tree)
# }
