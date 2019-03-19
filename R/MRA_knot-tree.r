library(dequer)


knot.tree.old = function(locs.tree, r, dim=2){

  M = get.M(locs.tree)
  Jm = get.Jm(locs.tree)
  knots = list()
  remaining = list(r=locs.tree[["r"]])
  exact = all(Jm==Jm[1]) && length(r)==1 && Jm[1]==r+1 && dim==1

  if( exact ) print("exact representation available!")

    for( ind in names(locs.tree) ){
      node.locs = locs.tree[[ind]]
      available = intersect(node.locs, remaining[[parent(ind)]])

      if(res(ind)==M ) knots[[ind]] = available       # if we are at the last resolution, everything is a knot
      else {                                          # if not at the last resolution:
        if( exact ){  #spectial case when we place knots at the split points in 1d
          children = Filter(function(s) startsWith(s, ind) && nchar(s)==nchar(ind)+1, names(locs.tree))
          knts = c()
          for( child in children[-1] ){
            locs.child.available = intersect(locs.tree[[child]], available)
            knts = c(knts, locs.child.available[1])
          }
          knots[[ind]] = knts
        } else {
          no_knots = getNKnt(r, ind)
          if(no_knots)  knots[[ind]] = available[1:min(no_knots, length(available))]
        }
      }

      remaining[[ind]] = setdiff(remaining[[parent(ind)]], knots[[ind]])
    }
  return(knots)
}






knot.tree = function(locs, mra.params){


  knots = new.env(hash=TRUE)
  indices = c()

  #knots = list()
  J = mra.params$J; M = mra.params$M; r = mra.params$r


  remaining = list()
  remaining[["r"]] = seq(n)


  #remaining = queue()
  #pushback( remaining, list("r", seq(n)) )

  while(length(remaining)>0){
    id = names(remaining)[1]
    #item = pop(remaining)
    #id = item[[1]]
    m = res(id)
    reg.inds = remaining[[id]]
    #reg.inds = item[[2]]

    indices = c(indices, id)

    if(m<M){
      r.eff = min(r[m+1], length(reg.inds))
      if(r.eff>0) {
        knots[[id]] = reg.inds[seq(r.eff)]
        reg.inds = reg.inds[-seq(r.eff)]
      }
      reg.locs = locs[reg.inds,]
      if(!is.matrix(reg.locs)) reg.locs = matrix(reg.locs, ncol=ncol(locs))

      if( length(reg.locs)==0 ) clusters = c()
      else{
        if( J[m+1]>nrow(reg.locs) ){
        clusters = seq(length(reg.locs))
        } else {
          clusters =  cluster.equal(reg.locs, K=J[m+1], dim.start=m%%2+1)
        }
      }
      for(child.no in 1:J[m+1]){
        child.id = paste(c(id, child.no), collapse="_")
        remaining[[child.id]] = reg.inds[clusters==child.no]
        #pushback(remaining, list(child.id, reg.inds[clusters==child.no]))
      }
    } else {
      knots[[id]] = reg.inds
    }
    remaining = remaining[-1]
  }
  knots$indices = indices
  return(knots)
}





getNKnt = function(r, ind){
  if( length(r)==1 ) {
    r
  } else {
    m = res(ind)
    if(m+1>length(r)){
      r[length(r)]
    } else
      r[m+1]
  }
}


# takes the list of knot indices and sorts them in the breadth-first manner
sortNodesBFS = function(knots){



}






# getNNmatrix = function(knot.tree,m){
getNNmatrix = function(knot.tree){

  indices = knot.tree$indices
  rm(indices, envir=knot.tree)

  neighbors = list()
  #fill out the list of neighbors for the root

  #find the root(s) of the tree
  min.length = min(sapply(indices, res))
  roots = which(sapply(indices, res)==min.length)


  for( root.no in roots ){
    root = indices[root.no]
    if( length(knot.tree[[root]])==0) next()
    root.ind = knot.tree[[root]][1]
    neighbors[[root.ind]] = c(root.ind); cond.set=c(root.ind); last.knot=root.ind
    for( knot in knot.tree[[root]][-1]){
      cond.set = c(knot, cond.set)
      neighbors[[knot]] = cond.set
      last.knot = knot
    }
  }


  # once the first knot is handled, fill out the list
  # for the remaining knots
  for( ind in indices[-roots] ){
    knots = knot.tree[[ind]]
    if(all(is.na(knots))) next()
    parent.knots = knot.tree[[parent(ind)]]
    last.knot = parent.knots[length(parent.knots)]
    cond.set = neighbors[[last.knot]]
    for( knot in knots) {
      cond.set = c(knot, cond.set)
      neighbors[[knot]] = cond.set
    }
  }

  list2matrix(neighbors)
  NNarray = list2matrix(neighbors)

  # #sometimes effective m is less than m. If so, we need to add columns
  # #to ensure NNarray has m+1 columns
  # if(ncol(NNarray)<m){
  #   Ncols = ncol(NNarray)
  #   padding = matrix(rep(NA, (m+1-Ncols)*nrow(NNarray)), ncol=m+1-Ncols)
  #   NNarray = cbind(NNarray, padding)
  # }

  return(NNarray)
}
