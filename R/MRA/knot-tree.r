knot.tree = function(locs.tree, r, dim=2){

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



getNNmatrix = function(knot.tree){

  neighbors = list()
  #fill out the list of neighbors for the root

  #find the root(s) of the tree
  min.length = min(sapply(names(knot.tree), res))
  roots = which(sapply(names(knot.tree), res)==min.length)

  for( root.no in roots ){
    root = names(knot.tree)[root.no]
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
  for( ind in names(knot.tree)[-roots] ){
    knots = knot.tree[[ind]]
    parent.knots = knot.tree[[parent(ind)]]
    last.knot = parent.knots[length(parent.knots)]
    cond.set = neighbors[[parent.knots[length(parent.knots)]]]
    for( knot in knots) {
      cond.set = c(knot, cond.set)
      neighbors[[knot]] = cond.set
    }
  }
  list2matrix(neighbors)
  NNarray = list2matrix(neighbors)
  return(NNarray)
}
