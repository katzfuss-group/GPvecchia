source('MRA/tree-methods.r')

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
  min.length = min(sapply(names(knot.tree), function(n) nchar(n)))
  roots = which(sapply(names(knot.tree), nchar)==min.length)

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


ord.knot.tree = function(knt.tree){
  ord = c()
  nlocs = max(sapply(knt.tree, function(node) max(node)))
  i=1
  for( node in names(knt.tree) ){
    node.locs = knt.tree[[node]]
    if(length(node.locs)==0) next
    n = min(length(node.locs), nlocs)
    ord = c(ord, node.locs)
    new.ind = seq(i, i+n-1)
    knt.tree[[node]] = new.ind
    i = i + n
  }
  return(list(knot.tree=knt.tree, ord=ord))
}
