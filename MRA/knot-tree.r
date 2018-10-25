source('MRA/tree-methods.r')
source('MRA/knot-tree-methods.r')


knot.tree = function(locs.tree, r){

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
