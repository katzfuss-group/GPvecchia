source('MRA/tree-methods.r')

getNNmatrix = function(knot.tree){

  neighbors = list()

  #fill out the list of neighbors for the root
  neighbors[[1]] = c(1); cond.set=c(1); last.knot=1

  for( knot in knot.tree[["r"]][-1]){
    cond.set = c(knot, cond.set)
    neighbors[[knot]] = cond.set
    last.knot = knot
  }

  # once the first knot is handled, fill out the list
  # for the remaining knots
  for( ind in names(knot.tree)[-1] ){
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
