source('MRA/tree-methods.r')
source('MRA/knot-tree-methods.r')


knot.tree = function(locs.tree, r, dim=2){

  M = get.M(locs.tree)
  Jm = get.Jm(locs.tree)

  knots = list()
  remaining = list(r=locs.tree[["r"]])
  exact = all(Jm==Jm[1]) && Jm[1]==r+1 && dim==1
  if( exact ) print("exact representation available!")

    for( ind in names(locs.tree) ){

      node.locs = locs.tree[[ind]]
      available = intersect(node.locs, remaining[[parent(ind)]])

      if(res(ind)==M ) knots[[ind]] = available       # if we are at the last resolution, everything is a knot

      else {                                          # if not at the last resolution:
        if( exact ){  #spectial case when we place knots at the split points in 1d
          #browser()
          children = Filter(function(s) startsWith(s, ind) && nchar(s)==nchar(ind)+1, names(locs.tree))
          knts = c()
          for( child in children[-1] ){
            locs.child.available = intersect(locs.tree[[child]], available)

            knts = c(knts, locs.child.available[1])
          }
          knots[[ind]] = knts
        } else {
          #browser()
          knots[[ind]] = available[1:min(r, length(available))]
        }
      }

      remaining[[ind]] = setdiff(remaining[[parent(ind)]], knots[[ind]])
    }
  return(knots)
}
