source('MRA/utility-functions.r')


parent = function(id){
  if( nchar(id) > 1){
    return(substr(id, 1, nchar(id)-1))
  } else {
    return(id)
  }
}

child.id = function(id){
  return(as.integer(substr(id, nchar(id), nchar(id))))
}

res = function(id){
  res = nchar(id)-1
  return(res)
}


get.M = function(tree){
  all.knots = names(tree)
  M = max(sapply(all.knots, function(x) nchar(x)))-1
  return(M)
}


get.Jm = function(tree, m=NULL){
  M = get.M(tree)
  all.knots = names(tree)
  lengths = sapply(all.knots, function(x) nchar(x))
  Js = rep(0, M)
  for( m in seq(2,max(M+1, 2))) {
    inds.m = all.knots[lengths==m]
    maxJ = max(sapply(inds.m, function(ind) as.integer(substr(ind, nchar(ind), nchar(ind) ))))
    Js[m-1] = maxJ
  }
  return(Js)
}
