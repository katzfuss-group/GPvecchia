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
