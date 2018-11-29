source('MRA/utility-functions.r')


parent = function(id){
  if( nchar(id) > 1){
    new_id = strsplit(id, "_")[[1]]
    return(new_id[1:length(new_id)-1])
  } else {
    return(id)
  }
}

child.id = function(id){
  new_id = strsplit(id, "_")[[1]]
  num = as.integer(new_id[length(new_id)])
  return(num)
}

res = function(id){
  res = length(strsplit(id, "_")[[1]])-1
  return(res)
}


get.M = function(tree){
  all.knots = names(tree)
  M = max(sapply(all.knots, function(x) length(strsplit(x,split="_")[[1]])))-1
  return(M)
}


get.Jm = function(tree, m=NULL){
  M = get.M(tree)
  all.knots = names(tree)
  lengths = sapply(all.knots, function(x) length(strsplit(x,split="_")[[1]]))
  Js = rep(0, M)
  for( m in seq(2,max(M+1, 2))) {
    inds.m = all.knots[lengths==m]
    maxJ = max(sapply(inds.m, function(ind) as.integer(strsplit(ind, split="_")[[1]][m])))
    Js[m-1] = maxJ
  }
  return(Js)
}
