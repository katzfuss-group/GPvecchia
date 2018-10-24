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


choose.M = function(n, m) {
  M=0
  while(2^(M+1)/(M+1) <= n/m) M=M+1
  
  ## for very small m:
  if(M+1>m) {
    M=m-1
    r=rep(1, M+1)
    J=c(1, rep(2, M-1), ceiling((n-sum(2^(0:(M-1))))/2^(M-1)))
  } else {
    J = c(1, rep(2,M))
    
    ## choose r based on m
    r=rep(ceiling(m/(M+1)),M+1)
    l=M
    while(sum(r)>m) {
      r[l] = r[l] - 1
      l=l-1
    }
  }
  
  if(sum(r)>m | sum(r*cumprod(J))<n) print ('ERROR')
  else return(M, r)
  
}


## generate a hierarchy of indices ordered by resolution and 
## lexicographically within each resolution
genInds = function(M, J=c(4)){
  
  genIndices = function(mLeft, J=c(4)){
    Jm = J[1]
    if(mLeft==1){
      return(sapply(seq(Jm), function(num) as.character(num)))
    } else {
      if( length(J)>1 ) J=J[2:length(J)]
      children = genIndices(mLeft-1, J=J)
      indices = list()
      for(i in 1:Jm){
        indices[[2*i]] = lapply(children, function(num) paste(c(as.character(i), num), collapse="") )
        indices[[2*i-1]] = list(i)
      }
      return(do.call(c, unlist(indices, recursive=FALSE)))
    }
  }
  
  inds = genIndices(M)
  inds = sapply(inds, function(ind) paste(c("r", ind), collapse=""))
  
  lengths = sapply(inds, function(s) nchar(s))
  ord = order(lengths)
  return(inds[ord])
}