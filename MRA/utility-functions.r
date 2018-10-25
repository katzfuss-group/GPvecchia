## Converts a list of numeric vectors to an n x m matrix,
## where n is the length of the list and m is the length
## of the longest vector. Rows wth shorter vectors are
## padded with `padding`.
list2matrix = function(L, padding=NA){
  a.width = max(sapply(L, function(x) length(x)))
  a.height = length(L)
  A = matrix(rep(padding, a.width*a.height), ncol=a.width)

  for( ind in 1:length(L) ){
    elem = L[[ind]]
    A[ind,1:length(elem)] = elem
  }
  print(A)
  return(A)
}


## generate a list containg a hierarchy of indices ordered
## by resolution and lexicographically within each resolution
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
