library(fields)
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
        indices[[2*i]] = lapply(children, function(num) paste(as.character(i), num, sep="_") )
        indices[[2*i-1]] = list(i)
      }
      return(do.call(c, unlist(indices, recursive=FALSE)))
    }
  }

  inds = genIndices(M,J)
  inds = sapply(inds, function(ind) paste("r", ind, sep="_"))

  lengths = sapply(inds, function(s) nchar(s))
  ord = order(lengths)
  return(inds[ord])
}



## plots ordered locations such that color intensity is decreasing with order
plot.locsord = function(locsord, col = "#000000", col2="#FFFFFF"){

  nlocs = length(locsord)/ncol(locsord)
  vals = seq(nlocs)
  #collist = grey.colors(nlocs, start=0.2, end=0.95)

  colf = function(f){
    color = rgb((1-f)*(t(col2rgb(col))) + f*(t(col2rgb(col2))), maxColorValue = 255)
    return(color)
  }

  collist = sapply(vals/nlocs, colf)
  fields::quilt.plot(locsord, vals, col=collist)
}
