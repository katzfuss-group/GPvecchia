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



get.mra.params2 = function(n,m, opts){

  params = list(m=opts$m)
  if(is.null(opts$J)) {
    J=2
    warning("J not specified. Setting J=2")
  } else J=opts$J
  if(is.null(opts$M) ){
    if(is.null(opts$r)) {
      M = 0
      while(J^(M+1)/(M+1) <= n/m ) M=M+1
      if(M+1>m){
        M = m-1
        r = rep(1,M+1)
        J = c(1, rep(2, M-1), ceiling((n-sum(2^(0:M-1)))/2^(M-1)))
      } else {
        J = c(1, rep(2,M))
        r = rep(ceiling(m/(M+1)), M+1)
        l = M
        while(sum(r)>m) {
          r[l]=r[l]-1
          l=l-1
        }
      }
    } else {
      r = opts$r
      if(length(r)>1) M=length(r)-1
      else if(length(J)>1) M=length(J)
      else M = floor((log(n/r)*(J-1)+1)/log(J))-1
    }
  } else if(is.null(opts$r)) {
    M = opts$M
    r = floor(m/M)
  } else {
    warning("M, r set for MRA. Parameter m will be overridden")
    M = opts$M; r = opts$r
  }
  params[['J']] = J; params[['M']] = M; params[['r']] = r

  return(params)
}



## figures out the parameters to pass to the function creating the NNmatrix with mra conditioning
## based on the mra.options list provided by the user
get.mra.params = function(n, m, opts){
  params = list(r=opts$r,M=opts$M, m=m, J=opts$J)

  if(!is.null(m) && 'J' %in% names(opts) && 'r' %in% names(opts) && 'M' %in% names(opts)) {

    if(length(opts$r)>opts$M) stop("The r vector has to have length at most M-1")
  }

  if(is.null(opts$J)) params$J=2

  if(is.null(opts$M) && 'J' %in% names(opts) && 'r' %in% names(opts) && length(opts$r)==1){
    M = floor(log((n/r)*(J-1) + 1)/log(J))
  } else if (!is.null(opts[['M']])){
    if(is.null(opts$r) && opts[['M']]==1 ){
      params$J=floor(2*n/m - 1)
      params$r=m %/% 2
    } else if(opts$M==1 && length(opts$r)==2 && opts$r[2]==0){
      params$J=n-opts$r[1]
    }
  }

  return(params)
}
