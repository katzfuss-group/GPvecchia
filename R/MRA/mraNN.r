source("R/MRA/domain-tree.r")
source("R/MRA/knot-tree.r")
source("R/MRA/tree-methods.r")
source("R/MRA/utility-functions.r")


choose.M = function(n, m) {

  M=1
  while(2^(M+1)/(M+1) <= n/m) M=M+1

  ## for very small m:
  if(M+1>m) {
    M=m-1
    r=rep(1,M+1)
    J=c(rep(2,M-1),ceiling((n-sum(2^(0:(M-1))))/2^(M-1)))
  } else{
    J=c(rep(2,M))

    ## choose r based on m
    r=rep(ceiling(m/(M+1)),M+1)
    l=M
    while(sum(r)>m) {
      r[l]=r[l]-1
      l=l-1
    }
  }

  ### check that choices are valid
  if(sum(r)>m | sum(r*cumprod(J))<n) print('ERROR')
  else return(list(M=M, r=r, J=J))
}


get.mra.params = function(n,m, opts){

  params = list(m=opts$m)
  if(is.null(opts$J)) {
    J=2
    warning("J not specified. Setting J=2")
  } else J=opts$J
  if(is.null(opts$M) ){
    if(is.null(opts$r)) {
      pars = choose.M(n,m)
      r = pars$r; M = pars$M; J = pars$J
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





findOrderedNN_mra = function(locs, m, mra.options){

  n = length(locs)/ncol(locs)
  mra.params = get.mra.params(n, m, mra.options)

  if((length(mra.params$r)==2 && mra.params$r[2]==0) || (mra.params$J!=2 && mra.params$J!=4)){
    if( mra.params$M>1) warning("When J is neither 2 nor 4 we always set M to 1 and use the Full scale approximation")
    ind.tree = domain.tree.FSA(locs, mra.params$r[1])
  } else if( mra.params[['J']]==2 ) ind.tree = domain.tree.J2(locs, mra.params)
  else if( mra.params[['J']]==4 ) ind.tree = domain.tree.J4(locs, mra.params)

  knt.tree = knot.tree(ind.tree, mra.params[['r']], dim=ncol(locs))
  plot.locs.tree(ind.tree, locs, knots=knt.tree)
  mat = getNNmatrix(knt.tree)
  print(paste("effective m is ", ncol(mat)-1, sep=""))

  return(mat)
}
