#source("GPvecchia/R/domain-tree.r")
#source("GPvecchia/R/knot-tree.r")
#source("GPvecchia/R/tree-methods.r")
#source("utility-functions.r")


choose.M = function(n, m) {

  M=1
  while(2^(M+1)/(M+1) <= n/m) M=M+1

  ## for very small m:
  if(M+1>m) {
    M=m-1
    r=rep(1,M+1)
    last.J = 2**ceiling(base::log((n-sum(2^(0:(M-1))))/2^(M-1),2))
    J=c(rep(2,max(M-1,0)),last.J)
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
  #if(sum(r)>m | sum(r*cumprod(J))<n) print('ERROR')
  if(sum(r)>m | sum(r*cumprod(c(1,J)))<n) print('ERROR')
  else return(list(M=M, r=r, J=J))
}



get.mra.params = function(n,opts,m){

  params = list(m=m)

  # plotting: yes/no
  if(is.null(opts$plots) || opts$plots==FALSE) params$plots=FALSE
  else params$plots=TRUE

  # set J first
  if(is.null(opts$J) && is.null(opts$M)) {
    if(opts$r[1]==0 && length(opts$r)==2) J = 2**ceiling(base::log(n/opts$r[2],2)) #needed for independent block
    else if(length(opts$r)==2 && opts$r[2]==1) J = 2**ceiling(base::log(n-opts$r[1],2)) #needed for low-rank
    else J=2
    warning("J not specified. Setting J=",J)
  } else J=opts$J

  # set M and r
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
    r = ceiling(m/(M+1))
    if(is.null(opts$J)){
      last.J = 2**(ceiling(base::log((n - r*(2^M-1))/r, 2)) - (M-1))
      J = c(rep(2, M-1), last.J)
    }
  } else {
    if( m>0 ) warning("M, r set for MRA. If parameter m was given, it will be overridden")
    M = opts$M; r = opts$r
  }
  params[['J']] = J; params[['M']] = M; params[['r']] = r
  return(params)
}





findOrderedNN_mra = function(locs, mra.options, m=-1){

  # make sure the low-rank case is covered as well

  n = length(locs)/ncol(locs)
  mra.params = get.mra.params(n, mra.options, m)
  print(mra.params)

  if(mra.params$J!=2 && mra.params$J!=4){
    if( mra.params$M>1) warning("When J is neither 2 nor 4 we always set M to 1 and use the Full scale approximation")
    if(length(mra.params$r)==2 && mra.params$r[2]==1) ind.tree = domain.tree.low.rank(locs, mra.params)
    else ind.tree = domain.tree.FSA(locs, mra.params)
  #} else if( all(mra.params[['J']]==2) ) ind.tree = domain.tree.J2(locs, mra.params)
  #else if( all(mra.params[['J']]==4) ) ind.tree = domain.tree.J4(locs, mra.params)
  #else stop(paste("J of the form c(", paste(mra.params$J, collapse=","), ") not supported. Try using a scalar J.", sep=""))
  } else if( any(mra.params[['J']]==2) ) ind.tree = domain.tree.J2(locs, mra.params)
  else if( all(mra.params[['J']]==4) ) ind.tree = domain.tree.J4(locs, mra.params)
  else stop(paste("J of the form c(", paste(mra.params$J, collapse=","), ") not supported. Try using a scalar J.", sep=""))

  knt.tree = knot.tree(ind.tree, mra.params[['r']], dim=ncol(locs))
  if(mra.params$plots==TRUE)  plot.locs.tree(ind.tree, locs, knots=knt.tree)
  mat = getNNmatrix(knt.tree)
  eff.m = ncol(mat)-1
  #print(paste("effective m: ", eff.m-1, sep=""))
  if(eff.m > 100) print(paste("Effective m is ", ncol(mat)-1, " which might slow down computations", sep=""))

  print(paste("MRA params: m=",eff.m, ", J=", paste(get.Jm(ind.tree), collapse=","), ", r=", paste(get.rm(knt.tree), collapse=","), ", M=", get.M(ind.tree), sep=""))

  return(mat)
}
