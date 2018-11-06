source('MRA/utility-functions.r')
source('MRA/tree-methods.r')


choose.M = function(n, m) {

  M=0
  while(2^(M+1)/(M+1) <= n/m) M=M+1

  ## for very small m:
  if(M+1>m) {

    M=m-1
    r=rep(1,M+1)
    J=c(1,rep(2,M-1),ceiling((n-sum(2^(0:(M-1))))/2^(M-1)))

  } else{

    J=c(1,rep(2,M))

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




  else return(list(M=M, r=r))

}
