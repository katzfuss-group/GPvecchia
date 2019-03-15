Rcpp::sourceCpp('~/GPvecchia/src/U_NZentries.cpp')

ichol = function(Sigma, S){

  n = dim(Sigma)[1]
  L = matrix(rep(0, n**2), ncol=n)
  for(i in 1:n){
    j=1
    while(j < i){
      if(S[i,j]!=0){
        L[i,j] = (Sigma[i,j] - sum(L[i,1:j]*L[j,1:j]))/L[j,j]
      } else {
        print(paste("Skipping entry (", i, ",",j,")", sep=""))
      }
      j = j+1
    }
    L[i,i] = sqrt(Sigma[i,i] - sum(L[i,1:i]**2))
  }

  return(L)

}




ichol.compare = function(Sigma, S){
  n = dim(Sigma)[1]
  Li = matrix(rep(0, n**2), ncol=n)
  L = matrix(rep(0, n**2), ncol=n)
  for(i in 1:n){
    j=1
    while(j < i){
      if(S[i,j]!=0) {
        Li[i,j] = (Sigma[i,j] - sum(Li[i,1:j]*Li[j,1:j]))/Li[j,j]
      } else {
        print(paste("Skipping entry (", i, ",",j,")", sep=""))
      }

      L[i,j] = (Sigma[i,j] - sum(L[i,1:j]*L[j,1:j]))/L[j,j]

      if(abs(L[i,j]-Li[i,j])>1e-6 & S[i,j]!=0) print("discrepancy here!")

      j = j+1
    }
    Li[i,i] = sqrt(Sigma[i,i] - sum(Li[i,1:i]**2))
    L[i,i] = sqrt(Sigma[i,i] - sum(L[i,1:i]**2))
    if(abs(L[i,i]-Li[i,i])>1e-6 & S[i,j]!=0) print("discrepancy here!")
  }

  return(list(Li, L))
}
