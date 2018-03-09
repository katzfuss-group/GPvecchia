
whichCondOnLatent <- function( NNarray, firstind.pred=nrow(NNarray)+1 ){
  
  # determine for which "neighbors" to condition on y (-> TRUE) vs z (-> FALSE)
  # current strategy: find "best" neighbor, 
  #      then only condition on y for that neighbor or its neighbors
  n <- nrow(NNarray)
  m <- ncol(NNarray)-1
  CondOnLatent <- matrix(NA,n,m+1)
  CondOnLatent[1,1] <- TRUE
  
  for(k in 2:n){
    latents=rep(0,m)
    for(ind in 2:(m+1)){
      l=NNarray[k,ind]
      if(!is.na(l) && l<firstind.pred) 
        latents[ind]=sum(is.element(NNarray[k,],NNarray[l,]*CondOnLatent[l,]))
    }
    ind=NNarray[k,(which(latents==max(latents)))[1]] # could also try max(which(latents==max(latents)))
    CondOnLatent[k, ] <- is.element( NNarray[k,], NNarray[ind,]*CondOnLatent[ind,] )
    CondOnLatent[k,NNarray[k,]>=firstind.pred] <- TRUE
    CondOnLatent[k,1] <- TRUE
    CondOnLatent[k,is.na(NNarray[k,])] <- NA
  }
  return(CondOnLatent)
  
}