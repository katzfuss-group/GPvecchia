## this function determines the sparsity structure of U
## similarly to a symbolic cholesky decomposition, it does not depend on parameter values
## hence, it only has to be carried out once for repeated evaluation of the likelihood

U_sparsity <- function( locs, NNarray, obs, Cond ){
  # locs are the locations
  # row i of NNarray gives the neighbors of location i
  # i-th element of observed is TRUE if loc i is observed

  nnp <- nrow(locs) # total number of locs (incl unobserved)
  n=sum(obs) # number of obs locs
  size=nnp+n # size (i.e., number of rows and columns in U)
  
  # Calculate how many nonzero entries we need to compute
  nentries <- sum( !is.na(NNarray) )
  
  # latent_map specifies which row in U corresponds to which latent variable
  # observed_map specifies which row in U corresponds to which obs variable
  cur <- 1
  latent_map <- rep(NA,n)
  observed_map <- rep(NA,n)
  for(k in 1:nnp){
    latent_map[k] <- cur
    cur <- cur+1
    if( obs[k] ){
      observed_map[k] <- cur
      cur <- cur+1
    }
  }  
  
  # pre-reorder everything
  revNNarray = t(apply(NNarray, 1,rev)) # column-reversed NNarray
  revCondOnLatent = t(apply(Cond, 1,rev)) # column-reversed CondOnLatent T/F
  
  # calculate the indexing of latent variables in U
  rowpointers = colindices = rep(NA,nentries)
  rowpointers[1] = colindices[1] = 1
  cur <- 0
  for( k in 1:nnp){
    inds    <- revNNarray[k,]
    inds0   <- inds[!is.na(inds)]
    n0 <- length(inds0)
    revCond <- revCondOnLatent[k,!is.na(inds)]
    
    cur_row <- latent_map[k]
    cur_cols <- rep(NA,n0)
    cur_cols[revCond]  <- latent_map[inds0[revCond]]
    cur_cols[!revCond] <- observed_map[inds0[!revCond]]
    
    # store pointers to appropriate rows and columns
    if( k > 1 ){
      rowpointers[cur + (1:n0)] <- as.integer(rep(cur_row,n0))
      colindices[cur + (1:n0)] <- as.integer(cur_cols)
    }
    cur=cur+n0
  }
  
  # separately, calculate indexing of obs variables in U
  Zrowpointers = Zcolindices = rep(NA,2*n)
  cur=0
  for(k in 1:nnp){
    if( obs[k] ){
      cur_row <- observed_map[k]
      cur_cols <- c( latent_map[k], observed_map[k] )
      Zrowpointers[cur + (1:2)] <- as.integer(rep(cur_row,2))
      Zcolindices[cur + (1:2)] <- as.integer(cur_cols)
      cur <- cur + 2
    }
  }         
  
  # combine pointers for latent and obs variables
  allrowpointers=c(rowpointers, Zrowpointers)
  allcolindices=c(colindices, Zcolindices)

  # number of cores to be used for creating U later
  n.cores=detectCores(all.tests = FALSE, logical = TRUE)
  
  return(list(revNNarray=revNNarray,revCond=revCondOnLatent,n.cores=n.cores,
        size=size,rowpointers=allrowpointers,colindices=allcolindices,y.ind=latent_map,observed_map=observed_map))
  
}

