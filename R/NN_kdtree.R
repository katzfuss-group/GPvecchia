library(FNN)

findOrderedNN_kdtree1 <- function(locs,m,mult=2){

    n <- nrow(locs)
    NNarray <- matrix(NA,n,m+1)
    NNarray[,1] <- 1:n
    NNarray[1:(mult*m+1),] <- findOrderedNN(locs[1:(mult*m+1),],m)

    data_inds <- 1:n
    query_inds <- (mult*m+2):n

    msearch <- mult*m

    while( length(query_inds) > 0 ){

        msearch <- min( max(query_inds), msearch )
        data_inds <- 1:max(query_inds)
        NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
        less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
        sum_less_than_k <- apply(less_than_k,1,sum)
        ind_less_than_k <- which(sum_less_than_k >= m+1)
        NN_less_than_k <- NN[ind_less_than_k,]

        NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))

        NNarray[ query_inds[ind_less_than_k], ] <- NN_m

        query_inds <- query_inds[-ind_less_than_k]
        #print(length(query_inds))
        msearch <- min( length(data_inds), msearch )


    }
    return(NNarray)
}


findOrderedNN_kdtree2 <- function(locs,m,mult=2){

    n <- nrow(locs)
    NNarray <- matrix(NA,n,m+1)
    NNarray[,1] <- 1:n
    maxval <- min(mult * m + 1, n)
    NNarray[1:(maxval),] <- findOrderedNN(locs[1:(maxval),,drop=FALSE],m)

    query_inds <- min(maxval + 1, n):n
    data_inds <- 1:n

    msearch <- m

    while( length(query_inds) > 0 ){

        msearch <- min( max(query_inds), 2*msearch )
        data_inds <- 1:max(query_inds)
        NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
        less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
        sum_less_than_k <- apply(less_than_k,1,sum)
        ind_less_than_k <- which(sum_less_than_k >= m+1)
        NN_less_than_k <- NN[ind_less_than_k,]

        NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))

        NNarray[ query_inds[ind_less_than_k], ] <- NN_m

        query_inds <- query_inds[-ind_less_than_k]
        # print(length(query_inds))

    }

    return(NNarray)
}


# naive nearest neighbor finder

findOrderedNN <- function( locs, m ){
  # find the m+1 nearest neighbors to locs[j,] in locs[1:j,]
  # by convention, this includes locs[j,], which is distance 0
  n <- dim(locs)[1]
  NNarray <- matrix(NA,n,m+1)
  for(j in 1:n ){
    distvec <- c(fields::rdist(locs[1:j,,drop=FALSE],locs[j,,drop=FALSE]) )
    NNarray[j,1:min(m+1,j)] <- order(distvec)[1:min(m+1,j)]
  }
  NNarray
}










