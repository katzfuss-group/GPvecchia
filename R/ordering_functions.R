#' Distance to specified point ordering
#'
#' Return the ordering of locations increasing in their
#' distance to some specified location
#'
#' @param locs A matrix of locations. Each row of \code{locs} contains a location, which can
#' be a point in Euclidean space R^d, a point in space-time R^d x T,
#' a longitude and latitude (in degrees) giving a point on the sphere,
#' or a longitude, latitude, and time giving a point in the sphere-time domain.
#' @param loc0 A vector containing a single location in R^d.
#' @param lonlat TRUE/FALSE whether locations are longitudes and latitudes.
#' @return A vector of indices giving the ordering, i.e.
#' the first element of this vector is the index of the location nearest to \code{loc0}.
#' @examples
#' n <- 100             # Number of locations
#' d <- 2               # dimension of domain
#' locs <- matrix( runif(n*d), n, d )
#' loc0 <- c(1/2,1/2)
#' ord <- order_dist_to_point(locs,loc0)
#' @export
order_dist_to_point <- function( locs, loc0, lonlat = FALSE ){

    if(lonlat){
        lon <- locs[,1]
        lat <- locs[,2]
        lonrad <- lon*2*pi/360
        latrad <- (lat+90)*2*pi/360
        x <- sin(latrad)*cos(lonrad)
        y <- sin(latrad)*sin(lonrad)
        z <- cos(latrad)
        locs <- cbind(x,y,z)
    }


    if (!requireNamespace("fields", quietly = TRUE)) {
        stop("fields package required for this function. Please install it.",
             call. = FALSE)
    }
    d <- ncol(locs)         # input dimension
    if(d != length(loc0)){
        stop("location in loc0 not in the same domain as the locations in locs")
    }
    loc0 <- matrix(c(loc0),1,d)             # order by distance to loc0
    distvec <- fields::rdist(locs,loc0)
    orderinds <- order(distvec)
    return(orderinds)
}

#' Middle-out ordering
#'
#' Return the ordering of locations increasing in their
#' distance to the average location
#'
#' @inheritParams order_dist_to_point
#'
#' @return A vector of indices giving the ordering, i.e.
#' the first element of this vector is the index of the location nearest the center.
#' @examples
#' n <- 100             # Number of locations
#' d <- 2               # dimension of domain
#' locs <- matrix( runif(n*d), n, d )
#' ord <- order_middleout(locs)
#' @export
order_middleout <- function( locs, lonlat = FALSE ){

    if(lonlat){
        lon <- locs[,1]
        lat <- locs[,2]
        lonrad <- lon*2*pi/360
        latrad <- (lat+90)*2*pi/360
        x <- sin(latrad)*cos(lonrad)
        y <- sin(latrad)*sin(lonrad)
        z <- cos(latrad)
        locs <- cbind(x,y,z)
    }

    d <- ncol(locs)
    loc0 <- matrix(colMeans(locs),1,d)
    orderinds <- order_dist_to_point(locs,loc0)
    return(orderinds)
}

#' Outside-in ordering
#'
#' Return the ordering of locations decreasing in their
#' distance to the average location.  Reverses middleout.
#'
#' @inheritParams order_dist_to_point
#'
#' @return A vector of indices giving the ordering, i.e.
#' the first element of this vector is the index of the location farthest from the center.
#' @examples
#' n <- 100             # Number of locations
#' d <- 2               # dimension of domain
#' locs <- matrix( runif(n*d), n, d )
#' ord <- order_outsidein(locs)
#' @export
order_outsidein <- function( locs, lonlat = FALSE ){
  orderinds_rev <- order_middleout(locs, lonlat)
  orderinds = rev(orderinds_rev)
  return(orderinds)
}




#' Sorted coordinate ordering
#'
#' Return the ordering of locations sorted along one of the
#' coordinates or the sum of multiple coordinates
#'
#' @param coordinate integer or vector of integers in {1,...,d}. If a single integer,
#' coordinates are ordered along that coordinate. If multiple integers,
#' coordinates are ordered according to the sum of specified coordinate values. For example,
#' when \code{d=2}, \code{coordinate = c(1,2)} orders from bottom left to top right.
#' @inheritParams order_dist_to_point
#' @return A vector of indices giving the ordering, i.e.
#' the first element of this vector is the index of the first location.
#' @examples
#' n <- 100             # Number of locations
#' d <- 2               # dimension of domain
#' locs <- matrix( runif(n*d), n, d )
#' ord1 <- order_coordinate(locs, 1 )
#' ord12 <- order_coordinate(locs, c(1,2) )
#' @export
order_coordinate <- function( locs, coordinate ){
    order(rowSums(locs[,coordinate,drop=FALSE]))
}



#' Maximum minimum distance ordering
#'
#' Return the indices of an exact maximum-minimum distance ordering.
#' The first point is chosen as the "center" point, minimizing L2 distance.
#' Dimensions d=2 and d=3 handled separately, dimensions d=1 and d>3 handled similarly.
#' Algorithm is exact and scales quasilinearly.
#'
#' @param locs Observation locations
#' @return A vector of indices giving the ordering, i.e.
#' the first element of this vector is the index of the first location.
#' @examples
#' n=100; locs <- cbind(runif(n),runif(n))
#' ord <- order_maxmin_exact(locs)
#'
#' @export
order_maxmin_exact<-function(locs){
  ord<-MaxMincpp(locs)
  return(ord)
}



## extension of the maxmin function, orders pred.locs last
# should be improved by extending MaxMincpp itself
#' @export
order_maxmin_exact_obs_pred<-function(locs, locs_pred){

  ord<-MaxMincpp(locs)
  ord_pred <-MaxMincpp(locs_pred)

  # The code below came from the method order_maxmin_obs_pred
  # The only difference is that the ord_pred (above) replaces
  # NN_pred (from knnx distances)
  locs_all = rbind(locs, locs_pred)

  n <- nrow(locs)
  m <- min( round(sqrt(n)), 200 )

  n_pred <- nrow(locs_pred)
  # next is to find 'ord_pred', a maxmin reordering of prediction locations
  NN <- FNN::get.knn( locs_all, k = m )$nn.index
  #NN_pred <- FNN::get.knnx( locs, locs_pred, k = 1 )$nn.dist
  # use ord, then order by NN_pred
  index_in_position <- c( ord, n + ord_pred, rep(NA,n_pred) )
  position_of_index <- order(index_in_position[1:(n+n_pred)])

  # move an index to the end if it is a
  # near neighbor of a previous location
  curlen <- n + n_pred
  nmoved <- 0
  for(j in (n+1):(n+2*n_pred) ){
    # nneigh tells us how many neighbors to look at
    # in order to decide whether the current point
    # has a previously ordered neighbor
    nneigh <- round( min(m,1*(n+n_pred)/(j-nmoved+1)) )
    neighbors <- NN[index_in_position[j],1:nneigh]
    if( min( position_of_index[neighbors], na.rm = TRUE ) < j ){
      nmoved <- nmoved+1
      curlen <- curlen + 1
      position_of_index[ index_in_position[j] ] <- curlen
      index_in_position[curlen] <- index_in_position[j]
      index_in_position[j] <- NA
    }
  }

  ord_pred <- index_in_position[ !is.na( index_in_position ) ][(n+1):(n+n_pred)] - n


  return(list(ord=ord, ord_pred =ord_pred))
}




## Maximum minimum distance ordering, observations then predictions
# Constrained so that observation locations ordered first, then prediction locations.
# order_maxmin_obs_pred <- function(locs, locs_pred, lonlat = FALSE,
#                                   space_time = FALSE, st_scale = NULL){
#
#   # FNN::get.knnx has strange behavior for exact matches
#   # so add a small amount of noise to each location
#   locs_all <- rbind( locs, locs_pred )
#   n_all <- nrow(locs_all)
#   ee <- min(apply( locs_all, 2, stats::sd ))
#   locs_all <- locs_all + matrix( ee*1e-4*stats::rnorm(n_all*ncol(locs_all)), n_all, ncol(locs_all) )
#
#   if(lonlat){ # convert lonlattime to xyztime or lonlat to xyz
#     lon <- locs_all[,1]
#     lat <- locs_all[,2]
#     lonrad <- lon*2*pi/360
#     latrad <- (lat+90)*2*pi/360
#     x <- sin(latrad)*cos(lonrad)
#     y <- sin(latrad)*sin(lonrad)
#     z <- cos(latrad)
#     if(space_time){
#       time <- locs_all[,3]
#       locs_all <- cbind(x,y,z,time)
#     } else {
#       locs_all <- cbind(x,y,z)
#     }
#   }
#
#   if(space_time){
#     d <- ncol(locs_all)-1
#     if( is.null(st_scale) ){
#       randinds <- sample(1:n_all, min(n_all,200) )
#       dvec <- c(fields::rdist( locs_all[randinds,1:d,drop=FALSE] ))
#       dvec <- dvec[ dvec > 0]
#       med1 <- mean(dvec)
#       dvec <- c(fields::rdist( locs_all[randinds, d + 1, drop=FALSE] ))
#       dvec <- dvec[ dvec > 0]
#       med2 <- mean(dvec)
#       st_scale <- c(med1,med2)
#     }
#     locs_all[ , 1:d] <- locs_all[ , 1:d]/st_scale[1]
#     locs_all[ , d+1] <- locs_all[ , d+1]/st_scale[2]
#   }
#
#   # get number of locs and redefine locs and locs_pred
#   n <- nrow(locs)
#   n_pred <- nrow(locs_pred)
#   locs <- locs_all[1:n,,drop=FALSE]
#   locs_pred <- locs_all[(n+1):(n+n_pred),,drop=FALSE]
#
#   m <- min( round(sqrt(n)), 200 )
#   # m is number of neighbors to search over
#   # get the past and future nearest neighbors
#   NN <- FNN::get.knn( locs, k = m )$nn.index
#   # pick a random ordering
#   index_in_position <- c( sample(n), rep(NA,n) )
#   position_of_index <- order(index_in_position[1:n])
#   # move an index to the end if it is a
#   # near neighbor of a previous location
#   curlen <- n
#   #curpos <- 1
#   nmoved <- 0
#   for(j in 2:(2*n) ){
#     nneigh <- round( min(m,n/(j-nmoved+1)) )
#     neighbors <- NN[index_in_position[j],1:nneigh]
#     if( min( position_of_index[neighbors], na.rm = TRUE ) < j ){
#       nmoved <- nmoved+1
#       curlen <- curlen + 1
#       position_of_index[ index_in_position[j] ] <- curlen
#       index_in_position[curlen] <- index_in_position[j]
#       index_in_position[j] <- NA
#     }
#   }
#   ord <- index_in_position[ !is.na( index_in_position ) ]
#
#   # now we have 'ord', a maxmin ordering of locs (observation locations)
#   # next is to find 'ord_pred', a maxmin reordering of prediction locations
#   NN <- FNN::get.knn( locs_all, k = m )$nn.index
#   NN_pred <- FNN::get.knnx( locs, locs_pred, k = 1 )$nn.dist
#   # use ord, then order by NN_pred
#   index_in_position <- c( ord, n + order(NN_pred,decreasing = TRUE), rep(NA,n_pred) )
#   position_of_index <- order(index_in_position[1:(n+n_pred)])
#   # move an index to the end if it is a
#   # near neighbor of a previous location
#   curlen <- n + n_pred
#   nmoved <- 0
#   for(j in (n+1):(n+2*n_pred) ){
#     # nneigh tells us how many neighbors to look at
#     # in order to decide whether the current point
#     # has a previously ordered neighbor
#     nneigh <- round( min(m,1*(n+n_pred)/(j-nmoved+1)) )
#     neighbors <- NN[index_in_position[j],1:nneigh]
#     if( min( position_of_index[neighbors], na.rm = TRUE ) < j ){
#       nmoved <- nmoved + 1
#       curlen <- curlen + 1
#       position_of_index[ index_in_position[j] ] <- curlen
#       index_in_position[curlen] <- index_in_position[j]
#       index_in_position[j] <- NA
#     }
#   }
#   ord_pred <- index_in_position[ !is.na( index_in_position ) ][(n+1):(n+n_pred)] - n
#   n+n_pred - nmoved
#
#
#   return(list(ord = ord, ord_pred = ord_pred))
# }
