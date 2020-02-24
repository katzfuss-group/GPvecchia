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




cluster.equal = function(locs, size, K=NULL, dim.start=2){

  n = nrow(locs)
  if(is.null(n)) return(c(1))
  if(!is.null(K)) size = ceiling(n/K)
  else K = n/size

  J = 2**ceiling(base::log(K,2))
  if(!K==J){
    warning(paste("the number of subregions fo the original domain is ", K, " but has to be a power of 2. Setting J=", J, sep=""))
  }
  regions = list(seq(1:nrow(locs)))
  for(power in 1:base::log(J,2)){
    new.regions = vector("list", 2**power)
    for( reg.id in 1:length(regions)) {
      region = regions[[reg.id]]
      if(ncol(locs)==2 && (power %% 2)==1) d=2*(dim.start==2)+1*(dim.start==1) else if(ncol(locs)==2) d=2*(dim.start==1)+1*(dim.start==2) else d=1
      locs.in.region = matrix(locs[region,], ncol=ncol(locs))
      cutoff = stats::quantile(locs.in.region[,d], 0.5)
      new.regions[[2*reg.id]] = region[which(locs.in.region[,d] > cutoff)]
      new.regions[[2*reg.id-1]] = region[which(locs.in.region[,d] < cutoff)]

      on.the.border = region[which(locs.in.region[,d]==cutoff)]
      len.diff = length(new.regions[[2*reg.id]]) - length(new.regions[[2*reg.id-1]])
      if(len.diff > 0 ){
        new.regions[[2*reg.id-1]] = c(new.regions[[2*reg.id-1]], on.the.border[1:len.diff])
        on.the.border = on.the.border[-(1:len.diff)]
      } else if(len.diff < 0 ){
        new.regions[[2*reg.id]] = c(new.regions[[2*reg.id]], on.the.border[1:-len.diff])
        on.the.border = on.the.border[-(1:-len.diff)]
      }
      if(length(on.the.border)>0){
        one.side = on.the.border[1:floor(length(on.the.border)/2)]
        other.side = on.the.border[-(1:floor(length(on.the.border)/2))]
        new.regions[[2*reg.id-1]] = c(new.regions[[2*reg.id-1]], one.side)
        new.regions[[2*reg.id]] = c(new.regions[[2*reg.id]], other.side)
      }
    }
    regions = new.regions
  }
  ## return data in a format consistent with kmeans clustering
  clusters = rep(0, nrow(locs))
  id = 1
  for(region in regions){
    clusters[region] = id
    id = id + 1
  }
  if( any(table(clusters)>size))  warning(paste("Something might be wrong with clustering. Max cluster size is ", max(table(clusters)), " and should be ", size, sep=""))
  return(clusters)
}


#' extract the required elements from the covariance matrix
#'
#' This function takes the entire covariance matrix and creates
#' a matrix of covariances based on the vecchia approximatino object
#' @param V the object returned by vecchia_specify
#' @param Sigma The full covariance matrix
#'
#' @return matrix of size n x (m+1) with only those elements that 
#' are used by the incomplete Cholesky decomposition
#' 
#' @export
getMatCov = function(V, Sigma){
  revNNarray = V$U.prep$revNNarray
  rows = c()
  cols = c()
  for(i in 1:nrow(revNNarray)){
    r = revNNarray[i,];
    newrows = rep(i, sum(!is.na(r)))
    newcols = r[!is.na(r)]
    rows = c(rows, newrows)
    cols = c(cols, newcols)
  }
  n = dim(Sigma)[1]
  inds = cbind(rows, cols)
  inds = as.vector(sapply(seq(nrow(inds)), function(r) inds[r,2]-1+n*(inds[r,1]-1)+1))
  
  Sig.sel = matrix(rep(NA, length(revNNarray)), nrow=ncol(revNNarray))
  inds_to_fill=which(!is.na(t(revNNarray)))
  Sigma.ord = Sigma[V$ord, V$ord]
  Sig.sel[inds_to_fill] = Sigma.ord[inds]
  Sig.sel = t(Sig.sel)
  return(Sig.sel)
}



