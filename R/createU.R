#' create the sparse triangular U matrix for specific parameters
#'
#' @param vecchia.approx: object returned by \code{\link{vecchia_specify}}
#' @param covparms: vector of covariance parameters
#' @param nuggets: nugget variances -- if a scalar is provided, variance is assumed constant
#' @param covmodel: covariance model. currently implemented:
#    matern: with covparms (var,range,smoothness)
#    esqe: exponential + squared exp with covparms (var1,range1,var2,range2)
#'
#' @return list containing the sparse upper triangular U,
#'     plus additional objects required for other functions
#' @examples
#' z=rnorm(9); locs=matrix(1:9,ncol=1); vecchia.approx=vecchia_specify(locs,m=5)
#' U.obj=createU(vecchia.approx,covparms=c(1,2,.5),nuggets=.2)
#' @export
createU <- function(vecchia.approx,covparms,nuggets,covmodel='matern') {

  n=sum(vecchia.approx$obs)
  size=vecchia.approx$U.prep$size
  latent = (1:size) %in% vecchia.approx$U.prep$y.ind
  ord=vecchia.approx$ord
  obs=vecchia.approx$obs

  # turn nugget into necessary format
  if(length(nuggets)==1) nuggets=rep(nuggets,n)
  nuggets.all=c(nuggets,rep(0,sum(latent)-n))
  if(vecchia.approx$cond.yz=='zy'){ ord.all=c(ord[1:n],ord+n) } else { ord.all=ord }
  nuggets.all.ord=nuggets.all[ord.all] # ordered nuggets for all locs
  nuggets.ord=nuggets.all[vecchia.approx$ord.z]# ordered nuggets for observed locs
  zero.nuggets=any(nuggets==0)
  n.obs = length(nuggets.ord)

  # cannot condition on observation with zero nugget
  if(zero.nuggets){
    zero.cond=which(vecchia.approx$U.prep$revNNarray %in% which(nuggets.ord==0))
    vecchia.approx$U.prep$revCond[zero.cond]=TRUE
  }

  # in the mra case we calculate the U matrix using incomplete Cholesky (ic0)
  if(vecchia.approx$conditioning=="mra"){
    inds = Filter(function(i) !is.na(i), as.vector(t(vecchia.approx$U.prep$revNNarray - 1)))
    ptrs = c(0, cumsum(apply(vecchia.approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))

    if(is.matrix(covmodel)){
      cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
      vals = createUcppM(ptrs, inds, cov.vals)
    } else if(is.function(covmodel)){
      ## function has to be of a certain form, specifically, it has to be able
      ## to take k pairs of locations and return a vector with distances which is
      ## of length k.
      f = function(r) rep(r[length(r)], length(which(!is.na(r))))
      inds1 = Filter(function(i) !is.na(i), as.vector(t(vecchia.approx$U.prep$revNNarray)))
      inds2 = unlist(apply(vecchia.approx$U.prep$revNNarray, 1, f))
      locs1 = vecchia.approx$locsord[inds1,]
      locs2 = vecchia.approx$locsord[inds2,]
      cov.vals = covmodel(locs1, locs2)
      vals = createUcppM(ptrs, inds, cov.vals)
    } else {
      vals = createUcpp(ptrs, inds, vecchia.approx$locsord, covparms)
    }

    Laux = Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
    Ulatent = Matrix::t(Matrix::solve(Laux, sparse=TRUE))

    N = nrow(vecchia.approx$U.prep$revNNarray)

    # build new matrix
    p.latent = c(0, cumsum(rep(2,n.obs)), rep(2*n.obs, N-n.obs)) + Ulatent@p
    p.obs = c(p.latent[2:(n.obs+1)] -2, rep(NA, N-n.obs+1))
    LZp = Filter(function(i) !is.na(i), c(rbind(p.latent,p.obs)))

    nuggets.inds = c(rbind(p.obs[1:n.obs], p.obs[1:n.obs]+1))+1
    nvals = length(Ulatent@x) + 2*n.obs

    LZvals = rep(0, nvals)
    LZvals[nuggets.inds] = c(rbind(-1/sqrt(nuggets.ord), 1/sqrt(nuggets.ord)))
    LZvals[-nuggets.inds] = Ulatent@x

    LZinds = rep(NA, nvals)
    LZinds[nuggets.inds] = seq(2*n.obs)-1

    # calculate indices for latent variables
    inds.lt.2nobs = which(Ulatent@i<n.obs)
    new.latent = rep(NA, nvals - length(nuggets.inds))
    new.latent[inds.lt.2nobs] = 2*Ulatent@i[inds.lt.2nobs]
    new.latent[-inds.lt.2nobs] = Ulatent@i[-inds.lt.2nobs]+n.obs
    LZinds[-nuggets.inds] = new.latent

    U = Matrix::t(Matrix::sparseMatrix(j=LZinds, p=LZp, x=LZvals, index1=FALSE))

  } else {

    if(is.matrix(covmodel)) U.entries=U_NZentries_mat(vecchia.approx$U.prep$n.cores, n, vecchia.approx$locsord,
                                                      vecchia.approx$U.prep$revNNarray, vecchia.approx$U.prep$revCond,
                                                      nuggets.all.ord, nuggets.ord, covmodel, covparms)
    else if(is.character(covmodel)) U.entries=U_NZentries(vecchia.approx$U.prep$n.cores, n, vecchia.approx$locsord,
                                                      vecchia.approx$U.prep$revNNarray, vecchia.approx$U.prep$revCond,
                                                      nuggets.all.ord, nuggets.ord, covmodel, covparms)
    else stop("argument 'covmodel' type not supported")


    # create sparse U matrix
    not.na=c(!is.na(apply(vecchia.approx$U.prep$revNNarray, 1,rev)))
    Lentries=c(t(U.entries$Lentries))[not.na]
    allLentries=c(Lentries, U.entries$Zentries)
    U=Matrix::sparseMatrix(i=vecchia.approx$U.prep$colindices,j=vecchia.approx$U.prep$rowpointers,
                  x=allLentries,dims=c(size,size))
  }



  # for zy ordering, remove rows/columns corresponding to dummy y's
  if(vecchia.approx$cond.yz=='zy') {
    dummy=2*(1:n)-1
    U=U[-dummy,-dummy]
    latent=latent[-dummy]
    obs=obs[-((1:n)+n)]
  }

  # remove rows/columns corresponding to zero nugget and store related info
  zero.nugg=list()
  if(zero.nuggets){

    # find and remove rows/columns corresponding to zero nugget
    inds.U=which(Matrix::diag(U)==Inf)
    cond.on=apply(U[,inds.U],2,function(x) min(which(x!=0)))
    U=U[-inds.U,-inds.U]

    # identify corresponding indices
    inds.z=which((1:size)[!latent] %in% inds.U)
    inds.locs=which((1:size)[latent] %in% cond.on)
    zero.nugg=list(inds.U=inds.U,inds.z=inds.z,inds.locs=inds.locs)

    # modify other quantities accordingly
    latent[cond.on]=FALSE
    latent=latent[-inds.U]
    ord=c(ord[-inds.locs],ord[inds.locs])
    obs=c(obs[-inds.locs],obs[inds.locs])

  }

  # return object
  U.obj=list(U=U,latent=latent,ord=ord,obs=obs,zero.nugg=zero.nugg,
             ord.pred=vecchia.approx$ord.pred,ord.z=vecchia.approx$ord.z,
             cond.yz=vecchia.approx$cond.yz)
  return(U.obj)

}
