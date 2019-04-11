# this function creates the sparse U matrix for specific parameters
# based on the output from specify.vecchia()
# covmodel: currently implemented
#    matern: with covparms (var,range,smoothness)
#    esqe: exponential + squared exp with covparms (var1,range1,var2,range2)

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

  # cannot condition on observation with zero nugget
  if(zero.nuggets){
    zero.cond=which(vecchia.approx$U.prep$revNNarray %in% which(nuggets.ord==0))
    vecchia.approx$U.prep$revCond[zero.cond]=TRUE
  }

  if(vecchia.approx$conditioning=="mra"){
    new = proc.time()
    inds = Filter(function(i) !is.na(i), as.vector(t(vecchia.approx$U.prep$revNNarray - 1)))
    ptrs = c(0, cumsum(apply(vecchia.approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))

    if(is.matrix(covmodel)){
      cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
      vals = createUcppM(ptrs, inds, cov.vals)
    } else if(is.function(covmodel)){
      ## function has to be of a certain form, specifically, it has to be able
      ## to take k pairs of locations and return a vector with distances which is
      ## of length k.

      # create a vector with pairs of locations
      #s1 = proc.time()
      f = function(r) rep(r[length(r)], length(which(!is.na(r))))
      inds1 = Filter(function(i) !is.na(i), as.vector(t(vecchia.approx$U.prep$revNNarray)))
      inds2 = unlist(apply(vecchia.approx$U.prep$revNNarray, 1, f))
      locs1 = vecchia.approx$locsord[inds1,]
      locs2 = vecchia.approx$locsord[inds2,]

      cov.vals = covmodel(locs1, locs2)
      vals = createUcppM(ptrs, inds, cov.vals)
      #print(proc.time()-s1)
    } else {
      #s2 = proc.time()
      vals = createUcpp(ptrs, inds, vecchia.approx$locsord)
      #print(proc.time() - s2)
    }

    Laux = sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
    Ulatent = t(solve(Laux, sparse=TRUE))

    N = nrow(vecchia.approx$U.prep$revNNarray)

    p1 = c(0, cumsum(rep(2,N))) + ptrs
    p2 = c(p1[-1] -2, NA)
    LZp = Filter(function(i) !is.na(i), c(rbind(p1,p2)))

    nuggets.inds = c(rbind(p2[1:N], p2[1:N]+1))+1
    nvals = length(vals) + 2*N
    LZvals = rep(0, nvals)
    LZvals[nuggets.inds] = c(rbind(-1/sqrt(nuggets.ord), 1/sqrt(nuggets.ord)))
    LZvals[-nuggets.inds] = Ulatent@x
    LZinds = Filter(function(i) !is.na(i), c(rbind(2*t(vecchia.approx$U.prep$revNNarray - 1), c(2*seq(N)-2), 2*seq(N)-1)))

    U = t(sparseMatrix(j=LZinds, p=LZp, x=LZvals, index1=FALSE))
    new = proc.time() - new
  } else {
    old = proc.time()
    # call Rcpp function to create the nonzero entries of U
    if(is.matrix(covmodel)) U.entries=U_NZentries_mat(vecchia.approx$U.prep$n.cores, n, vecchia.approx$locsord,
                                                      vecchia.approx$U.prep$revNNarray, vecchia.approx$U.prep$revCond,
                                                      nuggets.all.ord, nuggets.ord, covmodel[ord,ord], covparms)
    else if(is.character(covmodel)) U.entries=U_NZentries(vecchia.approx$U.prep$n.cores, n, vecchia.approx$locsord,
                                                      vecchia.approx$U.prep$revNNarray, vecchia.approx$U.prep$revCond,
                                                      nuggets.all.ord, nuggets.ord, covmodel, covparms)
    else stop("argument 'covmodel' type not supported")


    # create sparse U matrix
    not.na=c(!is.na(apply(vecchia.approx$U.prep$revNNarray, 1,rev)))
    Lentries=c(t(U.entries$Lentries))[not.na]
    allLentries=c(Lentries, U.entries$Zentries)
    U=sparseMatrix(i=vecchia.approx$U.prep$colindices,j=vecchia.approx$U.prep$rowpointers,
                  x=allLentries,dims=c(size,size))
    old = proc.time() - old

    print(new)
    print(old)

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
