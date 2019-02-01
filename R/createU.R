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
  nuggets.all=c(nuggets,rep(0,sum(!obs)))
  nuggets.all.ord=nuggets.all[ord] # ordered nuggets for all locs
  nuggets.ord=nuggets.all[vecchia.approx$ord.z]# ordered nuggets for observed locs
  zero.nuggets=any(nuggets==0)

  # cannot condition on observation with zero nugget
  if(zero.nuggets){
    zero.cond=which(vecchia.approx$U.prep$revNNarray %in% which(nuggets.ord==0))
    vecchia.approx$U.prep$revCond[zero.cond]=TRUE
  }

  # call Rcpp function to create the nonzero entries of U
  if(vecchia.approx$conditioning=='mra' && is.matrix(covmodel)) {

    print("special case")
    rows = c()
    cols = c()
    for(i in 1:nrow(vecchia.approx$U.prep$revNNarray)){
      r = vecchia.approx$U.prep$revNNarray[i,]; r = r[1:length(r)-1]
      newrows = rep(i, sum(!is.na(r)))
      newcols = r[!is.na(r)]
      #rows = c(rows, newrows, i, rev(newcols))
      #cols = c(cols, newcols, i, rev(newrows))
      rows = c(rows, newrows, i)
      cols = c(cols, newcols, i)
    }
    #inds = cbind(rows, cols)
    #inds = as.vector(sapply(seq(nrow(inds)), function(r) inds[r,1]-1+n*(inds[r,2]-1)+1))

    #values = t(cbind(covmodel, covmodel[,(ncol(covmodel)-1):1]))
    values = t(cbind(covmodel))

    P = as(n:1, "pMatrix")

    M = sparseMatrix(i=rows, j=cols, x=values[!is.na(values)], dims=c(n, n), symmetric=TRUE)
    S = P %*% M %*% P
    L = t(chol(S))
    U = P %*% solve(t(L))

    browser()

    m = ncol(vecchia.approx$U.prep$revCond)-1

    nnz_per_row = tabulate(U@i + 1)
    inds.s = (m+1)*(seq(n)-1)+1
    inds.e = inds.s + nnz_per_row-1
    inds = rbind(inds.s, inds.e)
    inds = split(inds, rep(1:ncol(inds), each=nrow(inds)))
    inds = Reduce(c, Map(f = function(v) seq(v[1], v[2]), inds))
    newU = matrix(rep(0, (m+1)*n), ncol=n)
    newU[inds] = t(U)@x
    newU = t(newU)

    Zentries = as.vector(t(matrix(c((-1)/sqrt(nuggets.ord), 1/sqrt(nuggets.ord)), ncol=2)))
    not.na=c(!is.na(apply(vecchia.approx$U.prep$revNNarray, 1,rev)))
    Lentries=c(t(newU))[not.na]
    allLentries=c(Lentries, Zentries)

  } else {
    if(is.matrix(covmodel)) U.entries=U_NZentries_full_mat(vecchia.approx$U.prep$n.cores, n, vecchia.approx$locsord,
                                       vecchia.approx$U.prep$revNNarray, vecchia.approx$U.prep$revCond,
                                       nuggets.all.ord, nuggets.ord, covmodel[ord,ord], covparms)

    else if(is.character(covmodel)) U.entries=U_NZentries(vecchia.approx$U.prep$n.cores, n, vecchia.approx$locsord,
                                                    vecchia.approx$U.prep$revNNarray, vecchia.approx$U.prep$revCond,
                                                    nuggets.all.ord, nuggets.ord, covmodel, covparms)
    else stop("argument 'covmodel' type not supported")

    not.na=c(!is.na(apply(vecchia.approx$U.prep$revNNarray, 1,rev)))
    Lentries=c(t(U.entries$Lentries))[not.na]
    allLentries=c(Lentries, U.entries$Zentries)

  }

  # create sparse U matrix
  U=sparseMatrix(i=vecchia.approx$U.prep$colindices,j=vecchia.approx$U.prep$rowpointers,
                x=allLentries,dims=c(size,size))

  S = Sigma.ord;# S[2,3]=S[3,2]=0
  Um = U[c(1, 3, 5), c(1, 3, 5)]
  print( Um %*% t(Um) %*% S )
  browser()


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
