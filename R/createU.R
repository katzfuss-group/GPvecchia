# this function creates the sparse U matrix for specific parameters
# based on the output from specify.vecchia()
# covmodel: currently implemented
#    matern: with covparms (var,range,smoothness)
#    esqe: exponential + squared exp with covparms (var1,range1,var2,range2)

createU <- function(vecchia.approx,covparms,nuggets,covmodel='matern') {

  n=sum(vecchia.approx$obs)
  size=vecchia.approx$U.prep$size
  if(length(nuggets)==1) nuggets=rep(nuggets,n)
  nuggets.all=c(nuggets,rep(0,sum(!vecchia.approx$obs)))
  nuggets.ord=nuggets.all[vecchia.approx$ord]

  # call Rcpp function to create the nonzero entries of U
  U.entries=U_NZentries(vecchia.approx$U.prep$n.cores,n,vecchia.approx$locsord,
          vecchia.approx$U.prep$revNNarray,vecchia.approx$U.prep$revCond,
          nuggets.ord,covmodel,covparms)

  not.na=c(!is.na(apply(vecchia.approx$U.prep$revNNarray, 1,rev)))
  Lentries=c(t(U.entries$Lentries))[not.na]
  allLentries=c(Lentries, U.entries$Zentries)

  U=sparseMatrix(i=vecchia.approx$U.prep$colindices,j=vecchia.approx$U.prep$rowpointers,
                x=allLentries,dims=c(size,size))
  return(U)

}
