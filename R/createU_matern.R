# this function creates the sparse U matrix for specific parameters
# based on the output from specify.vecchia()
# currently, only isotropic matern works

createU_matern <- function(vecchia.approx,covparms,nuggets) {

  n=sum(vecchia.approx$obs)
  size=vecchia.approx$U.prep$size
  if(length(nuggets==1)) nuggets=rep(nuggets,n)
  nuggets.all=c(nuggets,rep(0,sum(!vecchia.approx$obs)))
  nuggets.ord=nuggets.all[vecchia.approx$ord]
  
  # call Rcpp function to create the nonzero entries of U
  LmatZ=U_NZentries(vecchia.approx$U.prep$n.cores,n,vecchia.approx$locsord,
          vecchia.approx$U.prep$revNNarray,vecchia.approx$U.prep$revCond,
          nuggets.ord,covparms)
  
  not.na=c(!is.na(apply(vecchia.approx$U.prep$revNNarray, 1,rev)))
  Lentries=c(t(LmatZ$Lentries))[not.na]
  allLentries=c(Lentries, LmatZ$Zentries)
  
  U=sparseMatrix(i=vecchia.approx$U.prep$colindices,j=vecchia.approx$U.prep$rowpointers,
                x=allLentries,dims=c(size,size))
  return(U)
  
}
