## Test function
getL00 = function(vecchia.approx, covfun, locs){
  Sig.sel = getMatCov(vecchia.approx, covfun)
  inds = Filter(function(i) !is.na(i), as.vector(t(vecchia.approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(vecchia.approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(Sig.sel)))
  vals = ic0(ptrs, inds, cov.vals)
  Laux = Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  ro = order(vecchia.approx$ord)
  return(Laux[ro, ])
}

set.seed(1988)
spatial.dim = 2
n = 20**2
m = 30

## covariance parameters
sig2 = 1; range = .25; smooth = 0.5; 
covparms = c(sig2, range, smooth)
covfunErr = function(locs) MaternFun(fields::rdist(locs), covparms)
covfun = function(locs1, locs2) MaternFun(fields::rdist(locs1, locs2), covparms)
covfun.d = function(D) MaternFun(D, covparms)



## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 


## define Vecchia approximation
mra = vecchia_specify(locs, m, conditioning = 'mra')
L00 = getL00(mra, covfun.d, locs) / sqrt(sig2)

matCov = getMatCov(mra, covfun.d)
L = createL(mra, matCov)

expect_lt(max(abs(L - L00)), 1e-10)

expect_error(createL(mra, covfunErr),
             "The suplied covariance function has to have two arguments")
expect_error(createL(mra, "squaredExp"),
             "Argument covmodel has incorrect format")
expect_error(createL(mra, covfun),
             "The supplied covariance function is in a wrong format")
