## Test function
getL00 = function(vecchia.approx, covfun, locs){
  Sig.sel = GPvecchia::getMatCov(vecchia.approx, covfun)
  inds = Filter(function(i) !is.na(i), as.vector(t(vecchia.approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(vecchia.approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(Sig.sel)))
  vals = GPvecchia::ic0(ptrs, inds, cov.vals)
  Laux = Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  ro = order(vecchia.approx$ord)
  return(Laux[ro, ])
}

set.seed(1988)
spatial.dim = 2
n = 20**2
m = 30

## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside))

## covariance parameters
sig2 = 1; range = .25; smooth = 0.5; 
covparms = c(sig2, range, smooth)
covfunErr = function(locs) GPvecchia::MaternFun(fields::rdist(locs), covparms)
covfun = function(locs1, locs2) as.numeric(GPvecchia::MaternFun(matrix(fields::rdist.vec(locs1, locs2), ncol = 1), covparms))
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)

Sig0 = exp(-fields::rdist(locs)/range)

 


## define Vecchia approximation
vecchia = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm')
L00 = getL00(vecchia, covfun.d, locs) / sqrt(sig2)

matCov = getMatCov(vecchia, covfun.d)
L = createL(vecchia, covfun)



expect_lt(max(abs(L - L00)), 1e-10)
expect_lt(max(abs(Sig0 - L00 %*% Matrix::t(L00))), 1e-10)
expect_lt(max(abs(Sig0 - L %*% Matrix::t(L))), 1e-10)

expect_error(createL(vecchia, covfunErr),
             "The suplied covariance function has to have two arguments")
expect_error(createL(vecchia, "squaredExp"),
             "Argument covmodel has incorrect format")
