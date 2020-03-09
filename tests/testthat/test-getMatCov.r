######### set parameters #########
set.seed(1996)
n = 34**2
m = 50

## covariance parameters
sig2 = 0.5; range = .15; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) MaternFun(fields::rdist(locs),covparms)


## generate grid of pred.locs
grid = seq(0,1,length = sqrt(n))
locs = as.matrix(expand.grid(grid,grid)) 

## set initial state
Qmat = covfun(locs)

vecchia.approx = vecchia_specify(locs, m, conditioning = 'mra', verbose = TRUE)


########## filtering ##########
M1 = getMatCov(vecchia.approx, t(chol(Qmat)), factor = TRUE)
M2 = getMatCov(vecchia.approx, function(d) MaternFun(d, covparms))
M3 = getMatCov(vecchia.approx, Qmat)


test_that("all three getMatCov calls should give the same result", {
  expect_equal(sum(M2 - M3, na.rm = TRUE), 0)
  expect_equal(sum(M1 - M2, na.rm = TRUE), 0)
})