library(Matrix)
######### set parameters #########
set.seed(1996)
n = 100**2
m = 50

## covariance parameters
sig2 = 0.5; range = .15; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)
covfun = function(locs) covfun.d(fields::rdist(locs))

## generate grid of pred.locs
grid = seq(0,1,length = sqrt(n))
locs = as.matrix(expand.grid(grid,grid)) 

## set initial state
Qmat = covfun(locs)
Qc = t(chol(Qmat))
Qc.s = as(t(Qc), "dgCMatrix")

vecchia.approx = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')


########## filtering ##########
#M1 = getMatCov(vecchia.approx, Qc, factor = TRUE)
M2 = getMatCov(vecchia.approx, covfun.d)
#M3 = getMatCov(vecchia.approx, Qmat)
M4 = getMatCov(vecchia.approx, Qc.s, factor = TRUE)


test_that("all three getMatCov calls should give the same result", {
  #expect_equal(sum(M2 - M3, na.rm = TRUE), 0)
  #expect_equal(sum(M2 - M3, na.rm = TRUE), 0)
  #expect_equal(sum(M1 - M2, na.rm = TRUE), 0)
  expect_equal(sum(M2 - M4, na.rm = TRUE), 0)
})


#test_that("all three getMatCov calls should give the same result", {
#  expect_equal(sum(M2 - M3, na.rm = TRUE), 0)
#  expect_equal(sum(M1 - M2, na.rm = TRUE), 0)
#  expect_equal(sum(M1 - M4, na.rm = TRUE), 0)
#})