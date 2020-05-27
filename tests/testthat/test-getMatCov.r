library(Matrix)
######### set parameters #########
set.seed(1996)
n = 2**2
m = 3

## covariance parameters
sig21 = 0.5; range1 = .15; smooth1 = 0.5; 
covparms1 = c(sig21,range1,smooth1)
covfun1.d = function(D) GPvecchia::MaternFun(D, covparms1)
covfun1 = function(locs) covfun1.d(fields::rdist(locs))


sig22 = 0.5; range2 = .15; smooth2 = 0.5; 
covparms2 = c(sig22,range2,smooth2)
covfun2.d = function(D) GPvecchia::MaternFun(D, covparms2)
covfun2 = function(locs) covfun2.d(fields::rdist(locs))



## generate grid of pred.locs
grid = seq(0,1,length = sqrt(n))
locs = as.matrix(expand.grid(grid,grid)) - 0.5

## set initial state
Qmat = covfun1(locs)
Qc = t(chol(Qmat))
Qc.s = as(t(Qc), "dgCMatrix")
Vmat = covfun2(locs)



vecchia.approx = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')


########## filtering ##########
M1 = getMatCov(vecchia.approx, Qc, factor = TRUE)
M2 = getMatCov(vecchia.approx, covfun1.d)
M3 = getMatCov(vecchia.approx, Qmat)
M4 = getMatCov(vecchia.approx, Qc.s, factor = TRUE)
M5 = getMatCov(vecchia.approx, Qc %*% t(Qc) + Vmat)
M6 = getMatCov(vecchia.approx, Vmat)


test_that("all three getMatCov calls should give the same result", {
  expect_equal(sum(abs(M2 - M3), na.rm = TRUE), 0)
  expect_equal(sum(abs(M2 - M3), na.rm = TRUE), 0)
  expect_equal(sum(abs(M1 - M2), na.rm = TRUE), 0)
  expect_equal(sum(abs(M2 - M4), na.rm = TRUE), 0)
  expect_equal(sum(abs(M1 + M6 - M5), na.rm = TRUE), 0)
})
