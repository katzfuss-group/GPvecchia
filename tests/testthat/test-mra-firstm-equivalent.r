#library(GPvecchia)

#Generate 1D data
n = 10
locs=matrix(runif(n),ncol=1)
covparms=c(0.8,0.2,1.3) #covparms=c(sig2,range,smooth)
Om0 <- MaternFun(fields::rdist(locs),covparms)
y=as.numeric(t(chol(Om0))%*%rnorm(n))+1
nuggets=rep(.1,n)
z = rnorm(n,y,sqrt(nuggets))
# Approximate posterior
m=2
# Low Rank
lr.vecchia.approx = vecchia_specify(locs,m,conditioning = "firstm")
preds.lr = vecchia_prediction(z,lr.vecchia.approx,covparms,nuggets)

# Equivalent MRA
mra.vecchia.approx=vecchia_specify(locs,m,conditioning = "mra", mra.options = list(r=c(m,1)) )
preds.mra =vecchia_prediction(z,mra.vecchia.approx,covparms,nuggets)

loc_ord=order(locs)

# Plot showing discrepancy
test_that("MRA version of low rank is not equivalent", {
  expect_equal(sum(abs(preds.mra$mu.obs[loc_ord]-preds.lr$mu.obs[loc_ord])),0)
})

