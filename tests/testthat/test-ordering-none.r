#Generate 1D data
n = 10
locs=matrix(runif(n),ncol=1)
covparms=c(0.8,0.2,1.3)
Om0 <- MaternFun(fields::rdist(locs),covparms)
y=as.numeric(t(chol(Om0))%*%rnorm(n))+1
nuggets=rep(.1,n)
z = rnorm(n,y,sqrt(nuggets))

# Approximate posterior
m=2

# Equivalent MRA
mra.vecchia.approx=vecchia_specify(locs,m,ordering='none')#,conditioning = "mra",
                                   #mra.options = list(r=c(m,1)) )

#assert( "ordering is unchanged", sum(abs(mra.vecchia.approx$ord - seq(n)))==0)

test_that("ordering of locations is unchanged", {
          expect_equal(sum(abs(mra.vecchia.approx$ord - seq(n))), 0)
})
