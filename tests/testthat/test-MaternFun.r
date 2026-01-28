locs = matrix(runif(200), ncol=2)
D = fields::rdist(locs)


naiveMatern05 = function(D, covparms) {

    if(covparms[3]!=0.5) stop("smoothness parameter incorrect. It should be 0.5")
    scaledD = D/covparms[2]
    return(exp(-scaledD)*covparms[1])
}


naiveMatern15 = function(D, covparms) {
    if(covparms[3]!=1.5) stop("smoothness parameter incorrect. It should be 1.5")
    scaledD = D/covparms[2]
    return(covparms[1] * (1+ sqrt(3)*scaledD)*exp(-sqrt(3)*scaledD))
}


naiveMatern25 = function(D, covparms) {
    if(covparms[3]!=2.5) stop("smoothness parameter incorrect. It should be 2.5")
    scaledD = D/covparms[2]
    return(covparms[1] * (1 + sqrt(5)*scaledD + 5/3*(scaledD**2))*exp(-sqrt(5)*scaledD))
}



sig2 = 1
range = 0.2


D05 = naiveMatern05(D, c(sig2, range, 0.5)) - MaternFun(D, c(sig2, range, 0.5))
D15 = naiveMatern15(D, c(sig2, range, 1.5)) - MaternFun(D, c(sig2, range, 1.5))
D25 = naiveMatern25(D, c(sig2, range, 2.5)) - MaternFun(D, c(sig2, range, 2.5))


test_that("GPvecchia's Matern Function has the correct special form for smoothness 0.5, 1.5 and 2.5", {
    expect_lt(sum(abs(D05)), 1e-10)
    expect_lt(sum(abs(D15)), 1e-10)
    expect_lt(sum(abs(D25)), 1e-10)
})


## nus <- 3 * runif(10)
## diffs <- rep(0, length(nus))
## for (n in seq_along(nus)) {
##     diffs[n] <- MaternFunGen(D, c(sig2, range, nus[n])) - MaternFun(D, c(sig2, range, nus[n]))
## }


## test_that("The version using C++ bessel and gamma functions works the same as the boost version", {
##     expect_lt(sum(abs(diffs)), 1e-8)
## })
