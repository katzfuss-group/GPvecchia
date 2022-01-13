suppressPackageStartupMessages(library(RandomFields))
set.seed(1996)
n = 50**2
m = 50

sig2 = 1.0; range = 0.2; smooth = 1.5
a = 0.01

grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside))

x = RFsimulate(model = RMmatern(nu = smooth, scale = range, var = sig2),
                              x = locs[,1], y = locs[,2], spConform = FALSE)
y = rgamma(n, shape = a, rate = a*exp(-x))

mra = vecchia_specify(locs, m, conditioning = 'mra')
##the same thing happens when I don't use the  mra

test1 = function() {
    vecchia_laplace_likelihood(y, mra, likelihood_model="gamma",
                               covparms = c(sig2,range,smooth),
                               likparms = list(alpha=a))
}
expect_error(test1(), "Data invalid for likelihood type. Make sure that your data lies in the support of the likelihood function.")

