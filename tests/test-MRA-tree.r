rm(list=ls())

library(Rcpp)

setwd("/home/marcin/GPvecchia")
source("R/ordering_functions.R")
source("R/mraNN.r")
source("R/MRA_knot-tree.r")
source("R/MRA_tree-methods.r")
source("R/MRA_utility-functions.r")
sourceCpp("src/fastTree.cpp")

n=3**2
m=3
mra.options=NULL
mra.params = get.mra.params(n, mra.options, m)


# simulate locations
set.seed(10)
locs = cbind(runif(n),runif(n))

#knt.tree = knot.tree(locs, mra.params)
#mat = getNNmatrix(knt.tree)

NN0 = proc.time()
mat.new = generateNNarray(locs, mra.params[["J"]], mra.params[["M"]], mra.params[["r"]], m)
mat.new[mat.new==0]=NA

print(proc.time() - NN0)
