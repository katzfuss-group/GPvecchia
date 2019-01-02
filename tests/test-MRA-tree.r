rm(list=ls())

setwd("/home/marcin/GPvecchia")
source("R/ordering_functions.R")
source("R/MRA_utility-functions.r")
source("R/MRA_tree-methods.r")
source("R/MRA_knot-tree.r")
spatial.dim=2

r = c(20, 10, 4, 4)
m = sum(r)
M = length(r)-1
J = rep(2, M)

#r.tot = sum(cumprod(J)*r[-1])

n = 40

mra.params = get.mra.params(n, mra.options, m)


# simulate locations
set.seed(10)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
  ord = order_maxmin(locs)
  locs = locs[ord,]
}


NNarray = matrix(rep(NA, m*n), ncol=m)
knots = list()

inds = genInds(M,J)
remaining = list()
remaining[["r"]] = seq(n)

while(length(remaining)>0){
  id = names(remaining)[1]; m = res(id)
  reg.inds = remaining[[id]]

  r.eff = min(r, length(reg.inds))
  knots[[id]] = reg.inds[seq(r.eff)]

  reg.inds = reg.inds[-seq(r.eff)]; reg.locs = locs[reg.inds,]

  if(res(id)<M){

    clusters =  cluster.equal(reg.locs, K=J[m+1], dim.start=m%%2+1)

    for(child.no in 1:max(clusters)){
      child.id = paste(c(id, child.no), collapse="_")
      remaining[[child.id]] = reg.inds[clusters==child.no]
    }
  }

  remaining = remaining[-1]
}


NNmatrix = getNNmatrix(knots)
