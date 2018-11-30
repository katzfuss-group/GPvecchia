findOrderedNN_mra = function(locs, m, mra.options){

  n = length(locs)/ncol(locs)
  mra.params = get.mra.params(n, m, mra.options)

  if((length(mra.params$r)==2 && mra.params$r[2]==0) || (mra.params$J!=2 && mra.params$J!=4)){
    if( mra.params$M>1) warning("When J is neither 2 nor 4 we always set M to 1 and use the Full scale approximation")
    ind.tree = domain.tree.FSA(locs, mra.params$r[1])
  } else if( mra.params[['J']]==2 ) ind.tree = domain.tree.J2(locs, mra.params)
  else if( mra.params[['J']]==4 ) ind.tree = domain.tree.J4(locs, mra.params)

  knt.tree = knot.tree(ind.tree, mra.params[['r']], dim=ncol(locs))
  plot.locs.tree(ind.tree, locs, knots=knt.tree)

  mat = getNNmatrix(knt.tree)

  return(mat)
}



choose.M = function(n, m) {

  M=1
  while(2^(M+1)/(M+1) <= n/m) M=M+1

  ## for very small m:
  if(M+1>m) {
    M=m-1
    r=rep(1,M+1)
    J=c(1,rep(2,M-1),ceiling((n-sum(2^(0:(M-1))))/2^(M-1)))
  } else{
    J=c(1,rep(2,M))

    ## choose r based on m
    r=rep(ceiling(m/(M+1)),M+1)
    l=M
    while(sum(r)>m) {
      r[l]=r[l]-1
      l=l-1
    }
  }

  ### check that choices are valid
  if(sum(r)>m | sum(r*cumprod(J))<n) print('ERROR')
  else return(list(M=M, r=r, J=J))
}


get.mra.params = function(n,m, opts){

  params = list(m=opts$m)
  if(is.null(opts$J)) {
    J=2
    warning("J not specified. Setting J=2")
  } else J=opts$J
  if(is.null(opts$M) ){
    if(is.null(opts$r)) {
      pars = choose.M(n,m)
      # M = 0
      # while(J^(M+1)/(M+1) <= n/m ) M=M+1
      # if(M+1>m){
      #   M = m-1
      #   r = rep(1,M+1)
      #   J = c(1, rep(2, M-1), ceiling((n-sum(2^(0:M-1)))/2^(M-1)))
      # } else {
      #   J = c(1, rep(2,M))
      #   r = rep(ceiling(m/(M+1)), M+1)
      #   l = M
      #   while(sum(r)>m) {
      #     r[l]=r[l]-1
      #     l=l-1
      #   }
      # }
      #print(paste(r, " vs. ", pars$r, sep=""))
      r = pars$r
      M = pars$M
      J = pars$J
      print(J)
      #print(paste(r, " vs. ", pars$M, sep=""))
    } else {
      r = opts$r
      if(length(r)>1) M=length(r)-1
      else if(length(J)>1) M=length(J)
      else M = floor((log(n/r)*(J-1)+1)/log(J))-1
    }
  } else if(is.null(opts$r)) {
    M = opts$M
    r = floor(m/M)
  } else {
    warning("M, r set for MRA. Parameter m will be overridden")
    M = opts$M; r = opts$r
  }
  params[['J']] = J; params[['M']] = M; params[['r']] = r

  return(params)
}



########   domain tree functions  #########

domain.tree.FSA = function(locs, r){
  domain.tree=list(r=seq(nrow(locs)))
  for( idx in 1:(nrow(locs)-r) ){
    locs.idx = idx + r
    id = paste("r",idx,sep="_")
    domain.tree[[id]] = locs.idx
  }

  D = fields::rdist(locs[-(1:r),], locs[1:r,])
  knot.regions = apply(D, 2, which.min)
  for( knot.idx in seq(r,1)) {
    region = knot.regions[knot.idx]
    region.id = paste("r", region, sep="_")
    domain.tree[[region.id]] = c(knot.idx, domain.tree[[region.id]])

  }
  return(domain.tree)
}



domain.tree.J4 = function( locs, mra.options ){

  r = mra.options[['r']]; J = 4; n = length(locs)/ncol(locs)
  points = seq(n)

  M = floor(log((n/r)*(J-1) + 1)/log(J)) # not sure if -1 should be left in or out
  #M = floor(log((n/r)*(J-1) + 1)/log(J)) - 1
  if( M==0 ) stop(paste(c('ERROR: n=', n, ' points is not enough for J=4 and r=', r,' basis functions'), collapse=""))
  grid.tree = list(r=points)
  inds = genInds(M)

  for( ind in inds){

    if( child.id(ind)==1 ){

      par.inds = grid.tree[[parent(ind)]]
      par.locs = locs[par.inds,]

      if( ncol(locs)==2 ) {
        x_split = quantile(par.locs[,1], 0.5)
        y_split = quantile(par.locs[,2], 0.5)

        reg1 = par.inds[which(par.locs[,1]<=x_split & par.locs[,2]<=y_split)]
        reg2 = par.inds[which(par.locs[,1]> x_split & par.locs[,2]<=y_split)]
        reg3 = par.inds[which(par.locs[,1]<=x_split & par.locs[,2]> y_split)]
        reg4 = par.inds[which(par.locs[,1]> x_split & par.locs[,2]> y_split)]

        subregs = list(reg1, reg2, reg3, reg4)

      } else if( ncol(locs)==1 ) {

        x_split_1 = quantile(par.locs, 0.25)
        x_split_2 = quantile(par.locs, 0.5)
        x_split_3 = quantile(par.locs, 0.75)

        reg1 = par.inds[which(par.locs<= x_split_1)]
        reg2 = par.inds[which(par.locs>x_split_1 & par.locs<=x_split_2)]
        reg3 = par.inds[which(par.locs>x_split_2 & par.locs<=x_split_3)]
        reg4 = par.inds[which(par.locs>x_split_3)]

        subregs = list(reg1, reg2, reg3, reg4)
      }
    }

    grid.tree[[ind]] = subregs[[child.id(ind)]]
  }
  return(grid.tree)
}



domain.tree.J2 = function( locs, mra.options ){

  n = length(locs)/ncol(locs)
  points = seq(n)

  r = mra.options[['r']]
  J = mra.options[['J']]
  M = mra.options[['M']]
  grid.tree = list(r=points)
  inds = genInds(M,J=c(2))

  for( ind in inds ) {
    if( child.id(ind)==1 ){
      par.inds = grid.tree[[parent(ind)]]
      par.locs = locs[par.inds,]

      if( ncol(locs)==2 ){
        if( (res(ind) %% 2) == 1 ) {
          x_split = quantile(par.locs[,1],0.5)
          reg1 = par.inds[which(par.locs[,1]<=x_split)]
          reg2 = par.inds[which(par.locs[,1]>x_split)]
        } else {
          y_split = quantile(par.locs[,2],0.5)
          reg1 = par.inds[which(par.locs[,2]<=y_split)]
          reg2 = par.inds[which(par.locs[,2]>y_split)]
        }
      } else if(ncol(locs)==1 ) {
        x_split = quantile(par.locs, 0.5)
        reg1 = par.inds[which(par.locs<x_split)]
        reg2 = par.inds[which(par.locs>=x_split)]
      }
      subregs = list(reg1, reg2)
    }

    grid.tree[[ind]] = subregs[[child.id(ind)]]
  }
  return(grid.tree)
}

# domain.tree.FSA = function( locs, mra.options ){
#
#   J = mra.options[['J']]
#   n = length(locs)/ncol(locs)
#   points = seq(n)
#   M = 1
#   grid.tree = list(r=points)
#   centers = locs[1:J,]
#
#   D = fields::rdist(locs[-seq(J),], centers)
#   cents = apply(D, 1, which.min)
#   for( j in 1:J){
#     ind = paste("r", j, sep="")
#     grid.tree[[ind]] = c(j, J + which(cents==j))
#   }
#
#   return(grid.tree)
# }




#########   knot tree functions   #####

knot.tree = function(locs.tree, r, dim=2){

  M = get.M(locs.tree)
  Jm = get.Jm(locs.tree)
  knots = list()
  remaining = list(r=locs.tree[["r"]])
  exact = all(Jm==Jm[1]) && length(r)==1 && Jm[1]==r+1 && dim==1

  if( exact ) print("exact representation available!")

  for( ind in names(locs.tree) ){
    node.locs = locs.tree[[ind]]
    available = intersect(node.locs, remaining[[parent(ind)]])

    if(res(ind)==M ) knots[[ind]] = available       # if we are at the last resolution, everything is a knot

    else {                                          # if not at the last resolution:
      if( exact ){  #spectial case when we place knots at the split points in 1d
        children = Filter(function(s) startsWith(s, ind) && nchar(s)==nchar(ind)+1, names(locs.tree))
        knts = c()
        for( child in children[-1] ){
          locs.child.available = intersect(locs.tree[[child]], available)

          knts = c(knts, locs.child.available[1])
        }
        knots[[ind]] = knts
      } else {
        no_knots = getNKnt(r, ind)
        if(no_knots)  knots[[ind]] = available[1:min(no_knots, length(available))]
      }
    }

    remaining[[ind]] = setdiff(remaining[[parent(ind)]], knots[[ind]])
  }
  return(knots)
}




getNKnt = function(r, ind){
  if( length(r)==1 ) {
    r
  } else {
    m = res(ind)
    if(m+1>length(r)){
      r[length(r)]
    } else
      r[m+1]
  }
}



getNNmatrix = function(knot.tree){

  neighbors = list()
  #fill out the list of neighbors for the root

  #find the root(s) of the tree
  min.length = min(sapply(names(knot.tree), function(n) nchar(n)))
  roots = which(sapply(names(knot.tree), nchar)==min.length)

  for( root.no in roots ){
    root = names(knot.tree)[root.no]
    root.ind = knot.tree[[root]][1]
    neighbors[[root.ind]] = c(root.ind); cond.set=c(root.ind); last.knot=root.ind
    for( knot in knot.tree[[root]][-1]){
      cond.set = c(knot, cond.set)
      neighbors[[knot]] = cond.set
      last.knot = knot
    }
  }


  # once the first knot is handled, fill out the list
  # for the remaining knots
  for( ind in names(knot.tree)[-roots] ){
    knots = knot.tree[[ind]]
    parent.knots = knot.tree[[parent(ind)]]
    last.knot = parent.knots[length(parent.knots)]
    cond.set = neighbors[[parent.knots[length(parent.knots)]]]
    for( knot in knots) {
      cond.set = c(knot, cond.set)
      neighbors[[knot]] = cond.set
    }
  }
  list2matrix(neighbors)
  NNarray = list2matrix(neighbors)
  return(NNarray)
}





##########   tree methods    ########

parent = function(id){
  if( nchar(id) > 1){
    new_id = strsplit(id, "_")[[1]]
    assembly.elements = new_id[1:length(new_id)-1]
    return(paste(assembly.elements, collapse="_"))
  } else {
    return(id)
  }
}

child.id = function(id){
  new_id = strsplit(id, "_")[[1]]
  num = as.integer(new_id[length(new_id)])
  return(num)
}

res = function(id){
  res = length(strsplit(id, "_")[[1]])-1
  return(res)
}


get.M = function(tree){
  all.knots = names(tree)
  M = max(sapply(all.knots, function(x) length(strsplit(x,split="_")[[1]])))-1
  return(M)
}


get.Jm = function(tree, m=NULL){
  M = get.M(tree)
  all.knots = names(tree)
  lengths = sapply(all.knots, function(x) length(strsplit(x,split="_")[[1]]))
  Js = rep(0, M)
  for( m in seq(2,max(M+1, 2))) {
    inds.m = all.knots[lengths==m]
    maxJ = max(sapply(inds.m, function(ind) as.integer(strsplit(ind, split="_")[[1]][m])))
    Js[m-1] = maxJ
  }
  return(Js)
}


#### plotting methods ####

plot.locs.tree = function(locs.tree, locs, knots=NULL){
  if( ncol(locs)==2 ) plot.tree.2d(locs.tree, locs, knots)
  if( ncol(locs)==1 ) plot.tree.1d(locs.tree, locs, knots)


}



plot.tree.1d = function(locs.tree, locs, knots) {

  colors = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
  nc = length(colors)

  m = -1
  M = get.M(locs.tree)
  oldpar = par(mfrow=c(M+1, 1))

  xlim = c(min(locs), max(locs))

  for( ind in names(locs.tree) ){
    node.points = locs.tree[[ind]]
    if( length(node.points)==0) next
    if( res(ind) > m){
      m=res(ind)
      plot(locs[node.points], rep(0.1, length(node.points)), col=colors[1], pch=15,
           ylim=c(0, .5), xlim=xlim, yaxt="n",
           xlab="", ylab="", main=paste(c("resolution=", m), collapse=""))
    } else {
      points(locs[node.points], rep(0.1, length(node.points)), col=colors[1], pch=15)
    }

    if(!is.null(knots)){
      knot.points = knots[[ind]]
      points(locs[knot.points,1], rep(0.3, length(knot.points)), col="red", pch=4, cex=2)
    }


    colors = c(colors[2:nc], colors[1])
  }
  par(oldpar)
}


plot.tree.2d = function(locs.tree, locs, knots){

  colors = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
  nc = length(colors)

  m = -1
  for( ind in names(locs.tree) ){
    node.points = locs.tree[[ind]]
    if(length(node.points)==0) next
    if(res(ind)>m){
      m=res(ind)
      plot(locs[node.points,1], locs[node.points,2], col=colors[1], pch=16,
           xlim=c(0, 1), ylim=c(0, 1), xlab="", ylab="",
           main=paste(c("resolution=", m), collapse=""))
    } else {
      points(locs[node.points,1], locs[node.points,2], col=colors[1], pch=16)
    }

    if(!is.null(knots)){
      knot.points = knots[[ind]]
      points(locs[knot.points,1], locs[knot.points,2], col=colors[1], pch=4, cex=2)
    }

    colors = c(colors[2:nc], colors[1])
  }

}





########   utility functions   #######

## Converts a list of numeric vectors to an n x m matrix,
## where n is the length of the list and m is the length
## of the longest vector. Rows wth shorter vectors are
## padded with `padding`.
list2matrix = function(L, padding=NA){
  a.width = max(sapply(L, function(x) length(x)))
  a.height = length(L)
  A = matrix(rep(padding, a.width*a.height), ncol=a.width)

  for( ind in 1:length(L) ){
    elem = L[[ind]]
    A[ind,1:length(elem)] = elem
  }
  return(A)
}


## generate a list containg a hierarchy of indices ordered
## by resolution and lexicographically within each resolution
genInds = function(M, J=c(4)){

  genIndices = function(mLeft, J=c(4)){
    Jm = J[1]
    if(mLeft==1){
      return(sapply(seq(Jm), function(num) as.character(num)))
    } else {
      if( length(J)>1 ) J=J[2:length(J)]
      children = genIndices(mLeft-1, J=J)
      indices = list()
      for(i in 1:Jm){
        indices[[2*i]] = lapply(children, function(num) paste(as.character(i), num, sep="_") )
        indices[[2*i-1]] = list(i)
      }
      return(do.call(c, unlist(indices, recursive=FALSE)))
    }
  }

  inds = genIndices(M,J)
  inds = sapply(inds, function(ind) paste("r", ind, sep="_"))

  lengths = sapply(inds, function(s) nchar(s))
  ord = order(lengths)
  return(inds[ord])
}



## plots ordered locations such that color intensity is decreasing with order
plot.locsord = function(locsord, col = "#000000", col2="#FFFFFF"){

  nlocs = length(locsord)/ncol(locsord)
  vals = seq(nlocs)
  #collist = grey.colors(nlocs, start=0.2, end=0.95)

  colf = function(f){
    color = rgb((1-f)*(t(col2rgb(col))) + f*(t(col2rgb(col2))), maxColorValue = 255)
    return(color)
  }

  collist = sapply(vals/nlocs, colf)
  fields::quilt.plot(locsord, vals, col=collist)
}
