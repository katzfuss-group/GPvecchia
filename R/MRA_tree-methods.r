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

get.rm = function(knot.tree){
  M = get.M(knot.tree)
  all.knots = names(knot.tree)
  rs = rep(0, M)
  for( ind in all.knots ){
    m = res(ind)
    if(m==M) next
    r = length(knot.tree[[ind]])
    #print(paste("ind: ", ind, ", r: ", r, sep=""))
    if(rs[m+1]==0) rs[m+1]=r
    else if(rs[m+1]!=r) {
      warning("different number of knots for subregions of the same resolution")
    }
  }
  return(rs)
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
