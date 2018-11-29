plot.locs.tree = function(locs.tree, locs, knots=NULL){
  if( ncol(locs)==2 ) plot.tree.2d(locs.tree, locs, knots)
  if( ncol(locs)==1 ) plot.tree.1d(locs.tree, locs, knots)


}



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




plot.tree.1d = function(locs.tree, locs, knots) {

  colors = c('#e7298a','#377eb8','#4daf4a')#'#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
  nc = length(colors)

  m = 0
  M = get.M(locs.tree)
  oldpar = par(mfrow=c(M+1, 1))

  xlim = c(min(locs), max(locs))

  for( ind in names(locs.tree) ){
    node.points = locs.tree[[ind]]
    if( length(node.points)==0) next
    if( nchar(ind) > m){
      m=nchar(ind)
      plot(locs[node.points], rep(0.1, length(node.points)), col=colors[1], pch=15,
           ylim=c(0, .5), xlim=xlim, yaxt="n",
           xlab="", ylab="", main=paste(c("resolution=", m-1), collapse=""))
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
