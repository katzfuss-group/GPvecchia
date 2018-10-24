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




plot.locs.tree = function(locs.tree, locs, knots=NULL){
  
  colors = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
  nc = length(colors)
  
  m = 0
  for( ind in names(locs.tree) ){
    node.points = locs.tree[[ind]]
    if(length(node.points)==0) next
    if(nchar(ind)>m){
      m=nchar(ind)
      plot(locs[node.points,1], locs[node.points,2], col=colors[1], pch=16, 
           xlim=c(0, 1), ylim=c(0, 1), xlab="", ylab="", 
           main=paste(c("resolution=", m-1), collapse=""))
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