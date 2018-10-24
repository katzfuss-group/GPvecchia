rm(list=ls())
   
plot.locsord = function(locsord, col = "#000000", col2="#FFFFFF"){
  
  nlocs = length(locsord)/ncol(locsord)
  vals = seq(nlocs)
  #collist = grey.colors(nlocs, start=0.2, end=0.95)
  
  colf = function(f){
    color = rgb((1-f)*(t(col2rgb(col))) + f*(t(col2rgb(col2))), maxColorValue = 255)
    return(color)
  }
  
  collist = sapply(vals/nlocs, colf)
  print(collist)
  fields::quilt.plot(locsord, vals, col=collist)
}



#spatial.dim=2 # number of spatial dimensions
#n=21  # number of observed locs

# simulate locations
#set.seed(10)
#if(spatial.dim==1){
#  locs=matrix(runif(n),ncol=1)
#} else {
#  locs <- cbind(runif(n),runif(n))
#}

#plot.locsord(locs, col="#FF0000")