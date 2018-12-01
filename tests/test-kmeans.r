rm(list=ls())
#set.seed(1990)


plot.cl = function(clusters, centers=NULL, title=""){
  colors = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
  if(!is.null(centers)) {
    plot(centers[,1], centers[,2], col=colors[1:K], pch=as.character(seq(K)),
         xlim=c(0, 1), ylim=c(0,1), xlab="x", ylab="y", main=title)
  } else {
    plot(0, xlim=c(0, 1), ylim=c(0,1), xlab="x", ylab="y", main=title, pch="")
  }
  K = max(clusters)
  for(k in 1:K){
    cl = which(clusters==k)
    color.id = k %% length(colors)+1
    points(locs[cl,1], locs[cl,2], col=colors[color.id], pch=16)
  }
}


halve = function(splits){
  new.splits = c(splits[1])
  for(i in 2:length(splits)){
    new.splits = c(new.splits, (splits[i] - splits[i-1])/2, splits[i])
  }
  return(new.splits)
}


equal.kmeans3 = function(locs, size){

  K = n/size

  J = 2**ceiling(log(K,2))
  if(!K==J){
    warning(paste("the number of subregions fo the original domain is ", K, " but has to be a power of 2. Setting J=", J, sep=""))
  }

  regions = list(seq(1:nrow(locs)))
  for(power in 1:log(J,2)){
    new.regions = vector("list", 2**power)
    for( reg.id in 1:length(regions)) {
      region = regions[[reg.id]]
      if(ncol(locs)==2 && (power %% 2)==1) d=1 else d=2
      locs.in.region = matrix(locs[region,], ncol=ncol(locs))
      cutoff = quantile(locs.in.region[,d], 0.5)
      new.regions[[2*reg.id-1]] = region[which(locs.in.region[,d] <= cutoff)]
      new.regions[[2*reg.id]] = region[which(locs.in.region[,d] > cutoff)]
    }
    regions = new.regions
  }

  ## return data in a format consistent with kmeans clustering
  clusters = rep(0, nrow(locs))
  id = 1
  for(region in regions){
    clusters[region] = id
    id = id + 1
  }
  if( any(table(clusters)>size))  stop("something went wrong with clustering. Some clusters are too big")
  return(clusters)
}


n=200  # number of observed locs
size=7

# simulate locations
locs = cbind(runif(n),runif(n))



#clusters = equal.kmeans(locs, K, algo='marcin')
#clusters = equal.kmeans(locs, K, algo='kidzinski')
#clusters = equal.kmeans2(locs, K)
clusters = equal.kmeans3(locs, size)
plot.cl(clusters)
