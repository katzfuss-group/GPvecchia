# specify the vecchia approximation, prepare U
# this fct does not depend on parameter values
# only has to be run once before repeated likelihood evals

vecchia_specify=function(z,locs,m,ordering,cond.yz,locs.pred,ordering.pred,pred.cond) {
  ###  arguments:
  # locs: nxd matrix of obs locs
  # ordering: options are 'coord' or 'maxmin'
  # cond.yz: options are 'y', 'z', 'SGV', 'SGVT', and 'zy'
  # ordering.pred: options are 'obspred' or 'general'
  # pred.cond: prediction conditioning, options are 'general' or 'independent'
  
  spatial.dim=ncol(locs)
  n=nrow(locs)
  
  # default options
  if(missing(ordering)){ ordering = (if(spatial.dim==1) 'coord' else ord='maxmin') }
  if(missing(cond.yz)){ cond.yz = (if(missing(locs.pred)) 'SGV' else 'SGVT') }
  if(missing(pred.cond)){ pred.cond='general' }
  

  ### order locs and z

  if(missing(locs.pred)){  # no prediction
    
    if(ordering=='coord') ord=order_coordinate(locs) else {
      ord = order_maxmin(locs)
    }    
    zord=z[ord]
    locsord=locs[ord,,drop=FALSE]
    obs=rep(TRUE,n)
    ordering.pred='general'
    
  } else {    # prediction is desired
    
    n.p=nrow(locs.pred)
    locs.all=rbind(locs,locs.pred)
    observed.obspred=c(rep(TRUE,n),rep(FALSE,n.p))
    if(missing(ordering.pred))
      if(spatial.dim==1 & ordering=='coord') ordering.pred='general' else ordering.pred='obspred'  
    if(ordering.pred=='general'){
      if(ordering=='coord') ord=order_coordinate(locs.all) else {
        ord = order_maxmin(locs.all)
      }    
      ord.obs=ord[ord<=n]
    } else {
      if(ordering=='coord') {
        ord.obs=order_coordinate(locs)
        ord.pred=order_coordinate(locs.pred)
      } else {
        ord.obs = order_maxmin(locs)
        ord.pred = order_maxmin(locs.pred)
      }
      ord=c(ord.obs,ord.pred+n)
    }
    zord=z[ord.obs]
    locsord=locs.all[ord,,drop=FALSE]
    obs=observed.obspred[ord]
  }

  
  ### obtain nearest neighbors
  if(spatial.dim==1) {
    NNarray=findOrderedNN_kdtree2(locsord,m)
  } else NNarray <- find_ordered_nn(locsord,m)
  
  if(pred.cond=='independent'){
    if(ordering.pred=='obspred'){
      NNarray.pred <- array(dim=c(n.p,m+1))
      dist.predobs=rdist(locsord[n+(1:n.p),,drop=FALSE],locsord[1:n,,drop=FALSE])
      for(j in 1:n.p){
        dists=dist.predobs[j,]
        m.nearest.obs=sort(order(dists)[1:m],decreasing=TRUE)
        NNarray.pred[j,]=c(j+n,m.nearest.obs)
      }
      NNarray[n+(1:n.p),]=NNarray.pred
    } else print('indep. conditioning currently only implemented for obspred ordering')
  }
  
  
  ### conditioning on y or z?
  if(cond.yz=='SGV'){
    Cond=whichCondOnLatent(NNarray,firstind.pred=n+1)
  } else if(cond.yz=='SGVT'){
    Cond=rbind(whichCondOnLatent(NNarray[1:n,]),matrix(TRUE,nrow=n.p,ncol=m+1))
  } else if(cond.yz=='y'){ 
    Cond=matrix(NA,n+n.p,m+1); Cond[!is.na(NNarray)]=TRUE
  } else if(cond.yz=='zy'){ 
    Cond=(NNarray>n); Cond[,1]=TRUE  
  } else {  # cond.yz=='z'
    Cond=matrix(NA,n+n.p,m+1); Cond[!is.na(NNarray)]=FALSE; Cond[,1]=TRUE
  }
  
  
  ### determine the sparsity structure of U
  U.prep=U_sparsity( locsord, NNarray, obs, Cond )
  
  
  ### object that specifies the vecchia approximation
  vecchia.approx=list(zord=zord,locsord=locsord,obs=obs,ord=ord,ord.pred=ordering.pred,U.prep=U.prep)
  # U.prep has attributes: revNNarray,revCond,n.cores,size,rowpointers,colindices,y.ind)
  
  return(vecchia.approx)
  
}
