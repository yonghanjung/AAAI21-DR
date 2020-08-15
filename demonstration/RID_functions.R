fix_pred = function(pred){
  smallval = 1e-2
  pred[pred<0] = smallval; pred[pred>1] = 1-smallval
  return(pred)
}

mis_pred = function(predval,distortval){
  predval = predval + rnorm(length(predval),distortval,0.1)
  predval = fix_pred(predval)
  return(predval)
}

distortVar = function(myvar,seednum){
  set.seed(seednum)
  while(1){
    mysample = sample(c(0:(length(unique(myvar))-1)),size=length(myvar),replace=T)  
    if (length(unique(mysample)) == length(unique(myvar)) ){
      return(mysample)
    }
  }
}