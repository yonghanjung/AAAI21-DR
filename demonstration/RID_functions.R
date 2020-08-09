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
  if (length(unique(myvar)) == 2){
    mydistort = xor((myvar * rbinom(n=length(myvar),size=1,prob=0.5)),rbinom(n=length(myvar),size=1,prob=0.5))*1
  }else{
    myvar.unique = unique(myvar)[order(unique(myvar))]
    tmp = rep(sample(myvar.unique),floor(length(myvar)/length(myvar.unique)))
    if (length(tmp) < length(myvar)){
      tmp = c(tmp,rep(0,length(myvar)-length(tmp)))
    }
    mydistort = tmp 
  }
  return(mydistort)
}