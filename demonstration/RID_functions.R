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