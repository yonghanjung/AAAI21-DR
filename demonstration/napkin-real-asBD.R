source('RID_functions.R')
source('WERM_Heuristic.R')

asBDEstimator = function(OBS,mismode,seednum){
  ####################################################
  # Main 
  ####################################################
  W = OBS[,1] # High dim surrogate
  R = OBS[,2] # Cofounder 0-numCate
  X = OBS[,3]
  Y = OBS[,4]
  
  Wunique = unique(W)[order(unique(W))]
  Runique = unique(R)[order(unique(R))]
  Xunique = unique(X)[order(unique(X))]
  Yunique = unique(Y)[order(unique(Y))]

  DATA = data.frame(cbind(W,R,X,Y))
  DATA = subset(DATA,(is.na(W) == FALSE)&(is.na(R) == FALSE)&(is.na(X) == FALSE)&(is.na(Y) == FALSE))
  Ndata = nrow(DATA)
  
  Iy = (Y==1)*1
  
  model.Y = learnXG(inVar=data.matrix(data.frame(W=W,R=R,X=X)),labelval = Iy,binommode = 1,regval = rep(0,nrow(DATA)))
  
  Yx = rep(0,length(Xunique))
  idx = 1
  for (xval in Xunique){
    DATA.Xval = data.frame(cbind(W=W,R=R,X=rep(xval,nrow(OBS))))  
    pred.Yx = predict(model.Y,newdata=data.matrix(DATA.Xval),type='predict')
    Yx[idx] = mean(pred.Yx,na.rm=T)
    idx = idx + 1 
  }
  return(Yx)
  
}




