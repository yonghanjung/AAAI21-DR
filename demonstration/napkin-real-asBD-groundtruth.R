source('RID_functions.R')
source('WERM_Heuristic.R')

BDEstimator = function(DATA){
  ####################################################
  # Main 
  ####################################################
  K = DATA[,1]
  R = DATA[,3]
  TPR = DATA[,6]
  X = DATA[,4]
  Y = DATA[,5]
  
  Xunique = unique(X)[order(unique(X))]
  
  Iy = (Y==1)*1
  
  # model.Y.1 = learnXG(inVar=data.matrix(data.frame(TPR=TPR, X=X)),labelval = Iy,binommode = 1,regval = rep(0,nrow(DATA)))
  model.Y = learnXG(inVar=data.matrix(data.frame(X=X, K=K, R=R)),labelval = Iy,binommode = 1,regval = rep(0,nrow(DATA)))
  
  Yx = rep(0,length(Xunique))
  idx = 1
  for (xval in Xunique){
    DATA.Xval = data.frame(X=rep(xval,nrow(DATA)),K=K,R=R)
    # DATA.Xval = data.frame(cbind(W=W,R=R,X=rep(xval,nrow(DATA))))  
    pred.Yx = predict(model.Y,newdata=data.matrix(DATA.Xval),type='predict')
    Yx[idx] = mean(pred.Yx,na.rm=T)
    idx = idx + 1 
  }
  return(Yx)
  
  # W = DATA[,1] # High dim surrogate
  # R = DATA[,2] # Cofounder 0-numCate
  # X = DATA[,3]
  # Y = DATA[,4]
  # 
  # Wunique = unique(W)[order(unique(W))]
  # Runique = unique(R)[order(unique(R))]
  
  # Yunique = unique(Y)[order(unique(Y))]
  # 
  # DATA = data.frame(cbind(W,R,X,Y))
  # DATA = subset(DATA,(is.na(W) == FALSE)&(is.na(R) == FALSE)&(is.na(X) == FALSE)&(is.na(Y) == FALSE))
  # Ndata = nrow(DATA)
}




