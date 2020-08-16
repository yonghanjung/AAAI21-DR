source('RID_functions.R')
source('WERM_Heuristic.R')

returnUniqueFromVector = function(myvector){
  return(unique(myvector)[order(unique(myvector))])
}
myDivisionTruth = function(a,b){
  val = a/b
  if (is.nan(val)){
    val = 0
  }
  return(val)
}

BDEstimator = function(DATA){
  ####################################################
  # Main 
  ####################################################
  Ndata = nrow(DATA)
  
  PMB = DATA$PMB
  INT = DATA$INT
  KINK = DATA$KINK
  VTUB = DATA$VTUB
  SHNT = DATA$SHNT # X1
  VLNG = DATA$VLNG # X2
  Y = DATA$CCHL
  DATA = DATA[,c('PMB','INT','KINK','VTUB','SHNT','VLNG','CCHL')]
  
  PMBUnique = returnUniqueFromVector(PMB)
  INTUnique = returnUniqueFromVector(INT)
  KINKUnique = returnUniqueFromVector(KINK)
  VTUBUnique = returnUniqueFromVector(VTUB)
  SHNTUnique = returnUniqueFromVector(SHNT)
  VLNGUnique = returnUniqueFromVector(VLNG)
  YUnique = returnUniqueFromVector(Y)
  
  X1unique = SHNTUnique
  X2unique = VLNGUnique
  Iy = (Y==1)*1
  
  # model.Y.1 = learnXG(inVar=data.matrix(data.frame(TPR=TPR, X=X)),labelval = Iy,binommode = 1,regval = rep(0,nrow(DATA)))
  model.Y = learnXG(inVar=data.matrix(data.frame(PMB,INT,KINK,VTUB,SHNT,VLNG)),labelval = Iy,binommode = 1,regval = rep(0,nrow(DATA)))
  
  Yx = rep(0,length(X1unique)*length(X2unique))
  idx = 1
  for (x1val in X1unique){
    for (x2val in X2unique){
      DATA.Xval = data.frame(PMB,INT,KINK,VTUB,SHNT = rep(x1val,nrow(DATA)),VLNG = rep(x2val, nrow(DATA)))
      # DATA.Xval = data.frame(cbind(W=W,R=R,X=rep(xval,nrow(DATA))))  
      pred.Yx = predict(model.Y,newdata=data.matrix(DATA.Xval),type='predict')
      Yx[idx] = mean(pred.Yx,na.rm=T)
      idx = idx + 1   
    }
  }
  return(Yx)
}




