source('RID_functions.R')
source('WERM_Heuristic.R')
PlugInEstimator = function(OBS,mismode,seednum){
  W = OBS[,1] # High dim surrogate
  R = OBS[,2] # Cofounder 0-numCate
  X = OBS[,3]
  Y = OBS[,4]
  
  Wunique = unique(W)[order(unique(W))]
  Runique = unique(R)[order(unique(R))]
  Xunique = unique(X)[order(unique(X))]
  Yunique = unique(Y)[order(unique(Y))]
  
  yvalfix = 1 
  IyTrain = (Y == yvalfix)*1
  Rtrain = R 
  Xtrain = X 
  if (mismode == 1){
    IyTrain = distortVar(IyTrain,seednum)
    Xtrain = distortVar(Xtrain,seednum)
  }
  if (mismode == 2){
    Rtrain = distortVar(Rtrain,seednum)
  }
  
  # Setting
  DATA = cbind(W,R,X,Y)
  DATA = subset(DATA,(is.na(W) == FALSE)&(is.na(R) == FALSE)&(is.na(X) == FALSE)&(is.na(Y) == FALSE))
  DATA = data.frame(DATA)
  
  # Enumerate all possible values of column
  tmp = c()
  tmp = append(tmp,list(Wunique)) # W
  tmp = append(tmp, list(Runique)) # R
  tmp = append(tmp,list(Xunique)) # X
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('W','R','X')
  
  # Compute P(Y=1 | w,r,x)
  ExpYParam_Real = function(myallpossible,DATA,yval){
    modelY = learnXG(as.matrix(DATA[,c('W','R','X')]), IyTrain, rep(0,length(IyTrain)),binommode = 1)
    evalMat = as.matrix(myallpossible[,c('W','R','X')])
    predval = predict(modelY,newdata=evalMat,type='response')
    newcol = (ncol(myallpossible)+1)
    myallpossible[,newcol] = predval
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  # Compute P(x | w,r)
  ProbXParam_Real = function(myallpossible,DATA){
    modelX = learnXG(as.matrix(DATA[,c('W','R')]),Xtrain,rep(0,length(X)),binommode = 0)
    evalMat = as.matrix(myallpossible[,c('W','R')])
    predval = predict(modelX,newdata=evalMat,type='response')
    predval = t(matrix(predval,nrow=3))
    probX = rep(0,nrow(myallpossible))
    for (idx in 1:nrow(myallpossible)){
      xval = myallpossible$X[idx]
      probX[idx] = predval[idx,(xval+1)]
    }
    newcol = (ncol(myallpossible)+1)
    myallpossible[,newcol] = probX
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  # Compute P(r,x)
  ProbRXParam_Real = function(myallpossible,DATA){
    newcol = (ncol(myallpossible)+1)
    for (rval in Runique){
      for (xval in Xunique){
        filtered_DATA = subset(DATA,(X==xval & R==rval))
        prob_xr = nrow(filtered_DATA)/nrow(DATA)
        myallpossible[myallpossible$X==xval & myallpossible$R==rval,newcol] = prob_xr
      }
    }
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  
  Ytable = allpossible
  Ytable = ExpYParam_Real(Ytable,DATA,1)
  
  # Compute P(x | r,w )
  PxTable = allpossible 
  PxTable = ProbXParam_Real(PxTable,DATA)
  
  # Compute P(r,x)
  PrxTable = allpossible 
  PrxTable = ProbRXParam_Real(PrxTable,DATA)
  
  # Compute P(w)
  PwTable = allpossible  

  for (wval in Wunique){
    filtered_DATA_W = subset(DATA,W==wval)
    probval.W = nrow(filtered_DATA_W)/nrow(DATA)   
    PwTable[PwTable$W == wval,'prob'] = probval.W
  }
  allpossibleOrig = allpossible
  
  ComputeVal = allpossible
  ComputeVal$val1 = Ytable$prob * PxTable$prob * PwTable$prob
  ComputeVal$val2 = PxTable$prob * PwTable$prob
  
  ## Marginalizing over W
  tmp = c()
  tmp = append(tmp, list(Runique)) # R
  tmp = append(tmp,list(Xunique)) # X
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('R','X')
  ComputeVal.MarginW = allpossible
  ComputeVal.MarginW$val1 = 0
  ComputeVal.MarginW$val2 = 0
  
  for (rval in Runique){
    for(xval in Xunique){
      ComputeVal.MarginW[ComputeVal.MarginW$R==rval & ComputeVal.MarginW$X==xval,'val1'] = sum(ComputeVal[ComputeVal$X==xval & ComputeVal$R==rval,'val1'],na.rm=T)
      ComputeVal.MarginW[ComputeVal.MarginW$R==rval & ComputeVal.MarginW$X==xval,'val2'] = sum(ComputeVal[ComputeVal$X==xval & ComputeVal$R==rval,'val2'],na.rm=T)
    }
  }
  ComputeVal.MarginW$val3 = exp(log(ComputeVal.MarginW$val1) - log(ComputeVal.MarginW$val2))
  
  rList = rep(0,length(Runique))
  idx = 1 
  for (rval in Runique){
    rList[idx] = sum(PrxTable[PrxTable$W==0 & PrxTable$R==rval,'prob'])
    idx = idx + 1 
  }
  RFix = Runique[which.max(rList)]
  answer = ComputeVal.MarginW[ComputeVal.MarginW$R==RFix,'val3']
  return(answer)
}

