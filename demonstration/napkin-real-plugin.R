source('RID_functions.R')
source('WERM_Heuristic.R')
source('DRModule.R')

PlugInEstimator = function(OBS,mismode,seednum){
  # Compute P(Y=1 | w,r,x)
  ExpYParam_Real = function(myallpossible,DATA,mylambda){
    yvalfix = 1 
    IyTrain = (DATA$Y == yvalfix)*1
    Xtrain = DATA$X
    if (mismode == 1){
      IyTrain = distortVar(IyTrain,seednum)
      Xtrain = distortVar(Xtrain,seednum)
    }
    modelY = learnXG(as.matrix(DATA[,c('W','R','X')]), IyTrain, mylambda, binommode = 1)
    evalMat = as.matrix(myallpossible[,c('W','R','X')])
    predval = predict(modelY,newdata=evalMat,type='response')
    # Add noise 
    myN = nrow(DATA)
    if (mismode == 0){
      cvgrate = 4 
      predval = predval + rnorm(n=nrow(evalMat), mean = myN^(-1/cvgrate), sd = myN^(-1/cvgrate))    
      predval = fix_pred(predval)
    }
    if (mismode == 2){
      cvgrate = 2 
      predval = predval + rnorm(n=nrow(evalMat), mean = myN^(-1/cvgrate), sd = myN^(-1/cvgrate))    
      predval = fix_pred(predval)
    }
    
    newcol = (ncol(myallpossible)+1)
    myallpossible[,newcol] = predval
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  # Compute P(x | w,r)
  ProbXParam_Real = function(myallpossible,DATA,mylambda){
    Xtrain = DATA$X
    if (mismode == 1){
      Xtrain = distortVar(Xtrain,seednum)
    }
    modelX = learnXG(as.matrix(DATA[,c('W','R')]),Xtrain,mylambda,binommode = 0)
    evalMat = as.matrix(myallpossible[,c('W','R')])
    predval = predict(modelX,newdata=evalMat,type='response')
    myN = nrow(DATA)*2
    if (mismode == 0){
      cvgrate =  3
      predval = predval + rnorm(n=nrow(evalMat), mean = myN^(-1/cvgrate), sd = myN^(-1/cvgrate))    
      predval = fix_pred(predval)
    }
    if (mismode == 2){
      cvgrate =  2
      predval = predval + rnorm(n=nrow(evalMat), mean = myN^(-1/cvgrate), sd = myN^(-1/cvgrate))    
      predval = fix_pred(predval)
    }
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
  ProbRXParam_Real = function(myallpossible,DATA,mylambda){
    Xtrain = DATA$X
    Rtrain = DATA$R 
    if (mismode == 1){
      Xtrain = distortVar(Xtrain,seednum)
    }
    if (mismode == 2){
      Rtrain = distortVar(Rtrain,seednum)
    }
    newcol = (ncol(myallpossible)+1)
    for (rval in Runique){
      for (xval in Xunique){
        # filtered_DATA = subset(DATA,(X==xval & R==rval))
        prob_xr = sum(((Xtrain == xval) * (Rtrain ==rval))*1)/nrow(DATA)
        # prob_xr = nrow(filtered_DATA)/nrow(DATA)
        myallpossible[myallpossible$X==xval & myallpossible$R==rval,newcol] = prob_xr
      }
    }
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  ProbWParam_Real = function(myallpossible,DATA,mylambda){
    # Compute P(w)
    for (wval in Wunique){
      filtered_DATA_W = subset(DATA,W==wval)
      probval.W = nrow(filtered_DATA_W)/nrow(DATA)   
      myallpossible[myallpossible$W == wval,'prob'] = probval.W
    }
    return(myallpossible)
  }
  
  W = OBS[,1] # High dim surrogate
  R = OBS[,2] # Cofounder 0-numCate
  X = OBS[,3]
  Y = OBS[,4]
  
  Wunique = unique(W)[order(unique(W))]
  Runique = unique(R)[order(unique(R))]
  Xunique = unique(X)[order(unique(X))]
  Yunique = unique(Y)[order(unique(Y))]
  
  # Setting
  DATA = data.frame(cbind(W,R,X,Y))
  DATA = subset(DATA,(is.na(W) == FALSE)&(is.na(R) == FALSE)&(is.na(X) == FALSE)&(is.na(Y) == FALSE))
  Ndata = nrow(DATA)
  
  tmp = GoodSplit(DATA)
  DATA_Train = tmp[[1]]
  DATA_Eval = tmp[[2]]
  
  mylambda = rep(100/sqrt(nrow(DATA)),nrow(DATA)/2)
  
  # Enumerate all possible values of column
  tmp = c()
  tmp = append(tmp,list(Wunique)) # W
  tmp = append(tmp, list(Runique)) # R
  tmp = append(tmp,list(Xunique)) # X
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('W','R','X')
  
  Ytable1 = ExpYParam_Real(allpossible,DATA_Train,mylambda)
  Ytable2 = RunTryCatchProb_plugin(ExpYParam_Real,allpossible,DATA_Train,DATA_Eval,mylambda)
  Ytable = (Ytable1 + Ytable2)/2
  
  # Compute P(x | r,w )
  PxTable1 = ProbXParam_Real(allpossible,DATA_Train,mylambda)
  PxTable2 = RunTryCatchProb_plugin(ProbXParam_Real,allpossible,DATA_Train,DATA_Eval,mylambda)
  PxTable = (PxTable1 + PxTable2)/2 
  
  # Compute P(r,x)
  PrxTable1 = ProbRXParam_Real(allpossible,DATA_Train,mylambda)
  PrxTable2 = RunTryCatchProb_plugin(ProbRXParam_Real,allpossible,DATA_Train,DATA_Eval,mylambda)
  PrxTable = (PrxTable1+PrxTable2)/2  
  
  # Compute P(w)
  PwTable1 = ProbWParam_Real(allpossible,DATA_Train,mylambda)
  PwTable2 = RunTryCatchProb_plugin(ProbWParam_Real,allpossible,DATA_Train,DATA_Eval,mylambda)
  PwTable = (PwTable1+PwTable2)/2
  
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

