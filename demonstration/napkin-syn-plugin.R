source('RID_functions.R')
source('WERM_Heuristic.R')
source('DRModule.R')

highdim_reg_xgboost = function(OBS, outVarVector){
  ### P(W), P(W|X)
  for (d in 1:D){
    outVar = data.matrix(outVarVector[,d])
    if (d > 1){
      # P(r)
      inVar = data.matrix(outVarVector[,c(1:(d-1))])
      MyModel = xgboost(verbose = 0, data = data.matrix(inVar), label = outVar, nrounds = 20,max.depth=10,lambda=0,alpha=0, objective = "binary:logistic")
      list.Pr = c(list.Pr, list(MyModel))
    }else{
      inVar = data.matrix(rep(1,nrow(OBS)))
      MyModel = xgboost(verbose = 0, data = data.matrix(inVar), label = outVar, nrounds = 20,max.depth=10,lambda=0,alpha=0, objective = "binary:logistic")
      list.Pr = list(MyModel)
    }
  }
  return(list.Pr)
}

PlugInEstimator = function(OBS,mydim,mismode,seednum){
  # Compute P(Y=1 | w,r,x)
  ExpYParam_Real = function(myallpossible,myDATA,mylambda){
    W = myDATA[,c(1:mydim)] # High dim surrogate
    R = myDATA[,(mydim+1)] # Cofounder 0-numCate
    X = myDATA[,(mydim+2)]
    Y = myDATA[,(mydim+3)]
    
    yvalfix = 1 
    IyTrain = (myDATA$Y == yvalfix)*1
    Xtrain = myDATA$X
    if (mismode == 1){
      IyTrain = distortVar(IyTrain,seednum)
      Xtrain = distortVar(Xtrain,seednum)
    }
    modelY = learnXG(inVar=as.matrix(data.frame(W,R,X)), labelval=IyTrain, regval=mylambda, binommode = 1)
    evalMat = as.matrix(myallpossible)
    predval = predict(modelY,newdata=evalMat,type='response')
    # Add noise 
    myN = nrow(myDATA)
    if (mismode == 0){
      cvgrate = 3 
      predval = predval + rnorm(n=nrow(evalMat), mean = myN^(-1/cvgrate), sd = myN^(-1/cvgrate))    
      predval = fix_pred(predval)
    }
    # if (mismode == 2){
    #   cvgrate = 2 
    #   predval = predval + rnorm(n=nrow(evalMat), mean = myN^(-1/cvgrate), sd = myN^(-1/cvgrate))    
    #   predval = fix_pred(predval)
    # }
    
    newcol = (ncol(myallpossible)+1)
    myallpossible[,newcol] = predval
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  # Compute P(x | w,r)
  ProbXParam_Real = function(myallpossible,myDATA,mylambda){
    W = myDATA[,c(1:mydim)] # High dim surrogate
    R = myDATA[,(mydim+1)] # Cofounder 0-numCate
    X = myDATA[,(mydim+2)]
    Y = myDATA[,(mydim+3)]
    
    Xtrain = myDATA$X
    if (mismode == 1){
      Xtrain = distortVar(Xtrain,seednum)
    }
    modelX = learnXG(inVar=as.matrix(data.frame(W,R)), labelval=Xtrain, regval=mylambda, binommode = 1)
    evalMat = as.matrix(myallpossible[,c(1:(mydim+1))])
    predval = predict(modelX,newdata=evalMat,type='response')
    myN = nrow(myDATA)*2
    if (mismode == 0){
      cvgrate =  3
      predval = predval + rnorm(n=nrow(evalMat), mean = myN^(-1/cvgrate), sd = myN^(-1/cvgrate))    
      predval = fix_pred(predval)
    }
    # if (mismode == 2){
    #   cvgrate =  2
    #   predval = predval + rnorm(n=nrow(evalMat), mean = myN^(-1/cvgrate), sd = myN^(-1/cvgrate))    
    #   predval = fix_pred(predval)
    # }
    # predval = t(matrix(predval,nrow=3))
    probX = myallpossible$X*predval + (1-myallpossible$X)*(1-predval)
    # for (idx in 1:nrow(myallpossible)){
    #   xval = myallpossible$X[idx]
    #   probX[idx] = predval[idx,(xval+1)]
    # }
    newcol = (ncol(myallpossible)+1)
    myallpossible[,newcol] = probX
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  # Compute P(r,x)
  ProbRXParam_Real = function(myallpossible,myDATA,mylambda){
    W = myDATA[,c(1:mydim)] # High dim surrogate
    R = myDATA[,(mydim+1)] # Cofounder 0-numCate
    X = myDATA[,(mydim+2)]
    Y = myDATA[,(mydim+3)]
    
    Xtrain = myDATA$X
    Rtrain = myDATA$R 
    if (mismode == 1){
      Xtrain = distortVar(Xtrain,seednum)
    }
    if (mismode == 2){
      Rtrain = distortVar(Rtrain,seednum)
    }
    newcol = (ncol(myallpossible)+1)
    for (rval in Runique){
      for (xval in Xunique){
        # filtered_myDATA = subset(myDATA,(X==xval & R==rval))
        prob_xr = sum(((Xtrain == xval) * (Rtrain ==rval))*1)/nrow(myDATA)
        # prob_xr = nrow(filtered_myDATA)/nrow(myDATA)
        myallpossible[myallpossible$X==xval & myallpossible$R==rval,newcol] = prob_xr
      }
    }
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  ProbWParam_Real = function(myallpossible,myDATA,mylambda){
    W = myDATA[,c(1:mydim)] # High dim surrogate
    R = myDATA[,(mydim+1)] # Cofounder 0-numCate
    X = myDATA[,(mydim+2)]
    Y = myDATA[,(mydim+3)]
    
    ################################################################################
    # Learn P(w)
    ################################################################################
    model.W = highdim_reg_xgboost(myDATA,W)
    tmp = rep(1,nrow(allpossible))
    for (d in 1:D){
      if (d == 1){
        predInVar = data.matrix(rep(1,nrow(allpossible)))
      }else{
        predInVar = data.matrix(allpossible[,c(1:(d-1))])
      }
      predval = predict(model.W[[d]],newdata=predInVar,type='response') # P(Zd =1 | Z(d-1),Z(d-2),...,Z(1),X)
      # predval[predval < 0] = 1e-8
      # predval[predval > 1] = 1 - 1e-8
      resultval = predval * allpossible[d] + (1-predval) * (1-allpossible[d])
      tmp = tmp * resultval
    }
    myallpossible[,'prob'] = tmp
    # 
    # Wtable = allpossible
    # Wtable[,(ncol(allpossible)+1)] = tmp
    # colnames(Wtable)[ncol(Wtable)] = 'prob'
    return(myallpossible)
  }
  
  W = OBS[,c(1:mydim)] # High dim surrogate
  R = OBS[,(mydim+1)] # Cofounder 0-numCate
  X = OBS[,(mydim+2)]
  Y = OBS[,(mydim+3)]
  
  Wunique = c(0,1)
  Runique = unique(R)[order(unique(R))]
  Xunique = unique(X)[order(unique(X))]
  Yunique = unique(Y)[order(unique(Y))]
  
  # Setting
  tmp = GoodSplit(OBS)
  DATA_Train = tmp[[1]]
  DATA_Eval = tmp[[2]]
  
  mylambda = rep(100/sqrt(nrow(OBS)),nrow(OBS)/2)
  
  ################################################################################ 
  # Enumerate all possible values of column
  ################################################################################
  tmp = c()
  for (d in 1:D){
    tmp = append(tmp,list(Wunique)) # W 
  }
  Wname = paste("W",1:D,sep="")
  tmp = append(tmp,list(Runique)) # R 
  tmp = append(tmp,list(Xunique)) # X 
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c(Wname,'R','X')
  ################################################################################
  
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
  Yx = ComputeVal.MarginW[ComputeVal.MarginW$R==RFix,'val3']
  return(Yx)
}

