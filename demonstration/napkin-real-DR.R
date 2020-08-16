source('RID_functions.R')
source('WERM_Heuristic.R')

returnUnique = function(OBS){
  W = OBS[,1] # High dim surrogate
  R = OBS[,2] # Cofounder 0-numCate
  X = OBS[,3]
  Y = OBS[,4]
  Wunique = unique(W)[order(unique(W))]
  Runique = unique(R)[order(unique(R))]
  Xunique = unique(X)[order(unique(X))]
  Yunique = unique(Y)[order(unique(Y))]
  return(list(Wunique,Runique,Xunique,Yunique))
}

nonRandomSampleSplit = function(OBS,mysize){
  OBSuniqueList = returnUnique(OBS)
  # Enumerate all possible values of column
  tmp = c()
  for (idx in 1:length(OBSuniqueList)){
    tmp = append(tmp,list(OBSuniqueList[[idx]]))
  }
  allpossible = expand.grid(tmp)
  colnames(allpossible) = colnames(OBS)
  mycollect = c()
  for (rowidx in 1:nrow(allpossible)){
    filtered_OBS = subset(OBS, W == allpossible[rowidx,'W'] & R == allpossible[rowidx,'R'] & X == allpossible[rowidx,'X'] & Y == allpossible[rowidx,'Y'] )
    if (nrow(filtered_OBS) > 0){
      mycollect = rbind(mycollect,filtered_OBS[1,]) 
    }
  }
  if (nrow(mycollect) < mysize){
    for (rowidx in 1:nrow(OBS)){
      if (rowidx %in% as.numeric(rownames(mycollect)) == FALSE){
        mycollect = rbind(mycollect,OBS[rowidx,])  
      }
      if (nrow(mycollect) == mysize){
        break 
      }
    }
  }
  
  totalidx = c(1:nrow(OBS))
  mycollect2 = setdiff(totalidx,as.numeric(rownames(mycollect)))
  OBS_2 = OBS[mycollect2,]
  OBS_2_unique = returnUnique(OBS_2)
  OBS_1_unique = returnUnique(OBS_1)
  if (identical( OBS_2_unique[[length(OBS_2_unique)]], OBS_1_unique[[length(OBS_1_unique)]] ) == FALSE){
    while(1){
      splitidx = sample(c(1:nrow(OBS)),size=nrow(OBS)/2)
      OBS_2 = OBS[splitidx,]
      if (identical(OBS_2_unique[[length(OBS_2_unique)]], OBS_1_unique[[length(OBS_1_unique)]] )){
        break
      }
    }  
  }
  OBS_1 = mycollect
  rownames(OBS_1) = c(1:nrow(OBS_1))
  rownames(OBS_2) = c(1:nrow(OBS_2))
  return(list(OBS_1,OBS_2))
}

GoodSplit = function(OBS){
  totalidx = c(1:nrow(OBS))
  
  iteridx = 0
  iterMax = 10
  
  while(1){
    iteridx = iteridx + 1 
    
    splitidx_1 = sample(c(1:nrow(OBS)),size=nrow(OBS)/2)
    splitidx_2 = setdiff(totalidx,splitidx_1)
    OBS_1 = OBS[splitidx_1,]
    OBS_2 = OBS[splitidx_2,]
    
    Unique_1 = returnUnique(OBS_1)
    Unique_2 = returnUnique(OBS_2)
    stopSwitch = TRUE 
    if (identical(Unique_1,Unique_2) == FALSE){
      stopSwitch = FALSE 
    }
    # for (idx in 1:4){
    #   
    # }
    if (stopSwitch == TRUE){
      break
    }
    if (iteridx > iterMax){
      mytmp = nonRandomSampleSplit(OBS,mysize=nrow(OBS)/2)
      OBS_1 = mytmp[[1]]
      OBS_2 = mytmp[[2]]
    }
  }
  rownames(OBS_1) = c(1:nrow(OBS_1))
  rownames(OBS_2) = c(1:nrow(OBS_2))
  return(list(OBS_1,OBS_2))
}

returnRegularizer = function(myDATA,mismode){
  myregval = 100 
  
  if (mismode == 0){
    mylambda.Y = rep(myregval/sqrt(sqrt(nrow(myDATA))), nrow(myDATA))
    mylambda.X = rep(myregval/sqrt(sqrt(nrow(myDATA))), nrow(myDATA))
    mylambda.XY = rep(myregval/sqrt(sqrt(nrow(myDATA))), nrow(myDATA))
    mylambda.R = rep(myregval/sqrt(sqrt(nrow(myDATA))), nrow(myDATA))
  }
  if (mismode == 1){
    # mylambda.Y = rep(myregval/sqrt(sqrt(nrow(myDATA))), nrow(myDATA))
    # mylambda.X = rep(myregval/sqrt(sqrt(nrow(myDATA))), nrow(myDATA))
    # mylambda.XY = rep(myregval/sqrt(sqrt(nrow(myDATA))), nrow(myDATA))
    mylambda.Y = rep(0, nrow(myDATA))
    mylambda.X = rep(0, nrow(myDATA))
    mylambda.XY = rep(0, nrow(myDATA))
    mylambda.R = rep(0, nrow(myDATA))
  }
  if (mismode == 2){
    mylambda.Y = rep(0, nrow(myDATA))
    mylambda.X = rep(0, nrow(myDATA))
    mylambda.XY = rep(0, nrow(myDATA))
    # mylambda.R = rep(myregval/sqrt(sqrt(nrow(myDATA))), nrow(myDATA))
    mylambda.R = rep(0, nrow(myDATA))
  }
  return(list(mylambda.Y,mylambda.X,mylambda.XY,mylambda.R))
}

DREstimator = function(OBS,mismode,seednum){
  ####################################################
  # TrainModel 
  ####################################################
  TrainModel = function(DATA_Train, DATA_Eval, DATA, xfix, yfix, mismode){
    IyTrain = (DATA_Train$Y == yfix)*1
    Rtrain = DATA_Train$R 
    Xtrain = DATA_Train$X 
    
    if (mismode == 1){
      IyTrain = distortVar(IyTrain,seednum)
      Xtrain = distortVar(Xtrain,seednum)
    }
    if (mismode == 2){
      Rtrain = distortVar(Rtrain,seednum)
    }
    IxTrain = (Xtrain == xfix)*1
    IxyTrain = IyTrain*IxTrain
    
    mylambda = returnRegularizer(DATA,mismode)
    mylambda.Y = mylambda[[1]] 
    mylambda.X = mylambda[[2]] 
    mylambda.XY = mylambda[[3]] 
    mylambda.R = mylambda[[4]] 
    
    # regvallist = seq(0,10,by=0.2)
    # lambda.XY = learnHyperParam(regvallist=regvallist, invar=data.matrix(data.frame(R=DATA_Train$R, W=DATA_Train$W)), mylabel=IxyTrain, learningbinary=1, TFcontinuous=F)
    lambda.XY = mylambda.XY
    model.xy.RW = learnXG(inVar = data.matrix(data.frame(R=DATA_Train$R, W=DATA_Train$W)),labelval = IxyTrain, regval = lambda.XY, binommode = 1)
    
    # lambda.X = learnHyperParam(regvallist=regvallist, invar=data.matrix(data.frame(R=DATA_Train$R, W=DATA_Train$W)), mylabel=IxTrain, learningbinary=1, TFcontinuous=F)
    lambda.X = mylambda.X
    model.x.RW = learnXG(inVar = data.matrix(data.frame(R=DATA_Train$R, W=DATA_Train$W)),labelval = IxTrain, regval = lambda.X, binommode = 1)
    
    # lambda.R = learnHyperParam(regvallist=regvallist, invar = data.matrix(data.frame(W=DATA_Train$W)), mylabel=IxTrain, learningbinary=0, TFcontinuous=F)
    lambda.R = mylambda.R
    model.R.W = learnXG(inVar = data.matrix(data.frame(W=DATA_Train$W)), labelval = Rtrain, regval = lambda.R, binommode = 0)
    
    return(list(model.xy.RW, model.x.RW, model.R.W))
  }
  
  ####################################################
  # Choose R  
  ####################################################
  ChooseR = function(DATA,xfix){
    # Choose the fixed R 
    RProb = rep(length(Runique))
    idx = 1
    for (rval in Runique){
      filtered_DATA = subset(DATA,R==rval & X==xfix)
      RProb[idx] = nrow(filtered_DATA)/nrow(DATA)
      idx = idx + 1 
    } 
    rfix = Runique[which.max(RProb)]
    return(rfix)
  }
  
  Compute_UIF_M1 = function(DATA_Train, DATA_Eval, trainedlist, yfix, xfix){
    ############################################
    # M1 =  M[x,y | r, W]
    ############################################
    model.xy.RW = trainedlist[[1]]; model.x.RW = trainedlist[[2]]; model.R.W = trainedlist[[3]]
    
    rfix = ChooseR(DATA_Eval,xfix)
    IrEval = (DATA_Eval$R == rfix)*1
    IxEval = (DATA_Eval$X == xfix)*1
    IyEval = (DATA_Eval$Y == yfix)*1
    IxyEval = IxEval * IyEval
    
    ### Evaluate P(x,y|R,W)
    prob.xy.RW = predict(model.xy.RW,newdata=data.matrix(data.frame(R=DATA_Eval$R,W=DATA_Eval$W)),type='response')
    prob.xy.rW = predict(model.xy.RW,newdata=data.matrix(data.frame(R=rep(rfix,nrow(DATA_Eval)),W=DATA_Eval$W)),type='response')
    
    ### Evaluate P(R|W)
    pred.R.W = predict(model.R.W, newdata=data.matrix(data.frame(W=DATA_Eval$W)),type='response')
    pred.R.W  = t(matrix(pred.R.W,nrow=length(Runique)))
    prob.R.W = mapply(function(idx,riterval){
      return(pred.R.W[idx,(riterval+1)])
    },c(1:nrow(DATA_Eval)), DATA_Eval$R)
    
    UIF_M1 = (IrEval*(IxyEval - prob.xy.RW)/(prob.R.W)) + (prob.xy.rW) 
    return(UIF_M1)
  }
  
  Compute_UIF_M2 = function(DATA_Train, DATA_Eval, trainedlist, xfix){
    ############################################
    # M1 =  M[x,y | r, W]
    ############################################
    model.xy.RW = trainedlist[[1]]; model.x.RW = trainedlist[[2]]; model.R.W = trainedlist[[3]]
    
    rfix = ChooseR(DATA_Eval,xfix)
    IrEval = (DATA_Eval$R == rfix)*1
    IxEval = (DATA_Eval$X == xfix)*1

    ### Evaluate P(x|R,W)
    prob.x.RW = predict(model.x.RW,newdata=data.matrix(data.frame(R=DATA_Eval$R, W=DATA_Eval$W)),type='response')
    prob.x.rW = predict(model.x.RW,newdata=data.matrix(data.frame(R=rep(rfix,nrow(DATA_Eval)), W=DATA_Eval$W)),type='response')
    
    ### Evaluate P(R|W)
    pred.R.W = predict(model.R.W, newdata=data.matrix(data.frame(W=DATA_Eval$W)),type='response')
    pred.R.W  = t(matrix(pred.R.W,nrow=length(Runique)))
    prob.R.W = mapply(function(idx,riterval){
      return(pred.R.W[idx,(riterval+1)])
    },c(1:nrow(DATA_Eval)), DATA_Eval$R)
    
    UIF_M2 = (IrEval*(IxEval - prob.x.RW)/(prob.R.W)) + (prob.x.rW) 
    return(UIF_M2)
  }
  
  Compute_Yx = function(DATA_Train, DATA_Eval, trainedlist, yfix, xfix){
    UIF_M1 = Compute_UIF_M1(DATA_Train, DATA_Eval, trainedlist, yfix, xfix)
    UIF_M2 = Compute_UIF_M2(DATA_Train, DATA_Eval, trainedlist, xfix)
    
    prob.xy.dor = mean(UIF_M1, na.rm=T)
    prob.x.dor = mean(UIF_M2, na.rm=T)
    prob.y.dox = prob.xy.dor/prob.x.dor
    
    EIF_M1 = UIF_M1 - prob.xy.dor
    EIF_M2 = UIF_M2 - prob.x.dor
    UIF = (1/(prob.x.dor))*(UIF_M1 - EIF_M2*prob.y.dox)
    return(mean(UIF,na.rm=T))
  }
  
  ConductDoubleML = function(numIter,DATA,xfix,yfix,mismode){
    YxList = rep(0,numIter)
    tmp = GoodSplit(DATA)
    DATA_Train = tmp[[1]]
    DATA_Eval = tmp[[2]]
    for (iteridx in 1:numIter){
      trainedlist1 = TrainModel(DATA_Train, DATA_Eval, DATA, xfix, yfix, mismode)
      trainedlist2 = TrainModel(DATA_Eval, DATA_Train, DATA, xfix, yfix, mismode)
      Yxval = mean(Compute_Yx(DATA_Train, DATA_Eval, trainedlist1, yfix, xfix),
                   Compute_Yx(DATA_Eval,DATA_Train, trainedlist2, yfix, xfix),na.rm=T)
      YxList[iteridx] = Yxval
    }
    return(mean(YxList,na.rm=T))
  }
  
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
  
  tmp = GoodSplit(DATA)
  DATA_Train = tmp[[1]]
  DATA_Eval = tmp[[2]]
  
  yfix = 1 
  YxDR = rep(0,length(Xunique))
  idx = 1 
  for (xfix in Xunique){
    # YxDR[idx] = ConductDoubleML(numIter=1,DATA=DATA,xfix=xfix,yfix=yfix,mismode=mismode)
    trainedlist1 = TrainModel(DATA_Train, DATA_Eval, DATA, xfix, yfix, mismode)
    trainedlist2 = TrainModel(DATA_Eval, DATA_Train , DATA, xfix, yfix, mismode)
    # trainedlist3 = TrainModel(DATA, DATA , DATA, xfix, yfix, mismode)
    YxDR[idx] = mean(Compute_Yx(DATA_Train, DATA_Eval, trainedlist1, yfix, xfix),
                     Compute_Yx(DATA_Eval, DATA_Train , trainedlist2, yfix, xfix),
                     na.rm=T)
    # YxDR[idx] = mean(Compute_Yx(DATA, DATA, trainedlist3, yfix, xfix),
    #                  Compute_Yx(DATA, DATA , trainedlist3, yfix, xfix), 
    #                  na.rm=T)
    YxDR[idx] = max(YxDR[idx],0)
    YxDR[idx] = min(YxDR[idx],1)
    idx = idx + 1 
  }
  return(YxDR)
}




