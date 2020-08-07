source('RID_functions.R')
source('WERM_Heuristic.R')

returnUnique = function(OBS){
  X1 = OBS[,1] 
  Z = OBS[,2] 
  R = OBS[,3] 
  X2 = OBS[,4]  
  Y = OBS[,5]
  X1unique = unique(X1)[order(unique(X1))]
  Zunique = unique(Z)[order(unique(Z))]
  Runique = unique(R)[order(unique(R))]
  X2unique = unique(X2)[order(unique(X2))]
  Yunique = unique(Y)[order(unique(Y))]
  return(list(X1unique,Zunique,Runique,X2unique,Yunique))
}

GoodSplit = function(OBS){
  totalidx = c(1:nrow(OBS))
  
  while(1){
    splitidx_1 = sample(c(1:nrow(OBS)),size=nrow(OBS)/2)
    splitidx_2 = setdiff(totalidx,splitidx_1)
    OBS_1 = OBS[splitidx_1,]
    OBS_2 = OBS[splitidx_2,]
    
    Unique_1 = returnUnique(OBS_1)
    Unique_2 = returnUnique(OBS_2)
    stopSwitch = TRUE 
    for (idx in 1:ncol(OBS)){
      if (identical(Unique_1,Unique_2) == FALSE){
        stopSwitch = FALSE 
      }
    }
    if (stopSwitch == TRUE){
      break
    }
  }
  rownames(OBS_1) = c(1:nrow(OBS_1))
  rownames(OBS_2) = c(1:nrow(OBS_2))
  return(list(OBS_1,OBS_2))
}


DREstimator = function(OBS,distortval,mismode,seednum){
  set.seed(seednum)
  X1 = OBS[,1] 
  Z = OBS[,2] 
  R = OBS[,3] 
  X2 = OBS[,4]  
  Y = OBS[,5]
  
  X1unique = unique(X1)[order(unique(X1))]
  Zunique = unique(Z)[order(unique(Z))]
  Runique = unique(R)[order(unique(R))]
  X2unique = unique(X2)[order(unique(X2))]
  Yunique = unique(Y)[order(unique(Y))]
  
  yvalfix = 1 
  Iy = (Y==yvalfix)*1
  
  ############################################
  # Learn Models 
  ############################################
  TrainModel = function(DATA_Train, DATA_Eval, DATA){
    
    Iy.Train = (DATA_Train$Y==yvalfix)*1
    Iy.Test = (DATA_Eval$Y==yvalfix)*1
    
    regvallist = seq(0,10,by=0.2)
    lambda.Y = learnHyperParam(regvallist=regvallist, invar=data.matrix(data.frame(X1=DATA_Eval$X1, Z=DATA_Eval$Z, R=DATA_Eval$R, X2=DATA_Eval$X2)), mylabel=Iy.Test, learningbinary=1)
    # lambda.Y = rep(0,nrow(DATA_Train))
    model.Y = learnXG(inVar = data.matrix(data.frame(X1=DATA_Train$X1, Z=DATA_Train$Z, R=DATA_Train$R, X2=DATA_Train$X2)),labelval = Iy.Train, regval = lambda.Y,binommode = 1)
    
    # lambda.X2 = rep(1,nrow(DATA_Train))
    lambda.X2 = learnHyperParam(regvallist=regvallist, invar=data.matrix(data.frame(X1=DATA_Eval$X1, Z=DATA_Eval$Z)), mylabel=DATA_Eval$X2, learningbinary=0)
    model.X2.ZX1 = learnXG(inVar = data.matrix(data.frame(X1=DATA_Train$X1,Z=DATA_Train$Z)),labelval = DATA_Train$X2, regval = lambda.X2, binommode = 0)
  
    lambda.R = learnHyperParam(regvallist=regvallist, invar=data.matrix(data.frame(X1=DATA_Eval$X1)), mylabel=DATA_Eval$R, learningbinary=0)
    model.R.X1 = learnXG(inVar = data.matrix(data.frame(X1=DATA_Train$X1)),labelval = DATA_Train$R, regval = lambda.R,binommode = 0)
    
    return(list(model.Y,model.X2.ZX1,model.R.X1))
  }
  
  
  compute_UIF_M1 = function(DATA_Train,DATA_Eval,trainedlist,yval,rvalfix,x2valfix){
    ############################################
    # M1 =  M[y | (r,x2);(x1;z)]
    ############################################
    # rvalfix = 0 
    # x2valfix = 0
    Ir = (DATA_Eval$R==rvalfix)*1
    Ix2 = (DATA_Eval$X2==x2valfix)*1
    
    # mylist = TrainModel(DATA_Train,DATA_Eval)
    model.Y = trainedlist[[1]]; model.X2.ZX1 = trainedlist[[2]]; model.R.X1 = trainedlist[[3]]
    
    # Learn P(y | r,x2,X1,Z)
    pred.Y.rx2.X1Z = predict(model.Y,newdata=data.matrix(data.frame(X1=DATA_Eval$X1, Z=DATA_Eval$Z, R=rep(rvalfix,nrow(DATA_Eval)),X2=rep(x2valfix,nrow(DATA_Eval)))),type='response')
    # if (mismode == 1){
    #   pred.Y.rx2.X1Z = fix_pred(mis_pred(pred.Y.rx2.X1Z,distortval))
    # }
    # Learn P(y | R,X2,X1,Z)
    pred.Y.RX2.X1Z = predict(model.Y,newdata=data.matrix(data.frame(X1=DATA_Eval$X1, Z=DATA_Eval$Z, R=DATA_Eval$R, X2=DATA_Eval$X2)),type='response')
    # if (mismode == 1){
    #   pred.Y.RX2.X1Z = fix_pred(mis_pred(pred.Y.RX2.X1Z,distortval))
    # }
    # Learn P(R | X1)
    pred.R.X1 = predict(model.R.X1,newdata=data.matrix(data.frame(X1=DATA_Eval$X1)),reshape=TRUE)
    prob.R.X1 = mapply(function(idx,rval){
      pred.R.X1[idx,(rval+1)]
    },c(1:nrow(DATA_Eval)),DATA_Eval$R)
    # if (mismode == 1){
    #   prob.R.X1 = fix_pred(mis_pred(prob.R.X1,distortval))
    # }
    # Learn P(X2 | Z,X1)
    pred.X2.ZX1 = predict(model.X2.ZX1,newdata=data.matrix(data.frame(X1=DATA_Eval$X1,Z=DATA_Eval$Z)),reshape=TRUE)
    prob.X2.ZX1 = mapply(function(idx,x2val){
      pred.X2.ZX1[idx,(x2val+1)]
    },c(1:nrow(DATA_Eval)),DATA_Eval$X2)
    # if (mismode == 1){
    #   prob.X2.ZX1 = fix_pred(mis_pred(prob.X2.ZX1,distortval))
    # }
    UIF_M1 = pred.Y.rx2.X1Z + ((Ir*Ix2)/(prob.X2.ZX1*prob.R.X1))*(Iy - pred.Y.RX2.X1Z)
    # EIF_M1 = UIF_M1 - mean(UIF_M1)
    return(UIF_M1)
  }
  
  compute_UIF_M2 = function(DATA_Train, DATA_Eval, trainedlist, x1valfix, rvalfix){
    # X1 = DATA[,1] 
    # Z = DATA[,2] 
    # R = DATA[,3] 
    # X2 = DATA[,4]  
    # Y = DATA[,5]
    
    ############################################
    # M2 =  M[r | x1;]
    ############################################
    # x1valfix = 0 
    # rvalfix = 0 
    # mylist = TrainModel(DATA_Train,DATA_Eval)
    model.Y = trainedlist[[1]]; model.X2.ZX1 = trainedlist[[2]]; model.R.X1 = trainedlist[[3]]
    
    Ir = (DATA_Eval$R==rvalfix)*1
    Ix1 = (DATA_Eval$X1==x1valfix)*1
    
    prob.X1 = mean(DATA_Train$X1)*DATA_Eval$X1 + (1-mean(DATA_Train$X1))*(1-DATA_Eval$X1)
    
    pred.R.X1 = predict(model.R.X1,newdata=data.matrix(data.frame(X1=DATA_Eval$X1)),reshape=TRUE)
    prob.r.X1 = mapply(function(idx,rvalfix){
      pred.R.X1[idx,(rvalfix+1)]
    },c(1:nrow(DATA_Eval)),rep(rvalfix,nrow(DATA_Eval)))
    # if (mismode == 1){
    #   prob.r.X1 = fix_pred(mis_pred(prob.r.X1,distortval))
    # }
    prob.r.x1 = rep(predict(model.R.X1,newdata = data.matrix(data.frame(X1=x1valfix)),reshape=TRUE)[rvalfix+1],nrow(DATA_Eval))
    UIF_M2 = (Ix1/prob.X1)*(Ir-prob.r.X1)+prob.r.x1
    # EIF_M2 = UIF_M2 - mean(UIF_M2)
    return(UIF_M2)
  }
  
  compute_Yx = function(DATA_Train, DATA_Eval, trainedlist, yvalfix, x1valfix, x2valfix){
    ############################################
    # UIF =  \sum_{r}UIF_M1*mean(UIF_M2) + EIF_M2 * mean(UIF_M1)
    # M1 =  M[y | (r,x2);(x1;z)]
    # M2 =  M[r | x1;]
    ############################################
    # mylist = TrainModel(DATA_Train,DATA_Eval)
    # model.Y = mylist[[1]]; model.X2.ZX1 = mylist[[2]]; model.R.X1 = mylist[[3]]; pred.R.X1 = mylist[[4]]
    
    UIF = rep(0,nrow(DATA_Eval))
    for (rvalfix in Runique){
      UIF_M1 = compute_UIF_M1(DATA_Train=DATA_Train,DATA_Eval=DATA_Eval, trainedlist=trainedlist, yval=yvalfix,rvalfix=rvalfix,x2valfix=x2valfix)
      EIF_M1 = UIF_M1 - mean(UIF_M1,na.rm=T)
      
      UIF_M2 = compute_UIF_M2(DATA_Train=DATA_Train,DATA_Eval=DATA_Eval, trainedlist=trainedlist, x1valfix=x1valfix,rvalfix=rvalfix)
      EIF_M2 = UIF_M1 - mean(UIF_M1,na.rm=T)
      
      UIF = UIF + UIF_M1*mean(UIF_M2,na.rm=T) + EIF_M2*mean(UIF_M1,na.rm=T)
    }
    return(mean(UIF,na.rm=T))
  }
  
  ############################################
  # DATA setup 
  ############################################
  DATA = data.frame(X1,Z,R,X2,Y)
  DATA = subset(DATA,(is.na(X1) == FALSE)&(is.na(Z) == FALSE)&(is.na(R) == FALSE)&(is.na(X2) == FALSE)&(is.na(Y) == FALSE))
  Ndata = nrow(DATA)
  
  tmp = GoodSplit(DATA)
  DATA_Train = tmp[[1]]
  DATA_Eval = tmp[[2]]
  # DATA_Train = DATA
  # DATA_Eval = DATA
  
  trainedlist1 = TrainModel(DATA_Train=DATA_Train, DATA_Eval=DATA_Eval)
  trainedlist2 = TrainModel(DATA_Train=DATA_Eval, DATA_Eval=DATA_Train)
  
  Yx = rep(0,length(X1unique)*length(X2unique))
  idx = 1 
  for (x1valfix in X1unique){
    for (x2valfix in X2unique){
      Yx[idx] = mean(compute_Yx(DATA_Train=DATA_Train, DATA_Eval=DATA_Eval, trainedlist = trainedlist1, yvalfix=yvalfix,x1valfix=x1valfix,x2valfix=x2valfix),
                     compute_Yx(DATA_Train=DATA_Eval, DATA_Eval=DATA_Train, trainedlist = trainedlist2, yvalfix=yvalfix,x1valfix=x1valfix,x2valfix=x2valfix))
      # Yx[idx] = mean(compute_Yx(DATA_Train=DATA_Train, DATA_Eval=DATA_Eval, trainedlist = trainedlist1, yvalfix=yvalfix,x1valfix=x1valfix,x2valfix=x2valfix),
      #                compute_Yx(DATA_Train=DATA_Eval, DATA_Eval=DATA_Train, trainedlist = trainedlist2, yvalfix=yvalfix,x1valfix=x1valfix,x2valfix=x2valfix))
      # Yx[idx] = compute_Yx(DATA_Train=DATA_Train, DATA_Eval=DATA_Eval, trainedlist = trainedlist1, yvalfix=yvalfix,x1valfix=x1valfix,x2valfix=x2valfix)
      Yx[idx] = max(Yx[idx],0)
      Yx[idx] = min(Yx[idx],1)
      idx = idx + 1 
    }
  }
  return(Yx)
}




