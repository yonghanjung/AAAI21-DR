source('WERM_Heuristic.R')
source('RID_functions.R')
source('DRModule.R')

WERMEstimator = function(OBS,mismode,seednum){
  numRounds = 20 
  maxDepth = 20
  myverbose = F
  
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
  
  Train.P.Y = function(DATA_Train, DATA_Eval, mylambda){
    yvalfix = 1 
    X1Train = DATA_Train$X1
    ZTrain = DATA_Train$Z
    RTrain = DATA_Train$R
    X2Train = DATA_Train$X2
    Iy.Train = (DATA_Train$Y == yvalfix)*1
    if (mismode == 1){
      Iy.Train = distortVar(Iy.Train,seednum)
    }
    if (mismode == 2){
      RTrain = distortVar(RTrain,seednum)
      X2Train = distortVar(X2Train,seednum)
    }
    
    inVarTrain = data.frame(X1=X1Train,Z=ZTrain,R=RTrain,X2=X2Train)
    inVarEval = data.frame(X1=DATA_Eval$X1,Z=DATA_Eval$Z,R=DATA_Eval$R,X2=DATA_Eval$X2)
    
    model.Y = learnXG(inVar = data.matrix(inVarTrain), labelval = Iy.Train, regval = mylambda, binommode = 1)
    # model.Y = xgboost(verbose = 0, data = data.matrix(inVarTrain), label = Iy.Train, nrounds = numRounds,
    #                   maxdepth=maxDepth, lambda=mylambda, alpha=mylambda/2, objective="binary:logistic")  
    pred.Y = predict(model.Y,newdata=data.matrix(inVarEval),type='response')
    prob.Y = pred.Y*DATA_Eval$Y + (1-pred.Y)*(1-DATA_Eval$Y)
    return(prob.Y)
  }
  
  Train.P.R.X1 = function(DATA_Train, DATA_Eval, mylambda){
    yvalfix = 1 
    X1Train = DATA_Train$X1
    ZTrain = DATA_Train$Z
    RTrain = DATA_Train$R
    X2Train = DATA_Train$X2
    Iy.Train = (DATA_Train$Y == yvalfix)*1
    if (mismode == 1){
      Iy.Train = distortVar(Iy.Train,seednum)
    }
    if (mismode == 2){
      RTrain = distortVar(RTrain,seednum)
      X2Train = distortVar(X2Train,seednum)
    }
    
    inVarTrain = data.frame(X1=X1Train)
    inVarEval = data.frame(X1=DATA_Eval$X1)
  
    model.R.X1 =learnXG(data.matrix(inVarTrain),RTrain,mylambda,binommode=0)
    
    pred.R = t(matrix(predict(model.R.X1, newdata=as.matrix(inVarEval), type='response'), nrow=length(Runique)))
    prob.R.X1 = mapply(function(rowidx,rval){
      return(pred.R[rowidx,rval+1])
    },c(1:nrow(DATA_Eval)),DATA_Eval$R)
    # if (mismode == 2){
    #   prob.R.X1 = mis_pred(prob.R.X1)
    # }
    return(prob.R.X1)
  }
  
  Train.P.R.X2X1Z = function(DATA_Train, DATA_Eval, mylambda){
    yvalfix = 1 
    X1Train = DATA_Train$X1
    ZTrain = DATA_Train$Z
    RTrain = DATA_Train$R
    X2Train = DATA_Train$X2
    Iy.Train = (DATA_Train$Y == yvalfix)*1
    if (mismode == 1){
      Iy.Train = distortVar(Iy.Train,seednum)
    }
    if (mismode == 2){
      RTrain = distortVar(RTrain,seednum)
      X2Train = distortVar(X2Train,seednum)
    }
    
    inVarTrain = data.frame(X2=X2Train, X1=X1Train, Z=ZTrain)
    inVarEval = data.frame(X2=DATA_Eval$X2, X1=DATA_Eval$X1, Z=DATA_Eval$Z)
    
    model.R.X2X1Z =learnXG(data.matrix(inVarTrain),RTrain,mylambda,binommode=0)
    
    pred.R = t(matrix(predict(model.R.X2X1Z, newdata=as.matrix(inVarEval), type='response'), nrow=length(Runique)))
    prob.R.X2X1Z = mapply(function(rowidx,rval){
      return(pred.R[rowidx,rval+1])
    },c(1:nrow(DATA_Eval)),DATA_Eval$R)
    # if (mismode == 2){
    #   prob.R.X1 = mis_pred(prob.R.X1)
    # }
    return(prob.R.X2X1Z)
  }
  
  
  Train.P.X2.X1Z = function(DATA_Train, DATA_Eval, mylambda){
    yvalfix = 1 
    X1Train = DATA_Train$X1
    ZTrain = DATA_Train$Z
    RTrain = DATA_Train$R
    X2Train = DATA_Train$X2
    Iy.Train = (DATA_Train$Y == yvalfix)*1
    if (mismode == 1){
      Iy.Train = distortVar(Iy.Train,seednum)
    }
    if (mismode == 2){
      RTrain = distortVar(RTrain,seednum)
      X2Train = distortVar(X2Train,seednum)
    }
    
    inVarTrain = data.frame(X1=X1Train,Z=ZTrain)
    inVarEval = data.frame(X1=DATA_Eval$X1,Z=DATA_Eval$Z)
    
    model.X2.X1Z =learnXG(data.matrix(inVarTrain),X2Train,mylambda,binommode=0)
    
    pred.X2 = t(matrix(predict(model.X2.X1Z, newdata=as.matrix(inVarEval), type='response'), nrow=length(X2unique)))
    prob.X2.X1Z = mapply(function(rowidx,x2val){
      return(pred.X2[rowidx,x2val+1])
    },c(1:nrow(DATA_Eval)),DATA_Eval$X2)
    # if (mismode == 2){
    #   prob.X2.X1Z = mis_pred(prob.X2.X1Z)
    # }
    # 
    return(prob.X2.X1Z)
  }
  
  Train.P.R = function(DATA_Train, DATA_Eval, mylambda){
    yvalfix = 1 
    X1Train = DATA_Train$X1
    ZTrain = DATA_Train$Z
    RTrain = DATA_Train$R
    X2Train = DATA_Train$X2
    Iy.Train = (DATA_Train$Y == yvalfix)*1
    if (mismode == 1){
      Iy.Train = distortVar(Iy.Train,seednum)
    }
    if (mismode == 2){
      RTrain = distortVar(RTrain,seednum)
      X2Train = distortVar(X2Train,seednum)
    }
    
    prob.R.Array = rep(0,length(Runique))
    for (rval in Runique){
      prob.R.Array[rval+1] =  sum((RTrain==rval)*1) / nrow(DATA_Train) 
    }
    prob.R = mapply(function(rval){
      return(prob.R.Array[rval+1])
    },DATA_Eval$R)
    # if (mismode == 2){
    #   prob.R = mis_pred(prob.R)
    # }
    return(prob.R)
  }
  
  Train.P.X2 = function(DATA_Train, DATA_Eval, mylambda){
    yvalfix = 1 
    X1Train = DATA_Train$X1
    ZTrain = DATA_Train$Z
    RTrain = DATA_Train$R
    X2Train = DATA_Train$X2
    Iy.Train = (DATA_Train$Y == yvalfix)*1
    if (mismode == 1){
      Iy.Train = distortVar(Iy.Train,seednum)
    }
    if (mismode == 2){
      RTrain = distortVar(RTrain,seednum)
      X2Train = distortVar(X2Train,seednum)
    }
    
    prob.X2.Array = rep(0,length(X2unique))
    for (x2val in X2unique){
      prob.X2.Array[x2val+1] =  sum((X2Train==x2val)*1) / nrow(DATA_Train) 
    }
    prob.X2 = mapply(function(x2val){
      return(prob.X2.Array[x2val+1])
    },DATA_Eval$X2)
    # if (mismode == 2){
    #   prob.X2 = mis_pred(prob.X2)
    # }
    return(prob.X2)
  }
  
  Train.weighted.Y = function(DATA_Train, DATA_Eval, mylambda,SW_Y){
    yvalfix = 1 
    X1Train = DATA_Train$X1
    ZTrain = DATA_Train$Z
    RTrain = DATA_Train$R
    X2Train = DATA_Train$X2
    Iy.Train = (DATA_Train$Y == yvalfix)*1
    if (mismode == 1){
      Iy.Train = distortVar(Iy.Train,seednum)
    }
    if (mismode == 2){
      RTrain = distortVar(RTrain,seednum)
      X2Train = distortVar(X2Train,seednum)
    }
    
    inVarTrain = data.frame(R=RTrain,X2=X2Train)
    inVarEval = data.frame(R=DATA_Eval$R,X2=DATA_Eval$X2)
    
    # xgbMatrix = xgb.DMatrix(data.matrix(inVar_train), label=Iy.Train,weight =learned_W)
    # modelY_xgboost = xgboost(verbose=0, data=xgbMatrix,nrounds = numRounds,max.depth=maxDepth,lambda=mylambda, alpha=mylambda/2, objective = "binary:logistic")
    
    # Learn P^{SW_Y}(y|r,x2)
    Ybox = rep(0,nrow(DATA_Eval))
    bootstrap_iter = 1
    for (idx in 1:bootstrap_iter){
      # sampled_df = WERM_Sampler(DATA_Train,SW_Y)
      # Learn Pw(y|r,x2)
      
      xgbMatrix = xgb.DMatrix(data.matrix(inVarTrain), label=Iy.Train,weight = SW_Y)
      model.weighted.Y.rx2 = xgboost(verbose=0, data=xgbMatrix,nrounds = numRounds,max.depth=maxDepth,lambda=mylambda, alpha=mylambda/2)
      
      # model.weighted.Y.rx2 = learnXG(inVar = data.matrix(inVarTrain),labelval = Iy.Train, regval = mylambda, binommode = 1)
      # model.weighted.Y.rx2 = xgboost(verbose=0, data=data.matrix(inVarTrain), label = Iy.Train, nrounds = numRounds,max.depth=maxDepth,lambda=mylambda,alpha=mylambda/2,objective = "binary:logistic",weight = SW_Y)
      pred.weighted.Y.rx2 = predict(model.weighted.Y.rx2, newdata = data.matrix(inVarEval),type='response')
      Prob.weighted.Y.rx2 = DATA_Eval$Y*pred.weighted.Y.rx2 + (1-DATA_Eval$Y)*(1-pred.weighted.Y.rx2)
      Ybox = Ybox + Prob.weighted.Y.rx2
    }
    Prob.weighted.Y.rx2 = Ybox/bootstrap_iter
    return(Prob.weighted.Y.rx2)
  } 
  
  Train.SW = function(DATA_Train, DATA_Eval, mylambda,SW_importance_sampling){
    yvalfix = 1 
    X1Train = DATA_Train$X1
    ZTrain = DATA_Train$Z
    RTrain = DATA_Train$R
    X2Train = DATA_Train$X2
    Iy.Train = (DATA_Train$Y == yvalfix)*1
    YTrain = DATA_Train$Y
    if (mismode == 1){
      Iy.Train = distortVar(Iy.Train,seednum)
      YTrain = distortVar(YTrain,seednum)
    }
    if (mismode == 2){
      RTrain = distortVar(RTrain,seednum)
      X2Train = distortVar(X2Train,seednum)
    }
    
    inVarTrain = data.frame(X1=X1Train,Z=ZTrain,R=RTrain,X2=X2Train,Y=YTrain)
    inVarEval = data.frame(X1=DATA_Eval$X1,Z=DATA_Eval$Z,R=DATA_Eval$R,X2=DATA_Eval$X2,Y=DATA_Eval$Y)
    
    model_SW = xgboost(verbose = 0, data = data.matrix(inVarTrain), label = SW_importance_sampling, 
                       nrounds = numRounds, max.depth=maxDepth, lambda=mylambda,alpha=mylambda/2)
    
    pred.SW = predict(model_SW,newdata = data.matrix(inVarEval), type='response')
    
    # Handling if W contains the negative value. 
    if (length(pred.SW[pred.SW < 0 ]) > 0){
      pred.SW[pred.SW < 0 ] = runif(n=length(pred.SW[pred.SW < 0 ]),min=0,max=min(abs(pred.SW)))  
    }
    if (mismode == 2){
      distortval = 1
      sign_rv = 2*rbinom(n=length(pred.SW),size=1,prob=0.5)-1
      pred.SW = pred.SW + sign_rv*rnorm(length(pred.SW),distortval,1)
      pred.SW[pred.SW <= 0] = runif(n=length(pred.SW[pred.SW <= 0]),min=0.5,max=max(pred.SW))
    }
    return(pred.SW)
  }
  
  
  Train.Yx = function(DATA_Train, DATA_Eval, mylambda, x1val,x2val,learned_W){
    yvalfix = 1 
    X1Train = DATA_Train$X1
    ZTrain = DATA_Train$Z
    RTrain = DATA_Train$R
    X2Train = DATA_Train$X2
    Iy.Train = (DATA_Train$Y == yvalfix)*1
    if (mismode == 1){
      Iy.Train = distortVar(Iy.Train,seednum)
    }
    if (mismode == 2){
      RTrain = distortVar(RTrain,seednum)
      X2Train = distortVar(X2Train,seednum)
    }
    
    inVar_train=data.frame(X1=X1Train,X2=X2Train)
    inVar_eval=data.frame(X1=rep(x1val,nrow(DATA_Eval)),X2=rep(x2val,nrow(DATA_Eval)))
    
    xgbMatrix = xgb.DMatrix(data.matrix(inVar_train), label=Iy.Train, weight =learned_W)
    modelY_xgboost = xgboost(weight = learned_W,verbose=0, data=xgbMatrix,nrounds = numRounds,max.depth=maxDepth,lambda=mylambda, alpha=mylambda/2)
    
    predY = predict(modelY_xgboost,newdata=data.matrix(inVar_eval),type='response')
    Yx = mean(predY,na.rm=T)
    return(Yx)    
  }
 
  
  
  ############################################
  # DATA setup 
  ############################################
  DATA = data.frame(X1,Z,R,X2,Y)
  inVar = data.frame(X1,Z,R,X2)
  
  tmp = GoodSplit(DATA)
  DATA_Train = tmp[[1]]
  DATA_Eval = tmp[[2]]
  # if (nrow(DATA_Train) < nrow(DATA_Eval)){
  #   DATA_Eval = DATA_Eval[c(1:nrow(DATA_Train)),]
  # }
  # if (nrow(DATA_Eval) < nrow(DATA_Train)){
  #   DATA_Train = DATA_Train[c(1:nrow(DATA_Eval)),]
  # }
   
  
  mylambda = rep(100/sqrt(nrow(DATA)),nrow(DATA)/2)
  
  prob.Y.1 = Train.P.Y(DATA_Train,DATA_Eval,mylambda)
  prob.Y.2 = RunTryCatchProb_WERM(Train.P.Y,DATA_Train,DATA_Eval,mylambda) 
  prob.Y = (prob.Y.1 + prob.Y.2)/2
  
  
  prob.R.X1.1 = Train.P.R.X1(DATA_Train,DATA_Eval,mylambda)
  prob.R.X1.2 = RunTryCatchProb_WERM(Train.P.R.X1,DATA_Train,DATA_Eval,mylambda) 
  prob.R.X1 = (prob.R.X1.1 + prob.R.X1.2)/2
  
  
  prob.X2.X1Z.1 = Train.P.X2.X1Z(DATA_Train,DATA_Eval,mylambda)
  prob.X2.X1Z.2 = RunTryCatchProb_WERM(Train.P.X2.X1Z,DATA_Train,DATA_Eval,mylambda) 
  prob.X2.X1Z = (prob.R.X1.1 + prob.R.X1.2)/2
  
  
  prob.R.1 = Train.P.R(DATA_Train,DATA_Eval,mylambda)
  prob.R.2 = RunTryCatchProb_WERM(Train.P.R,DATA_Train,DATA_Eval,mylambda) 
  prob.R = (prob.R.1 + prob.R.2)/2
  
  
  
  prob.X2.1 = Train.P.X2(DATA_Train,DATA_Eval,mylambda)
  prob.X2.2 = RunTryCatchProb_WERM(Train.P.X2,DATA_Train,DATA_Eval,mylambda) 
  prob.X2 = (prob.X2.1+prob.X2.2)/2
  
  
  prob.R.X2X1Z.1 = Train.P.R.X2X1Z(DATA_Train,DATA_Eval,mylambda)
  prob.R.X2X1Z.2 = RunTryCatchProb_WERM(Train.P.R.X2X1Z,DATA_Train,DATA_Eval,mylambda) 
  prob.R.X2X1Z = (prob.R.X2X1Z.1+prob.R.X2X1Z.2)/2
  
  
  # Weight SW_Y  
  SW_Y = (prob.R * prob.X2)/(prob.R.X1 * prob.X2.X1Z)
  
  Prob.weighted.Y.rx2.1 = Train.weighted.Y(DATA_Train,DATA_Eval,mylambda,SW_Y)
  Prob.weighted.Y.rx2.2 = Train.weighted.Y(DATA_Eval,DATA_Train,mylambda,SW_Y)
  Prob.weighted.Y.rx2 = (Prob.weighted.Y.rx2.1+Prob.weighted.Y.rx2.2)/2
  
  SW_importance_sampling = (Prob.weighted.Y.rx2)/(prob.Y)
  # SW_importance_sampling = (Prob.weighted.Y.rx2)/(prob.Y)
  # if (mismode == 2){
  #   SW_importance_sampling = 2*SW_importance_sampling 
  # }
  # SW.1 = Train.SW(DATA_Train, DATA_Eval, mylambda,SW_importance_sampling)
  # SW.2 = Train.SW(DATA_Eval, DATA_Train, mylambda,SW_importance_sampling)
  # learned_W = SW.1
  learned_W = SW_importance_sampling
  if (myverbose){
    print(paste(mismode,"prob.Y",round(mean(prob.Y,na.rm=T),4),sep="-"))
    print(paste(mismode,"prob.R.X1",round(mean(prob.R.X1,na.rm=T),4),sep="-"))
    print(paste(mismode,"prob.X2.X1Z",round(mean(prob.X2.X1Z,na.rm=T),4),sep="-"))
    print(paste(mismode,"prob.R",round(mean(prob.R,na.rm=T),4),sep="-"))
    print(paste(mismode,"prob.X2",round(mean(prob.X2,na.rm=T),4),sep="-"))
    print(paste(mismode,"prob.R.X2X1Z",round(mean(prob.R.X2X1Z,na.rm=T),4),sep="-"))
    print(paste(mismode,"SW_Y",round(mean(SW_Y,na.rm=T),4),sep="-"))
    print(paste(mismode,"Prob.weighted.Y.rx2",round(mean(Prob.weighted.Y.rx2,na.rm=T),4),sep="-"))
    print(paste(mismode,"SW_importance_sampling",round(mean(SW_importance_sampling,na.rm=T),4),sep="-"))
    print(paste(mismode,"learned_W",round(mean(learned_W,na.rm=T),4),sep="-"))
  }
  
  YxWERM = rep(0,length(X1unique)*length(X2unique))
  idx = 1
  for (x1val in X1unique){
    for (x2val in X2unique){
      myresult1 = Train.Yx(DATA_Train, DATA_Eval, mylambda, x1val, x2val,learned_W)
      tryCatch(
        expr = {
          myresult2 = Train.Yx(DATA_Eval, DATA_Train , mylambda, x1val, x2val, learned_W)
        },
        error = function(e){
          print("Error in the Function")
        },
        finally = {
          myresult2 = Train.Yx(DATA_Train, DATA_Eval , mylambda, x1val, x2val, learned_W)
        }
      )
      YxWERM[idx] = (myresult1+myresult2)/2
      idx = idx + 1 
    }
  }
  
  return(YxWERM)
}





