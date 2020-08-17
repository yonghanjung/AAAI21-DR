source('RID_functions.R')
source('WERM_Heuristic.R')
source('DRModule.R')

WERMEstimator = function(OBS,mismode,seednum){
  numRounds = 10 
  maxDepth = 20
  
  Train.P.R.W = function(DATA_Train, DATA_Eval, mylambda){
    Rtrain = DATA_Train$R 
    if (mismode == 2){
      Rtrain = distortVar(Rtrain,seednum)
    }
    model.R.W = learnXG(inVar = data.matrix(data.frame(W=DATA_Train$W)),labelval = Rtrain, regval = mylambda, binommode = 0)  
    pred.R.W = predict(model.R.W,newdata=data.matrix(data.frame(W=DATA_Eval$W)),type='response')
    pred.R.W  = t(matrix(pred.R.W,nrow=length(Runique)))
    prob.R.W = rep(0,nrow(DATA_Eval))
    for (idx in 1:nrow(DATA_Eval)){
      prob.R.W[idx] = pred.R.W[idx,(DATA_Eval$R[idx]+1)]
    }
    return(prob.R.W)
  }
  
  Train.P.R = function(DATA_Train, DATA_Eval, mylambda){
    Rtrain = DATA_Train$R 
    if (mismode == 2){
      Rtrain = distortVar(Rtrain,seednum)
    }
    prob.R = mapply(function(rval){
      jointprob.R = sum((Rtrain == rval)*1)/nrow(DATA_Train)
      return(jointprob.R)
    },DATA_Train$R)
    return(prob.R)
  }
  
  Train.SW = function(DATA_Train, DATA_Eval, mylambda,SW_importance_sampling){
    Wtrain = DATA_Train$W 
    Rtrain = DATA_Train$R 
    if (mismode == 2){
      Rtrain = distortVar(Rtrain,seednum)
    }
    
    model_SW = xgboost(verbose = 0, data = data.matrix(data.frame(W=Wtrain,R=Rtrain)), label = SW_importance_sampling, 
                           nrounds = numRounds, max.depth=maxDepth, lambda=mylambda,alpha=mylambda/2)
    
    pred.SW = predict(model_SW,newdata = data.matrix(data.frame(W=DATA_Eval$W,R=DATA_Eval$R)), type='response')
    
    # Handling if W contains the negative value. 
    if (length(pred.SW[pred.SW < 0 ]) > 0){
      pred.SW[pred.SW < 0 ] = runif(n=length(pred.SW[pred.SW < 0 ]),min=0,max=min(abs(pred.SW)))  
    }
    return(pred.SW)
  }
  
  Train.Yx = function(DATA_Train, DATA_Eval, mylambda,xval,learned_W){
    yvalfix = 1 
    IyTrain = (DATA_Train$Y == yvalfix)*1
    Xtrain = DATA_Train$X 
    Rtrain = DATA_Train$R 
    if (mismode == 1){
      IyTrain = distortVar(IyTrain,seednum)
      Xtrain = distortVar(Xtrain,seednum)
    }
    if (mismode == 2){
      Rtrain = distortVar(Rtrain,seednum)
    }
    
    inVar_train=data.frame(X=Xtrain,R=Rtrain)
    inVar_eval=data.frame(X=rep(xval,nrow(DATA_Eval)),R=DATA_Eval$R)
    
    xgbMatrix = xgb.DMatrix(data.matrix(inVar_train), label=IyTrain)
    modelY_xgboost = xgboost(verbose=0, data=xgbMatrix,nrounds = numRounds,max.depth=maxDepth,lambda=mylambda, alpha=mylambda/2, objective = "binary:logistic", weight = learned_W)
    
    predY = predict(modelY_xgboost,newdata=data.matrix(inVar_eval),type='response')
    Yx = mean(predY)
    return(Yx)    
  }
  
  W = OBS[,1] # High dim surrogate
  R = OBS[,2] # Cofounder 0-numCate
  X = OBS[,3]
  Y = OBS[,4]
  Wunique = unique(W)[order(unique(W))]
  Runique = unique(R)[order(unique(R))]
  Xunique = unique(X)[order(unique(X))]
  Yunique = unique(Y)[order(unique(Y))]
  
  DATA = data.frame(W,R,X,Y)
  DATA = subset(DATA,(is.na(W) == FALSE)&(is.na(R) == FALSE)&(is.na(X) == FALSE)&(is.na(Y) == FALSE))
  Ndata = nrow(DATA)
  
  tmp = GoodSplit(DATA)
  DATA_Train = tmp[[1]]
  DATA_Eval = tmp[[2]]
  
  

  
  
  # Setting
  
  
  
  mylambda = rep(100/sqrt(nrow(DATA)),nrow(DATA))
  prob.R.W.1 = Train.P.R.W(DATA_Train,DATA_Eval,mylambda)
  prob.R.W.2 = Train.P.R.W(DATA_Eval,DATA_Train,mylambda)
  prob.R.W = (prob.R.W.1+prob.R.W.2)/2
  
  prob.R.1 = Train.P.R(DATA_Train,DATA_Eval,mylambda)
  prob.R.2 = Train.P.R(DATA_Eval,DATA_Train,mylambda)
  prob.R = (prob.R.1 + prob.R.2)/2
  
  # regvallist = seq(0,10,by=0.2)
  SW_importance_sampling = prob.R/prob.R.W
  learned_W1 = Train.SW(DATA_Train, DATA_Eval, mylambda,SW_importance_sampling)
  learned_W2 = Train.SW(DATA_Eval, DATA_Train, mylambda,SW_importance_sampling)
  learned_W = (learned_W1+learned_W2)/2
  
  lambda_h = mylambda
  YxWERM = rep(0,length(Xunique))
  idx = 1 
  for (xval in Xunique){
    # # Choose the fixed R 
    # RProb = rep(0,length(Xunique))
    # rval_idx = 1
    # for (rval in Runique){
    #   filtered_DATA = subset(DATA,R==rval & X==xval)
    #   RProb[rval_idx] = nrow(filtered_DATA)/nrow(DATA)
    #   rval_idx = rval_idx + 1 
    # } 
    # rfix = Runique[which.max(RProb)]
    rfix = 0
    myresult1 = Train.Yx(DATA_Train, DATA_Eval, mylambda, xval,learned_W)
    myresult2 = Train.Yx(DATA_Eval,DATA_Train , mylambda, xval,learned_W)
    YxWERM[idx] = (myresult1+myresult2)/2
    idx = idx + 1 
  }
  return(YxWERM)
}







