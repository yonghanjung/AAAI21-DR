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
    for (idx in 1:4){
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

DREstimator = function(OBS,mismode){
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
  
  tmp = GoodSplit(OBS)
  OBS_train = tmp[[1]]
  OBS_eval = tmp[[2]]
  
  ComputeDR = function(xfix,OBS_train,OBS_eval){
    # Choose the fixed R 
    RProb = rep(0,0,0)
    idx = 1
    for (rval in Runique){
      filtered_DATA = subset(OBS_eval,R==rval & X==xfix)
      RProb[idx] = nrow(filtered_DATA)/nrow(OBS_eval)
      idx = idx + 1 
    } 
    rfix = Runique[which.max(RProb)]
    yfix = 1 
    
    Ir.train = (OBS_train$R == rfix)*1
    Ix.train = (OBS_train$X == xfix)*1
    Iy.train = (OBS_train$Y == yfix)*1
    Ixy.train = Ix.train*Iy.train
    
    Ir.eval = (OBS_eval$R == rfix)*1
    Ix.eval = (OBS_eval$X == xfix)*1
    Iy.eval = (OBS_eval$Y == yfix)*1
    Ixy.eval = Ix.eval*Iy.eval
    
    if (mismode == 1){
      Iy.distorted.train = xor((Iy.train * rbinom(n=length(Iy.train),size=1,prob=0.5)),rbinom(n=length(Iy.train),size=1,prob=0.5))*1
      # Iy.distorted.eval = xor((Iy * rbinom(n=length(Iy.eval),size=1,prob=0.5)),rbinom(n=length(Iy.eval),size=1,prob=0.5))*1
      
      modeled_X.train = OBS_train$X 
      modeled_X.train = modeled_X.train*2
      modeled_X.train[modeled_X.train==4] = 1 

      Ix.distorted.train = (modeled_X.train == xfix)*1
      # Ix.distorted.eval = (modeled_X.eval == xfix)*1
      
      Ixy.distorted.train = Ix.distorted.train * Iy.distorted.train
      # Ixy.distorted.eval = Ix.distorted.eval * Iy.distorted.eval
    }
    
    ## Compute P(x,y|R,W)
    ### Train P(x,y|R,W)
    if (mismode == 1){
      model.xy.RW = learnXG(inVar = data.matrix(data.frame(R=OBS_train$R, W=OBS_train$W)),labelval = Ixy.distorted.train, regval = rep(0,nrow(OBS_train)),binommode = 1)
    }else{
      model.xy.RW = learnXG(inVar = data.matrix(data.frame(R=OBS_train$R, W=OBS_train$W)),labelval = Ixy.train, regval = rep(0,nrow(OBS_train)),binommode = 1)  
    }
    
    ### Evaluate P(x,y|R,W)
    prob.xy.RW = predict(model.xy.RW,newdata=data.matrix(data.frame(R=OBS_eval$R,W=OBS_eval$W)),type='response')
    prob.xy.rW = predict(model.xy.RW,newdata=data.matrix(data.frame(R=rep(rfix,nrow(OBS_eval)),W=OBS_eval$W)),type='response')
    
    
    ## Compute P(x|R,W)
    ### Train P(x|R,W)
    if (mismode == 1){
      model.x.RW = learnXG(inVar = data.matrix(data.frame(R=OBS_train$R, W=OBS_train$W)),labelval = Ix.distorted.train, regval = rep(0,nrow(OBS_train)),binommode = 1)
    }else{
      model.x.RW = learnXG(inVar = data.matrix(data.frame(R=OBS_train$R, W=OBS_train$W)),labelval = Ix.train, regval = rep(0,nrow(OBS_train)),binommode = 1)  
    }
    ### Evaluate P(x|R,W)
    prob.x.RW = predict(model.x.RW,newdata=data.matrix(data.frame(R=OBS_eval$R, W=OBS_eval$W)),type='response')
    prob.x.rW = predict(model.x.RW,newdata=data.matrix(data.frame(R=rep(rfix,nrow(OBS_eval)), W=OBS_eval$W)),type='response')
    
    ## Compute P(R|W)
    modeled_R.train = OBS_train$R 
    modeled_R.eval = OBS_eval$R 
    if (mismode == 2){
      modeled_R.train = 2* modeled_R.train
      modeled_R.train[modeled_R.train==4] = 1
    }
    ### Train P(R|W)
    model.R.W = learnXG(inVar = data.matrix(data.frame(W=OBS_train$W)),labelval = modeled_R.train, regval = rep(0,nrow(OBS_train)),binommode = 0)
    
    ### Evaluate P(R|W)
    pred.R.W = predict(model.R.W, newdata=data.matrix(data.frame(W=OBS_eval$W)),type='response')
    pred.R.W  = t(matrix(pred.R.W,nrow=length(Runique)))
    prob.R.W = rep(0,nrow(OBS_eval))
    for (idx in 1:nrow(OBS_eval)){
      prob.R.W[idx] = pred.R.W[idx,(OBS_eval$R[idx]+1)]
    }
    # if (mismode == 2){
    #   prob.R.W = fix_pred(mis_pred(prob.R.W,distortval))
    # }

    # Compute this by IPW 
    smallval = 0
    
    UIF_M1 = (Ir.eval*(Ixy.eval - prob.xy.RW)/(prob.R.W+smallval)) + prob.xy.rW 
    UIF_M2 = (Ir.eval*(Ix.eval - prob.x.RW)/(prob.R.W+smallval)) + (prob.x.rW)
    prob.xy.dor = mean(UIF_M1)
    prob.x.dor = mean(UIF_M2)
    prob.y.dox = prob.xy.dor/prob.x.dor
    
    EIF_M1 = UIF_M1 - mean(UIF_M1)
    EIF_M2 = UIF_M2 - mean(UIF_M2)
    
    UIF = (1/(prob.x.dor))*(UIF_M1 - EIF_M2*prob.y.dox)
    # UIF = (UIF_M1/prob.x.dor) - (EIF_M2/prob.x.dor)*(prob.xy.dor/prob.x.dor)
    return(mean(UIF))
  }
  
  YxDR = rep(0,length(Xunique))
  idx = 1 
  for (xval in Xunique){
    YxDR[idx] = mean(ComputeDR(xval,OBS_train = OBS_train, OBS_eval = OBS_eval),ComputeDR(xval,OBS_train = OBS_eval, OBS_eval = OBS_train))  
    YxDR[idx] = max(YxDR[idx],0)
    YxDR[idx] = min(YxDR[idx],1)
    idx = idx + 1 
  }
  return(YxDR)
}




