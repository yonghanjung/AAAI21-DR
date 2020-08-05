source('RID_functions.R')
source('WERM_Heuristic.R')
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
  # Construct ProbMat 
  # Enumerate all possible values of column
  tmp = c()
  tmp = append(tmp,list(Wunique)) # W
  tmp = append(tmp, list(Runique)) # R
  tmp = append(tmp,list(Xunique)) # X
  tmp = append(tmp,list(Yunique)) # Y
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('W','R','X','Y')
  
  # Compute the Joint prob 
  for (rowidx in 1:nrow(allpossible)){
    wval = allpossible[rowidx,'W']
    rval = allpossible[rowidx,'R']
    xval = allpossible[rowidx,'X']
    yval = allpossible[rowidx,'Y']
    filterDATA = subset(DATA,W==wval&R==rval&X==xval&Y==yval)
    allpossible[rowidx,'prob'] = nrow(filterDATA)/Ndata
  }
  
  ComputeDR = function(xfix){
    # Choose the fixed R 
    RProb = rep(0,0,0)
    idx = 1
    for (rval in Runique){
      filtered_DATA = subset(DATA,R==rval & X==xfix)
      RProb[idx] = nrow(filtered_DATA)/nrow(DATA)
      idx = idx + 1 
    } 
    rfix = Runique[which.max(RProb)]
    # IF_1 
    
    yfix = 1 
    Ir = (R == rfix)*1
    Ix = (X == xfix)*1
    Iy = (Y == yfix)*1
    Ixy = Ix*Iy
    
    if (mismode == 1){
      Iy.distorted = xor((Iy * rbinom(n=length(Iy),size=1,prob=0.5)),rbinom(n=length(Iy),size=1,prob=0.5))*1
      modeled_X = X 
      modeled_X = modeled_X*2
      modeled_X[modeled_X==4] = 1 
      Ix.distorted = (modeled_X == xfix)*1
      Ixy.distorted = Ix.distorted*Iy.distorted
    }
    
    ## Compute P(x,y|R,W)
    ### Compute P(x,y|R,W)
    if (mismode == 1){
      model.xy.RW = learnXG(inVar = data.matrix(data.frame(R,W)),labelval = Ixy.distorted, regval = rep(0,nrow(DATA)),binommode = 1)
    }else{
      model.xy.RW = learnXG(inVar = data.matrix(data.frame(R,W)),labelval = Ixy, regval = rep(0,nrow(DATA)),binommode = 1)  
    }
    prob.xy.RW = predict(model.xy.RW,newdata=data.matrix(data.frame(R,W)),type='response')
    # if (mismode == 1){
    #   prob.xy.RW = fix_pred(mis_pred(prob.xy.RW,distortval))
    # }
    ### Compute P(x|R,W)
    if (mismode == 1){
      model.x.RW = learnXG(inVar = data.matrix(data.frame(R,W)),labelval = Ix.distorted, regval = rep(0,nrow(DATA)),binommode = 1)
    }else{
      model.x.RW = learnXG(inVar = data.matrix(data.frame(R,W)),labelval = Ix, regval = rep(0,nrow(DATA)),binommode = 1)  
    }
    prob.x.RW = predict(model.x.RW,newdata=data.matrix(data.frame(R,W)),type='response')
    # if (mismode == 1){
    #   prob.x.RW = fix_pred(mis_pred(prob.x.RW,distortval))
    # }
    
    ### Compute P(x,y|r,W)
    # model.xy.rW = learnXG(inVar = data.matrix(data.frame(rep(rfix,nrow(OBS)),W)),labelval = Ixy, regval = rep(0,nrow(DATA)),binommode = 1)
    prob.xy.rW = predict(model.xy.RW,newdata=data.matrix(data.frame(R=rep(rfix,nrow(OBS)),W)),type='response')
    # if (mismode == 1){
    #   prob.xy.rW = fix_pred(mis_pred(prob.xy.rW,distortval))
    # }
    
    ### Compute P(x|r,W)
    # model.x.rW = learnXG(inVar = data.matrix(data.frame(rep(rfix,nrow(OBS)),W)),labelval = Ix, regval = rep(0,nrow(DATA)),binommode = 1)
    prob.x.rW = predict(model.x.RW,newdata=data.matrix(data.frame(R=rep(rfix,nrow(OBS)),W)),type='response')
    # if (mismode == 1){
    #   prob.x.rW = fix_pred(mis_pred(prob.x.rW,distortval))
    # }
    
    ### Compute P(R|W)
    modeled_R = R 
    if (mismode == 2){
      modeled_R = 2* modeled_R
      modeled_R[modeled_R==4] = 1
    }
    model.R.W = learnXG(inVar = data.matrix(data.frame(W)),labelval = modeled_R, regval = rep(0,nrow(DATA)),binommode = 0)
    pred.R.W = predict(model.R.W,newdata=data.matrix(data.frame(W)),type='response')
    pred.R.W  = t(matrix(pred.R.W,nrow=3))
    prob.R.W = rep(0,nrow(DATA))
    for (idx in 1:nrow(DATA)){
      prob.R.W[idx] = pred.R.W[idx,(R[idx]+1)]
    }
    # if (mismode == 2){
    #   prob.R.W = fix_pred(mis_pred(prob.R.W,distortval))
    # }

    # Compute this by IPW 
    smallval = 0
    prob.xy.dor = mean((Ir*(Ixy - prob.xy.RW)/(prob.R.W+smallval)) + prob.xy.rW)
    prob.x.dor = mean((Ir*(Ix - prob.x.RW)/(prob.R.W+smallval)) + (prob.x.rW))
    prob.y.dox = prob.xy.dor/prob.x.dor
    
    UIF_M1 = (Ir*(Ixy - prob.xy.RW)/(prob.R.W+smallval)) + prob.xy.rW 
    UIF_M2 = (Ir*(Ix - prob.x.RW)/(prob.R.W+smallval)) + (prob.x.rW)
    EIF_M1 = (Ir*(Ixy - prob.xy.RW)/(prob.R.W+smallval)) + prob.xy.rW - prob.xy.dor
    EIF_M2 = (Ir*(Ix - prob.x.RW)/(prob.R.W+smallval)) + (prob.x.rW - prob.x.dor)
    
    UIF = (1/(prob.x.dor))*(UIF_M1 - EIF_M2*prob.y.dox)
    # UIF = (UIF_M1/prob.x.dor) - (EIF_M2/prob.x.dor)*(prob.xy.dor/prob.x.dor)
    return(mean(UIF))
  }
  YxDR = rep(0,length(Xunique))
  idx = 1 
  for (xval in Xunique){
    YxDR[idx] = ComputeDR(xval)  
    idx = idx + 1 
  }
  return(YxDR)
}




