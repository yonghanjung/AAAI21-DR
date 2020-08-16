source('RID_functions.R')
source('WERM_Heuristic.R')
PlugInEstimator = function(OBS,mismode,seednum){
  X1 = OBS[,1] 
  Z = OBS[,2] 
  R = OBS[,3] 
  X2 = OBS[,4]  
  Y = OBS[,5]
  DATA = data.frame(X1,Z,R,X2,Y)
  
  X1unique = unique(X1)[order(unique(X1))]
  Zunique = unique(Z)[order(unique(Z))]
  Runique = unique(R)[order(unique(R))]
  X2unique = unique(X2)[order(unique(X2))]
  Yunique = unique(Y)[order(unique(Y))]
  
  IyTrain = Y 
  X2Train = X2
  Rtrain = R 
  if (mismode == 1){
    IyTrain = distortVar(IyTrain,seednum)
    # X2Train = distortVar(X2Train,seednum)
  }
  if (mismode == 2){
    Rtrain = distortVar(Rtrain,seednum)
    X2Train = distortVar(X2Train,seednum)
  }
  
  # Setting
  tmp = c()
  tmp = append(tmp,list(X1unique)) # X1 # RiskAversion: 0 highest 
  tmp = append(tmp, list(Zunique)) # Z
  tmp = append(tmp,list(Runique)) # R
  tmp = append(tmp,list(X2unique)) # X2 # Accident: 3 highest 
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('X1','Z','R','X2')
  
  # Causal Query 
  ## P(y|do(x)) = \sum_{r}P(r|x1)\sum_{x1',z}P(y|x1,z,r,x2)P(z,x1')
  ## Compute P(y|x1,z,r,x2)
  Expect.Y = function(myallpossible,DATA,yval){
    Iy = (DATA$Y == yval)*1
    modelY = learnXG(as.matrix(DATA[,c('X1','Z','R','X2')]),IyTrain,rep(0,length(Iy)),binommode = 1)
    evalMat = as.matrix(myallpossible[,c('X1','Z','R','X2')])
    predval = predict(modelY,newdata=evalMat,type='response')
    newcol = (ncol(myallpossible)+1)
    myallpossible[,newcol] = predval
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  ## Compute P(z|x1)
  Prob.Z.X1 = function(myallpossible,DATA){
    modelZ = learnXG(as.matrix(DATA[,c('X1')]),Z,rep(0,length(Z)),binommode = 0)
    evalMat = as.matrix(myallpossible[,c('X1')])
    predval = predict(modelZ,newdata=evalMat,type='response')
    predval = t(matrix(predval,nrow=length(Zunique)))
    probZ = rep(0,nrow(myallpossible))
    for (idx in 1:nrow(myallpossible)){
      zval = myallpossible$Z[idx]
      probZ[idx] = predval[idx,(zval+1)]
    }
    newcol = (ncol(myallpossible)+1)
    myallpossible[,newcol] = probZ
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  ## Compute P(r|x1)
  Prob.R.X1 = function(myallpossible,DATA){
    modelR = learnXG(as.matrix(DATA[,c('X1')]),Rtrain,rep(0,length(R)),binommode = 0)
    evalMat = as.matrix(myallpossible[,c('X1')])
    predval = predict(modelR,newdata=evalMat,type='response')
    predval = t(matrix(predval,nrow=length(Runique)))
    probR = rep(0,nrow(myallpossible))
    for (idx in 1:nrow(myallpossible)){
      rval = myallpossible$R[idx]
      probR[idx] = predval[idx,(rval+1)]
    }
    newcol = (ncol(myallpossible)+1)
    myallpossible[,newcol] = probR
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  Ytable = allpossible
  Ytable = Expect.Y(Ytable,DATA,1)
  
  # if (mismode == 1){
  #   Ytable[,'prob'] = fix_pred(mis_pred(Ytable[,'prob']))  
  # }
  
  
  Prob.Z.X1.Table = allpossible 
  Prob.Z.X1.Table = Prob.Z.X1(Prob.Z.X1.Table,DATA) 
  
  # if (mismode == 1){
  #   Prob.Z.X1.Table[,'prob'] = fix_pred(mis_pred(Prob.Z.X1.Table[,'prob']))  
  # }
  
  Prob.R.X1.Table = allpossible 
  Prob.R.X1.Table = Prob.R.X1(Prob.R.X1.Table,DATA) 
  
  # if (mismode == 1){
  #   Prob.R.X1.Table[,'prob'] = fix_pred(mis_pred(Prob.R.X1.Table[,'prob']))  
  # }
  
  idx = 1 
  Array.Prob.X1 = rep(0,length(X1unique))
  for (x1val in X1unique){
    Array.Prob.X1[idx] = nrow(subset(DATA,X1==x1val))/nrow(DATA)
    idx = idx + 1 
  }
  Prob.X1.Table = allpossible 
  Prob.X1.vector = mapply(function(x1val){
    return(Array.Prob.X1[x1val+1])
  }, Prob.X1.Table$X1)
  Prob.X1.Table[,'prob'] = Prob.X1.vector
  
  
  # Compute P(z,x1)
  Pzx1Table = allpossible 
  Pzx1Table[,'prob'] = Prob.X1.Table[,'prob'] * Prob.Z.X1.Table[,'prob']
  
  ComputeVal = allpossible
  ComputeVal$val1 = Ytable$prob * Pzx1Table$prob # P(y | x1,x2,r,z) * P(z,x1)
  ComputeVal$val2 = Prob.R.X1.Table$prob # P(r|x1)
  
  ## Marginalizing over X1,Z
  tmp = c()
  tmp = append(tmp,list(Runique)) # R
  tmp = append(tmp,list(X2unique)) # X2
  allpossible.Marginover.X1Z = expand.grid(tmp)
  colnames(allpossible.Marginover.X1Z) = c('R','X2')
  allpossible.Marginover.X1Z[,(ncol(allpossible.Marginover.X1Z)+1)] = 0
  colnames(allpossible.Marginover.X1Z)[ncol(allpossible.Marginover.X1Z)] = 'val1'
  
  for (rval in Runique){
    for(x2val in X2unique){
      allpossible.Marginover.X1Z[allpossible.Marginover.X1Z$R==rval & allpossible.Marginover.X1Z$X2==x2val,'val1'] = sum(ComputeVal[ComputeVal$X2==x2val & ComputeVal$R==rval,'val1'],na.rm=T)
    }
  }
  
  ## For all R,X1
  tmp = c()
  tmp = append(tmp,list(Runique)) # R
  tmp = append(tmp,list(X1unique)) # X1
  allpossible.X1R = expand.grid(tmp)
  colnames(allpossible.X1R) = c('R','X1')
  allpossible.X1R[,(ncol(allpossible.X1R)+1)] = 0
  colnames(allpossible.X1R)[ncol(allpossible.X1R)] = 'val1'
  
  for (rval in Runique){
    for(x1val in X1unique){
      allpossible.X1R[allpossible.X1R$R==rval & allpossible.X1R$X1==x1val,'val1'] = mean(ComputeVal[ComputeVal$X1==x1val & ComputeVal$R==rval,'val2'],na.rm=T)
    }
  }
  
  Yx = rep(0,length(X1unique)*length(X2unique))
  idx = 1 
  for (x1val in X1unique){
    for (x2val in X2unique){
      Yx[idx] = sum(allpossible.Marginover.X1Z[allpossible.Marginover.X1Z$X2==x2val,'val1'] * allpossible.X1R[allpossible.X1R$X1==x1val,'val1'])
      idx = idx + 1 
    }
  }
  return(Yx)
}

