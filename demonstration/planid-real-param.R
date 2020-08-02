
source('RID_functions.R')
source('WERM_Heuristic.R')

paramAdj = function(OBS,D,distortval,mismode){
  smallval = 1e-2 
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
  
  # Setting
  DATA = cbind(X1,Z,R,X2,Y)
  DATA = subset(DATA,(is.na(X1) == FALSE)&(is.na(Z) == FALSE)&(is.na(R) == FALSE)&(is.na(X2) == FALSE)&(is.na(Y) == FALSE))
  DATA = data.frame(DATA)
  
  # Enumerate all possible values of column
  tmp = c()
  tmp = append(tmp,list(X1unique)) # X1
  tmp = append(tmp,list(Zunique)) # Z
  tmp = append(tmp,list(Runique)) # R
  tmp = append(tmp,list(X2unique)) # X2
  tmp = append(tmp,list(Yunique)) # Y
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('X1','Z','R','X2','Y')
  
  # Enumerate all possible values of column
  tmp = c()
  tmp = append(tmp,list(c(0:1))) # X1
  tmp = append(tmp, list(c(0:1))) # Z
  tmp = append(tmp, list(c(0:1))) # R
  tmp = append(tmp,list(c(0,1))) # X2
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('X1','Z','R','X2')
  
  # Compute E[Y=1|x1,z,r,x2]
  ExpYParam = function(myallpossible,DATA,yval=1){
    Iy = (DATA$Y == yval)*1
    modelY = learnXG(as.matrix(DATA[,c('X1','Z','R','X2')]),Iy,rep(0,length(Iy)),binommode = 1)
    evalMat = as.matrix(myallpossible[,c('X1','Z','R','X2')])
    predval = predict(modelY,newdata=evalMat,type='response')
    newcol = (ncol(myallpossible)+1)
    myallpossible[,newcol] = predval
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  Ytable = ExpYParam(allpossible,DATA) # E[Y|X,Z,W]
  if (mismode == 1){
    Ytable$Y = fix_pred(mis_pred(Ytable$Y,distortval))
  }

  # Compute P(z | x1 )
  ProbZ.X1.Param_Real = function(myallpossible,DATA){
    model.Z.X1 = learnXG(as.matrix(DATA[,c('X1')]),DATA$Z,rep(0,length(Z)),binommode = 0)
    evalMat = as.matrix(myallpossible[,c('X1')])
    predval = predict(model.Z.X1,newdata=evalMat,type='response')
    predval = t(matrix(predval,nrow=3))
    prob.Z.X1 = rep(0,nrow(myallpossible))
    for (idx in 1:nrow(myallpossible)){
      zval = myallpossible$Z[idx]
      prob.Z.X1[idx] = predval[idx,(zval+1)]
    }
    newcol = (ncol(myallpossible)+1)
    myallpossible[,newcol] = prob.Z.X1
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  PzTable = ProbZ.X1.Param_Real(allpossible,DATA)
  
  PzTable = allpossible 
  formula.Z = as.formula('Z~X1')
  model.Z = bayesglm(formula.Z,family = binomial(), data=DATA)
  predZ = as.vector(fix_pred(predict.glm(model.Z,newdata = allpossible,type="response")))
  PzTable$prob = predZ*allpossible$Z + (1-predZ)*(1-allpossible$Z)
  if (mismode == 1){
    PzTable$prob = fix_pred(mis_pred(PzTable$prob,distortval))
  }
  
  # Compute P(x1)
  Px1Table = allpossible 
  predX1 = mean(OBS$X1)
  Px1Table$prob = predX1*allpossible$X1 + (1-predX1)*(1-allpossible$X1)
  if (mismode == 1){
    Px1Table$prob = fix_pred(mis_pred(Px1Table$prob,distortval))
  }
  
  # Compute P(r | x1 )
  PrTable = allpossible 
  formula.R = as.formula('R~X1')
  model.R = bayesglm(formula.R,family = binomial(), data=DATA)
  predR = as.vector(fix_pred(predict.glm(model.R,newdata = allpossible,type="response")))
  PrTable$prob = predR*allpossible$R + (1-predR)*(1-allpossible$R)
  if (mismode == 2){
    PrTable$prob = fix_pred(mis_pred(PrTable$prob,distortval))
  }
  
  ComputeVal = allpossible
  ComputeVal$val = Ytable$Y * PzTable$prob * Px1Table$prob
  
  # ComputeVal[ComputeVal$X1==x1val & ComputeVal$Z==z,'val']
  
  
  tmp = c()
  # tmp = append(tmp,list(c(0:1))) # X1
  # tmp = append(tmp, list(c(0:1))) # Z
  tmp = append(tmp, list(c(0:1))) # R
  tmp = append(tmp,list(c(0,1))) # X2
  allpossible_RX2 = expand.grid(tmp)
  colnames(allpossible_RX2) = c('R','X2')
  allpossible_RX2$val = 0
  
  for (x1val in (0:1)){
    for (zval in (0:1)){
      allpossible_RX2$val = allpossible_RX2$val + ComputeVal[ComputeVal$X1==x1val & ComputeVal$Z==zval,'val']
    }
  }
  
  Yx = rep(0,4)
  idx = 1
  for (x1val in (0:1)){
    for (x2val in (0:1)){
      Yx[idx] = mean(PrTable[PrTable$X1==x1val & PrTable$R==0,'prob'])*allpossible_RX2[allpossible_RX2$X2==x2val & allpossible_RX2$R==0,'val'] + 
        mean(PrTable[PrTable$X1==x1val & PrTable$R==1,'prob'])*allpossible_RX2[allpossible_RX2$X2==x2val & allpossible_RX2$R==1,'val']
      idx = idx + 1 
    }
  }
  for (idx in 1:4){
    Yx[idx] = min(Yx[idx],1)
    Yx[idx] = max(Yx[idx],0)
  }
  return(Yx)
}

