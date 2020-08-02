source('RID_functions.R')
source('WERM_Heuristic.R')

DREstimator = function(OBS,distortval,mismode){
  X1 = OBS[,1] 
  Z = OBS[,2] 
  R = OBS[,3] 
  X2 = OBS[,4]  
  Y = OBS[,5]
  
  yvalfix = 1 
  Iy = (Y==yvalfix)*1

  ############################################
  # DATA setup 
  ############################################
  DATA = data.frame(X1,Z,R,X2,Y)
  X1unique = unique(X1)[order(unique(X1))]
  Zunique = unique(Z)[order(unique(Z))]
  Runique = unique(R)[order(unique(R))]
  X2unique = unique(X2)[order(unique(X2))]
  Yunique = unique(Y)[order(unique(Y))]
  
  ############################################
  # Learn Models 
  ############################################
  model.Y = learnXG(inVar = data.matrix(data.frame(X1,Z,R,X2)),labelval = Y, regval = rep(0,nrow(DATA)),binommode = 1)
  model.X2.ZX1 = learnXG(inVar = data.matrix(data.frame(X1,Z)),labelval = X2, regval = rep(0,nrow(DATA)),binommode = 0)
  model.R.X1 = learnXG(inVar = data.matrix(data.frame(X1)),labelval = R, regval = rep(0,nrow(DATA)),binommode = 0)
  pred.R.X1 = predict(model.R.X1,newdata=data.matrix(data.frame(X1)),reshape=TRUE)
  
  compute_UIF_M1 = function(DATA,yval,rvalfix,x2valfix,mismode,distortval){
    X1 = DATA[,1] 
    Z = DATA[,2] 
    R = DATA[,3] 
    X2 = DATA[,4]  
    Y = DATA[,5]
    
    # yval = 1 
    # Iy = (Y==yval)*1
    
    ############################################
    # M1 =  M[y | (r,x2);(x1;z)]
    ############################################
    # rvalfix = 0 
    # x2valfix = 0
    Ir = (R==rvalfix)*1
    Ix2 = (X2==x2valfix)*1
    
    # Learn P(y | r,x2,X1,Z)
    pred.Y.rx2.X1Z = predict(model.Y,newdata=data.matrix(data.frame(X1,Z,R=rep(rvalfix,nrow(DATA)),X2=rep(x2valfix,nrow(DATA)))),type='response')
    if (mismode == 1){
      pred.Y.rx2.X1Z = fix_pred(mis_pred(pred.Y.rx2.X1Z,distortval))
    }
    # Learn P(y | R,X2,X1,Z)
    pred.Y.RX2.X1Z = predict(model.Y,newdata=data.matrix(data.frame(X1,Z,R,X2)),type='response')
    if (mismode == 1){
      pred.Y.RX2.X1Z = fix_pred(mis_pred(pred.Y.RX2.X1Z,distortval))
    }
    # Learn P(R | X1)
    pred.R.X1 = predict(model.R.X1,newdata=data.matrix(data.frame(X1)),reshape=TRUE)
    prob.R.X1 = mapply(function(idx,rval){
      pred.R.X1[idx,rval+1]
    },c(1:nrow(DATA)),R)
    if (mismode == 1){
      prob.R.X1 = fix_pred(mis_pred(prob.R.X1,distortval))
    }
    # Learn P(X2 | Z,X1)
    pred.X2.ZX1 = predict(model.X2.ZX1,newdata=data.matrix(data.frame(X1,Z)),reshape=TRUE)
    prob.X2.ZX1 = mapply(function(idx,x2val){
      pred.X2.ZX1[idx,x2val+1]
    },c(1:nrow(DATA)),X2)
    if (mismode == 1){
      prob.X2.ZX1 = fix_pred(mis_pred(prob.X2.ZX1,distortval))
    }
    UIF_M1 = pred.Y.rx2.X1Z + ((Ir*Ix2)/(prob.X2.ZX1*prob.R.X1))*(Iy - pred.Y.RX2.X1Z)
    # EIF_M1 = UIF_M1 - mean(UIF_M1)
    return(UIF_M1)
  }
  
  compute_UIF_M2 = function(DATA,x1valfix,rvalfix,mismode,distortval){
    X1 = DATA[,1] 
    Z = DATA[,2] 
    R = DATA[,3] 
    X2 = DATA[,4]  
    Y = DATA[,5]
    
    ############################################
    # M2 =  M[r | x1;]
    ############################################
    # x1valfix = 0 
    # rvalfix = 0 
    Ir = (R==rvalfix)*1
    Ix1 = (X1==x1valfix)*1
    
    prob.X1 = mean(X1)*X1 + (1-mean(X1))*(1-X1)
    prob.r.X1 = mapply(function(idx,rvalfix){
      pred.R.X1[idx,rvalfix+1]
    },c(1:nrow(DATA)),rep(rvalfix,nrow(DATA)))
    if (mismode == 1){
      prob.r.X1 = fix_pred(mis_pred(prob.r.X1,distortval))
    }
    prob.r.x1 = rep(predict(model.R.X1,newdata = data.matrix(data.frame(X1=x1valfix)),reshape=TRUE)[rvalfix+1],nrow(DATA))
    UIF_M2 = (Ix1/prob.X1)*(Ir-prob.r.X1)+prob.r.x1
    # EIF_M2 = UIF_M2 - mean(UIF_M2)
    return(UIF_M2)
  }
  
  compute_Yx = function(DATA,yvalfix,x1valfix,x2valfix,mismode,distortval){
    ############################################
    # UIF =  \sum_{r}UIF_M1*mean(UIF_M2) + EIF_M2 * mean(UIF_M1)
    # M1 =  M[y | (r,x2);(x1;z)]
    # M2 =  M[r | x1;]
    ############################################
    UIF = rep(0,nrow(DATA))
    for (rvalfix in Runique){
      UIF_M1 = compute_UIF_M1(DATA=DATA,yval=yvalfix,rvalfix=rvalfix,x2valfix=x2valfix,mismode,distortval)
      EIF_M1 = UIF_M1 - mean(UIF_M1,na.rm=T)
      
      UIF_M2 = compute_UIF_M2(DATA=DATA,x1valfix=x1valfix,rvalfix=rvalfix,mismode,distortval)
      EIF_M2 = UIF_M1 - mean(UIF_M1,na.rm=T)
      
      UIF = UIF + UIF_M1*mean(UIF_M2,na.rm=T) + EIF_M2*mean(UIF_M1,na.rm=T)
    }
    return(mean(UIF,na.rm=T))
  }
  
  Yx = rep(0,length(X1unique)*length(X2unique))
  idx = 1 
  for (x1valfix in X1unique){
    for (x2valfix in X2unique){
      Yx[idx] = compute_Yx(DATA=DATA,yvalfix=yvalfix,x1valfix=x1valfix,x2valfix=x2valfix,mismode,distortval)    
      idx = idx + 1 
    }
  }
  return(Yx)
}




